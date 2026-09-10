/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef REKS_GRADIENT_ENGINE_H
#define REKS_GRADIENT_ENGINE_H

/// @file reks_gradient_engine.h
/// @brief Assemble combined gradient g(N) and Hessian H(NxN) for TRAH solver
///
/// N = n_rot + sum_blocks block.count, where:
///   n_rot        = n_act*(n_act-1)/2
///   FON blocks   = one per (sector s, generation gen) with active geminals;
///                  block.count = sa_cassette.geminals_active(s, gen).size().
///
/// Components:
///   g = [g_orb; g_fon(block 0); g_fon(block 1); ...]
///   H block layout (block-diagonal across FON blocks):
///       [ H_oo     H_of(b0)  H_of(b1) ... ]
///       [ H_fo(b0) H_ff(b0)  0        ... ]
///       [ H_fo(b1) 0         H_ff(b1) ... ]
///
/// - g_orb:      orbital gradient from generalized Fock asymmetry
/// - g_fon(b):   block b's FON gradient via compute_weight_derivs[gen] accumulated
/// - H_oo:       orbital-orbital Hessian (diagonal regularized; exact 1e+2e for diagnostics)
/// - H_ff(b):    block b's FON x FON Hessian (diagonal via d2C, off-diag via mixed derivs)
/// - H_of(b):    orbital x block-b-FON cross-coupling

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_math.h"
#include "reks_reel.h"
#include "reks_si_types.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace reks {

/// Result of combined gradient+Hessian computation.
///
/// Layout of g (size N = n_rot + sum_blocks block.count):
///   [0 .. n_rot)                        kappa orbital rotations
///   [block.offset .. +block.count)      block b's (sector, generation) FON occupations
///
/// Same layout for H (row-major NxN). Cross-block FON x FON blocks are identically
/// zero: within a sector every config reads geminals of a single generation, and no
/// config reads another sector's geminals (block-diagonal Hessian).
struct CombinedGradient {
    std::vector<studio::FonBlock> fon_blocks;   ///< FON blocks, lexicographic (s, gen) order
    int N = 0;                                  ///< n_rot + sum_blocks block.count
    std::vector<double> g;      ///< Combined gradient (size N)
    std::vector<double> H;      ///< Combined Hessian (size N*N, row-major)
    /// Diagnostics (filled only when want_diagnostics), one slot per FON block.
    std::vector<double> hoo_evals_raw;  ///< H_oo eigenvalues
    std::vector<std::vector<double>> fon_hess_evals;  ///< [block] FON-Hessian eigenvalues
    std::vector<std::vector<double>> fon_hess_block;  ///< [block] FON-Hessian block (count^2)
};

/// Static methods assembling gradient and Hessian for TRAH from REKS quantities.
class REKSGradientEngine {
   public:
    /// Extract active-active orbital gradient from generalized Fock matrix.
    ///
    /// g_orb(i,j) = -2 * (F_gen[ai][aj] - F_gen[aj][ai])
    ///
    /// At convergence F_gen is symmetric (generalized Brillouin condition),
    /// so g_orb -> 0. The negative sign matches dE/d(kappa_{ij}).
    static std::vector<double> extract_orbital_gradient(SharedMatrix F_gen,
                                                        const std::vector<int>& active) {
        int n_act = static_cast<int>(active.size());
        int n_rot = n_act * (n_act - 1) / 2;
        std::vector<double> g_orb(n_rot);
        double** Fg = F_gen->pointer(0);

        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                int ai = active[i];
                int aj = active[j];
                g_orb[idx] = -2.0 * (Fg[ai][aj] - Fg[aj][ai]);
                ++idx;
            }
        }
        return g_orb;
    }

    /// View over the thread-local weight-derivative accumulators of one (s, gen)
    /// block. Row fi (stride-major) holds, over all microstates L,
    ///   dC[fi][L]  = sum_K w_K * dC_L^K  / dx_{g_fi}
    ///   d2C[fi][L] = sum_K w_K * d2C_L^K / dx_{g_fi}^2
    /// The next accumulate_weight_derivs call on this thread overwrites the storage;
    /// gen_ carries the sweep it was taken from and the row accessors throw if a newer
    /// sweep has run. [KKK]
    struct WeightDerivs {
        const double* dC  = nullptr;
        const double* d2C = nullptr;
        int n_fon  = 0;
        int stride = 0;
        std::uint64_t gen_ = 0;

        const double* dC_row(int fi)  const { check(); return dC  + static_cast<size_t>(fi) * stride; }
        const double* d2C_row(int fi) const { check(); return d2C + static_cast<size_t>(fi) * stride; }

       private:
        void check() const {
            if (gen_ != sweep_generation())
                throw std::runtime_error(
                    "REKSGradientEngine::WeightDerivs: read after a newer weight-derivative "
                    "sweep on the same thread");
        }
    };

    /// Counter bumped by every accumulate_weight_derivs call on this thread.
    static std::uint64_t& sweep_generation() {
        thread_local std::uint64_t g = 0;
        return g;
    }

    /// Weight-derivative sweep over run sector s's configs, one evaluation per
    /// active geminal of generation gen; Cassette::compute_weight_derivs fills
    /// dC and d2C together.
    static WeightDerivs accumulate_weight_derivs(int s, int gen,
                                                 const studio::Cassette& sa_cassette,
                                                 const studio::FONSnapshot& fons) {
        const auto& gem_active = sa_cassette.geminals_active(s, gen);
        const auto& positions  = sa_cassette.sector_config_positions(s);
        const int n_fon   = static_cast<int>(gem_active.size());
        const int stride  = sa_cassette.n_total_microstates();
        const int n_micro = sa_cassette.n_catalog_microstates();

        // dC_K/d2C_K are resized to the catalog count by compute_weight_derivs itself.
        thread_local std::vector<double> dC_all, d2C_all, dC_K, d2C_K;
        const size_t need = static_cast<size_t>(n_fon) * stride;
        if (dC_all.size()  < need) dC_all.resize(need);
        if (d2C_all.size() < need) d2C_all.resize(need);
        if (static_cast<int>(dC_K.size())  < n_micro) dC_K.resize(n_micro);
        if (static_cast<int>(d2C_K.size()) < n_micro) d2C_K.resize(n_micro);

        for (int fi = 0; fi < n_fon; ++fi) {
            double* dC  = dC_all.data()  + static_cast<size_t>(fi) * stride;
            double* d2C = d2C_all.data() + static_cast<size_t>(fi) * stride;
            for (int L : sa_cassette.microstates_active) { dC[L] = 0.0; d2C[L] = 0.0; }
            const int g = gem_active[fi];
            for (int pos : positions) {
                const int    K   = sa_cassette.K_indices[pos];
                const double w_K = sa_cassette.K_weights[pos];
                sa_cassette.compute_weight_derivs(gen, K, g, fons, dC_K, d2C_K);
                for (int L : sa_cassette.diag_deps(K).microstate_writes) {
                    dC[L]  += w_K * dC_K[L];
                    d2C[L] += w_K * d2C_K[L];
                }
            }
        }
        return WeightDerivs{dC_all.data(), d2C_all.data(), n_fon, stride, ++sweep_generation()};
    }

    /// Compute orbital x block-FON cross-Hessian H_of from a WeightDerivs sweep.
    ///
    /// H_{(p,q), i} = -2 * sum_L wd.dC_row(i)[L] * G_L(p,q)
    ///   where G_L(p,q) = n_alpha_Lp * Fa_L[p][q] + n_beta_Lp * Fb_L[p][q]
    ///                  - n_alpha_Lq * Fa_L[q][p] - n_beta_Lq * Fb_L[q][p]
    ///
    /// @return H_of as flat vector (n_rot * n_fon), row-major:
    ///         H_of[rot_idx * n_fon + fon_idx]
    /// F_MO_arows_*_L hold the active rows of the per-L MO Fock (row a_idx = MO
    /// active[a_idx]); only active x active elements are read here.
    static std::vector<double> cross_hessian_from(const WeightDerivs& wd,
                                                  const std::vector<SharedMatrix>& F_MO_arows_a_L,
                                                  const std::vector<SharedMatrix>& F_MO_arows_b_L,
                                                  const studio::Cassette& sa_cassette,
                                                  const std::vector<int>& active) {
        const int n_act   = static_cast<int>(active.size());
        const int n_rot   = n_act * (n_act - 1) / 2;
        const int n_fon   = wd.n_fon;

        std::vector<double> H_of(n_rot * n_fon, 0.0);
        for (int fi = 0; fi < n_fon; ++fi) {
            const double* dC_acc = wd.dC_row(fi);

            for (int L : sa_cassette.microstates_active) {
                if (std::abs(dC_acc[L]) < 1e-14) continue;
                double** Farows_a = F_MO_arows_a_L[L]->pointer(0);
                double** Farows_b = F_MO_arows_b_L[L]->pointer(0);
                const auto& m     = sa_cassette.microstate(L);
                const auto& alpha = m.alpha;
                const auto& beta  = m.beta;

                int idx = 0;
                for (int i = 0; i < n_act; ++i) {
                    for (int j = i + 1; j < n_act; ++j) {
                        const int ai = active[i];
                        const int aj = active[j];
                        const int na_i = alpha[i];
                        const int nb_i = beta[i];
                        const int na_j = alpha[j];
                        const int nb_j = beta[j];

                        const double G_L = na_i * Farows_a[i][aj] + nb_i * Farows_b[i][aj]
                                         - na_j * Farows_a[j][ai] - nb_j * Farows_b[j][ai];

                        H_of[idx * n_fon + fi] += -2.0 * dC_acc[L] * G_L;
                        ++idx;
                    }
                }
            }
        }
        return H_of;
    }

    /// Diagonal of the frozen-Fock 1e active-active orbital Hessian.
    ///
    /// H_oo[(i,j),(i,j)] ~ -2 * sum_L C_L *
    ///     [ (n_alpha_Li - n_alpha_Lj) * (Fa_L[ai][ai] - Fa_L[aj][aj])
    ///     + (n_beta_Li  - n_beta_Lj)  * (Fb_L[ai][ai] - Fb_L[aj][aj]) ]
    /// Sign and terms are the r1 == r2 case of compute_orbital_hessian_exact's A1e.
    static std::vector<double> compute_orbital_hessian_diag(const std::vector<SharedMatrix>& F_MO_arows_a_L,
                                                            const std::vector<SharedMatrix>& F_MO_arows_b_L,
                                                            const std::vector<double>& C_L,
                                                            const studio::Cassette& sa_cassette,
                                                            const std::vector<int>& active) {
        const int n_act = static_cast<int>(active.size());
        const int n_rot = n_act * (n_act - 1) / 2;

        std::vector<double> H_diag(n_rot, 0.0);

        for (int L : sa_cassette.microstates()) {
            if (std::abs(C_L[L]) < 1e-14) continue;

            double** Farows_a = F_MO_arows_a_L[L]->pointer(0);
            double** Farows_b = F_MO_arows_b_L[L]->pointer(0);
            const auto& m     = sa_cassette.microstate(L);
            const auto& alpha = m.alpha;
            const auto& beta  = m.beta;

            int idx = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    const int ai = active[i];
                    const int aj = active[j];
                    const int na_i = alpha[i];
                    const int nb_i = beta[i];
                    const int na_j = alpha[j];
                    const int nb_j = beta[j];

                    const double d_na = na_i - na_j;
                    const double d_nb = nb_i - nb_j;
                    const double d_ea = Farows_a[i][ai] - Farows_a[j][aj];
                    const double d_eb = Farows_b[i][ai] - Farows_b[j][aj];

                    H_diag[idx] += -2.0 * C_L[L] * (d_na * d_ea + d_nb * d_eb);
                    ++idx;
                }
            }
        }
        return H_diag;
    }

    /// Active-active orbital Hessian at fixed FON.
    ///
    /// H_oo[(i,j),(p,q)] = d(g_orb_ij) / d(kappa_pq) = -2*sym(A1e) + 4*sym(A2e),
    ///   sym(M)_{r1,r2} = 1/2 (M_{r1,r2} + M_{r2,r1})
    ///
    /// A1e_{(ij),(pq)} = sum_L C_L sum_sigma (n_sigma_iL - n_sigma_jL)
    ///   [ d_iq F_L^sigma[p,j] - d_ip F_L^sigma[q,j]
    ///   + d_jq F_L^sigma[i,p] - d_jp F_L^sigma[i,q] ]                  (F_L frozen)
    /// A2e_{(ij),(pq)} = sum_L C_L sum_{sigma,sigma'} (n_sigma_iL - n_sigma_jL)
    ///   (n_sigma'_pL - n_sigma'_qL)
    ///   * 1/2 [ (ij|qp) - d_{sigma sigma'} x_alpha (ip|qj) + (ij|pq) - d_{sigma sigma'} x_alpha (iq|pj) ]
    ///
    /// (ab|cd) is the chemist active-MO ERI eri.at(a,b,c,d); F_L^sigma[p,j] reads
    /// F_MO_arows_sigma_L[L] row p (active order), column active[j] (absolute MO).
    /// The 2e response is exact because the REKS 2-RDM is determinantal: the active
    /// Fock response under active rotation collapses to occupation-difference products.
    /// A1e/A2e are built non-symmetric, then symmetrized (the energy Hessian needs the
    /// symmetric part). x_alpha is the hybrid exact-exchange coefficient (1 for HF): the
    /// microstate Fock scales K by it, so the A2e exchange kernel does too; Coulomb is
    /// unscaled. The DFT XC-kernel density response (f_xc) is omitted -- the Hessian is
    /// approximate for DFT, the gradient stays exact.
    ///
    /// @param eri Active-MO ERI tiles; must cover every canonical pair (k<=l) over the
    ///            active orbitals, n_act*(n_act+1)/2 tiles.
    /// @return H_oo as flat n_rot*n_rot row-major (symmetric).
    static std::vector<double> compute_orbital_hessian_exact(
        const std::vector<SharedMatrix>& F_MO_arows_a_L,
        const std::vector<SharedMatrix>& F_MO_arows_b_L,
        const std::vector<double>& C_L,
        const ActiveEriTiles& eri,
        const studio::Cassette& sa_cassette,
        const std::vector<int>& active,
        double x_alpha) {

        const int n_act = static_cast<int>(active.size());
        const int n_rot = n_act * (n_act - 1) / 2;

        std::vector<std::pair<int, int>> rot_pairs(n_rot);
        {
            int idx = 0;
            for (int i = 0; i < n_act; ++i)
                for (int j = i + 1; j < n_act; ++j)
                    rot_pairs[idx++] = {i, j};
        }

        std::vector<double> A1e(n_rot * n_rot, 0.0);
        std::vector<double> A2e(n_rot * n_rot, 0.0);

        // Transcription of the A1e/A2e sums above. Outer: microstate L, weight CL=C_L[L],
        // occupations alpha[.]/beta[.] = n_{a,b}_L. r1=(i,j) and r2=(p,q) are the coupled
        // rotation pairs -> entry (r1,r2). Occupation differences become dn{a,b}_ij =
        // n_iL-n_jL and dn{a,b}_pq = n_pL-n_qL; the A1e Kronecker deltas d_ab become the
        // (i==q)/(i==p)/(j==q)/(j==p) guards; the A2e d_{sigma sigma'} is the same-spin vs
        // opposite-spin split (coul-exch vs coul). A1e/A2e accumulate per (r1,r2); the
        // final loop forms H_oo = -2*sym(A1e) + 4*sym(A2e).
        for (int L : sa_cassette.microstates()) {
            const double CL = C_L[L];
            if (std::abs(CL) < 1e-14) continue;

            double** Fa = F_MO_arows_a_L[L]->pointer(0);
            double** Fb = F_MO_arows_b_L[L]->pointer(0);
            const auto& alpha = sa_cassette.microstate(L).alpha;
            const auto& beta  = sa_cassette.microstate(L).beta;

            for (int r1 = 0; r1 < n_rot; ++r1) {
                const int i = rot_pairs[r1].first;
                const int j = rot_pairs[r1].second;
                const int dna_ij = static_cast<int>(alpha[i]) - static_cast<int>(alpha[j]);
                const int dnb_ij = static_cast<int>(beta[i])  - static_cast<int>(beta[j]);
                if (dna_ij == 0 && dnb_ij == 0) continue;

                for (int r2 = 0; r2 < n_rot; ++r2) {
                    const int p = rot_pairs[r2].first;
                    const int q = rot_pairs[r2].second;

                    // A1e: frozen-Fock 1e response (index rotation of g_orb).
                    double a1 = 0.0;
                    if (dna_ij != 0) {
                        double f = 0.0;
                        if (i == q) f += Fa[p][active[j]];
                        if (i == p) f -= Fa[q][active[j]];
                        if (j == q) f += Fa[i][active[p]];
                        if (j == p) f -= Fa[i][active[q]];
                        a1 += dna_ij * f;
                    }
                    if (dnb_ij != 0) {
                        double f = 0.0;
                        if (i == q) f += Fb[p][active[j]];
                        if (i == p) f -= Fb[q][active[j]];
                        if (j == q) f += Fb[i][active[p]];
                        if (j == p) f -= Fb[i][active[q]];
                        a1 += dnb_ij * f;
                    }
                    A1e[r1 * n_rot + r2] += CL * a1;

                    // A2e: 2e Hxc density response (determinantal 2-RDM).
                    const int dna_pq = static_cast<int>(alpha[p]) - static_cast<int>(alpha[q]);
                    const int dnb_pq = static_cast<int>(beta[p])  - static_cast<int>(beta[q]);
                    if ((dna_pq != 0) || (dnb_pq != 0)) {
                        // coul = 1/2[(ij|qp)+(ij|pq)]; exch = 1/2 x_alpha[(ip|qj)+(iq|pj)]
                        const double coul = 0.5 * (eri.at(i, j, q, p) + eri.at(i, j, p, q));
                        const double exch = 0.5 * x_alpha * (eri.at(i, p, q, j) + eri.at(i, q, p, j));
                        double a2 = 0.0;
                        a2 += dna_ij * dna_pq * (coul - exch);  // (alpha, alpha)
                        a2 += dna_ij * dnb_pq * coul;           // (alpha, beta)
                        a2 += dnb_ij * dna_pq * coul;           // (beta, alpha)
                        a2 += dnb_ij * dnb_pq * (coul - exch);  // (beta, beta)
                        A2e[r1 * n_rot + r2] += CL * a2;
                    }
                }
            }
        }

        std::vector<double> H_oo(n_rot * n_rot, 0.0);
        for (int r1 = 0; r1 < n_rot; ++r1) {
            for (int r2 = 0; r2 < n_rot; ++r2) {
                const double a1_sym = 0.5 * (A1e[r1 * n_rot + r2] + A1e[r2 * n_rot + r1]);
                const double a2_sym = 0.5 * (A2e[r1 * n_rot + r2] + A2e[r2 * n_rot + r1]);
                H_oo[r1 * n_rot + r2] = -2.0 * a1_sym + 4.0 * a2_sym;
            }
        }
        return H_oo;
    }

    /// g_fon[i] = sum_L wd.dC_row(i)[L] * E_L[L]
    static std::vector<double> fon_gradient_from(const WeightDerivs& wd,
                                                 const studio::Cassette& sa_cassette,
                                                 const std::vector<double>& E_L) {
        std::vector<double> g_fon(wd.n_fon, 0.0);
        for (int fi = 0; fi < wd.n_fon; ++fi) {
            const double* dC_acc = wd.dC_row(fi);
            double acc = 0.0;
            for (int L : sa_cassette.microstates_active)
                acc += dC_acc[L] * E_L[L];
            g_fon[fi] = acc;
        }
        return g_fon;
    }

    static std::vector<double> compute_fon_gradient(int s, int gen,
                                                    const studio::Cassette& sa_cassette,
                                                    const studio::FONSnapshot& fons,
                                                    const std::vector<double>& E_L) {
        return fon_gradient_from(accumulate_weight_derivs(s, gen, sa_cassette, fons),
                                 sa_cassette, E_L);
    }

    /// H_ff for block (s, gen), n_fon x n_fon row-major.
    /// H[i,i]        = sum_L wd.d2C_row(i)[L] * E_L[L]
    /// H[i,j] (i!=j) = sum_K w_K * sum_L (d2C^K / dx_g_i dx_g_j) * E_L[L]
    ///   via sa_cassette.compute_mixed_weight_derivs(gen, ...). K runs over run
    ///   sector s's configs only.
    static std::vector<double> fon_hessian_from(const WeightDerivs& wd, int s, int gen,
                                                const studio::Cassette& sa_cassette,
                                                const studio::FONSnapshot& fons,
                                                const std::vector<double>& E_L) {
        const auto& gem_active = sa_cassette.geminals_active(s, gen);
        const auto& positions  = sa_cassette.sector_config_positions(s);
        const int n_fon   = wd.n_fon;
        const int n_micro = static_cast<int>(E_L.size());

        thread_local std::vector<double> d2C_mix;
        if (static_cast<int>(d2C_mix.size()) < n_micro) d2C_mix.resize(n_micro);

        std::vector<double> H_ff(n_fon * n_fon, 0.0);
        for (int fi = 0; fi < n_fon; ++fi) {
            const double* d2C_acc = wd.d2C_row(fi);
            double acc = 0.0;
            for (int L : sa_cassette.microstates_active)
                acc += d2C_acc[L] * E_L[L];
            H_ff[fi * n_fon + fi] = acc;
        }

        for (int fi = 0; fi < n_fon; ++fi) {
            const int g_i = gem_active[fi];
            for (int fj = fi + 1; fj < n_fon; ++fj) {
                const int g_j = gem_active[fj];
                double acc = 0.0;
                for (int pos : positions) {
                    const int    K   = sa_cassette.K_indices[pos];
                    const double w_K = sa_cassette.K_weights[pos];
                    sa_cassette.compute_mixed_weight_derivs(gen, K, g_i, g_j, fons, d2C_mix);
                    for (int L : sa_cassette.diag_deps(K).microstate_writes)
                        acc += w_K * d2C_mix[L] * E_L[L];
                }
                H_ff[fi * n_fon + fj] = acc;
                H_ff[fj * n_fon + fi] = acc;
            }
        }
        return H_ff;
    }

    static std::vector<double> compute_fon_hessian(int s, int gen,
                                                   const studio::Cassette& sa_cassette,
                                                   const studio::FONSnapshot& fons,
                                                   const std::vector<double>& E_L) {
        return fon_hessian_from(accumulate_weight_derivs(s, gen, sa_cassette, fons),
                                s, gen, sa_cassette, fons, E_L);
    }

    /// Assemble the combined gradient and Hessian (layout: CombinedGradient).
    ///
    /// @param F_gen             Generalized Fock matrix
    /// @param reel              Reel (provides per-sector FON snapshots)
    /// @param sa_cassette       SA cassette (per-sector config sublists, active subsets)
    /// @param F_MO_arows_a_L   Per-microstate alpha MO Fock active rows
    /// @param F_MO_arows_b_L   Per-microstate beta MO Fock active rows
    /// @param E_L           Microstate energies
    /// @param C_L           Microstate weights
    /// @param active            Active MO indices (absolute, 0-based)
    /// @param want_diagnostics  Fill the raw Hessian arrays (hoo_evals_raw and the
    ///                          per-block fon_hess_evals/fon_hess_block).
    /// @param active_eri        Active-MO ERI tiles over every canonical pair (k<=l),
    ///                          needed for the diagnostic H_oo 2e response; null or
    ///                          empty leaves hoo_evals_raw empty.
    /// @param x_alpha           Hybrid exact-exchange coefficient.
    static CombinedGradient compute(SharedMatrix                       F_gen,
                                    studio::Reel&                        reel,
                                    const studio::Cassette&              sa_cassette,
                                    const std::vector<SharedMatrix>&   F_MO_arows_a_L,
                                    const std::vector<SharedMatrix>&   F_MO_arows_b_L,
                                    const std::vector<double>&         E_L,
                                    const std::vector<double>&         C_L,
                                    const std::vector<int>&            active,
                                    bool                               want_diagnostics,
                                    const ActiveEriTiles*              active_eri,
                                    double                             x_alpha) {
        const int n_act  = static_cast<int>(active.size());
        const int n_rot  = n_act * (n_act - 1) / 2;

        // FON blocks in lexicographic (s, gen) order, skipping empty ones; offsets
        // accumulate after the n_rot orbital rotations.
        CombinedGradient result;
        int acc_off = n_rot;
        for (int s = 0; s < sa_cassette.n_sectors(); ++s) {
            const int n_gen = sa_cassette.n_generations();
            for (int gen = 0; gen < n_gen; ++gen) {
                const int cnt =
                    static_cast<int>(sa_cassette.geminals_active(s, gen).size());
                if (cnt == 0) continue;
                result.fon_blocks.push_back(studio::FonBlock{s, gen, acc_off, cnt});
                acc_off += cnt;
            }
        }
        const int N = acc_off;
        result.N = N;
        result.g.assign(N, 0.0);
        result.H.assign(N * N, 0.0);
        result.fon_hess_evals.assign(result.fon_blocks.size(), {});
        result.fon_hess_block.assign(result.fon_blocks.size(), {});

        auto g_orb = extract_orbital_gradient(F_gen, active);
        for (int i = 0; i < n_rot; ++i) result.g[i] = g_orb[i];

        // H_oo diagonal: max(|h|, H_oo_min) over the frozen-Fock 1e diagonal. |h|
        // makes the solver diagonal positive whatever the curvature sign; the floor
        // bounds its magnitude below.
        auto H_oo_diag = compute_orbital_hessian_diag(F_MO_arows_a_L, F_MO_arows_b_L,
                                                      C_L, sa_cassette, active);
        constexpr double H_oo_min = 0.1;
        for (int i = 0; i < n_rot; ++i) {
            const double h = H_oo_diag[i];
            result.H[i * N + i] = std::max(std::abs(h), H_oo_min);
        }
        // Spectrum of the unregularized orbital Hessian (1e + 2e) for diagnostics.
        if (want_diagnostics && n_rot > 0 && active_eri && !active_eri->empty()) {
            auto H_oo_exact = compute_orbital_hessian_exact(F_MO_arows_a_L, F_MO_arows_b_L,
                                                            C_L, *active_eri, sa_cassette,
                                                            active, x_alpha);
            std::vector<double> evals, evecs;
            lapack_diagonalize(H_oo_exact.data(), n_rot, evals, evecs);
            result.hoo_evals_raw = std::move(evals);
        }

        // Per-block FON contributions: g_fon(b), H_ff(b) diagonal block, H_of(b)
        // orbital cross-coupling. Each block reads its own sector's FON snapshot.
        for (size_t bi = 0; bi < result.fon_blocks.size(); ++bi) {
            const studio::FonBlock& blk = result.fon_blocks[bi];
            const int nf  = blk.count;
            if (nf == 0) continue;
            const int off = blk.offset;
            const int gen = blk.gen;
            const studio::FONSnapshot& fons = reel.fon_state[blk.s];

            const WeightDerivs wd = accumulate_weight_derivs(blk.s, gen, sa_cassette, fons);

            auto g_fon = fon_gradient_from(wd, sa_cassette, E_L);
            for (int k = 0; k < nf; ++k) result.g[off + k] = g_fon[k];

            auto H_ff = fon_hessian_from(wd, blk.s, gen, sa_cassette, fons, E_L);
            for (int i = 0; i < nf; ++i)
                for (int j = 0; j < nf; ++j)
                    result.H[(off + i) * N + (off + j)] = H_ff[i * nf + j];

            auto H_of = cross_hessian_from(wd, F_MO_arows_a_L, F_MO_arows_b_L,
                                           sa_cassette, active);
            // Mirrored into the (off+k, i) transpose entry to keep H symmetric.
            for (int i = 0; i < n_rot; ++i)
                for (int k = 0; k < nf; ++k) {
                    const double v = H_of[i * nf + k];
                    result.H[i * N + (off + k)] = v;
                    result.H[(off + k) * N + i] = v;
                }

            // Diagnostics: per-block FON-Hessian spectrum + block.
            if (want_diagnostics) {
                std::vector<double> evals, evecs;
                lapack_diagonalize(H_ff.data(), nf, evals, evecs);  // leaves H_ff intact
                result.fon_hess_evals[bi] = std::move(evals);
                result.fon_hess_block[bi] = std::move(H_ff);
            }
        }

        return result;
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_GRADIENT_ENGINE_H
