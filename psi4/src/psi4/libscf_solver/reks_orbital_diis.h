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

#ifndef REKS_ORBITAL_DIIS_H
#define REKS_ORBITAL_DIIS_H

#include "reks_report_level.h"
#include "reks_cassette.h"
#include "reks_diis_config.h"
#include "reks_diis_core.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

namespace psi {
namespace reks {

/// @class OrbitalDIISEngine
/// @brief Orbital-rotation DIIS for REKS SCF (Ionova, I. V.; Carter, E. A. J. Chem. Phys. 1995,
///        102, 1251; r-GDIIS: Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472).
///
/// Extrapolates non-redundant skew-symmetric rotation generators kappa. The
/// error space has no virt-virt or core-core blocks (gauge-free).
///
/// C_new = C_old * exp(K), K antisym, K[i][j] = +kappa for i<j.
/// H_diag_pq is the diagonal one-electron orbital Hessian element.
///
/// Index layout of kappa vector (store and apply MUST agree on this layout):
///   block 0: core-active        size Ncore * n_act
///            idx = c*n_act + k                       c<Ncore, k<n_act
///   block 1: core-virt          size Ncore * n_virt
///            idx = off1 + c*n_virt + (v - first_virt)
///   block 2: active-virt        size n_act * n_virt
///            idx = off2 + k*n_virt + (v - first_virt)
///   block 3: active-active      size n_act*(n_act-1)/2 (i<j only)
///            idx = off3 + canonical_pair_index(i, j)
///
///   first_virt = active_mo_indices.back() + 1
///
/// Newton step:
///   A_pq   = F_gen[q][p] - F_gen[p][q]
///   B_pq   = -H_diag_pq/2 + gamma_pq + ipr_pq   // active-active
///   B_pq   = -H_diag_pq/2 + gamma_pq            // core-active
///   B_pq   = -H_diag_pq/2                       // core-virt and active-virt
///   denom  = -max(|B_pq|, B_floor)
///   e_pq   = A_pq / denom                       // Newton step in kappa
class OrbitalDIISEngine {
   public:
    OrbitalDIISEngine() : OrbitalDIISEngine(DiisConfig{}) {}
    explicit OrbitalDIISEngine(DiisConfig cfg)
        : config_(cfg),
          core_(cfg.max_vectors),
          chain_(assemble_conditioning_chain(cfg)) {}

    static int n_nonred(int Ncore, int n_act, int n_virt) {
        return Ncore * n_act + Ncore * n_virt + n_act * n_virt + n_act * (n_act - 1) / 2;
    }

    /// Start index of block b in the packed kappa vector.
    static int offset_block(int b, int Ncore, int n_act, int n_virt) {
        int o0 = 0;
        int o1 = o0 + Ncore * n_act;
        int o2 = o1 + Ncore * n_virt;
        int o3 = o2 + n_act * n_virt;
        switch (b) {
            case 0: return o0;
            case 1: return o1;
            case 2: return o2;
            case 3: return o3;
            default: return -1;
        }
    }

    /// K[q][p] = -K[p][q] over the strict upper triangle; tiled traversal keeps the
    /// transposed write within one tile's cache lines.
    static void mirror_skew_lower(double** K, int nmo) {
        constexpr int kTile = 64;
        for (int i0 = 0; i0 < nmo; i0 += kTile) {
            const int i1 = std::min(i0 + kTile, nmo);
            for (int j0 = i0; j0 < nmo; j0 += kTile) {
                const int j1 = std::min(j0 + kTile, nmo);
                for (int i = i0; i < i1; ++i)
                    for (int j = std::max(j0, i + 1); j < j1; ++j) K[j][i] = -K[i][j];
            }
        }
    }

    /// Build skew-symmetric K (nmo x nmo) from packed kappa vector.
    /// Inverse of the packing produced by compute_newton_rotation().
    static void unpack_kappa_to_K(const std::vector<double>& kappa,
                                  SharedMatrix K_out,
                                  const std::vector<int>& active_mo_indices,
                                  int Ncore, int nmo) {
        K_out->zero();
        double** K = K_out->pointer(0);
        const int n_act = static_cast<int>(active_mo_indices.size());
        const int first_virt = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;
        const int n_virt = nmo - first_virt;

        const int o0 = offset_block(0, Ncore, n_act, n_virt);
        const int o1 = offset_block(1, Ncore, n_act, n_virt);
        const int o2 = offset_block(2, Ncore, n_act, n_virt);
        const int o3 = offset_block(3, Ncore, n_act, n_virt);

        // block 0: core-active. Pair (c, ai) with c<ai by construction.
        for (int c = 0; c < Ncore; ++c) {
            for (int k = 0; k < n_act; ++k) {
                K[c][active_mo_indices[k]] = kappa[o0 + c * n_act + k];
            }
        }

        // block 1: core-virt
        for (int c = 0; c < Ncore; ++c) {
            for (int v = first_virt; v < nmo; ++v)
                K[c][v] = kappa[o1 + c * n_virt + (v - first_virt)];
        }

        // block 2: active-virt
        for (int k = 0; k < n_act; ++k) {
            int ai = active_mo_indices[k];
            for (int v = first_virt; v < nmo; ++v)
                K[ai][v] = kappa[o2 + k * n_virt + (v - first_virt)];
        }

        // block 3: active-active i<j pairs
        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            int ai = active_mo_indices[i];
            for (int j = i + 1; j < n_act; ++j) {
                K[ai][active_mo_indices[j]] = kappa[o3 + idx];
                ++idx;
            }
        }

        mirror_skew_lower(K, nmo);
    }

    /// Read the packed kappa vector back out of a skew-symmetric K (nmo x nmo).
    /// Inverse of unpack_kappa_to_K() on the non-redundant blocks; redundant entries of K
    /// are discarded.
    static void pack_K_to_kappa(SharedMatrix K_in, std::vector<double>& kappa_out,
                                const std::vector<int>& active_mo_indices,
                                int Ncore, int nmo) {
        double** K = K_in->pointer(0);
        const int n_act = static_cast<int>(active_mo_indices.size());
        const int first_virt = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;
        const int n_virt = nmo - first_virt;
        kappa_out.assign(n_nonred(Ncore, n_act, n_virt), 0.0);

        const int o0 = offset_block(0, Ncore, n_act, n_virt);
        const int o1 = offset_block(1, Ncore, n_act, n_virt);
        const int o2 = offset_block(2, Ncore, n_act, n_virt);
        const int o3 = offset_block(3, Ncore, n_act, n_virt);

        // block 0: core-active
        for (int c = 0; c < Ncore; ++c)
            for (int k = 0; k < n_act; ++k)
                kappa_out[o0 + c * n_act + k] = K[c][active_mo_indices[k]];

        // block 1: core-virt
        for (int c = 0; c < Ncore; ++c)
            for (int v = first_virt; v < nmo; ++v)
                kappa_out[o1 + c * n_virt + (v - first_virt)] = K[c][v];

        // block 2: active-virt
        for (int k = 0; k < n_act; ++k)
            for (int v = first_virt; v < nmo; ++v)
                kappa_out[o2 + k * n_virt + (v - first_virt)] = K[active_mo_indices[k]][v];

        // block 3: active-active i<j pairs
        int idx = 0;
        for (int i = 0; i < n_act; ++i)
            for (int j = i + 1; j < n_act; ++j)
                kappa_out[o3 + idx++] = K[active_mo_indices[i]][active_mo_indices[j]];
    }

    /// Compute Newton-step rotation in kappa basis (size = n_nonred).
    /// Layout matches unpack_kappa_to_K().
    ///
    /// @param F_gen           Asymmetric generalized Fock (nmo x nmo, MO basis)
    /// @param F_MO_diag_a_L  Per-microstate alpha MO Fock diagonal, [global L][nmo]
    /// @param F_MO_diag_b_L  Per-microstate beta  MO Fock diagonal, [global L][nmo]
    /// @param C_L         SA microstate weights, indexed by global L
    /// @param sa_cassette  Active pool (microstate occupations, iterates microstates_active)
    /// @param gamma_aa        Two-electron Hessian correction, active-active pairs
    /// @param gamma_ca        Two-electron Hessian correction, core-active pairs
    /// @param iteration       Iter number for tagging output
    /// @param narrate         False on probe calls, whose internals stay out of the report
    /// @param ipr_diag_aa     IPR Hessian diagonal, active-active pairs; length n_rot,
    ///                        zero-filled when IPR is off.
    /// @param g_packed        Non-null: receives the packed orbital gradient 2*A_pq.
    /// @param denom_packed    Non-null: receives the packed shaped denominator.
    ///
    /// g_packed non-null with denom_packed null selects gradient-only mode: only
    /// 2*A_pq is packed into g_packed; curvature, B, and the step are skipped and
    /// the returned kappa vector is empty.
    std::vector<double> compute_newton_rotation(
        SharedMatrix F_gen,
        const std::vector<std::vector<double>>& F_MO_diag_a_L,
        const std::vector<std::vector<double>>& F_MO_diag_b_L,
        const std::vector<double>& C_L,
        const studio::Cassette& sa_cassette,
        const std::vector<int>& active_mo_indices,
        int Ncore, int nmo,
        const std::vector<double>& gamma_aa = {},
        const std::vector<double>& gamma_ca = {},
        int iteration = -1,
        const std::vector<double>& ipr_diag_aa = {},
        std::vector<double>* g_packed = nullptr,
        std::vector<double>* denom_packed = nullptr,
        bool narrate = true) const {

        const int n_act = static_cast<int>(active_mo_indices.size());
        const int first_virt = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;
        const int n_virt = nmo - first_virt;

        const bool gradient_only = (g_packed != nullptr) && (denom_packed == nullptr);

        const int n_total = n_nonred(Ncore, n_act, n_virt);
        std::vector<double> kappa;
        if (!gradient_only) kappa.assign(n_total, 0.0);
        if (g_packed) g_packed->assign(n_total, 0.0);
        if (denom_packed) denom_packed->assign(n_total, 0.0);

        double** Fgen = F_gen->pointer(0);

        Composite1eProvider curvature;
        if (!gradient_only) {
            CurvatureInputs cin;
            cin.F_MO_diag_a_L = &F_MO_diag_a_L;
            cin.F_MO_diag_b_L = &F_MO_diag_b_L;
            cin.C_L = &C_L;
            cin.sa_cassette = &sa_cassette;
            cin.active_mo_indices = &active_mo_indices;
            cin.Ncore = Ncore;
            cin.nmo = nmo;
            cin.n_act = n_act;
            cin.first_virt = first_virt;
            cin.n_virt = n_virt;
            curvature.prepare(cin);
        }

        BAssemblyInputs corr;
        corr.gamma_aa = &gamma_aa;
        corr.gamma_ca = &gamma_ca;
        corr.ipr_diag_aa = &ipr_diag_aa;

        const DenominatorPolicy denom_policy(config_);
        const bool sat = config_.two_level_saturation;

        // Raw Newton step x = A/denom; two_level_saturation bounds it with the eigenvector
        // angle of a 2-level rotation:
        //   small x: 0.5*atan(2x) ~= x - 4x^3/3 + ...   (Newton recovers)
        //   large x: 0.5*atan(2x) -> sign(x) * pi/4      (saturates)
        // idx is this rotation's packed slot in g_packed/denom_packed.
        auto step_value = [&](double A, double B_total, int idx) -> double {
            const double denom = denom_policy.shape(B_total);
            if (g_packed) (*g_packed)[idx] = 2.0 * A;
            if (denom_packed) (*denom_packed)[idx] = denom;
            double raw = raw_step(A, denom);
            return sat ? 0.5 * std::atan(2.0 * raw) : raw;
        };

        if (narrate && report::reports(5)) {
            outfile->Printf("  [Kappa-DBG] iter=%d per-rotation Newton steps:\n", iteration);
            outfile->Printf("  %-4s %5s %5s  %12s  %12s  %12s  %12s  %12s  %12s  %s\n", "blk", "p",
                            "q", "A", "B_1e", "gamma+ipr", "B_total", "denom", "kappa", "sat?");
        }
        auto dbg_row = [&](const char* blk, int p, int q, double A, double H, double B_total,
                           double kap) {
            if (!narrate || !report::reports(5)) return;
            const double B_1e = -H / 2.0;
            const double denom = denom_policy.shape(B_total);
            outfile->Printf("  %-4s %5d %5d  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %s\n",
                            blk, p, q, A, B_1e, B_total - B_1e, B_total, denom, kap,
                            (sat && std::abs(kap) < std::abs(raw_step(A, denom))) ? "YES" : "no");
        };

        // Block 0: core-active
        const int o0 = offset_block(0, Ncore, n_act, n_virt);

        for (int c = 0; c < Ncore; ++c) {
            for (int k = 0; k < n_act; ++k) {
                int a = active_mo_indices[k];
                int idx = c * n_act + k;
                double A = Fgen[a][c] - Fgen[c][a];
                if (gradient_only) {
                    (*g_packed)[o0 + idx] = 2.0 * A;
                    continue;
                }
                double B_total = assemble_B(Block::CA, idx, curvature.ca(idx), corr);
                kappa[o0 + idx] = step_value(A, B_total, o0 + idx);
                dbg_row("ca", c, a, A, curvature.ca(idx), B_total, kappa[o0 + idx]);
            }
        }

        // Blocks 1 and 2: core-virt and active-virt
        const int o1 = offset_block(1, Ncore, n_act, n_virt);
        const int o2 = offset_block(2, Ncore, n_act, n_virt);

        // Block 1: core-virt
        {
            int idx = 0;
            for (int c = 0; c < Ncore; ++c) {
                for (int v = first_virt; v < nmo; ++v) {
                    // F_gen[v][c] = 0 (virtual rows zero) -> A = -F_gen[c][v]
                    double A = -Fgen[c][v];
                    if (gradient_only) {
                        (*g_packed)[o1 + idx] = 2.0 * A;
                        ++idx;
                        continue;
                    }
                    double B_total = assemble_B(Block::CV, idx, curvature.cv(idx), corr);
                    kappa[o1 + idx] = step_value(A, B_total, o1 + idx);
                    dbg_row("cv", c, v, A, curvature.cv(idx), B_total, kappa[o1 + idx]);
                    ++idx;
                }
            }
        }

        // Block 2: active-virt
        {
            int idx = 0;
            for (int k = 0; k < n_act; ++k) {
                int a = active_mo_indices[k];
                for (int v = first_virt; v < nmo; ++v) {
                    // F_gen[v][a] = 0 (virtual rows zero) -> A = -F_gen[a][v]
                    double A = -Fgen[a][v];
                    if (gradient_only) {
                        (*g_packed)[o2 + idx] = 2.0 * A;
                        ++idx;
                        continue;
                    }
                    double B_total = assemble_B(Block::AV, idx, curvature.av(idx), corr);
                    kappa[o2 + idx] = step_value(A, B_total, o2 + idx);
                    dbg_row("av", a, v, A, curvature.av(idx), B_total, kappa[o2 + idx]);
                    ++idx;
                }
            }
        }

        // Block 3: active-active i<j
        const int o3 = offset_block(3, Ncore, n_act, n_virt);
        const int n_rot_aa = n_act * (n_act - 1) / 2;

        {
            int idx = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    int ai = active_mo_indices[i];
                    int aj = active_mo_indices[j];
                    double A = Fgen[aj][ai] - Fgen[ai][aj];
                    if (gradient_only) {
                        (*g_packed)[o3 + idx] = 2.0 * A;
                        ++idx;
                        continue;
                    }
                    double B_total = assemble_B(Block::AA, idx, curvature.aa(idx), corr);
                    kappa[o3 + idx] = step_value(A, B_total, o3 + idx);
                    dbg_row("aa", ai, aj, A, curvature.aa(idx), B_total, kappa[o3 + idx]);
                    ++idx;
                }
            }
        }

        if (gradient_only) return {};

        if (narrate && report::reports(4)) {
            auto block_norm = [&](int off, int len) {
                double s = 0.0;
                for (int i = 0; i < len; ++i) s += kappa[off + i] * kappa[off + i];
                return std::sqrt(s);
            };
            outfile->Printf(
                "  [ORB-DIIS] iter=%d kappa-newton block norms: ca=%.3e cv=%.3e av=%.3e aa=%.3e\n",
                iteration,
                block_norm(o0, Ncore * n_act),
                block_norm(o1, Ncore * n_virt),
                block_norm(o2, n_act * n_virt),
                block_norm(o3, n_rot_aa));
        }

        return kappa;
    }

    /// cum_kappa = total rotation from the frame anchor Ca_ref; error = Newton step at this
    /// iterate. Mirrors the core store, keeping cum_kappa_ keyed by core slot.
    void store(std::vector<double> cum_kappa, std::vector<double> error) {
        DiisCore::StoreResult sr = core_.store(std::move(error));
        if (sr.wrote == static_cast<int>(cum_kappa_.size()))
            cum_kappa_.push_back(std::move(cum_kappa));
        else
            cum_kappa_[sr.wrote] = std::move(cum_kappa);
    }

    /// Removes the entry written by the last store(); no-op if no store is tracked.
    void drop_newest() {
        int dropped = core_.drop_newest();
        if (dropped >= 0 && dropped < static_cast<int>(cum_kappa_.size()))
            cum_kappa_.erase(cum_kappa_.begin() + dropped);
    }

    /// Pulay DIIS extrapolation: kappa^opt = sum_i c_i cum_kappa_i, c_i from the core solve
    /// (c_i = 0 on the slots the conditioning chain pruned).
    /// Returns an empty vector when the core declines to extrapolate.
    std::vector<double> extrapolate() {
        auto combine = [this](const std::vector<double>& c) -> std::vector<double> {
            const int n = core_.count();
            const int kappa_size = static_cast<int>(cum_kappa_[0].size());
            std::vector<double> kappa_opt(kappa_size, 0.0);
            for (int i = 0; i < n; ++i) {
                if (c[i] == 0.0) continue;
                C_DAXPY(kappa_size, c[i], cum_kappa_[i].data(), 1, kappa_opt.data(), 1);
            }
            return kappa_opt;
        };
        return core_.extrapolate(chain_, combine);
    }

    /// Rescales kappa_opt uniformly to max_i |kappa_opt[i]| <= max_step_kappa; true if the
    /// cap fired.
    bool apply_step_cap(std::vector<double>& kappa_opt) const {
        if (kappa_opt.empty()) return false;
        // C_IDAMAX (psi4 wrapper) returns a 0-based index, unlike Fortran IDAMAX.
        const double max_abs = std::abs(kappa_opt[C_IDAMAX(kappa_opt.size(), kappa_opt.data(), 1)]);
        if (max_abs > config_.max_step_kappa) {
            C_DSCAL(kappa_opt.size(), config_.max_step_kappa / max_abs, kappa_opt.data(), 1);
            return true;
        }
        return false;
    }

    /// Build Ca_new = Ca_ref * exp(K_cum) using Pade-4 with scaling/squaring.
    ///
    /// K is skew with an identically zero (virt,virt) block, K = [[A, B], [-B^T, 0]] over
    /// O = [0, first_virt) and V. Its range lies in span(e_O) + range(B), of dimension at
    /// most r = 2*no with no = |O| = Ncore + n_act. On that invariant subspace
    ///
    ///   exp(K) = I + U (exp(T) - I) U^T,   T = U^T K U   (r x r, skew)
    ///
    /// is exact, with U = [[I_no, 0], [0, Qb]] and Qb an orthonormal basis of range(B^T),
    /// nv x m. This turns the nmo x nmo exponential into an r x r one: for
    /// Ca_ref exp(K) = Ca_ref + (Ca_ref U)(exp(T) - I) U^T, the nmo^3 product becomes
    /// nso*r*nmo. r = no + m never exceeds nmo; at r = nmo (Ncore = 0, 2*n_act = nmo)
    /// U is square orthogonal and the expression remains exact.
    SharedMatrix apply_rotation(SharedMatrix Ca_ref,
                                const std::vector<double>& kappa,
                                const std::vector<int>& active_mo_indices,
                                int Ncore, int nmo) const {
        const int n_act = static_cast<int>(active_mo_indices.size());
        const int no = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;  // first_virt
        const int nv = nmo - no;

        auto K = std::make_shared<Matrix>("K_cum (orbital DIIS)", nmo, nmo);
        unpack_kappa_to_K(kappa, K, active_mo_indices, Ncore, nmo);
        double** Kp = K->pointer(0);

        // Qb = orthonormal basis of range(B^T), B = K[0:no, no:nmo] (B^T is nv x no).
        // Column-pivoted Householder QR: |R_kk| descending sets rank m via a relative
        // tolerance. Q is column-major, nv x no, ld = nv.
        std::vector<double> Q(static_cast<size_t>(nv) * no);
        for (int o = 0; o < no; ++o)
            for (int v = 0; v < nv; ++v) Q[static_cast<size_t>(o) * nv + v] = Kp[o][no + v];

        int m = 0;
        {
            std::vector<int> jpvt(no, 0);
            std::vector<double> tau(std::min(nv, no));
            double wkopt = 0.0;
            C_DGEQP3(nv, no, Q.data(), nv, jpvt.data(), tau.data(), &wkopt, -1);
            int lwork = std::max(1, static_cast<int>(wkopt));
            std::vector<double> work(lwork);
            C_DGEQP3(nv, no, Q.data(), nv, jpvt.data(), tau.data(), work.data(), lwork);

            const int kmax = std::min(nv, no);
            const double r00 = (kmax > 0) ? std::abs(Q[0]) : 0.0;
            const double tol = std::max(1e-12, 1e-12 * r00 * std::sqrt(static_cast<double>(nv)));
            while (m < kmax && std::abs(Q[static_cast<size_t>(m) * nv + m]) > tol) ++m;

            if (m > 0) {
                C_DORGQR(nv, m, m, Q.data(), nv, tau.data(), &wkopt, -1);
                lwork = std::max(1, static_cast<int>(wkopt));
                work.resize(lwork);
                C_DORGQR(nv, m, m, Q.data(), nv, tau.data(), work.data(), lwork);
            }
        }

        const int r = no + m;

        // T = U^T K U = [[A, B Qb], [-(B Qb)^T, 0]], with U = [[I_no, 0], [0, Qb]].
        auto T = std::make_shared<Matrix>("T (invariant subspace)", r, r);
        double** Tp = T->pointer(0);
        T->zero();
        for (int i = 0; i < no; ++i)
            for (int j = 0; j < no; ++j) Tp[i][j] = Kp[i][j];
        if (m > 0) {
            C_DGEMM('N', 'T', no, m, nv, 1.0, &Kp[0][no], nmo, Q.data(), nv, 0.0, &Tp[0][no], r);
            for (int i = 0; i < no; ++i)
                for (int j = 0; j < m; ++j) Tp[no + j][i] = -Tp[i][no + j];
        }

        T->expm(4, true);
        for (int i = 0; i < r; ++i) Tp[i][i] -= 1.0;  // exp(T) - I

        // W = Ca_ref U = [Ca_O | Ca_V Qb]   (nso x r, ld r)
        const int nso = Ca_ref->rowspi()[0];
        double** Cr = Ca_ref->pointer(0);
        std::vector<double> W(static_cast<size_t>(nso) * r);
        for (int mu = 0; mu < nso; ++mu)
            for (int o = 0; o < no; ++o) W[(size_t)mu * r + o] = Cr[mu][o];
        if (m > 0)
            C_DGEMM('N', 'T', nso, m, nv, 1.0, &Cr[0][no], nmo, Q.data(), nv, 0.0,
                    W.data() + no, r);

        // M = (exp(T) - I) U^T   (r x nmo, ld nmo): the first no columns land on O, the rest
        // go back through Qb^T.
        std::vector<double> M(static_cast<size_t>(r) * nmo, 0.0);
        for (int i = 0; i < r; ++i)
            for (int o = 0; o < no; ++o) M[(size_t)i * nmo + o] = Tp[i][o];
        if (m > 0)
            C_DGEMM('N', 'N', r, nv, m, 1.0, &Tp[0][no], r, Q.data(), nv, 0.0,
                    M.data() + no, nmo);

        // Ca_new = Ca_ref + W M
        auto Ca_new = Ca_ref->clone();
        Ca_new->set_name("Ca (orbital DIIS)");
        C_DGEMM('N', 'N', nso, nmo, r, 1.0, W.data(), r, M.data(), nmo, 1.0,
                Ca_new->pointer(0)[0], nmo);
        return Ca_new;
    }

    void reset() {
        core_.reset();
        cum_kappa_.clear();
    }

    int count() const { return core_.count(); }
    double last_max_error() const { return core_.last_max_error(); }
    double last_c_norm_sq() const { return core_.last_c_norm_sq(); }

    std::vector<double> last_coefficients() const { return core_.last_coefficients(); }

   private:
    DiisConfig config_;
    DiisCore core_;
    ConditioningChain chain_;
    std::vector<std::vector<double>> cum_kappa_;  ///< cumulative-kappa payload, keyed by core slot
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_ORBITAL_DIIS_H
