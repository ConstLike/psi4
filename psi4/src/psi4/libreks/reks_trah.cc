/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "reks_trah.h"

#include "reks_fon_derivs.h"
#include "reks_math.h"
#include "reks_objectives.h"
#include "reks_report.h"
#include "reks_report_level.h"
#include "reks_scf.h"
#include "reks_si_types.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace psi {
namespace reks {

using report::reports;

namespace {

/// g = [g_orb; g_fon(b)...], H = [[H_oo, H_of(b)...], [H_fo(b), H_ff(b)]] (row-major, N = n_rot + sum_b count):
/// kappa first, FON block b at [offset, offset + count); FON x FON entries of different blocks are zero.
struct CombinedGradient {
    std::vector<studio::FonBlock> fon_blocks;   ///< (s, gen) order
    std::vector<double> g;
    std::vector<double> H;
    /// Diagnostics, filled only when want_diagnostics.
    std::vector<double> hoo_evals_raw;  ///< H_oo eigenvalues
    std::vector<std::vector<double>> fon_hess_evals;  ///< [block] FON-Hessian eigenvalues
    std::vector<std::vector<double>> fon_hess_block;  ///< [block] FON-Hessian block (count^2)
};

/// g_orb(i,j) = -2 (F_gen[ai][aj] - F_gen[aj][ai]) = dE/dkappa_ij, zero at convergence (Brillouin).
std::vector<double> extract_orbital_gradient(SharedMatrix F_gen,
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

/// H_of[(pq) * n_fon + i] = -2 sum_L dC_L/dn_i G_L(p,q), G_L(p,q) = n^a_Lp Fa_L[p][q] + n^b_Lp Fb_L[p][q]
/// - n^a_Lq Fa_L[q][p] - n^b_Lq Fb_L[q][p], over the active rows of the per-L MO Fock.
std::vector<double> cross_hessian_from(const fon::WeightDerivs& wd,
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

/// Frozen-Fock 1e diagonal: H_oo[(ij),(ij)] = -2 sum_L C_L sum_s (n^s_Li - n^s_Lj)(F^s_L[ai][ai] - F^s_L[aj][aj]),
/// the r1 == r2 case of A1e in compute_orbital_hessian_exact.
std::vector<double> compute_orbital_hessian_diag(const std::vector<SharedMatrix>& F_MO_arows_a_L,
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

/// H_oo = -2 sym(A1e) + 4 sym(A2e) at fixed FON, dn^s_ij = n^s_iL - n^s_jL, (ab|cd) = eri.at(a,b,c,d):
/// A1e = sum_L C_L sum_s dn^s_ij [d_iq F^s_L[p,j] - d_ip F^s_L[q,j] + d_jq F^s_L[i,p] - d_jp F^s_L[i,q]],
/// A2e = sum_L C_L sum_ss' dn^s_ij dn^s'_pq [(ij|qp) + (ij|pq) - d_ss' x_alpha ((ip|qj) + (iq|pj))] / 2;
/// f_xc is omitted (approximate for DFT); eri covers every active pair k <= l.
std::vector<double> compute_orbital_hessian_exact(
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

    // A1e/A2e per microstate L and pair (r1, r2) = ((ij), (pq)): Kronecker deltas -> index guards,
    // d_ss' -> same-spin (Coulomb - exchange) vs opposite-spin (Coulomb) terms.
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

/// Gradient and Hessian in the CombinedGradient layout; want_diagnostics fills the raw Hessian arrays,
/// the exact H_oo among them only when active_eri (every active pair k <= l) is given.
CombinedGradient combined_gradient(SharedMatrix                       F_gen,
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

    // Non-empty FON blocks in (s, gen) order, offsets after the n_rot rotations.
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
    result.g.assign(N, 0.0);
    result.H.assign(N * N, 0.0);
    result.fon_hess_evals.assign(result.fon_blocks.size(), {});
    result.fon_hess_block.assign(result.fon_blocks.size(), {});

    auto g_orb = extract_orbital_gradient(F_gen, active);
    for (int i = 0; i < n_rot; ++i) result.g[i] = g_orb[i];

    // H_oo[i][i] = max(|h_i|, H_oo_min), h_i the frozen-Fock 1e diagonal.
    auto H_oo_diag = compute_orbital_hessian_diag(F_MO_arows_a_L, F_MO_arows_b_L,
                                                  C_L, sa_cassette, active);
    constexpr double H_oo_min = 0.1;
    for (int i = 0; i < n_rot; ++i) {
        const double h = H_oo_diag[i];
        result.H[i * N + i] = std::max(std::abs(h), H_oo_min);
    }
    // Spectrum of the unregularized H_oo (1e + 2e).
    if (want_diagnostics && n_rot > 0 && active_eri && !active_eri->empty()) {
        auto H_oo_exact = compute_orbital_hessian_exact(F_MO_arows_a_L, F_MO_arows_b_L,
                                                        C_L, *active_eri, sa_cassette,
                                                        active, x_alpha);
        std::vector<double> evals, evecs;
        lapack_diagonalize(H_oo_exact.data(), n_rot, evals, evecs);
        result.hoo_evals_raw = std::move(evals);
    }

    // Per block: g_fon(b), H_ff(b), H_of(b) from the block's sector FON snapshot.
    for (size_t bi = 0; bi < result.fon_blocks.size(); ++bi) {
        const studio::FonBlock& blk = result.fon_blocks[bi];
        const int nf  = blk.count;
        if (nf == 0) continue;
        const int off = blk.offset;
        const int gen = blk.gen;
        const studio::FONSnapshot& fons = reel.fon_state[blk.s];

        const fon::WeightDerivs wd = fon::accumulate_weight_derivs(blk.s, gen, sa_cassette, fons);

        auto g_fon = fon::fon_gradient_from(wd, sa_cassette, E_L);
        for (int k = 0; k < nf; ++k) result.g[off + k] = g_fon[k];

        auto H_ff = fon::fon_hessian_from(wd, blk.s, gen, sa_cassette, fons, E_L);
        for (int i = 0; i < nf; ++i)
            for (int j = 0; j < nf; ++j)
                result.H[(off + i) * N + (off + j)] = H_ff[i * nf + j];

        auto H_of = cross_hessian_from(wd, F_MO_arows_a_L, F_MO_arows_b_L,
                                       sa_cassette, active);
        for (int i = 0; i < n_rot; ++i)
            for (int k = 0; k < nf; ++k) {
                const double v = H_of[i * nf + k];
                result.H[i * N + (off + k)] = v;
                result.H[(off + k) * N + i] = v;
            }

        if (want_diagnostics) {
            std::vector<double> evals, evecs;
            lapack_diagonalize(H_ff.data(), nf, evals, evecs);  // leaves H_ff intact
            result.fon_hess_evals[bi] = std::move(evals);
            result.fon_hess_block[bi] = std::move(H_ff);
        }
    }

    return result;
}

}  // anonymous namespace

void TrahDriver::configure() {
    state_ = TRAHState();
    state_.trust_radius = 0.1;
}

void TrahDriver::resize_fon_blocks(const std::vector<studio::FonBlock>& fon_blocks) {
    // Seeded once; a block layout of the same size keeps the accumulated EMA.
    if (ahc_pred_fon_frac_smooth_.size() != fon_blocks.size()) {
        ahc_pred_fon_frac_smooth_.assign(fon_blocks.size(), 0.0);
        for (size_t b = 0; b < fon_blocks.size(); ++b)
            ahc_pred_fon_frac_smooth_[b] = (fon_blocks[b].gen == 0) ? 0.5 : 0.0;
    }
}

void TrahDriver::step(TrahStepContext& ctx) {
    timer_on("REKS: combined_step");

    int n_act = static_cast<int>(ctx.active.size());
    int n_rot = n_act * (n_act - 1) / 2;
    // Joint vector [n_rot | FON blocks in (s, gen) order], block b at [offset, offset + count).
    std::vector<reks::studio::FonBlock> blocks = ctx.fon_blocks;
    int N = n_rot;
    for (auto& blk : blocks) { blk.offset = N; N += blk.count; }
    // FON floor per generation: REKS_<G>_FON if >= 0, else 0.
    auto fon_lower = [&](int gen) {
        return ctx.fon_lower[gen] >= 0.0 ? ctx.fon_lower[gen] : 0.0;
    };
    // SA sector and its generation-0 FON count.
    const int sa_s = ctx.sa_cassette.sector();
    const int n_fon = static_cast<int>(ctx.sa_cassette.geminals_active_gen(0).size());
    bool read_guess = ctx.read_guess;

    if (reports(4) && ctx.iteration <= 2) {
        int n_sa = ctx.n_sa_microstates;
        std::string fon_counts;
        for (const auto& blk : blocks)
            fon_counts += " " + reks::studio::fon_channel_name(blk.gen, blk.s) +
                          "fon=" + std::to_string(blk.count);
        outfile->Printf("  [COMBINED_STEP] entry iter=%d  n_act=%d n_rot=%d%s n_sa=%d N=%d\n",
                        ctx.iteration, n_act, n_rot, fon_counts.c_str(), n_sa, N);
        outfile->Printf("  [COMBINED_STEP]   trah_state_.initialized=%d trust_radius=%.6f bfgs_update_count=%d\n",
                        state_.initialized ? 1 : 0, state_.trust_radius, bfgs_update_count_);
        outfile->Printf("  [COMBINED_STEP]   guess_Ca_=%s  gvb_diis_active=%d\n",
                        ctx.has_guess_Ca ? "SET" : "nullptr", ctx.ctrl_active ? 1 : 0);
        outfile->Printf("  [E_MICRO] iter=%d", ctx.iteration);
        { int shown = 0;
          for (int L : ctx.sa_cassette.microstates()) {
              if (shown++ >= 8) break;
              outfile->Printf(" E[%s]=%.6f",
                              reks::report::microstate_label(ctx.sa_cassette, L).c_str(), ctx.reel.E_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [C_MICRO] iter=%d", ctx.iteration);
        { int shown = 0;
          for (int L : ctx.sa_cassette.microstates()) {
              if (shown++ >= 8) break;
              outfile->Printf(" C[%s]=%.6f",
                              reks::report::microstate_label(ctx.sa_cassette, L).c_str(), ctx.reel.C_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [FON] iter=%d", ctx.iteration);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            const auto& gem = ctx.sa_cassette.geminals_active(blk.s, blk.gen);
            for (int k = 0; k < blk.count; ++k) {
                const auto& fon = ctx.reel.fon_state[blk.s].layers[blk.gen][gem[k]];
                outfile->Printf(" %s-FON%d=(%.8f, %.8f)", chan.c_str(), k, fon.p, fon.q);
            }
        }
        outfile->Printf("\n");
    }

    // rho = (E_prev - E_now) / pred_prev; skipped when pred_prev is negligible.
    if (state_.initialized) {
        double E_SA_actual = reks::scf::compute_E_SA(ctx.sa_cassette, ctx.reel.E_L, ctx.reel.C_L);

        constexpr double pred_threshold = 1e-10;
        if (state_.prev_predicted > pred_threshold) {
            double actual_decrease = state_.prev_energy - E_SA_actual;
            double rho = actual_decrease / state_.prev_predicted;
            state_.prev_rho = rho;
            if (reports(4)) {
                outfile->Printf("  [TRAH] iter=%d E_prev=%14.8f E_now=%14.8f actual=%12.8f pred=%12.8f rho=%8.4f %s\n",
                                ctx.iteration, state_.prev_energy, E_SA_actual, actual_decrease,
                                state_.prev_predicted, rho, (rho < 0.0 ? "E_UP" : "accept"));
                outfile->Printf("  [TRAH] iter=%d pred_decomp: orb=%+.6f", ctx.iteration,
                                state_.prev_pred_orb);
                for (size_t bi = 0; bi < blocks.size() && bi < state_.prev_pred_fon.size(); ++bi) {
                    const std::string lt = reks::studio::fon_channel_name(blocks[bi].gen, blocks[bi].s);
                    outfile->Printf(" %sfon=%+.6f cross_%s=%+.6f", lt.c_str(),
                                    state_.prev_pred_fon[bi], lt.c_str(),
                                    state_.prev_pred_cross[bi]);
                }
                outfile->Printf("\n");
            }
            reks::TrustRegionSolver::update_trust_radius(E_SA_actual, state_);

            if (actual_decrease >= 0) {
                // AHC state: EMA of rho (rho_s) and of the fractions |pred_b| / sum |pred|; the integral grows by
                // AHC_I_RATE_ (1 - rho_s) up to AHC_I_MAX_ at rho_s < 0.4, decays x0.7 at rho_s > 0.8.
                ahc_rho_smooth_ = (1.0 - AHC_BETA_) * ahc_rho_smooth_ + AHC_BETA_ * rho;
                double pred_total = std::abs(state_.prev_pred_orb) + 1e-15;
                for (size_t bi = 0; bi < blocks.size() && bi < state_.prev_pred_fon.size(); ++bi)
                    pred_total += std::abs(state_.prev_pred_fon[bi]);
                double orb_frac  = std::abs(state_.prev_pred_orb)  / pred_total;
                ahc_pred_orb_frac_smooth_  = (1.0 - AHC_BETA_) * ahc_pred_orb_frac_smooth_  + AHC_BETA_ * orb_frac;
                for (size_t bi = 0; bi < blocks.size(); ++bi) {
                    double fon_frac = (bi < state_.prev_pred_fon.size())
                        ? std::abs(state_.prev_pred_fon[bi]) / pred_total : 0.0;
                    ahc_pred_fon_frac_smooth_[bi] =
                        (1.0 - AHC_BETA_) * ahc_pred_fon_frac_smooth_[bi] + AHC_BETA_ * fon_frac;
                }
                if (ahc_rho_smooth_ < 0.4) {
                    ahc_integral_ += AHC_I_RATE_ * (1.0 - ahc_rho_smooth_);
                    ahc_integral_ = std::min(ahc_integral_, AHC_I_MAX_);
                } else if (ahc_rho_smooth_ > 0.8) {
                    ahc_integral_ *= 0.7;
                }
                consecutive_bad_rho_ = 0;
            } else {
                // Energy increase: AHC untouched; 2 in a row reset BFGS.
                consecutive_bad_rho_++;
                if (consecutive_bad_rho_ >= 2 && bfgs_update_count_ > 0) {
                    bfgs_initialized_ = false;
                    bfgs_update_count_ = 0;
                    bfgs_B_.clear();
                    bfgs_prev_step_.clear();
                    consecutive_bad_rho_ = 0;
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d RESET: 2 consecutive energy increases -> diagonal H_oo\n",
                                        ctx.iteration);
                    }
                }
            }
            // 2 consecutive kDIIS energy increases reset kDIIS.
            if (kdiis_used_prev_) {
                if (actual_decrease < 0) {
                    kdiis_energy_fail_++;
                    if (kdiis_energy_fail_ >= 2) {
                        kdiis_active_ = false;
                        kdiis_count_ = 0;
                        kdiis_oldest_ = 0;
                        kdiis_stag_count_ = 0;
                        kdiis_energy_fail_ = 0;
                        if (reports(4)) {
                            outfile->Printf("  [kDIIS] iter=%d RESET: 2 consecutive energy increases\n", ctx.iteration);
                        }
                    }
                } else {
                    kdiis_energy_fail_ = 0;
                }
            }
        } else {
            if (reports(4)) {
                outfile->Printf("  [TRAH] iter=%d pred=%.2e < threshold -- skipping TR update (DIIS-dominated)\n",
                                ctx.iteration, state_.prev_predicted);
            }
        }

        state_.prev_energy = E_SA_actual;
    }

    // First-call gen-0 FON seed: a READ guess is kept, else all FON = 1.0 (p = q = 1);
    // the CI-diagonal seed is computed for the debug print only.
    if (!state_.initialized) {
        const auto& sa_gem0 = ctx.sa_cassette.geminals_active_gen(0);
        if (read_guess && !read_boundary_trap_) {
            if (reports(4)) {
                outfile->Printf(
                    "  [TRAH] iter=%d CI-init SKIPPED: FON restored from READ guess, keeping",
                    ctx.iteration);
                for (int k = 0; k < n_fon; ++k) {
                    outfile->Printf(" %.10f", ctx.reel.fon_state[sa_s].layers[0][sa_gem0[k]].p);
                }
                outfile->Printf("\n");
            }
        } else if (read_guess && read_boundary_trap_) {
            // Boundary-trap escape (previous step was zero): reset FON to 1.0, drop the READ guess.
            for (int k = 0; k < n_fon; ++k) {
                ctx.reel.fon_state[sa_s].layers[0][sa_gem0[k]] = {1.0, 2.0 - 1.0};
            }
            read_boundary_trap_ = false;
            read_guess = false;
            ctx.clear_read_guess = true;
            if (reports(4)) {
                auto ci_fon = reks::fon_init_from_ci(ctx.sa_cassette.sector(), ctx.sa_cassette, ctx.reel.E_L,
                                                     ctx.iteration);
                outfile->Printf("  [TRAH] iter=%d BOUNDARY TRAP ESCAPE: CI-init FON (raw):", ctx.iteration);
                for (int k = 0; k < n_fon; ++k) outfile->Printf(" %.6f", ci_fon[k]);
                outfile->Printf(" -> all set to 1.0 (reset to explore full landscape)\n");
            }
        } else {
            for (int k = 0; k < n_fon; ++k) {
                ctx.reel.fon_state[sa_s].layers[0][sa_gem0[k]] = {1.0, 2.0 - 1.0};
            }
            if (reports(4)) {
                auto ci_fon = reks::fon_init_from_ci(ctx.sa_cassette.sector(), ctx.sa_cassette, ctx.reel.E_L,
                                                     ctx.iteration);
                outfile->Printf("  [TRAH] iter=%d CI-init FON (raw):", ctx.iteration);
                for (int k = 0; k < n_fon; ++k) outfile->Printf(" %.6f", ci_fon[k]);
                outfile->Printf(" -> all set to 1.0 (cold path)\n");
            }
        }
    }

    auto gh = combined_gradient(ctx.F_gen, ctx.reel, ctx.sa_cassette,
                                ctx.reel.F_MO_arows_a_L, ctx.reel.F_MO_arows_b_L,
                                ctx.reel.E_L, ctx.reel.C_L, ctx.active,
                                reports(4), ctx.hess_eri,
                                ctx.x_alpha);

    // Frozen (same-role inter-geminal) rotations: g_i = 0, row/column i of H zeroed, H_ii = 1.
    const std::vector<char>& aa_frozen = ctx.aa_frozen;
    auto is_frozen = [&](int i) { return i < static_cast<int>(aa_frozen.size()) && aa_frozen[i] != 0; };
    int n_rot_frozen = 0;
    for (int i = 0; i < n_rot; ++i) {
        if (!is_frozen(i)) continue;
        ++n_rot_frozen;
        gh.g[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            gh.H[i * N + j] = 0.0;
            gh.H[j * N + i] = 0.0;
        }
        gh.H[i * N + i] = 1.0;
    }
    if (reports(4) && n_rot_frozen > 0 && ctx.iteration <= 2) {
        outfile->Printf("  [COMBINED_STEP] iter=%d same-role projection: %d of %d rotations determined\n",
                        ctx.iteration, n_rot_frozen, n_rot);
    }

    // x_current = [kappa = 0 | FON values p].
    std::vector<double> x_current(N, 0.0);
    for (const auto& blk : blocks) {
        const auto& gem = ctx.sa_cassette.geminals_active(blk.s, blk.gen);
        for (int k = 0; k < blk.count; ++k)
            x_current[blk.offset + k] = ctx.reel.fon_state[blk.s].layers[blk.gen][gem[k]].p;
    }

    // Bounds: kappa in [-max_kappa, max_kappa]; FON in [fon_lower(gen), 2].
    double max_kappa = 0.3;
    std::vector<double> lower(N), upper(N);
    for (int i = 0; i < n_rot; ++i) {
        lower[i] = -max_kappa;
        upper[i] = max_kappa;
    }
    for (const auto& blk : blocks) {
        const double lo = fon_lower(blk.gen);
        for (int k = 0; k < blk.count; ++k) {
            lower[blk.offset + k] = lo;
            upper[blk.offset + k] = 2.0;
        }
    }

    // all_at_boundary: every FON pinned at a bound with the gradient pointing outward.
    // do_odb: ||g_orb|| fell by < 5% since the previous boundary iteration.
    bool all_at_boundary = (N > n_rot);
    bool do_odb = false;
    {
        for (size_t bi = 0; bi < blocks.size() && all_at_boundary; ++bi) {
            const auto& blk = blocks[bi];
            const auto& gem = ctx.sa_cassette.geminals_active(blk.s, blk.gen);
            for (int k = 0; k < blk.count && all_at_boundary; ++k) {
                int idx = blk.offset + k;
                double fon = ctx.reel.fon_state[blk.s].layers[blk.gen][gem[k]].p;
                double g_fon_k = gh.g[idx];
                bool at_upper = (fon >= upper[idx] - 1e-8) && (g_fon_k < 0);
                bool at_lower = (fon <= lower[idx] + 1e-8) && (g_fon_k > 0);
                if (!at_upper && !at_lower) all_at_boundary = false;
            }
        }

        if (all_at_boundary) {
            double g_orb_sq = 0.0;
            for (int i = 0; i < n_rot; ++i) g_orb_sq += gh.g[i] * gh.g[i];
            double g_orb_norm = std::sqrt(g_orb_sq);

            if (prev_gorb_at_boundary_ > 0.0) {
                do_odb = (g_orb_norm >= 0.95 * prev_gorb_at_boundary_);
            }
            prev_gorb_at_boundary_ = g_orb_norm;
        } else {
            prev_gorb_at_boundary_ = -1.0;
        }
    }

    if (do_odb) {
        kdiis_stag_count_++;
        if (reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag=%d/%d buf=%d/%d\n", ctx.iteration, kdiis_stag_count_, KDIIS_TRIGGER,
                            kdiis_count_, KDIIS_MAX_VEC);
        }
    } else if (!all_at_boundary) {
        if (kdiis_stag_count_ > 0 && reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag RESET: left boundary (was %d)\n", ctx.iteration, kdiis_stag_count_);
        }
        kdiis_stag_count_ = 0;
        if (kdiis_active_) {
            kdiis_active_ = false;
            kdiis_count_ = 0;
            kdiis_oldest_ = 0;
            kdiis_energy_fail_ = 0;
            if (reports(4)) {
                outfile->Printf("  [kDIIS] iter=%d deactivated: left boundary\n", ctx.iteration);
            }
        }
    } else {
        if (kdiis_stag_count_ > 0 && reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag RESET: g_orb improved (was %d)\n", ctx.iteration,
                            kdiis_stag_count_);
        }
        kdiis_stag_count_ = 0;
    }

    // AHC shift per block: H_diag += mean(H_diag) (P + I) frac_b, P = max(0, 1/rho_s - 1), I = ahc_integral_;
    // signed mean over free rotations, |.| mean over a FON block (H_ff may be negative).
    {
        double ahc_p = std::max(0.0, 1.0 / std::max(ahc_rho_smooth_, 0.01) - 1.0);
        double ahc_correction = ahc_p + ahc_integral_;
        if (ahc_correction > 1e-4) {
            double h_oo_mean = 0.0;
            int n_rot_free = 0;
            for (int i = 0; i < n_rot; ++i) {
                if (is_frozen(i)) continue;
                h_oo_mean += gh.H[i * N + i];
                ++n_rot_free;
            }
            if (n_rot_free > 0) h_oo_mean /= n_rot_free;

            double lambda_orb = h_oo_mean * ahc_correction * ahc_pred_orb_frac_smooth_;
            for (int i = 0; i < n_rot; ++i)
                if (!is_frozen(i)) gh.H[i * N + i] += lambda_orb;

            std::vector<double> lambda_fon(blocks.size(), 0.0);
            for (size_t bi = 0; bi < blocks.size(); ++bi) {
                const int nf = blocks[bi].count;
                if (nf == 0) continue;
                const int off = blocks[bi].offset;
                double h_mean = 0.0;
                for (int k = 0; k < nf; ++k)
                    h_mean += std::abs(gh.H[(off + k) * N + (off + k)]);
                h_mean /= nf;
                lambda_fon[bi] = h_mean * ahc_correction * ahc_pred_fon_frac_smooth_[bi];
                for (int k = 0; k < nf; ++k)
                    gh.H[(off + k) * N + (off + k)] += lambda_fon[bi];
            }

            if (reports(4)) {
                outfile->Printf(
                    "  [AHC] iter=%d rho_s=%.4f P=%.4f I=%.4f corr=%.4f lam_orb=%.6f",
                    ctx.iteration, ahc_rho_smooth_, ahc_p, ahc_integral_, ahc_correction, lambda_orb);
                for (size_t bi = 0; bi < blocks.size(); ++bi)
                    outfile->Printf(" lam_%sfon=%.6f",
                                    reks::studio::fon_channel_name(blocks[bi].gen, blocks[bi].s).c_str(),
                                    lambda_fon[bi]);
                outfile->Printf("\n");
            }
        }
    }

    // Powell-damped BFGS of H_oo below max|g_orb| < 0.25 (Chaban, G. et al. Theor. Chem. Acc. 1997, 97, 88);
    // H_ff and H_of stay exact.
    {
        double g_orb_max = 0.0;
        for (int i = 0; i < n_rot; ++i) g_orb_max = std::max(g_orb_max, std::abs(gh.g[i]));
        constexpr double BFGS_ACTIVATION_THRESHOLD = 0.25;
        if (g_orb_max >= BFGS_ACTIVATION_THRESHOLD) {
            bfgs_initialized_ = false;
            bfgs_update_count_ = 0;
            bfgs_B_.clear();
        }
    }
    if (n_rot > 0) {
        if (bfgs_initialized_ && !bfgs_prev_step_.empty()) {
            // s = prev kappa step, y = orbital gradient difference.
            std::vector<double> s_orb(n_rot), y_orb(n_rot);
            for (int i = 0; i < n_rot; ++i) {
                s_orb[i] = bfgs_prev_step_[i];
                y_orb[i] = gh.g[i] - bfgs_prev_g_[i];
            }

            double s_norm2 = 0.0;
            for (int i = 0; i < n_rot; ++i) s_norm2 += s_orb[i] * s_orb[i];

            // READ guess: B_oo diagonal *= prev_rho before the first update.
            if (bfgs_update_count_ == 0 && read_guess) {
                double scale = std::max(state_.prev_rho, 1.0);
                if (scale > 1.5) {
                    for (int i = 0; i < n_rot; ++i)
                        bfgs_B_[i * n_rot + i] *= scale;
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d pre-update H_oo scaled by rho=%.2f "
                                        "(READ curvature correction)\n", ctx.iteration, scale);
                    }
                }
            }

            if (std::sqrt(s_norm2) > 1e-8) {
                std::vector<double> Bs(n_rot, 0.0);
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) Bs[i] += bfgs_B_[i * n_rot + j] * s_orb[j];

                double sTBs = 0.0, yTs = 0.0;
                for (int i = 0; i < n_rot; ++i) {
                    sTBs += s_orb[i] * Bs[i];
                    yTs += y_orb[i] * s_orb[i];
                }

                // Powell damping: y <- theta y + (1 - theta) B s when y^T s < 0.2 s^T B s.
                double theta = 1.0;
                if (sTBs > 1e-15 && yTs < 0.2 * sTBs) {
                    theta = 0.8 * sTBs / (sTBs - yTs);
                    for (int i = 0; i < n_rot; ++i) y_orb[i] = theta * y_orb[i] + (1.0 - theta) * Bs[i];
                    yTs = 0.0;
                    for (int i = 0; i < n_rot; ++i) yTs += y_orb[i] * s_orb[i];
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d Powell damping: theta=%.4f yTs_new=%.6f\n",
                                        ctx.iteration, theta, yTs);
                    }
                }

                // B <- B - (Bs)(Bs)^T/(s^TBs) + (y)(y)^T/(y^Ts)
                if (sTBs > 1e-15 && yTs > 1e-15) {
                    for (int i = 0; i < n_rot; ++i)
                        for (int j = 0; j < n_rot; ++j)
                            bfgs_B_[i * n_rot + j] += -Bs[i] * Bs[j] / sTBs + y_orb[i] * y_orb[j] / yTs;
                    bfgs_update_count_++;

                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d update #%d: sTBs=%.6f yTs=%.6f theta=%.4f B_oo_diag=[",
                                        ctx.iteration, bfgs_update_count_, sTBs, yTs, theta);
                        for (int i = 0; i < n_rot; ++i) outfile->Printf(" %.4f", bfgs_B_[i * n_rot + i]);
                        outfile->Printf(" ]\n");
                    }
                } else if (reports(4)) {
                    outfile->Printf("  [BFGS] iter=%d update SKIPPED: sTBs=%.2e yTs=%.2e\n", ctx.iteration, sTBs, yTs);
                }
            } else if (reports(4)) {
                outfile->Printf("  [BFGS] iter=%d update SKIPPED: ||s_orb||=%.2e < 1e-8\n", ctx.iteration,
                                std::sqrt(s_norm2));
            }

            // BFGS H_oo replaces the diagonal one after the first successful update.
            if (bfgs_update_count_ > 0) {
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) gh.H[i * N + j] = bfgs_B_[i * n_rot + j];
            } else {
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) bfgs_B_[i * n_rot + j] = gh.H[i * N + j];
                constexpr double RHO_SCALE_THRESHOLD_SYNC = 2.0;
                if (read_guess && state_.prev_rho > RHO_SCALE_THRESHOLD_SYNC) {
                    double scale = state_.prev_rho;
                    for (int i = 0; i < n_rot; ++i) {
                        bfgs_B_[i * n_rot + i] *= scale;
                        gh.H[i * N + i] *= scale;
                    }
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d re-sync H_oo scaled by rho=%.2f "
                                        "(READ curvature correction)\n", ctx.iteration, scale);
                    }
                }
            }
        } else if (!bfgs_initialized_) {
            // B_oo <- H_oo; READ guess: diagonal *= prev_rho.
            bfgs_B_.resize(n_rot * n_rot);
            for (int i = 0; i < n_rot; ++i)
                for (int j = 0; j < n_rot; ++j) bfgs_B_[i * n_rot + j] = gh.H[i * N + j];

            constexpr double RHO_SCALE_THRESHOLD = 2.0;
            if (read_guess && state_.prev_rho > RHO_SCALE_THRESHOLD) {
                double scale = state_.prev_rho;
                for (int i = 0; i < n_rot; ++i) {
                    bfgs_B_[i * n_rot + i] *= scale;
                    gh.H[i * N + i] *= scale;
                }
                if (reports(4)) {
                    outfile->Printf("  [BFGS] iter=%d initialized B_oo(%dx%d), H_oo scaled by rho=%.2f "
                                    "(READ curvature correction)\n", ctx.iteration, n_rot, n_rot, scale);
                }
            } else if (reports(4)) {
                outfile->Printf("  [BFGS] iter=%d initialized B_oo(%dx%d) from diagonal Hessian\n",
                                ctx.iteration, n_rot, n_rot);
            }
            bfgs_initialized_ = true;
        }
    }

    // Gradient before KKT elimination.
    bfgs_prev_g_.assign(gh.g.begin(), gh.g.end());

    // Ring buffer of orbital gradients for kDIIS.
    if (all_at_boundary && n_rot > 0) {
        if (static_cast<int>(kdiis_g_history_.size()) < KDIIS_MAX_VEC) {
            kdiis_g_history_.resize(KDIIS_MAX_VEC);
        }
        int write_idx;
        if (kdiis_count_ < KDIIS_MAX_VEC) {
            write_idx = kdiis_count_;
            kdiis_count_++;
        } else {
            write_idx = kdiis_oldest_;
            kdiis_oldest_ = (kdiis_oldest_ + 1) % KDIIS_MAX_VEC;
        }
        kdiis_g_history_[write_idx].assign(gh.g.begin(), gh.g.begin() + n_rot);
        if (reports(5)) {
            double gn = 0;
            for (int r = 0; r < n_rot; ++r) gn += gh.g[r] * gh.g[r];
            outfile->Printf("  [kDIIS] iter=%d stored g_orb ||g||=%.6e buf=%d/%d\n", ctx.iteration, std::sqrt(gn),
                            kdiis_count_, KDIIS_MAX_VEC);
        }
    }

    // At boundary stagnation: zero the orbital x FON blocks of H.
    const int n_fon_total = N - n_rot;
    if (do_odb && n_fon_total > 0) {
        for (int i = 0; i < n_rot; ++i) {
            for (int c = n_rot; c < N; ++c) {
                gh.H[i * N + c] = 0.0;
                gh.H[c * N + i] = 0.0;
            }
        }
        if (reports(4)) {
            double g_orb_sq = 0.0;
            for (int i = 0; i < n_rot; ++i) g_orb_sq += gh.g[i] * gh.g[i];
            double g_fon_sq = 0.0;
            for (int c = n_rot; c < N; ++c) g_fon_sq += gh.g[c] * gh.g[c];
            outfile->Printf("  [AUTOPILOT] iter=%d BOUNDARY ODB: ||g_orb||=%.4e ||g_fon||=%.4e"
                            " (zeroed orbital-FON cross-coupling, all generations)\n",
                            ctx.iteration, std::sqrt(g_orb_sq), std::sqrt(g_fon_sq));
        }
    }

    // KKT: FON locked at a bound with the gradient pointing outward: g_i = 0, row/column i zeroed, H_ii = 1.
    for (const auto& blk : blocks) {
        const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
        for (int k = 0; k < blk.count; ++k) {
            int idx = blk.offset + k;
            double fon_k = x_current[idx];
            double g_fon_k = gh.g[idx];
            bool locked_upper = (fon_k >= upper[idx] - 1e-8) && (g_fon_k < 0);
            bool locked_lower = (fon_k <= lower[idx] + 1e-8) && (g_fon_k > 0);
            if (locked_upper || locked_lower) {
                gh.g[idx] = 0.0;
                for (int j = 0; j < N; ++j) {
                    gh.H[idx * N + j] = 0.0;
                    gh.H[j * N + idx] = 0.0;
                }
                gh.H[idx * N + idx] = 1.0;
                if (reports(4)) {
                    outfile->Printf("  [TRAH] iter=%d KKT-eliminate %sfon%d (f=%.6f, g=%+.6f)\n",
                                    ctx.iteration, chan.c_str(), k, fon_k, g_fon_k);
                }
            }
        }
    }

    auto step_result = reks::TrustRegionSolver::compute_step(gh.g, gh.H, lower, upper, x_current, state_);

    if (reports(4) && step_result.n_active_bounds > 0) {
        outfile->Printf("  [AUTOPILOT] iter=%d BOUNDS: %d vars clamped.", ctx.iteration, step_result.n_active_bounds);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            for (int k = 0; k < blk.count; ++k) {
                int idx = blk.offset + k;
                double fon = x_current[idx];
                double s_fon = step_result.step[idx];
                double fon_new = fon + s_fon;
                if (std::abs(s_fon) < 1e-12 && (fon >= 2.0 - 1e-8 || fon <= 1e-8)) {
                    outfile->Printf(" %sfon%d=%.4f(clamped,delta=0)", chan.c_str(), k, fon);
                } else {
                    outfile->Printf(" %sfon%d=%.4f->%.4f(delta=%+.6f)", chan.c_str(), k, fon, fon_new, s_fon);
                }
            }
        }
        double s_orb_sq = 0.0, s_fon_sq = 0.0;
        for (int i = 0; i < n_rot; ++i) s_orb_sq += step_result.step[i] * step_result.step[i];
        for (int c = n_rot; c < N; ++c) s_fon_sq += step_result.step[c] * step_result.step[c];
        double s_total = std::sqrt(s_orb_sq + s_fon_sq);
        outfile->Printf(" orb_frac=%.4f\n", s_total > 1e-15 ? std::sqrt(s_orb_sq) / s_total : 0.0);
    }

    // kDIIS on kappa at boundary stagnation (Fischer, T. H.; Almlof, J. J. Phys. Chem. 1992, 96, 9768).
    bool kdiis_used = false;
    bool aitken_fired_step = false;
    if (kdiis_stag_count_ >= KDIIS_TRIGGER && kdiis_count_ >= KDIIS_MIN_VEC && n_rot > 0) {
        kdiis_active_ = true;
        int m = kdiis_count_;

        // Error Gram B_ij = g_i^T g_j (m x m).
        std::vector<double> B_diis(m * m, 0.0);
        for (int i = 0; i < m; ++i) {
            const int ii = (kdiis_oldest_ + i) % KDIIS_MAX_VEC;
            for (int j = i; j < m; ++j) {
                const int jj = (kdiis_oldest_ + j) % KDIIS_MAX_VEC;
                const double dot = C_DDOT(n_rot, kdiis_g_history_[ii].data(), 1,
                                          kdiis_g_history_[jj].data(), 1);
                B_diis[i * m + j] = dot;
                B_diis[j * m + i] = dot;
            }
        }

        // Pulay coefficients c: min |sum_i c_i g_i|^2, sum_i c_i = 1.
        std::vector<double> diis_c;
        double c_norm_sq = 0.0;
        bool solve_ok = reks::diis_bordered_coefficients(B_diis.data(), m, diis_c, c_norm_sq);

        // Reject if any |c_i| > 10.
        if (solve_ok) {
            for (int i = 0; i < m; ++i) {
                if (std::abs(diis_c[i]) > 10.0) {
                    solve_ok = false;
                    break;
                }
            }
        }

        if (solve_ok) {
            // g_interp = sum_i c_i g_i.
            std::vector<double> g_interp(n_rot, 0.0);
            for (int i = 0; i < m; ++i) {
                int ii = (kdiis_oldest_ + i) % KDIIS_MAX_VEC;
                for (int r = 0; r < n_rot; ++r) {
                    g_interp[r] += diis_c[i] * kdiis_g_history_[ii][r];
                }
            }

            // H_oo kappa = -g_interp by Cholesky on the BFGS H_oo; diagonal fallback if not positive definite.
            std::vector<double> kappa_diis(n_rot, 0.0);
            if (bfgs_update_count_ > 0) {
                // C_DPOSV overwrites A_oo; row-major, so uplo 'U' reads the lower triangle.
                std::vector<double> A_oo(static_cast<size_t>(n_rot) * n_rot);
                for (int i = 0; i < n_rot; ++i)
                    C_DCOPY(n_rot, &gh.H[static_cast<size_t>(i) * N], 1,
                            &A_oo[static_cast<size_t>(i) * n_rot], 1);
                for (int i = 0; i < n_rot; ++i) kappa_diis[i] = -g_interp[i];

                const int info = C_DPOSV('U', n_rot, 1, A_oo.data(), n_rot,
                                         kappa_diis.data(), n_rot);
                if (info < 0)
                    throw PSIEXCEPTION("combined_step: C_DPOSV illegal argument, info=" +
                                       std::to_string(info));
                if (info > 0) {
                    for (int i = 0; i < n_rot; ++i) {
                        double h_ii = gh.H[i * N + i];
                        kappa_diis[i] = (std::abs(h_ii) > 1e-10) ? -g_interp[i] / h_ii : 0.0;
                    }
                }
            } else {
                for (int i = 0; i < n_rot; ++i) {
                    double h_ii = gh.H[i * N + i];
                    kappa_diis[i] = (std::abs(h_ii) > 1e-10) ? -g_interp[i] / h_ii : 0.0;
                }
            }

            for (int i = 0; i < n_rot; ++i) {
                kappa_diis[i] = std::max(-max_kappa, std::min(max_kappa, kappa_diis[i]));
            }

            double trs_kappa_norm = 0.0;
            for (int i = 0; i < n_rot; ++i) trs_kappa_norm += step_result.step[i] * step_result.step[i];
            trs_kappa_norm = std::sqrt(trs_kappa_norm);
            double trs_pred = step_result.predicted_decrease;

            for (int i = 0; i < n_rot; ++i) {
                step_result.step[i] = kappa_diis[i];
            }

            step_result.predicted_decrease =
                reks::TrustRegionSolver::predicted_decrease(gh.g.data(), gh.H.data(), step_result.step.data(), N);

            kdiis_used = true;

            if (reports(4)) {
                double kdiis_norm = 0.0;
                for (int i = 0; i < n_rot; ++i) kdiis_norm += kappa_diis[i] * kappa_diis[i];
                kdiis_norm = std::sqrt(kdiis_norm);
                double g_interp_norm = 0.0, g_orb_norm = 0.0;
                for (int i = 0; i < n_rot; ++i) {
                    g_interp_norm += g_interp[i] * g_interp[i];
                    g_orb_norm += gh.g[i] * gh.g[i];
                }
                g_interp_norm = std::sqrt(g_interp_norm);
                g_orb_norm = std::sqrt(g_orb_norm);
                outfile->Printf("  [kDIIS] iter=%d active: m=%d c=[", ctx.iteration, m);
                for (int i = 0; i < m; ++i) outfile->Printf(" %.3f", diis_c[i]);
                outfile->Printf(" ]\n");
                outfile->Printf("  [kDIIS] iter=%d ||g_interp||=%.6e ||g_orb||=%.6e ratio=%.4f\n", ctx.iteration,
                                g_interp_norm, g_orb_norm, g_orb_norm > 1e-15 ? g_interp_norm / g_orb_norm : 0.0);
                outfile->Printf(
                    "  [kDIIS] iter=%d ||kappa_diis||=%.6f ||kappa_trs||=%.6f "
                    "pred_diis=%.8f pred_trs=%.8f\n",
                    ctx.iteration, kdiis_norm, trs_kappa_norm, step_result.predicted_decrease, trs_pred);
            }
        } else {
            if (reports(4)) {
                outfile->Printf("  [kDIIS] iter=%d solve FAILED (m=%d) -- using TRS step\n", ctx.iteration, m);
            }
        }
    }

    // Aitken along one soft mode: for collinear steps (cos > 0.99) kappa_inf = kappa_n / (1 - rho),
    // rho = s_n . s_{n-1} / |s_{n-1}|^2.
    std::vector<double> cur_raw_orb;
    if (n_rot > 0) cur_raw_orb.assign(step_result.step.begin(), step_result.step.begin() + n_rot);
    if (!kdiis_used && n_rot > 0 &&
        static_cast<int>(prev_orb_step_.size()) == n_rot) {
        double dot = 0.0, np2 = 0.0, nn2 = 0.0;
        for (int i = 0; i < n_rot; ++i) {
            dot += cur_raw_orb[i] * prev_orb_step_[i];
            np2 += prev_orb_step_[i] * prev_orb_step_[i];
            nn2 += cur_raw_orb[i] * cur_raw_orb[i];
        }
        if (np2 > 1e-20 && nn2 > 1e-20) {
            double rho = dot / np2;
            double cosang = dot / std::sqrt(nn2 * np2);
            if (cosang > 0.99 && rho > 0.6 && rho < 0.999) {
                double factor = std::min(1.0 / (1.0 - rho), 50.0);
                for (int i = 0; i < n_rot; ++i) step_result.step[i] *= factor;
                aitken_fired_step = true;
                if (reports(4)) {
                    outfile->Printf("  [AITKEN] iter=%d cos=%.4f rho=%.4f factor=%.2f soft-mode step extrapolation\n",
                                    ctx.iteration, cosang, rho, factor);
                }
            }
        }
    }
    // Keeps the pre-extrapolation step; cleared after a kDIIS step.
    if (kdiis_used) {
        prev_orb_step_.clear();
    } else if (n_rot > 0) {
        prev_orb_step_ = cur_raw_orb;
    }

    // Intra-geminal rotation (i,j) with |eps_i - eps_j| = gap < GAP_THRESHOLD: kappa_ij *= gap / GAP_THRESHOLD.
    {
        constexpr double GAP_THRESHOLD = 0.1;
        double* eps = ctx.epsilon->pointer(0);
        bool any_damped = false;

        std::vector<int> mo_to_unit(n_act, -1);
        for (int u = 0; u < static_cast<int>(ctx.sa_cassette.geminals_active_gen(0).size()); ++u) {
            int g = ctx.sa_cassette.geminals_active_gen(0)[u];
            for (int orb : ctx.sa_cassette.geminal_templates()[g].orbitals) mo_to_unit[orb] = u;
        }

        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                bool is_intra = (mo_to_unit[i] >= 0 && mo_to_unit[i] == mo_to_unit[j]);
                if (is_intra) {
                    int mo_i = ctx.active[i];
                    int mo_j = ctx.active[j];
                    double gap = std::abs(eps[mo_i] - eps[mo_j]);

                    if (gap < GAP_THRESHOLD) {
                        double damp = gap / GAP_THRESHOLD;
                        double kappa_before = step_result.step[idx];
                        step_result.step[idx] *= damp;
                        any_damped = true;

                        if (reports(4)) {
                            outfile->Printf(
                                "  [DELOC] iter=%d damp pair(%d,%d) gap=%.6f damp=%.4f "
                                "kappa=%.6f->%.6f\n",
                                ctx.iteration, mo_i, mo_j, gap, damp,
                                kappa_before, step_result.step[idx]);
                        }
                    }
                }
                ++idx;
            }
        }

        if (any_damped) {
            step_result.predicted_decrease =
                reks::TrustRegionSolver::predicted_decrease(gh.g.data(), gh.H.data(), step_result.step.data(), N);
        }
    }

    // READ boundary trap: zero step on the first call; the next call re-seeds FON.
    if (read_guess && !state_.initialized) {
        double step_norm_trap = 0.0;
        for (int i = 0; i < N; ++i) step_norm_trap += step_result.step[i] * step_result.step[i];
        step_norm_trap = std::sqrt(step_norm_trap);
        if (step_norm_trap < 1e-10) {
            read_boundary_trap_ = true;
            state_.prev_energy = reks::scf::compute_E_SA(ctx.sa_cassette, ctx.reel.E_L, ctx.reel.C_L);
            if (reports(4)) {
                outfile->Printf(
                    "  [TRAH] iter=%d BOUNDARY TRAP DETECTED: ||step||=%.2e, "
                    "will restart FON via CI-init on next iteration\n",
                    ctx.iteration, step_norm_trap);
            }
            timer_off("REKS: combined_step");
            return;
        }
    }

    if (!state_.initialized) {
        state_.prev_energy = reks::scf::compute_E_SA(ctx.sa_cassette, ctx.reel.E_L, ctx.reel.C_L);
        state_.initialized = true;
    }
    // An Aitken step is off the quadratic model: zero pred so the next call skips the rho update.
    if (aitken_fired_step) step_result.predicted_decrease = 0.0;
    state_.prev_predicted = step_result.predicted_decrease;

    // Model -dE = -(g^T s + 0.5 s^T H s) split into the orbital term, one FON term per block and one
    // cross term s_orb^T H_of s_fon per block.
    {
        const auto& s = step_result.step;
        const auto& g = gh.g;
        const auto& H = gh.H;

        double pred_orb = 0.0;
        for (int i = 0; i < n_rot; ++i) {
            pred_orb += g[i] * s[i];
            for (int j = 0; j < n_rot; ++j) {
                pred_orb += 0.5 * s[i] * H[i * N + j] * s[j];
            }
        }

        std::vector<double> pred_fon(blocks.size(), 0.0), pred_cross(blocks.size(), 0.0);
        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            const int b = blocks[bi].offset;
            const int e = b + blocks[bi].count;
            for (int i = b; i < e; ++i) {
                pred_fon[bi] += g[i] * s[i];
                for (int j = b; j < e; ++j) {
                    pred_fon[bi] += 0.5 * s[i] * H[i * N + j] * s[j];
                }
            }
            // Printed only.
            if (reports(4)) {
                for (int i = 0; i < n_rot; ++i) {
                    for (int j = b; j < e; ++j) {
                        pred_cross[bi] += s[i] * H[i * N + j] * s[j];
                    }
                }
            }
        }

        // Sign flip: positive = decrease.
        state_.prev_pred_orb = -pred_orb;
        state_.prev_pred_fon.assign(blocks.size(), 0.0);
        state_.prev_pred_cross.assign(blocks.size(), 0.0);
        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            state_.prev_pred_fon[bi]   = -pred_fon[bi];
            state_.prev_pred_cross[bi] = -pred_cross[bi];
        }
        state_.prev_step_on_boundary = step_result.step_on_boundary;
    }

    // fon_post = clamp(fon_pre + step, [0, 2]) per block; the first entry near a bound stops at
    // fon_bnd_margin from it.
    constexpr double fon_bnd_margin = 1e-8;
    std::vector<std::vector<double>> fon_pre(blocks.size()), fon_post(blocks.size());
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        const auto& blk = blocks[bi];
        const auto& gem = ctx.sa_cassette.geminals_active(blk.s, blk.gen);
        fon_pre[bi].resize(blk.count);
        fon_post[bi].resize(blk.count);
        for (int k = 0; k < blk.count; ++k) {
            fon_pre[bi][k] = x_current[blk.offset + k];
            double new_fon = fon_pre[bi][k] + step_result.step[blk.offset + k];
            new_fon = std::max(0.0, std::min(2.0, new_fon));
            if (new_fon >= 2.0 - fon_bnd_margin && fon_pre[bi][k] < 2.0 - fon_bnd_margin) {
                new_fon = 2.0 - fon_bnd_margin;
            }
            if (new_fon <= fon_bnd_margin && fon_pre[bi][k] > fon_bnd_margin) {
                new_fon = fon_bnd_margin;
            }
            fon_post[bi][k] = new_fon;
            ctx.reel.fon_state[blk.s].layers[blk.gen][gem[k]] = {new_fon, 2.0 - new_fon};
        }
    }

    // Ca_ <- Ca_ exp(K), K skew and active x active only: rotates the n_act active columns.
    if (n_rot > 0) {
        int nso = ctx.Ca->rowspi()[0];

        if (!ctx.reel.trah_K)
            ctx.reel.trah_K = std::make_shared<Matrix>("trah K", n_act, n_act);
        if (!ctx.reel.trah_Ca_act)
            ctx.reel.trah_Ca_act = std::make_shared<Matrix>("trah Ca_act", nso, n_act);
        if (!ctx.reel.trah_Ca_act_rot)
            ctx.reel.trah_Ca_act_rot = std::make_shared<Matrix>("trah Ca_act_rot", nso, n_act);

        double** Kp = ctx.reel.trah_K->pointer(0);
        ctx.reel.trah_K->zero();

        int idx = 0;
        double max_kappa = 0.0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                double kappa_ij = step_result.step[idx];
                Kp[i][j] = kappa_ij;
                Kp[j][i] = -kappa_ij;
                if (std::abs(kappa_ij) > max_kappa) max_kappa = std::abs(kappa_ij);
                ++idx;
            }
        }

        ctx.reel.trah_K->expm(4, true);  // K_act -> U = exp(K_act)

        // Ca_act <- Ca_act U on the gathered active columns.
        double** Cap = ctx.Ca->pointer(0);
        double** Aa = ctx.reel.trah_Ca_act->pointer(0);
        for (int i = 0; i < n_act; ++i) {
            int ai = ctx.active[i];
            for (int mu = 0; mu < nso; ++mu) Aa[mu][i] = Cap[mu][ai];
        }

        double** Rp = ctx.reel.trah_Ca_act_rot->pointer(0);
        C_DGEMM('N', 'N', nso, n_act, n_act, 1.0, Aa[0], n_act, Kp[0], n_act, 0.0, Rp[0], n_act);

        for (int j = 0; j < n_act; ++j) {
            int aj = ctx.active[j];
            for (int mu = 0; mu < nso; ++mu) Cap[mu][aj] = Rp[mu][j];
        }

        if (reports(4)) {
            outfile->Printf("  [TRAH] iter=%d Applied orbital rotation: max|kappa|=%.6f\n", ctx.iteration, max_kappa);
        }
    }

    {
        bfgs_prev_step_.resize(N);
        for (int i = 0; i < n_rot; ++i) bfgs_prev_step_[i] = step_result.step[i];
        for (size_t bi = 0; bi < blocks.size(); ++bi)
            for (int k = 0; k < blocks[bi].count; ++k)
                bfgs_prev_step_[blocks[bi].offset + k] = fon_post[bi][k] - fon_pre[bi][k];
    }

    kdiis_used_prev_ = kdiis_used;

    if (reports(4)) {
        int it = ctx.iteration;

        outfile->Printf("  [TRAH] iter=%d", it);
        if (kdiis_used)
            outfile->Printf(" kDIIS=ON(stag=%d,buf=%d)", kdiis_stag_count_, kdiis_count_);
        else if (kdiis_stag_count_ > 0)
            outfile->Printf(" kDIIS=arming(%d/%d)", kdiis_stag_count_, KDIIS_TRIGGER);
        outfile->Printf("\n");

        double g_norm = 0.0;
        for (int i = 0; i < N; ++i) g_norm += gh.g[i] * gh.g[i];
        g_norm = std::sqrt(g_norm);

        double step_norm = 0.0;
        for (int i = 0; i < N; ++i) step_norm += step_result.step[i] * step_result.step[i];
        step_norm = std::sqrt(step_norm);

        double step_over_tr = (state_.trust_radius > 1e-15) ? step_norm / state_.trust_radius : 0.0;
        outfile->Printf("  [TRAH] iter=%d ||g||=%10.6f ||step||=%10.6f TR=%6.3f step/TR=%.4f mu=%10.6f pred=%12.8f", it,
                        g_norm, step_norm, state_.trust_radius, step_over_tr, step_result.mu,
                        step_result.predicted_decrease);
        if (step_result.n_active_bounds > 0) outfile->Printf(" bounds=%d", step_result.n_active_bounds);
        outfile->Printf("\n");

        // Step-gradient alignment per subspace (cos=-1 steepest descent, ~0 TR fighting Newton).
        {
            double dot_orb = 0.0, sn_orb = 0.0, gn_orb = 0.0;
            for (int i = 0; i < n_rot; ++i) {
                dot_orb += step_result.step[i] * gh.g[i];
                sn_orb  += step_result.step[i] * step_result.step[i];
                gn_orb  += gh.g[i] * gh.g[i];
            }
            double dot_total = 0.0;
            for (int i = 0; i < N; ++i) dot_total += step_result.step[i] * gh.g[i];
            double denom_orb   = std::sqrt(sn_orb) * std::sqrt(gn_orb);
            double denom_total = step_norm * g_norm;
            outfile->Printf("  [TRAH] iter=%d cos(step,g): total=%+.4f orb=%+.4f", it,
                            denom_total > 1e-15 ? dot_total / denom_total : 0.0,
                            denom_orb   > 1e-15 ? dot_orb   / denom_orb   : 0.0);
            for (const auto& blk : blocks) {
                if (blk.count == 0) continue;
                double dot = 0.0, sn = 0.0, gn = 0.0;
                for (int i = blk.offset; i < blk.offset + blk.count; ++i) {
                    dot += step_result.step[i] * gh.g[i];
                    sn  += step_result.step[i] * step_result.step[i];
                    gn  += gh.g[i] * gh.g[i];
                }
                double denom = std::sqrt(sn) * std::sqrt(gn);
                outfile->Printf(" %sfon=%+.4f",
                                reks::studio::fon_channel_name(blk.gen, blk.s).c_str(),
                                denom > 1e-15 ? dot / denom : 0.0);
            }
            outfile->Printf("\n");
        }

        std::vector<std::string> labels(N);
        {
            int idx = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    labels[idx] = "k(" + std::to_string(ctx.active[i]) + "," +
                                  std::to_string(ctx.active[j]) + ")";
                    ++idx;
                }
            }
            for (const auto& blk : blocks) {
                const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
                const auto& gem = ctx.sa_cassette.geminals_active(blk.s, blk.gen);
                for (int k = 0; k < blk.count; ++k) {
                    const auto& orbs = ctx.sa_cassette.geminal_templates()[gem[k]].orbitals;
                    labels[blk.offset + k] = chan + "(" +
                        std::to_string(ctx.active[orbs[0]]) + "," +
                        std::to_string(ctx.active[orbs[1]]) + ")";
                }
            }
        }

        int lw = 0;
        for (const auto& lb : labels) lw = std::max(lw, static_cast<int>(lb.size()));
        lw = std::max(lw, 10);

        outfile->Printf("  [TRAH] iter=%d g = [", it);
        for (int i = 0; i < N; ++i) {
            outfile->Printf(" %10.6f", gh.g[i]);
        }
        outfile->Printf(" ]\n");

        outfile->Printf("  [TRAH] iter=%d H_diag = [", it);
        for (int i = 0; i < N; ++i) {
            outfile->Printf(" %10.6f", gh.H[i * N + i]);
        }
        outfile->Printf(" ]\n");

        outfile->Printf("  [TRAH] iter=%d h_oo_intra:", it);
        for (int k = 0; k < n_fon; ++k) {
            const auto& orbs = ctx.sa_cassette.geminal_templates()[ctx.sa_cassette.geminals_active_gen(0)[k]].orbitals;
            int p = orbs[0], q = orbs[1];
            int idx = p * (2 * n_act - p - 1) / 2 + (q - p - 1);
            outfile->Printf(" unit%d=%.6f", k, gh.H[idx * N + idx]);
        }
        outfile->Printf("\n");

        if (bfgs_initialized_ && n_rot > 1) {
            outfile->Printf("  [BFGS_DIAG] iter=%d B_oo off-diag per row:", it);
            for (int i = 0; i < n_rot; ++i) {
                double max_off = 0.0;
                for (int j = 0; j < n_rot; ++j)
                    if (i != j) max_off = std::max(max_off, std::abs(bfgs_B_[i * n_rot + j]));
                outfile->Printf(" r%d=%.4f", i, max_off);
            }
            outfile->Printf("\n");

            outfile->Printf("  [BFGS_DIAG] iter=%d step_diag_only = [", it);
            for (int i = 0; i < n_rot; ++i) {
                double hii = gh.H[i * N + i];
                double s_diag = (std::abs(hii) > 1e-15) ? -gh.g[i] / hii : 0.0;
                outfile->Printf(" %+.6f", s_diag);
            }
            outfile->Printf(" ]\n");
            outfile->Printf("  [BFGS_DIAG] iter=%d step_actual   = [", it);
            for (int i = 0; i < n_rot; ++i) {
                outfile->Printf(" %+.6f", step_result.step[i]);
            }
            outfile->Printf(" ]\n");
        }

        outfile->Printf("  [TRAH] iter=%d h_of_cross:", it);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            for (int k = 0; k < blk.count; ++k) {
                int fon_idx = blk.offset + k;
                double max_cross = 0.0;
                for (int r = 0; r < n_rot; ++r) {
                    double v = std::abs(gh.H[r * N + fon_idx]);
                    if (v > max_cross) max_cross = v;
                }
                outfile->Printf(" %sfon%d_max=%.6f", chan.c_str(), k, max_cross);
            }
        }
        outfile->Printf("\n");

        if (!gh.hoo_evals_raw.empty()) {
            outfile->Printf("  [TRAH] iter=%d H_oo evals:", it);
            for (size_t i = 0; i < gh.hoo_evals_raw.size(); ++i) {
                outfile->Printf(" %.6f", gh.hoo_evals_raw[i]);
            }
            outfile->Printf("\n");
        }

        // Per (sector, generation) block: FON-Hessian eigenvalue spectrum, then the block matrix.
        for (size_t bi = 0; bi < gh.fon_blocks.size(); ++bi) {
            const auto& fblk = gh.fon_blocks[bi];
            const std::string chan = reks::studio::fon_channel_name(fblk.gen, fblk.s);
            if (!gh.fon_hess_evals[bi].empty()) {
                outfile->Printf("  [TRAH] iter=%d H_%sf evals:", it, chan.c_str());
                for (double ev : gh.fon_hess_evals[bi]) outfile->Printf(" %.6f", ev);
                outfile->Printf("\n");
            }
            if (gh.fon_hess_block[bi].empty()) continue;
            const int nf  = fblk.count;
            const int off = fblk.offset;
            const std::vector<double>& hblk = gh.fon_hess_block[bi];
            outfile->Printf("  [TRAH] iter=%d H_%sf: %*s", it, chan.c_str(), lw, "");
            for (int j = 0; j < nf; ++j) {
                outfile->Printf(" %*s", lw, labels[off + j].c_str());
            }
            outfile->Printf("\n");
            for (int i = 0; i < nf; ++i) {
                outfile->Printf("  [TRAH] iter=%d H_%sf: %*s", it, chan.c_str(), lw, labels[off + i].c_str());
                for (int j = 0; j < nf; ++j) {
                    outfile->Printf(" %*.6f", lw, hblk[i * nf + j]);
                }
                outfile->Printf("\n");
            }
        }

        if (N <= 25) {
            outfile->Printf("  [TRAH] iter=%d H: %*s", it, lw, "");
            for (int j = 0; j < N; ++j) {
                outfile->Printf(" %*s", lw, labels[j].c_str());
            }
            outfile->Printf("\n");
            for (int i = 0; i < N; ++i) {
                outfile->Printf("  [TRAH] iter=%d H: %*s", it, lw, labels[i].c_str());
                for (int j = 0; j < N; ++j) {
                    outfile->Printf(" %*.6f", lw, gh.H[i * N + j]);
                }
                outfile->Printf("\n");
            }
        }

        outfile->Printf("  [TRAH] iter=%d kappa:", it);
        for (int i = 0; i < n_rot; ++i) outfile->Printf(" %s=%+.6f", labels[i].c_str(), step_result.step[i]);
        outfile->Printf("\n");

        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            const auto& blk = blocks[bi];
            if (blk.count == 0) continue;
            outfile->Printf("  [TRAH] iter=%d %s-FON:", it,
                            reks::studio::fon_channel_name(blk.gen, blk.s).c_str());
            for (int k = 0; k < blk.count; ++k) {
                int idx = blk.offset + k;
                outfile->Printf(" %s=%.6f%+.6f->%.6f(delta=%+.6f)",
                                labels[idx].c_str(), fon_pre[bi][k], step_result.step[idx],
                                fon_post[bi][k], fon_post[bi][k] - fon_pre[bi][k]);
            }
            outfile->Printf("\n");
        }
    }

    timer_off("REKS: combined_step");
}

}  // namespace reks
}  // namespace psi
