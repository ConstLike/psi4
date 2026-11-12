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

/// DIIS on kappa = [core-act c n_act + k | core-virt c n_virt + v' | act-virt k n_virt + v' | act-act i < j], C exp(K);
/// v' = v - first_virt; e_pq = (F_gen[q][p] - F_gen[p][q]) / -max(|-H_pq/2 + gamma_pq|, B_floor), gamma for aa, ca.
/// Ionova, I. V.; Carter, E. A. J. Chem. Phys. 1995, 102, 1251. Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472.
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
    static void mirror_skew_lower(double** K, int n) {
        constexpr int kTile = 64;
        for (int i0 = 0; i0 < n; i0 += kTile) {
            const int i1 = std::min(i0 + kTile, n);
            for (int j0 = i0; j0 < n; j0 += kTile) {
                const int j1 = std::min(j0 + kTile, n);
                for (int i = i0; i < i1; ++i)
                    for (int j = std::max(j0, i + 1); j < j1; ++j) K[j][i] = -K[i][j];
            }
        }
    }

    /// Build rows [0, first_virt) of skew-symmetric K (first_virt x nmo) from packed kappa vector.
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

        mirror_skew_lower(K, first_virt);
    }

    /// Newton rotation e_pq in the kappa layout, H_diag from the per-L MO Fock diagonals and SA weights.
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
        int iteration = -1) const {

        const int n_act = static_cast<int>(active_mo_indices.size());
        const int first_virt = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;
        const int n_virt = nmo - first_virt;

        const int n_total = n_nonred(Ncore, n_act, n_virt);
        std::vector<double> kappa(n_total, 0.0);

        double** Fgen = F_gen->pointer(0);

        Composite1eProvider curvature;
        {
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

        const DenominatorPolicy denom_policy(config_);
        const bool sat = config_.two_level_saturation;

        // x = A / denom, saturated as the 2-level angle 0.5 atan(2x) (= x - 4x^3/3 small, +-pi/4 large).
        auto step_value = [&](double A, double B_total) -> double {
            const double denom = denom_policy.shape(B_total);
            double raw = raw_step(A, denom);
            return sat ? 0.5 * std::atan(2.0 * raw) : raw;
        };

        if (report::reports(5)) {
            outfile->Printf("  [Kappa-DBG] iter=%d per-rotation Newton steps:\n", iteration);
            outfile->Printf("  %-4s %5s %5s  %12s  %12s  %12s  %12s  %12s  %12s  %s\n", "blk", "p",
                            "q", "A", "B_1e", "gamma", "B_total", "denom", "kappa", "sat?");
        }
        auto dbg_row = [&](const char* blk, int p, int q, double A, double H, double B_total,
                           double kap) {
            if (!report::reports(5)) return;
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
                double B_total = assemble_B(Block::CA, idx, curvature.ca(idx), corr);
                kappa[o0 + idx] = step_value(A, B_total);
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
                    double B_total = assemble_B(Block::CV, idx, curvature.cv(idx), corr);
                    kappa[o1 + idx] = step_value(A, B_total);
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
                    double B_total = assemble_B(Block::AV, idx, curvature.av(idx), corr);
                    kappa[o2 + idx] = step_value(A, B_total);
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
                    double B_total = assemble_B(Block::AA, idx, curvature.aa(idx), corr);
                    kappa[o3 + idx] = step_value(A, B_total);
                    dbg_row("aa", ai, aj, A, curvature.aa(idx), B_total, kappa[o3 + idx]);
                    ++idx;
                }
            }
        }

        if (report::reports(4)) {
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

    /// kappa^opt = sum_i c_i cum_kappa_i (pruned slots c_i = 0); empty when the core does not extrapolate.
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

    /// Ca_out = Ca_ref exp(K) (no alias), K = [[A, B], [-B^T, 0]] over O = [0, no): exactly
    /// exp(K) = I + U (exp(T) - I) U^T, T = U^T K U (Pade-4), U = [[I, 0], [0, Qb]], Qb = basis of range(B^T).
    void apply_rotation(const SharedMatrix& Ca_ref, const SharedMatrix& Ca_out,
                        const std::vector<double>& kappa,
                        const std::vector<int>& active_mo_indices,
                        int Ncore, int nmo) const {
        const int n_act = static_cast<int>(active_mo_indices.size());
        const int no = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;  // first_virt
        const int nv = nmo - no;

        auto K = std::make_shared<Matrix>("K_cum (orbital DIIS)", no, nmo);
        unpack_kappa_to_K(kappa, K, active_mo_indices, Ncore, nmo);
        double** Kp = K->pointer(0);

        // Qb = basis of range(B^T), B = K[0:no, no:nmo], by pivoted QR (rank m at relative |R_kk| tolerance);
        // Q column-major nv x no.
        std::vector<double> Q(static_cast<size_t>(nv) * no);
        for (int o = 0; o < no; ++o)
            for (int v = 0; v < nv; ++v) Q[static_cast<size_t>(o) * nv + v] = Kp[o][no + v];

        int m = 0;
        if (nv > 0) {
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

        // Ca_out = Ca_ref + W M
        Ca_out->copy(Ca_ref);
        C_DGEMM('N', 'N', nso, nmo, r, 1.0, W.data(), r, M.data(), nmo, 1.0,
                Ca_out->pointer(0)[0], nmo);
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
