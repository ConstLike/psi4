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

#ifndef REKS_GVB_DIIS_H
#define REKS_GVB_DIIS_H

#include "reks_report_level.h"
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_diis_config.h"
#include "reks_diis_core.h"
#include "reks_diis_step.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <vector>
#include <memory>
#include <cmath>

namespace psi {
namespace reks {

/// @class GVBDIISEngine
/// @brief GVB-DIIS convergence accelerator for REKS SCF.
///
/// Implements the multishell Fock operator F^g and DIIS extrapolation
/// from Muller, R. P. et al. J. Chem. Phys. 1994, 100, 1226.
///
/// Each active orbital is its own shell. Core = 1 shell, virtual = 1 shell.
class GVBDIISEngine {
   public:
    GVBDIISEngine() : GVBDIISEngine(DiisConfig{}) {}
    explicit GVBDIISEngine(DiisConfig cfg)
        : config_(cfg),
          core_(cfg.max_vectors),
          chain_(assemble_conditioning_chain(cfg)) {}

    /// Build multishell Fock F^g in MO basis.
    ///
    /// Diagonal: F^g[i][i] = i * level_shift (artificial, 0-based MO index).
    /// Same-shell off-diag: 0.
    /// Off-diag pair (i,j), Eq. 3.15b': F^g[i][j] = (A_ij / denom_ij) * gap_ij,
    ///   A_ij     = F_gen[j][i] - F_gen[i][j]
    ///   B_ij     = B^{1e}_ij + gamma_ij + ipr_ij     active-active
    ///   B_ij     = B^{1e}_ij + gamma_ij              core-active
    ///   B_ij     = B^{1e}_ij                         core-virtual, active-virtual
    ///   denom_ij = -max(|B_ij|, B_floor)
    ///   gap_ij   = (j - i) * level_shift
    /// B^{1e} = -H_diag/2, H_diag from the per-microstate Fock diagonals throughout;
    /// ipr_ij = -ipr_diag_aa/2 (same -H/2 convention); gamma_ij = gamma_aa / gamma_ca as passed.
    ///
    /// @param F_gen           Asymmetric generalized Fock (nmo x nmo)
    /// @param F_MO_diag_a_L   Per-microstate alpha MO Fock diagonal, [global L][nmo]
    /// @param F_MO_diag_b_L   Per-microstate beta  MO Fock diagonal, [global L][nmo]
    /// @param C_L             SA microstate weights, indexed by global L
    /// @param sa_cassette     Supplies the microstate work list and their occupations
    /// @param gamma_aa        Two-electron gamma, active-active (empty = off); i<j pair index
    /// @param gamma_ca        Two-electron gamma, core-active (empty = off); [c*n_act + k]
    /// @param iteration       SCF iteration number for tagging output
    /// @param ipr_diag_aa     IPR Hessian diagonal, active-active; i<j pair index. Always
    ///                        length n_rot from the host, zero-filled when IPR is off.
    /// @return F^g in MO basis (nmo x nmo, symmetric)
    SharedMatrix build_multishell_fock(
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
        const std::vector<double>& ipr_diag_aa = {}) const {

        auto F_g = std::make_shared<Matrix>("F^g MO", nmo, nmo);
        double** Fg = F_g->pointer(0);
        double** Fgen = F_gen->pointer(0);

        int n_act = static_cast<int>(active_mo_indices.size());

        double ls = config_.level_shift;
        for (int i = 0; i < nmo; ++i) {
            Fg[i][i] = static_cast<double>(i) * ls;
        }

        int first_virt = (n_act > 0) ? active_mo_indices.back() + 1 : Ncore;
        int n_virt = nmo - first_virt;

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

        Composite1eProvider curvature;
        curvature.prepare(cin);

        BAssemblyInputs corr;
        corr.gamma_aa = &gamma_aa;
        corr.gamma_ca = &gamma_ca;
        corr.ipr_diag_aa = &ipr_diag_aa;

        const DenominatorPolicy denom_policy(config_);

        // Active-active off-diag (Eq. 3.15b').
        {
            if (report::reports(4)) {
                outfile->Printf("  [Fg-DBG] iter=%d active-active rotations:\n", iteration);
                outfile->Printf("  %4s %4s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %s\n",
                                "i", "j", "F_gen[j,i]", "F_gen[i,j]", "A_ij",
                                "B_1e", "gamma", "B_total", "denom", "clamped?");
            }
            int idx = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    int ai = active_mo_indices[i];
                    int aj = active_mo_indices[j];

                    double A = Fgen[aj][ai] - Fgen[ai][aj];

                    double B_1e = -curvature.aa(idx) / 2.0;
                    double B = assemble_B(Block::AA, idx, curvature.aa(idx), corr);

                    double denom = denom_policy.shape(B);
                    bool clamped = denom_policy.clamped(B);

                    double gap = static_cast<double>(aj - ai) * ls;
                    Fg[ai][aj] = raw_step(A, denom) * gap;
                    Fg[aj][ai] = Fg[ai][aj];

                    if (report::reports(4)) {
                        double gamma = (!gamma_aa.empty() && idx < static_cast<int>(gamma_aa.size()))
                                           ? gamma_aa[idx] : 0.0;
                        outfile->Printf("  %4d %4d  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %s\n",
                                        ai, aj, Fgen[aj][ai], Fgen[ai][aj], A,
                                        B_1e, gamma, B, denom, clamped ? "YES" : "no");
                        outfile->Printf("           gap=%.4f  Fg[%d,%d]=%.6e\n", gap, ai, aj, Fg[ai][aj]);
                    }

                    ++idx;
                }
            }
        }

        if (report::reports(5)) {
            outfile->Printf("  [Fg-DBG] iter=%d core-active rotations:\n", iteration);
        }
        // Core-active off-diag (Eq. 3.15b').
        for (int c = 0; c < Ncore; ++c) {
            for (int k = 0; k < n_act; ++k) {
                int a = active_mo_indices[k];
                int idx = c * n_act + k;

                double A = Fgen[a][c] - Fgen[c][a];

                double B_1e = -curvature.ca(idx) / 2.0;
                double B = assemble_B(Block::CA, idx, curvature.ca(idx), corr);

                double denom = denom_policy.shape(B);
                bool clamped = denom_policy.clamped(B);

                double gap = static_cast<double>(a - c) * ls;
                Fg[c][a] = raw_step(A, denom) * gap;
                Fg[a][c] = Fg[c][a];

                if (report::reports(5)) {
                    double gamma = (!gamma_ca.empty()) ? gamma_ca[idx] : 0.0;
                    outfile->Printf("    c=%d a=%d: A=%.6e B_1e=%.6e gamma=%.6e B=%.6e denom=%.6e Fg=%.6e %s\n",
                                    c, a, A, B_1e, gamma, B, denom, Fg[c][a], clamped ? "CLAMPED" : "");
                }
            }
        }

        // Core-virtual off-diag (Eq. 3.15b').
        {
            if (report::reports(5)) {
                outfile->Printf("  [Fg-DBG] iter=%d core-virtual rotations:\n", iteration);
            }
            int idx = 0;
            for (int c = 0; c < Ncore; ++c) {
                for (int v = first_virt; v < nmo; ++v) {
                    // A_cv = F_gen[v][c] - F_gen[c][v]; F_gen[v][c] = 0 for virtual rows
                    double A = -Fgen[c][v];

                    double B_1e = -curvature.cv(idx) / 2.0;
                    double B = assemble_B(Block::CV, idx, curvature.cv(idx), corr);

                    double denom = denom_policy.shape(B);
                    bool clamped = denom_policy.clamped(B);

                    double gap = static_cast<double>(v - c) * ls;
                    Fg[c][v] = raw_step(A, denom) * gap;
                    Fg[v][c] = Fg[c][v];

                    if (report::reports(5)) {
                        outfile->Printf("    c=%d v=%d: A=%.6e B_1e=%.6e denom=%.6e Fg=%.6e %s\n",
                                        c, v, A, B_1e, denom, Fg[c][v], clamped ? "CLAMPED" : "");
                    }
                    ++idx;
                }
            }
        }

        // Active-virtual off-diag (Eq. 3.15b').
        {
            if (report::reports(5)) {
                outfile->Printf("  [Fg-DBG] iter=%d active-virtual rotations:\n", iteration);
            }
            int idx = 0;
            for (int k = 0; k < n_act; ++k) {
                int a = active_mo_indices[k];
                for (int v = first_virt; v < nmo; ++v) {
                    // A_av = F_gen[v][a] - F_gen[a][v]; F_gen[v][a] = 0 for virtual rows
                    double A = -Fgen[a][v];

                    double B_1e = -curvature.av(idx) / 2.0;
                    double B = assemble_B(Block::AV, idx, curvature.av(idx), corr);

                    double denom = denom_policy.shape(B);
                    bool clamped = denom_policy.clamped(B);

                    double gap = static_cast<double>(v - a) * ls;
                    Fg[a][v] = raw_step(A, denom) * gap;
                    Fg[v][a] = Fg[a][v];

                    if (report::reports(5)) {
                        outfile->Printf("    a=%d v=%d: A=%.6e B_1e=%.6e denom=%.6e Fg=%.6e %s\n",
                                        a, v, A, B_1e, denom, Fg[a][v], clamped ? "CLAMPED" : "");
                    }
                    ++idx;
                }
            }
        }

        if (report::reports(4)) {
            outfile->Printf("  [Fg-DBG] iter=%d F_gen active block (asymmetric):\n", iteration);
            for (int i = 0; i < n_act; ++i) {
                int ai = active_mo_indices[i];
                outfile->Printf("    row %d:", ai);
                for (int j = 0; j < n_act; ++j) {
                    int aj = active_mo_indices[j];
                    outfile->Printf(" %12.6e", Fgen[ai][aj]);
                }
                outfile->Printf("\n");
            }

            outfile->Printf("  [Fg-DBG] iter=%d F^g active block (symmetric):\n", iteration);
            for (int i = 0; i < n_act; ++i) {
                int ai = active_mo_indices[i];
                outfile->Printf("    row %d:", ai);
                for (int j = 0; j < n_act; ++j) {
                    int aj = active_mo_indices[j];
                    outfile->Printf(" %12.6e", Fg[ai][aj]);
                }
                outfile->Printf("\n");
            }

            double fg_core_act = 0, fg_act_virt = 0, fg_core_virt = 0;
            for (int c = 0; c < Ncore; ++c)
                for (int k = 0; k < n_act; ++k) {
                    double v = Fg[c][active_mo_indices[k]];
                    fg_core_act += v * v;
                }
            for (int k = 0; k < n_act; ++k)
                for (int v = first_virt; v < nmo; ++v) {
                    double val = Fg[active_mo_indices[k]][v];
                    fg_act_virt += val * val;
                }
            for (int c = 0; c < Ncore; ++c)
                for (int v = first_virt; v < nmo; ++v) {
                    double val = Fg[c][v];
                    fg_core_virt += val * val;
                }
            outfile->Printf("  [Fg-DBG] iter=%d ||Fg|| by block: core-act=%.3e act-virt=%.3e core-virt=%.3e\n",
                            iteration, std::sqrt(fg_core_act), std::sqrt(fg_act_virt), std::sqrt(fg_core_virt));
        }

        return F_g;
    }

    /// Flatten STRICT upper triangle (i < j, NO diagonal) of symmetric matrix.
    static std::vector<double> flatten_upper_triangle(SharedMatrix M, int nmo) {
        int size = nmo * (nmo - 1) / 2;
        std::vector<double> v(size);
        double** Mp = M->pointer(0);
        int idx = 0;
        for (int i = 0; i < nmo; ++i) {
            for (int j = i + 1; j < nmo; ++j) {
                v[idx++] = Mp[i][j];
            }
        }
        return v;
    }

    /// Store F^g (AO) payload with its error vector; the payload mirrors the core ring on the
    /// slot the store wrote (append, or in-place overwrite once the ring is full).
    /// Takes ownership of F_g_AO and of error_AO; caller must not write through either afterward.
    void store(SharedMatrix F_g_AO, std::vector<double> error_AO) {
        DiisCore::StoreResult sr = core_.store(std::move(error_AO));
        if (sr.wrote == static_cast<int>(F_AO_.size()))
            F_AO_.push_back(std::move(F_g_AO));
        else
            F_AO_[sr.wrote] = std::move(F_g_AO);
    }

    /// Removes the entry written by the last store(); no-op if no store is tracked.
    void drop_newest() {
        int dropped = core_.drop_newest();
        if (dropped >= 0)
            F_AO_.erase(F_AO_.begin() + dropped);
    }

    /// DIIS extrapolation over the stored F^g (AO) payloads:
    ///   F^opt = sum_i c_i F^g_AO_i,  c slot-indexed (slots the chain pruned carry c_i = 0).
    /// Returns nullptr when the core declines to extrapolate: < 2 vectors, < 2 survivors left by
    /// the conditioning chain, or the coeff-norm guard still trips after the ladder floor.
    SharedMatrix extrapolate() {
        auto combine = [this](const std::vector<double>& c) -> SharedMatrix {
            const int n = core_.count();
            auto F_opt = std::make_shared<Matrix>("F^opt AO", F_AO_[0]->rowspi(),
                                                  F_AO_[0]->colspi());
            for (int i = 0; i < n; ++i)
                if (c[i] != 0.0) F_opt->axpy(c[i], F_AO_[i]);
            return F_opt;
        };
        return core_.extrapolate(chain_, combine);
    }

    /// Clear DIIS history and the payload ring.
    void reset() {
        core_.reset();
        F_AO_.clear();
    }

    int count() const { return core_.count(); }
    double last_max_error() const { return core_.last_max_error(); }
    double last_c_norm_sq() const { return core_.last_c_norm_sq(); }

    /// DIIS coefficients c_i from the last extrapolate() call, slot-indexed.
    std::vector<double> last_coefficients() const { return core_.last_coefficients(); }

   private:
    DiisConfig config_;
    DiisCore core_;
    ConditioningChain chain_;
    std::vector<SharedMatrix> F_AO_;       ///< F^g (AO) payload, keyed by core slot
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_GVB_DIIS_H
