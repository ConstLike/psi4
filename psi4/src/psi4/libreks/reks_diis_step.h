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

#ifndef REKS_DIIS_STEP_H
#define REKS_DIIS_STEP_H

#include "reks_cassette.h"
#include "reks_diis_config.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace psi {
namespace reks {

/// GVB-DIIS Newton step: h = -H_ii tables (CurvatureProvider), B = -h/2 + gamma (assemble_B),
/// denom = -max(|B|, B_floor) (DenominatorPolicy), step = A / denom (raw_step).

/// Index: AA i n_act - i(i+1)/2 + (j-i-1), CA c n_act + k, CV c n_virt + v', AV k n_virt + v', v' = v - first_virt;
/// occupations by local active slot, Fock diagonals by global MO.
enum class Block { AA, CA, CV, AV };

/// Pure-curvature inputs: excludes the additive gamma correction added in assemble_B.
struct CurvatureInputs {
    const std::vector<std::vector<double>>* F_MO_diag_a_L = nullptr;
    const std::vector<std::vector<double>>* F_MO_diag_b_L = nullptr;
    const std::vector<double>* C_L = nullptr;
    const studio::Cassette* sa_cassette = nullptr;
    const std::vector<int>* active_mo_indices = nullptr;
    int Ncore = 0;
    int nmo = 0;
    int n_act = 0;
    int first_virt = 0;
    int n_virt = 0;
};

/// prepare() batch-fills all four block tables once per iteration; accessors are
/// plain O(1) table reads.
class CurvatureProvider {
   public:
    virtual ~CurvatureProvider() = default;

    /// Deterministic; must leave the four tables sized for the block layouts.
    virtual void prepare(const CurvatureInputs& in) = 0;

    double aa(int idx) const { return h_aa_[idx]; }
    double ca(int idx) const { return h_ca_[idx]; }
    double cv(int idx) const { return h_cv_[idx]; }
    double av(int idx) const { return h_av_[idx]; }

   protected:
    std::vector<double> h_aa_, h_ca_, h_cv_, h_av_;
};

/// h_pq = 2 sum_L C_L sum_s w_s (F_s[p] - F_s[q]), w_s = n_i - n_j (aa), 1 - n_k (ca), 1 (cv), n_k (av), via
/// T = 2 sum_L C_L (Fa + Fb), Q_k with w = 1 - n_k, R_k = T - Q_k: h_cv = T[c] - T[v], h_ca = Q_k[c] - Q_k[a_k].
class Composite1eProvider final : public CurvatureProvider {
   public:
    void prepare(const CurvatureInputs& in) override {
        const int n_act = in.n_act;
        const int Ncore = in.Ncore;
        const int nmo = in.nmo;
        const int first_virt = in.first_virt;
        const int n_virt = in.n_virt;
        const auto& C_L = *in.C_L;
        const auto& sa_cassette = *in.sa_cassette;
        const auto& active_mo_indices = *in.active_mo_indices;

        h_aa_.assign(n_act * (n_act - 1) / 2, 0.0);
        h_ca_.assign(Ncore * n_act, 0.0);
        h_cv_.assign(Ncore * n_virt, 0.0);
        h_av_.assign(n_act * n_virt, 0.0);

        // T over all MOs; Q_k and R_k over the MOs their block reads (core + its own active
        // orbital for Q, virtuals + its own active orbital for R).
        const bool need_T = (Ncore > 0 && n_virt > 0);
        const bool need_Q = (Ncore > 0 && n_act > 0);
        const bool need_R = (n_virt > 0 && n_act > 0);
        std::vector<double> T(need_T ? nmo : 0, 0.0);
        std::vector<double> Q(need_Q ? static_cast<size_t>(n_act) * nmo : 0, 0.0);
        std::vector<double> R(need_R ? static_cast<size_t>(n_act) * nmo : 0, 0.0);

        for (int L : sa_cassette.microstates_active) {
            if (std::abs(C_L[L]) < 1e-14) continue;
            const double* Fdiag_a = (*in.F_MO_diag_a_L)[L].data();
            const double* Fdiag_b = (*in.F_MO_diag_b_L)[L].data();
            const auto& m = sa_cassette.microstate(L);
            const double w = 2.0 * C_L[L];

            int idx_aa = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    int ai = active_mo_indices[i];
                    int aj = active_mo_indices[j];
                    double d_na = m.alpha[i] - m.alpha[j];
                    double d_nb = m.beta[i] - m.beta[j];
                    double d_ea = Fdiag_a[ai] - Fdiag_a[aj];
                    double d_eb = Fdiag_b[ai] - Fdiag_b[aj];
                    h_aa_[idx_aa] += w * (d_na * d_ea + d_nb * d_eb);
                    ++idx_aa;
                }
            }

            // core (na=nb=1) minus virtual (na=nb=0) occupation diff is 1 on both spins,
            // absorbed into w.
            if (need_T)
                for (int p = 0; p < nmo; ++p) T[p] += w * (Fdiag_a[p] + Fdiag_b[p]);

            for (int k = 0; k < n_act; ++k) {
                const double na_k = m.alpha[k];
                const double nb_k = m.beta[k];
                if (need_Q) {
                    double* Qk = Q.data() + static_cast<size_t>(k) * nmo;
                    const double qa = w * (1.0 - na_k);
                    const double qb = w * (1.0 - nb_k);
                    for (int c = 0; c < Ncore; ++c) Qk[c] += qa * Fdiag_a[c] + qb * Fdiag_b[c];
                    const int a = active_mo_indices[k];
                    Qk[a] += qa * Fdiag_a[a] + qb * Fdiag_b[a];
                }
                if (need_R) {
                    double* Rk = R.data() + static_cast<size_t>(k) * nmo;
                    const double ra = w * na_k;
                    const double rb = w * nb_k;
                    for (int v = first_virt; v < nmo; ++v)
                        Rk[v] += ra * Fdiag_a[v] + rb * Fdiag_b[v];
                    const int a = active_mo_indices[k];
                    Rk[a] += ra * Fdiag_a[a] + rb * Fdiag_b[a];
                }
            }
        }

        for (int c = 0; c < Ncore; ++c)
            for (int k = 0; k < n_act; ++k) {
                const double* Qk = Q.data() + static_cast<size_t>(k) * nmo;
                h_ca_[c * n_act + k] = Qk[c] - Qk[active_mo_indices[k]];
            }

        for (int c = 0; c < Ncore; ++c)
            for (int v = first_virt; v < nmo; ++v)
                h_cv_[c * n_virt + (v - first_virt)] = T[c] - T[v];

        for (int k = 0; k < n_act; ++k) {
            const double* Rk = R.data() + static_cast<size_t>(k) * nmo;
            const double Ra = Rk[active_mo_indices[k]];
            for (int v = first_virt; v < nmo; ++v)
                h_av_[k * n_virt + (v - first_virt)] = Ra - Rk[v];
        }
    }
};

/// Additive B corrections. gamma_aa/gamma_ca may be null (no correction supplied).
struct BAssemblyInputs {
    const std::vector<double>* gamma_aa = nullptr;
    const std::vector<double>* gamma_ca = nullptr;
};

/// B = -h/2 plus the per-pair gamma correction; gamma enters unscaled.
inline double assemble_B(Block b, int idx, double h_diag, const BAssemblyInputs& corr) {
    const double B_1e = -h_diag / 2.0;
    switch (b) {
        case Block::AA: {
            const auto* g = corr.gamma_aa;
            double gamma = (g && idx < static_cast<int>(g->size())) ? (*g)[idx] : 0.0;
            return B_1e + gamma;
        }
        case Block::CA: {
            const auto* g = corr.gamma_ca;
            double gamma = (g && idx < static_cast<int>(g->size())) ? (*g)[idx] : 0.0;
            return B_1e + gamma;
        }
        default:
            return B_1e;
    }
}

/// denom = -max(|B|, B_floor): descent for either sign of B, |denom| >= B_floor.
class DenominatorPolicy {
   public:
    DenominatorPolicy() = default;
    explicit DenominatorPolicy(const DiisConfig& cfg) : b_floor_(cfg.B_floor) {}

    double shape(double B) const { return -std::max(std::abs(B), b_floor_); }
    /// True when shape() returns -b_floor_ instead of -|B| for this B.
    bool clamped(double B) const { return std::abs(B) < b_floor_; }

   private:
    double b_floor_ = DiisConfig{}.B_floor;
};

/// Newton step for one rotation.
inline double raw_step(double A, double denom) { return A / denom; }

/// Per-block Frobenius norms of the orbital gradient.
struct GradientBlocks {
    double aa = 0.0;
    double ca = 0.0;
    double cv = 0.0;
    double av = 0.0;
};

/// |g_orb| over aa, ca, cv, av, g[p,q] = -2 (F_gen[p,q] - F_gen[q,p]) (virtual rows zero); frozen[i<j] drops
/// an aa pair, blocks (optional) gets the per-block norms.
inline double orbital_gradient_norm(double* const* Fg, const std::vector<int>& active_mo_indices,
                                    int Ncore, int nmo, const std::vector<char>& frozen,
                                    GradientBlocks* blocks) {
    const int n_act = static_cast<int>(active_mo_indices.size());
    const int first_virt = active_mo_indices.back() + 1;
    double norm_sq = 0.0;
    double aa_sq = 0.0, ca_sq = 0.0, cv_sq = 0.0, av_sq = 0.0;

    int rot_idx = 0;
    for (int i = 0; i < n_act; ++i)
        for (int j = i + 1; j < n_act; ++j) {
            const bool is_frozen = rot_idx < static_cast<int>(frozen.size()) && frozen[rot_idx] != 0;
            ++rot_idx;
            if (is_frozen) continue;
            const int ai = active_mo_indices[i], aj = active_mo_indices[j];
            const double g = -2.0 * (Fg[ai][aj] - Fg[aj][ai]);
            norm_sq += g * g;
            aa_sq += g * g;
        }
    for (int c = 0; c < Ncore; ++c)
        for (int ai : active_mo_indices) {
            const double g = -2.0 * (Fg[c][ai] - Fg[ai][c]);
            norm_sq += g * g;
            ca_sq += g * g;
        }
    for (int c = 0; c < Ncore; ++c)
        for (int v = first_virt; v < nmo; ++v) {
            const double g = -2.0 * Fg[c][v];
            norm_sq += g * g;
            cv_sq += g * g;
        }
    for (int ai : active_mo_indices)
        for (int v = first_virt; v < nmo; ++v) {
            const double g = -2.0 * Fg[ai][v];
            norm_sq += g * g;
            av_sq += g * g;
        }
    if (blocks) {
        blocks->aa = std::sqrt(aa_sq);
        blocks->ca = std::sqrt(ca_sq);
        blocks->cv = std::sqrt(cv_sq);
        blocks->av = std::sqrt(av_sq);
    }
    return std::sqrt(norm_sq);
}

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_STEP_H
