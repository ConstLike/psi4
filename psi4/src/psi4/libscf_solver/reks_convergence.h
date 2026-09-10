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

#ifndef REKS_CONVERGENCE_H
#define REKS_CONVERGENCE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

namespace psi {
namespace reks {

/// Staircase level-shift of F_MO diag + MOM-style reorder of U columns (max |U[i][j]| picks slot).
class OrbitalGuard {
   public:
    /// Active GVB pair found inverted in pre_diag (F_p > F_q with p < q).
    struct InvertedPair {
        int p;
        int q;
        double F_p;
        double F_q;
    };

    const std::vector<InvertedPair>& inverted_pairs() const { return last_inverted_pairs_; }

    /// Mixing diagnostic: max_active_mixing = largest off-diag for any active row.
    struct MixingResult {
        double max_active_mixing = 0.0;
        int max_mixing_target = -1;
        int n_swaps = 0;
    };

    /// Apply staircase shift to F_MO active/virt diag in place; record the per-position shifts applied.
    inline double pre_diag(double** F_MO, int N, int Ncore, const std::vector<int>& active_indices,
                           const std::vector<std::pair<int, int>>& pairs, int iteration,
                           double user_shift = 0.0) {
        double shift = compute_adaptive_shift(iteration, user_shift);

        if (shift <= 0.0) return 0.0;

        int n_active = static_cast<int>(active_indices.size());
        int first_virt = active_indices.back() + 1;

        // Active staircase: slot k gets (k+1)*shift.
        last_all_shifts_.assign(N, 0.0);
        for (int k = 0; k < n_active; ++k) last_all_shifts_[active_indices[k]] = (k + 1) * shift;

        // GVB-pair inversion fix: if F_p > F_q within a pair, swap their shifts so the
        // more-bonding orbital still gets the smaller penalty after diag-sort.
        last_inverted_pairs_.clear();
        for (const auto& pr : pairs) {
            int p = pr.first;
            int q = pr.second;
            int mo_p = active_indices[p];
            int mo_q = active_indices[q];
            double F_p = F_MO[mo_p][mo_p];
            double F_q = F_MO[mo_q][mo_q];
            if (F_p > F_q) {
                std::swap(last_all_shifts_[mo_p], last_all_shifts_[mo_q]);
                last_inverted_pairs_.push_back({p, q, F_p, F_q});
            }
        }

        // Virtual block: flat shift one step above the highest active slot.
        double virt_shift = (n_active + 1) * shift;
        for (int v = first_virt; v < N; ++v) last_all_shifts_[v] = virt_shift;

        for (int i = 0; i < N; ++i) F_MO[i][i] += last_all_shifts_[i];

        return shift;
    }

    /// Phase-fixes U column signs, then MOM-reorders slots to max overlap in place; returns mixing diagnostics.
    inline MixingResult post_diag(double** U, double* eps, int N, int Ncore, const std::vector<int>& active_indices) {
        // Phase fix: make each diagonal positive so signs are stable across iters.
        for (int j = 0; j < N; ++j) {
            if (U[j][j] < 0.0) {
                for (int i = 0; i < N; ++i) {
                    U[i][j] = -U[i][j];
                }
            }
        }

        // MOM: for each slot i (core+active), assign column j with max |U[i][j]| (max overlap).
        // Gilbert, A. T. B.; Besley, N. A.; Gill, P. M. W. J. Phys. Chem. A 2008, 112, 13164.
        MixingResult result;
        int last_reorder = active_indices.empty() ? Ncore - 1 : active_indices.back();
        for (int i = 0; i <= last_reorder; ++i) {
            int best_j = i;
            double best_val = std::abs(U[i][i]);
            for (int j = i + 1; j < N; ++j) {
                double val = std::abs(U[i][j]);
                if (val > best_val) {
                    best_val = val;
                    best_j = j;
                }
            }
            if (best_j != i) {
                for (int k = 0; k < N; ++k) {
                    double tmp = U[k][i];
                    U[k][i] = U[k][best_j];
                    U[k][best_j] = tmp;
                }
                std::swap(eps[i], eps[best_j]);
                // last_all_shifts_ is by POSITION (burned into F_MO[i][i] in pre_diag);
                // NOT swapped here, so shift[i] stays aligned with the value F_MO[i][i] saw.

                ++result.n_swaps;
            }
            if (U[i][i] < 0.0) {
                for (int k = 0; k < N; ++k) U[k][i] = -U[k][i];
            }
        }

        // max_active_mixing = max_{a in active_indices} max_{j != a} |U[a][j]|
        if (active_indices.empty()) return result;
        for (int a : active_indices) {
            double max_mix = 0.0;
            int max_j = -1;
            for (int j = 0; j < N; ++j) {
                if (j == a) continue;
                if (std::abs(U[a][j]) > max_mix) {
                    max_mix = std::abs(U[a][j]);
                    max_j = j;
                }
            }
            if (max_mix > result.max_active_mixing) {
                result.max_active_mixing = max_mix;
                result.max_mixing_target = max_j;
            }
        }

        prev_max_mix_ = result.max_active_mixing;
        prev_n_swaps_ = result.n_swaps;

        return result;
    }

    /// Per-position shifts burned into F_MO by last pre_diag (post pair-swap, virt staircase).
    /// Valid until the next pre_diag.
    const std::vector<double>& applied_shifts() const { return last_all_shifts_; }

    /// LEVEL_SHIFT_ADAPT off: hold staircase at user_shift (no decay/inflate).
    void set_adaptive(bool enabled) { adaptive_enabled_ = enabled; }

   private:
    /// Adaptive shift: decay when stable (no swaps, low mix), inflate on scrambling.
    /// Capped at user_shift, floored at 0.01. iter<=3 holds at user_shift.
    double compute_adaptive_shift(int iter, double user_shift) {
        if (user_shift <= 0.0) {
            adaptive_shift_ = 0.0;
            adaptive_initialized_ = false;
            return 0.0;
        }

        if (!adaptive_enabled_) {
            adaptive_shift_ = user_shift;
            adaptive_initialized_ = true;
            return user_shift;
        }

        if (iter <= 3 || !adaptive_initialized_) {
            adaptive_shift_ = user_shift;
            adaptive_initialized_ = true;
            return adaptive_shift_;
        }

        bool stable = (prev_n_swaps_ == 0) && (prev_max_mix_ < 0.15);
        bool scrambled = (prev_n_swaps_ > 0) || (prev_max_mix_ > 0.30);

        if (stable) {
            double decay;
            if (prev_max_mix_ < 0.001)
                decay = 0.5;
            else if (prev_max_mix_ < 0.01)
                decay = 0.7;
            else if (prev_max_mix_ < 0.05)
                decay = 0.85;
            else
                decay = 0.95;

            adaptive_shift_ *= decay;
        } else if (scrambled) {
            adaptive_shift_ = std::min(adaptive_shift_ * 1.5, user_shift);
        }

        constexpr double shift_floor = 0.01;
        adaptive_shift_ = std::max(adaptive_shift_, shift_floor);
        adaptive_shift_ = std::min(adaptive_shift_, user_shift);

        return adaptive_shift_;
    }

    std::vector<double> last_all_shifts_;
    std::vector<InvertedPair> last_inverted_pairs_;

    double adaptive_shift_ = 0.0;
    bool adaptive_initialized_ = false;
    bool adaptive_enabled_ = true;
    double prev_max_mix_ = 1.0;
    int prev_n_swaps_ = 0;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_CONVERGENCE_H
