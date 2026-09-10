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

#ifndef REKS_DIIS_OUTCOME_H
#define REKS_DIIS_OUTCOME_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace psi {
namespace reks {

/// Outcome-state adjudicator: the energy-band departure verdicts (basin transition B vs
/// localization signature S) and the restart budget. E_sa and fon_disp (FON displacement from
/// the best iterate) enter through observe(); E_best and g_best are the on_capture() mirrors;
/// tau_E, theta_f and band_B come from configure().
///
///   band = max(tau_E, g_best^2 / (2 band_B))   [g^2/(2B): quadratic-model energy at curvature B]
///   B    = (E_sa <= E_best - band) AND (fon_disp >= theta_f)
///   S    = (E_sa <  E_best - band) AND (fon_disp <  theta_f)
///          AND (g_best^2 / (2 band_B) <= tau_E)          [stationary-grade best]
/// Both are false while best is unset.
class OutcomeGuard {
   public:
    void configure(double tau_E, double theta_f, double band_B, int latch_grace,
                   int restart_budget, double progress_eps, double landing_rtol,
                   int nonmonotone_m) {
        tau_e_ = tau_E;
        theta_f_ = theta_f;
        band_b_ = band_B;
        latch_grace_ = latch_grace;
        restart_budget_ = restart_budget;
        progress_eps_ = progress_eps;
        landing_rtol_ = landing_rtol;
        nonmonotone_m_ = std::max(1, nonmonotone_m);
    }

    void arm() {
        best_e_sa_ = 0.0;
        best_gorb_mirror_ = -1.0;
        basin_latch_ = false;
        latch_clear_streak_ = 0;
        b_fire_count_ = 0;
        restarts_since_improvement_ = 0;
        landing_set_ = false;
        landing_e_ = 0.0;
        landing_g_ = -1.0;
        have_cur_ = false;
        cur_e_sa_ = nan_();
        e_ring_.clear();
        cur_fon_disp_ = 0.0;
    }

    /// Must be called every iteration, null steps included: the latch-clear streak counts
    /// observations, and energy_descent() reads the ring of preceding ones.
    void observe(double E_sa, double fon_disp) {
        if (have_cur_) {
            e_ring_.push_back(cur_e_sa_);
            if (static_cast<int>(e_ring_.size()) > nonmonotone_m_) e_ring_.erase(e_ring_.begin());
        }
        cur_e_sa_ = E_sa;
        cur_fon_disp_ = fon_disp;
        have_cur_ = true;
        if (basin_latch_) {
            if (b_raw()) {
                latch_clear_streak_ = 0;
            } else if (++latch_clear_streak_ >= latch_grace_) {
                basin_latch_ = false;
                latch_clear_streak_ = 0;
            }
        }
    }

    void on_capture(double E_sa, double gorb) {
        // Reset the restart streak only on a relative gorb drop >= progress_eps_.
        if (best_gorb_mirror_ >= 0.0 && gorb <= (1.0 - progress_eps_) * best_gorb_mirror_)
            restarts_since_improvement_ = 0;
        best_e_sa_ = E_sa;
        best_gorb_mirror_ = gorb;
    }

    /// Precondition: called only after a basin transition (B) has fired.
    void on_best_invalidated() {
        best_e_sa_ = 0.0;
        best_gorb_mirror_ = -1.0;
        basin_latch_ = true;
        latch_clear_streak_ = 0;
        ++b_fire_count_;
        restarts_since_improvement_ = 0;
        landing_set_ = false;
        landing_e_ = 0.0;
        landing_g_ = -1.0;
    }

    void on_monitor_restart() { ++restarts_since_improvement_; }

    /// True when (E_best, gorb_best) matches the last recorded landing within relative
    /// tolerance landing_rtol:
    ///   |E_best - landing_e| <= rtol * max(1, |E_best|)
    ///   |gorb_best - landing_g| <= rtol * gorb_best
    bool landing_repeat(double E_best, double gorb_best) const {
        return landing_set_ &&
               std::abs(E_best - landing_e_) <= landing_rtol_ * std::max(1.0, std::abs(E_best)) &&
               std::abs(gorb_best - landing_g_) <= landing_rtol_ * gorb_best;
    }

    void record_landing(double E_best, double gorb_best) {
        landing_set_ = true;
        landing_e_ = E_best;
        landing_g_ = gorb_best;
    }

    bool budget_exhausted() const { return restarts_since_improvement_ >= restart_budget_; }

    bool basin_transition_fire() const { return b_raw() && !basin_latch_; }

    bool localization_signature() const { return s_raw(); }

    /// True on localization (S), or on a basin transition (B) recurring while the latch is
    /// still set.
    bool refuse_capture() const { return s_raw() || (basin_latch_ && b_raw()); }

    int b_fire_count() const { return b_fire_count_; }
    int restarts_since_improvement() const { return restarts_since_improvement_; }

    /// E_SA fell by more than tau_E against the nonmonotone reference min(E_SA over the last M
    /// observations) (Grippo, L. et al. SIAM J. Numer. Anal. 1986, 23, 707), M = nonmonotone_m.
    /// M = 1 is descent against the previous observation, under which a limit cycle of amplitude
    /// > tau_E reads as descent on every second iteration.
    bool energy_descent() const {
        if (!have_cur_ || e_ring_.empty()) return false;
        const double ref = *std::min_element(e_ring_.begin(), e_ring_.end());
        return cur_e_sa_ <= ref - tau_e_;
    }

   private:
    static double nan_() { return std::numeric_limits<double>::quiet_NaN(); }

    bool best_set() const { return best_gorb_mirror_ >= 0.0; }

    double band() const {
        return std::max(tau_e_, best_gorb_mirror_ * best_gorb_mirror_ / (2.0 * band_b_));
    }

    bool b_raw() const {
        return have_cur_ && best_set() && cur_e_sa_ <= best_e_sa_ - band() &&
               cur_fon_disp_ >= theta_f_;
    }

    bool s_raw() const {
        return have_cur_ && best_set() && cur_e_sa_ < best_e_sa_ - band() &&
               cur_fon_disp_ < theta_f_ &&
               best_gorb_mirror_ * best_gorb_mirror_ / (2.0 * band_b_) <= tau_e_;
    }

    double tau_e_ = 1e-7;
    double theta_f_ = 0.1;
    double band_b_ = 0.1;
    int    latch_grace_ = 3;
    int    restart_budget_ = 2;
    double progress_eps_ = 1e-2;
    double landing_rtol_ = 1e-12;
    int    nonmonotone_m_ = 3;   ///< depth M of the energy_descent() nonmonotone reference

    double best_e_sa_ = 0.0;
    double best_gorb_mirror_ = -1.0;  ///< < 0 marks best unset; best_e_sa_ is meaningless then
    bool   basin_latch_ = false;
    int    latch_clear_streak_ = 0;
    int    b_fire_count_ = 0;
    int    restarts_since_improvement_ = 0;
    bool   landing_set_ = false;   ///< a rewind-landing fingerprint is recorded
    double landing_e_ = 0.0;       ///< best E_sa at the last recorded landing
    double landing_g_ = -1.0;      ///< best gorb at the last recorded landing

    bool   have_cur_ = false;
    double cur_e_sa_ = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> e_ring_;  ///< the last nonmonotone_m_ observations preceding cur_e_sa_
    double cur_fon_disp_ = 0.0;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_OUTCOME_H
