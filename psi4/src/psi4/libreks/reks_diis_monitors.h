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

#ifndef REKS_DIIS_MONITORS_H
#define REKS_DIIS_MONITORS_H

#include <algorithm>
#include <memory>
#include <vector>

namespace psi {
namespace reks {

enum class MonitorVerdict { CONTINUE, REWIND_TO_BEST, RESTART_BASE_MAP };

/// One monitor's decision; reason must have static storage duration (a string literal).
struct MonitorAction {
    MonitorVerdict verdict = MonitorVerdict::CONTINUE;
    const char*    reason  = "";
};

struct TrajectoryView {
    int    iteration = 0;
    double gorb = 0.0;          ///< ||g_orb||_F at the current iterate (result of the previous step)
    double best_gorb = -1.0;    ///< lowest gorb seen this SCF (-1 = unset)
    bool   prev_did_extrap = false;
    bool   best_available = false;
    bool   cap_fired = false;   ///< ORB step cap fired on the previous step (ORB-only signal)
    int    last_rollback_iter = -100;  ///< iteration of the last executed rewind/restart/switch
                                       ///< (-100 = none; far in the past, so no suppression)
    bool   energy_descent = false;  ///< E_SA fell beyond tau_E against min(E_SA over the last M)
};

/// Monitor port: review() decides on the trajectory alone -- it never touches orbitals, FONs or
/// history; a monitor may update only its own detection state.
class IRestartMonitor {
   public:
    virtual ~IRestartMonitor() = default;
    virtual MonitorAction review(const TrajectoryView& t) = 0;
    virtual void arm() {}                       ///< re-init this monitor's detection state
};

/// Extrapolated step with gorb > rho max(gorb, last nonmonotone_m) -> REWIND_TO_BEST; streak_max such trips without
/// gorb <= (1 - progress_eps) best_gorb -> RESTART_BASE_MAP; silent for suppress_window after an action.
/// Toth, A.; Kelley, C. T. SIAM J. Numer. Anal. 2015, 53, 805. Pulay, P. J. Comput. Chem. 1982, 3, 556.
/// Grippo, L. et al. SIAM J. Numer. Anal. 1986, 23, 707.
class VerdictMonitor final : public IRestartMonitor {
   public:
    VerdictMonitor(bool enabled, double rho, int streak_max, int suppress_window, int nonmonotone_m,
                   double progress_eps)
        : enabled_(enabled), rho_(rho), streak_max_(streak_max),
          suppress_window_(suppress_window), nonmonotone_m_(nonmonotone_m),
          progress_eps_(progress_eps) {}

    void arm() override {
        rewind_streak_ = 0;
        ring_.clear();
    }

    MonitorAction review(const TrajectoryView& t) override {
        if (t.best_available && t.gorb <= (1.0 - progress_eps_) * t.best_gorb) rewind_streak_ = 0;
        double ref = -1.0;
        for (double g : ring_) ref = std::max(ref, g);
        const bool trip = enabled_ && t.prev_did_extrap &&
                          t.iteration > t.last_rollback_iter + suppress_window_ && ref > 0.0 &&
                          t.best_available && t.gorb > rho_ * ref;
        push_ring(t.gorb);
        if (!trip) return {};
        if (++rewind_streak_ >= streak_max_) {
            rewind_streak_ = 0;
            return {MonitorVerdict::RESTART_BASE_MAP, "rewind-streak"};
        }
        return {MonitorVerdict::REWIND_TO_BEST, "gorb increase beyond rho"};
    }

   private:
    void push_ring(double g) {
        ring_.push_back(g);
        if (static_cast<int>(ring_.size()) > nonmonotone_m_) ring_.erase(ring_.begin());
    }

    bool   enabled_;
    double rho_;
    int    streak_max_;
    int    suppress_window_;
    int    nonmonotone_m_;
    double progress_eps_;
    int    rewind_streak_ = 0;
    std::vector<double> ring_;
};

/// RESTART_BASE_MAP after `window` iterations without nonmonotone E_SA descent or gorb <= (1 - progress_eps)
/// best_gorb, or after `cap_streak` ORB step-cap fires (max |kappa| > max_step_kappa).
/// Toth, A.; Kelley, C. T. SIAM J. Numer. Anal. 2015, 53, 805. Pratapa, P. P.; Suryanarayana, P. Chem. Phys.
/// Lett. 2015, 635, 69. Zhang, J. et al. SIAM J. Optim. 2020, 30, 3170.
class CycleMonitor final : public IRestartMonitor {
   public:
    CycleMonitor(bool enabled, int window, int cap_streak, double progress_eps)
        : enabled_(enabled), window_(window), cap_streak_max_(cap_streak),
          progress_eps_(progress_eps) {}

    void arm() override {
        no_progress_count_ = 0;
        cap_streak_ = 0;
    }

    MonitorAction review(const TrajectoryView& t) override {
        if (!enabled_) return {};

        cap_streak_ = t.cap_fired ? cap_streak_ + 1 : 0;
        if (cap_streak_ >= cap_streak_max_) {
            cap_streak_ = 0;
            return {MonitorVerdict::RESTART_BASE_MAP, "cap-streak"};
        }

        // A descending iterate resets no_progress_count_ but not cap_streak_: the two
        // triggers count independent events.
        const bool progress = t.energy_descent || !t.best_available ||
                              t.gorb <= (1.0 - progress_eps_) * t.best_gorb;
        no_progress_count_ = progress ? 0 : no_progress_count_ + 1;
        if (no_progress_count_ >= window_) {
            no_progress_count_ = 0;
            return {MonitorVerdict::RESTART_BASE_MAP, "no-progress"};
        }
        return {};
    }

   private:
    bool   enabled_;
    int    window_;
    int    cap_streak_max_;
    double progress_eps_;
    int    no_progress_count_ = 0;
    int    cap_streak_ = 0;
};

class MonitorList {
   public:
    void add(std::unique_ptr<IRestartMonitor> m) { monitors_.push_back(std::move(m)); }
    void clear() { monitors_.clear(); }

    std::vector<MonitorAction> review(const TrajectoryView& t) {
        std::vector<MonitorAction> out;
        out.reserve(monitors_.size());
        for (auto& m : monitors_) out.push_back(m->review(t));
        return out;
    }
    void arm_all() {
        for (auto& m : monitors_) m->arm();
    }

   private:
    std::vector<std::unique_ptr<IRestartMonitor>> monitors_;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_MONITORS_H
