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

#ifndef REKS_DIIS_CONTROLLER_H
#define REKS_DIIS_CONTROLLER_H

#include "reks_diis_config.h"
#include "reks_diis_monitors.h"
#include "reks_diis_adapters.h"
#include "reks_diis_state.h"

#include <memory>
#include <vector>

namespace psi {
namespace reks {

/// Composed decision for one iteration: rewind and restart may co-fire, rewind
/// executing first (Toth-Kelley 2015 safeguarded restart).
///
/// Fields past `reason` are an already-made decision: the caller must execute them,
/// not re-derive them.
struct ResolvedAction {
    bool        rewind = false;
    bool        restart = false;
    const char* reason = "";

    /// The rewind retreats to best (quiescence ladder) rather than to the rejected iterate.
    bool        landed_on_best = false;
    bool        retain_history = false;       ///< rewind keeps the DIIS subspace (RETAIN rung)
    bool        record_landing = false;       ///< fingerprint this landing
    bool        count_toward_budget = false;  ///< charge the restart budget
    /// The rewind landed on an unchanged best without a co-fired restart. Distinct
    /// from count_toward_budget, which a co-fired restart also sets.
    bool        repeat_landing = false;
    bool        verdict_trip = false;         ///< counts as a monitor-verdict rewind
    ControllerState::Rung rung_target = ControllerState::Rung::ARMED;
    bool        set_rung = false;             ///< rung_target is meaningful
};

/// Fields are fixed before apply_quiescence() runs, independent of its output: no
/// ordering cycle.
struct QuiesceInputs {
    ControllerState::Rung rung = ControllerState::Rung::ARMED;
    bool quiesced = false;            ///< restart budget exhausted
    bool landing_repeat = false;      ///< best matches the last recorded landing fingerprint
    bool best_below_current = false;  ///< a best exists and is strictly below the current frame
};

/// Owns the restart monitors; routes each iteration's solve through the engine adapter selected
/// by formulation (CFM vs ORBITAL).
class DiisController {
   public:
    void configure(const DiisConfig& cfg) {
        monitors_.clear();
        monitors_.add(std::make_unique<VerdictMonitor>(
            cfg.mon_verdict, cfg.verdict_rho, cfg.verdict_rewind_streak, cfg.mon_suppress_window,
            cfg.mon_nonmonotone_m, cfg.mon_progress_eps));
        monitors_.add(std::make_unique<CycleMonitor>(cfg.mon_cycle, cfg.cycle_window,
                                                     cfg.cycle_cap_streak, cfg.mon_progress_eps));
    }

    /// Non-owning: both adapters must outlive the controller.
    void bind(CfmAdapter* cfm, OrbAdapter* orb) {
        cfm_ = cfm;
        orb_ = orb;
        active_ = cfm;
    }

    void select(bool orbital) {
        active_ = orbital ? static_cast<IEngineAdapter*>(orb_) : static_cast<IEngineAdapter*>(cfm_);
    }

    IEngineAdapter* active() const { return active_; }

    /// RESTART_BASE_MAP if any structural trigger fires; the check order sets only the reason label.
    MonitorAction pre_step_directive(bool swap_detected, bool ipr_guard_fires,
                                     bool basin_transition) const {
        if (basin_transition) return {MonitorVerdict::RESTART_BASE_MAP, "basin-transition"};
        if (swap_detected) return {MonitorVerdict::RESTART_BASE_MAP, "orbital-swap"};
        if (ipr_guard_fires) return {MonitorVerdict::RESTART_BASE_MAP, "ipr-basin"};
        return {};
    }

    StepResult step(const StepContext& ctx) { return active_->run_step(ctx); }

    void arm() { monitors_.arm_all(); }

    /// reason keeps the first verdict fired. A monitor-driven restart re-anchors from the best
    /// iterate, so it co-fires a rewind-to-best when best is strictly below the current frame
    /// (Toth-Kelley 2015 safeguarded restart).
    ResolvedAction review(const TrajectoryView& t) {
        ResolvedAction r;
        for (const MonitorAction& a : monitors_.review(t)) {
            switch (a.verdict) {
                case MonitorVerdict::REWIND_TO_BEST:
                    r.rewind = true;
                    if (r.reason[0] == '\0') r.reason = a.reason;
                    break;
                case MonitorVerdict::RESTART_BASE_MAP:
                    r.restart = true;
                    if (r.reason[0] == '\0') r.reason = a.reason;
                    break;
                case MonitorVerdict::CONTINUE:
                    break;
            }
        }
        if (r.restart && t.best_available && t.best_gorb >= 0.0 && t.best_gorb < t.gorb)
            r.rewind = true;
        return r;
    }

    /// Quiescence post-transform. Consumes the resolved DECISION (never the executed
    /// result), so the fold stays acyclic.
    ///
    /// Priority, in order:
    ///   FALLBACK  monitor verdicts are discarded outright until convergence or MAXITER
    ///   quiesced  past the restart budget, a verdict escalates the ladder instead of
    ///             acting: RETAIN on a fresh landing, FALLBACK once a landing recurs on
    ///             an unchanged best
    ///   otherwise the verdict acts: rewind (fingerprinted; a repeat landing without a
    ///             co-fired restart is charged to the budget) then restart
    ///
    /// Structural directives (basin transition, orbital swap, IPR guard) are NOT part of
    /// this transform -- they stay live on every rung.
    ResolvedAction apply_quiescence(ResolvedAction act, const QuiesceInputs& q) const {
        ResolvedAction out;
        out.reason = act.reason;

        if (q.rung == ControllerState::Rung::FALLBACK) return out;

        if (q.quiesced && (act.rewind || act.restart)) {
            if (q.rung == ControllerState::Rung::RETAIN && q.landing_repeat) {
                out.set_rung = true;
                out.rung_target = ControllerState::Rung::FALLBACK;
                return out;
            }
            out.set_rung = true;
            out.rung_target = ControllerState::Rung::RETAIN;
            // When best is not strictly better nothing executes, yet the rung still
            // advances: the ladder measures escalation, not work done.
            if (q.best_below_current) {
                out.rewind = true;
                out.retain_history = true;
                out.record_landing = true;
                out.landed_on_best = true;
            }
            return out;
        }

        if (act.rewind) {
            out.rewind = true;
            out.record_landing = true;
            out.verdict_trip = true;
            // A rewind-only landing on the same best replays a deterministic transient:
            // charge it. A co-fired restart is charged below instead.
            out.repeat_landing = !act.restart && q.landing_repeat;
            out.count_toward_budget = out.repeat_landing;
        }
        if (act.restart) {
            out.restart = true;
            out.count_toward_budget = true;
        }
        return out;
    }

   private:
    MonitorList     monitors_;
    CfmAdapter*     cfm_ = nullptr;
    OrbAdapter*     orb_ = nullptr;
    IEngineAdapter* active_ = nullptr;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_CONTROLLER_H
