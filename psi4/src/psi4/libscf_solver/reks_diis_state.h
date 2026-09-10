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

#ifndef REKS_DIIS_STATE_H
#define REKS_DIIS_STATE_H

namespace psi {
namespace reks {

/// GVB-DIIS controller state: two orthogonal axes plus the counters every
/// landing/restart path re-bases.
///
/// Both axes vary independently: the solver can be ACTIVE while the rung is
/// not ARMED.
///
///   activation  inactive -> active        (one-way, at the burn-in handoff)
///   rung        ARMED -> RETAIN -> FALLBACK (quiescence escalation, monotone
///               within an activation; reset to ARMED by activation and by a
///               basin invalidation)
///
/// Sole owner of last_rollback_iter and gate_streak.
///
/// PLANE ORDER (load-bearing, do not reorder): evidence (OutcomeGuard::observe)
/// runs before adjudication; adjudication before accounting.
class ControllerState {
   public:
    enum class Rung { ARMED, RETAIN, FALLBACK };

    bool is_active() const { return active_; }

    /// Burn-in handoff: switches the accelerator on and re-bases the counters.
    void activate() {
        active_ = true;
        rung_ = Rung::ARMED;
        last_rollback_iter_ = kNoRollback;
        gate_streak_ = 0;
    }

    Rung rung() const { return rung_; }
    void set_rung(Rung r) { rung_ = r; }

    /// Idempotent within an iteration: a rewind and a co-fired restart in the
    /// same iteration count once.
    void on_landing_or_restart(int iteration) {
        last_rollback_iter_ = iteration;
        gate_streak_ = 0;
    }

    int last_rollback_iter() const { return last_rollback_iter_; }

    int gate_streak() const { return gate_streak_; }
    void set_gate_streak(int s) { gate_streak_ = s; }

   private:
    static constexpr int kNoRollback = -100;

    bool active_ = false;
    Rung rung_ = Rung::ARMED;
    int last_rollback_iter_ = kNoRollback;
    int gate_streak_ = 0;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_STATE_H
