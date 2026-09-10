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

#ifndef REKS_DIIS_CONFIG_H
#define REKS_DIIS_CONFIG_H

namespace psi {
namespace reks {

/// Tunables for the REKS GVB-DIIS solver. Each field's source option,
/// when it has one, is noted in parentheses.
struct DiisConfig {
    int    max_vectors   = 10;      ///< ring buffer capacity (DIIS_MAX_VECS)

    double level_shift   = 1.0;     ///< F^g[i][i] = i * level_shift, 0-based MO index (REKS_GVB_LEVEL_SHIFT)
    double B_floor       = 0.1;     ///< Newton-step denominator floor: denom = -max(|B|, B_floor)

    // Conditioning chain.
    bool   cond_angle_filter = true;  ///< Pollock-Rebholz angle filter + Chupin stale-drop (REKS_DIIS_COND_ANGLE_FILTER)
    double angle_tol     = 0.1;       ///< c_s: min sin(angle to kept span) to keep a vector (REKS_DIIS_COND_ANGLE_TOL)
    double stale_delta   = 1e-4;      ///< Chupin stale-drop: drop error >= newest/sqrt(delta); 0 disables (REKS_DIIS_COND_STALE_DELTA)
    double kappa_bypass  = 1e9;       ///< decision cond(E) bypass boundary (raw-E SVD); above it -> engaged (REKS_DIIS_COND_KAPPA_BYPASS)
    double svd_rcond     = 1e-8;      ///< engaged SVD on F: drops sigma_k/sigma_max < svd_rcond (REKS_DIIS_COND_SVD_RCOND)
    double coeff_norm_max = 1e3;      ///< |c|^2 bound on every solve; trip -> escalation ladder (REKS_DIIS_COND_COEFF_NORM_MAX)
    bool   cond_tikhonov = true;      ///< Tikhonov ridge policy on the engaged SVD path (REKS_DIIS_COND_TIKHONOV)
    double tikhonov_scale = 0.0;      ///< engaged ridge magnitude; 0 -> delta = 0, no ridge (REKS_DIIS_COND_TIKHONOV_SCALE)

    // ORB step bounds.
    double max_step_kappa = 0.5;      ///< cap on |kappa^opt[i]|
    bool   two_level_saturation = true;  ///< per-rotation 2-level eigenvector-angle step bound

    // Restart monitors.
    bool   mon_verdict     = true;    ///< A9 post-hoc verdict (REKS_GVB_ROBUST_VERDICT)
    double verdict_rho     = 1.5;     ///< A9 trip ratio (Frobenius gorb): rewind if gorb > rho * max(gorb over last mon_nonmonotone_m iterates) (REKS_DIIS_MON_VERDICT_RHO)
    int    verdict_rewind_streak = 3; ///< consecutive A9 trips without sufficient decrease -> escalate-restart (REKS_DIIS_MON_VERDICT_STREAK)
    int    mon_suppress_window = 3;   ///< iterations after an executed rewind/restart/switch with no A9 verdict (REKS_DIIS_MON_SUPPRESS_WINDOW)
    int    mon_nonmonotone_m = 3;     ///< depth M of both nonmonotone references: max(gorb over last M iterates) and min(E_SA over last M observations) (REKS_DIIS_MON_NONMONOTONE_M)
    int    guard_latch_grace = 3;     ///< consecutive B-false iterations before the basin latch clears (REKS_DIIS_GUARD_LATCH_GRACE)
    double guard_band_B    = 0.1;     ///< curvature B in the adjudication band g_best^2 / (2 B) (REKS_DIIS_GUARD_BAND_B)
    double guard_landing_rtol = 1e-12;  ///< relative tolerance of the rewind-landing fingerprint match (REKS_DIIS_GUARD_LANDING_RTOL)
    double mon_progress_eps = 1e-2;   ///< eps in new <= (1 - eps) * old, a relative sufficient-decrease test (REKS_DIIS_MON_PROGRESS_EPS)
    double fon_branch_tol   = 0.1;    ///< theta_f: min FON displacement (vs best) of a branch transition (REKS_DIIS_FON_BRANCH_TOL)
    int    restart_budget   = 2;      ///< monitor restarts without sufficient improvement before quiesce (REKS_DIIS_RESTART_BUDGET)
    bool   mon_cycle       = true;    ///< no-progress / cap-streak restart monitor (REKS_DIIS_MON_CYCLE)
    int    cycle_window    = 12;      ///< no-progress window: iters without sufficient decrease before restart (REKS_DIIS_MON_CYCLE_WINDOW)
    int    cycle_cap_streak = 3;      ///< consecutive step-cap fires before restart (REKS_DIIS_MON_CYCLE_CAP_STREAK)
    double orb_base_map_damp = 0.6;   ///< ORB re-anchored base-map damping alpha (REKS_DIIS_ORB_BASE_MAP_DAMP)
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_CONFIG_H
