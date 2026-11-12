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

#include "reks_gvb_driver.h"

#include "reks_diis_monitors.h"
#include "reks_diis_step.h"
#include "reks_report_level.h"
#include "reks_scf.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace psi {
namespace reks {

using report::reports;

void GvbDiisDriver::configure(const GvbDriverConfig& cfg) {
    const DiisConfig& dcfg = cfg.diis;
    gvb_diis_ = reks::GVBDIISEngine(dcfg);
    gvb_diis_.reset();
    gvb_diis_start_ = cfg.start_iteration;
    diis_orb_base_map_damp_ = dcfg.orb_base_map_damp;
    diis_controller_.configure(dcfg);
    diis_controller_.bind(&cfm_adapter_, &orb_adapter_);
    gvb_d_conv_ = cfg.d_convergence;
    gvb_e_conv_ = cfg.e_convergence;
    gvb_scf_maxiter_ = cfg.maxiter;
    gvb_gate_streak_req_ = cfg.gate_streak;
    cfm_activation_gorb_ = cfg.cfm_activation_gorb;
    outcome_guard_.configure(cfg.e_convergence, dcfg.fon_branch_tol,
                             dcfg.guard_band_B, dcfg.guard_latch_grace, dcfg.restart_budget,
                             dcfg.mon_progress_eps, dcfg.guard_landing_rtol,
                             dcfg.mon_nonmonotone_m);

    // ORBITAL (default; Ionova, I. V.; Carter, E. A. J. Chem. Phys. 1995, 102, 1251;
    // Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472) vs CFM (Muller, R. P. et al.
    // J. Chem. Phys. 1994, 100, 1226).
    diis_formulation_ = cfg.orbital_formulation ? DIISFormulation::ORBITAL
                                                : DIISFormulation::CFM;
    orbital_runtime_active_ = false;
    orbital_cap_fired_this_iter_ = false;
    if (diis_formulation_ == DIISFormulation::ORBITAL) {
        if (reports(4)) {
            outfile->Printf("  GVB-DIIS formulation: ORBITAL (Ionova & Carter 1995; Sethio et al. 2024) "
                            "with limit-cycle re-anchored restart\n");
        }
        orbital_diis_ = reks::OrbitalDIISEngine(dcfg);
        orbital_diis_.reset();
        orbital_diis_cum_kappa_.clear();
        Ca_ref_.reset();
    }
}

// Rejected step = null step (Conn, A. R.; Gould, N. I. M.; Toint, P. L. Trust-Region
// Methods; SIAM: Philadelphia, 2000).
void GvbDiisDriver::execute_landed_rewind(GvbFockContext& ctx, GvbReportWitness witness,
                                          bool retain_history) {
    ctx.Ca->copy(gvb_best_Ca_);
    ctx.Cb->copy(gvb_best_Ca_);
    ctx.eps->copy(*gvb_best_eps_);
    reks::studio::restore_fons(ctx.reel, ctx.sa_cassette, ctx.fon_blocks, gvb_best_fons_);
    reks::StepContext step_ctx;
    step_ctx.Ca = &ctx.Ca;
    step_ctx.Ca_ref = &Ca_ref_;
    step_ctx.cum_kappa = &orbital_diis_cum_kappa_;
    step_ctx.Ncore = ctx.Ncore;
    // Frame guard: a frame-relative (ORB) history survives only in its own anchor frame,
    // compared by frame identity (epoch) and evaluated before restore().
    bool retain = retain_history;
    if (retain && gvb_best_engine_memento_.Ca_ref &&
        gvb_best_engine_memento_.epoch != diis_controller_.active()->epoch())
        retain = false;
    diis_controller_.active()->restore(gvb_best_engine_memento_, step_ctx);
    if (retain)
        diis_controller_.active()->on_rewind_retain(step_ctx);
    else
        diis_controller_.active()->on_rewind(step_ctx);
    ctrl_state_.on_landing_or_restart(ctx.iteration);
    ++gvb_landing_count_;
    gvb_null_step_ = true;
    gvb_landed_iter_ = ctx.iteration;
    gvb_report_witness_ = witness;
}

void GvbDiisDriver::clear_best() {
    gvb_best_gorb_ = -1.0;
    gvb_best_E_ = 0.0;
    gvb_best_iter_ = -1;
    gvb_best_Ca_.reset();
    gvb_best_eps_.reset();
    gvb_best_fons_.clear();
    gvb_best_engine_memento_ = reks::EngineMemento{};
    diis_controller_.arm();
}

void GvbDiisDriver::reset_to_activation(GvbFockContext& ctx) {
    ctrl_state_.activate();
    cfm_adapter_.reset();
    clear_best();
    gvb_prev_did_extrap_ = false;
    gvb_verdict_trips_ = 0;
    gvb_landing_count_ = 0;
    gvb_null_step_ = false;
    gvb_report_witness_ = GvbReportWitness::CURRENT;
    gvb_report_E_override_armed_ = false;
    outcome_guard_.arm();

    if (diis_formulation_ == DIISFormulation::ORBITAL) {
        int nmo_act = ctx.Ca->colspi()[0];
        int n_act = static_cast<int>(ctx.active.size());
        int first_virt = (n_act > 0) ? ctx.active.back() + 1 : ctx.Ncore;
        int n_virt = nmo_act - first_virt;
        orb_adapter_.reset();
        orbital_diis_cum_kappa_.assign(reks::OrbAdapter::n_nonred(ctx.Ncore, n_act, n_virt), 0.0);
        // Anchor: Ca_ref = Ca, cum_kappa = 0.
        reks::StepContext anchor_ctx;
        anchor_ctx.Ca = &ctx.Ca;
        anchor_ctx.Ca_ref = &Ca_ref_;
        anchor_ctx.cum_kappa = &orbital_diis_cum_kappa_;
        anchor_ctx.Ncore = ctx.Ncore;
        orb_adapter_.anchor_here(anchor_ctx);
        orbital_runtime_active_ = true;
        orbital_cap_fired_this_iter_ = false;
    }
}

void GvbDiisDriver::restart(GvbFockContext& ctx) {
    if (!ctrl_state_.is_active()) return;
    reset_to_activation(ctx);
    if (reports(4))
        outfile->Printf("  [GVB-DIIS] iter=%d RESTART: energy functional changed, activation reset\n",
                        ctx.iteration);
}

void GvbDiisDriver::set_terminal_report(double gorb) {
    gvb_report_gorb_ = gorb;
    gvb_report_E_override_armed_ = false;
}

void GvbDiisDriver::report_pair(GvbReportWitness witness) {
    switch (witness) {
        case GvbReportWitness::CURRENT:
        case GvbReportWitness::REJECTED:
            gvb_report_gorb_ = gvb_diis_gorb_norm_;
            gvb_report_E_override_armed_ = false;
            break;
        case GvbReportWitness::BEST:
            gvb_report_gorb_ = gvb_best_gorb_;
            gvb_report_E_override_armed_ = true;
            break;
    }
}

bool GvbDiisDriver::on_fock(GvbFockContext& ctx) {
    // Activation: ORBITAL at once, CFM at gorb < gorb_thresh (or without TRAH); F_gen holds the generalized
    // Fock once the CFM gate has built it.
    SharedMatrix F_gen;
    bool f_gen_current = false;
    if (!ctrl_state_.is_active() && ctx.iteration >= gvb_diis_start_) {
        bool gate_passes = false;
        double gorb_check = 0.0;
        if (diis_formulation_ == DIISFormulation::ORBITAL) {
            gate_passes = true;
        } else {
            F_gen = ctx.build_F_gen();
            f_gen_current = true;
            gorb_check = reks::orbital_gradient_norm(F_gen->pointer(0), ctx.active, ctx.Ncore,
                                                     ctx.nmo, ctx.aa_frozen, nullptr);
            const double gorb_thresh = cfm_activation_gorb_;
            gate_passes = (gorb_check < gorb_thresh || !ctx.use_trah);
            if (!gate_passes && reports(4)) {
                outfile->Printf("  [GVB-DIIS] iter=%d WAITING: gorb=%.3e > thresh=%.3e\n",
                                ctx.iteration, gorb_check, gorb_thresh);
            }
        }

        if (gate_passes) {
            reset_to_activation(ctx);
            if (diis_formulation_ == DIISFormulation::ORBITAL && reports(4)) {
                outfile->Printf("  [GVB-DIIS] iter=%d SWITCH: %s -> ORBITAL-DIIS (n_kappa=%d)\n",
                                ctx.iteration, ctx.use_trah ? "TRAH" : "burn-in",
                                static_cast<int>(orbital_diis_cum_kappa_.size()));
            } else if (diis_formulation_ == DIISFormulation::CFM && reports(4)) {
                outfile->Printf("  [GVB-DIIS] iter=%d SWITCH: %s -> GVB-DIIS (gorb=%.3e)\n",
                                ctx.iteration, ctx.use_trah ? "TRAH" : "burn-in", gorb_check);
            }
        }
    }

    if (!ctrl_state_.is_active()) return false;

    timer_on("REKS: gvb_diis_step");
    // Orbital swaps change MO ordering, invalidating the stored error vectors and the
    // cumulative kappa.
    const bool swap_detected = gvb_diis_had_swaps_;
    gvb_diis_had_swaps_ = false;
    if (swap_detected && reports(4))
        outfile->Printf("  [GVB-DIIS] iter=%d RESET: orbital swap detected\n", ctx.iteration);

    diis_controller_.select(orbital_runtime_active_);

    if (!f_gen_current) F_gen = ctx.build_F_gen();
    int nmo = ctx.nmo;

    // Orbital gradient at the current point (= the result of the previous step).
    gvb_diis_gorb_norm_ = reks::orbital_gradient_norm(F_gen->pointer(0), ctx.active, ctx.Ncore,
                                                      ctx.nmo, {}, nullptr);
    gvb_report_witness_ = GvbReportWitness::CURRENT;

    // Frame handles only: the verdict path passes step_ctx to restart(). The build inputs are
    // filled in the solve branch, which a null step skips.
    reks::StepContext step_ctx;
    step_ctx.Ca = &ctx.Ca;
    step_ctx.Ca_ref = &Ca_ref_;
    step_ctx.cum_kappa = &orbital_diis_cum_kappa_;
    step_ctx.orb_base_map_damp = diis_orb_base_map_damp_;

    // e_sa_cur = SA energy of the current frame;
    // fon_disp = max_{b,i} |fon_cur[b][i] - fon_best[b][i]|, 0 until a best exists.
    const double e_sa_cur = reks::scf::compute_E_SA(ctx.sa_cassette, ctx.reel.E_L, ctx.reel.C_L);
    double fon_disp = 0.0;
    if (gvb_best_gorb_ >= 0.0) {
        const auto fons_cur = reks::studio::capture_fons(ctx.reel, ctx.sa_cassette, ctx.fon_blocks);
        for (size_t b = 0; b < fons_cur.size() && b < gvb_best_fons_.size(); ++b)
            for (size_t i = 0; i < fons_cur[b].size() && i < gvb_best_fons_[b].size(); ++i)
                fon_disp =
                    std::max(fon_disp, std::abs(fons_cur[b][i] - gvb_best_fons_[b][i]));
    }
    outcome_guard_.observe(e_sa_cur, fon_disp);

    // Convergence refusal at a localization signature (FON-silent lower well vs a stationary best): land on
    // best, restart the accelerator, skip the monitor review.
    bool s_adjudicated = false;
    if (gvb_diis_gorb_norm_ / ctx.nmo < gvb_d_conv_ && outcome_guard_.localization_signature()) {
        const double gorb_refused = gvb_diis_gorb_norm_;
        execute_landed_rewind(ctx, GvbReportWitness::BEST);
        // Past the restart budget the refusal keeps the rewind but drops the accelerator
        // restart and the budget count.
        if (!outcome_guard_.budget_exhausted()) {
            diis_controller_.active()->restart(step_ctx);
            outcome_guard_.on_monitor_restart();
        }
        s_adjudicated = true;
        if (reports(4))
            outfile->Printf(
                "  [DIIS-CTRL] iter=%d CONVERGENCE REFUSED: orbital symmetry breaking detected "
                "(gorb %.3e, E_sa %.10f below stationary best %.10f, fon_disp %.3e)\n",
                ctx.iteration, gorb_refused, e_sa_cur, gvb_best_E_, fon_disp);
    }

    // Basin transition (FON branch changed): drop best, re-arm detection, restart the base map here
    // without rewind.
    bool basin_fire = false;
    if (!s_adjudicated && outcome_guard_.basin_transition_fire()) {
        basin_fire = true;
        clear_best();
        outcome_guard_.on_best_invalidated();
        ctrl_state_.set_rung(QuiesceLadder::ARMED);
        if (reports(4))
            outfile->Printf(
                "  [DIIS-CTRL] iter=%d BASIN TRANSITION #%d: FON branch moved (fon_disp %.3e), "
                "re-baseline from current frame\n",
                ctx.iteration, outcome_guard_.b_fire_count(), fon_disp);
    }

    // Verdicts: VerdictMonitor rewinds on |g_orb| growth beyond rho (a rewind streak restarts), CycleMonitor
    // restarts on a stall or ORB cap streak; rewinds compose before restarts.
    if (!s_adjudicated) {
        reks::TrajectoryView traj;
        traj.iteration = ctx.iteration;
        traj.gorb = gvb_diis_gorb_norm_;
        traj.best_gorb = gvb_best_gorb_;
        traj.prev_did_extrap = gvb_prev_did_extrap_;
        traj.best_available = (gvb_best_Ca_ != nullptr);
        traj.cap_fired = orbital_runtime_active_ && orbital_cap_fired_this_iter_;
        traj.last_rollback_iter = ctrl_state_.last_rollback_iter();
        traj.energy_descent = outcome_guard_.energy_descent();

        reks::MonitorAction pre = diis_controller_.pre_step_directive(swap_detected, basin_fire);
        const reks::ResolvedAction act = diis_controller_.review(traj);

        // Past the restart budget: RETAIN (rewind keeping the subspace minus its newest vector), then on a
        // repeated landing at the same best FALLBACK (verdicts ignored); basin/swap directives stay live.
        reks::QuiesceInputs qin;
        qin.rung = ctrl_state_.rung();
        qin.quiesced = outcome_guard_.budget_exhausted();
        qin.landing_repeat = outcome_guard_.landing_repeat(gvb_best_E_, gvb_best_gorb_);
        qin.best_below_current = (gvb_best_Ca_ != nullptr) && gvb_best_gorb_ >= 0.0 &&
                                 gvb_best_gorb_ < gvb_diis_gorb_norm_;

        const reks::ResolvedAction dec = diis_controller_.apply_quiescence(act, qin);

        const double gorb_bad = gvb_diis_gorb_norm_;

        if (dec.set_rung) ctrl_state_.set_rung(dec.rung_target);
        if (dec.rewind)
            execute_landed_rewind(
                ctx, dec.landed_on_best ? GvbReportWitness::BEST : GvbReportWitness::REJECTED,
                dec.retain_history);
        if (dec.record_landing) outcome_guard_.record_landing(gvb_best_E_, gvb_best_gorb_);
        if (dec.count_toward_budget) outcome_guard_.on_monitor_restart();
        if (dec.verdict_trip) ++gvb_verdict_trips_;

        if (dec.set_rung && dec.rung_target == QuiesceLadder::FALLBACK) {
            if (reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d QUIESCE FALLBACK: retain restart landed on an "
                    "unchanged best; monitor verdicts discarded until convergence or "
                    "MAXITER (%s)\n",
                    ctx.iteration, dec.reason);
        } else if (dec.set_rung && dec.rung_target == QuiesceLadder::RETAIN) {
            if (reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d QUIESCE RETAIN: restart budget exhausted (%d); "
                    "rewind to best keeping DIIS history (%s)\n",
                    ctx.iteration, outcome_guard_.restarts_since_improvement(), dec.reason);
        } else if (dec.rewind) {
            if (dec.repeat_landing && reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d REWIND_REPEAT: landing on unchanged best "
                    "(gorb %.3e) counted toward restart budget (%d)\n",
                    ctx.iteration, gvb_best_gorb_, outcome_guard_.restarts_since_improvement());
            if (reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d REWIND_TO_BEST (%s): gorb %.3e "
                    "-> best %.3e (null step)\n",
                    ctx.iteration, dec.reason, gorb_bad, gvb_best_gorb_);
        }

        if (dec.restart) {
            diis_controller_.active()->restart(step_ctx);
            ctrl_state_.on_landing_or_restart(ctx.iteration);
            if (reports(4))
                outfile->Printf("  [DIIS-CTRL] iter=%d RESTART_BASE_MAP (%s)\n", ctx.iteration,
                                dec.reason);
        }

        // Structural directives run regardless of the monitor verdict: idempotent history
        // flush + re-anchor from the current frame.
        if (pre.verdict == reks::MonitorVerdict::RESTART_BASE_MAP) {
            diis_controller_.active()->restart(step_ctx);
            ctrl_state_.on_landing_or_restart(ctx.iteration);
            ctrl_state_.set_gate_streak(0);
            if (reports(4))
                outfile->Printf("  [DIIS-CTRL] iter=%d RESTART_BASE_MAP (%s)\n", ctx.iteration,
                                pre.reason);
        }
    }

    // Terminal MAXITER landing: exhaustion without convergence ends on the best iterate in a
    // consistent frame.
    if (ctx.iteration >= gvb_scf_maxiter_ && gvb_diis_gorb_norm_ / ctx.nmo >= gvb_d_conv_ &&
        gvb_best_Ca_ && gvb_best_gorb_ >= 0.0) {
        if (gvb_null_step_) {
            gvb_report_witness_ = GvbReportWitness::BEST;
        } else if (gvb_best_gorb_ < gvb_diis_gorb_norm_) {
            const double gorb_last = gvb_diis_gorb_norm_;
            execute_landed_rewind(ctx, GvbReportWitness::BEST);
            if (reports(4))
                outfile->Printf("  [DIIS-CTRL] iter=%d MAXITER: landing on best (gorb %.3e -> %.3e)\n",
                                ctx.iteration, gorb_last, gvb_best_gorb_);
        }
    }

    // Best-iterate capture (lowest gorb seen this SCF): the rewind target. Skipped on a null
    // step -- the reported gorb belongs to the rejected iterate, not the restored frame.
    if (!gvb_null_step_ && !outcome_guard_.refuse_capture() &&
        (gvb_best_gorb_ < 0.0 || gvb_diis_gorb_norm_ < gvb_best_gorb_)) {
        gvb_best_gorb_ = gvb_diis_gorb_norm_;
        gvb_best_E_ = e_sa_cur;
        gvb_best_iter_ = ctx.iteration;
        if (!gvb_best_Ca_) gvb_best_Ca_ = ctx.Ca->clone();
        else gvb_best_Ca_->copy(ctx.Ca);
        if (!gvb_best_eps_) gvb_best_eps_ = std::make_shared<Vector>(*ctx.eps);
        else gvb_best_eps_->copy(*ctx.eps);
        gvb_best_fons_ = reks::studio::capture_fons(ctx.reel, ctx.sa_cassette, ctx.fon_blocks);
        reks::StepContext snap_ctx;
        snap_ctx.Ca_ref = &Ca_ref_;
        snap_ctx.cum_kappa = &orbital_diis_cum_kappa_;
        gvb_best_engine_memento_ = diis_controller_.active()->snapshot(snap_ctx);
        outcome_guard_.on_capture(gvb_best_E_, gvb_best_gorb_);
    }

    if (!gvb_null_step_) {
        // Adapter step: payload + error, store, extrapolate (else base map), apply (CFM diagonalizes, ORB
        // accumulates kappa).
        std::vector<double> gamma_aa, gamma_ca;
        reks::scf::compute_gamma_vectors(ctx.Ncore, static_cast<int>(ctx.active.size()), ctx.reel,
                                         ctx.sa_cassette, ctx.x_alpha, ctx.iteration, gamma_aa,
                                         gamma_ca);

        std::vector<double> fon_dbg;
        if (reports(4) && !orbital_runtime_active_)
            fon_dbg = reks::studio::fon_vector(ctx.reel, ctx.sa_cassette, ctx.sa_cassette.sector(), 0);

        step_ctx.F_gen             = F_gen;
        step_ctx.F_MO_diag_a_L     = &ctx.reel.F_MO_diag_a_L;
        step_ctx.F_MO_diag_b_L     = &ctx.reel.F_MO_diag_b_L;
        step_ctx.C_L               = &ctx.reel.C_L;
        step_ctx.sa_cassette       = &ctx.sa_cassette;
        step_ctx.active_mo_indices = &ctx.active;
        step_ctx.Ncore             = ctx.Ncore;
        step_ctx.nmo               = nmo;
        step_ctx.gamma_aa          = &gamma_aa;
        step_ctx.gamma_ca          = &gamma_ca;
        step_ctx.SC                = ctx.SC;
        step_ctx.F_reks_MO         = ctx.reel.F_reks_MO;
        step_ctx.reel_F_reks       = ctx.reel.F_reks;
        step_ctx.Fa                = &ctx.Fa;
        step_ctx.Fb                = &ctx.Fb;
        step_ctx.iteration         = ctx.iteration;
        step_ctx.gorb              = gvb_diis_gorb_norm_;
        step_ctx.n_fon_dbg         = fon_dbg.empty() ? nullptr : &fon_dbg;

        reks::StepResult sr = diis_controller_.step(step_ctx);
        gvb_prev_did_extrap_         = sr.extrap_ok;
        orbital_cap_fired_this_iter_ = sr.cap_fired;
    } else {
        // Null step: no solve.
        gvb_prev_did_extrap_ = false;
        orbital_cap_fired_this_iter_ = false;
    }
    report_pair(gvb_report_witness_);
    timer_off("REKS: gvb_diis_step");
    return true;
}

bool GvbDiisDriver::take_null_step() {
    if (!(ctrl_state_.is_active() && gvb_null_step_)) return false;
    gvb_null_step_ = false;
    return true;
}

void GvbDiisDriver::apply_orbital_rotation(const SharedMatrix& Ca_out, const std::vector<int>& active, int Ncore,
                                           int N) {
    orb_adapter_.apply_rotation(Ca_ref_, Ca_out, orbital_diis_cum_kappa_, active, Ncore, N);
}

double GvbDiisDriver::report_energy(std::map<std::string, double>& energies, int iteration) {
    // Stand-on-best: the best iterate's components and total replace this iteration's, so the
    // reported energy is the best's and its components sum to it; the gate keeps the current gorb.
    if (iteration == gvb_best_iter_) gvb_best_energies_ = energies;
    if (gvb_report_E_override_armed_) energies = gvb_best_energies_;
    const double E_report = energies["Total Energy"];

    // dE latch: shifts prev <- cur only on the first report_energy() call of a new iteration.
    if (gvb_gate_E_iter_ != iteration) {
        gvb_gate_dE_valid_ = (gvb_gate_E_iter_ >= 0);
        gvb_gate_E_prev_ = gvb_gate_E_cur_;
        gvb_gate_E_iter_ = iteration;
    }
    gvb_gate_E_cur_ = E_report;
    if (gvb_gate_dE_valid_) gvb_gate_dE_ = E_report - gvb_gate_E_prev_;

    return E_report;
}

void GvbDiisDriver::advance_gate(int iteration, double true_rms) {
    // Streak advances once per iteration; repeat calls this iteration are no-ops.
    if (gvb_gate_streak_iter_ != iteration) {
        gvb_gate_streak_iter_ = iteration;
        const bool joint = gvb_gate_dE_valid_ && std::abs(gvb_gate_dE_) < gvb_e_conv_ &&
                           true_rms < gvb_d_conv_;
        ctrl_state_.set_gate_streak(joint ? ctrl_state_.gate_streak() + 1 : 0);
    }
}

GvbDiisDriver::Quality GvbDiisDriver::quality(int nmo) const {
    Quality q;
    q.gorb_fro = gvb_report_gorb_;
    q.gorb_rms = q.gorb_fro / nmo;
    q.landings = gvb_landing_count_;
    q.verdict_trips = gvb_verdict_trips_;
    q.rung = static_cast<int>(ctrl_state_.rung());
    q.basin_fires = outcome_guard_.b_fire_count();
    q.grazing = (gvb_gate_dE_valid_ && std::abs(gvb_gate_dE_) > 0.1 * gvb_e_conv_) ||
                (q.gorb_rms > 0.1 * gvb_d_conv_);
    return q;
}

}  // namespace reks
}  // namespace psi
