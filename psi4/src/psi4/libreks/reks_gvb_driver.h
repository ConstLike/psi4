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

#ifndef REKS_GVB_DRIVER_H
#define REKS_GVB_DRIVER_H

/// GVB-DIIS driver: activation gate, accelerator step, best-iterate rewinds, joint (dE, rms) gate.

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_diis_adapters.h"
#include "reks_diis_config.h"
#include "reks_diis_controller.h"
#include "reks_diis_outcome.h"
#include "reks_diis_state.h"
#include "reks_gvb_diis.h"
#include "reks_orbital_diis.h"
#include "reks_reel.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace reks {

/// GVB-DIIS settings.
struct GvbDriverConfig {
    DiisConfig diis;
    int start_iteration = 5;            ///< REKS_GVB_DIIS_START
    double d_convergence = 1e-5;        ///< D_CONVERGENCE
    double e_convergence = 1e-6;        ///< E_CONVERGENCE
    int maxiter = 0;                    ///< MAXITER
    int gate_streak = 2;                ///< REKS_GVB_GATE_STREAK
    bool orbital_formulation = true;    ///< REKS_DIIS_FORMULATION == ORBITAL
    double cfm_activation_gorb = 1e-3;  ///< REKS_DIIS_CFM_ACTIVATION_GORB
};

/// Inputs of one GvbDiisDriver::on_fock call; the references are valid for that call only.
struct GvbFockContext {
    int iteration;
    bool use_trah;                                           ///< TRAH warmup precedes the activation
    const std::vector<int>& active;                          ///< active MO indices (absolute, 0-based)
    int Ncore;
    int nmo;
    const studio::Cassette& sa_cassette;
    /// E_L, C_L, F_MO_diag_*_L, F_reks(_MO) read; fon_state written.
    studio::Reel& reel;
    const std::vector<studio::FonBlock>& fon_blocks;         ///< (s, gen) blocks of the best-iterate FON memento
    SharedMatrix Ca;                                         ///< MO coefficients, rewound in place
    SharedMatrix Cb;
    SharedVector eps;                                        ///< orbital energies of Ca, rewound in place
    SharedMatrix Fa;                                         ///< SCF Fock targets, written in place
    SharedMatrix Fb;
    SharedMatrix SC;                                         ///< AO<-MO map S*Ca
    double x_alpha;                                          ///< hybrid exact-exchange coefficient
    std::function<SharedMatrix()> build_F_gen;               ///< generalized Fock at the current MOs, built on demand
    std::vector<char> aa_frozen;                             ///< aa pairs out of the gate gradient (i<j packing)
};

/// GVB-DIIS (Muller, R. P. et al. J. Chem. Phys. 1994, 100, 1226), ORBITAL engine by default (Ionova, I. V.;
/// Carter, E. A. J. Chem. Phys. 1995, 102, 1251; Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472).
class GvbDiisDriver {
   public:
    GvbDiisDriver() = default;
    GvbDiisDriver(const GvbDiisDriver&) = delete;
    GvbDiisDriver& operator=(const GvbDiisDriver&) = delete;

    /// Engines, controller, outcome guard and gate thresholds from cfg.
    void configure(const GvbDriverConfig& cfg);

    /// Accelerator activated (one-way).
    bool active() const { return ctrl_state_.is_active(); }
    /// ORBITAL formulation live (vs CFM fall-back).
    bool orbital_runtime_active() const { return orbital_runtime_active_; }
    /// REKS_DIIS_FORMULATION == ORBITAL.
    bool orbital_formulation() const { return diis_formulation_ == DIISFormulation::ORBITAL; }
    /// Live ORBITAL cumulative kappa (packed non-redundant pairs).
    const std::vector<double>& cum_kappa() const { return orbital_diis_cum_kappa_; }
    /// An orbital-column swap occurred this iteration.
    void note_swaps() { gvb_diis_had_swaps_ = true; }
    /// One-shot: a stand-on-best reported pair lives exactly one iteration.
    void disarm_report_override() { gvb_report_E_override_armed_ = false; }

    /// Activation gate + the active step on this iteration's Fock; false when the driver
    /// is not active after the gate (Fa/Fb untouched).
    bool on_fock(GvbFockContext& ctx);

    /// Active driver: the activation reset at ctx's frame; inactive: no-op.
    void restart(GvbFockContext& ctx);

    /// Consume the one-shot null step flagged by a rewind.
    bool take_null_step();

    /// Ca_out = Ca_ref * exp(K(cum_kappa)) over the non-redundant rotation pairs.
    void apply_orbital_rotation(const SharedMatrix& Ca_out, const std::vector<int>& active, int Ncore, int N);

    /// Reported energy and dE latch of the joint gate; under a stand-on-best override `energies` takes the
    /// best iterate's components.
    double report_energy(std::map<std::string, double>& energies, int iteration);

    /// Advance the joint (dE, rms) persistence streak once per iteration.
    void advance_gate(int iteration, double true_rms);

    /// Orbital gradient Frobenius norm ||g_orb||_F at the current Fock.
    double gorb() const { return gvb_diis_gorb_norm_; }

    /// True while the (dE, rms) gate streak is short of REKS_GVB_GATE_STREAK; false when inactive.
    bool gate_confirm_pending() const {
        return ctrl_state_.is_active() && ctrl_state_.gate_streak() < gvb_gate_streak_req_;
    }

    /// Iteration of the last executed landing; -1 = none.
    int landed_iter() const { return gvb_landed_iter_; }

    /// A landing executed on this iteration at or past MAXITER.
    bool terminal_landing(int iteration) const {
        return gvb_landed_iter_ == iteration && iteration >= gvb_scf_maxiter_;
    }

    /// ||g_orb||_F of the witnessed iterate (the reported one, not the gated one).
    double report_gorb() const { return gvb_report_gorb_; }

    /// Terminal landing: gorb of the re-evaluated landed frame becomes the reported gorb;
    /// the stand-on-best energy override is disarmed.
    void set_terminal_report(double gorb);

    /// Convergence quality of the finished run; grazing = final iterate within a decade of a gate.
    struct Quality {
        double gorb_fro;
        double gorb_rms;
        int landings;
        int verdict_trips;
        int rung;
        int basin_fires;
        bool grazing;
    };
    Quality quality(int nmo) const;

    int verdict_trips() const { return gvb_verdict_trips_; }
    double best_gorb() const { return gvb_best_gorb_; }

   private:
    /// Which single iterate the reported (E, gorb) pair describes: CURRENT = the live one,
    /// REJECTED (bad/bad) on a null-step rewind, BEST (best/best) on stand-on-best paths.
    enum class GvbReportWitness { CURRENT, REJECTED, BEST };

    /// REKS_DIIS_FORMULATION: ORBITAL (orbital-rotation DIIS) or CFM (composite Fock matrix).
    enum class DIISFormulation { CFM, ORBITAL };

    /// Anti-cycling ladder past the restart budget: ARMED -> RETAIN (rewind keeping the subspace, damped
    /// base map) -> FALLBACK (live DIIS, verdicts ignored); reset at activation and basin transition.
    using QuiesceLadder = ControllerState::Rung;

    /// Restore the best iterate, flush DIIS (retain_history: drop only the newest, ORB only on an identical
    /// anchor), flag the null step, record the landing iteration and witness.
    void execute_landed_rewind(GvbFockContext& ctx, GvbReportWitness witness,
                               bool retain_history = false);

    /// Reported gorb: best (BEST, arms the report_energy override) or current (CURRENT/REJECTED, disarms it).
    void report_pair(GvbReportWitness witness);

    /// Drop the best iterate and its engine memento, and re-arm the trajectory monitors.
    void clear_best();

    /// Activation state at ctx's frame: both controller axes, landing accounting, best iterate, outcome guard,
    /// DIIS subspaces; ORBITAL also re-anchors Ca_ref_ = Ca with cum_kappa = 0.
    void reset_to_activation(GvbFockContext& ctx);

    ControllerState ctrl_state_;                ///< Activation x quiescence axes + landing accounting
    GVBDIISEngine gvb_diis_;                    ///< DIIS engine (F^g storage + extrapolation)
    double gvb_diis_gorb_norm_ = 0.0;          ///< Orbital gradient Frobenius norm ||g_orb||_F (not rms)
    bool gvb_diis_had_swaps_ = false;          ///< Prev-iter orbital swaps occurred
    int gvb_diis_start_ = 5;                   ///< First iteration to use F^g (burn-in with F_reks)
    double cfm_activation_gorb_ = 1e-3;        ///< CFM activation threshold on ||g_orb||_F

    /// Best-iterate memento: the lowest-gorb iterate seen this SCF. Orbitals, orbital energies,
    /// FONs and engine frame form one snapshot, captured and restored together.
    SharedMatrix gvb_best_Ca_;
    SharedVector gvb_best_eps_;
    std::vector<std::vector<double>> gvb_best_fons_;  ///< [block][i] as studio::capture_fons()
    EngineMemento gvb_best_engine_memento_;
    double gvb_best_gorb_ = -1.0;                  ///< -1 = no best captured yet
    bool gvb_prev_did_extrap_ = false;
    int gvb_verdict_trips_ = 0;                    ///< Monitor-verdict rewinds this SCF
    int gvb_landing_count_ = 0;                    ///< Executed landings this SCF (any path)
    bool gvb_null_step_ = false;                   ///< One-shot: this iteration takes no step
    int gvb_landed_iter_ = -1;                     ///< iteration of the last executed landing; -1 = none
    double gvb_best_E_ = 0.0;                      ///< E_SA at the best iterate (valid when best is set)
    int gvb_best_iter_ = -1;                       ///< iteration of the best capture; -1 = none
    std::map<std::string, double> gvb_best_energies_;  ///< energy components of the best iterate

    GvbReportWitness gvb_report_witness_ = GvbReportWitness::CURRENT;
    bool gvb_report_E_override_armed_ = false;     ///< One-shot report_energy return override (BEST pair)
    double gvb_report_gorb_ = 0.0;                 ///< ||g_orb||_F of the witnessed iterate (reported, not gated)

    /// Streak of iterations with |dE| < E_CONVERGENCE and rms < D_CONVERGENCE (counter in ctrl_state_).
    int gvb_gate_streak_iter_ = -1;                ///< iteration of the last streak update
    int gvb_gate_streak_req_ = 2;                  ///< REKS_GVB_GATE_STREAK option: confirmation depth
    double gvb_gate_E_cur_ = 0.0;                  ///< report_energy return this iteration
    double gvb_gate_E_prev_ = 0.0;                 ///< report_energy return the previous iteration
    double gvb_gate_dE_ = 0.0;                     ///< gvb_gate_E_cur_ - gvb_gate_E_prev_
    bool gvb_gate_dE_valid_ = false;               ///< A previous-iteration energy exists
    int gvb_gate_E_iter_ = -1;                     ///< iteration of the last energy record

    DiisController diis_controller_;
    OutcomeGuard outcome_guard_;
    double gvb_d_conv_ = 1e-5;                     ///< D_CONVERGENCE option
    double gvb_e_conv_ = 1e-6;                     ///< E_CONVERGENCE option
    int gvb_scf_maxiter_ = 0;                      ///< MAXITER option
    double diis_orb_base_map_damp_ = 0.6;          ///< ORB re-anchored base-map damping alpha

    DIISFormulation diis_formulation_ = DIISFormulation::ORBITAL;
    OrbitalDIISEngine orbital_diis_;
    CfmAdapter cfm_adapter_{&gvb_diis_};
    OrbAdapter orb_adapter_{&orbital_diis_};
    /// ORBITAL frame: Ca = Ca_ref_ exp(K(cum_kappa)) from form_G to the step; set at activation, re-anchored
    /// (kappa = 0) on every engine restart.
    SharedMatrix Ca_ref_;
    std::vector<double> orbital_diis_cum_kappa_;   ///< Packed non-redundant pairs (size n_nonred)
    bool orbital_runtime_active_ = false;          ///< ORBITAL formulation live (vs CFM fall-back)
    /// Previous ORBITAL step's cap-fire flag.
    bool orbital_cap_fired_this_iter_ = false;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_GVB_DRIVER_H
