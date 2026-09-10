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

#ifndef REKS_H
#define REKS_H

#include "hf.h"
#include "reks_si_types.h"
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_reel.h"
#include "reks_rdm.h"
#include "reks_convergence.h"
#include "reks_trah_solver.h"
#include "reks_gradient_engine.h"
#include "reks_fon_solver.h"
#include "reks_gvb_diis.h"
#include "reks_orbital_diis.h"
#include "reks_diis_adapters.h"
#include "reks_diis_controller.h"
#include "reks_diis_outcome.h"
#include "reks_diis_state.h"
#include "reks_scf.h"
#include "psi4/libfock/v.h"
#include <vector>
#include <memory>
#include <utility>
#include <array>
#include <chrono>
#include <map>
#include <string>

namespace psi {
namespace reks { struct FONObjective; }
namespace scf {


/// @class REKS
/// @brief Restricted Ensemble-referenced Kohn-Sham (REKS) SCF engine.
///
/// SCF over the Catalog/Cassette/Reel active space; the SCF flow is
/// the same across all REKS(N,M) variants.
class REKS : public HF {
   protected:
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
    SharedMatrix wK_;

    /// Offset-ABI view of the loaded catalog blob (base pointer + header).
    reks::studio::Catalog catalog_{};

    /// Mutable runtime state owned by REKS for the Catalog path (FONs,
    /// per-microstate Focks/energies, lagrangians, base densities, ...).
    /// Sized from catalog_ + nsopi_/nso_.
    reks::studio::Reel reel_;

    /// SA cassette: the SCF ensemble pool (union over the run's declared sectors)
    /// plus the sector-less extra determinants at the tail; sole dispatcher over
    /// the global-L space, with per-sector config sublists for FON attribution.
    reks::studio::Cassette sa_cassette_;

    /// SI cassettes: one per user SI_REKS_CONFIGS entry, in order.
    std::vector<reks::studio::Cassette> si_cassettes_;

    /// Built once at init from the SA cassette and si_cassettes_ (see CallSheet).
    reks::studio::CallSheet call_sheet_;

    /// Joint-optimizer FON block table: one entry per (run sector s, generation gen)
    /// with gen in [0, sa_cassette_.n_generations()). The offset field is unset here
    /// (0); each optimizer pass lays out offsets after the orbital-rotation block.
    /// Rebuilt once per init from the SA cassette.
    std::vector<reks::studio::FonBlock> fon_blocks_;

    /// Blob sector -> run sector ordinal, -1 when a blob sector is not part of this
    /// run. Built once at init.
    std::vector<int> blob_to_run_sector_;

    /// Run sector ordinal s -> blob sector index (the inverse of blob_to_run_sector_
    /// over declared sectors). Size = number of declared run sectors.
    std::vector<int> run_to_blob_;

    /// Build fon_blocks_ from the SA cassette and (re)size persistent per-block
    /// optimizer state (AHC FON smoothing) seeded per (s, gen). Idempotent.
    void build_fon_blocks();

    /// Run sector ordinal owning global config index K (via the catalog's blob
    /// sector decode and blob_to_run_sector_). -1 if K is outside the run's sectors.
    [[nodiscard]] int run_sector_of_config(int K) const {
        const int b = catalog_.sector_of_config(K);
        return (b >= 0 && b < static_cast<int>(blob_to_run_sector_.size()))
            ? blob_to_run_sector_[b] : -1;
    }

    /// Catalog view bound to run sector s's blob manifold: sector-local config
    /// indices and names resolve within that sector's block, and n_generations()
    /// reflects it. Cheap value copy (the blob buffer is shared by pointer).
    [[nodiscard]] reks::studio::Catalog catalog_sector_view(int s) const {
        reks::studio::Catalog c = catalog_;
        c.active_sector = run_to_blob_[s];
        return c;
    }

    /// Assemble the symmetrized F_REKS_MO = sum_L C_L-weighted microstate Focks
    /// into reel_.F_reks_MO (symmetric).
    void assemble_F_reks_MO();

    /// Build the core density into reel_.D_core and the one-electron trace of
    /// every SA-referenced occupation pattern into reel_.base_density_e1.
    void build_base_densities();

    /// Da_ = reel_.D_core + sum_{i in [Ncore_, nalphapi_[0])} c_i c_i^T, the
    /// aufbau-proxy closed-shell density over the occupied columns of Ca_.
    /// Requires reel_.D_core current for the same Ca_.
    void form_Da_from_core();

    /// Build the SA per-L Focks (AO + MO) for the SA cassette's active microstate set.
    /// SI-only L's outside that set are left unfilled until post-SCF.
    void build_sa_focks();

    /// Alpha/beta base-density bitmask for microstate L: sum_i occ[i] << i over
    /// the active orbitals, read from sa_cassette_.microstate(L).
    [[nodiscard]] int base_idx(int L, bool alpha) const;

    int Ncore_ = -1;                       ///< Number of core (doubly occupied) orbitals
    std::vector<int> active_mo_indices_;

    bool si_computed_ = false;

    std::vector<reks::SIResult> si_results_;

    /// Primary (first) cassette result; empty while compute_si() has not run.
    [[nodiscard]] const reks::SIResult& primary_si_result() const;

    /// (ij|kl) tiles over the canonical union of the SI coupling pairs and the
    /// full off-diagonal active pair set. Rebuilt on every compute_si() run.
    reks::ActiveEriTiles active_eri_pool_;

    /// Persistent workspace for build_generalized_fock(). Aliased to its return
    /// value. Callers must not hold the SharedMatrix across another call.
    mutable SharedMatrix F_gen_workspace_;

    /// Per-iter Ca_subset("SO","OCC") cache, valid within one SCF iter (Ca_ does
    /// not change between fill and use).
    SharedMatrix C_occ_cache_;

    /// UKS XC for microstates whose D_alpha may differ from D_beta.
    std::shared_ptr<psi::UV> uv_potential_;

    /// RKS XC for the closed-shell SCF density Da_ (Db_ aliases Da_).
    std::shared_ptr<psi::RV> rv_potential_;

    double vv10_E_rv_ = 0.0;  ///< VV10 nonlocal correlation energy (RV functional), added to E_total

    /// Polarized (spin-unrestricted) version of the functional.
    std::shared_ptr<SuperFunctional> polarized_functional_;

    /// Stage timing (parallel to timer_on/off; gated by report level 4).
    mutable std::map<std::string, double> scf_iter_times_;
    mutable std::map<std::string, double> scf_total_times_;  ///< sum over all SCF iterations
    mutable std::map<std::string, double> post_scf_times_;   ///< post-SCF accumulator (compute_si etc.)
    mutable std::chrono::steady_clock::time_point iter_start_time_;
    mutable std::chrono::steady_clock::time_point scf_start_time_;
    mutable std::chrono::steady_clock::time_point post_scf_start_time_;  ///< post-SCF span anchor
    mutable int last_logged_iter_ = -1;
    mutable int n_iters_logged_ = 0;
    mutable bool scf_summary_printed_ = false;

    void log_iter_timings_(int iter, double iter_total_seconds) const;
    void log_scf_summary_() const;
    void log_post_scf_summary_(const char* phase_label) const;
    void log_jk_xc_summary_(const char* phase_label,
                            const std::map<std::string, double>& bucket,
                            double phase_total_seconds) const;

    SharedMatrix                      S_half_;           ///< S^{1/2}, built in reks_common_init
    std::vector<int>                  ao2atom_;          ///< AO -> atom center
    std::vector<std::vector<int>>     atom_to_ao_;       ///< inverse map
    int                               n_atoms_ipr_ = 0;  ///< molecule_->natom()

    double      lambda_ipr_ = 0.0;       ///< REKS_DELOC_IPR_PENALTY (constant over SCF)
    std::string ipr_method_ = "LOWDIN";  ///< LOWDIN | MULLIKEN

    /// Last-computed IPR diagnostics: total penalty, energy term, max gradient.
    mutable double last_ipr_total_ = 0.0;
    mutable double last_ipr_epen_  = 0.0;
    mutable double last_ipr_gmax_  = 0.0;

    /// Lowdin-population memo keyed on a clone of the active columns of Ca_:
    /// p_A/a_A reused while those columns compare equal, recomputed when Ca_ rotates.
    mutable SharedMatrix                          ipr_pop_Ca_act_;
    mutable std::vector<std::vector<double>>      ipr_pop_p_A_;
    mutable std::vector<std::vector<double>>      ipr_pop_a_A_;

    /// Rolling IPR history (2-slot): previous and previous-previous iteration.
    double ipr_prev_     = -1.0;
    double ipr_prevprev_ = -1.0;

    /// Streak fail-safe: guard disables itself when target basin is unreachable.
    int    ipr_guard_streak_           = 0;
    double ipr_guard_streak_start_ipr_ = -1.0;
    bool   ipr_guard_disabled_         = false;
    int    ipr_guard_streak_limit_     = 10;  ///< REKS_DELOC_IPR_GUARD_STREAK_LIMIT

    /// Build S^{1/2}, ao2atom, atom_to_ao. Throws on non-c1 or singular S.
    void build_ipr_static_data();

    /// Loewdin populations + pair overlaps for current Ca.
    ///   p_A[i][A] = sum_{mu in A} Ctil_{mu i}^2
    ///   a_A[pair_idx(i,j)][A] = sum_{mu in A} Ctil_{mu i} Ctil_{mu j}  (upper-triangle)
    /// LOWDIN uses Ctil = S^{1/2} * Ca; MULLIKEN uses the symmetric mixed-product form.
    void compute_loewdin_populations(
        std::vector<std::vector<double>>& p_A,
        std::vector<std::vector<double>>& a_A) const;

    /// Smallest active-MO participation ratio, PR_i = 1 / sum_A p_A(i)^2, over the
    /// populations of compute_loewdin_populations. -1.0 when no active space is set.
    double active_pr_min() const;

    /// Return E_pen = lambda * IPR_total. Populates last_ipr_{total,epen}.
    double compute_ipr_penalty_energy() const;

    /// Add IPR contribution to the MO-basis generalized Fock active-active
    /// block:  F_gen[p,q] += 2*lambda * sum_A p_A(p) * a_A(p,q).
    void add_ipr_contribution_to_F_gen(SharedMatrix F_gen) const;

    /// Add full IPR Hessian contribution (on- and off-diagonal) into hess
    /// active-active block. hess is row-major [N x N] with N = n_rot +
    /// sum_gen n_fon_gen[gen]; only [0..n_rot) x [0..n_rot) is touched.
    void add_ipr_contribution_to_hessian(
        std::vector<double>& hess, int N) const;

    /// Diagonal of the active-active IPR Hessian, one entry per rotation pair
    /// (p<q), in the same units as add_ipr_contribution_to_hessian.
    void ipr_hessian_diag_aa(std::vector<double>& diag_out) const;

    /// Per-iteration IPR diagnostic: summary line, plus n/m-FON line and, at report
    /// level 5, the per-pair gradient dump.
    void print_ipr_iteration_diagnostic(int iter, double g_orb_max) const;

    /// Advances the IPR history/streak; returns true on a basin-divide crossing.
    /// Never fires unless lambda_ipr > 0.
    bool ipr_basin_guard_prestep();

    /// Final summary after SCF convergence.
    void print_ipr_final_diagnostic() const;

    /// Per-generation FON lower bound, indexed by generation (0=n, 1=m, 2=u, ...).
    /// reks_fon_lower_[g] >= 0 pins active generation-g geminal FONs to
    /// [reks_fon_lower_[g], 2]; < 0 disables. Set from REKS_<G>_FON options,
    /// each defaulting to -1.
    std::array<double, reks::studio::kMaxGen> reks_fon_lower_{};


    std::vector<double> prev_fon_;       ///< prev-iter gen-0 FON, all sectors concatenated, for the printed df delta
    double prev_E_total_ = 0.0;          ///< prev-iter E_SA for dE / oscillation
    double prev_dE_ = 0.0;               ///< prev dE for sign-flip detection
    double last_computed_dE_ = 0.0;      ///< dE before prev_E_total_ is updated
    int energy_oscillation_count_ = 0;   ///< cumulative dE sign-flip count

    /// Adaptive level-shift + MOM-style reorder applied around the F_MO diagonalization.
    reks::OrbitalGuard orbital_guard_;

    /// TRAH combined solver (joint orbital+FON) toggle.
    bool use_trah_ = false;
    /// True only while GVB-DIIS has not activated (ctrl_state_ inactive); use_trah_ itself
    /// never resets after a TRAH warmup hands off.
    bool trah_live() const { return use_trah_ && !ctrl_state_.is_active(); }
    reks::TRAHState trah_state_;               ///< Trust region persistent state
    bool trah_kappa_disabled_ = false;         ///< Kappa rotations disabled after persistent cycle

    /// AHC (Adaptive Hessian Correction): PI-controller from rho feedback adds
    /// a diagonal shift to H_oo and to each active generation's H_ff[gen]; one
    /// EMA fraction per block, kept separate because block scales differ by
    /// ~O(10) at boundary.
    double ahc_rho_smooth_ = 1.0;                 ///< EMA of step quality rho
    double ahc_pred_orb_frac_smooth_ = 0.5;       ///< EMA orbital fraction of pred dE
    /// EMA FON fraction of pred dE, per FON block (parallel to fon_blocks_). gen-0
    /// blocks seed at 0.5, higher generations at 0.0 (inactive until they carry
    /// prediction weight). Sized and seeded in build_fon_blocks.
    std::vector<double> ahc_pred_fon_frac_smooth_;
    double ahc_integral_ = 0.0;                   ///< Accumulated I-term for AHC
    static constexpr double AHC_BETA_ = 0.3;      ///< EMA smoothing coefficient
    static constexpr double AHC_I_RATE_ = 0.15;   ///< I-term accumulation rate
    static constexpr double AHC_I_MAX_ = 20.0;    ///< I-term anti-windup cap

    /// KKT-boundary convergence state (PGC) and the previous orbital step.
    bool pgc_boundary_stationary_ = false;      ///< All FON at KKT boundary
    double pgc_gorb_norm_ = 1.0;               ///< ||g_orb|| at boundary
    double prev_gorb_at_boundary_ = -1.0;      ///< ||g_orb|| at previous boundary iter (-1 = not set)
    std::vector<double> prev_orb_step_;        ///< previous raw TRS orbital step (Aitken soft-mode extrapolation)

    /// READ guess with all FON pinned at boundary -> zero TRAH step (KKT-OK but
    /// wrong basin). On detection, CI-init restarts FON to escape the trap.
    bool read_boundary_trap_ = false;

    /// BFGS H_oo: rank-2 quasi-Newton update with Powell damping (PD-preserving),
    /// initialized from REKSGradientEngine diagonal.
    std::vector<double> bfgs_B_;           ///< H_oo Hessian approximation (n_rot*n_rot, row-major)
    std::vector<double> bfgs_prev_g_;      ///< Previous combined gradient (N)
    std::vector<double> bfgs_prev_step_;   ///< Previous applied step (N)
    bool bfgs_initialized_ = false;
    int bfgs_update_count_ = 0;            ///< Number of successful BFGS updates
    int consecutive_bad_rho_ = 0;          ///< Consecutive rho < 0 events (triggers BFGS reset at 2)

    /// kDIIS: orbital-gradient DIIS at boundary stagnation (FON locked, kappa
    /// converging slowly), SO/DIIS kappa-reset variant. Fischer, T. H.; Almlof, J.
    /// J. Phys. Chem. 1992, 96, 9768.
    std::vector<std::vector<double>> kdiis_g_history_;  ///< Ring buffer of orbital gradients (each n_rot)
    int kdiis_count_ = 0;              ///< Stored gradients in ring buffer (up to KDIIS_MAX_VEC)
    int kdiis_oldest_ = 0;             ///< Index of oldest entry in ring buffer
    bool kdiis_active_ = false;
    int kdiis_stag_count_ = 0;         ///< Consecutive boundary stagnation iterations
    int kdiis_energy_fail_ = 0;        ///< Consecutive energy increases while kDIIS active
    bool kdiis_used_prev_ = false;
    static constexpr int KDIIS_MAX_VEC = 6;    ///< Max gradient vectors in ring buffer
    static constexpr int KDIIS_MIN_VEC = 3;    ///< Min vectors before extrapolation
    static constexpr int KDIIS_TRIGGER = 3;    ///< Stagnation iterations before kDIIS activates

    /// TRAH combined FON+orbital step.
    void combined_step();

    /// Fixed base-sector ordinal of the SA cassette; the default sector for calls
    /// that take no explicit run-sector argument.
    int sa_sector() const { return sa_cassette_.sector(); }
    /// SA-context FON snapshot (reel_.fon_state[sa_sector()]).
    const reks::studio::FONSnapshot& sa_fons() const { return reel_.fon_state[sa_sector()]; }
    reks::studio::FONSnapshot&       sa_fons()       { return reel_.fon_state[sa_sector()]; }

    /// Per-block FON working vector (joint-position-indexed within (sector s,
    /// generation gen)): maps that block's active geminal pool to/from
    /// reel_.fon_state[s]. Position i <-> absolute geminal
    /// sa_cassette_.geminals_active(s, gen)[i].
    std::vector<double> fon_vector(int s, int gen) const {
        const auto& act = sa_cassette_.geminals_active(s, gen);
        std::vector<double> v;
        v.reserve(act.size());
        for (int g : act) v.push_back(reel_.fon_state[s].layers[gen][g].p);
        return v;
    }
    void set_fon_vector(int s, int gen, const std::vector<double>& fons) {
        const auto& act = sa_cassette_.geminals_active(s, gen);
        if (fons.size() != act.size())
            throw PSIEXCEPTION("set_fon_vector: size does not match active geminal pool");
        for (size_t i = 0; i < act.size(); ++i)
            reel_.fon_state[s].layers[gen][act[i]] = {fons[i], 2.0 - fons[i]};
    }
    /// Named wrappers for the SA sector's two lowest generations: n = gen 0, m = gen 1.
    std::vector<double> n_fon_vector() const { return fon_vector(sa_sector(), 0); }
    std::vector<double> m_fon_vector() const { return fon_vector(sa_sector(), 1); }

    /// FON snapshot indexed by FON block (parallel to fon_blocks_); entry b holds
    /// fon_vector(fon_blocks_[b].s, fon_blocks_[b].gen).
    std::vector<std::vector<double>> capture_all_fons() const {
        std::vector<std::vector<double>> out(fon_blocks_.size());
        for (size_t b = 0; b < fon_blocks_.size(); ++b)
            out[b] = fon_vector(fon_blocks_[b].s, fon_blocks_[b].gen);
        return out;
    }
    void restore_all_fons(const std::vector<std::vector<double>>& fons) {
        for (size_t b = 0; b < fons.size() && b < fon_blocks_.size(); ++b)
            set_fon_vector(fon_blocks_[b].s, fon_blocks_[b].gen, fons[b]);
    }

    /// GVB-DIIS (Muller, R. P. et al. J. Chem. Phys. 1994, 100, 1226).
    bool use_gvb_diis_ = false;                ///< Use GVB-DIIS convergence accelerator
    reks::ControllerState ctrl_state_;          ///< Activation x quiescence axes + landing accounting
    bool use_plain_ = false;                   ///< No accelerator: plain F_reks diag + FON Newton
    reks::GVBDIISEngine gvb_diis_;             ///< DIIS engine (F^g storage + extrapolation)
    SharedMatrix F_gen_;                       ///< Cached generalized Fock for reuse within the iter
    double gvb_diis_gorb_norm_ = 0.0;          ///< Orbital gradient Frobenius norm ||g_orb||_F (not rms)
    bool gvb_diis_had_swaps_ = false;          ///< Prev-iter orbital swaps occurred
    int gvb_diis_start_ = 5;                   ///< First iteration to use F^g (burn-in with F_reks)
    bool gvb_fon_initialized_ = false;         ///< First gvb_fon_step() call sets FON=1.0
    int gvb_fon_step_swaps_ = 0;               ///< Column swaps made by the last gvb_fon_step()

    /// Settings shared by both gvb_fon_step() micro-solves; block-local overrides are the
    /// per-generation lower bound and the history switch.
    reks::MicroSolverConfig fon_micro_cfg_;
    /// Distance from a FON bound inside which the gen-0 multi-start scan calls a minimum
    /// boundary-trapped and prefers the interior one.
    double fon_micro_boundary_skin_ = 0.05;

    /// Best-iterate memento: the lowest-gorb iterate seen this SCF. Orbitals, FONs and engine
    /// frame form one consistent snapshot, captured and restored together.
    SharedMatrix gvb_best_Ca_;
    std::vector<std::vector<double>> gvb_best_fons_;  ///< [block][i] as capture_all_fons()
    reks::EngineMemento gvb_best_engine_memento_;
    double gvb_best_gorb_ = -1.0;                  ///< -1 = no best captured yet
    bool gvb_prev_did_extrap_ = false;
    int gvb_verdict_trips_ = 0;                    ///< Monitor-verdict rewinds this SCF
    int gvb_landing_count_ = 0;                    ///< Executed landings this SCF (any path)
    bool gvb_null_step_ = false;                   ///< One-shot: this iteration takes no step
    double gvb_best_E_ = 0.0;                      ///< E_SA at the best iterate (valid when best is set)

    /// Which single iterate the reported (E, gorb) pair describes: CURRENT = the live one,
    /// REJECTED (bad/bad) on a null-step rewind, BEST (best/best) on stand-on-best paths.
    enum class GvbReportWitness { CURRENT, REJECTED, BEST };
    GvbReportWitness gvb_report_witness_ = GvbReportWitness::CURRENT;
    bool gvb_report_E_override_armed_ = false;     ///< One-shot compute_E return override (BEST pair)
    double gvb_report_E_override_ = 0.0;
    double gvb_report_gorb_ = 0.0;                 ///< ||g_orb||_F of the witnessed iterate (reported, not gated)

    /// Frobenius orbital gradient -> rms per matrix element (D_CONVERGENCE units).
    double gorb_rms(double gorb_fro) const { return gorb_fro / nmopi_[0]; }

    /// Persistence latch: streak of consecutive iterations where both
    /// |dE| < E_CONVERGENCE and rms < D_CONVERGENCE hold. The streak counter itself
    /// lives in ctrl_state_, which every landing/restart path clears.
    int gvb_gate_streak_iter_ = -1;                ///< iteration_ of the last streak update
    int gvb_gate_streak_req_ = 2;                  ///< REKS_GVB_GATE_STREAK option: confirmation depth
    double gvb_gate_E_cur_ = 0.0;                  ///< compute_E return this iteration
    double gvb_gate_E_prev_ = 0.0;                 ///< compute_E return the previous iteration
    double gvb_gate_dE_ = 0.0;                     ///< gvb_gate_E_cur_ - gvb_gate_E_prev_
    bool gvb_gate_dE_valid_ = false;               ///< A previous-iteration energy exists
    int gvb_gate_E_iter_ = -1;                     ///< iteration_ of the last energy record

    /// Restore the best iterate (orbitals, FONs, engine frame), flush the engine's DIIS history,
    /// flag the null step, record the landing iteration and the witness.
    /// retain_history keeps the subspace minus its newest vector instead of flushing; for the
    /// frame-relative ORB history it holds only when the best memento anchor bitwise-matches
    /// the current anchor.
    void execute_landed_rewind(GvbReportWitness witness, bool retain_history = false);

    /// Per-run convergence quality: wfn variables + the [QUALITY] trace line.
    void report_convergence_quality();
    /// Apply the witness to the reported pair: gvb_report_gorb_ takes gvb_best_gorb_ under BEST,
    /// which also arms the one-shot compute_E override with gvb_best_E_, and gvb_diis_gorb_norm_
    /// under CURRENT/REJECTED, which disarm it. The gated member is never written here.
    void gvb_report_pair(GvbReportWitness witness);

    reks::DiisController diis_controller_;
    reks::OutcomeGuard outcome_guard_;
    double gvb_d_conv_ = 1e-5;                     ///< D_CONVERGENCE option
    double gvb_e_conv_ = 1e-6;                     ///< E_CONVERGENCE option
    int gvb_scf_maxiter_ = 0;                      ///< MAXITER option
    /// Anti-cycling ladder past the restart budget: ARMED -> RETAIN (rewind to best keeping
    /// the DIIS subspace, damped base map) -> FALLBACK (live DIIS, monitor verdicts
    /// discarded). Reset at activation and on basin transition. Held in ctrl_state_.
    using QuiesceLadder = reks::ControllerState::Rung;
    double diis_orb_base_map_damp_ = 0.6;          ///< ORB re-anchored base-map damping alpha

    /// Orbital-rotation DIIS (Ionova, I. V.; Carter, E. A. J. Chem. Phys. 1995, 102, 1251;
    /// Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472); coexists with CFM (composite Fock
    /// matrix), selected by REKS_DIIS_FORMULATION = ORBITAL (default) or CFM.
    enum class DIISFormulation { CFM, ORBITAL };
    DIISFormulation diis_formulation_ = DIISFormulation::ORBITAL;
    reks::OrbitalDIISEngine orbital_diis_;
    reks::CfmAdapter cfm_adapter_{&gvb_diis_};
    reks::OrbAdapter orb_adapter_{&orbital_diis_};
    /// Live ORBITAL frame: Ca_ = Ca_ref_ * exp(K(orbital_diis_cum_kappa_)). Ca_ref_ is set
    /// at ORBITAL activation and re-anchored to Ca_ (kappa -> 0) on every engine restart.
    SharedMatrix Ca_ref_;
    std::vector<double> orbital_diis_cum_kappa_;   ///< Packed non-redundant pairs (size n_nonred)
    bool orbital_runtime_active_ = false;          ///< ORBITAL formulation live (vs CFM fall-back)
    /// Previous ORBITAL step's cap-fire flag.
    bool orbital_cap_fired_this_iter_ = false;

    /// Norm of the packed non-redundant generator carrying the orbitals from the previous
    /// GVB-DIIS iterate to the current one: K1 = (U - U^T)/2, U = Ca_prev^T S Ca_, masked to
    /// the ca/cv/av/aa blocks. Same definition on both formulations: CFM's diagonalization
    /// yields U directly; ORBITAL reads the first-order generator off the orbitals.
    double gvb_step_norm(SharedMatrix Ca_prev, const SharedMatrix& SC) const;

    SharedMatrix gvb_prev_Ca_;                     ///< previous iterate's orbitals
    double gvb_step_norm_ = -1.0;                  ///< this iterate's step norm (-1 = no sample)
    double gvb_prev_step_norm_ = -1.0;             ///< previous sample
    double gvb_prev_gorb_ = -1.0;                  ///< ||g_orb||_F of the previous sample

    /// Orbital-gradient norm over all rotation pairs (aa, ca, cv, av). A non-null blocks
    /// receives the per-block breakdown of the same pass.
    double compute_full_gradient_norm(SharedMatrix F_gen,
                                      reks::GradientBlocks* blocks = nullptr) const;

    /// Curvature inputs of the packed Newton denominator, in the OrbitalDIISEngine layout.
    /// Empty on a gradient-only probe, where the denominator is discarded.
    struct NKCurvature {
        std::vector<double> gamma_aa;
        std::vector<double> gamma_ca;
        std::vector<double> ipr_diag_aa;
    };

    /// Packed non-redundant rotation count in the OrbitalDIISEngine layout.
    int nk_n_rot() const;

    /// Bound state of every FON degree of freedom, ordered as capture_all_fons() flattened:
    /// -1 at the lower bound, +1 at the upper bound, 0 free.
    std::vector<int> nk_fon_active_set() const;

    /// Reduced orbital gradient at Ca = nk_Ca_anchor_ * exp(K(kappa)), packed in the
    /// compute_newton_rotation() layout; equivalent to a standard SCF Fock build, then
    /// restores Ca_, the FON state and the swap counters.
    /// A non-null denom_out also harvests the shaped Newton denominator, which requires the
    /// full-MO XC diagonal; a null one runs the cheaper gradient-only Fock build.
    /// With commit set the probe keeps its point: Ca_, the relaxed FON and reel_ are left at
    /// kappa.
    /// Returns false when a FON column swap or a FON active-set move makes the probe point
    /// incomparable with the base point.
    bool reduced_gradient_at(const std::vector<double>& kappa, const NKCurvature& curv,
                             std::vector<double>& g_out, std::vector<double>* denom_out,
                             bool commit = false);

    /// H.v by one-sided finite difference of reduced_gradient_at(); v is normalized here,
    /// g0 is the base-point gradient and is not recomputed. Ghat is g0 unpacked to a skew
    /// matrix, built once per base point because g0 is fixed for a whole Krylov solve.
    bool nk_Hv(const std::vector<double>& v, const std::vector<double>& g0, const SharedMatrix& Ghat,
               double h, const NKCurvature& curv, std::vector<double>& Hv_out);

    SharedMatrix nk_Ca_anchor_;                ///< Rotation anchor of the probe frame
    bool nk_probe_gradient_only_ = false;      ///< Probe Fock build narrows the MO diagonal

    /// Extent of the per-microstate MO Fock diagonal: true = all nmo, false = the active
    /// window [Ncore, Ncore + n_active) alone.
    bool need_full_mo_diag() const { return use_gvb_diis_ && !nk_probe_gradient_only_; }

    /// KKT multiplier of every FON pinned at a bound, zero for a free one, ordered as
    /// capture_all_fons() flattened. The multiplier is the unmasked FON gradient.
    std::vector<double> nk_fon_kkt_multipliers() const;

    /// miniTRAH episode at the current iterate: measures the finite-difference step, runs the
    /// outer Newton loop on g(kappa) = 0 in its own rotation frame anchored here, and on a
    /// non-zero accepted step lands Ca_, epsilon_a_ and reel_ at that point. Returns true
    /// when it landed.
    bool nk_run_episode(const NKCurvature& curv, const std::vector<double>& denom_live);

    bool nk_enabled_ = false;                  ///< REKS_GVB_NK
    int nk_trigger_ = 4;                       ///< REKS_GVB_NK_TRIGGER
    int nk_min_landings_ = 3;                  ///< REKS_GVB_NK_MIN_LANDINGS
    int nk_max_reject_ = 1;                    ///< REKS_GVB_NK_MAX_REJECT
    int nk_fon_micro_ = 1;                     ///< REKS_GVB_NK_FON_MICRO
    bool nk_require_active_fon_ = false;       ///< REKS_GVB_NK_REQUIRE_ACTIVE_FON
    int nk_max_grad_ = 150;                    ///< REKS_GVB_NK_MAX_GRAD
    int nk_max_macro_ = 10;                    ///< REKS_GVB_NK_MAX_MACRO
    int nk_max_inner_ = 20;                    ///< REKS_GVB_NK_MAX_INNER
    double nk_fd_h_ = 0.0;                     ///< REKS_GVB_NK_FD_H, 0 = measure at entry
    double nk_exit_factor_ = 0.5;              ///< REKS_GVB_NK_EXIT_FACTOR
    int nk_trigger_streak_ = 0;                ///< consecutive stall-signature iterations
    int nk_trigger_streak_max_ = 0;            ///< longest streak of this SCF
    bool nk_episode_done_ = false;             ///< an episode already ran this SCF
    int nk_grad_evals_ = 0;                    ///< gradient evaluations spent by the episode

    /// FON Newton solve with trust region + boundary guard (GVB-DIIS mode).
    void gvb_fon_step();

    /// Track E_SA across iterations: update last_computed_dE_, flag sign-flip oscillations.
    void print_sa_diagnostics();

    /// One-line per-iter diagnostic: "[REKS] iter=.. E_SA=.. dE=.." + FON + Fock gaps (grep-friendly).
    void print_iteration_diagnostics();

    /// READ-guess FON channel bypassing the reks_common_init reset of reel_.fon_state
    /// to FON=1.0: one gen-0 (n) matrix per run sector s. Presence of any sector's
    /// gen-0 guess + guess_Ca_ marks the loaded basis as trusted (resets skipped).
    std::vector<SharedMatrix> guess_fon_;

    /// READ guesses for the higher FON generations (gen >= 1: m, u, ...), per run
    /// sector s: guess_fon_upper_[s][gen] holds that generation's saved FON, null
    /// when absent.
    std::vector<std::array<SharedMatrix, reks::studio::kMaxGen>> guess_fon_upper_;

    /// True when any run sector carries a gen-0 READ FON guess (trusted-basis sentinel).
    [[nodiscard]] bool has_read_fon_guess() const {
        for (const auto& m : guess_fon_) if (m) return true;
        return false;
    }
    /// Clear every sector's READ FON guess (boundary-trap escape).
    void clear_read_fon_guess() {
        for (auto& m : guess_fon_) m.reset();
    }

    /// Allocates SCF matrices with beta aliased to alpha (REKS is restricted);
    /// sets the report level from REKS_REPORT_LEVEL, then dispatches to subclass_init().
    void reks_common_init();

    /// iter<=0 form_C on a trusted READ guess: no-op on Ca_, copies F_MO diagonal
    /// into epsilon_a_.
    void form_C_read_guess_update();

    /// Validate reel_ array sizes against the catalog (base-density count,
    /// SA microstate count); throws on mismatch. Allocates nothing.
    void allocate_reks_matrices();

    /// Print per-pool memory footprint of REKS-owned matrices.
    /// `tag` is appended to the header.
    void print_memory_footprint(const std::string& tag) const;

    /// Compute microstate energies E_L
    void compute_sa_energies();

    /// Compute the per-microstate weighting factors C_L.
    void compute_weighting_factors();

    /// Build each SI cassette's Hamiltonian from the converged SCF orbitals into
    /// si_results_. No-op if already run (si_computed_).
    void compute_si();

    /// Result of state-property computation: permanent/transition dipoles + osc strengths.
    struct SIProperties {
        /// Every member is sized over the reported states n_rep = SIResult::n_report():
        /// n_pairs = n_rep*(n_rep-1)/2, and si_pair_index takes n_rep as its bound.
        std::vector<std::array<double, 3>> perm_dipole;   ///< [n_rep], atomic units
        std::vector<std::array<double, 3>> trans_dipole;  ///< [n_pairs], atomic units
        std::vector<double> osc_str;                      ///< [n_pairs]
        /// Adiabatic 1-RDM tensor, block-row stride n_rep: block [I*n_rep+J]*N_act^2.
        /// Diagonal (I=J) = permanent state density; off-diagonal (I!=J) = transition
        /// 1-RDM. Empty when no diabatic RDM.
        std::vector<double> rho_adiab;
    };

    /// Cassette-independent dipole data: nuclear + core dipole and the active-block
    /// MO dipole over the frozen orbitals (Ca_fock_transform).
    struct DipoleIntegralCache {
        std::array<double, 3> mu_nuc  = {0.0, 0.0, 0.0};
        std::array<double, 3> mu_core = {0.0, 0.0, 0.0};
        std::array<std::vector<double>, 3> d_act;   ///< [axis][p*N_act+q]
    };

    /// Build the cassette-independent dipole integrals.
    DipoleIntegralCache build_dipole_integral_cache() const;

    /// Pure: permanent + transition dipoles + oscillator strengths for the given
    /// SI cassette and its diagonalized result.
    SIProperties compute_state_properties(const reks::studio::Cassette& cassette,
                                          const reks::SIResult& sir,
                                          const reks::rdm::DiabaticRdm& rho_diab,
                                          const DipoleIntegralCache& dip) const;

    /// Always-on MO-consistency block (Ca vs Ca_fock_transform; window around active).
    void print_mo_consistency() const;

    /// Debug-only F_REKS_MO block (off-diag occ/act -> 0 at convergence).
    void print_f_reks_mo_debug() const;

    /// `==> METHOD State Properties (K=N) <==` section banner. Label reflects the
    /// SI cassette's own pool dimension (not the variant catalog default).
    void print_state_properties_header(const reks::studio::Cassette& cassette, int K) const;

    /// Print the SI Hamiltonian + Overlap matrices.
    void print_si_hamiltonian_overlap(const reks::studio::Cassette& cassette,
                                      const reks::SIResult& sir) const;

    /// Per-active-orbital occupations plus orbital energies and the geminal pair
    /// partitioning of every scheme the cassette spans. One FON row per (active
    /// generation, scheme) the cassette populates, labelled by generation_letter.
    void print_fon_table(const reks::studio::Cassette& cassette) const;

    /// Trans dipoles + oscillator strengths (K labels module).
    void print_state_dipoles(const SIProperties& props, int K,
                             const reks::SIResult& sir) const;

    /// Register SI-module data as wfn variables under "<BASE> K=<K>" (and bare
    /// base when primary_alias). Run sector s qualifies the names for s > 0
    /// ("REKS SECTOR s ..."); s == 0 keeps the byte-stable single-sector names.
    /// Arrays over the config axis (coefficients, Hamiltonian, overlap, diabatic
    /// 1-RDM) are exported from REKS_REPORT_LEVEL 3 up and carry no bare alias.
    void register_si_wfn_variables(int K, int s,
                                   const reks::SIResult& sir,
                                   const SIProperties& props,
                                   bool primary_alias,
                                   const reks::rdm::DiabaticRdm& rho_diab);

    /// Pure; returns empty when SI data is missing or the catalog carries no diabatic RDM.
    /// rho_adiab is the adiabatic 1-RDM tensor from rdm::adiabatic.
    reks::StateNaturalOrbitals compute_state_natural(
        const reks::studio::Cassette& cassette, const reks::SIResult& sir,
        const std::vector<double>& rho_adiab) const;

    /// Tables A (full SI eigenvector) + B (per-state FONs / NOs) + integrity block.
    void print_state_natural(const reks::StateNaturalOrbitals& nos, int K,
                             const reks::studio::Cassette& cassette,
                             const reks::SIResult& sir,
                             const reks::rdm::DiabaticRdm& rho_diab) const;

    /// Register NO data as wfn variables: " K=<K>" per module, " STATE=<k>" per state;
    /// also bare aliases when primary_alias. Run sector s qualifies the names for
    /// s > 0 ("REKS SECTOR s ..."); s == 0 keeps the byte-stable single-sector names.
    /// Exports every state nos carries, i.e. the reported states.
    void register_natural_wfn_variables(
        int K, int s,
        const reks::StateNaturalOrbitals& nos,
        bool primary_alias);

    void setup_potential() override;  ///< Create UV potential instead of RV
    void save_density_and_energy() override;
    double compute_initial_E() override;
    void damping_update(double damp) override;
    /// REKS-specific guess: READ propagates full Ca (incl. active virtuals)
    /// without touching nalphapi_/nalpha_; all other types delegate to HF::guess().
    void guess() override;
    /// S>0: symmetric occupation proxy nalphapi_ = nbetapi_ = Ncore_ + n_active_orbitals()
    /// (Ms=0 densities are spin-symmetric); S=0 delegates to the HF aufbau.
    void find_occupation() override;
    /// Makes C a full S-orthonormal basis on irrep h: S-metric MGS over the loaded cols
    /// 0..n_loaded-1 (dropping linearly dependent ones, survivors compacted to the front
    /// in order), then seeds the remaining cols from X. Returns the surviving loaded count.
    static int complete_orbital_basis(SharedMatrix C, SharedMatrix S, SharedMatrix X,
                                      int n_loaded, int h);
    void form_D() override;
    void form_G() override;
    void form_V() override;           ///< Use UV for spin-polarized XC
    void form_F() override;
    SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso) override;  ///< PGC: use ||g_orb|| for convergence at KKT boundary
    void form_C(double shift = 0.0) override;
    double compute_E() override;
    void finalize() override;         ///< Clean up UV potential
    void print_orbitals() override;   ///< 3-section print: Core / Active (with FON) / Virtual

    /// Global 0-based energy rank of every MO (irrep 0; REKS is C1).
    std::vector<int> mo_energy_rank() const;

    /// Energy-rank symmetry label ("5A","10A",...) for each active slot
    /// k = 0..n_active_orbitals()-1.
    std::vector<std::string> active_orbital_labels() const;

    /// (ij|kl) MO ERI tiles for the canonical pairs (k<=l) in rank1_pairs.
    /// Precondition: rank1_pairs canonical and strictly sorted.
    /// stage_bucket accumulates the JK wall time under key "jk_compute_active_eri".
    [[nodiscard]] reks::ActiveEriTiles compute_active_mo_eri(
        const std::vector<int>& active_mo_indices,
        const std::vector<std::pair<int,int>>& rank1_pairs,
        std::map<std::string, double>* stage_bucket);

    /// Two-electron Hessian diagonals gamma_aa and gamma_ca, 2 sum_L C_L A2e(L,pq)
    /// over the pair ERIs (pp|kk) and (pk|pk). Those are the MO diagonals of the
    /// J_act/K_act Fock slots, built this iteration by build_sa_focks_MO and kept on
    /// reel_.slot_diag; a wcombine JK leaves no plain K to read and throws.
    void compute_gamma_vectors(
        std::vector<double>& gamma_aa,
        std::vector<double>& gamma_ca) const;

    /// Asymmetric generalized Fock for orbital gradient G[i,a] = 2*(F_gen[i,a]-F_gen[a,i]):
    ///   F_gen[p,q] = sum_L C_L * (n_alpha_L[p] * Fa_L[p,q] + n_beta_L[p] * Fb_L[p,q]).
    /// Asymmetric because row occupation weights differ across core/active/virt.
    SharedMatrix build_generalized_fock() const;

    /// Build reel_.F_acc_MO = transform of Sum_L C_L (Fa_AO_L + Fb_AO_L), the
    /// C_L-weighted accumulated MO Fock. core_rows_only limits the update to the
    /// first Ncore_ rows; other rows keep their previous values.
    void build_F_acc_MO(bool core_rows_only = false);

   public:
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
         Options& options, std::shared_ptr<PSIO> psio);
    ~REKS() override;

    /// Create C1 symmetry deep copy
    std::shared_ptr<REKS> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    std::shared_ptr<VBase> V_potential() const override { return rv_potential_; }

    /// Seed the generation-0 (n) active-space FON from a saved wavefunction
    /// (READ-guess path). Matrix must have sa_cassette_.geminals_active_gen(0).size()
    /// entries; values within 1e-8 of {0,2} are clamped inside to avoid boundary
    /// traps. When unset (e.g. RHF->REKS cross-read or pre-FON .npy), form_C iter=0
    /// diagonalizes without an active-active freeze.
    void guess_fon(SharedMatrix fon, int s = 0) {
        if (s >= static_cast<int>(guess_fon_.size())) guess_fon_.resize(s + 1);
        guess_fon_[s] = fon;
    }

    /// Seed a higher-generation (gen >= 1: m, u, ...) active FON of run sector s from
    /// a saved wavefunction. Size must equal sa_cassette_.geminals_active_gen(gen).size().
    void guess_fon_upper(int gen, SharedMatrix fon, int s = 0) {
        if (s >= static_cast<int>(guess_fon_upper_.size())) guess_fon_upper_.resize(s + 1);
        guess_fon_upper_[s][gen] = fon;
    }

    /// Total microstate count (Catalog basis).
    [[nodiscard]] int n_microstates() const {
        return sa_cassette_.n_catalog_microstates();
    }

    /// SA microstate count = cardinality of sa_cassette_'s F_MO build set.
    [[nodiscard]] int n_sa_microstates() const {
        return static_cast<int>(sa_cassette_.microstates().size());
    }

    [[nodiscard]] double get_microstate_energy(int L) const {
        return reel_.E_L[L];
    }

    /// Get f_interp value for n-geminal g of run sector s.
    [[nodiscard]] double get_f_value(int g = 0, int s = 0) const {
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        const auto& layer = reel_.fon_state[s].layers[0];
        if (g < 0 || g >= static_cast<int>(layer.size())) return 0.0;
        const auto& fon = layer[g];
        return reks::f_interp(fon.p * fon.q);
    }

    /// Get SI Hamiltonian matrix (pre-diagonalization)
    [[nodiscard]] SharedMatrix SI_hamiltonian() const;
    /// Get SI eigenvalues (adiabatic state energies)
    [[nodiscard]] SharedVector SI_energies() const;
    /// Get SI eigenvectors (rows = states)
    [[nodiscard]] SharedMatrix SI_coefficients() const;
    [[nodiscard]] SharedMatrix SI_overlap() const;

    /// Bonding-orbital FON of geminal g in generation gen (0=n, 1=m, 2=u, ...) of
    /// run sector s.
    [[nodiscard]] double get_fon_p(int g, int gen = 0, int s = 0) const {
        if (gen < 0 || gen >= reks::studio::kMaxGen) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        const auto& layer = reel_.fon_state[s].layers[gen];
        return (g >= 0 && g < static_cast<int>(layer.size())) ? layer[g].p : 0.0;
    }

    /// Antibonding-orbital FON of geminal g in generation gen (0=n, 1=m, 2=u, ...) of
    /// run sector s.
    [[nodiscard]] double get_fon_q(int g, int gen = 0, int s = 0) const {
        if (gen < 0 || gen >= reks::studio::kMaxGen) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        const auto& layer = reel_.fon_state[s].layers[gen];
        return (g >= 0 && g < static_cast<int>(layer.size())) ? layer[g].q : 0.0;
    }

    /// Per-spin MO occupations (alpha=beta) with run sector s's n-FON written onto
    /// its active n-geminal orbitals.
    [[nodiscard]] SharedVector fon_occupation(int s = 0) const;

    [[nodiscard]] double get_lagrangian(int idx) const {
        const auto& l = reel_.lagrangians;
        return (idx >= 0 && idx < static_cast<int>(l.size())) ? l[idx] : 0.0;
    }

    /// Configuration energy E_K = sum_L C_L(K, fons) * E_L[L], evaluated against run
    /// sector s's FON snapshot; C_L, config_microstate_writes and the microstate
    /// count are catalog-global.
    [[nodiscard]] double get_config_energy(int K, int s = 0) const {
        if (K < 0 || K >= sa_cassette_.n_configs()) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        std::vector<double> C_K(sa_cassette_.n_catalog_microstates(), 0.0);
        sa_cassette_.compute_C_L(K, reel_.fon_state[s], C_K);
        double E = 0.0;
        for (int L : sa_cassette_.config_microstate_writes(K))
            E += C_K[L] * reel_.E_L[L];
        return E;
    }

    /// Public per-sector FON vector accessors (active bonding p values of run sector s,
    /// generation gen). set mutates the live snapshot without re-running SCF.
    [[nodiscard]] std::vector<double> get_fon_vector(int s, int gen) const {
        return fon_vector(s, gen);
    }
    void put_fon_vector(int s, int gen, const std::vector<double>& fons) {
        set_fon_vector(s, gen, fons);
    }

    /// Total geminal slots in the catalog (generation-independent).
    [[nodiscard]] int n_pairs() const { return sa_cassette_.n_geminals(); }

    /// Scheme tag of geminal slot g (0 = SCF-optimized, >0 = post-SCF).
    [[nodiscard]] int get_n_scheme(int g) const {
        return (g >= 0 && g < sa_cassette_.n_geminals())
            ? sa_cassette_.geminal_templates()[g].scheme : -1;
    }

    /// Total configuration count (variant-global codegen constant).
    [[nodiscard]] int n_configs() const { return sa_cassette_.n_configs(); }

    /// Number of Lagrangian elements (variant-global).
    [[nodiscard]] int n_lagrangians() const {
        return sa_cassette_.n_lagrangian_pairs();
    }

    /// Dimension of the SI Hamiltonian (variant config count, HS excluded).
    [[nodiscard]] int n_si_states() const { return sa_cassette_.n_configs(); }

    /// Number of doubly-occupied core orbitals (MOs below the active space).
    [[nodiscard]] int get_Ncore() const { return Ncore_; }

    [[nodiscard]] int n_active_orbitals() const {
        return sa_cassette_.n_active_orbitals();
    }

    /// MO indices of active orbitals (length n_active_orbitals).
    [[nodiscard]] std::vector<int> get_active_mo_indices() const {
        return active_mo_indices_;
    }

    /// Active accelerators this SCF iter for `@REKS iter ...` status string;
    /// labels in {"CFM-GVB-DIIS","ORB-GVB-DIIS","TRAH","PLAIN","kDIIS"}, empty at iter<=0.
    [[nodiscard]] std::vector<std::string> get_iter_accel_labels() const;

    /// True while the GVB-DIIS joint (dE, rms) gate still needs further consecutive iterations
    /// to reach REKS_GVB_GATE_STREAK; false whenever the controller is inactive.
    [[nodiscard]] bool gvb_gate_confirm_pending() const {
        return ctrl_state_.is_active() && ctrl_state_.gate_streak() < gvb_gate_streak_req_;
    }

};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
