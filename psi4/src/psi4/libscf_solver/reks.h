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

#include "rhf.h"
#include "reks_active_space.h"
#include <vector>
#include <memory>

namespace psi {
namespace scf {

/// FON optimization phase for adaptive solver (FON_SOLVER=4)
enum class FONPhase {
    STABLE,         ///< Hessian positive definite, well away from transition
    PRE_TRANSITION, ///< Hessian still positive but eigenvalues decreasing (approaching transition)
    TRANSITION      ///< Hessian indefinite (at or past saddle point)
};

/**
 * @class REKS
 * @brief Restricted Ensemble Kohn-Sham (REKS) SCF engine
 *
 * This class implements the REKS SCF procedure using a pluggable
 * active space definition (REKSActiveSpace). The SCF algorithm is
 * the same for all REKS variants (2,2), (4,4), etc.
 *
 * Architecture:
 * - Level 1: reks_math.h - pure math functions (f_interp, etc.)
 * - Level 2: reks_active_space.h - active space interface + implementations
 * - Level 3: reks.h/cc - SCF engine (this class)
 *
 * Current Implementation: REKS(2,2) via REKS22Space
 * Future: REKS(4,4) via REKS44Space
 *
 * References:
 * - Filatov, M. "Note on SA-REKS, SSR, CP-REKS implementation" (2024)
 * - Filatov, M.; Shaik, S. Chem. Phys. Lett. 1999, 304, 429-437
 * - Filatov, M. et al. J. Chem. Phys. 147, 064104 (2017)
 */
class REKS : public RHF {
   protected:
    // ========================================================================
    // Active Space (polymorphic)
    // ========================================================================

    /// Active space definition (REKS22Space, REKS44Space, etc.)
    std::unique_ptr<reks::REKSActiveSpace> active_space_;

    // ========================================================================
    // Core Orbital Information
    // ========================================================================

    int Ncore_ = -1;                       ///< Number of core (doubly occupied) orbitals
    std::vector<int> active_mo_indices_;   ///< Active orbital MO indices (generalized from active_r_, active_s_)

    // ========================================================================
    // Base Density Matrices (generalized for REKS(N,M))
    // ========================================================================

    /// Base densities indexed by occupation bitmask
    /// For REKS(2,2): 4 densities (patterns 0-3)
    /// For REKS(4,4): 16 densities (patterns 0-15)
    std::vector<SharedMatrix> base_densities_;

    // ========================================================================
    // Microstate Data (dynamically sized)
    // ========================================================================

    std::vector<SharedMatrix> D_alpha_micro_;  ///< Alpha density for each microstate
    std::vector<SharedMatrix> D_beta_micro_;   ///< Beta density for each microstate
    std::vector<SharedMatrix> F_alpha_micro_;  ///< Alpha Fock for each microstate
    std::vector<SharedMatrix> F_beta_micro_;   ///< Beta Fock for each microstate
    std::vector<double> E_micro_;              ///< Energy for each microstate
    std::vector<double> C_L_;                  ///< Weighting factors for each microstate

    // ========================================================================
    // Coupling Fock Matrix
    // ========================================================================

    SharedMatrix F_reks_;          ///< Assembled coupling Fock matrix (AO basis)
    SharedMatrix F_reks_MO_;       ///< F_reks in MO basis
    double Wrs_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for SI (REKS(2,2): r,s; REKS(4,4): a,d)
    double Wbc_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for (b,c) pair (REKS(4,4) only)
    double Wac_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for (a,c) inter-pair (REKS(4,4) 9SI only)
    double Wbd_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for (b,d) inter-pair (REKS(4,4) 9SI only)

    // ========================================================================
    // Pre-allocated Work Matrices (for HPC efficiency)
    // ========================================================================

    /// Work matrices for C columns for each base density pattern
    std::vector<SharedMatrix> C_base_;     ///< Generalized from C_D00_, C_D10_, etc.
    SharedMatrix J_work_;                   ///< Work matrix: temporary for J_total
    SharedMatrix Temp_work_;                ///< Work matrix: temporary for AO-to-MO transforms
    std::vector<SharedMatrix> F_alpha_MO_;  ///< Pre-allocated MO-basis alpha Fock
    std::vector<SharedMatrix> F_beta_MO_;   ///< Pre-allocated MO-basis beta Fock

    // ========================================================================
    // XC (Exchange-Correlation) Support
    // ========================================================================

    std::vector<SharedMatrix> D_total_micro_;  ///< Total density (alpha+beta) for each microstate
    std::vector<SharedMatrix> V_xc_micro_;     ///< XC potential for each microstate
    std::vector<double> E_xc_micro_;           ///< XC energy for each microstate

    // ========================================================================
    // Debug Level
    // ========================================================================

    int reks_debug_ = 0;           ///< 0=none, 1=energies, 2=matrices, 3=all

    // ========================================================================
    // Orbital Localization for Active Space
    // ========================================================================

    std::string localization_type_;    ///< "NONE" or "BOYS"
    bool localization_done_ = false;   ///< Track if localization has been applied
    int localization_freeze_iter_ = 0; ///< Freeze localization after this iteration (0=never freeze)

    /// Apply Boys localization to active space orbitals
    void localize_active_space();

    /// Reorder active orbitals for REKS(4,4) GVB pairs based on centroids
    void reorder_active_orbitals_for_gvb_pairs();

    // ========================================================================
    // FON Optimization Control (Two-Phase SCF)
    // ========================================================================
    //
    // Phase 1: Freeze FON for initial iterations to let orbitals converge.
    //          DIIS works correctly with consistent C_L weights.
    // Phase 2: Enable FON optimization with step limiting.
    //
    // This is the standard approach for multi-configurational SCF methods
    // (CASSCF, RASSCF) and ensures proper DIIS convergence.

    /// Number of SCF iterations to freeze FON at initial value (Phase 1)
    /// At n_r = 1.0 (symmetric), all weights are non-negative, so no freeze is needed.
    /// GAMESS optimizes FON from iteration 1, so we match that behavior.
    int fon_freeze_iters_ = 0;

    /// Maximum FON change per SCF iteration (Phase 2 damping)
    /// Set to 2.0 to effectively disable (was 0.3)
    /// This allows FON to quickly escape from saddle point n=1.0
    double fon_max_delta_ = 2.0;

    /// Previous FON values for delta limiting
    double prev_n_r_ = 1.0;        ///< REKS(2,2): previous n_r (starts symmetric)
    double prev_n_a_ = 1.0;        ///< REKS(4,4): previous n_a (starts symmetric)
    double prev_n_b_ = 1.0;        ///< REKS(4,4): previous n_b (starts symmetric)

    /// Flag: reset DIIS when transitioning from Phase 1 to Phase 2
    bool fon_phase_transition_ = false;

    /// Enable line search with golden ratio for REKS(4,4) FON optimization.
    /// Prevents oscillations when Newton-Raphson step increases energy.
    bool reks_line_search_ = true;

    /// Threshold for "shake it up" - if FON is within this of 1.0, reset to 0.5.
    /// Prevents FON locking at the degenerate point n=1.0. Set to 0.0 to disable.
    double reks_shake_threshold_ = 1e-6;

    /// Coupling threshold for REKS(4,4). If coupling element is below this,
    /// fix corresponding FON at boundary (2.0). Set to 0.0 to disable (default).
    double reks_coupling_threshold_ = 0.0;

    /// Initial FON value for n_a (REKS(4,4) geometry scan continuation).
    /// Set to -1.0 to use HAM-DIAG initialization (default).
    /// For geometry scans, set to previous geometry's FON for smooth PES.
    double reks_initial_na_ = -1.0;

    /// Initial FON value for n_b (REKS(4,4) geometry scan continuation).
    /// Set to -1.0 to use HAM-DIAG initialization (default).
    double reks_initial_nb_ = -1.0;

    /// FON history for oscillation detection (last 3 values)
    double prev_prev_n_a_ = 1.0;   ///< n_a from 2 iterations ago
    double prev_prev_n_b_ = 1.0;   ///< n_b from 2 iterations ago

    /// Counter for consecutive oscillations
    int fon_oscillation_count_ = 0;

    /// Flag: FONs are frozen due to oscillation detection
    bool fon_oscillation_frozen_ = false;

    /// Flag: REKS(4,4) FONs have been initialized from Hamiltonian diagonalization
    bool fon_initialized_44_ = false;

    // ========================================================================
    // Trust-Region FON Optimization (REKS(4,4))
    // ========================================================================
    //
    // Trust-Region Augmented Hessian (TRAH) method for robust FON optimization.
    // Handles indefinite Hessian at saddle points and near boundaries.
    //
    // Reference: J. Chem. Phys. 156, 204104 (2022)

    /// Flag: use Trust-Region method instead of Newton-Raphson for REKS(4,4)
    bool use_trust_region_44_ = false;

    /// Current trust region radius (persistent across SCF iterations)
    double trust_radius_44_ = reks::constants::TR_DELTA_INIT;

    // ========================================================================
    // Multi-Start FON Optimization
    // ========================================================================

    /// Flag: enable multi-start FON optimization
    bool reks_multi_start_ = false;

    /// Branch selection criterion: "ENERGY", "DIABATIC", or "AUTO"
    std::string reks_branch_criterion_ = "AUTO";

    /// Energy tolerance for AUTO branch selection (Ha)
    double reks_energy_tolerance_ = 0.01;

    /// Report all solutions (debug)
    bool reks_report_all_solutions_ = false;

    /// FON solver level: 1=Fast(NR), 2=Balanced(TR), 3=Robust(TR+MS), 4=Adaptive
    int fon_solver_level_ = 2;

    // ========================================================================
    // Adaptive FON Solver (FON_SOLVER=4)
    // ========================================================================

    /// Previous Hessian determinant for trend detection
    double det_H_prev_ = 1.0;

    /// Flag: FON values should be symmetric (n_a = n_b) based on molecular symmetry
    bool is_symmetric_fon_ = true;

    /// Flag: enforce FON symmetry from user option
    bool reks_enforce_symmetry_ = false;

    // ========================================================================
    // multi_scf() Support: Precomputed J/K Matrices
    // ========================================================================

    /// Precomputed J matrices for base densities (from multi_scf shared JK)
    std::vector<SharedMatrix> precomputed_J_base_;

    /// Precomputed K matrices for base densities (from multi_scf shared JK)
    std::vector<SharedMatrix> precomputed_K_base_;

    /// Precomputed wK matrices for base densities (from multi_scf shared JK, LRC functionals)
    std::vector<SharedMatrix> precomputed_wK_base_;

    /// Flag: use precomputed J/K instead of calling jk_->compute()
    bool use_precomputed_jk_base_ = false;

    // ========================================================================
    // Protected Methods
    // ========================================================================

    /// REKS-specific initialization
    void reks_common_init();

    /// Allocate all REKS-specific matrices
    void allocate_reks_matrices();

    /// Build base density matrices D00, D10, D01, D11
    void build_base_densities();

    /// Map base densities to microstate alpha/beta densities
    void build_microstate_densities();

    /// Build UHF-like Fock matrices for each microstate
    void build_microstate_focks();

    /// Compute microstate energies E_L
    void compute_microstate_energies();

    /// Newton-Raphson FON optimization (dispatcher)
    void rex_solver();

    /// REKS(2,2) FON optimization: single variable n_r
    void rex_solver_22();

    /// REKS(4,4) FON optimization: two variables (n_a, n_b)
    void rex_solver_44();

    /// Trust-Region FON optimization for REKS(4,4)
    /// Handles indefinite Hessian robustly via secular equation
    void trust_region_fon_optimizer_44();

    /// Adaptive Trust-Region + Multi-Start for REKS(4,4) (FON_SOLVER=3)
    /// Uses Trust-Region by default, triggers Multi-Start only when needed
    /// (indefinite Hessian, first iteration, or large FON change)
    void trust_region_fon_optimizer_44_adaptive();

    /// Adaptive FON solver for REKS(4,4) (FON_SOLVER=4)
    /// Continuation-first approach with early transition detection.
    /// Auto-probes for diradical solutions when Hessian indicates transition.
    void adaptive_fon_solver_44();

    /// Detect current FON optimization phase based on Hessian analysis
    /// @param grad Gradient vector [dE/dn_a, dE/dn_b]
    /// @param hess Hessian matrix (2x2 as flat vector [H_aa, H_ab, H_ba, H_bb])
    /// @return Current phase: STABLE, PRE_TRANSITION, or TRANSITION
    FONPhase detect_fon_phase(const std::vector<double>& grad, const std::vector<double>& hess);

    /// Multi-start FON optimization for REKS(4,4)
    /// Tries multiple starting points and selects best based on criterion
    void rex_solver_44_multi_start();

    /// Optimize FON from a single starting point using Trust-Region
    /// @param start Starting point (n_a, n_b)
    /// @param start_id ID of this starting point (for tracking)
    /// @return Converged FON solution
    reks::FONSolution optimize_from_start_44(const reks::FONStartPoint& start, int start_id);

    /// Select best solution from list based on criterion
    reks::FONSolution select_solution(
        const std::vector<reks::FONSolution>& solutions,
        double prev_na, double prev_nb,
        const std::string& criterion,
        double energy_tol);

    /// Solve Trust-Region subproblem: min g·p + ½p·H·p s.t. ||p|| ≤ Δ
    /// @param g Gradient vector [g_a, g_b]
    /// @param H Hessian matrix (2x2)
    /// @param Delta Trust region radius
    /// @return Step vector [δn_a, δn_b]
    std::vector<double> solve_trust_region_subproblem_2d(
        const std::vector<double>& g, const std::vector<std::vector<double>>& H, double Delta);

    /// Solve secular equation for boundary-constrained step
    /// When Newton step is outside trust region or Hessian is indefinite
    /// @param g1, g2 Gradient in eigenspace
    /// @param lam1, lam2 Hessian eigenvalues
    /// @param v1, v2 Hessian eigenvectors (2-element each)
    /// @param Delta Trust region radius
    /// @return Step vector in original space [δn_a, δn_b]
    std::vector<double> solve_boundary_constrained_2d(
        double g1, double g2, double lam1, double lam2,
        const std::vector<double>& v1, const std::vector<double>& v2, double Delta);

    /// Initialize FON from 4×4 CI Hamiltonian diagonalization (REKS(4,4) only)
    /// This provides a physically motivated initial guess instead of saddle point n=1.0
    void initialize_fon_from_hamiltonian_44();

    /// Compute weighting factors C_L (delegates to active_space_)
    void compute_weighting_factors();

    /// Add microstate L contribution to coupling Fock matrix (dispatcher)
    void fock_micro_to_macro(int L, double Cl, double** Fa, double** Fb);

    /// REKS(2,2) version of coupling Fock assembly
    void fock_micro_to_macro_22(int L, double Cl, double** Fa, double** Fb);

    /// REKS(4,4) version of coupling Fock assembly
    void fock_micro_to_macro_44(int L, double Cl, double** Fa, double** Fb);

    /// Debug output: header
    void print_reks_header() const;
    /// Debug output: iteration info
    void print_reks_iteration(int iter, double E_tot, double dE, double dD) const;
    /// Debug output: microstate energies
    void print_microstate_energies() const;
    /// Debug output: FON info
    void print_fon_info() const;
    /// Print SI-SA state energies (called after SCF convergence)
    void print_SI_energies() const;

    // ========================================================================
    // SCF Method Overrides
    // ========================================================================

    void form_D() override;
    void form_G() override;
    // Note: form_V() NOT overridden - use RHF::form_V() directly
    void form_F() override;
    void form_C(double shift = 0.0) override;
    double compute_E() override;

   public:
    // ========================================================================
    // Constructors
    // ========================================================================

    /**
     * @brief Construct REKS with default REKS(2,2) active space
     * @param ref_wfn Reference wavefunction
     * @param functional DFT functional (or HF if nullptr)
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);

    /**
     * @brief Construct REKS with explicit options and PSIO
     * @param ref_wfn Reference wavefunction
     * @param functional DFT functional
     * @param options Psi4 options object
     * @param psio PSIO object for file I/O
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
         Options& options, std::shared_ptr<PSIO> psio);

    /// Destructor
    ~REKS() override;

    // ========================================================================
    // Public Interface
    // ========================================================================

    // ========================================================================
    // multi_scf() Interface (standalone, no longer overrides)
    // ========================================================================

    /// Number of electronic states for multi_scf JK batching
    /// REKS requires n_base + 1 J/K matrices (all base densities + RHF density)
    int n_states() const;

    /// Returns orbital matrices for multi_scf shared JK computation
    /// Returns all C_base_ matrices + C_occ (RHF occupied orbitals)
    std::vector<SharedMatrix> get_orbital_matrices() const;

    /// Sets precomputed J/K matrices from multi_scf shared JK
    /// @param J_list: J matrices for each base density + RHF
    /// @param K_list: K matrices for each base density + RHF
    /// @param wK_list: long-range K matrices (empty for HF)
    void set_jk_matrices(const std::vector<SharedMatrix>& J_list,
                         const std::vector<SharedMatrix>& K_list,
                         const std::vector<SharedMatrix>& wK_list);

    /// Create C1 symmetry deep copy
    std::shared_ptr<REKS> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    // ========================================================================
    // Accessors for Testing & Analysis
    // ========================================================================

    /// Get active space definition
    [[nodiscard]] const reks::REKSActiveSpace& active_space() const {
        return *active_space_;
    }

    /// Get number of microstates
    [[nodiscard]] int n_microstates() const {
        return active_space_->n_microstates();
    }

    /// Get microstate energy for index L
    [[nodiscard]] double get_microstate_energy(int L) const {
        return E_micro_[L];
    }

    /// Get weighting factor for microstate L
    [[nodiscard]] double get_C_L(int L) const {
        return C_L_[L];
    }

    // ========================================================================
    // FON Accessors
    // ========================================================================

    // --- REKS(2,2) Accessors (pair 0: r,s) ---

    /// Get FON for orbital r (REKS(2,2), or n_a for REKS(4,4))
    [[nodiscard]] double get_n_r() const {
        return active_space_->pair(0).fon_p;
    }

    /// Get FON for orbital s (REKS(2,2), or n_d for REKS(4,4))
    [[nodiscard]] double get_n_s() const {
        return active_space_->pair(0).fon_q;
    }

    /// Get f_interp value at current FON (only valid for REKS(2,2))
    [[nodiscard]] double get_f_value() const {
        const auto& p = active_space_->pair(0);
        return reks::f_interp(p.fon_p * p.fon_q);
    }

    // --- REKS(4,4) Accessors (pair 0: a,d and pair 1: b,c) ---

    /// Get FON for orbital a (pair 0, REKS(4,4))
    [[nodiscard]] double get_n_a() const {
        return active_space_->pair(0).fon_p;
    }

    /// Get FON for orbital d (pair 0, REKS(4,4))
    [[nodiscard]] double get_n_d() const {
        return active_space_->pair(0).fon_q;
    }

    /// Get FON for orbital b (pair 1, REKS(4,4) only)
    /// @return n_b if REKS(4,4), 0.0 otherwise
    [[nodiscard]] double get_n_b() const {
        if (active_space_->n_pairs() >= 2) {
            return active_space_->pair(1).fon_p;
        }
        return 0.0;
    }

    /// Get FON for orbital c (pair 1, REKS(4,4) only)
    /// @return n_c if REKS(4,4), 0.0 otherwise
    [[nodiscard]] double get_n_c() const {
        if (active_space_->n_pairs() >= 2) {
            return active_space_->pair(1).fon_q;
        }
        return 0.0;
    }

    // ========================================================================
    // Response Lagrangian Accessors
    // ========================================================================

    /// Get W_ad Lagrange multiplier (off-diagonal SI coupling for a,d pair)
    /// For REKS(2,2) this is Wrs
    [[nodiscard]] double get_Wrs() const {
        return Wrs_lagr_;
    }

    /// Get W_ad Lagrange multiplier (same as Wrs, alternative name for REKS(4,4))
    [[nodiscard]] double get_W_ad() const {
        return Wrs_lagr_;
    }

    /// Get W_bc Lagrange multiplier (off-diagonal SI coupling for b,c pair)
    /// Only valid for REKS(4,4)
    [[nodiscard]] double get_Wbc() const {
        return Wbc_lagr_;
    }

    /// Get W_bc Lagrange multiplier (alternative name)
    [[nodiscard]] double get_W_bc() const {
        return Wbc_lagr_;
    }

    /// Get W_ac Lagrange multiplier (inter-pair SI coupling for a,c)
    /// Only valid for REKS(4,4) 9SI
    [[nodiscard]] double get_W_ac() const {
        return Wac_lagr_;
    }

    /// Get W_bd Lagrange multiplier (inter-pair SI coupling for b,d)
    /// Only valid for REKS(4,4) 9SI
    [[nodiscard]] double get_W_bd() const {
        return Wbd_lagr_;
    }

    // ========================================================================
    // SI-SA-REKS State Energies
    // ========================================================================

    /// Get SI state energy for state index (0 = ground, 1 = S1, etc.)
    /// @param state State index (0 to n_si_states-1)
    /// @param n_si_states Number of SI states (2 or 3 for 3SI, up to 9 for 9SI)
    /// @return State energy in Hartree
    [[nodiscard]] double get_SI_energy(int state, int n_si_states = 2) const {
        // For REKS(4,4)+, use compute_SI_energies_44 with both Lagrangians
        // Check n_pairs() >= 2 to detect REKS(4,4) or higher
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, n_si_states);
            if (state >= 0 && state < static_cast<int>(si.energies.size())) {
                return si.energies[state];
            }
            return 0.0;
        }
        // For REKS(2,2), use original method
        auto si = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, n_si_states);
        if (state >= 0 && state < static_cast<int>(si.energies.size())) {
            return si.energies[state];
        }
        return 0.0;
    }

    /// Get 9SI state energy for REKS(4,4)
    /// Uses all 4 Lagrangians (W_ad, W_bc, W_ac, W_bd) for full 9x9 Hamiltonian
    /// @param state State index (0 to n_si_states-1)
    /// @param n_si_states Number of SI states (up to 9)
    /// @return State energy in Hartree
    [[nodiscard]] double get_SI_energy_9x9(int state, int n_si_states = 9) const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_9x9(
                E_micro_, Wrs_lagr_, Wbc_lagr_, Wac_lagr_, Wbd_lagr_, n_si_states);
            if (state >= 0 && state < static_cast<int>(si.energies.size())) {
                return si.energies[state];
            }
        }
        return 0.0;
    }

    /// Get PPS diagonal energy (H_11)
    [[nodiscard]] double get_E_PPS() const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);
            return si.E_PPS;
        }
        auto si = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, 2);
        return si.E_PPS;
    }

    /// Get OSS diagonal energy (H_22)
    /// For REKS(4,4), returns E_OSS1
    [[nodiscard]] double get_E_OSS() const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);
            return si.E_OSS;  // E_OSS1 for REKS(4,4)
        }
        auto si = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, 2);
        return si.E_OSS;
    }

    /// Get DES diagonal energy (H_33)
    /// For REKS(4,4), returns E_OSS2
    [[nodiscard]] double get_E_DES() const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);
            return si.E_DES;  // E_OSS2 for REKS(4,4)
        }
        auto si = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, 3);
        return si.E_DES;
    }

    // --- REKS(4,4) Specific State Energy Accessors ---

    /// Get E_OSS1 diagonal energy (REKS(4,4) only)
    [[nodiscard]] double get_E_OSS1() const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);
            return si.E_OSS;
        }
        return 0.0;
    }

    /// Get E_OSS2 diagonal energy (REKS(4,4) only)
    [[nodiscard]] double get_E_OSS2() const {
        if (active_space_->n_pairs() >= 2) {
            auto si = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);
            return si.E_DES;
        }
        return 0.0;
    }

    /// Get triplet energy (from microstate 3 for REKS(2,2))
    [[nodiscard]] double get_E_triplet() const {
        if (E_micro_.size() > 3) {
            return E_micro_[3];
        }
        return 0.0;
    }

    // --- REKS(4,4) Additional Configuration Energies ---

    /// Get DOSS (Double Open-Shell Singlet) configuration energy
    /// Only valid for REKS(4,4). Uses fixed FONs: n'_a=n'_d=n'_b=n'_c=1
    [[nodiscard]] double get_E_DOSS() const {
        if (active_space_->n_pairs() >= 2) {
            // Cast to REKS44Space to access DOSS energy function
            auto* space44 = dynamic_cast<const reks::REKS44Space*>(active_space_.get());
            if (space44) {
                return space44->compute_energy_DOSS(E_micro_);
            }
        }
        return 0.0;
    }

    /// Get DSPS (Double Single-Pair Singlet) configuration energy
    /// Only valid for REKS(4,4). Uses fixed FONs: n'_a=n'_d=n'_b=n'_c=1
    [[nodiscard]] double get_E_DSPS() const {
        if (active_space_->n_pairs() >= 2) {
            // Cast to REKS44Space to access DSPS energy function
            auto* space44 = dynamic_cast<const reks::REKS44Space*>(active_space_.get());
            if (space44) {
                return space44->compute_energy_DSPS(E_micro_);
            }
        }
        return 0.0;
    }

    /// Get OSS3 (Inter-pair Open-Shell Singlet) configuration energy
    /// Only valid for REKS(4,4). INTER-PAIR pairing: (a,c) OSS, (b,d) GVB
    [[nodiscard]] double get_E_OSS3() const {
        if (active_space_->n_pairs() >= 2) {
            return active_space_->compute_energy_OSS3(E_micro_);
        }
        return 0.0;
    }

    /// Get OSS4 (Inter-pair Open-Shell Singlet) configuration energy
    /// Only valid for REKS(4,4). INTER-PAIR pairing: (b,d) OSS, (a,c) GVB
    [[nodiscard]] double get_E_OSS4() const {
        if (active_space_->n_pairs() >= 2) {
            return active_space_->compute_energy_OSS4(E_micro_);
        }
        return 0.0;
    }

    /// Get DES1 (Doubly Excited Singlet type 1) configuration energy
    /// Only valid for REKS(4,4). Standard pairing, doubly excited
    [[nodiscard]] double get_E_DES1() const {
        if (active_space_->n_pairs() >= 2) {
            auto* space44 = dynamic_cast<const reks::REKS44Space*>(active_space_.get());
            if (space44) {
                return space44->compute_energy_DES1(E_micro_);
            }
        }
        return 0.0;
    }

    /// Get DES2 (Doubly Excited Singlet type 2) configuration energy
    /// Only valid for REKS(4,4). Inter-pair pairing, doubly excited
    [[nodiscard]] double get_E_DES2() const {
        if (active_space_->n_pairs() >= 2) {
            auto* space44 = dynamic_cast<const reks::REKS44Space*>(active_space_.get());
            if (space44) {
                return space44->compute_energy_DES2(E_micro_);
            }
        }
        return 0.0;
    }
};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
