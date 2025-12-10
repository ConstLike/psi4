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

/**
 * @file reks.cc
 * @brief Implementation of Restricted Ensemble Kohn-Sham (REKS) method.
 *
 * Level 3: SCF Engine using pluggable REKSActiveSpace.
 *
 * Architecture:
 * - Level 1: reks_math.h - pure math functions (f_interp, GVBPair, Microstate)
 * - Level 2: reks_active_space.h/cc - abstract interface + REKS22Space, REKS44Space
 * - Level 3: reks.h/cc - SCF engine (this file)
 *
 * Current Implementation: REKS(2,2) via REKS22Space
 * Future: REKS(4,4) via REKS44Space
 *
 * Reference: Filatov, M. "Note on SA-REKS, SSR, CP-REKS implementation" (2024)
 */

#include "reks.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"  // For VBase (XC potential)
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"  // For C_DGEMM, C_DGER (BLAS)
#include "psi4/libmints/local.h"  // For Boys localization
#include "psi4/libmints/mintshelper.h"  // For dipole integrals

#include <cmath>
#include <algorithm>
#include <array>

namespace psi {
namespace scf {

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : RHF(ref_wfn, func) {
    reks_common_init();
}

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
           Options& options, std::shared_ptr<PSIO> psio)
    : RHF(ref_wfn, func, options, psio) {
    reks_common_init();
}

REKS::~REKS() {}

void REKS::reks_common_init() {
    // Update wavefunction name to REKS
    name_ = "REKS";

    // Use existing DEBUG option for REKS debug output
    reks_debug_ = options_.get_int("DEBUG");

    // Get localization type for active space orbitals
    localization_type_ = options_.get_str("REKS_LOCALIZATION");
    localization_done_ = false;
    localization_freeze_iter_ = options_.get_int("REKS_LOCALIZATION_FREEZE");

    // Get FON optimizer options
    reks_line_search_ = options_.get_bool("REKS_LINE_SEARCH");
    reks_shake_threshold_ = options_.get_double("REKS_SHAKE_THRESHOLD");
    reks_coupling_threshold_ = options_.get_double("REKS_COUPLING_THRESHOLD");
    fon_max_delta_ = options_.get_double("REKS_FON_MAX_DELTA");
    reks_initial_na_ = options_.get_double("REKS_INITIAL_NA");
    reks_initial_nb_ = options_.get_double("REKS_INITIAL_NB");
    use_trust_region_44_ = options_.get_bool("REKS_USE_TRUST_REGION");

    // Multi-Start FON Optimization options
    reks_multi_start_ = options_.get_bool("REKS_MULTI_START");
    reks_branch_criterion_ = options_.get_str("REKS_BRANCH_CRITERION");
    reks_energy_tolerance_ = options_.get_double("REKS_ENERGY_TOLERANCE");
    reks_report_all_solutions_ = options_.get_bool("REKS_REPORT_ALL_SOLUTIONS");

    // FON Solver Level (simplified interface)
    fon_solver_level_ = options_.get_int("FON_SOLVER");
    if (fon_solver_level_ < 1) fon_solver_level_ = 1;
    if (fon_solver_level_ > 4) fon_solver_level_ = 4;

    // Adaptive solver options
    reks_enforce_symmetry_ = options_.get_bool("REKS_ENFORCE_SYMMETRY");

    // Set trust region and multi-start flags based on FON_SOLVER level
    // Level 1: NR only  -> no trust region, no multi-start
    // Level 2: TR only  -> trust region, no multi-start, localization_freeze=10
    // Level 3: TR + MS  -> trust region + multi-start, localization_freeze=10
    // Level 4: Adaptive -> trust region + continuation, no multi-start (uses probing)
    if (fon_solver_level_ >= 2) {
        use_trust_region_44_ = true;
        // Auto-enable localization freeze for FON_SOLVER>=2 if not explicitly set
        // This prevents SCF oscillations from repeated localization
        if (localization_freeze_iter_ == 0) {
            localization_freeze_iter_ = 10;
        }
    }
    if (fon_solver_level_ >= 3) {
        reks_multi_start_ = true;
    }
    if (fon_solver_level_ >= 4) {
        // Adaptive solver: uses multi-start + coupling-based direction
        // Each point finds the global minimum independently
        reks_multi_start_ = true;
    }

    // Get REKS active space type from options (default: 2,2)
    // REKS_PAIRS = 1 -> REKS(2,2), REKS_PAIRS = 2 -> REKS(4,4)
    int n_pairs = options_.get_int("REKS_PAIRS");
    if (n_pairs < 1) n_pairs = 1;  // Default to REKS(2,2)

    // Compute number of core (doubly occupied) orbitals
    // For REKS(N,M): Ncore = nalpha - n_pairs
    Ncore_ = nalpha_ - n_pairs;

    // Create active space based on number of pairs
    double w_PPS = 0.5;  // Default SA weight for PPS
    double w_OSS = 0.5;  // Default SA weight for OSS

    if (n_pairs == 1) {
        // REKS(2,2): 2 electrons in 2 orbitals
        active_space_ = reks::REKSActiveSpace::create_2_2(w_PPS, w_OSS);
    } else if (n_pairs == 2) {
        // REKS(4,4): 4 electrons in 4 orbitals
        // For REKS(4,4) with 3-state averaging: w_PPS = 1/3
        w_PPS = 1.0 / 3.0;
        active_space_ = reks::REKSActiveSpace::create_4_4(w_PPS);
    } else {
        throw PSIEXCEPTION("REKS: Unsupported REKS_PAIRS value: " + std::to_string(n_pairs)
                          + ". Supported values: 1 (REKS(2,2)), 2 (REKS(4,4))");
    }

    // Get active orbital MO indices from active_space_ (generalized)
    active_mo_indices_ = active_space_->get_active_mo_indices(Ncore_);

    // Initialize vectors for microstate energies and weights
    int n_micro = n_microstates();
    E_micro_.assign(n_micro, 0.0);
    C_L_.assign(n_micro, 0.0);

    // Allocate REKS-specific matrices
    allocate_reks_matrices();

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  ============================================================\n");
        outfile->Printf("                REKS(%d,%d) Initialization\n",
                        active_space_->n_electrons(), active_space_->n_orbitals());
        outfile->Printf("  ------------------------------------------------------------\n");
        outfile->Printf("  Active space: %d electrons in %d orbitals\n",
                        active_space_->n_electrons(), active_space_->n_orbitals());
        outfile->Printf("  Number of GVB pairs:   %d\n", active_space_->n_pairs());
        outfile->Printf("  Number of microstates: %d\n", n_micro);
        outfile->Printf("  Number of base densities: %d\n", active_space_->n_base_densities());
        outfile->Printf("  Core orbitals (Ncore): %3d\n", Ncore_);
        outfile->Printf("  Active MO indices:     ");
        for (int idx : active_mo_indices_) {
            outfile->Printf("%d ", idx);
        }
        outfile->Printf("\n");
        outfile->Printf("  ------------------------------------------------------------\n");
        if (n_pairs == 1) {
            outfile->Printf("  Initial FON: n_r = %.6f, n_s = %.6f\n", get_n_r(), get_n_s());
        } else {
            outfile->Printf("  Initial FONs: n_a=%.4f, n_d=%.4f, n_b=%.4f, n_c=%.4f\n",
                            active_space_->pair(0).fon_p, active_space_->pair(0).fon_q,
                            active_space_->pair(1).fon_p, active_space_->pair(1).fon_q);
        }
        outfile->Printf("  SA weight w_PPS = %.4f\n", active_space_->w_PPS());
        outfile->Printf("  Filatov delta parameter: %.4f\n", reks::DELTA);
        outfile->Printf("  ------------------------------------------------------------\n");
        outfile->Printf("  f(0.5) = %.8f  (should be 1.0)\n", reks::f_interp(0.5));
        outfile->Printf("  f(0.25) = %.8f\n", reks::f_interp(0.25));
        outfile->Printf("  f(0.0) = %.8f  (should be 0.0)\n", reks::f_interp(0.0));
        outfile->Printf("  ============================================================\n\n");
    }

    // NOTE: Ca_ is allocated here but EMPTY (zeros).
    // Actual orbitals come from guess() in compute_energy().
    // Debug output for densities will appear in SCF iterations.
}

void REKS::allocate_reks_matrices() {
    int nso = nsopi_[0];  // Assuming C1 symmetry
    int n_micro = n_microstates();
    int n_base = active_space_->n_base_densities();
    int n_active = active_space_->n_orbitals();

    // Allocate base density matrices (generalized for REKS(N,M))
    // For REKS(2,2): 4 densities (patterns 0-3)
    // For REKS(4,4): 16 densities (patterns 0-15)
    base_densities_.resize(n_base);
    for (int p = 0; p < n_base; ++p) {
        std::string name = "D_base_" + std::to_string(p);
        base_densities_[p] = std::make_shared<Matrix>(name, nsopi_, nsopi_);
    }

    // Allocate microstate density and Fock matrices (dynamically sized)
    D_alpha_micro_.resize(n_micro);
    D_beta_micro_.resize(n_micro);
    F_alpha_micro_.resize(n_micro);
    F_beta_micro_.resize(n_micro);

    for (int L = 0; L < n_micro; ++L) {
        std::string suffix = " L=" + std::to_string(L);
        D_alpha_micro_[L] = std::make_shared<Matrix>("D_alpha" + suffix, nsopi_, nsopi_);
        D_beta_micro_[L] = std::make_shared<Matrix>("D_beta" + suffix, nsopi_, nsopi_);
        F_alpha_micro_[L] = std::make_shared<Matrix>("F_alpha" + suffix, nsopi_, nsopi_);
        F_beta_micro_[L] = std::make_shared<Matrix>("F_beta" + suffix, nsopi_, nsopi_);
    }

    // Allocate coupling Fock matrix
    F_reks_ = std::make_shared<Matrix>("F_REKS (coupling)", nsopi_, nsopi_);

    // Pre-allocate work matrices for JK computation (HPC optimization)
    // These are reused every SCF iteration instead of being reallocated
    C_base_.resize(n_base);
    for (int p = 0; p < n_base; ++p) {
        // Number of occupied orbitals for pattern p = Ncore + popcount(p)
        int n_occ = Ncore_ + __builtin_popcount(p);
        n_occ = std::max(1, n_occ);  // At least 1 column to avoid empty matrix
        std::string name = "C_base_" + std::to_string(p);
        C_base_[p] = std::make_shared<Matrix>(name, nso, n_occ);
    }
    J_work_ = std::make_shared<Matrix>("J_work", nsopi_, nsopi_);

    // Pre-allocate work matrices for AO-to-MO transforms in form_F()
    Temp_work_ = std::make_shared<Matrix>("Temp_work", nso, nso);
    F_alpha_MO_.resize(n_micro);
    F_beta_MO_.resize(n_micro);

    for (int L = 0; L < n_micro; ++L) {
        std::string suffix = " L=" + std::to_string(L);
        F_alpha_MO_[L] = std::make_shared<Matrix>("F_alpha_MO" + suffix, nso, nso);
        F_beta_MO_[L] = std::make_shared<Matrix>("F_beta_MO" + suffix, nso, nso);
    }

    // Allocate XC work matrices (for DFT support)
    D_total_micro_.resize(n_micro);
    V_xc_micro_.resize(n_micro);
    E_xc_micro_.assign(n_micro, 0.0);

    for (int L = 0; L < n_micro; ++L) {
        std::string suffix = " L=" + std::to_string(L);
        D_total_micro_[L] = std::make_shared<Matrix>("D_total" + suffix, nsopi_, nsopi_);
        V_xc_micro_[L] = std::make_shared<Matrix>("V_xc" + suffix, nsopi_, nsopi_);
    }
}

// ---------------------------------------------------------------------------
// Density Matrix Construction (Filatov 2024, Section 2)
// ---------------------------------------------------------------------------

void REKS::build_base_densities() {
    // Builds base density matrices for all occupation patterns from current orbitals Ca_
    // For REKS(2,2): patterns 0-3 (D00, D01, D10, D11)
    // For REKS(4,4): patterns 0-15
    //
    // Pattern bitmask: bit i set = active orbital i occupied

    // Debug: Print C matrix after GUESS (iteration 1 only)
    if (iteration_ == 1 && reks_debug_ >= 2) {
        outfile->Printf("\n  === GUESS C Matrix (Psi4) ===\n");
        outfile->Printf("  First 5 orbitals, first 5 basis functions:\n");
        double** Cp = Ca_->pointer(0);
        int nprint = std::min(5, nsopi_[0]);
        int nmo_print = std::min(5, nsopi_[0]);
        for (int mu = 0; mu < nprint; ++mu) {
            outfile->Printf("  mu=%d:", mu);
            for (int i = 0; i < nmo_print; ++i) {
                outfile->Printf(" %12.8f", Cp[mu][i]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("  ======================================\n");
    }

    int nso = nsopi_[0];  // Assuming C1 symmetry
    int nmo = Ca_->colspi()[0];
    double** Cp = Ca_->pointer(0);
    int n_active = active_space_->n_orbitals();
    int n_base = active_space_->n_base_densities();

    // Build base density for pattern 0 (core only, no active orbitals occupied)
    base_densities_[0]->zero();
    if (Ncore_ > 0) {
        double** Dp = base_densities_[0]->pointer(0);
        C_DGEMM('N', 'T', nso, nso, Ncore_, 1.0, Cp[0], nmo, Cp[0], nmo, 0.0, Dp[0], nso);
    }

    // Build remaining base densities by adding active orbitals based on pattern bits
    for (int p = 1; p < n_base; ++p) {
        // Start with core density
        base_densities_[p]->copy(base_densities_[0]);
        double** Dp = base_densities_[p]->pointer(0);

        // Add each active orbital that is occupied in pattern p
        for (int i = 0; i < n_active; ++i) {
            if (p & (1 << i)) {
                int mo_idx = active_mo_indices_[i];
                C_DGER(nso, nso, 1.0, &Cp[0][mo_idx], nmo, &Cp[0][mo_idx], nmo, Dp[0], nso);
            }
        }
    }

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Base Densities (Generalized) ===\n");
        for (int p = 0; p < n_base; ++p) {
            int n_occ = Ncore_ + __builtin_popcount(p);
            outfile->Printf("  Pattern %d (n_occ=%d): tr(D*S) = %8.4f\n",
                           p, n_occ, base_densities_[p]->vector_dot(S_));
        }
    }
}

void REKS::build_microstate_densities() {
    // Map base densities to microstate alpha/beta densities using active_space_ interface
    // This is now generalized for any REKS(N,M)

    int n_micro = n_microstates();
    for (int L = 0; L < n_micro; ++L) {
        int alpha_idx = active_space_->get_alpha_base_idx(L);
        int beta_idx = active_space_->get_beta_base_idx(L);

        D_alpha_micro_[L]->copy(base_densities_[alpha_idx]);
        D_beta_micro_[L]->copy(base_densities_[beta_idx]);
    }

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Microstate Densities ===\n");
        int n_active = active_space_->n_orbitals();
        for (int L = 0; L < n_micro; ++L) {
            outfile->Printf("  L=%d: base_idx(alpha=%d, beta=%d), tr(D^alpha)=%10.6f, tr(D^beta)=%10.6f\n",
                L, active_space_->get_alpha_base_idx(L), active_space_->get_beta_base_idx(L),
                D_alpha_micro_[L]->trace(), D_beta_micro_[L]->trace());

            // Print active orbital occupations
            outfile->Printf("        occupations: ");
            for (int i = 0; i < n_active; ++i) {
                outfile->Printf("orb%d(a=%d,b=%d) ", i,
                    active_space_->alpha_occ(L, i), active_space_->beta_occ(L, i));
            }
            outfile->Printf("\n");
        }
    }
}

void REKS::compute_weighting_factors() {
    // Compute microstate weights from active space
    active_space_->compute_weights(C_L_);

    if (reks_debug_ >= 2) {
        double n_r = get_n_r();
        double n_s = get_n_s();
        double f = reks::f_interp(n_r * n_s);

        outfile->Printf("\n  === REKS Weighting Factors C_L ===\n");
        outfile->Printf("  n_r = %.6f, n_s = %.6f, f(n_r*n_s) = %.6f\n", n_r, n_s, f);
        outfile->Printf("  w_PPS = %.4f, w_OSS = %.4f\n",
                        active_space_->w_PPS(), active_space_->w_OSS());

        int n_micro = n_microstates();
        double sum = 0.0;
        for (int L = 0; L < n_micro; ++L) {
            outfile->Printf("  C_L[%d] = %12.8f\n", L, C_L_[L]);
            sum += C_L_[L];
        }
        outfile->Printf("  Sum C_L = %12.8f (should be ~1.0)\n", sum);
    }
}

// ---------------------------------------------------------------------------
// Build Microstate Fock Matrices (Filatov 2024, Section 2)
// ---------------------------------------------------------------------------

void REKS::build_microstate_focks() {
    // Build UHF-like Fock matrices for each microstate.
    // Single JK batch for all unique base densities + RHF density for G_
    //
    // For REKS(2,2): 4 base densities (patterns 0-3)
    // For REKS(4,4): 16 base densities (patterns 0-15)

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  build_microstate_focks: starting...\n");
    }

    int nso = nsopi_[0];
    int nmo = Ca_->colspi()[0];
    int n_base = active_space_->n_base_densities();
    int n_active = active_space_->n_orbitals();
    int n_micro = n_microstates();

    // Use pre-allocated work matrices (allocated once in allocate_reks_matrices)
    auto C_occ = Ca_subset("SO", "OCC");  // RHF occupied orbitals for G_

    double** Cp = Ca_->pointer(0);

    // Populate C_base_ matrices for JK computation
    for (int p = 0; p < n_base; ++p) {
        double** Cbp = C_base_[p]->pointer(0);
        int n_occ_in_pattern = __builtin_popcount(p);
        int col = 0;

        // Copy core orbitals first
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_; ++k) {
                Cbp[mu][col + k] = Cp[mu][k];
            }
        }
        col = Ncore_;

        // Add active orbitals based on pattern bits
        for (int i = 0; i < n_active; ++i) {
            if (p & (1 << i)) {
                int mo_idx = active_mo_indices_[i];
                for (int mu = 0; mu < nso; ++mu) {
                    Cbp[mu][col] = Cp[mu][mo_idx];
                }
                col++;
            }
        }

        // Zero remaining columns if any (shouldn't happen if allocation is correct)
        int expected_cols = Ncore_ + n_occ_in_pattern;
        expected_cols = std::max(1, expected_cols);
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = col; k < expected_cols; ++k) {
                Cbp[mu][k] = 0.0;
            }
        }
    }

    // Handle special case: Ncore=0 with pattern 0 (empty density)
    if (Ncore_ == 0) {
        C_base_[0]->zero();
    }

    // Batch JK computation: n_base + 1 densities in one call
    // [0..n_base-1]: base densities
    // [n_base]: RHF occupied for G_

    // Support for multi_scf(): use precomputed J/K if available
    std::vector<SharedMatrix> J_vec, K_vec, wK_vec;

    if (use_precomputed_jk_base_ && !precomputed_J_base_.empty()) {
        // Use precomputed J/K from multi_scf shared JK
        if (reks_debug_ >= 1) {
            outfile->Printf("  build_microstate_focks: using precomputed J/K from multi_scf\n");
        }
        J_vec = precomputed_J_base_;
        K_vec = precomputed_K_base_;

        // Use precomputed wK for LRC functionals
        if (functional_->is_x_lrc() && !precomputed_wK_base_.empty()) {
            wK_vec = precomputed_wK_base_;
        }

        // Reset flag for next iteration
        use_precomputed_jk_base_ = false;
    } else {
        // Normal path: compute J/K using JK builder
        std::vector<SharedMatrix>& C_left = jk_->C_left();
        std::vector<SharedMatrix>& C_right = jk_->C_right();
        C_left.clear();
        C_right.clear();

        for (int p = 0; p < n_base; ++p) {
            C_left.push_back(C_base_[p]);
        }
        C_left.push_back(C_occ);  // Last one for RHF G_

        jk_->compute();

        // Copy references to results
        const std::vector<SharedMatrix>& J_ref = jk_->J();
        const std::vector<SharedMatrix>& K_ref = jk_->K();
        J_vec.assign(J_ref.begin(), J_ref.end());
        K_vec.assign(K_ref.begin(), K_ref.end());

        // Get wK for LRC (long-range corrected) functionals
        if (functional_->is_x_lrc()) {
            const std::vector<SharedMatrix>& wK_ref = jk_->wK();
            wK_vec.assign(wK_ref.begin(), wK_ref.end());
        }
    }

    // Handle special case: pattern 0 with Ncore=0 may have garbage
    if (Ncore_ == 0) {
        J_vec[0]->zero();
        K_vec[0]->zero();
        if (!wK_vec.empty()) wK_vec[0]->zero();
    }

    // HF exchange fraction: 0.0 for pure DFT, 0.5 for BHHLYP, 1.0 for HF
    double alpha = functional_->x_alpha();

    // Build RHF J_, K_, G_ for DIIS and compute_E() compatibility
    // RHF convention: G_ = 2*J - alpha*K + V_xc (if DFT)
    J_ = J_vec[n_base];
    K_ = K_vec[n_base];

    // Compute V_xc for RHF density (for G_ and DIIS compatibility)
    if (functional_->needs_xc() && potential_) {
        potential_->set_D({Da_});
        potential_->compute_V({Va_});
        G_->copy(Va_);
    } else {
        G_->zero();
    }
    G_->axpy(2.0, J_);

    // Exchange contribution - match RHF logic exactly (rhf.cc:243-257)
    // When LRC + wcombine, JK combines alpha*K + beta*wK into wK, don't add K separately
    if (functional_->is_x_hybrid() && !(functional_->is_x_lrc() && jk_->get_wcombine())) {
        G_->axpy(-alpha, K_);
    }

    // Long-range exchange contribution for range-separated functionals
    if (functional_->is_x_lrc() && !wK_vec.empty()) {
        double beta = functional_->x_beta();
        if (jk_->get_wcombine()) {
            // wK contains combined alpha*K + beta*wK
            G_->axpy(-1.0, wK_vec[n_base]);
        } else {
            // Add LR exchange separately
            G_->axpy(-beta, wK_vec[n_base]);
        }
    }

    // Build Fock matrices for each microstate
    // F^sigma_L = H + J_total_L - alpha*K^sigma_L + V_xc_L
    // where J_total = J^alpha + J^beta, K^sigma uses the sigma density

    // XC support: compute V_xc for each microstate if functional needs XC
    // Each microstate gets its own XC potential based on its total density
    bool needs_xc = functional_->needs_xc();

    if (needs_xc && potential_) {
        if (reks_debug_ >= 2) {
            outfile->Printf("  Computing XC for %d microstates...\n", n_micro);
        }

        // First, compute V_xc for the RHF density Da_ (for G_ and finalize compatibility)
        // This ensures RHO_A etc are populated for finalize_energy()
        potential_->set_D({Da_});
        potential_->compute_V({Va_});

        // Compute XC for each microstate
        for (int L = 0; L < n_micro; ++L) {
            // Build total density for microstate L: D_total = D_alpha + D_beta
            D_total_micro_[L]->copy(D_alpha_micro_[L]);
            D_total_micro_[L]->add(D_beta_micro_[L]);

            // Zero V_xc before computing
            V_xc_micro_[L]->zero();

            // IMPORTANT: RKS grid expects "alpha density" and multiplies by 2 internally
            // to get total density. So we pass D_total/2 to get the correct XC.
            auto D_half = std::make_shared<Matrix>("D_half", D_total_micro_[L]->rowspi(), D_total_micro_[L]->colspi());
            D_half->copy(D_total_micro_[L]);
            D_half->scale(0.5);

            // Set density for XC computation (RKS-style expects half-density)
            potential_->set_D({D_half});

            // Compute V_xc
            potential_->compute_V({V_xc_micro_[L]});

            // Get XC energy from quadrature
            E_xc_micro_[L] = potential_->quadrature_values()["FUNCTIONAL"];

            if (reks_debug_ >= 2) {
                outfile->Printf("    L=%d: tr(D)=%.4f, E_xc=%.10f\n",
                               L, D_total_micro_[L]->trace(), E_xc_micro_[L]);
            }
        }
    }

    for (int L = 0; L < n_micro; ++L) {
        int alpha_idx = active_space_->get_alpha_base_idx(L);
        int beta_idx = active_space_->get_beta_base_idx(L);

        // J_total = J_alpha + J_beta
        J_work_->copy(J_vec[alpha_idx]);
        J_work_->add(J_vec[beta_idx]);

        // F^sigma = H + J_total
        F_alpha_micro_[L]->copy(H_);
        F_alpha_micro_[L]->add(J_work_);
        F_beta_micro_[L]->copy(H_);
        F_beta_micro_[L]->add(J_work_);

        // Exchange contribution - must match RHF logic (rhf.cc:243-257)
        // When LRC + wcombine, K is already combined into wK, don't add K separately
        if (functional_->is_x_hybrid()) {
            if (!(functional_->is_x_lrc() && jk_->get_wcombine())) {
                F_alpha_micro_[L]->axpy(-alpha, K_vec[alpha_idx]);
                F_beta_micro_[L]->axpy(-alpha, K_vec[beta_idx]);
            }
        }

        // Long-range exchange contribution for range-separated functionals
        if (functional_->is_x_lrc() && !wK_vec.empty()) {
            double beta = functional_->x_beta();
            if (jk_->get_wcombine()) {
                // wK contains combined alpha*K + beta*wK
                F_alpha_micro_[L]->axpy(-1.0, wK_vec[alpha_idx]);
                F_beta_micro_[L]->axpy(-1.0, wK_vec[beta_idx]);
            } else {
                // Add LR exchange separately
                F_alpha_micro_[L]->axpy(-beta, wK_vec[alpha_idx]);
                F_beta_micro_[L]->axpy(-beta, wK_vec[beta_idx]);
            }
        }

        // Add XC potential to both alpha and beta Fock matrices
        if (needs_xc && potential_) {
            F_alpha_micro_[L]->add(V_xc_micro_[L]);
            F_beta_micro_[L]->add(V_xc_micro_[L]);
        }
    }

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Microstate Fock Matrices (Generalized) ===\n");
        outfile->Printf("  Exchange scaling alpha = %.4f\n", alpha);
        outfile->Printf("  Number of base densities in JK batch: %d\n", n_base);

        // Compare with RHF Fock for validation
        Fa_->copy(H_);
        Fa_->add(G_);
        double tr_rhf = Fa_->trace();
        outfile->Printf("  RHF Fock: tr(Fa_RHF) = %15.10f\n", tr_rhf);

        for (int L = 0; L < n_micro; ++L) {
            double tr_alpha = F_alpha_micro_[L]->trace();
            double tr_beta = F_beta_micro_[L]->trace();

            int alpha_idx = active_space_->get_alpha_base_idx(L);
            int beta_idx = active_space_->get_beta_base_idx(L);

            outfile->Printf("  L=%d: base_idx(alpha=%d, beta=%d), tr(F^alpha)=%15.10f, tr(F^beta)=%15.10f\n",
                           L, alpha_idx, beta_idx, tr_alpha, tr_beta);
        }
    }
}

// ---------------------------------------------------------------------------
// Compute Microstate Energies (Filatov 2024, Eq. 5)
// ---------------------------------------------------------------------------

void REKS::compute_microstate_energies() {
    // E_L = tr(D^alpha_L * H) + tr(D^beta_L * H)
    //     + 0.5 * [tr(D^alpha_L * G^alpha_L) + tr(D^beta_L * G^beta_L)]
    //     + E_nuc + E_xc_correction
    // where G^sigma = F^sigma - H = J - alpha*K^sigma + V_xc
    //
    // XC correction: The Fock matrix energy includes 0.5*tr(D*V_xc), but
    // the actual XC energy is E_xc from quadrature. So we need to correct:
    // E_xc_correction = E_xc - 0.5*tr(D_total*V_xc)

    // Use nuclearrep_ directly - energies_["Nuclear"] is not set in REKS
    // (unlike RHF::compute_E() which sets it)
    double E_nuc = nuclearrep_;
    bool needs_xc = functional_->needs_xc();
    double alpha = functional_->x_alpha();

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  Microstate energy computation:\n");
        outfile->Printf("  E_nuc = %.10f, needs_xc = %d, alpha = %.4f\n", E_nuc, needs_xc, alpha);
    }

    int n_micro = n_microstates();
    for (int L = 0; L < n_micro; ++L) {
        SharedMatrix Da = D_alpha_micro_[L];
        SharedMatrix Db = D_beta_micro_[L];
        SharedMatrix Fa = F_alpha_micro_[L];
        SharedMatrix Fb = F_beta_micro_[L];

        // One-electron energy: tr(D^alpha * H) + tr(D^beta * H)
        double E_1e = Da->vector_dot(H_) + Db->vector_dot(H_);

        // Two-electron energy from Fock matrices:
        // E_2e = 0.5 * [tr(D^alpha * (F^alpha - H)) + tr(D^beta * (F^beta - H))]
        //      = 0.5 * [tr(D^alpha * F^alpha) + tr(D^beta * F^beta) - E_1e]
        double E_Fa = Da->vector_dot(Fa);
        double E_Fb = Db->vector_dot(Fb);
        double E_2e = 0.5 * (E_Fa + E_Fb - E_1e);

        E_micro_[L] = E_1e + E_2e + E_nuc;

        // XC energy correction for DFT
        // The Fock matrix trick gives: 0.5*tr(D*V_xc) instead of E_xc
        // Correction = E_xc - 0.5*tr(D_total*V_xc)
        if (needs_xc && potential_) {
            double tr_D_Vxc = D_total_micro_[L]->vector_dot(V_xc_micro_[L]);
            double xc_correction = E_xc_micro_[L] - 0.5 * tr_D_Vxc;

            if (reks_debug_ >= 2) {
                outfile->Printf("  L=%d: E_1e=%.10f, E_2e=%.10f, E_nuc=%.10f\n", L, E_1e, E_2e, E_nuc);
                outfile->Printf("        E_micro(before XC)=%.10f\n", E_micro_[L]);
                outfile->Printf("        E_xc=%.10f, tr(D*Vxc)=%.10f, 0.5*tr=%.10f\n",
                               E_xc_micro_[L], tr_D_Vxc, 0.5*tr_D_Vxc);
                outfile->Printf("        xc_correction=%.10f\n", xc_correction);
            }

            E_micro_[L] += xc_correction;

            if (reks_debug_ >= 2) {
                outfile->Printf("        E_micro(after XC)=%.10f\n", E_micro_[L]);
            }
        }
    }

    if (reks_debug_ >= 2) {
        outfile->Printf("  Microstate energies: E[0]=%.10f E[1]=%.10f E[2]=%.10f E[3]=%.10f\n",
                        E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
    }

    if (reks_debug_ >= 1) {
        print_microstate_energies();
    }
}

// ---------------------------------------------------------------------------
// RexSolver - Newton-Raphson FON Optimization
// ---------------------------------------------------------------------------

void REKS::rex_solver() {
    // Optimize FON variables to minimize E_SA = sum_L FACT_L * C_L(FON) * E_L
    //
    // For REKS(2,2): 1 variable (n_r), constraint n_s = 2 - n_r
    // For REKS(4,4): 2 variables (n_a, n_b), constraints n_d = 2 - n_a, n_c = 2 - n_b
    //
    // Uses Newton-Raphson optimization with gradient and Hessian from active_space_.

    // Skip if E_micro_ not yet computed (first SCF iteration with SAD guess)
    if (std::abs(E_micro_[0]) < reks::constants::ENERGY_THRESHOLD &&
        std::abs(E_micro_[1]) < reks::constants::ENERGY_THRESHOLD) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  RexSolver: skipping (E_micro not yet computed)\n");
        }
        return;
    }

    int n_pairs = active_space_->n_pairs();

    if (n_pairs == 1) {
        // REKS(2,2): single FON variable
        rex_solver_22();
    } else if (n_pairs == 2) {
        // REKS(4,4): two FON variables
        // FON_SOLVER levels: 1=Fast(NR), 2=Balanced(TR), 3=Robust(TR+MS)
        if (reks_debug_ >= 1 && iteration_ == 1) {
            outfile->Printf("\n  ================================================================================\n");
            outfile->Printf("  REKS(4,4) FON Solver Configuration\n");
            outfile->Printf("  ================================================================================\n");
            outfile->Printf("  FON_SOLVER level: %d", fon_solver_level_);
            if (fon_solver_level_ == 1) outfile->Printf(" (Fast: Newton-Raphson)\n");
            else if (fon_solver_level_ == 2) outfile->Printf(" (Balanced: Trust-Region)\n");
            else if (fon_solver_level_ == 3) outfile->Printf(" (Robust: Trust-Region + Multi-Start)\n");
            else outfile->Printf(" (Adaptive: Continuation + Early-Transition)\n");
            if (reks_enforce_symmetry_) {
                outfile->Printf("  Symmetry enforcement: ENABLED (n_a = n_b)\n");
            }
            outfile->Printf("  ================================================================================\n");
        }

        if (fon_solver_level_ >= 4) {
            // FON_SOLVER=4: Adaptive solver
            // Continuation-first + early transition detection + symmetry enforcement
            adaptive_fon_solver_44();
        } else if (fon_solver_level_ >= 3 && reks_multi_start_) {
            // FON_SOLVER=3: Adaptive Multi-Start
            // Use Trust-Region by default, trigger multi-start only when needed
            // (first iteration, indefinite Hessian, or large FON change)
            trust_region_fon_optimizer_44_adaptive();
        } else if (use_trust_region_44_) {
            // FON_SOLVER=2: Trust-Region Augmented Hessian method
            trust_region_fon_optimizer_44();
        } else {
            // FON_SOLVER=1: Original Newton-Raphson with coupling-based direction
            rex_solver_44();
        }
    } else {
        throw PSIEXCEPTION("RexSolver: unsupported number of GVB pairs: " + std::to_string(n_pairs));
    }
}

void REKS::rex_solver_22() {
    // REKS(2,2) FON optimization: single variable n_r
    //
    // TWO-PHASE SCF FOR STABLE CONVERGENCE:
    //
    // Phase 1 (iter < fon_freeze_iters_):
    //   - Freeze FON at initial value (e.g., n_r = 1.0)
    //   - Let orbitals converge to this fixed electronic state
    //   - DIIS works correctly with consistent C_L weights
    //
    // Phase 2 (iter >= fon_freeze_iters_):
    //   - Enable FON optimization
    //   - Limit FON change per SCF iteration (fon_max_delta_)
    //   - Ensures Fock matrix continuity for DIIS
    //
    // This approach is used in CASSCF and other multi-configurational methods.

    double n_r_current = get_n_r();

    // =========================================================================
    // Phase 1: Freeze FON for initial iterations
    // =========================================================================
    if (iteration_ < fon_freeze_iters_) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === RexSolver(2,2): Phase 1 (FON frozen) ===\n");
            outfile->Printf("  Iteration %d < %d: keeping n_r = %.6f (frozen)\n",
                           iteration_, fon_freeze_iters_, n_r_current);
        }
        // Update previous FON for Phase 2 tracking
        prev_n_r_ = n_r_current;
        compute_weighting_factors();
        return;
    }

    // Mark phase transition (for potential DIIS reset)
    if (iteration_ == fon_freeze_iters_) {
        fon_phase_transition_ = true;
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === Phase transition: enabling FON optimization ===\n");
        }
    }

    // =========================================================================
    // Phase 2: Optimize FON with step limiting
    // =========================================================================

    const int max_iter = reks::constants::FON_MAX_ITER;
    const double tol = reks::constants::FON_TOL;
    const double max_step = reks::constants::FON_MAX_STEP;

    double w_PPS = active_space_->w_PPS();
    double w_OSS = active_space_->w_OSS();

    // Lambda to compute E_SA at given n_r
    auto compute_E_SA = [&](double nr) -> double {
        double ns = 2.0 - nr;
        double x = nr * ns;
        double f = reks::f_interp(x);
        double C0 = w_PPS * nr / 2.0;
        double C1 = w_PPS * ns / 2.0;
        double C2 = w_OSS - 0.5 * w_PPS * f;
        double C3 = 0.5 * w_PPS * f - 0.5 * w_OSS;
        return C0 * E_micro_[0] + C1 * E_micro_[1]
             + 2.0 * C2 * E_micro_[2] + 2.0 * C3 * E_micro_[3];
    };

    // Multi-start optimization to find the global minimum
    // Try multiple starting points to avoid getting stuck in local minima.
    // This is essential for finding correct multi-configurational solutions.
    std::vector<double> start_points = {n_r_current, 0.5, 1.0, 1.5};

    double best_n_r = n_r_current;
    double best_E_SA = compute_E_SA(n_r_current);

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver(2,2): Phase 2 (FON optimization) ===\n");
        outfile->Printf("  E_micro: [0]=%.8f, [1]=%.8f, [2]=%.8f, [3]=%.8f\n",
                       E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
        outfile->Printf("  prev_n_r = %.6f, fon_max_delta = %.4f\n",
                       prev_n_r_, fon_max_delta_);
    }

    for (double n_r_start : start_points) {
        double n_r = n_r_start;
        double n_s = 2.0 - n_r;

        if (reks_debug_ >= 1) {
            outfile->Printf("  Starting from: n_r=%10.6f, n_s=%10.6f\n", n_r, n_s);
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            n_s = 2.0 - n_r;

            // Interpolating function and derivatives at x = n_r*n_s
            double x = n_r * n_s;
            double f = reks::f_interp(x);
            double df_dx = reks::df_interp(x);
            double d2f_dx2 = reks::d2f_interp(x);

            // Chain rule: df/dn_r = (df/dx) * (dx/dn_r) = (df/dx) * (n_s - n_r)
            double dx_dnr = n_s - n_r;
            double df_dnr = df_dx * dx_dnr;
            double d2f_dnr2 = d2f_dx2 * dx_dnr * dx_dnr - 2.0 * df_dx;

            // Compute C_L at current n_r
            double C0 = w_PPS * n_r / 2.0;
            double C1 = w_PPS * n_s / 2.0;
            double C2 = w_OSS - 0.5 * w_PPS * f;
            double C3 = 0.5 * w_PPS * f - 0.5 * w_OSS;

            double E_SA_current = C0 * E_micro_[0] + C1 * E_micro_[1]
                                + 2.0 * C2 * E_micro_[2] + 2.0 * C3 * E_micro_[3];

            // Energy gradient and Hessian
            double dE_dnr = (w_PPS / 2.0) * (E_micro_[0] - E_micro_[1])
                          + w_PPS * df_dnr * (E_micro_[3] - E_micro_[2]);
            double d2E_dnr2 = w_PPS * d2f_dnr2 * (E_micro_[3] - E_micro_[2]);

            // Newton-Raphson step
            double delta;
            if (std::abs(d2E_dnr2) < reks::constants::HESSIAN_THRESHOLD) {
                delta = -std::copysign(max_step, dE_dnr);
            } else {
                delta = -dE_dnr / d2E_dnr2;
            }

            // Damping
            if (std::abs(delta) > max_step) {
                delta = std::copysign(max_step, delta);
            }

            // Update n_r with bounds [0.0, 2.0]
            double n_r_new = std::clamp(n_r + delta, 0.0, 2.0);
            delta = n_r_new - n_r;
            n_r = n_r_new;

            if (reks_debug_ >= 2) {
                outfile->Printf("    iter %2d: n_r=%.6f E=%.10f dE=%.2e d2E=%.2e delta=%.6f\n",
                               iter, n_r, E_SA_current, dE_dnr, d2E_dnr2, delta);
            }

            // Convergence check
            if (std::abs(delta) < tol) {
                break;
            }
        }

        // Check if this starting point gave better result
        double E_SA_final = compute_E_SA(n_r);
        if (E_SA_final < best_E_SA) {
            best_E_SA = E_SA_final;
            best_n_r = n_r;
        }

        if (reks_debug_ >= 1) {
            outfile->Printf("  Result: n_r=%.8f, E_SA=%.10f\n", n_r, E_SA_final);
        }
    }

    // =========================================================================
    // Apply FON change limiting for DIIS stability
    // =========================================================================
    // Limit the change from previous SCF iteration to fon_max_delta_
    // This ensures Fock matrices remain similar between DIIS iterations

    double delta_from_prev = best_n_r - prev_n_r_;

    if (std::abs(delta_from_prev) > fon_max_delta_) {
        double limited_n_r = prev_n_r_ + std::copysign(fon_max_delta_, delta_from_prev);
        limited_n_r = std::clamp(limited_n_r, 0.0, 2.0);

        if (reks_debug_ >= 1) {
            outfile->Printf("  FON change limited: %.6f -> %.6f (delta %.4f -> %.4f)\n",
                           best_n_r, limited_n_r, delta_from_prev, limited_n_r - prev_n_r_);
        }

        best_n_r = limited_n_r;
        best_E_SA = compute_E_SA(best_n_r);
    }

    // Update FON and weights
    active_space_->set_pair_fon(0, best_n_r);
    compute_weighting_factors();

    // Store current n_r for next iteration's delta limiting
    prev_n_r_ = best_n_r;

    if (reks_debug_ >= 1) {
        outfile->Printf("  Final: n_r = %12.8f, n_s = %12.8f, E_SA = %20.12f\n",
                        get_n_r(), get_n_s(), best_E_SA);
    }
}

void REKS::initialize_fon_from_hamiltonian_44() {
    // =========================================================================
    // 4×4 CI Hamiltonian Diagonalization for Initial FON Guess
    // =========================================================================
    //
    // Based on c-code rexslv4x4() lines 340-361 (Filatov's GAMESS implementation)
    //
    // Instead of starting at the saddle point n_a = n_b = 1.0, we diagonalize
    // the 4×4 CI Hamiltonian to get physically motivated initial FONs.
    //
    // Hamiltonian basis states:
    //   |1⟩ = |aābb̄⟩  (E_micro[0] = eab)
    //   |2⟩ = |aācc̄⟩  (E_micro[1] = eac)
    //   |3⟩ = |bb̄dd̄⟩  (E_micro[2] = ebd)
    //   |4⟩ = |cc̄dd̄⟩  (E_micro[3] = ecd)
    //
    // Off-diagonal coupling elements:
    //   dad: coupling between (a,d) pair (derived from E_micro[4-7])
    //   dbc: coupling between (b,c) pair (derived from E_micro[8-11])

    if (E_micro_.size() < 12) {
        outfile->Printf("  [HAM-DIAG] Not enough microstates (need 12), skipping\n");
        return;
    }

    // Extract closed-shell microstate energies (diagonal elements)
    double eab = E_micro_[0];  // |aābb̄|
    double eac = E_micro_[1];  // |aācc̄|
    double ebd = E_micro_[2];  // |bb̄dd̄|
    double ecd = E_micro_[3];  // |cc̄dd̄|

    // Compute coupling elements from singlet-triplet splitting
    // Formula derived from weight analysis:
    //   PPS weights: C_4=C_5=-f_ad/2 (singlet), C_6=C_7=+f_ad/2 (triplet)
    //   Coupling contribution: -f_ad * dad = -f_ad/2*(E[4]+E[5]) + f_ad/2*(E[6]+E[7])
    //   Therefore: dad = (E[4]+E[5]-E[6]-E[7])/2
    double dad = (E_micro_[4] + E_micro_[5] - E_micro_[6] - E_micro_[7]) / 2.0;
    double dbc = (E_micro_[8] + E_micro_[9] - E_micro_[10] - E_micro_[11]) / 2.0;

    // Coupling threshold check (from c-code lines 365-366)
    const double coupling_threshold = 1.0e-8;
    bool fix_na = (std::fabs(dad) <= coupling_threshold);
    bool fix_nb = (std::fabs(dbc) <= coupling_threshold);

    if (fix_na && fix_nb) {
        // Both couplings negligible → closed-shell solution
        active_space_->set_pair_fon(0, 2.0);  // n_a = 2.0
        active_space_->set_pair_fon(1, 2.0);  // n_b = 2.0
        outfile->Printf("  [HAM-DIAG] Weak coupling (dad=%.2e, dbc=%.2e): n_a=n_b=2.0\n", dad, dbc);
        compute_weighting_factors();
        return;
    }

    // Build 4×4 Hamiltonian (column-major for LAPACK)
    // Matrix structure (from c-code):
    //  | eab | dbc | dad |  0  |
    //  | dbc | eac |  0  | dad |
    //  | dad |  0  | ebd | dbc |
    //  |  0  | dad | dbc | ecd |
    std::vector<double> H(16, 0.0);
    // Column-major indexing: H[i + j*4] = H[i][j]
    H[0]  = eab;   // H[0,0]
    H[5]  = eac;   // H[1,1]
    H[10] = ebd;   // H[2,2]
    H[15] = ecd;   // H[3,3]
    H[1]  = dbc;   // H[1,0] = H[0,1]
    H[4]  = dbc;   // H[0,1]
    H[2]  = dad;   // H[2,0] = H[0,2]
    H[8]  = dad;   // H[0,2]
    H[7]  = dad;   // H[3,1] = H[1,3]
    H[13] = dad;   // H[1,3]
    H[11] = dbc;   // H[3,2] = H[2,3]
    H[14] = dbc;   // H[2,3]

    // Eigenvalue decomposition using LAPACK dsyev
    int n = 4;
    std::vector<double> w(4);      // eigenvalues
    std::vector<double> work(16);
    int lwork = 16;
    int info;

    // C_DSYEV: computes eigenvalues and eigenvectors of symmetric matrix
    // 'V' = compute eigenvectors, 'L' = use lower triangle
    // Returns info: 0 = success, >0 = failed to converge
    info = C_DSYEV('V', 'L', n, H.data(), n, w.data(), work.data(), lwork);

    if (info != 0) {
        outfile->Printf("  [HAM-DIAG] dsyev failed, info=%d, using n_a=n_b=1.0\n", info);
        active_space_->set_pair_fon(0, 1.0);
        active_space_->set_pair_fon(1, 1.0);
        compute_weighting_factors();
        return;
    }

    // Extract FONs from ground state eigenvector (column 0 after dsyev)
    // Column-major: v[i] = H[i + 0*4] = H[i]
    double v0 = H[0], v1 = H[1], v2 = H[2], v3 = H[3];

    // From c-code (lines 359-360):
    // x1 = v0² + v1²  (FON for pair a,d in x∈[0,1] coordinates)
    // x2 = v0² + v2²  (FON for pair b,c in x∈[0,1] coordinates)
    // Convert to n∈[0,2]: n = 2*x
    double x1 = v0*v0 + v1*v1;
    double x2 = v0*v0 + v2*v2;
    double n_a = 2.0 * x1;
    double n_b = 2.0 * x2;

    // Apply coupling threshold fixes
    if (fix_na) n_a = 2.0;
    if (fix_nb) n_b = 2.0;

    // Boundary protection (from c-code line 215)
    const double boundary_eps = 1.0e-9;
    n_a = std::clamp(n_a, boundary_eps, 2.0 - boundary_eps);
    n_b = std::clamp(n_b, boundary_eps, 2.0 - boundary_eps);

    // Set FONs and update weights
    active_space_->set_pair_fon(0, n_a);
    active_space_->set_pair_fon(1, n_b);
    compute_weighting_factors();

    // Update previous FON values for delta limiting
    prev_n_a_ = n_a;
    prev_n_b_ = n_b;

    outfile->Printf("\n  [HAM-DIAG] Initial FON from 4×4 CI Hamiltonian diagonalization:\n");
    outfile->Printf("    E_micro[0-3] (closed-shell): %.6f %.6f %.6f %.6f\n",
                   E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
    outfile->Printf("    Coupling elements: dad=%.6f dbc=%.6f\n", dad, dbc);
    outfile->Printf("    Eigenvalues: %.6f %.6f %.6f %.6f\n", w[0], w[1], w[2], w[3]);
    outfile->Printf("    Eigenvector[0]: [%.4f, %.4f, %.4f, %.4f]\n", v0, v1, v2, v3);
    outfile->Printf("    Initial FON: n_a=%.6f n_b=%.6f\n", n_a, n_b);
}

void REKS::rex_solver_44() {
    // REKS(4,4) FON optimization: two variables (n_a, n_b)
    //
    // Uses 4x4 CI Hamiltonian diagonalization for initial guess (every call)
    // followed by Newton-Raphson refinement (like c-code rexslv4x4).
    //
    // Uses 2D Newton-Raphson with gradient [dE/dn_a, dE/dn_b]
    // and 2x2 Hessian [d²E/dn_a², d²E/dn_a dn_b; d²E/dn_b dn_a, d²E/dn_b²]
    //
    // Constraints: n_d = 2 - n_a, n_c = 2 - n_b

    // Get current FONs
    double n_a = active_space_->pair(0).fon_p;
    double n_b = active_space_->pair(1).fon_p;

    // =========================================================================
    // Hamiltonian Diagonalization for Initial FON Guess (ONCE after localization)
    // =========================================================================
    // HAM-DIAG provides physically motivated initial guess from 4x4 CI Hamiltonian.
    // Done ONCE after orbital localization (iteration >= 2).
    //
    // NOTE: C-code does HAM-DIAG every call, but in PSI4 this destabilizes SCF
    // because orbitals change significantly between iterations for problematic
    // geometries. So we do HAM-DIAG once, then use previous FON as starting point.

    if (!fon_initialized_44_) {
        // Before localization (iteration_ < 2): keep FONs frozen at n=1.0
        if (iteration_ < 2) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  [HAM-DIAG] Waiting for localization (iteration %d < 2), keeping n_a=n_b=1.0\n", iteration_);
            }
            active_space_->set_pair_fon(0, 1.0);
            active_space_->set_pair_fon(1, 1.0);
            compute_weighting_factors();
            return;
        }

        // After localization (iteration_ >= 2): initialize FON
        // Check if user provided initial FON values (for geometry scan continuation)
        if (reks_initial_na_ >= 0.0 && reks_initial_nb_ >= 0.0) {
            // User-provided initial FON values for smooth PES scans
            active_space_->set_pair_fon(0, reks_initial_na_);
            active_space_->set_pair_fon(1, reks_initial_nb_);
            n_a = reks_initial_na_;
            n_b = reks_initial_nb_;
            if (reks_debug_ >= 1) {
                outfile->Printf("  [FON-INIT] Using user-provided initial FON: n_a=%.6f, n_b=%.6f\n",
                               reks_initial_na_, reks_initial_nb_);
            }
        } else {
            // Default: HAM-DIAG initialization
            initialize_fon_from_hamiltonian_44();
            n_a = active_space_->pair(0).fon_p;
            n_b = active_space_->pair(1).fon_p;
        }
        fon_initialized_44_ = true;
    }

    // =========================================================================
    // "Shake it up" - Prevent FON locking at degenerate point n=1.0
    // =========================================================================
    // At n_a = n_d = 1.0 (or n_b = n_c = 1.0), the gradient vanishes and
    // the optimizer gets stuck. This is from Filatov's original GAMESS code.
    if (reks_shake_threshold_ > 0.0) {
        const double shake_target = 0.5;  // Reset FON to this value
        bool shaken = false;

        if (std::abs(n_a - 1.0) < reks_shake_threshold_) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  Shake it up: n_a was %.8f (near 1.0), resetting to %.4f\n",
                               n_a, shake_target);
            }
            n_a = shake_target;
            active_space_->set_pair_fon(0, n_a);
            shaken = true;
        }
        if (std::abs(n_b - 1.0) < reks_shake_threshold_) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  Shake it up: n_b was %.8f (near 1.0), resetting to %.4f\n",
                               n_b, shake_target);
            }
            n_b = shake_target;
            active_space_->set_pair_fon(1, n_b);
            shaken = true;
        }
        if (shaken) {
            compute_weighting_factors();
        }
    }

    // =========================================================================
    // Coupling check (optional) - Fix FON if coupling is negligible
    // =========================================================================
    // If coupling element Δad or Δbc is close to zero, fix corresponding FON.
    // This is from Filatov's original code: if(fabs(dad)<=1.e-8) *x1 = one;
    if (reks_coupling_threshold_ > 0.0 && static_cast<int>(E_micro_.size()) >= 12) {
        // Approximate coupling from microstate energies
        // dad ~ |E_coupling_ad - E_closed_shell|
        // dbc ~ |E_coupling_bc - E_closed_shell|
        double E_closed = (E_micro_[0] + E_micro_[1] + E_micro_[2] + E_micro_[3]) / 4.0;
        double E_coup_ad = (E_micro_[4] + E_micro_[5] + E_micro_[6] + E_micro_[7]) / 4.0;
        double E_coup_bc = (E_micro_[8] + E_micro_[9] + E_micro_[10] + E_micro_[11]) / 4.0;

        double dad = std::abs(E_coup_ad - E_closed);
        double dbc = std::abs(E_coup_bc - E_closed);

        if (dad < reks_coupling_threshold_) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  Coupling dad=%.2e < threshold: fixing n_a=2.0\n", dad);
            }
            n_a = 2.0;
            active_space_->set_pair_fon(0, n_a);
            compute_weighting_factors();
        }
        if (dbc < reks_coupling_threshold_) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  Coupling dbc=%.2e < threshold: fixing n_b=2.0\n", dbc);
            }
            n_b = 2.0;
            active_space_->set_pair_fon(1, n_b);
            compute_weighting_factors();
        }
    }

    // =========================================================================
    // FON Optimization (Newton-Raphson)
    // =========================================================================
    //
    // Uses Newton-Raphson for 2D optimization of (n_a, n_b) with:
    // - Step size limiting (max_step = 0.2)
    // - Boundary handling (0 <= n <= 2)
    // - Oscillation detection for gradient descent fallback

    // Skip optimization if FONs are frozen due to oscillation detection
    if (fon_oscillation_frozen_) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === RexSolver(4,4): FON frozen (oscillation), skipping optimization ===\n");
            outfile->Printf("  Frozen FONs: n_a=%.6f, n_b=%.6f\n", n_a, n_b);
        }
        compute_weighting_factors();
        return;
    }

    const int max_iter = reks::constants::FON_MAX_ITER;
    const double tol = reks::constants::FON_TOL;
    const double osc_threshold = reks::constants::OSCILLATION_THRESHOLD;  // 0.02

    // Track if saddle escape was used (to bypass fon_max_delta limiting)
    bool used_saddle_escape_a = false;
    bool used_saddle_escape_b = false;

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver(4,4): FON Optimization ===\n");
        outfile->Printf("  Initial FONs: n_a=%.6f, n_b=%.6f\n", n_a, n_b);
        outfile->Printf("  prev_n_a=%.6f, prev_n_b=%.6f, fon_max_delta=%.4f\n",
                       prev_n_a_, prev_n_b_, fon_max_delta_);
        outfile->Printf("  E_micro: ");
        for (int L = 0; L < 12 && L < static_cast<int>(E_micro_.size()); ++L) {
            outfile->Printf("[%d]=%.8f ", L, E_micro_[L]);
        }
        outfile->Printf("\n");
    }

    // =========================================================================
    // Coupling Elements: Physical Compass for FON Direction
    // =========================================================================
    // From Filatov's GVB-PP Hamiltonian (c-code reks_solver.cpp:340-361):
    //   dad = (E[4] + E[5] - E[6] - E[7]) / 2  -> (a,d) pair coupling
    //   dbc = (E[8] + E[9] - E[10] - E[11]) / 2  -> (b,c) pair coupling
    //
    // Physical meaning:
    //   dad < 0: exchange stabilization (singlet lower than triplet) in (a,d)
    //   |dad| large: strong correlation -> n_a should tend towards 1.0
    //   |dad| small: weak correlation -> n_a can be near 2.0 (closed-shell)
    //
    double coupling_dad = 0.0;
    double coupling_dbc = 0.0;
    if (E_micro_.size() >= 12) {
        coupling_dad = (E_micro_[4] + E_micro_[5] - E_micro_[6] - E_micro_[7]) / 2.0;
        coupling_dbc = (E_micro_[8] + E_micro_[9] - E_micro_[10] - E_micro_[11]) / 2.0;
    }

    // Threshold for "weak coupling" - below this, closed-shell (n=2) may be valid
    const double coupling_threshold = 1.0e-3;  // Hartree, ~0.6 kcal/mol

    // Determine preferred direction for each FON based on coupling strength
    // prefer_diradical = true means prefer n → 1.0 (equal mixture)
    // prefer_diradical = false means n → 2.0 (closed-shell) may be acceptable
    bool prefer_diradical_a = (std::abs(coupling_dad) > coupling_threshold);
    bool prefer_diradical_b = (std::abs(coupling_dbc) > coupling_threshold);

    if (reks_debug_ >= 1) {
        outfile->Printf("  Coupling elements: dad=%.6f dbc=%.6f (threshold=%.1e)\n",
                       coupling_dad, coupling_dbc, coupling_threshold);
        outfile->Printf("  Preferred direction: n_a→%s, n_b→%s\n",
                       prefer_diradical_a ? "1.0 (diradical)" : "2.0 (closed-shell allowed)",
                       prefer_diradical_b ? "1.0 (diradical)" : "2.0 (closed-shell allowed)");
    }

    // Lambda to compute E_SA at current FON values
    auto compute_E_SA = [&]() -> double {
        std::vector<double> C_L;
        active_space_->compute_weights(C_L);
        double E_SA = 0.0;
        int n_micro = active_space_->n_microstates();
        for (int L = 0; L < n_micro; ++L) {
            double fact = active_space_->get_symmetry_factor(L);
            E_SA += fact * C_L[L] * E_micro_[L];
        }
        return E_SA;
    };

    // Lambda to compute E_PPS (PPS diagonal energy) - the TRUE target function for FON optimization
    // Use PSI4's existing method which computes E_PPS from weight-based formula
    auto compute_E_PPS = [&]() -> double {
        return active_space_->compute_energy_PPS(E_micro_);
    };

    // Oscillation detection: track last FON values
    double prev_n_a = n_a, prev_n_b = n_b;
    double prev2_n_a = -1.0, prev2_n_b = -1.0;
    int oscillation_count = 0;

    // Maximum step size for FON optimization
    const double max_step = reks::constants::FON_MAX_STEP;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute gradient [dE_PPS/dn_a, dE_PPS/dn_b] - use PPS gradient for consistency
        // with E_PPS line search (this is what c-code does - optimize GVB-PP energy)
        std::vector<double> grad = active_space_->compute_gradient_PPS(E_micro_);

        // --- Boundary-aware convergence check ---
        // At a boundary, if the gradient points outward, we're at the constrained optimum.
        // Effective gradient considers only interior-pointing components.
        double eff_grad_na = grad[0];
        double eff_grad_nb = grad[1];

        // If at upper bound (n=2) and gradient wants to increase (negative gradient), set to 0
        if (n_a >= 2.0 - 1e-10 && grad[0] < 0.0) eff_grad_na = 0.0;
        if (n_b >= 2.0 - 1e-10 && grad[1] < 0.0) eff_grad_nb = 0.0;
        // If at lower bound (n=0) and gradient wants to decrease (positive gradient), set to 0
        if (n_a <= 1e-10 && grad[0] > 0.0) eff_grad_na = 0.0;
        if (n_b <= 1e-10 && grad[1] > 0.0) eff_grad_nb = 0.0;

        // Check convergence using effective gradient
        double eff_grad_norm = std::sqrt(eff_grad_na*eff_grad_na + eff_grad_nb*eff_grad_nb);
        if (eff_grad_norm < tol * 10.0) {
            if (reks_debug_ >= 2) {
                double E_PPS = compute_E_PPS();
                outfile->Printf("    iter %2d: n_a=%.6f n_b=%.6f E_PPS=%.10f |eff_grad|=%.2e CONVERGED (boundary)\n",
                               iter, n_a, n_b, E_PPS, eff_grad_norm);
            }
            break;
        }

        // Compute 2x2 Hessian [H_aa, H_ab; H_ba, H_bb] - use PPS Hessian for consistency
        std::vector<double> hess = active_space_->compute_hessian_PPS(E_micro_);

        double H_aa = hess[0];
        double H_ab = hess[1];
        double H_bb = hess[3];

        // Compute determinant for 2x2 matrix inversion
        double det = H_aa * H_bb - H_ab * H_ab;

        double delta_na, delta_nb;

        // Check if at boundary where Hessian becomes singular.
        // At n_a=2 or n_a=0, x = n_a*n_d = 0, so df/dx = 0 and H_aa = 0.
        // In this case, we should do 1D optimization in the other direction.
        //
        // However, we only "stay at boundary" if the gradient points INTO the boundary.
        // If gradient points AWAY from boundary (wants to leave), we should allow that step.
        bool at_upper_a = (n_a >= 2.0 - 1e-8);
        bool at_lower_a = (n_a <= 1e-8);
        bool at_upper_b = (n_b >= 2.0 - 1e-8);
        bool at_lower_b = (n_b <= 1e-8);

        // Check if constrained at boundary (can't move in that direction)
        // At upper bound: constrained if gradient is negative (wants to increase)
        // At lower bound: constrained if gradient is positive (wants to decrease)
        bool constrained_a = (at_upper_a && grad[0] < 0.0) || (at_lower_a && grad[0] > 0.0);
        bool constrained_b = (at_upper_b && grad[1] < 0.0) || (at_lower_b && grad[1] > 0.0);

        // Check if Hessian is singular in a direction (happens at boundary with x=0)
        bool singular_a = std::abs(H_aa) < reks::constants::HESSIAN_THRESHOLD;
        bool singular_b = std::abs(H_bb) < reks::constants::HESSIAN_THRESHOLD;

        if ((constrained_a || singular_a) && !singular_b && std::abs(H_bb) > reks::constants::HESSIAN_THRESHOLD) {
            // n_a direction blocked or singular: do 1D Newton-Raphson in n_b direction only
            delta_na = 0.0;
            delta_nb = -grad[1] / H_bb;  // 1D Newton step

            // =========================================================================
            // SADDLE POINT ESCAPE: When n_a ≈ 1.0 (where n_a = n_d) and H_aa is singular,
            // we're at a saddle point. The gradient is zero by symmetry but we need to
            // escape. Add a perturbation that the line search will validate.
            // =========================================================================
            const double saddle_threshold = 0.15;  // Detection threshold around n_a = 1.0
            const double saddle_perturbation = 0.2;  // Perturbation size to escape
            if (singular_a && !constrained_a && std::abs(n_a - 1.0) < saddle_threshold) {
                // At saddle point: n_a ≈ 1.0 with H_aa ≈ 0
                // Try moving n_a toward 2.0 first (bonding character)
                // The line search with two-pass will find the correct direction
                delta_na = saddle_perturbation;
                used_saddle_escape_a = true;  // Mark for bypassing fon_max_delta
                if (reks_debug_ >= 1) {
                    outfile->Printf("    [Saddle escape: n_a=%.4f near 1.0, H_aa=%.4e singular, adding delta_na=%.4f]\n",
                                   n_a, H_aa, delta_na);
                }
            }

            if (reks_debug_ >= 2) {
                outfile->Printf("    [1D Newton in n_b: H_bb=%.4e, grad_b=%.4e, constrained_a=%d, singular_a=%d]\n",
                               H_bb, grad[1], constrained_a, singular_a);
            }
        } else if ((constrained_b || singular_b) && !singular_a && std::abs(H_aa) > reks::constants::HESSIAN_THRESHOLD) {
            // n_b direction blocked or singular: do 1D Newton-Raphson in n_a direction only
            delta_na = -grad[0] / H_aa;  // 1D Newton step
            delta_nb = 0.0;

            // SADDLE POINT ESCAPE for n_b
            const double saddle_threshold = 0.15;
            const double saddle_perturbation = 0.2;
            if (singular_b && !constrained_b && std::abs(n_b - 1.0) < saddle_threshold) {
                delta_nb = saddle_perturbation;
                used_saddle_escape_b = true;  // Mark for bypassing fon_max_delta
                if (reks_debug_ >= 1) {
                    outfile->Printf("    [Saddle escape: n_b=%.4f near 1.0, H_bb=%.4e singular, adding delta_nb=%.4f]\n",
                                   n_b, H_bb, delta_nb);
                }
            }

            if (reks_debug_ >= 2) {
                outfile->Printf("    [1D Newton in n_a: H_aa=%.4e, grad_a=%.4e, constrained_b=%d, singular_b=%d]\n",
                               H_aa, grad[0], constrained_b, singular_b);
            }
        } else if (std::abs(det) < reks::constants::HESSIAN_THRESHOLD || det < 0.0) {
            // =========================================================================
            // Indefinite or singular Hessian: Coupling-Based Direction (Phase 1)
            // =========================================================================
            // When Hessian is indefinite, we're at/near a saddle point. The gradient
            // may point to either local minimum (diradical n≈1 or closed-shell n≈2).
            //
            // KEY INSIGHT: Coupling elements tell us which minimum is PHYSICALLY correct:
            // - Strong coupling (|dad| > threshold): diradical solution is correct
            // - Weak coupling (|dad| < threshold): closed-shell may be acceptable
            //
            // Strategy:
            // 1. If strong coupling: force movement towards n=1.0 regardless of gradient
            // 2. If weak coupling: follow gradient direction (may go to n=2.0)
            //
            // This prevents getting trapped in the WRONG minimum.

            delta_na = 0.0;
            delta_nb = 0.0;

            const double grad_step = 0.05;  // Larger step for faster escape from saddle

            // Handle n_a based on coupling strength
            if (prefer_diradical_a) {
                // Strong (a,d) coupling: force towards n_a = 1.0
                if (n_a > 1.0) {
                    delta_na = -grad_step;  // Decrease towards 1.0
                } else if (n_a < 1.0) {
                    delta_na = grad_step;   // Increase towards 1.0
                }
            } else {
                // Weak (a,d) coupling: follow gradient direction
                // gradient > 0 means dE/dn > 0, so decreasing n lowers E
                // gradient < 0 means dE/dn < 0, so increasing n lowers E
                if (grad[0] > 0.0) {
                    delta_na = -grad_step;  // Move in gradient descent direction
                } else if (grad[0] < 0.0) {
                    delta_na = grad_step;
                }
            }

            // Handle n_b based on coupling strength
            if (prefer_diradical_b) {
                // Strong (b,c) coupling: force towards n_b = 1.0
                if (n_b > 1.0) {
                    delta_nb = -grad_step;  // Decrease towards 1.0
                } else if (n_b < 1.0) {
                    delta_nb = grad_step;   // Increase towards 1.0
                }
            } else {
                // Weak (b,c) coupling: follow gradient direction
                if (grad[1] > 0.0) {
                    delta_nb = -grad_step;
                } else if (grad[1] < 0.0) {
                    delta_nb = grad_step;
                }
            }

            if (reks_debug_ >= 2) {
                outfile->Printf("    [Hessian indef: det=%.4e, grad=(%.4e,%.4e)]\n",
                               det, grad[0], grad[1]);
                outfile->Printf("    [Coupling-based: prefer_birad_a=%d prefer_birad_b=%d, delta=(%.4f,%.4f)]\n",
                               prefer_diradical_a, prefer_diradical_b, delta_na, delta_nb);
            }
        } else {
            // Newton-Raphson: delta = -H^{-1} * grad
            // For 2x2: [a b; b d]^{-1} = (1/det) * [d -b; -b a]
            delta_na = -(H_bb * grad[0] - H_ab * grad[1]) / det;
            delta_nb = -(-H_ab * grad[0] + H_aa * grad[1]) / det;
        }

        // =====================================================================
        // Line Search with Golden Ratio (from Filatov's original GAMESS code)
        // =====================================================================
        // Instead of blindly applying the Newton step, check if it actually
        // decreases the energy. If not, reduce the step size.

        double n_a_new, n_b_new;

        if (reks_line_search_) {
            const double golden = reks::constants::GOLDEN_RATIO;
            const int max_ls_iter = reks::constants::LINE_SEARCH_MAX_ITER;

            // Apply initial damping to limit step size
            if (std::abs(delta_na) > max_step) {
                delta_na = std::copysign(max_step, delta_na);
            }
            if (std::abs(delta_nb) > max_step) {
                delta_nb = std::copysign(max_step, delta_nb);
            }

            // Energy before the step
            // CRITICAL: Use E_PPS (GVB-PP energy) as target function, NOT E_SA!
            // C-code (reks_solver.cpp) optimizes PPS energy in line 493-503.
            // SA-REKS minimum can be at closed-shell (n=2.0), but PPS minimum
            // is at the correlated solution (n~1.5 for stretched geometries).
            double E_PPS_old = compute_E_PPS();

            double factor = 1.0;
            int pass = 0;  // 0 = positive direction, 1 = negative direction
            bool step_accepted = false;

            for (int ls_iter = 0; ls_iter < max_ls_iter; ++ls_iter) {
                // Try step with current factor
                double n_a_try = std::clamp(n_a + delta_na * factor, 0.0, 2.0);
                double n_b_try = std::clamp(n_b + delta_nb * factor, 0.0, 2.0);

                // Temporarily set FONs to compute energy
                active_space_->set_pair_fon(0, n_a_try);
                active_space_->set_pair_fon(1, n_b_try);
                compute_weighting_factors();

                double E_PPS_new = compute_E_PPS();

                if (E_PPS_new < E_PPS_old) {
                    // Energy decreased - accept this step
                    n_a_new = n_a_try;
                    n_b_new = n_b_try;
                    step_accepted = true;
                    if (reks_debug_ >= 2) {
                        outfile->Printf("    [Line search: factor=%.4f, E_PPS_old=%.10f, E_PPS_new=%.10f, accepted]\n",
                                       factor, E_PPS_old, E_PPS_new);
                    }
                    break;
                }

                // Energy increased - reduce step size
                factor /= (golden * golden);

                // Check if step is too small
                double step_norm = std::abs(factor) * std::sqrt(delta_na*delta_na + delta_nb*delta_nb);
                if (step_norm < tol) {
                    // Try negative direction (two-pass line search)
                    if (pass == 0) {
                        factor = -1.0;
                        pass = 1;
                        if (reks_debug_ >= 2) {
                            outfile->Printf("    [Line search: trying negative direction]\n");
                        }
                    } else {
                        // Both directions failed - take zero step
                        if (reks_debug_ >= 2) {
                            outfile->Printf("    [Line search: both directions failed, keeping current FON]\n");
                        }
                        n_a_new = n_a;
                        n_b_new = n_b;
                        step_accepted = true;
                        break;
                    }
                }

                // Two-pass: after half the iterations in positive direction, try negative
                if (ls_iter + 1 >= max_ls_iter / 2 && pass == 0) {
                    factor = -1.0;
                    pass = 1;
                    if (reks_debug_ >= 2) {
                        outfile->Printf("    [Line search: switching to negative direction]\n");
                    }
                }
            }

            // If loop completed without acceptance, keep current FON
            if (!step_accepted) {
                n_a_new = n_a;
                n_b_new = n_b;
                if (reks_debug_ >= 2) {
                    outfile->Printf("    [Line search: max iterations reached, keeping current FON]\n");
                }
            }

            // Restore FONs to the accepted values (may have changed during line search)
            active_space_->set_pair_fon(0, n_a_new);
            active_space_->set_pair_fon(1, n_b_new);
            compute_weighting_factors();

        } else {
            // Original simple damping (no line search)
            if (std::abs(delta_na) > max_step) {
                delta_na = std::copysign(max_step, delta_na);
            }
            if (std::abs(delta_nb) > max_step) {
                delta_nb = std::copysign(max_step, delta_nb);
            }

            // Update FONs with bounds [0.0, 2.0]
            n_a_new = std::clamp(n_a + delta_na, 0.0, 2.0);
            n_b_new = std::clamp(n_b + delta_nb, 0.0, 2.0);
        }

        // Actual deltas after bounds enforcement
        double actual_delta_na = n_a_new - n_a;
        double actual_delta_nb = n_b_new - n_b;

        // =====================================================================
        // Oscillation Detection
        // =====================================================================
        // Check for oscillation ONLY when using gradient descent (not Newton)
        bool using_gradient_descent = (std::abs(det) < reks::constants::HESSIAN_THRESHOLD || det < 0.0) &&
                                       std::abs(actual_delta_na) > tol && std::abs(actual_delta_nb) > tol;

        if (using_gradient_descent && prev2_n_a >= 0.0) {
            double diff_a = std::abs(n_a_new - prev2_n_a);
            double diff_b = std::abs(n_b_new - prev2_n_b);
            if (diff_a < osc_threshold && diff_b < osc_threshold) {
                oscillation_count++;
                if (oscillation_count >= 3) {
                    // Oscillation detected! Take midpoint
                    n_a_new = (n_a_new + prev_n_a) / 2.0;
                    n_b_new = (n_b_new + prev_n_b) / 2.0;
                    if (reks_debug_ >= 1) {
                        outfile->Printf("    Oscillation detected! Midpoint: n_a=%.6f n_b=%.6f\n",
                                       n_a_new, n_b_new);
                    }
                    oscillation_count = 0;
                }
            } else {
                oscillation_count = 0;
            }
        }

        // Update history
        prev2_n_a = prev_n_a;
        prev2_n_b = prev_n_b;
        prev_n_a = n_a;
        prev_n_b = n_b;

        n_a = n_a_new;
        n_b = n_b_new;

        // Update active space FONs
        active_space_->set_pair_fon(0, n_a);
        active_space_->set_pair_fon(1, n_b);

        // Recompute weights
        compute_weighting_factors();

        if (reks_debug_ >= 2) {
            double E_PPS = compute_E_PPS();
            outfile->Printf("    iter %2d: n_a=%.6f n_b=%.6f E_PPS=%.10f |grad|=%.2e |eff_grad|=%.2e det=%.2e\n",
                           iter, n_a, n_b, E_PPS,
                           std::sqrt(grad[0]*grad[0] + grad[1]*grad[1]), eff_grad_norm, det);
        }

        // Convergence check on actual step
        if (std::abs(actual_delta_na) < tol && std::abs(actual_delta_nb) < tol) {
            break;
        }
    }

    // Save previous values for later use
    double saved_prev_n_a = prev_n_a_;
    double saved_prev_n_b = prev_n_b_;

    // =========================================================================
    // Apply FON change limiting for DIIS stability
    // =========================================================================
    // Limit the change from previous SCF iteration to fon_max_delta_
    // This ensures Fock matrices remain similar between DIIS iterations

    double delta_a_from_prev = n_a - saved_prev_n_a;
    double delta_b_from_prev = n_b - saved_prev_n_b;

    bool limited = false;
    if (!fon_oscillation_frozen_) {  // Don't limit if frozen
        // Don't limit if saddle escape was used (needs full perturbation to escape)
        if (std::abs(delta_a_from_prev) > fon_max_delta_ && !used_saddle_escape_a) {
            n_a = saved_prev_n_a + std::copysign(fon_max_delta_, delta_a_from_prev);
            n_a = std::clamp(n_a, 0.0, 2.0);
            limited = true;
        }
        if (std::abs(delta_b_from_prev) > fon_max_delta_ && !used_saddle_escape_b) {
            n_b = saved_prev_n_b + std::copysign(fon_max_delta_, delta_b_from_prev);
            n_b = std::clamp(n_b, 0.0, 2.0);
            limited = true;
        }
    }

    if (limited) {
        // Update active space with limited FONs
        active_space_->set_pair_fon(0, n_a);
        active_space_->set_pair_fon(1, n_b);
        compute_weighting_factors();

        if (reks_debug_ >= 1) {
            outfile->Printf("  FON change limited: n_a: %.6f->%.6f (delta %.4f), n_b: %.6f->%.6f (delta %.4f)\n",
                           saved_prev_n_a, n_a, n_a - saved_prev_n_a, saved_prev_n_b, n_b, n_b - saved_prev_n_b);
        }
    }

    // =========================================================================
    // Orbital Swap for Convergence (from Filatov's original GAMESS code)
    // =========================================================================
    // "Flipped occupations (a stronger occupied orb is above) are bad for convergence!"
    // When x < 0.5 (n < 1.0), the virtual orbital is more occupied than the occupied one.
    // Swap them to put the more occupied orbital in the "occupied" position.
    const double swap_tol = 1.0 - 1e-4;  // n < 1.0 means swap needed

    if (n_a < swap_tol) {
        // Swap orbitals a and d
        // Active orbital indices: a=0, b=1, c=2, d=3 (relative to Ncore_ start)
        int a_idx = Ncore_;
        int d_idx = Ncore_ + 3;

        // Swap columns in Ca_
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, a_idx);
            Ca_->set(i, a_idx, Ca_->get(i, d_idx));
            Ca_->set(i, d_idx, tmp);
        }

        // Transform FON: n_a -> 2 - n_a
        n_a = 2.0 - n_a;
        active_space_->set_pair_fon(0, n_a);

        if (reks_debug_ >= 1) {
            outfile->Printf("  [Orbital swap] n_a < 1.0, swapped orbitals %d<->%d, n_a -> %.6f\n",
                           a_idx + 1, d_idx + 1, n_a);
        }
    }

    if (n_b < swap_tol) {
        // Swap orbitals b and c
        int b_idx = Ncore_ + 1;
        int c_idx = Ncore_ + 2;

        // Swap columns in Ca_
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, b_idx);
            Ca_->set(i, b_idx, Ca_->get(i, c_idx));
            Ca_->set(i, c_idx, tmp);
        }

        // Transform FON: n_b -> 2 - n_b
        n_b = 2.0 - n_b;
        active_space_->set_pair_fon(1, n_b);

        if (reks_debug_ >= 1) {
            outfile->Printf("  [Orbital swap] n_b < 1.0, swapped orbitals %d<->%d, n_b -> %.6f\n",
                           b_idx + 1, c_idx + 1, n_b);
        }
    }

    // Recompute weights after potential swap
    compute_weighting_factors();

    // =========================================================================
    // FON Oscillation Detection and Handling (AFTER limiting and swap)
    // =========================================================================
    // Detect 2-cycle oscillations: if current LIMITED FON is close to 2-iterations-ago
    // but different from previous iteration, we have an oscillation.
    // When detected, freeze FONs at the average value.

    const double oscillation_threshold = 0.025;  // Slightly larger than fon_max_delta

    if (iteration_ >= 5 && !fon_oscillation_frozen_) {  // Start detection after warmup
        // Check if current (limited) FON is similar to 2-iterations-ago
        double diff_a_2ago = std::abs(n_a - prev_prev_n_a_);
        double diff_b_2ago = std::abs(n_b - prev_prev_n_b_);
        double diff_a_1ago = std::abs(n_a - saved_prev_n_a);
        double diff_b_1ago = std::abs(n_b - saved_prev_n_b);

        bool oscillation_detected = false;
        // True 2-cycle oscillation:
        // - diff(2-ago) is very small (current ≈ 2-ago)
        // - diff(2-ago) << diff(1-ago) (not just monotonic increase)
        // For monotonic: diff(2-ago) = 2*diff(1-ago)
        // For oscillation: diff(2-ago) << diff(1-ago)
        if (diff_a_2ago < oscillation_threshold && diff_a_1ago > oscillation_threshold * 0.3 &&
            diff_a_2ago < diff_a_1ago * 0.5) {  // Key: diff_2ago much smaller than diff_1ago
            oscillation_detected = true;
        }
        if (diff_b_2ago < oscillation_threshold && diff_b_1ago > oscillation_threshold * 0.3 &&
            diff_b_2ago < diff_b_1ago * 0.5) {
            oscillation_detected = true;
        }

        if (oscillation_detected) {
            fon_oscillation_count_++;
            if (reks_debug_ >= 1) {
                outfile->Printf("  FON oscillation detected (count=%d): n_a 2-ago=%.4f, 1-ago=%.4f, now=%.4f\n",
                               fon_oscillation_count_, prev_prev_n_a_, saved_prev_n_a, n_a);
                outfile->Printf("                                      n_b 2-ago=%.4f, 1-ago=%.4f, now=%.4f\n",
                               prev_prev_n_b_, saved_prev_n_b, n_b);
            }

            // If oscillation persists for 3 consecutive iterations, freeze at average
            if (fon_oscillation_count_ >= 3) {
                double avg_n_a = (n_a + saved_prev_n_a + prev_prev_n_a_) / 3.0;
                double avg_n_b = (n_b + saved_prev_n_b + prev_prev_n_b_) / 3.0;
                n_a = avg_n_a;
                n_b = avg_n_b;
                fon_oscillation_frozen_ = true;
                // Update active space with frozen FONs
                active_space_->set_pair_fon(0, n_a);
                active_space_->set_pair_fon(1, n_b);
                compute_weighting_factors();
                if (reks_debug_ >= 1) {
                    outfile->Printf("  FON FROZEN at average: n_a=%.6f, n_b=%.6f\n", n_a, n_b);
                }
            }
        } else {
            fon_oscillation_count_ = 0;  // Reset counter if no oscillation
        }
    }

    // Update history for next iteration
    prev_prev_n_a_ = saved_prev_n_a;
    prev_prev_n_b_ = saved_prev_n_b;
    prev_n_a_ = n_a;
    prev_n_b_ = n_b;

    // =========================================================================
    // OSS1 FON Optimization: Optimize n'_a (oss1_fon_a_) for OSS1 configuration
    // =========================================================================

    double oss1_fon = active_space_->get_oss1_fon_a();

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver(4,4): OSS1 FON Optimization ===\n");
        outfile->Printf("  Initial n'_a=%.6f\n", oss1_fon);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        double grad = active_space_->compute_gradient_OSS1(E_micro_);
        double hess = active_space_->compute_hessian_OSS1(E_micro_);

        // Boundary-aware convergence check
        double eff_grad = grad;
        if (oss1_fon >= 2.0 - 1e-10 && grad < 0.0) eff_grad = 0.0;
        if (oss1_fon <= 1e-10 && grad > 0.0) eff_grad = 0.0;

        if (std::abs(eff_grad) < tol * 10.0) {
            if (reks_debug_ >= 2) {
                outfile->Printf("    OSS1 iter %2d: n'_a=%.6f |eff_grad|=%.2e CONVERGED\n",
                               iter, oss1_fon, std::abs(eff_grad));
            }
            break;
        }

        double delta;
        if (std::abs(hess) < reks::constants::HESSIAN_THRESHOLD || hess > 0.0) {
            // Near-singular or positive (wrong curvature) Hessian: use gradient descent
            delta = -std::copysign(max_step * 0.5, eff_grad);
            if (std::abs(eff_grad) < 1e-10) delta = 0.0;
        } else {
            // Newton-Raphson: delta = -grad/hess
            delta = -grad / hess;
        }

        // Damping
        if (std::abs(delta) > max_step) {
            delta = std::copysign(max_step, delta);
        }

        double oss1_fon_new = std::clamp(oss1_fon + delta, 0.0, 2.0);
        delta = oss1_fon_new - oss1_fon;
        oss1_fon = oss1_fon_new;

        active_space_->set_oss1_fon_a(oss1_fon);

        if (reks_debug_ >= 2) {
            double E_OSS1 = active_space_->compute_energy_OSS1(E_micro_);
            outfile->Printf("    OSS1 iter %2d: n'_a=%.6f E_OSS1=%.10f |grad|=%.2e hess=%.2e\n",
                           iter, oss1_fon, E_OSS1, std::abs(grad), hess);
        }

        if (std::abs(delta) < tol) break;
    }

    // =========================================================================
    // OSS2 FON Optimization: Optimize n'_b (oss2_fon_b_) for OSS2 configuration
    // =========================================================================
    double oss2_fon = active_space_->get_oss2_fon_b();

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver(4,4): OSS2 FON Optimization ===\n");
        outfile->Printf("  Initial n'_b=%.6f\n", oss2_fon);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        double grad = active_space_->compute_gradient_OSS2(E_micro_);
        double hess = active_space_->compute_hessian_OSS2(E_micro_);

        // Boundary-aware convergence check
        double eff_grad = grad;
        if (oss2_fon >= 2.0 - 1e-10 && grad < 0.0) eff_grad = 0.0;
        if (oss2_fon <= 1e-10 && grad > 0.0) eff_grad = 0.0;

        if (std::abs(eff_grad) < tol * 10.0) {
            if (reks_debug_ >= 2) {
                outfile->Printf("    OSS2 iter %2d: n'_b=%.6f |eff_grad|=%.2e CONVERGED\n",
                               iter, oss2_fon, std::abs(eff_grad));
            }
            break;
        }

        double delta;
        if (std::abs(hess) < reks::constants::HESSIAN_THRESHOLD || hess > 0.0) {
            // Near-singular or positive (wrong curvature) Hessian: use gradient descent
            delta = -std::copysign(max_step * 0.5, eff_grad);
            if (std::abs(eff_grad) < 1e-10) delta = 0.0;
        } else {
            // Newton-Raphson: delta = -grad/hess
            delta = -grad / hess;
        }

        // Damping
        if (std::abs(delta) > max_step) {
            delta = std::copysign(max_step, delta);
        }

        double oss2_fon_new = std::clamp(oss2_fon + delta, 0.0, 2.0);
        delta = oss2_fon_new - oss2_fon;
        oss2_fon = oss2_fon_new;

        active_space_->set_oss2_fon_b(oss2_fon);

        if (reks_debug_ >= 2) {
            double E_OSS2 = active_space_->compute_energy_OSS2(E_micro_);
            outfile->Printf("    OSS2 iter %2d: n'_b=%.6f E_OSS2=%.10f |grad|=%.2e hess=%.2e\n",
                           iter, oss2_fon, E_OSS2, std::abs(grad), hess);
        }

        if (std::abs(delta) < tol) break;
    }

    // =========================================================================
    // Recompute combined SA weights after all FON optimizations
    // =========================================================================
    compute_weighting_factors();

    if (reks_debug_ >= 1) {
        double E_PPS = active_space_->compute_energy_PPS(E_micro_);
        double E_OSS1 = active_space_->compute_energy_OSS1(E_micro_);
        double E_OSS2 = active_space_->compute_energy_OSS2(E_micro_);
        double w_pps = active_space_->w_PPS();
        double w_oss = active_space_->w_OSS();
        double E_SA = w_pps * E_PPS + w_oss * E_OSS1 + w_oss * E_OSS2;

        double n_d = active_space_->pair(0).fon_q;
        double n_c = active_space_->pair(1).fon_q;
        outfile->Printf("\n  === 3SA-REKS(4,4) Final Results ===\n");
        outfile->Printf("  PPS FONs:  n_a=%12.8f, n_d=%12.8f, n_b=%12.8f, n_c=%12.8f\n",
                        n_a, n_d, n_b, n_c);
        outfile->Printf("  OSS1 FON:  n'_a=%12.8f (n'_d=%12.8f)\n",
                        oss1_fon, 2.0 - oss1_fon);
        outfile->Printf("  OSS2 FON:  n'_b=%12.8f (n'_c=%12.8f)\n",
                        oss2_fon, 2.0 - oss2_fon);
        outfile->Printf("  E_PPS  = %20.12f\n", E_PPS);
        outfile->Printf("  E_OSS1 = %20.12f\n", E_OSS1);
        outfile->Printf("  E_OSS2 = %20.12f\n", E_OSS2);
        outfile->Printf("  E_SA   = %20.12f  (w_PPS=%.4f, w_OSS=%.4f)\n", E_SA, w_pps, w_oss);
    }
}

// ===========================================================================
// Trust-Region FON Optimization Methods
// ===========================================================================
// Trust-Region Augmented Hessian (TRAH) method for robust FON optimization.
// Handles indefinite Hessian at saddle points via secular equation.
//
// Reference: J. Chem. Phys. 156, 204104 (2022)
// ===========================================================================

std::vector<double> REKS::solve_boundary_constrained_2d(
    double g1, double g2, double lam1, double lam2,
    const std::vector<double>& v1, const std::vector<double>& v2, double Delta)
{
    // Solve secular equation for boundary-constrained step:
    //   ||p(μ)||² = Δ²  where p_i(μ) = -g_i / (λ_i + μ)
    //
    // The secular equation is: φ(μ) = Σ g_i²/(λ_i+μ)² - Δ² = 0
    //
    // Newton iteration on φ(μ):
    //   φ'(μ) = -2 Σ g_i²/(λ_i+μ)³
    //   μ_new = μ - φ(μ)/φ'(μ)

    using namespace reks::constants;

    // Regularization: if both gradients are tiny, return zero step
    const double g_norm_sq = g1*g1 + g2*g2;
    if (g_norm_sq < 1e-20) {
        return {0.0, 0.0};
    }

    // Find lower bound for μ: must satisfy λ_i + μ > 0 for all i
    // For indefinite Hessian, at least one λ < 0, so μ_lower > -min(λ)
    double lam_min = std::min(lam1, lam2);
    double mu_lower = std::max(0.0, -lam_min) + 1e-10;

    // Initial guess for μ
    double mu = mu_lower + 0.1;

    // Newton iteration to solve secular equation
    for (int iter = 0; iter < TR_MAX_SECULAR_ITER; ++iter) {
        double d1 = lam1 + mu;
        double d2 = lam2 + mu;

        // Safety: avoid division by very small numbers
        if (std::abs(d1) < 1e-15) d1 = std::copysign(1e-15, d1);
        if (std::abs(d2) < 1e-15) d2 = std::copysign(1e-15, d2);

        // φ(μ) = g1²/d1² + g2²/d2² - Δ²
        double phi = (g1*g1)/(d1*d1) + (g2*g2)/(d2*d2) - Delta*Delta;

        // Check convergence
        if (std::abs(phi) < TR_SECULAR_TOL) {
            break;
        }

        // φ'(μ) = -2*(g1²/d1³ + g2²/d2³)
        double dphi = -2.0 * ((g1*g1)/(d1*d1*d1) + (g2*g2)/(d2*d2*d2));

        // Avoid division by zero in Newton step
        if (std::abs(dphi) < 1e-20) {
            break;
        }

        // Newton step
        double dmu = -phi / dphi;
        mu += dmu;

        // Keep μ above the pole
        if (mu < mu_lower) {
            mu = mu_lower;
        }
    }

    // Compute step in eigenspace
    double d1 = lam1 + mu;
    double d2 = lam2 + mu;
    if (std::abs(d1) < 1e-15) d1 = std::copysign(1e-15, d1);
    if (std::abs(d2) < 1e-15) d2 = std::copysign(1e-15, d2);

    double p_eigen1 = -g1 / d1;
    double p_eigen2 = -g2 / d2;

    // Transform back to original space: p = V * p_eigen
    // where V = [v1, v2] (columns are eigenvectors)
    double delta_na = v1[0] * p_eigen1 + v2[0] * p_eigen2;
    double delta_nb = v1[1] * p_eigen1 + v2[1] * p_eigen2;

    return {delta_na, delta_nb};
}

std::vector<double> REKS::solve_trust_region_subproblem_2d(
    const std::vector<double>& g, const std::vector<std::vector<double>>& H, double Delta)
{
    // Solve Trust-Region subproblem for 2D case:
    //   minimize m(p) = g·p + ½ p·H·p
    //   subject to ||p|| ≤ Δ
    //
    // Cases:
    // 1. H positive definite & Newton step inside trust region -> use Newton
    // 2. H positive definite & Newton step outside -> scale to boundary
    // 3. H indefinite (det < 0) -> solve on boundary via secular equation
    // 4. H singular -> gradient descent direction

    using namespace reks::constants;

    // Extract Hessian elements
    double H11 = H[0][0], H12 = H[0][1], H22 = H[1][1];
    double g1 = g[0], g2 = g[1];

    // Compute eigendecomposition of 2x2 symmetric Hessian
    // H = V * Λ * V^T where Λ = diag(lam1, lam2)
    double trace = H11 + H22;
    double det = H11 * H22 - H12 * H12;

    // Eigenvalues: λ = (trace ± sqrt(trace² - 4*det)) / 2
    double discriminant = trace * trace - 4.0 * det;
    if (discriminant < 0.0) discriminant = 0.0;  // Numerical safeguard
    double sqrt_disc = std::sqrt(discriminant);

    double lam1 = 0.5 * (trace - sqrt_disc);  // Smaller eigenvalue
    double lam2 = 0.5 * (trace + sqrt_disc);  // Larger eigenvalue

    // Eigenvectors for 2x2 symmetric matrix
    std::vector<double> v1(2), v2(2);
    if (std::abs(H12) > 1e-14) {
        // Off-diagonal is nonzero: standard eigenvector formula
        v1[0] = lam1 - H22;
        v1[1] = H12;
        double norm1 = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
        if (norm1 > NORM_THRESHOLD) {
            v1[0] /= norm1;
            v1[1] /= norm1;
        } else {
            v1[0] = 1.0; v1[1] = 0.0;
        }

        v2[0] = lam2 - H22;
        v2[1] = H12;
        double norm2 = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
        if (norm2 > NORM_THRESHOLD) {
            v2[0] /= norm2;
            v2[1] /= norm2;
        } else {
            v2[0] = 0.0; v2[1] = 1.0;
        }
    } else {
        // Diagonal matrix: eigenvectors are coordinate axes
        if (H11 <= H22) {
            v1[0] = 1.0; v1[1] = 0.0;
            v2[0] = 0.0; v2[1] = 1.0;
        } else {
            v1[0] = 0.0; v1[1] = 1.0;
            v2[0] = 1.0; v2[1] = 0.0;
            std::swap(lam1, lam2);
        }
    }

    // Transform gradient to eigenspace: g_eigen = V^T * g
    double g_eigen1 = v1[0]*g1 + v1[1]*g2;
    double g_eigen2 = v2[0]*g1 + v2[1]*g2;

    // Case 3: Indefinite Hessian (det < 0, meaning lam1 < 0 < lam2)
    if (det < -HESSIAN_THRESHOLD) {
        // Must solve on the boundary using secular equation
        return solve_boundary_constrained_2d(g_eigen1, g_eigen2, lam1, lam2, v1, v2, Delta);
    }

    // Case 4: Singular or near-singular Hessian
    if (std::abs(det) < HESSIAN_THRESHOLD || lam1 < HESSIAN_THRESHOLD) {
        // Use gradient descent direction, scaled to trust radius
        double g_norm = std::sqrt(g1*g1 + g2*g2);
        if (g_norm < 1e-14) {
            return {0.0, 0.0};
        }
        return {-Delta * g1 / g_norm, -Delta * g2 / g_norm};
    }

    // Cases 1 & 2: H is positive definite
    // Compute Newton step: p_N = -H^{-1} * g
    double inv_det = 1.0 / det;
    double p_newton_a = -inv_det * (H22 * g1 - H12 * g2);
    double p_newton_b = -inv_det * (-H12 * g1 + H11 * g2);

    double newton_norm = std::sqrt(p_newton_a*p_newton_a + p_newton_b*p_newton_b);

    // Case 1: Newton step is inside trust region
    if (newton_norm <= Delta) {
        return {p_newton_a, p_newton_b};
    }

    // Case 2: Newton step is outside trust region
    // Scale to the boundary: p = (Δ/||p_N||) * p_N
    // Note: For strictly convex case, this is not optimal (should use secular equation)
    // but it's a reasonable approximation that avoids the complexity
    return solve_boundary_constrained_2d(g_eigen1, g_eigen2, lam1, lam2, v1, v2, Delta);
}

void REKS::trust_region_fon_optimizer_44() {
    // Trust-Region FON optimization for REKS(4,4)
    //
    // This method replaces the Newton-Raphson optimizer with a Trust-Region
    // approach that handles indefinite Hessians robustly via the secular equation.
    //
    // Algorithm:
    // 1. Compute gradient and Hessian of E_PPS w.r.t. (n_a, n_b)
    // 2. Solve Trust-Region subproblem: min m(p) = g·p + ½p·H·p  s.t. ||p|| ≤ Δ
    // 3. Evaluate actual vs predicted reduction
    // 4. Accept/reject step and update trust radius
    //
    // Reference: Trust-Region Augmented Hessian (TRAH)

    using namespace reks::constants;

    // Skip if FONs are frozen
    if (fon_oscillation_frozen_) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === Trust-Region REKS(4,4): FON frozen, skipping ===\n");
        }
        compute_weighting_factors();
        return;
    }

    // Get current FONs
    double n_a = active_space_->pair(0).fon_p;
    double n_b = active_space_->pair(1).fon_p;

    // =========================================================================
    // FON Initialization (same logic as rex_solver_44)
    // =========================================================================
    if (!fon_initialized_44_) {
        // Before localization (iteration_ < 2): keep FONs frozen at n=1.0
        if (iteration_ < 2) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  [TR] Waiting for localization (iteration %d < 2), keeping n_a=n_b=1.0\n", iteration_);
            }
            active_space_->set_pair_fon(0, 1.0);
            active_space_->set_pair_fon(1, 1.0);
            compute_weighting_factors();
            return;
        }

        // After localization: initialize FON
        // Priority: (1) User-provided, (2) FON history prediction, (3) HAM-DIAG
        bool initialized = false;

        // (1) User-provided initial FON values
        if (reks_initial_na_ >= 0.0 && reks_initial_nb_ >= 0.0) {
            active_space_->set_pair_fon(0, reks_initial_na_);
            active_space_->set_pair_fon(1, reks_initial_nb_);
            n_a = reks_initial_na_;
            n_b = reks_initial_nb_;
            initialized = true;
            if (reks_debug_ >= 1) {
                outfile->Printf("  [TR] Using user-provided initial FON: n_a=%.6f, n_b=%.6f\n",
                               reks_initial_na_, reks_initial_nb_);
            }
        }

        // (2) Default: HAM-DIAG initialization
        if (!initialized) {
            initialize_fon_from_hamiltonian_44();
            n_a = active_space_->pair(0).fon_p;
            n_b = active_space_->pair(1).fon_p;
        }
        fon_initialized_44_ = true;
    }

    // Lambda to compute E_PPS at given FON values
    auto compute_E_PPS_at = [&](double na, double nb) -> double {
        double old_na = active_space_->pair(0).fon_p;
        double old_nb = active_space_->pair(1).fon_p;
        active_space_->set_pair_fon(0, na);
        active_space_->set_pair_fon(1, nb);
        double E = active_space_->compute_energy_PPS(E_micro_);
        active_space_->set_pair_fon(0, old_na);
        active_space_->set_pair_fon(1, old_nb);
        return E;
    };

    // Trust-Region parameters (persistent radius from member variable)
    double Delta = trust_radius_44_;

    // Store initial FON values for summary
    double n_a_initial = n_a;
    double n_b_initial = n_b;
    int fon_iters_used = 0;
    bool hessian_was_indefinite = false;

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  ================================================================================\n");
        outfile->Printf("  REKS(4,4) FON Optimizer (Trust-Region) - SCF Iteration %d\n", iteration_);
        outfile->Printf("  ================================================================================\n");
        outfile->Printf("  Initial: n_a = %.6f, n_b = %.6f (n_d = %.6f, n_c = %.6f)\n",
                       n_a, n_b, 2.0 - n_a, 2.0 - n_b);
        outfile->Printf("  Trust radius: Δ = %.4f\n", Delta);
    }

    double E_old = compute_E_PPS_at(n_a, n_b);

    if (reks_debug_ >= 1) {
        outfile->Printf("  Initial E_PPS = %.10f Ha\n\n", E_old);
    }

    for (int iter = 0; iter < FON_MAX_ITER; ++iter) {
        // 1. Compute gradient [dE/dn_a, dE/dn_b]
        std::vector<double> grad = active_space_->compute_gradient_PPS(E_micro_);

        // Boundary-aware effective gradient
        double eff_grad_na = grad[0];
        double eff_grad_nb = grad[1];

        // At upper bound (n≈2) with negative gradient: constrained
        if (n_a >= 2.0 - 1e-10 && grad[0] < 0.0) eff_grad_na = 0.0;
        if (n_b >= 2.0 - 1e-10 && grad[1] < 0.0) eff_grad_nb = 0.0;
        // At lower bound (n≈0) with positive gradient: constrained
        if (n_a <= 1e-10 && grad[0] > 0.0) eff_grad_na = 0.0;
        if (n_b <= 1e-10 && grad[1] > 0.0) eff_grad_nb = 0.0;

        double eff_grad_norm = std::sqrt(eff_grad_na*eff_grad_na + eff_grad_nb*eff_grad_nb);

        // 2. Check convergence
        if (eff_grad_norm < FON_TOL * 10.0) {
            fon_iters_used = iter;
            if (reks_debug_ >= 1) {
                outfile->Printf("  --- FON iter %2d: CONVERGED ---\n", iter);
                outfile->Printf("    n = (%.6f, %.6f), |eff_grad| = %.2e < %.2e\n",
                               n_a, n_b, eff_grad_norm, FON_TOL * 10.0);
            }
            break;
        }

        // 3. Compute Hessian [H_aa, H_ab; H_ab, H_bb]
        std::vector<double> hess_flat = active_space_->compute_hessian_PPS(E_micro_);
        std::vector<std::vector<double>> H = {
            {hess_flat[0], hess_flat[1]},
            {hess_flat[2], hess_flat[3]}
        };

        // Compute Hessian determinant and eigenvalues for diagnostics
        double det_H = H[0][0] * H[1][1] - H[0][1] * H[1][0];
        double trace_H = H[0][0] + H[1][1];
        // Eigenvalues: λ = (trace ± sqrt(trace² - 4*det)) / 2
        double discriminant = trace_H * trace_H - 4.0 * det_H;
        double lambda1 = 0.0, lambda2 = 0.0;
        if (discriminant >= 0.0) {
            lambda1 = 0.5 * (trace_H - std::sqrt(discriminant));
            lambda2 = 0.5 * (trace_H + std::sqrt(discriminant));
        }
        bool is_positive_definite = (lambda1 > 0.0 && lambda2 > 0.0);
        if (!is_positive_definite) hessian_was_indefinite = true;

        if (reks_debug_ >= 1) {
            outfile->Printf("  --- FON iter %2d ---\n", iter);
            outfile->Printf("    n = (%.6f, %.6f), E_PPS = %.10f\n", n_a, n_b, E_old);
            outfile->Printf("    grad = [%+.6e, %+.6e], |g| = %.2e\n", grad[0], grad[1],
                           std::sqrt(grad[0]*grad[0] + grad[1]*grad[1]));
            outfile->Printf("    Hess = [[%+.6f, %+.6f], [%+.6f, %+.6f]]\n",
                           H[0][0], H[0][1], H[1][0], H[1][1]);
            outfile->Printf("    det(H) = %+.6e, eigenvalues = (%.4f, %.4f) %s\n",
                           det_H, lambda1, lambda2,
                           is_positive_definite ? "[positive definite]" : "[INDEFINITE]");
        }

        // 4. Solve Trust-Region subproblem
        std::vector<double> p = solve_trust_region_subproblem_2d({eff_grad_na, eff_grad_nb}, H, Delta);

        // 5. Trial point with boundary clamping
        double n_a_trial = std::clamp(n_a + p[0], 0.0, 2.0);
        double n_b_trial = std::clamp(n_b + p[1], 0.0, 2.0);

        // Actual step (may be smaller due to clamping)
        double actual_p0 = n_a_trial - n_a;
        double actual_p1 = n_b_trial - n_b;
        double step_norm = std::sqrt(actual_p0*actual_p0 + actual_p1*actual_p1);

        // 6. Compute actual vs predicted reduction
        double E_new = compute_E_PPS_at(n_a_trial, n_b_trial);
        double actual_reduction = E_old - E_new;

        // Predicted reduction from quadratic model: -g·p - ½ p·H·p
        double predicted_reduction = -(eff_grad_na * actual_p0 + eff_grad_nb * actual_p1)
                                    - 0.5 * (H[0][0]*actual_p0*actual_p0 + 2*H[0][1]*actual_p0*actual_p1 + H[1][1]*actual_p1*actual_p1);

        // Ratio ρ = actual / predicted
        double rho = (std::abs(predicted_reduction) > 1e-15) ? actual_reduction / predicted_reduction : 0.0;

        if (reks_debug_ >= 1) {
            outfile->Printf("    step p = (%+.6f, %+.6f), |p| = %.4f, Δ = %.4f\n",
                           actual_p0, actual_p1, step_norm, Delta);
            outfile->Printf("    trial n = (%.6f, %.6f), E_trial = %.10f\n", n_a_trial, n_b_trial, E_new);
            outfile->Printf("    actual_red = %+.6e, pred_red = %+.6e, ρ = %.4f\n",
                           actual_reduction, predicted_reduction, rho);
        }

        // 7. Accept/reject step and update trust radius
        const char* step_status = "";
        if (rho > TR_ETA_VERY_GOOD) {
            // Very good step: accept and increase trust radius
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
            Delta = std::min(2.0 * Delta, TR_DELTA_MAX);
            step_status = "VERY_GOOD: accepted, Δ increased";
        } else if (rho > TR_ETA_GOOD) {
            // Good step: accept, keep trust radius
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
            step_status = "GOOD: accepted";
        } else if (rho > TR_ETA_ACCEPT) {
            // Acceptable step: accept but reduce trust radius
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
            Delta = std::max(0.5 * step_norm, TR_DELTA_MIN);
            step_status = "ACCEPTABLE: accepted, Δ reduced";
        } else {
            // Bad step: reject and shrink trust radius
            Delta = std::max(0.25 * step_norm, TR_DELTA_MIN);
            step_status = "REJECTED: Δ shrunk";
        }

        fon_iters_used = iter + 1;

        if (reks_debug_ >= 1) {
            outfile->Printf("    --> %s, new Δ = %.4f\n\n", step_status, Delta);
        }

        // Check for tiny trust radius (algorithm stalled)
        if (Delta < TR_DELTA_MIN * 10.0 && step_norm < FON_TOL) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  Trust radius too small, stopping\n");
            }
            break;
        }
    }

    // Save trust radius for next SCF iteration
    trust_radius_44_ = Delta;

    // Update active space FONs
    active_space_->set_pair_fon(0, n_a);
    active_space_->set_pair_fon(1, n_b);

    // Orbital swap if needed (same logic as rex_solver_44)
    const double swap_tol = 1.0 - 1e-4;
    if (n_a < swap_tol) {
        int a_idx = Ncore_;
        int d_idx = Ncore_ + 3;
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, a_idx);
            Ca_->set(i, a_idx, Ca_->get(i, d_idx));
            Ca_->set(i, d_idx, tmp);
        }
        n_a = 2.0 - n_a;
        active_space_->set_pair_fon(0, n_a);
        if (reks_debug_ >= 1) {
            outfile->Printf("  [TR] Orbital swap a<->d, n_a -> %.6f\n", n_a);
        }
    }
    if (n_b < swap_tol) {
        int b_idx = Ncore_ + 1;
        int c_idx = Ncore_ + 2;
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, b_idx);
            Ca_->set(i, b_idx, Ca_->get(i, c_idx));
            Ca_->set(i, c_idx, tmp);
        }
        n_b = 2.0 - n_b;
        active_space_->set_pair_fon(1, n_b);
        if (reks_debug_ >= 1) {
            outfile->Printf("  [TR] Orbital swap b<->c, n_b -> %.6f\n", n_b);
        }
    }

    // Update prev FON history
    prev_n_a_ = n_a;
    prev_n_b_ = n_b;

    // Recompute weights
    compute_weighting_factors();

    if (reks_debug_ >= 1) {
        double E_PPS = active_space_->compute_energy_PPS(E_micro_);
        outfile->Printf("  --------------------------------------------------------------------------------\n");
        outfile->Printf("  FON Optimization Summary:\n");
        outfile->Printf("    FON iterations:    %d\n", fon_iters_used);
        outfile->Printf("    Hessian status:    %s\n",
                       hessian_was_indefinite ? "INDEFINITE (at some iteration)" : "positive definite");
        outfile->Printf("    Initial n:         (%.6f, %.6f)\n", n_a_initial, n_b_initial);
        outfile->Printf("    Final n:           (%.6f, %.6f)\n", n_a, n_b);
        outfile->Printf("    Change Δn:         (%+.6f, %+.6f)\n", n_a - n_a_initial, n_b - n_b_initial);
        outfile->Printf("    Final E_PPS:       %.10f Ha\n", E_PPS);
        outfile->Printf("    Trust radius:      %.4f\n", trust_radius_44_);
        outfile->Printf("  ================================================================================\n\n");
    }
}

// ---------------------------------------------------------------------------
// Adaptive Trust-Region + Multi-Start (FON_SOLVER=3)
// ---------------------------------------------------------------------------

void REKS::trust_region_fon_optimizer_44_adaptive() {
    // FON_SOLVER=3: Adaptive Multi-Start
    //
    // Strategy:
    // 1. Always use Trust-Region as the primary optimizer
    // 2. Check Hessian after TR optimization
    // 3. Trigger Multi-Start ONLY when needed:
    //    - Hessian is indefinite (det < 0) -> multiple minima possible
    //    - First SCF iteration -> no prior state to trust
    //    - Large FON change (|Δn| > 0.5) -> possible state switching
    //
    // This gives robustness without unnecessary computational cost.

    using namespace reks::constants;

    // Step 1: Run Trust-Region optimization (same as FON_SOLVER=2)
    trust_region_fon_optimizer_44();

    // Step 2: Check if multi-start is needed
    double n_a = active_space_->pair(0).fon_p;
    double n_b = active_space_->pair(1).fon_p;
    double E_current = active_space_->compute_energy_PPS(E_micro_);

    // Compute Hessian to check for multiple minima
    std::vector<double> hess = active_space_->compute_hessian_PPS(E_micro_);
    double det = hess[0] * hess[3] - hess[1] * hess[2];
    double trace = hess[0] + hess[3];
    double discriminant = trace * trace - 4.0 * det;
    double lambda1 = (trace + std::sqrt(std::max(0.0, discriminant))) / 2.0;
    double lambda2 = (trace - std::sqrt(std::max(0.0, discriminant))) / 2.0;

    // Check if multi-start needed
    bool indefinite_hessian = (det < -1e-6);  // Negative determinant = saddle point
    bool first_iteration = (iteration_ <= 2 && !fon_initialized_44_);
    bool large_fon_change = false;
    if (prev_n_a_ >= 0.0) {
        double delta_na = std::abs(n_a - prev_n_a_);
        double delta_nb = std::abs(n_b - prev_n_b_);
        large_fon_change = (delta_na > 0.5 || delta_nb > 0.5);
    }

    bool need_multi_start = indefinite_hessian || first_iteration || large_fon_change;

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === Adaptive Multi-Start Check (FON_SOLVER=3) ===\n");
        outfile->Printf("    Hessian: det=%.6f, λ₁=%.4f, λ₂=%.4f\n", det, lambda1, lambda2);
        outfile->Printf("    Indefinite Hessian: %s\n", indefinite_hessian ? "YES" : "NO");
        outfile->Printf("    First iteration:    %s\n", first_iteration ? "YES" : "NO");
        outfile->Printf("    Large FON change:   %s\n", large_fon_change ? "YES" : "NO");
        outfile->Printf("    Multi-Start trigger: %s\n", need_multi_start ? "TRIGGERED" : "SKIPPED");
    }

    if (!need_multi_start) {
        // Trust-Region solution is good enough
        prev_n_a_ = n_a;
        prev_n_b_ = n_b;
        return;
    }

    // Step 3: Try alternative starting points
    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === Multi-Start Alternative Solutions ===\n");
    }

    // Current solution (from Trust-Region)
    reks::FONSolution current;
    current.n_a = n_a;
    current.n_b = n_b;
    current.E_PPS = E_current;
    current.converged = true;
    current.classify();

    std::vector<reks::FONSolution> solutions;
    solutions.push_back(current);

    // Try diradical starting point (n=1.0, 1.0)
    if (n_a > 1.2 || n_b > 1.2) {  // Only if current is not already diradical
        reks::FONStartPoint birad_start(1.0, 1.0, "DIRADICAL");
        reks::FONSolution birad = optimize_from_start_44(birad_start, 1);
        if (birad.converged) {
            birad.classify();
            solutions.push_back(birad);
        }
    }

    // Try closed-shell starting point (n=2.0, 2.0)
    if (n_a < 1.8 || n_b < 1.8) {  // Only if current is not already closed-shell
        reks::FONStartPoint closed_start(2.0, 2.0, "CLOSED_SHELL");
        reks::FONSolution closed = optimize_from_start_44(closed_start, 2);
        if (closed.converged) {
            closed.classify();
            solutions.push_back(closed);
        }
    }

    // Remove duplicates
    std::vector<reks::FONSolution> unique;
    for (const auto& sol : solutions) {
        bool is_dup = false;
        for (const auto& u : unique) {
            if (sol.is_duplicate_of(u)) {
                is_dup = true;
                break;
            }
        }
        if (!is_dup) {
            unique.push_back(sol);
        }
    }

    if (reks_debug_ >= 1) {
        outfile->Printf("    Solutions found: %zu unique\n", unique.size());
        for (size_t i = 0; i < unique.size(); ++i) {
            const char* char_str = unique[i].character == reks::FONCharacter::DIRADICAL ? "DIRADICAL" :
                                   unique[i].character == reks::FONCharacter::CLOSED_SHELL ? "CLOSED_SHELL" :
                                   "MIXED";
            outfile->Printf("      [%zu] E=%.10f, n=(%.4f, %.4f), %s\n",
                           i, unique[i].E_PPS, unique[i].n_a, unique[i].n_b, char_str);
        }
    }

    // Select best solution (use criterion from options)
    reks::FONSolution selected = select_solution(unique, prev_n_a_, prev_n_b_,
                                                  reks_branch_criterion_, reks_energy_tolerance_);

    // Check if a better solution was found
    if (selected.E_PPS < E_current - 1e-8) {
        if (reks_debug_ >= 1) {
            outfile->Printf("    Multi-Start found LOWER energy solution!\n");
            outfile->Printf("      Previous: E=%.10f, n=(%.4f, %.4f)\n", E_current, n_a, n_b);
            outfile->Printf("      New:      E=%.10f, n=(%.4f, %.4f)\n",
                           selected.E_PPS, selected.n_a, selected.n_b);
        }

        // Apply better solution
        active_space_->set_pair_fon(0, selected.n_a);
        active_space_->set_pair_fon(1, selected.n_b);
        compute_weighting_factors();
    }

    // Update previous FON for next geometry
    prev_n_a_ = selected.n_a;
    prev_n_b_ = selected.n_b;
}

// ---------------------------------------------------------------------------
// Adaptive FON Solver for REKS(4,4) (FON_SOLVER=4)
// ---------------------------------------------------------------------------
// This solver uses continuation-first approach with early transition detection.
// Key features:
// - Uses FON history prediction when available (smooth PES)
// - Probes for diradical solutions when Hessian indicates transition
// - Enforces symmetry (n_a = n_b) for symmetric molecules
// - No multi-start: probing preserves continuity better

FONPhase REKS::detect_fon_phase(const std::vector<double>& grad, const std::vector<double>& hess) {
    // Detect current FON optimization phase based on Hessian analysis
    //
    // Phase detection logic:
    // - STABLE: det(H) > 0, λ_min > threshold -> well away from transition
    // - PRE_TRANSITION: det(H) > 0 but decreasing OR λ_min < threshold -> approaching transition
    // - TRANSITION: det(H) <= 0 -> at or past saddle point
    //
    // Returns phase enum for adaptive strategy selection.

    const double LAMBDA_MIN_THRESHOLD = 0.1;  // Ha^-1

    // Compute Hessian eigenvalues
    double det_H = hess[0] * hess[3] - hess[1] * hess[2];
    double trace_H = hess[0] + hess[3];
    double discriminant = trace_H * trace_H - 4.0 * det_H;

    double lambda_min = 0.0;
    if (discriminant >= 0.0) {
        lambda_min = 0.5 * (trace_H - std::sqrt(discriminant));
    }

    // Check if det(H) is decreasing (approaching transition)
    bool det_decreasing = (det_H > 0.0 && det_H < det_H_prev_ * 0.9);

    // Update previous determinant for next call
    det_H_prev_ = det_H;

    if (det_H <= 0.0) {
        return FONPhase::TRANSITION;
    } else if (lambda_min < LAMBDA_MIN_THRESHOLD || det_decreasing) {
        return FONPhase::PRE_TRANSITION;
    } else {
        return FONPhase::STABLE;
    }
}

void REKS::adaptive_fon_solver_44() {
    // Adaptive FON Solver for REKS(4,4) (FON_SOLVER=4)
    //
    // Key algorithm:
    // 1. Check FON history for prediction (continuation-first)
    // 2. Detect phase: STABLE, PRE_TRANSITION, or TRANSITION
    // 3. For PRE_TRANSITION: probe diradical solution to catch early transition
    // 4. Use Trust-Region for optimization
    // 5. Enforce symmetry if enabled (n_a = n_b for symmetric molecules)
    // 6. Quality control: reject large jumps, prefer history prediction
    //
    // Goal: Smooth E(R) and FON(R) curves for PES scans

    using namespace reks::constants;

    // Skip if FONs are frozen
    if (fon_oscillation_frozen_) {
        compute_weighting_factors();
        return;
    }

    // Constants for adaptive solver
    const double STEP_PROBE = 0.1;      // Probing step towards diradical
    const double JUMP_ALERT = 0.3;      // Warn if FON changes this much from prediction
    const double JUMP_CRITICAL = 0.5;   // Reject changes larger than this

    // Get current FONs
    double n_a = active_space_->pair(0).fon_p;
    double n_b = active_space_->pair(1).fon_p;

    // =========================================================================
    // FON Initialization
    // =========================================================================
    if (!fon_initialized_44_) {
        // Before localization: keep FONs frozen
        if (iteration_ < 2) {
            active_space_->set_pair_fon(0, 1.0);
            active_space_->set_pair_fon(1, 1.0);
            compute_weighting_factors();
            return;
        }

        // Initialize FON
        bool initialized = false;

        // (1) User-provided initial values
        if (reks_initial_na_ >= 0.0 && reks_initial_nb_ >= 0.0) {
            n_a = reks_initial_na_;
            n_b = reks_initial_nb_;
            active_space_->set_pair_fon(0, n_a);
            active_space_->set_pair_fon(1, n_b);
            initialized = true;
            if (reks_debug_ >= 1) {
                outfile->Printf("  [Adaptive] Using user-provided FON: n_a=%.6f, n_b=%.6f\n", n_a, n_b);
            }
        }

        // (2) Multi-start search: ALWAYS try multiple starting points
        // Each calculation finds the global minimum independently
        if (!initialized) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  [Adaptive] Running multi-start search\n");
            }

            // Start with HAM-DIAG initialization
            initialize_fon_from_hamiltonian_44();
            n_a = active_space_->pair(0).fon_p;
            n_b = active_space_->pair(1).fon_p;

            // Optimize from HAM-DIAG starting point (using Trust-Region)
            active_space_->set_pair_fon(0, n_a);
            active_space_->set_pair_fon(1, n_b);
            reks::FONStartPoint ham_start(n_a, n_b, "HAM_DIAG");
            reks::FONSolution ham_sol = optimize_from_start_44(ham_start, 0);

            // Also try DIRADICAL starting point (n=1.0, 1.0)
            reks::FONStartPoint birad_start(1.0, 1.0, "DIRADICAL");
            reks::FONSolution birad_sol = optimize_from_start_44(birad_start, 1);

            // Also try CLOSED_SHELL starting point (n=2.0, 2.0)
            reks::FONStartPoint closed_start(2.0, 2.0, "CLOSED_SHELL");
            reks::FONSolution closed_sol = optimize_from_start_44(closed_start, 2);

            // Collect converged solutions
            std::vector<reks::FONSolution> solutions;
            if (ham_sol.converged) solutions.push_back(ham_sol);
            if (birad_sol.converged) solutions.push_back(birad_sol);
            if (closed_sol.converged) solutions.push_back(closed_sol);

            // Helper to get name from start_id
            auto get_name = [](int start_id) -> const char* {
                switch(start_id) {
                    case 0: return "HAM_DIAG";
                    case 1: return "DIRADICAL";
                    case 2: return "CLOSED_SHELL";
                    default: return "UNKNOWN";
                }
            };

            if (reks_debug_ >= 1) {
                outfile->Printf("  [Adaptive] Multi-start results:\n");
                for (const auto& sol : solutions) {
                    outfile->Printf("    %s: E=%.8f, n_a=%.4f, n_b=%.4f\n",
                                   get_name(sol.start_id), sol.E_PPS, sol.n_a, sol.n_b);
                }
            }

            // Select the lowest energy solution
            if (!solutions.empty()) {
                auto best = std::min_element(solutions.begin(), solutions.end(),
                                            [](const reks::FONSolution& a, const reks::FONSolution& b) {
                                                return a.E_PPS < b.E_PPS;
                                            });
                n_a = best->n_a;
                n_b = best->n_b;
                if (reks_debug_ >= 1) {
                    outfile->Printf("  [Adaptive] Selected %s: E=%.8f, n_a=%.4f, n_b=%.4f\n",
                                   get_name(best->start_id), best->E_PPS, n_a, n_b);
                }
            } else {
                // Fallback to HAM-DIAG if nothing converged
                initialize_fon_from_hamiltonian_44();
                n_a = active_space_->pair(0).fon_p;
                n_b = active_space_->pair(1).fon_p;
            }

            active_space_->set_pair_fon(0, n_a);
            active_space_->set_pair_fon(1, n_b);
            initialized = true;
        }
        fon_initialized_44_ = true;
    }

    // Lambda to compute E_PPS at given FON values
    auto compute_E_PPS_at = [&](double na, double nb) -> double {
        double old_na = active_space_->pair(0).fon_p;
        double old_nb = active_space_->pair(1).fon_p;
        active_space_->set_pair_fon(0, na);
        active_space_->set_pair_fon(1, nb);
        double E = active_space_->compute_energy_PPS(E_micro_);
        active_space_->set_pair_fon(0, old_na);
        active_space_->set_pair_fon(1, old_nb);
        return E;
    };

    // =========================================================================
    // Trust-Region Optimization with Coupling-Based Direction
    // =========================================================================
    double Delta = trust_radius_44_;
    double E_old = compute_E_PPS_at(n_a, n_b);

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  ================================================================================\n");
        outfile->Printf("  REKS(4,4) FON Optimizer (Adaptive) - SCF Iteration %d\n", iteration_);
        outfile->Printf("  ================================================================================\n");
        outfile->Printf("  Initial: n_a = %.6f, n_b = %.6f, E_PPS = %.10f\n", n_a, n_b, E_old);
    }

    // =========================================================================
    // Compute coupling elements for diradical detection
    // =========================================================================
    // Coupling elements from microstate energies (same as FON_SOLVER=1)
    // dad = coupling between microstates 4,5 and 6,7 for pair (a,d)
    // dbc = coupling between microstates 8,9 and 10,11 for pair (b,c)
    double coupling_dad = (E_micro_[4] + E_micro_[5] - E_micro_[6] - E_micro_[7]) / 2.0;
    double coupling_dbc = (E_micro_[8] + E_micro_[9] - E_micro_[10] - E_micro_[11]) / 2.0;

    const double COUPLING_THRESHOLD = 1e-3;  // Ha - threshold for "strong" coupling
    bool prefer_diradical_a = (std::abs(coupling_dad) > COUPLING_THRESHOLD);
    bool prefer_diradical_b = (std::abs(coupling_dbc) > COUPLING_THRESHOLD);

    if (reks_debug_ >= 1) {
        outfile->Printf("  Coupling: dad=%.6f (|>%.1e? %s), dbc=%.6f (|>%.1e? %s)\n",
                       coupling_dad, COUPLING_THRESHOLD, prefer_diradical_a ? "YES" : "no",
                       coupling_dbc, COUPLING_THRESHOLD, prefer_diradical_b ? "YES" : "no");
    }

    for (int iter = 0; iter < FON_MAX_ITER; ++iter) {
        // 1. Compute gradient and Hessian
        std::vector<double> grad = active_space_->compute_gradient_PPS(E_micro_);
        std::vector<double> hess = active_space_->compute_hessian_PPS(E_micro_);

        // Boundary-aware effective gradient
        // KEY FIX: At upper bound (n≈2), allow gradient if coupling prefers diradical
        double eff_grad_na = grad[0];
        double eff_grad_nb = grad[1];

        // At lower bound (n≈0): zero out if gradient pushes further down
        if (n_a <= 1e-10 && grad[0] > 0.0) eff_grad_na = 0.0;
        if (n_b <= 1e-10 && grad[1] > 0.0) eff_grad_nb = 0.0;

        // At upper bound (n≈2): only constrain if coupling is weak (closed-shell is OK)
        // If coupling is strong, we WANT to move towards diradical
        if (n_a >= 2.0 - 1e-10 && grad[0] < 0.0 && !prefer_diradical_a) {
            eff_grad_na = 0.0;
        }
        if (n_b >= 2.0 - 1e-10 && grad[1] < 0.0 && !prefer_diradical_b) {
            eff_grad_nb = 0.0;
        }

        double eff_grad_norm = std::sqrt(eff_grad_na*eff_grad_na + eff_grad_nb*eff_grad_nb);

        // 2. Check convergence
        if (eff_grad_norm < FON_TOL * 10.0) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  --- FON iter %2d: CONVERGED |g|=%.2e ---\n", iter, eff_grad_norm);
            }
            break;
        }

        // 3. Detect phase and adapt strategy
        FONPhase phase = detect_fon_phase(grad, hess);

        if (reks_debug_ >= 2) {
            const char* phase_str = (phase == FONPhase::STABLE) ? "STABLE" :
                                   (phase == FONPhase::PRE_TRANSITION) ? "PRE_TRANSITION" : "TRANSITION";
            outfile->Printf("  --- FON iter %2d: phase=%s, |g|=%.2e ---\n", iter, phase_str, eff_grad_norm);
        }

        // 4. COUPLING-BASED DIRADICAL PROBING
        // If at closed-shell (n≈2) but coupling is strong, probe diradical direction
        if ((n_a > 1.9 && prefer_diradical_a) || (n_b > 1.9 && prefer_diradical_b)) {
            // Probe if diradical direction lowers energy
            double n_a_probe = prefer_diradical_a && n_a > 1.9 ? std::max(n_a - STEP_PROBE, 1.0) : n_a;
            double n_b_probe = prefer_diradical_b && n_b > 1.9 ? std::max(n_b - STEP_PROBE, 1.0) : n_b;
            double E_probe = compute_E_PPS_at(n_a_probe, n_b_probe);

            if (E_probe < E_old) {
                // Diradical direction is favorable - take the step
                n_a = n_a_probe;
                n_b = n_b_probe;
                E_old = E_probe;
                active_space_->set_pair_fon(0, n_a);
                active_space_->set_pair_fon(1, n_b);
                if (reks_debug_ >= 1) {
                    outfile->Printf("  [Adaptive] Coupling-probe accepted: n_a=%.4f, n_b=%.4f, dE=%.6f\n",
                                   n_a, n_b, E_probe - E_old);
                }
                // Update coupling after step
                coupling_dad = (E_micro_[4] + E_micro_[5] - E_micro_[6] - E_micro_[7]) / 2.0;
                coupling_dbc = (E_micro_[8] + E_micro_[9] - E_micro_[10] - E_micro_[11]) / 2.0;
                prefer_diradical_a = (std::abs(coupling_dad) > COUPLING_THRESHOLD);
                prefer_diradical_b = (std::abs(coupling_dbc) > COUPLING_THRESHOLD);
                continue;  // Skip trust-region step, already took a step
            }
        }

        // 5. PRE_TRANSITION: Additional probe even without strong coupling
        if (phase == FONPhase::PRE_TRANSITION && n_a > 1.9 && n_b > 1.9) {
            // We're in closed-shell but approaching transition
            // Probe if diradical direction lowers energy
            double n_probe = std::max(n_a - STEP_PROBE, 1.0);
            double E_probe = compute_E_PPS_at(n_probe, n_probe);

            if (E_probe < E_old) {
                // Diradical direction is favorable - take the step
                n_a = n_probe;
                n_b = n_probe;
                E_old = E_probe;
                active_space_->set_pair_fon(0, n_a);
                active_space_->set_pair_fon(1, n_b);
                if (reks_debug_ >= 1) {
                    outfile->Printf("  [Adaptive] PRE_TRANSITION: probe accepted, n -> %.4f\n", n_a);
                }
                continue;  // Skip trust-region step, already took a step
            }
        }

        // 5. Compute step: Trust-Region OR Coupling-Based Direction
        double n_a_trial, n_b_trial;
        double det = hess[0] * hess[3] - hess[1] * hess[2];

        // Use coupling-based direction when:
        // 1. Hessian is indefinite (det < 0), OR
        // 2. Coupling is strong AND we're far from diradical optimum
        //    This prevents Trust-Region from converging to wrong minimum
        //
        // Note: threshold 1.28 is chosen so that coupling-based direction continues
        // to push towards diradical until we're close to the correct solution (~1.26)
        bool use_coupling_direction = (det < -1e-8);
        if (!use_coupling_direction) {
            // Check if coupling is strong and we're not yet at diradical
            const double BIRAD_THRESHOLD = 1.28;  // Close to typical diradical FON
            if ((prefer_diradical_a && n_a > BIRAD_THRESHOLD) ||
                (prefer_diradical_b && n_b > BIRAD_THRESHOLD)) {
                use_coupling_direction = true;
            }
        }

        if (use_coupling_direction) {
            // =====================================================================
            // Coupling-Based Direction
            // =====================================================================
            // KEY INSIGHT (from FON_SOLVER=1): When Hessian is indefinite OR coupling
            // is strong, use physics (coupling) to guide direction:
            // - Strong coupling → force towards n=1.0 (diradical)
            // - Weak coupling → follow gradient (may go to closed-shell)
            const double grad_step = 0.05;

            double delta_na = 0.0, delta_nb = 0.0;

            if (prefer_diradical_a) {
                // Strong coupling: force towards n_a = 1.0
                delta_na = (n_a > 1.0) ? -grad_step : (n_a < 1.0 ? grad_step : 0.0);
            } else {
                // Weak coupling: follow gradient
                delta_na = (grad[0] > 0.0) ? -grad_step : (grad[0] < 0.0 ? grad_step : 0.0);
            }

            if (prefer_diradical_b) {
                delta_nb = (n_b > 1.0) ? -grad_step : (n_b < 1.0 ? grad_step : 0.0);
            } else {
                delta_nb = (grad[1] > 0.0) ? -grad_step : (grad[1] < 0.0 ? grad_step : 0.0);
            }

            n_a_trial = std::clamp(n_a + delta_na, 0.0, 2.0);
            n_b_trial = std::clamp(n_b + delta_nb, 0.0, 2.0);

            if (reks_debug_ >= 2) {
                outfile->Printf("  --- Coupling-based direction: det=%.2e, prefer_birad=(%d,%d), delta=(%.3f,%.3f)\n",
                               det, prefer_diradical_a, prefer_diradical_b, delta_na, delta_nb);
            }
        } else {
            // Standard Trust-Region step
            std::vector<std::vector<double>> H = {
                {hess[0], hess[1]},
                {hess[2], hess[3]}
            };
            std::vector<double> p = solve_trust_region_subproblem_2d({eff_grad_na, eff_grad_nb}, H, Delta);

            n_a_trial = std::clamp(n_a + p[0], 0.0, 2.0);
            n_b_trial = std::clamp(n_b + p[1], 0.0, 2.0);
        }

        double E_new = compute_E_PPS_at(n_a_trial, n_b_trial);
        double actual_reduction = E_old - E_new;

        double actual_p0 = n_a_trial - n_a;
        double actual_p1 = n_b_trial - n_b;

        // Accept/reject logic
        if (use_coupling_direction) {
            // Coupling-based direction: always accept if energy doesn't increase
            // (no predicted_reduction calculation needed - we trust coupling)
            if (actual_reduction > -1e-12) {  // Accept if not increasing energy
                n_a = n_a_trial;
                n_b = n_b_trial;
                E_old = E_new;
            }
            // Don't adjust trust radius for coupling-based steps
        } else {
            // Standard Trust-Region with ratio test
            double predicted_reduction = -(eff_grad_na * actual_p0 + eff_grad_nb * actual_p1)
                                        - 0.5 * (hess[0]*actual_p0*actual_p0 + 2*hess[1]*actual_p0*actual_p1 + hess[3]*actual_p1*actual_p1);

            double rho = (std::abs(predicted_reduction) > 1e-15) ? actual_reduction / predicted_reduction : 0.0;

            if (rho > TR_ETA_ACCEPT) {
                n_a = n_a_trial;
                n_b = n_b_trial;
                E_old = E_new;
                if (rho > TR_ETA_VERY_GOOD) {
                    Delta = std::min(2.0 * Delta, TR_DELTA_MAX);
                } else if (rho < TR_ETA_GOOD) {
                    double step_norm = std::sqrt(actual_p0*actual_p0 + actual_p1*actual_p1);
                    Delta = std::max(0.5 * step_norm, TR_DELTA_MIN);
                }
            } else {
                double step_norm = std::sqrt(actual_p0*actual_p0 + actual_p1*actual_p1);
                Delta = std::max(0.25 * step_norm, TR_DELTA_MIN);
            }
        }

        if (Delta < TR_DELTA_MIN * 10.0) break;
    }

    // No post-optimization correction - coupling-based direction is handled within TR loop

    // Update FON values
    active_space_->set_pair_fon(0, n_a);
    active_space_->set_pair_fon(1, n_b);

    // =========================================================================
    // Symmetry Enforcement
    // =========================================================================
    if (reks_enforce_symmetry_) {
        double n_avg = (n_a + n_b) / 2.0;
        n_a = n_avg;
        n_b = n_avg;
        if (reks_debug_ >= 1) {
            outfile->Printf("  [Adaptive] Symmetry enforced: n_a = n_b = %.6f\n", n_avg);
        }
    }

    // =========================================================================
    // Finalize
    // =========================================================================
    trust_radius_44_ = Delta;
    active_space_->set_pair_fon(0, n_a);
    active_space_->set_pair_fon(1, n_b);

    // Orbital swap if needed (keep n in [1, 2] range)
    const double swap_tol = 1.0 - 1e-4;
    if (n_a < swap_tol) {
        int a_idx = Ncore_;
        int d_idx = Ncore_ + 3;
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, a_idx);
            Ca_->set(i, a_idx, Ca_->get(i, d_idx));
            Ca_->set(i, d_idx, tmp);
        }
        n_a = 2.0 - n_a;
        active_space_->set_pair_fon(0, n_a);
    }
    if (n_b < swap_tol) {
        int b_idx = Ncore_ + 1;
        int c_idx = Ncore_ + 2;
        for (int i = 0; i < nso_; ++i) {
            double tmp = Ca_->get(i, b_idx);
            Ca_->set(i, b_idx, Ca_->get(i, c_idx));
            Ca_->set(i, c_idx, tmp);
        }
        n_b = 2.0 - n_b;
        active_space_->set_pair_fon(1, n_b);
    }

    // Update history
    prev_n_a_ = n_a;
    prev_n_b_ = n_b;

    // Recompute weights
    compute_weighting_factors();

    if (reks_debug_ >= 1) {
        double E_PPS = active_space_->compute_energy_PPS(E_micro_);
        outfile->Printf("  --------------------------------------------------------------------------------\n");
        outfile->Printf("  Adaptive FON Summary: n_a = %.6f, n_b = %.6f, E_PPS = %.10f\n", n_a, n_b, E_PPS);
        outfile->Printf("  ================================================================================\n");
    }
}

// ---------------------------------------------------------------------------
// FockMicro2Macro - Coupling Operator Assembly (Filatov 2024, Eq.11-15)
// ---------------------------------------------------------------------------

void REKS::fock_micro_to_macro(int L, double Cl, double** Fa, double** Fb) {
    // Add contribution from microstate L to coupling Fock matrix F_reks_
    // Different orbital blocks have different effective Fock operators.
    //
    // Block structure:
    //   - Core-core, virtual-virtual, core-virtual: F_c = 0.5*C_L*(F_alpha + F_beta)
    //   - Active diagonal: weighted by n^σ_i / f_i
    //   - Core-active coupling: weighted by (1-n^σ_i) / (1-f_i)
    //   - Active-active coupling (i,j): sign(f_i - f_j) weighted
    //
    // Generalized for REKS(N,M) with N/2 active orbitals.

    int n_active = active_space_->n_orbitals();

    if (n_active == 2) {
        fock_micro_to_macro_22(L, Cl, Fa, Fb);
    } else if (n_active == 4) {
        fock_micro_to_macro_44(L, Cl, Fa, Fb);
    } else {
        throw PSIEXCEPTION("fock_micro_to_macro: unsupported number of active orbitals: "
                          + std::to_string(n_active));
    }
}

void REKS::fock_micro_to_macro_22(int L, double Cl, double** Fa, double** Fb) {
    // REKS(2,2) version: 2 active orbitals (r, s)
    //
    // Original specialized implementation for maximum efficiency.

    double** Freks = F_reks_->pointer(0);

    int active_r = active_mo_indices_[0];
    int active_s = active_mo_indices_[1];

    // Get microstate occupations
    int nra = active_space_->alpha_occ(L, 0);
    int nrb = active_space_->beta_occ(L, 0);
    int nsa = active_space_->alpha_occ(L, 1);
    int nsb = active_space_->beta_occ(L, 1);

    // Effective FONs
    double fr = active_space_->get_effective_fon(0);
    double fs = active_space_->get_effective_fon(1);

    double Wc = 0.5 * Cl;

    // Core-core block
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = i; j < Ncore_; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Virtual-virtual block
    for (int i = active_s + 1; i < nsopi_[0]; ++i) {
        for (int j = i; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Core-virtual block
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = active_s + 1; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Active diagonal
    double Wr_a = (fr > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * nra / fr : 0.0;
    double Wr_b = (fr > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * nrb / fr : 0.0;
    Freks[active_r][active_r] += Wr_a * Fa[active_r][active_r] + Wr_b * Fb[active_r][active_r];

    double Ws_a = (fs > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * nsa / fs : 0.0;
    double Ws_b = (fs > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * nsb / fs : 0.0;
    Freks[active_s][active_s] += Ws_a * Fa[active_s][active_s] + Ws_b * Fb[active_s][active_s];

    // Core-active coupling
    double omfr = 1.0 - fr;
    double Wcr_a = (omfr > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * (1 - nra) / omfr : 0.0;
    double Wcr_b = (omfr > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * (1 - nrb) / omfr : 0.0;

    double omfs = 1.0 - fs;
    double Wcs_a = (omfs > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * (1 - nsa) / omfs : 0.0;
    double Wcs_b = (omfs > reks::constants::FON_THRESHOLD) ? 0.5 * Cl * (1 - nsb) / omfs : 0.0;

    for (int c = 0; c < Ncore_; ++c) {
        Freks[c][active_r] += Wcr_a * Fa[c][active_r] + Wcr_b * Fb[c][active_r];
        Freks[c][active_s] += Wcs_a * Fa[c][active_s] + Wcs_b * Fb[c][active_s];
    }

    // Active-virtual coupling
    for (int v = active_s + 1; v < nsopi_[0]; ++v) {
        Freks[active_r][v] += Wr_a * Fa[active_r][v] + Wr_b * Fb[active_r][v];
        Freks[active_s][v] += Ws_a * Fa[active_s][v] + Ws_b * Fb[active_s][v];
    }

    // Active-active coupling (r,s)
    int signrs = (fr > fs) ? 1 : ((fr < fs) ? -1 : 0);
    double Wrs_a = Cl * (nra - nsa) * signrs;
    double Wrs_b = Cl * (nrb - nsb) * signrs;
    Freks[active_r][active_s] += Wrs_a * Fa[active_r][active_s] + Wrs_b * Fb[active_r][active_s];

    // Accumulate Lagrange multiplier for SI
    Wrs_lagr_ += Cl * (nra * Fa[active_r][active_s] + nrb * Fb[active_r][active_s]);
}

void REKS::fock_micro_to_macro_44(int L, double Cl, double** Fa, double** Fb) {
    // REKS(4,4) version: 4 active orbitals (a, b, c, d)
    //
    // Generalized from REKS(2,2) to handle 4 active orbitals.
    // Active orbital indices: a=0, b=1, c=2, d=3

    double** Freks = F_reks_->pointer(0);
    int n_active = 4;

    // Get MO indices for all active orbitals
    std::vector<int> act(n_active);
    std::vector<int> n_alpha(n_active), n_beta(n_active);
    std::vector<double> f_eff(n_active), W_diag_a(n_active), W_diag_b(n_active);

    int last_active = 0;
    for (int i = 0; i < n_active; ++i) {
        act[i] = active_mo_indices_[i];
        n_alpha[i] = active_space_->alpha_occ(L, i);
        n_beta[i] = active_space_->beta_occ(L, i);
        f_eff[i] = active_space_->get_effective_fon(i);

        // Diagonal weight: 0.5 * C_L * n^sigma / f
        W_diag_a[i] = (f_eff[i] > reks::constants::FON_THRESHOLD) ?
                      0.5 * Cl * n_alpha[i] / f_eff[i] : 0.0;
        W_diag_b[i] = (f_eff[i] > reks::constants::FON_THRESHOLD) ?
                      0.5 * Cl * n_beta[i] / f_eff[i] : 0.0;

        if (act[i] > last_active) last_active = act[i];
    }

    double Wc = 0.5 * Cl;  // Core block weight

    // Core-core block (i,j < Ncore)
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = i; j < Ncore_; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Virtual-virtual block (i,j > last_active)
    for (int i = last_active + 1; i < nsopi_[0]; ++i) {
        for (int j = i; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Core-virtual block
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = last_active + 1; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Active diagonal blocks
    for (int i = 0; i < n_active; ++i) {
        int ai = act[i];
        Freks[ai][ai] += W_diag_a[i] * Fa[ai][ai] + W_diag_b[i] * Fb[ai][ai];
    }

    // Core-active coupling for each active orbital
    for (int i = 0; i < n_active; ++i) {
        int ai = act[i];
        double omf = 1.0 - f_eff[i];
        double Wca_a = (omf > reks::constants::FON_THRESHOLD) ?
                       0.5 * Cl * (1 - n_alpha[i]) / omf : 0.0;
        double Wca_b = (omf > reks::constants::FON_THRESHOLD) ?
                       0.5 * Cl * (1 - n_beta[i]) / omf : 0.0;

        for (int c = 0; c < Ncore_; ++c) {
            Freks[c][ai] += Wca_a * Fa[c][ai] + Wca_b * Fb[c][ai];
        }
    }

    // Active-virtual coupling for each active orbital
    for (int i = 0; i < n_active; ++i) {
        int ai = act[i];
        for (int v = last_active + 1; v < nsopi_[0]; ++v) {
            Freks[ai][v] += W_diag_a[i] * Fa[ai][v] + W_diag_b[i] * Fb[ai][v];
        }
    }

    // Active-active coupling for all pairs (i,j) with i < j
    for (int i = 0; i < n_active; ++i) {
        for (int j = i + 1; j < n_active; ++j) {
            int ai = act[i];
            int aj = act[j];
            int sign_ij = (f_eff[i] > f_eff[j]) ? 1 : ((f_eff[i] < f_eff[j]) ? -1 : 0);
            double Wij_a = Cl * (n_alpha[i] - n_alpha[j]) * sign_ij;
            double Wij_b = Cl * (n_beta[i] - n_beta[j]) * sign_ij;
            Freks[ai][aj] += Wij_a * Fa[ai][aj] + Wij_b * Fb[ai][aj];
        }
    }

    // Accumulate Lagrange multipliers for SI (REKS(4,4) has two coupling pairs)
    //
    // For REKS(4,4), there are two primary coupling pairs:
    // - W_ad (stored in Wrs_lagr_): Lagrangian for (a,d) pair
    // - W_bc (stored in Wbc_lagr_): Lagrangian for (b,c) pair
    //
    // The Lagrange multiplier is:
    //   W_pq = sum_L C_L * (n_alpha_p * F_alpha[p][q] + n_beta_p * F_beta[p][q])
    //
    // where p,q are the active orbitals in the coupling pair.

    // (a,d) pair: orbitals 0 and 3
    int a_idx = act[0];  // orbital a
    int d_idx = act[3];  // orbital d
    Wrs_lagr_ += Cl * (n_alpha[0] * Fa[a_idx][d_idx] + n_beta[0] * Fb[a_idx][d_idx]);

    // (b,c) pair: orbitals 1 and 2
    int b_idx = act[1];  // orbital b
    int c_idx = act[2];  // orbital c
    Wbc_lagr_ += Cl * (n_alpha[1] * Fa[b_idx][c_idx] + n_beta[1] * Fb[b_idx][c_idx]);

    // Inter-pair couplings for 9SI:
    // (a,c) pair: orbitals 0 and 2
    Wac_lagr_ += Cl * (n_alpha[0] * Fa[a_idx][c_idx] + n_beta[0] * Fb[a_idx][c_idx]);

    // (b,d) pair: orbitals 1 and 3
    Wbd_lagr_ += Cl * (n_alpha[1] * Fa[b_idx][d_idx] + n_beta[1] * Fb[b_idx][d_idx]);

    if (reks_debug_ >= 3) {
        outfile->Printf("  fock_micro_to_macro_44: L=%d Cl=%.6f\n", L, Cl);
        for (int i = 0; i < n_active; ++i) {
            outfile->Printf("    orb %d: act=%d n_a=%d n_b=%d f=%.6f\n",
                           i, act[i], n_alpha[i], n_beta[i], f_eff[i]);
        }
        outfile->Printf("    W_ad contribution: %.10f\n",
                       Cl * (n_alpha[0] * Fa[a_idx][d_idx] + n_beta[0] * Fb[a_idx][d_idx]));
        outfile->Printf("    W_bc contribution: %.10f\n",
                       Cl * (n_alpha[1] * Fa[b_idx][c_idx] + n_beta[1] * Fb[b_idx][c_idx]));
    }
}

// ---------------------------------------------------------------------------
// SCF Method Overrides
// ---------------------------------------------------------------------------

void REKS::form_D() {
    // 1. Build standard RHF density Da_ (needed for DIIS, convergence checks)
    RHF::form_D();

    // 1.5. Apply orbital localization
    //      Re-localization "anchors" orbitals to physical bonds.
    //      However, for tight geometries (e.g., compressed square H4), repeated
    //      localization can cause SCF oscillations. Use REKS_LOCALIZATION_FREEZE
    //      to stop re-localizing after a certain iteration.
    if (localization_type_ != "NONE") {
        bool should_localize = true;
        if (localization_freeze_iter_ > 0 && iteration_ > localization_freeze_iter_) {
            should_localize = false;
            if (reks_debug_ >= 2 && iteration_ == localization_freeze_iter_ + 1) {
                outfile->Printf("  [REKS] Localization frozen after iteration %d\n",
                               localization_freeze_iter_);
            }
        }
        if (should_localize) {
            localize_active_space();
        }
    }

    // 2. Build REKS-specific base densities from current Ca_
    build_base_densities();

    // 3. Map to microstate densities
    build_microstate_densities();

    // 4. Compute weighting factors (depend on current FON)
    compute_weighting_factors();

    // 5. Monitor active orbitals
    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === Active Orbitals (after form_D) ===\n");

        // Get active MO indices
        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];

        // Overlap <r|s> = C_r^T · S · C_s
        int nso = nsopi_[0];
        double** Cp = Ca_->pointer(0);
        double** Sp = S_->pointer(0);
        double overlap_rs = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                overlap_rs += Cp[mu][active_r] * Sp[mu][nu] * Cp[nu][active_s];
            }
        }
        outfile->Printf("  <r|s> = %12.8f (should be ~0)\n", overlap_rs);

        // Norms <r|r>, <s|s>
        double norm_r = 0.0, norm_s = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                norm_r += Cp[mu][active_r] * Sp[mu][nu] * Cp[nu][active_r];
                norm_s += Cp[mu][active_s] * Sp[mu][nu] * Cp[nu][active_s];
            }
        }
        outfile->Printf("  <r|r> = %12.8f, <s|s> = %12.8f (should be 1.0)\n", norm_r, norm_s);
    }
}

void REKS::form_G() {
    // For SAD guess iteration (iteration_ <= 0): delegate to RHF::form_G()
    // This ensures proper XC initialization (form_V() gets called)
    // Note: sad_ remains true after SAD completes, so we check iteration_ too
    if (sad_ && iteration_ <= 0) {
        RHF::form_G();
        return;
    }

    // REKS completely overrides form_G() - builds microstate Focks + RHF G_ in one JK batch.

    // 1. Build microstate Fock matrices (single JK batch for D00, D10, D01, D11, Da_)
    build_microstate_focks();

    // 2. Compute microstate energies E_L
    compute_microstate_energies();

    // 3. Optimize FON using Newton-Raphson (RexSolver)
    rex_solver();

    // G_ is set inside build_microstate_focks() from Da_ JK results
}

void REKS::form_F() {
    // For SAD guess iteration: use standard RHF Fock matrix
    if (sad_ && iteration_ <= 0) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  form_F: SAD iteration, delegating to RHF\n");
        }
        RHF::form_F();
        return;
    }

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  form_F: REKS iteration %d\n", iteration_);
    }

    // Build REKS coupling Fock matrix F_reks_ from microstate Fock matrices.
    // Procedure:
    //   1. Build F^sigma_L in AO basis
    //   2. Transform to MO: F^sigma_L_MO = C^T * F^sigma_L_AO * C
    //   3. Apply fock_micro_to_macro in MO basis
    //   4. Transform F_reks back to AO for DIIS

    // Zero F_reks_ and Lagrange multipliers
    F_reks_->zero();
    Wrs_lagr_ = 0.0;  // W_ad for REKS(4,4), Wrs for REKS(2,2)
    Wbc_lagr_ = 0.0;  // W_bc for REKS(4,4) only
    Wac_lagr_ = 0.0;  // W_ac for REKS(4,4) 9SI only
    Wbd_lagr_ = 0.0;  // W_bd for REKS(4,4) 9SI only

    int N = nsopi_[0];
    int nmo = Ca_->colspi()[0];
    double** Cp = Ca_->pointer(0);
    double** Tp = Temp_work_->pointer(0);

    // Step 1: Transform microstate Fock matrices from AO to MO basis
    // F_MO = C^T * F_AO * C  (using BLAS DGEMM, reusing pre-allocated matrices)
    // Two-step: Temp = F_AO * C, then F_MO = C^T * Temp
    int n_micro = n_microstates();
    for (int L = 0; L < n_micro; ++L) {
        // Alpha: F_alpha_MO_[L] = C^T * F_alpha_micro_[L] * C
        double** Fa_AO = F_alpha_micro_[L]->pointer(0);
        double** Fa_MO = F_alpha_MO_[L]->pointer(0);
        C_DGEMM('N', 'N', N, N, N, 1.0, Fa_AO[0], N, Cp[0], nmo, 0.0, Tp[0], N);
        C_DGEMM('T', 'N', N, N, N, 1.0, Cp[0], nmo, Tp[0], N, 0.0, Fa_MO[0], N);

        // Beta: F_beta_MO_[L] = C^T * F_beta_micro_[L] * C
        double** Fb_AO = F_beta_micro_[L]->pointer(0);
        double** Fb_MO = F_beta_MO_[L]->pointer(0);
        C_DGEMM('N', 'N', N, N, N, 1.0, Fb_AO[0], N, Cp[0], nmo, 0.0, Tp[0], N);
        C_DGEMM('T', 'N', N, N, N, 1.0, Cp[0], nmo, Tp[0], N, 0.0, Fb_MO[0], N);
    }

    // Step 2: Assemble F_reks_ in MO basis from all 4 microstates
    // Note: L=2 (open-shell) represents L=3 and L=4 in paper - factor 2 already in C_L_[2]
    //       L=3 (triplet-like) represents L=5 and L=6 - factor 2 already in C_L_[3]
    for (int L = 0; L < n_micro; ++L) {
        double** Fa_MO = F_alpha_MO_[L]->pointer(0);
        double** Fb_MO = F_beta_MO_[L]->pointer(0);

        // Scale C_L by symmetry factor for degenerate microstates
        double fact = active_space_->get_symmetry_factor(L);
        double Cl_scaled = fact * C_L_[L];

        fock_micro_to_macro(L, Cl_scaled, Fa_MO, Fb_MO);
    }

    // Step 3: Symmetrize F_reks_ in MO basis (copy upper triangle to lower triangle)
    double** Freks = F_reks_->pointer(0);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Freks[j][i] = Freks[i][j];
        }
    }

    // At this point F_reks_ is in MO basis - save for debug output
    auto F_reks_MO = F_reks_->clone();

    // Store F_reks_MO for MO-basis orbital update in form_C()
    F_reks_MO_ = F_reks_->clone();
    F_reks_MO_->set_name("F_REKS_MO");

    // Step 4: Transform F_reks from MO back to AO basis for DIIS
    // F_AO = S * C * F_MO * C^T * S
    // This ensures [F,P] = F*D*S - S*D*F -> 0 at convergence
    auto SC = linalg::doublet(S_, Ca_, false, false);  // S * C
    auto F_reks_AO = linalg::triplet(SC, F_reks_, SC, false, false, true);  // SC * F_MO * SC^T

    F_reks_->copy(F_reks_AO);

    // Set Fa_ = F_reks_ (now in AO basis)
    Fa_->copy(F_reks_);

    // For RHF: Fb_ = Fa_
    Fb_->copy(Fa_);

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Coupling Fock Matrix ===\n");
        outfile->Printf("  tr(F_reks_AO) = %15.10f\n", F_reks_->trace());
        outfile->Printf("  tr(F_reks_MO) = %15.10f\n", F_reks_MO->trace());
        outfile->Printf("  Wrs_lagr      = %15.10f\n", Wrs_lagr_);

        // Get active MO indices for debug output
        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];

        // Print MO-basis F_reks diagonal
        double** Fmo = F_reks_MO->pointer(0);
        outfile->Printf("  REKS Fock diagonal: ");
        for (int i = 0; i < std::min(8, N); ++i) {
            outfile->Printf("F(%d)=%10.6f ", i, Fmo[i][i]);
        }
        outfile->Printf("\n");
        outfile->Printf("  REKS Fock active: F(r,r)=%10.6f F(s,s)=%10.6f F(r,s)=%10.6f\n",
                       Fmo[active_r][active_r], Fmo[active_s][active_s], Fmo[active_r][active_s]);

        // Check symmetry of MO-basis F_reks
        double max_asym = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double asym = std::abs(Fmo[i][j] - Fmo[j][i]);
                if (asym > max_asym) max_asym = asym;
            }
        }
        outfile->Printf("  Max MO asymmetry: %.3e\n", max_asym);
    }
}

void REKS::form_C(double shift) {
    // For initial guess (iteration_ <= 0): use standard RHF diagonalization
    // This handles both SAD (iteration_ = -1) and HUCKEL/CORE/GWH (iteration_ = 0) guesses
    // where Ca_ is either not yet set or needs to be built from Fa_ via AO-basis diagonalization
    if (iteration_ <= 0) {
        RHF::form_C(shift);
        return;
    }

    // MO-basis orbital update with DIIS support:
    // Convert DIIS-extrapolated Fa_ back to MO basis: F_MO = C^T * Fa_ * C
    // C is S-orthonormal (C^T*S*C = I), so this correctly recovers F_MO

    int N = nsopi_[0];

    // Convert Fa_ (possibly DIIS-extrapolated) back to MO basis: F_MO = C^T * Fa_ * C
    auto F_MO_diis = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);

    // Diagonalize F_MO: F_MO * U = U * eps
    auto U = std::make_shared<Matrix>("Eigenvectors", N, N);
    auto eps = std::make_shared<Vector>("Eigenvalues", N);

    F_MO_diis->diagonalize(U, eps);

    // Fix orbital phase: ensure diagonal elements of U are positive
    double** Up = U->pointer(0);
    for (int j = 0; j < N; ++j) {
        if (Up[j][j] < 0.0) {
            // Flip sign of this column
            for (int i = 0; i < N; ++i) {
                Up[i][j] = -Up[i][j];
            }
        }
    }

    // Update orbitals: C_new = C_old * U
    auto C_old = Ca_->clone();
    auto C_new = linalg::doublet(C_old, U, false, false);
    Ca_->copy(C_new);
    Cb_->copy(Ca_);  // RHF: Cb = Ca

    // Store orbital energies (eigenvalues of F_MO)
    double* eps_p = eps->pointer();
    double* eps_a = epsilon_a_->pointer(0);
    for (int i = 0; i < N; ++i) {
        eps_a[i] = eps_p[i];
    }
    epsilon_b_->copy(*epsilon_a_);

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === MO-basis Orbital Update ===\n");
        // Get active MO indices for debug output
        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];
        outfile->Printf("  Eigenvalues: ");
        for (int i = std::max(0, active_r - 2); i < std::min(N, active_s + 3); ++i) {
            outfile->Printf("eps(%d)=%10.6f ", i, eps_p[i]);
        }
        outfile->Printf("\n");
    }
}

double REKS::compute_E() {
    // For SAD guess iteration: use standard RHF energy
    if (sad_ && iteration_ <= 0) {
        return RHF::compute_E();
    }

    // Compute state-averaged REKS energy: E_SA = sum_L FACT * C_L * E_L
    // Symmetry factor accounts for degenerate microstates (e.g., REKS(2,2): L=2,3 have factor 2)
    double E_SA = 0.0;
    int n_micro = n_microstates();
    for (int L = 0; L < n_micro; ++L) {
        double fact = active_space_->get_symmetry_factor(L);
        E_SA += fact * C_L_[L] * E_micro_[L];
    }

    energies_["Total Energy"] = E_SA;

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS State-Averaged Energy ===\n");
        outfile->Printf("  E_SA = %20.12f\n", E_SA);
    }

    // Print SI-SA state energies (always print after SCF)
    // Note: This runs every iteration, but output is useful for final state
    if (reks_debug_ >= 1) {
        print_SI_energies();
    }

    return E_SA;
}

// ---------------------------------------------------------------------------
// Debug Output Functions
// ---------------------------------------------------------------------------

void REKS::print_microstate_energies() const {
    outfile->Printf("\n  === REKS Microstate Energies ===\n");
    outfile->Printf("  E_1 (L=0, r-up r-down)      = %20.12f\n", E_micro_[0]);
    outfile->Printf("  E_2 (L=1, s-up s-down)      = %20.12f\n", E_micro_[1]);
    outfile->Printf("  E_3 (L=2, r-up s-down)      = %20.12f\n", E_micro_[2]);
    outfile->Printf("  E_5 (L=3, r-up s-up/r-dn s-dn)= %20.12f\n", E_micro_[3]);
}

void REKS::print_fon_info() const {
    double n_r = get_n_r();
    double n_s = get_n_s();
    double f = reks::f_interp(n_r * n_s);
    outfile->Printf("\n  === FON Analysis ===\n");
    outfile->Printf("  n_r = %12.8f, n_s = %12.8f\n", n_r, n_s);
    outfile->Printf("  n_r + n_s = %12.8f (should be 2.0)\n", n_r + n_s);
    outfile->Printf("  f(n_r*n_s) = %12.8f\n", f);
}

void REKS::print_SI_energies() const {
    // Compute SI-SA state energies after SCF convergence
    //
    // REKS(2,2):
    // - 2SI-2SA: Ground state (S0) and first excited singlet (S1) energies
    // - 3SI-2SA: Also includes doubly excited state (S2)
    //
    // REKS(4,4):
    // - 3SI-3SA: S0, S1 (from OSS1), S2 (from OSS2) energies

    // Skip if energies are not computed yet
    if (E_micro_.empty() || std::abs(E_micro_[0]) < reks::constants::ENERGY_THRESHOLD) {
        return;
    }

    // Detect REKS(4,4) by checking number of pairs (>= 2 means REKS(4,4) or higher)
    bool is_reks44 = (active_space_->n_pairs() >= 2);

    outfile->Printf("\n  ============================================================\n");
    if (is_reks44) {
        outfile->Printf("                 3SI-3SA-REKS(4,4) State Energies\n");
    } else {
        outfile->Printf("                    SI-SA-REKS State Energies\n");
    }
    outfile->Printf("  ------------------------------------------------------------\n");

    if (is_reks44) {
        // ================================================================
        // REKS(4,4): Use compute_SI_energies_44 with both Lagrangians
        // ================================================================
        auto si3 = active_space_->compute_SI_energies_44(E_micro_, Wrs_lagr_, Wbc_lagr_, 3);

        // Print diagonal elements of SI Hamiltonian
        outfile->Printf("\n  Diagonal Configuration Energies:\n");
        outfile->Printf("    E_PPS  (Perfectly Paired Singlet) = %20.12f Ha\n", si3.E_PPS);
        outfile->Printf("    E_OSS1 (Open-Shell Singlet 1)     = %20.12f Ha\n", si3.E_OSS);
        outfile->Printf("    E_OSS2 (Open-Shell Singlet 2)     = %20.12f Ha\n", si3.E_DES);

        // Print off-diagonal couplings
        outfile->Printf("\n  Off-diagonal Coupling (Response Lagrangians):\n");
        outfile->Printf("    W_ad (a,d pair Lagrangian)        = %20.12f Ha\n", Wrs_lagr_);
        outfile->Printf("    W_bc (b,c pair Lagrangian)        = %20.12f Ha\n", Wbc_lagr_);
        outfile->Printf("    H_01 (PPS-OSS1)                   = %20.12f Ha\n", si3.H_12);
        outfile->Printf("    H_02 (PPS-OSS2)                   = %20.12f Ha\n", si3.H_23);

        // Print 3SI-3SA eigenvalues
        outfile->Printf("\n  3SI-3SA-REKS(4,4) Adiabatic State Energies:\n");
        outfile->Printf("    %-20s %22s %14s %14s %14s\n", "", "E (Ha)", "C_PPS", "C_OSS1", "C_OSS2");
        outfile->Printf("    %s\n", std::string(84, '-').c_str());
        for (int i = 0; i < 3; ++i) {
            outfile->Printf("    SSR State S%d         %22.12f %14.8f %14.8f %14.8f\n",
                            i, si3.energies[i],
                            si3.coeffs[i * 3 + 0], si3.coeffs[i * 3 + 1], si3.coeffs[i * 3 + 2]);
        }

        // Print excitation energies
        if (si3.energies.size() >= 3) {
            double exc01_eV = (si3.energies[1] - si3.energies[0]) * 27.211386;
            double exc02_eV = (si3.energies[2] - si3.energies[0]) * 27.211386;
            outfile->Printf("\n  Excitation Energies:\n");
            outfile->Printf("    S0 -> S1: Delta E = %12.6f Ha = %12.6f eV\n",
                            si3.energies[1] - si3.energies[0], exc01_eV);
            outfile->Printf("    S0 -> S2: Delta E = %12.6f Ha = %12.6f eV\n",
                            si3.energies[2] - si3.energies[0], exc02_eV);
        }
    } else {
        // ================================================================
        // REKS(2,2): Use original compute_SI_energies with single Wrs
        // ================================================================
        // Compute 2SI-2SA energies
        auto si2 = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, 2);

        // Print diagonal elements of SI Hamiltonian
        outfile->Printf("\n  Diagonal State Energies:\n");
        outfile->Printf("    E_PPS (Perfectly Paired Singlet)  = %20.12f Ha\n", si2.E_PPS);
        outfile->Printf("    E_OSS (Open-Shell Singlet)        = %20.12f Ha\n", si2.E_OSS);
        outfile->Printf("    Triplet                           = %20.12f Ha\n", E_micro_[3]);

        // Print off-diagonal coupling
        outfile->Printf("\n  Off-diagonal Coupling:\n");
        outfile->Printf("    Wrs (Lagrange multiplier)         = %20.12f Ha\n", Wrs_lagr_);
        outfile->Printf("    H_12 = Wrs*(sqrt(n_r)-sqrt(n_s))*sqrt(2) = %20.12f Ha\n", si2.H_12);

        // Print 2SI-2SA eigenvalues
        outfile->Printf("\n  2SI-2SA-REKS State Energies:\n");
        outfile->Printf("    %-20s %22s %14s %14s\n", "", "E (Ha)", "C_PPS", "C_OSS");
        outfile->Printf("    %s\n", std::string(70, '-').c_str());
        for (int i = 0; i < 2; ++i) {
            outfile->Printf("    SSR State S%d         %22.12f %14.8f %14.8f\n",
                            i, si2.energies[i],
                            si2.coeffs[i * 2 + 0], si2.coeffs[i * 2 + 1]);
        }

        // Print excitation energy
        if (si2.energies.size() >= 2) {
            double excitation_eV = (si2.energies[1] - si2.energies[0]) * 27.211386;  // Ha to eV
            outfile->Printf("\n  S0 -> S1 Excitation Energy:\n");
            outfile->Printf("    Delta E = %20.12f Ha = %12.6f eV\n",
                            si2.energies[1] - si2.energies[0], excitation_eV);
        }

        // Compute 3SI-2SA energies
        auto si3 = active_space_->compute_SI_energies(E_micro_, Wrs_lagr_, 3);

        outfile->Printf("\n  3SI-2SA-REKS State Energies:\n");
        outfile->Printf("    E_DES (Doubly Excited Singlet)    = %20.12f Ha\n", si3.E_DES);
        outfile->Printf("    H_23 = Wrs*(sqrt(n_r)+sqrt(n_s))*sqrt(2) = %20.12f Ha\n", si3.H_23);
        outfile->Printf("\n    %-20s %22s %14s %14s %14s\n", "", "E (Ha)", "C_PPS", "C_OSS", "C_DES");
        outfile->Printf("    %s\n", std::string(84, '-').c_str());
        for (int i = 0; i < 3; ++i) {
            outfile->Printf("    SSR State S%d         %22.12f %14.8f %14.8f %14.8f\n",
                            i, si3.energies[i],
                            si3.coeffs[i * 3 + 0], si3.coeffs[i * 3 + 1], si3.coeffs[i * 3 + 2]);
        }
    }

    outfile->Printf("  ============================================================\n\n");
}

// ============================================================================
// multi_scf() Interface Implementation
// ============================================================================

int REKS::n_states() const {
    // REKS needs n_base + 1 J/K matrices:
    // - n_base matrices for base density patterns (0..n_base-1)
    // - 1 matrix for RHF density (used for G_ and DIIS)
    if (!active_space_) {
        return 1;  // Not initialized yet, return default
    }
    return active_space_->n_base_densities() + 1;
}

std::vector<SharedMatrix> REKS::get_orbital_matrices() const {
    // Return all C matrices needed for JK computation:
    // C_base_[0..n_base-1] + C_occ (RHF occupied orbitals)
    //
    // Note: C_base_ matrices are populated in build_microstate_focks()
    // which runs AFTER get_orbital_matrices() is called by multi_scf.
    // So we need to populate them here based on current Ca_.

    if (!active_space_) {
        // Not initialized - return RHF default
        return {Ca_subset("SO", "OCC")};
    }

    int n_base = active_space_->n_base_densities();
    int n_active = active_space_->n_orbitals();
    int nso = nsopi_[0];

    std::vector<SharedMatrix> result;
    result.reserve(n_base + 1);

    double** Cp = Ca_->pointer(0);

    // Populate C matrices for each base density pattern
    for (int p = 0; p < n_base; ++p) {
        int n_occ_in_pattern = __builtin_popcount(p);
        int n_cols = std::max(1, Ncore_ + n_occ_in_pattern);

        auto C_p = std::make_shared<Matrix>("C_base_" + std::to_string(p), nso, n_cols);
        double** C_ptr = C_p->pointer(0);

        int col = 0;

        // Copy core orbitals first
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_; ++k) {
                C_ptr[mu][col + k] = Cp[mu][k];
            }
        }
        col = Ncore_;

        // Add active orbitals based on pattern bits
        for (int i = 0; i < n_active; ++i) {
            if (p & (1 << i)) {
                int mo_idx = active_mo_indices_[i];
                for (int mu = 0; mu < nso; ++mu) {
                    C_ptr[mu][col] = Cp[mu][mo_idx];
                }
                col++;
            }
        }

        // Handle special case: Ncore=0 with pattern 0 (empty density)
        if (Ncore_ == 0 && p == 0) {
            C_p->zero();
        }

        result.push_back(C_p);
    }

    // Add RHF occupied orbitals (last entry)
    result.push_back(Ca_subset("SO", "OCC"));

    return result;
}

void REKS::set_jk_matrices(const std::vector<SharedMatrix>& J_list,
                           const std::vector<SharedMatrix>& K_list,
                           const std::vector<SharedMatrix>& wK_list) {
    // Store precomputed J/K matrices for use in build_microstate_focks()
    //
    // Expected layout (matching get_orbital_matrices()):
    // - J_list[0..n_base-1]: J matrices for base density patterns
    // - J_list[n_base]: J matrix for RHF density
    // Same for K_list.

    if (!active_space_) {
        // Not initialized - call base class
        HF::set_jk_matrices(J_list, K_list, wK_list);
        return;
    }

    int n_base = active_space_->n_base_densities();
    int expected_size = n_base + 1;

    if (static_cast<int>(J_list.size()) != expected_size ||
        static_cast<int>(K_list.size()) != expected_size) {
        throw PSIEXCEPTION("REKS::set_jk_matrices: expected " + std::to_string(expected_size) +
                           " J/K matrices, got " + std::to_string(J_list.size()));
    }

    // Store precomputed matrices
    precomputed_J_base_ = J_list;
    precomputed_K_base_ = K_list;
    use_precomputed_jk_base_ = true;

    // Store wK for base densities (LRC functionals)
    precomputed_wK_base_ = wK_list;

    // Also set base class precomputed matrices for RHF parts
    // (last entry is RHF density)
    precomputed_J_ = {J_list[n_base]};
    precomputed_K_ = {K_list[n_base]};
    if (!wK_list.empty() && n_base < static_cast<int>(wK_list.size())) {
        precomputed_wK_ = {wK_list[n_base]};
    }
    use_precomputed_jk_ = true;

    if (reks_debug_ >= 1) {
        outfile->Printf("  REKS::set_jk_matrices: received %d J/K matrices\n", expected_size);
    }
}

std::shared_ptr<REKS> REKS::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    // Create base wavefunction copy
    auto wfn = Wavefunction::c1_deep_copy(basis);

    // Create REKS wavefunction from the copy
    auto reks_wfn = std::make_shared<REKS>(wfn, functional_, wfn->options(), wfn->psio());

    // Copy matrices that REKS/RHF initializes
    // Include only persistent matrices (some are deleted in finalize())
    if (Ca_) {
        reks_wfn->Ca_ = Ca_subset("AO", "ALL");
        reks_wfn->Cb_ = reks_wfn->Ca_;
    }
    if (Da_) {
        reks_wfn->Da_ = Da_subset("AO");
        reks_wfn->Db_ = reks_wfn->Da_;
    }
    if (Fa_) {
        reks_wfn->Fa_ = Fa_subset("AO");
        reks_wfn->Fb_ = reks_wfn->Fa_;
    }
    if (epsilon_a_) {
        reks_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nalphapi_, "AO", "ALL");
        reks_wfn->epsilon_b_ = reks_wfn->epsilon_a_;
    }

    // H_ and X_ are reset in HF constructor, copy them over
    auto SO2AO = aotoso()->transpose();
    if (H_) reks_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) reks_wfn->X_->remove_symmetry(X_, SO2AO);

    return reks_wfn;
}

// ============================================================================
// Orbital Localization for Active Space
// ============================================================================

void REKS::localize_active_space() {
    if (localization_type_ == "NONE") return;
    // Note: We re-localize every iteration to maintain locality after SCF rotations

    int n_active = active_space_->n_orbitals();
    int nso = Ca_->rowspi()[0];

    // Only print header on first call
    static bool first_call = true;
    if (first_call) {
        outfile->Printf("\n  === REKS Active Space Localization ===\n");
        outfile->Printf("  Method: %s\n", localization_type_.c_str());
        outfile->Printf("  Active orbitals: %d (indices:", n_active);
        for (int i = 0; i < n_active; i++) {
            outfile->Printf(" %d", active_mo_indices_[i]);
        }
        outfile->Printf(")\n");
        first_call = false;
    }

    // Extract active orbitals into separate matrix
    auto C_active = std::make_shared<Matrix>("C_active", nso, n_active);
    double** Cp = Ca_->pointer();
    double** Cap = C_active->pointer();
    for (int i = 0; i < n_active; i++) {
        int mo_idx = active_mo_indices_[i];
        for (int mu = 0; mu < nso; mu++) {
            Cap[mu][i] = Cp[mu][mo_idx];
        }
    }

    // Build localizer (suppress output)
    auto localizer = Localizer::build(localization_type_, basisset_, C_active);
    localizer->set_convergence(1e-8);
    localizer->set_maxiter(50);
    localizer->set_print(0);  // Suppress Boys localizer output

    // Localize
    localizer->localize();

    auto L_active = localizer->L();
    double** Lp = L_active->pointer();

    // For REKS(4,4): Localize occupied (first 2) and virtual (last 2) SEPARATELY
    // Then pair by spatial overlap to get bond orbitals
    if (n_active == 4) {
        // Compute dipole integrals
        auto mints = std::make_shared<MintsHelper>(basisset_);
        std::vector<SharedMatrix> dipole = mints->ao_dipole();
        SharedMatrix dipole_z = dipole[2];  // Z-component

        // Separate localization for occupied and virtual
        // Extract first 2 columns (occupied) and last 2 (virtual)
        auto C_occ = std::make_shared<Matrix>("C_occ", nso, 2);
        auto C_vir = std::make_shared<Matrix>("C_vir", nso, 2);
        double** Cocc = C_occ->pointer();
        double** Cvir = C_vir->pointer();
        double** Ctmp = C_active->pointer();

        for (int mu = 0; mu < nso; mu++) {
            Cocc[mu][0] = Ctmp[mu][0];
            Cocc[mu][1] = Ctmp[mu][1];
            Cvir[mu][0] = Ctmp[mu][2];
            Cvir[mu][1] = Ctmp[mu][3];
        }

        // Localize occupied using user-specified method
        auto loc_occ = Localizer::build(localization_type_, basisset_, C_occ);
        loc_occ->set_convergence(1e-8);
        loc_occ->set_maxiter(50);
        loc_occ->set_print(0);
        loc_occ->localize();
        auto L_occ = loc_occ->L();
        double** Locc = L_occ->pointer();

        // Localize virtual using user-specified method
        auto loc_vir = Localizer::build(localization_type_, basisset_, C_vir);
        loc_vir->set_convergence(1e-8);
        loc_vir->set_maxiter(50);
        loc_vir->set_print(0);
        loc_vir->localize();
        auto L_vir = loc_vir->L();
        double** Lvir = L_vir->pointer();

        // Compute 3D centroids for each localized orbital (x, y, z)
        // This enables proper localization for 2D geometries (square, rhombus)
        SharedMatrix dipole_x = dipole[0];
        SharedMatrix dipole_y = dipole[1];

        std::vector<std::array<double, 3>> occ_cents_3d(2), vir_cents_3d(2);
        for (int i = 0; i < 2; i++) {
            double x_occ = 0.0, y_occ = 0.0, z_occ = 0.0;
            double x_vir = 0.0, y_vir = 0.0, z_vir = 0.0;
            for (int mu = 0; mu < nso; mu++) {
                for (int nu = 0; nu < nso; nu++) {
                    double dx = dipole_x->get(mu, nu);
                    double dy = dipole_y->get(mu, nu);
                    double dz = dipole_z->get(mu, nu);
                    x_occ += Locc[mu][i] * dx * Locc[nu][i];
                    y_occ += Locc[mu][i] * dy * Locc[nu][i];
                    z_occ += Locc[mu][i] * dz * Locc[nu][i];
                    x_vir += Lvir[mu][i] * dx * Lvir[nu][i];
                    y_vir += Lvir[mu][i] * dy * Lvir[nu][i];
                    z_vir += Lvir[mu][i] * dz * Lvir[nu][i];
                }
            }
            occ_cents_3d[i] = {x_occ, y_occ, z_occ};
            vir_cents_3d[i] = {x_vir, y_vir, z_vir};
        }

        // Sort localized orbitals by position to ensure deterministic ordering
        // Use principal axis: find axis with largest spread
        double spread_x = std::abs(occ_cents_3d[0][0] - occ_cents_3d[1][0])
                        + std::abs(vir_cents_3d[0][0] - vir_cents_3d[1][0]);
        double spread_y = std::abs(occ_cents_3d[0][1] - occ_cents_3d[1][1])
                        + std::abs(vir_cents_3d[0][1] - vir_cents_3d[1][1]);
        double spread_z = std::abs(occ_cents_3d[0][2] - occ_cents_3d[1][2])
                        + std::abs(vir_cents_3d[0][2] - vir_cents_3d[1][2]);

        int sort_axis = 0;  // Default: x-axis
        if (spread_y > spread_x && spread_y > spread_z) sort_axis = 1;
        if (spread_z > spread_x && spread_z > spread_y) sort_axis = 2;

        // Sort occupied orbitals: ensure occ[0] has smaller coordinate on sort_axis
        if (occ_cents_3d[0][sort_axis] > occ_cents_3d[1][sort_axis]) {
            // Swap columns in L_occ
            for (int mu = 0; mu < nso; mu++) {
                std::swap(Locc[mu][0], Locc[mu][1]);
            }
            std::swap(occ_cents_3d[0], occ_cents_3d[1]);
        }

        // Sort virtual orbitals: ensure vir[0] has smaller coordinate on sort_axis
        if (vir_cents_3d[0][sort_axis] > vir_cents_3d[1][sort_axis]) {
            // Swap columns in L_vir
            for (int mu = 0; mu < nso; mu++) {
                std::swap(Lvir[mu][0], Lvir[mu][1]);
            }
            std::swap(vir_cents_3d[0], vir_cents_3d[1]);
        }

        if (reks_debug_ >= 2) {
            outfile->Printf("  REKS(4,4) 3D localization centroids:\n");
            outfile->Printf("    Occ 0: (%.4f, %.4f, %.4f)\n",
                            occ_cents_3d[0][0], occ_cents_3d[0][1], occ_cents_3d[0][2]);
            outfile->Printf("    Occ 1: (%.4f, %.4f, %.4f)\n",
                            occ_cents_3d[1][0], occ_cents_3d[1][1], occ_cents_3d[1][2]);
            outfile->Printf("    Vir 0: (%.4f, %.4f, %.4f)\n",
                            vir_cents_3d[0][0], vir_cents_3d[0][1], vir_cents_3d[0][2]);
            outfile->Printf("    Vir 1: (%.4f, %.4f, %.4f)\n",
                            vir_cents_3d[1][0], vir_cents_3d[1][1], vir_cents_3d[1][2]);
        }

        // Match bonding and antibonding orbitals by 3D proximity
        // For each occupied orbital, find the closest virtual orbital
        // Distance: d = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
        double dist_00 = 0.0, dist_01 = 0.0, dist_10 = 0.0, dist_11 = 0.0;
        for (int k = 0; k < 3; k++) {
            dist_00 += std::pow(occ_cents_3d[0][k] - vir_cents_3d[0][k], 2);
            dist_01 += std::pow(occ_cents_3d[0][k] - vir_cents_3d[1][k], 2);
            dist_10 += std::pow(occ_cents_3d[1][k] - vir_cents_3d[0][k], 2);
            dist_11 += std::pow(occ_cents_3d[1][k] - vir_cents_3d[1][k], 2);
        }
        dist_00 = std::sqrt(dist_00);
        dist_01 = std::sqrt(dist_01);
        dist_10 = std::sqrt(dist_10);
        dist_11 = std::sqrt(dist_11);

        if (reks_debug_ >= 2) {
            outfile->Printf("  Distance matrix (occ->vir):\n");
            outfile->Printf("    occ0->vir0: %.4f, occ0->vir1: %.4f\n", dist_00, dist_01);
            outfile->Printf("    occ1->vir0: %.4f, occ1->vir1: %.4f\n", dist_10, dist_11);
        }

        // Optimal pairing: minimize sum of distances
        // Option A: (occ0,vir0) + (occ1,vir1) -> total = dist_00 + dist_11
        // Option B: (occ0,vir1) + (occ1,vir0) -> total = dist_01 + dist_10
        int occ_a, occ_b, vir_c, vir_d;
        if (dist_00 + dist_11 <= dist_01 + dist_10) {
            // Option A: occ0 pairs with vir0, occ1 pairs with vir1
            occ_a = 0; vir_d = 0;  // Pair 0: (a,d)
            occ_b = 1; vir_c = 1;  // Pair 1: (b,c)
            if (reks_debug_ >= 2) {
                outfile->Printf("  Pairing: (occ0,vir0)+(occ1,vir1), total dist=%.4f\n", dist_00 + dist_11);
            }
        } else {
            // Option B: occ0 pairs with vir1, occ1 pairs with vir0
            occ_a = 0; vir_d = 1;  // Pair 0: (a,d)
            occ_b = 1; vir_c = 0;  // Pair 1: (b,c)
            if (reks_debug_ >= 2) {
                outfile->Printf("  Pairing: (occ0,vir1)+(occ1,vir0), total dist=%.4f\n", dist_01 + dist_10);
            }
        }

        // Insert into Ca_
        // REKS(4,4) ordering: a, b, c, d at positions 0, 1, 2, 3
        // GVB pair 0: (a, d) - bonding and antibonding
        // GVB pair 1: (b, c) - bonding and antibonding
        int mo_a = active_mo_indices_[0];
        int mo_b = active_mo_indices_[1];
        int mo_c = active_mo_indices_[2];
        int mo_d = active_mo_indices_[3];

        for (int mu = 0; mu < nso; mu++) {
            Cp[mu][mo_a] = Locc[mu][occ_a];  // a = bonding orbital 0
            Cp[mu][mo_b] = Locc[mu][occ_b];  // b = bonding orbital 1
            Cp[mu][mo_c] = Lvir[mu][vir_c];  // c = antibonding for b
            Cp[mu][mo_d] = Lvir[mu][vir_d];  // d = antibonding for a
        }
    } else {
        // For REKS(2,2) or other: just insert directly
        for (int i = 0; i < n_active; i++) {
            int mo_idx = active_mo_indices_[i];
            for (int mu = 0; mu < nso; mu++) {
                Cp[mu][mo_idx] = Lp[mu][i];
            }
        }
    }

    localization_done_ = true;
}

void REKS::reorder_active_orbitals_for_gvb_pairs() {
    // Compute orbital centroids using ALL THREE dipole components (x, y, z)
    // This enables proper localization for 2D geometries (square, rhombus)
    // where z-centroid alone cannot distinguish orbitals

    auto mints = std::make_shared<MintsHelper>(basisset_);
    std::vector<SharedMatrix> dipole = mints->ao_dipole();
    SharedMatrix dipole_x = dipole[0];  // X-component
    SharedMatrix dipole_y = dipole[1];  // Y-component
    SharedMatrix dipole_z = dipole[2];  // Z-component

    int n_active = 4;
    int nso = Ca_->rowspi()[0];
    double** Cp = Ca_->pointer();

    // Compute 3D centroid (x, y, z) for each active orbital
    std::vector<std::array<double, 3>> centroids_3d(n_active);
    for (int i = 0; i < n_active; i++) {
        int mo_idx = active_mo_indices_[i];
        double x_cent = 0.0, y_cent = 0.0, z_cent = 0.0;
        for (int mu = 0; mu < nso; mu++) {
            for (int nu = 0; nu < nso; nu++) {
                double c_mu = Cp[mu][mo_idx];
                double c_nu = Cp[nu][mo_idx];
                x_cent += c_mu * dipole_x->get(mu, nu) * c_nu;
                y_cent += c_mu * dipole_y->get(mu, nu) * c_nu;
                z_cent += c_mu * dipole_z->get(mu, nu) * c_nu;
            }
        }
        centroids_3d[i] = {x_cent, y_cent, z_cent};
        // Always print centroid info for debugging 3D localization
        outfile->Printf("  Orbital %d (MO %d): centroid = (%.4f, %.4f, %.4f)\n",
                        i, mo_idx, x_cent, y_cent, z_cent);
    }

    // Compute distance matrix between all orbital pairs
    std::vector<std::vector<double>> dist(n_active, std::vector<double>(n_active, 0.0));
    for (int i = 0; i < n_active; i++) {
        for (int j = i + 1; j < n_active; j++) {
            double dx = centroids_3d[j][0] - centroids_3d[i][0];
            double dy = centroids_3d[j][1] - centroids_3d[i][1];
            double dz = centroids_3d[j][2] - centroids_3d[i][2];
            dist[i][j] = dist[j][i] = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    }

    if (reks_debug_ >= 1) {
        outfile->Printf("  Distance matrix:\n");
        for (int i = 0; i < n_active; i++) {
            outfile->Printf("    ");
            for (int j = 0; j < n_active; j++) {
                outfile->Printf("%.3f ", dist[i][j]);
            }
            outfile->Printf("\n");
        }
    }

    // Find GVB pairs based on spatial proximity
    // Strategy: find optimal pairing that minimizes total intra-pair distance
    // For 4 orbitals, there are only 3 possible pairings:
    //   Option A: {0,1} + {2,3}
    //   Option B: {0,2} + {1,3}
    //   Option C: {0,3} + {1,2}
    // Choose the one with minimum sum of intra-pair distances

    // Define all possible pairings
    int pairings[3][4] = {
        {0, 1, 2, 3},  // Option A: (0,1) and (2,3)
        {0, 2, 1, 3},  // Option B: (0,2) and (1,3)
        {0, 3, 1, 2}   // Option C: (0,3) and (1,2)
    };

    double min_total_dist = 1e10;
    int best_pairing = 0;
    for (int p = 0; p < 3; p++) {
        double total = dist[pairings[p][0]][pairings[p][1]] +
                       dist[pairings[p][2]][pairings[p][3]];
        if (total < min_total_dist) {
            min_total_dist = total;
            best_pairing = p;
        }
    }

    int pair1_i = pairings[best_pairing][0];
    int pair1_j = pairings[best_pairing][1];
    int pair2_i = pairings[best_pairing][2];
    int pair2_j = pairings[best_pairing][3];

    // Always print pairing info
    outfile->Printf("  Optimal pairing (total dist=%.4f):\n", min_total_dist);
    outfile->Printf("    GVB pair 0 (a,d): orbitals %d and %d (dist=%.4f)\n",
                    pair1_i, pair1_j, dist[pair1_i][pair1_j]);
    outfile->Printf("    GVB pair 1 (b,c): orbitals %d and %d (dist=%.4f)\n",
                    pair2_i, pair2_j, dist[pair2_i][pair2_j]);

    // REKS(4,4) orbital ordering: a, b, c, d in positions 0, 1, 2, 3
    // GVB pair 0: (a, d) - positions 0 and 3
    // GVB pair 1: (b, c) - positions 1 and 2
    //
    // Assign:
    // a = pair1_i (position 0)
    // d = pair1_j (position 3, pairs with a)
    // b = pair2_i (position 1)
    // c = pair2_j (position 2, pairs with b)

    std::vector<int> new_order = {
        pair1_i,  // a
        pair2_i,  // b
        pair2_j,  // c
        pair1_j   // d
    };

    if (reks_debug_ >= 1) {
        outfile->Printf("  Reordering: old -> new mapping:\n");
        for (int i = 0; i < n_active; i++) {
            int orb = new_order[i];
            outfile->Printf("    position %d: orbital %d (centroid = %.4f, %.4f, %.4f)\n",
                            i, orb, centroids_3d[orb][0], centroids_3d[orb][1], centroids_3d[orb][2]);
        }
    }

    // Check if reordering is needed
    bool needs_reorder = false;
    for (int i = 0; i < n_active; i++) {
        if (new_order[i] != i) {
            needs_reorder = true;
            break;
        }
    }

    if (!needs_reorder) {
        outfile->Printf("  Orbitals already in correct order.\n");
        return;
    }

    // Create temporary storage for columns
    std::vector<std::vector<double>> temp_cols(n_active, std::vector<double>(nso));
    for (int i = 0; i < n_active; i++) {
        int mo_idx = active_mo_indices_[new_order[i]];
        for (int mu = 0; mu < nso; mu++) {
            temp_cols[i][mu] = Cp[mu][mo_idx];
        }
    }

    // Write back in new order
    for (int i = 0; i < n_active; i++) {
        int mo_idx = active_mo_indices_[i];
        for (int mu = 0; mu < nso; mu++) {
            Cp[mu][mo_idx] = temp_cols[i][mu];
        }
    }

    outfile->Printf("  Orbitals reordered for GVB pairs: (a,d)=(0,3), (b,c)=(1,2)\n");
}

// ---------------------------------------------------------------------------
// Multi-Start FON Optimization for REKS(4,4)
// ---------------------------------------------------------------------------

void REKS::rex_solver_44_multi_start() {
    // Multi-start FON optimization:
    // 1. Generate starting points (diradical, closed-shell, continuation, HAM-DIAG)
    // 2. Optimize from each start using Trust-Region
    // 3. Select best solution based on criterion (ENERGY/DIABATIC/AUTO)

    using namespace reks;

    // Skip if FONs are frozen
    if (fon_oscillation_frozen_) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === Multi-Start REKS(4,4): FON frozen, skipping ===\n");
        }
        compute_weighting_factors();
        return;
    }

    // Wait for localization (same as trust_region_fon_optimizer_44)
    if (!fon_initialized_44_ && iteration_ < 2) {
        if (reks_debug_ >= 1) {
            outfile->Printf("  [MS] Waiting for localization (iteration %d < 2)\n", iteration_);
        }
        active_space_->set_pair_fon(0, 1.0);
        active_space_->set_pair_fon(1, 1.0);
        compute_weighting_factors();
        return;
    }

    // Get HAM-DIAG suggestion
    initialize_fon_from_hamiltonian_44();
    double ham_na = active_space_->pair(0).fon_p;
    double ham_nb = active_space_->pair(1).fon_p;

    // Generate starting points
    MultiStartGenerator generator;
    auto starts = generator.generate(prev_n_a_, prev_n_b_, ham_na, ham_nb);

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === Multi-Start REKS(4,4): FON Optimization ===\n");
        outfile->Printf("  Starting points: %zu\n", starts.size());
        for (size_t i = 0; i < starts.size(); ++i) {
            outfile->Printf("    [%zu] %s: (%.4f, %.4f)\n",
                           i, starts[i].name.c_str(), starts[i].n_a, starts[i].n_b);
        }
    }

    // Optimize from each starting point
    std::vector<FONSolution> solutions;
    for (size_t i = 0; i < starts.size(); ++i) {
        FONSolution sol = optimize_from_start_44(starts[i], static_cast<int>(i));
        if (sol.converged) {
            sol.classify();
            solutions.push_back(sol);
        }
    }

    if (solutions.empty()) {
        // No solutions converged - fall back to HAM-DIAG
        outfile->Printf("  [MS] Warning: No solutions converged, using HAM-DIAG\n");
        active_space_->set_pair_fon(0, ham_na);
        active_space_->set_pair_fon(1, ham_nb);
        compute_weighting_factors();
        fon_initialized_44_ = true;
        return;
    }

    // Remove duplicates
    std::vector<FONSolution> unique;
    for (const auto& sol : solutions) {
        bool is_dup = false;
        for (const auto& u : unique) {
            if (sol.is_duplicate_of(u)) {
                is_dup = true;
                break;
            }
        }
        if (!is_dup) {
            unique.push_back(sol);
        }
    }

    if (reks_debug_ >= 1 || reks_report_all_solutions_) {
        outfile->Printf("  Converged solutions: %zu total, %zu unique\n",
                       solutions.size(), unique.size());
        for (size_t i = 0; i < unique.size(); ++i) {
            const char* char_str = unique[i].character == FONCharacter::DIRADICAL ? "DIRADICAL" :
                                   unique[i].character == FONCharacter::CLOSED_SHELL ? "CLOSED_SHELL" :
                                   "MIXED";
            outfile->Printf("    [%zu] E=%.10f, n=(%.4f, %.4f), %s\n",
                           i, unique[i].E_PPS, unique[i].n_a, unique[i].n_b, char_str);
        }
    }

    // Select best solution
    FONSolution selected = select_solution(unique, prev_n_a_, prev_n_b_,
                                           reks_branch_criterion_, reks_energy_tolerance_);

    // Apply selected solution
    active_space_->set_pair_fon(0, selected.n_a);
    active_space_->set_pair_fon(1, selected.n_b);
    fon_initialized_44_ = true;

    compute_weighting_factors();

    if (reks_debug_ >= 1) {
        const char* char_str = selected.character == FONCharacter::DIRADICAL ? "DIRADICAL" :
                               selected.character == FONCharacter::CLOSED_SHELL ? "CLOSED_SHELL" :
                               "MIXED";
        outfile->Printf("  Selected: E=%.10f, n=(%.4f, %.4f), %s [%s criterion]\n",
                       selected.E_PPS, selected.n_a, selected.n_b,
                       char_str, reks_branch_criterion_.c_str());
    }
}

reks::FONSolution REKS::optimize_from_start_44(const reks::FONStartPoint& start, int start_id) {
    // Optimize FON from a single starting point using Trust-Region
    // Returns converged FONSolution

    using namespace reks::constants;

    reks::FONSolution result;
    result.n_a = start.n_a;
    result.n_b = start.n_b;
    result.start_id = start_id;
    result.converged = false;

    // Set starting FON
    active_space_->set_pair_fon(0, start.n_a);
    active_space_->set_pair_fon(1, start.n_b);

    double n_a = start.n_a;
    double n_b = start.n_b;

    // Lambda to compute E_PPS at given FON values
    auto compute_E_PPS_at = [&](double na, double nb) -> double {
        double old_na = active_space_->pair(0).fon_p;
        double old_nb = active_space_->pair(1).fon_p;
        active_space_->set_pair_fon(0, na);
        active_space_->set_pair_fon(1, nb);
        double E = active_space_->compute_energy_PPS(E_micro_);
        active_space_->set_pair_fon(0, old_na);
        active_space_->set_pair_fon(1, old_nb);
        return E;
    };

    // Trust-Region optimization (simplified, fewer iterations)
    double Delta = TR_DELTA_INIT;
    double E_old = compute_E_PPS_at(n_a, n_b);

    for (int iter = 0; iter < 50; ++iter) {  // Reduced iterations for multi-start
        // Compute gradient
        active_space_->set_pair_fon(0, n_a);
        active_space_->set_pair_fon(1, n_b);
        std::vector<double> grad = active_space_->compute_gradient_PPS(E_micro_);

        // Boundary-aware effective gradient
        double eff_grad_na = grad[0];
        double eff_grad_nb = grad[1];
        if (n_a >= 2.0 - 1e-10 && grad[0] < 0.0) eff_grad_na = 0.0;
        if (n_b >= 2.0 - 1e-10 && grad[1] < 0.0) eff_grad_nb = 0.0;
        if (n_a <= 1e-10 && grad[0] > 0.0) eff_grad_na = 0.0;
        if (n_b <= 1e-10 && grad[1] > 0.0) eff_grad_nb = 0.0;

        double grad_norm = std::sqrt(eff_grad_na * eff_grad_na + eff_grad_nb * eff_grad_nb);

        // Check convergence
        if (grad_norm < FON_TOL) {
            result.n_a = n_a;
            result.n_b = n_b;
            result.E_PPS = compute_E_PPS_at(n_a, n_b);
            result.converged = true;
            break;
        }

        // Compute Hessian (flat -> 2x2)
        std::vector<double> hess_flat = active_space_->compute_hessian_PPS(E_micro_);
        std::vector<std::vector<double>> Hess = {
            {hess_flat[0], hess_flat[1]},
            {hess_flat[2], hess_flat[3]}
        };

        // Solve Trust-Region subproblem
        std::vector<double> eff_grad = {eff_grad_na, eff_grad_nb};
        std::vector<double> step = solve_trust_region_subproblem_2d(eff_grad, Hess, Delta);

        // Trial point with clamping
        double n_a_trial = std::max(0.0, std::min(2.0, n_a + step[0]));
        double n_b_trial = std::max(0.0, std::min(2.0, n_b + step[1]));

        // Compute actual vs predicted reduction
        double E_new = compute_E_PPS_at(n_a_trial, n_b_trial);
        double actual = E_old - E_new;
        double predicted = -(eff_grad_na * step[0] + eff_grad_nb * step[1]) -
                          0.5 * (step[0] * (Hess[0][0] * step[0] + Hess[0][1] * step[1]) +
                                 step[1] * (Hess[1][0] * step[0] + Hess[1][1] * step[1]));

        if (std::abs(predicted) < 1e-15) predicted = 1e-15;
        double rho = actual / predicted;

        // Accept/reject and update trust radius
        double step_norm = std::sqrt(step[0] * step[0] + step[1] * step[1]);

        if (rho > TR_ETA_VERY_GOOD) {
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
            Delta = std::min(2.0 * Delta, TR_DELTA_MAX);
        } else if (rho > TR_ETA_GOOD) {
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
        } else if (rho > TR_ETA_ACCEPT) {
            n_a = n_a_trial;
            n_b = n_b_trial;
            E_old = E_new;
            Delta = std::max(0.5 * step_norm, TR_DELTA_MIN);
        } else {
            Delta = std::max(0.25 * step_norm, TR_DELTA_MIN);
        }

        // Check for trust radius collapse
        if (Delta < TR_DELTA_MIN * 10.0) {
            result.n_a = n_a;
            result.n_b = n_b;
            result.E_PPS = E_old;
            result.converged = true;
            break;
        }
    }

    // Final result
    if (!result.converged) {
        result.n_a = n_a;
        result.n_b = n_b;
        result.E_PPS = E_old;
        result.converged = true;  // Accept even if not fully converged
    }

    return result;
}

reks::FONSolution REKS::select_solution(
    const std::vector<reks::FONSolution>& solutions,
    double prev_na, double prev_nb,
    const std::string& criterion,
    double energy_tol)
{
    // Select best solution based on criterion

    if (solutions.empty()) {
        reks::FONSolution empty;
        return empty;
    }

    if (solutions.size() == 1) {
        return solutions[0];
    }

    // Find lowest energy solution
    const reks::FONSolution* lowest = &solutions[0];
    for (const auto& sol : solutions) {
        if (sol.E_PPS < lowest->E_PPS) {
            lowest = &sol;
        }
    }

    // Find diabatic solution (closest to previous)
    const reks::FONSolution* diabatic = &solutions[0];
    if (prev_na >= 0.0 && prev_nb >= 0.0) {
        double min_dist = solutions[0].diabatic_distance(prev_na, prev_nb);
        for (const auto& sol : solutions) {
            double dist = sol.diabatic_distance(prev_na, prev_nb);
            if (dist < min_dist) {
                min_dist = dist;
                diabatic = &sol;
            }
        }
    }

    // Selection based on criterion
    if (criterion == "ENERGY") {
        return *lowest;
    } else if (criterion == "DIABATIC") {
        return *diabatic;
    } else {  // AUTO
        // Use diabatic if within energy tolerance of lowest
        double gap = diabatic->E_PPS - lowest->E_PPS;
        if (gap <= energy_tol) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  [AUTO] Diabatic within tolerance (gap=%.6f Ha ≤ %.6f Ha)\n",
                               gap, energy_tol);
            }
            return *diabatic;
        } else {
            if (reks_debug_ >= 1) {
                outfile->Printf("  [AUTO] Using lowest energy (gap=%.6f Ha > %.6f Ha)\n",
                               gap, energy_tol);
            }
            return *lowest;
        }
    }
}

}  // namespace scf
}  // namespace psi
