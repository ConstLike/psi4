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
        rex_solver_44();
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

void REKS::rex_solver_44() {
    // REKS(4,4) FON optimization: two variables (n_a, n_b)
    //
    // TWO-PHASE SCF FOR STABLE CONVERGENCE (same as REKS(2,2)):
    //
    // Phase 1: Freeze FON for initial iterations
    // Phase 2: Enable FON optimization with step limiting
    //
    // Uses 2D Newton-Raphson with gradient [dE/dn_a, dE/dn_b]
    // and 2x2 Hessian [d²E/dn_a², d²E/dn_a dn_b; d²E/dn_b dn_a, d²E/dn_b²]
    //
    // Constraints: n_d = 2 - n_a, n_c = 2 - n_b

    // Get current FONs
    double n_a = active_space_->pair(0).fon_p;
    double n_b = active_space_->pair(1).fon_p;

    // =========================================================================
    // Phase 1: Freeze FON for initial iterations
    // =========================================================================
    if (iteration_ < fon_freeze_iters_) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === RexSolver(4,4): Phase 1 (FON frozen) ===\n");
            outfile->Printf("  Iteration %d < %d: keeping n_a=%.6f, n_b=%.6f (frozen)\n",
                           iteration_, fon_freeze_iters_, n_a, n_b);
        }
        // Update previous FONs for Phase 2 tracking
        prev_n_a_ = n_a;
        prev_n_b_ = n_b;
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
    double max_step = reks::constants::FON_MAX_STEP;

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver(4,4): Phase 2 (2D FON Optimization) ===\n");
        outfile->Printf("  Initial FONs: n_a=%.6f, n_b=%.6f\n", n_a, n_b);
        outfile->Printf("  prev_n_a=%.6f, prev_n_b=%.6f, fon_max_delta=%.4f\n",
                       prev_n_a_, prev_n_b_, fon_max_delta_);
        outfile->Printf("  E_micro: ");
        for (int L = 0; L < 12 && L < static_cast<int>(E_micro_.size()); ++L) {
            outfile->Printf("[%d]=%.8f ", L, E_micro_[L]);
        }
        outfile->Printf("\n");
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

    // Oscillation detection: track last two FON values
    double prev_n_a = n_a, prev_n_b = n_b;
    double prev2_n_a = -1.0, prev2_n_b = -1.0;
    int oscillation_count = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute gradient [dE/dn_a, dE/dn_b]
        std::vector<double> grad = active_space_->compute_energy_gradient(E_micro_, C_L_);

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
                double E_SA = compute_E_SA();
                outfile->Printf("    iter %2d: n_a=%.6f n_b=%.6f E=%.10f |eff_grad|=%.2e CONVERGED (boundary)\n",
                               iter, n_a, n_b, E_SA, eff_grad_norm);
            }
            break;
        }

        // Compute 2x2 Hessian [H_aa, H_ab; H_ba, H_bb]
        std::vector<double> hess = active_space_->compute_energy_hessian(E_micro_, C_L_);

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
            if (reks_debug_ >= 2) {
                outfile->Printf("    [1D Newton in n_b: H_bb=%.4e, grad_b=%.4e, constrained_a=%d, singular_a=%d]\n",
                               H_bb, grad[1], constrained_a, singular_a);
            }
        } else if ((constrained_b || singular_b) && !singular_a && std::abs(H_aa) > reks::constants::HESSIAN_THRESHOLD) {
            // n_b direction blocked or singular: do 1D Newton-Raphson in n_a direction only
            delta_na = -grad[0] / H_aa;  // 1D Newton step
            delta_nb = 0.0;
            if (reks_debug_ >= 2) {
                outfile->Printf("    [1D Newton in n_a: H_aa=%.4e, grad_a=%.4e, constrained_b=%d, singular_b=%d]\n",
                               H_aa, grad[0], constrained_b, singular_b);
            }
        } else if (std::abs(det) < reks::constants::HESSIAN_THRESHOLD || det < 0.0) {
            // Both at boundary or indefinite Hessian: use gradient descent
            // But use adaptive step based on gradient magnitude
            double grad_norm = std::sqrt(grad[0]*grad[0] + grad[1]*grad[1]);
            double step_scale = (grad_norm > 0.1) ? max_step : max_step * 0.5;

            // Scale by inverse gradient to approximate Newton behavior
            delta_na = (std::abs(eff_grad_na) > 1e-10) ? -step_scale * std::copysign(1.0, eff_grad_na) : 0.0;
            delta_nb = (std::abs(eff_grad_nb) > 1e-10) ? -step_scale * std::copysign(1.0, eff_grad_nb) : 0.0;

            if (reks_debug_ >= 2) {
                outfile->Printf("    [Gradient descent: det=%.4e, grad_norm=%.4e]\n", det, grad_norm);
            }
        } else {
            // Newton-Raphson: delta = -H^{-1} * grad
            // For 2x2: [a b; b d]^{-1} = (1/det) * [d -b; -b a]
            delta_na = -(H_bb * grad[0] - H_ab * grad[1]) / det;
            delta_nb = -(-H_ab * grad[0] + H_aa * grad[1]) / det;
        }

        // Damping: limit step size
        if (std::abs(delta_na) > max_step) {
            delta_na = std::copysign(max_step, delta_na);
        }
        if (std::abs(delta_nb) > max_step) {
            delta_nb = std::copysign(max_step, delta_nb);
        }

        // Update FONs with bounds [0.0, 2.0]
        double n_a_new = std::clamp(n_a + delta_na, 0.0, 2.0);
        double n_b_new = std::clamp(n_b + delta_nb, 0.0, 2.0);

        // Actual deltas after bounds enforcement
        delta_na = n_a_new - n_a;
        delta_nb = n_b_new - n_b;

        // Check for oscillation ONLY when using gradient descent (not Newton)
        // Newton methods naturally converge quadratically and shouldn't trigger oscillation detection
        // We identify gradient descent when det < 0 and both directions are being updated
        bool using_gradient_descent = (std::abs(det) < reks::constants::HESSIAN_THRESHOLD || det < 0.0) &&
                                       std::abs(delta_na) > tol && std::abs(delta_nb) > tol;

        if (using_gradient_descent && prev2_n_a >= 0.0) {
            double diff_a = std::abs(n_a_new - prev2_n_a);
            double diff_b = std::abs(n_b_new - prev2_n_b);
            // Very tight threshold: only trigger when bouncing in small region
            if (diff_a < 0.02 && diff_b < 0.02) {
                oscillation_count++;
                if (oscillation_count >= 3) {  // Require 3 oscillations before intervening
                    // Oscillation detected! Take midpoint of oscillating values
                    n_a_new = (n_a_new + prev_n_a) / 2.0;
                    n_b_new = (n_b_new + prev_n_b) / 2.0;
                    if (reks_debug_ >= 1) {
                        outfile->Printf("    Oscillation detected (count=%d)! Using midpoint: n_a=%.6f n_b=%.6f\n",
                                       oscillation_count, n_a_new, n_b_new);
                    }
                    // Also reduce step size for future iterations
                    max_step *= 0.5;
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
            double E_SA = compute_E_SA();
            outfile->Printf("    iter %2d: n_a=%.6f n_b=%.6f E=%.10f |grad|=%.2e |eff_grad|=%.2e det=%.2e\n",
                           iter, n_a, n_b, E_SA,
                           std::sqrt(grad[0]*grad[0] + grad[1]*grad[1]), eff_grad_norm, det);
        }

        // Convergence check on actual step
        if (std::abs(delta_na) < tol && std::abs(delta_nb) < tol) {
            break;
        }
    }

    // =========================================================================
    // Apply FON change limiting for DIIS stability
    // =========================================================================
    // Limit the change from previous SCF iteration to fon_max_delta_
    // This ensures Fock matrices remain similar between DIIS iterations

    double delta_a_from_prev = n_a - prev_n_a_;
    double delta_b_from_prev = n_b - prev_n_b_;

    bool limited = false;
    if (std::abs(delta_a_from_prev) > fon_max_delta_) {
        n_a = prev_n_a_ + std::copysign(fon_max_delta_, delta_a_from_prev);
        n_a = std::clamp(n_a, 0.0, 2.0);
        limited = true;
    }
    if (std::abs(delta_b_from_prev) > fon_max_delta_) {
        n_b = prev_n_b_ + std::copysign(fon_max_delta_, delta_b_from_prev);
        n_b = std::clamp(n_b, 0.0, 2.0);
        limited = true;
    }

    if (limited) {
        // Update active space with limited FONs
        active_space_->set_pair_fon(0, n_a);
        active_space_->set_pair_fon(1, n_b);
        compute_weighting_factors();

        if (reks_debug_ >= 1) {
            outfile->Printf("  FON change limited: n_a: %.6f->%.6f (delta %.4f), n_b: %.6f->%.6f (delta %.4f)\n",
                           prev_n_a_, n_a, n_a - prev_n_a_, prev_n_b_, n_b, n_b - prev_n_b_);
        }
    }

    // Store current FONs for next iteration's delta limiting
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

    // 1.5. Apply orbital localization on every REKS iteration (after first)
    //      This keeps active orbitals localized despite SCF rotations
    if (iteration_ >= 1 && localization_type_ != "NONE") {
        localize_active_space();
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

        // Localize occupied
        auto loc_occ = Localizer::build("BOYS", basisset_, C_occ);
        loc_occ->set_convergence(1e-8);
        loc_occ->set_maxiter(50);
        loc_occ->set_print(0);
        loc_occ->localize();
        auto L_occ = loc_occ->L();
        double** Locc = L_occ->pointer();

        // Localize virtual
        auto loc_vir = Localizer::build("BOYS", basisset_, C_vir);
        loc_vir->set_convergence(1e-8);
        loc_vir->set_maxiter(50);
        loc_vir->set_print(0);
        loc_vir->localize();
        auto L_vir = loc_vir->L();
        double** Lvir = L_vir->pointer();

        // Compute z-centroid for each localized orbital
        std::vector<std::pair<double, int>> occ_cents(2), vir_cents(2);
        for (int i = 0; i < 2; i++) {
            double z_occ = 0.0, z_vir = 0.0;
            for (int mu = 0; mu < nso; mu++) {
                for (int nu = 0; nu < nso; nu++) {
                    double dz = dipole_z->get(mu, nu);
                    z_occ += Locc[mu][i] * dz * Locc[nu][i];
                    z_vir += Lvir[mu][i] * dz * Lvir[nu][i];
                }
            }
            occ_cents[i] = {z_occ, i};
            vir_cents[i] = {z_vir, i};
        }

        // Sort by z-coordinate
        std::sort(occ_cents.begin(), occ_cents.end());
        std::sort(vir_cents.begin(), vir_cents.end());

        if (reks_debug_ >= 1) {
            outfile->Printf("  REKS(4,4) localization centroids:\n");
            outfile->Printf("    Occ left  z=%.4f (idx %d) -> MO a\n", occ_cents[0].first, occ_cents[0].second);
            outfile->Printf("    Occ right z=%.4f (idx %d) -> MO b\n", occ_cents[1].first, occ_cents[1].second);
            outfile->Printf("    Vir left  z=%.4f (idx %d) -> MO d\n", vir_cents[0].first, vir_cents[0].second);
            outfile->Printf("    Vir right z=%.4f (idx %d) -> MO c\n", vir_cents[1].first, vir_cents[1].second);
        }

        // Now we have:
        // occ_cents[0] = left bond (σ12 at z ~ H1-H2 midpoint)
        // occ_cents[1] = right bond (σ34 at z ~ H3-H4 midpoint)
        // vir_cents[0] = left antibond (σ*12)
        // vir_cents[1] = right antibond (σ*34)
        //
        // REKS(4,4) ordering: a, b, c, d
        // Pair 0 (a,d): left bond (σ12) and left antibond (σ*12)
        // Pair 1 (b,c): right bond (σ34) and right antibond (σ*34)
        //
        // Assignment:
        // MO 0 (a) = occ left (σ12)
        // MO 1 (b) = occ right (σ34)
        // MO 2 (c) = vir right (σ*34) - pairs with b
        // MO 3 (d) = vir left (σ*12) - pairs with a

        // Insert into Ca_
        int mo_a = active_mo_indices_[0];
        int mo_b = active_mo_indices_[1];
        int mo_c = active_mo_indices_[2];
        int mo_d = active_mo_indices_[3];

        for (int mu = 0; mu < nso; mu++) {
            Cp[mu][mo_a] = Locc[mu][occ_cents[0].second];  // a = left bond
            Cp[mu][mo_b] = Locc[mu][occ_cents[1].second];  // b = right bond
            Cp[mu][mo_c] = Lvir[mu][vir_cents[1].second];  // c = right antibond
            Cp[mu][mo_d] = Lvir[mu][vir_cents[0].second];  // d = left antibond
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
    // Compute orbital centroids using dipole integrals
    // For linear H4 along z-axis, sort by z-centroid

    auto mints = std::make_shared<MintsHelper>(basisset_);
    std::vector<SharedMatrix> dipole = mints->ao_dipole();
    SharedMatrix dipole_z = dipole[2];  // Z-component

    int n_active = 4;
    int nso = Ca_->rowspi()[0];
    double** Cp = Ca_->pointer();

    // Compute z-centroid for each active orbital
    std::vector<std::pair<double, int>> centroids(n_active);
    for (int i = 0; i < n_active; i++) {
        int mo_idx = active_mo_indices_[i];
        double z_cent = 0.0;
        for (int mu = 0; mu < nso; mu++) {
            for (int nu = 0; nu < nso; nu++) {
                z_cent += Cp[mu][mo_idx] * dipole_z->get(mu, nu) * Cp[nu][mo_idx];
            }
        }
        centroids[i] = {z_cent, i};
        if (reks_debug_ >= 1) {
            outfile->Printf("  Orbital %d (MO %d): z-centroid = %.4f\n",
                            i, mo_idx, z_cent);
        }
    }

    // Sort by z-coordinate
    std::sort(centroids.begin(), centroids.end());

    // REKS(4,4) orbital ordering: a, b, c, d in positions 0, 1, 2, 3
    // GVB pair 0: (a, d) - positions 0 and 3
    // GVB pair 1: (b, c) - positions 1 and 2
    //
    // For linear H4 with Boys localization, we get 4 atom-centered orbitals.
    // After sorting by z-centroid:
    // centroids[0] = H1 (leftmost)
    // centroids[1] = H2
    // centroids[2] = H3
    // centroids[3] = H4 (rightmost)
    //
    // For symmetric GVB pairs, we need:
    // - GVB pair 0: H1 + H2 (adjacent atoms, left bond)
    // - GVB pair 1: H3 + H4 (adjacent atoms, right bond)
    //
    // Correct assignment:
    // a = H1 (centroids[0]) - position 0
    // d = H2 (centroids[1]) - position 3 (pair with a)
    // b = H3 (centroids[2]) - position 1
    // c = H4 (centroids[3]) - position 2 (pair with b)
    //
    // Final order: [H1, H3, H4, H2] for positions [a, b, c, d]

    std::vector<int> new_order = {
        centroids[0].second,  // a = H1 (leftmost)
        centroids[2].second,  // b = H3 (right-left)
        centroids[3].second,  // c = H4 (rightmost)
        centroids[1].second   // d = H2 (left-right, pairs with a)
    };

    if (reks_debug_ >= 1) {
        outfile->Printf("  Reordering: old -> new mapping:\n");
        for (int i = 0; i < n_active; i++) {
            outfile->Printf("    position %d: orbital %d (z=%.4f)\n",
                            i, new_order[i], centroids[i].first);
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

}  // namespace scf
}  // namespace psi
