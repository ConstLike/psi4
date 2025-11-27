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

    // Compute number of core (doubly occupied) orbitals
    // For REKS(2,2): Ncore = nalpha - 1 (HOMO becomes active)
    Ncore_ = nalpha_ - 1;

    // Create REKS(2,2) active space with default SA weights
    // Initial FONs: n_r=2.0, n_s=0.0 (closed-shell limit)
    double w_PPS = 0.5;
    double w_OSS = 0.5;
    active_space_ = reks::REKSActiveSpace::create_2_2(w_PPS, w_OSS);

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
        outfile->Printf("  Initial FON: n_r = %.6f, n_s = %.6f\n", get_n_r(), get_n_s());
        outfile->Printf("  SA weights:  w_PPS = %.4f, w_OSS = %.4f\n",
                        active_space_->w_PPS(), active_space_->w_OSS());
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

    std::vector<SharedMatrix>& C_left = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();

    for (int p = 0; p < n_base; ++p) {
        C_left.push_back(C_base_[p]);
    }
    C_left.push_back(C_occ);  // Last one for RHF G_

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // Handle special case: pattern 0 with Ncore=0 may have garbage
    if (Ncore_ == 0) {
        J[0]->zero();
        K[0]->zero();
    }

    // HF exchange fraction: 0.0 for pure DFT, 0.5 for BHHLYP, 1.0 for HF
    double alpha = functional_->x_alpha();

    // Build RHF J_, K_, G_ for DIIS and compute_E() compatibility
    // RHF convention: G_ = 2*J - alpha*K + V_xc (if DFT)
    J_ = J[n_base];
    K_ = K[n_base];

    // Compute V_xc for RHF density (for G_ and DIIS compatibility)
    if (functional_->needs_xc() && potential_) {
        potential_->set_D({Da_});
        potential_->compute_V({Va_});
        G_->copy(Va_);
    } else {
        G_->zero();
    }
    G_->axpy(2.0, J_);
    if (functional_->is_x_hybrid()) {
        G_->axpy(-alpha, K_);
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
        J_work_->copy(J[alpha_idx]);
        J_work_->add(J[beta_idx]);

        // F^alpha = H + J_total - alpha*K^alpha
        F_alpha_micro_[L]->copy(H_);
        F_alpha_micro_[L]->add(J_work_);
        F_alpha_micro_[L]->axpy(-alpha, K[alpha_idx]);

        // F^beta = H + J_total - alpha*K^beta
        F_beta_micro_[L]->copy(H_);
        F_beta_micro_[L]->add(J_work_);
        F_beta_micro_[L]->axpy(-alpha, K[beta_idx]);

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

    double E_nuc = energies_["Nuclear"];
    bool needs_xc = functional_->needs_xc();
    double alpha = functional_->x_alpha();

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === compute_microstate_energies DEBUG ===\n");
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
        outfile->Printf("  DEBUG_PSI4_E_MICRO: E[0]=%.10f E[1]=%.10f E[2]=%.10f E[3]=%.10f\n",
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
    // Optimize n_r to minimize E_SA = sum_L C_L(n_r) * E_L
    // Constraint: n_s = 2 - n_r
    //
    // Uses multi-start optimization to avoid getting stuck at boundaries.
    // At n_r=2 (n_s=0), the gradient df/dn_r = 0, which can trap the optimizer.

    // Skip if E_micro_ not yet computed (first SCF iteration with SAD guess)
    if (std::abs(E_micro_[0]) < 1e-10 && std::abs(E_micro_[1]) < 1e-10) {
        if (reks_debug_ >= 1) {
            outfile->Printf("\n  === RexSolver: skipping (E_micro not yet computed) ===\n");
        }
        return;
    }

    const int max_iter = 50;
    const double tol = 1e-10;
    const double max_step = 0.2;

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

    // Multi-start: try starting from current n_r AND from n_r=1.0 (symmetric point)
    double n_r_current = get_n_r();
    std::vector<double> start_points = {n_r_current};

    // If at boundary, also try symmetric point
    if (n_r_current > 1.9 || n_r_current < 0.1) {
        start_points.push_back(1.0);  // Symmetric point n_r = n_s = 1
    }

    double best_n_r = n_r_current;
    double best_E_SA = compute_E_SA(n_r_current);

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver: FON Optimization ===\n");
        outfile->Printf("  E_micro: [0]=%.8f, [1]=%.8f, [2]=%.8f, [3]=%.8f\n",
                       E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
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
            // since x = n_r*n_s = n_r*(2-n_r) and dx/dn_r = 2 - 2*n_r = n_s - n_r
            double dx_dnr = n_s - n_r;  // = 2 - 2*n_r
            double df_dnr = df_dx * dx_dnr;
            // d²f/dn_r² = d²f/dx² * (dx/dn_r)² + df/dx * d²x/dn_r²
            //           = d²f/dx² * (n_s - n_r)² + df/dx * (-2)
            double d2f_dnr2 = d2f_dx2 * dx_dnr * dx_dnr - 2.0 * df_dx;

            // Compute C_L at current n_r
            double C0 = w_PPS * n_r / 2.0;
            double C1 = w_PPS * n_s / 2.0;
            double C2 = w_OSS - 0.5 * w_PPS * f;
            double C3 = 0.5 * w_PPS * f - 0.5 * w_OSS;

            double E_SA_current = C0 * E_micro_[0] + C1 * E_micro_[1]
                                + 2.0 * C2 * E_micro_[2] + 2.0 * C3 * E_micro_[3];

            // Energy gradient dE_SA/dn_r
            // dC0/dn_r = w_PPS/2, dC1/dn_r = -w_PPS/2
            // dC2/dn_r = -0.5*w_PPS*df/dn_r, dC3/dn_r = 0.5*w_PPS*df/dn_r
            double dE_dnr = (w_PPS / 2.0) * (E_micro_[0] - E_micro_[1])
                          + w_PPS * df_dnr * (E_micro_[3] - E_micro_[2]);

            // Energy Hessian
            double d2E_dnr2 = w_PPS * d2f_dnr2 * (E_micro_[3] - E_micro_[2]);

            // Newton-Raphson step (or gradient descent if Hessian ~ 0)
            double delta;
            if (std::abs(d2E_dnr2) < 1e-10) {
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

    // Use the best n_r found
    active_space_->set_pair_fon(0, best_n_r);
    compute_weighting_factors();

    if (reks_debug_ >= 1) {
        outfile->Printf("  Best: n_r = %12.8f, n_s = %12.8f, E_SA = %20.12f\n",
                        get_n_r(), get_n_s(), best_E_SA);
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
    //   - Active r diagonal: weighted by n^σ_r / f_r
    //   - Active s diagonal: weighted by n^σ_s / f_s
    //   - Core-active coupling: weighted by (1-n^σ_q) / (1-f_q)
    //   - Active-active (r,s): sign(f_r - f_s) weighted
    //
    // NOTE: Currently specialized for REKS(2,2). Generalization for REKS(4,4)+
    // requires extending to multiple active orbital pairs.

    double** Freks = F_reks_->pointer(0);

    // Local variables for active orbital indices (from generalized structure)
    // For REKS(2,2): active_r = active_mo_indices_[0], active_s = active_mo_indices_[1]
    int active_r = active_mo_indices_[0];
    int active_s = active_mo_indices_[1];

    // Get microstate occupations from active_space_
    int nra = active_space_->alpha_occ(L, 0);  // n_r^alpha
    int nrb = active_space_->beta_occ(L, 0);   // n_r^beta
    int nsa = active_space_->alpha_occ(L, 1);  // n_s^alpha
    int nsb = active_space_->beta_occ(L, 1);   // n_s^beta

    // Get current FONs and SA weights from active_space_
    double n_r = get_n_r();
    double n_s = get_n_s();
    double w_PPS = active_space_->w_PPS();
    double w_OSS = active_space_->w_OSS();

    // Use generalized effective FON (delegates to active_space_)
    double fr = active_space_->get_effective_fon(0);  // f_r
    double fs = active_space_->get_effective_fon(1);  // f_s

    double Wc = 0.5 * Cl;  // Core block weight

    // Core-core block (i,j < Ncore)
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = i; j < Ncore_; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // Virtual-virtual block (i,j > active_s)
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

    // Active orbital r diagonal
    // Weight: 0.5 * C_L * n^sigma_r / f_r
    double Wr_a = (fr > 1e-10) ? 0.5 * Cl * nra / fr : 0.0;
    double Wr_b = (fr > 1e-10) ? 0.5 * Cl * nrb / fr : 0.0;
    Freks[active_r][active_r] += Wr_a * Fa[active_r][active_r] + Wr_b * Fb[active_r][active_r];

    // Active orbital s diagonal
    // Weight: 0.5 * C_L * n^sigma_s / f_s
    double Ws_a = (fs > 1e-10) ? 0.5 * Cl * nsa / fs : 0.0;
    double Ws_b = (fs > 1e-10) ? 0.5 * Cl * nsb / fs : 0.0;
    double contrib_ss = Ws_a * Fa[active_s][active_s] + Ws_b * Fb[active_s][active_s];
    Freks[active_s][active_s] += contrib_ss;

    // Debug: Trace F(s,s) contributions in MO basis
    if (reks_debug_ >= 2) {
        outfile->Printf("  DEBUG_FM2FE_MO L=%d: Cl=%.6f fr=%.6f fs=%.6f nsa=%d nsb=%d\n",
                       L, Cl, fr, fs, nsa, nsb);
        outfile->Printf("    Fa_MO(s,s)=%.6f Fb_MO(s,s)=%.6f Ws_a=%.6f Ws_b=%.6f\n",
                       Fa[active_s][active_s], Fb[active_s][active_s], Ws_a, Ws_b);
        outfile->Printf("    contrib_ss=%.6f cumul_F_MO(s,s)=%.6f\n",
                       contrib_ss, Freks[active_s][active_s]);
    }

    // Core-active coupling (core to r)
    // Weight: 0.5 * C_L * (1 - n^sigma_r) / (1 - f_r)
    double omfr = 1.0 - fr;  // 1 - f_r
    double Wcr_a = (omfr > 1e-10) ? 0.5 * Cl * (1 - nra) / omfr : 0.0;
    double Wcr_b = (omfr > 1e-10) ? 0.5 * Cl * (1 - nrb) / omfr : 0.0;

    // Core-active coupling (core to s)
    double omfs = 1.0 - fs;  // 1 - f_s
    double Wcs_a = (omfs > 1e-10) ? 0.5 * Cl * (1 - nsa) / omfs : 0.0;
    double Wcs_b = (omfs > 1e-10) ? 0.5 * Cl * (1 - nsb) / omfs : 0.0;

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
    // Weight: C_L * (n^sigma_r - n^sigma_s) * sign(f_r - f_s)
    int signrs = (fr > fs) ? 1 : ((fr < fs) ? -1 : 0);
    double Wrs_a = Cl * (nra - nsa) * signrs;
    double Wrs_b = Cl * (nrb - nsb) * signrs;
    Freks[active_r][active_s] += Wrs_a * Fa[active_r][active_s] + Wrs_b * Fb[active_r][active_s];

    // Accumulate Lagrange multiplier for SI (state interaction)
    Wrs_lagr_ += Cl * (nra * Fa[active_r][active_s] + nrb * Fb[active_r][active_s]);

    // Note: Symmetrization will be done ONCE in form_F() after all microstates are accumulated
}

// ---------------------------------------------------------------------------
// SCF Method Overrides
// ---------------------------------------------------------------------------

void REKS::form_D() {
    // 1. Build standard RHF density Da_ (needed for DIIS, convergence checks)
    RHF::form_D();

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

    // Zero F_reks_ and Lagrange multiplier
    F_reks_->zero();
    Wrs_lagr_ = 0.0;

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
        outfile->Printf("  DEBUG_PSI4_FREKS_MO_DIAG: ");
        for (int i = 0; i < std::min(8, N); ++i) {
            outfile->Printf("F(%d)=%10.6f ", i, Fmo[i][i]);
        }
        outfile->Printf("\n");
        outfile->Printf("  DEBUG_PSI4_FREKS_MO_ACTIVE: F(r,r)=%10.6f F(s,s)=%10.6f F(r,s)=%10.6f\n",
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
    // For SAD guess iteration: use standard RHF diagonalization
    if (sad_ && iteration_ <= 0) {
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
    // This provides:
    // - 2SI-2SA: Ground state (S0) and first excited singlet (S1) energies
    // - 3SI-2SA: Also includes doubly excited state (S2)

    // Skip if energies are not computed yet
    if (E_micro_.empty() || std::abs(E_micro_[0]) < 1e-10) {
        return;
    }

    outfile->Printf("\n  ============================================================\n");
    outfile->Printf("                    SI-SA-REKS State Energies\n");
    outfile->Printf("  ------------------------------------------------------------\n");

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

    outfile->Printf("  ============================================================\n\n");
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

}  // namespace scf
}  // namespace psi
