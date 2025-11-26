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
    // Initial FONs: n_r=2.0, n_s=0.0 (closed-shell limit like GAMESS)
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

    // Pre-allocate work matrices for AO→MO transforms in form_F() (HPC optimization)
    Temp_work_ = std::make_shared<Matrix>("Temp_work", nso, nso);
    F_alpha_MO_.resize(n_micro);
    F_beta_MO_.resize(n_micro);

    for (int L = 0; L < n_micro; ++L) {
        std::string suffix = " L=" + std::to_string(L);
        F_alpha_MO_[L] = std::make_shared<Matrix>("F_alpha_MO" + suffix, nso, nso);
        F_beta_MO_[L] = std::make_shared<Matrix>("F_beta_MO" + suffix, nso, nso);
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
    // Delegate to active_space_ (GAMESS REXCM function)
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

    double alpha = functional_->is_x_hybrid() ? functional_->x_alpha() : 1.0;

    // Build RHF J_, K_, G_ for DIIS and compute_E() compatibility
    // RHF convention: G_ = 2*J - alpha*K
    J_ = J[n_base];
    K_ = K[n_base];

    G_->zero();
    G_->axpy(2.0, J_);
    if (functional_->is_x_hybrid()) {
        G_->axpy(-alpha, K_);
    }

    // Build Fock matrices for each microstate
    // F^sigma_L = H + J_total_L - alpha*K^sigma_L
    // where J_total = J^alpha + J^beta, K^sigma uses the sigma density

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
    //     + E_nuc
    // where G^sigma = F^sigma - H = J - alpha*K^sigma

    double E_nuc = energies_["Nuclear"];

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
// RexSolver - Newton-Raphson FON Optimization (Filatov 2024)
// ---------------------------------------------------------------------------

void REKS::rex_solver() {
    // Optimize n_r to minimize E_SA = sum_L C_L(n_r) * E_L
    // Constraint: n_s = 2 - n_r
    //
    // For SA-REKS energy (w_PPS, w_OSS weights):
    //   E_SA = C_0*E_0 + C_1*E_1 + C_2*E_2 + C_3*E_3
    // where C_L depend on n_r through f(n_r/2)
    //
    // Gradient and Hessian w.r.t. n_r are computed analytically.

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

    // Get current FONs and SA weights from active_space_
    double n_r = get_n_r();
    double n_s = get_n_s();
    double w_PPS = active_space_->w_PPS();
    double w_OSS = active_space_->w_OSS();

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver: FON Optimization ===\n");
        outfile->Printf("  Initial: n_r=%10.6f, n_s=%10.6f\n", n_r, n_s);
        outfile->Printf("  E_micro: [0]=%.8f, [1]=%.8f, [2]=%.8f, [3]=%.8f\n",
                       E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
        if (reks_debug_ >= 2) {
            outfile->Printf("\n  Iter   n_r      n_s      f(x)      C_L[0]    C_L[1]    C_L[2]    C_L[3]   Sum(C_L)    E_SA         dE/dn    d2E/dn2   delta\n");
            outfile->Printf("  ---- -------- -------- -------- --------- --------- --------- --------- -------- ------------ --------- --------- --------\n");
        }
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        n_s = 2.0 - n_r;
        active_space_->set_pair_fon(0, n_r);  // Updates both fon_p and fon_q

        // Interpolating function and derivatives at x = n_r*n_s
        // GAMESS uses f(n_r*n_s), not f(n_r/2)
        double x = n_r * n_s;  // x = n_r*n_s in [0, 1]
        double f = reks::f_interp(x);
        double df_dx = reks::df_interp(x);
        double d2f_dx2 = reks::d2f_interp(x);

        // Chain rule: df/dn_r = (df/dx) * (dx/dn_r) = (df/dx) * n_s
        // (since x = n_r*n_s and dx/dn_r = n_s)
        double df_dnr = df_dx * n_s;
        // d²f/dn_r² = d/dn_r[df/dx * n_s] = d²f/dx² * n_s * dx/dn_r + df/dx * dn_s/dn_r
        //           = d²f/dx² * n_s² + df/dx * (-1)  (since n_s = 2-n_r)
        double d2f_dnr2 = d2f_dx2 * n_s * n_s - df_dx;

        // Compute C_L at current n_r (before update)
        double C_L_current[4];
        C_L_current[0] = w_PPS * n_r / 2.0;
        C_L_current[1] = w_PPS * n_s / 2.0;
        C_L_current[2] = w_OSS - 0.5 * w_PPS * f;
        C_L_current[3] = 0.5 * w_PPS * f - 0.5 * w_OSS;
        double sum_CL = C_L_current[0] + C_L_current[1] + C_L_current[2] + C_L_current[3];

        // GAMESS applies FACT=2 for L>=2 in REXEM2EE
        double E_SA_current = C_L_current[0] * E_micro_[0] + C_L_current[1] * E_micro_[1]
                            + 2.0 * C_L_current[2] * E_micro_[2] + 2.0 * C_L_current[3] * E_micro_[3];

        // Energy gradient dE_SA/dn_r
        // C_L[0] = w_PPS * n_r/2             -> dC_0/dn_r = w_PPS/2
        // C_L[1] = w_PPS * n_s/2             -> dC_1/dn_r = -w_PPS/2
        // C_L[2] = w_OSS - 0.5*w_PPS*f       -> dC_2/dn_r = -0.5*w_PPS*df/dn_r
        // C_L[3] = 0.5*w_PPS*f - 0.5*w_OSS   -> dC_3/dn_r = 0.5*w_PPS*df/dn_r
        //
        // dE/dn_r = sum_L FACT * (dC_L/dn_r) * E_L, FACT=2 for L>=2
        double dE_dnr = (w_PPS / 2.0) * (E_micro_[0] - E_micro_[1])
                      + 2.0 * 0.5 * w_PPS * df_dnr * (E_micro_[3] - E_micro_[2]);

        // Energy Hessian d2E_SA/dn_r2
        // Only C_L[2] and C_L[3] have second derivatives (through f)
        // d2C_2/dn_r2 = -0.5*w_PPS * d2f/dn_r2
        // d2C_3/dn_r2 = 0.5*w_PPS * d2f/dn_r2
        // Apply FACT=2 for L>=2
        double d2E_dnr2 = 2.0 * 0.5 * w_PPS * d2f_dnr2 * (E_micro_[3] - E_micro_[2]);

        // Newton-Raphson step (or gradient descent if Hessian ~ 0)
        double delta;
        if (std::abs(d2E_dnr2) < 1e-10) {
            // Hessian near zero (saddle point at x=0.5) - use gradient descent
            delta = -std::copysign(max_step, dE_dnr);
            if (reks_debug_ >= 2) {
                outfile->Printf("    iter %2d: zero Hessian, gradient descent\n", iter);
            }
        } else {
            // Standard Newton-Raphson
            delta = -dE_dnr / d2E_dnr2;
        }

        // Damping for large steps
        if (std::abs(delta) > max_step) {
            delta = std::copysign(max_step, delta);
        }

        // Print before update
        if (reks_debug_ >= 2) {
            outfile->Printf("  %4d %8.5f %8.5f %8.5f %9.6f %9.6f %9.6f %9.6f %8.5f %12.8f %9.2e %9.2e",
                           iter, n_r, n_s, f,
                           C_L_current[0], C_L_current[1], C_L_current[2], C_L_current[3],
                           sum_CL, E_SA_current, dE_dnr, d2E_dnr2);
        }

        // Update n_r with bounds [0.0, 2.0] (same as GAMESS x2 ∈ [0, 1])
        double n_r_new = n_r + delta;
        n_r_new = std::clamp(n_r_new, 0.0, 2.0);
        delta = n_r_new - n_r;
        n_r = n_r_new;

        if (reks_debug_ >= 2) {
            outfile->Printf(" %8.5f\n", delta);
        }

        // Convergence check
        if (std::abs(delta) < tol) {
            if (reks_debug_ >= 1) {
                outfile->Printf("  RexSolver converged in %d iterations\n", iter + 1);
            }
            break;
        }
    }

    // Update active_space_ with final FON
    active_space_->set_pair_fon(0, n_r);

    // Compute SA-REKS energy with optimized FON
    compute_weighting_factors();  // Update C_L with new n_r

    // GAMESS REXEM2EE multiplies by 2 for L>=3 (1-based), i.e. L>=2 in 0-based
    double E_SA = C_L_[0] * E_micro_[0] + C_L_[1] * E_micro_[1]
                + 2.0 * C_L_[2] * E_micro_[2] + 2.0 * C_L_[3] * E_micro_[3];

    if (reks_debug_ >= 1) {
        outfile->Printf("  Final: n_r = %12.8f, n_s = %12.8f, f(x) = %12.8f\n",
                        get_n_r(), get_n_s(), reks::f_interp(get_n_r() * get_n_s()));
        outfile->Printf("  E_SA = %20.12f\n", E_SA);
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
    //   - Core-core, virtual-virtual, core-virtual: F_c = 0.5*C_L*(F^α + F^β)
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
    // For SAD guess iteration: Ca_ doesn't exist yet, delegate to RHF
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
        RHF::form_F();
        return;
    }

    // Build REKS coupling Fock matrix F_reks_ from microstate Fock matrices.
    // The REKS coupling formulas (Filatov 2024) are designed for MO basis.
    // GAMESS procedure:
    //   1. Build F^sigma_L in AO basis
    //   2. Transform to MO: F^sigma_L_MO = C^T * F^sigma_L_AO * C
    //   3. Apply fock_micro_to_macro in MO basis (active_mo_indices_ are MO indices)
    //   4. Transform F_reks back to AO: F_reks_AO = C * F_reks_MO * C^T

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
    // Note: L=2 (open-shell) represents L=3 and L=4 in paper → factor 2 already in C_L_[2]
    //       L=3 (triplet-like) represents L=5 and L=6 → factor 2 already in C_L_[3]
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

    // Store F_reks_MO for GAMESS-style orbital update in form_C()
    F_reks_MO_ = F_reks_->clone();
    F_reks_MO_->set_name("F_REKS_MO");

    // Step 4: Transform F_reks from MO back to AO basis for DIIS
    // GAMESS formula (reks.src lines 1363-1367): F_AO = S * C * F_MO * C^T * S
    // The commutator [F,P] = F*D*S - S*D*F goes to zero at convergence
    // only if F_AO is constructed this way.
    // temp = S * C, then F_AO = temp * F_MO * temp^T
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

        // Print MO-basis F_reks diagonal - this is what GAMESS prints!
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

    // GAMESS-style orbital update with DIIS support (reks.src lines 1380-1420):
    //
    // After form_F(), Psi4's DIIS may have extrapolated Fa_ (which was set to F_AO = S*C*F_MO*C^T*S).
    // GAMESS converts F_AO back to MO basis (lines 1380-1386) before diagonalization.
    //
    // The back-transformation is: F_MO = C^T * F_AO * C
    // This works because C is S-orthonormal: C^T * S * C = I
    // So: C^T * (S*C*F_MO*C^T*S) * C = (C^T*S*C) * F_MO * (C^T*S*C) = I * F_MO * I = F_MO
    //
    // After DIIS extrapolation, F_AO may have changed, so we get:
    // F_MO_after_DIIS = C^T * Fa_diis * C

    int N = nsopi_[0];

    // Convert Fa_ (possibly DIIS-extrapolated) back to MO basis: F_MO = C^T * Fa_ * C
    auto F_MO_diis = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);

    // Diagonalize F_MO: F_MO * U = U * ε
    auto U = std::make_shared<Matrix>("Eigenvectors", N, N);
    auto eps = std::make_shared<Vector>("Eigenvalues", N);

    F_MO_diis->diagonalize(U, eps);

    // Fix orbital phase (GAMESS lines 1413-1417): ensure diagonal elements of U are positive
    double** Up = U->pointer(0);
    for (int j = 0; j < N; ++j) {
        if (Up[j][j] < 0.0) {
            // Flip sign of this column
            for (int i = 0; i < N; ++i) {
                Up[i][j] = -Up[i][j];
            }
        }
    }

    // Update orbitals: C_new = C_old * U (GAMESS lines 1419-1420)
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
        outfile->Printf("\n  === GAMESS-style Orbital Update ===\n");
        // Get active MO indices for debug output
        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];
        outfile->Printf("  Eigenvalues: ");
        for (int i = std::max(0, active_r - 2); i < std::min(N, active_s + 3); ++i) {
            outfile->Printf("ε(%d)=%10.6f ", i, eps_p[i]);
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
    //
    // Reference: GAMESS reks.src lines 1784-1850

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
