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
 * Current Implementation Status:
 * - Phase 1: REKS class structure inheriting from RHF
 * - All calculations delegate to RHF (functionally equivalent to RHF)
 * - Foundation for REKS-specific features in future phases
 *
 * Future Development Phases:
 * - Phase 2: Active space definition and fractional occupations
 * - Phase 3: Ensemble energy functional
 * - Phase 4: SI-SA-REKS (state-interaction state-averaged REKS)
 * - Phase 5: Analytical gradients
 */

#include "reks.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfunctional/superfunctional.h"

#include <cmath>
#include <algorithm>

namespace psi {
namespace scf {

// Define static constexpr members (required for ODR-use in C++14/17)
constexpr int REKS::nr_alpha_[4];
constexpr int REKS::nr_beta_[4];
constexpr int REKS::ns_alpha_[4];
constexpr int REKS::ns_beta_[4];

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

    // Compute active space indices
    // For REKS(2,2): 2 electrons in 2 orbitals at HOMO/LUMO
    Ncore_ = nalpha_ - 1;  // All doubly occupied except HOMO
    active_r_ = Ncore_;    // HOMO
    active_s_ = Ncore_ + 1; // LUMO (will be fractionally occupied)

    // Initialize FONs (start at closed-shell limit like GAMESS)
    // GAMESS: DNR=1.0, DNS=0.0 → n_r=2.0, n_s=0.0
    n_r_ = 2.0;
    n_s_ = 0.0;

    // Initialize SA weights (default: equal weighting)
    w_PPS_ = 0.5;
    w_OSS_ = 0.5;

    // Initialize microstate energies
    E_micro_.fill(0.0);
    C_L_.fill(0.0);

    // Allocate REKS-specific matrices
    allocate_reks_matrices();

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  ╔══════════════════════════════════════════════════════════╗\n");
        outfile->Printf("  ║              REKS(2,2) Initialization                     ║\n");
        outfile->Printf("  ╠══════════════════════════════════════════════════════════╣\n");
        outfile->Printf("  ║  Active space: %d electrons in %d orbitals               ║\n",
                        n_active_electrons_, n_active_orbitals_);
        outfile->Printf("  ║  Core orbitals (Ncore): %3d                              ║\n", Ncore_);
        outfile->Printf("  ║  Active orbital r:      %3d (HOMO)                       ║\n", active_r_);
        outfile->Printf("  ║  Active orbital s:      %3d (LUMO)                       ║\n", active_s_);
        outfile->Printf("  ╠══════════════════════════════════════════════════════════╣\n");
        outfile->Printf("  ║  Initial FON: n_r = %.6f, n_s = %.6f              ║\n", n_r_, n_s_);
        outfile->Printf("  ║  SA weights:  w_PPS = %.4f, w_OSS = %.4f               ║\n", w_PPS_, w_OSS_);
        outfile->Printf("  ║  Filatov delta parameter: %.4f                          ║\n", DELTA_);
        outfile->Printf("  ╠══════════════════════════════════════════════════════════╣\n");
        outfile->Printf("  ║  f(0.5) = %.8f  (should be 1.0)                    ║\n", f_interp(0.5));
        outfile->Printf("  ║  f(0.25) = %.8f                                    ║\n", f_interp(0.25));
        outfile->Printf("  ║  f(0.0) = %.8f  (should be 0.0)                    ║\n", f_interp(0.0));
        outfile->Printf("  ╚══════════════════════════════════════════════════════════╝\n\n");
    }

    // NOTE: Ca_ is allocated here but EMPTY (zeros).
    // Actual orbitals come from guess() in compute_energy().
    // Debug output for densities will appear in SCF iterations.
}

void REKS::allocate_reks_matrices() {
    // Allocate base density matrices
    D00_ = std::make_shared<Matrix>("D00 (core)", nsopi_, nsopi_);
    D10_ = std::make_shared<Matrix>("D10 (core+r)", nsopi_, nsopi_);
    D01_ = std::make_shared<Matrix>("D01 (core+s)", nsopi_, nsopi_);
    D11_ = std::make_shared<Matrix>("D11 (core+r+s)", nsopi_, nsopi_);

    // Allocate microstate density and Fock matrices
    for (int L = 0; L < N_MICRO_; ++L) {
        std::string suffix = " L=" + std::to_string(L);
        D_alpha_micro_[L] = std::make_shared<Matrix>("D_alpha" + suffix, nsopi_, nsopi_);
        D_beta_micro_[L] = std::make_shared<Matrix>("D_beta" + suffix, nsopi_, nsopi_);
        F_alpha_micro_[L] = std::make_shared<Matrix>("F_alpha" + suffix, nsopi_, nsopi_);
        F_beta_micro_[L] = std::make_shared<Matrix>("F_beta" + suffix, nsopi_, nsopi_);
    }

    // Allocate coupling Fock matrix
    F_reks_ = std::make_shared<Matrix>("F_REKS (coupling)", nsopi_, nsopi_);
}

// ═══════════════════════════════════════════════════════════════════════════
// Interpolating Function f(x) and Derivatives (Filatov 2024, Eq. 4)
// ═══════════════════════════════════════════════════════════════════════════

double REKS::f_interp(double x) const {
    // GAMESS REXCONVF function (used for x = n_r*n_s)
    // f(x) = x^(1 - 0.5*(x+δ)/(1+δ))
    // Domain: 0 <= x <= 1, where x = n_r*n_s
    // At x = 0: f = 0 (no coupling)
    // At x = 1: f = 1 (maximum coupling, n_r=n_s=1)

    if (x <= 0.0) return 0.0;
    if (x >= 1.0) {
        // At x=1: TMP = 0.5, so f = 1^0.5 = 1
        return 1.0;
    }

    double tmp = 0.5 * (x + DELTA_) / (1.0 + DELTA_);
    double exponent = 1.0 - tmp;
    return std::pow(x, exponent);
}

double REKS::df_interp(double x) const {
    // df/dx for f(x) = x^(1-tmp) where tmp = 0.5*(x+δ)/(1+δ)
    // Using logarithmic derivative: ln(f) = (1-tmp)*ln(x)
    // d[ln(f)]/dx = (1-tmp)/x - 0.5*ln(x)/(1+δ)
    // df/dx = f(x) * [(1-tmp)/x - 0.5*ln(x)/(1+δ)]

    if (x <= 1e-10) return 0.0;
    if (x >= 1.0 - 1e-10) return 0.0;  // At x=1, derivative is 0

    double tmp = 0.5 * (x + DELTA_) / (1.0 + DELTA_);
    double f_val = std::pow(x, 1.0 - tmp);

    double term1 = (1.0 - tmp) / x;
    double term2 = -0.5 * std::log(x) / (1.0 + DELTA_);

    return f_val * (term1 + term2);
}

double REKS::d2f_interp(double x) const {
    // d²f/dx² - numerical differentiation for robustness
    const double h = 1e-5;
    double fp = df_interp(x + h);
    double fm = df_interp(x - h);
    return (fp - fm) / (2.0 * h);
}

// ═══════════════════════════════════════════════════════════════════════════
// Density Matrix Construction (Filatov 2024, Section 2)
// ═══════════════════════════════════════════════════════════════════════════

void REKS::build_base_densities() {
    // Builds D00, D10, D01, D11 from current orbitals Ca_
    // D00 = core only (Ncore doubly occupied)
    // D10 = core + orbital r singly occupied
    // D01 = core + orbital s singly occupied
    // D11 = core + r + s both singly occupied

    int nso = nsopi_[0];  // Assuming C1 symmetry

    // D00: core electrons only
    D00_->zero();
    if (Ncore_ > 0) {
        double** D00p = D00_->pointer(0);
        double** Cp = Ca_->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                double val = 0.0;
                for (int k = 0; k < Ncore_; ++k) {
                    val += Cp[mu][k] * Cp[nu][k];
                }
                D00p[mu][nu] = val;
            }
        }
    }

    // D10: core + orbital r
    D10_->copy(D00_);
    {
        double** D10p = D10_->pointer(0);
        double** Cp = Ca_->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                D10p[mu][nu] += Cp[mu][active_r_] * Cp[nu][active_r_];
            }
        }
    }

    // D01: core + orbital s
    D01_->copy(D00_);
    {
        double** D01p = D01_->pointer(0);
        double** Cp = Ca_->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                D01p[mu][nu] += Cp[mu][active_s_] * Cp[nu][active_s_];
            }
        }
    }

    // D11: core + r + s (both singly occupied)
    D11_->copy(D10_);
    {
        double** D11p = D11_->pointer(0);
        double** Cp = Ca_->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                D11p[mu][nu] += Cp[mu][active_s_] * Cp[nu][active_s_];
            }
        }
    }

    if (reks_debug_ >= 2) {
        // Note: tr(D) != N_electrons in non-orthogonal basis. Use tr(D*S) for electrons.
        outfile->Printf("\n  === REKS Base Densities ===\n");
        outfile->Printf("  D00 (core):    tr(D*S) = %8.4f, N_elec = %d\n",
                        D00_->vector_dot(S_), Ncore_);
        outfile->Printf("  D10 (core+r):  tr(D*S) = %8.4f, N_elec = %d\n",
                        D10_->vector_dot(S_), Ncore_ + 1);
        outfile->Printf("  D01 (core+s):  tr(D*S) = %8.4f, N_elec = %d\n",
                        D01_->vector_dot(S_), Ncore_ + 1);
        outfile->Printf("  D11 (core+rs): tr(D*S) = %8.4f, N_elec = %d\n",
                        D11_->vector_dot(S_), Ncore_ + 2);
    }
}

void REKS::build_microstate_densities() {
    // Map base densities to microstate α/β densities
    // Only 4 unique microstates needed (L=3≡L=4, L=5≡L=6)
    //
    // L=0 (paper L=1): D^α = D10, D^β = D10  (closed-shell n_r)
    // L=1 (paper L=2): D^α = D01, D^β = D01  (closed-shell n_s)
    // L=2 (paper L=3): D^α = D10, D^β = D01  (open-shell r↑s↓)
    // L=3 (paper L=5): D^α = D11, D^β = D00  (triplet-like)

    D_alpha_micro_[0]->copy(D10_);
    D_beta_micro_[0]->copy(D10_);

    D_alpha_micro_[1]->copy(D01_);
    D_beta_micro_[1]->copy(D01_);

    D_alpha_micro_[2]->copy(D10_);
    D_beta_micro_[2]->copy(D01_);

    D_alpha_micro_[3]->copy(D11_);
    D_beta_micro_[3]->copy(D00_);

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Microstate Densities ===\n");
        for (int L = 0; L < N_MICRO_; ++L) {
            outfile->Printf("  L=%d: tr(D^α) = %10.6f, tr(D^β) = %10.6f, "
                           "n_r^α=%d, n_r^β=%d, n_s^α=%d, n_s^β=%d\n",
                L, D_alpha_micro_[L]->trace(), D_beta_micro_[L]->trace(),
                nr_alpha_[L], nr_beta_[L], ns_alpha_[L], ns_beta_[L]);
        }
    }
}

void REKS::compute_weighting_factors() {
    // Compute C_L weights for SA-REKS (GAMESS REXCM function)
    // CRITICAL: f_interp argument is n_r*n_s, NOT n_r/2!
    // GAMESS: TMP = 4*DNR*DNS = 4*(n_r/2)*(n_s/2) = n_r*n_s
    double f = f_interp(n_r_ * n_s_);  // f(n_r*n_s), NOT f(n_r/2)!

    // C_L for SA-REKS (STATAVG branch, wpps≠1.0)
    // L=0: WPPS*DNR = w_PPS*(n_r/2)
    C_L_[0] = w_PPS_ * n_r_ / 2.0;

    // L=1: WPPS*DNS = w_PPS*(n_s/2)
    C_L_[1] = w_PPS_ * n_s_ / 2.0;

    // L=2 (open, ×2): WOSS - 0.5*F = w_OSS - 0.5*w_PPS*f(n_r*n_s)
    C_L_[2] = w_OSS_ - 0.5 * w_PPS_ * f;

    // L=3 (triplet, ×2): 0.5*F - 0.5*WOSS = 0.5*w_PPS*f(n_r*n_s) - 0.5*w_OSS
    C_L_[3] = 0.5 * w_PPS_ * f - 0.5 * w_OSS_;

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Weighting Factors C_L ===\n");
        outfile->Printf("  n_r = %.6f, n_s = %.6f, f(n_r*n_s) = %.6f\n", n_r_, n_s_, f);
        outfile->Printf("  w_PPS = %.4f, w_OSS = %.4f\n", w_PPS_, w_OSS_);
        for (int L = 0; L < N_MICRO_; ++L) {
            outfile->Printf("  C_L[%d] = %12.8f\n", L, C_L_[L]);
        }
        double sum = C_L_[0] + C_L_[1] + C_L_[2] + C_L_[3];
        outfile->Printf("  Sum C_L = %12.8f (should be ~1.0)\n", sum);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Step 7: Build Microstate Fock Matrices (Filatov 2024, Section 2)
// ═══════════════════════════════════════════════════════════════════════════

void REKS::build_microstate_focks() {
    // Build UHF-like Fock matrices for each of the 4 unique microstates.
    // Single JK batch: D00, D10, D01, D11 + Da_ (for RHF G_ compatibility)
    //
    // Microstate density mapping:
    //   L=0: D^α = D10, D^β = D10  (closed-shell r)
    //   L=1: D^α = D01, D^β = D01  (closed-shell s)
    //   L=2: D^α = D10, D^β = D01  (open-shell r↑s↓)
    //   L=3: D^α = D11, D^β = D00  (triplet-like)

    int nso = nsopi_[0];

    // Build C matrices for the 4 unique base densities + RHF occupied
    auto C_D00 = std::make_shared<Matrix>("C for D00", nso, Ncore_);  // Fixed: was std::max(1, Ncore_)
    auto C_D10 = std::make_shared<Matrix>("C for D10", nso, Ncore_ + 1);
    auto C_D01 = std::make_shared<Matrix>("C for D01", nso, Ncore_ + 1);
    auto C_D11 = std::make_shared<Matrix>("C for D11", nso, Ncore_ + 2);
    auto C_occ = Ca_subset("SO", "OCC");  // RHF occupied orbitals for G_

    double** Cp = Ca_->pointer(0);

    // C_D00: columns [0..Ncore-1] - core only
    if (Ncore_ > 0) {
        double** C00p = C_D00->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_; ++k) {
                C00p[mu][k] = Cp[mu][k];
            }
        }
    }

    // C_D10: columns [0..Ncore] - core + r (contiguous)
    {
        double** C10p = C_D10->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_ + 1; ++k) {
                C10p[mu][k] = Cp[mu][k];
            }
        }
    }

    // C_D01: columns [0..Ncore-1, s] - core + s (NOT contiguous, s is at Ncore+1)
    {
        double** C01p = C_D01->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_; ++k) {
                C01p[mu][k] = Cp[mu][k];
            }
            // s orbital at position Ncore+1 in Ca_, placed at Ncore in C_D01
            C01p[mu][Ncore_] = Cp[mu][active_s_];
        }
    }

    // C_D11: columns [0..Ncore+1] - core + r + s (contiguous)
    {
        double** C11p = C_D11->pointer(0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int k = 0; k < Ncore_ + 2; ++k) {
                C11p[mu][k] = Cp[mu][k];
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Batch JK computation: 5 densities in one call
    // ═══════════════════════════════════════════════════════════════════════
    //   [0] C_D00 → J00, K00
    //   [1] C_D10 → J10, K10
    //   [2] C_D01 → J01, K01
    //   [3] C_D11 → J11, K11
    //   [4] C_occ → J_rhf, K_rhf (for RHF G_ matrix)

    std::vector<SharedMatrix>& C_left = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();

    if (Ncore_ > 0) {
        C_left.push_back(C_D00);  // [0]
    } else {
        auto C_empty = std::make_shared<Matrix>("C empty", nso, 0);
        C_left.push_back(C_empty);  // [0] D00 = 0 for H2
    }
    C_left.push_back(C_D10);  // [1]
    C_left.push_back(C_D01);  // [2]
    C_left.push_back(C_D11);  // [3]
    C_left.push_back(C_occ);  // [4] RHF occupied for G_

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // J00=J[0], J10=J[1], J01=J[2], J11=J[3], J_rhf=J[4]
    // K00=K[0], K10=K[1], K01=K[2], K11=K[3], K_rhf=K[4]

    // For H2 (Ncore=0), J[0]/K[0] may contain garbage from empty C matrix
    // Explicitly zero them to avoid contaminating microstate energies
    if (Ncore_ == 0) {
        J[0]->zero();
        K[0]->zero();
    }

    double alpha = functional_->is_x_hybrid() ? functional_->x_alpha() : 1.0;

    // ═══════════════════════════════════════════════════════════════════════
    // Build RHF J_, K_, G_ for DIIS and compute_E() compatibility
    // RHF convention: G_ = 2*J - α*K
    // ═══════════════════════════════════════════════════════════════════════
    J_ = J[4];
    K_ = K[4];

    G_->zero();
    G_->axpy(2.0, J_);
    if (functional_->is_x_hybrid()) {
        G_->axpy(-alpha, K_);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Build Fock matrices for each microstate
    // F^σ_L = H + J_total_L - α*K^σ_L
    // ═══════════════════════════════════════════════════════════════════════

    // L=0 (closed-shell r): D^α = D^β = D10
    // J_total = 2*J10, K^α = K^β = K10
    {
        auto J_total = J[1]->clone();
        J_total->scale(2.0);

        F_alpha_micro_[0]->copy(H_);
        F_alpha_micro_[0]->add(J_total);
        F_alpha_micro_[0]->axpy(-alpha, K[1]);

        F_beta_micro_[0]->copy(F_alpha_micro_[0]);  // Same for closed-shell
    }

    // L=1 (closed-shell s): D^α = D^β = D01
    // J_total = 2*J01, K^α = K^β = K01
    {
        auto J_total = J[2]->clone();
        J_total->scale(2.0);

        F_alpha_micro_[1]->copy(H_);
        F_alpha_micro_[1]->add(J_total);
        F_alpha_micro_[1]->axpy(-alpha, K[2]);

        F_beta_micro_[1]->copy(F_alpha_micro_[1]);  // Same for closed-shell
    }

    // L=2 (open-shell r↑s↓): D^α = D10, D^β = D01
    // J_total = J10 + J01, K^α = K10, K^β = K01
    {
        auto J_total = J[1]->clone();
        J_total->add(J[2]);

        F_alpha_micro_[2]->copy(H_);
        F_alpha_micro_[2]->add(J_total);
        F_alpha_micro_[2]->axpy(-alpha, K[1]);

        F_beta_micro_[2]->copy(H_);
        F_beta_micro_[2]->add(J_total);
        F_beta_micro_[2]->axpy(-alpha, K[2]);
    }

    // L=3 (triplet-like): D^α = D11, D^β = D00
    // J_total = J11 + J00, K^α = K11, K^β = K00
    {
        auto J_total = J[3]->clone();
        J_total->add(J[0]);

        F_alpha_micro_[3]->copy(H_);
        F_alpha_micro_[3]->add(J_total);
        F_alpha_micro_[3]->axpy(-alpha, K[3]);

        F_beta_micro_[3]->copy(H_);
        F_beta_micro_[3]->add(J_total);
        F_beta_micro_[3]->axpy(-alpha, K[0]);
    }

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Microstate Fock Matrices ===\n");
        outfile->Printf("  Exchange scaling α = %.4f\n", alpha);

        // Compare with RHF Fock for validation
        Fa_->copy(H_);
        Fa_->add(G_);
        double tr_rhf = Fa_->trace();
        outfile->Printf("  RHF Fock: tr(Fa_RHF) = %15.10f\n", tr_rhf);

        for (int L = 0; L < N_MICRO_; ++L) {
            double tr_alpha = F_alpha_micro_[L]->trace();
            double tr_beta = F_beta_micro_[L]->trace();

            // Check active orbital diagonal elements
            double** Fa = F_alpha_micro_[L]->pointer(0);
            double** Fb = F_beta_micro_[L]->pointer(0);
            double Frr_a = Fa[active_r_][active_r_];
            double Fss_a = Fa[active_s_][active_s_];
            double Frs_a = Fa[active_r_][active_s_];
            double Frr_b = Fb[active_r_][active_r_];
            double Fss_b = Fb[active_s_][active_s_];
            double Frs_b = Fb[active_r_][active_s_];

            outfile->Printf("  L=%d: tr(F^α)=%15.10f, tr(F^β)=%15.10f\n", L, tr_alpha, tr_beta);
            outfile->Printf("        F^α[r,r]=%12.6f, F^α[s,s]=%12.6f, F^α[r,s]=%12.6f\n",
                           Frr_a, Fss_a, Frs_a);
            outfile->Printf("        F^β[r,r]=%12.6f, F^β[s,s]=%12.6f, F^β[r,s]=%12.6f\n",
                           Frr_b, Fss_b, Frs_b);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Step 8: Compute Microstate Energies (Filatov 2024, Eq. 5)
// ═══════════════════════════════════════════════════════════════════════════

void REKS::compute_microstate_energies() {
    // E_L = tr(D^α_L · H) + tr(D^β_L · H)
    //     + 0.5 * [tr(D^α_L · G^α_L) + tr(D^β_L · G^β_L)]
    //     + E_nuc
    // where G^σ = F^σ - H = J - α*K^σ

    double E_nuc = energies_["Nuclear"];

    for (int L = 0; L < N_MICRO_; ++L) {
        SharedMatrix Da = D_alpha_micro_[L];
        SharedMatrix Db = D_beta_micro_[L];
        SharedMatrix Fa = F_alpha_micro_[L];
        SharedMatrix Fb = F_beta_micro_[L];

        // One-electron energy: tr(D^α · H) + tr(D^β · H)
        double E_1e = Da->vector_dot(H_) + Db->vector_dot(H_);

        // Two-electron energy from Fock matrices:
        // E_2e = 0.5 * [tr(D^α · (F^α - H)) + tr(D^β · (F^β - H))]
        //      = 0.5 * [tr(D^α · F^α) + tr(D^β · F^β) - E_1e]
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

// ═══════════════════════════════════════════════════════════════════════════
// Step 9: RexSolver - Newton-Raphson FON Optimization (Filatov 2024)
// ═══════════════════════════════════════════════════════════════════════════

void REKS::rex_solver() {
    // Optimize n_r to minimize E_SA = Σ C_L(n_r) * E_L
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

    if (reks_debug_ >= 1) {
        outfile->Printf("\n  === RexSolver: FON Optimization ===\n");
        outfile->Printf("  Initial: n_r=%10.6f, n_s=%10.6f\n", n_r_, n_s_);
        outfile->Printf("  E_micro: [0]=%.8f, [1]=%.8f, [2]=%.8f, [3]=%.8f\n",
                       E_micro_[0], E_micro_[1], E_micro_[2], E_micro_[3]);
        if (reks_debug_ >= 2) {
            outfile->Printf("\n  Iter   n_r      n_s      f(x)      C_L[0]    C_L[1]    C_L[2]    C_L[3]   Sum(C_L)    E_SA         dE/dn    d2E/dn2   delta\n");
            outfile->Printf("  ---- -------- -------- -------- --------- --------- --------- --------- -------- ------------ --------- --------- --------\n");
        }
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        n_s_ = 2.0 - n_r_;

        // Interpolating function and derivatives at x = n_r*n_s
        // CRITICAL: GAMESS uses f(n_r*n_s), NOT f(n_r/2)!
        double x = n_r_ * n_s_;  // x = n_r*n_s ∈ [0, 1]
        double f = f_interp(x);
        double df_dx = df_interp(x);
        double d2f_dx2 = d2f_interp(x);

        // Chain rule: df/dn_r = (df/dx) * (dx/dn_r) = (df/dx) * n_s
        // (since x = n_r*n_s and dx/dn_r = n_s)
        double df_dnr = df_dx * n_s_;
        // d²f/dn_r² = d/dn_r[df/dx * n_s] = d²f/dx² * n_s * dx/dn_r + df/dx * dn_s/dn_r
        //           = d²f/dx² * n_s² + df/dx * (-1)  (since n_s = 2-n_r)
        double d2f_dnr2 = d2f_dx2 * n_s_ * n_s_ - df_dx;

        // Compute C_L at current n_r (before update)
        double C_L_current[4];
        C_L_current[0] = w_PPS_ * n_r_ / 2.0;
        C_L_current[1] = w_PPS_ * n_s_ / 2.0;
        C_L_current[2] = w_OSS_ - 0.5 * w_PPS_ * f;
        C_L_current[3] = 0.5 * w_PPS_ * f - 0.5 * w_OSS_;
        double sum_CL = C_L_current[0] + C_L_current[1] + C_L_current[2] + C_L_current[3];

        // CRITICAL FIX: Apply FACT=2 for L>=2 as in GAMESS REXEM2EE
        double E_SA_current = C_L_current[0] * E_micro_[0] + C_L_current[1] * E_micro_[1]
                            + 2.0 * C_L_current[2] * E_micro_[2] + 2.0 * C_L_current[3] * E_micro_[3];

        // ─────────── Energy gradient dE_SA/dn_r ───────────
        // C_L[0] = w_PPS * n_r/2             → dC_0/dn_r = w_PPS/2
        // C_L[1] = w_PPS * n_s/2             → dC_1/dn_r = -w_PPS/2
        // C_L[2] = w_OSS - 0.5*w_PPS*f       → dC_2/dn_r = -0.5*w_PPS*df/dn_r
        // C_L[3] = 0.5*w_PPS*f - 0.5*w_OSS   → dC_3/dn_r = 0.5*w_PPS*df/dn_r
        //
        // dE/dn_r = Σ FACT * (dC_L/dn_r) * E_L, FACT=2 for L>=2
        double dE_dnr = (w_PPS_ / 2.0) * (E_micro_[0] - E_micro_[1])
                      + 2.0 * 0.5 * w_PPS_ * df_dnr * (E_micro_[3] - E_micro_[2]);

        // ─────────── Energy Hessian d²E_SA/dn_r² ───────────
        // Only C_L[2] and C_L[3] have second derivatives (through f)
        // d²C_2/dn_r² = -0.5*w_PPS * d²f/dn_r²
        // d²C_3/dn_r² = 0.5*w_PPS * d²f/dn_r²
        // Apply FACT=2 for L>=2
        double d2E_dnr2 = 2.0 * 0.5 * w_PPS_ * d2f_dnr2 * (E_micro_[3] - E_micro_[2]);

        // Newton-Raphson step (or gradient descent if Hessian ≈ 0)
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
                           iter, n_r_, n_s_, f,
                           C_L_current[0], C_L_current[1], C_L_current[2], C_L_current[3],
                           sum_CL, E_SA_current, dE_dnr, d2E_dnr2);
        }

        // Update n_r with bounds [0.0, 2.0] (same as GAMESS x2 ∈ [0, 1])
        double n_r_new = n_r_ + delta;
        n_r_new = std::clamp(n_r_new, 0.0, 2.0);
        delta = n_r_new - n_r_;
        n_r_ = n_r_new;

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

    n_s_ = 2.0 - n_r_;

    // Compute SA-REKS energy with optimized FON
    compute_weighting_factors();  // Update C_L with new n_r

    // CRITICAL FIX: GAMESS REXEM2EE multiplies by 2 for L>=3 (1-based)
    // In Psi4 (0-based): multiply by 2 for L>=2
    // GAMESS code: IF(I.GE.3) FACT=TWO; REXEM2EE = REXEM2EE + FACT*CM(I)*EM(I)
    double E_SA = C_L_[0] * E_micro_[0] + C_L_[1] * E_micro_[1]
                + 2.0 * C_L_[2] * E_micro_[2] + 2.0 * C_L_[3] * E_micro_[3];

    if (reks_debug_ >= 1) {
        outfile->Printf("  Final: n_r = %12.8f, n_s = %12.8f, f(x) = %12.8f\n",
                        n_r_, n_s_, f_interp(n_r_ / 2.0));
        outfile->Printf("  E_SA = %20.12f\n", E_SA);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Step 10: FockMicro2Macro - Coupling Operator Assembly (Filatov 2024, Eq.11-15)
// ═══════════════════════════════════════════════════════════════════════════

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

    double** Freks = F_reks_->pointer(0);

    int nra = nr_alpha_[L];
    int nrb = nr_beta_[L];
    int nsa = ns_alpha_[L];
    int nsb = ns_beta_[L];

    // GAMESS formula (reks.src line 1178-1179): FR = DNR*WPPS + 0.5*WOSS
    // where DNR = n_r/2, DNS = n_s/2 (normalized to [0,1])
    // So: FR = (n_r/2)*WPPS + 0.5*WOSS
    double fr = (n_r_ / 2.0) * w_PPS_ + 0.5 * w_OSS_;
    double fs = (n_s_ / 2.0) * w_PPS_ + 0.5 * w_OSS_;

    double Wc = 0.5 * Cl;  // Core block weight

    // === Core-core block (i,j < Ncore) ===
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = i; j < Ncore_; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // === Virtual-virtual block (i,j > active_s) ===
    for (int i = active_s_ + 1; i < nsopi_[0]; ++i) {
        for (int j = i; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // === Core-virtual block ===
    for (int i = 0; i < Ncore_; ++i) {
        for (int j = active_s_ + 1; j < nsopi_[0]; ++j) {
            Freks[i][j] += Wc * (Fa[i][j] + Fb[i][j]);
        }
    }

    // === Active orbital r diagonal ===
    // Weight: 0.5 * C_L * n^σ_r / f_r
    double Wr_a = (fr > 1e-10) ? 0.5 * Cl * nra / fr : 0.0;
    double Wr_b = (fr > 1e-10) ? 0.5 * Cl * nrb / fr : 0.0;
    Freks[active_r_][active_r_] += Wr_a * Fa[active_r_][active_r_] + Wr_b * Fb[active_r_][active_r_];

    // === Active orbital s diagonal ===
    // Weight: 0.5 * C_L * n^σ_s / f_s
    double Ws_a = (fs > 1e-10) ? 0.5 * Cl * nsa / fs : 0.0;
    double Ws_b = (fs > 1e-10) ? 0.5 * Cl * nsb / fs : 0.0;
    Freks[active_s_][active_s_] += Ws_a * Fa[active_s_][active_s_] + Ws_b * Fb[active_s_][active_s_];

    // === Core-active coupling (core to r) ===
    // Weight: 0.5 * C_L * (1 - n^σ_r) / (1 - f_r)
    double omfr = 1.0 - fr;  // 1 - f_r
    double Wcr_a = (omfr > 1e-10) ? 0.5 * Cl * (1 - nra) / omfr : 0.0;
    double Wcr_b = (omfr > 1e-10) ? 0.5 * Cl * (1 - nrb) / omfr : 0.0;

    // === Core-active coupling (core to s) ===
    double omfs = 1.0 - fs;  // 1 - f_s
    double Wcs_a = (omfs > 1e-10) ? 0.5 * Cl * (1 - nsa) / omfs : 0.0;
    double Wcs_b = (omfs > 1e-10) ? 0.5 * Cl * (1 - nsb) / omfs : 0.0;

    for (int c = 0; c < Ncore_; ++c) {
        Freks[c][active_r_] += Wcr_a * Fa[c][active_r_] + Wcr_b * Fb[c][active_r_];
        Freks[c][active_s_] += Wcs_a * Fa[c][active_s_] + Wcs_b * Fb[c][active_s_];
    }

    // === Active-virtual coupling ===
    for (int v = active_s_ + 1; v < nsopi_[0]; ++v) {
        Freks[active_r_][v] += Wr_a * Fa[active_r_][v] + Wr_b * Fb[active_r_][v];
        Freks[active_s_][v] += Ws_a * Fa[active_s_][v] + Ws_b * Fb[active_s_][v];
    }

    // === Active-active coupling (r,s) ===
    // Weight: C_L * (n^σ_r - n^σ_s) * sign(f_r - f_s)
    int signrs = (fr > fs) ? 1 : ((fr < fs) ? -1 : 0);
    double Wrs_a = Cl * (nra - nsa) * signrs;
    double Wrs_b = Cl * (nrb - nsb) * signrs;
    Freks[active_r_][active_s_] += Wrs_a * Fa[active_r_][active_s_] + Wrs_b * Fb[active_r_][active_s_];

    // Accumulate Lagrange multiplier for SI (state interaction)
    Wrs_lagr_ += Cl * (nra * Fa[active_r_][active_s_] + nrb * Fb[active_r_][active_s_]);

    // Note: Symmetrization will be done ONCE in form_F() after all microstates are accumulated
}

// ═══════════════════════════════════════════════════════════════════════════
// SCF Method Overrides
// ═══════════════════════════════════════════════════════════════════════════

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

        // Overlap <r|s> = C_r^T · S · C_s
        int nso = nsopi_[0];
        double** Cp = Ca_->pointer(0);
        double** Sp = S_->pointer(0);
        double overlap_rs = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                overlap_rs += Cp[mu][active_r_] * Sp[mu][nu] * Cp[nu][active_s_];
            }
        }
        outfile->Printf("  <r|s> = %12.8f (should be ~0)\n", overlap_rs);

        // Norms <r|r>, <s|s>
        double norm_r = 0.0, norm_s = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                norm_r += Cp[mu][active_r_] * Sp[mu][nu] * Cp[nu][active_r_];
                norm_s += Cp[mu][active_s_] * Sp[mu][nu] * Cp[nu][active_s_];
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
    // Then set Fa_ = F_reks_ for orbital optimization via diagonalization.
    // Formula from Filatov 2024, pages 7-9: FockMicro2Macro function.

    // Zero F_reks_ and Lagrange multiplier
    F_reks_->zero();
    Wrs_lagr_ = 0.0;

    int N = nsopi_[0];

    // Assemble F_reks_ from all 4 microstates
    // Note: L=2 (open-shell) represents L=3 and L=4 in paper → factor 2 already in C_L_[2]
    //       L=3 (triplet-like) represents L=5 and L=6 → factor 2 already in C_L_[3]
    for (int L = 0; L < N_MICRO_; ++L) {
        double** Fa = F_alpha_micro_[L]->pointer(0);
        double** Fb = F_beta_micro_[L]->pointer(0);

        // CRITICAL FIX: GAMESS multiplies C_L by 2 for L>=3 (1-based)
        // In Psi4 (0-based): multiply by 2 for L>=2
        // GAMESS code: IF(L.GE.3) CMTMP = TWO*CMTMP
        double Cl_scaled = (L >= 2) ? 2.0 * C_L_[L] : C_L_[L];

        fock_micro_to_macro(L, Cl_scaled, Fa, Fb);
    }

    // Symmetrize F_reks_ (copy upper triangle to lower triangle)
    // This must be done ONCE after all microstates are accumulated
    double** Freks = F_reks_->pointer(0);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Freks[j][i] = Freks[i][j];
        }
    }

    // Set Fa_ = F_reks_ (F_reks_ IS the full Fock matrix, already contains H)
    Fa_->copy(F_reks_);

    // DIIS still uses Fa_
    Fb_->copy(Fa_);

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS Coupling Fock Matrix ===\n");
        outfile->Printf("  tr(F_reks) = %15.10f\n", F_reks_->trace());
        outfile->Printf("  Wrs_lagr   = %15.10f\n", Wrs_lagr_);

        // Check active orbital elements
        double** Fr = Freks;
        outfile->Printf("  F_reks[r,r] = %12.6f, F_reks[s,s] = %12.6f\n",
                       Fr[active_r_][active_r_], Fr[active_s_][active_s_]);
        outfile->Printf("  F_reks[r,s] = %12.6f\n", Fr[active_r_][active_s_]);

        // Check symmetry
        double max_asym = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double asym = std::abs(Fr[i][j] - Fr[j][i]);
                if (asym > max_asym) max_asym = asym;
            }
        }
        outfile->Printf("  Max asymmetry: %.3e\n", max_asym);
    }
}

double REKS::compute_E() {
    // For SAD guess iteration: use standard RHF energy
    if (sad_ && iteration_ <= 0) {
        return RHF::compute_E();
    }

    // Compute state-averaged REKS energy: E_SA = Σ FACT * C_L * E_L
    // CRITICAL FIX: GAMESS REXEM2EE multiplies by 2 for L>=3 (1-based)
    // In Psi4 (0-based): multiply by 2 for L>=2
    // GAMESS code: IF(I.GE.3) FACT=TWO; REXEM2EE = REXEM2EE + FACT*CM(I)*EM(I)
    double E_SA = 0.0;
    for (int L = 0; L < N_MICRO_; ++L) {
        double fact = (L >= 2) ? 2.0 : 1.0;
        E_SA += fact * C_L_[L] * E_micro_[L];
    }

    energies_["Total Energy"] = E_SA;

    if (reks_debug_ >= 2) {
        outfile->Printf("\n  === REKS State-Averaged Energy ===\n");
        outfile->Printf("  E_SA = %20.12f\n", E_SA);
    }

    return E_SA;
}

// ═══════════════════════════════════════════════════════════════════════════
// Debug Output Functions
// ═══════════════════════════════════════════════════════════════════════════

void REKS::print_microstate_energies() const {
    outfile->Printf("\n  === REKS Microstate Energies ===\n");
    outfile->Printf("  E_1 (L=0, r↑r↓)     = %20.12f\n", E_micro_[0]);
    outfile->Printf("  E_2 (L=1, s↑s↓)     = %20.12f\n", E_micro_[1]);
    outfile->Printf("  E_3 (L=2, r↑s↓)     = %20.12f\n", E_micro_[2]);
    outfile->Printf("  E_5 (L=3, r↑s↑/r↓s↓)= %20.12f\n", E_micro_[3]);
}

void REKS::print_fon_info() const {
    double f = f_interp(n_r_ / 2.0);
    outfile->Printf("\n  === FON Analysis ===\n");
    outfile->Printf("  n_r = %12.8f, n_s = %12.8f\n", n_r_, n_s_);
    outfile->Printf("  n_r + n_s = %12.8f (should be 2.0)\n", n_r_ + n_s_);
    outfile->Printf("  f(n_r/2) = %12.8f\n", f);
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
