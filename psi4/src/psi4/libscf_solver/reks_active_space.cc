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
 * @file reks_active_space.cc
 * @brief Implementation of REKS active space classes
 *
 * Contains implementations for:
 * - REKSActiveSpace factory methods
 * - REKS22Space (2 electrons in 2 orbitals)
 * - REKS44Space placeholder (4 electrons in 4 orbitals)
 */

#include "reks_active_space.h"
#include <stdexcept>
#include <cmath>

namespace psi {
namespace reks {

// ============================================================================
// Factory Methods
// ============================================================================

std::unique_ptr<REKSActiveSpace> REKSActiveSpace::create_2_2(double w_pps, double w_oss) {
    return std::make_unique<REKS22Space>(w_pps, w_oss);
}

std::unique_ptr<REKSActiveSpace> REKSActiveSpace::create_4_4(double w_pps) {
    return std::make_unique<REKS44Space>(w_pps);
}

// ============================================================================
// REKS22Space Implementation
// ============================================================================

REKS22Space::REKS22Space(double w_pps, double w_oss)
    : w_pps_(w_pps), w_oss_(w_oss) {
    // Initialize GVB pair with closed-shell starting point
    // n_r = 2.0, n_s = 0.0 (like GAMESS default)
    pair_ = GVBPair(0, 1, 2.0);

    // Initialize microstate occupation patterns
    init_microstates();
}

void REKS22Space::init_microstates() {
    // 4 unique microstates for REKS(2,2)
    // Using spin symmetry: L=3=L=4, L=5=L=6
    //
    // Orbital indices: r=0, s=1
    //
    // From Filatov 2024, Table 2:
    //   L=0 (paper L=1): r doubly occupied  -> alpha={1,0}, beta={1,0}
    //   L=1 (paper L=2): s doubly occupied  -> alpha={0,1}, beta={0,1}
    //   L=2 (paper L=3): r-alpha, s-beta    -> alpha={1,0}, beta={0,1}
    //   L=3 (paper L=5): triplet-like       -> alpha={1,1}, beta={0,0}

    // L=0: r doubly occupied (closed-shell r)
    microstates_[0] = Microstate({1, 0}, {1, 0});

    // L=1: s doubly occupied (closed-shell s)
    microstates_[1] = Microstate({0, 1}, {0, 1});

    // L=2: r-alpha, s-beta (open-shell)
    microstates_[2] = Microstate({1, 0}, {0, 1});

    // L=3: triplet-like component (alpha on both r and s)
    microstates_[3] = Microstate({1, 1}, {0, 0});
}

void REKS22Space::compute_weights(std::vector<double>& C_L) const {
    // Resize if needed
    if (C_L.size() != 4) {
        C_L.resize(4);
    }

    // Get current FONs
    double n_r = pair_.fon_p;
    double n_s = pair_.fon_q;

    // Compute f(n_r * n_s) - the interpolating function
    double f = f_interp(n_r * n_s);

    // SA-REKS weight formulas (GAMESS REXCM function, STATAVG branch)
    //
    // L=0: C_0 = w_PPS * n_r/2
    C_L[0] = w_pps_ * n_r / 2.0;

    // L=1: C_1 = w_PPS * n_s/2
    C_L[1] = w_pps_ * n_s / 2.0;

    // L=2: C_2 = w_OSS - 0.5 * w_PPS * f(n_r*n_s)
    C_L[2] = w_oss_ - 0.5 * w_pps_ * f;

    // L=3: C_3 = 0.5 * w_PPS * f(n_r*n_s) - 0.5 * w_OSS
    C_L[3] = 0.5 * w_pps_ * f - 0.5 * w_oss_;
}

void REKS22Space::compute_weight_derivs(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // Resize if needed
    if (dC_dfon.size() != 4) {
        dC_dfon.resize(4);
    }
    if (d2C_dfon2.size() != 4) {
        d2C_dfon2.resize(4);
    }

    // Get current FONs
    double n_r = pair_.fon_p;
    double n_s = pair_.fon_q;

    // Compute f and derivatives at x = n_r * n_s
    double x = n_r * n_s;
    double f = f_interp(x);
    double df_dx = df_interp(x);
    double d2f_dx2 = d2f_interp(x);

    // Chain rule: x = n_r * n_s, dx/dn_r = n_s (since n_s = 2 - n_r)
    // df/dn_r = df/dx * dx/dn_r = df/dx * n_s
    double df_dnr = df_dx * n_s;

    // d²f/dn_r² = d/dn_r[df/dx * n_s]
    //           = d²f/dx² * (dx/dn_r)² + df/dx * d(n_s)/dn_r
    //           = d²f/dx² * n_s² + df/dx * (-1)
    // Since dn_s/dn_r = -1
    double d2f_dnr2 = d2f_dx2 * n_s * n_s - df_dx;

    // First derivatives dC_L/dn_r:
    //
    // C_0 = w_PPS * n_r/2             -> dC_0/dn_r = w_PPS/2
    dC_dfon[0] = w_pps_ / 2.0;

    // C_1 = w_PPS * n_s/2             -> dC_1/dn_r = w_PPS * (-1/2) = -w_PPS/2
    dC_dfon[1] = -w_pps_ / 2.0;

    // C_2 = w_OSS - 0.5*w_PPS*f       -> dC_2/dn_r = -0.5*w_PPS*df/dn_r
    dC_dfon[2] = -0.5 * w_pps_ * df_dnr;

    // C_3 = 0.5*w_PPS*f - 0.5*w_OSS   -> dC_3/dn_r = 0.5*w_PPS*df/dn_r
    dC_dfon[3] = 0.5 * w_pps_ * df_dnr;

    // Second derivatives d²C_L/dn_r²:
    //
    // C_0 and C_1 are linear in n_r, so second derivative = 0
    d2C_dfon2[0] = 0.0;
    d2C_dfon2[1] = 0.0;

    // C_2: d²C_2/dn_r² = -0.5*w_PPS * d²f/dn_r²
    d2C_dfon2[2] = -0.5 * w_pps_ * d2f_dnr2;

    // C_3: d²C_3/dn_r² = 0.5*w_PPS * d²f/dn_r²
    d2C_dfon2[3] = 0.5 * w_pps_ * d2f_dnr2;
}

// ============================================================================
// REKS44Space Implementation (Placeholder)
// ============================================================================

REKS44Space::REKS44Space(double w_pps) : w_pps_(w_pps) {
    // Initialize two GVB pairs
    // Pair 0: (a, d) -> orbitals 0 and 3
    // Pair 1: (b, c) -> orbitals 1 and 2
    pairs_[0] = GVBPair(0, 3, 1.0);  // n_a = n_d = 1.0
    pairs_[1] = GVBPair(1, 2, 1.0);  // n_b = n_c = 1.0

    init_microstates();
}

void REKS44Space::init_microstates() {
    // REKS(4,4) has 12 unique microstates (Table IV in Filatov 2017)
    // This is a placeholder - full implementation requires careful
    // enumeration of all configurations

    // For now, just create empty microstates
    // Full implementation will follow paper's Table IV

    microstates_.resize(12);
    for (int L = 0; L < 12; ++L) {
        microstates_[L] = Microstate(4);  // 4 active orbitals
    }

    // TODO: Populate microstates according to Filatov 2017 Table IV
}

int REKS44Space::n_microstates() const {
    return static_cast<int>(microstates_.size());
}

void REKS44Space::compute_weights(std::vector<double>& C_L) const {
    // Resize to 12 microstates
    if (C_L.size() != microstates_.size()) {
        C_L.resize(microstates_.size());
    }

    // Get FONs for both pairs
    double n_a = pairs_[0].fon_p;
    double n_d = pairs_[0].fon_q;
    double n_b = pairs_[1].fon_p;
    double n_c = pairs_[1].fon_q;

    // From Filatov 2017, Eq. 3a-3c:
    //
    // Closed-shell weights (Eq. 3a):
    // C_1 = n_a*n_b/4,  C_2 = n_a*n_c/4,  C_3 = n_b*n_d/4,  C_4 = n_c*n_d/4
    //
    // Coupling weights (Eq. 3b-3c):
    // C_5 = C_6 = -C_7 = -C_8 = -f(n_a*n_d)/2
    // C_9 = C_10 = -C_11 = -C_12 = -f(n_b*n_c)/2

    // NOTE: This is simplified - full implementation needs to match
    // exact microstate ordering from Table IV

    double f_ad = f_interp(n_a * n_d);
    double f_bc = f_interp(n_b * n_c);

    // Placeholder: set all weights to zero for now
    for (auto& c : C_L) {
        c = 0.0;
    }

    // TODO: Implement full weight computation following Filatov 2017
    // This requires careful matching of microstate indices

    // Basic closed-shell weights (indices 0-3 in our ordering)
    if (C_L.size() >= 4) {
        C_L[0] = n_a * n_b / 4.0;
        C_L[1] = n_a * n_c / 4.0;
        C_L[2] = n_b * n_d / 4.0;
        C_L[3] = n_c * n_d / 4.0;
    }

    // Coupling weights (indices 4-11)
    if (C_L.size() >= 12) {
        C_L[4] = -0.5 * f_ad;
        C_L[5] = -0.5 * f_ad;
        C_L[6] = 0.5 * f_ad;
        C_L[7] = 0.5 * f_ad;
        C_L[8] = -0.5 * f_bc;
        C_L[9] = -0.5 * f_bc;
        C_L[10] = 0.5 * f_bc;
        C_L[11] = 0.5 * f_bc;
    }
}

void REKS44Space::compute_weight_derivs(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // Placeholder implementation
    int n = n_microstates();
    if (dC_dfon.size() != static_cast<size_t>(n)) {
        dC_dfon.resize(n, 0.0);
    }
    if (d2C_dfon2.size() != static_cast<size_t>(n)) {
        d2C_dfon2.resize(n, 0.0);
    }

    // TODO: Implement full derivative computation for REKS(4,4)
    // This requires derivatives w.r.t. both FON variables
}

}  // namespace reks
}  // namespace psi
