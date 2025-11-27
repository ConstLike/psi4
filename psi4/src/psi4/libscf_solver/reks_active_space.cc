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
#include <algorithm>

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
    // Initialize GVB pair at closed-shell limit: n_r=2.0, n_s=0.0
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

    // SA-REKS weight formulas for state-averaged ensemble
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

int REKS22Space::get_alpha_base_idx(int L) const {
    // Convert occupation pattern to base density index using bitmask encoding.
    // index = sum(alpha[i] * 2^i) where i is orbital index.
    // Patterns: {0,0}->0, {1,0}->1, {0,1}->2, {1,1}->3
    const auto& micro = microstates_[L];
    return micro.alpha[0] + 2 * micro.alpha[1];
}

int REKS22Space::get_beta_base_idx(int L) const {
    // Same logic as alpha
    const auto& micro = microstates_[L];
    return micro.beta[0] + 2 * micro.beta[1];
}

std::vector<double> REKS22Space::compute_energy_gradient(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // For REKS(2,2): single FON variable n_r
    // dE_SA/dn_r = sum_L { FACT_L * dC_L/dn_r * E_L }
    //
    // Where FACT_L accounts for spin symmetry:
    //   L=0,1: FACT=1.0 (closed-shell)
    //   L=2,3: FACT=2.0 (open-shell with spin partner)

    std::vector<double> dC_dfon(4), d2C_dfon2(4);
    compute_weight_derivs(dC_dfon, d2C_dfon2);

    double gradient = 0.0;
    for (int L = 0; L < 4; ++L) {
        double FACT = (L >= 2) ? 2.0 : 1.0;
        gradient += FACT * dC_dfon[L] * E_micro[L];
    }

    return {gradient};
}

std::vector<double> REKS22Space::compute_energy_hessian(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // For REKS(2,2): single FON variable n_r
    // d²E_SA/dn_r² = sum_L { FACT_L * d²C_L/dn_r² * E_L }

    std::vector<double> dC_dfon(4), d2C_dfon2(4);
    compute_weight_derivs(dC_dfon, d2C_dfon2);

    double hessian = 0.0;
    for (int L = 0; L < 4; ++L) {
        double FACT = (L >= 2) ? 2.0 : 1.0;
        hessian += FACT * d2C_dfon2[L] * E_micro[L];
    }

    return {hessian};
}

double REKS22Space::get_effective_fon(int orbital_idx) const {
    // Effective FON for Fock matrix coupling: f_i = (n_i/2)*w_PPS + 0.5*w_OSS
    double n_i = (orbital_idx == 0) ? pair_.fon_p : pair_.fon_q;
    return (n_i / 2.0) * w_pps_ + 0.5 * w_oss_;
}

REKSActiveSpace::SIResult REKS22Space::compute_SI_energies(
    const std::vector<double>& E_micro,
    double Wrs,
    int n_si_states) const {

    // SI-SA-REKS: build and diagonalize state interaction Hamiltonian
    // Couples PPS, OSS, and optionally DES states

    SIResult result;
    result.n_states = n_si_states;

    double n_r = pair_.fon_p;
    double n_s = pair_.fon_q;

    // Diagonal energies computed with pure PPS weights (w_PPS=1, w_OSS=0)
    // C_L^{PPS} = {n_r/2, n_s/2, -f/2, f/2}
    double f = f_interp(n_r * n_s);
    double C_pps[4] = {n_r / 2.0, n_s / 2.0, -f / 2.0, f / 2.0};

    // E_PPS = sum_L C_L^{PPS} * FACT_L * E_L
    result.E_PPS = C_pps[0] * E_micro[0] + C_pps[1] * E_micro[1]
                 + 2.0 * C_pps[2] * E_micro[2] + 2.0 * C_pps[3] * E_micro[3];

    // E_OSS = 2*E_L[2] - E_L[3]
    result.E_OSS = 2.0 * E_micro[2] - E_micro[3];

    // E_DES: swap n_r/n_s weights and negate open-shell terms
    result.E_DES = C_pps[1] * E_micro[0] + C_pps[0] * E_micro[1]
                 - 2.0 * C_pps[2] * E_micro[2] - 2.0 * C_pps[3] * E_micro[3];

    // Off-diagonal coupling elements
    // H_12 = Wrs * (sqrt(n_r) - sqrt(n_s)) * sqrt(2)
    result.H_12 = Wrs * (std::sqrt(n_r) - std::sqrt(n_s)) * std::sqrt(2.0);

    // H_23 = Wrs * (sqrt(n_r) + sqrt(n_s)) * sqrt(2)
    result.H_23 = Wrs * (std::sqrt(n_r) + std::sqrt(n_s)) * std::sqrt(2.0);

    // Step 3: Build and diagonalize SI Hamiltonian
    if (n_si_states == 2) {
        // 2SI-2SA: 2x2 symmetric matrix
        // H = | E_PPS   H_12  |
        //     | H_12    E_OSS |

        double H11 = result.E_PPS;
        double H22 = result.E_OSS;
        double H12 = result.H_12;

        // Analytic diagonalization of 2x2 symmetric matrix
        double trace = H11 + H22;
        double det = H11 * H22 - H12 * H12;
        double disc = std::sqrt(trace * trace / 4.0 - det);

        double E0 = trace / 2.0 - disc;  // Lower eigenvalue
        double E1 = trace / 2.0 + disc;  // Upper eigenvalue

        result.energies = {E0, E1};

        // Eigenvectors: for eigenvalue E, (H - E*I) v = 0
        // v1 = (H12, E - H11) normalized, or (E - H22, H12) normalized
        //
        // For E0: v0 = (H12, E0 - H11) / norm
        double v0_pps = H12;
        double v0_oss = E0 - H11;
        double norm0 = std::sqrt(v0_pps * v0_pps + v0_oss * v0_oss);
        if (norm0 > constants::NORM_THRESHOLD) {
            v0_pps /= norm0;
            v0_oss /= norm0;
        } else {
            // Degenerate case: E0 = H11
            v0_pps = 1.0;
            v0_oss = 0.0;
        }

        // For E1: v1 = (H12, E1 - H11) / norm
        double v1_pps = H12;
        double v1_oss = E1 - H11;
        double norm1 = std::sqrt(v1_pps * v1_pps + v1_oss * v1_oss);
        if (norm1 > constants::NORM_THRESHOLD) {
            v1_pps /= norm1;
            v1_oss /= norm1;
        } else {
            v1_pps = 0.0;
            v1_oss = 1.0;
        }

        // Store coefficients: coeffs[state * n_states + component]
        // State 0 (ground): coeffs[0] = C_PPS, coeffs[1] = C_OSS
        // State 1 (excited): coeffs[2] = C_PPS, coeffs[3] = C_OSS
        result.coeffs = {v0_pps, v0_oss, v1_pps, v1_oss};

    } else if (n_si_states == 3) {
        // 3SI-2SA: 3x3 symmetric matrix
        // H = | E_PPS   H_12    0    |
        //     | H_12    E_OSS   H_23 |
        //     |   0     H_23  E_DES  |
        //
        // Use explicit 3x3 eigenvalue solver (Cardano's formula)
        // For simplicity, we use iterative Jacobi diagonalization

        double H[3][3] = {
            {result.E_PPS, result.H_12, 0.0},
            {result.H_12, result.E_OSS, result.H_23},
            {0.0, result.H_23, result.E_DES}
        };

        // Initialize eigenvectors to identity
        double V[3][3] = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        // Jacobi iteration for 3x3 symmetric matrix
        const int max_iter = constants::FON_MAX_ITER;
        const double tol = constants::NORM_THRESHOLD;

        for (int iter = 0; iter < max_iter; ++iter) {
            // Find largest off-diagonal element
            double max_off = 0.0;
            int p = 0, q = 1;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    if (std::abs(H[i][j]) > max_off) {
                        max_off = std::abs(H[i][j]);
                        p = i;
                        q = j;
                    }
                }
            }

            if (max_off < tol) break;

            // Compute rotation angle
            double theta = (H[q][q] - H[p][p]) / (2.0 * H[p][q]);
            double t = (theta >= 0) ?
                       1.0 / (theta + std::sqrt(1.0 + theta * theta)) :
                      -1.0 / (-theta + std::sqrt(1.0 + theta * theta));
            double c = 1.0 / std::sqrt(1.0 + t * t);
            double s = t * c;

            // Apply Jacobi rotation to H
            double H_pp = H[p][p];
            double H_qq = H[q][q];
            double H_pq = H[p][q];

            H[p][p] = c * c * H_pp - 2.0 * s * c * H_pq + s * s * H_qq;
            H[q][q] = s * s * H_pp + 2.0 * s * c * H_pq + c * c * H_qq;
            H[p][q] = H[q][p] = 0.0;

            for (int k = 0; k < 3; ++k) {
                if (k != p && k != q) {
                    double H_kp = H[k][p];
                    double H_kq = H[k][q];
                    H[k][p] = H[p][k] = c * H_kp - s * H_kq;
                    H[k][q] = H[q][k] = s * H_kp + c * H_kq;
                }
            }

            // Apply rotation to eigenvectors
            for (int k = 0; k < 3; ++k) {
                double V_kp = V[k][p];
                double V_kq = V[k][q];
                V[k][p] = c * V_kp - s * V_kq;
                V[k][q] = s * V_kp + c * V_kq;
            }
        }

        // Extract eigenvalues and sort them
        std::vector<std::pair<double, int>> eig_pairs = {
            {H[0][0], 0}, {H[1][1], 1}, {H[2][2], 2}
        };
        std::sort(eig_pairs.begin(), eig_pairs.end());

        result.energies.resize(3);
        result.coeffs.resize(9);
        for (int i = 0; i < 3; ++i) {
            result.energies[i] = eig_pairs[i].first;
            int col = eig_pairs[i].second;
            for (int j = 0; j < 3; ++j) {
                result.coeffs[i * 3 + j] = V[j][col];
            }
        }
    }

    return result;
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

int REKS44Space::get_alpha_base_idx(int L) const {
    // Convert alpha occupation pattern to bitmask index
    // For REKS(4,4): 4 active orbitals, pattern = sum(alpha[i] << i)
    const auto& micro = microstates_[L];
    int pattern = 0;
    for (int i = 0; i < 4; ++i) {
        pattern += micro.alpha[i] << i;
    }
    return pattern;
}

int REKS44Space::get_beta_base_idx(int L) const {
    const auto& micro = microstates_[L];
    int pattern = 0;
    for (int i = 0; i < 4; ++i) {
        pattern += micro.beta[i] << i;
    }
    return pattern;
}

double REKS44Space::get_symmetry_factor(int L) const {
    // TODO: Implement proper symmetry factors for REKS(4,4)
    // This requires analysis of microstate spin symmetry from Table IV
    // For now, return 1.0 (placeholder)
    return 1.0;
}

std::vector<double> REKS44Space::compute_energy_gradient(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // TODO: Implement multi-variable gradient for REKS(4,4)
    // Returns {dE/dn_a, dE/dn_b}
    // For now, return zeros (placeholder)
    return {0.0, 0.0};
}

std::vector<double> REKS44Space::compute_energy_hessian(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // TODO: Implement 2x2 Hessian for REKS(4,4)
    // Returns {d²E/dn_a², d²E/dn_a dn_b, d²E/dn_b dn_a, d²E/dn_b²}
    // For now, return zeros (placeholder)
    return {0.0, 0.0, 0.0, 0.0};
}

double REKS44Space::get_effective_fon(int orbital_idx) const {
    // Effective FON for Fock coupling
    // For REKS(4,4): orbitals ordered as (a, b, c, d) = (0, 1, 2, 3)
    // Pair 0: (a, d), Pair 1: (b, c)
    //
    // f_i = (n_i/2)*w_PPS + 0.5*w_OSS
    // where w_OSS = (1 - w_PPS) / 2 for REKS(4,4)

    double n_i = 0.0;
    switch (orbital_idx) {
        case 0: n_i = pairs_[0].fon_p; break;  // n_a
        case 1: n_i = pairs_[1].fon_p; break;  // n_b
        case 2: n_i = pairs_[1].fon_q; break;  // n_c
        case 3: n_i = pairs_[0].fon_q; break;  // n_d
        default: n_i = 0.0;
    }

    double w_oss = (1.0 - w_pps_) / 2.0;
    return (n_i / 2.0) * w_pps_ + 0.5 * w_oss;
}

REKSActiveSpace::SIResult REKS44Space::compute_SI_energies(
    const std::vector<double>& E_micro,
    double Wrs,
    int n_si_states) const {

    // TODO: Implement SI-SA for REKS(4,4)
    // This requires computing 5-state SI Hamiltonian
    // For now, return placeholder result

    SIResult result;
    result.n_states = n_si_states;
    result.E_PPS = 0.0;
    result.E_OSS = 0.0;
    result.E_DES = 0.0;
    result.H_12 = 0.0;
    result.H_23 = 0.0;
    result.energies.resize(n_si_states, 0.0);
    result.coeffs.resize(n_si_states * n_si_states, 0.0);

    return result;
}

}  // namespace reks
}  // namespace psi
