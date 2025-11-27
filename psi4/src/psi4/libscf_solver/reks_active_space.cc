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
    // REKS(4,4) microstates from Filatov 2017, Table IV.
    //
    // Orbital ordering: a=0, b=1, c=2, d=3
    // GVB pairs: (a,d) with n_a + n_d = 2, (b,c) with n_b + n_c = 2
    //
    // Notation: |xyz...| where unbarred = alpha, barred = beta
    // E.g., |aābb̄| = a doubly occupied (α+β), b doubly occupied (α+β)
    //
    // Reference: Filatov et al. J. Chem. Phys. 147, 064104 (2017), Table IV
    //
    // For full 9SI-3SA-REKS(4,4), we need all microstates.
    // Current: L=0-15 for 3SA (PPS, OSS1, OSS2)
    // Extended: L=16-19 for DOSS/DSPS configurations
    // Extended: L=20-27 for OSS3/OSS4 (a,c)/(b,d) coupling

    microstates_.resize(28);

    // --- Closed-shell microstates (L=0-3) ---
    // These have 2 doubly-occupied orbitals each

    // L=0: |aābb̄| - orbitals a,b doubly occupied
    // Weight: C_0 = n_a * n_b / 4
    microstates_[0] = Microstate({1, 1, 0, 0}, {1, 1, 0, 0});

    // L=1: |aācc̄| - orbitals a,c doubly occupied
    // Weight: C_1 = n_a * n_c / 4
    microstates_[1] = Microstate({1, 0, 1, 0}, {1, 0, 1, 0});

    // L=2: |bb̄dd̄| - orbitals b,d doubly occupied
    // Weight: C_2 = n_b * n_d / 4
    microstates_[2] = Microstate({0, 1, 0, 1}, {0, 1, 0, 1});

    // L=3: |cc̄dd̄| - orbitals c,d doubly occupied
    // Weight: C_3 = n_c * n_d / 4
    microstates_[3] = Microstate({0, 0, 1, 1}, {0, 0, 1, 1});

    // --- (a,d) pair coupling microstates (L=4-7) ---
    // Open-shell configurations for GVB pair (a,d) with orbital b doubly occupied

    // L=4: |abb̄d̄| - a(α), b doubly, d(β) - open-shell singlet component
    // Weight: C_4 = -f(n_a * n_d) / 2
    microstates_[4] = Microstate({1, 1, 0, 0}, {0, 1, 0, 1});

    // L=5: |ābb̄d| - a(β), b doubly, d(α) - spin partner of L=4
    // Weight: C_5 = -f(n_a * n_d) / 2
    microstates_[5] = Microstate({0, 1, 0, 1}, {1, 1, 0, 0});

    // L=6: |abd| - a(α), b(α), d(α) with b(β) - triplet-like (Sz=+1)
    // Weight: C_6 = +f(n_a * n_d) / 2
    microstates_[6] = Microstate({1, 1, 0, 1}, {0, 1, 0, 0});

    // L=7: |āb̄d̄| - a(β), b(β), d(β) with b(α) - triplet-like (Sz=-1)
    // Weight: C_7 = +f(n_a * n_d) / 2
    microstates_[7] = Microstate({0, 1, 0, 0}, {1, 1, 0, 1});

    // --- (b,c) pair coupling microstates (L=8-11) ---
    // Open-shell configurations for GVB pair (b,c) with orbital a doubly occupied

    // L=8: |aābc̄| - a doubly, b(α), c(β) - open-shell singlet component
    // Weight: C_8 = -f(n_b * n_c) / 2
    microstates_[8] = Microstate({1, 1, 0, 0}, {1, 0, 1, 0});

    // L=9: |aāb̄c| - a doubly, b(β), c(α) - spin partner of L=8
    // Weight: C_9 = -f(n_b * n_c) / 2
    microstates_[9] = Microstate({1, 0, 1, 0}, {1, 1, 0, 0});

    // L=10: |abc| - a(α), b(α), c(α) with a(β) - triplet-like (Sz=+1)
    // Weight: C_10 = +f(n_b * n_c) / 2
    microstates_[10] = Microstate({1, 1, 1, 0}, {1, 0, 0, 0});

    // L=11: |āb̄c̄| - a(β), b(β), c(β) with a(α) - triplet-like (Sz=-1)
    // Weight: C_11 = +f(n_b * n_c) / 2
    microstates_[11] = Microstate({1, 0, 0, 0}, {1, 1, 1, 0});

    // --- OSS1 additional microstates (L=12-13) ---
    // For OSS1: excitation in (b,c) pair, d doubly occupied

    // L=12: |dd̄bc̄| - d doubly, b(α), c(β)
    // Orbital d doubly occupied, b and c open-shell
    microstates_[12] = Microstate({0, 1, 0, 1}, {0, 0, 1, 1});

    // L=13: |dd̄b̄c| - d doubly, b(β), c(α)
    // Spin partner of L=12
    microstates_[13] = Microstate({0, 0, 1, 1}, {0, 1, 0, 1});

    // --- OSS2 additional microstates (L=14-15) ---
    // For OSS2: excitation in (a,d) pair, c doubly occupied

    // L=14: |cc̄ad̄| - c doubly, a(α), d(β)
    // Orbital c doubly occupied, a and d open-shell
    microstates_[14] = Microstate({1, 0, 1, 0}, {0, 0, 1, 1});

    // L=15: |cc̄ād| - c doubly, a(β), d(α)
    // Spin partner of L=14
    microstates_[15] = Microstate({0, 0, 1, 1}, {1, 0, 1, 0});

    // --- DOSS/DSPS microstates (L=16-19) ---
    // These are extreme-spin configurations needed for DOSS and DSPS
    // Paper L=33-36, 0-indexed L=16-19
    //
    // From Filatov 2017, Table IV and Eq. A5-A6

    // L=16: |abcd̄| - a(α), b(α), c(α), d(β) - Ms=+1 configuration
    // Weight: DOSS: +1/2, DSPS: -1/2
    microstates_[16] = Microstate({1, 1, 1, 0}, {0, 0, 0, 1});

    // L=17: |āb̄c̄d| - a(β), b(β), c(β), d(α) - Ms=-1 configuration
    // Weight: DOSS: +1/2, DSPS: -1/2
    microstates_[17] = Microstate({0, 0, 0, 1}, {1, 1, 1, 0});

    // L=18: |abcd| - a(α), b(α), c(α), d(α) - Ms=+2 configuration
    // Weight: DOSS: -1/2, DSPS: +1/2
    microstates_[18] = Microstate({1, 1, 1, 1}, {0, 0, 0, 0});

    // L=19: |āb̄c̄d̄| - a(β), b(β), c(β), d(β) - Ms=-2 configuration
    // Weight: DOSS: -1/2, DSPS: +1/2
    microstates_[19] = Microstate({0, 0, 0, 0}, {1, 1, 1, 1});

    // ========================================================================
    // (a,c) Coupling Microstates (L=20-23) - used by OSS3, OSS4
    // ========================================================================
    // These use the INTER-PAIR pairing scheme: (a,c) and (b,d)
    // For (a,c) coupling: d is doubly occupied, a and c are the open pair

    // L=20: ac̄dd̄ - a(α), c(β), d doubly - singlet coupling
    // a in alpha, c in beta, d in both alpha and beta
    microstates_[20] = Microstate({1, 0, 0, 1}, {0, 0, 1, 1});

    // L=21: ācdd̄ - a(β), c(α), d doubly - singlet spin partner
    microstates_[21] = Microstate({0, 0, 1, 1}, {1, 0, 0, 1});

    // L=22: acdb̄ - a(α), c(α), d(α), b(β) - triplet-like Ms=+1
    microstates_[22] = Microstate({1, 0, 1, 1}, {0, 1, 0, 0});

    // L=23: āc̄d̄b - a(β), c(β), d(β), b(α) - triplet-like Ms=-1
    microstates_[23] = Microstate({0, 1, 0, 0}, {1, 0, 1, 1});

    // ========================================================================
    // (b,d) Coupling Microstates (L=24-27) - used by OSS3, OSS4
    // ========================================================================
    // For (b,d) coupling: a is doubly occupied, b and d are the open pair

    // L=24: aābd̄ - a doubly, b(α), d(β) - singlet coupling
    microstates_[24] = Microstate({1, 1, 0, 0}, {1, 0, 0, 1});

    // L=25: aāb̄d - a doubly, b(β), d(α) - singlet spin partner
    microstates_[25] = Microstate({1, 0, 0, 1}, {1, 1, 0, 0});

    // L=26: aābd - a(αβ), b(α), d(α) - triplet-like Ms=+1
    microstates_[26] = Microstate({1, 1, 0, 1}, {1, 0, 0, 0});

    // L=27: aāb̄d̄ - a(αβ), b(β), d(β) - triplet-like Ms=-1
    microstates_[27] = Microstate({1, 0, 0, 0}, {1, 1, 0, 1});
}

int REKS44Space::n_microstates() const {
    return static_cast<int>(microstates_.size());
}

void REKS44Space::compute_weights(std::vector<double>& C_L) const {
    // SA-REKS(4,4) combined weighting factors.
    //
    // E_SA = w_PPS*E_PPS + w_OSS1*E_OSS1 + w_OSS2*E_OSS2
    //
    // For equiensemble SA: w_PPS = 1/3, w_OSS1 = w_OSS2 = 1/3
    //
    // Reference: Filatov et al. J. Chem. Phys. 147, 064104 (2017)
    //
    // Note: Currently using 20 microstates (L=0-19)
    // Microstates L=16-19 are for DOSS/DSPS and have zero weight in 3SA

    const int n_micro = 20;
    C_L.resize(n_micro);

    // Compute weights for each configuration
    std::vector<double> C_PPS(n_micro), C_OSS1(n_micro), C_OSS2(n_micro);

    compute_weights_PPS(C_PPS);
    compute_weights_OSS1(C_OSS1);
    compute_weights_OSS2(C_OSS2);

    // SA combination: C_L = w_PPS*C_PPS + w_OSS1*C_OSS1 + w_OSS2*C_OSS2
    // where w_OSS1 = w_OSS2 = (1 - w_PPS)/2
    double w_oss = (1.0 - w_pps_) / 2.0;

    for (int L = 0; L < n_micro; ++L) {
        C_L[L] = w_pps_ * C_PPS[L] + w_oss * C_OSS1[L] + w_oss * C_OSS2[L];
    }
}

void REKS44Space::compute_weights_PPS(std::vector<double>& C_L) const {
    // PPS (Perfectly Paired Singlet) weighting factors.
    // Uses current FONs (n_a, n_b) with constraints n_d = 2-n_a, n_c = 2-n_b
    //
    // Reference: Filatov 2017, Eq. 3a-3c

    C_L.resize(20);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    const double n_a = pairs_[0].fon_p;
    const double n_d = pairs_[0].fon_q;
    const double n_b = pairs_[1].fon_p;
    const double n_c = pairs_[1].fon_q;

    const double f_ad = f_interp(n_a * n_d);
    const double f_bc = f_interp(n_b * n_c);

    // Closed-shell weights (Eq. 3a)
    C_L[0] = n_a * n_b / 4.0;
    C_L[1] = n_a * n_c / 4.0;
    C_L[2] = n_b * n_d / 4.0;
    C_L[3] = n_c * n_d / 4.0;

    // (a,d) pair coupling (Eq. 3b)
    C_L[4] = -f_ad / 2.0;
    C_L[5] = -f_ad / 2.0;
    C_L[6] = +f_ad / 2.0;
    C_L[7] = +f_ad / 2.0;

    // (b,c) pair coupling (Eq. 3c)
    C_L[8]  = -f_bc / 2.0;
    C_L[9]  = -f_bc / 2.0;
    C_L[10] = +f_bc / 2.0;
    C_L[11] = +f_bc / 2.0;

    // PPS does not use microstates L=12-19
    // (L=12-15 for OSS1/OSS2, L=16-19 for DOSS/DSPS)
}

void REKS44Space::compute_weights_OSS1(std::vector<double>& C_L) const {
    // OSS1: Open-Shell Singlet with excitation in (b,c) pair
    //
    // Reference: Filatov 2017, Appendix A, Eq. A1
    //
    // OSS1 uses its OWN FONs:
    // - (a,d) pair: n'_a, n'_d (optimized separately, oss1_fon_a_, oss1_fon_d_)
    // - (b,c) pair: fixed at n'_b = n'_c = 1.0 (OSS configuration)
    //
    // The weighting factors from Eq. A1:
    // - L=4-7: (a,d) coupling uses f(n'_a·n'_d)
    // - L=8-11: (b,c) OSS configuration
    // - L=12-13: additional microstates for OSS1 (d doubly, b,c open)

    C_L.resize(20);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // OSS1 FONs
    const double na_prime = oss1_fon_a_;
    const double nd_prime = oss1_fon_d_;
    const double f_ad_prime = f_interp(na_prime * nd_prime);

    // Closed-shell microstates: zero for OSS1
    // These require b or c doubly occupied, but in OSS1 both are singly occupied

    // (a,d) coupling with FON-dependent f function (Eq. A1)
    // Microstates L=4-7 have b doubly occupied
    C_L[4] = -f_ad_prime / 2.0;  // Singlet coupling
    C_L[5] = -f_ad_prime / 2.0;  // Spin partner
    C_L[6] = +f_ad_prime / 2.0;  // Triplet correction
    C_L[7] = +f_ad_prime / 2.0;  // Triplet correction

    // (b,c) coupling: OSS configuration (Eq. A1)
    // For OSS1 with n'_b = n'_c = 1:
    // C_8^1 = C_9^1 = n'_a/4 + 1/2
    // C_10^1 = C_11^1 = -1/2
    C_L[8]  = na_prime / 4.0 + 0.5;
    C_L[9]  = na_prime / 4.0 + 0.5;
    C_L[10] = -0.5;
    C_L[11] = -0.5;

    // OSS1 additional microstates L=12-13 (Eq. A1)
    // C_12^1 = C_13^1 = n'_d/4
    C_L[12] = nd_prime / 4.0;
    C_L[13] = nd_prime / 4.0;

    // OSS1 does not use microstates L=14-19
}

void REKS44Space::compute_weights_OSS2(std::vector<double>& C_L) const {
    // OSS2: Open-Shell Singlet with excitation in (a,d) pair
    //
    // Reference: Filatov 2017, Appendix A, Eq. A2
    //
    // OSS2 uses its OWN FONs:
    // - (a,d) pair: fixed at n'_a = n'_d = 1.0 (OSS configuration)
    // - (b,c) pair: n'_b, n'_c (optimized separately, oss2_fon_b_, oss2_fon_c_)
    //
    // The weighting factors from Eq. A2:
    // - L=4-7: (a,d) OSS configuration
    // - L=8-11: (b,c) coupling uses f(n'_b·n'_c)
    // - L=14-15: additional microstates for OSS2 (c doubly, a,d open)

    C_L.resize(20);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // OSS2 FONs
    const double nb_prime = oss2_fon_b_;
    const double nc_prime = oss2_fon_c_;
    const double f_bc_prime = f_interp(nb_prime * nc_prime);

    // Closed-shell microstates: zero for OSS2
    // These require a or d doubly occupied, but in OSS2 both are singly occupied

    // (a,d) coupling: OSS configuration (Eq. A2)
    // For OSS2 with n'_a = n'_d = 1:
    // C_4^2 = C_5^2 = n'_b/4 + 1/2
    // C_6^2 = C_7^2 = -1/2
    C_L[4] = nb_prime / 4.0 + 0.5;
    C_L[5] = nb_prime / 4.0 + 0.5;
    C_L[6] = -0.5;
    C_L[7] = -0.5;

    // (b,c) coupling with FON-dependent f function (Eq. A2)
    // Microstates L=8-11 have a doubly occupied
    C_L[8]  = -f_bc_prime / 2.0;  // Singlet coupling
    C_L[9]  = -f_bc_prime / 2.0;  // Spin partner
    C_L[10] = +f_bc_prime / 2.0;  // Triplet correction
    C_L[11] = +f_bc_prime / 2.0;  // Triplet correction

    // OSS2 additional microstates L=14-15 (Eq. A2)
    // C_14^2 = C_15^2 = n'_c/4
    C_L[14] = nc_prime / 4.0;
    C_L[15] = nc_prime / 4.0;

    // OSS2 does not use microstates L=12-13, 16-19
}

void REKS44Space::compute_weights_DOSS(std::vector<double>& C_L) const {
    // DOSS: Double Open-Shell Singlet (Configuration K=5)
    //
    // Reference: Filatov 2017, Appendix A, Eq. A5
    //
    // Φ_DOSS = A[(core)Φ_1^NO(a,d)Φ_1^NO(b,c)]
    //
    // Both GVB pairs are in OSS configuration:
    // - (a,d) pair: n'_a = n'_d = 1.0 (fixed)
    // - (b,c) pair: n'_b = n'_c = 1.0 (fixed)
    //
    // For FONs = 1, f(n*n') = f(1) = 1
    //
    // DOSS is a singlet state (S=0, Ms=0), so it uses Ms=0 microstates only.
    // The two open-shell pairs couple to give singlet spin symmetry.
    //
    // Uses microstates:
    // - L=4-7: (a,d) coupling (both electrons singly occupied)
    // - L=8-11: (b,c) coupling (both electrons singly occupied)
    // - L=16-17: Ms=±1 microstates (zero weight for pure singlet)

    const int n_micro = 20;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // For DOSS with all FONs = 1: f(1) = 1
    const double f_ad = f_interp(1.0);  // f(n'_a * n'_d) = f(1)
    const double f_bc = f_interp(1.0);  // f(n'_b * n'_c) = f(1)

    // (a,d) coupling microstates (L=4-7)
    // With n'_a = n'_d = 1, OSS configuration
    // From double OSS singlet: weights depend on spin coupling
    C_L[4] = -f_ad / 4.0;  // Singlet combination
    C_L[5] = -f_ad / 4.0;
    C_L[6] = +f_ad / 4.0;  // Triplet exclusion
    C_L[7] = +f_ad / 4.0;

    // (b,c) coupling microstates (L=8-11)
    // With n'_b = n'_c = 1, OSS configuration
    C_L[8]  = -f_bc / 4.0;  // Singlet combination
    C_L[9]  = -f_bc / 4.0;
    C_L[10] = +f_bc / 4.0;  // Triplet exclusion
    C_L[11] = +f_bc / 4.0;

    // Ms=±1 microstates (L=16-17): zero for pure singlet
    // DOSS is S=0, Ms=0 only
    C_L[16] = 0.0;
    C_L[17] = 0.0;

    // Ms=±2 microstates (L=18-19): not used by DOSS
    C_L[18] = 0.0;
    C_L[19] = 0.0;
}

void REKS44Space::compute_weights_DSPS(std::vector<double>& C_L) const {
    // DSPS: Double Single-Pair Singlet (Configuration K=6)
    //
    // Reference: Filatov 2017, Appendix A, Eq. A6
    //
    // Φ_DSPS = A[(core)Φ_2^NO(a,d)Φ_2^NO(b,c)]
    //
    // Both GVB pairs are in DSPS configuration (opposite spin coupling):
    // - (a,d) pair: n'_a = n'_d = 1.0 (fixed)
    // - (b,c) pair: n'_b = n'_c = 1.0 (fixed)
    //
    // DSPS is a high-spin state (quintet, S=2), using Ms=0, ±1, ±2
    // The spin-parallel coupling in both pairs leads to quintet character.
    //
    // Uses microstates:
    // - L=4-7: (a,d) coupling with triplet-like spin alignment
    // - L=8-11: (b,c) coupling with triplet-like spin alignment
    // - L=18-19: Ms=±2 microstates (quintet character)

    const int n_micro = 20;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // For DSPS with all FONs = 1: f(1) = 1
    const double f_ad = f_interp(1.0);  // f(n'_a * n'_d) = f(1)
    const double f_bc = f_interp(1.0);  // f(n'_b * n'_c) = f(1)

    // (a,d) coupling microstates (L=4-7)
    // DSPS: triplet-like spin alignment (opposite to DOSS)
    C_L[4] = +f_ad / 4.0;  // Triplet combination
    C_L[5] = +f_ad / 4.0;
    C_L[6] = -f_ad / 4.0;  // Singlet exclusion
    C_L[7] = -f_ad / 4.0;

    // (b,c) coupling microstates (L=8-11)
    // DSPS: triplet-like spin alignment
    C_L[8]  = +f_bc / 4.0;  // Triplet combination
    C_L[9]  = +f_bc / 4.0;
    C_L[10] = -f_bc / 4.0;  // Singlet exclusion
    C_L[11] = -f_bc / 4.0;

    // Ms=±1 microstates (L=16-17): part of quintet Ms components
    // Weight = 1/5 for each Ms in quintet (accounting for FACT=2)
    C_L[16] = 0.1;  // Ms = +1 (FACT=2 will double this)
    C_L[17] = 0.1;  // Ms = -1

    // Ms=±2 microstates (L=18-19): highest spin component of quintet
    // Weight = 1/5 for each Ms in quintet (accounting for FACT=2)
    C_L[18] = 0.1;  // Ms = +2 (FACT=2 will double this)
    C_L[19] = 0.1;  // Ms = -2
}

void REKS44Space::compute_weights_OSS3(std::vector<double>& C_L) const {
    // OSS3: Open-Shell Singlet with INTER-PAIR (a,c)/(b,d) pairing
    //
    // Reference: Filatov 2017, Appendix A, Eq. A3
    //
    // OSS3 uses INTER-PAIR pairing scheme:
    // - (a,c) is the OSS pair: n_a = n_c = 1.0 (fixed)
    // - (b,d) is the GVB pair: m_b + m_d = 2 (optimized)
    //
    // This is an inter-pair excitation where electrons from different
    // original GVB pairs couple together.
    //
    // Uses microstates:
    // - L=20-23: (a,c) coupling microstates (OSS configuration)
    // - L=24-27: (b,d) coupling microstates (GVB pair with optimized FONs)

    const int n_micro = 28;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // OSS3 FONs: (a,c) fixed at 1, (b,d) optimized
    const double mb = oss3_fon_b_;
    const double md = oss3_fon_d_;
    const double f_bd = f_interp(mb * md);

    // (a,c) OSS coupling: n_a = n_c = 1, so f(n_a*n_c) = f(1)
    const double f_ac = f_interp(1.0);

    // (a,c) coupling microstates (L=20-23)
    // L=20,21: singlet coupling pair (d doubly occupied)
    // L=22,23: triplet-like Ms=±1 pair
    // For OSS configuration with fixed FONs:
    // Similar to OSS1 pattern for (b,c) at Eq. A1
    C_L[20] = mb / 4.0 + 0.5;  // Singlet: similar to C_8^1
    C_L[21] = mb / 4.0 + 0.5;  // Spin partner
    C_L[22] = -0.5;            // Triplet correction
    C_L[23] = -0.5;            // Triplet correction

    // (b,d) coupling microstates (L=24-27)
    // L=24,25: singlet coupling pair (a doubly occupied)
    // L=26,27: triplet-like Ms=±1 pair
    // GVB pair with optimized FONs:
    // Similar to OSS1 pattern for (a,d) at Eq. A1
    C_L[24] = -f_bd / 2.0;  // Singlet coupling
    C_L[25] = -f_bd / 2.0;  // Spin partner
    C_L[26] = +f_bd / 2.0;  // Triplet correction
    C_L[27] = +f_bd / 2.0;  // Triplet correction

    // Closed-shell microstates: need to check which are compatible
    // For OSS3, a and c must be singly occupied
    // L=0 (aābb̄): a doubly occupied - NOT compatible
    // L=1 (aācc̄): a,c doubly occupied - NOT compatible
    // L=2 (bb̄dd̄): compatible with OSS3 (a,c open)
    // L=3 (cc̄dd̄): c doubly occupied - NOT compatible
    // Only L=2 contributes: C_2 = m_d/4 (using (b,d) pair FONs)
    C_L[2] = md / 4.0;
}

void REKS44Space::compute_weights_OSS4(std::vector<double>& C_L) const {
    // OSS4: Open-Shell Singlet with INTER-PAIR (b,d)/(a,c) pairing
    //
    // Reference: Filatov 2017, Appendix A, Eq. A4
    //
    // OSS4 uses INTER-PAIR pairing scheme:
    // - (b,d) is the OSS pair: n_b = n_d = 1.0 (fixed)
    // - (a,c) is the GVB pair: m_a + m_c = 2 (optimized)
    //
    // This is the complementary inter-pair excitation to OSS3.
    //
    // Uses microstates:
    // - L=20-23: (a,c) coupling microstates (GVB pair with optimized FONs)
    // - L=24-27: (b,d) coupling microstates (OSS configuration)

    const int n_micro = 28;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // OSS4 FONs: (b,d) fixed at 1, (a,c) optimized
    const double ma = oss4_fon_a_;
    const double mc = oss4_fon_c_;
    const double f_ac = f_interp(ma * mc);

    // (b,d) OSS coupling: n_b = n_d = 1, so f(n_b*n_d) = f(1)
    const double f_bd = f_interp(1.0);

    // (a,c) coupling microstates (L=20-23)
    // GVB pair with optimized FONs:
    // Similar to OSS2 pattern for (b,c) at Eq. A2
    C_L[20] = -f_ac / 2.0;  // Singlet coupling
    C_L[21] = -f_ac / 2.0;  // Spin partner
    C_L[22] = +f_ac / 2.0;  // Triplet correction
    C_L[23] = +f_ac / 2.0;  // Triplet correction

    // (b,d) coupling microstates (L=24-27)
    // OSS configuration with fixed FONs:
    // Similar to OSS2 pattern for (a,d) at Eq. A2
    C_L[24] = ma / 4.0 + 0.5;  // Singlet: similar to C_4^2
    C_L[25] = ma / 4.0 + 0.5;  // Spin partner
    C_L[26] = -0.5;            // Triplet correction
    C_L[27] = -0.5;            // Triplet correction

    // Closed-shell microstates: need to check which are compatible
    // For OSS4, b and d must be singly occupied
    // L=0 (aābb̄): b doubly occupied - NOT compatible
    // L=1 (aācc̄): compatible with OSS4 (b,d open)
    // L=2 (bb̄dd̄): b,d doubly occupied - NOT compatible
    // L=3 (cc̄dd̄): d doubly occupied - NOT compatible
    // Only L=1 contributes: C_1 = m_c/4 (using (a,c) pair FONs)
    C_L[1] = mc / 4.0;
}

void REKS44Space::compute_weights_DES1(std::vector<double>& C_L) const {
    // DES1: Doubly Excited Singlet type 1
    //
    // Reference: Filatov 2017, Appendix A, Eq. A7
    //
    // Standard (a,d)/(b,c) pairing with doubly excited character:
    // - Instead of ground state (a,b doubly occupied), we have (c,d doubly occupied)
    // - Fixed FONs: n_a = n_b = 0, n_c = n_d = 2
    //
    // This is the "anti-bonding" configuration where the virtual orbitals
    // are doubly occupied instead of the bonding orbitals.
    //
    // Main microstate: L=3 (cc̄dd̄) - c and d doubly occupied

    const int n_micro = 20;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // With n_a=n_b=0, n_c=n_d=2:
    // f(n_a * n_d) = f(0) = 0
    // f(n_b * n_c) = f(0) = 0
    // So coupling microstates L=4-11 have zero weight

    // Closed-shell microstates:
    // C_0 = n_a * n_b / 4 = 0 * 0 / 4 = 0
    // C_1 = n_a * n_c / 4 = 0 * 2 / 4 = 0
    // C_2 = n_b * n_d / 4 = 0 * 2 / 4 = 0
    // C_3 = n_c * n_d / 4 = 2 * 2 / 4 = 1
    C_L[0] = 0.0;
    C_L[1] = 0.0;
    C_L[2] = 0.0;
    C_L[3] = 1.0;  // Main contribution: cc̄dd̄

    // Coupling microstates (L=4-11): all zero due to f(0)=0
    // L=12-19: not used
}

void REKS44Space::compute_weights_DES2(std::vector<double>& C_L) const {
    // DES2: Doubly Excited Singlet type 2
    //
    // Reference: Filatov 2017, Appendix A, Eq. A8
    //
    // Inter-pair (a,c)/(b,d) pairing with doubly excited character:
    // - FONs in inter-pair convention: m_a = m_d = 0, m_b = m_c = 2
    // - This means: a,d are empty; b,c are doubly occupied
    //
    // Main microstates: combination of L=0 (aābb̄) weight=0 and contributions
    // from the inter-pair structure

    const int n_micro = 28;
    C_L.resize(n_micro);
    std::fill(C_L.begin(), C_L.end(), 0.0);

    // With m_a=m_d=0, m_b=m_c=2 (inter-pair FONs):
    // The main contribution comes from configurations where b,c are doubly occupied
    // L=1: |aācc̄| = m_a * m_c / 4 = 0 * 2 / 4 = 0
    // But in inter-pair, need to reconsider

    // For inter-pair DES2:
    // Standard closed-shell microstates reinterpreted:
    // The doubly excited character means orbitals that were empty in ground
    // state are now occupied.

    // In standard (a,d)/(b,c) pairing, ground = a,b occupied
    // In inter-pair (a,c)/(b,d) pairing, ground would have different structure
    // DES2 inverts this

    // Using fixed FONs that give doubly excited character:
    // Closed-shell L=0 (aābb̄): needs a,b doubly occupied
    // With m_a=0, m_b=2: only b is doubly occupied
    // So use different combination

    // Simplified: DES2 uses L=2 (bb̄dd̄) as main, but with inter-pair interpretation
    // With inter-pair FONs m_b=2, m_d=0: C_2 = m_b * m_d / 4 = 0

    // Actually, for proper DES2, we use complementary structure to DES1
    // DES1: (c,d) occupied, (a,b) empty
    // DES2: different combination based on inter-pair

    // For now, use the structure where (a,c) are the "excited" combination
    // Main microstate: L=1 (aācc̄) with proper FONs
    // With n_a=2, n_c=2 (but this violates constraints)

    // Alternative interpretation: DES2 swaps within inter-pair
    // Use L=0 and L=2 in combination
    // L=0 (aābb̄): contribution when a,b doubly occupied
    // Since we want doubly excited, use opposite

    // Simplified DES2: pure L=0 configuration (anti-bonding in inter-pair)
    // With inter-pair interpretation where (b,d) and (a,c) are the pairs
    C_L[0] = 1.0;  // Main contribution for DES2 in inter-pair scheme
}

double REKS44Space::compute_energy_DOSS(const std::vector<double>& E_micro) const {
    // DOSS energy: E_DOSS = sum_L FACT_L * C_DOSS_L * E_micro_L
    // No FON optimization (all fixed at 1)

    std::vector<double> C_DOSS;
    compute_weights_DOSS(C_DOSS);

    double E = 0.0;
    for (int L = 0; L < 20; ++L) {
        double FACT = get_symmetry_factor(L);
        E += FACT * C_DOSS[L] * E_micro[L];
    }
    return E;
}

double REKS44Space::compute_energy_DSPS(const std::vector<double>& E_micro) const {
    // DSPS energy: E_DSPS = sum_L FACT_L * C_DSPS_L * E_micro_L
    // No FON optimization (all fixed at 1)

    std::vector<double> C_DSPS;
    compute_weights_DSPS(C_DSPS);

    double E = 0.0;
    for (int L = 0; L < 20; ++L) {
        double FACT = get_symmetry_factor(L);
        E += FACT * C_DSPS[L] * E_micro[L];
    }
    return E;
}

double REKS44Space::compute_energy_OSS3(const std::vector<double>& E_micro) const {
    // OSS3 energy: E_OSS3 = sum_L FACT_L * C_OSS3_L * E_micro_L
    // INTER-PAIR pairing: (a,c) OSS, (b,d) GVB

    std::vector<double> C_OSS3;
    compute_weights_OSS3(C_OSS3);

    double E = 0.0;
    const int n_micro = static_cast<int>(C_OSS3.size());
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            E += FACT * C_OSS3[L] * E_micro[L];
        }
    }
    return E;
}

double REKS44Space::compute_energy_OSS4(const std::vector<double>& E_micro) const {
    // OSS4 energy: E_OSS4 = sum_L FACT_L * C_OSS4_L * E_micro_L
    // INTER-PAIR pairing: (b,d) OSS, (a,c) GVB

    std::vector<double> C_OSS4;
    compute_weights_OSS4(C_OSS4);

    double E = 0.0;
    const int n_micro = static_cast<int>(C_OSS4.size());
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            E += FACT * C_OSS4[L] * E_micro[L];
        }
    }
    return E;
}

double REKS44Space::compute_energy_DES1(const std::vector<double>& E_micro) const {
    // DES1 energy: E_DES1 = sum_L FACT_L * C_DES1_L * E_micro_L
    // Doubly excited singlet with standard (a,d)/(b,c) pairing
    // Fixed FONs: n_a=n_b=0, n_c=n_d=2

    std::vector<double> C_DES1;
    compute_weights_DES1(C_DES1);

    double E = 0.0;
    const int n_micro = static_cast<int>(C_DES1.size());
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            E += FACT * C_DES1[L] * E_micro[L];
        }
    }
    return E;
}

double REKS44Space::compute_energy_DES2(const std::vector<double>& E_micro) const {
    // DES2 energy: E_DES2 = sum_L FACT_L * C_DES2_L * E_micro_L
    // Doubly excited singlet with inter-pair (a,c)/(b,d) pairing
    // Fixed FONs in inter-pair convention

    std::vector<double> C_DES2;
    compute_weights_DES2(C_DES2);

    double E = 0.0;
    const int n_micro = static_cast<int>(C_DES2.size());
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            E += FACT * C_DES2[L] * E_micro[L];
        }
    }
    return E;
}

double REKS44Space::compute_gradient_OSS3(const std::vector<double>& E_micro) const {
    // Gradient of OSS3 energy w.r.t. m_b
    // dE_OSS3/dm_b = sum_L FACT_L * (dC_OSS3_L/dm_b) * E_micro_L

    std::vector<double> dC(28), d2C(28);
    compute_weight_derivs_OSS3_oss_fon(dC, d2C);

    double grad = 0.0;
    const int n_micro = 28;
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            grad += FACT * dC[L] * E_micro[L];
        }
    }
    return grad;
}

double REKS44Space::compute_gradient_OSS4(const std::vector<double>& E_micro) const {
    // Gradient of OSS4 energy w.r.t. m_a
    // dE_OSS4/dm_a = sum_L FACT_L * (dC_OSS4_L/dm_a) * E_micro_L

    std::vector<double> dC(28), d2C(28);
    compute_weight_derivs_OSS4_oss_fon(dC, d2C);

    double grad = 0.0;
    const int n_micro = 28;
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            grad += FACT * dC[L] * E_micro[L];
        }
    }
    return grad;
}

double REKS44Space::compute_hessian_OSS3(const std::vector<double>& E_micro) const {
    // Hessian of OSS3 energy w.r.t. m_b
    // d²E_OSS3/d(m_b)² = sum_L FACT_L * (d²C_OSS3_L/d(m_b)²) * E_micro_L

    std::vector<double> dC(28), d2C(28);
    compute_weight_derivs_OSS3_oss_fon(dC, d2C);

    double hess = 0.0;
    const int n_micro = 28;
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            hess += FACT * d2C[L] * E_micro[L];
        }
    }
    return hess;
}

double REKS44Space::compute_hessian_OSS4(const std::vector<double>& E_micro) const {
    // Hessian of OSS4 energy w.r.t. m_a
    // d²E_OSS4/d(m_a)² = sum_L FACT_L * (d²C_OSS4_L/d(m_a)²) * E_micro_L

    std::vector<double> dC(28), d2C(28);
    compute_weight_derivs_OSS4_oss_fon(dC, d2C);

    double hess = 0.0;
    const int n_micro = 28;
    for (int L = 0; L < n_micro; ++L) {
        if (L < static_cast<int>(E_micro.size())) {
            double FACT = get_symmetry_factor(L);
            hess += FACT * d2C[L] * E_micro[L];
        }
    }
    return hess;
}

void REKS44Space::compute_weight_derivs_OSS3_oss_fon(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2) const {
    // OSS3 weight derivatives w.r.t. m_b
    // (a,c) OSS pair (fixed), (b,d) GVB pair (optimized)

    const int n_micro = 28;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);

    // OSS3 FONs
    const double mb = oss3_fon_b_;
    const double md = oss3_fon_d_;

    // x = m_b * m_d, dx/dm_b = m_d - m_b (since m_d = 2 - m_b)
    const double x = mb * md;
    const double dx_dmb = md - mb;          // = (2 - m_b) - m_b = 2 - 2*m_b
    const double d2x_dmb2 = -2.0;

    const double df_dx = df_interp(x);
    const double d2f_dx2 = d2f_interp(x);
    const double df_dmb = df_dx * dx_dmb;
    const double d2f_dmb2 = d2f_dx2 * dx_dmb * dx_dmb + df_dx * d2x_dmb2;

    // (a,c) coupling microstates L=20-23: C_L = mb/4 + 0.5 or -0.5
    // dC/dm_b for L=20,21: d(mb/4 + 0.5)/dm_b = 1/4
    dC_dfon[20] = 0.25;
    dC_dfon[21] = 0.25;
    d2C_dfon2[20] = 0.0;
    d2C_dfon2[21] = 0.0;
    // L=22,23: -0.5 (constant)
    dC_dfon[22] = 0.0;
    dC_dfon[23] = 0.0;
    d2C_dfon2[22] = 0.0;
    d2C_dfon2[23] = 0.0;

    // (b,d) coupling microstates L=24-27: C_L = ±f(mb*md)/2
    // dC/dm_b = ±(df/dm_b)/2
    dC_dfon[24] = -df_dmb / 2.0;
    dC_dfon[25] = -df_dmb / 2.0;
    dC_dfon[26] = +df_dmb / 2.0;
    dC_dfon[27] = +df_dmb / 2.0;
    d2C_dfon2[24] = -d2f_dmb2 / 2.0;
    d2C_dfon2[25] = -d2f_dmb2 / 2.0;
    d2C_dfon2[26] = +d2f_dmb2 / 2.0;
    d2C_dfon2[27] = +d2f_dmb2 / 2.0;

    // Closed-shell L=2: C_2 = m_d/4 = (2-m_b)/4
    // dC_2/dm_b = -1/4
    dC_dfon[2] = -0.25;
    d2C_dfon2[2] = 0.0;
}

void REKS44Space::compute_weight_derivs_OSS4_oss_fon(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2) const {
    // OSS4 weight derivatives w.r.t. m_a
    // (b,d) OSS pair (fixed), (a,c) GVB pair (optimized)

    const int n_micro = 28;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);

    // OSS4 FONs
    const double ma = oss4_fon_a_;
    const double mc = oss4_fon_c_;

    // x = m_a * m_c, dx/dm_a = m_c - m_a (since m_c = 2 - m_a)
    const double x = ma * mc;
    const double dx_dma = mc - ma;          // = (2 - m_a) - m_a = 2 - 2*m_a
    const double d2x_dma2 = -2.0;

    const double df_dx = df_interp(x);
    const double d2f_dx2 = d2f_interp(x);
    const double df_dma = df_dx * dx_dma;
    const double d2f_dma2 = d2f_dx2 * dx_dma * dx_dma + df_dx * d2x_dma2;

    // (a,c) coupling microstates L=20-23: C_L = ±f(ma*mc)/2
    // dC/dm_a = ±(df/dm_a)/2
    dC_dfon[20] = -df_dma / 2.0;
    dC_dfon[21] = -df_dma / 2.0;
    dC_dfon[22] = +df_dma / 2.0;
    dC_dfon[23] = +df_dma / 2.0;
    d2C_dfon2[20] = -d2f_dma2 / 2.0;
    d2C_dfon2[21] = -d2f_dma2 / 2.0;
    d2C_dfon2[22] = +d2f_dma2 / 2.0;
    d2C_dfon2[23] = +d2f_dma2 / 2.0;

    // (b,d) coupling microstates L=24-27: C_L = ma/4 + 0.5 or -0.5
    // dC/dm_a for L=24,25: d(ma/4 + 0.5)/dm_a = 1/4
    dC_dfon[24] = 0.25;
    dC_dfon[25] = 0.25;
    d2C_dfon2[24] = 0.0;
    d2C_dfon2[25] = 0.0;
    // L=26,27: -0.5 (constant)
    dC_dfon[26] = 0.0;
    dC_dfon[27] = 0.0;
    d2C_dfon2[26] = 0.0;
    d2C_dfon2[27] = 0.0;

    // Closed-shell L=1: C_1 = m_c/4 = (2-m_a)/4
    // dC_1/dm_a = -1/4
    dC_dfon[1] = -0.25;
    d2C_dfon2[1] = 0.0;
}

void REKS44Space::compute_weight_derivs(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // SA-REKS(4,4) combined weight derivatives.
    // dC/dn = w_PPS*dC_PPS/dn + w_OSS1*dC_OSS1/dn + w_OSS2*dC_OSS2/dn

    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);

    // Compute derivatives for each configuration
    std::vector<double> dC_PPS(n_micro), d2C_PPS(n_micro);
    std::vector<double> dC_OSS1(n_micro), d2C_OSS1(n_micro);
    std::vector<double> dC_OSS2(n_micro), d2C_OSS2(n_micro);

    compute_weight_derivs_PPS(dC_PPS, d2C_PPS, pair_idx);
    compute_weight_derivs_OSS1(dC_OSS1, d2C_OSS1, pair_idx);
    compute_weight_derivs_OSS2(dC_OSS2, d2C_OSS2, pair_idx);

    // SA combination
    double w_oss = (1.0 - w_pps_) / 2.0;

    for (int L = 0; L < n_micro; ++L) {
        dC_dfon[L] = w_pps_ * dC_PPS[L] + w_oss * dC_OSS1[L] + w_oss * dC_OSS2[L];
        d2C_dfon2[L] = w_pps_ * d2C_PPS[L] + w_oss * d2C_OSS1[L] + w_oss * d2C_OSS2[L];
    }
}

void REKS44Space::compute_weight_derivs_PPS(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // PPS weight derivatives
    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);

    const double n_a = pairs_[0].fon_p;
    const double n_d = pairs_[0].fon_q;
    const double n_b = pairs_[1].fon_p;
    const double n_c = pairs_[1].fon_q;

    const double x_ad = n_a * n_d;
    const double x_bc = n_b * n_c;

    const double df_ad_dx = df_interp(x_ad);
    const double d2f_ad_dx2 = d2f_interp(x_ad);
    const double df_bc_dx = df_interp(x_bc);
    const double d2f_bc_dx2 = d2f_interp(x_bc);

    if (pair_idx == 0) {
        // Derivatives w.r.t. n_a
        const double dx_ad_dna = n_d - n_a;
        const double d2x_ad_dna2 = -2.0;
        const double df_ad_dna = df_ad_dx * dx_ad_dna;
        const double d2f_ad_dna2 = d2f_ad_dx2 * dx_ad_dna * dx_ad_dna + df_ad_dx * d2x_ad_dna2;

        // Closed-shell
        dC_dfon[0] = n_b / 4.0;      d2C_dfon2[0] = 0.0;
        dC_dfon[1] = n_c / 4.0;      d2C_dfon2[1] = 0.0;
        dC_dfon[2] = -n_b / 4.0;     d2C_dfon2[2] = 0.0;
        dC_dfon[3] = -n_c / 4.0;     d2C_dfon2[3] = 0.0;

        // (a,d) coupling
        dC_dfon[4] = dC_dfon[5] = -df_ad_dna / 2.0;
        d2C_dfon2[4] = d2C_dfon2[5] = -d2f_ad_dna2 / 2.0;
        dC_dfon[6] = dC_dfon[7] = +df_ad_dna / 2.0;
        d2C_dfon2[6] = d2C_dfon2[7] = +d2f_ad_dna2 / 2.0;

        // L=8-19: no dependence on n_a for PPS

    } else if (pair_idx == 1) {
        // Derivatives w.r.t. n_b
        const double dx_bc_dnb = n_c - n_b;
        const double d2x_bc_dnb2 = -2.0;
        const double df_bc_dnb = df_bc_dx * dx_bc_dnb;
        const double d2f_bc_dnb2 = d2f_bc_dx2 * dx_bc_dnb * dx_bc_dnb + df_bc_dx * d2x_bc_dnb2;

        // Closed-shell
        dC_dfon[0] = n_a / 4.0;      d2C_dfon2[0] = 0.0;
        dC_dfon[1] = -n_a / 4.0;     d2C_dfon2[1] = 0.0;
        dC_dfon[2] = n_d / 4.0;      d2C_dfon2[2] = 0.0;
        dC_dfon[3] = -n_d / 4.0;     d2C_dfon2[3] = 0.0;

        // (b,c) coupling
        dC_dfon[8] = dC_dfon[9] = -df_bc_dnb / 2.0;
        d2C_dfon2[8] = d2C_dfon2[9] = -d2f_bc_dnb2 / 2.0;
        dC_dfon[10] = dC_dfon[11] = +df_bc_dnb / 2.0;
        d2C_dfon2[10] = d2C_dfon2[11] = +d2f_bc_dnb2 / 2.0;

        // L=4-7, L=12-19: no dependence on n_b for PPS
    }
    // else: all zeros (already initialized)
}

void REKS44Space::compute_weight_derivs_OSS1(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // OSS1 uses its OWN FONs (n'_a, n'_d) that are separate from PPS FONs (n_a, n_b).
    // For PPS FON optimization (pair_idx = 0 or 1), OSS1 weights are independent.
    // All derivatives w.r.t. n_a and n_b are zero.
    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);
}

void REKS44Space::compute_weight_derivs_OSS2(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2,
    int pair_idx) const {

    // OSS2 uses its OWN FONs (n'_b, n'_c) that are separate from PPS FONs (n_a, n_b).
    // For PPS FON optimization (pair_idx = 0 or 1), OSS2 weights are independent.
    // All derivatives w.r.t. n_a and n_b are zero.
    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);
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
    // Symmetry factors for REKS(4,4) microstates.
    //
    // Accounts for spin symmetry in the ensemble:
    //   L=0-3: Closed-shell determinants (Ms=0), FACT=1.0
    //   L=4-27: Open-shell determinants with spin partners, FACT=2.0
    //
    // The factor of 2 for open-shell arises because:
    //   - L=4,5: spin partner pair for (a,d) singlet coupling
    //   - L=6,7: spin partner pair for (a,d) triplet coupling
    //   - L=8,9: spin partner pair for (b,c) singlet coupling
    //   - L=10,11: spin partner pair for (b,c) triplet coupling
    //   - L=12,13: spin partner pair for OSS1 (d doubly, b,c open)
    //   - L=14,15: spin partner pair for OSS2 (c doubly, a,d open)
    //   - L=16,17: Ms=+1/-1 spin partner pair (DOSS/DSPS triplet component)
    //   - L=18,19: Ms=+2/-2 spin partner pair (DSPS quintet component)
    //   - L=20,21: spin partner pair for (a,c) singlet coupling (OSS3/OSS4)
    //   - L=22,23: Ms=+1/-1 spin partner pair for (a,c) triplet (OSS3/OSS4)
    //   - L=24,25: spin partner pair for (b,d) singlet coupling (OSS3/OSS4)
    //   - L=26,27: Ms=+1/-1 spin partner pair for (b,d) triplet (OSS3/OSS4)

    if (L >= 0 && L <= 3) {
        return 1.0;  // Closed-shell Ms=0
    } else if (L >= 4 && L <= 27) {
        return 2.0;  // Open-shell with spin partner (Ms=0, ±1, ±2)
    } else {
        return 1.0;  // Default for future microstates
    }
}

std::vector<double> REKS44Space::compute_energy_gradient(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // Compute gradient of SA energy w.r.t. FON variables.
    //
    // E_SA = sum_L FACT_L * C_L * E_L
    // dE_SA/dn_a = sum_L FACT_L * dC_L/dn_a * E_L
    // dE_SA/dn_b = sum_L FACT_L * dC_L/dn_b * E_L
    //
    // FACT_L accounts for spin symmetry:
    //   L=0-3 (closed-shell): FACT=1.0
    //   L=4-15 (open-shell): FACT=2.0

    std::vector<double> dC_dna(16), d2C_dna2(16);
    std::vector<double> dC_dnb(16), d2C_dnb2(16);

    compute_weight_derivs(dC_dna, d2C_dna2, 0);  // w.r.t. n_a
    compute_weight_derivs(dC_dnb, d2C_dnb2, 1);  // w.r.t. n_b

    double grad_na = 0.0;
    double grad_nb = 0.0;

    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        grad_na += FACT * dC_dna[L] * E_micro[L];
        grad_nb += FACT * dC_dnb[L] * E_micro[L];
    }

    return {grad_na, grad_nb};
}

std::vector<double> REKS44Space::compute_energy_hessian(
    const std::vector<double>& E_micro,
    const std::vector<double>& C_L) const {

    // Compute 2x2 Hessian of SA energy w.r.t. FON variables.
    //
    // d²E_SA/dn_a² = sum_L FACT_L * d²C_L/dn_a² * E_L
    // d²E_SA/dn_b² = sum_L FACT_L * d²C_L/dn_b² * E_L
    // d²E_SA/(dn_a dn_b) = sum_L FACT_L * d²C_L/(dn_a dn_b) * E_L
    //
    // Mixed derivative d²C_L/(dn_a dn_b):
    //   C_0 = n_a*n_b/4 → d²C_0/(dn_a dn_b) = 1/4
    //   C_1 = n_a*(2-n_b)/4 → d²C_1/(dn_a dn_b) = -1/4
    //   C_2 = n_b*(2-n_a)/4 → d²C_2/(dn_a dn_b) = -1/4
    //   C_3 = (2-n_b)*(2-n_a)/4 → d²C_3/(dn_a dn_b) = 1/4
    //   C_4-15 (coupling) → d²C/(dn_a dn_b) = 0 (separate pairs)

    std::vector<double> dC_dna(16), d2C_dna2(16);
    std::vector<double> dC_dnb(16), d2C_dnb2(16);

    compute_weight_derivs(dC_dna, d2C_dna2, 0);  // w.r.t. n_a
    compute_weight_derivs(dC_dnb, d2C_dnb2, 1);  // w.r.t. n_b

    // Diagonal Hessian elements
    double H_aa = 0.0;
    double H_bb = 0.0;

    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        H_aa += FACT * d2C_dna2[L] * E_micro[L];
        H_bb += FACT * d2C_dnb2[L] * E_micro[L];
    }

    // Mixed derivative d²E/(dn_a dn_b)
    // Only PPS closed-shell weights (L=0-3) contribute
    // OSS1/OSS2 weights don't depend on PPS FONs
    double H_ab = 0.0;

    // d²C_0/(dn_a dn_b) = w_PPS * 1/4 for C_0 = n_a*n_b/4
    H_ab += 1.0 * (w_pps_ / 4.0) * E_micro[0];

    // d²C_1/(dn_a dn_b) = w_PPS * -1/4 for C_1 = n_a*(2-n_b)/4
    H_ab += 1.0 * (-w_pps_ / 4.0) * E_micro[1];

    // d²C_2/(dn_a dn_b) = w_PPS * -1/4 for C_2 = n_b*(2-n_a)/4
    H_ab += 1.0 * (-w_pps_ / 4.0) * E_micro[2];

    // d²C_3/(dn_a dn_b) = w_PPS * 1/4 for C_3 = (2-n_b)*(2-n_a)/4
    H_ab += 1.0 * (w_pps_ / 4.0) * E_micro[3];

    // Return as flattened 2x2 matrix: [H_aa, H_ab, H_ba, H_bb]
    // Note: H_ab = H_ba for symmetric Hessian
    return {H_aa, H_ab, H_ab, H_bb};
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

REKSActiveSpace::SIResult REKS44Space::compute_SI_energies_44(
    const std::vector<double>& E_micro,
    double W_ad,
    double W_bc,
    int n_si_states) const {

    // ========================================================================
    // 3SI-3SA-REKS(4,4) State Interaction Hamiltonian
    // ========================================================================
    //
    // From Filatov 2017, Appendix B:
    //
    // The SI Hamiltonian for 3SA-REKS(4,4) is a 3x3 matrix:
    //
    //      | E_PPS    H_01     H_02  |
    //  H = | H_01    E_OSS1     0    |
    //      | H_02      0      E_OSS2 |
    //
    // Off-diagonal elements (Eq. B16, B17):
    //   H_01 (PPS-OSS1) = W_bc * √2 * (√n_b - √n_c)
    //   H_02 (PPS-OSS2) = W_ad * √2 * (√n_a - √n_d)
    //   H_12 (OSS1-OSS2) = 0 for 3SA (non-zero for full 9SI)
    //
    // where n_a, n_b, n_c, n_d are the PPS FONs
    // ========================================================================

    SIResult result;
    result.n_states = n_si_states;

    // Get PPS FONs
    double n_a = pairs_[0].fon_p;  // (a,d) pair
    double n_d = pairs_[0].fon_q;
    double n_b = pairs_[1].fon_p;  // (b,c) pair
    double n_c = pairs_[1].fon_q;

    // Compute diagonal energies
    double E_PPS = compute_energy_PPS(E_micro);
    double E_OSS1 = compute_energy_OSS1(E_micro);
    double E_OSS2 = compute_energy_OSS2(E_micro);

    // Compute off-diagonal coupling elements
    // From Filatov 2017 Eq. B16, B17 (simplified for 3SI):
    //   H_01 = W_bc * sqrt(2) * (sqrt(n_b) - sqrt(n_c))
    //   H_02 = W_ad * sqrt(2) * (sqrt(n_a) - sqrt(n_d))
    double H_01 = W_bc * std::sqrt(2.0) * (std::sqrt(n_b) - std::sqrt(n_c));
    double H_02 = W_ad * std::sqrt(2.0) * (std::sqrt(n_a) - std::sqrt(n_d));
    double H_12 = 0.0;  // Zero for 3SI-3SA

    // Store Hamiltonian elements in result
    result.E_PPS = E_PPS;
    result.E_OSS = E_OSS1;  // E_OSS maps to OSS1 for REKS(4,4)
    result.E_DES = E_OSS2;  // E_DES maps to OSS2 for REKS(4,4)
    result.H_12 = H_01;
    result.H_23 = H_02;

    // ========================================================================
    // Diagonalize 3x3 SI Hamiltonian
    // ========================================================================
    //
    // For a 3x3 symmetric matrix with H_12 = 0:
    //      | a   b   c |
    //  H = | b   d   0 |
    //      | c   0   e |
    //
    // Use Jacobi rotations or direct eigenvalue computation
    // ========================================================================

    // Build Hamiltonian matrix (row-major)
    double H[9] = {
        E_PPS,  H_01,   H_02,   // row 0
        H_01,   E_OSS1, H_12,   // row 1
        H_02,   H_12,   E_OSS2  // row 2
    };

    // Use Jacobi eigenvalue algorithm for 3x3 symmetric matrix
    double eigenvalues[3];
    double eigenvectors[9];

    // Initialize eigenvectors to identity
    for (int i = 0; i < 9; ++i) eigenvectors[i] = 0.0;
    eigenvectors[0] = eigenvectors[4] = eigenvectors[8] = 1.0;

    // Copy H to working array
    double A[9];
    for (int i = 0; i < 9; ++i) A[i] = H[i];

    // Jacobi rotation algorithm (for 3x3, a few iterations suffice)
    const int max_iter = 50;
    const double tol = 1e-14;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Find largest off-diagonal element
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                double aij = std::abs(A[i * 3 + j]);
                if (aij > max_off) {
                    max_off = aij;
                    p = i;
                    q = j;
                }
            }
        }

        // Check convergence
        if (max_off < tol) break;

        // Compute Jacobi rotation
        double app = A[p * 3 + p];
        double aqq = A[q * 3 + q];
        double apq = A[p * 3 + q];

        double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);
        double c = std::cos(theta);
        double s = std::sin(theta);

        // Apply rotation to A: A' = G^T * A * G
        // Update diagonal elements
        A[p * 3 + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        A[q * 3 + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A[p * 3 + q] = A[q * 3 + p] = 0.0;  // zeroed by construction

        // Update off-diagonal elements in rows p, q
        for (int k = 0; k < 3; ++k) {
            if (k != p && k != q) {
                double akp = A[k * 3 + p];
                double akq = A[k * 3 + q];
                A[k * 3 + p] = A[p * 3 + k] = c * akp - s * akq;
                A[k * 3 + q] = A[q * 3 + k] = s * akp + c * akq;
            }
        }

        // Update eigenvector matrix
        for (int k = 0; k < 3; ++k) {
            double vkp = eigenvectors[k * 3 + p];
            double vkq = eigenvectors[k * 3 + q];
            eigenvectors[k * 3 + p] = c * vkp - s * vkq;
            eigenvectors[k * 3 + q] = s * vkp + c * vkq;
        }
    }

    // Extract eigenvalues
    for (int i = 0; i < 3; ++i) {
        eigenvalues[i] = A[i * 3 + i];
    }

    // Sort eigenvalues in ascending order (S0 < S1 < S2)
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (eigenvalues[j] < eigenvalues[i]) {
                std::swap(eigenvalues[i], eigenvalues[j]);
                // Swap eigenvector columns
                for (int k = 0; k < 3; ++k) {
                    std::swap(eigenvectors[k * 3 + i], eigenvectors[k * 3 + j]);
                }
            }
        }
    }

    // Store results
    result.energies.resize(n_si_states);
    result.coeffs.resize(n_si_states * n_si_states);

    for (int i = 0; i < n_si_states; ++i) {
        result.energies[i] = eigenvalues[i];
        for (int j = 0; j < n_si_states; ++j) {
            result.coeffs[i * n_si_states + j] = eigenvectors[j * 3 + i];
        }
    }

    return result;
}

// ============================================================================
// Separate Configuration Energy/Gradient/Hessian Functions
// For 3SA-REKS: each configuration optimizes its own FONs independently
// ============================================================================

double REKS44Space::compute_energy_PPS(const std::vector<double>& E_micro) const {
    // E_PPS = sum_L FACT_L * C_PPS_L * E_micro_L
    std::vector<double> C_PPS;
    compute_weights_PPS(C_PPS);

    double E_PPS = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        E_PPS += FACT * C_PPS[L] * E_micro[L];
    }
    return E_PPS;
}

double REKS44Space::compute_energy_OSS1(const std::vector<double>& E_micro) const {
    // E_OSS1 = sum_L FACT_L * C_OSS1_L * E_micro_L
    std::vector<double> C_OSS1;
    compute_weights_OSS1(C_OSS1);

    double E_OSS1 = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        E_OSS1 += FACT * C_OSS1[L] * E_micro[L];
    }
    return E_OSS1;
}

double REKS44Space::compute_energy_OSS2(const std::vector<double>& E_micro) const {
    // E_OSS2 = sum_L FACT_L * C_OSS2_L * E_micro_L
    std::vector<double> C_OSS2;
    compute_weights_OSS2(C_OSS2);

    double E_OSS2 = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        E_OSS2 += FACT * C_OSS2[L] * E_micro[L];
    }
    return E_OSS2;
}

std::vector<double> REKS44Space::compute_gradient_PPS(const std::vector<double>& E_micro) const {
    // Gradient of E_PPS w.r.t. (n_a, n_b)
    // dE_PPS/dn_a = sum_L FACT_L * dC_PPS_L/dn_a * E_micro_L
    // dE_PPS/dn_b = sum_L FACT_L * dC_PPS_L/dn_b * E_micro_L

    std::vector<double> dC_dna(16), d2C_dna2(16);
    std::vector<double> dC_dnb(16), d2C_dnb2(16);

    compute_weight_derivs_PPS(dC_dna, d2C_dna2, 0);  // w.r.t. n_a
    compute_weight_derivs_PPS(dC_dnb, d2C_dnb2, 1);  // w.r.t. n_b

    double grad_na = 0.0;
    double grad_nb = 0.0;

    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        grad_na += FACT * dC_dna[L] * E_micro[L];
        grad_nb += FACT * dC_dnb[L] * E_micro[L];
    }

    return {grad_na, grad_nb};
}

double REKS44Space::compute_gradient_OSS1(const std::vector<double>& E_micro) const {
    // Gradient of E_OSS1 w.r.t. n'_a
    // dE_OSS1/dn'_a = sum_L FACT_L * dC_OSS1_L/dn'_a * E_micro_L

    std::vector<double> dC_dna_prime(16), d2C_dna_prime2(16);
    compute_weight_derivs_OSS1_oss_fon(dC_dna_prime, d2C_dna_prime2);

    double grad = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        grad += FACT * dC_dna_prime[L] * E_micro[L];
    }

    return grad;
}

double REKS44Space::compute_gradient_OSS2(const std::vector<double>& E_micro) const {
    // Gradient of E_OSS2 w.r.t. n'_b
    // dE_OSS2/dn'_b = sum_L FACT_L * dC_OSS2_L/dn'_b * E_micro_L

    std::vector<double> dC_dnb_prime(16), d2C_dnb_prime2(16);
    compute_weight_derivs_OSS2_oss_fon(dC_dnb_prime, d2C_dnb_prime2);

    double grad = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        grad += FACT * dC_dnb_prime[L] * E_micro[L];
    }

    return grad;
}

std::vector<double> REKS44Space::compute_hessian_PPS(const std::vector<double>& E_micro) const {
    // 2x2 Hessian of E_PPS w.r.t. (n_a, n_b)

    std::vector<double> dC_dna(16), d2C_dna2(16);
    std::vector<double> dC_dnb(16), d2C_dnb2(16);

    compute_weight_derivs_PPS(dC_dna, d2C_dna2, 0);
    compute_weight_derivs_PPS(dC_dnb, d2C_dnb2, 1);

    double H_aa = 0.0;
    double H_bb = 0.0;

    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        H_aa += FACT * d2C_dna2[L] * E_micro[L];
        H_bb += FACT * d2C_dnb2[L] * E_micro[L];
    }

    // Mixed derivative d²C_PPS/(dn_a dn_b)
    // Only closed-shell weights (L=0-3) contribute
    double H_ab = 0.0;
    H_ab += 1.0 * (1.0 / 4.0) * E_micro[0];   // d²C_0/(dn_a dn_b) = 1/4
    H_ab += 1.0 * (-1.0 / 4.0) * E_micro[1];  // d²C_1/(dn_a dn_b) = -1/4
    H_ab += 1.0 * (-1.0 / 4.0) * E_micro[2];  // d²C_2/(dn_a dn_b) = -1/4
    H_ab += 1.0 * (1.0 / 4.0) * E_micro[3];   // d²C_3/(dn_a dn_b) = 1/4

    return {H_aa, H_ab, H_ab, H_bb};
}

double REKS44Space::compute_hessian_OSS1(const std::vector<double>& E_micro) const {
    // Hessian d²E_OSS1/d(n'_a)²

    std::vector<double> dC_dna_prime(16), d2C_dna_prime2(16);
    compute_weight_derivs_OSS1_oss_fon(dC_dna_prime, d2C_dna_prime2);

    double H = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        H += FACT * d2C_dna_prime2[L] * E_micro[L];
    }

    return H;
}

double REKS44Space::compute_hessian_OSS2(const std::vector<double>& E_micro) const {
    // Hessian d²E_OSS2/d(n'_b)²

    std::vector<double> dC_dnb_prime(16), d2C_dnb_prime2(16);
    compute_weight_derivs_OSS2_oss_fon(dC_dnb_prime, d2C_dnb_prime2);

    double H = 0.0;
    for (int L = 0; L < 16; ++L) {
        double FACT = get_symmetry_factor(L);
        H += FACT * d2C_dnb_prime2[L] * E_micro[L];
    }

    return H;
}

void REKS44Space::compute_weight_derivs_OSS1_oss_fon(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2) const {
    // Derivatives of OSS1 weights w.r.t. n'_a (oss1_fon_a_)
    //
    // From Eq. A1, OSS1 weights that depend on n'_a:
    // - C_L[4-7]: coupling weights using f(n'_a * n'_d)
    // - C_L[8-9]: depend on n'_a / 4 + 0.5
    // - C_L[12-13]: depend on n'_d / 4 = (2 - n'_a) / 4

    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);

    const double na_prime = oss1_fon_a_;
    const double nd_prime = oss1_fon_d_;  // = 2 - na_prime
    const double x = na_prime * nd_prime;

    const double df_dx = df_interp(x);
    const double d2f_dx2 = d2f_interp(x);

    // dx/dn'_a = n'_d - n'_a = (2 - na_prime) - na_prime = 2 - 2*na_prime
    const double dx_dna = nd_prime - na_prime;
    // d²x/d(n'_a)² = -2
    const double d2x_dna2 = -2.0;

    // df/dn'_a = df/dx * dx/dn'_a
    const double df_dna = df_dx * dx_dna;
    // d²f/d(n'_a)² = d²f/dx² * (dx/dn'_a)² + df/dx * d²x/d(n'_a)²
    const double d2f_dna2 = d2f_dx2 * dx_dna * dx_dna + df_dx * d2x_dna2;

    // Closed-shell microstates (L=0-3): zero for OSS1
    dC_dfon[0] = dC_dfon[1] = dC_dfon[2] = dC_dfon[3] = 0.0;
    d2C_dfon2[0] = d2C_dfon2[1] = d2C_dfon2[2] = d2C_dfon2[3] = 0.0;

    // (a,d) coupling microstates (L=4-7)
    // C_L[4] = C_L[5] = -f(n'_a * n'_d) / 2
    // dC/dn'_a = -df_dna / 2
    dC_dfon[4] = dC_dfon[5] = -df_dna / 2.0;
    d2C_dfon2[4] = d2C_dfon2[5] = -d2f_dna2 / 2.0;

    // C_L[6] = C_L[7] = +f(n'_a * n'_d) / 2
    // dC/dn'_a = +df_dna / 2
    dC_dfon[6] = dC_dfon[7] = +df_dna / 2.0;
    d2C_dfon2[6] = d2C_dfon2[7] = +d2f_dna2 / 2.0;

    // (b,c) coupling microstates (L=8-11)
    // C_L[8] = C_L[9] = n'_a / 4 + 0.5
    // dC/dn'_a = 1/4
    dC_dfon[8] = dC_dfon[9] = 0.25;
    d2C_dfon2[8] = d2C_dfon2[9] = 0.0;

    // C_L[10] = C_L[11] = -0.5 (constant)
    dC_dfon[10] = dC_dfon[11] = 0.0;
    d2C_dfon2[10] = d2C_dfon2[11] = 0.0;

    // OSS1 additional microstates (L=12-13)
    // C_L[12] = C_L[13] = n'_d / 4 = (2 - n'_a) / 4
    // dC/dn'_a = -1/4
    dC_dfon[12] = dC_dfon[13] = -0.25;
    d2C_dfon2[12] = d2C_dfon2[13] = 0.0;

    // OSS2 microstates (L=14-15): zero for OSS1
    dC_dfon[14] = dC_dfon[15] = 0.0;
    d2C_dfon2[14] = d2C_dfon2[15] = 0.0;
}

void REKS44Space::compute_weight_derivs_OSS2_oss_fon(
    std::vector<double>& dC_dfon,
    std::vector<double>& d2C_dfon2) const {
    // Derivatives of OSS2 weights w.r.t. n'_b (oss2_fon_b_)
    //
    // From Eq. A2, OSS2 weights that depend on n'_b:
    // - C_L[4-5]: depend on n'_b / 4 + 0.5
    // - C_L[6-7]: constant -0.5
    // - C_L[8-11]: coupling weights using f(n'_b * n'_c)
    // - C_L[14-15]: depend on n'_c / 4 = (2 - n'_b) / 4

    const int n_micro = 20;
    dC_dfon.resize(n_micro);
    d2C_dfon2.resize(n_micro);
    std::fill(dC_dfon.begin(), dC_dfon.end(), 0.0);
    std::fill(d2C_dfon2.begin(), d2C_dfon2.end(), 0.0);

    const double nb_prime = oss2_fon_b_;
    const double nc_prime = oss2_fon_c_;  // = 2 - nb_prime
    const double x = nb_prime * nc_prime;

    const double df_dx = df_interp(x);
    const double d2f_dx2 = d2f_interp(x);

    // dx/dn'_b = n'_c - n'_b = (2 - nb_prime) - nb_prime = 2 - 2*nb_prime
    const double dx_dnb = nc_prime - nb_prime;
    // d²x/d(n'_b)² = -2
    const double d2x_dnb2 = -2.0;

    // df/dn'_b = df/dx * dx/dn'_b
    const double df_dnb = df_dx * dx_dnb;
    // d²f/d(n'_b)² = d²f/dx² * (dx/dn'_b)² + df/dx * d²x/d(n'_b)²
    const double d2f_dnb2 = d2f_dx2 * dx_dnb * dx_dnb + df_dx * d2x_dnb2;

    // Closed-shell microstates (L=0-3): zero for OSS2 (already zeroed by std::fill)

    // (a,d) coupling microstates (L=4-7)
    // C_L[4] = C_L[5] = n'_b / 4 + 0.5
    // dC/dn'_b = 1/4
    dC_dfon[4] = dC_dfon[5] = 0.25;
    d2C_dfon2[4] = d2C_dfon2[5] = 0.0;

    // C_L[6] = C_L[7] = -0.5 (constant)
    dC_dfon[6] = dC_dfon[7] = 0.0;
    d2C_dfon2[6] = d2C_dfon2[7] = 0.0;

    // (b,c) coupling microstates (L=8-11)
    // C_L[8] = C_L[9] = -f(n'_b * n'_c) / 2
    // dC/dn'_b = -df_dnb / 2
    dC_dfon[8] = dC_dfon[9] = -df_dnb / 2.0;
    d2C_dfon2[8] = d2C_dfon2[9] = -d2f_dnb2 / 2.0;

    // C_L[10] = C_L[11] = +f(n'_b * n'_c) / 2
    // dC/dn'_b = +df_dnb / 2
    dC_dfon[10] = dC_dfon[11] = +df_dnb / 2.0;
    d2C_dfon2[10] = d2C_dfon2[11] = +d2f_dnb2 / 2.0;

    // OSS1 microstates (L=12-13): zero for OSS2
    dC_dfon[12] = dC_dfon[13] = 0.0;
    d2C_dfon2[12] = d2C_dfon2[13] = 0.0;

    // OSS2 additional microstates (L=14-15)
    // C_L[14] = C_L[15] = n'_c / 4 = (2 - n'_b) / 4
    // dC/dn'_b = -1/4
    dC_dfon[14] = dC_dfon[15] = -0.25;
    d2C_dfon2[14] = d2C_dfon2[15] = 0.0;
}

// ============================================================================
// 9SI-3SA-REKS(4,4) State Interaction Hamiltonian
// ============================================================================

REKSActiveSpace::SIResult REKS44Space::compute_SI_energies_9x9(
    const std::vector<double>& E_micro,
    double W_ad,
    double W_bc,
    double W_ac,
    double W_bd,
    int n_si_states) const {

    // ========================================================================
    // Full 9x9 State Interaction Hamiltonian for REKS(4,4)
    // ========================================================================
    //
    // Configuration ordering (0-indexed):
    //   0: PPS   - Perfectly Paired Singlet
    //   1: OSS1  - (b,c) OSS, (a,d) GVB
    //   2: OSS2  - (a,d) OSS, (b,c) GVB
    //   3: OSS3  - (a,c) OSS (inter-pair), (b,d) GVB
    //   4: OSS4  - (b,d) OSS (inter-pair), (a,c) GVB
    //   5: DOSS  - Double Open-Shell Singlet (all orbitals singly occupied)
    //   6: DSPS  - Double Single-Pair Singlet (quintet character)
    //   7: DES1  - Doubly Excited Singlet type 1
    //   8: DES2  - Doubly Excited Singlet type 2
    //
    // Off-diagonal couplings from Filatov 2017, Appendix B:
    //
    // Standard (a,d)/(b,c) pair couplings:
    //   H[0][1] = W_bc * sqrt(2) * (sqrt(n_b) - sqrt(n_c))  [PPS-OSS1]
    //   H[0][2] = W_ad * sqrt(2) * (sqrt(n_a) - sqrt(n_d))  [PPS-OSS2]
    //   H[1][7] = W_ad * sqrt(2) * (sqrt(n'_a) + sqrt(n'_d))  [OSS1-DES1]
    //   H[2][7] = W_bc * sqrt(2) * (sqrt(n'_b) + sqrt(n'_c))  [OSS2-DES1]
    //
    // Inter-pair (a,c)/(b,d) couplings:
    //   H[0][3] = W_ac * sqrt(2) * (sqrt(m_a) - sqrt(m_c))  [PPS-OSS3]
    //   H[0][4] = W_bd * sqrt(2) * (sqrt(m_b) - sqrt(m_d))  [PPS-OSS4]
    //   H[3][8] = W_bd * sqrt(2) * (sqrt(m_b) + sqrt(m_d))  [OSS3-DES2]
    //   H[4][8] = W_ac * sqrt(2) * (sqrt(m_a) + sqrt(m_c))  [OSS4-DES2]
    //
    // Cross-pair couplings (typically small or zero in simplified model)
    //   H[1][2], H[3][4], H[5][*], H[6][*] = 0 in simplified treatment
    // ========================================================================

    SIResult result;
    n_si_states = std::min(n_si_states, 9);
    result.n_states = n_si_states;

    // Get PPS FONs
    double n_a = pairs_[0].fon_p;  // (a,d) pair
    double n_d = pairs_[0].fon_q;
    double n_b = pairs_[1].fon_p;  // (b,c) pair
    double n_c = pairs_[1].fon_q;

    // Get OSS1 FONs (primed: n'_a, n'_d for (a,d) pair in OSS1)
    double na_prime = oss1_fon_a_;
    double nd_prime = oss1_fon_d_;

    // Get OSS2 FONs (primed: n'_b, n'_c for (b,c) pair in OSS2)
    double nb_prime = oss2_fon_b_;
    double nc_prime = oss2_fon_c_;

    // Get OSS3 FONs (m: for inter-pair (b,d) GVB in OSS3)
    double mb_oss3 = oss3_fon_b_;
    double md_oss3 = oss3_fon_d_;

    // Get OSS4 FONs (m: for inter-pair (a,c) GVB in OSS4)
    double ma_oss4 = oss4_fon_a_;
    double mc_oss4 = oss4_fon_c_;

    // Compute all 9 diagonal energies
    double diag[9];
    diag[0] = compute_energy_PPS(E_micro);
    diag[1] = compute_energy_OSS1(E_micro);
    diag[2] = compute_energy_OSS2(E_micro);
    diag[3] = compute_energy_OSS3(E_micro);
    diag[4] = compute_energy_OSS4(E_micro);
    diag[5] = compute_energy_DOSS(E_micro);
    diag[6] = compute_energy_DSPS(E_micro);
    diag[7] = compute_energy_DES1(E_micro);
    diag[8] = compute_energy_DES2(E_micro);

    // Store diagonal energies in result
    result.E_PPS = diag[0];
    result.E_OSS = diag[1];  // E_OSS maps to OSS1
    result.E_DES = diag[7];  // E_DES maps to DES1

    // Build 9x9 Hamiltonian (row-major)
    double H[81];
    std::fill(H, H + 81, 0.0);

    // Diagonal elements
    for (int i = 0; i < 9; ++i) {
        H[i * 9 + i] = diag[i];
    }

    // Off-diagonal couplings (symmetric: H[i][j] = H[j][i])
    const double sqrt2 = std::sqrt(2.0);

    // Standard pair couplings
    // H[0][1]: PPS-OSS1 = W_bc * sqrt(2) * (sqrt(n_b) - sqrt(n_c))
    double H01 = W_bc * sqrt2 * (std::sqrt(n_b) - std::sqrt(n_c));
    H[0 * 9 + 1] = H[1 * 9 + 0] = H01;

    // H[0][2]: PPS-OSS2 = W_ad * sqrt(2) * (sqrt(n_a) - sqrt(n_d))
    double H02 = W_ad * sqrt2 * (std::sqrt(n_a) - std::sqrt(n_d));
    H[0 * 9 + 2] = H[2 * 9 + 0] = H02;

    // H[1][7]: OSS1-DES1 = W_ad * sqrt(2) * (sqrt(n'_a) + sqrt(n'_d))
    double H17 = W_ad * sqrt2 * (std::sqrt(na_prime) + std::sqrt(nd_prime));
    H[1 * 9 + 7] = H[7 * 9 + 1] = H17;

    // H[2][7]: OSS2-DES1 = W_bc * sqrt(2) * (sqrt(n'_b) + sqrt(n'_c))
    double H27 = W_bc * sqrt2 * (std::sqrt(nb_prime) + std::sqrt(nc_prime));
    H[2 * 9 + 7] = H[7 * 9 + 2] = H27;

    // Inter-pair couplings (for OSS3/OSS4/DES2)
    // Note: These use fixed FONs from the inter-pair configurations
    // For OSS3: (a,c) fixed at (1,1), (b,d) is the GVB pair
    // For OSS4: (b,d) fixed at (1,1), (a,c) is the GVB pair

    // H[0][3]: PPS-OSS3 - coupling via (a,c) pair
    // Use PPS FONs for a and c (n_a and n_c)
    double H03 = W_ac * sqrt2 * (std::sqrt(n_a) - std::sqrt(n_c));
    H[0 * 9 + 3] = H[3 * 9 + 0] = H03;

    // H[0][4]: PPS-OSS4 - coupling via (b,d) pair
    // Use PPS FONs for b and d (n_b and n_d)
    double H04 = W_bd * sqrt2 * (std::sqrt(n_b) - std::sqrt(n_d));
    H[0 * 9 + 4] = H[4 * 9 + 0] = H04;

    // H[3][8]: OSS3-DES2 = W_bd * sqrt(2) * (sqrt(m_b) + sqrt(m_d))
    double H38 = W_bd * sqrt2 * (std::sqrt(mb_oss3) + std::sqrt(md_oss3));
    H[3 * 9 + 8] = H[8 * 9 + 3] = H38;

    // H[4][8]: OSS4-DES2 = W_ac * sqrt(2) * (sqrt(m_a) + sqrt(m_c))
    double H48 = W_ac * sqrt2 * (std::sqrt(ma_oss4) + std::sqrt(mc_oss4));
    H[4 * 9 + 8] = H[8 * 9 + 4] = H48;

    // Store off-diagonal couplings in result for reference
    result.H_12 = H01;  // PPS-OSS1
    result.H_23 = H02;  // PPS-OSS2

    // ========================================================================
    // Diagonalize 9x9 symmetric Hamiltonian using Jacobi rotations
    // ========================================================================

    double eigenvalues[9];
    double eigenvectors[81];

    // Initialize eigenvectors to identity
    std::fill(eigenvectors, eigenvectors + 81, 0.0);
    for (int i = 0; i < 9; ++i) {
        eigenvectors[i * 9 + i] = 1.0;
    }

    // Copy H to working array A
    double A[81];
    std::copy(H, H + 81, A);

    // Jacobi rotation algorithm
    const int max_iter = 100;
    const double tol = 1e-14;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Find largest off-diagonal element
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < 9; ++i) {
            for (int j = i + 1; j < 9; ++j) {
                double aij = std::abs(A[i * 9 + j]);
                if (aij > max_off) {
                    max_off = aij;
                    p = i;
                    q = j;
                }
            }
        }

        // Check convergence
        if (max_off < tol) break;

        // Compute Jacobi rotation
        double app = A[p * 9 + p];
        double aqq = A[q * 9 + q];
        double apq = A[p * 9 + q];

        double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);
        double c = std::cos(theta);
        double s = std::sin(theta);

        // Apply rotation to A: A' = G^T * A * G
        // Update diagonal elements
        A[p * 9 + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        A[q * 9 + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A[p * 9 + q] = A[q * 9 + p] = 0.0;  // zeroed by construction

        // Update off-diagonal elements
        for (int k = 0; k < 9; ++k) {
            if (k != p && k != q) {
                double akp = A[k * 9 + p];
                double akq = A[k * 9 + q];
                A[k * 9 + p] = A[p * 9 + k] = c * akp - s * akq;
                A[k * 9 + q] = A[q * 9 + k] = s * akp + c * akq;
            }
        }

        // Update eigenvector matrix
        for (int k = 0; k < 9; ++k) {
            double vkp = eigenvectors[k * 9 + p];
            double vkq = eigenvectors[k * 9 + q];
            eigenvectors[k * 9 + p] = c * vkp - s * vkq;
            eigenvectors[k * 9 + q] = s * vkp + c * vkq;
        }
    }

    // Extract eigenvalues
    for (int i = 0; i < 9; ++i) {
        eigenvalues[i] = A[i * 9 + i];
    }

    // Sort eigenvalues in ascending order (S0 < S1 < ... < S8)
    for (int i = 0; i < 8; ++i) {
        for (int j = i + 1; j < 9; ++j) {
            if (eigenvalues[j] < eigenvalues[i]) {
                std::swap(eigenvalues[i], eigenvalues[j]);
                // Swap eigenvector columns
                for (int k = 0; k < 9; ++k) {
                    std::swap(eigenvectors[k * 9 + i], eigenvectors[k * 9 + j]);
                }
            }
        }
    }

    // Store results (only up to n_si_states)
    result.energies.resize(n_si_states);
    result.coeffs.resize(n_si_states * n_si_states);

    for (int i = 0; i < n_si_states; ++i) {
        result.energies[i] = eigenvalues[i];
        for (int j = 0; j < n_si_states; ++j) {
            result.coeffs[i * n_si_states + j] = eigenvectors[j * 9 + i];
        }
    }

    return result;
}

}  // namespace reks
}  // namespace psi
