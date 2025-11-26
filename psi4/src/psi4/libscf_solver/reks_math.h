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

#ifndef REKS_MATH_H
#define REKS_MATH_H

/**
 * @file reks_math.h
 * @brief Level 1: Pure mathematical functions for REKS(N,M)
 *
 * This header-only file contains:
 * - Filatov interpolating function f(x) and derivatives (universal for all REKS)
 * - GVBPair struct for orbital pair with FON constraint
 * - Microstate struct for Slater determinant occupation pattern
 *
 * All functions are inline for maximum performance (no virtual calls).
 *
 * References:
 * - Filatov, M. et al. J. Chem. Phys. 147, 064104 (2017) - Eq. 5
 * - Filatov, M.; Shaik, S. Chem. Phys. Lett. 1999, 304, 429-437
 */

#include <vector>
#include <cmath>
#include <string>
#include <sstream>

namespace psi {
namespace reks {

// ============================================================================
// Constants
// ============================================================================

/// Filatov interpolation parameter delta - universal for all REKS variants
/// From: Filatov et al. J. Chem. Phys. 147, 064104 (2017), Eq. 5
constexpr double DELTA = 0.4;

// ============================================================================
// Interpolating Function f(x) and Derivatives
// ============================================================================

/**
 * @brief Filatov interpolating function f(x)
 *
 * From Eq. 5 of Filatov et al. J. Chem. Phys. 147, 064104 (2017):
 *   f(x) = x^(1 - (x + delta) / (2 * (1 + delta)))
 *
 * where x = n_p * n_q is the product of FONs for a GVB pair.
 *
 * This function smoothly interpolates between:
 * - Strongly correlated limit (x -> 1): f(1) = 1
 * - Weakly correlated limit (x -> 0): f(0) = 0
 *
 * @param x Product of FONs: n_p * n_q (range: 0 to 1)
 * @return Interpolated value (range: 0 to 1)
 */
inline double f_interp(double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    // exponent = 1 - (x + delta) / (2 * (1 + delta))
    const double c = 0.5 / (1.0 + DELTA);  // = 1 / (2 * (1 + delta))
    const double exponent = 1.0 - c * (x + DELTA);

    return std::pow(x, exponent);
}

/**
 * @brief First derivative df/dx of interpolating function
 *
 * Using chain rule on f(x) = x^g(x) where g(x) = 1 - c*(x + delta):
 *   df/dx = f(x) * (g(x)/x - c*ln(x))
 *         = f(x) * (exponent/x - c*ln(x))
 *
 * @param x Product of FONs: n_p * n_q
 * @return First derivative df/dx
 */
inline double df_interp(double x) {
    // Handle boundary cases
    if (x <= 1e-14) return 0.0;  // Avoid log(0) and division by zero
    if (x >= 1.0) return 0.0;    // Derivative is zero at saturation

    const double c = 0.5 / (1.0 + DELTA);
    const double exponent = 1.0 - c * (x + DELTA);
    const double f = std::pow(x, exponent);
    const double log_x = std::log(x);

    // df/dx = f * (exponent/x - c*log(x))
    return f * (exponent / x - c * log_x);
}

/**
 * @brief Second derivative d²f/dx² of interpolating function
 *
 * Differentiating df/dx = f * (exponent/x - c*ln(x)):
 *   d²f/dx² = (df/dx) * (exponent/x - c*ln(x))
 *           + f * (-exponent/x² - c/x)
 *
 * Let term1 = exponent/x - c*ln(x), then:
 *   d²f/dx² = f * (term1² - exponent/x² - c/x)
 *
 * @param x Product of FONs: n_p * n_q
 * @return Second derivative d²f/dx²
 */
inline double d2f_interp(double x) {
    // Handle boundary cases
    if (x <= 1e-14) return 0.0;
    if (x >= 1.0) return 0.0;

    const double c = 0.5 / (1.0 + DELTA);
    const double exponent = 1.0 - c * (x + DELTA);
    const double f = std::pow(x, exponent);
    const double log_x = std::log(x);

    const double term1 = exponent / x - c * log_x;
    const double term2 = -exponent / (x * x) - c / x;

    return f * (term1 * term1 + term2);
}

// ============================================================================
// GVB Pair Structure
// ============================================================================

/**
 * @brief GVB orbital pair with FON constraint
 *
 * Represents a pair of orbitals (p, q) in the active space with
 * fractional occupation numbers satisfying: fon_p + fon_q = 2.0
 *
 * For REKS(2,2): Single pair (r, s) with n_r + n_s = 2
 * For REKS(4,4): Two pairs (a, d) and (b, c):
 *   - n_a + n_d = 2
 *   - n_b + n_c = 2
 */
struct GVBPair {
    int p;           ///< First orbital index in active space (0-based)
    int q;           ///< Second orbital index in active space (0-based)
    double fon_p;    ///< FON of orbital p
    double fon_q;    ///< FON of orbital q = 2 - fon_p

    /**
     * @brief Construct GVB pair with initial FON
     * @param p_ First orbital index
     * @param q_ Second orbital index
     * @param n_p Initial FON for orbital p (default: 1.0 for equal occupancy)
     */
    GVBPair(int p_ = 0, int q_ = 1, double n_p = 1.0)
        : p(p_), q(q_), fon_p(n_p), fon_q(2.0 - n_p) {}

    /**
     * @brief Set FON for orbital p (automatically updates q)
     * @param n_p New FON for orbital p
     */
    void set_fon(double n_p) {
        fon_p = n_p;
        fon_q = 2.0 - n_p;
    }

    /**
     * @brief Get product of FONs (argument for f_interp)
     * @return fon_p * fon_q
     */
    double product() const {
        return fon_p * fon_q;
    }

    /**
     * @brief Get f_interp value for this pair's FONs
     * @return f_interp(fon_p * fon_q)
     */
    double f_value() const {
        return f_interp(fon_p * fon_q);
    }
};

// ============================================================================
// Microstate Structure
// ============================================================================

/**
 * @brief Single Slater determinant occupation pattern
 *
 * Represents integer occupations (0 or 1) for alpha and beta electrons
 * in each active orbital. The full density is built from core orbitals
 * (always doubly occupied) plus this active space pattern.
 *
 * For REKS(2,2), 4 unique microstates (using symmetry):
 *   L=0: alpha={1,0}, beta={1,0}  (r doubly occupied)
 *   L=1: alpha={0,1}, beta={0,1}  (s doubly occupied)
 *   L=2: alpha={1,0}, beta={0,1}  (r-alpha, s-beta)
 *   L=3: alpha={1,1}, beta={0,0}  (triplet-like component)
 */
struct Microstate {
    std::vector<int> alpha;  ///< Alpha spin occupations (0 or 1) per active orbital
    std::vector<int> beta;   ///< Beta spin occupations (0 or 1) per active orbital

    /**
     * @brief Default constructor
     */
    Microstate() = default;

    /**
     * @brief Construct with given size
     * @param n_orbitals Number of active orbitals
     */
    explicit Microstate(int n_orbitals)
        : alpha(n_orbitals, 0), beta(n_orbitals, 0) {}

    /**
     * @brief Construct with explicit occupations
     * @param a Alpha occupations
     * @param b Beta occupations
     */
    Microstate(std::vector<int> a, std::vector<int> b)
        : alpha(std::move(a)), beta(std::move(b)) {}

    /**
     * @brief Number of active orbitals
     */
    int n_orbitals() const { return static_cast<int>(alpha.size()); }

    /**
     * @brief Total alpha electrons in active space
     */
    int n_alpha() const {
        int n = 0;
        for (int a : alpha) n += a;
        return n;
    }

    /**
     * @brief Total beta electrons in active space
     */
    int n_beta() const {
        int n = 0;
        for (int b : beta) n += b;
        return n;
    }

    /**
     * @brief Total electrons in active space
     */
    int n_electrons() const {
        return n_alpha() + n_beta();
    }

    /**
     * @brief Spin projection Sz = (n_alpha - n_beta) / 2
     */
    double Sz() const {
        return 0.5 * (n_alpha() - n_beta());
    }

    /**
     * @brief Check if closed-shell (alpha == beta)
     */
    bool is_closed_shell() const {
        return alpha == beta;
    }

    /**
     * @brief Total occupation of orbital i (0, 1, or 2)
     */
    int occupation(int i) const {
        return alpha[i] + beta[i];
    }

    /**
     * @brief String representation for debugging
     */
    std::string to_string() const {
        std::ostringstream oss;
        oss << "alpha={";
        for (size_t i = 0; i < alpha.size(); i++) {
            if (i > 0) oss << ",";
            oss << alpha[i];
        }
        oss << "}, beta={";
        for (size_t i = 0; i < beta.size(); i++) {
            if (i > 0) oss << ",";
            oss << beta[i];
        }
        oss << "}";
        return oss.str();
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_MATH_H