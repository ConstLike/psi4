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
#include <algorithm>
#include <utility>
#include <limits>

namespace psi {
namespace reks {

// ============================================================================
// Constants
// ============================================================================

/// Filatov interpolation parameter delta for all REKS variants.
/// Reference: Filatov et al. J. Chem. Phys. 147, 064104 (2017), Eq. 5
constexpr double DELTA = 0.4;

/// Numerical thresholds for FON optimization and convergence checks.
namespace constants {
    constexpr double ENERGY_THRESHOLD = 1e-10;   ///< Energy comparison threshold
    constexpr double FON_THRESHOLD = 1e-10;      ///< FON value threshold
    constexpr double HESSIAN_THRESHOLD = 1e-10;  ///< Hessian singularity threshold
    constexpr double NORM_THRESHOLD = 1e-14;     ///< Eigenvector normalization threshold
    constexpr int FON_MAX_ITER = 50;             ///< Maximum FON optimization iterations
    constexpr double FON_TOL = 1e-10;            ///< FON convergence tolerance
    constexpr double FON_MAX_STEP = 0.2;         ///< Maximum FON update step size

    // Oscillation detection threshold for FON optimizer
    constexpr double OSCILLATION_THRESHOLD = 0.02;  ///< Oscillation detection threshold

    // Line search constants (from Filatov's original GAMESS implementation)
    constexpr double GOLDEN_RATIO = 1.6180339887498949;  ///< Golden ratio φ = (1+√5)/2
    constexpr int LINE_SEARCH_MAX_ITER = 20;             ///< Maximum line search iterations

    // ========================================================================
    // Trust-Region constants for FON optimization
    // Reference: Trust-Region Augmented Hessian (TRAH) method
    // ========================================================================
    constexpr double TR_DELTA_INIT = 0.2;      ///< Initial trust region radius
    constexpr double TR_DELTA_MIN = 1e-10;     ///< Minimum trust region radius
    constexpr double TR_DELTA_MAX = 1.0;       ///< Maximum trust region radius
    constexpr double TR_ETA_ACCEPT = 0.0;      ///< Threshold for step acceptance (ρ > η_accept)
    constexpr double TR_ETA_GOOD = 0.25;       ///< Threshold for good step (radius unchanged)
    constexpr double TR_ETA_VERY_GOOD = 0.75;  ///< Threshold for very good step (radius increases)
    constexpr int TR_MAX_SECULAR_ITER = 50;    ///< Max iterations for secular equation solver
    constexpr double TR_SECULAR_TOL = 1e-12;   ///< Tolerance for secular equation convergence
}

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

// ============================================================================
// FON History for Geometry Scan Continuation
// ============================================================================

/**
 * @brief Single entry in FON history for geometry scans
 *
 * Stores converged FON values at a specific geometry parameter (R).
 * Used for extrapolation to provide smooth initial guesses.
 */
struct FONHistoryEntry {
    double R;         ///< Geometry parameter (scan coordinate)
    double n_a;       ///< Converged FON for pair (a,d)
    double n_b;       ///< Converged FON for pair (b,c)
    double E;         ///< Total energy
    int scf_iters;    ///< SCF iterations (diagnostic)

    FONHistoryEntry() : R(0.0), n_a(1.0), n_b(1.0), E(0.0), scf_iters(0) {}

    FONHistoryEntry(double r, double na, double nb, double e, int iters = 0)
        : R(r), n_a(na), n_b(nb), E(e), scf_iters(iters) {}
};

/**
 * @brief FON history manager for smooth PES scans
 *
 * Maintains a history of converged FON values across geometries.
 * Provides extrapolation methods to predict FON at new geometries,
 * enabling smooth transitions and diabatic tracking.
 *
 * Usage:
 *   1. After SCF convergence: history.add_point(R, n_a, n_b, E)
 *   2. At new geometry: auto [pred_na, pred_nb] = history.predict(R_new)
 *   3. Use predicted values as initial guess for FON optimization
 *
 * The class is designed to be used as a static singleton for persistence
 * across PSI4 energy() calls during geometry scans.
 */
class FONHistory {
public:
    /**
     * @brief Add a converged point to history
     * @param R Geometry parameter (scan coordinate)
     * @param n_a Converged FON for pair (a,d)
     * @param n_b Converged FON for pair (b,c)
     * @param E Total energy
     * @param iters SCF iterations (optional diagnostic)
     */
    void add_point(double R, double n_a, double n_b, double E, int iters = 0) {
        // Check if point already exists (update it)
        for (auto& entry : points_) {
            if (std::abs(entry.R - R) < 1e-10) {
                entry.n_a = n_a;
                entry.n_b = n_b;
                entry.E = E;
                entry.scf_iters = iters;
                return;
            }
        }

        // Add new point
        points_.emplace_back(R, n_a, n_b, E, iters);

        // Sort by R
        std::sort(points_.begin(), points_.end(),
                  [](const FONHistoryEntry& a, const FONHistoryEntry& b) {
                      return a.R < b.R;
                  });

        // Trim if too large
        while (static_cast<int>(points_.size()) > max_size_) {
            points_.erase(points_.begin());  // Remove oldest (smallest R)
        }
    }

    /**
     * @brief Predict FON at new geometry using extrapolation
     * @param R_new New geometry parameter
     * @param order Extrapolation order: 0=nearest, 1=linear, 2=quadratic
     * @return Pair (n_a_predicted, n_b_predicted)
     */
    std::pair<double, double> predict(double R_new, int order = 2) const {
        if (points_.empty()) {
            return {1.0, 1.0};  // Default: equal mixture
        }

        if (points_.size() == 1) {
            return {points_[0].n_a, points_[0].n_b};
        }

        if (order == 0 || points_.size() < 2) {
            // Nearest neighbor
            return extrapolate_nearest(R_new);
        } else if (order == 1 || points_.size() < 3) {
            // Linear extrapolation
            return extrapolate_linear(R_new);
        } else {
            // Quadratic extrapolation
            return extrapolate_quadratic(R_new);
        }
    }

    /**
     * @brief Clear all history
     */
    void clear() {
        points_.clear();
    }

    /**
     * @brief Number of points in history
     */
    int size() const {
        return static_cast<int>(points_.size());
    }

    /**
     * @brief Check if history is empty
     */
    bool empty() const {
        return points_.empty();
    }

    /**
     * @brief Get last (most recent) entry
     */
    const FONHistoryEntry& last() const {
        return points_.back();
    }

    /**
     * @brief Get entry at index
     */
    const FONHistoryEntry& at(int i) const {
        return points_.at(i);
    }

    /**
     * @brief Set maximum history size
     */
    void set_max_size(int n) {
        max_size_ = n;
    }

private:
    std::vector<FONHistoryEntry> points_;  ///< History storage
    int max_size_ = 100;                    ///< Maximum entries to keep

    /**
     * @brief Find indices of n nearest points to R
     */
    std::vector<int> find_nearest(double R, int n) const {
        std::vector<std::pair<double, int>> distances;
        for (int i = 0; i < static_cast<int>(points_.size()); ++i) {
            distances.emplace_back(std::abs(points_[i].R - R), i);
        }
        std::sort(distances.begin(), distances.end());

        std::vector<int> indices;
        for (int i = 0; i < n && i < static_cast<int>(distances.size()); ++i) {
            indices.push_back(distances[i].second);
        }
        return indices;
    }

    /**
     * @brief Nearest neighbor extrapolation
     */
    std::pair<double, double> extrapolate_nearest(double R_new) const {
        auto idx = find_nearest(R_new, 1);
        if (idx.empty()) return {1.0, 1.0};
        return {points_[idx[0]].n_a, points_[idx[0]].n_b};
    }

    /**
     * @brief Linear extrapolation using 2 nearest points
     */
    std::pair<double, double> extrapolate_linear(double R_new) const {
        auto idx = find_nearest(R_new, 2);
        if (idx.size() < 2) return extrapolate_nearest(R_new);

        // Sort by R
        if (points_[idx[0]].R > points_[idx[1]].R) std::swap(idx[0], idx[1]);

        double R0 = points_[idx[0]].R, R1 = points_[idx[1]].R;
        double dR = R1 - R0;
        if (std::abs(dR) < 1e-14) return extrapolate_nearest(R_new);

        // Linear interpolation/extrapolation
        double t = (R_new - R0) / dR;
        double n_a = points_[idx[0]].n_a + t * (points_[idx[1]].n_a - points_[idx[0]].n_a);
        double n_b = points_[idx[0]].n_b + t * (points_[idx[1]].n_b - points_[idx[0]].n_b);

        // Clamp to valid range
        n_a = std::max(0.0, std::min(2.0, n_a));
        n_b = std::max(0.0, std::min(2.0, n_b));

        return {n_a, n_b};
    }

    /**
     * @brief Quadratic extrapolation using 3 nearest points (Lagrange)
     */
    std::pair<double, double> extrapolate_quadratic(double R_new) const {
        auto idx = find_nearest(R_new, 3);
        if (idx.size() < 3) return extrapolate_linear(R_new);

        // Sort by R
        std::sort(idx.begin(), idx.end(),
                  [this](int a, int b) { return points_[a].R < points_[b].R; });

        double R0 = points_[idx[0]].R, R1 = points_[idx[1]].R, R2 = points_[idx[2]].R;

        // Lagrange coefficients
        double denom0 = (R0 - R1) * (R0 - R2);
        double denom1 = (R1 - R0) * (R1 - R2);
        double denom2 = (R2 - R0) * (R2 - R1);

        // Avoid division by zero
        if (std::abs(denom0) < 1e-20 || std::abs(denom1) < 1e-20 || std::abs(denom2) < 1e-20) {
            return extrapolate_linear(R_new);
        }

        double L0 = ((R_new - R1) * (R_new - R2)) / denom0;
        double L1 = ((R_new - R0) * (R_new - R2)) / denom1;
        double L2 = ((R_new - R0) * (R_new - R1)) / denom2;

        double n_a = L0 * points_[idx[0]].n_a + L1 * points_[idx[1]].n_a + L2 * points_[idx[2]].n_a;
        double n_b = L0 * points_[idx[0]].n_b + L1 * points_[idx[1]].n_b + L2 * points_[idx[2]].n_b;

        // Clamp to valid range
        n_a = std::max(0.0, std::min(2.0, n_a));
        n_b = std::max(0.0, std::min(2.0, n_b));

        return {n_a, n_b};
    }
};

// ============================================================================
// Multi-Start FON Optimization Support
// ============================================================================

/**
 * @brief Characterization of FON solution type
 */
enum class FONCharacter {
    DIRADICAL,     ///< n_a, n_b ≈ 1.0 (strong correlation)
    CLOSED_SHELL,  ///< n_a, n_b ≈ 2.0 (weak correlation)
    MIXED          ///< Intermediate values
};

/**
 * @brief Represents a converged FON solution from multi-start optimization
 */
struct FONSolution {
    double n_a = 1.0;     ///< FON for pair (a,d)
    double n_b = 1.0;     ///< FON for pair (b,c)
    double E_PPS = 0.0;   ///< Energy at this solution

    FONCharacter character = FONCharacter::MIXED;
    int start_id = -1;    ///< Which starting point led here
    bool converged = false;

    /**
     * @brief Check if this is a duplicate of another solution
     * @param other Another solution
     * @param tol Tolerance for FON comparison
     */
    bool is_duplicate_of(const FONSolution& other, double tol = 0.01) const {
        return std::abs(n_a - other.n_a) < tol && std::abs(n_b - other.n_b) < tol;
    }

    /**
     * @brief Compute distance to previous FON values (for diabatic tracking)
     */
    double diabatic_distance(double prev_na, double prev_nb) const {
        return std::sqrt(std::pow(n_a - prev_na, 2) + std::pow(n_b - prev_nb, 2));
    }

    /**
     * @brief Classify the character based on FON values
     */
    void classify() {
        double avg_n = (n_a + n_b) / 2.0;
        if (avg_n < 1.3) {
            character = FONCharacter::DIRADICAL;
        } else if (avg_n > 1.7) {
            character = FONCharacter::CLOSED_SHELL;
        } else {
            character = FONCharacter::MIXED;
        }
    }
};

/**
 * @brief Branch selection criterion for multi-start optimization
 */
enum class BranchCriterion {
    ENERGY,   ///< Always select lowest energy (adiabatic)
    DIABATIC, ///< Select closest to previous geometry
    AUTO      ///< Diabatic if within tolerance, else energy
};

/**
 * @brief Starting point for multi-start optimization
 */
struct FONStartPoint {
    double n_a;
    double n_b;
    std::string name;

    FONStartPoint(double na, double nb, const std::string& n)
        : n_a(na), n_b(nb), name(n) {}
};

/**
 * @brief Generate starting points for multi-start FON optimization
 */
class MultiStartGenerator {
public:
    /**
     * @brief Generate all starting points
     * @param prev_na Previous geometry FON (or -1 if none)
     * @param prev_nb Previous geometry FON (or -1 if none)
     * @param ham_na HAM-DIAG suggested FON
     * @param ham_nb HAM-DIAG suggested FON
     * @return Vector of starting points
     */
    std::vector<FONStartPoint> generate(
        double prev_na, double prev_nb,
        double ham_na, double ham_nb) const
    {
        std::vector<FONStartPoint> starts;

        // 1. Canonical points - always included
        starts.emplace_back(1.0, 1.0, "DIRADICAL");
        starts.emplace_back(2.0, 2.0, "CLOSED_SHELL");

        // 2. HAM-DIAG suggestion (if different from canonical)
        if (std::abs(ham_na - 1.0) > 0.1 || std::abs(ham_nb - 1.0) > 0.1) {
            if (std::abs(ham_na - 2.0) > 0.1 || std::abs(ham_nb - 2.0) > 0.1) {
                starts.emplace_back(ham_na, ham_nb, "HAM_DIAG");
            }
        }

        // 3. Continuation from previous geometry (if available)
        if (prev_na >= 0.0 && prev_nb >= 0.0) {
            starts.emplace_back(prev_na, prev_nb, "CONTINUATION");
        }

        return starts;
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_MATH_H