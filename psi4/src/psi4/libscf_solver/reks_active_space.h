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

#ifndef REKS_ACTIVE_SPACE_H
#define REKS_ACTIVE_SPACE_H

/**
 * @file reks_active_space.h
 * @brief Level 2: Abstract interface for REKS(N,M) active space
 *
 * Defines the REKSActiveSpace abstract class that encapsulates:
 * - Active space parameters (N electrons, M orbitals)
 * - GVB pair structure and FON management
 * - Microstate occupation patterns
 * - Weight computation (the key abstraction for different REKS variants)
 *
 * Concrete implementations:
 * - REKS22Space: REKS(2,2) - 2 electrons in 2 orbitals
 * - REKS44Space: REKS(4,4) - 4 electrons in 4 orbitals (future)
 */

#include "reks_math.h"
#include <memory>
#include <vector>
#include <array>

namespace psi {
namespace reks {

// ============================================================================
// Abstract Base Class: REKSActiveSpace
// ============================================================================

/**
 * @class REKSActiveSpace
 * @brief Abstract interface defining REKS(N,M) active space structure
 *
 * This class provides the interface that the REKS SCF engine uses to
 * interact with specific active space implementations. The key abstraction
 * is weight computation: different (N,M) setups have different formulas
 * for computing C_L weights, but the SCF algorithm remains the same.
 */
class REKSActiveSpace {
public:
    virtual ~REKSActiveSpace() = default;

    // ========================================================================
    // Active Space Parameters
    // ========================================================================

    /// Number of active electrons (N in REKS(N,M))
    virtual int n_electrons() const = 0;

    /// Number of active orbitals (M in REKS(N,M))
    virtual int n_orbitals() const = 0;

    /// Number of GVB pairs (typically N/2)
    virtual int n_pairs() const = 0;

    /// Number of unique microstates
    virtual int n_microstates() const = 0;

    // ========================================================================
    // GVB Pair Management
    // ========================================================================

    /// Get GVB pair by index (0-based)
    virtual const GVBPair& pair(int i) const = 0;

    /// Set FON for pair i (fon_p; fon_q is automatically 2-fon_p)
    virtual void set_pair_fon(int i, double fon_p) = 0;

    // ========================================================================
    // Microstate Access
    // ========================================================================

    /// Get microstate by index (0-based)
    virtual const Microstate& microstate(int L) const = 0;

    // ========================================================================
    // Weight Computation (THE KEY ABSTRACTION!)
    // ========================================================================

    /**
     * @brief Compute C_L weights for all microstates
     *
     * The weights depend on current FONs and state-averaging parameters.
     * For REKS(2,2) with SA weights (w_PPS, w_OSS):
     *   C_0 = w_PPS * n_r/2
     *   C_1 = w_PPS * n_s/2
     *   C_2 = w_OSS - 0.5*w_PPS*f(n_r*n_s)
     *   C_3 = 0.5*w_PPS*f(n_r*n_s) - 0.5*w_OSS
     *
     * @param[out] C_L Vector of weights (resized if needed)
     */
    virtual void compute_weights(std::vector<double>& C_L) const = 0;

    /**
     * @brief Compute derivatives of C_L weights w.r.t. FON
     *
     * Required for Newton-Raphson FON optimization.
     *
     * @param[out] dC_dfon First derivatives dC_L/d(fon_p) for each L
     * @param[out] d2C_dfon2 Second derivatives d²C_L/d(fon_p)² for each L
     * @param pair_idx Index of the GVB pair (for REKS(2,2), always 0)
     */
    virtual void compute_weight_derivs(
        std::vector<double>& dC_dfon,
        std::vector<double>& d2C_dfon2,
        int pair_idx = 0) const = 0;

    // ========================================================================
    // State Averaging Parameters
    // ========================================================================

    /// Weight for PPS (Perfectly Paired Singlet) state
    virtual double w_PPS() const = 0;

    /// Weight for OSS (Open-Shell Singlet) state
    virtual double w_OSS() const = 0;

    /// Set state-averaging weights
    virtual void set_sa_weights(double w_pps, double w_oss) = 0;

    // ========================================================================
    // Microstate Occupation Access (for density matrix construction)
    // ========================================================================

    /**
     * @brief Get alpha occupation for orbital in microstate
     * @param L Microstate index
     * @param orbital Active orbital index (0-based)
     * @return 0 or 1
     */
    virtual int alpha_occ(int L, int orbital) const = 0;

    /**
     * @brief Get beta occupation for orbital in microstate
     * @param L Microstate index
     * @param orbital Active orbital index (0-based)
     * @return 0 or 1
     */
    virtual int beta_occ(int L, int orbital) const = 0;

    // ========================================================================
    // Orbital Index Management (NEW for REKS(N,M) generalization)
    // ========================================================================

    /**
     * @brief Get active orbital indices in MO space
     * @param Ncore Number of core (doubly occupied) orbitals
     * @return Vector of active MO indices (0-based)
     *
     * For REKS(2,2): returns {Ncore, Ncore+1}
     * For REKS(4,4): returns {Ncore, Ncore+1, Ncore+2, Ncore+3}
     */
    [[nodiscard]] virtual std::vector<int> get_active_mo_indices(int Ncore) const = 0;

    // ========================================================================
    // Base Density Pattern Enumeration (for efficient JK batching)
    // ========================================================================

    /**
     * @brief Get unique base density occupation patterns
     * @return Vector of bitmask patterns
     *
     * Pattern: bitmask where bit i = occupation of active orbital i
     * REKS(2,2): {0b00, 0b01, 0b10, 0b11} = {0, 1, 2, 3}
     * REKS(4,4): {0b0000, ..., 0b1111} = {0, ..., 15}
     */
    [[nodiscard]] virtual std::vector<int> get_base_density_patterns() const = 0;

    /**
     * @brief Number of unique base densities
     * @return 2^M for REKS(N,M)
     */
    [[nodiscard]] virtual int n_base_densities() const = 0;

    // ========================================================================
    // Microstate-to-BaseDensity Mapping
    // ========================================================================

    /**
     * @brief Get base density index for alpha density of microstate L
     * @param L Microstate index
     * @return Index into base_density_patterns array
     *
     * Converts alpha occupation pattern to bitmask index
     */
    [[nodiscard]] virtual int get_alpha_base_idx(int L) const = 0;

    /**
     * @brief Get base density index for beta density of microstate L
     * @param L Microstate index
     * @return Index into base_density_patterns array
     */
    [[nodiscard]] virtual int get_beta_base_idx(int L) const = 0;

    // ========================================================================
    // Energy Weighting (Symmetry Factors)
    // ========================================================================

    /**
     * @brief Get symmetry factor (FACT) for microstate L
     * @param L Microstate index
     * @return Multiplier for energy contribution (1.0 or 2.0)
     *
     * Accounts for spin-paired microstates counted once.
     * For REKS(2,2): FACT=1.0 for L=0,1; FACT=2.0 for L=2,3
     */
    [[nodiscard]] virtual double get_symmetry_factor(int L) const = 0;

    // ========================================================================
    // Multi-variable FON Optimization
    // ========================================================================

    /**
     * @brief Number of independent FON variables
     * @return Number of GVB pairs (1 for REKS(2,2), 2 for REKS(4,4), etc.)
     */
    [[nodiscard]] virtual int n_fon_variables() const = 0;

    /**
     * @brief Get all FON values as vector
     * @return Vector of fon_p values for each pair
     */
    [[nodiscard]] virtual std::vector<double> get_fon_vector() const = 0;

    /**
     * @brief Set FON values from vector
     * @param fons Vector of fon_p values (one per pair)
     */
    virtual void set_fon_vector(const std::vector<double>& fons) = 0;

    /**
     * @brief Compute gradient of E_SA w.r.t. FON variables
     * @param E_micro Vector of microstate energies
     * @param C_L Vector of current weights
     * @return Gradient vector of length n_fon_variables()
     *
     * For REKS(2,2): returns {dE/dn_r}
     * For REKS(4,4): returns {dE/dn_a, dE/dn_b}
     */
    [[nodiscard]] virtual std::vector<double> compute_energy_gradient(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const = 0;

    /**
     * @brief Compute Hessian of E_SA w.r.t. FON variables
     * @param E_micro Vector of microstate energies
     * @param C_L Vector of current weights
     * @return Hessian matrix (row-major, n_fon_variables x n_fon_variables)
     *
     * For REKS(2,2): returns {d²E/dn_r²}
     * For REKS(4,4): returns {d²E/dn_a², d²E/dn_a dn_b, d²E/dn_b dn_a, d²E/dn_b²}
     */
    [[nodiscard]] virtual std::vector<double> compute_energy_hessian(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const = 0;

    // ========================================================================
    // Fock Coupling Weights
    // ========================================================================

    /**
     * @brief Get effective FON for active orbital
     * @param orbital_idx Active orbital index (0-based within active space)
     * @return Effective occupation for Fock coupling
     *
     * f_i = (n_i/2)*w_PPS + 0.5*w_OSS for SA-REKS
     *
     * For REKS(2,2): orbital_idx=0 -> f_r, orbital_idx=1 -> f_s
     */
    [[nodiscard]] virtual double get_effective_fon(int orbital_idx) const = 0;

    // ========================================================================
    // SI-SA (State Interaction) Methods
    // ========================================================================

    /**
     * @brief SI Hamiltonian result structure
     *
     * For 2SI-2SA-REKS: 2x2 matrix between PPS and OSS
     * For 3SI-2SA-REKS: 3x3 matrix including DES
     */
    struct SIResult {
        int n_states;                    ///< Number of states (2 or 3)
        std::vector<double> energies;    ///< Eigenvalues (state energies)
        std::vector<double> coeffs;      ///< Eigenvectors (row-major: coeffs[state*n_states + component])
        double E_PPS;                    ///< PPS diagonal element
        double E_OSS;                    ///< OSS diagonal element
        double E_DES;                    ///< DES diagonal element (3SI only)
        double H_12;                     ///< PPS-OSS coupling
        double H_23;                     ///< OSS-DES coupling (3SI only)
    };

    /**
     * @brief Compute SI Hamiltonian and solve eigenvalue problem
     *
     * @param E_micro Microstate energies from SCF
     * @param Wrs Lagrange multiplier (off-diagonal coupling)
     * @param n_si_states Number of SI states (2 for 2SI, 3 for 3SI)
     * @return SIResult with eigenvalues and eigenvectors
     *
     * 2SI-2SA Hamiltonian:
     *   H_11 = E_PPS = sum_L C_L^{PPS} * FACT * E_L  (C_L computed with w_PPS=1, w_OSS=0)
     *   H_22 = E_OSS = 2*E_L[2] - E_L[3]
     *   H_12 = Wrs * (sqrt(n_r) - sqrt(n_s)) * sqrt(2)
     *
     * 3SI-2SA adds:
     *   H_33 = E_DES
     *   H_23 = Wrs * (sqrt(n_r) + sqrt(n_s)) * sqrt(2)
     *   H_13 = 0
     */
    [[nodiscard]] virtual SIResult compute_SI_energies(
        const std::vector<double>& E_micro,
        double Wrs,
        int n_si_states = 2) const = 0;

    // ========================================================================
    // Factory Methods
    // ========================================================================

    /// Create REKS(2,2) active space
    static std::unique_ptr<REKSActiveSpace> create_2_2(double w_pps = 0.5, double w_oss = 0.5);

    /// Create REKS(4,4) active space (future implementation)
    static std::unique_ptr<REKSActiveSpace> create_4_4(double w_pps = 1.0/3.0);
};

// ============================================================================
// REKS(2,2) Implementation
// ============================================================================

/**
 * @class REKS22Space
 * @brief REKS(2,2) active space: 2 electrons in 2 orbitals
 *
 * Single GVB pair (r, s) with FON constraint: n_r + n_s = 2.
 * 4 unique microstates (using spin symmetry).
 *
 * Microstate structure:
 *   L=0: r doubly occupied      alpha={1,0}, beta={1,0}
 *   L=1: s doubly occupied      alpha={0,1}, beta={0,1}
 *   L=2: r-alpha, s-beta        alpha={1,0}, beta={0,1}
 *   L=3: triplet-like           alpha={1,1}, beta={0,0}
 */
class REKS22Space final : public REKSActiveSpace {
public:
    /**
     * @brief Construct REKS(2,2) active space
     * @param w_pps Weight for PPS state (default: 0.5)
     * @param w_oss Weight for OSS state (default: 0.5)
     */
    explicit REKS22Space(double w_pps = 0.5, double w_oss = 0.5);

    // ========================================================================
    // REKSActiveSpace Interface Implementation
    // ========================================================================

    int n_electrons() const override { return 2; }
    int n_orbitals() const override { return 2; }
    int n_pairs() const override { return 1; }
    int n_microstates() const override { return 4; }

    const GVBPair& pair(int i) const override { return pair_; }
    void set_pair_fon(int i, double fon_p) override { pair_.set_fon(fon_p); }

    const Microstate& microstate(int L) const override { return microstates_[L]; }

    void compute_weights(std::vector<double>& C_L) const override;
    void compute_weight_derivs(
        std::vector<double>& dC_dfon,
        std::vector<double>& d2C_dfon2,
        int pair_idx = 0) const override;

    double w_PPS() const override { return w_pps_; }
    double w_OSS() const override { return w_oss_; }
    void set_sa_weights(double w_pps, double w_oss) override {
        w_pps_ = w_pps;
        w_oss_ = w_oss;
    }

    int alpha_occ(int L, int orbital) const override {
        return microstates_[L].alpha[orbital];
    }
    int beta_occ(int L, int orbital) const override {
        return microstates_[L].beta[orbital];
    }

    // ========================================================================
    // NEW Interface Methods Implementation
    // ========================================================================

    [[nodiscard]] std::vector<int> get_active_mo_indices(int Ncore) const override {
        return {Ncore, Ncore + 1};
    }

    [[nodiscard]] std::vector<int> get_base_density_patterns() const override {
        // 4 patterns: 0b00=0 (none), 0b01=1 (s only), 0b10=2 (r only), 0b11=3 (r and s)
        return {0, 1, 2, 3};
    }

    [[nodiscard]] int n_base_densities() const override { return 4; }  // 2^2

    [[nodiscard]] int get_alpha_base_idx(int L) const override;
    [[nodiscard]] int get_beta_base_idx(int L) const override;

    [[nodiscard]] double get_symmetry_factor(int L) const override {
        // L=0,1: closed-shell, FACT=1.0
        // L=2,3: open-shell, FACT=2.0 (spin symmetry)
        return (L >= 2) ? 2.0 : 1.0;
    }

    [[nodiscard]] int n_fon_variables() const override { return 1; }

    [[nodiscard]] std::vector<double> get_fon_vector() const override {
        return {pair_.fon_p};
    }

    void set_fon_vector(const std::vector<double>& fons) override {
        if (!fons.empty()) {
            pair_.set_fon(fons[0]);
        }
    }

    [[nodiscard]] std::vector<double> compute_energy_gradient(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const override;

    [[nodiscard]] std::vector<double> compute_energy_hessian(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const override;

    [[nodiscard]] double get_effective_fon(int orbital_idx) const override;

    [[nodiscard]] SIResult compute_SI_energies(
        const std::vector<double>& E_micro,
        double Wrs,
        int n_si_states = 2) const override;

    // ========================================================================
    // REKS(2,2) Specific Accessors
    // ========================================================================

    /// Get FON for orbital r (same as pair_.fon_p)
    double n_r() const { return pair_.fon_p; }

    /// Get FON for orbital s (same as pair_.fon_q)
    double n_s() const { return pair_.fon_q; }

    /// Set FON for orbital r (n_s is automatically 2 - n_r)
    void set_n_r(double n) { pair_.set_fon(n); }

    /// Get f_interp(n_r * n_s)
    double f_value() const { return f_interp(pair_.fon_p * pair_.fon_q); }

private:
    GVBPair pair_{0, 1, 1.0};           ///< Single GVB pair (r=0, s=1)
    double w_pps_;                       ///< SA weight for PPS
    double w_oss_;                       ///< SA weight for OSS
    std::array<Microstate, 4> microstates_;  ///< 4 unique microstates

    /// Initialize microstate occupation patterns
    void init_microstates();
};

// ============================================================================
// REKS(4,4) Implementation (placeholder for future)
// ============================================================================

/**
 * @class REKS44Space
 * @brief REKS(4,4) active space: 4 electrons in 4 orbitals (future)
 *
 * Two GVB pairs: (a, d) and (b, c)
 * Constraints: n_a + n_d = 2, n_b + n_c = 2
 *
 * Based on: Filatov et al. J. Chem. Phys. 147, 064104 (2017)
 */
class REKS44Space final : public REKSActiveSpace {
public:
    explicit REKS44Space(double w_pps = 1.0/3.0);

    int n_electrons() const override { return 4; }
    int n_orbitals() const override { return 4; }
    int n_pairs() const override { return 2; }
    int n_microstates() const override;

    const GVBPair& pair(int i) const override { return pairs_[i]; }
    void set_pair_fon(int i, double fon_p) override { pairs_[i].set_fon(fon_p); }

    const Microstate& microstate(int L) const override { return microstates_[L]; }

    void compute_weights(std::vector<double>& C_L) const override;
    void compute_weight_derivs(
        std::vector<double>& dC_dfon,
        std::vector<double>& d2C_dfon2,
        int pair_idx = 0) const override;

    double w_PPS() const override { return w_pps_; }
    double w_OSS() const override { return (1.0 - w_pps_) / 2.0; }
    void set_sa_weights(double w_pps, double w_oss) override { w_pps_ = w_pps; }

    int alpha_occ(int L, int orbital) const override {
        return microstates_[L].alpha[orbital];
    }
    int beta_occ(int L, int orbital) const override {
        return microstates_[L].beta[orbital];
    }

    // ========================================================================
    // NEW Interface Methods Implementation
    // ========================================================================

    [[nodiscard]] std::vector<int> get_active_mo_indices(int Ncore) const override {
        return {Ncore, Ncore + 1, Ncore + 2, Ncore + 3};
    }

    [[nodiscard]] std::vector<int> get_base_density_patterns() const override {
        // 16 patterns: 0b0000=0 to 0b1111=15
        std::vector<int> patterns;
        patterns.reserve(16);
        for (int p = 0; p < 16; ++p) {
            patterns.push_back(p);
        }
        return patterns;
    }

    [[nodiscard]] int n_base_densities() const override { return 16; }  // 2^4

    [[nodiscard]] int get_alpha_base_idx(int L) const override;
    [[nodiscard]] int get_beta_base_idx(int L) const override;

    [[nodiscard]] double get_symmetry_factor(int L) const override;

    [[nodiscard]] int n_fon_variables() const override { return 2; }

    [[nodiscard]] std::vector<double> get_fon_vector() const override {
        return {pairs_[0].fon_p, pairs_[1].fon_p};
    }

    void set_fon_vector(const std::vector<double>& fons) override {
        if (fons.size() >= 2) {
            pairs_[0].set_fon(fons[0]);
            pairs_[1].set_fon(fons[1]);
        }
    }

    [[nodiscard]] std::vector<double> compute_energy_gradient(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const override;

    [[nodiscard]] std::vector<double> compute_energy_hessian(
        const std::vector<double>& E_micro,
        const std::vector<double>& C_L) const override;

    [[nodiscard]] double get_effective_fon(int orbital_idx) const override;

    [[nodiscard]] SIResult compute_SI_energies(
        const std::vector<double>& E_micro,
        double Wrs,
        int n_si_states = 2) const override;

private:
    std::array<GVBPair, 2> pairs_;       ///< Two GVB pairs: (a,d) and (b,c)
    double w_pps_;
    std::vector<Microstate> microstates_;

    void init_microstates();
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_ACTIVE_SPACE_H
