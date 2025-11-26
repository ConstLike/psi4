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

private:
    std::array<GVBPair, 2> pairs_;       ///< Two GVB pairs: (a,d) and (b,c)
    double w_pps_;
    std::vector<Microstate> microstates_;

    void init_microstates();
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_ACTIVE_SPACE_H
