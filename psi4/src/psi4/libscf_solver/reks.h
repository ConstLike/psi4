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

#ifndef REKS_H
#define REKS_H

#include "rhf.h"
#include "reks_active_space.h"
#include <vector>
#include <memory>

namespace psi {
namespace scf {

/**
 * @class REKS
 * @brief Restricted Ensemble Kohn-Sham (REKS) SCF engine
 *
 * This class implements the REKS SCF procedure using a pluggable
 * active space definition (REKSActiveSpace). The SCF algorithm is
 * the same for all REKS variants (2,2), (4,4), etc.
 *
 * Architecture:
 * - Level 1: reks_math.h - pure math functions (f_interp, etc.)
 * - Level 2: reks_active_space.h - active space interface + implementations
 * - Level 3: reks.h/cc - SCF engine (this class)
 *
 * Current Implementation: REKS(2,2) via REKS22Space
 * Future: REKS(4,4) via REKS44Space
 *
 * References:
 * - Filatov, M. "Note on SA-REKS, SSR, CP-REKS implementation" (2024)
 * - Filatov, M.; Shaik, S. Chem. Phys. Lett. 1999, 304, 429-437
 * - Filatov, M. et al. J. Chem. Phys. 147, 064104 (2017)
 */
class REKS : public RHF {
   protected:
    // ========================================================================
    // Active Space (polymorphic)
    // ========================================================================

    /// Active space definition (REKS22Space, REKS44Space, etc.)
    std::unique_ptr<reks::REKSActiveSpace> active_space_;

    // ========================================================================
    // Core Orbital Information
    // ========================================================================

    int Ncore_ = -1;                       ///< Number of core (doubly occupied) orbitals
    std::vector<int> active_mo_indices_;   ///< Active orbital MO indices (generalized from active_r_, active_s_)

    // ========================================================================
    // Base Density Matrices (generalized for REKS(N,M))
    // ========================================================================

    /// Base densities indexed by occupation bitmask
    /// For REKS(2,2): 4 densities (patterns 0-3)
    /// For REKS(4,4): 16 densities (patterns 0-15)
    std::vector<SharedMatrix> base_densities_;

    // ========================================================================
    // Microstate Data (dynamically sized)
    // ========================================================================

    std::vector<SharedMatrix> D_alpha_micro_;  ///< Alpha density for each microstate
    std::vector<SharedMatrix> D_beta_micro_;   ///< Beta density for each microstate
    std::vector<SharedMatrix> F_alpha_micro_;  ///< Alpha Fock for each microstate
    std::vector<SharedMatrix> F_beta_micro_;   ///< Beta Fock for each microstate
    std::vector<double> E_micro_;              ///< Energy for each microstate
    std::vector<double> C_L_;                  ///< Weighting factors for each microstate

    // ========================================================================
    // Coupling Fock Matrix
    // ========================================================================

    SharedMatrix F_reks_;          ///< Assembled coupling Fock matrix (AO basis)
    SharedMatrix F_reks_MO_;       ///< F_reks in MO basis
    double Wrs_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for SI

    // ========================================================================
    // Pre-allocated Work Matrices (for HPC efficiency)
    // ========================================================================

    /// Work matrices for C columns for each base density pattern
    std::vector<SharedMatrix> C_base_;     ///< Generalized from C_D00_, C_D10_, etc.
    SharedMatrix J_work_;                   ///< Work matrix: temporary for J_total
    SharedMatrix Temp_work_;                ///< Work matrix: temporary for AO→MO transforms
    std::vector<SharedMatrix> F_alpha_MO_;  ///< Pre-allocated MO-basis alpha Fock
    std::vector<SharedMatrix> F_beta_MO_;   ///< Pre-allocated MO-basis beta Fock

    // ========================================================================
    // Debug Level
    // ========================================================================

    int reks_debug_ = 0;           ///< 0=none, 1=energies, 2=matrices, 3=all

    // ========================================================================
    // Protected Methods
    // ========================================================================

    /// REKS-specific initialization
    void reks_common_init();

    /// Allocate all REKS-specific matrices
    void allocate_reks_matrices();

    /// Build base density matrices D00, D10, D01, D11
    void build_base_densities();

    /// Map base densities to microstate alpha/beta densities
    void build_microstate_densities();

    /// Build UHF-like Fock matrices for each microstate
    void build_microstate_focks();

    /// Compute microstate energies E_L
    void compute_microstate_energies();

    /// Newton-Raphson FON optimization
    void rex_solver();

    /// Compute weighting factors C_L (delegates to active_space_)
    void compute_weighting_factors();

    /// Add microstate L contribution to coupling Fock matrix
    void fock_micro_to_macro(int L, double Cl, double** Fa, double** Fb);

    /// Debug output: header
    void print_reks_header() const;
    /// Debug output: iteration info
    void print_reks_iteration(int iter, double E_tot, double dE, double dD) const;
    /// Debug output: microstate energies
    void print_microstate_energies() const;
    /// Debug output: FON info
    void print_fon_info() const;
    /// Print SI-SA state energies (called after SCF convergence)
    void print_SI_energies() const;

    // ========================================================================
    // SCF Method Overrides
    // ========================================================================

    void form_D() override;
    void form_G() override;
    void form_F() override;
    void form_C(double shift = 0.0) override;
    double compute_E() override;

   public:
    // ========================================================================
    // Constructors
    // ========================================================================

    /**
     * @brief Construct REKS with default REKS(2,2) active space
     * @param ref_wfn Reference wavefunction
     * @param functional DFT functional (or HF if nullptr)
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);

    /**
     * @brief Construct REKS with explicit options and PSIO
     * @param ref_wfn Reference wavefunction
     * @param functional DFT functional
     * @param options Psi4 options object
     * @param psio PSIO object for file I/O
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
         Options& options, std::shared_ptr<PSIO> psio);

    /// Destructor
    ~REKS() override;

    // ========================================================================
    // Public Interface
    // ========================================================================

    /// Number of electronic states
    int n_states() const override { return 1; }

    /// Create C1 symmetry deep copy
    std::shared_ptr<REKS> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    // ========================================================================
    // Accessors for Testing & Analysis
    // ========================================================================

    /// Get active space definition
    [[nodiscard]] const reks::REKSActiveSpace& active_space() const {
        return *active_space_;
    }

    /// Get number of microstates
    [[nodiscard]] int n_microstates() const {
        return active_space_->n_microstates();
    }

    /// Get microstate energy for index L
    [[nodiscard]] double get_microstate_energy(int L) const {
        return E_micro_[L];
    }

    /// Get weighting factor for microstate L
    [[nodiscard]] double get_C_L(int L) const {
        return C_L_[L];
    }

    // --- REKS(2,2) Specific Accessors (backward compatibility) ---

    /// Get FON for orbital r (only valid for REKS(2,2))
    [[nodiscard]] double get_n_r() const {
        return active_space_->pair(0).fon_p;
    }

    /// Get FON for orbital s (only valid for REKS(2,2))
    [[nodiscard]] double get_n_s() const {
        return active_space_->pair(0).fon_q;
    }

    /// Get f_interp value at current FON (only valid for REKS(2,2))
    [[nodiscard]] double get_f_value() const {
        const auto& p = active_space_->pair(0);
        return reks::f_interp(p.fon_p * p.fon_q);
    }
};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
