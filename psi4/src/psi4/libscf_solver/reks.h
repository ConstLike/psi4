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
#include <array>

namespace psi {
namespace scf {

/**
 * @class REKS
 * @brief Restricted Ensemble Kohn-Sham (REKS) method implementation.
 *
 * REKS is a multi-configurational DFT method that treats static correlation
 * through ensemble averaging of fractionally occupied Kohn-Sham determinants.
 *
 * Current Implementation: REKS(2,2) - 2 electrons in 2 active orbitals
 * Based on: Filatov, M. "Note on the implementation of SA-REKS, SSR, and CP-REKS" (2024)
 *
 * Key features:
 * - 4 unique microstates (from 6 total, using symmetry: L=3=L=4, L=5=L=6)
 * - Fractional occupation numbers (FONs) with constraint: n_r + n_s = 2
 * - Newton-Raphson FON optimization (RexSolver)
 * - Coupling operator technique for orbital optimization
 * - State-averaging (SA-REKS) for ground and excited states
 *
 * References:
 * - Filatov, M.; Shaik, S. Chem. Phys. Lett. 1999, 304, 429-437
 * - Filatov, M. WIREs Comput. Mol. Sci. 2015, 5, 146-167
 * - Filatov, M. "Note on SA-REKS, SSR, CP-REKS implementation" 2024
 */
class REKS : public RHF {
   protected:
    // --- Active Space Parameters ---
    int n_active_electrons_ = 2;   ///< Number of active electrons (2 for REKS(2,2))
    int n_active_orbitals_ = 2;    ///< Number of active orbitals (2 for REKS(2,2))
    int Ncore_ = -1;               ///< Number of core orbitals (computed from nalpha - 1)
    int active_r_ = -1;            ///< Index of active orbital r (= Ncore)
    int active_s_ = -1;            ///< Index of active orbital s (= Ncore + 1)

    // --- Fractional Occupation Numbers ---
    double n_r_ = 1.0;             ///< FON for orbital r, constraint: n_r + n_s = 2
    double n_s_ = 1.0;             ///< FON for orbital s = 2 - n_r
    static constexpr double DELTA_ = 0.4;  ///< Filatov interpolation parameter

    // --- State Averaging Weights ---
    double w_PPS_ = 0.5;           ///< Weight for PPS (perfectly paired singlet) state
    double w_OSS_ = 0.5;           ///< Weight for OSS (open-shell singlet) state

    // --- Base Density Matrices ---
    SharedMatrix D00_;             ///< Core only density
    SharedMatrix D10_;             ///< Core + r singly occupied
    SharedMatrix D01_;             ///< Core + s singly occupied
    SharedMatrix D11_;             ///< Core + r + s singly occupied

    // --- Microstate Data (4 unique microstates) ---
    static constexpr int N_MICRO_ = 4;  ///< Only 4 needed (L=3=L=4, L=5=L=6 by symmetry)

    std::array<SharedMatrix, N_MICRO_> D_alpha_micro_;  ///< D^alpha for L=1,2,3,5
    std::array<SharedMatrix, N_MICRO_> D_beta_micro_;   ///< D^beta for L=1,2,3,5
    std::array<SharedMatrix, N_MICRO_> F_alpha_micro_;  ///< F^alpha for L=1,2,3,5
    std::array<SharedMatrix, N_MICRO_> F_beta_micro_;   ///< F^beta for L=1,2,3,5
    std::array<double, N_MICRO_> E_micro_;              ///< E_L for L=1,2,3,5

    // --- Coupling Fock Matrix ---
    SharedMatrix F_reks_;          ///< Assembled coupling Fock matrix (AO basis)
    SharedMatrix F_reks_MO_;       ///< F_reks in MO basis (for GAMESS-style orbital update)
    double Wrs_lagr_ = 0.0;        ///< Off-diagonal Lagrange multiplier for SI

    // --- Weighting Factors ---
    std::array<double, N_MICRO_> C_L_;  ///< Weighting factors for SA-REKS

    // --- Microstate Occupation Tables ---
    /// From Filatov 2024 Table 2 (for 4 unique microstates: L=1,2,3,5)
    static constexpr int nr_alpha_[4] = {1, 0, 1, 1};  ///< n^alpha_r for L=1,2,3,5
    static constexpr int nr_beta_[4]  = {1, 0, 0, 0};  ///< n^beta_r for L=1,2,3,5
    static constexpr int ns_alpha_[4] = {0, 1, 0, 1};  ///< n^alpha_s for L=1,2,3,5
    static constexpr int ns_beta_[4]  = {0, 1, 1, 0};  ///< n^beta_s for L=1,2,3,5

    // --- Debug Level ---
    int reks_debug_ = 0;           ///< 0=none, 1=energies, 2=matrices, 3=all

    // --- Protected Methods ---
    /// REKS-specific initialization
    void reks_common_init();

    /// Allocate all REKS-specific matrices
    void allocate_reks_matrices();

    /// Interpolating function f(x) - Eq. (4) Filatov 2024
    double f_interp(double x) const;
    /// First derivative df/dx
    double df_interp(double x) const;
    /// Second derivative d2f/dx2
    double d2f_interp(double x) const;

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

    /// Compute weighting factors C_L for SA/PPS/OSS
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

    // --- SCF Method Overrides ---
    /// Build REKS density matrices (D00, D10, D01, D11 -> microstate densities)
    void form_D() override;
    /// Build microstate Fock matrices and compute microstate energies
    void form_G() override;
    /// Build REKS coupling Fock matrix from microstate Fock matrices
    /// NOTE: Also updates orbitals (Ca_, epsilon_a_) GAMESS-style via MO-basis diagonalization
    void form_F() override;
    /// REKS uses GAMESS-style orbital update in form_F(), so form_C does nothing
    void form_C(double shift = 0.0) override;
    /// Compute state-averaged REKS energy
    double compute_E() override;

   public:
    /**
     * @brief Construct REKS wavefunction from reference wavefunction and functional.
     * @param ref_wfn Reference wavefunction (provides basis, molecule, etc.)
     * @param functional DFT functional (or HF if nullptr)
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);

    /**
     * @brief Construct REKS wavefunction with explicit options and PSIO.
     * @param ref_wfn Reference wavefunction
     * @param functional DFT functional
     * @param options Psi4 options object
     * @param psio PSIO object for file I/O
     */
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
         Options& options, std::shared_ptr<PSIO> psio);

    /// Destructor
    ~REKS() override;

    /**
     * @brief Number of electronic states handled by this method.
     * @return Number of states (1 for single-state REKS, N for SA-REKS)
     */
    int n_states() const override { return 1; }

    /**
     * @brief Create a C1 symmetry deep copy of this wavefunction.
     * @param basis Basis set for the new wavefunction
     * @return New REKS wavefunction in C1 symmetry
     */
    std::shared_ptr<REKS> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    // --- Accessors for Testing & Analysis ---
    /// Get FON for orbital r
    double get_n_r() const { return n_r_; }
    /// Get FON for orbital s
    double get_n_s() const { return n_s_; }
    /// Get microstate energy for index L (0-3)
    double get_microstate_energy(int L) const { return E_micro_[L]; }
    /// Get interpolating function value at current FON
    double get_f_value() const { return f_interp(n_r_ / 2.0); }
};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
