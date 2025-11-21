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

namespace psi {
namespace scf {

/**
 * @class REKS
 * @brief Restricted Ensemble Kohn-Sham (REKS) method implementation.
 *
 * REKS is a multi-configurational DFT method that treats static correlation
 * through ensemble averaging of fractionally occupied Kohn-Sham determinants.
 * This class currently inherits from RHF and delegates all functionality to it,
 * serving as a foundation for future REKS-specific implementations.
 *
 * Design Philosophy:
 * - Clean separation from RHF to allow independent REKS development
 * - All RHF methods are available but can be overridden for REKS-specific behavior
 * - Prepared for multi-state extensions (SA-REKS)
 *
 * Future Extensions:
 * - Fractional occupation numbers for active orbitals
 * - Ensemble energy expression with state-averaging
 * - Analytical gradients for geometry optimization
 * - State-specific properties
 *
 * References:
 * - Filatov, M.; Shaik, S. Chem. Phys. Lett. 1999, 304, 429-437
 * - Filatov, M. WIREs Comput. Mol. Sci. 2015, 5, 146-167
 */
class REKS : public RHF {
   protected:
    /// REKS-specific initialization
    void reks_common_init();

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
     *
     * Currently returns 1 (single-state), will be extended for SA-REKS.
     */
    int n_states() const override { return 1; }

    /**
     * @brief Create a C1 symmetry deep copy of this wavefunction.
     * @param basis Basis set for the new wavefunction
     * @return New REKS wavefunction in C1 symmetry
     */
    std::shared_ptr<REKS> c1_deep_copy(std::shared_ptr<BasisSet> basis);
};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
