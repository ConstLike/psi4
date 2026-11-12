/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef REKS_PROPERTIES_H
#define REKS_PROPERTIES_H

/// SI state properties: dipoles, natural orbitals and their wavefunction variables.

#include "reks_report.h"
#include "reks_si_types.h"
#include "reks_cassette.h"
#include "reks_rdm.h"

#include <array>
#include <memory>
#include <vector>

namespace psi {

class Wavefunction;
class Molecule;
class BasisSet;

namespace reks {
namespace properties {

/// Result of state-property computation: permanent/transition dipoles + osc strengths.
struct SIProperties {
    /// Every member is sized over the reported states n_rep = SIResult::n_report():
    /// n_pairs = n_rep*(n_rep-1)/2, and si_pair_index takes n_rep as its bound.
    std::vector<std::array<double, 3>> perm_dipole;   ///< [n_rep], atomic units
    std::vector<std::array<double, 3>> trans_dipole;  ///< [n_pairs], atomic units
    std::vector<double> osc_str;                      ///< [n_pairs]
    /// Adiabatic 1-RDMs, block (I n_rep + J) N_act^2: state density (I = J), transition (I != J); empty
    /// without a diabatic RDM.
    std::vector<double> rho_adiab;
};

/// Cassette-independent dipole data: nuclear + core dipole and the active-block
/// MO dipole over the frozen orbitals (Ca_fock_transform).
struct DipoleIntegralCache {
    std::array<double, 3> mu_nuc  = {0.0, 0.0, 0.0};
    std::array<double, 3> mu_core = {0.0, 0.0, 0.0};
    std::array<std::vector<double>, 3> d_act;   ///< [axis][p*N_act+q]
};

/// Build the cassette-independent dipole integrals.
DipoleIntegralCache build_dipole_integral_cache(const report::ReportFrame& frame,
                                                std::shared_ptr<Molecule> molecule,
                                                std::shared_ptr<BasisSet> basisset);

/// Permanent + transition dipoles + oscillator strengths of the SI states in sir.
SIProperties compute_state_properties(const report::ReportFrame& frame,
                                      const studio::Cassette& cassette,
                                      const SIResult& sir,
                                      const rdm::DiabaticRdm& rho_diab,
                                      const DipoleIntegralCache& dip);

/// Prints the transition-dipole and oscillator-strength tables.
void print_state_dipoles(const SIProperties& props, const SIResult& sir);

/// SI wfn variables "<BASE> K=<K>" (bare too if primary_alias), sector-qualified for s > 0; config-axis
/// arrays only from REKS_REPORT_LEVEL 3, without alias.
void register_si_wfn_variables(Wavefunction& wfn, const report::ReportFrame& frame,
                               int K, int s,
                               const SIResult& sir,
                               const SIProperties& props,
                               bool primary_alias,
                               const rdm::DiabaticRdm& rho_diab);

/// Empty when SI data is missing or the catalog carries no diabatic RDM.
StateNaturalOrbitals compute_state_natural(const studio::Catalog& catalog,
                                           const studio::Cassette& cassette,
                                           const SIResult& sir,
                                           const std::vector<double>& rho_adiab);

/// Tables A (full SI eigenvector) + B (per-state FONs / NOs) + integrity block.
/// `spin2` is the cassette sector's 2S.
void print_state_natural(const report::ReportFrame& frame,
                         const StateNaturalOrbitals& nos, int K,
                         const studio::Cassette& cassette,
                         const SIResult& sir,
                         const rdm::DiabaticRdm& rho_diab, int spin2);

/// Natural-orbital wfn variables " K=<K>" and " STATE=<k>" of the reported states (bare too if
/// primary_alias), sector-qualified for s > 0.
void register_natural_wfn_variables(Wavefunction& wfn, const report::ReportFrame& frame,
                                    int K, int s,
                                    const StateNaturalOrbitals& nos,
                                    bool primary_alias);

}  // namespace properties
}  // namespace reks
}  // namespace psi

#endif  // REKS_PROPERTIES_H
