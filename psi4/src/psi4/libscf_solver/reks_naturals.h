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

#pragma once

#include "reks_si_types.h"
#include "reks_cassette.h"

namespace psi {
namespace reks {
namespace naturals {

/// Natural occupations (descending, in [0,2]) and orbitals (columns, active-MO
/// basis) of each reported state k in [0, sir.n_report()), diagonalizing rho^(k),
/// the diagonal block of rho_adiab (rdm::adiabatic).
/// Returns an empty result when sir.n_states < 1, n_report() < 1, or
/// si_cassette.n_active_orbitals() == 0.
/// Throws unless rho_adiab is sized n_report^2 * n_active^2.
reks::StateNaturalOrbitals
compute(const std::vector<double>& rho_adiab,
        const studio::Cassette&  si_cassette,
        const reks::SIResult& sir);

}  // namespace naturals
}  // namespace reks
}  // namespace psi
