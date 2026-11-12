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

#pragma once

#include "reks_arch.h"
#include "reks_cassette.h"

#include <vector>

namespace psi {
namespace reks {
namespace fon_relax {

/// Per-SI-cassette FON output.
/// layers[gen] is the relaxed FON vector of generation gen (n, m, u, ...).
struct Result {
    std::array<std::vector<studio::GeminalFON>, studio::kMaxGen> layers;
};

/// Box-constrained FON minimum of the geminals the SI cassette adds to sector sa_sector's SA set; SA geminals
/// keep their converged FONs.
Result optimize(const studio::Cassette&      si_cassette,
                           const studio::Cassette&      sa_cassette,
                           int                        sa_sector,
                           const studio::FONSnapshot&   fons_initial,
                           const std::vector<double>& E_L);

}  // namespace fon_relax
}  // namespace reks
}  // namespace psi
