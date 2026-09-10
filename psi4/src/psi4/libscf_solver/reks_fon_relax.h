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

#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_report_level.h"

#include <vector>

namespace psi {
namespace reks {
namespace fon_relax {

/// Per-SI-cassette FON output.
/// layers[gen] is the relaxed FON vector of generation gen (n, m, u, ...).
struct Result {
    std::array<std::vector<studio::GeminalFON>, studio::kMaxGen> layers;
};

/// State-specific FON relaxation for one si_cassette (all its configs in run
/// sector sa_sector). Relaxes every geminal this SI cassette adds on top of the
/// SCF pool for that sector (SI.active \ SA.geminals_active[sa_sector]) to a
/// box-constrained local energy minimum; geminals in that sector's SA active
/// set keep their SA-converged FON (re-relaxing destroys the SA minimum).
/// sa_sector: run sector ordinal whose SA active set is frozen.
/// E_L: per-microstate energy vector, global-L indexed.
Result optimize(const studio::Cassette&      si_cassette,
                           const studio::Cassette&      sa_cassette,
                           int                        sa_sector,
                           const studio::FONSnapshot&   fons_initial,
                           const std::vector<double>& E_L);

}  // namespace fon_relax
}  // namespace reks
}  // namespace psi
