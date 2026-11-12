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
#include "reks_fon_solver.h"

#include <memory>
#include <vector>

namespace psi {
namespace reks {

/// GVB-PP CI FONs in [eps, 2 - eps] of sector s's active n-geminals; all 1.0 for a diradical, an unresolved
/// structure or a failed eigensolve.
std::vector<double> fon_init_from_ci(int s,
                                     const studio::Cassette&      sa_cassette,
                                     const std::vector<double>& E_L,
                                     int iteration, bool narrate = true);

/// Sector s's SA energy over the FONs of generation gen (x[i] -> snap layer of geminals_active(s, gen)[i]);
/// sa_cassette and E_L are held by reference.
FONObjective
build_geminal_objective(int                                s,
                        int                                gen,
                        const studio::Cassette&              sa_cassette,
                        std::shared_ptr<studio::FONSnapshot> snap,
                        const std::vector<double>&         E_L);

}  // namespace reks
}  // namespace psi
