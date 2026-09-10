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

#include "reks_report_level.h"
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_fon_solver.h"

#include <memory>
#include <vector>

namespace psi {
namespace reks {

/// CI-based FON initialization over run sector s's microstates. Every active
/// n-geminal is a 2-electron GVB pair (no SOMO). Returns a vector sized
/// sa_cassette.geminals_active(s, 0).size() (entry i = FON of active geminal
/// sa_cassette.geminals_active(s, 0)[i]) in [eps, 2-eps]; diradical (all
/// couplings weak), unresolved closed-shell GVB structure and a failed
/// eigensolve return all-1.0.
/// narrate = false suppresses report output.
std::vector<double> fon_init_from_ci(int s,
                                     const studio::Cassette&      sa_cassette,
                                     const std::vector<double>& E_L,
                                     int iteration, bool narrate = true);

/// Multi-dim FON objective for run sector s, generation gen (n=0, m=1, u=2, ...).
/// ndim = sa_cassette.geminals_active(s, gen).size(); a trial x[i] writes into
/// snap->layers[gen][sa_cassette.geminals_active(s, gen)[i]]. The energy is sector
/// s's partial SA energy; snap must carry sector s's own FON state.
/// Precondition: sa_cassette, E_L must outlive the returned objective (captured
/// by reference); snap must outlive it (shared).
FONObjective
build_geminal_objective(int                                s,
                        int                                gen,
                        const studio::Cassette&              sa_cassette,
                        std::shared_ptr<studio::FONSnapshot> snap,
                        const std::vector<double>&         E_L);

}  // namespace reks
}  // namespace psi
