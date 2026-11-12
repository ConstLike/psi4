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

#ifndef REKS_FON_STEP_H
#define REKS_FON_STEP_H

/// @file reks_fon_step.h
/// @brief FON Newton relaxation step: multi-start seed, per-block joint NR, pair swaps.

#include <array>
#include <vector>

#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_fon_solver.h"
#include "reks_reel.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace psi {
namespace reks {

/// Inputs of one FonStep::relax call; the references are valid for that call only.
struct FonStepContext {
    int iteration;
    const std::vector<int>& active;                          ///< active MO indices (absolute, 0-based)
    const studio::Cassette& sa_cassette;
    studio::Reel& reel;                                      ///< E_L read; fon_state written
    const std::vector<studio::FonBlock>& fon_blocks;         ///< (s, gen) blocks; debug snapshot of all FONs
    const std::array<double, studio::kMaxGen>& fon_lower;    ///< REKS_<G>_FON per generation, < 0 = unset
    SharedMatrix Ca;                                         ///< MO coefficients; pair columns swapped in place
    SharedVector epsilon;                                    ///< orbital energies; null skips the fgap inv-swap
    bool swapped = false;                                    ///< OUT: an orbital-column swap occurred
};

/// FON relaxation: first NR per sector from {(1,...,1), CI} keeping the lowest interior E_SA, then NR per
/// (sector, generation) from the current snapshot.
class FonStep {
   public:
    /// Micro-solver settings and the boundary skin of the gen-0 multi-start scan.
    void configure(const MicroSolverConfig& cfg, double boundary_skin);

    /// One relaxation of the FON snapshot from the computed microstate energies.
    void relax(FonStepContext& ctx);

   private:
    bool initialized_ = false;         ///< First relax() call seeds the FON

    /// Settings shared by both relax() micro-solves; block-local overrides are the
    /// per-generation lower bound and the history switch.
    MicroSolverConfig cfg_;
    /// Distance from a FON bound inside which the gen-0 multi-start scan calls a minimum
    /// boundary-trapped and prefers the interior one.
    double boundary_skin_ = 0.05;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_FON_STEP_H
