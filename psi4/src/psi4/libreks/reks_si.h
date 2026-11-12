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

#include "reks_si_types.h"
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_reel.h"

#include <vector>

namespace psi {
namespace reks {
namespace si {

/// SI H: H_ii = sum_L C_L(K_i) E_L, H_ij from the pair rows with reel.fock_aa and active_eri (packs = this
/// cassette's selector rows); B50 scaling when hf_exchange_fraction < 1.
void build_hamiltonian(const studio::Reel&                 reel,
                            const studio::Cassette&             si_cassette,
                            const studio::SelectorPacks&      packs,
                            const std::vector<double>&        E_L,
                            const std::vector<double>&        lagrangians,
                            double                            hf_exchange_fraction,
                            std::vector<double>&              H,
                            int&                              n_out,
                            const ActiveEriTiles*             active_eri = nullptr);

/// SI overlap: S_ii = 1, S_ij from the pure-FON s channel of the pair (zero without it).
void build_overlap(const studio::Reel&      reel,
                        const studio::Cassette&  si_cassette,
                        BlockedMatrix&         S);

/// H c = E S c: standard eigensolver if S = I (to 1e-12), else canonical orthogonalization; H and S move into
/// the result.
reks::SIResult diagonalize(const studio::Cassette&  si_cassette,
                                          std::vector<double>    H,
                                          BlockedMatrix          S,
                                          int                    n);

}  // namespace si
}  // namespace reks
}  // namespace psi
