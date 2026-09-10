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
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_reel.h"
#include "psi4/libmints/typedefs.h"

#include <vector>

namespace psi {
namespace reks {
namespace si {

/// Build the n x n SI Hamiltonian over si_cassette.K_indices via the cassette's
/// capability formulas. Diagonal H[i,i] = sum_L C_L(K_i) * E_L[L];
/// off-diagonal: coupling tables from the cassette's pair rows and fon_pool() evaluated via
/// reel.fock_aa (the active block of the per-microstate Fock spin difference) and
/// active_eri; zero when no pair entry.
/// `packs` carries the fock/eri selector rows this cassette reaches, in its sector's
/// numbering (decode_selector_packs over selector_demand of the same cassette).
/// B50 scaling fires when hf_exchange_fraction is below 1. H is resized to n*n;
/// n_out returns si_cassette.K_indices.size().
void build_hamiltonian(const studio::Reel&                 reel,
                            const studio::Cassette&             si_cassette,
                            const studio::SelectorPacks&      packs,
                            const std::vector<double>&        E_L,
                            const std::vector<double>&        lagrangians,
                            double                            hf_exchange_fraction,
                            std::vector<double>&              H,
                            int&                              n_out,
                            const ActiveEriTiles*             active_eri = nullptr);

/// Build the SI overlap over si_cassette.K_indices from the pure-FON s channel of each
/// coupled pair; S[i,i] = 1 by normalization. Only pairs carrying that channel are
/// nonzero. n_out returns si_cassette.K_indices.size().
void build_overlap(const studio::Reel&      reel,
                        const studio::Cassette&  si_cassette,
                        BlockedMatrix&         S,
                        int&                   n_out);

/// Solve the n x n SI eigenproblem H c = E S c. Branches to standard LAPACK
/// diagonalization when S == I (off-diagonal magnitude below 1e-12) and to
/// rank-revealing Lowdin generalized diagonalization otherwise.
///
/// H and S are consumed and handed to the returned SIResult in .hamiltonian / .overlap.
reks::SIResult diagonalize(const studio::Cassette&  si_cassette,
                                          std::vector<double>    H,
                                          BlockedMatrix          S,
                                          int                    n);

}  // namespace si
}  // namespace reks
}  // namespace psi
