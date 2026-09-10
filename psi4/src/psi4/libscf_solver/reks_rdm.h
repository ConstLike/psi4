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
#include "reks_reel.h"

#include <vector>

namespace psi {
namespace reks {
namespace rdm {

/// Diabatic 1-RDM in the active orbital basis, CSR over the configuration rows of
/// si_cassette.K_indices. Row K lists the cells the catalog declares for that row,
/// ascending in (L, pq) with pq = p*n_active + q. Blocks the catalog does not
/// declare are identically zero and carry no entry; a declared cell is kept even
/// when its value vanishes at the current FONs.
/// Diagonal blocks hold p == q only:
///   rho^(KK)_{pp} = sum_M C_M^(K)(fons) * (alpha[p] + beta[p]),
/// p != q inside K = L being identically zero in REKS.
struct DiabaticRdm {
    int                 n_si     = 0;
    int                 n_active = 0;
    std::vector<size_t> row_ptr;   ///< [n_si + 1]
    std::vector<int>    col_L;     ///< [nnz] partner configuration, local index
    std::vector<int>    col_pq;    ///< [nnz] p*n_active + q
    std::vector<double> val;       ///< [nnz]

    size_t nnz() const { return val.size(); }
    bool   empty() const { return n_si < 1; }
};

/// Build the diabatic 1-RDM of si_cassette at the Reel's current FONs.
DiabaticRdm diabatic(const studio::Reel&      reel,
                     const studio::Cassette&  si_cassette);

/// Adiabatic cross-state 1-RDM tensor in the active orbital basis. Layout
///   rho_adiab[((I*n_rep + J)*n_active + p)*n_active + q]
///   rho_adiab[I,J,p,q] = sum_{K,L} c_I^K c_J^L rho_diab[K,L,p,q]
/// where sir.coeffs is row-major c_I^K = coeffs[I*sir.n_states + K]. I,J run over
/// the reported states [0, sir.n_report()); n_rep = sir.n_report() is the block
/// stride. The null-space suffix and the states above the report cap carry no
/// transform and are not stored. K,L run over the full configuration space.
/// Diagonal blocks I == J are the per-state active-block 1-RDM rho^(I).
/// Returns empty vector when sir.n_states < 1.
std::vector<double> adiabatic(const DiabaticRdm&       rho_diab,
                              const studio::Cassette&  si_cassette,
                              const reks::SIResult& sir);

}  // namespace rdm
}  // namespace reks
}  // namespace psi
