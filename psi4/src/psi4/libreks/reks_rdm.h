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
#include "reks_cassette.h"
#include "reks_reel.h"

#include <vector>

namespace psi {
namespace reks {
namespace rdm {

/// Diabatic active 1-RDM, CSR over config rows K: the catalog-declared cells ascending in (L, p n_active + q);
/// diagonal blocks only p == q, rho^(KK)_pp = sum_M C_M^(K)(fons) (alpha_p + beta_p).
struct DiabaticRdm {
    int                 n_si     = 0;
    int                 n_active = 0;
    std::vector<size_t> row_ptr;   ///< [n_si + 1]
    std::vector<int>    col_L;     ///< [nnz] partner configuration, local index
    std::vector<int>    col_pq;    ///< [nnz] p*n_active + q
    std::vector<double> val;       ///< [nnz]
};

/// Build the diabatic 1-RDM of si_cassette at the Reel's current FONs.
DiabaticRdm diabatic(const studio::Reel&      reel,
                     const studio::Cassette&  si_cassette);

/// rho_adiab[((I n_rep + J) n_active + p) n_active + q] = sum_{K,L} c_I^K c_J^L rho_diab[K,L,p,q] over the
/// reported states I, J < n_rep (I = J: state 1-RDM); empty if sir.n_states < 1.
std::vector<double> adiabatic(const DiabaticRdm&       rho_diab,
                              const studio::Cassette&  si_cassette,
                              const reks::SIResult& sir);

}  // namespace rdm
}  // namespace reks
}  // namespace psi
