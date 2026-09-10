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

#include "reks_naturals.h"

#include "reks_math.h"

#include "psi4/libpsi4util/exception.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace psi {
namespace reks {
namespace naturals {

reks::StateNaturalOrbitals
compute(const std::vector<double>& rho_adiab,
        const studio::Cassette&  si_cassette,
        const reks::SIResult& sir) {
    reks::StateNaturalOrbitals out;
    if (sir.n_states < 1) return out;

    const int n_st = sir.n_report();   // reported states; rho_adiab block-row stride
    const int N    = si_cassette.n_active_orbitals();
    const int N2   = N * N;
    if (n_st < 1 || N2 == 0) return out;

    const size_t need = static_cast<size_t>(n_st) * n_st * N2;
    if (rho_adiab.size() != need) {
        throw PSIEXCEPTION(
            "naturals::compute: rho_adiab must be sized n_report^2 * n_active^2 (expected " +
            std::to_string(need) + ", got " + std::to_string(rho_adiab.size()) + ")");
    }

    out.n_states = n_st;
    out.n_active = N;
    out.occupations.assign(static_cast<size_t>(n_st) * N,  0.0);
    out.orbitals  .assign(static_cast<size_t>(n_st) * N2, 0.0);
    out.rdm       .assign(static_cast<size_t>(n_st) * N2, 0.0);

    // rho^(k) = rho_adiab[k,k], per the block layout in rdm::adiabatic.
    for (int k = 0; k < n_st; ++k) {
        const double* src = &rho_adiab[(static_cast<size_t>(k) * n_st + k) * N2];
        std::copy(src, src + N2, &out.rdm[static_cast<size_t>(k) * N2]);
    }

    // Natural-orbital eigenproblem rho^(k) v_alpha = n_alpha v_alpha: n_alpha is the
    // natural occupation, v_alpha the natural orbital (active MO basis).
    // lapack_diagonalize returns n_alpha ascending; reversed here for descending occupations.
    std::vector<double> evals;
    std::vector<double> evecs;
    for (int k = 0; k < n_st; ++k) {
        lapack_diagonalize(&out.rdm[static_cast<size_t>(k) * N2], N, evals, evecs);
        for (int alpha = 0; alpha < N; ++alpha) {
            const int src = N - 1 - alpha;
            out.occupations[static_cast<size_t>(k) * N + alpha] = evals[src];
            for (int p = 0; p < N; ++p) {
                out.orbitals[static_cast<size_t>(k) * N2 + p * N + alpha]
                    = evecs[p * N + src];
            }
        }
    }

    return out;
}

}  // namespace naturals
}  // namespace reks
}  // namespace psi
