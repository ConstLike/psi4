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

#include "reks_rdm.h"

#include "reks_studio_eval.h"
#include "reks_math.h"

#include "psi4/libpsi4util/exception.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>

namespace psi {
namespace reks {
namespace rdm {

DiabaticRdm diabatic(const studio::Reel&      reel,
                     const studio::Cassette&  si_cassette) {
    DiabaticRdm out;
    out.n_si     = static_cast<int>(si_cassette.K_indices.size());
    out.n_active = si_cassette.n_active_orbitals();
    if (out.n_si < 1) return out;

    const int n_si  = out.n_si;
    const int N_act = out.n_active;
    const int n_mic = si_cassette.n_catalog_microstates();

    const studio::FONSnapshot& fons = reel.fon_state[si_cassette.sector()];

    // FON factors depend only on the snapshot, not on the block; cache once.
    std::vector<double> fon_cache;
    studio::build_fon_cache(si_cassette.fon_pool(), si_cassette.n_fon_pool(),
                            si_cassette.fon_slot_pool(), si_cassette.fon_slot_idx(),
                            fons, &f_interp, fon_cache);

    // blk_ptr becomes the CSR row-pointer of a local-row adjacency: for each local
    // row i (K_indices[i]), count the rdm_elems rows whose partner key_j also maps
    // to a local index, once per direction (i->j and j->i).
    const studio::RdmRow* rdm_rows = si_cassette.rdm_elems();
    const int*            row_ptr  = si_cassette.rdm_row_ptr();

    std::vector<size_t> blk_ptr(static_cast<size_t>(n_si) + 1, 0);
    size_t cell_bound = 0;
    if (rdm_rows != nullptr && row_ptr != nullptr) {
        for (int i = 0; i < n_si; ++i) {
            const int K = si_cassette.K_indices[i];
            for (int r = row_ptr[K]; r < row_ptr[K + 1]; ++r) {
                if (rdm_rows[r].key_j == K) continue;
                const int j = si_cassette.local_index(rdm_rows[r].key_j);
                if (j < 0) continue;
                ++blk_ptr[static_cast<size_t>(i) + 1];
                ++blk_ptr[static_cast<size_t>(j) + 1];
                cell_bound += static_cast<size_t>(rdm_rows[r].cell_n);
            }
        }
    }
    for (int i = 0; i < n_si; ++i) blk_ptr[i + 1] += blk_ptr[i];

    // Entry packs local partner j (high bits) with rdm_elems row r (low bits); an
    // ascending sort per row slice orders the partners and carries r along.
    std::vector<uint64_t> blk(blk_ptr[n_si]);
    if (rdm_rows != nullptr && row_ptr != nullptr) {
        std::vector<size_t> cur(blk_ptr.begin(), blk_ptr.end() - 1);
        for (int i = 0; i < n_si; ++i) {
            const int K = si_cassette.K_indices[i];
            for (int r = row_ptr[K]; r < row_ptr[K + 1]; ++r) {
                if (rdm_rows[r].key_j == K) continue;
                const int j = si_cassette.local_index(rdm_rows[r].key_j);
                if (j < 0) continue;
                const uint64_t row = static_cast<uint64_t>(r);
                blk[cur[i]++] = (static_cast<uint64_t>(j) << 32) | row;
                blk[cur[j]++] = (static_cast<uint64_t>(i) << 32) | row;
            }
        }
    }
    for (int i = 0; i < n_si; ++i)
        std::sort(blk.begin() + blk_ptr[i], blk.begin() + blk_ptr[i + 1]);

    // Upper bound on emitted entries: 2 positions per cell (compute_RDM_cells'
    // p != q mirrors) times 2 directions per off-diagonal block (4*cell_bound),
    // plus N_act diagonal entries per local row.
    const size_t bound = 4 * cell_bound + static_cast<size_t>(n_si) * N_act;
    out.row_ptr.assign(static_cast<size_t>(n_si) + 1, 0);
    out.col_L.reserve(bound);
    out.col_pq.reserve(bound);
    out.val.reserve(bound);

    std::vector<double> C_K(static_cast<size_t>(n_mic));
    std::vector<studio::RdmCell> cells;

    // The diagonal block of row i: occupations of the active orbitals in config K,
    // over the microstates config_microstate_writes(K) names.
    auto emit_diagonal = [&](int i) {
        const int K = si_cassette.K_indices[i];
        si_cassette.compute_C_L(K, fons, C_K, fon_cache.data());
        for (int p = 0; p < N_act; ++p) {
            double n_p = 0.0;
            for (int M : si_cassette.config_microstate_writes(K)) {
                const auto& m = si_cassette.microstate(M);
                n_p += C_K[M] * static_cast<double>(m.alpha[p] + m.beta[p]);
            }
            out.col_L.push_back(i);
            out.col_pq.push_back(p * N_act + p);
            out.val.push_back(n_p);
        }
    };

    for (int i = 0; i < n_si; ++i) {
        bool diagonal_done = false;
        for (size_t z = blk_ptr[i]; z < blk_ptr[i + 1]; ++z) {
            const int L = static_cast<int>(blk[z] >> 32);
            const int r = static_cast<int>(blk[z] & 0xffffffffu);
            if (!diagonal_done && L > i) {
                emit_diagonal(i);
                diagonal_done = true;
            }
            const studio::RdmRow& br = rdm_rows[r];
            if (cells.size() < static_cast<size_t>(2 * br.cell_n))
                cells.resize(static_cast<size_t>(2 * br.cell_n));
            const int n_cells = si_cassette.compute_RDM_cells(br.cell_off, br.cell_n, fons,
                                                             N_act, cells.data(),
                                                             fon_cache.data());
            std::sort(cells.begin(), cells.begin() + n_cells);
            for (int c = 0; c < n_cells; ++c) {
                out.col_L.push_back(L);
                out.col_pq.push_back(cells[c].pq);
                out.val.push_back(cells[c].val);
            }
        }
        if (!diagonal_done) emit_diagonal(i);
        out.row_ptr[i + 1] = out.val.size();
    }

    return out;
}

std::vector<double> adiabatic(const DiabaticRdm&       rho_diab,
                              const studio::Cassette&  si_cassette,
                              const reks::SIResult& sir) {
    if (sir.n_states < 1) return {};

    const int    n_si  = sir.n_states;
    const int    n_rep = sir.n_report();
    const int    N_act = si_cassette.n_active_orbitals();
    const size_t N2    = static_cast<size_t>(N_act) * static_cast<size_t>(N_act);

    if (static_cast<int>(si_cassette.K_indices.size()) != n_si) {
        throw PSIEXCEPTION(
            "rdm::adiabatic: si_cassette.K_indices size (" +
            std::to_string(si_cassette.K_indices.size()) +
            ") must match sir.n_states (" + std::to_string(n_si) + ")");
    }

    // Precondition: n_physical <= n_states; coeffs holds the presented states as
    // rows 0..n_show-1, and n_rep <= n_show bounds the read below.
    const size_t n_coeff_rows = static_cast<size_t>(sir.n_show(n_si));
    if (sir.n_physical > n_si ||
        sir.coeffs.size() != n_coeff_rows * static_cast<size_t>(n_si)) {
        throw PSIEXCEPTION(
            "rdm::adiabatic: expects n_physical <= n_states and coeffs sized "
            "n_show(n_states)*n_states (got n_physical=" + std::to_string(sir.n_physical) +
            ", n_states=" + std::to_string(n_si) +
            ", coeffs=" + std::to_string(sir.coeffs.size()) + ")");
    }

    // Factors the two-sided transform (reks_rdm.h) into two contractions over the
    // flattened (L,pq) / pq trailing axes:
    //   T[L][I][pq]         = sum_K c[I][K] rho_diab[K][L][pq]   (over the CSR rows of rho_diab)
    //   rho_adiab[J][I][pq] = sum_L c[J][L] T[L][I][pq]          (one GEMM, T as n_si x n_rep*N2)
    // rho_diab[K,L] == rho_diab[L,K], so the GEMM's [J][I] output equals the wanted [I][J].
    const size_t out_row = static_cast<size_t>(n_rep) * N2;
    std::vector<double> rho_adiab(static_cast<size_t>(n_rep) * out_row, 0.0);
    if (n_rep < 1 || N2 == 0) return rho_adiab;

    double* coeffs = const_cast<double*>(sir.coeffs.data());  // C_DGEMM needs non-const; not modified.

    if (rho_diab.n_si != n_si || rho_diab.n_active != N_act) {
        throw PSIEXCEPTION(
            "rdm::adiabatic: rho_diab is (" + std::to_string(rho_diab.n_si) + ", " +
            std::to_string(rho_diab.n_active) + ") but the cassette is (" +
            std::to_string(n_si) + ", " + std::to_string(N_act) + ")");
    }

    // State I writes only the [I*N2, (I+1)*N2) stripe of every T row: the loop
    // over I below is race-free.
    std::vector<double> T(static_cast<size_t>(n_si) * out_row, 0.0);
#pragma omp parallel for schedule(static)
    for (int I = 0; I < n_rep; ++I) {
        const double* ci = coeffs + static_cast<size_t>(I) * n_si;
        const size_t base = static_cast<size_t>(I) * N2;
        for (int K = 0; K < n_si; ++K) {
            const double a = ci[K];
            if (a == 0.0) continue;
            for (size_t z = rho_diab.row_ptr[K]; z < rho_diab.row_ptr[K + 1]; ++z)
                T[static_cast<size_t>(rho_diab.col_L[z]) * out_row + base +
                  rho_diab.col_pq[z]] += a * rho_diab.val[z];
        }
    }

    C_DGEMM('N', 'N', n_rep, static_cast<int>(out_row), n_si, 1.0, coeffs, n_si,
            T.data(), static_cast<int>(out_row), 0.0, rho_adiab.data(),
            static_cast<int>(out_row));
    return rho_adiab;
}

}  // namespace rdm
}  // namespace reks
}  // namespace psi
