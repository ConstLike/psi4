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

#include "reks_si.h"
#include "reks_si_types.h"
#include "reks_math.h"
#include "reks_report_level.h"
#include "reks_studio_eval.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace psi {
namespace reks {
namespace si {

using report::reports;

namespace {

/// Connected components of the SI configs under the pair-overlap graph (edges are pairs
/// with nonzero n_s). Numbered by the first index they claim.
std::vector<std::vector<int>> overlap_blocks(const studio::Cassette& si_cassette) {
    const int n = static_cast<int>(si_cassette.K_indices.size());
    std::vector<int> parent(n);
    for (int i = 0; i < n; ++i) parent[i] = i;
    auto find = [&parent](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };

    const int* row_ptr = si_cassette.pair_row_ptr();
    const int* key_j   = si_cassette.pair_key_j();
    if (row_ptr != nullptr && key_j != nullptr) {
        for (int i = 0; i < n; ++i) {
            const int K_i = si_cassette.K_indices[i];
            for (long long p = row_ptr[K_i]; p < row_ptr[K_i + 1]; ++p) {
                const int j = si_cassette.local_index(key_j[p]);
                if (j < 0 || si_cassette.pair_ref(p).sh->n_s == 0) continue;
                const int ri = find(i), rj = find(j);
                if (ri != rj) parent[ri] = rj;
            }
        }
    }

    std::vector<std::vector<int>> comps;
    std::vector<int> comp_of(n, -1);
    for (int i = 0; i < n; ++i) {
        const int r = find(i);
        if (comp_of[r] < 0) {
            comp_of[r] = static_cast<int>(comps.size());
            comps.emplace_back();
        }
        comps[comp_of[r]].push_back(i);
    }
    return comps;
}

}  // namespace

void build_hamiltonian(const studio::Reel&                 reel,
                            const studio::Cassette&             si_cassette,
                            const studio::SelectorPacks&      packs,
                            const std::vector<double>&        E_L,
                            const std::vector<double>&        lagrangians,
                            double                            hf_exchange_fraction,
                            std::vector<double>&              H,
                            int&                              n_out,
                            const ActiveEriTiles*             active_eri) {
    const int n = static_cast<int>(si_cassette.K_indices.size());
    n_out = n;
    H.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);

    const studio::FONSnapshot& fons = reel.fon_state[si_cassette.sector()];

    // C_K holds C_L(K_i); compute_C_L sizes it on first touch per thread.
    thread_local std::vector<double> C_K;

    // FON factors depend only on the snapshot, not on the config or pair; cache once.
    std::vector<double> fon_cache;
    studio::build_fon_cache(si_cassette.fon_pool(), si_cassette.n_fon_pool(),
                            si_cassette.fon_slot_pool(), si_cassette.fon_slot_idx(),
                            fons, &f_interp, fon_cache);

    double* const H_data = H.data();
    const size_t stride = static_cast<size_t>(n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        const int K_i = si_cassette.K_indices[i];
        si_cassette.compute_C_L(K_i, fons, C_K, fon_cache.data());
        double E_K = 0.0;
        for (int L : si_cassette.config_microstate_writes(K_i))
            E_K += C_K[L] * E_L[L];
        H_data[static_cast<size_t>(i) * stride + i] = E_K;
    }

    // Off-diagonals: uncoupled (i,j) keep 0.0 from assign(); B50 damping uses the
    // diagonal filled above and the pair's has_scaling flag.
    const bool b50_active = hf_exchange_fraction < 1.0 - 1e-10;
    constexpr double b50_p1 = 0.619;
    constexpr double b50_p2 = 3.27;

    const int* row_ptr = si_cassette.pair_row_ptr();
    const int* key_j   = si_cassette.pair_key_j();
    thread_local std::vector<int> win_scratch;

    // Row lengths vary widely across the catalog; dynamic scheduling balances load.
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < n; ++i) {
        const int K_i = si_cassette.K_indices[i];
        double* const H_row_i = H_data + static_cast<size_t>(i) * stride;
        for (long long p = row_ptr[K_i]; p < row_ptr[K_i + 1]; ++p) {
            const int j = si_cassette.local_index(key_j[p]);
            if (j < 0) continue;

            const studio::coupling::PairRef pr = si_cassette.pair_ref(p);
            double H_ij = studio::coupling::eval_pair_H(
                si_cassette.pair_values(p, win_scratch), si_cassette.fon_pool(),
                si_cassette.fon_idx(),
                si_cassette.fon_slot_pool(), si_cassette.fon_slot_idx(),
                si_cassette.coeff_pool(), si_cassette.e_val_pool(),
                si_cassette.fock_idxpart(), si_cassette.fock_valpart(),
                packs.fock.data(),
                si_cassette.eri_idxpart(), si_cassette.eri_valpart(),
                packs.eri.data(),
                si_cassette.lagr_row_pool(), fons, E_L,
                reel.fock_aa_row, reel.fock_aa, si_cassette.n_active_orbitals(),
                *active_eri, lagrangians, fon_cache.data());

            // H_ii, H_jj: written by the earlier loop, already final here.
            if (b50_active && pr.sh->has_scaling) {
                const double dE = H_row_i[i] - H_data[static_cast<size_t>(j) * stride + j];
                const double dE4 = dE * dE * dE * dE;
                H_ij *= b50_p1 * std::exp(-b50_p2 * dE4);
            }

            H_row_i[j] = H_ij;
            H_data[static_cast<size_t>(j) * stride + i] = H_ij;
        }
    }
}

void build_overlap(const studio::Reel&      reel,
                        const studio::Cassette&  si_cassette,
                        BlockedMatrix&         S,
                        int&                   n_out) {
    const int n = static_cast<int>(si_cassette.K_indices.size());
    n_out = n;

    S.reset(n, overlap_blocks(si_cassette));

    const studio::FONSnapshot& fons = reel.fon_state[si_cassette.sector()];

    for (int i = 0; i < n; ++i) {
        const int c = S.block_of[i];
        S.tile(c)[static_cast<size_t>(S.pos_of[i]) * S.dim(c) + S.pos_of[i]] = 1.0;
    }

    const int* row_ptr = si_cassette.pair_row_ptr();
    const int* key_j   = si_cassette.pair_key_j();
    thread_local std::vector<int> win_scratch;

    // FON factors depend only on the snapshot, not on the pair; cache once.
    std::vector<double> fon_cache;
    studio::build_fon_cache(si_cassette.fon_pool(), si_cassette.n_fon_pool(),
                            si_cassette.fon_slot_pool(), si_cassette.fon_slot_idx(),
                            fons, &f_interp, fon_cache);

    // One visit per pair; writes are disjoint across iterations.
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < n; ++i) {
        const int K_i = si_cassette.K_indices[i];
        const int c = S.block_of[i];
        const int b = S.dim(c);
        double* const tile = S.tile(c);
        for (long long p = row_ptr[K_i]; p < row_ptr[K_i + 1]; ++p) {
            const int j = si_cassette.local_index(key_j[p]);
            if (j < 0) continue;
            // n_s == 0: sum is exactly the tile's preset 0.0, skip. Otherwise i and j share
            // this overlap block (overlap_blocks), so pos_of[j] indexes within this tile.
            if (si_cassette.pair_ref(p).sh->n_s == 0) continue;
            const double S_ij = studio::coupling::eval_pair_S(
                si_cassette.pair_values(p, win_scratch).s_idx, si_cassette.fon_pool(),
                si_cassette.fon_idx(), si_cassette.fon_slot_pool(),
                si_cassette.fon_slot_idx(), si_cassette.coeff_pool(),
                si_cassette.s_row_pool(), fons, fon_cache.data());
            tile[static_cast<size_t>(S.pos_of[i]) * b + S.pos_of[j]] = S_ij;
            tile[static_cast<size_t>(S.pos_of[j]) * b + S.pos_of[i]] = S_ij;
        }
    }
}

reks::SIResult diagonalize(const studio::Cassette&  si_cassette,
                                          std::vector<double>    H,
                                          BlockedMatrix          S,
                                          int                    n) {
    reks::SIResult result;
    result.n_states    = n;
    result.report_cap  = si_cassette.report_states;

    // Off-block entries are structural zeros; the sweep covers only the tiles. Block sizes
    // vary widely across the catalog; dynamic scheduling balances load.
    double max_offdiag = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(max : max_offdiag)
    for (int c = 0; c < S.n_blocks(); ++c) {
        const int b = S.dim(c);
        const double* const tile = S.tile(c);
        double blk_max = 0.0;
        for (int r = 0; r < b; ++r)
            for (int s = r + 1; s < b; ++s)
                blk_max = std::max(blk_max, std::abs(tile[static_cast<size_t>(r) * b + s]));
        max_offdiag = std::max(max_offdiag, blk_max);
    }

    // generalized_diagonalize (reks_math.h) builds only the n_show presented columns;
    // `evec_stride` holds however many evecs carries: n for lapack_diagonalize, n_show
    // otherwise.
    const int n_show = result.n_show(n);
    std::vector<double> evals, evecs;
    int evec_stride = n;
    if (max_offdiag < 1e-12) {
        if (reports(2))
            outfile->Printf("\n  SI Eigenproblem: standard (S = I), dimension %d\n", n);
        timer_on("REKS: si_diagonalize");
        lapack_diagonalize(H.data(), n, evals, evecs);
        timer_off("REKS: si_diagonalize");
        result.n_physical = n;  // S == I: full rank, every state physical.
    } else {
        // S is block-diagonal under the catalog's overlap graph: each block diagonalizes
        // independently; the Lowdin transform keeps that block structure.
        const std::vector<std::vector<int>>& comps = S.blocks;
        std::vector<double> s_evals;
        BlockedEigenvectors s_evecs;
        timer_on("REKS: si_diagonalize_overlap");
        lapack_diagonalize_blocked(S, s_evals, s_evecs);
        timer_off("REKS: si_diagonalize_overlap");

        constexpr double sigma = 1e-10;
        constexpr double mu    = 10.0;
        double s_min = *std::min_element(s_evals.begin(), s_evals.end());
        double s_max = *std::max_element(s_evals.begin(), s_evals.end());
        int n_near_null = 0;
        for (int i = 0; i < n; ++i) {
            if (s_evals[i] < sigma) ++n_near_null;
        }

        if (reports(2)) {
            outfile->Printf("\n  SI Eigenproblem: generalized (S != I), dimension %d\n", n);
            outfile->Printf("  Regularization: rank-revealing Lowdin, null-space truncation + "
                            "sentinel penalty (sigma = %.0e, mu = %.0e, E_penalty = %.0f)\n",
                            sigma, mu, mu / sigma);
        }
        if (reports(4)) {
            outfile->Printf("  Overlap S eigenvalues:\n");
            for (int i = 0; i < n; ++i) {
                const char* flag = (s_evals[i] < sigma) ? "  [NULL-SPACE]" : "";
                outfile->Printf("    s[%2d] = %14.6e%s\n", i, s_evals[i], flag);
            }
        }
        if (reports(2)) {
            outfile->Printf("  Physical: %d / %d states, null-space: %d (s < sigma)\n",
                            n - n_near_null, n, n_near_null);
            outfile->Printf("  s_min = %.6e, s_max = %.6e\n", s_min, s_max);
        }

        // n_physical = rank of S (eigenvalues >= sigma).
        result.n_physical = n - n_near_null;

        timer_on("REKS: si_diagonalize");
        generalized_diagonalize(H.data(), n, s_evals, s_evecs, comps, sigma, mu, n_show,
                                evals, evecs);
        timer_off("REKS: si_diagonalize");
        evec_stride = n_show;
    }

    result.energies = std::move(evals);

    // Invariant: physical roots occupy the prefix [0, n_physical); the suffix carries
    // the mu/sigma sentinel.
    {
        constexpr double kSentinelFloor = 1e10;  // physical SI energies ~ -2 Ha; sentinel = mu/sigma = 1e11.
        for (int i = 0; i < result.n_physical; ++i)
            if (!(result.energies[i] < kSentinelFloor))
                throw std::runtime_error("reks::si::diagonalize: physical state " +
                    std::to_string(i) + " carries sentinel energy -- ordering invariant broken");
        for (int i = result.n_physical; i < n; ++i)
            if (!(result.energies[i] > kSentinelFloor))
                throw std::runtime_error("reks::si::diagonalize: null-space state " +
                    std::to_string(i) + " lacks sentinel energy -- ordering invariant broken");
    }

    // evecs holds columns (LAPACK layout); result.coeffs holds states as rows, built up
    // to n_show.
    result.coeffs.resize(static_cast<size_t>(n_show) * n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_show; ++i) {
        double* const c_row_i = result.coeffs.data() + static_cast<size_t>(i) * n;
        const double* const ev = evecs.data() + i;
        for (int j = 0; j < n; ++j) c_row_i[j] = ev[static_cast<size_t>(j) * evec_stride];
    }

    if (reports(2)) {
        // |c^T S c| ~ 1: physical eigenvector; ~ 0: null-space artifact.
        outfile->Printf("  SI eigenvalues and composition (|w_j| > 5%%, w_j = c_j*(Sc)_j):\n");
        // ScAll[i*n+j] = sum_k S[j][k] c_i[k] = (S c_i)_j, for each row i result.coeffs carries.
        // S is block-diagonal: each block is gathered, multiplied by its tile in one GEMM,
        // and scattered back.
        std::vector<double> ScAll(static_cast<size_t>(n_show) * n, 0.0);
        if (n_show > 0) {
            std::vector<double> c_blk, sc_blk;
            for (int c = 0; c < S.n_blocks(); ++c) {
                const std::vector<int>& idx = S.blocks[c];
                const int b = static_cast<int>(idx.size());
                c_blk.assign(static_cast<size_t>(n_show) * b, 0.0);
                sc_blk.assign(static_cast<size_t>(n_show) * b, 0.0);
                for (int i = 0; i < n_show; ++i)
                    for (int r = 0; r < b; ++r)
                        c_blk[static_cast<size_t>(i) * b + r] =
                            result.coeffs[static_cast<size_t>(i) * n + idx[r]];
                C_DGEMM('N', 'T', n_show, b, b, 1.0, c_blk.data(), b,
                        const_cast<double*>(S.tile(c)), b, 0.0, sc_blk.data(), b);
                for (int i = 0; i < n_show; ++i)
                    for (int r = 0; r < b; ++r)
                        ScAll[static_cast<size_t>(i) * n + idx[r]] =
                            sc_blk[static_cast<size_t>(i) * b + r];
            }
        }

        for (int i = 0; i < n_show; ++i) {
            const double* Sc = ScAll.data() + static_cast<size_t>(i) * n;
            const double cSc = C_DDOT(n, &result.coeffs[static_cast<size_t>(i) * n], 1, Sc, 1);

            const char* flag = (std::abs(cSc) < 0.5) ? "  [null-space]" : "";
            outfile->Printf("    E[%2d] = %20.12f%s  ", i, result.energies[i], flag);

            if (std::abs(cSc) < 0.5) {
                outfile->Printf("[c^T S c = %.2e]\n", cSc);
                continue;
            }

            std::vector<std::pair<double, int>> contribs;
            const double* const c_row_i = result.coeffs.data() + static_cast<size_t>(i) * n;
            for (int j = 0; j < n; ++j) {
                double w = c_row_i[j] * Sc[j];
                if (std::abs(w) > 0.05) contribs.push_back({w, j});
            }
            std::sort(contribs.begin(), contribs.end(),
                      [](const auto& a, const auto& b) {
                          return std::abs(a.first) > std::abs(b.first);
                      });
            for (const auto& [w, j] : contribs) {
                const int K = (j >= 0 && j < static_cast<int>(si_cassette.K_indices.size()))
                                  ? si_cassette.K_indices[j] : j;
                const char* name = si_cassette.si_config_name(K);
                std::string fallback;
                if (!name) {
                    fallback = "cfg" + std::to_string(K);
                    name = fallback.c_str();
                }
                outfile->Printf(" %s(%+3.0f%%)", name, w * 100.0);
            }
            outfile->Printf("\n");
        }
    }

    // Neither solver mutates H or S in place.
    result.hamiltonian = std::move(H);
    result.overlap     = std::move(S);

    return result;
}

}  // namespace si
}  // namespace reks
}  // namespace psi
