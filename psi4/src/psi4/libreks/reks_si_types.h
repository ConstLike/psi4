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

#include <algorithm>
#include <utility>
#include <vector>

namespace psi {
namespace reks {

/// Max columns per block when printing wide SI matrices.
inline constexpr int kSIBlockColumns = 18;

/// Position of pair (I < J) in the upper-triangular packing over n states.
inline int si_pair_index(int I, int J, int n) { return I * n - I * (I + 1) / 2 + (J - I - 1); }

/// Symmetric matrix block-diagonal under a permutation: block c = b_c x b_c row-major tile at data[off[c]]
/// over indices blocks[c], zero elsewhere.
struct BlockedMatrix {
    int n = 0;
    std::vector<std::vector<int>> blocks;  ///< partition of [0, n), ascending inside each
    std::vector<int> block_of;             ///< n entries: index -> its block
    std::vector<int> pos_of;               ///< n entries: index -> its position in that block
    std::vector<double> data;              ///< tiles back to back
    std::vector<size_t> off;               ///< n_blocks + 1 offsets into data

    int n_blocks() const { return static_cast<int>(blocks.size()); }
    int dim(int c) const { return static_cast<int>(blocks[c].size()); }
    double* tile(int c) { return data.data() + off[c]; }
    const double* tile(int c) const { return data.data() + off[c]; }
    bool empty() const { return n == 0; }

    /// Element (i, j); zero whenever the two indices sit in different blocks.
    double at(int i, int j) const {
        const int c = block_of[i];
        if (block_of[j] != c) return 0.0;
        return data[off[c] + static_cast<size_t>(pos_of[i]) * blocks[c].size() + pos_of[j]];
    }

    /// Lay the tiles out over `parts` and zero them.
    void reset(int dim_n, std::vector<std::vector<int>> parts) {
        n = dim_n;
        blocks = std::move(parts);
        block_of.assign(n, 0);
        pos_of.assign(n, 0);
        off.assign(blocks.size() + 1, 0);
        for (size_t c = 0; c < blocks.size(); ++c) {
            off[c + 1] = off[c] + blocks[c].size() * blocks[c].size();
            for (size_t r = 0; r < blocks[c].size(); ++r) {
                block_of[blocks[c][r]] = static_cast<int>(c);
                pos_of[blocks[c][r]]   = static_cast<int>(r);
            }
        }
        data.assign(off.back(), 0.0);
    }
};

/// SI Hamiltonian result. Square n_states x n_states, one row/col per SI config.
struct SIResult {
    int n_states = 0;                ///< Config-space dimension (= number of SI configs); stride of coeffs/H/S
    int n_physical = 0;              ///< Number of physical eigenstates (rank of S); prefix [0, n_physical)
                                     ///< is physical, rest are null-space sentinels; equals n_states iff S == I.
    std::vector<double> energies;    ///< Eigenvalues (adiabatic state energies) [n_states]
    std::vector<double> coeffs;      ///< Eigenvectors, row-major [n_show(n_states) * n_states]
                                     ///< (row = presented state, col = config)
    std::vector<double> hamiltonian; ///< H matrix, row-major [n_states * n_states]
    BlockedMatrix overlap;           ///< S; nonzero only inside the catalog's overlap blocks

    int report_cap = 0;              ///< Cap on presented states; 0 = no cap

    /// Number of states to present out of `full`, honouring the cap.
    int n_show(int full) const { return report_cap > 0 ? std::min(full, report_cap) : full; }

    /// State-axis size: the physical roots under the cap (n_states is the config axis).
    int n_report() const { return n_show(n_physical); }
};

/// Natural orbitals of the reported states (active basis): occupations[k N + a] descending,
/// orbitals[k N^2 + p N + a] (column a), rdm[k N^2 + p N + q] = rho^(k).
struct StateNaturalOrbitals {
    int n_states = 0;
    int n_active = 0;
    std::vector<double> occupations;
    std::vector<double> orbitals;
    std::vector<double> rdm;
};

/// Active (ij|kl) as one N x N row-major tile per registered pair k <= l.
class ActiveEriTiles {
    int N_ = 0;
    /// Canonical (k<=l) keys, lexicographically sorted; parallel to the tiles.
    std::vector<std::pair<int,int>> keys_;
    /// Tile p occupies [p*N*N, (p+1)*N*N), row-major N x N.
    std::vector<double> tiles_;
    /// tile_of_[k*N + l] == tile index of the pair, symmetric in (k,l); -1 when unregistered.
    std::vector<int> tile_of_;

    /// First element of the tile covering (k,l).
    size_t tile_base_(int k, int l) const {
        return static_cast<size_t>(tile_of_[static_cast<size_t>(k) * N_ + l]) *
               static_cast<size_t>(N_) * N_;
    }

   public:
    ActiveEriTiles() = default;
    ActiveEriTiles(int N, std::vector<std::pair<int,int>> canonical_keys);

    bool empty() const { return tiles_.empty(); }

    /// (ij|kl), k and l in either order; unchecked.
    [[nodiscard]] double at(int i, int j, int k, int l) const {
        return tiles_[tile_base_(k, l) + static_cast<size_t>(i) * N_ + j];
    }

    /// Mutable tile of (k,l), N*N row-major.
    double* tile_data(int k, int l) { return tiles_.data() + tile_base_(k, l); }
};

}  // namespace reks
}  // namespace psi
