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

#include "reks_si_types.h"

#include <stdexcept>
#include <string>

namespace psi {
namespace reks {

ActiveEriTiles::ActiveEriTiles(int N, std::vector<std::pair<int,int>> canonical_keys)
    : N_(N), keys_(std::move(canonical_keys)) {
    if (N_ < 0) {
        throw std::invalid_argument("ActiveEriTiles: N must be non-negative");
    }
    for (size_t p = 0; p < keys_.size(); ++p) {
        const auto& kv = keys_[p];
        if (kv.first < 0 || kv.second < kv.first || kv.second >= N_) {
            throw std::invalid_argument(
                "ActiveEriTiles: invalid pair (" + std::to_string(kv.first) + "," +
                std::to_string(kv.second) + ") for N=" + std::to_string(N_) +
                "; expected 0 <= k <= l < N");
        }
        if (p > 0 && !(keys_[p - 1] < kv)) {
            throw std::invalid_argument(
                "ActiveEriTiles: keys must be lexicographically strictly sorted "
                "(no duplicates)");
        }
    }
    tiles_.assign(keys_.size() * static_cast<size_t>(N_) * N_, 0.0);

    tile_of_.assign(static_cast<size_t>(N_) * N_, -1);
    for (size_t p = 0; p < keys_.size(); ++p) {
        const auto& kv = keys_[p];
        tile_of_[static_cast<size_t>(kv.first) * N_ + kv.second] = static_cast<int>(p);
        tile_of_[static_cast<size_t>(kv.second) * N_ + kv.first] = static_cast<int>(p);
    }
}

}  // namespace reks
}  // namespace psi
