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

namespace psi {
namespace reks {
namespace report {

/// REKS_REPORT_LEVEL, the verbosity of the REKS report: 1 fatal diagnostics only,
/// 2 the run report, 3 extended plus the full microstate list, 4 per-iteration traces,
/// pool maps and SI matrices wider than kSIMatrixReportDim, 5 matrices over the config,
/// microstate and orbital axes.
inline constexpr int kMinLevel = 1;
inline constexpr int kMaxLevel = 5;

namespace detail {
/// Set once before any parallel region reads it; no synchronization needed.
inline int current = 2;
}  // namespace detail

inline void set_level(int lvl) { detail::current = lvl; }

inline int level() { return detail::current; }

inline bool reports(int lvl) { return detail::current >= lvl; }

}  // namespace report
}  // namespace reks
}  // namespace psi
