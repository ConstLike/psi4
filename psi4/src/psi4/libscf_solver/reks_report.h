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

#include "reks_report_level.h"
#include "reks_si_types.h"
#include "reks_arch.h"
#include "reks_cassette.h"

#include <string>
#include <vector>

namespace psi {
namespace reks {
namespace report {

/// Lowest-energy microstates listed in the run report; the full list needs level 3.
inline constexpr int kMicrostateReportRows = 15;

/// Largest SI matrix printed in the run report; wider ones need level 4. An nsi x nsi
/// matrix costs ceil(nsi/kSIBlockColumns) * (nsi+1) lines.
inline constexpr int kSIMatrixReportDim = 150;

/// REKS input templates for the Python API and the PSithon (.dat) front-ends.
std::string input_usage_examples();

/// Print detail to outfile in a framed banner, then throw PsiException; never
/// returns. show_usage appends input_usage_examples(); set false when the fix is
/// unrelated to the active space / config pool (symmetry, scalar option ranges).
[[noreturn]] void input_error(const std::string& detail, bool show_usage = true);

/// Consume the next table tag "[tagN]" and advance the counter.
std::string next_table_tag();

/// Restart the table-tag counter at 1.
void reset_table_tags();

/// Letter for active-orbital index p (0-based): 'a' + p. Uniform across variants.
std::string orbital_label(int p);

/// Determinant string for global microstate L (catalog then runtime extras),
/// e.g. "|a,b_|". Returns "" when L is outside [0, n_total_microstates).
std::string microstate_label(const studio::Cassette& sa_cassette, int L);

/// Catalog SI config display name for global config index K, else "cfg<K>".
std::string si_config_name(const studio::Cassette& cassette, int K);

/// Compound method label "<K>SI-<L>SA-REKS(N,M)", e.g. "18SI-3SA-REKS(4,4)".
/// n_si = SI cassette dimension; n_sa = SA cassette configs + runtime extra
/// determinants (extras carry an SA weight and count as SA states).
/// Drops the SI prefix when n_si == 0 or n_si == n_sa.
std::string method_label(const studio::Cassette& sa_cassette, int n_si);

/// Multiplicity name for total spin given as 2S (mult = 2S + 1): "singlet",
/// "triplet", ...; "2S+1=<mult>" beyond the named range.
std::string spin_name(int spin2);

/// Spin tag "2S=<2S> (<Multiplicity>)", multiplicity name capitalized,
/// e.g. "2S=2 (Triplet)".
std::string spin_tag(int spin2);

/// Adiabatic SI eigenstate energies (S0..S_{n-1}), one per physical state.
/// No-op when sir.n_physical < 1.
void print_adiabatic_energies(const reks::SIResult& sir);

/// S0 -> Si excitation energies (Ha + eV). No-op when sir.n_physical < 2.
void print_excitation_energies(const reks::SIResult& sir);

/// Microstate energies sorted ascending. E_L is the global microstate energy
/// array (indexed by global L); the k-th row reads E_L[studio_L_indices[k]].
/// Level 2 lists the kMicrostateReportRows lowest and counts the rest, level 3 lists
/// the whole pool.
void print_microstate_energies(const studio::Cassette&      sa_cassette,
                               const std::vector<int>&    studio_L_indices,
                               const std::vector<double>& E_L);

/// Joint Lagrangian eps_pq + (pq|st) table over canonical orbital pairs.
/// ERI cells render as "-" placeholders when pair_eris is empty or its active
/// count differs from cassette.n_active_orbitals() (the eps_pq column always prints).
void print_lagrangian_eri_table(const studio::Cassette&               cassette,
                                const std::vector<double>&          lagrangians,
                                const ActiveEriTiles&               pair_eris);

/// Static masthead: title, method line, authors, institution.
void print_studio_banner();

/// "REKS Studio Setup": variant (N,M,scheme,spin), active space + core, active
/// MO indices. `spin` is the total spin S; `n_core` the doubly-occupied core.
void print_studio_setup(const studio::Cassette& sa, int spin, int n_core,
                        const std::vector<int>& active_mo_indices);

/// "REKS Sector Map": one line per run sector with its 2S, SA config count
/// ("SA-empty" for 0), SI cassette count ("-" for 0), and FON generations.
/// All vectors are indexed by run sector ordinal.
void print_sector_map(const std::vector<int>& sector_2S,
                      const std::vector<int>& sa_counts,
                      const std::vector<int>& si_counts,
                      const std::vector<int>& n_generations);

/// "Pairing schemes": each scheme's geminal orbital pairs followed by the
/// canonical GVB(2,2) wavefunction legend. `show_triplet` adds the ms=0 triplet
/// geminal Psi_To legend.
void print_pairing_schemes(const studio::Cassette& sa, bool show_triplet);

/// "SA pool": the K|config|weight|definition table over the SCF pool, runtime
/// extra determinants as inline rows, weight sum, and the E[SA-REKS] equation.
void print_sa_pool(const studio::Cassette& sa_cassette);

/// "=> SA-REKS Energy Decomposition <=": one row per config (and per runtime extra
/// determinant) with weight, E(config) = sum_L C_L^{K} E_L, and w * E(config);
/// the contributions sum to E[SA-REKS] = e_sa. `e_config` is parallel to
/// sa_cassette.K_indices, `e_extra` parallel to the extra determinants. When e_pen or
/// e_vv10 is nonzero, appends IPR-penalty / VV10 / Total-Energy reconciliation rows.
void print_sa_energy_decomposition(const studio::Cassette&    sa_cassette,
                                   const std::vector<double>& e_config,
                                   const std::vector<double>& e_extra,
                                   double e_sa, double e_pen, double e_vv10);

/// "SI pool": per SI cassette, a header (index, pool dim, method label via the
/// cassette's OWN dim) and its K|config|definition table. A cassette whose config
/// set equals the SA pool is skipped; the whole section is suppressed when no
/// genuine SI cassette remains. `sector_2S` (indexed by run sector) tags each
/// cassette header with its 2S and multiplicity in multi-sector runs.
void print_si_pool(const std::vector<studio::Cassette>& si_cassettes,
                   const studio::Cassette& sa_cassette,
                   const std::vector<int>& sector_2S);

/// Print (nsi x nsi) row-major matrix `M`. When K_indices.size() == nsi, row/col
/// labels are si_config_name(cassette, K_indices[i]); otherwise (including the
/// empty default) labels are si_config_name(cassette, i) for i in [0..nsi).
/// Reports from level 2, or from level 4 when nsi exceeds kSIMatrixReportDim; below
/// the level it needs, a wide matrix leaves one line naming the level that shows it.
void print_si_matrix(const studio::Cassette&      cassette,
                     const char*                title,
                     const std::vector<double>& M,
                     int                        nsi,
                     const std::vector<int>&    K_indices = {});

/// Same table for a matrix kept by its blocks; entries outside a block print as zero.
void print_si_matrix(const studio::Cassette&   cassette,
                     const char*             title,
                     const BlockedMatrix&    M,
                     int                     nsi,
                     const std::vector<int>& K_indices = {});

}  // namespace report
}  // namespace reks
}  // namespace psi
