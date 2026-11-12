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

#include "reks_report_level.h"
#include "reks_si_types.h"
#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_reel.h"
#include "reks_convergence.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

#include <string>
#include <tuple>
#include <vector>

namespace psi {

class Wavefunction;

namespace reks {
namespace report {

/// Lowest-energy microstates listed in the run report; the full list needs level 3.
inline constexpr int kMicrostateReportRows = 15;

/// Largest SI matrix printed in the run report; wider ones need level 4. An nsi x nsi
/// matrix costs ceil(nsi/kSIBlockColumns) * (nsi+1) lines.
inline constexpr int kSIMatrixReportDim = 150;

/// Const snapshot of the SCF state behind one report call.
struct ReportFrame {
    const studio::Catalog&  catalog;
    const studio::Cassette& sa_cassette;
    const studio::Reel&     reel;
    const std::vector<int>& active;      ///< active MO indices (absolute, 0-based)
    int                     Ncore;
    int                     nmo;         ///< irrep 0 (REKS is C1)
    SharedMatrix            Ca;          ///< SCF MO coefficients; reel.Ca_fock_transform is the Fock-basis set
    SharedVector            epsilon;     ///< orbital energies, may be null
    std::string             irrep;       ///< irrep label of irrep 0
};

/// Variable name of run sector s: s = 0 keeps it, s > 0 gives "REKS SECTOR s <name without leading 'REKS '>"
/// ("REKS FON" -> "REKS SECTOR 1 FON").
std::string sector_psivar_name(int s, const std::string& name);

/// Global 0-based energy rank of every MO (irrep 0; REKS is C1).
std::vector<int> mo_energy_rank(const SharedVector& epsilon, int nmo);

/// Energy-rank symmetry label ("5A","10A",...) for each active slot
/// k = 0..n_active_orbitals()-1.
std::vector<std::string> active_orbital_labels(const ReportFrame& frame);

/// "Active Fractionally Occupied" rows: rank, energy and FON per (sector, generation, scheme) of each
/// active orbital; `active` = (energy, rank + 1, slot) sorted by slot.
void print_active_orbitals(const ReportFrame& frame, const std::vector<int>& rank,
                           const std::vector<std::tuple<double, int, int>>& active);

/// Active occupations, orbital energies and geminal pairs of every scheme; one FON row per populated
/// (generation, scheme).
void print_fon_table(const ReportFrame& frame, const studio::Cassette& cassette);

/// Per-pool memory footprint of REKS-owned matrices; `tag` is appended to the header.
void print_memory_footprint(const std::string& tag, const ReportFrame& frame,
                            const Wavefunction& wfn);

/// Always-on MO-consistency block (Ca vs Ca_fock_transform; window around active).
void print_mo_consistency(const ReportFrame& frame);

/// Debug-only F_REKS_MO block (off-diag occ/act -> 0 at convergence).
void print_f_reks_mo_debug(const ReportFrame& frame);

/// Print the SI Hamiltonian + Overlap matrices.
void print_si_hamiltonian_overlap(const ReportFrame& frame, const studio::Cassette& cassette,
                                  const SIResult& sir);

/// `==> METHOD State Properties (K=N) <==` banner, labelled by the SI cassette's own dimension and 2S.
void print_state_properties_header(const ReportFrame& frame, const studio::Cassette& cassette,
                                   int K, int spin2);

/// Track E_SA across iterations: update trace.last_computed_dE, flag sign-flip oscillations.
void print_sa_diagnostics(const ReportFrame& frame, int iteration, ConvergenceTrace& trace);

/// One-line per-iter diagnostic: "[REKS] iter=.. E_SA=.. dE=.." + FON + Fock gaps (grep-friendly).
void print_iteration_diagnostics(const ReportFrame& frame, int iteration, ConvergenceTrace& trace);

/// REKS input templates for the Python API and the PSithon (.dat) front-ends.
std::string input_usage_examples();

/// Print detail in a framed banner (with input_usage_examples() if show_usage) and throw PsiException.
[[noreturn]] void input_error(const std::string& detail, bool show_usage = true);

/// Consume the next table tag "[tagN]" and advance the counter.
std::string next_table_tag();

/// Restart the table-tag counter at 1.
void reset_table_tags();

/// Letter for active-orbital index p (0-based): 'a' + p; "?" if p < 0.
std::string orbital_label(int p);

/// Determinant string for global microstate L (catalog then runtime extras),
/// e.g. "|a,b_|". Returns "" when L is outside [0, n_total_microstates).
std::string microstate_label(const studio::Cassette& sa_cassette, int L);

/// Catalog SI config display name for global config index K.
std::string si_config_name(const studio::Cassette& cassette, int K);

/// "<n_si>SI-<n_sa>SA-REKS(N,M)", n_sa counting extra determinants; no SI prefix if n_si is 0 or n_sa.
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

/// Microstate energies E_L[studio_L_indices[k]] ascending; level 2 lists the kMicrostateReportRows
/// lowest, level 3 all.
void print_microstate_energies(const studio::Cassette&      sa_cassette,
                               const std::vector<int>&    studio_L_indices,
                               const std::vector<double>& E_L);

/// eps_pq and (pq|st) over orbital pairs; pair_eris covers every pair p < q.
void print_lagrangian_eri_table(const studio::Cassette&               cassette,
                                const std::vector<double>&          lagrangians,
                                const ActiveEriTiles&               pair_eris);

/// Static masthead: title, method line, authors, institution.
void print_studio_banner();

/// "REKS Studio Setup": variant (N,M,scheme,spin), active space + core, active
/// MO indices. `spin` is the total spin S; `n_core` the doubly-occupied core.
void print_studio_setup(const studio::Cassette& sa, int spin, int n_core,
                        const std::vector<int>& active_mo_indices);

/// "REKS Sector Map": per run sector its 2S, SA configs, SI cassettes and FON generations.
void print_sector_map(const std::vector<int>& sector_2S,
                      const std::vector<int>& sa_counts,
                      const std::vector<int>& si_counts,
                      const std::vector<int>& n_generations);

/// "Pairing schemes": geminal pairs per scheme and the GVB(2,2) legend (+ triplet Psi_To if show_triplet).
void print_pairing_schemes(const studio::Cassette& sa, bool show_triplet);

/// "SA pool": the K|config|weight|definition table over the SCF pool, runtime
/// extra determinants as inline rows, weight sum, and the E[SA-REKS] equation.
void print_sa_pool(const studio::Cassette& sa_cassette);

/// SA energy decomposition: per config and extra determinant w, E_K = sum_L C_L^K E_L and w E_K,
/// summing to e_sa.
void print_sa_energy_decomposition(const studio::Cassette&    sa_cassette,
                                   const std::vector<double>& e_config,
                                   const std::vector<double>& e_extra,
                                   double e_sa);

/// "SI pool": header and K|config|definition table per SI cassette that differs from the SA pool;
/// multi-sector headers carry 2S and multiplicity.
void print_si_pool(const std::vector<studio::Cassette>& si_cassettes,
                   const studio::Cassette& sa_cassette,
                   const std::vector<int>& sector_2S);

/// nsi x nsi matrix labelled by K_indices (or 0..nsi-1); level 2, or level 4 above kSIMatrixReportDim
/// (one line naming that level below it).
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
