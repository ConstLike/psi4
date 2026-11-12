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

#include "reks_cassette.h"

#include "reks_studio_eval.h"
#include "reks_math.h"
#include "psi4/libpsi4util/process.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <utility>

namespace psi {
namespace reks {
namespace studio {

namespace {

template <class T>
std::vector<T> sorted_unique_from_set(const std::set<T>& s) {
    return std::vector<T>(s.begin(), s.end());
}

// Ids stamped in a dense byte map, read out ascending; ids are bounded by the
// catalog (geminals by n_geminals, deltas by n_delta_specs).
std::vector<int> sorted_unique_from_stamps(const std::vector<uint8_t>& stamps) {
    std::vector<int> out;
    for (size_t i = 0; i < stamps.size(); ++i)
        if (stamps[i]) out.push_back(static_cast<int>(i));
    return out;
}

void merge_stamps(std::vector<uint8_t>& dst, const std::vector<uint8_t>& src) {
    for (size_t i = 0; i < dst.size(); ++i) dst[i] |= src[i];
}

}  // namespace

Cassette make_cassette(Catalog                 catalog,
                       std::vector<int>          K_indices,
                       std::vector<double>       K_weights,
                       CassetteRole              role,
                       std::vector<Microstate>   extra_microstates,
                       std::vector<double>       extra_weights) {
    Cassette c;
    c.catalog_   = catalog;
    c.K_indices = std::move(K_indices);
    c.K_weights = std::move(K_weights);
    c.extra_microstates = std::move(extra_microstates);
    c.extra_weights     = std::move(extra_weights);

    if (!c.K_indices.empty()) {
        const auto mm = std::minmax_element(c.K_indices.begin(), c.K_indices.end());
        c.K_lo_ = *mm.first;
        c.local_of_.assign(static_cast<size_t>(*mm.second - c.K_lo_ + 1), -1);
        for (int i = 0; i < static_cast<int>(c.K_indices.size()); ++i)
            c.local_of_[c.K_indices[i] - c.K_lo_] = i;
    }

    // ms_seen[L] stamps microstate L active, dense in [0, n_microstates); gem_seen[gen][g]
    // and delta_seen[d] stamp the same way, dense in [0, n_geminals) and [0, n_delta_specs).
    const size_t n_ms_stamps    = static_cast<size_t>(catalog.n_microstates());
    const size_t n_gem_stamps   = static_cast<size_t>(catalog.data->n_geminals);
    const size_t n_delta_stamps = static_cast<size_t>(catalog.n_delta_specs());
    std::vector<uint8_t> ms_seen(n_ms_stamps, 0);
    std::array<std::vector<uint8_t>, kMaxGen> gem_seen;
    for (auto& stamps : gem_seen) stamps.assign(n_gem_stamps, 0);
    std::vector<uint8_t> delta_seen(n_delta_stamps, 0);

    // Diagonal contributions: every K in K_indices participates as H[K,K].
    for (int K : c.K_indices) {
        const DiagDepsView d = catalog.diag_deps_view(K);
        for (int L : d.microstate_writes) ms_seen[L] = 1;
        for (int gen = 0; gen < kMaxGen; ++gen)
            for (int g : d.geminal_reads[gen]) gem_seen[gen][g] = 1;
    }

    // SI: microstates and deltas of every pair H[K_i,K_j] in the cassette, one pass over the CSR rows of the
    // lower keys; all configs belong to the bound sector.
    if (role == CassetteRole::SI) {
        const int n = static_cast<int>(c.K_indices.size());

        // Rows split over threads (dynamic: lengths vary > 10x), per-thread maps folded at the end.
        const int nthread = Process::environment.get_n_threads();
#pragma omp parallel num_threads(nthread)
        {
            std::vector<uint8_t> ms_local(n_ms_stamps, 0);
            std::array<std::vector<uint8_t>, kMaxGen> gem_local;
            for (auto& stamps : gem_local) stamps.assign(n_gem_stamps, 0);
            std::vector<uint8_t> delta_local(n_delta_stamps, 0);
            // Dependency windows are varint-coded in the blob; one buffer serves every pair.
            std::vector<int> deps_scratch;
            // Selector reach of this thread's pairs, folded out of the same walk.
            long long fock_local = 0, eri_local = 0;

#pragma omp for schedule(dynamic, 8)
            for (int i = 0; i < n; ++i) {
                const auto [q_lo, q_hi] = catalog.pair_row(c.K_indices[i]);
                for (long long q = q_lo; q < q_hi; ++q) {
                    if (c.local_index(catalog.pair_partner(q)) < 0) continue;
                    const PairDepsView pd =
                        catalog.pair_deps(q, deps_scratch, fock_local, eri_local);
                    for (int L : pd.fock_footprint) ms_local[L] = 1;
                    for (int gen = 0; gen < kMaxGen; ++gen)
                        for (int g : pd.geminal_reads[gen]) gem_local[gen][g] = 1;
                    for (int d : pd.delta_reads) delta_local[d] = 1;
                }
            }

#pragma omp critical
            {
                merge_stamps(ms_seen, ms_local);
                for (int gen = 0; gen < kMaxGen; ++gen)
                    merge_stamps(gem_seen[gen], gem_local[gen]);
                merge_stamps(delta_seen, delta_local);
                c.selector_demand_.fock_rows =
                    std::max(c.selector_demand_.fock_rows, fock_local);
                c.selector_demand_.eri_rows =
                    std::max(c.selector_demand_.eri_rows, eri_local);
            }
        }
    }

    c.microstates_active.clear();
    for (size_t L = 0; L < ms_seen.size(); ++L)
        if (ms_seen[L]) c.microstates_active.push_back(static_cast<int>(L));
    for (int gen = 0; gen < kMaxGen; ++gen)
        c.geminals_active_pool_[gen] = sorted_unique_from_stamps(gem_seen[gen]);
    c.deltas_active      = sorted_unique_from_stamps(delta_seen);

    // Extra L = n_catalog + e exceeds every catalog L, so push_back keeps
    // microstates_active sorted-unique without re-sorting.
    const int n_catalog = catalog.n_microstates();
    const int n_extra   = static_cast<int>(c.extra_microstates.size());
    for (int e = 0; e < n_extra; ++e)
        c.microstates_active.push_back(n_catalog + e);

    c.active_L.assign(n_catalog + n_extra, 0);
    for (int L : c.microstates_active) c.active_L[L] = 1;

    c.w_configs.assign(catalog.n_configs(), 0.0);
    if (!c.K_weights.empty()) {
        for (size_t i = 0; i < c.K_indices.size(); ++i)
            c.w_configs[c.K_indices[i]] = c.K_weights[i];
    }

    // Default single sector: sector 0 owns every config, in K_indices order.
    c.sector_config_pos.assign(1, std::vector<int>(c.K_indices.size()));
    std::iota(c.sector_config_pos[0].begin(), c.sector_config_pos[0].end(), 0);

    // Single-sector per-sector view = the whole-pool aggregate.
    c.sector_gem_active_.assign(1, c.geminals_active_pool_);

    return c;
}

SelectorDemand selector_demand(const Cassette& cassette) {
    return cassette.selector_demand_;
}

void Cassette::install_sector_partition(std::vector<std::vector<int>> positions,
                                        const std::vector<int>&      n_gen_per_sector) {
    const int n_sec = static_cast<int>(positions.size());
    sector_config_pos = std::move(positions);

    sector_gem_active_.assign(n_sec, std::array<std::vector<int>, kMaxGen>{});
    for (int s = 0; s < n_sec; ++s) {
        std::array<std::set<int>, kMaxGen> gem_sets;
        for (int pos : sector_config_pos[s]) {
            const DiagDepsView d = catalog_.diag_deps_view(K_indices[pos]);
            for (int gen = 0; gen < kMaxGen; ++gen)
                gem_sets[gen].insert(d.geminal_reads[gen].begin(),
                                     d.geminal_reads[gen].end());
        }
        for (int gen = 0; gen < kMaxGen; ++gen)
            sector_gem_active_[s][gen] = sorted_unique_from_set(gem_sets[gen]);
    }

    n_generations_agg_ = 0;
    for (int ng : n_gen_per_sector) n_generations_agg_ = std::max(n_generations_agg_, ng);
}

void Cassette::compute_C_L(int K, const FONSnapshot& fons,
                           std::vector<double>& C_L,
                           const double* fon_cache) const {
    if (static_cast<int>(C_L.size()) < catalog_.n_microstates())
        C_L.resize(catalog_.n_microstates());
    for (int L : config_microstate_writes(K)) C_L[L] = 0.0;
    const Catalog::WeightBlock w = catalog_.weight_block(K);
    if (!w.row) return;
    EvalCtx ctx{ ::psi::reks::f_interp, fon_cache };
    eval_weight(w.bp, w.row->cell_off, w.row->cell_n, fons, ctx, C_L);
}

void Cassette::compute_weight_derivs(int gen, int K, int g,
                                     const FONSnapshot& fons,
                                     std::vector<double>& dC,
                                     std::vector<double>& d2C) const {
    for (int L : config_microstate_writes(K)) { dC[L] = 0.0; d2C[L] = 0.0; }
    const Catalog::WeightBlock w = catalog_.weight_block(K);
    if (!w.row) return;
    DerivCtx dctx{ ::psi::reks::f_interp, ::psi::reks::df_interp,
                   ::psi::reks::d2f_interp };
    eval_weight_deriv(w.bp, w.row->cell_off, w.row->cell_n, fons, dctx, gen, g, dC, d2C);
}

void Cassette::compute_mixed_weight_derivs(int gen, int K, int g_i, int g_j,
                                           const FONSnapshot& fons,
                                           std::vector<double>& d2C_mixed) const {
    for (int L : config_microstate_writes(K)) d2C_mixed[L] = 0.0;
    if (g_i == g_j) return;
    const Catalog::WeightBlock w = catalog_.weight_block(K);
    if (!w.row) return;
    DerivCtx dctx{ ::psi::reks::f_interp, ::psi::reks::df_interp,
                   ::psi::reks::d2f_interp };
    eval_weight_mixed_deriv(w.bp, w.row->cell_off, w.row->cell_n, fons, dctx, gen, g_i, g_j,
                            d2C_mixed);
}

int Cassette::compute_RDM_cells(int64_t cell_off, int cell_n, const FONSnapshot& fons,
                                int N_active, RdmCell* out, const double* fon_cache) const {
    const BlockPools bp = catalog_.block_pools(catalog_.base);
    EvalCtx ctx{ ::psi::reks::f_interp, fon_cache };
    return eval_rdm_cells(bp, cell_off, cell_n, fons, ctx, N_active, out);
}

IntSpan Cassette::config_microstate_writes(int K) const {
    return catalog_.diag_deps_view(K).microstate_writes;
}

DiagDepsView Cassette::diag_deps(int K) const {
    return catalog_.diag_deps_view(K);
}

const Microstate& Cassette::microstate(int L) const {
    const int n_catalog = catalog_.n_microstates();
    return L < n_catalog ? catalog_.microstates()[L]
                         : extra_microstates[L - n_catalog];
}

const LagrangianPair& Cassette::lagrangian_pair(int p) const {
    return catalog_.lagrangian_pairs(catalog_.base)[p];
}

std::vector<LagrangianPair> Cassette::lagrangian_pairs_view() const {
    return reks::studio::lagrangian_pairs_view(catalog_);
}

const GeminalTemplate* Cassette::geminal_templates() const {
    return catalog_.geminal_templates(catalog_.base);
}

int Cassette::n_catalog_microstates() const { return catalog_.n_microstates(); }
int Cassette::n_active_orbitals()     const { return catalog_.data->n_active_orbitals; }
int Cassette::n_geminals()            const { return catalog_.data->n_geminals; }
int Cassette::n_generations()         const {
    return n_generations_agg_ >= 0 ? n_generations_agg_ : catalog_.sector_n_generations();
}
int Cassette::n_lagrangian_pairs()    const { return catalog_.data->n_lagrangian_pairs; }
int Cassette::n_configs()             const { return catalog_.n_configs(); }
int Cassette::n_electrons()           const { return catalog_.data->n_electrons; }
int Cassette::scheme()                const { return catalog_.data->scheme; }
int Cassette::n_orbitals_per_geminal() const { return catalog_.data->n_orbitals_per_geminal; }

const FonFactor* Cassette::fon_pool() const { return catalog_.fon_pool(catalog_.base); }
int Cassette::n_fon_pool() const { return catalog_.data->n_fon_pool; }

int Cassette::n_delta_specs() const { return catalog_.n_delta_specs(); }
const DeltaSpec& Cassette::delta_spec(int i) const { return catalog_.delta_spec(i); }

const char* Cassette::si_config_name(int K) const {
    return catalog_.si_config_name(K);
}

const char* Cassette::si_config_def(int K) const {
    return catalog_.si_config_def(K);
}

CallSheet build_call_sheet(const Cassette&              sa_cassette,
                       const std::vector<Cassette>& si_cassettes) {
    const int n_active = sa_cassette.n_active_orbitals();

    CallSheet plan;
    plan.sa_active = sa_cassette.microstates_active;

    // Each cassette already folds its K2b_Fock delta F_MO carriers into
    // microstates_active, so the plain union covers every dereferenced microstate.
    std::set<int> union_set(plan.sa_active.begin(), plan.sa_active.end());
    for (const Cassette& e : si_cassettes)
        union_set.insert(e.microstates_active.begin(),
                         e.microstates_active.end());

    plan.si_sa_active = sorted_unique_from_set(union_set);

    std::set<int> bm_set;
    for (int L : plan.si_sa_active) {
        const auto& m = sa_cassette.microstate(L);
        bm_set.insert(base_density_index(m, n_active, /*beta=*/false));
        bm_set.insert(base_density_index(m, n_active, /*beta=*/true));
    }
    plan.referenced_bitmask = sorted_unique_from_set(bm_set);

    std::set<int> sa_bm_set;
    for (int L : plan.sa_active) {
        const auto& m = sa_cassette.microstate(L);
        sa_bm_set.insert(base_density_index(m, n_active, /*beta=*/false));
        sa_bm_set.insert(base_density_index(m, n_active, /*beta=*/true));
    }
    plan.sa_referenced_bitmask = sorted_unique_from_set(sa_bm_set);
    for (int p : plan.referenced_bitmask)
        if (!sa_bm_set.count(p)) plan.si_only_bitmask.push_back(p);

    std::set<std::pair<int,int>> eri_set;
    for (const Cassette& e : si_cassettes) {
        for (const std::pair<int,int>& p : collect_eri_pairs(e)) eri_set.insert(p);
    }
    plan.eri_pairs  = std::vector<std::pair<int,int>>(eri_set.begin(), eri_set.end());

    return plan;
}

}  // namespace studio
}  // namespace reks
}  // namespace psi
