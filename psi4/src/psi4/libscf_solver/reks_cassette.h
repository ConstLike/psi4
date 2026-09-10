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

#include "reks_arch.h"

#include <cstdint>
#include <vector>

namespace psi {
namespace reks {
namespace studio {

/// SA = SCF ensemble (diagonal H[K,K] only); SI = state-interaction pool (off-diagonal H[K_i,K_j]).
enum class CassetteRole { SA, SI };

struct Cassette {
    /// Run sector ordinal this cassette belongs to: index into the run's
    /// declared sector list.
    int                        run_sector = 0;

    /// Cap on states presented for this pool (print + psivars); 0 = no cap.
    int                        report_states = 0;

    std::vector<int>           K_indices;
    std::vector<double>        K_weights;     ///< SCF ensemble weights; empty for an SI pool

    /// Dependency union of microstates over K_indices (global L); extra
    /// determinants, if any, are appended at the tail.
    std::vector<int>           microstates_active;
    /// Active geminal slots per FON generation, unioned over the WHOLE pool:
    /// geminals_active_pool_[gen] are the slots any config varies in generation
    /// gen (gen=0 n, 1 m, 2 u, ...). Only [0, n_generations()) are populated.
    /// The aggregate view; per-sector slots live in sector_gem_active_.
    std::array<std::vector<int>, kMaxGen>  geminals_active_pool_;
    /// Delta-integral indices (into delta_spec) referenced by this pool's SI
    /// pairs; empty for SA pools.
    std::vector<DeltaKey>      deltas_active;

    std::vector<std::int8_t>   active_L;               ///< global-L -> is-active in this pool; size n_total_microstates
    std::vector<double>        w_configs;             ///< size n_configs

    /// Runtime SA-only extra determinants: occupation records and constant
    /// weights for global L in [n_catalog_microstates, n_total_microstates).
    /// Empty for SI pools and for the default SA pool.
    std::vector<Microstate>    extra_microstates;
    std::vector<double>        extra_weights;

    /// Per-run-sector partition of the pool's configs: sector_config_pos[s] holds
    /// the positions into K_indices/K_weights whose config belongs to run sector s.
    /// Installed via install_sector_partition; a declared sector may hold no positions.
    std::vector<std::vector<int>>  sector_config_pos;

    /// Per-run-sector active geminal slots: sector_gem_active_[s][gen] are the
    /// slots sector s's configs vary in generation gen. Empty for a declared
    /// sector with no SA configs. Installed via install_sector_partition.
    std::vector<std::array<std::vector<int>, kMaxGen>>  sector_gem_active_;

    /// max over declared sectors of n_generations; -1 = read the catalog's
    /// active-sector count instead.
    int  n_generations_agg_ = -1;

    /// Global config id -> position in K_indices, over the span the ids cover:
    /// local_of_[K - K_lo_], -1 for a config outside this cassette. Built in
    /// make_cassette.
    int                K_lo_ = 0;
    std::vector<int>   local_of_;

    /// Position of global config K in K_indices, -1 when K is not in this cassette.
    int local_index(int K) const {
        const int off = K - K_lo_;
        if (off < 0 || off >= static_cast<int>(local_of_.size())) return -1;
        return local_of_[off];
    }

    /// C_L[L] = weight of microstate L in config K, as a function of fons.
    /// Defined on config_microstate_writes(K); entries outside that window keep
    /// whatever the buffer held. The buffer grows to the catalog microstate count
    /// if it is shorter. A null fon_cache evaluates each FON factor directly.
    void compute_C_L(int K, const FONSnapshot& fons,
                     std::vector<double>& C_L,
                     const double* fon_cache = nullptr) const;
    /// dC_L/dx_g and d2C_L/dx_g^2 within generation gen (x = that generation's
    /// FON). Same output window and sizing rule as compute_C_L.
    void compute_weight_derivs(int gen, int K, int g, const FONSnapshot& fons,
                               std::vector<double>& dC,
                               std::vector<double>& d2C) const;
    /// d2C_L/dx_{g_i} dx_{g_j} within generation gen, g_i != g_j (0 when equal).
    /// Same output window and sizing rule as compute_C_L.
    void compute_mixed_weight_derivs(int gen, int K, int g_i, int g_j,
                                     const FONSnapshot& fons,
                                     std::vector<double>& d2C_mixed) const;

    /// Diagonal H[K,K] dependencies of config K: the microstate-write indices
    /// (global L) and, via diag_deps, also the per-generation geminal-read slots.
    /// Both views point into the mapped blob (process lifetime).
    IntSpan      config_microstate_writes(int K) const;
    DiagDepsView diag_deps(int K) const;
    /// Catalog record for L < n_catalog_microstates, runtime extra determinant beyond.
    const Microstate&       microstate(int L) const;
    const LagrangianPair&   lagrangian_pair(int p) const;
    std::vector<LagrangianPair> lagrangian_pairs_view() const;
    /// Geminal templates: one table shared across all FON generations (geometry
    /// is generation-independent). Length n_geminals().
    const GeminalTemplate*  geminal_templates() const;

    /// Catalog dimensions. n_catalog_microstates sizes global-L bookkeeping
    /// buffers; it is a sizing dimension, not a work-loop bound.
    int n_catalog_microstates() const;
    int n_active_orbitals()     const;
    int n_geminals()            const;   ///< total geminal slots in the catalog, generation-independent
    int n_generations()         const;   ///< FON generations = active pairs = M/2
    int n_lagrangian_pairs()    const;
    int n_configs()             const;
    int n_electrons()           const;
    int scheme()                const;
    int n_orbitals_per_geminal() const;

    /// Off-diagonal SI-Hamiltonian coupling tables: pair coupling records and the
    /// FON-factor pools they reference. Mapped blob root for KIND-1 offset resolution.
    const uint8_t* base() const { return catalog_.base; }

    /// Partner key of every pair, in CSR order; pair p of row K lies in
    /// [pair_row_ptr()[K], pair_row_ptr()[K+1]), with K_j ascending inside a row.
    const int*                    pair_key_j()   const { return catalog_.pair_key_j(); }
    const int*                    pair_row_ptr() const { return catalog_.pair_row_ptr(); }
    /// Join pair p to its palette row; the counts and the damping flag live there.
    coupling::PairRef pair_ref(long long p) const { return catalog_.pair_ref(p); }
    /// Value channels of pair p, decoded from the blob into `scratch`; the spans stay
    /// valid until the next call on the same scratch.
    coupling::PairValues pair_values(long long p, std::vector<int>& scratch) const {
        return catalog_.pair_values(p, scratch);
    }
    /// fon_idx()[k] and fon_slot_idx()[k] are row ids into fon_pool()/fon_slot_pool().
    /// fon_pool, fon_slot_pool and fon_slot_idx number rows over the whole catalog;
    /// fon_idx() alone rebases to the catalog's bound sector.
    const FonFactor*              fon_pool() const;
    int                           n_fon_pool() const;
    const int*                    fon_idx() const { return catalog_.fon_idx(catalog_.base); }
    const FonSlot* fon_slot_pool() const { return catalog_.fon_slot_pool(catalog_.base); }
    const int*     fon_slot_idx()  const { return catalog_.fon_slot_idx(catalog_.base); }

    /// Global per-family DISTINCT-row pools, nullptr when the variant has no terms
    /// of that family. FOCK and ERI are column-decomposed into idx/val sub-pools
    /// plus a per-row pack.
    const coupling::FockIdxPart* fock_idxpart() const { return catalog_.fock_idxpart(catalog_.base); }
    const coupling::ValRow*      fock_valpart() const { return catalog_.fock_valpart(catalog_.base); }
    const coupling::EriIdxPart*  eri_idxpart()  const { return catalog_.eri_idxpart(catalog_.base); }
    const coupling::ValRow*      eri_valpart()  const { return catalog_.eri_valpart(catalog_.base); }
    const coupling::LagrTerm* lagr_row_pool() const { return catalog_.lagr_row_pool(catalog_.base); }
    const coupling::ValRow*      e_val_pool()   const { return catalog_.e_val_pool(catalog_.base); }
    const coupling::ValRow*      s_row_pool()   const { return catalog_.s_row_pool(catalog_.base); }
    /// Coefficient palette, numbered over the whole catalog (not per-family).
    const Coeff*                 coeff_pool()   const { return catalog_.coeff_pool(catalog_.base); }

    /// RDM block index: one row per live (key_i <= key_j) block, sorted on the
    /// key pair. Blocks absent from it are identically zero.
    const RdmRow* rdm_elems()   const { return catalog_.rdm_elems(catalog_.base); }
    const int*    rdm_row_ptr() const { return catalog_.rdm_row_ptr(catalog_.base); }
    int           n_rdm_elems() const { return catalog_.data->n_rdm_elems; }
    /// Evaluate one rdm_elems row's cells into positions p*N_active + q and the
    /// p != q mirrors; `out` must hold 2 * cell_n. Returns how many were written.
    /// A null fon_cache evaluates each FON factor directly.
    int compute_RDM_cells(int cell_off, int cell_n, const FONSnapshot& fons, int N_active,
                          RdmCell* out, const double* fon_cache = nullptr) const;

    int sector() const { return run_sector; }

    /// Whole-pool active geminals of generation gen (aggregate over sectors).
    const std::vector<int>& geminals_active_gen(int gen) const {
        return geminals_active_pool_[gen];
    }
    /// Active geminals of run sector s in generation gen (per-sector view).
    const std::vector<int>& geminals_active(int s, int gen) const {
        return sector_gem_active_[s][gen];
    }
    /// Total active FON degrees of freedom = sum over generations of the active
    /// slot count.
    int n_active_geminals() const {
        int total = 0;
        for (int gen = 0; gen < n_generations(); ++gen)
            total += static_cast<int>(geminals_active_pool_[gen].size());
        return total;
    }

    /// Install the multi-sector partition: positions[s] = positions into K_indices
    /// belonging to run sector s (empty allowed); n_gen_per_sector[s] = that
    /// sector's generation count. Rebuilds the per-sector active-geminal tables;
    /// the whole-pool aggregate is untouched.
    void install_sector_partition(std::vector<std::vector<int>> positions,
                                  const std::vector<int>&      n_gen_per_sector);

    /// Off-diagonal delta-integral specs; n_delta_specs() == 0 for variants with
    /// no off-diagonal delta integrals.
    int               n_delta_specs()      const;
    const DeltaSpec&  delta_spec(int i)    const;

    /// SI config display name for config K; nullptr if the variant has no names
    /// or K is out of range.
    const char*       si_config_name(int K) const;

    /// SI config wavefunction def string for config K; nullptr if the variant
    /// has no def strings or K is out of range.
    const char*       si_config_def(int K) const;

    /// Config index for a display name (case-insensitive); throws if absent.
    int               config_index_by_name(const std::string& name) const;

    /// The active-microstate work list (global L): catalog dependency union then
    /// the extra determinants at the tail.
    const std::vector<int>& microstates()  const { return microstates_active; }
    int                     n_microstates_active() const {
        return static_cast<int>(microstates_active.size());
    }
    bool active(int L) const { return active_L[L] != 0; }

    int n_extra_microstates() const {
        return static_cast<int>(extra_microstates.size());
    }
    int n_total_microstates() const {
        return n_catalog_microstates() + n_extra_microstates();
    }

    /// Declared run-sector count of this pool (1 at construction;
    /// install_sector_partition may install more).
    int n_sectors() const { return static_cast<int>(sector_config_pos.size()); }
    /// Positions into K_indices/K_weights whose config belongs to run sector s.
    const std::vector<int>& sector_config_positions(int s) const {
        return sector_config_pos[s];
    }

   private:
    Catalog catalog_{};
    /// Function of the blob and K_indices alone; fixed for the cassette's lifetime.
    SelectorDemand selector_demand_{};

    friend Cassette make_cassette(Catalog, std::vector<int>, std::vector<double>,
                                  CassetteRole, std::vector<Microstate>,
                                  std::vector<double>);
    friend SelectorDemand selector_demand(const Cassette&);
};

/// How far this cassette's pairs reach into each selector column, in rows counted from
/// its sector's first row. Settled when the cassette was built.
SelectorDemand selector_demand(const Cassette& cassette);

/// Builds a Cassette over K_indices from catalog. SA role unions each config's
/// diagonal H[K,K] dependencies; SI role adds every pairwise off-diagonal
/// H[K_i,K_j] dependency. extra_microstates/extra_weights append runtime
/// SA-only determinants beyond the catalog (empty for SI pools).
Cassette make_cassette(Catalog                 catalog,
                       std::vector<int>          K_indices,
                       std::vector<double>       K_weights,
                       CassetteRole              role,
                       std::vector<Microstate>   extra_microstates = {},
                       std::vector<double>       extra_weights     = {});

/// Run-wide active set. Indices are GLOBAL catalog L; no remap.
///
///   sa_active          = sa_cassette.microstates_active
///   si_sa_active       = sort_unique(sa_active U_e si_cassettes[e].microstates_active)
///   referenced_bitmask = sort_unique{ alpha_bitmask(L), beta_bitmask(L)
///                                     : L in si_sa_active }       (<= 2^M)
///
/// si_sa_active is fixed at construction.
struct CallSheet {
    std::vector<int>  sa_active;
    std::vector<int>  si_sa_active;
    std::vector<int>  referenced_bitmask;
    /// sa_referenced_bitmask: bitmask patterns over sa_active;
    /// si_only_bitmask = referenced_bitmask \ sa_referenced_bitmask.
    std::vector<int>  sa_referenced_bitmask;
    std::vector<int>  si_only_bitmask;
    /// Active-orbital ERI pairs: union over SI cassettes, canonical k<=l, sorted-unique.
    std::vector<std::pair<int,int>>  eri_pairs;
};

/// Builds the CallSheet. Pure: depends only on its arguments, does not mutate them.
CallSheet build_call_sheet(const Cassette&              sa_cassette,
                       const std::vector<Cassette>& si_cassettes);

}  // namespace studio
}  // namespace reks
}  // namespace psi
