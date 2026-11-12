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

#include "reks_arch.h"

#include <cstdint>
#include <vector>

namespace psi {
namespace reks {
namespace studio {

/// SA = SCF ensemble (diagonal H[K,K] only); SI = state-interaction pool (off-diagonal H[K_i,K_j]).
enum class CassetteRole { SA, SI };

struct Cassette {
    /// Index into the run's declared sector list this cassette belongs to.
    int                        run_sector = 0;

    /// Cap on states presented for this pool (print + psivars); 0 = no cap.
    int                        report_states = 0;

    std::vector<int>           K_indices;
    std::vector<double>        K_weights;     ///< SCF ensemble weights; empty for an SI pool

    /// Dependency union of microstates over K_indices (global L); extra
    /// determinants, if any, are appended at the tail.
    std::vector<int>           microstates_active;
    /// Geminal slots varied by any config of the pool, per generation [0, n_generations()).
    std::array<std::vector<int>, kMaxGen>  geminals_active_pool_;
    /// Delta-integral indices (into delta_spec) referenced by this pool's SI
    /// pairs; empty for SA pools.
    std::vector<DeltaKey>      deltas_active;

    std::vector<std::int8_t>   active_L;               ///< global-L -> is-active in this pool; size n_total_microstates
    std::vector<double>        w_configs;             ///< size n_configs

    /// SA extra determinants, global L in [n_catalog_microstates, n_total_microstates), constant weights.
    std::vector<Microstate>    extra_microstates;
    std::vector<double>        extra_weights;

    /// Positions into K_indices/K_weights of run sector s's configs (may be empty).
    std::vector<std::vector<int>>  sector_config_pos;

    /// Geminal slots varied by sector s's configs in generation gen.
    std::vector<std::array<std::vector<int>, kMaxGen>>  sector_gem_active_;

    /// max over declared sectors of n_generations; -1 = the bound sector's count.
    int  n_generations_agg_ = -1;

    /// Position of global config K in K_indices = local_of_[K - K_lo_], -1 outside the cassette.
    int                K_lo_ = 0;
    std::vector<int>   local_of_;

    /// Position of global config K in K_indices, -1 when K is not in this cassette.
    int local_index(int K) const {
        const int off = K - K_lo_;
        if (off < 0 || off >= static_cast<int>(local_of_.size())) return -1;
        return local_of_[off];
    }

    /// C_L = weight of L in config K at fons, written on config_microstate_writes(K) only (buffer grown to
    /// the microstate count); null fon_cache evaluates factors directly.
    void compute_C_L(int K, const FONSnapshot& fons,
                     std::vector<double>& C_L,
                     const double* fon_cache = nullptr) const;
    /// dC_L/dx_g and d2C_L/dx_g^2 within generation gen (x = that generation's
    /// FON). Same output window as compute_C_L; dC, d2C hold >= n_catalog_microstates() entries.
    void compute_weight_derivs(int gen, int K, int g, const FONSnapshot& fons,
                               std::vector<double>& dC,
                               std::vector<double>& d2C) const;
    /// d2C_L/dx_{g_i} dx_{g_j} within generation gen, g_i != g_j (0 when equal).
    /// Same output window as compute_C_L; d2C_mixed holds >= n_catalog_microstates() entries.
    void compute_mixed_weight_derivs(int gen, int K, int g_i, int g_j,
                                     const FONSnapshot& fons,
                                     std::vector<double>& d2C_mixed) const;

    /// Microstates (global L) of H[K,K]; diag_deps adds the geminal reads. Views into the mapped blob.
    IntSpan      config_microstate_writes(int K) const;
    DiagDepsView diag_deps(int K) const;
    /// Catalog record for L < n_catalog_microstates, runtime extra determinant beyond.
    const Microstate&       microstate(int L) const;
    const LagrangianPair&   lagrangian_pair(int p) const;
    std::vector<LagrangianPair> lagrangian_pairs_view() const;
    /// Geminal templates: one table shared across all FON generations (geometry
    /// is generation-independent). Length n_geminals().
    const GeminalTemplate*  geminal_templates() const;

    /// Catalog dimensions; n_catalog_microstates sizes global-L buffers.
    int n_catalog_microstates() const;
    int n_active_orbitals()     const;
    int n_geminals()            const;   ///< total geminal slots in the catalog, generation-independent
    int n_generations()         const;   ///< FON generations = active pairs = M/2
    int n_lagrangian_pairs()    const;
    int n_configs()             const;
    int n_electrons()           const;
    int scheme()                const;
    int n_orbitals_per_geminal() const;

    /// Off-diagonal SI coupling tables: pair records and the FON-factor pools they reference.
    /// Pairs of config K as [first, last) in CSR order, partners ascending.
    std::pair<int64_t, int64_t> pair_row(int K) const { return catalog_.pair_row(K); }
    /// Config K of the partner of pair p.
    int pair_partner(int64_t p) const { return catalog_.pair_partner(p); }
    /// Join pair p to its palette row; the counts and the damping flag live there.
    coupling::PairRef pair_ref(long long p) const { return catalog_.pair_ref(p); }
    /// Value channels of pair p, decoded from the blob into `scratch`; the spans stay
    /// valid until the next call on the same scratch.
    coupling::PairValues pair_values(long long p, std::vector<int64_t>& scratch) const {
        return catalog_.pair_values(p, scratch);
    }
    /// fon_idx/fon_slot_idx: catalog-wide row ids into fon_pool/fon_slot_pool; fon_win: absolute fon_idx
    /// offsets of the value rows, rebased to the bound sector.
    const FonFactor*              fon_pool() const;
    int                           n_fon_pool() const;
    const int*                    fon_idx() const { return catalog_.fon_idx(catalog_.base); }
    const int64_t*                fon_win() const { return catalog_.fon_win(catalog_.base); }
    const FonSlot* fon_slot_pool() const { return catalog_.fon_slot_pool(catalog_.base); }
    const int*     fon_slot_idx()  const { return catalog_.fon_slot_idx(catalog_.base); }

    /// Distinct-row pools per family (null if absent); FOCK and ERI split into idx/val parts and a row pack.
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
    /// rdm_elems rows of config K as [first, last), partners ascending.
    std::pair<int64_t, int64_t> rdm_row(int K) const { return catalog_.rdm_row(K); }
    /// Config K of the partner of rdm_elems row r.
    int rdm_partner(int64_t r) const { return catalog_.rdm_partner(r); }
    /// One rdm_elems row into positions p*N_active + q and their p != q mirrors (`out` holds 2 cell_n);
    /// returns the count written.
    int compute_RDM_cells(int64_t cell_off, int cell_n, const FONSnapshot& fons, int N_active,
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

    /// Multi-sector partition (positions and generation count per run sector) and its per-sector geminal tables.
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

    /// Active-microstate work list (global L).
    const std::vector<int>& microstates()  const { return microstates_active; }
    bool active(int L) const { return active_L[L] != 0; }

    int n_extra_microstates() const {
        return static_cast<int>(extra_microstates.size());
    }
    int n_total_microstates() const {
        return n_catalog_microstates() + n_extra_microstates();
    }

    /// Declared run-sector count.
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
/// its sector's first row.
SelectorDemand selector_demand(const Cassette& cassette);

/// Cassette over K_indices: SA unions the H[K,K] dependencies, SI adds those of H[K_i,K_j]; extras are SA-only.
Cassette make_cassette(Catalog                 catalog,
                       std::vector<int>          K_indices,
                       std::vector<double>       K_weights,
                       CassetteRole              role,
                       std::vector<Microstate>   extra_microstates = {},
                       std::vector<double>       extra_weights     = {});

/// Run-wide sets in global L: si_sa_active = SA U all SI microstates (fixed at construction),
/// referenced_bitmask = their alpha/beta bitmasks (sorted, unique).
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

/// Builds the CallSheet.
CallSheet build_call_sheet(const Cassette&              sa_cassette,
                       const std::vector<Cassette>& si_cassettes);

}  // namespace studio
}  // namespace reks
}  // namespace psi
