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

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace psi {
namespace reks {
namespace studio {

struct Cassette;

/// Runtime bound on FON generations (openness orders): n, m, u, v, w, x, y, z; covers
/// REKS up to (16,16). Not an ABI constant: the blob sizes its records by its own
/// n_generations, and the loader rejects a catalog that exceeds this bound.
inline constexpr int kMaxGen = 8;

struct GeminalFON {
    double p = 1.0;   ///< occupation of orbital p of the geminal pair
    double q = 1.0;   ///< occupation of orbital q of the geminal pair
};

struct FONSnapshot {
    /// Generation-indexed FON storage: layers[gen][g] is the FON of geminal slot
    /// g in generation gen (gen = openness order: 0=n PPS, 1=m OSS, 2=u DOSS, ...).
    std::array<std::vector<GeminalFON>, kMaxGen> layers;
};

/// One contiguous FON block of the joint optimizer vector, ordered lexicographically
/// by (s, gen).
struct FonBlock {
    int s;        ///< run sector ordinal
    int gen;      ///< FON generation
    int offset;   ///< block's start index in the joint vector
    int count;    ///< number of FON DOF in this block (0 if none active)
};

/// Index into the delta_specs pool (see DeltaSpec).
using DeltaKey = int;

/// Decoded diagonal-block dependency record: (byte_offset, count) windows into a
/// shared int pool. The blob carries one gem slot per generation it declares; the
/// loader widens it to kMaxGen, leaving the unused tail empty. Resolved to a
/// DiagDepsView via the blob base (see below).
struct DiagDeps {
    uint64_t mw_off, mw_n;                        ///< microstate_writes window
    struct { uint64_t off, n; } gem[kMaxGen];     ///< [gen] -> geminal_reads window
};

/// Non-owning view of a run of ints.
struct IntSpan {
    const int* data = nullptr;
    int        n    = 0;
    constexpr const int* begin() const noexcept { return data; }
    constexpr const int* end()   const noexcept { return data + n; }
    constexpr bool       empty() const noexcept { return n == 0; }
    constexpr int        size()  const noexcept { return n; }
    constexpr int        operator[](int i) const noexcept { return data[i]; }
};

/// Off-diagonal dependency record: lengths of a pair's three dependency
/// windows, read from the pair's PairShape row. Widened to kMaxGen like DiagDeps.
struct PairDeps {
    uint32_t fock_n;              ///< fock_footprint window length
    uint32_t gem_n[kMaxGen];      ///< [gen] -> geminal_reads window length
    uint32_t delta_n;             ///< delta_reads window length
};

/// Base-resolved view of a DiagDeps record; IntSpans point into the mapped blob.
struct DiagDepsView {
    IntSpan microstate_writes;
    IntSpan geminal_reads[kMaxGen];
};

/// Decoded view of a PairDeps record; IntSpans point into the scratch buffer
/// decode_pair_deps was given, valid until its next call on that buffer.
struct PairDepsView {
    IntSpan fock_footprint;            ///< microstate indices L whose F_MO[L] is read
    IntSpan geminal_reads[kMaxGen];    ///< [gen] -> slots read; index == FON generation
    IntSpan delta_reads;
};

inline DiagDepsView resolve_diag_deps(const uint8_t* base, const DiagDeps& d) {
    auto span = [base](uint64_t off, uint64_t n) {
        return IntSpan{ n ? reinterpret_cast<const int*>(base + off) : nullptr,
                        static_cast<int>(n) };
    };
    DiagDepsView v;
    v.microstate_writes = span(d.mw_off, d.mw_n);
    for (int g = 0; g < kMaxGen; ++g) v.geminal_reads[g] = span(d.gem[g].off, d.gem[g].n);
    return v;
}

/// Runtime bound on active orbitals; covers REKS up to (12,12). Not an ABI constant:
/// the blob stores n_active_orbitals occupations per string, and the loader rejects
/// a catalog that exceeds this bound.
inline constexpr int kMaxActiveOrb = 12;

/// One decoded Slater microstate over the active space: alpha[i]/beta[i] in {0,1}
/// is the occupation of active orbital i (i < n_active); the tail up to
/// kMaxActiveOrb is zero.
struct Microstate {
    std::array<int8_t, kMaxActiveOrb> alpha{};
    std::array<int8_t, kMaxActiveOrb> beta{};
};

/// Occupation bitmask of a microstate's alpha (beta=false) or beta string:
/// sum_i occ[i] << i over the n_active active orbitals.
inline int base_density_index(const Microstate& m, int n_active, bool beta) {
    int idx = 0;
    for (int i = 0; i < n_active; ++i)
        idx += static_cast<int>(beta ? m.beta[i] : m.alpha[i]) << i;
    return idx;
}

/// One geminal slot: scheme tags its coupling group (0 = canonical SCF pairs);
/// orbitals is the (p, q) active-orbital index pair.
struct GeminalTemplate {
    int                scheme = 0;
    std::array<int, 2> orbitals{};
};

struct LagrangianPair {
    int orb_from = 0;
    int orb_to   = 0;
};

// Off-diagonal integral forms. Equation/appendix numbers refer to
// Filatov, M. et al. J. Chem. Phys. 2017, 147, 064104.
enum class DeltaKind : int {
    K2b_Fock = 0,   ///< Fock-difference form, Eqs. B2 + B38-B49
    K2b_ERI  = 1,   ///< ERI fallback when no clean (L1, L2) microstate pair exists
    K2c_ERI  = 2    ///< 4-distinct-orbital ERI
};

/// One off-diagonal delta; kind (DeltaKind) selects between the Fock-difference
/// fields (fock_orbs, L1, L2, sign) and the ERI fields (eri_idx).
struct DeltaSpec {
    DeltaKind   kind;
    int         fock_orbs[2];  ///< (p_active, q_active) for K2b_Fock; else {-1, -1}
    int         L1;            ///< clean Ms=0 microstate (K2b_Fock); else -1
    int         L2;            ///< clean Ms=+1 microstate (K2b_Fock); else -1
    int         sign;          ///< +1 if common in alpha of L1; -1 if in beta; 0 for ERI kinds
    int         eri_idx[4];    ///< (i1, i2, i3, i4) active-orbital indices for ERI kinds; else {-1, -1, -1, -1}
};

/// One spin manifold of the merged (N,M) catalog. Configs of sector b occupy the
/// global block [config_start, config_start + config_count). sa_default holds GLOBAL
/// config indices. Sectors are ordered ascending spin2; sector 0 is the lowest 2S.
///
/// Pools that grow with the config pairs are numbered per sector, from zero at
/// each sector's own first row; add the matching base_* for an absolute index.
/// Origins ascend with the sector index, and sector 0 starts every pool at 0.
struct SectorEntry {
    int      spin2;          ///< total spin as 2S for this manifold
    int      config_start;   ///< global config index of the first config in this sector
    int      config_count;   ///< number of configs in this sector
    int      n_generations;  ///< FON generations of this sector's configs
    int      n_sa_default;
    uint64_t off_sa_default;      ///< default SA pool (GLOBAL config indices)
    uint64_t base_e_val;          ///< row origin in e_val_pool
    uint64_t base_fock_idxpart;
    uint64_t base_fock_valpart;
    uint64_t base_fock_pack;      ///< ROW origin in fock_pack (2 ints per row)
    uint64_t base_eri_idxpart;
    uint64_t base_eri_valpart;
    uint64_t base_eri_pack;       ///< ROW origin in eri_pack (2 ints per row)
    uint64_t base_lagr_row;
    uint64_t base_s_row;
    uint64_t base_fon_idx;        ///< origin of the FON windows the value rows carry
};

// Value-formula IR: each element (H off-diagonal, 1-RDM cell, overlap, config weight)
// is a signed sum of terms, each coeff * product(FON factors) * optional payload.

enum FonKind : int {
    kFonLinear    = 0,   // * (use_q ? q : p)
    kFonSqrtGroup = 1,   // * sqrt( prod over radicand window of (use_q ? q : p) )
    kFonFinterp   = 2    // * f_interp(p * q)
};

struct FonSlot {         // one geminal read: fons.layers[gen][g].p (or .q)
    int gen;             // FON generation (n=0, m=1, u=2, v=3, w=4, x=5, y=6, z=7)
    int g;               // global geminal slot
    int use_q;           // read .q instead of .p (LINEAR / SQRT radicand only)
};

struct FonFactor {       // one multiplicative factor of a term
    int kind;            // FonKind
    int slot;            // kFonLinear / kFonFinterp: index into the FonSlot pool
    int rad_off;         // kFonSqrtGroup: window [rad_off, rad_off+rad_n) into FonSlot pool
    int rad_n;           // kFonSqrtGroup: radicand length (>=1); 0 for non-group kinds
};

struct Term {            // coeff * prod(fon[fon_off .. fon_off+fon_n])
    double coeff;
    int    fon_off;      // window into the FonFactor pool
    int    fon_n;
};

// One entry of the catalog's coefficient palette.
struct Coeff { double coeff; };

// Block / vector output: RDM (cells over (p,q)) AND weights (cells over (config, L)).
struct Cell   { int p, q; int term_off, term_n; };   // RDM: out[p*N+q]; weight: p=L, q=-1

// Interned (p, q) part of a Cell.
struct CellPq  { int p, q; };

// One row per config, in config order (row index == config K).
struct WeightRow { int cell_off, cell_n; };

// A 1-RDM row keyed (K_i, K_j); K_i is the CSR row (see Catalog::rdm_row_ptr).
struct RdmRow { int key_j; int cell_off, cell_n; };

// One evaluated RDM cell at flat position pq = p*N + q; operator< sorts a
// block into ascending (p, q) order.
struct RdmCell {
    int    pq;
    double val;
    bool operator<(const RdmCell& o) const { return pq < o.pq; }
};

// Pools shared by every block row; fon_idx is the DIAG_FON_IDX block row selector
// (distinct from the coupling FON_IDX). cell_win and term_win are CSR boundary
// arrays: window w spans [win[w], win[w+1]).
struct BlockPools {
    const CellPq*    cell_pq;
    const int*       cell_win;
    const int*       cell_pack;
    const Coeff*     coeff_pool;
    const int*       term_win;
    const int*       term_pack;
    const FonFactor* fon_pool;
    const FonSlot*   fon_slot_pool;
    const int*       fon_idx;
    const int*       fon_slot_idx;
};

// Reassemble Cell c from the CELL_POOL sub-pools: (p, q) from cell_pq and the
// Term window from cell_win, selected by the (pq_id, win_id) pack entry.
inline Cell load_cell(const BlockPools& bp, int c) {
    const CellPq& pq = bp.cell_pq[bp.cell_pack[2 * c]];
    const int     w  = bp.cell_pack[2 * c + 1];
    return Cell{pq.p, pq.q, bp.cell_win[w], bp.cell_win[w + 1] - bp.cell_win[w]};
}

// Reassemble Term t: the coefficient from the palette and the FON window from
// term_win, selected by the (coeff_id, win_id) pack entry.
inline Term load_term(const BlockPools& bp, int t) {
    const int w = bp.term_pack[2 * t + 1];
    return Term{bp.coeff_pool[bp.term_pack[2 * t]].coeff,
                bp.term_win[w], bp.term_win[w + 1] - bp.term_win[w]};
}

// Off-diagonal H_ij = <K_i|H|K_j> coupling IR: 5-channel decomposition
// (energy, Fock, Lagrangian, ERI, overlap-FON).
namespace coupling {

// The value part every coupling channel shares: `word` packs the row's coefficient
// index in the palette, the length of its FON window and the window's start.
struct ValRow { int64_t word; };

// coeff * prod(fons) * (F^alpha_L - F^beta_L)[p, q]
struct FockTerm{ int L; int p; int q;
                 double coeff; int64_t fon_off; int fon_n; };
// coeff * prod(fons) * W[w];  w = lagrangian_pairs index. The pair index puts the
// packed word at offset 4, so the record carries it unaligned in 12 bytes.
#pragma pack(push, 4)
struct LagrTerm{ int w; int64_t word; };
#pragma pack(pop)
// coeff * prod(fons) * (p,q|r,s)   spatial Mulliken active-orbital ERI.
struct EriTerm { int p; int q; int r; int s;
                 double coeff; int64_t fon_off; int fon_n; };

// Index part of a FockTerm's (L, p, q); joined to its ValRow via a
// (idxpart_id, valpart_id) pack entry.
struct FockIdxPart { int L; int p; int q; };

// Index part of an EriTerm's (p, q, r, s); joined to its ValRow via a
// (idxpart_id, valpart_id) pack entry.
struct EriIdxPart { int p; int q; int r; int s; };

/// Count tuple a pair's windows follow; pairs with identical counts share one
/// row via pair_shape_id.
struct PairShape {
    int has_scaling;          // element carries an empirical damping header
    int n_e;                  // n_e (L, val_id) pairs; 2*n_e ints in the window
    int n_fock, n_lagr, n_eri, n_s;
    PairDeps deps;            // SI active-set footprint: fock/geminal/delta reads
};

/// A pair joined to its palette row. K_i is not stored: a pair is always reached
/// through pair_row_ptr, whose row names it.
struct PairRef {
    int              K_j;     ///< partner config, energy-catalog K-order, K_i < K_j
    const PairShape* sh;      ///< palette row carrying the counts and the damping flag
};

/// The five value channels of one pair, decoded out of its varint block.
struct PairValues {
    IntSpan e_pack;           ///< 2*n_e ints: (L, val_id) pairs into e_val_pool
    IntSpan fock_idx;         ///< selector rows into the fock pack
    IntSpan lagr_idx;         ///< rows of lagr_row_pool
    IntSpan eri_idx;          ///< selector rows into the eri pack
    IntSpan s_idx;            ///< rows of s_row_pool
};

}  // namespace coupling

/// Advance `off` past a pair's five value channels to its dependency windows,
/// raising `fock_rows`/`eri_rows` to the rows this pair reaches in each selector
/// column. The delta runs are signed and not monotone, so the reach is the running
/// maximum over a window, not its last value; no channel is materialised.
void scan_pair_values(const uint8_t* base, uint64_t& off,
                      const coupling::PairShape& sh,
                      long long& fock_rows, long long& eri_rows);

/// Decode one pair's five value channels straight from the blob, `off` starting at its
/// block and ending at its dependency windows. No resident pool carries them.
/// `scratch` is reused across calls; the returned spans point into it and stay
/// valid until the next call on the same scratch.
coupling::PairValues decode_pair_values(const uint8_t* base, uint64_t& off,
                                        const coupling::PairShape& sh,
                                        std::vector<int>& scratch);

/// Decode one pair's three dependency windows, `off` starting at them (past the value
/// channels) and ending at the block end. Same scratch contract as decode_pair_values.
PairDepsView decode_pair_deps(const uint8_t* base, uint64_t& off,
                              const coupling::PairShape& sh, std::vector<int>& scratch);

// Offset-ABI catalog header (POD, trivially-copyable): counts are int64, and
// byte offsets are ABSOLUTE from the blob base (KIND-1) with their byte length
// beside them as uint64. Resolve offsets via the Catalog accessors below.
// Counts reuse the scalar dimensions where one exists; pools without one carry
// their own n_* field.
struct CatalogData {
    int64_t n_microstates;
    int64_t n_configs;
    int64_t n_geminals;              ///< slot g names a (scheme, pair), generation-independent
    int64_t n_generations;           ///< FON generations, max over sectors (per-sector count in SectorEntry)
    int64_t n_lagrangian_pairs;
    int64_t n_active_orbitals;
    int64_t n_electrons;             ///< active electrons; N==M variants have n_electrons == n_active_orbitals
    int64_t scheme;                  ///< variant scheme selector (0 = default; multi-scheme variants use higher values)
    int64_t n_orbitals_per_geminal;

    uint64_t string_table_off;   ///< byte offset of the string table (KIND-1); labels relative to it
    uint64_t has_si_config_defs; ///< 1 if si_config_defs present, else 0

    uint64_t off_microstates;        // count n_microstates
    uint64_t off_geminal_templates;  // count n_geminals
    uint64_t off_lagrangian_pairs;   // count n_lagrangian_pairs
    uint64_t off_diag_deps;          // count n_configs (DiagDeps records)
    uint64_t off_si_config_names;    // count n_configs (u32[] of relative string offsets)
    uint64_t off_si_config_defs;     // count n_configs (u32[]); valid iff has_si_config_defs

    int64_t n_delta_specs;   uint64_t off_delta_specs;

    // Pair section: column arrays over the n_pairs emitted config pairs, a palette
    // of their count tuples (PairShape), and one varint block per pair in
    // pair_win_data holding all of that pair's windows back to back.
    // n_pair_rows == n_configs + 1.
    int64_t n_pairs;
    int64_t n_pair_rows;     uint64_t off_pair_row_ptr;   // int32[n_pair_rows], CSR over K_i
    uint64_t off_pair_key_j;                              // int32[n_pairs]
    uint64_t off_pair_shape;                              // uint32[n_pairs], index into pair_shapes
    uint64_t off_pair_win_off;                            // uint64[n_pairs], absolute byte offsets
    int64_t n_pair_shapes;   uint64_t off_pair_shapes;    // (8 + n_generations) int32 per shape
    uint64_t off_pair_win_data;  uint64_t pair_win_bytes;

    int64_t n_fon_pool;      uint64_t off_fon_pool;     // FON_VALUE_POOL (coupling + block fon_idx select)
    int64_t n_fon_idx;       uint64_t off_fon_idx;      // FON_IDX coupling row selector
    int64_t n_fock_idxpart;  uint64_t off_fock_idxpart;
    int64_t n_fock_valpart;  uint64_t off_fock_valpart;
    // Selector rows (idxpart id, valpart id), stored as two delta+zigzag varint
    // columns of n_*_pack_rows values spanning *_pack_bytes from off_*_pack.
    int64_t n_fock_pack_rows; uint64_t off_fock_pack;  int64_t fock_pack_bytes;
    int64_t n_eri_idxpart;   uint64_t off_eri_idxpart;
    int64_t n_eri_valpart;   uint64_t off_eri_valpart;
    int64_t n_eri_pack_rows; uint64_t off_eri_pack;    int64_t eri_pack_bytes;
    int64_t n_lagr_pool;     uint64_t off_lagr_row_pool;
    int64_t n_e_val_pool;    uint64_t off_e_val_pool;
    int64_t n_s_row_pool;    uint64_t off_s_row_pool;
    int64_t n_weight_elems;  uint64_t off_weight_elems;   // one WeightRow per config, in config order
    int64_t n_rdm_elems;
    int64_t n_rdm_rows;      uint64_t off_rdm_row_ptr;    // int32[n_rdm_rows], CSR over K_i
    uint64_t off_rdm_elems;
    int64_t n_fon_slot_pool; uint64_t off_fon_slot_pool;
    int64_t n_fon_slot_idx;  uint64_t off_fon_slot_idx;

    // Block-evaluator (DIAG) pools. n_cell_win and n_term_win count WINDOWS;
    // each boundary array holds one more entry than that.
    int64_t n_diag_fon_idx;  uint64_t off_diag_fon_idx;  // DIAG_FON_IDX block row selector
    int64_t n_cell_pq;       uint64_t off_cell_pq;
    int64_t n_cell_win;      uint64_t off_cell_win;
    int64_t n_cell_pack;     uint64_t off_cell_pack;
    int64_t n_coeff_pool;    uint64_t off_coeff_pool;    // coefficient palette, whole-catalog
    int64_t n_term_win;      uint64_t off_term_win;
    int64_t n_term_pack;     uint64_t off_term_pack;

    int64_t n_sectors;       uint64_t off_sectors;  // SectorEntry table; one per spin manifold
};

struct Catalog {
    const uint8_t*     base = nullptr;   ///< mapped blob root (KIND-1 offset resolution)
    const CatalogData* data = nullptr;   ///< catalog header; owned copy, not a blob alias
    int                active_sector = -1;  ///< run's selected blob sector; -1 = unbound (whole catalog)

    // The three pools the blob sizes by its own dimensions, widened at load into
    // fixed-layout arrays owned alongside the blob buffer (process lifetime).
    // Everything else stays in the mapping and is read where it is used.
    const Microstate*             microstates_ = nullptr;
    const DiagDeps*               diag_deps_   = nullptr;
    const coupling::PairShape*    pair_shapes_ = nullptr;

    template <class T>
    const T* pool_(const uint8_t* b, uint64_t off, int64_t count) const {
        return count > 0 ? reinterpret_cast<const T*>(b + off) : nullptr;
    }

    /// Sector table entry of the bound sector; the ten base_* origins live here.
    const SectorEntry& bound_sector(const uint8_t* b) const { return sectors(b)[active_sector]; }

    /// Pool base advanced to the bound sector's first row: a zero-based id from
    /// that sector indexes the result directly. `stride` is the elements one
    /// row occupies (2 for the interleaved selector packs).
    template <class T>
    const T* sector_pool_(const T* pool, uint64_t base_row, int64_t stride = 1) const {
        return pool != nullptr ? pool + base_row * stride : nullptr;
    }

    const Microstate*             microstates() const { return microstates_; }

    /// Partner key of every pair, in CSR order; the row it lies in names K_i.
    const int* pair_key_j() const {
        if (base == nullptr || data == nullptr || data->n_pairs <= 0) return nullptr;
        return reinterpret_cast<const int*>(base + data->off_pair_key_j);
    }
    /// Palette id of every pair, in CSR order.
    const uint32_t* pair_shape_id() const {
        if (base == nullptr || data == nullptr || data->n_pairs <= 0) return nullptr;
        return reinterpret_cast<const uint32_t*>(base + data->off_pair_shape);
    }
    /// Absolute byte offset of each pair's varint block in pair_win_data.
    const uint64_t* pair_win_off() const {
        if (base == nullptr || data == nullptr || data->n_pairs <= 0) return nullptr;
        return reinterpret_cast<const uint64_t*>(base + data->off_pair_win_off);
    }

    /// Join pair p to its palette row.
    coupling::PairRef pair_ref(long long p) const {
        return coupling::PairRef{ pair_key_j()[p], pair_shapes_ + pair_shape_id()[p] };
    }

    /// Value channels of pair p, decoded on demand (see decode_pair_values).
    coupling::PairValues pair_values(long long p, std::vector<int>& scratch) const {
        uint64_t off = pair_win_off()[p];
        return decode_pair_values(base, off, pair_shapes_[pair_shape_id()[p]], scratch);
    }

    /// Dependency windows of pair p, decoded on demand (see decode_pair_deps).
    /// Walking past the pair's value channels (see scan_pair_values) also
    /// raises `fock_rows`/`eri_rows` to this pair's selector reach.
    PairDepsView pair_deps(long long p, std::vector<int>& scratch,
                           long long& fock_rows, long long& eri_rows) const {
        const coupling::PairShape& sh = pair_shapes_[pair_shape_id()[p]];
        uint64_t off = pair_win_off()[p];
        scan_pair_values(base, off, sh, fock_rows, eri_rows);
        return decode_pair_deps(base, off, sh, scratch);
    }

    /// CSR row pointer over K_i into the pair columns: pair p of row K lies in
    /// [pair_row_ptr()[K], pair_row_ptr()[K+1]), with K_j ascending inside a row.
    /// Null when the blob carries no row pointer.
    const int* pair_row_ptr() const {
        if (base == nullptr || data == nullptr || data->n_pair_rows <= 0) return nullptr;
        return reinterpret_cast<const int*>(base + data->off_pair_row_ptr);
    }

    const GeminalTemplate* geminal_templates(const uint8_t* b) const { return pool_<GeminalTemplate>(b, data->off_geminal_templates, data->n_geminals); }
    const LagrangianPair*  lagrangian_pairs(const uint8_t* b)  const { return pool_<LagrangianPair>(b, data->off_lagrangian_pairs, data->n_lagrangian_pairs); }
    const DeltaSpec*       delta_specs(const uint8_t* b)       const { return pool_<DeltaSpec>(b, data->off_delta_specs, data->n_delta_specs); }
    const FonFactor*       fon_pool(const uint8_t* b)          const { return pool_<FonFactor>(b, data->off_fon_pool, data->n_fon_pool); }

    // Pools numbered per sector: the pair windows' channel ids, the selector
    // rows' two column ids, and the value rows' FON window offsets all count
    // from zero at the emitting sector; each accessor below rebases to the
    // bound sector's first row.
    const int*             fon_idx(const uint8_t* b)           const { return sector_pool_(pool_<int>(b, data->off_fon_idx, data->n_fon_idx), bound_sector(b).base_fon_idx); }
    const coupling::FockIdxPart* fock_idxpart(const uint8_t* b) const { return sector_pool_(pool_<coupling::FockIdxPart>(b, data->off_fock_idxpart, data->n_fock_idxpart), bound_sector(b).base_fock_idxpart); }
    const coupling::ValRow* fock_valpart(const uint8_t* b)     const { return sector_pool_(pool_<coupling::ValRow>(b, data->off_fock_valpart, data->n_fock_valpart), bound_sector(b).base_fock_valpart); }
    const coupling::EriIdxPart* eri_idxpart(const uint8_t* b)  const { return sector_pool_(pool_<coupling::EriIdxPart>(b, data->off_eri_idxpart, data->n_eri_idxpart), bound_sector(b).base_eri_idxpart); }
    const coupling::ValRow* eri_valpart(const uint8_t* b)      const { return sector_pool_(pool_<coupling::ValRow>(b, data->off_eri_valpart, data->n_eri_valpart), bound_sector(b).base_eri_valpart); }
    const coupling::LagrTerm* lagr_row_pool(const uint8_t* b)  const { return sector_pool_(pool_<coupling::LagrTerm>(b, data->off_lagr_row_pool, data->n_lagr_pool), bound_sector(b).base_lagr_row); }
    const coupling::ValRow* e_val_pool(const uint8_t* b)       const { return sector_pool_(pool_<coupling::ValRow>(b, data->off_e_val_pool, data->n_e_val_pool), bound_sector(b).base_e_val); }
    const coupling::ValRow* s_row_pool(const uint8_t* b)       const { return sector_pool_(pool_<coupling::ValRow>(b, data->off_s_row_pool, data->n_s_row_pool), bound_sector(b).base_s_row); }

    /// Coefficient palette, numbered over the whole catalog (no sector origin).
    const Coeff*           coeff_pool(const uint8_t* b)        const { return pool_<Coeff>(b, data->off_coeff_pool, data->n_coeff_pool); }

    const WeightRow*       weight_elems(const uint8_t* b)      const { return pool_<WeightRow>(b, data->off_weight_elems, data->n_weight_elems); }
    const RdmRow*          rdm_elems(const uint8_t* b)         const { return pool_<RdmRow>(b, data->off_rdm_elems, data->n_rdm_elems); }
    /// CSR row pointer over K_i into rdm_elems: the rows of config K are
    /// [rdm_row_ptr()[K], rdm_row_ptr()[K+1]), with K_j ascending inside a row.
    const int*             rdm_row_ptr(const uint8_t* b)       const { return pool_<int>(b, data->off_rdm_row_ptr, data->n_rdm_rows); }
    const FonSlot*         fon_slot_pool(const uint8_t* b)     const { return pool_<FonSlot>(b, data->off_fon_slot_pool, data->n_fon_slot_pool); }
    const int*             fon_slot_idx(const uint8_t* b)      const { return pool_<int>(b, data->off_fon_slot_idx, data->n_fon_slot_idx); }
    const int*             diag_fon_idx(const uint8_t* b)      const { return pool_<int>(b, data->off_diag_fon_idx, data->n_diag_fon_idx); }
    const CellPq*          cell_pq(const uint8_t* b)           const { return pool_<CellPq>(b, data->off_cell_pq, data->n_cell_pq); }
    const int*             cell_win(const uint8_t* b)          const { return pool_<int>(b, data->off_cell_win, data->n_cell_win); }
    const int*             cell_pack(const uint8_t* b)         const { return pool_<int>(b, data->off_cell_pack, data->n_cell_pack); }
    const int*             term_win(const uint8_t* b)          const { return pool_<int>(b, data->off_term_win, data->n_term_win); }
    const int*             term_pack(const uint8_t* b)         const { return pool_<int>(b, data->off_term_pack, data->n_term_pack); }
    const SectorEntry*     sectors(const uint8_t* b)           const { return pool_<SectorEntry>(b, data->off_sectors, data->n_sectors); }

    int n_sectors_available() const { return static_cast<int>(data->n_sectors); }
    // b is a blob sector index in [0, n_sectors_available()).
    const SectorEntry& sector_entry(int b) const { return sectors(base)[b]; }
    /// GLOBAL config index K -> its blob sector index; -1 if K is out of range.
    int sector_of_config(int K) const {
        const SectorEntry* se = sectors(base);
        for (int b = 0; b < n_sectors_available(); ++b)
            if (K >= se[b].config_start && K < se[b].config_start + se[b].config_count) return b;
        return -1;
    }
    /// Default SA pool of blob sector sec (GLOBAL config indices).
    const int* sector_sa_default(const uint8_t* b, int sec) const {
        const SectorEntry& se = sectors(b)[sec];
        return se.n_sa_default > 0 ? reinterpret_cast<const int*>(b + se.off_sa_default) : nullptr;
    }

    // Bound-sector convenience: read the run's selected sector (active_sector).
    int sector_config_start()  const { return sector_entry(active_sector).config_start; }
    int sector_config_count()  const { return sector_entry(active_sector).config_count; }
    int sector_n_generations() const { return sector_entry(active_sector).n_generations; }
    int sector_n_sa_default()  const { return sector_entry(active_sector).n_sa_default; }

    BlockPools block_pools(const uint8_t* b) const {
        return BlockPools{ cell_pq(b), cell_win(b), cell_pack(b), coeff_pool(b),
                           term_win(b), term_pack(b), fon_pool(b), fon_slot_pool(b),
                           diag_fon_idx(b), fon_slot_idx(b) };
    }
    DiagDepsView diag_deps_view(const uint8_t* b, int K) const {
        return resolve_diag_deps(b, diag_deps_[K]);
    }

    // String accessors (RELATIVE scheme): rel is a byte offset RELATIVE to string_table_off.
    const char* string_at(uint32_t rel) const {
        return reinterpret_cast<const char*>(base + data->string_table_off + rel);
    }
    const char* si_config_name(int K) const {
        if (!base || K < 0 || K >= data->n_configs) return nullptr;
        const uint32_t* idx = reinterpret_cast<const uint32_t*>(base + data->off_si_config_names);
        return string_at(idx[K]);
    }
    const char* si_config_def(int K) const {
        if (!base || !data->has_si_config_defs || K < 0 || K >= data->n_configs) return nullptr;
        const uint32_t* idx = reinterpret_cast<const uint32_t*>(base + data->off_si_config_defs);
        return string_at(idx[K]);
    }
};

/// How far into the fock/eri selector columns, from the bound sector's first
/// row, a caller's pairs need decoded.
struct SelectorDemand {
    long long fock_rows = 0;
    long long eri_rows  = 0;
};

/// Decoded fock/eri selector rows -- (idxpart id, valpart id) per row -- up to
/// the requested demand. Row 0 is the bound sector's first row, matching the
/// pair windows' numbering.
struct SelectorPacks {
    std::vector<int> fock;   ///< 2 ints per row
    std::vector<int> eri;
};

/// Decode each selector column of cat's bound sector as far as `demand`
/// reaches, with every stored id range-checked against that sector's
/// sub-pools; only the demanded rows are returned.
SelectorPacks decode_selector_packs(const Catalog& cat, const SelectorDemand& demand);

/// Loads (or returns the cached) merged reks{N}{M}.bin catalog and binds it to the
/// spin manifold 2S = spin. Throws if the blob is absent, fails validation, or has
/// no matching sector.
Catalog make_catalog(int N, int M, int spin);

/// Loads (or returns the cached) merged reks{N}{M}.bin catalog bound to blob
/// sector 0, the lowest available 2S by sector-table order. Throws if the blob
/// is absent or fails validation.
Catalog make_catalog(int N, int M);

/// Display letter for a FON generation index: 0->'n', 1->'m', 2->'u', 3->'v',
/// 4->'w', 5->'x', 6->'y', 7->'z'. Throws if gen is outside [0, kMaxGen).
char generation_letter(int gen);

/// Display name of a FON channel: generation_letter(gen), then the scheme
/// digit if scheme >= 0, then "t<run_sector>" if run_sector > 0.
/// E.g. gen=1 -> "m"; gen=1,run_sector=1 -> "mt1";
/// gen=2,run_sector=2,scheme=1 -> "u1t2".
std::string fon_channel_name(int gen, int run_sector, int scheme = -1);

/// Wavefunction array-variable base name for a FON generation: gen 0 -> "REKS FON",
/// gen >= 1 -> "REKS <LETTER>_FON" (M_FON, U_FON, V_FON, W_FON).
std::string fon_array_name(int gen);

/// True iff a merged reks{N}{M}.bin is installed (filename probe; the spin
/// manifold is resolved and checked in make_catalog).
bool variant_supported(int N, int M);

/// Installed active spaces formatted as "[N,M], ...".
std::string supported_variants_str();

/// Default SI si_cassettes: single si_cassette holding the bound sector's config block,
/// [config_start, config_start + config_count) in GLOBAL config indices.
std::vector<std::vector<int>> default_si_indices(const Catalog& s);

/// si_config_names display name -> global config index K (case-insensitive),
/// or -1 if the name is not in the catalog.
int find_config_index(const Catalog& s, const std::string& name);

/// si_config_names joined as "NAME1, NAME2, ..." in catalog order.
std::string config_names_csv(const Catalog& s);

/// Config type of a display name: the name minus its trailing running index
/// (_?[0-9]+). "PPS1"->"PPS", "DOSS_SPS_0"->"DOSS_SPS".
std::string config_type_of(const std::string& name);

/// Distinct config types of the bound sector's block, joined "TYPE1, TYPE2, ..."
/// in first-occurrence catalog order.
std::string config_types_csv(const Catalog& s);

/// Bulk config-selection token -> the bound sector's matching GLOBAL config
/// indices, in catalog order. Tokens (case-insensitive): "full" selects the
/// whole sector block; "sps" the spin-projected base; a pattern holding '*' or
/// '?' the names it matches; a config-type name selects every config of that type;
/// "single-scheme[:N]" selects the configs of pairing scheme N; "mixed-scheme[:N]"
/// selects the SSR state set on scheme N (its single-excitation configs plus every
/// open-shell-singlet config). The spin-projected configs belong to no single
/// scheme and are the manifold's base: every scheme token carries them whole and
/// none repeats them. A scheme token
/// takes one index or several ("mixed-scheme-1-2-3" = their union), each behind
/// ':', '-', '_' or ',', and defaults to scheme 0; it throws when an index is
/// malformed, repeated, out of range, or the selection is empty. Empty when the
/// token is neither "full" nor a type present in the sector.
std::vector<int> expand_config_token(const Catalog& s, const std::string& token);

/// find_config_index, throwing PsiException with the valid-name list on miss.
int config_index_by_name(const Catalog& s, const std::string& name);

/// Uniform weights vector of length n; each entry = 1.0 / n. Returns empty
/// vector when n <= 0.
std::vector<double> uniform_weights(int n);

/// True iff the catalog provides a diabatic RDM (n_rdm_elems > 0).
bool has_diabatic_rdm(const Catalog& cat);

/// Active-orbital ERI pairs (eri_idx[2], eri_idx[3]) of the cassette's K2b_ERI and
/// K2c_ERI deltas, canonical (k<=l) and sorted lexicographically. Empty when the
/// cassette references neither.
std::vector<std::pair<int,int>> collect_eri_pairs(const Cassette&   cassette);

/// Owning-vector copy of the catalog's lagrangian_pairs (length n_lagrangian_pairs).
std::vector<LagrangianPair> lagrangian_pairs_view(const Catalog& cat);

// Compile-time ABI lock: fixed-stride struct geometry + semantic constants must equal the
// reks_layout.json single source (drives the runtime layout_hash; a drift breaks the build).
// Microstate, DiagDeps and PairShape are absent: the blob sizes them by its own dimensions
// and the loader decodes them, so their C++ geometry is not ABI.
static_assert(sizeof(FonFactor) == 16 && alignof(FonFactor) == 4, "FonFactor");
static_assert(sizeof(FonSlot) == 12 && alignof(FonSlot) == 4, "FonSlot");
static_assert(sizeof(Coeff) == 8 && alignof(Coeff) == 8, "Coeff");
static_assert(sizeof(CellPq) == 8 && alignof(CellPq) == 4, "CellPq");
static_assert(sizeof(GeminalTemplate) == 12 && alignof(GeminalTemplate) == 4, "GeminalTemplate");
static_assert(sizeof(LagrangianPair) == 8 && alignof(LagrangianPair) == 4, "LagrangianPair");
static_assert(sizeof(coupling::FockIdxPart) == 12 && alignof(coupling::FockIdxPart) == 4, "FockIdxPart");
static_assert(sizeof(coupling::EriIdxPart) == 16 && alignof(coupling::EriIdxPart) == 4, "EriIdxPart");
static_assert(sizeof(coupling::ValRow) == 8 && alignof(coupling::ValRow) == 8, "ValRow");
static_assert(sizeof(coupling::LagrTerm) == 12 && alignof(coupling::LagrTerm) == 4, "LagrTerm");
static_assert(sizeof(WeightRow) == 8 && alignof(WeightRow) == 4, "WeightRow");
static_assert(sizeof(RdmRow) == 12 && alignof(RdmRow) == 4, "RdmRow");
static_assert(sizeof(DeltaSpec) == 40 && alignof(DeltaSpec) == 4, "DeltaSpec");
static_assert(sizeof(SectorEntry) == 112 && alignof(SectorEntry) == 8, "SectorEntry");
static_assert(offsetof(SectorEntry, off_sa_default) == 24, "SectorEntry.off_sa_default");
// make_catalog memcpy's the header out of the blob by sizeof, so its size/layout is ABI.
static_assert(sizeof(CatalogData) == 632 && alignof(CatalogData) == 8, "CatalogData");
static_assert(std::is_trivially_copyable<CatalogData>::value, "CatalogData trivially copyable");
static_assert(kFonLinear == 0 && kFonSqrtGroup == 1 && kFonFinterp == 2, "FonKind values");
static_assert(static_cast<int>(DeltaKind::K2b_Fock) == 0 &&
              static_cast<int>(DeltaKind::K2b_ERI) == 1 &&
              static_cast<int>(DeltaKind::K2c_ERI) == 2, "DeltaKind values");

}  // namespace studio
}  // namespace reks
}  // namespace psi
