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

/// Runtime bound on FON generations n, m, u, v, w, x, y, z (REKS up to (16,16)); not ABI, the loader
/// rejects a catalog above it.
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

/// Diagonal-block dependency record: (byte_offset, count) windows into a shared int pool, gem slots
/// widened to kMaxGen at load (unused tail empty).
struct DiagDeps {
    uint64_t mw_off, mw_n;                        ///< microstate_writes window
    struct { uint64_t off, n; } gem[kMaxGen];     ///< [gen] -> geminal_reads window
};

/// Non-owning view of a run of integers.
template <class T>
struct Span {
    const T* data = nullptr;
    int      n    = 0;
    constexpr const T* begin() const noexcept { return data; }
    constexpr const T* end()   const noexcept { return data + n; }
    constexpr bool     empty() const noexcept { return n == 0; }
    constexpr int      size()  const noexcept { return n; }
    constexpr T        operator[](int i) const noexcept { return data[i]; }
};
/// Ids into bounded tables (microstates, geminal slots, deltas).
using IntSpan = Span<int>;
/// Ids into the pools that grow with the config pairs.
using IdSpan = Span<int64_t>;

/// Off-diagonal dependency record: lengths of a pair's three dependency
/// windows, read from the pair's PairShape row. Widened to kMaxGen like DiagDeps.
struct PairDeps {
    uint32_t fock_n;              ///< fock_footprint window length
    uint32_t gem_n[kMaxGen];      ///< [gen] -> geminal_reads window length
    uint32_t delta_n;             ///< delta_reads window length
};

/// Base-resolved view of a DiagDeps record.
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

/// Runtime bound on active orbitals (REKS up to (12,12)); not ABI, the loader rejects a catalog above it.
inline constexpr int kMaxActiveOrb = 12;

/// Slater microstate: alpha[i], beta[i] in {0,1} occupy active orbital i < n_active; the tail is zero.
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

/// Spin manifold of a blob (ascending spin2), configs [config_start, config_start + config_count); its
/// pair pools count from zero, base_* gives the absolute row.
struct SectorEntry {
    int      spin2;          ///< total spin as 2S for this manifold
    int      config_start;   ///< blob config index of the first config in this sector
    int      config_count;   ///< number of configs in this sector
    int      n_generations;  ///< FON generations of this sector's configs
    int      n_sa_default;
    uint64_t off_sa_default;      ///< default SA pool (blob config indices)
    uint64_t base_e_val;          ///< row origin in e_val_pool
    uint64_t base_fock_idxpart;
    uint64_t base_fock_valpart;
    uint64_t base_fock_pack;      ///< ROW origin in fock_pack (2 ids per row)
    uint64_t base_eri_idxpart;
    uint64_t base_eri_valpart;
    uint64_t base_eri_pack;       ///< ROW origin in eri_pack (2 ids per row)
    uint64_t base_lagr_row;
    uint64_t base_s_row;
    uint64_t base_fon_idx;        ///< origin in fon_idx; fon_win entries already include it
    uint64_t base_fon_win;        ///< origin in fon_win of the value rows' window ids
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
    double  coeff;
    int64_t fon_off;     // window into the FonFactor pool
    int     fon_n;
};

// One entry of the catalog's coefficient palette.
struct Coeff { double coeff; };

// Block / vector output: RDM (cells over (p,q)) AND weights (cells over (config, L)).
struct Cell { int p, q; int64_t term_off; int term_n; };   // RDM: out[p*N+q]; weight: p=L, q=-1

// Interned (p, q) part of a Cell.
struct CellPq  { int p, q; };

#pragma pack(push, 4)
// One cell: its term window in cell_win and its (p, q) in cell_pq.
struct CellSel { int64_t win; int pq; };
// One diagonal term: its FON window in term_win and its coefficient in coeff_pool.
struct TermSel { int64_t win; int coeff; };
// One row per config, in config order (row index == config K).
struct WeightRow { int64_t cell_off; int cell_n; };
#pragma pack(pop)

// A 1-RDM row keyed (K_i, K_j); K_i is the CSR row (see Catalog::rdm_row_ptr).
struct RdmRow { int64_t cell_off; int key_j; int cell_n; };

// One evaluated RDM cell at flat position pq = p*N + q; operator< sorts a
// block into ascending (p, q) order.
struct RdmCell {
    int    pq;
    double val;
    bool operator<(const RdmCell& o) const { return pq < o.pq; }
};

// Pools shared by every block row; fon_idx selects DIAG_FON_IDX rows, cell_win/term_win are CSR
// boundaries (window w = [win[w], win[w+1])).
struct BlockPools {
    const CellPq*    cell_pq;
    const int64_t*   cell_win;
    const CellSel*   cell_pack;
    const Coeff*     coeff_pool;
    const int64_t*   term_win;
    const TermSel*   term_pack;
    const FonFactor* fon_pool;
    const FonSlot*   fon_slot_pool;
    const int*       fon_idx;
    const int*       fon_slot_idx;
};

// Reassemble Cell c: (p, q) from cell_pq and its term window from cell_win.
inline Cell load_cell(const BlockPools& bp, int64_t c) {
    const CellSel& s = bp.cell_pack[c];
    const CellPq&  pq = bp.cell_pq[s.pq];
    return Cell{pq.p, pq.q, bp.cell_win[s.win],
                static_cast<int>(bp.cell_win[s.win + 1] - bp.cell_win[s.win])};
}

// Reassemble Term t: the coefficient from the palette and the FON window from term_win.
inline Term load_term(const BlockPools& bp, int64_t t) {
    const TermSel& s = bp.term_pack[t];
    return Term{bp.coeff_pool[s.coeff].coeff, bp.term_win[s.win],
                static_cast<int>(bp.term_win[s.win + 1] - bp.term_win[s.win])};
}

// Off-diagonal H_ij = <K_i|H|K_j> coupling IR: 5-channel decomposition
// (energy, Fock, Lagrangian, ERI, overlap-FON).
namespace coupling {

// The value part every coupling channel shares: the row's FON window id
// (sector-local, into fon_win) and its coefficient id in the palette.
#pragma pack(push, 4)
struct ValRow { int64_t win; int coeff; };
#pragma pack(pop)

// coeff * prod(fons) * (F^alpha_L - F^beta_L)[p, q]
struct FockTerm{ int L; int p; int q;
                 double coeff; int64_t fon_off; int fon_n; };
// coeff * prod(fons) * W[w];  w = lagrangian_pairs index.
struct LagrTerm{ int64_t win; int coeff; int w; };
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
    IdSpan e_pack;            ///< 2*n_e ids: (L, val_id) pairs into e_val_pool
    IdSpan fock_idx;          ///< selector rows into the fock pack
    IdSpan lagr_idx;          ///< rows of lagr_row_pool
    IdSpan eri_idx;           ///< selector rows into the eri pack
    IdSpan s_idx;             ///< rows of s_row_pool
};

}  // namespace coupling

/// Advance `off` past a pair's value channels, raising `fock_rows`/`eri_rows` to the running maximum
/// selector row (signed deltas are not monotone); nothing is materialised.
void scan_pair_values(const uint8_t* base, uint64_t& off,
                      const coupling::PairShape& sh,
                      long long& fock_rows, long long& eri_rows);

/// Decode a pair's five value channels from the blob, advancing `off` to its dependency windows; L maps
/// through `l_global`, the spans live in `scratch` until its next use.
coupling::PairValues decode_pair_values(const uint8_t* base, uint64_t& off,
                                        const coupling::PairShape& sh, const int* l_global,
                                        std::vector<int64_t>& scratch);

/// Decode a pair's three dependency windows up to the block end; L maps through `l_global`, delta ids
/// shift by `delta_shift`, scratch as in decode_pair_values.
PairDepsView decode_pair_deps(const uint8_t* base, uint64_t& off,
                              const coupling::PairShape& sh, const int* l_global,
                              int delta_shift, std::vector<int>& scratch);

// Offset-ABI catalog header (POD): int64 counts, uint64 byte offsets absolute from the blob base with
// their lengths; resolved through the Catalog accessors.
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

    // Pair section: columns over n_pairs config pairs, a PairShape palette, one varint block per pair in
    // pair_win_data; n_pair_rows == n_configs + 1.
    int64_t n_pairs;
    int64_t n_pair_rows;     uint64_t off_pair_row_ptr;   // int64[n_pair_rows], CSR over K_i
    uint64_t off_pair_key_j;                              // int32[n_pairs]
    uint64_t off_pair_shape;                              // uint32[n_pairs], index into pair_shapes
    uint64_t off_pair_win_off;                            // uint64[n_pairs], absolute byte offsets
    int64_t n_pair_shapes;   uint64_t off_pair_shapes;    // (8 + n_generations) int32 per shape
    uint64_t off_pair_win_data;  uint64_t pair_win_bytes;

    int64_t n_fon_pool;      uint64_t off_fon_pool;     // FON_VALUE_POOL (coupling + block fon_idx select)
    int64_t n_fon_idx;       uint64_t off_fon_idx;      // FON_IDX coupling row selector
    // Coupling FON windows: int64[n_fon_win + 1] absolute offsets into fon_idx;
    // window g spans [fon_win[g], fon_win[g+1]).
    int64_t n_fon_win;       uint64_t off_fon_win;
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
    int64_t n_rdm_rows;      uint64_t off_rdm_row_ptr;    // int64[n_rdm_rows], CSR over K_i
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

/// One mapped blob of a shelf with its config/microstate/delta maps onto shelf numbering; renumbered or
/// self-sized pools are decoded at load, the rest stays mapped.
struct CatalogFile {
    const uint8_t* base = nullptr;   ///< mapped blob root (KIND-1 offset resolution)
    CatalogData    hdr{};            ///< owned header copy
    int            k_shift     = 0;  ///< shelf K = file K + k_shift
    int            delta_shift = 0;  ///< shelf delta id = file delta id + delta_shift
    std::vector<int>                   l_global;      ///< file L -> shelf L
    std::vector<coupling::PairShape>   pair_shapes;   ///< shape palette, widened to kMaxGen
    std::vector<CellPq>                cell_pq;       ///< weight cells (q == -1) carry shelf L in p
    std::vector<coupling::FockIdxPart> fock_idxpart;  ///< L is shelf L
};

/// One spin manifold of a shelf; sectors ascend in 2S. Its configs occupy the shelf
/// block [config_start, config_start + config_count).
struct SectorInfo {
    int              spin2         = 0;
    int              config_start  = 0;
    int              config_count  = 0;
    int              n_generations = 0;
    std::vector<int> sa_default;       ///< default SA pool, shelf K
    int              file          = 0;  ///< shelf file holding the sector
    int              file_sector   = 0;  ///< its row in that file's sector table
};

/// Every catalog blob of one (N,M), numbered as one catalog.
struct CatalogShelf;

/// Shelf view bound to one sector: `base`/`data` are its file, ids come out in shelf numbering.
struct Catalog {
    const CatalogShelf* shelf = nullptr;
    const CatalogFile*  file  = nullptr;  ///< file of the bound sector
    const uint8_t*     base = nullptr;    ///< file->base
    const CatalogData* data = nullptr;    ///< &file->hdr
    int                active_sector = -1;  ///< bound shelf sector; -1 = unbound
    int                file_sector   = -1;  ///< the bound sector's row in its file
    const Microstate*  microstates_  = nullptr;  ///< shelf microstate table

    template <class T>
    const T* pool_(const uint8_t* b, uint64_t off, int64_t count) const {
        return count > 0 ? reinterpret_cast<const T*>(b + off) : nullptr;
    }

    /// This view rebound to shelf sector b.
    Catalog bound_to(int b) const;

    /// Sector table entry of the bound sector in its file; the base_* origins live here.
    const SectorEntry& bound_sector(const uint8_t* b) const { return sectors(b)[file_sector]; }

    /// Pool base at the bound sector's first row; `stride` = elements per row (2 for selector packs).
    template <class T>
    const T* sector_pool_(const T* pool, uint64_t base_row, int64_t stride = 1) const {
        return pool != nullptr ? pool + base_row * stride : nullptr;
    }

    /// Shelf microstate table, indexed by shelf L.
    const Microstate*             microstates() const { return microstates_; }

    /// Pairs of config K (shelf K, in the bound file) as [first, last) in CSR order,
    /// partners ascending; empty when the blob carries no pairs.
    std::pair<int64_t, int64_t> pair_row(int K) const {
        const int64_t* rp = pool_<int64_t>(base, data->off_pair_row_ptr, data->n_pair_rows);
        if (rp == nullptr || data->n_pairs <= 0) return {0, 0};
        const int k = K - file->k_shift;
        return {rp[k], rp[k + 1]};
    }
    /// Shelf K of the partner config of pair p.
    int pair_partner(int64_t p) const {
        return reinterpret_cast<const int*>(base + data->off_pair_key_j)[p] + file->k_shift;
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
        return coupling::PairRef{ pair_partner(p), file->pair_shapes.data() + pair_shape_id()[p] };
    }

    /// Value channels of pair p, decoded on demand (see decode_pair_values).
    coupling::PairValues pair_values(long long p, std::vector<int64_t>& scratch) const {
        uint64_t off = pair_win_off()[p];
        return decode_pair_values(base, off, file->pair_shapes[pair_shape_id()[p]],
                                  file->l_global.data(), scratch);
    }

    /// Dependency windows of pair p, decoded on demand; raises `fock_rows`/`eri_rows` to its selector reach.
    PairDepsView pair_deps(long long p, std::vector<int>& scratch,
                           long long& fock_rows, long long& eri_rows) const {
        const coupling::PairShape& sh = file->pair_shapes[pair_shape_id()[p]];
        uint64_t off = pair_win_off()[p];
        scan_pair_values(base, off, sh, fock_rows, eri_rows);
        return decode_pair_deps(base, off, sh, file->l_global.data(), file->delta_shift,
                                scratch);
    }

    const GeminalTemplate* geminal_templates(const uint8_t* b) const {
        return pool_<GeminalTemplate>(b, data->off_geminal_templates, data->n_geminals);
    }
    const LagrangianPair* lagrangian_pairs(const uint8_t* b) const {
        return pool_<LagrangianPair>(b, data->off_lagrangian_pairs, data->n_lagrangian_pairs);
    }
    const FonFactor* fon_pool(const uint8_t* b) const {
        return pool_<FonFactor>(b, data->off_fon_pool, data->n_fon_pool);
    }

    /// Coupling FON factor ids; the fon_win offsets into it are absolute.
    const int* fon_idx(const uint8_t* b) const {
        return pool_<int>(b, data->off_fon_idx, data->n_fon_idx);
    }

    // Per-sector pools (channel, selector column and FON window ids count from zero per sector), rebased
    // to the bound sector's first row.
    const int64_t* fon_win(const uint8_t* b) const {
        return sector_pool_(pool_<int64_t>(b, data->off_fon_win, data->n_fon_win),
                            bound_sector(b).base_fon_win);
    }
    const coupling::FockIdxPart* fock_idxpart(const uint8_t* b) const {
        return sector_pool_(file->fock_idxpart.empty() ? nullptr : file->fock_idxpart.data(),
                            bound_sector(b).base_fock_idxpart);
    }
    const coupling::ValRow* fock_valpart(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::ValRow>(b, data->off_fock_valpart, data->n_fock_valpart),
                            bound_sector(b).base_fock_valpart);
    }
    const coupling::EriIdxPart* eri_idxpart(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::EriIdxPart>(b, data->off_eri_idxpart, data->n_eri_idxpart),
                            bound_sector(b).base_eri_idxpart);
    }
    const coupling::ValRow* eri_valpart(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::ValRow>(b, data->off_eri_valpart, data->n_eri_valpart),
                            bound_sector(b).base_eri_valpart);
    }
    const coupling::LagrTerm* lagr_row_pool(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::LagrTerm>(b, data->off_lagr_row_pool, data->n_lagr_pool),
                            bound_sector(b).base_lagr_row);
    }
    const coupling::ValRow* e_val_pool(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::ValRow>(b, data->off_e_val_pool, data->n_e_val_pool),
                            bound_sector(b).base_e_val);
    }
    const coupling::ValRow* s_row_pool(const uint8_t* b) const {
        return sector_pool_(pool_<coupling::ValRow>(b, data->off_s_row_pool, data->n_s_row_pool),
                            bound_sector(b).base_s_row);
    }

    /// Coefficient palette, numbered over the whole catalog (no sector origin).
    const Coeff* coeff_pool(const uint8_t* b) const {
        return pool_<Coeff>(b, data->off_coeff_pool, data->n_coeff_pool);
    }

    const WeightRow* weight_elems(const uint8_t* b) const {
        return pool_<WeightRow>(b, data->off_weight_elems, data->n_weight_elems);
    }
    const RdmRow* rdm_elems(const uint8_t* b) const {
        return pool_<RdmRow>(b, data->off_rdm_elems, data->n_rdm_elems);
    }
    /// rdm_elems rows of config K (shelf K, in the bound file) as [first, last),
    /// partners ascending; empty when the blob carries no RDM.
    std::pair<int64_t, int64_t> rdm_row(int K) const {
        const int64_t* rp = pool_<int64_t>(base, data->off_rdm_row_ptr, data->n_rdm_rows);
        if (rp == nullptr || data->n_rdm_elems <= 0) return {0, 0};
        const int k = K - file->k_shift;
        return {rp[k], rp[k + 1]};
    }
    /// Shelf K of the partner config of rdm_elems row r.
    int rdm_partner(int64_t r) const { return rdm_elems(base)[r].key_j + file->k_shift; }
    const FonSlot* fon_slot_pool(const uint8_t* b) const {
        return pool_<FonSlot>(b, data->off_fon_slot_pool, data->n_fon_slot_pool);
    }
    const int* fon_slot_idx(const uint8_t* b) const {
        return pool_<int>(b, data->off_fon_slot_idx, data->n_fon_slot_idx);
    }
    const int* diag_fon_idx(const uint8_t* b) const {
        return pool_<int>(b, data->off_diag_fon_idx, data->n_diag_fon_idx);
    }
    const CellPq* cell_pq(const uint8_t*) const {
        return file->cell_pq.empty() ? nullptr : file->cell_pq.data();
    }
    const int64_t* cell_win(const uint8_t* b) const {
        return pool_<int64_t>(b, data->off_cell_win, data->n_cell_win);
    }
    const CellSel* cell_pack(const uint8_t* b) const {
        return pool_<CellSel>(b, data->off_cell_pack, data->n_cell_pack);
    }
    const int64_t* term_win(const uint8_t* b) const {
        return pool_<int64_t>(b, data->off_term_win, data->n_term_win);
    }
    const TermSel* term_pack(const uint8_t* b) const {
        return pool_<TermSel>(b, data->off_term_pack, data->n_term_pack);
    }
    const SectorEntry* sectors(const uint8_t* b) const {
        return pool_<SectorEntry>(b, data->off_sectors, data->n_sectors);
    }

    // Shelf-wide dimensions and sector table; b is a shelf sector in
    // [0, n_sectors_available()), K a shelf config.
    int n_sectors_available() const;
    const SectorInfo& sector_entry(int b) const;
    /// Shelf sector owning config K; -1 if K is out of range.
    int sector_of_config(int K) const;
    /// Default SA pool of shelf sector b (shelf K).
    const std::vector<int>& sector_sa_default(int b) const;
    int n_microstates() const;
    int n_configs() const;
    int n_delta_specs() const;
    /// Delta spec of shelf delta id d; its microstates are shelf L.
    const DeltaSpec& delta_spec(int d) const;
    /// True iff some file of the shelf carries a diabatic RDM.
    bool has_rdm() const;

    // Bound-sector (active_sector) entries.
    int sector_config_start()  const { return sector_entry(active_sector).config_start; }
    int sector_config_count()  const { return sector_entry(active_sector).config_count; }
    int sector_n_generations() const { return sector_entry(active_sector).n_generations; }

    BlockPools block_pools(const uint8_t* b) const {
        return BlockPools{ cell_pq(b), cell_win(b), cell_pack(b), coeff_pool(b),
                           term_win(b), term_pack(b), fon_pool(b), fon_slot_pool(b),
                           diag_fon_idx(b), fon_slot_idx(b) };
    }
    /// Diagonal dependencies of config K; microstate writes are shelf L.
    DiagDepsView diag_deps_view(int K) const;

    /// A config's weight row, paired with the block pools it indexes into.
    struct WeightBlock {
        BlockPools       bp{};
        const WeightRow* row = nullptr;
    };
    /// Weight row of config K with the pools of the file holding it; row is null when
    /// that file carries no weights or K is out of range.
    WeightBlock weight_block(int K) const;

    /// Display name / wavefunction def of config K; nullptr when absent.
    const char* si_config_name(int K) const;
    const char* si_config_def(int K) const;
};

/// How far into the fock/eri selector columns, from the bound sector's first
/// row, a caller's pairs need decoded.
struct SelectorDemand {
    long long fock_rows = 0;
    long long eri_rows  = 0;
};

/// Fock/ERI selector rows (idxpart id, valpart id) up to the demand; row 0 = the bound sector's first.
struct SelectorPacks {
    std::vector<int64_t> fock;   ///< 2 ids per row
    std::vector<int64_t> eri;
};

/// Decode each selector column of cat's bound sector as far as `demand`
/// reaches; only the demanded rows are returned.
SelectorPacks decode_selector_packs(const Catalog& cat, const SelectorDemand& demand);

/// Cached shelf of (N,M) (reks{N}{M}.bin and/or reks{N}{M}_2spin{2S}.bin) bound to shelf
/// sector 0, the lowest available 2S; throws if a blob is invalid or the blobs disagree.
Catalog make_catalog(int N, int M);

/// Display letter for a FON generation index: 0->'n', 1->'m', 2->'u', 3->'v',
/// 4->'w', 5->'x', 6->'y', 7->'z'. Throws if gen is outside [0, kMaxGen).
char generation_letter(int gen);

/// FON channel name: generation letter, scheme digit if scheme >= 0, "t<run_sector>" if run_sector > 0
/// (gen 2, run_sector 2, scheme 1 -> "u1t2").
std::string fon_channel_name(int gen, int run_sector, int scheme = -1);

/// Wavefunction array-variable base name for a FON generation: gen 0 -> "REKS FON",
/// gen >= 1 -> "REKS <LETTER>_FON" (M_FON, U_FON, V_FON, W_FON).
std::string fon_array_name(int gen);

/// True iff reks{N}{M}.bin or some reks{N}{M}_2spin{2S}.bin is installed (filename
/// probe; the spin manifolds are resolved and checked in make_catalog).
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

/// Bulk token (grammar: option SA_REKS_CONFIGS) -> the bound sector's global config indices in catalog
/// order; empty for an unknown type, throws on a malformed scheme token.
std::vector<int> expand_config_token(const Catalog& s, const std::string& token);

/// Uniform weights vector of length n; each entry = 1.0 / n. Returns empty
/// vector when n <= 0.
std::vector<double> uniform_weights(int n);

/// True iff the catalog provides a diabatic RDM (some file has n_rdm_elems > 0).
bool has_diabatic_rdm(const Catalog& cat);

/// Active ERI pairs (eri_idx[2], eri_idx[3]) of the cassette's K2b/K2c_ERI deltas, k <= l, sorted.
std::vector<std::pair<int,int>> collect_eri_pairs(const Cassette&   cassette);

/// Owning-vector copy of the catalog's lagrangian_pairs (length n_lagrangian_pairs).
std::vector<LagrangianPair> lagrangian_pairs_view(const Catalog& cat);

// ABI lock: fixed-stride struct geometry and constants == reks_layout.json (source of layout_hash);
// Microstate, DiagDeps and PairShape are decoded at load and are not ABI.
static_assert(sizeof(FonFactor) == 16 && alignof(FonFactor) == 4, "FonFactor");
static_assert(sizeof(FonSlot) == 12 && alignof(FonSlot) == 4, "FonSlot");
static_assert(sizeof(Coeff) == 8 && alignof(Coeff) == 8, "Coeff");
static_assert(sizeof(CellPq) == 8 && alignof(CellPq) == 4, "CellPq");
static_assert(sizeof(GeminalTemplate) == 12 && alignof(GeminalTemplate) == 4, "GeminalTemplate");
static_assert(sizeof(LagrangianPair) == 8 && alignof(LagrangianPair) == 4, "LagrangianPair");
static_assert(sizeof(coupling::FockIdxPart) == 12 && alignof(coupling::FockIdxPart) == 4, "FockIdxPart");
static_assert(sizeof(coupling::EriIdxPart) == 16 && alignof(coupling::EriIdxPart) == 4, "EriIdxPart");
static_assert(sizeof(coupling::ValRow) == 12 && alignof(coupling::ValRow) == 4, "ValRow");
static_assert(sizeof(coupling::LagrTerm) == 16 && alignof(coupling::LagrTerm) == 8, "LagrTerm");
static_assert(sizeof(CellSel) == 12 && alignof(CellSel) == 4, "CellSel");
static_assert(sizeof(TermSel) == 12 && alignof(TermSel) == 4, "TermSel");
static_assert(sizeof(WeightRow) == 12 && alignof(WeightRow) == 4, "WeightRow");
static_assert(sizeof(RdmRow) == 16 && alignof(RdmRow) == 8, "RdmRow");
static_assert(sizeof(DeltaSpec) == 40 && alignof(DeltaSpec) == 4, "DeltaSpec");
static_assert(sizeof(SectorEntry) == 120 && alignof(SectorEntry) == 8, "SectorEntry");
static_assert(offsetof(SectorEntry, off_sa_default) == 24, "SectorEntry.off_sa_default");
// make_catalog memcpy's the header out of the blob by sizeof, so its size/layout is ABI.
static_assert(sizeof(CatalogData) == 648 && alignof(CatalogData) == 8, "CatalogData");
static_assert(std::is_trivially_copyable<CatalogData>::value, "CatalogData trivially copyable");
static_assert(kFonLinear == 0 && kFonSqrtGroup == 1 && kFonFinterp == 2, "FonKind values");
static_assert(static_cast<int>(DeltaKind::K2b_Fock) == 0 &&
              static_cast<int>(DeltaKind::K2b_ERI) == 1 &&
              static_cast<int>(DeltaKind::K2c_ERI) == 2, "DeltaKind values");

}  // namespace studio
}  // namespace reks
}  // namespace psi
