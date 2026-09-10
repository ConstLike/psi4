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

#include "reks_arch.h"
#include "reks_layout_hash.h"  // generated: constexpr uint64_t REKS_LAYOUT_HASH

#include "reks_cassette.h"
#include "reks_report_level.h"
#include "reks_studio_eval.h"
#include "reks_si_types.h"
#include "reks_math.h"

#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <atomic>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace psi {
namespace reks {
namespace studio {

namespace coupling {

// The packed value-record word a ValRow (and a LagrTerm) carries: the low
// REKS_COEFF_BITS are the row's index into the coefficient palette, the next
// REKS_FON_N_BITS the length of its FON window, the next REKS_FON_OFF_BITS the
// window's start in fon_idx, and the sign bit stays clear.
static_assert(REKS_COEFF_BITS + REKS_FON_N_BITS + REKS_FON_OFF_BITS <= 63,
              "a value record's word leaves the sign bit clear");
constexpr int64_t kCoeffMask = (int64_t{1} << REKS_COEFF_BITS) - 1;
constexpr int64_t kFonNMask  = (int64_t{1} << REKS_FON_N_BITS) - 1;

inline int word_coeff_id(int64_t word) { return static_cast<int>(word & kCoeffMask); }
inline int word_fon_n(int64_t word) {
    return static_cast<int>((word >> REKS_COEFF_BITS) & kFonNMask);
}
inline int64_t word_fon_off(int64_t word) {
    return word >> (REKS_COEFF_BITS + REKS_FON_N_BITS);
}

// Reassemble a Fock row from its two sub-pools, selected by the (idxpart_id,
// valpart_id) pack entry, with the coefficient resolved through the palette.
inline FockTerm gather_fock(const FockIdxPart* idxpart, const ValRow* valpart,
                            const int* pack, const Coeff* coeffs, int t) {
    const FockIdxPart& ip = idxpart[pack[2 * t]];
    const ValRow&      vp = valpart[pack[2 * t + 1]];
    return FockTerm{ip.L, ip.p, ip.q, coeffs[word_coeff_id(vp.word)].coeff,
                    word_fon_off(vp.word), word_fon_n(vp.word)};
}

inline EriTerm gather_eri(const EriIdxPart* idxpart, const ValRow* valpart,
                          const int* pack, const Coeff* coeffs, int t) {
    const EriIdxPart& ip = idxpart[pack[2 * t]];
    const ValRow&     vp = valpart[pack[2 * t + 1]];
    return EriTerm{ip.p, ip.q, ip.r, ip.s, coeffs[word_coeff_id(vp.word)].coeff,
                   word_fon_off(vp.word), word_fon_n(vp.word)};
}

}  // namespace coupling

namespace {

std::filesystem::path catalog_dir() {
    return std::filesystem::path(Process::environment.get_datadir()) / "reks_catalog";
}
std::filesystem::path catalog_blob_path(int N, int M) {
    return catalog_dir() / ("reks" + std::to_string(N) + std::to_string(M) + ".bin");
}

}  // namespace

// Must not map or validate the blob here.
bool variant_supported(int N, int M) {
    std::error_code ec;
    return std::filesystem::exists(catalog_blob_path(N, M), ec);
}

// List installed active spaces by scanning reks_catalog/ for reks{N}{M}.bin.
// N == M, so the active-space digits split in half.
std::string supported_variants_str() {
    std::error_code ec;
    std::vector<std::string> found;
    for (const auto& entry : std::filesystem::directory_iterator(catalog_dir(), ec)) {
        const std::string name = entry.path().filename().string();
        if (name.rfind("reks", 0) != 0) continue;
        const auto bin_pos = name.rfind(".bin");
        if (bin_pos == std::string::npos || bin_pos <= 4) continue;
        if (bin_pos + 4 != name.size()) continue;  // a sidecar beside the blob is not a blob
        const std::string nm = name.substr(4, bin_pos - 4);
        if (nm.size() < 2 || nm.size() % 2 != 0) continue;
        if (nm.find_first_not_of("0123456789") != std::string::npos) continue;
        found.push_back("[" + nm.substr(0, nm.size() / 2) + "," + nm.substr(nm.size() / 2) + "]");
    }
    std::sort(found.begin(), found.end());
    std::string out;
    for (size_t i = 0; i < found.size(); ++i) {
        if (i) out += ", ";
        out += found[i];
    }
    return out.empty() ? "(none found in " + catalog_dir().string() + ")" : out;
}

namespace {

// Blob preamble constants; the loader validates them and fails closed on mismatch.
constexpr uint32_t kBlobMagic     = 0x31534B52u;  // 'RKS1'
constexpr uint32_t kBlobVersion   = 10u;
constexpr uint32_t kEndianTag     = 0x04030201u;
// Layout hash from reks_layout.json (build-generated into reks_layout_hash.h).
constexpr uint64_t kReksLayoutHash = REKS_LAYOUT_HASH;
// magic+version+endian+reserved (u32*4) + hash (u64) + blob_size (u64). The reserved
// word puts both 8-byte values on their own alignment and leaves the header 8-aligned.
constexpr size_t   kPreambleSize   = 32;

void validate_preamble(const uint8_t* base, size_t size) {
    // Reject a truncated file before any header read goes out of bounds: the
    // preamble + the fixed CatalogData header must both fit.
    if (size < kPreambleSize + sizeof(CatalogData))
        throw PSIEXCEPTION("REKS catalog blob: file too small (truncated)");
    uint32_t magic, version, endian;
    uint64_t hash, blob_size;
    std::memcpy(&magic,     base + 0,  4);
    std::memcpy(&version,   base + 4,  4);
    std::memcpy(&endian,    base + 8,  4);
    std::memcpy(&hash,      base + 16, 8);
    std::memcpy(&blob_size, base + 24, 8);
    if (magic != kBlobMagic)       throw PSIEXCEPTION("REKS catalog blob: bad magic");
    if (version != kBlobVersion)   throw PSIEXCEPTION("REKS catalog blob: format version mismatch");
    if (endian != kEndianTag)      throw PSIEXCEPTION("REKS catalog blob: endian mismatch");
    if (hash != kReksLayoutHash)   throw PSIEXCEPTION("REKS catalog blob: layout hash mismatch");
    if (blob_size != size)         throw PSIEXCEPTION("REKS catalog blob: size mismatch");
}


// Byte stride of the three record pools, sized from the header's own dimensions:
// one occupation byte per active orbital per spin string, one window slot per FON
// generation, one geminal count per FON generation. Mirrors reks_layout.json
// Microstate/DiagDeps/PairShape.
size_t microstate_stride(const CatalogData* h) { return 2u * h->n_active_orbitals; }
size_t diag_deps_stride(const CatalogData* h)  { return 16u + 16u * h->n_generations; }
size_t pair_shape_stride(const CatalogData* h) { return 32u + 4u * h->n_generations; }

// Byte cap of one varint element: five seven-bit groups cover int32 after zigzag.
// One LEB128 zigzag varint: seven value bits per byte, low group first, bit 7 set on
// every byte but the last. Unzigzags via v = (acc>>1) ^ -(acc&1); `off` advances past
// the element.
long long varint_read(const uint8_t* data, uint64_t& off) {
    uint64_t acc = 0;
    int shift = 0;
    for (;;) {
        const uint8_t byte = data[off++];
        acc |= static_cast<uint64_t>(byte & 0x7Fu) << shift;
        if (byte < 0x80u) break;
        shift += 7;
    }
    return static_cast<long long>((acc >> 1) ^ (~(acc & 1ull) + 1ull));
}

// `n` values of one varint window appended to `out`; `delta` decodes first
// differences, the mode every ascending window is written in.
void varint_window(const uint8_t* data, uint64_t& off, long long n, bool delta,
                   std::vector<int>& out) {
    int prev = 0;
    for (long long k = 0; k < n; ++k) {
        int v = static_cast<int>(varint_read(data, off));
        if (delta) { v += prev; prev = v; }
        out.push_back(v);
    }
}

// Advances `off` past `n` varint elements, parsed but not delta-folded or stored.
void varint_skip(const uint8_t* data, uint64_t& off, long long n) {
    for (long long k = 0; k < n; ++k) varint_read(data, off);
}


// Read-only mmap of a catalog blob file; every accessor resolves offsets off
// the base pointer.
class MappedBlob {
    const uint8_t* p_ = nullptr;
    size_t         n_ = 0;

  public:
    MappedBlob() = default;
    MappedBlob(const MappedBlob&) = delete;
    MappedBlob& operator=(const MappedBlob&) = delete;
    MappedBlob(MappedBlob&& o) noexcept : p_(o.p_), n_(o.n_) { o.p_ = nullptr; o.n_ = 0; }
    MappedBlob& operator=(MappedBlob&& o) noexcept {
        if (this != &o) {
            reset();
            p_ = o.p_; n_ = o.n_;
            o.p_ = nullptr; o.n_ = 0;
        }
        return *this;
    }
    ~MappedBlob() { reset(); }

    void reset() {
        if (p_ != nullptr) ::munmap(const_cast<uint8_t*>(p_), n_);
        p_ = nullptr;
        n_ = 0;
    }

    void map(const std::string& path) {
        reset();
        const int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0)
            throw PSIEXCEPTION("REKS catalog blob not found: " + path +
                               " (expected an installed reks_catalog/*.bin in the Psi4 data dir)");
        struct stat st {};
        if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
            ::close(fd);
            throw PSIEXCEPTION("REKS catalog blob: cannot size " + path);
        }
        void* addr = ::mmap(nullptr, static_cast<size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
        ::close(fd);  // the mapping keeps its own reference
        if (addr == MAP_FAILED)
            throw PSIEXCEPTION("REKS catalog blob: mmap failed: " + path);
        p_ = static_cast<const uint8_t*>(addr);
        n_ = static_cast<size_t>(st.st_size);
    }

    const uint8_t* data() const { return p_; }
    size_t         size() const { return n_; }
};

struct CatalogRegistryEntry {
    MappedBlob           buf;
    CatalogData          hdr;
    // Widened copies of the three pools the blob stores at its own (variable) stride.
    std::vector<Microstate>             microstates;
    std::vector<DiagDeps>               diag_deps;
    std::vector<coupling::PairShape>    pair_shapes;
};

// One spin string of a microstate record into its fixed-width occupation array; the
// tail up to kMaxActiveOrb stays zero.
void copy_occupations(std::array<int8_t, kMaxActiveOrb>& dst, const uint8_t* src, int n) {
    if (n < 0 || n > kMaxActiveOrb)
        throw PSIEXCEPTION("REKS catalog blob: microstate string of " + std::to_string(n) +
                           " orbitals exceeds kMaxActiveOrb=" + std::to_string(kMaxActiveOrb) +
                           ", this build's bound");
    for (int k = 0; k < n; ++k) dst[k] = static_cast<int8_t>(src[k]);
}

// Elements the five value channels occupy: e_pack is plain-coded at two ints per
// element, the other four are one each.
long long pair_value_elements(const coupling::PairShape& sh) {
    return 2ll * sh.n_e + sh.n_fock + sh.n_lagr + sh.n_eri + sh.n_s;
}

// Widen the three variable-stride pools into the fixed-layout records the runtime
// indexes. The blob carries n_active_orbitals occupations per spin string and one
// window slot per FON generation it declares; the tail up to kMaxActiveOrb / kMaxGen
// stays zero, which reads as an empty window everywhere it is looped over.
void decode_records(CatalogRegistryEntry& e) {
    const uint8_t*     base = e.buf.data();
    const CatalogData* h    = &e.hdr;
    const int n_orb = static_cast<int>(h->n_active_orbitals);
    const int n_gen = static_cast<int>(h->n_generations);
    const size_t ms_stride = microstate_stride(h);
    const size_t dd_stride = diag_deps_stride(h);

    e.microstates.assign(h->n_microstates, Microstate{});
    for (long long i = 0; i < h->n_microstates; ++i) {
        const uint8_t* r = base + h->off_microstates + i * ms_stride;
        Microstate& m = e.microstates[i];
        copy_occupations(m.alpha, r, n_orb);
        copy_occupations(m.beta, r + n_orb, n_orb);
    }

    e.diag_deps.assign(h->n_configs, DiagDeps{});
    for (long long K = 0; K < h->n_configs; ++K) {
        const uint8_t* r = base + h->off_diag_deps + K * dd_stride;
        DiagDeps& d = e.diag_deps[K];
        std::memcpy(&d.mw_off, r, 2 * sizeof(uint64_t));
        for (int g = 0; g < n_gen; ++g)
            std::memcpy(&d.gem[g], r + 16 + 16 * g, 2 * sizeof(uint64_t));
    }

    // Shape palette: (has_scaling, six window counts, one geminal count per
    // generation, the delta count).
    const size_t shape_stride = pair_shape_stride(h);
    const auto* shapes = h->n_pair_shapes > 0
                             ? reinterpret_cast<const int*>(base + h->off_pair_shapes)
                             : nullptr;
    e.pair_shapes.assign(static_cast<size_t>(h->n_pair_shapes), coupling::PairShape{});
    for (long long s = 0; s < h->n_pair_shapes; ++s) {
        const int* row = reinterpret_cast<const int*>(
            reinterpret_cast<const uint8_t*>(shapes) + s * shape_stride);
        // Widen the blob's n_gen geminal counts to kMaxGen; the tail reads as empty.
        coupling::PairShape& ps = e.pair_shapes[s];
        ps.has_scaling = row[0];
        ps.n_e    = row[1];
        ps.n_fock = row[2];
        ps.n_lagr = row[3];
        ps.n_eri  = row[4];
        ps.n_s    = row[5];
        ps.deps.fock_n = static_cast<uint32_t>(row[6]);
        for (int g = 0; g < n_gen; ++g)
            ps.deps.gem_n[g] = static_cast<uint32_t>(row[7 + g]);
        ps.deps.delta_n = static_cast<uint32_t>(row[7 + n_gen]);
    }

}


// Unbound Catalog over a registry entry; empty decoded pools null out, the same
// sentinel the blob-resolved pools use.
Catalog bind_catalog(const CatalogRegistryEntry& e) {
    Catalog cat;
    cat.base         = e.buf.data();
    cat.data         = &e.hdr;
    cat.microstates_ = e.microstates.empty() ? nullptr : e.microstates.data();
    cat.diag_deps_   = e.diag_deps.empty()   ? nullptr : e.diag_deps.data();
    cat.pair_shapes_ = e.pair_shapes.empty() ? nullptr : e.pair_shapes.data();
    return cat;
}






// One byte chunk of a pack region, pulled forward onto the first element starting
// inside it. `count` elements begin in [begin, next chunk's begin) and `first` is where
// they land on the region's element axis. The region's bytes alone fix all three.
struct PackChunk {
    uint64_t  begin = 0;
    long long count = 0;
    long long first = 0;
};

// What one band adds to a chunk: `sum[c]` folds the deltas of the chunk's elements in
// column c, `org[c]` is the column's running value where the chunk opens.
struct PackBandState {
    long long sum[2]  = {0, 0};
    long long org[2]  = {0, 0};
    bool      decoded = false;
};

// A varint closes on a byte with bit 7 clear, so element boundaries are recoverable
// from any byte offset: an element starts at off0 and after every such byte. That makes
// the census possible, and with it the column-1 origin without decoding column 0.
//
// Chunk table of one pack region, built on first use and kept for the process. The
// region's bytes fix the table, so every band decoded from that region reuses it: the
// sweep reads the whole region, which the bands do not. The sweep also carries the
// region's whole-extent checks, and a region that fails one is never tabled.
const std::vector<PackChunk>& pack_census(const uint8_t* base, uint64_t off0,
                                          long long nbytes) {
    // Nodes are pointer-stable and never erased, so a table stays valid once returned.
    static std::mutex mtx;
    static std::map<std::pair<const uint8_t*, uint64_t>, std::vector<PackChunk>> tables;

    std::lock_guard<std::mutex> lock(mtx);
    const auto key = std::make_pair(base, off0);
    const auto known = tables.find(key);
    if (known != tables.end()) return known->second;

    const uint64_t end = off0 + static_cast<uint64_t>(nbytes);
    if (nbytes <= 0) return tables.emplace(key, std::vector<PackChunk>{}).first->second;

    // Chunks sized for load balance across the pool, not for cache: the sweep below reads
    // every byte, the decode after it touches only the chunks the bands reach.
    constexpr long long kChunkBytes = 1 << 18;
    const int n_chunks = static_cast<int>(
        std::min<long long>(std::max<long long>(1, nbytes / kChunkBytes), 8192));
    std::vector<PackChunk> ch(n_chunks);
    ch[0].begin = off0;
    for (int k = 1; k < n_chunks; ++k) {
        const uint64_t at = off0 + static_cast<uint64_t>((nbytes * k) / n_chunks);
        uint64_t j = at - 1;
        while (j < end && base[j] >= 0x80u) ++j;
        ch[k].begin = std::max(ch[k - 1].begin, (j < end) ? j + 1 : end);
    }

    // Terminator bytes count the elements a chunk holds.
#pragma omp parallel for schedule(static)
    for (int k = 0; k < n_chunks; ++k) {
        const uint64_t hi = (k + 1 < n_chunks) ? ch[k + 1].begin : end;
        long long cnt = 0;
        for (uint64_t j = ch[k].begin; j < hi; ++j)
            if (base[j] < 0x80u) ++cnt;
        ch[k].count = cnt;
    }
    long long total = 0;
    for (int k = 0; k < n_chunks; ++k) {
        ch[k].first = total;
        total += ch[k].count;
    }
    return tables.emplace(key, std::move(ch)).first->second;
}

// Rows [lo, lo + n) of one selector column, decoded into `out` as interleaved
// (idxpart id, valpart id) pairs. The region holds the full idxpart column then the
// full valpart column, each a delta chain from row 0. Rows [lo, lo+n) are the sector's
// own, so each wanted band is a PREFIX of its column and the chain reaching it is
// self-contained: rows past lo+n carry nothing the sector reads.
void decode_pack_rows(const uint8_t* base, long long rows, uint64_t off0, long long nbytes,
                      long long lo, long long n, std::vector<int>& out, const char* what) {
    out.clear();
    if (rows <= 0) return;

    const std::vector<PackChunk>& ch = pack_census(base, off0, nbytes);
    if (n == 0) return;
    out.assign(2 * static_cast<size_t>(n), 0);

    const int n_chunks = static_cast<int>(ch.size());
    std::vector<PackBandState> st(n_chunks);
    const long long col_lo[2]  = {0, rows};
    const long long want_hi[2] = {lo + n, rows + lo + n};

    // Only the chunks the two bands reach are decoded; values come out relative to the
    // chunk, which the fold below lifts onto the column.
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < n_chunks; ++k) {
        const long long g0 = ch[k].first;
        const long long g1 = g0 + ch[k].count;
        if (g1 == g0) continue;
        if (!(g0 < want_hi[0] || (g1 > col_lo[1] && g0 < want_hi[1]))) continue;
        st[k].decoded = true;
        uint64_t off = ch[k].begin;
        long long run[2] = {0, 0};
        for (long long g = g0; g < g1; ++g) {
            const int c = (g < rows) ? 0 : 1;
            const long long t = g - col_lo[c];
            run[c] += varint_read(base, off);
            if (t < lo || t >= lo + n) continue;
            out[2 * (t - lo) + c] = static_cast<int>(run[c]);
        }
        st[k].sum[0] = run[0];
        st[k].sum[1] = run[1];
    }

    // Chunks past a band are never decoded and fold in as zero, and no decoded chunk
    // follows them in that column, so each column's running value stays exact.
    long long acc[2] = {0, 0};
    for (int k = 0; k < n_chunks; ++k) {
        st[k].org[0] = acc[0];
        st[k].org[1] = acc[1];
        acc[0] += st[k].sum[0];
        acc[1] += st[k].sum[1];
    }

    // The chunk-relative value each row carries is lifted onto its column.
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < n_chunks; ++k) {
        if (!st[k].decoded) continue;
        const long long g0 = ch[k].first;
        const long long g1 = g0 + ch[k].count;
        for (int c = 0; c < 2; ++c) {
            // t is this chunk's rows in column c, clipped to the decoded band [lo, lo+n).
            const long long t_lo = std::max(std::max(g0 - col_lo[c], 0LL), lo);
            const long long t_hi = std::min(std::min(g1 - col_lo[c], rows), lo + n);
            for (long long t = t_lo; t < t_hi; ++t)
                out[2 * (t - lo) + c] = static_cast<int>(out[2 * (t - lo) + c] + st[k].org[c]);
        }
    }
}

}  // namespace

SelectorPacks decode_selector_packs(const Catalog& cat, const SelectorDemand& demand) {
    const uint8_t*     base = cat.base;
    const CatalogData* h    = cat.data;
    const SectorEntry* se = reinterpret_cast<const SectorEntry*>(base + h->off_sectors);
    const int b = cat.active_sector;

    // fock_pack and eri_pack are independent byte regions, decoded one after the
    // other so each decode_pack_rows call spreads across the whole thread pool.
    SelectorPacks p;
    decode_pack_rows(base, h->n_fock_pack_rows, h->off_fock_pack, h->fock_pack_bytes,
                     static_cast<long long>(se[b].base_fock_pack), demand.fock_rows,
                     p.fock, "fock_pack");
    decode_pack_rows(base, h->n_eri_pack_rows, h->off_eri_pack, h->eri_pack_bytes,
                     static_cast<long long>(se[b].base_eri_pack), demand.eri_rows,
                     p.eri, "eri_pack");
    return p;
}

void scan_pair_values(const uint8_t* base, uint64_t& off,
                      const coupling::PairShape& sh,
                      long long& fock_rows, long long& eri_rows) {
    // Channel order is reks_layout.json pair_windows: e_pack (plain, two ints per
    // element), then the four delta-coded runs. Only fock_idx and eri_idx are folded;
    // the rest advance `off` and are dropped.
    varint_skip(base, off, 2ll * sh.n_e);
    long long run = 0;
    for (int t = 0; t < sh.n_fock; ++t) {
        run += varint_read(base, off);
        fock_rows = std::max(fock_rows, run + 1);
    }
    varint_skip(base, off, sh.n_lagr);
    run = 0;
    for (int t = 0; t < sh.n_eri; ++t) {
        run += varint_read(base, off);
        eri_rows = std::max(eri_rows, run + 1);
    }
    varint_skip(base, off, sh.n_s);
}

coupling::PairValues decode_pair_values(const uint8_t* base, uint64_t& off,
                                        const coupling::PairShape& sh,
                                        std::vector<int>& scratch) {
    scratch.clear();
    scratch.reserve(static_cast<size_t>(pair_value_elements(sh)));
    // Channel order is the one reks_layout.json pair_windows fixes; e_pack is
    // plain-coded and holds two ints per element, every other window is a delta-coded
    // ascending run.
    varint_window(base, off, 2 * sh.n_e, false, scratch);
    varint_window(base, off, sh.n_fock, true, scratch);
    varint_window(base, off, sh.n_lagr, true, scratch);
    varint_window(base, off, sh.n_eri, true, scratch);
    varint_window(base, off, sh.n_s, true, scratch);

    size_t at = 0;
    auto take = [&](int n) {
        const IntSpan sp{ n ? scratch.data() + at : nullptr, n };
        at += static_cast<size_t>(n);
        return sp;
    };
    coupling::PairValues v;
    v.e_pack   = take(2 * sh.n_e);
    v.fock_idx = take(sh.n_fock);
    v.lagr_idx = take(sh.n_lagr);
    v.eri_idx  = take(sh.n_eri);
    v.s_idx    = take(sh.n_s);
    return v;
}

PairDepsView decode_pair_deps(const uint8_t* base, uint64_t& off,
                              const coupling::PairShape& sh, std::vector<int>& scratch) {
    scratch.clear();
    varint_window(base, off, sh.deps.fock_n, true, scratch);
    for (int g = 0; g < kMaxGen; ++g)
        varint_window(base, off, sh.deps.gem_n[g], true, scratch);
    varint_window(base, off, sh.deps.delta_n, true, scratch);

    size_t at = 0;
    auto take = [&](uint32_t n) {
        const IntSpan sp{ n ? scratch.data() + at : nullptr, static_cast<int>(n) };
        at += n;
        return sp;
    };
    PairDepsView v;
    v.fock_footprint = take(sh.deps.fock_n);
    for (int g = 0; g < kMaxGen; ++g) v.geminal_reads[g] = take(sh.deps.gem_n[g]);
    v.delta_reads = take(sh.deps.delta_n);
    return v;
}

char generation_letter(int gen) {
    static const char letters[kMaxGen] = {'n', 'm', 'u', 'v', 'w', 'x', 'y', 'z'};
    if (gen < 0 || gen >= kMaxGen)
        throw PSIEXCEPTION("generation_letter: gen " + std::to_string(gen) +
                           " out of range [0," + std::to_string(kMaxGen) + ")");
    return letters[gen];
}

std::string fon_array_name(int gen) {
    if (gen == 0) return "REKS FON";
    const char up = static_cast<char>(generation_letter(gen) - 'a' + 'A');
    return std::string("REKS ") + up + "_FON";
}

std::string fon_channel_name(int gen, int run_sector, int scheme) {
    std::string s(1, generation_letter(gen));
    if (scheme >= 0) s += std::to_string(scheme);
    // run_sector 0 adds no suffix: single-sector display stays byte-stable.
    if (run_sector > 0) s += "t" + std::to_string(run_sector);
    return s;
}

// Blob sector index whose spin2 == spin, or throw with the available manifolds.
static int resolve_sector(const Catalog& cat, const std::filesystem::path& path, int spin) {
    const SectorEntry* se = cat.sectors(cat.base);
    std::string avail;
    for (int b = 0; b < cat.data->n_sectors; ++b) {
        if (se[b].spin2 == spin) return b;
        if (b) avail += ", ";
        avail += std::to_string(se[b].spin2);
    }
    throw PSIEXCEPTION("REKS catalog blob " + path.string() + " has no 2S=" +
                       std::to_string(spin) + " spin manifold (available 2S: " + avail + ").");
}

// Loads (or returns the cached) merged (N,M) blob as an unbound Catalog
// (active_sector == -1); a sector must be bound before use.
static Catalog load_catalog_blob(int N, int M) {
    // std::map nodes are pointer-stable; buf.data()/&hdr stay valid for process
    // lifetime. Buffer is immutable after init; reads are lock-free. One merged blob
    // per (N,M) carries every spin sector.
    static std::mutex mtx;
    static std::map<std::pair<int, int>, CatalogRegistryEntry> registry;

    std::lock_guard<std::mutex> lock(mtx);
    const auto key = std::make_pair(N, M);
    auto it = registry.find(key);
    const std::filesystem::path path = catalog_blob_path(N, M);
    if (it == registry.end()) {
        CatalogRegistryEntry entry;
        // mmap's base is page-aligned, so a pool offset that is itself 8-aligned
        // reaches an 8-aligned address.
        entry.buf.map(path.string());

        // Fail-closed: reject any magic/version/endian/layout-hash/size mismatch.
        validate_preamble(entry.buf.data(), entry.buf.size());
        std::memcpy(&entry.hdr, entry.buf.data() + kPreambleSize, sizeof(CatalogData));
        // The filename is the only binding of bytes to (N,M); cross-check the
        // active-space identity against the header.
        if (entry.hdr.n_electrons != N)
            throw PSIEXCEPTION("REKS catalog blob: electron-count mismatch for " + path.string() +
                               " (header n_electrons=" + std::to_string(entry.hdr.n_electrons) +
                               ", requested N=" + std::to_string(N) + ")");
        if (entry.hdr.n_active_orbitals != M)
            throw PSIEXCEPTION("REKS catalog blob: active-space mismatch for " + path.string() +
                               " (header n_active_orbitals=" + std::to_string(entry.hdr.n_active_orbitals) +
                               ", requested M=" + std::to_string(M) + ")");
        it = registry.emplace(key, std::move(entry)).first;
        decode_records(it->second);
    }

    return bind_catalog(it->second);
}

Catalog make_catalog(int N, int M, int spin) {
    Catalog cat = load_catalog_blob(N, M);
    cat.active_sector = resolve_sector(cat, catalog_blob_path(N, M), spin);
    return cat;
}

Catalog make_catalog(int N, int M) {
    Catalog cat = load_catalog_blob(N, M);
    cat.active_sector = 0;
    return cat;
}

std::vector<std::vector<int>> default_si_indices(const Catalog& s) {
    const int start = s.sector_config_start();
    const int count = s.sector_config_count();
    std::vector<int> Ks;
    Ks.reserve(count);
    for (int i = 0; i < count; ++i) Ks.push_back(start + i);
    return { std::move(Ks) };
}

std::vector<double> uniform_weights(int n) {
    if (n <= 0) return {};
    return std::vector<double>(static_cast<size_t>(n), 1.0 / static_cast<double>(n));
}

// Global config range to search for names: the bound sector's block, or the whole
// catalog when unbound (active_sector < 0). Returns [begin, end).
static std::pair<int, int> name_scan_range(const Catalog& s) {
    if (s.active_sector < 0) return {0, s.data->n_configs};
    const int start = s.sector_config_start();
    return {start, start + s.sector_config_count()};
}

int find_config_index(const Catalog& s, const std::string& name) {
    std::string upper = name;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    const auto [begin, end] = name_scan_range(s);
    for (int K = begin; K < end; ++K) {
        const char* nm = s.si_config_name(K);
        if (nm && upper == nm) return K;
    }
    return -1;
}

std::string config_names_csv(const Catalog& s) {
    const auto [begin, end] = name_scan_range(s);
    std::string out;
    for (int K = begin; K < end; ++K) {
        if (K > begin) out += ", ";
        const char* nm = s.si_config_name(K);
        out += nm ? nm : "";
    }
    return out;
}

std::string config_type_of(const std::string& name) {
    size_t end = name.size();
    while (end > 0 && std::isdigit(static_cast<unsigned char>(name[end - 1]))) --end;
    if (end > 0 && end < name.size() && name[end - 1] == '_') --end;
    return name.substr(0, end);
}

std::string config_types_csv(const Catalog& s) {
    const auto [begin, end] = name_scan_range(s);
    std::vector<std::string> seen;
    for (int K = begin; K < end; ++K) {
        const char* nm = s.si_config_name(K);
        if (!nm) continue;
        std::string ty = config_type_of(nm);
        if (std::find(seen.begin(), seen.end(), ty) == seen.end()) seen.push_back(ty);
    }
    std::string out;
    for (size_t i = 0; i < seen.size(); ++i) {
        if (i) out += ", ";
        out += seen[i];
    }
    return out;
}

// Geminal-multiplicity markers of the catalog generator's two-orbital naming
// convention: the marker counts the geminals of its kind, "PPS" names the
// perfectly-paired remainder, a run of "T" the sector's triplet background and
// "SPS" the spin-projection tag.
struct MarkerCount {
    const char* marker;
    int         count;
};
static constexpr MarkerCount kOssMarkers[] = {
    {"OSS", 1}, {"DOSS", 2}, {"TOSS", 3}, {"QOSS", 4}, {"QUOSS", 5},
    {"SOSS", 6}, {"SPOSS", 7}, {"OOSS", 8}, {"NOSS", 9}};
static constexpr MarkerCount kDesMarkers[] = {
    {"DES", 1}, {"DDES", 2}, {"TDES", 3}, {"QDES", 4}, {"QUDES", 5},
    {"SDES", 6}, {"SPDES", 7}, {"ODES", 8}, {"NDES", 9}};

// Geminal counts a config type name carries: Psi_1 (open-shell singlet) and
// Psi_2 (doubly excited closed singlet) multiplicities, plus the SPS tag. A
// spin-projected config is built from one parent expansion and carries no Psi_2
// count, so sps forces n_psi2 to zero.
struct ConfigSignature {
    int  n_psi1 = 0;
    int  n_psi2 = 0;
    bool sps    = false;
};

// Signature of a config type name; false when a part of the name matches no
// marker of the two-orbital convention.
static bool config_signature_of(const std::string& type_name, ConfigSignature& sig) {
    sig = ConfigSignature{};
    std::string up = type_name;
    std::transform(up.begin(), up.end(), up.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    for (size_t pos = 0; pos <= up.size();) {
        const size_t sep  = up.find('_', pos);
        const size_t stop = (sep == std::string::npos) ? up.size() : sep;
        const std::string part = up.substr(pos, stop - pos);
        pos = stop + 1;
        if (part.empty()) return false;
        if (part == "SPS") {
            sig.sps = true;
            continue;
        }
        if (part == "PPS") continue;
        if (part.find_first_not_of('T') == std::string::npos) continue;
        bool matched = false;
        for (const MarkerCount& m : kOssMarkers)
            if (part == m.marker) {
                sig.n_psi1 = m.count;
                matched    = true;
                break;
            }
        if (matched) continue;
        for (const MarkerCount& m : kDesMarkers)
            if (part == m.marker) {
                sig.n_psi2 = m.count;
                matched    = true;
                break;
            }
        if (!matched) return false;
    }
    if (sig.sps) sig.n_psi2 = 0;
    return true;
}

// Pairing schemes of the catalog, numbered 0 .. n-1 over the geminal slots.
static int catalog_n_schemes(const Catalog& s) {
    const GeminalTemplate* gems = s.geminal_templates(s.base);
    if (gems == nullptr) return 0;
    int hi = -1;
    for (int64_t g = 0; g < s.data->n_geminals; ++g) hi = std::max(hi, gems[g].scheme);
    return hi + 1;
}

// Pairing scheme of a config: the scheme of any geminal slot it reads, the slots
// of one config sharing a scheme. Spin-projected configs read no slot and belong
// to no single scheme -- they span all of them -- and get -1; a one-scheme catalog
// leaves nothing to resolve, so everything there is scheme 0.
static int config_scheme_of(const Catalog& s, int K) {
    const GeminalTemplate* gems = s.geminal_templates(s.base);
    if (gems == nullptr) return 0;
    const DiagDepsView dv = s.diag_deps_view(s.base, K);
    for (int gen = 0; gen < kMaxGen; ++gen)
        for (int slot : dv.geminal_reads[gen])
            if (slot >= 0 && slot < s.data->n_geminals) return gems[slot].scheme;
    return catalog_n_schemes(s) == 1 ? 0 : -1;
}

// Name of the scanned block: "the catalog" when unbound, else the active
// sector's 2S manifold.
static std::string scan_range_label(const Catalog& s) {
    if (s.active_sector < 0) return "the catalog";
    return "the 2S=" + std::to_string(s.sector_entry(s.active_sector).spin2) + " manifold";
}

// Which configs of the listed schemes a scheme token keeps. The spin-projected
// family belongs to no single scheme and is the manifold's base: both modes take
// it whole and neither repeats it, which also keeps a one-scheme selection the
// same size whichever scheme it names.
enum class SchemeMode {
    Whole,  ///< single-scheme: the configs of the schemes, plus the base
    Mixed   ///< mixed-scheme: their single excitations, plus every open-shell singlet of the manifold
};

// The scheme tokens, spelled with '-' between the words.
struct SchemeToken {
    const char* base;
    SchemeMode  mode;
};
static constexpr SchemeToken kSchemeTokens[] = {{"MIXED-SCHEME", SchemeMode::Mixed},
                                                {"SINGLE-SCHEME", SchemeMode::Whole}};

// Scheme-token selection over the scanned block, taking the union over the
// listed schemes:
//   single-scheme:N = the configs of scheme N, plus the spin-projected base
//   mixed-scheme:N  = the single-excitation configs of scheme N, every
//                     open-shell-singlet config and the base
static std::vector<int> expand_scheme_token(const Catalog& s, SchemeMode mode,
                                            const std::vector<int>& schemes,
                                            const std::string& token) {
    if (s.data->n_orbitals_per_geminal != 2)
        throw PSIEXCEPTION("REKS config token '" + token + "' needs two-orbital geminals; "
                           "this catalog carries " +
                           std::to_string(s.data->n_orbitals_per_geminal) +
                           " orbitals per geminal and names its configs by sub-pair.");
    const int n_schemes = catalog_n_schemes(s);
    for (int scheme : schemes)
        if (scheme < 0 || scheme >= n_schemes)
            throw PSIEXCEPTION("REKS config token '" + token + "': scheme " +
                               std::to_string(scheme) + " is outside [0, " +
                               std::to_string(n_schemes) + ") of this catalog.");
    const auto [begin, end] = name_scan_range(s);
    std::vector<int> out;
    for (int K = begin; K < end; ++K) {
        const char* nm = s.si_config_name(K);
        if (nm == nullptr) continue;
        const std::string type = config_type_of(nm);
        ConfigSignature sig;
        if (!config_signature_of(type, sig))
            throw PSIEXCEPTION("REKS config token '" + token + "': config type '" + type +
                               "' carries no geminal-multiplicity marker.");
        const int scheme = config_scheme_of(s, K);
        if (scheme < 0 && !sig.sps)
            throw PSIEXCEPTION("REKS config token '" + token + "': config '" +
                               std::string(nm) +
                               "' reads no geminal and is not spin-projected, so its "
                               "pairing scheme cannot be resolved.");
        const bool in_scheme = std::find(schemes.begin(), schemes.end(), scheme) != schemes.end();
        const bool single_ex = sig.n_psi1 + sig.n_psi2 <= 1;
        const bool open_sing = sig.n_psi1 >= 1 && sig.n_psi2 == 0;
        const bool keep = (mode == SchemeMode::Mixed)
                              ? ((in_scheme && single_ex) || open_sing || sig.sps)
                              : (in_scheme || sig.sps);
        if (keep) out.push_back(K);
    }
    if (out.empty())
        throw PSIEXCEPTION("REKS config token '" + token + "' selects no configuration of " +
                           scan_range_label(s) +
                           ". An empty list declares a manifold without configs.");
    return out;
}

// Case-insensitive match of a whole name against a pattern whose '*' stands for
// any run of characters and '?' for one.
static bool glob_match(const char* pat, const char* str) {
    const char* star = nullptr;
    const char* mark = nullptr;
    auto upper = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    while (*str != '\0') {
        if (*pat == '?' || (*pat != '\0' && upper(*pat) == upper(*str))) {
            ++pat;
            ++str;
        } else if (*pat == '*') {
            star = pat++;
            mark = str;
        } else if (star != nullptr) {
            pat = star + 1;
            str = ++mark;
        } else {
            return false;
        }
    }
    while (*pat == '*') ++pat;
    return *pat == '\0';
}

std::vector<int> expand_config_token(const Catalog& s, const std::string& token) {
    std::string up = token;
    std::transform(up.begin(), up.end(), up.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    // The spin-projected base, and name patterns.
    const bool is_pattern =
        up.find('*') != std::string::npos || up.find('?') != std::string::npos;
    if (up == "SPS" || is_pattern) {
        const auto [pbegin, pend] = name_scan_range(s);
        std::vector<int> hits;
        for (int K = pbegin; K < pend; ++K) {
            const char* nm = s.si_config_name(K);
            if (nm == nullptr) continue;
            if (is_pattern) {
                if (glob_match(up.c_str(), nm)) hits.push_back(K);
                continue;
            }
            ConfigSignature sig;
            if (config_signature_of(config_type_of(nm), sig) && sig.sps) hits.push_back(K);
        }
        if (hits.empty())
            throw PSIEXCEPTION("REKS config token '" + token + "' matches no configuration of " +
                               scan_range_label(s) + ".");
        return hits;
    }

    // A scheme token is its base followed by the scheme indices, each behind
    // ':', '-', '_' or ',' -- all four equivalent. Config type names carry '_',
    // so the separators are folded on a copy. No index means scheme 0.
    std::string norm = up;
    for (char& c : norm)
        if (c == '_' || c == ':' || c == ',') c = '-';
    for (const SchemeToken& st : kSchemeTokens) {
        const std::string base(st.base);
        if (norm.rfind(base, 0) != 0) continue;
        std::vector<int> schemes;
        for (size_t pos = base.size(); pos < norm.size();) {
            const size_t stop   = norm.find('-', pos + 1);
            const size_t argend = (stop == std::string::npos) ? norm.size() : stop;
            const std::string arg = norm.substr(pos + 1, argend - pos - 1);
            if (norm[pos] != '-' || arg.empty() || arg.size() > 9 ||
                arg.find_first_not_of("0123456789") != std::string::npos)
                throw PSIEXCEPTION("REKS config token '" + token + "': '" + base +
                                   "' takes scheme indices (non-negative integers) behind "
                                   "':', '-', '_' or ',', and nothing else.");
            const int scheme = std::stoi(arg);
            if (std::find(schemes.begin(), schemes.end(), scheme) != schemes.end())
                throw PSIEXCEPTION("REKS config token '" + token + "': scheme " +
                                   std::to_string(scheme) + " is listed twice.");
            schemes.push_back(scheme);
            pos = argend;
        }
        if (schemes.empty()) schemes.push_back(0);
        return expand_scheme_token(s, st.mode, schemes, token);
    }
    const auto [begin, end] = name_scan_range(s);
    std::vector<int> out;
    if (up == "FULL") {
        out.reserve(end - begin);
        for (int K = begin; K < end; ++K) out.push_back(K);
        return out;
    }
    for (int K = begin; K < end; ++K) {
        const char* nm = s.si_config_name(K);
        if (nm && config_type_of(nm) == up) out.push_back(K);
    }
    return out;
}

int config_index_by_name(const Catalog& s, const std::string& name) {
    const int K = find_config_index(s, name);
    if (K < 0)
        throw PSIEXCEPTION("Unknown REKS config name '" + name +
                           "'. Valid names: " + config_names_csv(s) + ".");
    return K;
}

bool has_diabatic_rdm(const Catalog& cat) {
    return cat.data && cat.data->n_rdm_elems > 0;
}

std::vector<std::pair<int,int>> collect_eri_pairs(const Cassette&  cassette) {
    std::set<std::pair<int,int>> pairs;
    const int n_specs = cassette.n_delta_specs();
    for (DeltaKey d : cassette.deltas_active) {
        if (d < 0 || d >= n_specs) continue;
        const DeltaSpec& s = cassette.delta_spec(d);
        // K2b_Fock carries the {-1,-1,-1,-1} sentinel; both ERI forms read a tile.
        if (s.kind == DeltaKind::K2b_Fock) continue;
        int k = s.eri_idx[2];
        int l = s.eri_idx[3];
        if (k > l) std::swap(k, l);
        pairs.emplace(k, l);
    }
    return {pairs.begin(), pairs.end()};
}

std::vector<LagrangianPair> lagrangian_pairs_view(const Catalog& cat) {
    const LagrangianPair* lp = cat.lagrangian_pairs(cat.base);
    return std::vector<LagrangianPair>(lp, lp + cat.data->n_lagrangian_pairs);
}

namespace coupling {

// H_{KL} for one SI coupling pair: sum over e/fock/lagr/eri term classes,
//   each term c_t * prod(FONs_t) * X_t, X_t in {E_L, F_pq, lagrangian, (pq|rs)}.
double eval_pair_H(const PairValues& win,
                   const FonFactor* fon_pool,
                   const int* fon_idx,
                   const FonSlot* slot_pool,
                   const int* slot_idx,
                   const Coeff* coeff_pool,
                   const ValRow* e_val_pool,
                   const FockIdxPart* fock_idxpart,
                   const ValRow* fock_valpart,
                   const int* fock_pack,
                   const EriIdxPart* eri_idxpart,
                   const ValRow* eri_valpart,
                   const int* eri_pack,
                   const LagrTerm* lagr_row_pool,
                   const FONSnapshot& fons,
                   const std::vector<double>& E_L,
                   const std::vector<int>& fock_aa_row,
                   const std::vector<double>& fock_aa,
                   int n_active,
                   const reks::ActiveEriTiles& active_eri,
                   const std::vector<double>& lagrangians,
                   const double* fon_cache) {
    double v = 0.0;
    const int n_e = win.e_pack.n / 2;   // e_pack holds (L, val_id) per element
    for (int t = 0; t < n_e; ++t) {
        const int L = win.e_pack[2 * t];
        const ValRow& e = e_val_pool[win.e_pack[2 * t + 1]];
        v += coeff_pool[word_coeff_id(e.word)].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, word_fon_off(e.word), word_fon_n(e.word), fons, &f_interp, fon_cache)
               * E_L[L];
    }
    for (int t = 0; t < win.fock_idx.n; ++t) {
        const FockTerm f = gather_fock(fock_idxpart, fock_valpart, fock_pack, coeff_pool,
                                       win.fock_idx[t]);
        // Active block of F_MO_a - F_MO_b: row = f.p, column = f.q, both in active order.
        const double fval = fock_aa[static_cast<size_t>(fock_aa_row[f.L]) * n_active * n_active
                                    + static_cast<size_t>(f.p) * n_active + f.q];
        v += f.coeff * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, f.fon_off, f.fon_n, fons, &f_interp, fon_cache) * fval;
    }
    for (int t = 0; t < win.lagr_idx.n; ++t) {
        const LagrTerm& w = lagr_row_pool[win.lagr_idx[t]];
        v += coeff_pool[word_coeff_id(w.word)].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, word_fon_off(w.word), word_fon_n(w.word), fons, &f_interp, fon_cache)
               * lagrangians[w.w];
    }
    for (int t = 0; t < win.eri_idx.n; ++t) {
        const EriTerm g = gather_eri(eri_idxpart, eri_valpart, eri_pack, coeff_pool,
                                     win.eri_idx[t]);
        v += g.coeff * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, g.fon_off, g.fon_n, fons, &f_interp, fon_cache)
               * active_eri.at(g.p, g.q, g.r, g.s);
    }
    return v;
}

// S_{KL} for one SI overlap pair: sum over pure-FON terms, each term c_t * prod(FONs_t).
double eval_pair_S(const IntSpan& s_idx, const FonFactor* fon_pool,
                   const int* fon_idx, const FonSlot* slot_pool, const int* slot_idx,
                   const Coeff* coeff_pool, const ValRow* s_row_pool,
                   const FONSnapshot& fons, const double* fon_cache) {
    double v = 0.0;
    for (int t = 0; t < s_idx.n; ++t) {
        const ValRow& s = s_row_pool[s_idx[t]];
        v += coeff_pool[word_coeff_id(s.word)].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, word_fon_off(s.word), word_fon_n(s.word), fons, &f_interp, fon_cache);
    }
    return v;
}

}  // namespace coupling

}  // namespace studio
}  // namespace reks
}  // namespace psi
