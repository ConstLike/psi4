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

#include "reks_arch.h"
#include "reks_layout_hash.h"  // generated: constexpr uint64_t REKS_LAYOUT_HASH

#include "reks_cassette.h"
#include "reks_studio_eval.h"
#include "reks_si_types.h"
#include "reks_math.h"

#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _MSC_VER
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace psi {
namespace reks {
namespace studio {

namespace coupling {

// Span of window `win` in fon_idx; fon_win is rebased to the bound sector.
struct FonWindow {
    int64_t off;
    int     n;
};
inline FonWindow fon_window(const int64_t* fon_win, int64_t win) {
    return FonWindow{fon_win[win], static_cast<int>(fon_win[win + 1] - fon_win[win])};
}

// Reassemble a Fock row from its two sub-pools, selected by the (idxpart_id,
// valpart_id) pack entry, with the coefficient resolved through the palette.
inline FockTerm gather_fock(const FockIdxPart* idxpart, const ValRow* valpart,
                            const int64_t* pack, const Coeff* coeffs,
                            const int64_t* fon_win, int64_t t) {
    const FockIdxPart& ip = idxpart[pack[2 * t]];
    const ValRow&      vp = valpart[pack[2 * t + 1]];
    const FonWindow    fw = fon_window(fon_win, vp.win);
    return FockTerm{ip.L, ip.p, ip.q, coeffs[vp.coeff].coeff, fw.off, fw.n};
}

inline EriTerm gather_eri(const EriIdxPart* idxpart, const ValRow* valpart,
                          const int64_t* pack, const Coeff* coeffs,
                          const int64_t* fon_win, int64_t t) {
    const EriIdxPart& ip = idxpart[pack[2 * t]];
    const ValRow&     vp = valpart[pack[2 * t + 1]];
    const FonWindow   fw = fon_window(fon_win, vp.win);
    return EriTerm{ip.p, ip.q, ip.r, ip.s, coeffs[vp.coeff].coeff, fw.off, fw.n};
}

}  // namespace coupling

namespace {

// $PSI4_REKS_CATALOG if set and non-empty, else <PSIDATADIR>/reks_catalog.
std::filesystem::path catalog_dir() {
    const char* env = std::getenv("PSI4_REKS_CATALOG");
    if (env != nullptr && *env != '\0') return std::filesystem::path(env);
    return std::filesystem::path(Process::environment.get_datadir()) / "reks_catalog";
}

// Active-space digits and spin tag of a blob name: "reks{N}{M}.bin" -> empty tag,
// "reks{N}{M}_2spin{2S}.bin" -> 2S; false for any other name.
bool parse_blob_name(const std::string& name, std::string& nm, std::string& spin_tag) {
    if (name.size() < 8 || name.rfind("reks", 0) != 0) return false;
    if (name.compare(name.size() - 4, 4, ".bin") != 0) return false;
    std::string stem = name.substr(4, name.size() - 8);
    spin_tag.clear();
    const auto tag = stem.find("_2spin");
    if (tag != std::string::npos) {
        spin_tag = stem.substr(tag + 6);
        if (spin_tag.empty() || spin_tag.find_first_not_of("0123456789") != std::string::npos)
            return false;
        stem.resize(tag);
    }
    if (stem.size() < 2 || stem.size() % 2 != 0) return false;
    if (stem.find_first_not_of("0123456789") != std::string::npos) return false;
    nm = stem;
    return true;
}

// Blobs of (N,M) in reks_catalog/ with their spin tags, in path order.
std::vector<std::pair<std::filesystem::path, std::string>> catalog_blob_paths(int N, int M) {
    const std::string want = std::to_string(N) + std::to_string(M);
    std::vector<std::pair<std::filesystem::path, std::string>> out;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(catalog_dir(), ec)) {
        std::string nm, tag;
        if (parse_blob_name(entry.path().filename().string(), nm, tag) && nm == want)
            out.emplace_back(entry.path(), tag);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

// Must not map or validate the blobs here.
bool variant_supported(int N, int M) { return !catalog_blob_paths(N, M).empty(); }

// List installed active spaces by scanning reks_catalog/ for reks{N}{M}.bin and
// reks{N}{M}_2spin{2S}.bin. N == M, so the active-space digits split in half.
std::string supported_variants_str() {
    std::error_code ec;
    std::set<std::string> found;
    for (const auto& entry : std::filesystem::directory_iterator(catalog_dir(), ec)) {
        std::string nm, tag;
        if (!parse_blob_name(entry.path().filename().string(), nm, tag)) continue;
        found.insert("[" + nm.substr(0, nm.size() / 2) + "," + nm.substr(nm.size() / 2) + "]");
    }
    std::string out;
    for (const std::string& v : found) {
        if (!out.empty()) out += ", ";
        out += v;
    }
    return (out.empty() ? std::string("none") : out) + " (catalog directory " + catalog_dir().string() +
           "; PSI4_REKS_CATALOG overrides it)";
}

namespace {

// Blob preamble constants; the loader validates them and fails closed on mismatch.
constexpr uint32_t kBlobMagic     = 0x31534B52u;  // 'RKS1'
constexpr uint32_t kBlobVersion   = 11u;
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

// Byte strides of the Microstate/DiagDeps/PairShape pools from the header: one byte per active orbital
// per spin string, one window slot and one geminal count per FON generation.
size_t microstate_stride(const CatalogData* h) { return 2u * h->n_active_orbitals; }
size_t diag_deps_stride(const CatalogData* h)  { return 16u + 16u * h->n_generations; }
size_t pair_shape_stride(const CatalogData* h) { return 32u + 4u * h->n_generations; }

// LEB128 zigzag varint (7 bits per byte, low first, bit 7 = continue): v = (acc >> 1) ^ -(acc & 1);
// `off` advances past the element.
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

// Append `n` values of one varint window to `out`; `delta` decodes first differences (ascending
// windows). T = int64_t for pair-scaling pool ids, int for bounded tables.
template <class T>
void varint_window(const uint8_t* data, uint64_t& off, long long n, bool delta,
                   std::vector<T>& out) {
    long long prev = 0;
    for (long long k = 0; k < n; ++k) {
        long long v = varint_read(data, off);
        if (delta) { v += prev; prev = v; }
        out.push_back(static_cast<T>(v));
    }
}

// Advances `off` past `n` varint elements, parsed but not delta-folded or stored.
void varint_skip(const uint8_t* data, uint64_t& off, long long n) {
    for (long long k = 0; k < n; ++k) varint_read(data, off);
}

}  // namespace

// Read-only mapping of a catalog blob (mmap on POSIX, MapViewOfFile on Windows).
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
#ifdef _MSC_VER
        if (p_ != nullptr) ::UnmapViewOfFile(p_);
#else
        if (p_ != nullptr) ::munmap(const_cast<uint8_t*>(p_), n_);
#endif
        p_ = nullptr;
        n_ = 0;
    }

    void map(const std::string& path) {
        reset();
#ifdef _MSC_VER
        const HANDLE file = ::CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                                          OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (file == INVALID_HANDLE_VALUE)
            throw PSIEXCEPTION("REKS catalog blob: cannot open " + path);
        LARGE_INTEGER size {};
        if (!::GetFileSizeEx(file, &size) || size.QuadPart <= 0) {
            ::CloseHandle(file);
            throw PSIEXCEPTION("REKS catalog blob: cannot size " + path);
        }
        const HANDLE mapping = ::CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
        ::CloseHandle(file);  // the mapping keeps its own reference
        if (mapping == nullptr)
            throw PSIEXCEPTION("REKS catalog blob: mmap failed: " + path);
        void* addr = ::MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
        ::CloseHandle(mapping);  // the view keeps its own reference
        if (addr == nullptr)
            throw PSIEXCEPTION("REKS catalog blob: mmap failed: " + path);
        p_ = static_cast<const uint8_t*>(addr);
        n_ = static_cast<size_t>(size.QuadPart);
#else
        const int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0)
            throw PSIEXCEPTION("REKS catalog blob: cannot open " + path);
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
#endif
    }

    const uint8_t* data() const { return p_; }
    size_t         size() const { return n_; }
};

namespace {

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

// Widen the variable-stride pools into fixed-stride records in file numbering (shape palette into `f`,
// microstates and diagonal records into the buffers); tails up to kMaxActiveOrb / kMaxGen stay zero.
void decode_records(CatalogFile& f, std::vector<Microstate>& microstates,
                    std::vector<DiagDeps>& diag_deps) {
    const uint8_t*     base = f.base;
    const CatalogData* h    = &f.hdr;
    const int n_orb = static_cast<int>(h->n_active_orbitals);
    const int n_gen = static_cast<int>(h->n_generations);
    const size_t ms_stride = microstate_stride(h);
    const size_t dd_stride = diag_deps_stride(h);

    microstates.assign(h->n_microstates, Microstate{});
    for (long long i = 0; i < h->n_microstates; ++i) {
        const uint8_t* r = base + h->off_microstates + i * ms_stride;
        Microstate& m = microstates[i];
        copy_occupations(m.alpha, r, n_orb);
        copy_occupations(m.beta, r + n_orb, n_orb);
    }

    diag_deps.assign(h->n_configs, DiagDeps{});
    for (long long K = 0; K < h->n_configs; ++K) {
        const uint8_t* r = base + h->off_diag_deps + K * dd_stride;
        DiagDeps& d = diag_deps[K];
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
    f.pair_shapes.assign(static_cast<size_t>(h->n_pair_shapes), coupling::PairShape{});
    for (long long s = 0; s < h->n_pair_shapes; ++s) {
        const int* row = reinterpret_cast<const int*>(
            reinterpret_cast<const uint8_t*>(shapes) + s * shape_stride);
        // Widen the blob's n_gen geminal counts to kMaxGen; the tail reads as empty.
        coupling::PairShape& ps = f.pair_shapes[s];
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

// Sector table of a blob.
const SectorEntry* file_sectors(const CatalogFile& f) {
    return reinterpret_cast<const SectorEntry*>(f.base + f.hdr.off_sectors);
}

// Spin manifolds a blob holds, as "0, 2, 4".
std::string blob_spins_csv(const CatalogFile& f) {
    const SectorEntry* se = file_sectors(f);
    std::string out;
    for (int64_t b = 0; b < f.hdr.n_sectors; ++b) {
        if (b) out += ", ";
        out += std::to_string(se[b].spin2);
    }
    return out;
}

// Map one blob of (N,M), check it, and decode its variable-stride records. A blob
// named for one spin manifold (non-empty spin_tag) must hold exactly that sector.
void open_blob(const std::filesystem::path& path, const std::string& spin_tag, int N, int M,
               MappedBlob& map, CatalogFile& f, std::vector<Microstate>& microstates,
               std::vector<DiagDeps>& diag_deps) {
    // mmap's base is page-aligned, so a pool offset that is itself 8-aligned
    // reaches an 8-aligned address.
    map.map(path.string());

    // Fail-closed: reject any magic/version/endian/layout-hash/size mismatch.
    validate_preamble(map.data(), map.size());
    f.base = map.data();
    std::memcpy(&f.hdr, f.base + kPreambleSize, sizeof(CatalogData));
    // The filename is the only binding of bytes to (N,M) and to a spin manifold;
    // cross-check both against the header.
    if (f.hdr.n_electrons != N)
        throw PSIEXCEPTION("REKS catalog blob: electron-count mismatch for " + path.string() +
                           " (header n_electrons=" + std::to_string(f.hdr.n_electrons) +
                           ", requested N=" + std::to_string(N) + ")");
    if (f.hdr.n_active_orbitals != M)
        throw PSIEXCEPTION("REKS catalog blob: active-space mismatch for " + path.string() +
                           " (header n_active_orbitals=" + std::to_string(f.hdr.n_active_orbitals) +
                           ", requested M=" + std::to_string(M) + ")");
    if (!spin_tag.empty()) {
        if (f.hdr.n_sectors != 1 || std::to_string(file_sectors(f)[0].spin2) != spin_tag)
            throw PSIEXCEPTION("REKS catalog blob " + path.string() + " is named for the 2S=" +
                               spin_tag + " manifold but holds 2S = " + blob_spins_csv(f));
    }
    decode_records(f, microstates, diag_deps);
}

// Byte chunk of a pack region aligned to its first element start: `count` elements begin in
// [begin, next begin), `first` is the index of the first one.
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

// Chunk table of a pack region, built once per process: an element starts at off0 and after every byte
// with bit 7 clear, so chunks (and the column-1 origin) are found without decoding.
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
    const int nthread = Process::environment.get_n_threads();
#pragma omp parallel for schedule(static) num_threads(nthread)
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

// Rows [lo, lo + n) of a selector region (idxpart column, then valpart column, each a delta chain from
// row 0) into `out` as interleaved (idxpart, valpart) pairs; rows past lo + n are not decoded.
void decode_pack_rows(const uint8_t* base, long long rows, uint64_t off0, long long nbytes,
                      long long lo, long long n, std::vector<int64_t>& out) {
    out.clear();
    if (rows <= 0) return;

    const std::vector<PackChunk>& ch = pack_census(base, off0, nbytes);
    if (n == 0) return;
    out.assign(2 * static_cast<size_t>(n), 0);

    const int n_chunks = static_cast<int>(ch.size());
    std::vector<PackBandState> st(n_chunks);
    const long long col_lo[2]  = {0, rows};
    const long long want_hi[2] = {lo + n, rows + lo + n};
    const int nthread = Process::environment.get_n_threads();

    // Only the chunks the two bands reach are decoded; values come out relative to the
    // chunk, which the fold below lifts onto the column.
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread)
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
            out[2 * (t - lo) + c] = run[c];
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
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread)
    for (int k = 0; k < n_chunks; ++k) {
        if (!st[k].decoded) continue;
        const long long g0 = ch[k].first;
        const long long g1 = g0 + ch[k].count;
        for (int c = 0; c < 2; ++c) {
            // t is this chunk's rows in column c, clipped to the decoded band [lo, lo+n).
            const long long t_lo = std::max(std::max(g0 - col_lo[c], 0LL), lo);
            const long long t_hi = std::min(std::min(g1 - col_lo[c], rows), lo + n);
            for (long long t = t_lo; t < t_hi; ++t)
                out[2 * (t - lo) + c] += st[k].org[c];
        }
    }
}

}  // namespace

SelectorPacks decode_selector_packs(const Catalog& cat, const SelectorDemand& demand) {
    const uint8_t*     base = cat.base;
    const CatalogData* h    = cat.data;
    const SectorEntry& se   = cat.bound_sector(base);

    SelectorPacks p;
    decode_pack_rows(base, h->n_fock_pack_rows, h->off_fock_pack, h->fock_pack_bytes,
                     static_cast<long long>(se.base_fock_pack), demand.fock_rows, p.fock);
    decode_pack_rows(base, h->n_eri_pack_rows, h->off_eri_pack, h->eri_pack_bytes,
                     static_cast<long long>(se.base_eri_pack), demand.eri_rows, p.eri);
    return p;
}

void scan_pair_values(const uint8_t* base, uint64_t& off,
                      const coupling::PairShape& sh,
                      long long& fock_rows, long long& eri_rows) {
    // pair_windows order: e_pack (plain, 2 ints per element), four delta runs; only fock_idx and eri_idx
    // are folded.
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
                                        const coupling::PairShape& sh, const int* l_global,
                                        std::vector<int64_t>& scratch) {
    scratch.clear();
    scratch.reserve(static_cast<size_t>(pair_value_elements(sh)));
    // pair_windows order: e_pack (plain, 2 ints per element), then ascending delta runs.
    varint_window(base, off, 2 * sh.n_e, false, scratch);
    for (int t = 0; t < sh.n_e; ++t) scratch[2 * t] = l_global[scratch[2 * t]];
    varint_window(base, off, sh.n_fock, true, scratch);
    varint_window(base, off, sh.n_lagr, true, scratch);
    varint_window(base, off, sh.n_eri, true, scratch);
    varint_window(base, off, sh.n_s, true, scratch);

    size_t at = 0;
    auto take = [&](int n) {
        const IdSpan sp{ n ? scratch.data() + at : nullptr, n };
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
                              const coupling::PairShape& sh, const int* l_global,
                              int delta_shift, std::vector<int>& scratch) {
    scratch.clear();
    varint_window(base, off, sh.deps.fock_n, true, scratch);
    for (int& L : scratch) L = l_global[L];
    for (int g = 0; g < kMaxGen; ++g)
        varint_window(base, off, sh.deps.gem_n[g], true, scratch);
    const size_t delta_at = scratch.size();
    varint_window(base, off, sh.deps.delta_n, true, scratch);
    for (size_t k = delta_at; k < scratch.size(); ++k) scratch[k] += delta_shift;

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
    if (run_sector > 0) s += "t" + std::to_string(run_sector);
    return s;
}

// All blobs of one (N,M) in merged-blob numbering: K in sector blocks ascending 2S, L first-seen over
// the files in that order, delta ids file after file.
struct CatalogShelf {
    std::vector<MappedBlob>   maps;          ///< one mapping per file
    std::vector<CatalogFile>  files;         ///< ascending lowest 2S
    std::vector<SectorInfo>   sectors;       ///< ascending 2S
    std::vector<Microstate>   microstates;   ///< shelf L
    std::vector<int>          mw_pool;       ///< diagonal microstate writes, shelf L
    std::vector<DiagDepsView> diag;          ///< per shelf K
    std::vector<DeltaSpec>    delta_specs;   ///< shelf delta ids, shelf L
    bool                    has_rdm       = false;
};

namespace {

// Identity of a microstate across blobs: its alpha and beta occupation bitmasks.
uint32_t microstate_key(const Microstate& m, int n_active) {
    return static_cast<uint32_t>(base_density_index(m, n_active, false)) |
           (static_cast<uint32_t>(base_density_index(m, n_active, true)) << n_active);
}

// Throws unless blob f describes the same active space as blob f0: the dimensions
// every sector shares, the geminal slots and the lagrangian pairs.
void check_same_active_space(const CatalogFile& f0, const std::filesystem::path& p0,
                             const CatalogFile& f, const std::filesystem::path& p) {
    auto fail = [&](const std::string& what) {
        throw PSIEXCEPTION("REKS catalog blobs " + p0.string() + " and " + p.string() +
                           " disagree on " + what + "; they are not one catalog.");
    };
    const CatalogData& a = f0.hdr;
    const CatalogData& b = f.hdr;
    if (a.n_geminals != b.n_geminals || a.n_lagrangian_pairs != b.n_lagrangian_pairs ||
        a.n_orbitals_per_geminal != b.n_orbitals_per_geminal || a.scheme != b.scheme)
        fail("the geminal layout");
    if (a.n_geminals > 0 &&
        std::memcmp(f0.base + a.off_geminal_templates, f.base + b.off_geminal_templates,
                    static_cast<size_t>(a.n_geminals) * sizeof(GeminalTemplate)) != 0)
        fail("the geminal slots");
    if (a.n_lagrangian_pairs > 0 &&
        std::memcmp(f0.base + a.off_lagrangian_pairs, f.base + b.off_lagrangian_pairs,
                    static_cast<size_t>(a.n_lagrangian_pairs) * sizeof(LagrangianPair)) != 0)
        fail("the lagrangian pairs");
}

// Shelf of (N,M) from its installed blobs; throws on a repeated or interleaved 2S manifold or on
// different active spaces.
CatalogShelf build_shelf(int N, int M) {
    const auto paths = catalog_blob_paths(N, M);
    const int n_files = static_cast<int>(paths.size());

    std::vector<MappedBlob>              maps(n_files);
    std::vector<CatalogFile>             files(n_files);
    std::vector<std::vector<Microstate>> file_ms(n_files);
    std::vector<std::vector<DiagDeps>>   file_dd(n_files);
    for (int i = 0; i < n_files; ++i)
        open_blob(paths[i].first, paths[i].second, N, M, maps[i], files[i], file_ms[i],
                  file_dd[i]);

    for (int i = 0; i < n_files; ++i) {
        if (files[i].hdr.n_sectors <= 0)
            throw PSIEXCEPTION("REKS catalog blob " + paths[i].first.string() +
                               " holds no spin manifold");
        if (i > 0) check_same_active_space(files[0], paths[0].first, files[i], paths[i].first);
    }

    // Files in ascending order of their lowest manifold; a file's table ascends in 2S.
    std::vector<int> order(n_files);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
        return file_sectors(files[x])[0].spin2 < file_sectors(files[y])[0].spin2;
    });

    CatalogShelf shelf;
    for (int i : order) {
        shelf.maps.push_back(std::move(maps[i]));
        shelf.files.push_back(std::move(files[i]));
    }

    // Shelf sectors ascend in 2S and walk the files in order, each file's table in turn.
    struct Row { int spin2, file, file_sector; };
    std::vector<Row> rows;
    for (int f = 0; f < n_files; ++f)
        for (int s = 0; s < shelf.files[f].hdr.n_sectors; ++s)
            rows.push_back({file_sectors(shelf.files[f])[s].spin2, f, s});
    std::stable_sort(rows.begin(), rows.end(),
                     [](const Row& x, const Row& y) { return x.spin2 < y.spin2; });
    for (size_t k = 1; k < rows.size(); ++k) {
        const Row& a = rows[k - 1];
        const Row& b = rows[k];
        if (a.spin2 == b.spin2)
            throw PSIEXCEPTION("REKS catalog: the 2S=" + std::to_string(b.spin2) +
                               " manifold is in both " + paths[order[a.file]].first.string() +
                               " and " + paths[order[b.file]].first.string() + ".");
        const bool next_in_file = a.file == b.file && b.file_sector == a.file_sector + 1;
        const bool next_file    = b.file == a.file + 1 && b.file_sector == 0;
        if (!next_in_file && !next_file)
            throw PSIEXCEPTION("REKS catalog: the spin manifolds of " +
                               paths[order[a.file]].first.string() + " and " +
                               paths[order[b.file]].first.string() +
                               " interleave; each blob must hold a run of consecutive 2S.");
    }

    // Numbering, file after file: K and delta ids by offset, L by first-seen union.
    std::unordered_map<uint32_t, int> l_of_key;
    int k_total = 0;
    int d_total = 0;
    for (int f = 0; f < n_files; ++f) {
        CatalogFile& cf = shelf.files[f];
        const std::vector<Microstate>& ms = file_ms[order[f]];
        const int n_act = static_cast<int>(cf.hdr.n_active_orbitals);
        cf.k_shift     = k_total;
        cf.delta_shift = d_total;
        cf.l_global.resize(ms.size());
        for (size_t L = 0; L < ms.size(); ++L) {
            const auto ins = l_of_key.emplace(microstate_key(ms[L], n_act),
                                              static_cast<int>(shelf.microstates.size()));
            if (ins.second) shelf.microstates.push_back(ms[L]);
            cf.l_global[L] = ins.first->second;
        }
        k_total += static_cast<int>(cf.hdr.n_configs);
        d_total += static_cast<int>(cf.hdr.n_delta_specs);

        // Pools whose rows name a microstate, owned with it in shelf L. Weight cells
        // are the cell palette's q == -1 entries, p their microstate.
        const CellPq* pq = reinterpret_cast<const CellPq*>(cf.base + cf.hdr.off_cell_pq);
        if (cf.hdr.n_cell_pq > 0) cf.cell_pq.assign(pq, pq + cf.hdr.n_cell_pq);
        for (CellPq& c : cf.cell_pq)
            if (c.q == -1) c.p = cf.l_global[c.p];
        const auto* fi = reinterpret_cast<const coupling::FockIdxPart*>(
            cf.base + cf.hdr.off_fock_idxpart);
        if (cf.hdr.n_fock_idxpart > 0) cf.fock_idxpart.assign(fi, fi + cf.hdr.n_fock_idxpart);
        for (coupling::FockIdxPart& x : cf.fock_idxpart) x.L = cf.l_global[x.L];
        const auto* ds = reinterpret_cast<const DeltaSpec*>(cf.base + cf.hdr.off_delta_specs);
        for (int64_t d = 0; d < cf.hdr.n_delta_specs; ++d) {
            DeltaSpec s = ds[d];
            if (s.L1 >= 0) s.L1 = cf.l_global[s.L1];
            if (s.L2 >= 0) s.L2 = cf.l_global[s.L2];
            shelf.delta_specs.push_back(s);
        }

        shelf.has_rdm = shelf.has_rdm || cf.hdr.n_rdm_elems > 0;
    }

    // Diagonal records per shelf K; microstate writes are copied out in shelf L,
    // into a pool sized up front so the views stay put.
    size_t n_mw = 0;
    for (const auto& dd : file_dd)
        for (const DiagDeps& d : dd) n_mw += d.mw_n;
    shelf.mw_pool.reserve(n_mw);
    shelf.diag.reserve(static_cast<size_t>(k_total));
    for (int f = 0; f < n_files; ++f) {
        const CatalogFile& cf = shelf.files[f];
        for (const DiagDeps& d : file_dd[order[f]]) {
            DiagDepsView v = resolve_diag_deps(cf.base, d);
            const size_t at = shelf.mw_pool.size();
            for (int L : v.microstate_writes) shelf.mw_pool.push_back(cf.l_global[L]);
            v.microstate_writes = IntSpan{ v.microstate_writes.n ? shelf.mw_pool.data() + at
                                                                 : nullptr,
                                           v.microstate_writes.n };
            shelf.diag.push_back(v);
        }
    }

    // The sector blocks tile the shelf's configs in sector order.
    int k_next = 0;
    for (const Row& r : rows) {
        const CatalogFile& cf = shelf.files[r.file];
        const SectorEntry& se = file_sectors(cf)[r.file_sector];
        if (cf.k_shift + se.config_start != k_next)
            throw PSIEXCEPTION("REKS catalog blob " + paths[order[r.file]].first.string() +
                               ": the 2S=" + std::to_string(se.spin2) +
                               " config block does not follow the one before it.");
        k_next += se.config_count;
        SectorInfo si;
        si.spin2         = se.spin2;
        si.config_start  = cf.k_shift + se.config_start;
        si.config_count  = se.config_count;
        si.n_generations = se.n_generations;
        const int* sa = reinterpret_cast<const int*>(cf.base + se.off_sa_default);
        for (int i = 0; i < se.n_sa_default; ++i) si.sa_default.push_back(sa[i] + cf.k_shift);
        si.file        = r.file;
        si.file_sector = r.file_sector;
        shelf.sectors.push_back(std::move(si));
    }
    if (k_next != k_total)
        throw PSIEXCEPTION("REKS catalog: the spin manifolds cover " + std::to_string(k_next) +
                           " of the " + std::to_string(k_total) + " configs the blobs hold.");
    return shelf;
}

// Cached shelf of (N,M), immutable and pointer-stable for the process lifetime.
const CatalogShelf& load_shelf(int N, int M) {
    static std::mutex mtx;
    static std::map<std::pair<int, int>, CatalogShelf> registry;

    std::lock_guard<std::mutex> lock(mtx);
    const auto key = std::make_pair(N, M);
    auto it = registry.find(key);
    if (it == registry.end()) it = registry.emplace(key, build_shelf(N, M)).first;
    return it->second;
}

// Unbound view of the shelf of (N,M).
Catalog shelf_catalog(int N, int M) {
    const CatalogShelf& shelf = load_shelf(N, M);
    Catalog cat;
    cat.shelf        = &shelf;
    cat.microstates_ = shelf.microstates.empty() ? nullptr : shelf.microstates.data();
    return cat;
}

}  // namespace

Catalog Catalog::bound_to(int b) const {
    const SectorInfo& si = shelf->sectors[b];
    Catalog c = *this;
    c.active_sector = b;
    c.file_sector   = si.file_sector;
    c.file          = &shelf->files[si.file];
    c.base          = c.file->base;
    c.data          = &c.file->hdr;
    return c;
}

int Catalog::n_sectors_available() const { return static_cast<int>(shelf->sectors.size()); }
const SectorInfo& Catalog::sector_entry(int b) const { return shelf->sectors[b]; }

int Catalog::sector_of_config(int K) const {
    for (int b = 0; b < n_sectors_available(); ++b) {
        const SectorInfo& si = shelf->sectors[b];
        if (K >= si.config_start && K < si.config_start + si.config_count) return b;
    }
    return -1;
}

const std::vector<int>& Catalog::sector_sa_default(int b) const {
    return shelf->sectors[b].sa_default;
}

int Catalog::n_microstates() const { return static_cast<int>(shelf->microstates.size()); }
int Catalog::n_configs() const { return static_cast<int>(shelf->diag.size()); }
int Catalog::n_delta_specs() const { return static_cast<int>(shelf->delta_specs.size()); }
const DeltaSpec& Catalog::delta_spec(int d) const { return shelf->delta_specs[d]; }
bool Catalog::has_rdm() const { return shelf->has_rdm; }
DiagDepsView Catalog::diag_deps_view(int K) const { return shelf->diag[K]; }

Catalog::WeightBlock Catalog::weight_block(int K) const {
    const int b = sector_of_config(K);
    if (b < 0) return {};
    const Catalog c = bound_to(b);
    const WeightRow* rows = c.weight_elems(c.base);
    const int k = K - c.file->k_shift;
    if (rows == nullptr || k >= c.data->n_weight_elems) return {};
    return {c.block_pools(c.base), rows + k};
}

const char* Catalog::si_config_name(int K) const {
    const int b = sector_of_config(K);
    if (b < 0) return nullptr;
    const CatalogFile& f = shelf->files[shelf->sectors[b].file];
    const uint32_t* idx = reinterpret_cast<const uint32_t*>(f.base + f.hdr.off_si_config_names);
    return reinterpret_cast<const char*>(f.base + f.hdr.string_table_off + idx[K - f.k_shift]);
}

const char* Catalog::si_config_def(int K) const {
    const int b = sector_of_config(K);
    if (b < 0) return nullptr;
    const CatalogFile& f = shelf->files[shelf->sectors[b].file];
    if (!f.hdr.has_si_config_defs) return nullptr;
    const uint32_t* idx = reinterpret_cast<const uint32_t*>(f.base + f.hdr.off_si_config_defs);
    return reinterpret_cast<const char*>(f.base + f.hdr.string_table_off + idx[K - f.k_shift]);
}

Catalog make_catalog(int N, int M) { return shelf_catalog(N, M).bound_to(0); }

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

// Global config range of the bound sector's block, [begin, end).
static std::pair<int, int> name_scan_range(const Catalog& s) {
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

// Geminal-count markers of config type names; "PPS" = perfectly paired remainder, a run of "T" = triplet
// background, "SPS" = spin-projection tag.
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

// Psi_1 (open-shell singlet) and Psi_2 (doubly excited singlet) counts of a config type name and the
// SPS tag; sps forces n_psi2 = 0.
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

// Pairing scheme of a config = scheme of any geminal slot it reads; -1 for spin-projected configs
// (no slot), 0 in a one-scheme catalog.
static int config_scheme_of(const Catalog& s, int K) {
    const GeminalTemplate* gems = s.geminal_templates(s.base);
    if (gems == nullptr) return 0;
    const DiagDepsView dv = s.diag_deps_view(K);
    for (int gen = 0; gen < kMaxGen; ++gen)
        for (int slot : dv.geminal_reads[gen])
            if (slot >= 0 && slot < s.data->n_geminals) return gems[slot].scheme;
    return catalog_n_schemes(s) == 1 ? 0 : -1;
}

// Name of the bound sector's 2S manifold.
static std::string scan_range_label(const Catalog& s) {
    return "the 2S=" + std::to_string(s.sector_entry(s.active_sector).spin2) + " manifold";
}

// Configs a scheme token keeps; both modes take the spin-projected base whole, once.
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

// Union over the listed schemes N: single = configs of N; mixed = single-excitation configs of N and
// every open-shell-singlet config; both add the spin-projected base.
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

// Value of a non-negative integer of at most 9 decimal digits, false for any other string.
static bool parse_index(const std::string& digits, int& value) {
    if (digits.empty() || digits.size() > 9 || digits.find_first_not_of("0123456789") != std::string::npos)
        return false;
    value = std::stoi(digits);
    return true;
}

// Cycle type of a scheme class, false otherwise: T = (2), D<k> = k disjoint transpositions,
// C<l>[X<l>...] = cycles of those lengths (k, l >= 2), ascending.
static bool scheme_class_of(const std::string& arg, std::vector<int>& cycle) {
    auto count = [](const std::string& digits, int& value) { return parse_index(digits, value) && value >= 2; };
    cycle.clear();
    if (arg == "T") {
        cycle = {2};
        return true;
    }
    if (arg.size() < 2) return false;
    const std::string body = arg.substr(1);
    int k = 0;
    if (arg[0] == 'D') {
        if (!count(body, k)) return false;
        cycle.assign(k, 2);
        return true;
    }
    if (arg[0] != 'C') return false;
    for (size_t pos = 0; pos <= body.size();) {
        const size_t stop = std::min(body.find('X', pos), body.size());
        if (!count(body.substr(pos, stop - pos), k)) return false;
        cycle.push_back(k);
        pos = stop + 1;
    }
    std::sort(cycle.begin(), cycle.end());
    return true;
}

// Cycle type of each scheme s vs scheme 0 from the geminal templates (b = lower orbital, a_s(b) its
// partner): ascending cycle lengths >= 2 of tau_s : a_0(b) -> a_s(b).
static std::map<int, std::vector<int>> scheme_cycle_types(const Catalog& s,
                                                          const std::string& token) {
    std::map<int, std::map<int, int>> partners;
    if (const GeminalTemplate* gems = s.geminal_templates(s.base))
        for (int64_t g = 0; g < s.data->n_geminals; ++g) {
            const auto [lo, hi] = std::minmax(gems[g].orbitals[0], gems[g].orbitals[1]);
            partners[gems[g].scheme][lo] = hi;
        }
    const auto ref = partners.find(0);
    if (ref == partners.end())
        throw PSIEXCEPTION("REKS config token '" + token + "': the catalog carries no geminal "
                           "of scheme 0, the reference of the scheme classes.");
    std::map<int, std::vector<int>> types;
    for (const auto& [scheme, own] : partners) {
        std::map<int, int> tau;
        for (const auto& [b, a0] : ref->second) {
            const auto it = own.find(b);
            if (it == own.end() || own.size() != ref->second.size())
                throw PSIEXCEPTION("REKS config token '" + token + "': scheme " +
                                   std::to_string(scheme) +
                                   " pairs other bonding orbitals than scheme 0.");
            tau[a0] = it->second;
        }
        std::vector<int> lengths;
        std::set<int> seen;
        for (const auto& [start, unused] : tau) {
            int length = 0;
            for (int a = start; seen.insert(a).second; ++length) {
                const auto next = tau.find(a);
                if (next == tau.end())
                    throw PSIEXCEPTION("REKS config token '" + token + "': scheme " +
                                       std::to_string(scheme) +
                                       " pairs other antibonding orbitals than scheme 0.");
                a = next->second;
            }
            if (length >= 2) lengths.push_back(length);
        }
        std::sort(lengths.begin(), lengths.end());
        types[scheme] = std::move(lengths);
    }
    return types;
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

    // Scheme token = base, then indices and classes behind ':', '-', '_' or ',' (folded on a copy, type
    // names carry '_'); a class adds its schemes in catalog order, no argument means scheme 0.
    std::string norm = up;
    for (char& c : norm)
        if (c == '_' || c == ':' || c == ',') c = '-';
    for (const SchemeToken& st : kSchemeTokens) {
        const std::string base(st.base);
        if (norm.rfind(base, 0) != 0) continue;
        std::vector<int> schemes;
        std::map<int, std::vector<int>> types;
        for (size_t pos = base.size(); pos < norm.size();) {
            const size_t stop   = norm.find('-', pos + 1);
            const size_t argend = (stop == std::string::npos) ? norm.size() : stop;
            const std::string arg = norm.substr(pos + 1, argend - pos - 1);
            std::vector<int> picked;
            std::vector<int> cycle;
            int index = 0;
            if (norm[pos] == '-' && parse_index(arg, index)) {
                picked.push_back(index);
            } else if (norm[pos] == '-' && scheme_class_of(arg, cycle)) {
                if (types.empty()) types = scheme_cycle_types(s, token);
                for (const auto& [scheme, type] : types)
                    if (type == cycle) picked.push_back(scheme);
                if (picked.empty())
                    throw PSIEXCEPTION("REKS config token '" + token + "': scheme class '" + arg +
                                       "' has no scheme among the " +
                                       std::to_string(types.size()) + " of this catalog.");
            } else {
                throw PSIEXCEPTION("REKS config token '" + token + "': '" + base +
                                   "' takes scheme indices (non-negative integers) and scheme "
                                   "classes (T, D<k>, C<l>[X<l>...]) behind ':', '-', '_' or ',', "
                                   "and nothing else.");
            }
            for (int scheme : picked) {
                if (std::find(schemes.begin(), schemes.end(), scheme) != schemes.end())
                    throw PSIEXCEPTION("REKS config token '" + token + "': scheme " +
                                       std::to_string(scheme) + " is reached twice.");
                schemes.push_back(scheme);
            }
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

bool has_diabatic_rdm(const Catalog& cat) {
    return cat.shelf && cat.has_rdm();
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
                   const int64_t* fon_win,
                   const FonSlot* slot_pool,
                   const int* slot_idx,
                   const Coeff* coeff_pool,
                   const ValRow* e_val_pool,
                   const FockIdxPart* fock_idxpart,
                   const ValRow* fock_valpart,
                   const int64_t* fock_pack,
                   const EriIdxPart* eri_idxpart,
                   const ValRow* eri_valpart,
                   const int64_t* eri_pack,
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
        const int64_t L = win.e_pack[2 * t];
        const ValRow& e = e_val_pool[win.e_pack[2 * t + 1]];
        const FonWindow fw = fon_window(fon_win, e.win);
        v += coeff_pool[e.coeff].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, fw.off, fw.n, fons,
                                     &f_interp, fon_cache)
               * E_L[L];
    }
    for (int t = 0; t < win.fock_idx.n; ++t) {
        const FockTerm f = gather_fock(fock_idxpart, fock_valpart, fock_pack, coeff_pool,
                                       fon_win, win.fock_idx[t]);
        // Active block of F_MO_a - F_MO_b: row = f.p, column = f.q, both in active order.
        const double fval = fock_aa[static_cast<size_t>(fock_aa_row[f.L]) * n_active * n_active
                                    + static_cast<size_t>(f.p) * n_active + f.q];
        v += f.coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, f.fon_off, f.fon_n, fons,
                                     &f_interp, fon_cache)
               * fval;
    }
    for (int t = 0; t < win.lagr_idx.n; ++t) {
        const LagrTerm& w = lagr_row_pool[win.lagr_idx[t]];
        const FonWindow fw = fon_window(fon_win, w.win);
        v += coeff_pool[w.coeff].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, fw.off, fw.n, fons,
                                     &f_interp, fon_cache)
               * lagrangians[w.w];
    }
    for (int t = 0; t < win.eri_idx.n; ++t) {
        const EriTerm g = gather_eri(eri_idxpart, eri_valpart, eri_pack, coeff_pool,
                                     fon_win, win.eri_idx[t]);
        v += g.coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, g.fon_off, g.fon_n, fons,
                                     &f_interp, fon_cache)
               * active_eri.at(g.p, g.q, g.r, g.s);
    }
    return v;
}

// S_{KL} for one SI overlap pair: sum over pure-FON terms, each term c_t * prod(FONs_t).
double eval_pair_S(const IdSpan& s_idx, const FonFactor* fon_pool,
                   const int* fon_idx, const int64_t* fon_win, const FonSlot* slot_pool,
                   const int* slot_idx, const Coeff* coeff_pool, const ValRow* s_row_pool,
                   const FONSnapshot& fons, const double* fon_cache) {
    double v = 0.0;
    for (int t = 0; t < s_idx.n; ++t) {
        const ValRow& s = s_row_pool[s_idx[t]];
        const FonWindow fw = fon_window(fon_win, s.win);
        v += coeff_pool[s.coeff].coeff
               * studio::fon_product(fon_pool, slot_pool, fon_idx, slot_idx, fw.off, fw.n, fons,
                                     &f_interp, fon_cache);
    }
    return v;
}

}  // namespace coupling

}  // namespace studio
}  // namespace reks
}  // namespace psi
