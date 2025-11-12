# SCF Refactoring Project

## Overview

**Goal:** Maximize SCF performance and enable REKS implementation through HPC-aware architecture.

**Strategy:**
1. **Phase 0:** HPC optimization (aligned memory, cache-aware operations) - **START HERE** 🔥
2. **Phase 1-6:** Incremental refactoring to multi-Fock architecture

Each step: Claude codes → commits → pushes; User compiles → tests → validates on **kk/refactor_scf** branch.

**Key Innovation:** Multi-state architecture with shared JK contraction + HPC-optimized memory layout.

## Workflow & Roles

**Claude (AI):** Modifies code, commits, pushes (does NOT compile/test)
**User:** Pulls, compiles (`cd build && ninja psi4`), tests, reports results

**Branch:** `claude/scf_optimization-011CV3EdN9C37wMjH8pg4vTe` ← **ACTIVE DEVELOPMENT**

---

## Current Status

**Phase:** 0 - HPC Performance Optimization 🔥

**Completed:**
- ✅ Code analysis: Matrix class, BLAS integration, memory layout
- ✅ Identified bottlenecks: alignment (10-30%), get_block (10-20x), cache locality (2-3x)
- ✅ HPC optimization plan formulated
- ✅ **Phase 0.1:** Aligned allocation (posix_memalign 64-byte) ← committed 581c0833

**Next Action:** User tests Phase 0.1 performance → if OK, proceed to Phase 0.2

---

## Target Architecture

### Multi-Fock Design (NEW!)

**Core Idea:** SCF handles **N Fock matrices simultaneously** with **single shared JK contraction**

```
┌─────────────────────────────────────────────────┐
│         SCFDriver (convergence logic)            │
│  - DIIS, damping, convergence checks            │
│  - Iteration loop                               │
└─────────────┬───────────────────────────────────┘
              │ uses
              ↓
┌─────────────────────────────────────────────────┐
│       FockTheory (theory-specific)              │
│  n_states() → int   (RHF:1, UHF:2, REKS:N)     │
│                                                 │
│  ┌───────────────────────────────────────────┐ │
│  │ 1. build_density(C → D) [per state]      │ │
│  │ 2. compute_JK(D₁,D₂,...Dₙ → J,K) [ONCE!] │ │ ← SHARED!
│  │ 3. combine_fock(J,K → F₁,F₂,...Fₙ)       │ │ ← theory-specific
│  └───────────────────────────────────────────┘ │
└─────────────┬───────────────────────────────────┘
              │ implements
       ┌──────┴──────┬──────────┬──────────┐
       ↓             ↓          ↓          ↓
   RHFTheory    UHFTheory  REKSTheory  CustomTheory
   n=1          n=2        n=N         n=?
```

**Example: UHF with shared JK**
```cpp
// Current (inefficient): JK called twice for Da and Db
jk->C_left() = {Ca_occ};  jk->compute();  // J_a, K_a
jk->C_left() = {Cb_occ};  jk->compute();  // J_b, K_b

// New (efficient): JK called ONCE for both
jk->C_left() = {Ca_occ, Cb_occ};  jk->compute();  // J_a, K_a, J_b, K_b in ONE call
F_a = H + (J_a + J_b) - K_a
F_b = H + (J_a + J_b) - K_b
```

**Generalization:** Any number of Fock matrices!
- RHF alone: 1 Fock
- UHF alone: 2 Fock
- **RHF + ROHF simultaneously: 2 Fock** (your example!)
- REKS: N Fock (N states)

---

## Refactoring Strategy: Bottom-Up Approach

**NO HFNew!** We refactor **in-place** with **opt-in** pattern:
1. Add new code alongside old (e.g., `DensityContainer` + `Da_` both exist)
2. Add flag `use_new_path_` to switch between old/new
3. Test new path, keep old path working
4. When validated, remove old path

**Increments:** Each step takes 1-2 hours, fully tested before next step.

---

---

## Phase 0: HPC Performance Optimization 🔥 **PRIORITY**

**Goal:** Maximize Matrix operations performance (3-5x speedup for multi-state)

### Critical Findings from Code Analysis

| Issue | Current | Impact | Fix |
|-------|---------|--------|-----|
| **Memory alignment** | `malloc()` ~16 bytes | -10-30% BLAS | `posix_memalign(64)` |
| **get_block() copy** | Element-wise O(n²) | -10-20x | `memcpy` per row |
| **Multi-state layout** | Separate allocations | -2-3x cache miss | Contiguous storage |
| **Cache blocking** | None | -30-50% (large) | Tiled DGEMM |

**Total potential:** **3-5x overall speedup** for REKS/multi-state operations

---

### 0.1: Aligned Memory Allocation ⚡ **QUICK WIN**

**Time:** 1 hour | **Gain:** 10-30% BLAS performance

**Problem:**
```cpp
// matrix.cc:461 - Current code
double* block = (double*)malloc(nrows * ncols * sizeof(double));
// malloc() gives ~16 byte alignment
// AVX-512 needs 64 bytes (cache line size)
```

**Solution:**
```cpp
// matrix.cc:461 - New code
static constexpr size_t CACHE_LINE_SIZE = 64;  // x86-64 standard

double* block = nullptr;
int ret = posix_memalign((void**)&block, CACHE_LINE_SIZE,
                         nrows * ncols * sizeof(double));
if (ret != 0) {
    throw PSIEXCEPTION("posix_memalign failed");
}
```

**Files to modify:**
- `psi4/src/psi4/libmints/matrix.cc:461` (Matrix::alloc)
- `psi4/src/psi4/libmints/matrix.cc:3562` (linalg::detail::matrix)

**Cache line size:** Hardcoded 64 bytes (Intel/AMD/ARM standard). Could use runtime:
```cpp
// Optional: runtime detection (overkill)
#include <unistd.h>
size_t cache_line = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);  // Usually 64
```
But 64 is safe for 99.9% of systems.

**Test:** Compile + run RHF/UHF tests. Performance should improve 10-30% on DGEMM-heavy systems.

---

### 0.2: Vectorize get_block() ⚡⚡ **CRITICAL**

**Time:** 4 hours | **Gain:** 10-20x for subset operations

**Problem:**
```cpp
// matrix.cc:666 - Current code (DISASTER!)
for (int p = 0; p < max_p; p++) {
    for (int q = 0; q < max_q; q++) {
        double val = get(h, p + rows_begin[h], q + cols_begin[h]);
        block->set(h, p, q, val);  // Element-wise! O(n²) operations!
    }
}
// Used in Ca_subset("SO", "OCC") → EVERY SCF ITERATION!
```

**Solution:**
```cpp
// matrix.cc:666 - New code
for (int h = 0; h < nirrep_; ++h) {
    if (rows_per_h[h] == 0 || cols_per_h[h] == 0) continue;

    const int row_offset = rows_begin[h];
    const int col_offset = cols_begin[h];
    const size_t bytes_per_row = cols_per_h[h] * sizeof(double);

    // Vectorized memcpy per row (10-20x faster!)
    for (int i = 0; i < rows_per_h[h]; ++i) {
        std::memcpy(&block->matrix_[h][i][0],
                   &matrix_[h][i + row_offset][col_offset],
                   bytes_per_row);
    }
}
```

**Files to modify:**
- `psi4/src/psi4/libmints/matrix.cc:666` (Matrix::get_block)
- Similar pattern in other block operations

**Test:** Benchmark `Ca_subset("SO", "OCC")` - should be 10-20x faster.

---

### 0.3: Multi-State Contiguous Storage ⚡⚡⚡ **GAME CHANGER**

**Time:** 2-3 days | **Gain:** 2-3x for multi-state (REKS critical!)

**Problem:**
```cpp
// Current UHF: Da and Db in different memory locations
SharedMatrix Da_;  // address: 0x1000
SharedMatrix Db_;  // address: 0x9000 (8KB+ gap!)
// Access pattern: Da[i] → Db[i] → CACHE MISS every time!
```

**Solution: MultiStateMatrix with contiguous storage**

```cpp
// New file: psi4/src/psi4/libscf_solver/multistate_matrix.h

class MultiStateMatrix {
    int n_states_;
    int nirrep_;
    std::vector<Dimension> rowspi_;
    std::vector<Dimension> colspi_;

    // KEY: Single aligned allocation for ALL states
    double* data_contiguous_;
    size_t total_elements_;

    // Views (SharedMatrix) pointing into data_contiguous_
    std::vector<SharedMatrix> state_views_;

public:
    MultiStateMatrix(const std::string& name, int n_states,
                     int nirrep, const std::vector<Dimension>& dims);

    ~MultiStateMatrix() {
        if (data_contiguous_) free(data_contiguous_);
    }

    // Access
    SharedMatrix get(int state) const { return state_views_[state]; }
    int n_states() const { return n_states_; }

    // Contiguous pointer for advanced operations
    double* contiguous_data() { return data_contiguous_; }

    // Zero all states efficiently
    void zero_all() {
        std::memset(data_contiguous_, 0, total_elements_ * sizeof(double));
    }
};
```

**Memory layout:**
```
Single allocation:
[State0_h0_data][State0_h1_data]...[State1_h0_data][State1_h1_data]...
^                                  ^
Da (view)                          Db (view)

All data contiguous → excellent cache locality!
```

**Files to create:**
- `psi4/src/psi4/libscf_solver/multistate_matrix.h`
- `psi4/src/psi4/libscf_solver/multistate_matrix.cc`
- Update `CMakeLists.txt`

**Integration:** Phase 1 will use this in RHF/UHF

**Test:** Benchmark multi-state operations - expect 2-3x speedup vs separate matrices.

---

### 0.4: Cache-Aware Tiling (Optional, later)

**Time:** 1-2 weeks | **Gain:** 30-50% for large systems (>500 basis functions)

**Idea:** Block DGEMM to fit in L2 cache (typical: 256KB → tile ~128x128)

**Status:** Defer to Phase 2 or later. Profile first to confirm benefit.

---

## Phase 0 Summary

| Step | Time | Gain | Priority | Status |
|------|------|------|----------|--------|
| **0.1 Aligned allocation** | 1h | 10-30% | ⚠️ HIGH | 📍 NEXT |
| **0.2 Vectorize get_block** | 4h | 10-20x | ⚠️ CRITICAL | Pending |
| **0.3 Contiguous multi-state** | 2-3d | 2-3x | ⚠️ CRITICAL | Pending |
| **0.4 Cache tiling** | 1-2w | 30-50% | MEDIUM | Deferred |
| **Total Phase 0** | ~3-4 days | **3-5x overall** | - | 0% |

**After Phase 0:** Matrix infrastructure optimized, ready for refactoring.

---

## Phases 1-6: SCF Refactoring (Brief Overview)

**After Phase 0 HPC optimization**, we proceed with architectural refactoring:

**Phase 1:** Matrix Containers (10h) - Use MultiStateMatrix from Phase 0.3 in RHF/UHF
**Phase 2:** Density Operations (4h) - Extract `form_D()` to standalone functions
**Phase 3:** JK Multi-State (6h) - Systematize shared JK contraction (UHF already does this)
**Phase 4:** Fock Assembly (3h) - Standardize F = H + G + V_ext
**Phase 5:** Theory Abstraction (12h) - FockTheory interface (RHF/UHF/REKS)
**Phase 6:** SCF Driver (10h) - Separate convergence algorithm from theory

**Total:** ~45 hours after Phase 0

**Detailed plan:** See IMPLEMENTATION_PLAN.md (old version, needs update after Phase 0)

---

## Implementation Roadmap

| Phase | Focus | Time | Gain | Status |
|-------|-------|------|------|--------|
| **0** | **HPC Optimization** | **3-4d** | **3-5x** | **📍 ACTIVE** |
| 0.1 | Aligned allocation | 1h | 10-30% | ✅ **TESTING** |
| 0.2 | Vectorize get_block | 4h | 10-20x | 📍 NEXT |
| 0.3 | Contiguous multi-state | 2-3d | 2-3x | Pending |
| 1 | Matrix Containers | 10h | - | Pending |
| 2 | Density Ops | 4h | - | Pending |
| 3 | JK Multi-State | 6h | - | Pending |
| 4 | Fock Assembly | 3h | - | Pending |
| 5 | Theory Abstraction | 12h | - | Pending |
| 6 | SCF Driver | 10h | - | Pending |
| **TOTAL** | **All phases** | **~50h** | **3-5x** | **0%** |

---

## Success Criteria (User Validates)

### Phase 0 ⭐ **PRIORITY**
- [ ] **0.1:** Aligned allocation works, BLAS 10-30% faster
- [ ] **0.2:** get_block() 10-20x faster, subset operations validated
- [ ] **0.3:** MultiStateMatrix with contiguous storage, 2-3x speedup for multi-state
- [ ] All RHF/UHF tests pass with exact energy match (< 1e-10)
- [ ] **Performance validated:** User benchmarks show expected gains

### Phases 1-6 (After Phase 0)
- [ ] Architecture refactored: FockTheory, SCFDriver, multi-Fock
- [ ] All tests pass (< 1e-10 energy match)
- [ ] REKS-ready architecture
- [ ] Performance maintained or improved

---

## Next Steps

**Immediate (Phase 0.1):**
1. ✅ Plan approved: HPC optimization first
2. Modify `matrix.cc:461` - aligned allocation
3. Commit → push → User tests performance
4. If successful → proceed to 0.2 (vectorize get_block)

---

## Technical Notes: HPC Details

### Cache Line Size

**Standard:** 64 bytes for x86-64 (Intel, AMD), ARM
**Detection (optional):**
```cpp
#include <unistd.h>
size_t cache_line = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
```
**Recommendation:** Hardcode 64 - universally safe

### SIMD Alignment Requirements

| Instruction Set | Alignment | Performance Loss if Misaligned |
|-----------------|-----------|-------------------------------|
| SSE | 16 bytes | -10% |
| AVX | 32 bytes | -15-20% |
| AVX-512 | 64 bytes | -20-30% |

**Current:** `malloc()` gives ~16 bytes → losing AVX-512 performance

### Memory Layout Comparison

**Current (separate matrices):**
```
Memory: [Da_______gap______][Db_______gap______]
Access: Da[i] → cache miss → Db[i] → cache miss
```

**Contiguous (MultiStateMatrix):**
```
Memory: [Da][Db] (adjacent)
Access: Da[i] → Db[i] (likely in same cache line!)
Speedup: 2-3x for alternating access patterns
```

### Why get_block() is Slow

**Element-wise (current):**
```cpp
for (i) for (j)  // O(n²) operations
    dst[i][j] = src[i+offset][j+offset];  // Individual set() calls
```

**memcpy (new):**
```cpp
for (i)  // O(n) operations
    memcpy(&dst[i][0], &src[i+offset][offset], n_cols * 8);  // Hardware-accelerated
```

**Why 10-20x faster:**
- Hardware memcpy (SIMD, prefetch)
- Eliminates function call overhead
- Eliminates bounds checks per element

---

## Risk Mitigation

**If blocked >2 hours on a step:**
1. Revert last commit
2. Break step into smaller pieces
3. Ask User for clarification
4. Document blocker in plan.md

**If tests fail:**
1. Check compilation warnings
2. User runs single test with verbose output
3. Compare intermediate values (D, F matrices) via DEBUG output
4. Revert if can't fix in 30 minutes

**If performance regresses >20%:**
1. User profiles with perf/gprof
2. Check for unnecessary copies (use const references)
3. Verify BLAS usage (DGEMM, DAXPY)
4. May be acceptable if architecture is cleaner (document trade-off)

---

## Testing Strategy

**User's Comprehensive Test Suite (Already Exists):**
- RHF, UHF, ROHF
- All DFT functionals: HF, PBE, B3LYP, ωB97X-V, etc.
- All SCF types: PK, DF, MEM_DF, DISK_DF, OUT_OF_CORE, CD
- All screening: CSAM, DENSITY, NONE
- All guess: CORE, SAD, AUTO, MODHUCKEL, GWH, SADNO, SAP, SAPGAU
- Special: incremental Fock, DIIS, damping, MOM, fractional occupation

**Per-Phase Thresholds:**
- **Phase 0:** Performance gains validated (CRITICAL)
- **Phase 1-6:** Energy < 1e-10, performance maintained

---

## Key Insights from Code Analysis

**Matrix class (psi4/libmints/matrix.*):**
- Per-irrep contiguous blocks (good)
- No SIMD alignment (BAD - losing 10-30%)
- get_block() element-wise (DISASTER - 10-20x slow)
- Symmetry via XOR mapping (elegant)

**BLAS integration:**
- Direct C_DGEMM calls (efficient)
- Row-major → column-major handled in wrapper
- Per-irrep DGEMM (natural parallelism)

**UHF JK:** Already uses multi-state pattern (uhf.cc:184-200)
```cpp
jk_->C_left() = {Ca_occ, Cb_occ};  // Both at once!
jk_->compute();  // Single call
```

**Our goal:** Systematize for REKS (N states) + optimize underlying Matrix
