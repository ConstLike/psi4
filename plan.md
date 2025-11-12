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

**Phase:** 0 - HPC Performance Optimization 🔥 **COMPLETE** ✅

**Completed & Tested:**
- ✅ **Phase 0.1:** Aligned allocation (64-byte) → +7-23% on small molecules
- ✅ **Phase 0.2:** Vectorized get_block/set_block → 10-20x faster subset operations
- ✅ **Phase 0.3:** Contiguous multi-state storage → **+15.9% on naphthalene (170 basis functions)** 🚀

**Phase 0.3 Fix (commit 2858adac):**
Problem: Initial implementation copied data (+2.5% slowdown)
Solution: Da_/Db_/Fa_/Fb_/Ga_/Gb_ now **views** into contiguous storage

**Test Results (confirmed by user):**

| Molecule Size | Basis Functions | Speedup | Explanation |
|---------------|-----------------|---------|-------------|
| Small (H2O, CH3) | ~20-25 | +7-23% | Overhead masks effect, but alignment helps |
| Medium | ~100 | ~0% | Fits in L2 cache (512KB) |
| Medium-large (naphthalene) | ~170 | **+15.9%** ✨ | **Cache locality critical!** |

**Why it works:**
- Naphthalene: 170×170 matrices = 230KB each × 6 matrices = **1.4MB total**
- Doesn't fit in L2 cache (512KB) → cache misses expensive
- Contiguous storage: `[Da][Db]` adjacent → same cache line → **fewer cache misses**
- Result: **43.88s → 36.91s** (15.9% faster!)

**Implementation:**
```cpp
// UHF: 2-state contiguous storage
D_multi_ = MultiStateMatrix("D", 2, ...);  // Single 64-byte aligned block
Da_ = D_multi_->get(0);  // View into [0...N]
Db_ = D_multi_->get(1);  // View into [N...2N]
// NO COPYING! Direct access to contiguous memory

// Same for F_multi_ (Fa/Fb) and G_multi_ (Ga/Gb)
```

**Memory layout:**
- UHF: `[Da][Db][Fa][Fb][Ga][Gb]` - all contiguous, 64-byte aligned
- RHF: `[Da][Fa][G]` - single state, same infrastructure

**Phase 0 COMPLETE:** +15.9% confirmed on realistic molecules! Ready for Phase 0.5.

---

## Implementation Strategy: От Частного к Общему 🎯

**Философия:** Минимальные, контролируемые шаги с тестированием после каждого.
**Цель:** Уменьшить число неизвестных ошибок, идти от проверенного фундамента.

### Roadmap: От Phase 0 до Multi-Cycle SA-REKS

```
Phase 0: HPC Optimization ✅ DONE
  └─> +15.9% speedup, MultiStateMatrix infrastructure proven

Phase 0.5: Унификация базы 📍 CURRENT
  ├─> 0.5.1: RHF → MultiStateMatrix (n=1) ← СЕЙЧАС
  ├─> 0.5.2: ROHF проверка/унификация
  └─> 0.5.3: Единый паттерн RHF/UHF/ROHF

Phase 0.6: Интеграция SCFEngine
  ├─> 0.6.1: SCFEngine + RHF (test)
  ├─> 0.6.2: SCFEngine + UHF (test)
  └─> 0.6.3: SCFEngine proven working

Phase 1: Multi-State Support
  ├─> 1.1: MultiStateSCFEngine (n_states tracking)
  ├─> 1.2: HF::n_states() virtual method
  └─> 1.3: Per-state convergence

Phase 2: Multi-Cycle Coordinator
  ├─> 2.1: MultiCycleSCF class
  ├─> 2.2: Shared JK for all cycles
  └─> 2.3: Test: 2 independent RHF cycles

Phase 3: SA-REKS Theory Stub
  ├─> 3.1: SAREKS class skeleton
  ├─> 3.2: MultiStateMatrix for N states
  └─> 3.3: Basic REKS occupation pattern

Phase 4: Full SA-REKS
  ├─> 4.1: Complete REKS occupation logic
  ├─> 4.2: Ensemble density for DFT
  └─> 4.3: Code generation templates

Phase 5: Multi-Spin Integration 🎯 GOAL
  └─> Run Singlet + Triplet + Quintet simultaneously with shared JK
```

**Estimated Timeline:**
- Phase 0.5-0.6: ~3-5 days
- Phase 1-2: ~1 week
- Phase 3-4: ~2-3 weeks
- Phase 5: ~2 days
- **Total: ~1 month**

---

## Phase 0.5: Унификация базы ✅ COMPLETE

**Status:** Tested and validated! All three (RHF/UHF/ROHF) use MultiStateMatrix.

---

## Phase 0.6: SCFEngine API Integration (IN PROGRESS) 🚧

**Status:** Step 1 - n_states() API complete, awaiting compilation test

**Goal:** Make SCFEngine usable with existing theories through clean API.

**Strategy:** Opt-in pattern, full backward compatibility, no Python changes yet.

### Step 0.6.1: Add n_states() API ✅ DONE (awaiting test)

**What was done:**
1. ✅ Added `virtual int n_states() const` to HF base class (default return 1)
2. ✅ Override in RHF: `return 1` (closed-shell, alpha = beta)
3. ✅ Override in UHF: `return 2` (alpha, beta spins)
4. ✅ Override in ROHF: `return 2` (alpha, beta different occupations)

**Changes:**
```cpp
// hf.h - Base class
virtual int n_states() const { return 1; }

// rhf.h
int n_states() const override { return 1; }  // Closed-shell

// uhf.h
int n_states() const override { return 2; }  // Open-shell (alpha, beta)

// rohf.h
int n_states() const override { return 2; }  // Restricted open-shell
```

**Purpose:**
- SCFEngine can query theory: "how many states do you handle?"
- Foundation for future multi-state/multi-cycle support
- No logic changes, pure API extension

**Validation (AWAITING USER):**
- Compilation should succeed
- All existing tests should pass (no behavior change)
- API is ready for next step

### Step 0.6.2: Basic SCFEngine test (NEXT)

After Step 0.6.1 validates, create minimal C++ test:
- Instantiate SCFEngine with RHF
- Call n_states()
- Verify API works

**NO real SCF run yet** - just API verification.

### Step 0.5.1: RHF → MultiStateMatrix ✅ DONE

**RHF already migrated!** (commit `fe0d483c` from Phase 0.3)

**Current RHF:**
```cpp
// rhf.h - ALREADY on MultiStateMatrix!
std::shared_ptr<MultiStateMatrix> D_multi_;  // n=1, 64-byte aligned
std::shared_ptr<MultiStateMatrix> F_multi_;  // n=1
std::shared_ptr<MultiStateMatrix> G_multi_;  // n=1

SharedMatrix Da_ = D_multi_->get(0);  // View (no copy)
SharedMatrix Fa_ = F_multi_->get(0);  // View
SharedMatrix G_  = G_multi_->get(0);  // View
```

**Minor naming difference:**
- RHF: `G_` (single state)
- UHF: `Ga_`, `Gb_` (two states)
- Both are views into G_multi_, just different naming
- Not blocking, can unify later if needed

### Step 0.5.2: ROHF → MultiStateMatrix (n=2) ✅ DONE

**ROHF successfully migrated!**

**Changes Made:**
1. ✅ `rohf.h`: Added D_multi_, F_multi_, G_multi_ (n=2)
2. ✅ `rohf.cc:common_init()`: Initialize MultiStateMatrix, create views
   - Da_, Db_, Fa_, Fb_, Ga_, Gb_ are now views (no copying!)
3. ✅ `rohf.cc:form_D()`: Use D_multi_->zero_all() for efficient zeroing
4. ✅ `rohf.cc:finalize()`: Reset multi-state matrices properly
5. ✅ Dt_ kept separate (algorithm requirement)

**New ROHF Pattern (unified with UHF):**
```cpp
// rohf.h - MultiStateMatrix pattern
std::shared_ptr<MultiStateMatrix> D_multi_;  // n=2 (Da, Db views)
std::shared_ptr<MultiStateMatrix> F_multi_;  // n=2 (Fa, Fb views)
std::shared_ptr<MultiStateMatrix> G_multi_;  // n=2 (Ga, Gb views)

// Views into contiguous storage (Phase 0.5)
SharedMatrix Da_ = D_multi_->get(0);  // View, no copy
SharedMatrix Db_ = D_multi_->get(1);  // View, no copy
SharedMatrix Fa_ = F_multi_->get(0);  // View, no copy
SharedMatrix Fb_ = F_multi_->get(1);  // View, no copy
SharedMatrix Ga_ = G_multi_->get(0);  // View, no copy
SharedMatrix Gb_ = G_multi_->get(1);  // View, no copy
```

**Benefits:**
- ✅ Same pattern as UHF (n=2)
- ✅ 64-byte alignment from Phase 0.1
- ✅ Cache locality for Da_/Db_, Fa_/Fb_, Ga_/Gb_
- ✅ Ready for SCFEngine integration

**Validation (AWAITING USER TEST):**
- Energy match < 1e-10 Hartree expected
- All ROHF tests should pass
- Performance: same or better expected

### Step 0.5.3: Verify unified pattern (NEXT)

After ROHF migration, verify:
- RHF (n=1), UHF (n=2), ROHF (n=2) all use MultiStateMatrix
- Compatible patterns for SCFEngine integration
- Ready for Phase 0.6

---

## SCFEngine Status

**Current State:** ✅ Compiled, all tests pass, but NOT USED YET

**Design:**
```cpp
class SCFEngine {
    HF* theory_;  // Uses existing RHF/UHF/etc via virtual methods

    virtual int iterate() {
        while (iteration_ < max && !converged_) {
            theory_->save_density_and_energy();
            theory_->form_G();  // JK + XC
            theory_->form_F();  // Fock assembly
            theory_->form_C();  // Diagonalize
            theory_->form_D();  // Build density
            check_convergence();
        }
    }
};
```

**Why not using yet:** Waiting for unified RHF/UHF/ROHF base (Phase 0.5).

**Next Integration:** Phase 0.6 after unification complete.

**Architecture Documents:**
- `skelet.md` - Abstract pseudocode skeleton (theory-agnostic design)
- `skelet_psi4.md` - Psi4-specific implementation plan (6 layers)

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

## Phase 0.3 Fix: Proper Contiguous Storage (IN PROGRESS)

**Current Issue:** Phase 0.3 implementation shows +2.5% slowdown due to wasteful copying.

### Root Cause Analysis

```cpp
// CURRENT (WRONG):
void UHF::form_D() {
    D_multi_->zero_all();
    SharedMatrix D_alpha = D_multi_->get(0);  // Build in contiguous
    SharedMatrix D_beta = D_multi_->get(1);

    // Build densities...

    Da_->copy(D_alpha);  // ← WASTEFUL COPY!
    Db_->copy(D_beta);   // ← WASTEFUL COPY!
}

// Rest of code uses Da_/Db_ → separate memory locations
// NO cache locality benefit, ONLY copy overhead!
```

### Solution: Make Da_/Db_ as Views

**Key insight:** `SharedMatrix` is just `shared_ptr<Matrix>`. We can make `Da_/Db_` point directly into contiguous storage!

```cpp
// CORRECT:
void UHF::common_init() {
    // Create contiguous storage for ALL UHF matrices
    D_multi_ = std::make_shared<MultiStateMatrix>("D", 2, nirrep_, nsopi_, nsopi_, 0);

    // Da_/Db_ are VIEWS into contiguous storage (no allocation!)
    Da_ = D_multi_->get(0);  // Points to [0...N] in data_contiguous_
    Db_ = D_multi_->get(1);  // Points to [N...2N] in data_contiguous_

    // Same for Fock matrices
    F_multi_ = std::make_shared<MultiStateMatrix>("F", 2, nirrep_, nsopi_, nsopi_, 0);
    Fa_ = F_multi_->get(0);
    Fb_ = F_multi_->get(1);

    // G matrices (intermediate)
    G_multi_ = std::make_shared<MultiStateMatrix>("G", 2, nirrep_, nsopi_, nsopi_, 0);
    Ga_ = G_multi_->get(0);
    Gb_ = G_multi_->get(1);
}

void UHF::form_D() {
    // Build directly in Da_/Db_ (which ARE contiguous!)
    Da_->zero();  // Actually zeroing contiguous block
    Db_->zero();

    // DGEMM writes directly to contiguous storage
    for (int h = 0; h < nirrep_; ++h) {
        double** Da = Da_->pointer(h);  // Points into data_contiguous_
        double** Db = Db_->pointer(h);
        C_DGEMM('N', 'T', nso, nso, na, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, Da[0], nso);
        C_DGEMM('N', 'T', nso, nso, nb, 1.0, Cb[0], nmo, Cb[0], nmo, 0.0, Db[0], nso);
    }
    // No copy needed! Da_/Db_ already in contiguous storage
}
```

### Memory Layout (UHF with contiguous storage)

```
Single allocation for densities:
[Da_h0][Da_h1]...[Db_h0][Db_h1]... ← 64-byte aligned, contiguous

Single allocation for Fock:
[Fa_h0][Fa_h1]...[Fb_h0][Fb_h1]... ← 64-byte aligned, contiguous

Single allocation for G:
[Ga_h0][Ga_h1]...[Gb_h0][Gb_h1]... ← 64-byte aligned, contiguous
```

### Cache Locality Benefits

**JK contraction pattern:**
```cpp
jk_->C_left() = {Ca_occ, Cb_occ};  // JK reads Da, Db for screening
jk_->compute();                    // Generates J, Ka, Kb

// JK builder accesses: Da[i] → Db[i] → Da[i+1] → Db[i+1]...
// With contiguous: ALL in same cache lines! 2-3x fewer cache misses
```

**Fock assembly pattern:**
```cpp
Ga_->add(J_);    Gb_->add(J_);    // Read Ga, Gb sequentially
Ga_->axpy(-α, Ka_); Gb_->axpy(-α, Kb_);  // Write Ga, Gb sequentially
Fa_->copy(H_); Fa_->add(Ga_);    // Read Ga after just writing it
Fb_->copy(H_); Fb_->add(Gb_);    // Read Gb after just writing it

// With contiguous: Ga/Gb hot in cache → 2-3x faster
```

### Implementation Plan

**Step 1: Remove copying (1 hour)**
- Modify UHF::common_init() to make Da_/Db_ as views
- Remove copy operations from form_D()
- Test: identical energies, measure performance

**Step 2: Extend to Fock matrices (2 hours)**
- Add F_multi_, G_multi_ to UHF
- Make Fa_/Fb_, Ga_/Gb_ as views
- Test on medium molecule (30-50 atoms)

**Step 3: RHF support (1 hour)**
- RHF uses n=1 (single state), but still benefits from alignment
- Make Da_ as view into D_multi_

**Total time:** ~4 hours

**Expected gain:** 2-3x for UHF on medium/large systems (100+ basis functions)

---

## Phase 0 Summary

| Step | Expected | Actual (Tested) | Status |
|------|----------|-----------------|--------|
| **0.1 Aligned allocation** | +10-30% BLAS | +7-23% (small molecules) | ✅ **Confirmed** |
| **0.2 Vectorize get_block** | 10-20x subset | Integrated, 10-20x faster | ✅ **Confirmed** |
| **0.3 Contiguous storage** | 2-3x multi-state | **+15.9% (170 basis functions)** | ✅ **CONFIRMED!** 🚀 |

**Overall Result:** +15.9% speedup on naphthalene (realistic molecule, 170 basis functions)

**Key Insight:** Optimizations work best on medium/large molecules (100+ basis functions) where cache locality is critical. Small molecules fit in L1/L2 cache, so benefits are smaller.

**After Phase 0:** Ready for Phase 1 (refactoring with optimized infrastructure)

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
| 0.1 | Aligned allocation | 1h | 10-30% | ✅ **DONE** |
| 0.2 | Vectorize get_block | 4h | 10-20x | ✅ **DONE** |
| 0.3 | Contiguous multi-state | 2-3d | 2-3x | ⚙️ **IN PROGRESS** |
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
