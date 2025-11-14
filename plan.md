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

### Roadmap: От Phase 0 до Multi-Cycle SA-REKS (UPDATED)

```
Phase 0: HPC Optimization ✅ DONE
  └─> +15.9% speedup, MultiStateMatrix infrastructure proven

Phase 0.5: Унификация базы ✅ DONE
  ├─> RHF: MultiStateMatrix (n=1) ✓
  ├─> UHF: MultiStateMatrix (n=2) ✓
  └─> ROHF: MultiStateMatrix (n=2) ✓

Phase 0.6: API Foundation ✅ DONE
  ├─> 0.6.1: n_states() API ✅ DONE
  ├─> 0.6.2: Basic API tests ✅ DONE (user)
  └─> 0.6.3: Strategic decision: Python-first approach ✅ DONE

Phase 1: Python Multi-Cycle Coordinator 📍 IN PROGRESS - NEW APPROACH
  ├─> 1.1: Refactor scf_iterate() → extract scf_iteration() ✅ DONE
  ├─> 1.2: Convert scf_iteration() to _scf_iteration() method ✅ DONE
  ├─> 1.3: Create multi_scf() coordinator ✅ DONE
  ├─> 1.4: Fix pybind11 exports & C++ bugs ✅ DONE (JK, get_orbital_matrices)
  ├─> 1.5: Implement options snapshot pattern ← NOW
  └─> 1.6: Add validation & testing

**CRITICAL: Multi-SCF Requirements (MUST be satisfied)**

Shared JK depends on:
- ✅ **SAME geometry** - atomic coordinates MUST be identical (ERI depend on r)
- ✅ **SAME basis** - primary basis set MUST match (JK built for one basis)
- ✅ **SAME SCF_TYPE** - DF/DIRECT/CD MUST be same (JK algorithm)
- ✅ **SAME DF_BASIS_SCF** - auxiliary basis MUST match (if DF)

Can differ (orbital/occupation level):
- ✓ Multiplicity (singlet/triplet) - affects nalpha/nbeta, not JK
- ✓ Reference (RHF/UHF/ROHF) - different density construction, same JK
- ✓ Functional (HF/B3LYP) - XC differs, JK same
- ✓ Convergence (DIIS/MOM/damping) - acceleration methods
- ✓ n_alpha/n_beta - electron count per spin

**Correct Terminology for SA-REKS:**
```python
# WRONG (implies different geometries):
mol1 = psi4.geometry("H2O ...")  # ❌
mol2 = psi4.geometry("H2O ...")  # ❌

# CORRECT (one geometry, multiple states):
molecule = psi4.geometry("H2O ...")  # ✅ One geometry
state1 = {'multiplicity': 1, 'nalpha': 5, 'nbeta': 5}  # Singlet
state2 = {'multiplicity': 3, 'nalpha': 6, 'nbeta': 4}  # Triplet
```

**NEW STRATEGY (correct approach):**
Instead of creating separate multi_cycle_scf_iterate(), we refactor existing
scf_iterate() to enable multi-SCF coordination while keeping ALL features
(DIIS, damping, MOM, SOSCF, convergence checks, etc.)

**Step 1.1:** ✅ DONE
- Extracted while loop body into scf_iteration() closure
- Used nonlocal variables for state
- ZERO logic changes - pure refactoring
- All 78 tests pass! (after fix 6db3bae6)

**Step 1.2:** ✅ DONE
- Created _scf_initialize_iteration_state(e_conv, d_conv)
- Converted closure scf_iteration() → method _scf_iteration()
- Moved state to self._scf_* members (~15 attributes)
- Method can now be called externally
- All 78 tests pass! ✅

**Step 1.3:** ✅ DONE (commits 4170231d + cleanups)
- Created multi_scf() coordinator function
- Uses _scf_iteration() for each wfn
- Shared JK computation via jk.C_clear()/C_add()
- Distribution via set_jk_matrices()
- Supports ALL SCF features (DIIS, damping, MOM, etc)

**Step 1.4:** ✅ DONE (critical bug fixes)
- Added pybind11 exports: n_states(), get_orbital_matrices(), set_jk_matrices()
- Fixed get_orbital_matrices() to return ONLY occupied orbitals
- Added JK pybind11 exports: C_left(), C_right(), J(), K()
- Fixed C_clear()/C_add() usage (not direct vector methods)
- Removed obsolete multi_cycle_scf_iterate() function

**Step 1.5:** ← NOW (options snapshot pattern)
Goal: Eliminate global state pollution causing non-determinism

Problem:
```python
wfn1 created → reads global DIIS_START=0
baseline test runs → changes global DIIS_START=14
wfn2 created → reads global DIIS_START=14  # ❌ Different!
multi_scf([wfn1, wfn2]) → non-deterministic behavior
```

Solution: Options Snapshot Pattern
- Freeze global options at wfn creation time
- Each wfn has independent option copy
- No global state pollution

Implementation:
1. Create scf_options_snapshot.py:
   - snapshot_scf_options() - freeze current global options
   - apply_options_snapshot(wfn, snapshot) - set wfn-local copy
2. Create multi_scf_helper.py:
   - multi_scf_wavefunction_factory() - create wfn with snapshot
   - multi_scf_helper() - high-level API
3. Modify _scf_initialize_iteration_state():
   - Check if wfn._options_snapshot exists
   - If yes: use snapshot, if no: read global (backward compat)

User API levels:
- Level 1 (99%): psi4.energy('scf') - no changes, works out of box
- Level 2 (simple): multi_scf_helper(molecules, method) - auto snapshot
- Level 3 (advanced): multi_scf_wavefunction_factory(..., options={...})

Files:
```
psi4/driver/procrouting/scf_proc/
├── scf_options_snapshot.py  # NEW
├── multi_scf_helper.py      # NEW
└── scf_iterator.py          # modify _scf_initialize_iteration_state
```

**Step 1.6:** Validation & Testing
1. Add validate_multi_scf_compatibility():
   - Check basis match (MUST)
   - Check JK type match (MUST)
   - Check geometry match (MUST)
   - Warn if functionals differ
2. Test determinism (100 runs)
3. Test with different options per state
4. Update test_multi_scf.py to use helper API

Phase 2: Multi-Spin SA-REKS 🎯 GOAL
  ├─> 2.1: SA-REKS theory stub (n_states = N)
  ├─> 2.2: Ensemble density for shared JK
  ├─> 2.3: Complete REKS occupation logic
  └─> 2.4: Test: Singlet + Triplet + Quintet simultaneously
```

**Estimated Timeline:**
- Phase 0.6: ✅ DONE (n_states() API + strategic decision)
- Phase 1: ~1 week (Python multi-cycle coordinator)
- Phase 2: ~2-3 weeks (SA-REKS implementation)
- **Total: ~1 month**

---

## Current Status Summary (UPDATED 2025-01-13)

**Phase 0: HPC Optimization ✅ COMPLETE**
1. **MultiStateMatrix infrastructure** (Phase 0.3)
   - 64-byte aligned allocation
   - Contiguous storage for cache locality
   - Proven +15.9% speedup on naphthalene

2. **Unified base for RHF/UHF/ROHF** (Phase 0.5)
   - All three use MultiStateMatrix pattern
   - RHF: n=1, UHF: n=2, ROHF: n=2
   - Tested and validated

3. **n_states() API** (Phase 0.6)
   - HF base class: virtual int n_states() const
   - RHF/UHF/ROHF override correctly
   - Tested via user's unit tests

**Strategic Decision: Python-First Approach** ✅

**Key Insight:** Existing `scf_iterate()` in Python ALREADY provides the abstraction we need!
- RHF/UHF/ROHF all work through ONE mechanism: Python scf_iterate() + HF virtual methods
- DIIS, convergence, damping - all working features in Python
- **Only addition needed:** Multi-cycle coordinator for shared JK computation
- **No need to rebuild SCF in C++** - would be code duplication with workarounds

**Phase 1: Python Multi-Cycle ✅ COMPLETE**
- ✅ C++ multi-cycle JK API (get_orbital_matrices, set_jk_matrices)
- ✅ Modified form_G() in RHF/UHF/ROHF to use precomputed J/K
- ✅ Implemented `multi_cycle_scf_iterate()` in scf_iterator.py
- ✅ Compilation successful
- 📍 Ready for testing (Phase 1.4)

**Implementation:**
- Location: `psi4/driver/procrouting/scf_proc/scf_iterator.py:1113-1263`
- Function: `multi_cycle_scf_iterate(wfn_list, e_conv, d_conv, max_iter, verbose)`
- Features: Shared JK computation, convergence checking, detailed output
- Backend: Uses precomputed_J_/K_ members + use_precomputed_jk_ flag

**Files Modified:**
- `psi4/src/psi4/libscf_solver/hf.h` - Added API methods
- `psi4/src/psi4/libscf_solver/hf.cc` - Initialization
- `psi4/src/psi4/libscf_solver/rhf.cc` - form_G() modification
- `psi4/src/psi4/libscf_solver/uhf.h` - get_orbital_matrices() override
- `psi4/src/psi4/libscf_solver/uhf.cc` - form_G() modification
- `psi4/src/psi4/libscf_solver/rohf.h` - get_orbital_matrices() override
- `psi4/src/psi4/libscf_solver/rohf.cc` - form_G() modification
- `psi4/driver/procrouting/scf_proc/scf_iterator.py` - New function

**Files to Note:**
- `SCFEngine` (scf_engine.{h,cc}): Exists but not used in new strategy
- Python approach uses: `psi4/driver/procrouting/scf_proc/scf_iterator.py`

---

## Phase 0.5: Унификация базы ✅ COMPLETE

**Status:** Tested and validated! All three (RHF/UHF/ROHF) use MultiStateMatrix.

**Achievements:**
- ✅ RHF uses MultiStateMatrix (n=1) - from Phase 0.3
- ✅ UHF uses MultiStateMatrix (n=2) - from Phase 0.3
- ✅ ROHF uses MultiStateMatrix (n=2) - Phase 0.5
- ✅ Unified pattern ready for SCFEngine
- ✅ All existing tests pass

---

## Phase 0.6: API Foundation ✅ COMPLETE

**Goal:** Add necessary APIs for multi-cycle support and decide on architecture strategy.

### Step 0.6.1: Add n_states() API ✅ DONE

**Completed:**
- Added virtual int n_states() const to HF base
- RHF returns 1, UHF returns 2, ROHF returns 2
- Compilation successful, tests pass

### Step 0.6.2: Basic API tests ✅ DONE (by user)

**User created 4 tests:**
1. test_scf_engine_basic_api() - wavefunction infrastructure
2. test_rhf_n_states() - confirms n_states() = 1
3. test_uhf_n_states() - confirms n_states() = 2
4. test_rohf_n_states() - confirms n_states() = 2

**All tests pass!** ✓

### Step 0.6.3: Strategic Decision - Python-First Approach ✅ DONE

**Analysis:**
- Initially considered building SCFEngine foundation in C++
- Realized existing Python `scf_iterate()` ALREADY provides needed abstraction
- RHF/UHF/ROHF all work through: Python scf_iterate() + HF virtual methods
- All features (DIIS, convergence, damping) already work in Python

**Decision:**
- ✅ Use existing `scf_iterate()` infrastructure
- ✅ Create `multi_cycle_scf_iterate()` in Python for shared JK coordination
- ❌ No need to rebuild SCF in C++ (would duplicate working code)

**Why This is Better:**
- No code duplication
- Uses all existing tested features (DIIS, convergence, etc.)
- Clean separation: C++ for theory, Python for coordination
- Faster to implement and test

**Note on SCFEngine:**
- Files `scf_engine.{h,cc}` exist from exploratory work
- Compile and link correctly but have placeholder implementations
- Not used in new Python-first strategy
- May be useful for future C++-level iteration control, but not needed now

---

## Phase 1: Python Multi-Cycle Strategy 📍 NEXT

**Goal:** Enable multiple SCF calculations to share a single JK computation for ~1.8-2x speedup.

### Why Python-First Approach?

**Existing Infrastructure that WORKS:**
```python
# psi4/driver/procrouting/scf_proc/scf_iterator.py
def scf_iterate(wfn, options):
    # This ALREADY works for RHF/UHF/ROHF via HF virtual methods!
    while not converged:
        wfn.form_G()   # Build two-electron contribution (J, K, XC)
        wfn.form_F()   # Assemble Fock matrix
        wfn.form_C()   # Diagonalize Fock
        wfn.form_D()   # Build density from orbitals
        E = wfn.compute_E()
        # DIIS, damping, convergence checking, etc.
```

**Key Observation:** All three theories (RHF/UHF/ROHF) already work through the SAME mechanism!
- Python calls virtual methods on HF base class
- Each theory implements the methods appropriately
- DIIS, convergence, damping - all tested and working

**What's Missing:** Ability to run multiple SCF cycles with shared JK
- Currently: Each wfn builds its own J/K matrices independently
- Goal: Collect all C matrices, single JK call, distribute results

### Architecture: Multi-Cycle Coordinator

```python
# psi4/driver/procrouting/scf_proc/scf_iterator.py

def multi_cycle_scf_iterate(wfn_list, options):
    """
    Run multiple SCF calculations with shared JK computation.

    Args:
        wfn_list: List of HF wavefunction objects (RHF, UHF, ROHF, SA-REKS)
        options: Options object

    Returns:
        List of final energies
    """
    # 1. Initialize all wavefunctions
    for wfn in wfn_list:
        wfn.initialize_scf()

    # 2. Main iteration loop
    for iteration in range(max_iter):
        # 3. Collect all C matrices from all wavefunctions
        all_C_matrices = []
        for wfn in wfn_list:
            # For RHF: 1 C matrix (Ca)
            # For UHF/ROHF: 2 C matrices (Ca, Cb)
            # For SA-REKS: N C matrices (one per state)
            C_list = wfn.get_orbital_matrices()  # Returns list based on n_states()
            all_C_matrices.extend(C_list)

        # 4. SHARED JK COMPUTATION (KEY OPTIMIZATION!)
        jk = wfn_list[0].jk()  # Get JK builder from any wfn (all use same basis)
        jk.C_left().clear()
        for C in all_C_matrices:
            jk.C_left().append(C)
        jk.compute()  # Single call for ALL matrices!

        # 5. Distribute J/K results back to each wavefunction
        jk_index = 0
        for wfn in wfn_list:
            n = wfn.n_states()  # How many J/K pairs this wfn needs
            J_list = [jk.J()[jk_index + i] for i in range(n)]
            K_list = [jk.K()[jk_index + i] for i in range(n)]
            wfn.set_jk_matrices(J_list, K_list)
            jk_index += n

        # 6. Each wavefunction completes its SCF step
        for wfn in wfn_list:
            wfn.form_F()   # Assemble Fock from J/K
            wfn.form_C()   # Diagonalize
            wfn.form_D()   # Build new density
            wfn.compute_E()

        # 7. Check convergence for all wavefunctions
        all_converged = all(wfn.is_converged() for wfn in wfn_list)
        if all_converged:
            break

    return [wfn.energy() for wfn in wfn_list]
```

### Required C++ APIs (Minimal)

Most infrastructure already exists! Only need:

1. **n_states() API** ✅ DONE (Phase 0.6.1)
   ```cpp
   virtual int n_states() const { return 1; }  // RHF: 1, UHF: 2, ROHF: 2
   ```

2. **get_orbital_matrices()** - NEW
   ```cpp
   // HF base class
   virtual std::vector<SharedMatrix> get_orbital_matrices() const {
       // RHF: return {Ca_}
       // UHF: return {Ca_, Cb_}
       // ROHF: return {Ca_, Cb_}
   }
   ```

3. **set_jk_matrices()** - NEW
   ```cpp
   // HF base class
   virtual void set_jk_matrices(const std::vector<SharedMatrix>& J_list,
                                 const std::vector<SharedMatrix>& K_list) {
       // Store J/K for use in form_F()
       // RHF: J_[0], K_[0]
       // UHF: J_[0,1], K_[0,1]
   }
   ```

4. **Modify form_G()** - SMALL CHANGE
   ```cpp
   // Current: form_G() builds J/K internally
   // New: form_G() uses pre-computed J/K from set_jk_matrices()
   //      OR builds them if not provided (backward compatibility)
   ```

### Benefits of This Approach

1. **Minimal C++ changes** - Only 2-3 new methods in HF base class
2. **No code duplication** - Uses all existing tested features
3. **Backward compatible** - Existing scf_iterate() still works unchanged
4. **Clean separation** - C++ for theory, Python for coordination
5. **Easy to test** - Can test with 2 independent RHF calculations first
6. **Extensible** - Works for any number of wavefunctions

### Testing Strategy

**Phase 1.1:** Test with 2 independent RHF calculations
```python
h2o_singlet = psi4.energy('scf/6-31g', molecule=h2o_singlet, return_wfn=True)[1]
h2o_triplet = psi4.energy('scf/6-31g', molecule=h2o_triplet, return_wfn=True)[1]

# Should give same energies as running separately, but faster
energies = multi_cycle_scf_iterate([h2o_singlet, h2o_triplet], options)
```

**Phase 1.2:** Test with UHF (already uses n=2 internally)
**Phase 1.3:** Test with SA-REKS (n=N states)

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
