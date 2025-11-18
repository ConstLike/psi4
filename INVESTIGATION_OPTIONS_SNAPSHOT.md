# Investigation: Options Snapshot Necessity

**Date**: 2025-11-18
**Task**: Determine if options snapshot system is necessary or can be simplified
**Method**: Code analysis of snapshot implementation and usage patterns

---

## Executive Summary

**VERDICT**: ⚠️ **PARTIALLY JUSTIFIED**

Options snapshot solves a **REAL but RARE problem** (global state pollution between wfn creation).

**Current implementation** is more complex than necessary:
- 34 call sites using `get_option_from_snapshot()`
- External module with 148 lines
- Partial caching in `_scf_initialize_iteration_state()`

**RECOMMENDATION**: Keep the concept but **SIMPLIFY** implementation.

---

## The Problem

### Stated Issue (from `scf_options_snapshot.py` lines 7-11):

```python
# Problem:
wfn1 created → reads global DIIS_START=0
baseline test runs → changes global DIIS_START=14
wfn2 created → reads global DIIS_START=14  # Different!
multi_scf([wfn1, wfn2]) → non-deterministic behavior
```

### Is This a Real Problem?

**YES**, but RARE:

**Scenario 1 (typical, 99% of cases)**:
```python
psi4.set_options({'DIIS_START': 1})
wfn1 = create_wfn(...)  # Reads DIIS_START=1
wfn2 = create_wfn(...)  # Reads DIIS_START=1
multi_scf([wfn1, wfn2])  # ← NO problem (same options)
```

**Scenario 2 (edge case, <1%)**:
```python
psi4.set_options({'DIIS_START': 1})
wfn1 = create_wfn(...)  # Reads DIIS_START=1

# Someone/something changes global options
psi4.set_options({'DIIS_START': 14})

wfn2 = create_wfn(...)  # Reads DIIS_START=14 (DIFFERENT!)
multi_scf([wfn1, wfn2])  # ← PROBLEM!
```

**When can this happen?**
1. **Interactive sessions**: User creates wfn, then fiddles with options, then creates more wfn
2. **Test suites**: Previous test modifies global options, affects next test
3. **Library usage**: Calling code changes options between wfn creations

**Conclusion**: Problem is REAL for robustness, but won't affect typical usage.

---

## Current Implementation

### Architecture

**3 functions in `scf_options_snapshot.py`**:
1. `snapshot_scf_options()` - Freeze 35+ options into dict
2. `apply_options_snapshot(wfn, snapshot)` - Store dict in wfn._options_snapshot
3. `get_option_from_snapshot(wfn, name)` - Read from snapshot or fallback to global

### Usage in `scf_iterator.py`

**Snapshot creation** (`multi_scf()` line 1628-1635):
```python
# Snapshot global options
options_snapshot = snapshot_scf_options()

# Apply to each wfn
for wfn in wfn_list:
    apply_options_snapshot(wfn, options_snapshot)
```

**Option reads** (34 call sites):
```python
diis_start = get_option_from_snapshot(self, 'DIIS_START')
soscf_conv = get_option_from_snapshot(self, 'SOSCF_CONV')
# ... 32 more calls
```

### Metrics

| Component | Lines | Complexity |
|-----------|-------|------------|
| `scf_options_snapshot.py` | 148 | LOW-MED |
| Call sites in `scf_iterator.py` | 34 | LOW each |
| Total impact | ~180 lines | MEDIUM |

---

## Analysis

### What Options Are Snapshotted?

**35 options** across categories:

1. **DIIS** (8 options): DIIS, DIIS_START, DIIS_MIN_VECS, DIIS_MAX_VECS, DIIS_RMS_ERROR, + AEDIIS options
2. **Damping** (2): DAMPING_PERCENTAGE, DAMPING_CONVERGENCE
3. **SOSCF** (6): SOSCF, SOSCF_START_CONVERGENCE, SOSCF_CONV, SOSCF_MIN/MAX_ITER, SOSCF_PRINT
4. **MOM** (2): MOM_START, MOM_OCC
5. **FRAC** (4): FRAC_START, FRAC_OCC, FRAC_VAL, FRAC_RENORMALIZE
6. **Convergence** (4): MAXITER, E_CONVERGENCE, D_CONVERGENCE, FAIL_ON_MAXITER
7. **Other** (6): LEVEL_SHIFT, LEVEL_SHIFT_CUTOFF, PRINT, COSX_MAXITER_FINAL

**Critical for correctness**:
- DIIS_START - when DIIS kicks in
- E_CONVERGENCE, D_CONVERGENCE - convergence thresholds
- MAXITER - iteration limit

**Less critical**:
- PRINT - verbosity (doesn't affect results)
- SOSCF options - if enabled, affects algorithm
- MOM/FRAC - advanced features

---

## Pros and Cons

### Pros ✅

1. **Robustness**: Handles global state pollution gracefully
2. **Determinism**: Guarantees all wfn use same convergence parameters
3. **Backward compatible**: Falls back to global if no snapshot (single-cycle SCF works)
4. **Clean API**: Simple `snapshot → apply → use` pattern
5. **Solves real problem**: Even if rare, prevents hard-to-debug issues

### Cons ❌

1. **Complexity**: 34 call sites need `get_option_from_snapshot()` instead of `core.get_option()`
2. **External module**: 148-line module for what could be simpler
3. **Partial caching**: Some options cached in `_scf_initialize_iteration_state()`, others not
4. **Inconsistency**: Some code uses `core.get_option()`, some uses `get_option_from_snapshot()`
5. **Runtime overhead**: Dict lookup + fallback for every option read (negligible, but exists)

---

## Alternative Solutions

### Option 1: Documentation Only ❌

**Approach**:
```python
def multi_scf(wfn_list):
    """
    ...

    Important: All wavefunctions must be created with the same global options.
    Do not modify global options between wfn creation!
    """
```

**Pros**: Zero code
**Cons**: User error-prone, doesn't actually prevent problem

**Verdict**: ❌ Too fragile

---

### Option 2: Validation Only ⚠️

**Approach**:
```python
def validate_options_match(wfn_list):
    ref_diis_start = ref_wfn.diis_start_
    for wfn in wfn_list[1:]:
        if wfn.diis_start_ != ref_diis_start:
            raise ValidationError("DIIS_START mismatch!")
    # ... check other options
```

**Pros**:
- Simpler than snapshot (no `get_option_from_snapshot` calls)
- Catches problem early with clear error

**Cons**:
- Requires options to be cached in wfn at creation (C++ changes?)
- Still complex validation logic
- Doesn't FIX problem, just detects it

**Verdict**: ⚠️ Better than nothing, but doesn't solve problem

---

### Option 3: Context Manager ⚠️

**Approach**:
```python
with freeze_global_options():
    wfn1 = create_wfn(...)
    wfn2 = create_wfn(...)
    energies = multi_scf([wfn1, wfn2])
```

**Pros**:
- Pythonic, explicit
- Prevents option changes during critical section

**Cons**:
- Doesn't help if wfn already created with different options
- Requires user cooperation (must use context manager)

**Verdict**: ⚠️ Good for prevention, but not robust against pre-existing wfn

---

### Option 4: Current Snapshot ✅

**Approach**: (as implemented)

**Pros**:
- Handles ALL scenarios robustly
- Works regardless of when/how wfn created
- Python-only fix (no C++ changes)

**Cons**:
- Complex implementation (34 call sites)

**Verdict**: ✅ Most robust, but could be simpler

---

## Simplification Opportunities

### Current Partial Caching

`_scf_initialize_iteration_state()` (lines 273-308) already caches some options:

```python
self._scf_e_conv = e_conv
self._scf_d_conv = d_conv
self._scf_is_dfjk = ...
self._scf_verbose = get_option_from_snapshot(self, 'PRINT')
self._scf_soscf_enabled = ...
self._scf_damping_enabled = ...
# ... etc
```

**Problem**: Not all options are cached, so still need 34 `get_option_from_snapshot()` calls elsewhere.

### Proposed Simplification 1: Full Caching

**Approach**: Cache ALL snapshot options as wfn._scf_* members in `_scf_initialize_iteration_state()`:

```python
def _scf_initialize_iteration_state(self, e_conv, d_conv):
    # Cache ALL options from snapshot
    self._scf_diis_start = get_option_from_snapshot(self, 'DIIS_START')
    self._scf_diis_max_vecs = get_option_from_snapshot(self, 'DIIS_MAX_VECS')
    self._scf_soscf_conv = get_option_from_snapshot(self, 'SOSCF_CONV')
    # ... etc for all 35 options

# Then in other code:
# OLD: diis_start = get_option_from_snapshot(self, 'DIIS_START')
# NEW: diis_start = self._scf_diis_start
```

**Pros**:
- Reduces 34 call sites to ~1 (all in one place)
- Faster (direct member access vs dict lookup)
- Clearer (explicit caching)

**Cons**:
- More members on wfn object
- Still need snapshot infrastructure

**Impact**:
- Remove: 34 scattered `get_option_from_snapshot()` calls
- Add: 35 cache assignments in `_scf_initialize_iteration_state()`
- **Net**: Much cleaner!

---

### Proposed Simplification 2: Snapshot at wfn Creation

**Approach**: Apply snapshot when wfn is CREATED, not in `multi_scf()`:

```python
# In wavefunction factory or initialization
def create_wfn(...):
    wfn = core.RHF(...)

    # Immediately snapshot and cache options
    wfn._scf_diis_start = core.get_option('SCF', 'DIIS_START')
    wfn._scf_diis_max_vecs = core.get_option('SCF', 'DIIS_MAX_VECS')
    # ... etc

    return wfn
```

**Pros**:
- Options frozen at creation (most intuitive)
- No need for `apply_options_snapshot()` in `multi_scf()`
- Even simpler

**Cons**:
- Requires finding all wfn creation points
- May need C++ changes for robust implementation

**Impact**:
- Remove: `snapshot_scf_options()` and `apply_options_snapshot()` calls from `multi_scf()`
- Add: Option caching at wfn creation
- **Net**: Cleaner separation of concerns

---

### Proposed Simplification 3: Hybrid - Lazy Snapshot

**Approach**: Snapshot on first option read:

```python
def get_option(wfn, option_name):
    if not hasattr(wfn, '_scf_options_cache'):
        wfn._scf_options_cache = {}

    if option_name not in wfn._scf_options_cache:
        # First read - cache it
        wfn._scf_options_cache[option_name] = core.get_option('SCF', option_name)

    return wfn._scf_options_cache[option_name]
```

**Pros**:
- No explicit snapshot call needed
- Only caches options actually used
- Automatic

**Cons**:
- Implicit behavior (less clear when snapshot happens)
- Still 34 call sites (just different function name)

**Impact**:
- Similar to current, but more automatic

---

## Recommended Approach

### Short-term (1-2 months): ✅ Keep current snapshot

**Reason**:
- Already working
- Solves real problem
- Not performance-critical
- Risk of breaking existing functionality

**Minor cleanup**:
- Add comment explaining WHY snapshot is needed (link to this document)
- Consider moving snapshot earlier (before wfn naming, etc.)

---

### Long-term (3-6 months): ⚠️ Implement Simplification 1 (Full Caching)

**Approach**:
1. Cache ALL 35 options as `self._scf_*` members in `_scf_initialize_iteration_state()`
2. Replace 34 `get_option_from_snapshot()` calls with direct member access
3. Keep `snapshot_scf_options()` and `apply_options_snapshot()` infrastructure

**Benefits**:
- **Cleaner code**: Single location for caching (one function)
- **Faster**: Direct member access vs dict lookup
- **More explicit**: Clear what options are cached
- **Easier to maintain**: All caching logic in one place

**Effort**: ~2-3 hours
- Add 35 cache assignments to `_scf_initialize_iteration_state()`
- Replace 34 `get_option_from_snapshot()` calls with `self._scf_*`
- Test to verify no regressions

**Example**:
```python
# In _scf_initialize_iteration_state():
def _scf_initialize_iteration_state(self, e_conv, d_conv):
    # ... existing code ...

    # Cache all snapshot options (ONCE)
    self._scf_diis_start = get_option_from_snapshot(self, 'DIIS_START')
    self._scf_diis_min_vecs = get_option_from_snapshot(self, 'DIIS_MIN_VECS')
    self._scf_diis_max_vecs = get_option_from_snapshot(self, 'DIIS_MAX_VECS')
    self._scf_soscf_conv = get_option_from_snapshot(self, 'SOSCF_CONV')
    # ... all 35 options ...

# Everywhere else:
# BEFORE: nmicro = self.soscf_update(get_option_from_snapshot(self, 'SOSCF_CONV'), ...)
# AFTER:  nmicro = self.soscf_update(self._scf_soscf_conv, ...)
```

---

## Verdict

### Is options snapshot necessary?

**YES** ✅ - Solves real problem (global state pollution)

### Is current implementation optimal?

**NO** ⚠️ - Can be simplified with full caching

### Should we remove it?

**NO** ❌ - Problem is real, snapshot is best solution

### Should we simplify it?

**YES** ✅ - Long-term refactoring recommended

---

## Final Recommendation

**Keep options snapshot** but **simplify implementation** via full caching:

1. **Short-term** (now): Keep as-is, add documentation
2. **Long-term** (3-6 months): Implement Simplification 1 (full caching)

**Impact after simplification**:
- From: 34 scattered `get_option_from_snapshot()` calls
- To: 35 cached assignments in ONE function
- **Net**: Much cleaner, easier to maintain

---

## Conclusion

Options snapshot is **JUSTIFIED** but **OVERCOMPLICATED** in current form.

**Philosophy**:
> "The problem is real (global state pollution), the solution is correct (snapshot),
> but the implementation could be simpler (full caching instead of scattered calls)."

Not a high-priority simplification, but a good refactoring opportunity for the future.

---

## Code Locations

**Snapshot module**: `scf_options_snapshot.py` (148 lines)
**Usage in multi_scf**: `scf_iterator.py:1628-1635`
**Call sites**: 34 locations in `scf_iterator.py`
**Partial caching**: `_scf_initialize_iteration_state()` lines 273-308

---

## Testing Requirements

If simplifying with full caching:

```python
def test_options_snapshot_prevents_pollution():
    """Verify snapshot prevents global option pollution"""
    psi4.set_options({'DIIS_START': 1})
    wfn1 = create_wfn(...)

    # Change global (simulating pollution)
    psi4.set_options({'DIIS_START': 14})
    wfn2 = create_wfn(...)

    # Apply snapshot
    snapshot = snapshot_scf_options()  # Should be DIIS_START=1 (before pollution)
    apply_options_snapshot(wfn1, snapshot)
    apply_options_snapshot(wfn2, snapshot)

    # Both should use DIIS_START=1, not 14
    energies = multi_scf([wfn1, wfn2])
    # Verify both used same DIIS_START (check logs or internal state)
```

**Investigation complete!** ⚠️
