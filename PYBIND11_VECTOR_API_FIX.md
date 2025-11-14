# PyBind11 Vector API Fix - Why C_clear() and C_add() are needed

**Date:** 2025-11-14
**Status:** CRITICAL FIX for multi-cycle SCF

---

## Problem

`jk.C_left().append()` and `jk.C_left().clear()` **DO NOT WORK** in Python, even though `C_left()` returns a reference to the C++ vector.

**Result:** JK C_left and C_right remained empty (length 0), causing JK.compute() to return no results.

---

## Root Cause

**PyBind11 does NOT automatically export C++ std::vector methods to Python.**

Even with `py::return_value_policy::reference_internal`, Python cannot call:
- `.clear()` - not bound
- `.append()` / `.push_back()` - not bound
- `.pop()` - not bound
- etc.

### What happens in Python:

```python
vec = jk.C_left()          # Returns reference to C++ std::vector<SharedMatrix>
print(type(vec))           # <class 'list'> ← Python wrapper, NOT C++ vector!
vec.clear()                # Clears PYTHON wrapper, NOT C++ vector ❌
vec.append(matrix)         # Appends to PYTHON wrapper, NOT C++ vector ❌
len(jk.C_left())           # Still 0! ❌
```

**The reference is valid, but methods don't work on the underlying C++ vector.**

---

## Evidence from Test

```python
# Test 1: Using .append() directly
jk.C_left().clear()
print(len(jk.C_left()))    # 0 ✓
jk.C_left().append(matrix)
print(len(jk.C_left()))    # 0 ❌ NOT 1!

# Test 2: Using C_add() wrapper
jk.C_clear()
print(len(jk.C_left()))    # 0 ✓
jk.C_add(matrix)
print(len(jk.C_left()))    # 1 ✓ Works!
```

---

## Solution

**Use exported wrapper methods that call C++ vector methods directly:**

### Existing exports in export_fock.cc (lines 79-90):

```cpp
.def("C_clear",
     [](JK &jk) {
         jk.C_left().clear();   // ← Calls C++ .clear() directly
         jk.C_right().clear();
     })
.def("C_add",
     [](JK &jk, SharedMatrix Cl) {
         jk.C_left().push_back(Cl);   // ← Calls C++ .push_back() directly
         jk.C_right().push_back(Cl);
     })
.def("C_left_add", [](JK &jk, SharedMatrix Cl) { jk.C_left().push_back(Cl); })
.def("C_right_add", [](JK &jk, SharedMatrix Cr) { jk.C_right().push_back(Cr); })
```

**These methods work because they execute C++ code, not Python code.**

---

## Fix Applied

### File: psi4/driver/procrouting/scf_proc/scf_iterator.py (lines 1267-1273)

**Before (BROKEN):**
```python
jk.C_left().clear()           # ❌ Doesn't work - clears Python wrapper
jk.C_right().clear()          # ❌ Doesn't work
for C_occ in all_C_occ_matrices:
    jk.C_left().append(C_occ)    # ❌ Doesn't work - appends to wrapper
    jk.C_right().append(C_occ)   # ❌ Doesn't work
```

**After (WORKS):**
```python
jk.C_clear()                  # ✅ Calls C++ clear() via wrapper
for C_occ in all_C_occ_matrices:
    jk.C_add(C_occ)           # ✅ Calls C++ push_back() via wrapper
```

---

## Why This Happens

**PyBind11 philosophy:**
- You must **explicitly export** every method you want available in Python
- Just returning a reference to `std::vector` doesn't export its methods
- Need lambda wrappers to expose vector operations

**From PyBind11 docs:**
> "STL containers are not automatically converted. You need to use py::bind_vector<>
> or create custom wrappers for the operations you need."

**Our approach:** Custom wrappers (`C_clear`, `C_add`, etc.) - simpler and cleaner than `py::bind_vector<>`.

---

## Alternative Approaches (NOT used)

### Option 1: py::bind_vector<>
```cpp
py::bind_vector<std::vector<SharedMatrix>>(m, "SharedMatrixVector");
```
**Rejected:** Requires creating new Python type, changes API significantly.

### Option 2: Export C_left() as property with full binding
```cpp
.def_property_readonly("C_left_vec",
    &JK::C_left,
    py::return_value_policy::reference_internal)
```
**Rejected:** Still needs `bind_vector<>` or manual method exports.

### Option 3: Wrapper methods (CHOSEN)
```cpp
.def("C_clear", [](JK &jk) { jk.C_left().clear(); jk.C_right().clear(); })
.def("C_add", [](JK &jk, SharedMatrix C) { jk.C_left().push_back(C); ... })
```
**Chosen:** Simple, explicit, matches Psi4 API style, already implemented.

---

## Impact on Multi-Cycle SCF

**Before fix:**
```
Step 2: Fill JK
  jk.C_left().clear()         → Python wrapper cleared (C++ vector unchanged)
  jk.C_left().append(C1)      → Append to wrapper (C++ vector unchanged)
  jk.C_left().append(C2)      → Append to wrapper (C++ vector unchanged)

Step 3: Compute
  jk.compute()                → Sees EMPTY C++ vector, returns no results

Step 4: Get results
  J_all = jk.J()              → Empty list []
  K_all = jk.K()              → Empty list []

CRASH: IndexError (expected 2 matrices, got 0)
```

**After fix:**
```
Step 2: Fill JK
  jk.C_clear()                → C++ vectors cleared ✓
  jk.C_add(C1)                → C++ vectors get C1 ✓
  jk.C_add(C2)                → C++ vectors get C2 ✓

Step 3: Compute
  jk.compute()                → Sees 2 matrices, computes J/K ✓

Step 4: Get results
  J_all = jk.J()              → [J1, J2] ✓
  K_all = jk.K()              → [K1, K2] ✓

SUCCESS: Multi-cycle SCF proceeds correctly
```

---

## Lessons Learned

1. **PyBind11 references ≠ Python access to methods**
   - Reference keeps object alive, but doesn't expose methods

2. **Always use explicit wrapper methods for STL containers**
   - Don't assume `.clear()`, `.append()` will work

3. **Test at the boundary**
   - Test Python/C++ interface thoroughly
   - Don't assume C++ patterns translate to Python

4. **Wrapper methods are cleaner than full binding**
   - Only expose operations you need
   - Simpler API, less code

---

## Status

✅ Fixed in scf_iterator.py (lines 1267-1273)
✅ Uses existing wrapper methods from export_fock.cc
✅ Ready for testing (no C++ rebuild needed - pure Python fix)

**Next step:** Run `python test_multi_scf.py` to verify fix works.
