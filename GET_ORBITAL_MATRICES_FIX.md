# get_orbital_matrices() Fix - Occupied vs Full Orbitals

**Date:** 2025-11-14
**Status:** CRITICAL FIX for multi-cycle SCF

---

## Problem

`get_orbital_matrices()` was returning **FULL orbital matrices** instead of **ONLY OCCUPIED orbitals**.

**Result:** JK compute failed because JK expects only occupied orbitals, causing empty J/K lists.

---

## Root Cause

**In normal SCF (form_G()):** ALL classes use **ONLY occupied orbitals**
- **RHF:** `Ca_subset("SO", "OCC")` - only occupied (line rhf.cc:220)
- **UHF:** `Ca_subset("SO", "OCC")`, `Cb_subset("SO", "OCC")` (line uhf.cc:223-224)
- **ROHF:** `get_block()` for docc/socc - only occupied (line rohf.cc:971-976)

**In multi-cycle SCF (get_orbital_matrices()):** Returned **FULL matrices**
- **HF (base):** `Ca_` - full matrix (nso x nmo) ❌
- **UHF:** `{Ca_, Cb_}` - full matrices ❌
- **ROHF:** `{Ca_, Cb_}` - full matrices ❌

**Mismatch:** JK expects occupied orbitals, but got full matrices!

---

## Evidence from Test

**Molecule:** H₂O, cc-pVDZ
- nso = 24 (basis functions)
- nmo = 24 (molecular orbitals)
- nocc = 5 (occupied orbitals, 10 electrons / 2)

**Before fix:**
```python
get_orbital_matrices()[0].shape = (24, 24)  # ❌ WRONG - full matrix
```

**After fix:**
```python
get_orbital_matrices()[0].shape = (24, 5)   # ✅ CORRECT - only occupied
```

---

## Fix Applied

### File: psi4/src/psi4/libscf_solver/hf.h (lines 333-343)

**Before:**
```cpp
virtual std::vector<SharedMatrix> get_orbital_matrices() const {
    return {Ca_};  // ❌ Full matrix
}
```

**After:**
```cpp
virtual std::vector<SharedMatrix> get_orbital_matrices() const {
    return {Ca_subset("SO", "OCC")};  // ✅ Only occupied
}
```

### File: psi4/src/psi4/libscf_solver/uhf.h (lines 86-90)

**Before:**
```cpp
std::vector<SharedMatrix> get_orbital_matrices() const override {
    return {Ca_, Cb_};  // ❌ Full matrices
}
```

**After:**
```cpp
std::vector<SharedMatrix> get_orbital_matrices() const override {
    return {Ca_subset("SO", "OCC"), Cb_subset("SO", "OCC")};  // ✅ Only occupied
}
```

### File: psi4/src/psi4/libscf_solver/rohf.h (lines 91-99)

**Before:**
```cpp
std::vector<SharedMatrix> get_orbital_matrices() const override {
    return {Ca_, Cb_};  // ❌ Full matrices
}
```

**After:**
```cpp
std::vector<SharedMatrix> get_orbital_matrices() const override {
    // ROHF needs separate docc and socc blocks (same as form_G())
    Dimension dim_zero(nirrep_, "Zero Dim");
    SharedMatrix Cdocc = Ca_->get_block({dim_zero, nsopi_}, {dim_zero, nbetapi_});
    SharedMatrix Csocc = Ca_->get_block({dim_zero, nsopi_}, {nbetapi_, nalphapi_});
    return {Cdocc, Csocc};  // ✅ Only occupied (docc + socc)
}
```

---

## Why This Matters

**JK algorithm expects only occupied orbitals:**
1. Forms density: `D = C_occ * C_occ^T`
2. Computes J: `J[μν] = (μν|ρσ) D[ρσ]`
3. Computes K: `K[μν] = (μρ|νσ) D[ρσ]`

**If full C matrix is passed:**
- Density includes virtual orbitals (wrong!)
- J/K computation fails or produces garbage
- Multi-cycle SCF breaks

**This is exactly what normal SCF does:**
- RHF: `jk_->C_left().push_back(Ca_subset("SO", "OCC"))`
- UHF: `jk_->C_left().push_back(Ca_subset("SO", "OCC"))` + beta
- ROHF: `jk_->C_left().push_back(Cdocc)` + `Csocc`

**Multi-cycle SCF now does the same thing!**

---

## Timeline

1. **Phase 0.6:** Added `get_orbital_matrices()` API
2. **Initial implementation:** Returned full Ca_ (copy-paste from other methods)
3. **Testing revealed:** JK compute returned empty lists
4. **Root cause:** Full matrices instead of occupied
5. **Fix:** Changed to return only occupied (matching form_G() logic)

---

## Status

✅ Fixed in all three classes (HF base, UHF, ROHF)
✅ Logic matches normal SCF form_G() exactly
✅ Ready for testing after C++ rebuild

**Next step:** Rebuild C++ and test `python test_multi_scf.py`
