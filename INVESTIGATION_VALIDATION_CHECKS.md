# Investigation: Validation MUST Checks (do_wK, alpha, beta)

**Date**: 2025-11-18
**Task**: Verify which validation checks are truly MUST for shared JK
**Method**: Deep code analysis of C++ implementation (JK, MemDFJK, DFHelper)

---

## Executive Summary

**FINDINGS**:

| Check | MUST Match? | Justification |
|-------|-------------|---------------|
| **omega** | ✅ **YES** | Used in erf(ωr)/r operator for wK integrals |
| **do_wK** | ✅ **YES** | Capacity baked at build - cannot change after |
| **alpha** | ❌ **NO** | Used ONLY if wcombine=TRUE (default=FALSE) |
| **beta** | ❌ **NO** | Used ONLY if wcombine=TRUE (default=FALSE) |

**ACTION REQUIRED**: Remove alpha/beta MUST checks from validation function!

---

## Investigation Method

Traced through C++ code to understand:
1. When parameters are used (build vs runtime)
2. Whether parameters affect INTEGRALS or just Fock assembly
3. Which parameters can be reconfigured vs baked-in

**Files analyzed**:
- `/home/user/psi4/psi4/src/psi4/libfock/jk.cc` - JK base class
- `/home/user/psi4/psi4/src/psi4/libfock/MemDFJK.cc` - DF JK implementation
- `/home/user/psi4/psi4/src/psi4/lib3index/dfhelper.cc` - 3-index integral engine

---

## Finding 1: omega ✅ MUST Match

### Code Evidence

**Usage in DFHelper** (`dfhelper.cc:595`):
```cpp
weri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
```

omega is the **range-separation parameter** passed to erf_eri() to create long-range integral engine.

**Physics**:
```
wK_μν = Σ_ρσ P_ρσ (μρ|erf(ωr)/r|νσ)
```

Different omega → different erf(ωr) → **different integrals** → WRONG RESULTS!

**Reconfigured in multi_scf?** NO
- Set in `initialize_jk()` line 116: `jk.set_omega(...)`
- NOT called in `_multi_scf_inner()`
- wfn[1:] skip initialize_jk() due to idempotency

### Verdict: ✅ **MUST match** (JUSTIFIED)

---

## Finding 2: do_wK ✅ MUST Match

### Code Evidence

**At JK build** (`jk.cc` common_init):
```cpp
do_wK_ = false;  // Default
```

**At initialize** (`MemDFJK.cc:99`):
```cpp
dfh_->set_do_wK(do_wK_);
// ...
dfh_->initialize();
```

**In DFHelper initialize** (`dfhelper.cc:191-193`):
```cpp
if (do_wK_) {
    prepare_AO_wK_core();  // Allocate memory, prepare wK integrals
} else {
    // wK NOT prepared!
}
```

**Critical**: `prepare_AO_wK_core()` is called **ONLY IF do_wK_=True** at initialize time!

### What happens if do_wK doesn't match?

**Scenario**:
- wfn[0]: non-LRC → JK built with `do_wK=False` → wK NOT prepared
- wfn[1]: LRC functional → needs wK → **ERROR!** (wK not prepared)

**Result**: JK cannot compute wK for wfn[1] because memory and integrals not prepared!

### Can set_do_wK() change this after build?

**NO!** Testing code path:

1. `JK.build(..., do_wK=False)` → creates JK
2. `jk.initialize()` → calls `prepare_AO_wK_core()` ONLY if do_wK_=True
3. If do_wK_=False → wK integral storage NOT allocated
4. Later `jk.set_do_wK(True)` → changes flag but NOT storage!
5. `jk.compute()` → tries to compute wK → **CRASH** (no storage!)

### Verdict: ✅ **MUST match** (JUSTIFIED)

**Note**: Current check `wfn_is_lrc != ref_is_lrc` is correct.
Alternative (more permissive): Allow if `ref_is_lrc or not any(wfn.is_lrc() for wfn in wfn_list[1:])` (i.e., ref has capability even if others don't need it). But strict equality is safer.

---

## Finding 3: alpha ❌ CAN Differ

### Code Evidence

**wcombine flag** (`jk.cc` common_init):
```cpp
wcombine_ = false;  // DEFAULT!
```

**In DFHelper prepare** (`dfhelper.cc:655-676`):
```cpp
if ( wcombine_ ) {
    // computes (Q|mn) and (Q|w|mn)
    compute_sparse_pQq_blocking_p_symm_abw(start, stop, M1p, M2p, eri, weri);
    // ^^ Uses alpha/beta here!
} else {
    // computes just (Q|mn)
    compute_sparse_pQq_blocking_p_symm(start, stop, M1p, eri);
    // ^^ NO alpha/beta!
}
```

**In compute_sparse_pQq_blocking_p_symm_abw** (`dfhelper.cc:1387`):
```cpp
just_Mp[ind1] = buffer[rank][...];  // "Just" K integral
param_Mp[ind1] = omega_alpha_ * buffer[...] + omega_beta_ * wbuffer[...];
// ^^ Parametrized: alpha*K + beta*wK
```

**Two branches**:

1. **wcombine=FALSE (default)**:
   - Calls `compute_sparse_pQq_blocking_p_symm()` (NO alpha/beta)
   - Computes K and wK separately
   - Returns "pure" K and wK matrices
   - alpha/beta used **LATER** in Fock assembly (not in integrals!)
   - **Different alpha/beta → same integrals → OK!**

2. **wcombine=TRUE**:
   - Calls `compute_sparse_pQq_blocking_p_symm_abw()` (WITH alpha/beta)
   - Computes combined matrix: alpha×K + beta×wK
   - Returns single parametrized matrix
   - alpha/beta **baked into integrals**
   - **Different alpha/beta → different integrals → WRONG!**

### Which mode does Psi4 use?

**Checked**: wcombine is NOT set in `/home/user/psi4/psi4/driver/procrouting/scf_proc/`
**Result**: Uses **default wcombine=FALSE**

### Verdict: ❌ **CAN differ** (NOT MUST in default mode)

**Recommendation**: Remove alpha MUST check from validation!

---

## Finding 4: beta ❌ CAN Differ

### Code Evidence

Same as alpha (see Finding 3).

**Usage**:
```cpp
param_Mp[ind1] = omega_alpha_ * buffer[...] + omega_beta_ * wbuffer[...];
```

Only used if wcombine=TRUE.
Default wcombine=FALSE → beta NOT used in integrals.

### Verdict: ❌ **CAN differ** (NOT MUST in default mode)

**Recommendation**: Remove beta MUST check from validation!

---

## Python Binding Documentation

From `/home/user/psi4/psi4/src/export_fock.cc`:

```cpp
.def("set_omega_alpha", &JK::set_omega_alpha,
     "Weight for HF exchange term in range-separated DFT", "alpha"_a)
.def("set_omega_beta", &JK::set_omega_beta,
     "Weight for dampened exchange term in range-separated DFT", "beta"_a)
```

Documentation says "**Weight**" - i.e., **coefficients**, not integral parameters!

This aligns with finding that alpha/beta are mixing coefficients, not integral operators (in default wcombine=FALSE mode).

---

## Physical Interpretation

### RSH Functional Energy

Typical range-separated hybrid (RSH) functional:

```
E_xc = a_SR * E_x^HF,SR(ω) + a_LR * E_x^HF,LR(ω) + E_xc^DFT
```

where:
- E_x^HF,SR - short-range HF exchange with erfc(ωr)/r
- E_x^HF,LR - long-range HF exchange with erf(ωr)/r
- a_SR, a_LR - mixing coefficients (alpha, beta)

### Integral Computation

**Standard approach** (wcombine=FALSE):
1. Compute K (standard exchange, 1/r)
2. Compute wK (long-range exchange, erf(ωr)/r)
3. Compute K^SR = K - wK (using identity: erfc + erf = 1)
4. Form Fock: F = H + J + a_SR×K^SR + a_LR×wK + V_xc

alpha/beta used in **step 4** (Fock assembly), not in steps 1-2 (integrals).

**Optimized approach** (wcombine=TRUE):
1. Compute combined: K_combined = alpha×K + beta×wK directly
2. Form Fock: F = H + J + K_combined + V_xc

alpha/beta used in **step 1** (integral computation).

**Psi4 default**: Standard approach (wcombine=FALSE)

---

## Reconfiguration in multi_scf

Checked `_multi_scf_inner()` lines 1777-1778:

```python
jk.set_do_J(True)
jk.set_do_K(True)
```

**NOT called**:
- `jk.set_do_wK()` - NOT reconfigured!
- `jk.set_omega()` - NOT reconfigured!
- `jk.set_omega_alpha()` - NOT reconfigured!
- `jk.set_omega_beta()` - NOT reconfigured!

**Reason**: These are baked into integrals at initialize time.
Only do_J and do_K are runtime flags (always set to True for safety).

---

## Summary of Defaults

From `jk.cc` common_init():
```cpp
do_J_ = true;
do_K_ = true;
do_wK_ = false;
wcombine_ = false;    // ← KEY!
lr_symmetric_ = false;
omega_ = 0.0;
omega_alpha_ = 1.0;
omega_beta_ = 0.0;
```

**wcombine=FALSE is default** → alpha/beta NOT used in integrals!

---

## Recommendations

### Remove from validation (NOT MUST):
1. ❌ **alpha check** (lines 1462-1486) - Only needed if wcombine=TRUE
2. ❌ **beta check** (lines 1488-1519) - Only needed if wcombine=TRUE

### Keep in validation (MUST):
1. ✅ **omega check** (lines 1429-1460) - Always affects wK integrals
2. ✅ **do_wK check** (lines 1373-1408) - Capacity baked at build

### Optional Enhancement:

Add wcombine support in future:
```python
# If wcombine is ever enabled:
if any_wfn_uses_wcombine():
    # THEN check alpha/beta
    validate_alpha_beta_match()
```

But for current code (wcombine=FALSE always), alpha/beta checks are **UNNECESSARY**.

---

## Testing Recommendations

### Test 1: Verify alpha/beta CAN differ (should PASS)
```python
mol = psi4.geometry(...)

# Two wavefunctions with SAME omega but DIFFERENT alpha
wfn1 = create_wfn_with_functional(omega=0.3, alpha=0.2, beta=0.8)
wfn2 = create_wfn_with_functional(omega=0.3, alpha=0.3, beta=0.7)

# Should NOT raise error (alpha/beta used in Fock, not integrals)
energies = multi_scf([wfn1, wfn2])
```

### Test 2: Verify omega MUST match (should FAIL)
```python
wfn1 = create_wfn_with_functional(omega=0.3, ...)
wfn2 = create_wfn_with_functional(omega=0.4, ...)

# Should raise ValidationError
with pytest.raises(ValidationError, match="omega"):
    multi_scf([wfn1, wfn2])
```

### Test 3: Verify do_wK MUST match (should FAIL)
```python
wfn1 = create_wfn(functional="HF")      # non-LRC, do_wK=False
wfn2 = create_wfn(functional="wB97X")   # LRC, do_wK=True

# Should raise ValidationError
with pytest.raises(ValidationError, match="LRC capability"):
    multi_scf([wfn1, wfn2])
```

---

## Conclusion

**Current validation is TOO STRICT!**

Checks that should be REMOVED:
- alpha MUST match → **NOT needed** (wcombine=FALSE default)
- beta MUST match → **NOT needed** (wcombine=FALSE default)

Checks that are CORRECT:
- omega MUST match → **Justified** (affects wK integrals)
- do_wK MUST match → **Justified** (capacity baked at build)

**Validation function should be simplified from 317 lines → ~200 lines** by removing unnecessary alpha/beta checks.

---

## Code Locations

**Validation function**: `scf_iterator.py:1203-1519` (317 lines)

**Lines to remove**:
- Alpha check: lines 1462-1486 (25 lines)
- Beta check: lines 1488-1519 (32 lines)
- Total reduction: ~57 lines

**Lines to keep**:
- omega check: lines 1429-1460 (32 lines) ✅
- do_wK check: lines 1373-1408 (36 lines) ✅

---

## Final Verdict

| Parameter | Current Check | Correct? | Action |
|-----------|---------------|----------|--------|
| omega | MUST match | ✅ YES | Keep |
| do_wK | MUST match | ✅ YES | Keep |
| alpha | MUST match | ❌ NO | **Remove** |
| beta | MUST match | ❌ NO | **Remove** |

**Investigation complete!** ✅
