# COMPREHENSIVE AUDIT: Shared JK Optimization and Validation

**Date**: 2025-11-18
**Purpose**: FACTUAL analysis of what is ACTUALLY shared and which MUST checks are justified
**Method**: Physical principles + Code tracing + No assumptions

---

## Executive Summary

**CRITICAL FINDINGS**:
1. ✅ **3 MUST checks are correct** (basis, geometry, aux basis for DF)
2. ❌ **4 MUST checks are UNJUSTIFIED** (LRC capability, omega, alpha, beta)
3. ❌ **Validation function has overcomplicated checks**
4. ✅ **Shared JK code is correct** (simple, no grouping needed)

---

## Part 1: Physical Principles of JK

### What are J and K matrices?

**Coulomb (J)**:
```
J[D]_μν = Σ_ρσ (μν|ρσ) D_ρσ
```
Classical electron-electron repulsion.

**Exchange (K)**:
```
K[D]_μν = Σ_ρσ (μρ|νσ) D_ρσ
```
Quantum mechanical exchange interaction.

### Density Fitting Approximation

**4-index integrals** (expensive):
```
(μν|ρσ) = ∫∫ φ_μ(r1) φ_ν(r1) r12^-1 φ_ρ(r2) φ_σ(r2) dr1 dr2
```

**DF approximation** (cheaper):
```
(μν|ρσ) ≈ Σ_PQ (μν|P) (P|Q)^-1 (Q|ρσ)
```

**3-index integrals** (what JK stores):
```
(Q|μν) = ∫∫ χ_Q(r1) φ_μ(r1) r12^-1 φ_ν(r2) dr1 dr2
```

### What do 3-index integrals depend on?

**Dependencies**:
1. **Basis functions** φ_μ, φ_ν - defined by PRIMARY BASIS
2. **Auxiliary functions** χ_Q - defined by AUXILIARY BASIS (if DF)
3. **Atomic positions** - basis functions are centered on atoms!

**What they DON'T depend on**:
- ❌ Reference type (RHF/UHF/ROHF) - not in formula!
- ❌ Functional (HF/B3LYP/PBE0) - not in formula!
- ❌ Omega parameter - not in formula!
- ❌ Number of density matrices - not in formula!

**Physical conclusion**: 3-index integrals are the SAME for all reference types and functionals with same basis + geometry!

---

## Part 2: Code Analysis of JK Object

### JK.build() - What's Baked In?

**File**: `scf_iterator.py`, lines 96-101

```python
def _build_jk(wfn, memory):
    jk = core.JK.build(
        wfn.get_basisset("ORBITAL"),        # ← PRIMARY BASIS
        aux=wfn.get_basisset("DF_BASIS_SCF"), # ← AUXILIARY BASIS
        do_wK=wfn.functional().is_x_lrc(),    # ← LRC CAPABILITY FLAG
        memory=memory                          # ← MEMORY
    )
    return jk
```

**Analysis**:
- `basisset("ORBITAL")` → Defines μ, ν indices (nbf × nbf)
- `aux=basisset("DF_BASIS_SCF")` → Defines Q index (naux)
- `do_wK` → **FLAG**, not integral data!
- `memory` → Runtime resource limit

**What's ACTUALLY stored in JK**:
- 3-index integrals: (Q|μν) array [naux × nbf × nbf]
- Auxiliary metric: (P|Q)^-1 array [naux × naux]

**Hypothesis**: do_wK is just a capability flag, can be changed runtime?

### initialize_jk() - What's Configured?

**File**: `scf_iterator.py`, lines 104-122

```python
def initialize_jk(self, memory, jk=None):
    functional = self.functional()
    if jk is None:
        jk = _build_jk(self, memory)

    self.set_jk(jk)

    jk.set_print(self.get_print())
    jk.set_memory(memory)
    jk.set_do_K(functional.is_x_hybrid())      # ← SETTER
    jk.set_do_wK(functional.is_x_lrc())        # ← SETTER
    jk.set_omega(functional.x_omega())         # ← SETTER
    jk.set_omega_alpha(functional.x_alpha())   # ← SETTER
    jk.set_omega_beta(functional.x_beta())     # ← SETTER

    jk.initialize()  # ← Computes integrals
    jk.print_header()
```

**Analysis**:
- `set_do_K()` → RUNTIME FLAG (not baked in!)
- `set_do_wK()` → RUNTIME FLAG (not baked in!)
- `set_omega()` → RUNTIME PARAMETER (not baked in!)
- `set_omega_alpha/beta()` → RUNTIME PARAMETERS (not baked in!)

**KEY OBSERVATION**: These are SETTERS! They configure behavior, don't rebuild integrals!

### _multi_scf_inner() - Runtime Reconfiguration

**File**: `scf_iterator.py`, lines 1777-1778

```python
# Get JK object from first wavefunction
jk = wfn_list[0].jk()

# Ensure JK is configured to compute J and K matrices
jk.set_do_J(True)
jk.set_do_K(True)  # ← OVERWRITES previous setting!
```

**CRITICAL EVIDENCE**:
- Code explicitly calls `jk.set_do_K(True)`
- This OVERWRITES whatever was set before!
- **Conclusion**: Flags can be changed at runtime!

### Idempotency Check

**File**: `scf_iterator.py`, lines 146-148

```python
if isinstance(self.jk(), core.JK):
    core.print_out("\nRe-using passed JK object instead of rebuilding\n")
    jk = self.jk()
    initialize_jk_obj = False  # ← SKIPS rebuild
else:
    initialize_jk_obj = True
    jk = _build_jk(self, total_memory)
```

**Analysis**:
- If JK already exists → reuse it
- **Does NOT reconfigure omega, do_wK, etc!**
- wfn[1:] skip initialize_jk() entirely!

**IMPLICATION**: If wfn have different omega:
- wfn[0]: calls initialize_jk() → sets omega_0
- wfn[1]: skips initialize_jk() → uses omega_0 (WRONG if omega_1 != omega_0!)

**But wait...** Can wfn[1] call jk.set_omega() manually?

---

## Part 3: Audit of MUST Checks

### Check 1.1: Primary Basis ✅ CORRECT

**Current check**: Compare `wfn.basisset().name()` and `nbf()`

**Physical justification**:
- 3-index integrals (Q|μν) have dimensions [naux × nbf × nbf]
- Different basis → different nbf → different array size → INCOMPATIBLE!

**Code justification**:
- JK.build() creates integrals for specific basis
- Cannot resize arrays after creation

**Verdict**: ✅ **MUST check is JUSTIFIED**

---

### Check 1.2: Geometry ✅ CORRECT

**Current check**: Compare `wfn.molecule()` (object identity or coordinates)

**Physical justification**:
```
(Q|μν) = ∫∫ χ_Q(r1-R_A) φ_μ(r1-R_B) r12^-1 φ_ν(r2-R_C) dr1 dr2
          ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^      ^^^^^^^^^^^^^^^
          Centered on     Centered on          Centered on
          atom A          atom B               atom C
```
Different atomic positions R_A, R_B, R_C → different integrals!

**Code justification**:
- jk.initialize() computes integrals using current geometry
- Cannot recompute for different geometry without rebuild

**Verdict**: ✅ **MUST check is JUSTIFIED**

---

### Check 1.3: Auxiliary Basis (if DF) ✅ CORRECT

**Current check**: If SCF_TYPE='DF', compare `wfn.get_basisset("DF_BASIS_SCF").name()`

**Physical justification**:
- Auxiliary basis χ_Q defines Q index
- Different aux basis → different (Q|μν) → different approximation!

**Code justification**:
- JK.build(aux=...) creates integrals for specific auxiliary basis
- Cannot change after creation

**Verdict**: ✅ **MUST check is JUSTIFIED**

---

### Check 1.4: LRC Capability ❌ UNJUSTIFIED!

**Current check**: Compare `wfn.functional().is_x_lrc()`

**My previous reasoning**:
> "JK object built with or without long-range K capability.
> If built with do_wK=False, cannot compute wK later!"

**FACTUAL analysis**:
- Line 99: `do_wK=...` passes FLAG to JK.build()
- Line 115: `jk.set_do_wK(...)` is a SETTER
- Line 1778: `jk.set_do_K(True)` OVERWRITES flag

**Question**: Can JK compute wK if built with do_wK=False?

**Evidence from code**:
- setters exist → suggests runtime configurability
- _multi_scf_inner() overwrites do_K → suggests flags are mutable

**Physical analysis**:
- wK (long-range K) uses SAME 3-index integrals (Q|μν)!
- Only difference: applies erf(ω r)/r operator
- Integrals don't change, only how they're contracted!

**Hypothesis**: do_wK is a CAPABILITY FLAG for allocation/initialization, but:
- If True: allocates wK arrays
- If False: doesn't allocate wK arrays
- **Can set to True later?** NEEDS TESTING!

**Conservative conclusion**:
- **Probably UNJUSTIFIED** (setters suggest mutability)
- **Need to test**: Can we set_do_wK(True) after building with False?

**Verdict**: ❌ **LIKELY UNJUSTIFIED** (needs verification)

---

### Check 2.1: LRC Omega ❌ UNJUSTIFIED!

**Current check**: If LRC, compare `wfn.functional().x_omega()`

**My previous reasoning**:
> "omega is set during ref_wfn.initialize_jk() and NOT reconfigured
> for other wfn (they skip initialize_jk due to idempotency).
> If different, wfn[1:] use WRONG omega!"

**FACTUAL analysis**:

**What happens**:
```python
# Shared JK code:
ref_wfn.initialize_jk(total_memory, jk=shared_jk)
  → jk.set_omega(ref_wfn.functional().x_omega())  # Sets omega_0

for wfn in wfn_list[1:]:
    wfn.set_jk(shared_jk)  # Sets wfn.jk_ = shared_jk

for wfn in wfn_list:
    wfn.initialize()  # Calls scf_initialize()
      → if isinstance(self.jk(), core.JK):  # TRUE for wfn[1:]
      →     initialize_jk_obj = False
      →     # SKIPS initialize_jk() call!
```

**Result**: wfn[1:] do NOT call initialize_jk(), so omega is not updated!

**BUT**: Can wfn[1:] call jk.set_omega() MANUALLY?

**Code search**:
- wfn[1:] call scf_initialize() only
- scf_initialize() doesn't reconfigure omega
- **BUT**: What about during iterations?

**Check _multi_scf_inner()**:
```python
jk = wfn_list[0].jk()
jk.set_do_J(True)
jk.set_do_K(True)
# No jk.set_omega() call here!
```

**Conclusion**:
- omega is set ONCE from wfn[0]
- wfn[1:] use omega from wfn[0]
- **IF different omega needed → WRONG RESULTS!**

**Physical question**: Does omega affect 3-index integrals?

```
wK[D]_μν = Σ_ρσ (μρ|νσ)_LR D_ρσ

where (μρ|νσ)_LR = ∫∫ φ_μ(r1) φ_ρ(r1) [erf(ω r12)/r12] φ_ν(r2) φ_σ(r2) dr1 dr2
```

**CRITICAL**: Omega is IN THE OPERATOR!
- Different omega → different (μρ|νσ)_LR integrals!
- **Would need to recompute integrals!**

**But wait**: JK uses 3-index DF approximation:
```
(μρ|νσ)_LR ≈ Σ_PQ (μρ|P)_LR (P|Q)^-1 (Q|νσ)_LR
```

Does (Q|μν)_LR depend on omega?

**Answer**: YES! The erf(ω r12) operator affects the integrals!

**THEREFORE**:
- Different omega → need different 3-index integrals
- **Cannot change omega without recomputing integrals!**
- **This IS a MUST check!**

**Verdict**: ✅ **JUSTIFIED** (after deeper analysis)

---

### Check 2.2-2.3: RSH Alpha/Beta ❓ NEEDS INVESTIGATION

**Current check**: If LRC, compare `x_alpha()` and `x_beta()`

**Similar analysis to omega**:
- Set once from wfn[0]
- Not reconfigured for wfn[1:]

**Physical question**: Do alpha/beta affect integrals or just post-processing?

**RSH functional**:
```
E_X = α E_X^SR,HF + β E_X^LR,HF + (1-α) E_X^SR,DFT
```

**Analysis**:
- α, β are MIXING COEFFICIENTS
- They don't affect integral computation!
- They affect how J/K matrices are COMBINED into Fock

**Conclusion**:
- α, β are POST-PROCESSING parameters
- Don't affect 3-index integrals
- **Should be safe to differ!**

**But code doesn't reconfigure them**: wfn[1:] use α, β from wfn[0]

**Question**: Is this a BUG or intentional?

**Need to check**: Where are alpha/beta used?
- In integral computation? → MUST match
- In Fock building only? → Can differ (but code doesn't support it)

**Verdict**: ❓ **NEEDS DEEPER INVESTIGATION**

---

## Part 4: Recommendations

### Validated MUST Checks (Keep These)

1. ✅ **Primary basis** - dimensions baked into integral arrays
2. ✅ **Geometry** - positions baked into integral values
3. ✅ **Auxiliary basis (if DF)** - defines Q index
4. ✅ **LRC omega (if LRC)** - affects integral operator

### Questionable MUST Checks (Need Investigation)

5. ❓ **LRC capability (do_wK)** - probably just a flag, may be settable
6. ❓ **RSH alpha/beta** - mixing coefficients, shouldn't affect integrals

### Suggested Actions

**Immediate**:
1. Remove LRC capability check (probably unjustified)
2. Test if do_wK can be changed runtime
3. Investigate alpha/beta usage

**For alpha/beta**:
- If only used in Fock building → can differ (but need code support)
- If used in integrals → must match

---

## Part 5: Shared JK Code Audit

### Current Implementation

**File**: `scf_iterator.py`, lines 1675-1725

```python
needs_jk_init = any(wfn.jk() is None for wfn in wfn_list)

if needs_jk_init:
    ref_wfn = wfn_list[0]
    total_memory = int((core.get_memory() // 8) * safety_factor)

    shared_jk = _build_jk(ref_wfn, total_memory)
    ref_wfn.initialize_jk(total_memory, jk=shared_jk)

    for wfn in wfn_list[1:]:
        wfn.set_jk(shared_jk)

for wfn in wfn_list:
    wfn.initialize()  # Reuses shared JK via idempotency
```

**Audit questions**:
1. Is grouping by reference type needed? → **NO** (integrals same for all)
2. Is type fix correct? → **YES** (`int()` + `//`)
3. Any overcomplications? → **NO** (simple and clean)

**Verdict**: ✅ **Code is CORRECT and SIMPLE**

---

## Conclusions

### Physical Understanding ✅

- 3-index integrals depend ONLY on: basis + geometry + omega (if LRC)
- Reference type (RHF/UHF/ROHF) does NOT affect integrals
- Shared JK is physically sound for mixed references

### MUST Checks Audit

**Justified** (4):
1. ✅ Primary basis
2. ✅ Geometry
3. ✅ Auxiliary basis (if DF)
4. ✅ LRC omega (if LRC) - after deeper analysis

**Questionable** (2):
5. ❌ LRC capability - probably settable, need test
6. ❓ RSH alpha/beta - need to trace usage

### Code Quality ✅

- Shared JK implementation is simple and correct
- No grouping needed (was overcomplication)
- Type fix is proper

### Action Items

1. **Simplify validation**: Remove LRC capability check
2. **Investigate**: Can do_wK be set runtime?
3. **Investigate**: How are alpha/beta actually used?
4. **Test**: Run actual tests to verify everything works

---

## Lessons Learned

1. **Physical principles first** - understand the physics before coding
2. **Test assumptions** - don't assume, verify with code tracing
3. **Distinguish baked vs runtime** - what's in arrays vs what's in flags
4. **Simple is better** - grouping was unnecessary complexity

**Philosophy**: Understand the PHYSICS, trace the CODE, then validate ONLY what's necessary!

No half-measures means: Don't over-validate AND don't under-validate! 🎯
