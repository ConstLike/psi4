# Shared JK Compatibility Requirements

**Date**: 2025-11-18
**Purpose**: Define EXACT requirements for wavefunction compatibility in multi_scf()
**Method**: **FACTUAL code analysis** (not assumptions!)

---

## Methodology

This document is based on **tracing the actual code paths** to identify what parameters are used by shared JK. Only parameters that affect shared components are required to match.

**Principle**: Check ONLY what is SHARED, nothing more, nothing less.

---

## Code Analysis: Shared JK Creation

### _build_jk() - Line 96-101

```python
def _build_jk(wfn, memory):
    jk = core.JK.build(wfn.get_basisset("ORBITAL"),       # ← PARAMETER 1: Primary basis
                       aux=wfn.get_basisset("DF_BASIS_SCF"),  # ← PARAMETER 2: Auxiliary basis
                       do_wK=wfn.functional().is_x_lrc(),     # ← PARAMETER 3: LRC flag
                       memory=memory)
    return jk
```

**What's built into JK object**:
1. **Primary basis set** (`"ORBITAL"`) - defines AO basis functions
2. **Auxiliary basis set** (`"DF_BASIS_SCF"`) - for density fitting (if DF)
3. **LRC capability** (`do_wK`) - whether JK can compute long-range K

**Consequence**: These parameters are BAKED INTO the JK object structure. Cannot be changed after build!

### initialize_jk() - Lines 104-122

```python
def initialize_jk(self, memory, jk=None):
    functional = self.functional()
    if jk is None:
        jk = _build_jk(self, memory)

    self.set_jk(jk)

    jk.set_print(self.get_print())
    jk.set_memory(memory)
    jk.set_do_K(functional.is_x_hybrid())         # ← PARAMETER 4: Hybrid flag
    jk.set_do_wK(functional.is_x_lrc())           # ← PARAMETER 5: LRC flag (runtime)
    jk.set_omega(functional.x_omega())            # ← PARAMETER 6: LRC omega
    jk.set_omega_alpha(functional.x_alpha())      # ← PARAMETER 7: LRC alpha
    jk.set_omega_beta(functional.x_beta())        # ← PARAMETER 8: LRC beta

    jk.initialize()  # ← COMPUTES 3-INDEX INTEGRALS
    jk.print_header()
```

**What's configured in JK**:
4. **do_K** - whether to compute K matrices
5. **do_wK** - whether to compute long-range K (runtime check)
6. **omega** - LRC range-separation parameter
7. **omega_alpha** - RSH (range-separated hybrid) parameter
8. **omega_beta** - RSH parameter

**jk.initialize()** - what does it compute?
- **3-index integrals (Q|μν)** based on:
  - Primary basis (from build)
  - Auxiliary basis (from build, if DF)
  - **Geometry** (atomic coordinates!) ← CRITICAL DEPENDENCY

---

## Code Analysis: Shared JK Usage

### multi_scf() Shared JK Pattern - Lines 1367-1392

```python
if needs_jk_init:
    ref_wfn = wfn_list[0]  # ← Uses FIRST wavefunction as reference!

    # Step 1: Build JK from ref_wfn parameters
    shared_jk = _build_jk(ref_wfn, total_memory)
    # Uses: ref_wfn.get_basisset("ORBITAL")
    # Uses: ref_wfn.get_basisset("DF_BASIS_SCF")
    # Uses: ref_wfn.functional().is_x_lrc()

    # Step 2: Configure JK from ref_wfn functional
    ref_wfn.initialize_jk(total_memory, jk=shared_jk)
    # Sets: jk.omega = ref_wfn.functional().x_omega()
    # Sets: jk.omega_alpha = ref_wfn.functional().x_alpha()
    # Sets: jk.omega_beta = ref_wfn.functional().x_beta()
    # Computes: 3-index integrals for ref_wfn geometry!

    # Step 3: Share JK with other wavefunctions
    for wfn in wfn_list[1:]:
        wfn.set_jk(shared_jk)  # ← Just sets pointer!

# Step 4: Initialize all wavefunctions
for wfn in wfn_list:
    wfn.initialize()  # ← Calls scf_initialize()
```

### scf_initialize() Idempotency - Lines 146-148

```python
if isinstance(self.jk(), core.JK):  # ← JK already exists!
    core.print_out("\nRe-using passed JK object instead of rebuilding\n")
    jk = self.jk()
    initialize_jk_obj = False  # ← SKIPS initialize_jk() call!
else:
    initialize_jk_obj = True
    jk = _build_jk(self, total_memory)

# ... later (line 192-193) ...
if initialize_jk_obj:  # ← False for wfn_list[1:]!
    self.initialize_jk(self.memory_jk_, jk=jk)  # ← SKIPPED!
```

**CRITICAL OBSERVATION**:
- wfn_list[0]: Calls `initialize_jk()` → configures omega, alpha, beta
- wfn_list[1:]: **SKIPS** `initialize_jk()` → uses omega from wfn[0]!

**Consequence**: If wfn have different omega, wfn[1:] will use **WRONG omega**!

### _multi_scf_inner() JK Reconfiguration - Lines 1444-1445

```python
jk = wfn_list[0].jk()  # Shared JK

# Reconfigure for all wfn
jk.set_do_J(True)   # ← Always compute J
jk.set_do_K(True)   # ← Always compute K (even for pure GGA)
# NOTE: omega is NOT reconfigured here!
```

**Observation**:
- do_J, do_K are **overwritten** → safe if different
- omega, omega_alpha, omega_beta are **NOT overwritten** → MUST be same!

---

## MUST Match Requirements

Based on code analysis, these parameters MUST match for all wfn in multi_scf():

### Category 1: JK Structure (Baked into JK.build())

These cannot be changed after JK is built:

#### 1.1 Primary Basis Set ✅ CRITICAL

**Code location**: `_build_jk()` line 97
```python
wfn.get_basisset("ORBITAL")
```

**Why**: JK object is built with specific basis set dimensions (nbf, nshell).
3-index integrals (Q|μν) are computed for this specific basis.

**Check**:
```python
ref_basis = wfn_list[0].basisset()
for wfn in wfn_list[1:]:
    if wfn.basisset().name() != ref_basis.name():
        raise Error("Basis set mismatch")
    if wfn.basisset().nbf() != ref_basis.nbf():
        raise Error("Basis function count mismatch")
```

**Consequence if mismatch**: Segfault or garbage integrals!

#### 1.2 Geometry (Molecule) ✅ CRITICAL

**Code location**: `jk.initialize()` line 121 (computes integrals)

**Why**: 3-index integrals (Q|μν) depend on atomic positions:
```
(Q|μν) = ∫∫ χ_Q(r1) χ_μ(r1) r12^-1 χ_ν(r2) dr1 dr2
```
Basis functions χ are centered on atoms → integrals change if geometry changes!

**Check**:
```python
ref_mol = wfn_list[0].molecule()
for wfn in wfn_list[1:]:
    if wfn.molecule() is not ref_mol:  # Object identity
        # Must be SAME molecule object!
        raise Error("Different molecule objects")
    # Could also check geometry hash if molecules are cloned
```

**Consequence if mismatch**: Wrong integrals, wrong energy, wrong everything!

#### 1.3 Auxiliary Basis (for DF only) ✅ CRITICAL if SCF_TYPE='DF'

**Code location**: `_build_jk()` line 98
```python
aux=wfn.get_basisset("DF_BASIS_SCF")
```

**Why**: For density fitting, JK uses auxiliary basis to approximate 4-index integrals:
```
(μν|ρσ) ≈ (μν|P) (P|Q)^-1 (Q|ρσ)
```
Different auxiliary basis → different approximation!

**Check**:
```python
if core.get_global_option('SCF_TYPE') == 'DF':
    ref_aux = wfn_list[0].get_basisset("DF_BASIS_SCF")
    for wfn in wfn_list[1:]:
        aux = wfn.get_basisset("DF_BASIS_SCF")
        if aux.name() != ref_aux.name():
            raise Error("DF auxiliary basis mismatch")
```

**Consequence if mismatch**: Different JK results for different wfn!

#### 1.4 LRC Capability ✅ CRITICAL

**Code location**: `_build_jk()` line 99
```python
do_wK=wfn.functional().is_x_lrc()
```

**Why**: JK object built with or without long-range K capability.
If built with `do_wK=False`, cannot compute wK later!

**Check**:
```python
ref_is_lrc = wfn_list[0].functional().is_x_lrc()
for wfn in wfn_list[1:]:
    if wfn.functional().is_x_lrc() != ref_is_lrc:
        raise Error("LRC capability mismatch: "
                    "all wfn must be LRC or all non-LRC")
```

**Consequence if mismatch**:
- If ref is non-LRC but wfn[i] is LRC → wK not computed → wrong energy!
- If ref is LRC but wfn[i] is non-LRC → overhead (wK computed but not used)

### Category 2: JK Configuration (Set in initialize_jk())

These are configured during initialize_jk() and NOT reconfigured later:

#### 2.1 LRC Omega Parameter ✅ CRITICAL if is_x_lrc()

**Code location**: `initialize_jk()` line 116
```python
jk.set_omega(functional.x_omega())
```

**Why**: Range-separation parameter for LRC functionals:
```
1/r = erf(ω r)/r + erfc(ω r)/r
      \___ LR ___/   \____ SR ____/
```
Different omega → different long-range/short-range splitting!

**CRITICAL**: Only set during ref_wfn.initialize_jk(), NOT updated for other wfn!

**Check**:
```python
if ref_is_lrc:
    ref_omega = wfn_list[0].functional().x_omega()
    for wfn in wfn_list[1:]:
        if abs(wfn.functional().x_omega() - ref_omega) > 1e-10:
            raise Error(f"LRC omega mismatch: {wfn.functional().x_omega()} "
                        f"vs {ref_omega}")
```

**Consequence if mismatch**: wfn[1:] use WRONG omega → wrong energy!

#### 2.2 RSH Alpha Parameter ✅ CRITICAL if is_x_lrc()

**Code location**: `initialize_jk()` line 118
```python
jk.set_omega_alpha(functional.x_alpha())
```

**Why**: Range-separated hybrid (RSH) parameter:
```
E_X = α E_X^SR,HF + β E_X^LR,HF + (1-α) E_X^SR,DFT
```

**Check**:
```python
if ref_is_lrc:
    ref_alpha = wfn_list[0].functional().x_alpha()
    for wfn in wfn_list[1:]:
        if abs(wfn.functional().x_alpha() - ref_alpha) > 1e-10:
            raise Error("RSH alpha mismatch")
```

#### 2.3 RSH Beta Parameter ✅ CRITICAL if is_x_lrc()

**Code location**: `initialize_jk()` line 119
```python
jk.set_omega_beta(functional.x_beta())
```

**Why**: Long-range HF exchange fraction.

**Check**: Same as alpha.

### Category 3: Runtime Configuration (Set in _multi_scf_inner())

These are reconfigured in _multi_scf_inner(), so differences are OK:

#### 3.1 do_K Flag ✅ CAN DIFFER (overwritten)

**Code location**: `_multi_scf_inner()` line 1445
```python
jk.set_do_K(True)  # ← Always True, regardless of functional!
```

**Why**: Multi-SCF always computes K (even for pure GGA where it's not used).

**Check**: NOT NEEDED - will be overwritten to True.

**Note**: If one wfn is pure GGA (do_K=False) and another is hybrid (do_K=True),
the shared JK will compute K for both (slight overhead for GGA, but safe).

#### 3.2 do_J Flag ✅ CAN DIFFER (overwritten)

**Code location**: `_multi_scf_inner()` line 1444
```python
jk.set_do_J(True)  # ← Always True
```

**Why**: All functionals need J.

**Check**: NOT NEEDED.

### Category 4: Algorithm Selection (Protected by options_snapshot)

#### 4.1 SCF_TYPE ✅ Already Protected

**Code location**: `multi_scf()` lines 1305-1312
```python
options_snapshot = snapshot_scf_options()
for wfn in wfn_list:
    apply_options_snapshot(wfn, options_snapshot)
```

**Why**: Options snapshot ensures all wfn read same SCF_TYPE.

**Check**: NOT NEEDED - already protected by snapshot.

**Paranoid check** (optional):
```python
ref_scf_type = core.get_global_option('SCF_TYPE')
# All wfn will use this from snapshot
```

---

## CAN Differ (Safe to be Different)

These parameters do NOT affect shared JK and can differ between wfn:

### 1. Multiplicity / Occupation ✅ OK

**Why**: Different occupation → different C matrices → different densities.
But JK computes J/K for ANY density, doesn't care about occupation.

**Example**: RHF (singlet) + UHF (triplet) → OK!

### 2. Reference Type ✅ OK

**Why**: RHF (n_states=1), UHF (n_states=2), ROHF (n_states=2) return different
number of C matrices, but JK processes any list of C matrices.

**Example**: RHF + UHF + ROHF → OK! (verified in SHARED_JK_VERIFICATION_RHF_UHF_ROHF.md)

### 3. Charge ✅ OK

**Why**: Different charge → different occupation, but doesn't affect integrals.

**Example**: Neutral + cation → OK!

### 4. Non-LRC XC Functional ✅ OK (if all non-LRC)

**Why**: For non-LRC functionals, only J and K differ. XC contribution is computed
separately per-wfn, doesn't go through shared JK.

**Example**:
- B3LYP + PBE0 → OK (both hybrid, non-LRC, different XC but same J/K needs)
- HF + B3LYP → OK (both need K, HF has 100% HF exchange, B3LYP has 20%)

**Constraint**: Both must have is_x_lrc() = False!

### 5. Hybrid Fraction ✅ OK

**Why**: do_K is overwritten to True in _multi_scf_inner().
Different HF exchange fractions don't affect J/K computation, only how they're combined.

**Example**: HF (100% HF) + B3LYP (20% HF) → OK!

### 6. Convergence Settings ✅ OK

**Why**: DIIS, damping, SOSCF, etc. are per-wfn settings, don't affect shared JK.

**Example**: wfn1 with DIIS_START=1, wfn2 with DIIS_START=5 → OK!

### 7. Per-Wfn Memory Settings ✅ OK

**Why**: Total memory for JK is set once, per-wfn settings don't affect shared JK.

---

## Summary: Validation Checklist

### MUST Match (CRITICAL - will cause wrong results or crashes):

| Parameter | Check Method | If DF Only? | If LRC Only? |
|-----------|--------------|-------------|--------------|
| **Primary basis** | `wfn.basisset().name()` | Always | Always |
| **Geometry** | `wfn.molecule() is ref_mol` | Always | Always |
| **Auxiliary basis** | `wfn.get_basisset("DF_BASIS_SCF").name()` | ✅ Yes | No |
| **LRC capability** | `wfn.functional().is_x_lrc()` | No | Always |
| **LRC omega** | `wfn.functional().x_omega()` | No | ✅ Yes |
| **RSH alpha** | `wfn.functional().x_alpha()` | No | ✅ Yes |
| **RSH beta** | `wfn.functional().x_beta()` | No | ✅ Yes |

### CAN Differ (Safe):

- Multiplicity
- Reference type (RHF/UHF/ROHF)
- Charge
- Non-LRC XC functional
- Hybrid fraction (if do_K reconfigured, which we do)
- Convergence settings
- Per-wfn memory

### Already Protected (No Check Needed):

- SCF_TYPE (options snapshot)
- do_J, do_K (overwritten in _multi_scf_inner)

---

## Implementation Priority

**HIGH PRIORITY** (2-3 hours):
1. Implement validation function with these checks
2. Call at start of multi_scf() before any initialization
3. Clear error messages indicating WHAT doesn't match and WHY it matters

**Example error message**:
```
ValidationError: Wavefunction 2 has different basis set:
  Expected: cc-pVDZ (from wavefunction 0)
  Got: aug-cc-pVDZ

  Reason: All wavefunctions must use the same primary basis set
  because they share a single JK object with 3-index integrals
  computed for a specific basis. Different basis sets would
  require different integrals.

  Solution: Use the same basis set for all wavefunctions in multi_scf().
```

---

## Conclusion

This analysis is based on **ACTUAL code tracing**, not assumptions.

**Key Insights**:
1. Shared JK is built from wfn[0] parameters → all wfn must match those parameters
2. wfn[1:] skip initialize_jk() due to idempotency → configuration from wfn[0] is used
3. Some parameters (do_J, do_K) are reconfigured → safe if different
4. Most critical: basis, geometry, LRC parameters (omega)

**No half-measures**: Every MUST requirement has factual code justification! 🎯
