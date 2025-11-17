# Multi-SCF Input Design - Snapshot Pattern

## Проблема

**Текущий flow:**
1. User sets global options: `psi4.set_options({'diis_start': 0})`
2. wfn1 создается
3. Baseline test runs → меняет global options на `diis_start=14`
4. wfn2 создается
5. multi_scf([wfn1, wfn2]) calls scf_iterate() → **ЧИТАЕТ РАЗНЫЕ ОПЦИИ!**

**Результат:** wfn1.diis_start_=0, wfn2.diis_start_=14 → non-determinism

## Архитектурное решение

### 1. Snapshot Pattern (Freeze options ONCE per wfn)

```python
# File: psi4/driver/procrouting/scf_proc/scf_options_snapshot.py (NEW FILE)

def snapshot_scf_options():
    """
    Capture current global SCF options into a frozen dictionary.

    This allows each wavefunction to have its own independent copy
    of options, preventing pollution from later option changes.

    Returns:
        dict: Frozen snapshot of all SCF-relevant options
    """
    snapshot = {}

    # DIIS options
    snapshot['DIIS_START'] = core.get_option('SCF', 'DIIS_START')
    snapshot['DIIS_MIN_VECS'] = core.get_option('SCF', 'DIIS_MIN_VECS')
    snapshot['DIIS_MAX_VECS'] = core.get_option('SCF', 'DIIS_MAX_VECS')
    snapshot['DIIS_RMS_ERROR'] = core.get_option('SCF', 'DIIS_RMS_ERROR')

    # Convergence acceleration
    snapshot['MOM_START'] = core.get_option('SCF', 'MOM_START')
    snapshot['FRAC_START'] = core.get_option('SCF', 'FRAC_START')
    snapshot['FRAC_OCC'] = core.get_option('SCF', 'FRAC_OCC')
    snapshot['FRAC_VAL'] = core.get_option('SCF', 'FRAC_VAL')
    snapshot['DAMPING_PERCENTAGE'] = core.get_option('SCF', 'DAMPING_PERCENTAGE')
    snapshot['SOSCF_START_CONVERGENCE'] = core.get_option('SCF', 'SOSCF_START_CONVERGENCE')

    # Convergence thresholds
    snapshot['E_CONVERGENCE'] = core.get_option('SCF', 'E_CONVERGENCE')
    snapshot['D_CONVERGENCE'] = core.get_option('SCF', 'D_CONVERGENCE')
    snapshot['MAXITER'] = core.get_option('SCF', 'MAXITER')

    # Other SCF options
    snapshot['PRINT'] = core.get_option('SCF', 'PRINT')
    snapshot['REFERENCE'] = core.get_option('SCF', 'REFERENCE')

    return snapshot


def apply_options_snapshot(wfn, snapshot):
    """
    Apply frozen options snapshot to wavefunction.

    This sets wavefunction-local attributes that scf_iterate()
    will use INSTEAD of reading from global options.

    Args:
        wfn: Wavefunction object (RHF/UHF/ROHF)
        snapshot: dict from snapshot_scf_options()
    """
    # Store snapshot on wfn for later use
    wfn._options_snapshot = snapshot

    # Apply immediately (for backward compatibility)
    # These will be used by scf_iterate() if present
    wfn.diis_start_ = snapshot['DIIS_START']
    wfn.diis_min_vecs_ = snapshot['DIIS_MIN_VECS']
    wfn.diis_max_vecs_ = snapshot['DIIS_MAX_VECS']

    # ... etc for other options


def multi_scf_wavefunction_factory(name, molecule, reference='RHF', options=None):
    """
    Create wavefunction with frozen options for multi-SCF.

    This is the multi-SCF equivalent of scf_wavefunction_factory().
    It creates wfn with OPTIONS SNAPSHOT to prevent pollution.

    Args:
        name: Method name ('hf', 'b3lyp', etc)
        molecule: Molecule object
        reference: 'RHF', 'UHF', 'ROHF', 'CUHF'
        options: dict of options to override global (optional)

    Returns:
        Wavefunction with frozen options snapshot

    Example:
        # Create 2 wfn with different DIIS settings
        wfn1 = multi_scf_wavefunction_factory('hf', mol, options={'DIIS_START': 0})
        wfn2 = multi_scf_wavefunction_factory('hf', mol, options={'DIIS_START': 5})

        # Or use current global options for all
        wfn1 = multi_scf_wavefunction_factory('hf', mol1)
        wfn2 = multi_scf_wavefunction_factory('hf', mol2)
    """
    # Snapshot current global options BEFORE any changes
    if options is None:
        snapshot = snapshot_scf_options()
    else:
        # Merge user options with global
        snapshot = snapshot_scf_options()
        snapshot.update(options)

    # Create base wavefunction
    ref_wfn = core.Wavefunction.build(molecule, core.get_global_option('BASIS'))

    # Build functional
    from psi4.driver.procrouting.proc import build_functional_and_disp
    superfunc, _disp_functor = build_functional_and_disp(name, restricted=(reference in ["RKS", "RHF"]))

    # Create wavefunction
    core.prepare_options_for_module("SCF")
    if reference in ["RHF", "RKS"]:
        wfn = core.RHF(ref_wfn, superfunc)
    elif reference == "ROHF":
        wfn = core.ROHF(ref_wfn, superfunc)
    elif reference in ["UHF", "UKS"]:
        wfn = core.UHF(ref_wfn, superfunc)
    elif reference == "CUHF":
        wfn = core.CUHF(ref_wfn, superfunc)
    else:
        raise ValidationError(f"Unknown reference: {reference}")

    # Set auxiliary basis
    if core.get_global_option("SCF_TYPE") in ["DF", "MEM_DF", "DISK_DF"]:
        aux_basis = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SCF",
                                         core.get_option("SCF", "DF_BASIS_SCF"),
                                         "JKFIT", core.get_global_option('BASIS'),
                                         puream=wfn.basisset().has_puream())
        wfn.set_basisset("DF_BASIS_SCF", aux_basis)

    # CRITICAL: Apply options snapshot BEFORE initialize()
    apply_options_snapshot(wfn, snapshot)

    # Initialize wavefunction
    wfn.initialize()

    return wfn
```

### 2. User API (3 levels of convenience)

#### Level 1: Automatic (99% users) - API UNCHANGED

```python
# Existing code works without modification
psi4.set_options({'basis': 'cc-pvdz', 'scf_type': 'df'})
E = psi4.energy('scf')  # Internally: scf_helper() → multi_scf_helper([wfn]) with N=1
                         # User sees: No changes, same API, same performance
```

#### Level 2: Simple multi-SCF (same options for all)

```python
from psi4.driver.procrouting.scf_proc.multi_scf_helper import multi_scf_helper

# Define molecules
mol1 = psi4.geometry("H2O ...")
mol2 = psi4.geometry("H2O ...")

# Set common options
psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'df',
    'diis_start': 0,
    'e_convergence': 1e-8
})

# Create and run multi-SCF (options snapshotted automatically)
energies = multi_scf_helper(
    molecules=[mol1, mol2],
    method='hf'
)
```

#### Level 3: Advanced (different options per wfn)

```python
from psi4.driver.procrouting.scf_proc.multi_scf_helper import (
    multi_scf_wavefunction_factory,
    multi_scf
)

# Common base options
psi4.set_options({'basis': 'cc-pvdz', 'scf_type': 'df'})

# Create wfn with individual options
wfn1 = multi_scf_wavefunction_factory(
    'hf', mol1,
    options={'diis_start': 0, 'mom_start': 5}  # Aggressive
)

wfn2 = multi_scf_wavefunction_factory(
    'hf', mol2,
    options={'diis_start': 5, 'damping_percentage': 20}  # Conservative
)

# Run multi-SCF
energies = multi_scf([wfn1, wfn2])
```

### 3. Unified Architecture (CRITICAL DESIGN DECISION)

**IMPORTANT:** ALL SCF calculations go through multi-SCF coordinator!

- Single-cycle SCF → `scf_helper()` → `multi_scf_helper([wfn])` with N=1
- Multi-state SCF → `multi_scf_helper([wfn1, wfn2, ...])` with N>1
- **NO separate code paths** → single unified implementation
- User API unchanged: `psi4.energy('scf')` works out of box
- Internally: ALL use same coordinator with options snapshot pattern

## Implementation Plan

### Phase A: Infrastructure (unified architecture)

1. Create `psi4/driver/procrouting/scf_proc/scf_options_snapshot.py`
   - `snapshot_scf_options()` function
   - `apply_options_snapshot()` function

2. Create `psi4/driver/procrouting/scf_proc/multi_scf_helper.py`
   - `multi_scf_wavefunction_factory()` - creates wfn with snapshot
   - `multi_scf_helper()` - high-level API

3. Modify `scf_helper()` in `proc.py`
   - Always calls `multi_scf_helper([wfn])` for N=1 case
   - User API unchanged, internally unified

4. Modify `scf_iterator.py` `_scf_initialize_iteration_state()`
   - Read from `wfn._options_snapshot` (ALWAYS present)
   - No fallback to global options (snapshot always set by multi_scf_wavefunction_factory)

### Phase B: Testing

1. Test single-cycle SCF (N=1) works through unified coordinator
2. Test multi-SCF (N=2) with same options for all states
3. Test multi-SCF (N=2) with different options per state
4. Test that baseline pollution is eliminated (determinism across 100 runs)
5. Verify performance: single-cycle same speed, multi-state ~2x faster

### Phase C: Documentation

1. Update multi-SCF examples
2. Document snapshot pattern rationale
3. Add "Best Practices" guide

## Benefits

✅ **No changes for 99% users** - `psi4.energy('scf')` API unchanged
✅ **Unified architecture** - single code path for N=1 and N>1 (easier maintenance)
✅ **Eliminates global state pollution** - each wfn has frozen options
✅ **Flexible for advanced users** - can override per-state options
✅ **No code duplication** - single coordinator handles all cases
✅ **Forward compatible with SA-REKS** - snapshot per state
✅ **Testable** - can verify snapshot isolation and determinism
✅ **Performance** - single-cycle unchanged, multi-state ~2x faster

## File Structure

```
psi4/driver/procrouting/scf_proc/
├── scf_iterator.py              # Existing (minor mod to read snapshot)
├── scf_options_snapshot.py      # NEW - snapshot/apply functions
├── multi_scf_helper.py          # NEW - high-level multi-SCF API
└── __init__.py                  # Export new functions
```

## Next Steps

1. Implement `scf_options_snapshot.py` (30 min)
2. Implement `multi_scf_helper.py` (1 hour)
3. Modify `scf_iterator.py` to check snapshot (15 min)
4. Test with current `test_multi_scf.py` (30 min)
5. Verify determinism across 100 runs

Total: ~2.5 hours to implement complete solution
