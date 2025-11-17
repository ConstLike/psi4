# Single Source of Truth: Unified SCF Architecture

## Motivation

**Problem**: До этого изменения в Psi4 было **два параллельных пути** для SCF:

1. **Single SCF**: `psi4.energy('hf')` → `scf_compute_energy()` → `iterations()`
2. **Multi-SCF**: `multi_scf([wfn1, wfn2])` → `_multi_scf_inner()`

Это создавало риски:
- ❌ Дублирование кода (DF_SCF_GUESS, DIIS, damping, etc.)
- ❌ Соблазн развивать один путь и забывать про другой
- ❌ Разное поведение для single vs multi-SCF
- ❌ Двойная работа при добавлении новых фич

## Solution: Single Source of Truth

**Архитектурное решение**: Single SCF = частный случай multi-SCF

```python
# Before (two separate paths)
def scf_compute_energy(self):
    if DF_SCF_GUESS and DIRECT:
        # DF iterations
        self.initialize()
        self.iterations()
        # DIRECT iterations
        self.initialize_jk()
    else:
        self.initialize()
    self.iterations()  # Main loop
    return self.finalize_energy()

def multi_scf(wfn_list):
    if DF_SCF_GUESS and DIRECT:
        _multi_scf_inner(wfn_list)  # DF
        _multi_scf_inner(wfn_list)  # DIRECT
    else:
        _multi_scf_inner(wfn_list)
    # Duplicate logic!

# After (single source of truth)
def scf_compute_energy(self):
    """Single SCF is just multi_scf([self])"""
    energies = multi_scf([self])
    return self.finalize_energy()

def multi_scf(wfn_list):
    """Universal SCF coordinator"""
    if DF_SCF_GUESS and DIRECT:
        _multi_scf_inner(wfn_list)  # DF
        _multi_scf_inner(wfn_list)  # DIRECT
    else:
        _multi_scf_inner(wfn_list)
    return energies
```

## Architecture Flow

```
psi4.energy('hf')
    ↓
scf_helper()
    ↓
wfn.compute_energy()
    ↓
scf_compute_energy(self)  ← WRAPPER
    ↓
multi_scf([self])  ← SINGLE SOURCE OF TRUTH
    ↓
_multi_scf_inner([self])
    ↓
wfn._scf_iteration() × N iterations
```

**Key insight**: `psi4.energy('hf')` теперь использует `multi_scf([wfn])` с одной wfn!

## Benefits

### ✅ Code Deduplication
- **Before**: DF_SCF_GUESS implemented twice (scf_compute_energy + multi_scf)
- **After**: DF_SCF_GUESS implemented ONCE (multi_scf only)

### ✅ Automatic Feature Propagation
Любая фича добавленная в `multi_scf()` автоматически работает в single SCF:
- DIIS
- Damping
- SOSCF
- MOM
- FRAC
- DF_SCF_GUESS
- Convergence acceleration
- Future features!

### ✅ Identical Behavior
Single и multi-SCF теперь **гарантированно** используют одинаковую логику:
- Same iteration counts
- Same convergence behavior
- Same energy
- Same ALL features

### ✅ Maintainability
- **One function to maintain**: `multi_scf()` + `_multi_scf_inner()`
- **No code duplication**: Features implemented once
- **No divergence risk**: Can't forget to update one path

## Implementation Details

### scf_compute_energy() (wrapper)
```python
def scf_compute_energy(self):
    """
    Single source of truth: All SCF calculations flow through multi_scf().
    """
    self.iteration_energies = []  # Backward compatibility

    try:
        energies = multi_scf([self], verbose=True)
        scf_energy = energies[0]
    except SCFConvergenceError as e:
        if core.get_option("SCF", "FAIL_ON_MAXITER"):
            raise e
        else:
            scf_energy = self.get_energies("Total Energy")

    return self.finalize_energy()
```

### multi_scf() (universal coordinator)
```python
def multi_scf(wfn_list, e_conv=None, d_conv=None, max_iter=None, verbose=True):
    """
    Universal SCF coordinator for 1 to N wavefunctions.

    Handles:
    - Single SCF (len(wfn_list) == 1)
    - Multi-SCF (len(wfn_list) > 1)
    - DF_SCF_GUESS for DIRECT
    - All convergence features
    """
    # ... validation ...

    if DF_SCF_GUESS and SCF_TYPE == 'DIRECT':
        # Phase 1: DF pre-iterations
        _multi_scf_inner(wfn_list)
        # Phase 2: DIRECT final iterations
        _multi_scf_inner(wfn_list)
    else:
        _multi_scf_inner(wfn_list)

    return energies
```

### _multi_scf_inner() (iteration loop)
```python
def _multi_scf_inner(wfn_list, e_conv, d_conv, max_iter, verbose):
    """
    Inner SCF iteration loop (works for 1 to N wfn).
    """
    # Initialize
    jk = wfn_list[0].jk()

    for iteration in range(max_iter):
        # Collect C matrices from ALL wfn
        all_C = [wfn.get_orbital_matrices() for wfn in wfn_list]

        # Shared JK computation
        jk.compute()

        # Distribute J/K to each wfn
        for wfn in wfn_list:
            wfn.set_jk_matrices(J, K, wK)
            wfn._scf_iteration()

        # Check convergence
        if all_converged:
            break

    return energies
```

## Testing

### Before: Two separate tests needed
```python
def test_single_scf_df_guess():
    """Test DF_SCF_GUESS for single SCF"""
    # ...

def test_multi_scf_df_guess():
    """Test DF_SCF_GUESS for multi-SCF"""
    # ...
```

### After: One test covers both!
```python
def test_df_scf_guess():
    """Test DF_SCF_GUESS (works for both single and multi)"""
    # Single SCF
    e_single = psi4.energy('hf')  # Uses multi_scf([wfn])

    # Multi-SCF
    energies = multi_scf([wfn1, wfn2])

    # Both use same code path!
```

## Performance

**No overhead** for single SCF:
- `multi_scf([wfn])` with one wfn has zero overhead vs old path
- Same number of iterations
- Same JK computation
- Same convergence

**Benefit** for multi-SCF:
- All single-SCF optimizations automatically apply!

## Backward Compatibility

✅ **100% backward compatible**:
- `psi4.energy('hf')` still works identically
- `wfn.compute_energy()` still works
- All tests pass
- No API changes

## Future Work

Now that we have single source of truth, future improvements are **automatic**:

1. **Threading**: Parallelize wfn._scf_iteration() → works for single AND multi
2. **New convergence algorithms**: Implement in multi_scf → works everywhere
3. **Performance optimizations**: Shared JK pre-initialization → works for all
4. **New features**: MOM improvements, better DIIS → automatic propagation

## Philosophy

> "Make it work, make it right, make it fast"
> - Kent Beck

**This change is "make it right"**:
- ✅ Eliminates code duplication
- ✅ Ensures correctness through single implementation
- ✅ Makes future development easier
- ✅ Reduces maintenance burden

**No half-measures**: We don't maintain two parallel SCF implementations!

## Conclusion

**Before**:
```
Single SCF path ──────┐
                      ├─→ Features (DIIS, damping, etc.)
Multi-SCF path  ──────┘
    ↑ Risk of divergence!
```

**After**:
```
Single SCF (1 wfn) ──┐
                     ├─→ multi_scf() ─→ Features (DIIS, damping, etc.)
Multi-SCF (N wfn) ───┘
    ↑ Single source of truth!
```

This is the **right architecture** for long-term maintainability! 🎯
