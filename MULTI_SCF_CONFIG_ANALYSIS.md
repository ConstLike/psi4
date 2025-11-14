# Multi-SCF Configuration Analysis

## Концепция

Для multi-SCF нужно поддерживать:
1. **Одновременный запуск РАЗНЫХ типов**: RHF + ROHF, RHF + UHF, etc
2. **Одновременный запуск с РАЗНЫМИ опциями**: RHF с MOM + RHF без MOM, разные DIIS settings, etc
3. **Общие ресурсы**: shared JK computation

## Классификация параметров

### ОБЯЗАТЕЛЬНО одинаковые (для shared JK)

| Параметр | Почему | Проверка |
|----------|--------|----------|
| `BASIS` | JK создается для одного basis set | `wfn.basisset()` |
| `SCF_TYPE` | JK тип (DF/DIRECT/CD/etc) | `core.get_option('SCF', 'SCF_TYPE')` |
| `DF_BASIS_SCF` | Auxiliary basis для DF-JK | Если SCF_TYPE=DF |
| JK cutoff | Screening threshold | `jk.cutoff()` |
| omega | Для LRC functionals | `jk.omega()` (если используется) |

**Validation**: multi_scf() должен проверить что все wfn compatible.

### МОЖЕТ быть разным (per-wavefunction)

| Параметр | Где | Категория |
|----------|-----|-----------|
| **Type** | RHF/UHF/ROHF/CUHF | Wavefunction type |
| **Functional** | HF vs B3LYP vs PBE | Exchange-correlation |
| **DIIS settings** | | Convergence acceleration |
| - `DIIS_START` | Line 321 | When to start DIIS |
| - `DIIS_MIN_VECS` | DIIS code | Min vectors for extrapolation |
| - `DIIS_MAX_VECS` | DIIS code | Max vectors in buffer |
| **MOM** | Line 320 | Maximum Overlap Method |
| **Damping** | Line 322 | Density damping |
| **SOSCF** | Line 323 | Second-order SCF |
| **FRAC** | Line 324 | Fractional occupation |
| **Level shift** | Line 332 | Virtual orbital shifting |
| **Convergence** | | |
| - `E_CONVERGENCE` | Energy threshold | Per-wfn |
| - `D_CONVERGENCE` | Density threshold | Per-wfn |
| - `MAXITER` | Max iterations | Per-wfn |
| **Guess** | | Initial orbitals |
| - `GUESS` | SAD/CORE/READ/etc | Can differ |
| - `GUESS_MIX` | Line 334 | Orbital mixing |
| **Occupation** | | |
| - `nalpha` / `nbeta` | Electron count | Different molecules! |
| - `multiplicity` | Spin state | Singlet vs triplet |

### Зависит от контекста (может быть общим или разным)

| Параметр | Общий? | Разный? | Примечания |
|----------|--------|---------|------------|
| `PRINT` | ✓ | ✓ | Можно общий verbose для всех |
| `SOSCF_START_CONVERGENCE` | ✓ | ✓ | Может зависеть от wfn type |
| `DAMPING_PERCENTAGE` | | ✓ | Обычно per-wfn |

## Архитектурные решения

### Option 1: Pre-configured wavefunctions (CURRENT)

```python
# User creates each wfn with desired options
psi4.set_options({'diis_start': 0, 'mom_start': 5})
wfn1 = psi4.core.RHF(...)
wfn1.initialize()

psi4.set_options({'diis_start': 2, 'mom_start': 0})
wfn2 = psi4.core.RHF(...)
wfn2.initialize()

# Problem: wfn read options at DIFFERENT times → pollution
energies = multi_scf([wfn1, wfn2])
```

**Проблема**: Глобальное состояние. Baseline test меняет опции → последующие wfn affected.

### Option 2: Options as arguments (BETTER)

```python
# Each wfn gets explicit option dict
opts1 = {'diis_start': 0, 'mom_start': 5}
opts2 = {'diis_start': 2, 'mom_start': 0}

energies = multi_scf(
    [wfn1, wfn2],
    options=[opts1, opts2],  # Per-wfn options
    e_conv=1e-8,             # Common convergence
    d_conv=1e-6
)
```

**Плюс**: Explicit, no global state pollution
**Минус**: Need to refactor scf_iterate() to accept options dict

### Option 3: Isolated option contexts (CLEANEST)

```python
# Create isolated option contexts
with psi4.option_context({'diis_start': 0, 'mom_start': 5}):
    wfn1 = psi4.core.RHF(...)
    wfn1.initialize()

with psi4.option_context({'diis_start': 2, 'mom_start': 0}):
    wfn2 = psi4.core.RHF(...)
    wfn2.initialize()

energies = multi_scf([wfn1, wfn2])
```

**Плюс**: Clean isolation, no pollution
**Минус**: Need to implement option_context() (context manager)

### Option 4: Wavefunction-level options (BEST LONG-TERM)

```python
# Each wfn stores its own options internally
wfn1 = psi4.core.RHF(...)
wfn1.set_options({'diis_start': 0, 'mom_start': 5})
wfn1.initialize()

wfn2 = psi4.core.RHF(...)
wfn2.set_options({'diis_start': 2, 'mom_start': 0})
wfn2.initialize()

energies = multi_scf([wfn1, wfn2])
```

**Плюс**: Complete encapsulation, no global state
**Минус**: Major refactoring of HF class to store options

## Immediate Fix (Pragmatic)

**Goal**: Make multi_scf() work NOW without major refactoring

**Solution**: Detect option mismatches and WARN user, but DON'T force sync (let them differ)

```python
def multi_scf(wfn_list, ...):
    # Validate REQUIRED compatibility
    basis_set = wfn_list[0].basisset()
    for i, wfn in enumerate(wfn_list[1:], 1):
        if not wfn.basisset().equiv(basis_set):
            raise ValidationError(f"wfn[{i}] has incompatible basis set")

    # Check JK type compatibility
    jk_type = type(wfn_list[0].jk())
    for i, wfn in enumerate(wfn_list[1:], 1):
        if type(wfn.jk()) != jk_type:
            raise ValidationError(f"wfn[{i}] has incompatible JK type")

    # DETECT but DON'T fix option mismatches
    if wfn_list[0].diis_start_ != wfn_list[1].diis_start_:
        core.print_out("  ⚠ WARNING: wavefunctions have different DIIS settings\n")
        core.print_out("  This may be intentional or due to polluted global options\n")
        # Let them proceed - user may WANT different settings!

    # ... rest of multi_scf()
```

## Long-term Solution (Phase 2-3)

1. **Implement wavefunction-level options** (Option 4 above)
2. **Deprecate global options** for per-wfn settings
3. **Document multi-SCF configuration** patterns

## Example Use Cases

### Case 1: Same molecule, different spin states (SA-REKS goal)

```python
# Singlet state
wfn_singlet = psi4.core.RHF(mol, functional)
wfn_singlet.set_multiplicity(1)

# Triplet state
wfn_triplet = psi4.core.RHF(mol, functional)
wfn_triplet.set_multiplicity(3)

energies = multi_scf([wfn_singlet, wfn_triplet])
```

### Case 2: Different molecules (different sizes OK!)

```python
# H2O (10 electrons)
wfn_h2o = psi4.core.RHF(mol_h2o, functional)

# NH3 (18 electrons)
wfn_nh3 = psi4.core.RHF(mol_nh3, functional)

# Different nalpha/nbeta is OK! JK handles arbitrary C matrix sizes
energies = multi_scf([wfn_h2o, wfn_nh3])
```

### Case 3: Different convergence strategies

```python
# Aggressive DIIS (fast but risky)
wfn1 = psi4.core.RHF(...)
wfn1.set_options({'diis_start': 0, 'damping_percentage': 0})

# Conservative (slow but stable)
wfn2 = psi4.core.RHF(...)
wfn2.set_options({'diis_start': 5, 'damping_percentage': 20})

energies = multi_scf([wfn1, wfn2])
```

## References

- JK documentation: `psi4/src/psi4/libfock/jk.h` lines 90-169
- scf_iterate options: `psi4/driver/procrouting/scf_proc/scf_iterator.py` lines 314-340
- Phase 0.6 universal API: `n_states()`, `get_orbital_matrices()`, `set_jk_matrices()`
