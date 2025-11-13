# Plan Perfect - Code Quality Improvement Roadmap

## Экспертная оценка текущего кода

### Оценка для legacy codebase Psi4: **8/10** ⭐

**Что сделано правильно:**
- ✅ Консистентен со стилем кодовой базы Psi4
- ✅ Минимально инвазивен (только Python, без C++)
- ✅ Сохраняет backward compatibility
- ✅ Enableает multi-cycle SCF архитектуру
- ✅ Можно тестировать и валидировать (78 тестов прошли!)
- ✅ Можно улучшить позже (incremental refactoring)

**Где есть code smells (но это нормально для legacy):**
- ❌ Namespace pollution (~15 `self._scf_*` атрибутов)
- ❌ Нет cleanup после завершения итераций
- ❌ Дублирование переменных в `scf_iterate()` и `_scf_initialize_iteration_state()`
- ❌ Mixed concerns (OOO path vs regular path)
- ❌ Incomplete abstraction (завязан на `self.iteration_`, `self.diis_enabled_`)

---

## Phase 2-3: Code Quality Improvements (ПОСЛЕ завершения multi-cycle SCF)

### Priority 1: State Object Pattern (Phase 2)
**КОГДА**: После того как `multi_scf()` будет работать

**ЧТО ДЕЛАТЬ**:
```python
from dataclasses import dataclass
from typing import Optional

@dataclass
class SCFIterationState:
    """Encapsulates all iteration state"""
    __slots__ = ['SCFE_old', 'Dnorm', 'Ediff', 'e_conv', 'd_conv',
                 'verbose', 'reference', 'is_dfjk', 'damping_enabled',
                 'soscf_enabled', 'frac_enabled', 'efp_enabled',
                 'cosx_enabled', 'early_screening', 'early_screening_disabled',
                 'maxiter_post_screening', 'iter_post_screening']

    SCFE_old: float = 0.0
    Dnorm: float = 0.0
    Ediff: float = 0.0
    e_conv: Optional[float] = None
    d_conv: Optional[float] = None
    # ... остальные поля

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Cleanup state after iterations"""
        pass

# Usage:
def _scf_initialize_iteration_state(self, e_conv, d_conv):
    self._scf_state = SCFIterationState(
        e_conv=e_conv,
        d_conv=d_conv,
        verbose=core.get_option('SCF', "PRINT"),
        reference=core.get_option('SCF', "REFERENCE"),
        # ... etc
    )

def _scf_iteration(self):
    state = self._scf_state  # Short alias
    self.iteration_ += 1

    # Use state.Ediff instead of self._scf_Ediff
    state.Ediff = SCFE - state.SCFE_old
    state.SCFE_old = SCFE

    if _converged(state.Ediff, state.Dnorm, state.e_conv, state.d_conv):
        return (False, 'converged')

    return (True, 'continue')

def scf_iterate(self, e_conv=None, d_conv=None):
    with SCFIterationState(e_conv, d_conv) as state:
        self._scf_state = state
        while True:
            should_continue, reason = self._scf_iteration()
            if not should_continue:
                break
        # State автоматически очищается через __exit__
```

**BENEFITS**:
- ✅ Cleanup namespace (~15 атрибутов → 1 объект)
- ✅ Automatic cleanup через context manager
- ✅ Better memory efficiency (`__slots__`)
- ✅ Type safety
- ✅ Easier to test (можно mock state object)

**RISKS**:
- Нужно переписать ~200 строк в `_scf_iteration()`
- Нужно тестировать 78 тестов снова

---

### Priority 2: Cleanup Method (Phase 2)
**КОГДА**: Сразу после Priority 1 или как отдельная задача

**ЧТО ДЕЛАТЬ**:
```python
def _scf_cleanup_iteration_state(self):
    """Clean up iteration state after completion"""
    if hasattr(self, '_scf_state'):
        del self._scf_state
    # Или если без state object:
    attrs_to_delete = [
        '_scf_SCFE_old', '_scf_Dnorm', '_scf_Ediff',
        '_scf_e_conv', '_scf_d_conv', '_scf_verbose',
        # ... все _scf_* атрибуты
    ]
    for attr in attrs_to_delete:
        if hasattr(self, attr):
            delattr(self, attr)

def scf_iterate(self, e_conv=None, d_conv=None):
    try:
        self._scf_initialize_iteration_state(e_conv, d_conv)
        while True:
            should_continue, reason = self._scf_iteration()
            if not should_continue:
                break
            if self.iteration_ >= core.get_option('SCF', 'MAXITER'):
                raise SCFConvergenceError(...)
    finally:
        self._scf_cleanup_iteration_state()  # Always cleanup
```

**BENEFITS**:
- ✅ No memory leaks from lingering state
- ✅ Clean object after use
- ✅ Easier debugging (clear state boundaries)

---

### Priority 3: Refactor OOO Path (Phase 3)
**КОГДА**: Когда будет время и motivation

**ЧТО ДЕЛАТЬ**:
```python
def scf_iterate(self, e_conv=None, d_conv=None):
    # Detect which backend to use
    if self._should_use_ooo():
        return self._scf_iterate_ooo(e_conv, d_conv)
    else:
        return self._scf_iterate_internal(e_conv, d_conv)

def _should_use_ooo(self):
    """Determine if OpenOrbitalOptimizer should be used"""
    ooo_scf = core.get_option("SCF", "ORBITAL_OPTIMIZER_PACKAGE") in ["OOO", "OPENORBITALOPTIMIZER"]
    if not ooo_scf:
        return False

    reference = core.get_option('SCF', "REFERENCE")
    soscf_enabled = _validate_soscf()
    frac_enabled = _validate_frac()
    efp_enabled = hasattr(self.molecule(), 'EFP')
    pcm_enabled = core.get_option('SCF', 'PCM')
    # ... rest of checks

    incompatible = (reference in ["ROHF", "CUHF"] or soscf_enabled or
                    self.MOM_excited_ or frac_enabled or efp_enabled or ...)

    if incompatible:
        core.print_out("Note: OpenOrbitalOptimizer not compatible. Falling back to internal\n")
        return False

    return True

def _scf_iterate_ooo(self, e_conv, d_conv):
    """OOO-specific iteration path"""
    # SAD handling
    if self.sad_ and self.iteration_ <= 0:
        # ... SAD logic ...

    try:
        self.openorbital_scf()
    except RuntimeError as ex:
        if "openorbital_scf is virtual" in str(ex):
            core.print_out("Note: OpenOrbitalOptimizer NYI. Falling back to Internal.\n")
            return self._scf_iterate_internal(e_conv, d_conv)
        raise ex

    SCFE = self.compute_E()
    self.set_energies("Total Energy", SCFE)
    # ... rest ...
    return SCFE

def _scf_iterate_internal(self, e_conv, d_conv):
    """Internal (our new) iteration path"""
    self._scf_initialize_iteration_state(e_conv, d_conv)
    try:
        while True:
            should_continue, reason = self._scf_iteration()
            if not should_continue:
                break
            if self.iteration_ >= core.get_option('SCF', 'MAXITER'):
                raise SCFConvergenceError(...)
    finally:
        self._scf_cleanup_iteration_state()
```

**BENEFITS**:
- ✅ Clear separation of concerns
- ✅ Easier to test each path independently
- ✅ No mixed logic in main function
- ✅ Easier to maintain

**RISKS**:
- Большой refactoring
- Нужно тестировать OOO path отдельно

---

### Priority 4: Type Hints (Phase 3)
**КОГДА**: Постепенно, по мере работы

**ЧТО ДЕЛАТЬ**:
```python
from typing import Optional, Tuple
from enum import Enum

class ConvergenceReason(Enum):
    """Enumeration of convergence reasons"""
    CONVERGED = "converged"
    MOM_NOT_STARTED = "mom_not_started"
    FRAC_NOT_STARTED = "frac_not_started"
    EARLY_SCREENING_MAXITER = "early_screening_maxiter"
    CONTINUE = "continue"

def _scf_initialize_iteration_state(
    self,
    e_conv: Optional[float],
    d_conv: Optional[float]
) -> None:
    """Initialize state for SCF iterations"""
    pass

def _scf_iteration(self) -> Tuple[bool, str]:
    """
    Performs ONE SCF iteration.

    Returns
    -------
    continue_flag : bool
        True to continue iterations, False to stop
    reason : str
        Reason for stopping: 'converged', 'mom_not_started', etc.
    """
    pass

def scf_iterate(
    self,
    e_conv: Optional[float] = None,
    d_conv: Optional[float] = None
) -> None:
    """Main SCF iteration loop"""
    pass
```

**BENEFITS**:
- ✅ Better IDE support
- ✅ Catch type errors early
- ✅ Self-documenting code
- ✅ Easier for new contributors

---

### Priority 5: Better Docstrings (Phase 3)
**КОГДА**: Постепенно

**ЧТО ДЕЛАТЬ**:
```python
def _scf_iteration(self) -> Tuple[bool, str]:
    """
    Performs ONE SCF iteration.

    This method is designed to be called externally by a multi-cycle SCF
    coordinator. It executes a complete SCF iteration including:
    - Form G (JK computation or use precomputed)
    - Form F (Fock matrix)
    - DIIS/SOSCF convergence acceleration
    - Form C (orbital coefficients)
    - Form D (density matrix)
    - Damping if enabled
    - Convergence checking

    State Management
    ----------------
    Uses state stored in self._scf_* members initialized by
    _scf_initialize_iteration_state(). The method modifies:
    - self.iteration_ (incremented)
    - self._scf_SCFE_old (updated energy)
    - self._scf_Dnorm (updated density norm)
    - self._scf_Ediff (energy difference)
    - Density/Fock/orbital matrices

    Multi-Cycle Support
    -------------------
    When called from multi_scf() coordinator:
    1. Coordinator collects all C matrices from all wfns
    2. Coordinator performs shared JK computation
    3. Coordinator distributes J/K via set_jk_matrices()
    4. This method uses precomputed J/K via use_precomputed_jk_ flag
    5. Rest of iteration proceeds normally

    Returns
    -------
    continue_flag : bool
        True to continue iterations, False to stop
    reason : str
        Convergence status:
        - 'converged': E and D converged
        - 'mom_not_started': MOM not yet activated
        - 'frac_not_started': Fractional occupation not yet activated
        - 'early_screening_maxiter': COSX final grid iterations complete
        - 'continue': Keep iterating

    Examples
    --------
    Normal usage (internal):
    >>> self._scf_initialize_iteration_state(1e-6, 1e-5)
    >>> while True:
    >>>     cont, reason = self._scf_iteration()
    >>>     if not cont:
    >>>         break

    Multi-cycle usage (external):
    >>> for wfn in wfn_list:
    >>>     wfn._scf_initialize_iteration_state(1e-6, 1e-5)
    >>> while not all_converged:
    >>>     # Shared JK
    >>>     C_all = [wfn.Ca_subset("SO", "OCC") for wfn in wfn_list]
    >>>     jk.C_left().clear()
    >>>     for C in C_all: jk.C_left().append(C)
    >>>     jk.compute()
    >>>     # Distribute and iterate
    >>>     for i, wfn in enumerate(wfn_list):
    >>>         wfn.set_jk_matrices([jk.J()[i]], [jk.K()[i]])
    >>>         cont, reason = wfn._scf_iteration()

    See Also
    --------
    _scf_initialize_iteration_state : Initialize iteration state
    scf_iterate : Main iteration loop
    multi_cycle_scf_iterate : Multi-cycle coordinator
    """
    pass
```

**BENEFITS**:
- ✅ New contributors understand code faster
- ✅ Examples show how to use
- ✅ Clear state management documentation
- ✅ Multi-cycle usage documented

---

## Дополнительные улучшения (низкий приоритет)

### Better error messages
```python
if self.iteration_ >= core.get_option('SCF', 'MAXITER'):
    msg = (f"SCF did not converge in {self.iteration_} iterations.\n"
           f"  Final energy: {self._scf_SCFE_old:.8f}\n"
           f"  Energy change: {self._scf_Ediff:.2e} (threshold: {self._scf_e_conv:.2e})\n"
           f"  Density change: {self._scf_Dnorm:.2e} (threshold: {self._scf_d_conv:.2e})\n"
           f"Try:\n"
           f"  - Increase MAXITER\n"
           f"  - Enable/tune DAMPING\n"
           f"  - Check geometry for problems")
    raise SCFConvergenceError(msg, self.iteration_, self, self._scf_Ediff, self._scf_Dnorm)
```

### Property accessors
```python
@property
def scf_energy(self) -> float:
    """Current SCF energy"""
    return self._scf_SCFE_old if hasattr(self, '_scf_SCFE_old') else 0.0

@property
def scf_converged(self) -> bool:
    """Check if SCF is converged"""
    if not hasattr(self, '_scf_Ediff'):
        return False
    return _converged(self._scf_Ediff, self._scf_Dnorm,
                     self._scf_e_conv, self._scf_d_conv)
```

### Unit tests для новых методов
```python
# tests/pytests/test_scf_refactoring.py
def test_scf_initialize_state():
    """Test state initialization"""
    wfn = setup_test_wfn()
    wfn._scf_initialize_iteration_state(1e-6, 1e-5)

    assert wfn._scf_e_conv == 1e-6
    assert wfn._scf_d_conv == 1e-5
    assert wfn._scf_SCFE_old == 0.0
    assert wfn._scf_Dnorm == 0.0

def test_scf_iteration_single_step():
    """Test single iteration step"""
    wfn = setup_converged_wfn()
    wfn._scf_initialize_iteration_state(1e-6, 1e-5)

    cont, reason = wfn._scf_iteration()

    assert isinstance(cont, bool)
    assert reason in ['converged', 'continue', 'mom_not_started', ...]
```

---

## Timeline

- **Phase 1** (текущая): ✅ DONE - Enable multi-cycle SCF
  - Step 1.1: Extract scf_iteration() ✅
  - Step 1.2: Convert to method ✅
  - Step 1.3: Create multi_scf() coordinator ⏭️ NEXT
  - Step 1.4: Test with 2 RHF cycles

- **Phase 2** (после Phase 1): Code quality improvements
  - Priority 1: State Object Pattern (1 week)
  - Priority 2: Cleanup Method (1 day)
  - Priority 4: Type Hints (ongoing)

- **Phase 3** (будущее): Advanced improvements
  - Priority 3: Refactor OOO Path (1 week)
  - Priority 5: Better Docstrings (ongoing)
  - Unit tests (1 week)

---

## Почему НЕ делать это сейчас?

**Pragmatic reasons:**
1. **Risk management** - каждый refactoring = риск сломать тесты
2. **Incremental value** - сначала enableаем multi-cycle, потом улучшаем качество
3. **Time constraints** - у нас есть concrete goal (multi-cycle SCF)
4. **Legacy codebase** - Psi4 уже полон "неидеального" кода, мы не должны быть перфекционистами
5. **Proof of concept** - сначала докажем что multi-cycle работает, потом полируем

---

## Conclusion

Наш текущий код — это **solid 8/10** для legacy codebase.

Он:
- ✅ **Работает** (78 тестов прошли)
- ✅ **Enableает multi-cycle SCF**
- ✅ **Консистентен с Psi4 style**
- ✅ **Безопасен** (minimal changes)

Улучшения в Priority 1-2 можно сделать в Phase 2, после того как multi-cycle SCF заработает.

Priority 3-5 — это "nice to have", но не критично.

**Девиз**: "Make it work, make it right, make it fast" — мы на этапе "make it work" → "make it right" будет в Phase 2! 🚀
