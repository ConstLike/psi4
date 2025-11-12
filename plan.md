# SCF Refactoring Project

## Overview

**Goal:** Refactor Psi4 SCF to separate theory (Fock formation) from convergence algorithm, enabling clean implementation of advanced methods like REKS (Restricted Ensemble Kohn-Sham).

**Strategy:** **Bottom-up incremental refactoring** - from small, safe changes to complete architecture.
Each step: Claude modifies code → commits → pushes; User compiles → tests → validates.

**Key Innovation:** Multi-Fock architecture with **shared JK contraction** for efficiency.

## Workflow & Roles

**Claude (AI) - "The Hands":**
- Modifies code (C++, Python, CMake)
- Commits and pushes changes
- **DOES NOT:** compile, run tests, create test files

**User - "The Tester/Validator":**
- Pulls changes on chc2 server
- Compiles: `cd build && ninja psi4`
- Runs comprehensive test suite (covers ALL SCF features: RHF/UHF/ROHF, all DFT functionals, all SCF types, all options)
- Reports: "✅ Success" or "❌ Error: [details]"

**Key:** User's test suite is production-ready and comprehensive. Trust user's validation completely.

---

## Current Status

**Phase:** Planning & Analysis ✅

**Branch:** `claude/improve-scf-code-refactor-011CUySHMdVr9aeEhprxf6EJ`

**Completed:**
- ✅ Phase 0.0: DEBUG output added to RHF, UHF, ROHF (validated by user)
- ✅ Code analysis: current structure understood
- ✅ New refactoring plan formulated (bottom-up, multi-Fock)

**Next Action:** Begin Phase 1 - Matrix Containers

---

## Target Architecture

### Multi-Fock Design (NEW!)

**Core Idea:** SCF handles **N Fock matrices simultaneously** with **single shared JK contraction**

```
┌─────────────────────────────────────────────────┐
│         SCFDriver (convergence logic)            │
│  - DIIS, damping, convergence checks            │
│  - Iteration loop                               │
└─────────────┬───────────────────────────────────┘
              │ uses
              ↓
┌─────────────────────────────────────────────────┐
│       FockTheory (theory-specific)              │
│  n_states() → int   (RHF:1, UHF:2, REKS:N)     │
│                                                 │
│  ┌───────────────────────────────────────────┐ │
│  │ 1. build_density(C → D) [per state]      │ │
│  │ 2. compute_JK(D₁,D₂,...Dₙ → J,K) [ONCE!] │ │ ← SHARED!
│  │ 3. combine_fock(J,K → F₁,F₂,...Fₙ)       │ │ ← theory-specific
│  └───────────────────────────────────────────┘ │
└─────────────┬───────────────────────────────────┘
              │ implements
       ┌──────┴──────┬──────────┬──────────┐
       ↓             ↓          ↓          ↓
   RHFTheory    UHFTheory  REKSTheory  CustomTheory
   n=1          n=2        n=N         n=?
```

**Example: UHF with shared JK**
```cpp
// Current (inefficient): JK called twice for Da and Db
jk->C_left() = {Ca_occ};  jk->compute();  // J_a, K_a
jk->C_left() = {Cb_occ};  jk->compute();  // J_b, K_b

// New (efficient): JK called ONCE for both
jk->C_left() = {Ca_occ, Cb_occ};  jk->compute();  // J_a, K_a, J_b, K_b in ONE call
F_a = H + (J_a + J_b) - K_a
F_b = H + (J_a + J_b) - K_b
```

**Generalization:** Any number of Fock matrices!
- RHF alone: 1 Fock
- UHF alone: 2 Fock
- **RHF + ROHF simultaneously: 2 Fock** (your example!)
- REKS: N Fock (N states)

---

## Refactoring Strategy: Bottom-Up Approach

**NO HFNew!** We refactor **in-place** with **opt-in** pattern:
1. Add new code alongside old (e.g., `DensityContainer` + `Da_` both exist)
2. Add flag `use_new_path_` to switch between old/new
3. Test new path, keep old path working
4. When validated, remove old path

**Increments:** Each step takes 1-2 hours, fully tested before next step.

---

## Phases Overview

### Phase 1: Matrix Containers (Encapsulation)
**Goal:** Wrap ALL SCF matrices in multi-state containers

**1.1 Create `MultiStateMatrix` class**
- Wrapper for N matrices: `states_[i] → SharedMatrix` for i=0..N-1
- Methods: `get(i)`, `set(i, M)`, `n_states()`, `zero_all()`, `copy_all_from()`
- Location: `psi4/src/psi4/libscf_solver/multistate_matrix.h`
- Test: compile only, not used yet

**1.2 Add `MultiStateMatrix` for density in RHF (opt-in)**
- Add member: `std::shared_ptr<MultiStateMatrix> density_msm_;` (n_states=1)
- Add flag: `bool use_multistate_density_ = false;`
- In `form_D()`: if flag set, sync `density_msm_->get(0)` ↔ `Da_`
- Test: both paths work

**1.3 Add `MultiStateMatrix` for Fock in RHF (opt-in)**
- Add member: `std::shared_ptr<MultiStateMatrix> fock_msm_;` (n_states=1)
- In `form_F()`: if flag set, sync `fock_msm_->get(0)` ↔ `Fa_`
- Test: both paths work

**1.4 Add `MultiStateMatrix` for C, epsilon in RHF (opt-in)**
- Orbitals: `orbitals_msm_` → contains `Ca_`
- Energies: `epsilon_msm_` → contains `epsilon_a_`
- Test: both paths work

**1.5 Switch RHF to `MultiStateMatrix` fully**
- Remove opt-in flags, use only `*_msm_` members
- Create aliases: `SharedMatrix& Da_ = density_msm_->get(0);` for compatibility
- Test: all RHF tests pass

**1.6 Switch UHF to `MultiStateMatrix` (n_states=2)**
- `density_msm_(2)`: states[0]=Da, states[1]=Db
- `fock_msm_(2)`: states[0]=Fa, states[1]=Fb
- `orbitals_msm_(2)`: states[0]=Ca, states[1]=Cb
- Test: all UHF tests pass

**Time:** ~8-10 hours | **Risk:** LOW (pure wrapping, no logic change)

---

### Phase 2: Density Operations (Extraction)
**Goal:** Extract density formation into standalone functions

**2.1 Create `density_builder.h/cc`**
- Function: `build_density_from_orbitals(C, nocc_pi, D_out)`
- Uses DGEMM: `D = C_occ × C_occ^T`
- Location: `psi4/src/psi4/libscf_solver/density_builder.h`
- Test: unit test (standalone)

**2.2 Use `build_density_from_orbitals` in RHF**
```cpp
void RHF::form_D() {
    density_ops::build_density_from_orbitals(
        orbitals_msm_->get(0), nalphapi_, density_msm_->get(0));
}
```
- Test: all RHF tests pass

**2.3 Use `build_density_from_orbitals` in UHF**
```cpp
void UHF::form_D() {
    for (int s = 0; s < 2; ++s) {
        density_ops::build_density_from_orbitals(
            orbitals_msm_->get(s),
            s == 0 ? nalphapi_ : nbetapi_,
            density_msm_->get(s));
    }
}
```
- Test: all UHF tests pass

**Time:** ~4 hours | **Risk:** LOW (simple extraction)

---

### Phase 3: JK Operations (Multi-State Contraction)
**Goal:** Single JK call for ALL densities, store results in `JKContainer`

**3.1 Create `JKContainer` class**
- Stores J[i], K[i], wK[i] for i=0..N-1
- Methods: `get_J(i)`, `get_K(i)`, `get_wK(i)`, `n_states()`
- Location: `psi4/src/psi4/libscf_solver/jk_container.h`

**3.2 Create `compute_jk_multistate()` function**
```cpp
void compute_jk_multistate(
    std::shared_ptr<JK> jk,
    const std::vector<SharedMatrix>& C_occ_list,  // [C_a_occ, C_b_occ, ...]
    JKContainer& jk_results  // output: J[i], K[i], wK[i]
) {
    jk->C_left().clear();
    for (const auto& C : C_occ_list) {
        jk->C_left().push_back(C);
    }
    jk->compute();  // ONE CALL FOR ALL!

    for (int i = 0; i < C_occ_list.size(); ++i) {
        jk_results.set_J(jk->J()[i], i);
        jk_results.set_K(jk->K()[i], i);
        if (jk->wK().size() > 0) {
            jk_results.set_wK(jk->wK()[i], i);
        }
    }
}
```
- Test: unit test

**3.3 Use in RHF::form_G() (opt-in)**
```cpp
void RHF::form_G() {
    auto C_occ = Ca_subset("SO", "OCC");
    JKContainer jk_results(1);  // 1 state
    jk_ops::compute_jk_multistate(jk_, {C_occ}, jk_results);

    G_->zero();
    if (functional_->needs_xc()) { G_->add(Va_); }
    G_->axpy(2.0, jk_results.get_J(0));
    G_->axpy(-alpha, jk_results.get_K(0));
}
```
- Test: RHF passes

**3.4 Use in UHF::form_G() - KEY OPTIMIZATION!**
```cpp
void UHF::form_G() {
    auto Ca_occ = Ca_subset("SO", "OCC");
    auto Cb_occ = Cb_subset("SO", "OCC");

    JKContainer jk_results(2);  // 2 states
    jk_ops::compute_jk_multistate(jk_, {Ca_occ, Cb_occ}, jk_results);  // ONE CALL!

    // J_total = J_a + J_b
    auto J_total = jk_results.get_J(0)->clone();
    J_total->add(jk_results.get_J(1));

    // F_a = H + J_total - K_a
    // F_b = H + J_total - K_b
    Ga_->copy(Va_);
    Ga_->add(J_total);
    Ga_->axpy(-alpha, jk_results.get_K(0));

    Gb_->copy(Vb_);
    Gb_->add(J_total);
    Gb_->axpy(-alpha, jk_results.get_K(1));
}
```
- Test: UHF passes, **verify performance gain!**

**Time:** ~6 hours | **Risk:** MEDIUM (core performance optimization)

---

### Phase 4: Fock Operations (Assembly)
**Goal:** Standardize Fock = H + G + V_ext construction

**4.1 Create `build_fock_from_components()` function**
```cpp
void build_fock_from_components(
    const SharedMatrix& H,
    const SharedMatrix& G,
    const std::vector<SharedMatrix>& external_pots,
    SharedMatrix& F_out
) {
    F_out->copy(H);
    F_out->add(G);
    for (const auto& V : external_pots) { F_out->add(V); }
}
```
- Location: `psi4/src/psi4/libscf_solver/fock_builder.h`

**4.2 Use in RHF::form_F()**
```cpp
void RHF::form_F() {
    fock_ops::build_fock_from_components(H_, G_, external_potentials_, Fa_);
}
```

**4.3 Use in UHF::form_F()**
```cpp
void UHF::form_F() {
    fock_ops::build_fock_from_components(H_, Ga_, external_potentials_, Fa_);
    fock_ops::build_fock_from_components(H_, Gb_, external_potentials_, Fb_);
}
```
- Test: RHF, UHF pass

**Time:** ~3 hours | **Risk:** LOW

---

### Phase 5: Theory Abstraction (FockTheory Interface)
**Goal:** Encapsulate theory-specific logic (RHF/UHF/REKS) behind interface

**5.1 Create `FockTheory` abstract base class**
```cpp
class FockTheory {
public:
    virtual ~FockTheory() = default;

    virtual int n_states() const = 0;  // 1, 2, N

    virtual void build_all_densities(
        const MultiStateMatrix& orbitals,
        const std::vector<Dimension>& nocc,
        MultiStateMatrix& densities_out
    ) = 0;

    virtual void build_all_focks(
        std::shared_ptr<JK> jk,
        std::shared_ptr<VBase> potential,
        const MultiStateMatrix& densities,
        const MultiStateMatrix& orbitals,
        const SharedMatrix& H,
        const std::vector<SharedMatrix>& external_pots,
        MultiStateMatrix& focks_out
    ) = 0;
};
```

**5.2 Implement `RHFTheory : public FockTheory`**
```cpp
class RHFTheory : public FockTheory {
    int n_states() const override { return 1; }

    void build_all_densities(...) override {
        density_ops::build_density_from_orbitals(
            orbitals.get(0), nocc[0], densities_out.get(0));
    }

    void build_all_focks(...) override {
        // 1. Compute JK (1 state)
        auto C_occ = extract_occupied(orbitals.get(0), nocc[0]);
        JKContainer jk_results(1);
        jk_ops::compute_jk_multistate(jk, {C_occ}, jk_results);

        // 2. Build G = V_xc + 2J - αK
        auto G = focks_out.get(0);
        G->zero();
        if (potential) { /* add V_xc */ }
        G->axpy(2.0, jk_results.get_J(0));
        G->axpy(-alpha, jk_results.get_K(0));

        // 3. F = H + G + V_ext
        fock_ops::build_fock_from_components(H, G, external_pots, focks_out.get(0));
    }
};
```

**5.3 Implement `UHFTheory : public FockTheory`**
```cpp
class UHFTheory : public FockTheory {
    int n_states() const override { return 2; }

    void build_all_densities(...) override {
        for (int s = 0; s < 2; ++s) {
            density_ops::build_density_from_orbitals(
                orbitals.get(s), nocc[s], densities_out.get(s));
        }
    }

    void build_all_focks(...) override {
        // 1. Compute JK (2 states) - ONE CALL!
        auto Ca_occ = extract_occupied(orbitals.get(0), nocc[0]);
        auto Cb_occ = extract_occupied(orbitals.get(1), nocc[1]);
        JKContainer jk_results(2);
        jk_ops::compute_jk_multistate(jk, {Ca_occ, Cb_occ}, jk_results);

        // 2. J_total = J_a + J_b
        auto J_total = jk_results.get_J(0)->clone();
        J_total->add(jk_results.get_J(1));

        // 3. Build Ga, Gb
        // Ga = V_xc_a + J_total - αK_a
        // Gb = V_xc_b + J_total - αK_b
        // ...

        // 4. Fa = H + Ga, Fb = H + Gb
        fock_ops::build_fock_from_components(H, Ga, external_pots, focks_out.get(0));
        fock_ops::build_fock_from_components(H, Gb, external_pots, focks_out.get(1));
    }
};
```

**5.4 Integrate `FockTheory` into RHF (opt-in)**
```cpp
// In rhf.h
protected:
    std::shared_ptr<FockTheory> fock_theory_;
    bool use_fock_theory_ = false;

// In rhf.cc constructor
fock_theory_ = std::make_shared<RHFTheory>();

// In form_D(), form_F()
if (use_fock_theory_) {
    fock_theory_->build_all_densities(...);
    fock_theory_->build_all_focks(...);
} else {
    // old path
}
```
- Test: both paths work

**5.5 Switch RHF to `FockTheory` fully**
- Remove old `form_D()`, `form_G()`, `form_F()` code
- Keep only `fock_theory_->build_all_*()` calls
- Test: all RHF tests pass

**5.6 Switch UHF to `FockTheory` fully**
- Test: all UHF tests pass

**Time:** ~12 hours | **Risk:** MEDIUM-HIGH (major refactor)

---

### Phase 6: SCF Driver (Algorithm Separation)
**Goal:** Extract convergence logic from HF into standalone driver

**6.1 Create `SCFDriver` class**
```cpp
class SCFDriver {
    std::shared_ptr<FockTheory> theory_;
    std::shared_ptr<DIISManager> diis_;
    double E_convergence_;
    double D_convergence_;

public:
    double iterate_to_convergence() {
        for (int iter = 0; iter < max_iter_; ++iter) {
            // 1. Build densities
            theory_->build_all_densities(...);

            // 2. Build Focks
            theory_->build_all_focks(...);

            // 3. DIIS extrapolation
            if (diis_) { diis_->extrapolate(...); }

            // 4. Diagonalize
            for (int s = 0; s < theory_->n_states(); ++s) {
                diagonalize_fock(focks_.get(s), orbitals_.get(s), epsilon_.get(s));
            }

            // 5. Compute energy
            double E = theory_->compute_energy(...);

            // 6. Check convergence
            if (converged(E, D)) break;
        }
    }
};
```

**6.2 Integrate into RHF (opt-in)**
**6.3 Switch RHF to SCFDriver fully**
**6.4 Switch UHF to SCFDriver fully**

**Time:** ~10 hours | **Risk:** HIGH (complete separation)

---

## Implementation Roadmap

| Phase | Focus | Steps | Time | Risk | Status |
|-------|-------|-------|------|------|--------|
| 0 | Analysis | - | 2h | - | ✅ DONE |
| 1 | Matrix Containers | 6 | 10h | LOW | 📍 NEXT |
| 2 | Density Ops | 3 | 4h | LOW | Pending |
| 3 | JK Ops (Multi-State!) | 4 | 6h | MED | Pending |
| 4 | Fock Ops | 3 | 3h | LOW | Pending |
| 5 | Theory Abstraction | 6 | 12h | MED-HIGH | Pending |
| 6 | SCF Driver | 4 | 10h | HIGH | Pending |
| **Total** | **RHF+UHF refactored** | **26** | **~45h** | - | 0% |

---

## Success Criteria (User Validates)

### Phase 1
- [ ] `MultiStateMatrix` compiles
- [ ] RHF with `*_msm_` passes all tests
- [ ] UHF with `*_msm_` (n=2) passes all tests

### Phase 2
- [ ] Density formation extracted to standalone function
- [ ] RHF, UHF use `density_ops::build_density_from_orbitals()`
- [ ] All tests pass

### Phase 3 ⭐ **KEY PHASE**
- [ ] `JKContainer` works
- [ ] UHF uses single JK call for both α and β
- [ ] **Performance:** UHF JK time reduced by ~40-50%
- [ ] All tests pass

### Phase 4
- [ ] Fock assembly standardized
- [ ] All tests pass

### Phase 5
- [ ] `FockTheory` interface complete
- [ ] `RHFTheory`, `UHFTheory` implemented
- [ ] RHF, UHF use `FockTheory`
- [ ] All tests pass, exact energy match (< 1e-10)

### Phase 6
- [ ] `SCFDriver` functional
- [ ] RHF, UHF use `SCFDriver`
- [ ] All tests pass
- [ ] **Architecture ready for REKS!**

---

## Next Steps

**Immediate:**
1. Review this plan with User
2. Adjust if needed
3. Begin Phase 1.1: Create `MultiStateMatrix`

**Questions for User:**
1. Agree with bottom-up approach (no HFNew)?
2. Agree with multi-Fock architecture?
3. Phase 3 (multi-state JK) is critical - concerns?
4. Any missing aspects?

---

## Technical Notes

### Why Multi-State JK is Efficient

**Current UHF:**
```cpp
// Two separate JK calls
jk->C_left() = {Ca_occ};  jk->compute();  // ~T seconds
jk->C_left() = {Cb_occ};  jk->compute();  // ~T seconds
// Total: 2T
```

**New UHF:**
```cpp
// Single JK call
jk->C_left() = {Ca_occ, Cb_occ};  jk->compute();  // ~T seconds (not 2T!)
// Total: T  (2x faster!)
```

**Why?** JK internally:
- Builds Schwarz screening once
- Loops over shell quartets once
- Contracts with multiple densities in inner loop (vectorized)

**Real-world check (current UHF):**
```cpp
// In uhf.cc line 184-200
std::vector<SharedMatrix>& C = jk_->C_left();
C.clear();
C.push_back(Ca_subset("SO", "OCC"));
C.push_back(Cb_subset("SO", "OCC"));  // Already supports multi-state!
jk_->compute();  // ONE call
```
✅ **Good news:** Current UHF already does this! But we'll systematize it for REKS.

---

## Risk Mitigation

**If blocked >2 hours on a step:**
1. Revert last commit
2. Break step into smaller pieces
3. Ask User for clarification
4. Document blocker in plan.md

**If tests fail:**
1. Check compilation warnings
2. User runs single test with verbose output
3. Compare intermediate values (D, F matrices) via DEBUG output
4. Revert if can't fix in 30 minutes

**If performance regresses >20%:**
1. User profiles with perf/gprof
2. Check for unnecessary copies (use const references)
3. Verify BLAS usage (DGEMM, DAXPY)
4. May be acceptable if architecture is cleaner (document trade-off)

---

## Testing Strategy

**User's Comprehensive Test Suite (Already Exists):**
- RHF, UHF, ROHF
- All DFT functionals: HF, PBE, B3LYP, ωB97X-V, etc.
- All SCF types: PK, DF, MEM_DF, DISK_DF, OUT_OF_CORE, CD
- All screening: CSAM, DENSITY, NONE
- All guess: CORE, SAD, AUTO, MODHUCKEL, GWH, SADNO, SAP, SAPGAU
- Special: incremental Fock, DIIS, damping, MOM, fractional occupation

**Per-Phase Thresholds:**
- Phase 1: Compilation success, no crashes
- Phase 2: Compilation + all tests pass
- Phase 3: Energy deviation < 1e-10 (exact match)
- Phase 4: Energy deviation < 1e-10
- Phase 5: Energy deviation < 1e-10 + exact iteration match
- Phase 6: Full equivalence + performance acceptable

---

## Future: REKS Implementation

After Phase 6, the architecture will support REKS naturally:

```cpp
class REKSTheory : public FockTheory {
    int n_states() const override { return n_ensemble_; }  // 2-6 typically

    void build_all_focks(...) override {
        // 1. Build sub-densities D_I, D_J from ensemble
        // 2. Compute JK for all sub-densities (ONE call)
        JKContainer jk_results(n_states_);
        jk_ops::compute_jk_multistate(jk, {D_I, D_J, ...}, jk_results);

        // 3. Build state-averaged Fock or multiple Fock matrices
        // F = Σ w_K * (H + J_total - K_K)
        // ...
    }
};
```

**Estimated time:** ~15-20 hours after Phase 6 complete.

---

## Appendix: Current Code Structure

### RHF key methods (rhf.cc)
- `form_D()` line 352: D = C × C^T via DGEMM ✅ already optimal
- `form_G()` line 187: JK contraction + DFT
- `form_F()` line 277: F = H + G + V_ext
- `form_C()` line 317: Diagonalize F
- `compute_E()` line 398: Energy calculation

### UHF key methods (uhf.cc)
- `form_D()` line 427: Da, Db separately
- `form_G()` line 184: **Already multi-state JK!** (α and β in one call)
- `form_F()` line 156: Fa, Fb separately

### Key observation:
✅ UHF already uses multi-state JK contraction efficiently
❌ But architecture is not generalized for N states (REKS)
→ Our refactoring systematizes this pattern
