# SCF Refactoring: Implementation Roadmap

## Executive Summary

**Goal:** Refactor Psi4 SCF to separate theory (Fock formation) from algorithm (convergence), enabling clean REKS implementation.

**Strategy:** Incremental refactoring with continuous testing. Each micro-step: I modify code → commit → push; User pulls → compiles → tests → reports results.

**Roles:**
- **Claude (AI)**: Modifies code, commits changes (the "hands")
- **User**: Compiles on chc2, runs comprehensive test suite, analyzes results, reports errors (the "tester")

**Timeline:** ~20-25 micro-steps over phases 0-3.

---

## Workflow (CRITICAL - READ FIRST!)

### Division of Labor

**Claude (AI) - "The Hands":**
- Reads existing code
- Modifies code directly (C++, Python, CMake)
- Commits changes with descriptive messages
- Pushes to branch `claude/improve-scf-code-refactor-011CUySHMdVr9aeEhprxf6EJ`
- **DOES NOT**: compile, run tests, create test files, validate results

**User - "The Tester/Validator":**
- Pulls latest changes from Claude's branch
- Compiles on chc2 server: `cd build && ninja psi4`
- Runs comprehensive test suite (already exists!)
  - Covers: RHF, UHF, ROHF
  - All DFT functionals (HF, PBE, B3LYP, etc.)
  - All SCF types (PK, DF, MEM_DF, DISK_DF, OUT_OF_CORE, CD)
  - All screening methods, guess methods
  - Special features (incremental Fock, DIIS, damping, MOM, fractional occupation)
- Analyzes results
- Reports back: "✅ Compilation success, all tests pass" OR "❌ Error in file X, line Y"

### Per-Step Workflow

For EACH micro-step in this plan:

1. **Claude**: Modifies code according to step description
2. **Claude**: Commits with message from plan
3. **Claude**: Pushes to remote
4. **Claude**: Informs user "Ready for testing"
5. **User**: Pulls, compiles, tests
6. **User**: Reports results
7. **If errors**: Claude fixes → repeat from step 1
8. **If success**: Move to next step

### Key Principle

User's test suite is **comprehensive and production-ready**. No need to create new tests. Trust user's validation completely.

---

## Architecture Overview (Target)

```
┌─────────────────────────────────────────────────────┐
│                   SCFDriver                         │
│  (convergence loop, algorithm-agnostic)             │
│  - iteration loop                                   │
│  - convergence check                                │
└──────────────┬──────────────────────────────────────┘
               │ uses
               ↓
┌─────────────────────────────────────────────────────┐
│                FockBuilder (abstract)                │
│  Interface:                                         │
│  - n_states() → int                                 │
│  - build_density()                                  │
│  - apply_kernel()  (D → F via ERI/DFT)             │
│  - diagonalize_fock()                               │
│  - compute_energy()                                 │
│  - get_fock(state), get_density(state), ...        │
└──────────────┬──────────────────────────────────────┘
               │ implements
        ┌──────┴──────┬─────────┬──────────┐
        ↓             ↓         ↓          ↓
   RHFFockBuilder  UHF     REKS       RKS
```

**Key principles:**
- One Fock per "state" (RHF: 1, UHF: 2, REKS: N or 1 SA)
- `n_states()` determines API multiplicity
- Theory encapsulated in FockBuilder
- Algorithm (DIIS, damping, MOM) in separate accelerators/modifiers
- Backward compatible via HF adapter

---

## Phase 0: Setup & Baseline (3 micro-steps)

**Goal:** Establish baseline and create sandbox (HFNew) for refactoring.

**Note:** User already has comprehensive test suite covering all SCF features (RHF, UHF, ROHF, all DFT functionals, all SCF types, all screening methods, all guess methods, special features). No need to create tests.

### Step 0.0: Add DEBUG output to all SCF variants ✅ COMPLETED
**What:** Added comprehensive DEBUG output to RHF, UHF, ROHF
**Status:** ✅ Completed and validated by user
**Result:** All SCF code working correctly, DEBUG output functional

---

### Step 0.1: Copy HF class to HFNew (exact copy)
**What:** Duplicate `hf.h/cc` → `hf_new.h/cc` with zero changes
**Why:** Create sandbox for refactoring without touching original
**Files:**
- Create `psi4/src/psi4/libscf_solver/hf_new.h`
- Create `psi4/src/psi4/libscf_solver/hf_new.cc`
- Update CMakeLists.txt
**Changes:**
```cpp
// hf_new.h
class HFNew : public Wavefunction {  // Exact copy of HF
    // ... all methods identical
};
```
**Workflow:**
1. I create HFNew files → commit → push
2. User pulls → compiles → reports compilation result

**Acceptance (User confirms):**
- [ ] Code compiles successfully
- [ ] HFNew class exists but not used anywhere
- [ ] Original HF unchanged
**Time:** ~30 minutes
**Commit:** "Create HFNew class as exact copy of HF for refactoring"

---

### Step 0.2: Wire HFNew to Python (parallel to HF)
**What:** Make HFNew callable from Python without breaking HF
**Why:** Enable testing of new implementation
**Files:**
- `psi4/src/psi4/libscf_solver/export_scf.cc` (add HFNew binding)
- `psi4/driver/procrouting/scf_proc.py` (add 'scf_new' method)
**Changes:**
```python
# scf_proc.py
def run_scf_new(name, **kwargs):
    # Call HFNew instead of HF
    ref_wfn = core.HFNew(ref_wfn, functional)
    return ref_wfn.compute_energy()
```
**Workflow:**
1. I wire HFNew to Python → commit → push
2. User pulls → compiles → runs test suite comparing 'scf' vs 'scf_new' → reports results

**Acceptance (User confirms):**
- [ ] Code compiles successfully
- [ ] `psi4.energy('scf_new')` works
- [ ] Gives identical results to 'scf' (< 1e-10 Hartree)
- [ ] All tests pass for both 'scf' and 'scf_new'
**Time:** ~2 hours
**Commit:** "Add scf_new method using HFNew (identical to scf)"

---

**Phase 0 Complete Checklist:**
- [x] DEBUG output added to RHF, UHF, ROHF (Step 0.0)
- [ ] HFNew class created (exact copy of HF) (Step 0.1)
- [ ] scf_new callable from Python (Step 0.2)
- [ ] User validates: scf_new ≡ scf (< 1e-10)

**Time:** ~3 hours remaining (Step 0.0 done)
**Risk:** LOW - no actual refactoring yet, just setup

---

## Phase 1: Extract Core Interfaces (8 micro-steps)

**Goal:** Create FockBuilder interface and RHFFockBuilder without breaking anything.

### Step 1.1: Create FockBuilder abstract interface (empty)
**What:** Define FockBuilder.h with pure virtual methods
**Why:** Establish API contract before implementation
**Files:** Create `psi4/src/psi4/libscf_solver/fock_builder.h`
**Code:**
```cpp
class FockBuilder {
public:
    virtual ~FockBuilder() = default;
    
    virtual void initialize(...) = 0;
    virtual int n_states() const = 0;
    virtual void build_density() = 0;
    virtual void apply_kernel() = 0;
    virtual void diagonalize_fock() = 0;
    virtual double compute_energy() const = 0;
    
    virtual SharedMatrix get_fock(int state = 0) = 0;
    virtual SharedMatrix get_density(int state = 0) = 0;
    virtual SharedMatrix get_orbitals(int state = 0) = 0;
    virtual SharedMatrix get_overlap() const = 0;
    
    virtual void set_fock(SharedMatrix F, int state = 0) = 0;
};
```
**Test:**
```bash
ninja psi4  # Should compile
# No functional tests - interface only
```
**Acceptance:**
- [ ] FockBuilder.h compiles
- [ ] No implementations yet
- [ ] No breakage of existing code
**Time:** ~1 hour
**Commit:** "Add FockBuilder abstract interface"

---

### Step 1.2: Identify HFNew methods mapping to FockBuilder
**What:** Document which HFNew methods correspond to which FockBuilder methods
**Why:** Plan the extraction
**Artifact:** Create doc comment in HFNew:
```cpp
// HFNew method mapping to FockBuilder:
// form_D()      → build_density()
// form_G()      → apply_kernel() [partial: ERI contraction]
// form_F()      → apply_kernel() [partial: F = H + G]
// form_C()      → diagonalize_fock()
// compute_E()   → compute_energy()
```
**Test:** Documentation only, no code changes
**Acceptance:**
- [ ] All major HFNew methods mapped
- [ ] Documented in code comments
**Time:** ~30 minutes
**Commit:** "Document HFNew to FockBuilder method mapping"

---

### Step 1.3: Extract build_density from HFNew
**What:** Create RHFFockBuilder::build_density() by copying from HFNew
**Why:** First concrete method extraction
**Files:** Create `psi4/src/psi4/libscf_solver/rhf_fock_builder.h/cc`
**Code:**
```cpp
class RHFFockBuilder : public FockBuilder {
public:
    void build_density() override {
        // Copy implementation from HFNew::form_D()
        // Use DGEMM as planned
        D_prev_->copy(D_);
        
        auto Cp = C_->pointer();
        auto Dp = D_->pointer();
        C_DGEMM('N', 'T', nbf_, nbf_, nocc_,
                1.0, Cp[0], nbf_,
                     Cp[0], nbf_,
                0.0, Dp[0], nbf_);
    }
    
    int n_states() const override { return 1; }
    
    // Other methods: throw std::runtime_error("Not implemented yet")
};
```
**Test:**
```cpp
// Unit test
RHFFockBuilder builder;
builder.initialize(...);
builder.build_density();
// Check D = C × C^T
```
**Acceptance:**
- [ ] Compiles
- [ ] Unit test for build_density passes
- [ ] HFNew still works (not using builder yet)
**Time:** ~2 hours
**Commit:** "Implement RHFFockBuilder::build_density with DGEMM"

---

### Step 1.4: Extract apply_kernel from HFNew (RHF)
**What:** Implement RHFFockBuilder::apply_kernel()
**Why:** Second major method
**Code:**
```cpp
void RHFFockBuilder::apply_kernel() override {
    // Copy from HFNew::form_G() + form_F()
    jk_->C_left().clear();
    jk_->C_left().push_back(C_);
    jk_->compute();
    
    F_->copy(H_);
    F_->axpy(2.0, jk_->J()[0]);
    F_->axpy(-1.0, jk_->K()[0]);
}
```
**Test:**
```cpp
// Unit test
builder.build_density();
builder.apply_kernel();
// Check F = H + 2J - K
```
**Acceptance:**
- [ ] apply_kernel unit test passes
- [ ] Fock matrix correct
**Time:** ~2 hours
**Commit:** "Implement RHFFockBuilder::apply_kernel"

---

### Step 1.5: Extract diagonalize_fock from HFNew
**What:** Implement RHFFockBuilder::diagonalize_fock()
**Code:**
```cpp
void RHFFockBuilder::diagonalize_fock() override {
    // Copy from HFNew::form_C()
    F_->diagonalize(C_, epsilon_, S_);
}
```
**Test:**
```cpp
builder.diagonalize_fock();
// Check FC = SCε
```
**Acceptance:**
- [ ] Diagonalization test passes
**Time:** ~1 hour
**Commit:** "Implement RHFFockBuilder::diagonalize_fock"

---

### Step 1.6: Extract compute_energy from HFNew
**What:** Implement RHFFockBuilder::compute_energy()
**Code:**
```cpp
double RHFFockBuilder::compute_energy() const override {
    // Copy from HFNew::compute_E()
    double E_elec = D_->vector_dot(H_);
    E_elec += D_->vector_dot(F_);
    return E_elec;  // Nuclear added by driver
}
```
**Test:**
```cpp
double E = builder.compute_energy();
// Check against expected
```
**Acceptance:**
- [ ] Energy calculation correct
**Time:** ~1 hour
**Commit:** "Implement RHFFockBuilder::compute_energy"

---

### Step 1.7: Create SCFDriverNew using RHFFockBuilder
**What:** Create minimal SCF loop using FockBuilder
**Files:** Create `psi4/src/psi4/libscf_solver/scf_driver.h/cc`
**Code:**
```cpp
class SCFDriver {
public:
    SCFDriver(std::shared_ptr<FockBuilder> builder) 
        : builder_(builder) {}
    
    double solve() {
        int iter = 0;
        bool converged = false;
        double E_prev = 0.0;
        
        while (!converged && iter < max_iter_) {
            builder_->build_density();
            builder_->apply_kernel();
            builder_->diagonalize_fock();
            
            double E = builder_->compute_energy() + E_nuc_;
            double dE = E - E_prev;
            
            converged = (std::abs(dE) < e_thresh_);
            
            E_prev = E;
            iter++;
        }
        
        return E_prev;
    }
    
private:
    std::shared_ptr<FockBuilder> builder_;
    double E_nuc_;
    double e_thresh_ = 1e-8;
    int max_iter_ = 100;
};
```
**Test:**
```cpp
auto builder = std::make_shared<RHFFockBuilder>();
builder->initialize(...);

SCFDriver driver(builder);
double E = driver.solve();
// Compare with HFNew result
```
**Acceptance:**
- [ ] SCFDriver compiles
- [ ] Simple RHF calculation works
- [ ] Energy matches HFNew to 1e-10
**Time:** ~3 hours
**Commit:** "Add SCFDriver using FockBuilder (basic loop, no DIIS)"

---

### Step 1.8: Wire SCFDriver to scf_new
**What:** Make scf_new use SCFDriver + RHFFockBuilder
**Why:** End-to-end test of new architecture
**Changes:**
```cpp
// In HFNew::compute_energy()
if (use_new_driver_) {
    auto builder = std::make_shared<RHFFockBuilder>();
    builder->initialize(basis_, jk_, H_, S_, nocc_);
    
    auto driver = std::make_shared<SCFDriver>(builder);
    return driver->solve();
} else {
    // Old path (original HF code)
    return old_scf_loop();
}
```
**Test:**
```python
E_new = psi4.energy('scf_new')  # Uses SCFDriver
E_old = psi4.energy('scf')      # Uses original HF
assert abs(E_new - E_old) < 1e-8  # May differ slightly (no DIIS yet)
```
**Acceptance:**
- [ ] scf_new uses SCFDriver
- [ ] Energy within 1e-8 of original (no DIIS okay for now)
- [ ] All baseline tests pass (with relaxed tolerance)
**Time:** ~2 hours
**Commit:** "Wire scf_new to use SCFDriver+RHFFockBuilder"

---

**Phase 1 Complete Checklist:**
- [ ] FockBuilder interface defined
- [ ] RHFFockBuilder implements all methods
- [ ] SCFDriver basic loop works
- [ ] scf_new uses new architecture
- [ ] Baseline tests pass (within 1e-8)

**Time:** ~12 hours
**Risk:** LOW-MEDIUM - isolated changes, old path still works

---

## Phase 2: Add Accelerators & Converge to Baseline (7 micro-steps)

**Goal:** Add DIIS, achieve exact equivalence with original SCF.

### Step 2.1: Extract DIIS as DIISAccelerator
**What:** Create standalone DIIS class
**Files:** Create `psi4/src/psi4/libscf_solver/diis_accelerator.h/cc`
**Code:**
```cpp
class DIISAccelerator {
public:
    void accelerate(FockBuilder* builder, int iteration) {
        int n = builder->n_states();
        
        for (int s = 0; s < n; ++s) {
            auto F = builder->get_fock(s);
            auto D = builder->get_density(s);
            auto S = builder->get_overlap();
            
            auto error = compute_error(F, D, S);
            
            fock_history_[s].push_back(F);
            error_history_[s].push_back(error);
            
            if (fock_history_[s].size() >= min_vectors_) {
                auto F_extrap = extrapolate(s);
                builder->set_fock(F_extrap, s);
            }
        }
    }
    
private:
    std::vector<std::deque<SharedMatrix>> fock_history_;
    std::vector<std::deque<SharedMatrix>> error_history_;
    
    SharedMatrix compute_error(SharedMatrix F, SharedMatrix D, SharedMatrix S);
    SharedMatrix extrapolate(int state);
};
```
**Test:**
```cpp
// Unit test DIIS
DIISAccelerator diis;
// Run 5 iterations, check error decreases
```
**Acceptance:**
- [ ] DIIS compiles
- [ ] Unit test shows error reduction
**Time:** ~4 hours
**Commit:** "Add DIISAccelerator (standalone DIIS implementation)"

---

### Step 2.2: Integrate DIIS into SCFDriver
**What:** Add accelerator support to SCFDriver
**Code:**
```cpp
class SCFDriver {
public:
    void add_accelerator(std::shared_ptr<SCFAccelerator> acc) {
        accelerators_.push_back(acc);
    }
    
    double solve() {
        while (!converged) {
            builder_->build_density();
            builder_->apply_kernel();
            
            // Apply accelerators
            for (auto& acc : accelerators_) {
                acc->accelerate(builder_.get(), iter);
            }
            
            builder_->diagonalize_fock();
            // ...
        }
    }
    
private:
    std::vector<std::shared_ptr<SCFAccelerator>> accelerators_;
};
```
**Test:**
```cpp
auto driver = std::make_shared<SCFDriver>(builder);
auto diis = std::make_shared<DIISAccelerator>(8);
driver->add_accelerator(diis);
double E = driver->solve();
// Should converge in ~10 iterations like original
```
**Acceptance:**
- [ ] DIIS integration works
- [ ] Convergence in similar iterations as original
- [ ] Energy matches original to 1e-10
**Time:** ~2 hours
**Commit:** "Integrate DIIS into SCFDriver"

---

### Step 2.3: Add damping accelerator
**What:** Implement DampingAccelerator
**Time:** ~2 hours
**Commit:** "Add DampingAccelerator"

---

### Step 2.4: Add density convergence check
**What:** Implement RMS(ΔD) in convergence criteria
**Code:**
```cpp
double rms_D = builder_->compute_rms_density_change();
bool converged = (std::abs(dE) < e_thresh_) && (rms_D < d_thresh_);
```
**Time:** ~1 hour
**Commit:** "Add density convergence criterion to SCFDriver"

---

### Step 2.5: Match iteration printing to original
**What:** Make scf_new output identical to scf
**Time:** ~1 hour
**Commit:** "Match SCFDriver iteration output to original HF"

---

### Step 2.6: Achieve exact equivalence (< 1e-10)
**What:** Fine-tune until all baseline tests match to 1e-10
**Why:** Prove perfect equivalence
**Test:**
```python
for test in baseline_tests:
    E_old = psi4.energy('scf', **test)
    E_new = psi4.energy('scf_new', **test)
    assert abs(E_old - E_new) < 1e-10
```
**Acceptance:**
- [ ] All baseline tests < 1e-10 deviation
- [ ] Iteration counts match
- [ ] Output format matches
**Time:** ~3 hours (may need debugging)
**Commit:** "Achieve exact equivalence with original SCF (< 1e-10)"

---

### Step 2.7: Performance comparison
**What:** Benchmark scf vs scf_new
**Test:**
```python
import time
t0 = time.time()
E_old = psi4.energy('scf')
t_old = time.time() - t0

t0 = time.time()
E_new = psi4.energy('scf_new')
t_new = time.time() - t0

print(f"Old: {t_old:.2f}s, New: {t_new:.2f}s, Ratio: {t_new/t_old:.2f}")
```
**Acceptance:**
- [ ] scf_new within 1.2x of scf performance
- [ ] Document any regressions
**Time:** ~2 hours
**Commit:** "Add performance benchmarks (scf vs scf_new)"

---

**Phase 2 Complete Checklist:**
- [ ] DIIS accelerator working
- [ ] Exact equivalence achieved (< 1e-10)
- [ ] Performance acceptable (< 1.2x slowdown)
- [ ] All baseline tests pass

**Time:** ~15 hours
**Risk:** MEDIUM - DIIS is complex, may need debugging

---

## Phase 3: Add UHF Support (5 micro-steps)

**Goal:** Implement UHFFockBuilder and validate.

### Step 3.1: Create UHFFockBuilder skeleton
**What:** UHFFockBuilder with n_states() = 2
**Files:** Create `psi4/src/psi4/libscf_solver/uhf_fock_builder.h/cc`
**Code:**
```cpp
class UHFFockBuilder : public FockBuilder {
public:
    int n_states() const override { return 2; }
    
    SharedMatrix get_fock(int s) override {
        return (s == 0) ? F_alpha_ : F_beta_;
    }
    
    // Implement all methods for alpha/beta
};
```
**Time:** ~3 hours
**Commit:** "Add UHFFockBuilder skeleton (n_states=2)"

---

### Step 3.2: Implement UHF build_density
**Time:** ~2 hours
**Commit:** "Implement UHFFockBuilder::build_density"

---

### Step 3.3: Implement UHF apply_kernel
**What:** Form F_α and F_β with J_total
**Code:**
```cpp
void apply_kernel() override {
    jk_->C_left() = {C_alpha_, C_beta_};
    jk_->compute();
    
    auto J_total = jk_->J()[0]->clone();
    J_total->add(jk_->J()[1]);
    
    F_alpha_ = H_ + J_total - jk_->K()[0];
    F_beta_  = H_ + J_total - jk_->K()[1];
}
```
**Time:** ~2 hours
**Commit:** "Implement UHFFockBuilder::apply_kernel"

---

### Step 3.4: Wire UHF to scf_new
**What:** Detect reference=uhf and use UHFFockBuilder
**Time:** ~1 hour
**Commit:** "Add UHF support to scf_new"

---

### Step 3.5: Validate UHF baseline
**What:** Test UHF against original
**Test:**
```python
E_old = psi4.energy('scf', reference='uhf')
E_new = psi4.energy('scf_new', reference='uhf')
assert abs(E_old - E_new) < 1e-10
```
**Acceptance:**
- [ ] UHF baseline tests pass
- [ ] Exact equivalence achieved
**Time:** ~2 hours
**Commit:** "Validate UHF equivalence with baseline"

---

**Phase 3 Complete Checklist:**
- [ ] UHFFockBuilder working
- [ ] UHF baseline tests pass (< 1e-10)
- [ ] Both RHF and UHF work in scf_new

**Time:** ~10 hours
**Risk:** MEDIUM

---

## Phase 4: REKS Implementation (Future)

*(To be detailed after Phase 3 completion)*

**Overview:**
- Step 4.1: Create REKSFockBuilder skeleton
- Step 4.2: Implement sub-density formation
- Step 4.3: Implement SA-Fock formation
- Step 4.4: Add REKS-specific tests
- Step 4.5: Validate against literature

**Time:** ~20 hours
**Risk:** HIGH (new functionality)

---

## Continuous Practices Throughout

### After EVERY micro-step:

**Claude does:**
1. **Modify code** according to step requirements
2. **Commit:**
   ```bash
   git add -A
   git commit -m "<descriptive message from plan>"
   git push
   ```
3. **Inform user:** "Step X.Y complete, ready for testing"

**User does:**
1. **Pull changes:**
   ```bash
   git pull origin claude/improve-scf-code-refactor-011CUySHMdVr9aeEhprxf6EJ
   ```
2. **Compile:**
   ```bash
   cd build && ninja psi4
   ```
3. **Run test suite** (user's comprehensive tests)
4. **Report results:** "✅ Success" or "❌ Error: [details]"

**Acceptance thresholds per phase:**
- Phase 0-1: Compilation success, no crashes
- Phase 1: Energy deviation < 1e-8
- Phase 2-4: Energy deviation < 1e-10

---

## Success Metrics (User Validates)

### Phase 0:
- [x] DEBUG output working (Step 0.0) ✅
- [ ] HFNew compiles (Step 0.1)
- [ ] scf_new callable and equivalent to scf (Step 0.2)

### Phase 1:
- [ ] FockBuilder interface complete
- [ ] RHFFockBuilder working
- [ ] SCFDriver basic loop functional
- [ ] User confirms: energy deviation < 1e-8

### Phase 2:
- [ ] DIIS accelerator working
- [ ] User confirms: exact equivalence < 1e-10
- [ ] User confirms: performance < 1.2x slowdown

### Phase 3:
- [ ] UHFFockBuilder working
- [ ] User confirms: both RHF and UHF < 1e-10

### Phase 4:
- [ ] REKS compiles and runs
- [ ] User confirms: matches literature benchmarks

---

## Risk Mitigation

**If stuck >2 hours on a micro-step:**
1. Revert commit
2. Break into smaller steps
3. Ask for help / review architecture
4. Document blocker in ISSUES.md

**If tests fail after change:**
1. Check compilation warnings
2. Run single test with debug output
3. Compare intermediate values (D, F, energy)
4. Revert if can't fix in 30min

**If performance regresses >20%:**
1. Profile with perf/gprof
2. Check for unnecessary copies
3. Verify BLAS calls
4. Document in PERF.md

---

## Estimated Timeline

| Phase | Micro-steps | Time | Cumulative |
|-------|-------------|------|------------|
| 0 | 5 | 6h | 6h |
| 1 | 8 | 12h | 18h |
| 2 | 7 | 15h | 33h |
| 3 | 5 | 10h | 43h |
| **Total Phase 0-3** | **25** | **43h** | - |
| 4 (Future) | ~10 | ~20h | ~63h |

**Realistic:** ~1-2 weeks for Phase 0-3 (RHF+UHF equivalence)
**Stretch:** Add 1-2 weeks for Phase 4 (REKS)

---

## Next Steps

**Current Status:**
- ✅ Phase 0.0 complete (DEBUG output added to RHF, UHF, ROHF)
- ✅ Branch: `claude/improve-scf-code-refactor-011CUySHMdVr9aeEhprxf6EJ`
- 📍 **Next:** Phase 0, Step 0.1 - Copy HF to HFNew

**For Claude:**
1. Read workflow section carefully
2. Begin Phase 0, Step 0.1
3. Follow micro-steps sequentially (do NOT skip!)
4. Commit after each step with message from plan
5. Inform user when ready for testing
6. Wait for user validation before proceeding

**For User:**
1. After Claude informs "ready for testing"
2. Pull changes: `git pull origin claude/improve-scf-code-refactor-011CUySHMdVr9aeEhprxf6EJ`
3. Compile: `cd build && ninja psi4`
4. Run your test suite
5. Report results to Claude
6. Proceed to next step if success

**Good luck! 🚀**
