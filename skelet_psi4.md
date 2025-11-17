# Multi-Cycle SCF Architecture - Psi4 Implementation Plan

**Based on:** `skelet.md` (abstract design)
**Target:** Psi4 codebase with minimal disruption
**Strategy:** "Onion layers" - incremental, testable steps from core to full wrapper

---

## Current Psi4 Architecture (Layer 0 - BASELINE)

### What We Have:

**Main Loop:** Python (`psi4/driver/procrouting/scf_proc/scf_iterator.py`)
```python
def scf_iterate(wfn, ..., maxiter=50):
    while iteration < maxiter:
        wfn.form_G()       # Virtual: theory-specific
        wfn.form_F()       # Virtual: theory-specific

        # DIIS (Python-side)
        diis_manager.add_entry(...)
        diis_manager.extrapolate(...)

        wfn.form_C()       # Virtual: diagonalize
        wfn.form_D()       # Virtual: build density
        E = wfn.compute_E()  # Virtual: energy

        # Check convergence
        if converged: break
```

**Class Hierarchy:**
```
HF (abstract base) - /libscf_solver/hf.{h,cc}
├── RHF - rhf.{h,cc}
├── UHF - uhf.{h,cc}
├── ROHF - rohf.{h,cc}
└── CUHF - cuhf.{h,cc}
```

**Virtual Methods (Theory Interface Already Exists!):**
```cpp
class HF : public Wavefunction {
protected:
    // Theory-specific operations (ALREADY VIRTUAL!)
    virtual void form_D() = 0;      // Build density
    virtual void form_G() = 0;      // JK + XC
    virtual void form_F() = 0;      // Fock assembly
    virtual void form_C(double shift = 0.0) = 0;  // Diagonalize
    virtual double compute_E() = 0;  // Energy

    // Common infrastructure
    std::shared_ptr<JK> jk_;
    py::object diis_manager_;
    int iteration_;
    bool converged_;
};
```

**What We Already Have from Phase 0.3:**
```cpp
// UHF with MultiStateMatrix (contiguous storage for n=2)
class UHF : public HF {
    std::shared_ptr<MultiStateMatrix> D_multi_;  // n=2
    std::shared_ptr<MultiStateMatrix> F_multi_;
    std::shared_ptr<MultiStateMatrix> G_multi_;

    SharedMatrix Da_ = D_multi_->get(0);  // Views!
    SharedMatrix Db_ = D_multi_->get(1);
};
```

### Key Insight:

**Psi4 ALREADY HAS the theory abstraction we need!**
- Virtual methods = `SCFTheory` interface (in skeleton)
- HF base class = partial `SCFCycle` (missing multi-state coordination)
- Python loop = partial `MultiCycleSCF` (missing shared JK)

**We don't need to refactor the whole thing - just ADD layers on top!**

---

## "Onion Layer" Implementation Plan

### 🎯 Strategy:
- Each layer ADDS functionality without breaking existing code
- Each layer is testable with existing tests
- Each layer brings us closer to the skeleton design
- Final layer = full multi-cycle SA-REKS

---

## Layer 1: Extract Common SCF Logic (Foundation)

**Goal:** Create `SCFEngine` base that coordinates virtual methods
**Time:** 2-3 days
**Risk:** Low (non-breaking, opt-in)

### What to Create:

**File:** `psi4/src/psi4/libscf_solver/scf_engine.h`
```cpp
namespace psi { namespace scf {

// Base class for SCF iteration engines
class SCFEngine {
protected:
    // Reference to theory implementation
    HF* theory_;

    // Iteration control
    int max_iterations_;
    double e_convergence_;
    double d_convergence_;

    // State tracking
    int iteration_;
    bool converged_;
    std::vector<double> energy_history_;

public:
    SCFEngine(HF* theory, Options& options)
        : theory_(theory),
          max_iterations_(options.get_int("MAXITER")),
          e_convergence_(options.get_double("E_CONVERGENCE")),
          d_convergence_(options.get_double("D_CONVERGENCE")),
          iteration_(0),
          converged_(false) {}

    virtual ~SCFEngine() = default;

    // Main iteration method (default implementation)
    virtual int iterate() {
        while (iteration_ < max_iterations_ && !converged_) {
            iterate_step();
            check_convergence();
            iteration_++;
        }
        return iteration_;
    }

protected:
    // Single iteration (calls virtual methods on theory_)
    virtual void iterate_step() {
        theory_->form_G();
        theory_->form_F();
        theory_->form_C();
        theory_->form_D();
    }

    // Convergence check (theory-independent)
    virtual void check_convergence() {
        double E_new = theory_->compute_E();

        if (iteration_ > 0) {
            double delta_E = std::abs(E_new - energy_history_.back());
            // RMS_D computed from densities...

            if (delta_E < e_convergence_ /* && RMS_D < d_convergence_ */) {
                converged_ = true;
            }
        }

        energy_history_.push_back(E_new);
    }

public:
    bool is_converged() const { return converged_; }
    int get_iteration() const { return iteration_; }
};

}} // namespace psi::scf
```

### Integration with Existing Code:

**No changes to RHF/UHF!** Just add optional usage:

```cpp
// In Python or C++ driver:
auto engine = std::make_shared<SCFEngine>(rhf.get(), options);
int iters = engine->iterate();
// OR keep using existing scf_iterate() in Python
```

### Test:
- Run existing RHF/UHF tests → should pass (no code change)
- Add test using new `SCFEngine` → should match old results

**Outcome:** Proof-of-concept that we can wrap existing virtual methods in engine.

---

## Layer 2: Multi-State Support in SCFEngine

**Goal:** Allow engine to handle N states (extend MultiStateMatrix)
**Time:** 2-3 days
**Risk:** Low (RHF/UHF still use n=1,2; ready for n>2)

### What to Add:

**File:** `psi4/src/psi4/libscf_solver/scf_engine.h` (extend Layer 1)
```cpp
class MultiStateSCFEngine : public SCFEngine {
protected:
    int n_states_;

    // Per-state tracking
    std::vector<bool> state_converged_;
    std::vector<double> state_energies_;
    std::vector<double> state_delta_E_;

public:
    MultiStateSCFEngine(HF* theory, int n_states, Options& options)
        : SCFEngine(theory, options),
          n_states_(n_states),
          state_converged_(n_states, false),
          state_energies_(n_states, 0.0),
          state_delta_E_(n_states, 1e10) {}

    // Override for multi-state
    void check_convergence() override {
        // Check each state independently
        // For RHF/UHF: n_states = 1 or 2 (backward compatible)
        // For SAREKS: n_states = N
    }
};
```

### Extend RHF/UHF to Report n_states:

**Add to HF base:**
```cpp
class HF : public Wavefunction {
public:
    // New: theory reports how many states it handles
    virtual int n_states() const { return 1; }  // Default: RHF
};

class UHF : public HF {
public:
    int n_states() const override { return 2; }  // Alpha, Beta
};
```

### Test:
- RHF with `MultiStateSCFEngine(rhf, n=1)` → works
- UHF with `MultiStateSCFEngine(uhf, n=2)` → works
- Prepare for n>2 (no theory yet, but engine ready)

**Outcome:** Engine can handle arbitrary N states (tested with n=1,2).

---

## Layer 3: Cycle Coordinator (Multi-Cycle Support)

**Goal:** Manage multiple independent cycles with shared JK
**Time:** 3-4 days
**Risk:** Medium (new concept, but testable with 1 cycle first)

### What to Create:

**File:** `psi4/src/psi4/libscf_solver/multi_cycle_scf.h`
```cpp
class MultiCycleSCF {
private:
    // Cycles (each is an independent theory)
    std::vector<std::shared_ptr<MultiStateSCFEngine>> cycles_;

    // Shared JK for all cycles
    std::shared_ptr<JK> shared_jk_;

    // Global settings
    int max_iterations_;

public:
    // Add a cycle
    void add_cycle(std::shared_ptr<MultiStateSCFEngine> cycle) {
        cycles_.push_back(cycle);
    }

    // Main loop
    void run() {
        for (int iter = 0; iter < max_iterations_; ++iter) {
            // Phase 1: Collect all orbitals from all cycles
            std::vector<SharedMatrix> all_C_occ;
            for (auto& cycle : cycles_) {
                auto C = cycle->get_orbitals();
                all_C_occ.insert(all_C_occ.end(), C.begin(), C.end());
            }

            // Phase 2: Shared JK (ONE call for ALL cycles!)
            shared_jk_->C_left() = all_C_occ;
            shared_jk_->compute();

            const auto& J = shared_jk_->J();
            const auto& K = shared_jk_->K();

            // Phase 3: Distribute JK results and iterate each cycle
            int j_offset = 0;
            for (auto& cycle : cycles_) {
                int n = cycle->n_states();

                // Extract JK for this cycle
                std::vector<SharedMatrix> J_cycle(J.begin() + j_offset,
                                                   J.begin() + j_offset + n);
                std::vector<SharedMatrix> K_cycle(K.begin() + j_offset,
                                                   K.begin() + j_offset + n);

                // Set JK for this cycle
                cycle->set_JK(J_cycle, K_cycle);

                // Iterate this cycle (form_F, form_C, form_D, etc.)
                cycle->iterate_step();

                j_offset += n;
            }

            // Phase 4: Check global convergence
            if (all_converged()) break;
        }
    }

private:
    bool all_converged() const {
        for (auto& cycle : cycles_) {
            if (!cycle->is_converged()) return false;
        }
        return true;
    }
};
```

### Modification to MultiStateSCFEngine:

Add methods for JK injection:
```cpp
class MultiStateSCFEngine : public SCFEngine {
public:
    // For multi-cycle coordination
    std::vector<SharedMatrix> get_orbitals() const;
    void set_JK(const std::vector<SharedMatrix>& J,
                const std::vector<SharedMatrix>& K);

protected:
    // Modified iterate_step (skip form_G, use injected JK)
    void iterate_step() override {
        // Use pre-computed J/K from coordinator
        theory_->form_F();  // Uses J/K set by coordinator
        theory_->form_C();
        theory_->form_D();
    }
};
```

### Test:

**Test 1: Single cycle (backward compatible)**
```cpp
auto rhf = std::make_shared<RHF>(...);
auto engine = std::make_shared<MultiStateSCFEngine>(rhf.get(), 1, options);

MultiCycleSCF coordinator;
coordinator.add_cycle(engine);
coordinator.run();
// Should match existing RHF results
```

**Test 2: Two independent RHF cycles (proof of concept)**
```cpp
auto rhf1 = std::make_shared<RHF>(...);  // Molecule 1
auto rhf2 = std::make_shared<RHF>(...);  // Molecule 2

auto engine1 = std::make_shared<MultiStateSCFEngine>(rhf1.get(), 1, options);
auto engine2 = std::make_shared<MultiStateSCFEngine>(rhf2.get(), 1, options);

MultiCycleSCF coordinator;
coordinator.add_cycle(engine1);
coordinator.add_cycle(engine2);
coordinator.run();
// Both converge, shared JK used
```

**Outcome:** Multi-cycle infrastructure works, tested with existing theories.

---

## Layer 4: SA-REKS Theory Stub (Minimal Implementation)

**Goal:** Create simplest possible SA-REKS that fits in the framework
**Time:** 3-4 days
**Risk:** Medium (new theory, but follows existing pattern)

### What to Create:

**File:** `psi4/src/psi4/libscf_solver/sareks.h`
```cpp
class SAREKS : public HF {
private:
    int n_states_;
    int n_electrons_;
    int n_active_orbitals_;
    int spin_;  // 0=singlet, 1=triplet, etc.
    std::vector<double> weights_;

    // MultiStateMatrix for N states (extend Phase 0.3!)
    std::shared_ptr<MultiStateMatrix> D_multi_;
    std::shared_ptr<MultiStateMatrix> F_multi_;
    std::shared_ptr<MultiStateMatrix> G_multi_;

    // Views into contiguous storage (Phase 0.3 pattern!)
    std::vector<SharedMatrix> D_states_;  // [D[0], D[1], ..., D[n_states-1]]
    std::vector<SharedMatrix> F_states_;
    std::vector<SharedMatrix> G_states_;

public:
    SAREKS(int n_electrons, int n_orbitals, int n_states, int spin,
           const std::vector<double>& weights, SharedWavefunction ref, ...)
        : HF(ref, ...),
          n_states_(n_states),
          n_electrons_(n_electrons),
          n_active_orbitals_(n_orbitals),
          spin_(spin),
          weights_(weights) {

        // Phase 0.3: Contiguous storage for N states
        D_multi_ = std::make_shared<MultiStateMatrix>("D", n_states, nirrep_, nsopi_, nsopi_, 0);
        F_multi_ = std::make_shared<MultiStateMatrix>("F", n_states, nirrep_, nsopi_, nsopi_, 0);
        G_multi_ = std::make_shared<MultiStateMatrix>("G", n_states, nirrep_, nsopi_, nsopi_, 0);

        // Views (no copying!)
        for (int i = 0; i < n_states; ++i) {
            D_states_.push_back(D_multi_->get(i));
            F_states_.push_back(F_multi_->get(i));
            G_states_.push_back(G_multi_->get(i));
        }
    }

    // Report n_states to engine
    int n_states() const override { return n_states_; }

    // Theory-specific implementations
    void form_D() override {
        // Build N densities from orbitals (REKS occupation logic)
        for (int i = 0; i < n_states_; ++i) {
            // Simplified: proper REKS occupation pattern goes here
            build_density_for_state(i);
        }
    }

    void form_G() override {
        // Called by coordinator with pre-computed J/K for all states
        // Build G[i] = J[i] + V_xc - K[i]
        // V_xc from ensemble density
    }

    void form_F() override {
        // F[i] = H + G[i]
        for (int i = 0; i < n_states_; ++i) {
            F_states_[i]->copy(H_);
            F_states_[i]->add(G_states_[i]);
        }
    }

    void form_C() override {
        // Diagonalize EACH Fock separately
        for (int i = 0; i < n_states_; ++i) {
            diagonalize_F(F_states_[i], C_states_[i], epsilon_states_[i]);
        }
    }

    double compute_E() override {
        // State-averaged energy
        double E_total = 0.0;
        for (int i = 0; i < n_states_; ++i) {
            double E_i = compute_state_energy(i);
            E_total += weights_[i] * E_i;
        }
        return E_total + nuclearrep_;
    }
};
```

### Test:

**Minimal SA-REKS(2,2,S=0):** 2 electrons, 2 orbitals, singlet (should match RHF!)
```python
sareks = SAREKS(n_elec=2, n_orb=2, n_states=1, spin=0, weights=[1.0], ...)
engine = MultiStateSCFEngine(sareks, n_states=1, ...)

coordinator = MultiCycleSCF()
coordinator.add_cycle(engine)
coordinator.run()

# Should give RHF-like energy (sanity check)
```

**Outcome:** SA-REKS skeleton in place, ready for real REKS logic.

---

## Layer 5: Full SA-REKS Implementation (Generated Code)

**Goal:** Complete REKS occupation logic, support N states, all spins
**Time:** 1-2 weeks (with code generation)
**Risk:** High (complex REKS math, but framework proven)

### What to Implement:

1. **REKS occupation pattern** (state-specific, fractional)
2. **Proper ensemble density** for DFT
3. **State-specific convergence** criteria
4. **All spin cases** (singlet, triplet, quintet, ...)

### Code Generation Pattern:

**Template:** `sareks_template.cpp.j2` (Jinja2 or similar)
```cpp
// Generated for SA-REKS({{n_elec}},{{n_orb}},S={{spin}})
void SAREKS_{{n_elec}}e_{{n_orb}}orb_S{{spin}}::assign_occupations() {
    {% for state in range(n_states) %}
    // State {{state}}: occupation pattern
    occupations_[{{state}}] = { {{occupation_pattern[state]}} };
    {% endfor %}
}
```

**Generator:** `generate_sareks.py`
```python
def generate_sareks_theory(n_elec, n_orb, spin_configs):
    template = load_template("sareks_template.cpp.j2")

    for spin in spin_configs:
        occupation_pattern = compute_reks_occupations(n_elec, n_orb, spin)
        code = template.render(
            n_elec=n_elec,
            n_orb=n_orb,
            spin=spin,
            n_states=len(spin_configs[spin]),
            occupation_pattern=occupation_pattern
        )
        write_file(f"sareks_{n_elec}e_{n_orb}orb_S{spin}.cc", code)
```

### Test:
- SA-REKS literature benchmarks
- Compare with reference implementations

**Outcome:** Full SA-REKS working in multi-cycle framework.

---

## Layer 6: Multi-Spin-State Calculations (Final Goal!)

**Goal:** Run Singlet + Triplet + Quintet SA-REKS simultaneously
**Time:** 1-2 days (framework already done!)
**Risk:** Low (just combining existing pieces)

### Usage Example:

```python
# SA-REKS(10,6,S=0): 5 singlet states
singlet = SAREKS(n_elec=10, n_orb=6, n_states=5, spin=0,
                 weights=[0.4, 0.3, 0.2, 0.05, 0.05], ...)
singlet_engine = MultiStateSCFEngine(singlet, n_states=5, ...)

# SA-REKS(10,6,S=1): 3 triplet states
triplet = SAREKS(n_elec=10, n_orb=6, n_states=3, spin=1,
                 weights=[0.6, 0.3, 0.1], ...)
triplet_engine = MultiStateSCFEngine(triplet, n_states=3, ...)

# SA-REKS(10,6,S=2): 1 quintet state
quintet = SAREKS(n_elec=10, n_orb=6, n_states=1, spin=2,
                 weights=[1.0], ...)
quintet_engine = MultiStateSCFEngine(quintet, n_states=1, ...)

# Run all together!
coordinator = MultiCycleSCF()
coordinator.add_cycle(singlet_engine)
coordinator.add_cycle(triplet_engine)
coordinator.add_cycle(quintet_engine)
coordinator.run()

# Shared JK for 5+3+1=9 states! (~1.8x speedup)

print(f"Singlet energy: {singlet_engine.get_energy()}")
print(f"Triplet energy: {triplet_engine.get_energy()}")
print(f"Quintet energy: {quintet_engine.get_energy()}")
```

**Outcome:** Full multi-cycle, multi-spin SA-REKS working!

---

## Mapping to Abstract Skeleton

| Skeleton Component | Psi4 Implementation | Layer |
|--------------------|---------------------|-------|
| `SCFTheory` (interface) | `HF` virtual methods | **Layer 0** (exists!) |
| `RHFTheory`, `UHFTheory` | `RHF`, `UHF` classes | **Layer 0** (exists!) |
| `SCFCycle` (engine) | `SCFEngine` → `MultiStateSCFEngine` | **Layers 1-2** |
| `MultiCycleSCF` (coordinator) | `MultiCycleSCF` class | **Layer 3** |
| `SAREKSTheory` | `SAREKS` class | **Layers 4-5** |
| Multi-spin usage | Python driver | **Layer 6** |

---

## Timeline Summary

| Layer | Component | Time | Cumulative |
|-------|-----------|------|------------|
| 1 | `SCFEngine` (foundation) | 2-3 days | 3 days |
| 2 | Multi-state support | 2-3 days | 6 days |
| 3 | `MultiCycleSCF` coordinator | 3-4 days | 10 days |
| 4 | SA-REKS stub | 3-4 days | 14 days |
| 5 | Full SA-REKS (generated) | 1-2 weeks | ~1 month |
| 6 | Multi-spin integration | 1-2 days | ~1 month |

**Total:** ~1 month for full implementation

---

## Testing Strategy (Per Layer)

### Layer 1: SCFEngine
- Run all existing RHF/UHF/ROHF tests → must pass
- Compare `SCFEngine` results vs existing Python loop → must match

### Layer 2: MultiStateSCFEngine
- RHF with n=1 → match Layer 1
- UHF with n=2 → match Layer 1
- Ready for n>2 (no theory yet)

### Layer 3: MultiCycleSCF
- Single cycle → match Layer 2
- Two independent RHF → both converge
- Measure JK time: 2 cycles < 2× single cycle (proof of sharing)

### Layer 4: SA-REKS stub
- Minimal case (should match RHF/UHF)
- Convergence test (does it converge at all?)

### Layer 5: Full SA-REKS
- Literature benchmarks
- Compare with reference implementations

### Layer 6: Multi-spin
- Independent runs: Singlet, Triplet, Quintet separately
- Combined run: all together → same energies
- Performance: combined < 3× slowest single (proof of sharing)

---

## Risk Mitigation

### If a layer fails:
1. **Roll back to previous layer** (each layer is working!)
2. **Debug in isolation** (each layer is testable)
3. **Ask for help** with theory-specific issues (REKS math)

### If performance regresses:
1. **Profile with perf/gprof**
2. **Check MultiStateMatrix usage** (Phase 0.3 should help!)
3. **Verify shared JK is actually working**

---

## Key Advantages of This Approach

1. ✅ **Non-breaking:** Each layer adds, doesn't change existing code
2. ✅ **Testable:** Each layer passes existing tests
3. ✅ **Incremental:** Can stop at any layer and still have working code
4. ✅ **Reuses Phase 0.3:** MultiStateMatrix already proven (+15.9%!)
5. ✅ **Follows Psi4 patterns:** Virtual methods, Python integration
6. ✅ **Prepares for generation:** Clean abstractions ready for templates

---

## Next Steps

1. **Review this plan** - Does the layering make sense?
2. **Adjust timeline** - More/less time per layer?
3. **Start Layer 1** - Create `SCFEngine` foundation
4. **Test each layer** - User compiles, tests, validates
5. **Iterate** - Adjust based on findings

**Ready to start Layer 1?** 🚀
