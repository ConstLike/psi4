# Multi-Cycle SCF Architecture Skeleton

**Goal:** Universal SCF engine supporting RHF → ROHF → SA-REKS(N,M,S) with multiple independent cycles.

**Design Philosophy:**
- SCF engine = **abstract algorithm** (iteration, convergence, DIIS)
- Theory = **concrete implementation** (how to build D, F, check convergence)
- Clear separation of concerns → easy code generation

---

## Core Abstraction Layers

```
┌─────────────────────────────────────────────────────┐
│         MultiCycleSCF (Coordinator)                  │
│  - Manages N independent SCF cycles                  │
│  - Coordinates shared resources (JK, grid)           │
│  - Handles inter-cycle dependencies (optional)       │
└──────────────┬──────────────────────────────────────┘
               │ manages
               ↓
┌─────────────────────────────────────────────────────┐
│           SCFCycle (Single Cycle)                    │
│  - Generic SCF iteration algorithm                   │
│  - Uses SCFTheory for theory-specific operations     │
│  - Handles convergence, DIIS, damping                │
└──────────────┬──────────────────────────────────────┘
               │ uses
               ↓
┌─────────────────────────────────────────────────────┐
│          SCFTheory (Abstract Interface)              │
│  - Defines what SCF needs from theory                │
│  - Theory-specific: n_states, build_D, build_F, etc.│
└──────────────┬──────────────────────────────────────┘
               │ implements
       ┌───────┴────────┬──────────────┬──────────────┐
       ↓                ↓              ↓              ↓
   RHFTheory      UHFTheory      SAREKSTheory   CustomTheory
   (n=1)          (n=2)          (n=N, spin=S)  (user-defined)
```

---

## 1. SCFTheory Interface (Abstract Base)

**What SCF needs from any theory:**

```pseudocode
ABSTRACT CLASS SCFTheory:

    // ========== CONFIGURATION ==========
    ABSTRACT METHOD n_states() -> int
        // How many Fock matrices does this theory need?
        // RHF: 1, UHF: 2, SA-REKS(5 states): 5

    ABSTRACT METHOD spin_type() -> enum {RESTRICTED, UNRESTRICTED, ENSEMBLE}
        // What kind of spin treatment?

    ABSTRACT METHOD state_labels() -> list[string]
        // ["S0", "S1", "S2"] or ["T1", "T2"] or ["Alpha", "Beta"]

    // ========== DENSITY BUILDING ==========
    ABSTRACT METHOD build_densities(orbitals, occupations) -> list[Matrix]
        // Build N density matrices from orbitals
        // Input: orbitals = [C[0], C[1], ..., C[n_states-1]]
        //        occupations = theory-specific occupation pattern
        // Output: [D[0], D[1], ..., D[n_states-1]]

    ABSTRACT METHOD ensemble_density(densities, weights) -> Matrix
        // For DFT: build ensemble density for V_xc
        // Input: densities = [D[0], ..., D[N-1]]
        //        weights = [w[0], ..., w[N-1]]  (theory-specific)
        // Output: D_ensemble = Σ w[i] * D[i]
        // RHF/UHF: return single density or None

    // ========== FOCK BUILDING ==========
    ABSTRACT METHOD build_focks(H, J_list, K_list, V_xc) -> list[Matrix]
        // Build N Fock matrices from JK results
        // Input: H = core Hamiltonian
        //        J_list = [J[0], ..., J[n_states-1]]
        //        K_list = [K[0], ..., K[n_states-1]]
        //        V_xc = DFT potential (or None for HF)
        // Output: [F[0], F[1], ..., F[n_states-1]]

        // Examples:
        // RHF:     F[0] = H + 2*J[0] - K[0]
        // UHF:     F[0] = H + (J[0]+J[1]) - K[0]
        //          F[1] = H + (J[0]+J[1]) - K[1]
        // SAREKS:  F[i] = H + J[i] + V_xc - K[i]  (state-specific J!)

    // ========== ORBITAL UPDATE ==========
    ABSTRACT METHOD diagonalize_focks(focks) -> (orbitals, eigenvalues)
        // Diagonalize N Fock matrices
        // Input: [F[0], F[1], ..., F[n_states-1]]
        // Output: orbitals = [C[0], C[1], ..., C[n_states-1]]
        //         eigenvalues = [eps[0], eps[1], ..., eps[n_states-1]]

    ABSTRACT METHOD assign_occupations(eigenvalues) -> occupations
        // Determine occupation pattern for next iteration
        // Theory-specific logic (Aufbau, fractional, excited-state, etc.)

    // ========== CONVERGENCE ==========
    ABSTRACT METHOD compute_energy(densities, H, focks) -> float
        // Compute total energy for this theory
        // Theory-specific formula

    ABSTRACT METHOD convergence_criteria() -> dict
        // Return theory-specific convergence thresholds
        // Example: {"delta_E": 1e-8, "RMS_D": 1e-6, "max_gradient": 1e-5}

    ABSTRACT METHOD check_convergence(old_state, new_state) -> bool
        // Check if THIS theory has converged
        // old_state/new_state contain energies, densities, gradients, etc.

END CLASS
```

---

## 2. SCFCycle (Single Cycle Engine)

**Generic SCF iteration logic - theory-independent:**

```pseudocode
CLASS SCFCycle:

    // ========== MEMBERS ==========
    theory: SCFTheory            // Abstract theory interface
    cycle_id: string             // "singlet", "triplet", etc.

    // State variables
    orbitals: list[Matrix]       // Current orbitals [C[0], ..., C[n-1]]
    densities: list[Matrix]      // Current densities [D[0], ..., D[n-1]]
    focks: list[Matrix]          // Current Fock matrices [F[0], ..., F[n-1]]
    energies: list[float]        // Per-state energies

    // Convergence tracking
    iteration: int
    converged: bool
    convergence_history: list[dict]

    // DIIS for THIS cycle
    diis_manager: DIISManager

    // ========== CONSTRUCTOR ==========
    METHOD __init__(theory: SCFTheory, cycle_id: string):
        self.theory = theory
        self.cycle_id = cycle_id
        self.converged = False
        self.iteration = 0

        // Allocate storage based on theory
        n = theory.n_states()
        self.orbitals = allocate_list(n)
        self.densities = allocate_list(n)
        self.focks = allocate_list(n)
        self.energies = allocate_list(n)

        // Create DIIS for N states (shared coefficients)
        self.diis_manager = DIISManager(n_matrices=n)

    // ========== INITIALIZATION ==========
    METHOD initialize(guess_orbitals):
        // Set initial guess
        self.orbitals = guess_orbitals

        // Build initial densities
        occupations = self.theory.assign_occupations(None)  // Initial occ
        self.densities = self.theory.build_densities(self.orbitals, occupations)

    // ========== CORE SCF STEP ==========
    METHOD iterate_step(H, J_list, K_list, V_xc):
        // One SCF iteration (called by MultiCycleSCF)

        // 1. Build Fock matrices
        self.focks = self.theory.build_focks(H, J_list, K_list, V_xc)

        // 2. Compute error vectors (for DIIS)
        errors = []
        FOR each F in self.focks:
            error = compute_error_vector(F, self.densities[i])
            errors.append(error)

        // 3. DIIS extrapolation (if enabled)
        IF self.iteration >= diis_start:
            self.diis_manager.add_entry(self.focks, errors)
            self.focks = self.diis_manager.extrapolate(self.focks)

        // 4. Diagonalize Fock matrices
        self.orbitals, eigenvalues = self.theory.diagonalize_focks(self.focks)

        // 5. Update occupations
        occupations = self.theory.assign_occupations(eigenvalues)

        // 6. Build new densities
        old_densities = self.densities
        self.densities = self.theory.build_densities(self.orbitals, occupations)

        // 7. Compute energy
        old_energy = self.energies[-1] IF self.energies ELSE 0.0
        new_energy = self.theory.compute_energy(self.densities, H, self.focks)
        self.energies.append(new_energy)

        // 8. Check convergence
        state = {
            "energy": new_energy,
            "delta_E": new_energy - old_energy,
            "densities": self.densities,
            "old_densities": old_densities,
            "errors": errors
        }
        self.converged = self.theory.check_convergence(old_state, state)

        self.iteration += 1

    // ========== ACCESSORS ==========
    METHOD get_densities() -> list[Matrix]:
        RETURN self.densities

    METHOD get_orbitals() -> list[Matrix]:
        RETURN self.orbitals

    METHOD is_converged() -> bool:
        RETURN self.converged

    METHOD get_energy() -> float:
        RETURN self.energies[-1]

END CLASS
```

---

## 3. MultiCycleSCF (Coordinator)

**Manages multiple independent cycles with shared resources:**

```pseudocode
CLASS MultiCycleSCF:

    // ========== MEMBERS ==========
    cycles: list[SCFCycle]         // List of independent cycles
    shared_resources: SharedResources  // JK builder, DFT grid, etc.

    // Global settings
    max_iterations: int
    global_convergence_policy: enum {ALL, ANY, CUSTOM}

    // ========== CONSTRUCTOR ==========
    METHOD __init__(cycles: list[SCFCycle], shared_resources):
        self.cycles = cycles
        self.shared_resources = shared_resources

    // ========== MAIN LOOP ==========
    METHOD run():
        PRINT "Starting Multi-Cycle SCF"
        PRINT "  Number of cycles: ", len(self.cycles)
        FOR cycle in self.cycles:
            PRINT "    - ", cycle.cycle_id, ": ", cycle.theory.n_states(), " states"

        WHILE NOT global_converged() AND iteration < max_iterations:

            // ===== PHASE 1: Collect densities from ALL cycles =====
            all_densities = []
            FOR cycle in self.cycles:
                IF NOT cycle.is_converged():  // Skip converged cycles
                    densities = cycle.get_densities()
                    all_densities.extend(densities)

            // ===== PHASE 2: Shared JK computation (ONE call!) =====
            all_orbitals = extract_orbitals_from_densities(all_densities)
            J_list, K_list = self.shared_resources.compute_JK(all_orbitals)

            // ===== PHASE 3: DFT potential (per-cycle ensemble) =====
            dft_potentials = {}  // Map cycle_id -> V_xc
            FOR cycle in self.cycles:
                IF NOT cycle.is_converged():
                    // Build ensemble density for THIS cycle
                    densities = cycle.get_densities()
                    weights = cycle.theory.get_ensemble_weights()
                    D_ens = cycle.theory.ensemble_density(densities, weights)

                    IF D_ens is not None:  // DFT calculation
                        V_xc = self.shared_resources.compute_V_xc(D_ens)
                        dft_potentials[cycle.cycle_id] = V_xc
                    ELSE:
                        dft_potentials[cycle.cycle_id] = None

            // ===== PHASE 4: Distribute results and iterate each cycle =====
            j_offset = 0
            k_offset = 0
            FOR cycle in self.cycles:
                IF NOT cycle.is_converged():
                    n = cycle.theory.n_states()

                    // Extract JK for THIS cycle
                    J_cycle = J_list[j_offset : j_offset+n]
                    K_cycle = K_list[k_offset : k_offset+n]
                    V_xc = dft_potentials[cycle.cycle_id]

                    // Run one iteration for THIS cycle
                    cycle.iterate_step(H, J_cycle, K_cycle, V_xc)

                    j_offset += n
                    k_offset += n

                    // Print progress
                    PRINT "  Cycle ", cycle.cycle_id, " iter ", cycle.iteration,
                          " E=", cycle.get_energy(),
                          " converged=", cycle.is_converged()

            // ===== PHASE 5: Check global convergence =====
            IF global_converged():
                PRINT "All cycles converged!"
                BREAK

        // Print final results
        PRINT "\nFinal Energies:"
        FOR cycle in self.cycles:
            PRINT "  ", cycle.cycle_id, ": ", cycle.get_energy()

    // ========== CONVERGENCE POLICY ==========
    METHOD global_converged() -> bool:
        IF global_convergence_policy == ALL:
            // All cycles must converge
            RETURN all(cycle.is_converged() for cycle in self.cycles)

        ELIF global_convergence_policy == ANY:
            // At least one cycle converged
            RETURN any(cycle.is_converged() for cycle in self.cycles)

        ELSE:  // CUSTOM
            // User-defined logic (e.g., weighted convergence)
            RETURN custom_convergence_check()

END CLASS
```

---

## 4. Concrete Theory Implementations

### 4.1 RHF Theory

```pseudocode
CLASS RHFTheory IMPLEMENTS SCFTheory:

    METHOD n_states() -> int:
        RETURN 1

    METHOD spin_type() -> enum:
        RETURN RESTRICTED

    METHOD state_labels() -> list[string]:
        RETURN ["RHF"]

    METHOD build_densities(orbitals, occupations) -> list[Matrix]:
        C = orbitals[0]
        n_occ = occupations[0]
        D = C[:, 0:n_occ] @ C[:, 0:n_occ].T
        RETURN [D]

    METHOD ensemble_density(densities, weights) -> Matrix:
        RETURN None  // No ensemble for RHF

    METHOD build_focks(H, J_list, K_list, V_xc) -> list[Matrix]:
        J = J_list[0]
        K = K_list[0]
        F = H + 2*J - K
        IF V_xc is not None:
            F += V_xc
        RETURN [F]

    METHOD diagonalize_focks(focks) -> (orbitals, eigenvalues):
        F = focks[0]
        C, eps = diagonalize(F)
        RETURN [C], [eps]

    METHOD assign_occupations(eigenvalues) -> occupations:
        // Aufbau principle
        n_elec = get_n_electrons()
        RETURN [n_elec / 2]  // Doubly occupied

    METHOD compute_energy(densities, H, focks) -> float:
        D = densities[0]
        F = focks[0]
        E = trace(D @ (H + F))
        RETURN E + nuclear_repulsion

    METHOD check_convergence(old_state, new_state) -> bool:
        delta_E = abs(new_state["delta_E"])
        RMS_D = compute_RMS(new_state["densities"][0] - old_state["densities"][0])
        RETURN delta_E < 1e-8 AND RMS_D < 1e-6

END CLASS
```

### 4.2 UHF Theory

```pseudocode
CLASS UHFTheory IMPLEMENTS SCFTheory:

    METHOD n_states() -> int:
        RETURN 2  // Alpha, Beta

    METHOD spin_type() -> enum:
        RETURN UNRESTRICTED

    METHOD state_labels() -> list[string]:
        RETURN ["Alpha", "Beta"]

    METHOD build_densities(orbitals, occupations) -> list[Matrix]:
        C_a = orbitals[0]
        C_b = orbitals[1]
        n_a = occupations[0]
        n_b = occupations[1]

        D_a = C_a[:, 0:n_a] @ C_a[:, 0:n_a].T
        D_b = C_b[:, 0:n_b] @ C_b[:, 0:n_b].T
        RETURN [D_a, D_b]

    METHOD ensemble_density(densities, weights) -> Matrix:
        RETURN None  // No ensemble for UHF (or return total D_a + D_b)

    METHOD build_focks(H, J_list, K_list, V_xc) -> list[Matrix]:
        J_a = J_list[0]
        J_b = J_list[1]
        K_a = K_list[0]
        K_b = K_list[1]

        J_total = J_a + J_b
        F_a = H + J_total - K_a
        F_b = H + J_total - K_b

        IF V_xc is not None:
            V_a, V_b = V_xc  // Spin-polarized DFT
            F_a += V_a
            F_b += V_b

        RETURN [F_a, F_b]

    METHOD diagonalize_focks(focks) -> (orbitals, eigenvalues):
        F_a = focks[0]
        F_b = focks[1]
        C_a, eps_a = diagonalize(F_a)
        C_b, eps_b = diagonalize(F_b)
        RETURN [C_a, C_b], [eps_a, eps_b]

    METHOD assign_occupations(eigenvalues) -> occupations:
        n_alpha = get_n_alpha()
        n_beta = get_n_beta()
        RETURN [n_alpha, n_beta]

    METHOD compute_energy(densities, H, focks) -> float:
        D_a = densities[0]
        D_b = densities[1]
        F_a = focks[0]
        F_b = focks[1]
        E = 0.5 * (trace(D_a @ (H + F_a)) + trace(D_b @ (H + F_b)))
        RETURN E + nuclear_repulsion

    METHOD check_convergence(old_state, new_state) -> bool:
        delta_E = abs(new_state["delta_E"])
        RMS_D_a = compute_RMS(new_state["densities"][0] - old_state["densities"][0])
        RMS_D_b = compute_RMS(new_state["densities"][1] - old_state["densities"][1])
        RMS_D = max(RMS_D_a, RMS_D_b)
        RETURN delta_E < 1e-8 AND RMS_D < 1e-6

END CLASS
```

### 4.3 SA-REKS Theory

```pseudocode
CLASS SAREKSTheory IMPLEMENTS SCFTheory:

    // Configuration
    n_electrons: int       // N
    n_orbitals: int        // M
    n_states_config: int   // Number of states
    spin: int              // 0=singlet, 1=triplet, 2=quintet, etc.
    weights: list[float]   // State-averaging weights

    METHOD __init__(n_electrons, n_orbitals, n_states, spin, weights):
        self.n_electrons = n_electrons
        self.n_orbitals = n_orbitals
        self.n_states_config = n_states
        self.spin = spin
        self.weights = weights

    METHOD n_states() -> int:
        RETURN self.n_states_config

    METHOD spin_type() -> enum:
        RETURN ENSEMBLE

    METHOD state_labels() -> list[string]:
        IF self.spin == 0:
            RETURN ["S0", "S1", "S2", ...][:self.n_states_config]
        ELIF self.spin == 1:
            RETURN ["T1", "T2", "T3", ...][:self.n_states_config]
        ELIF self.spin == 2:
            RETURN ["Q1", "Q2", ...][:self.n_states_config]

    METHOD build_densities(orbitals, occupations) -> list[Matrix]:
        // Build N state-specific densities
        densities = []
        FOR i in range(self.n_states_config):
            C_i = orbitals[i]
            occ_i = occupations[i]  // State-specific fractional occupations
            D_i = build_density_from_fractional_occ(C_i, occ_i)
            densities.append(D_i)
        RETURN densities

    METHOD ensemble_density(densities, weights) -> Matrix:
        // KEY: State-averaged density for DFT
        D_ensemble = zeros_like(densities[0])
        FOR i in range(len(densities)):
            D_ensemble += weights[i] * densities[i]
        RETURN D_ensemble

    METHOD build_focks(H, J_list, K_list, V_xc) -> list[Matrix]:
        // Build N Fock matrices
        focks = []
        FOR i in range(self.n_states_config):
            J_i = J_list[i]
            K_i = K_list[i]
            F_i = H + J_i - K_i

            IF V_xc is not None:
                F_i += V_xc  // Shared V_xc from ensemble density!

            focks.append(F_i)
        RETURN focks

    METHOD diagonalize_focks(focks) -> (orbitals, eigenvalues):
        // Diagonalize EACH Fock separately
        orbitals = []
        eigenvalues = []
        FOR i in range(self.n_states_config):
            C_i, eps_i = diagonalize(focks[i])
            orbitals.append(C_i)
            eigenvalues.append(eps_i)
        RETURN orbitals, eigenvalues

    METHOD assign_occupations(eigenvalues) -> occupations:
        // State-specific fractional occupations (REKS logic)
        // This is theory-specific and can be GENERATED from template
        occupations = []
        FOR i in range(self.n_states_config):
            occ_i = compute_reks_occupations(
                eigenvalues[i],
                self.n_electrons,
                self.n_orbitals,
                self.spin,
                state_index=i
            )
            occupations.append(occ_i)
        RETURN occupations

    METHOD compute_energy(densities, H, focks) -> float:
        // State-averaged energy
        E_total = 0.0
        FOR i in range(self.n_states_config):
            D_i = densities[i]
            F_i = focks[i]
            E_i = trace(D_i @ (H + F_i))
            E_total += self.weights[i] * E_i
        RETURN E_total + nuclear_repulsion

    METHOD check_convergence(old_state, new_state) -> bool:
        // Check convergence for ensemble
        delta_E = abs(new_state["delta_E"])

        // Max RMS across all states
        RMS_max = 0.0
        FOR i in range(self.n_states_config):
            RMS_i = compute_RMS(new_state["densities"][i] - old_state["densities"][i])
            RMS_max = max(RMS_max, RMS_i)

        RETURN delta_E < 1e-8 AND RMS_max < 1e-6

    METHOD get_ensemble_weights() -> list[float]:
        RETURN self.weights

END CLASS
```

---

## 5. Usage Examples

### Example 1: Simple RHF

```pseudocode
// Create theory
theory = RHFTheory()

// Create cycle
cycle = SCFCycle(theory, "RHF")
cycle.initialize(guess_orbitals)

// Create coordinator with one cycle
coordinator = MultiCycleSCF([cycle], shared_resources)

// Run
coordinator.run()
```

### Example 2: UHF + ROHF (two independent cycles)

```pseudocode
// Create theories
uhf_theory = UHFTheory()
rohf_theory = ROHFTheory()

// Create cycles
uhf_cycle = SCFCycle(uhf_theory, "UHF")
rohf_cycle = SCFCycle(rohf_theory, "ROHF")

// Initialize both
uhf_cycle.initialize(guess_uhf)
rohf_cycle.initialize(guess_rohf)

// Create coordinator with TWO cycles
coordinator = MultiCycleSCF([uhf_cycle, rohf_cycle], shared_resources)

// Run - shared JK for both!
coordinator.run()
```

### Example 3: SA-REKS with multiple spin states

```pseudocode
// SA-REKS(10,6,S=0): 10 electrons, 6 active orbitals, 5 singlet states
singlet_theory = SAREKSTheory(
    n_electrons=10,
    n_orbitals=6,
    n_states=5,
    spin=0,
    weights=[0.4, 0.3, 0.2, 0.05, 0.05]
)
singlet_cycle = SCFCycle(singlet_theory, "Singlet-REKS")

// SA-REKS(10,6,S=1): 10 electrons, 6 active orbitals, 3 triplet states
triplet_theory = SAREKSTheory(
    n_electrons=10,
    n_orbitals=6,
    n_states=3,
    spin=1,
    weights=[0.6, 0.3, 0.1]
)
triplet_cycle = SCFCycle(triplet_theory, "Triplet-REKS")

// SA-REKS(10,6,S=2): 10 electrons, 6 active orbitals, 1 quintet state
quintet_theory = SAREKSTheory(
    n_electrons=10,
    n_orbitals=6,
    n_states=1,
    spin=2,
    weights=[1.0]
)
quintet_cycle = SCFCycle(quintet_theory, "Quintet-REKS")

// Initialize all
singlet_cycle.initialize(guess_singlet)
triplet_cycle.initialize(guess_triplet)
quintet_cycle.initialize(guess_quintet)

// Create coordinator with THREE independent cycles!
coordinator = MultiCycleSCF(
    [singlet_cycle, triplet_cycle, quintet_cycle],
    shared_resources
)

// Run - ONE JK call for ALL (5+3+1=9 states)!
coordinator.run()

// Results:
//   Singlet-REKS: E = -100.234 Ha
//   Triplet-REKS: E = -99.876 Ha
//   Quintet-REKS: E = -99.345 Ha
```

---

## 6. Code Generation Pattern

**For code generator:**

```pseudocode
// Template for SA-REKS(N, M, S)
FUNCTION generate_sareks_theory(N, M, states_per_spin):

    code = """
    CLASS SAREKSTheory_{N}e_{M}orb IMPLEMENTS SCFTheory:

        METHOD n_states() -> int:
            RETURN {total_states}

        METHOD assign_occupations(eigenvalues) -> occupations:
            // GENERATED occupation pattern for N={N}, M={M}
            {generated_occupation_logic}

        // ... rest of methods follow template
    """

    // Fill in template with specific N, M, S configuration
    generated_occupation_logic = generate_occupation_pattern(N, M, states_per_spin)
    total_states = sum(states_per_spin.values())

    RETURN code

// Example usage:
config = {
    "spin_0": 5,  // 5 singlet states
    "spin_1": 3,  // 3 triplet states
    "spin_2": 1   // 1 quintet state
}

theory_code = generate_sareks_theory(N=10, M=6, states_per_spin=config)
```

---

## 7. Key Architectural Decisions

### What is ABSTRACT (Theory-specific):
1. ✅ Number of states (`n_states()`)
2. ✅ Density building from orbitals (`build_densities()`)
3. ✅ Fock matrix construction (`build_focks()`)
4. ✅ Occupation assignment (`assign_occupations()`)
5. ✅ Energy computation (`compute_energy()`)
6. ✅ Convergence criteria (`check_convergence()`)
7. ✅ Ensemble density (if needed) (`ensemble_density()`)

### What is COMMON (SCF engine):
1. ✅ Iteration loop (Phase 1-5 in MultiCycleSCF.run())
2. ✅ DIIS management (DIISManager, shared per cycle)
3. ✅ Convergence tracking (iteration count, history)
4. ✅ JK coordination (collect orbitals → compute → distribute)
5. ✅ DFT grid management (compute V_xc per cycle)
6. ✅ Inter-cycle orchestration (run multiple cycles together)

### Benefits of this architecture:
1. ✅ **Clear separation:** SCF logic vs theory logic
2. ✅ **Easy to generate:** SAREKSTheory template + fill parameters
3. ✅ **Extensible:** Add new theory by implementing interface
4. ✅ **Testable:** Each theory can be tested independently
5. ✅ **HPC-optimized:** Shared JK for all cycles (1.8-2x speedup)
6. ✅ **Flexible:** Can run any combination of theories together

---

## Next Steps

1. **Review this skeleton** - Is the abstraction clear?
2. **Agree on interface** - Any missing methods in SCFTheory?
3. **Map to Psi4** - Align with actual Psi4 classes/functions
4. **Implement prototype** - Start with RHF/UHF theories
5. **Code generator** - Create template for SA-REKS(N,M,S)

**Ready to proceed?** 🚀
