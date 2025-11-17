# Root Cause Analysis - 3x Verification
## Bug: UHF +9 Extra Iterations with ROHF

### Observed Data (FACTS)
```
Independent UHF: 7 iterations → E = -74.362669 Ha
multi_scf([UHF only]): 8 iterations (+1 vs independent)
multi_scf([UHF, ROHF]): 16 iterations (+9 vs independent)
```

**Grace pattern debug output:**
```
Iter 6: ROHF CONVERGED, will freeze on next iteration
Iter 7: GRACE ITERATION, freezing 2 CONVERGED C matrices
Iter 8+: FROZEN, using 2 frozen C matrices
```

**Energy timeline:**
```
Iter    UHF Energy         ROHF Energy        Status
1       -73.608110        -73.608110         ----/----
2       -74.349002        -74.348128         ----/----
3       -74.362430        -74.361534         ----/----
4       -74.362606        -74.361561         ----/----
5       -74.361053 ←DIVERGE  -74.361562       ----/----
6       -74.359090        -74.361562 ←CONV   ----/CONV
7       -74.359455        -74.361562         ----/CONV ← Grace
8-15    oscillating       -74.361562         ----/CONV
16      -74.362669 ←CONV  -74.361562         CONV/CONV
```

### First Analysis - When Does Problem Occur?

**Critical observation:** UHF diverges on iteration 5, BEFORE ROHF converges (iter 6).

Grace pattern cannot prevent this because:
- Grace starts at iter 7 (after ROHF convergence)
- Damage occurs at iter 5 (before ROHF convergence)

**Question 1:** Why does UHF diverge on iteration 5?
- Iter 4: UHF = -74.362606 Ha (reasonable)
- Iter 5: UHF = -74.361053 Ha (jump of 0.0015 Ha!)
- This is BEFORE any freezing occurs

**Question 2:** Is there physical coupling between UHF and ROHF?

Need to verify: Does JK compute create physical coupling, or is it just computational batching?

### Second Analysis - Understanding JK Coupling

**Code review (lines 1379-1388):**
```python
jk.C_clear()
for C_occ in all_C_occ_matrices:  # C from ALL wfn
    jk.C_add(C_occ)
jk.compute()  # Single call
```

**Code review (lines 1409-1416):**
```python
for i, wfn in enumerate(wfn_list):
    J_list = [J_all[jk_index + j] for j in range(n_states)]
    K_list = [K_all[jk_index + j] for j in range(n_states)]
    wfn.set_jk_matrices(J_list, K_list)
    jk_index += n_states  # Each wfn gets its OWN J/K
```

**CRITICAL QUESTION:** When JK object computes with multiple C matrices, does it:
- A) Compute J[i] from C[i] ONLY (independent)
- B) Compute J[i] from ALL C matrices (physical coupling)

**Need to verify JK.compute() behavior!**

If (A): No physical coupling → +9 iterations is a BUG in our code
If (B): Physical coupling → +9 iterations may be EXPECTED

### Third Analysis - Comparing with multi_scf([UHF only])

**Key data point:** multi_scf([UHF only]) = 8 iterations (+1 vs independent)

This tells us:
1. multi_scf framework adds +1 iteration overhead (acceptable)
2. The remaining +8 iterations (16 - 8 = 8) are specifically from ROHF coupling

**Breakdown:**
- Framework overhead: +1 iteration
- ROHF coupling effect: +8 iterations
- Total: +9 iterations

**Question:** What causes the +8 from ROHF coupling if there's no physical coupling?

### Hypothesis Testing

**Hypothesis 1:** Physical coupling through JK
- If J[UHF] includes contributions from D[ROHF], this is physical coupling
- This would create a coupled optimization surface
- DIIS may fail when D[ROHF] changes

**Test:** Check JK.compute() documentation/code to verify behavior

**Hypothesis 2:** Numerical/algorithmic interference
- Even without physical coupling, sharing JK object may create artifacts
- Convergence checks, DIIS management, etc. may interfere

**Test:** Run multi_scf([UHF, ROHF]) but with COMPLETELY separate JK objects

**Hypothesis 3:** Grace pattern timing issue
- Freezing after convergence is too late
- Need to freeze EARLIER (before convergence)

**Test:** Try "early freeze" when wfn is "close" to convergence (e.g., ΔE < 10*threshold)

### Action Items - Before Implementation

**VERIFY (do not assume!):**

1. ✓ Grace pattern freezes CONVERGED C (not pre-converged)
   - Verified: Iteration 7 gets C from AFTER form_C() on iter 6

2. ? How does JK.compute() work with multiple C matrices?
   - **MUST READ**: JK source code to understand compute() behavior
   - Does it couple densities or compute independently?

3. ? Why does multi_scf([UHF only]) add +1 iteration?
   - Need to identify specific difference vs independent SCF
   - Options? Initialization? Convergence checks?

4. ? Can we reproduce with different molecules?
   - Test on H2O, NH3, etc.
   - Is this OH-specific or general?

5. ? What if we run multi_scf([ROHF, UHF]) (reversed order)?
   - Does UHF still diverge?
   - Or is it order-dependent?

### Next Steps (IN ORDER)

**Step 1:** READ JK.compute() source code
- File: psi4/src/psi4/libfock/jk.cc (or similar)
- Understand: Does J[i] depend on C[i] only, or ALL C matrices?
- Answer definitively: Is there physical coupling?

**Step 2:** If NO physical coupling:
- Debug why UHF behaves differently when paired with ROHF
- Check for algorithmic interference
- May be a BUG in our multi_scf implementation

**Step 3:** If YES physical coupling:
- This is fundamental to multi-cycle approach
- May need to accept extra iterations
- OR: Decouple after convergence (run separately)

**Step 4:** Test alternative approaches
- Try resetting DIIS at grace iteration
- Try early freeze (before full convergence)
- Try micro-iterations

**DO NOT implement anything until Step 1-2 are complete!**

### Critical Questions to Answer

1. **Does JK.compute([C1, C2, ...]) create physical coupling?**
   - YES: J[i] depends on ALL densities → explains divergence
   - NO: J[i] depends on C[i] only → this is a bug

2. **Why does multi_scf([UHF only]) add +1 iteration?**
   - Framework overhead
   - Options difference
   - Other?

3. **Can we eliminate coupling while keeping batching?**
   - Compute J/K independently but batch for efficiency?
   - Or is coupling fundamental to shared JK?

### Conclusion of Analysis

**Before implementing ANY solution:**
1. Read JK source code (Step 1)
2. Understand coupling mechanism (Step 2)
3. Identify root cause definitively (Step 3)
4. Design solution based on root cause (Step 4)
5. Verify solution logic 2x before coding (Step 5)
6. Test on multiple molecules (Step 6)

**DO NOT:**
- Guess at solutions
- Implement without understanding root cause
- Skip verification steps

This is a complex coupling issue. Need rigorous analysis before solutions.
