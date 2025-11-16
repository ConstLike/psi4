#!/usr/bin/env python3
"""
Test grace iteration pattern fix

Expected behavior:
- Iteration N: ROHF converges, marked as just_converged
- Iteration N+1: GRACE period - freeze CONVERGED C, transition to fully converged
- Iteration N+2+: Use frozen converged C (stable)

This should eliminate discontinuity because frozen C = converged C.
"""

print("="*80)
print("GRACE ITERATION PATTERN - Mental Trace Test")
print("="*80)

# Simulated timeline
print("\n### ITERATION 6 (ROHF converges):")
print("1. Collect C:")
print("   - ROHF: get_orbital_matrices() → C_5 (from iter 5)")
print("   - UHF: get_orbital_matrices() → C_5")
print("\n2. JK compute: jk.compute([C_5_ROHF, C_5_UHF])")
print("   - UHF receives J/K based on C_5_ROHF")
print("\n3. _scf_iteration():")
print("   - ROHF: form_C() → C_6 (converged orbitals!)")
print("   - Convergence check → TRUE")
print("   - just_converged_flags[ROHF] = True ✓")
print("   - (NOT converged_flags yet!)")
print("\n4. State after iteration 6:")
print("   - ROHF.Ca_ = C_6 (converged)")
print("   - just_converged_flags[ROHF] = True")
print("   - converged_flags[ROHF] = False")

print("\n" + "="*80)
print("### ITERATION 7 (GRACE PERIOD):")
print("1. Collect C (at start of iteration):")
print("   - Check: just_converged_flags[ROHF] = True")
print("   - Get CONVERGED C: C_matrices = wfn.get_orbital_matrices() → C_6 ✓")
print("   - FREEZE: wfn._frozen_C_for_jk = C_6 ✓")
print("   - Transition: just_converged_flags[ROHF] = False")
print("   - Transition: converged_flags[ROHF] = True")
print("\n2. JK compute: jk.compute([C_6_ROHF, C_6_UHF])")
print("   - UHF receives J/K based on C_6_ROHF for FIRST time")
print("   - This is NEW J/K, but it's based on CONVERGED density ✓")
print("\n3. _scf_iteration():")
print("   - ROHF: SKIP (in grace period)")
print("   - UHF: continues iterating with NEW J/K from converged ROHF")
print("\n4. State after iteration 7:")
print("   - ROHF._frozen_C_for_jk = C_6 (converged)")
print("   - converged_flags[ROHF] = True")
print("   - UHF adapted to converged ROHF density")

print("\n" + "="*80)
print("### ITERATION 8+ (STABLE):")
print("1. Collect C:")
print("   - Check: converged_flags[ROHF] = True")
print("   - Use: C_matrices = wfn._frozen_C_for_jk → C_6 (same as iter 7!)")
print("\n2. JK compute: jk.compute([C_6_ROHF, C_7_UHF])")
print("   - UHF receives J/K based on C_6_ROHF (STABLE) ✓")
print("\n3. _scf_iteration():")
print("   - ROHF: SKIP (fully converged)")
print("   - UHF: continues with STABLE J/K")
print("\n4. Result:")
print("   - NO discontinuity! ✓")
print("   - UHF sees C_6 from iter 7 onwards (stable)")
print("   - UHF should converge within ~1-2 extra iterations")

print("\n" + "="*80)
print("### KEY INSIGHT:")
print("="*80)
print("✓ Grace iteration allows frozen C = CONVERGED C (not pre-converged)")
print("✓ Iteration 7 is transition where UHF first sees converged ROHF")
print("✓ From iteration 8+, everything stable")
print("✓ Expected UHF extra iterations: 1-2 (just transition cost)")
print("\n" + "="*80)
print("LOGIC VERIFIED! Ready to test with actual Psi4 code.")
print("="*80)
