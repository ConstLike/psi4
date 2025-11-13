#!/usr/bin/env python
"""
Test script for multi_scf() coordinator - Step 1.4

This script tests the NEW multi-cycle SCF architecture by converging
two independent RHF calculations simultaneously with shared JK computation.

Expected behavior:
- Both RHF cycles converge correctly
- Energies match independent runs (within tolerance)
- Performance gain ~1.5-2x (shared JK saves time)
"""

import psi4
import time
from psi4.driver.procrouting.scf_proc.scf_iterator import multi_scf

# Set memory and output
psi4.set_memory('2 GB')
psi4.core.set_output_file('test_multi_scf.out', False)

print("\n" + "="*70)
print("  Multi-Cycle SCF Test - Step 1.4")
print("  Testing: 2 independent H2O RHF calculations")
print("="*70 + "\n")

# Define two identical H2O molecules (just for testing)
# In real SA-REKS, these would be different spin states
mol1_str = """
O  0.000000  0.000000  0.117083
H  0.000000  0.755453 -0.468333
H  0.000000 -0.755453 -0.468333
symmetry c1
"""

mol2_str = """
O  0.000000  0.000000  0.117083
H  0.000000  0.755453 -0.468333
H  0.000000 -0.755453 -0.468333
symmetry c1
"""

# Create molecules
mol1 = psi4.geometry(mol1_str)
mol2 = psi4.geometry(mol2_str)

# Set basis
basis = psi4.core.BasisSet.build(mol1, "ORBITAL", "cc-pvdz")

# Set options
psi4.set_options({
    'scf_type': 'df',
    'e_convergence': 1e-8,
    'd_convergence': 1e-6,
    'maxiter': 50,
    'print': 1
})

print("Basis: cc-pVDZ")
print("SCF Type: DF")
print("E convergence: 1e-8")
print("D convergence: 1e-6\n")

# Test 1: Independent runs (baseline)
print("-" * 70)
print("Test 1: Independent SCF runs (baseline)")
print("-" * 70 + "\n")

t0 = time.time()

print("Running SCF for molecule 1...")
E1_independent = psi4.energy('scf', molecule=mol1)
print(f"  Energy 1: {E1_independent:.12f}\n")

print("Running SCF for molecule 2...")
E2_independent = psi4.energy('scf', molecule=mol2)
print(f"  Energy 2: {E2_independent:.12f}\n")

t_independent = time.time() - t0
print(f"Independent runs total time: {t_independent:.3f} seconds\n")

# Test 2: Multi-cycle SCF (shared JK)
print("-" * 70)
print("Test 2: Multi-cycle SCF (shared JK)")
print("-" * 70 + "\n")

# Create fresh wavefunctions
wfn1 = psi4.core.RHF(psi4.core.Wavefunction.build(mol1, basis), psi4.core.get_global_option('REFERENCE'))
wfn2 = psi4.core.RHF(psi4.core.Wavefunction.build(mol2, basis), psi4.core.get_global_option('REFERENCE'))

# Initialize both
print("Initializing wavefunctions...")
wfn1.initialize()
wfn2.initialize()

# Run multi-cycle SCF
print("\nRunning multi_scf() with NEW architecture...\n")
print("NOTE: Requires C++ recompilation for Phase 0.6 pybind11 exports\n")
t0 = time.time()

energies_multi = multi_scf([wfn1, wfn2], e_conv=1e-8, d_conv=1e-6, verbose=True)

t_multi = time.time() - t0

print(f"\nMulti-cycle SCF total time: {t_multi:.3f} seconds\n")

# Test 3: Verify results
print("-" * 70)
print("Test 3: Verification")
print("-" * 70 + "\n")

E1_multi = energies_multi[0]
E2_multi = energies_multi[1]

diff1 = abs(E1_multi - E1_independent)
diff2 = abs(E2_multi - E2_independent)

print(f"Energy 1 (independent): {E1_independent:.12f}")
print(f"Energy 1 (multi-cycle): {E1_multi:.12f}")
print(f"Difference:             {diff1:.2e}\n")

print(f"Energy 2 (independent): {E2_independent:.12f}")
print(f"Energy 2 (multi-cycle): {E2_multi:.12f}")
print(f"Difference:             {diff2:.2e}\n")

# Performance comparison
speedup = t_independent / t_multi if t_multi > 0 else 0.0
print(f"Performance:")
print(f"  Independent: {t_independent:.3f} s")
print(f"  Multi-cycle: {t_multi:.3f} s")
print(f"  Speedup:     {speedup:.2f}x\n")

# Pass/fail criteria
TOLERANCE = 1e-8
pass1 = diff1 < TOLERANCE
pass2 = diff2 < TOLERANCE
pass_overall = pass1 and pass2

print("-" * 70)
if pass_overall:
    print("✓ TEST PASSED - Energies match within tolerance")
else:
    print("✗ TEST FAILED - Energy mismatch exceeds tolerance")
print("-" * 70 + "\n")

# Summary
print("\n" + "="*70)
print("  Summary")
print("="*70)
print(f"  Status: {'PASS ✓' if pass_overall else 'FAIL ✗'}")
print(f"  Speedup: {speedup:.2f}x")
print(f"  Energy agreement: {max(diff1, diff2):.2e}")
print("="*70 + "\n")

if not pass_overall:
    raise Exception("Multi-cycle SCF test FAILED - energy mismatch")

print("Multi-cycle SCF test completed successfully!\n")
