#!/usr/bin/env python3
"""
Test script for C matrix freeze pattern fix (Step 1.5.1)

This tests the bug where UHF required +8 extra iterations when running
with ROHF in multi_scf() due to discontinuous J/K from converged wfn.

Expected behavior AFTER fix:
- Both UHF and ROHF should converge within reasonable iterations
- UHF should NOT require +8 extra iterations after ROHF converges
"""

import psi4
import sys

# Setup
psi4.set_output_file('test_freeze_pattern.out', False)
psi4.set_memory('500 MB')
psi4.set_num_threads(1)

# OH radical molecule (from bug report)
mol = psi4.geometry('''
    0 2
    O  0.0  0.0  0.0
    H  0.0  0.0  0.97
    symmetry c1
''')

# Options
psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'df',
    'e_convergence': 1e-8,
    'd_convergence': 1e-6,
    'maxiter': 50,
    'print': 2,  # Enable debug output
})

print("="*80)
print("Testing C matrix freeze pattern fix")
print("="*80)

# Test 1: Independent UHF (baseline)
print("\n" + "="*80)
print("TEST 1: Independent UHF (baseline)")
print("="*80)

psi4.set_options({'reference': 'uhf'})
e_uhf_independent = psi4.energy('scf')
uhf_iterations_independent = psi4.core.variable('SCF ITERATIONS')

print(f"\nUHF (independent): {uhf_iterations_independent} iterations")
print(f"Energy: {e_uhf_independent:.14f} Ha")

# Test 2: Independent ROHF (baseline)
print("\n" + "="*80)
print("TEST 2: Independent ROHF (baseline)")
print("="*80)

psi4.set_options({'reference': 'rohf'})
e_rohf_independent = psi4.energy('scf')
rohf_iterations_independent = psi4.core.variable('SCF ITERATIONS')

print(f"\nROHF (independent): {rohf_iterations_independent} iterations")
print(f"Energy: {e_rohf_independent:.14f} Ha")

# Test 3: Multi-SCF (UHF + ROHF) - THE MAIN TEST
print("\n" + "="*80)
print("TEST 3: Multi-SCF (UHF + ROHF) - Testing freeze pattern")
print("="*80)

# Import multi_scf
sys.path.insert(0, '/home/user/psi4/psi4/driver/procrouting/scf_proc')
from scf_iterator import multi_scf, create_wavefunction

# Create wavefunctions
wfn_uhf = create_wavefunction(mol, 'HF', 'UHF')
wfn_rohf = create_wavefunction(mol, 'HF', 'ROHF')

# Run multi_scf
print("\nRunning multi_scf([UHF, ROHF])...")
print("Debug output will show when C matrices are frozen/active:\n")

energies = multi_scf([wfn_uhf, wfn_rohf], e_conv=1e-8, d_conv=1e-6)

# Get iteration counts
uhf_iterations_multi = wfn_uhf.iteration_
rohf_iterations_multi = wfn_rohf.iteration_

print("\n" + "="*80)
print("RESULTS SUMMARY")
print("="*80)

print(f"\nIndependent UHF: {uhf_iterations_independent} iterations, E = {e_uhf_independent:.14f} Ha")
print(f"Independent ROHF: {rohf_iterations_independent} iterations, E = {e_rohf_independent:.14f} Ha")
print(f"\nMulti-SCF UHF: {uhf_iterations_multi} iterations, E = {energies[0]:.14f} Ha")
print(f"Multi-SCF ROHF: {rohf_iterations_multi} iterations, E = {energies[1]:.14f} Ha")

# Analysis
extra_uhf_iterations = uhf_iterations_multi - uhf_iterations_independent

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)

print(f"\nUHF extra iterations in multi_scf: {extra_uhf_iterations}")

# Energy check
energy_match_uhf = abs(energies[0] - e_uhf_independent) < 1e-8
energy_match_rohf = abs(energies[1] - e_rohf_independent) < 1e-8

print(f"\nEnergy match UHF: {'✓ PASS' if energy_match_uhf else '✗ FAIL'}")
print(f"Energy match ROHF: {'✓ PASS' if energy_match_rohf else '✗ FAIL'}")

# Iteration check (allow up to 2 extra iterations as reasonable overhead)
iteration_ok = extra_uhf_iterations <= 2

print(f"\nIteration overhead acceptable (≤2): {'✓ PASS' if iteration_ok else '✗ FAIL'}")

if extra_uhf_iterations >= 8:
    print("\n⚠️  WARNING: UHF required +8 or more extra iterations!")
    print("   This suggests the freeze pattern is NOT working correctly.")
elif extra_uhf_iterations <= 2:
    print("\n✓ SUCCESS: Freeze pattern is working!")
    print(f"  UHF only needed +{extra_uhf_iterations} extra iterations (acceptable overhead).")
else:
    print(f"\n⚠️  MARGINAL: UHF needed +{extra_uhf_iterations} extra iterations.")
    print("   Expected ≤2, but less than the bug's +8.")

print("\n" + "="*80)

# Exit with appropriate status
if energy_match_uhf and energy_match_rohf and iteration_ok:
    print("OVERALL: ✓ ALL TESTS PASSED")
    sys.exit(0)
else:
    print("OVERALL: ✗ SOME TESTS FAILED")
    sys.exit(1)
