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

import sys
import os

# Ensure we can import psi4 from installed location, not source tree
# This allows running from any directory
if 'psi4' in sys.modules:
    del sys.modules['psi4']

# Add installed psi4 to path if needed
psi4_install = os.path.join(os.path.dirname(__file__), '..', 'psi4_install', 'lib')
if os.path.exists(psi4_install) and psi4_install not in sys.path:
    sys.path.insert(0, psi4_install)

import psi4
import time
from psi4.driver.procrouting.scf_proc.scf_iterator import multi_scf
from psi4.driver.procrouting.proc import build_functional_and_disp

# Enable debug mode for verbose output
DEBUG = True

def debug_print(msg):
    """Print debug message if DEBUG is enabled"""
    if DEBUG:
        print(f"[DEBUG] {msg}")

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
mol1.update_geometry()

mol2 = psi4.geometry(mol2_str)
mol2.update_geometry()

# Set options
psi4.set_options({
    'basis': 'cc-pvdz',
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
ref_wfn1 = psi4.core.Wavefunction.build(mol1, psi4.core.get_global_option('BASIS'))
ref_wfn2 = psi4.core.Wavefunction.build(mol2, psi4.core.get_global_option('BASIS'))

# Build superfunctional for HF
psi4.core.prepare_options_for_module("SCF")
superfunc, _ = build_functional_and_disp('HF', restricted=True)

wfn1 = psi4.core.RHF(ref_wfn1, superfunc)
wfn2 = psi4.core.RHF(ref_wfn2, superfunc)

# Set DF auxiliary basis (THIS WAS MISSING!)
aux_basis1 = psi4.core.BasisSet.build(wfn1.molecule(), "DF_BASIS_SCF",
                                     psi4.core.get_option("SCF", "DF_BASIS_SCF"),
                                     "JKFIT", psi4.core.get_global_option('BASIS'),
                                     puream=wfn1.basisset().has_puream())
wfn1.set_basisset("DF_BASIS_SCF", aux_basis1)

aux_basis2 = psi4.core.BasisSet.build(wfn2.molecule(), "DF_BASIS_SCF",
                                     psi4.core.get_option("SCF", "DF_BASIS_SCF"),
                                     "JKFIT", psi4.core.get_global_option('BASIS'),
                                     puream=wfn2.basisset().has_puream())
wfn2.set_basisset("DF_BASIS_SCF", aux_basis2)

# Initialize both
print("Initializing wavefunctions...")
wfn1.initialize()
wfn2.initialize()

# DEBUG: Check wavefunction state after initialization
debug_print("=" * 60)
debug_print("WAVEFUNCTION STATE AFTER INITIALIZATION")
debug_print("=" * 60)
debug_print(f"wfn1.Ca() shape: {wfn1.Ca().shape}")
debug_print(f"wfn1.nalpha(): {wfn1.nalpha()}")
debug_print(f"wfn1.nbeta(): {wfn1.nbeta()}")
debug_print(f"wfn1.nso(): {wfn1.nso()}")
debug_print(f"wfn1.nmo(): {wfn1.nmo()}")
debug_print(f"wfn1.Ca_subset('SO', 'OCC') shape: {wfn1.Ca_subset('SO', 'OCC').shape}")

orb_mats1 = wfn1.get_orbital_matrices()
debug_print(f"\nwfn1.get_orbital_matrices() returned {len(orb_mats1)} matrices")
for i, mat in enumerate(orb_mats1):
    debug_print(f"  Matrix {i} shape: {mat.shape}")

orb_mats2 = wfn2.get_orbital_matrices()
debug_print(f"\nwfn2.get_orbital_matrices() returned {len(orb_mats2)} matrices")
for i, mat in enumerate(orb_mats2):
    debug_print(f"  Matrix {i} shape: {mat.shape}")

# DEBUG: Check JK object
jk = wfn1.jk()
debug_print(f"\nJK object: {jk}")
debug_print(f"JK type: {type(jk)}")
debug_print(f"JK C_left() length before clear: {len(jk.C_left())}")
debug_print(f"JK C_right() length before clear: {len(jk.C_right())}")

# DEBUG: Test if C_clear() method exists
debug_print(f"\nChecking available JK methods:")
debug_print(f"  hasattr(jk, 'C_clear'): {hasattr(jk, 'C_clear')}")
debug_print(f"  hasattr(jk, 'C_add'): {hasattr(jk, 'C_add')}")
debug_print(f"  hasattr(jk, 'C_left_add'): {hasattr(jk, 'C_left_add')}")
debug_print(f"  hasattr(jk, 'C_right_add'): {hasattr(jk, 'C_right_add')}")
debug_print("=" * 60)

# Run multi-cycle SCF
print("\nRunning multi_scf() with NEW architecture...\n")
t0 = time.time()

try:
    energies_multi = multi_scf([wfn1, wfn2], e_conv=1e-8, d_conv=1e-6, verbose=True)
except Exception as e:
    print("\n" + "=" * 70)
    print("ERROR IN multi_scf():")
    print("=" * 70)
    print(f"Exception type: {type(e).__name__}")
    print(f"Exception message: {str(e)}")

    # DEBUG: Check JK state when error occurred
    debug_print("\n" + "=" * 60)
    debug_print("JK STATE WHEN ERROR OCCURRED")
    debug_print("=" * 60)
    debug_print(f"JK C_left() length: {len(jk.C_left())}")
    debug_print(f"JK C_right() length: {len(jk.C_right())}")
    debug_print(f"JK J() length: {len(jk.J())}")
    debug_print(f"JK K() length: {len(jk.K())}")

    # Test manual JK operation
    debug_print("\n" + "=" * 60)
    debug_print("TESTING MANUAL JK OPERATIONS")
    debug_print("=" * 60)

    # Test 1: Using .append() (known to fail)
    debug_print("\nTest 1: Using jk.C_left().append() (expected to FAIL)")
    jk.C_left().clear()
    debug_print(f"  After clear: C_left length = {len(jk.C_left())}")
    jk.C_left().append(orb_mats1[0])
    debug_print(f"  After append: C_left length = {len(jk.C_left())} (should be 1, but is 0)")

    # Test 2: Using C_add() (should work)
    if hasattr(jk, 'C_clear') and hasattr(jk, 'C_add'):
        debug_print("\nTest 2: Using jk.C_clear() and jk.C_add() (should WORK)")
        jk.C_clear()
        debug_print(f"  After C_clear(): C_left length = {len(jk.C_left())}")
        jk.C_add(orb_mats1[0])
        debug_print(f"  After C_add(): C_left length = {len(jk.C_left())} (should be 1)")
        jk.C_add(orb_mats2[0])
        debug_print(f"  After C_add(): C_left length = {len(jk.C_left())} (should be 2)")

        if len(jk.C_left()) == 2:
            debug_print("\n  ✓ C_add() works! Matrices added successfully.")
            debug_print("  → scf_iterator.py needs to use jk.C_clear() and jk.C_add()")
        else:
            debug_print("\n  ✗ C_add() also failed!")

    debug_print("=" * 60)
    print("\nTest aborted due to error.")
    sys.exit(1)

t_multi = time.time() - t0

print(f"\nMulti-cycle SCF total time: {t_multi:.3f} seconds\n")

# DEBUG: Success case - show what worked
debug_print("=" * 60)
debug_print("SUCCESS - multi_scf() completed")
debug_print("=" * 60)
debug_print(f"Returned energies: {energies_multi}")
debug_print(f"Energy type: {type(energies_multi)}")
debug_print(f"Number of energies: {len(energies_multi)}")
debug_print("=" * 60)

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
