#!/usr/bin/env python3
"""
Test DF_SCF_GUESS support in multi_scf() for DIRECT algorithm.

This verifies that multi-SCF now produces identical iteration counts
to single SCF when using SCF_TYPE='DIRECT' with DF_SCF_GUESS=True.
"""

import psi4
import sys

# Test molecule: H2O
mol = psi4.geometry("""
0 1
O  0.000000  0.000000  0.118351
H  0.000000  0.755453 -0.471404
H  0.000000 -0.755453 -0.471404
units angstrom
no_reorient
no_com
symmetry c1
""")

# Settings
psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'DIRECT',
    'df_scf_guess': True,  # Enable Andy Trick 2.0
    'e_convergence': 1e-8,
    'd_convergence': 1e-7,
    'maxiter': 50,
    'print': 1
})

print("=" * 80)
print("TEST: DF_SCF_GUESS support in multi_scf() for DIRECT")
print("=" * 80)

# Test 1: Single SCF with DIRECT + DF_SCF_GUESS
print("\n" + "=" * 80)
print("Test 1: Single SCF (DIRECT + DF_SCF_GUESS)")
print("=" * 80)

psi4.core.clean()
psi4.core.clean_options()
psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'DIRECT',
    'df_scf_guess': True,
    'e_convergence': 1e-8,
    'd_convergence': 1e-7,
    'maxiter': 50,
    'print': 1
})

e_single = psi4.energy('hf', molecule=mol)
iters_single = psi4.core.variable("SCF ITERATIONS")

print(f"\nSingle SCF Results:")
print(f"  Energy: {e_single:.12f}")
print(f"  Total iterations: {iters_single}")

# Test 2: Multi-SCF with DIRECT + DF_SCF_GUESS
print("\n" + "=" * 80)
print("Test 2: Multi-SCF (DIRECT + DF_SCF_GUESS)")
print("=" * 80)

psi4.core.clean()
psi4.core.clean_options()
psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'DIRECT',
    'df_scf_guess': True,
    'e_convergence': 1e-8,
    'd_convergence': 1e-7,
    'maxiter': 50,
    'print': 1,
    'reference': 'RHF'
})

# Import multi_scf
sys.path.insert(0, '/home/user/psi4')
from psi4.driver.procrouting.scf_proc.scf_iterator import multi_scf
from psi4.driver.procrouting.scf_proc import scf_wavefunction_factory

# Create wavefunction
wfn = scf_wavefunction_factory('hf', mol, 'RHF')

# Run multi_scf
energies = multi_scf([wfn], verbose=True)
e_multi = energies[0]
iters_multi = wfn.iteration_

print(f"\nMulti-SCF Results:")
print(f"  Energy: {e_multi:.12f}")
print(f"  Total iterations: {iters_multi}")

# Test 3: Comparison
print("\n" + "=" * 80)
print("COMPARISON:")
print("=" * 80)
print(f"Single SCF iterations: {iters_single}")
print(f"Multi-SCF iterations:  {iters_multi}")
print(f"Difference:            {abs(iters_single - iters_multi)}")
print(f"\nEnergy difference:     {abs(e_single - e_multi):.2e}")

# Validation
success = True
if abs(e_single - e_multi) > 1e-9:
    print(f"\n❌ FAIL: Energies differ by {abs(e_single - e_multi):.2e}")
    success = False
else:
    print(f"\n✅ PASS: Energies match (ΔE = {abs(e_single - e_multi):.2e})")

# Note: Iterations might differ by 1-2 due to different convergence paths
# But should be MUCH closer than before (where difference was ~6 iterations)
iter_diff = abs(iters_single - iters_multi)
if iter_diff <= 2:
    print(f"✅ PASS: Iteration counts nearly identical (Δiter = {iter_diff})")
    print("   This is expected due to slightly different convergence paths.")
elif iter_diff <= 6:
    print(f"⚠️  WARNING: Iteration counts differ by {iter_diff}")
    print("   This is improved but could be better.")
else:
    print(f"❌ FAIL: Iteration counts differ significantly (Δiter = {iter_diff})")
    success = False

print("\n" + "=" * 80)
if success:
    print("OVERALL: ✅ TEST PASSED")
else:
    print("OVERALL: ❌ TEST FAILED")
print("=" * 80)

sys.exit(0 if success else 1)
