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
import numpy as np
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

# Create COMPLETELY FRESH molecules (don't reuse mol1/mol2 from baseline test)
mol1_multi = psi4.geometry(mol1_str)
mol1_multi.update_geometry()

mol2_multi = psi4.geometry(mol2_str)
mol2_multi.update_geometry()

# Create fresh wavefunctions
ref_wfn1 = psi4.core.Wavefunction.build(mol1_multi, psi4.core.get_global_option('BASIS'))
ref_wfn2 = psi4.core.Wavefunction.build(mol2_multi, psi4.core.get_global_option('BASIS'))

# Build SEPARATE superfunctional for each wavefunction (don't share!)
psi4.core.prepare_options_for_module("SCF")
superfunc1, _ = build_functional_and_disp('HF', restricted=True)
superfunc2, _ = build_functional_and_disp('HF', restricted=True)

wfn1 = psi4.core.RHF(ref_wfn1, superfunc1)
wfn2 = psi4.core.RHF(ref_wfn2, superfunc2)

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

# CRITICAL FIX: Ensure DIIS settings are IDENTICAL for both wavefunctions
# The problem is that wfn1.diis_start_ and wfn2.diis_start_ may differ
# because they read from global options at different times (possibly modified by baseline test)
# We must manually synchronize them to avoid non-deterministic behavior
print(f"Synchronizing DIIS settings: wfn1.diis_start_={wfn1.diis_start_}, wfn2.diis_start_={wfn2.diis_start_}")
if wfn1.diis_start_ != wfn2.diis_start_:
    print(f"WARNING: DIIS_START mismatch! Setting both to {wfn1.diis_start_}")
    wfn2.diis_start_ = wfn1.diis_start_
if wfn1.diis_enabled_ != wfn2.diis_enabled_:
    print(f"WARNING: DIIS_ENABLED mismatch! Setting both to {wfn1.diis_enabled_}")
    wfn2.diis_enabled_ = wfn1.diis_enabled_

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

# CRITICAL: Check if DIIS is enabled for both (this could cause non-determinism!)
debug_print(f"\nCHECKING DIIS SETTINGS:")
debug_print(f"  wfn1.diis_enabled_: {wfn1.diis_enabled_}")
debug_print(f"  wfn2.diis_enabled_: {wfn2.diis_enabled_}")
debug_print(f"  wfn1.diis_start_: {wfn1.diis_start_}")
debug_print(f"  wfn2.diis_start_: {wfn2.diis_start_}")
debug_print(f"  ARE DIIS SETTINGS SAME? {wfn1.diis_enabled_ == wfn2.diis_enabled_}")
if wfn1.diis_enabled_ != wfn2.diis_enabled_:
    debug_print(f"  *** WARNING: WFN1 and WFN2 have DIFFERENT DIIS settings! ***")
    debug_print(f"  *** This will cause non-deterministic behavior! ***")

orb_mats1 = wfn1.get_orbital_matrices()
debug_print(f"\nwfn1.get_orbital_matrices() returned {len(orb_mats1)} matrices")
for i, mat in enumerate(orb_mats1):
    debug_print(f"  Matrix {i} shape: {mat.shape}")

orb_mats2 = wfn2.get_orbital_matrices()
debug_print(f"\nwfn2.get_orbital_matrices() returned {len(orb_mats2)} matrices")
for i, mat in enumerate(orb_mats2):
    debug_print(f"  Matrix {i} shape: {mat.shape}")

# DEBUG: Check JK object
jk1 = wfn1.jk()
jk2 = wfn2.jk()
debug_print(f"\nJK object 1: {jk1} (id: {id(jk1)})")
debug_print(f"JK object 2: {jk2} (id: {id(jk2)})")
debug_print(f"JK objects are same? {jk1 is jk2}")
debug_print(f"JK type: {type(jk1)}")

# CRITICAL FIX: Must use wfn1's JK, NOT wfn2's, because only wfn1's JK has been properly initialized
# with auxiliary basis sets. Using wfn2's JK would cause inconsistent behavior.
jk = jk1
debug_print(f"\nUsing JK from wfn1 for shared computation (wfn2's JK not initialized)")
debug_print(f"JK C_left() length before clear: {len(jk.C_left())}")
debug_print(f"JK C_right() length before clear: {len(jk.C_right())}")

# DEBUG: Test if C_clear() method exists
debug_print(f"\nChecking available JK methods:")
debug_print(f"  hasattr(jk, 'C_clear'): {hasattr(jk, 'C_clear')}")
debug_print(f"  hasattr(jk, 'C_add'): {hasattr(jk, 'C_add')}")
debug_print(f"  hasattr(jk, 'C_left_add'): {hasattr(jk, 'C_left_add')}")
debug_print(f"  hasattr(jk, 'C_right_add'): {hasattr(jk, 'C_right_add')}")
debug_print("=" * 60)

# Run multi-cycle SCF with iteration monitoring
print("\nRunning multi-cycle SCF with shared JK computation...\n")

# Manual multi-SCF loop with detailed monitoring
t0 = time.time()
max_iter = 20  # Run 20 iterations to see convergence pattern
e_conv = 1e-8
d_conv = 1e-6

print(f"\nConvergence thresholds: E={e_conv:.2e}, D={d_conv:.2e}")
print(f"Running {max_iter} iterations...\n")
print("Note: Both WFN1 and WFN2 start with SAME initial guess (SAD),")
print("      so they should converge to SAME result (identical molecules).\n")

# Initialize iteration state
wfn1._scf_initialize_iteration_state(e_conv, d_conv)
wfn2._scf_initialize_iteration_state(e_conv, d_conv)

# CRITICAL: Reset iteration counters to ensure both start from same state
# iteration_ may have been incremented by initialize() or baseline test
print(f"Resetting iteration counters: wfn1.iteration_={wfn1.iteration_}, wfn2.iteration_={wfn2.iteration_}")
wfn1.iteration_ = 0
wfn2.iteration_ = 0

converged_1 = False
converged_2 = False
E_prev_1 = 0.0
E_prev_2 = 0.0

print(f"{'Iter':>4} | {'E1 (Ha)':>16} | {'ΔE1':>10} | {'ΔD1':>10} | {'E2 (Ha)':>16} | {'ΔE2':>10} | {'ΔD2':>10} | {'Status':>15}")
print("-" * 120)

try:
    for iteration in range(1, max_iter + 1):
        # Collect C matrices from both wavefunctions
        C_occ_list = []
        wfn_state_counts = []

        for wfn in [wfn1, wfn2]:
            C_matrices = wfn.get_orbital_matrices()
            C_occ_list.extend(C_matrices)
            wfn_state_counts.append(len(C_matrices))

        # DEBUG: Check C matrix identity for iteration 2
        if iteration == 2 and DEBUG:
            debug_print(f"\n  === C MATRICES BEFORE JK ===")
            debug_print(f"  C_occ_list[0][0,0] = {C_occ_list[0].get(0,0):.12f} (from wfn1)")
            debug_print(f"  C_occ_list[1][0,0] = {C_occ_list[1].get(0,0):.12f} (from wfn2)")
            debug_print(f"  id(C_occ_list[0]) = {id(C_occ_list[0])}")
            debug_print(f"  id(C_occ_list[1]) = {id(C_occ_list[1])}")

        # Shared JK computation
        jk.C_clear()
        for idx, C_occ in enumerate(C_occ_list):
            jk.C_add(C_occ)
            if iteration == 2 and DEBUG:
                debug_print(f"  Added C[{idx}] to JK: C[0,0]={C_occ.get(0,0):.12f}")
        jk.compute()

        # Get J/K matrices
        J_all = jk.J()
        K_all = jk.K()

        # DEBUG: Check JK results for first iteration
        if iteration == 1 and DEBUG:
            debug_print(f"\n  JK computation results:")
            debug_print(f"    len(C_occ_list) = {len(C_occ_list)}")
            debug_print(f"    len(J_all) = {len(J_all)}")
            debug_print(f"    len(K_all) = {len(K_all)}")
            debug_print(f"    J_all[0] id: {id(J_all[0])}, shape: {J_all[0].shape}")
            debug_print(f"    J_all[1] id: {id(J_all[1])}, shape: {J_all[1].shape}")
            debug_print(f"    J_all[0] is J_all[1]: {J_all[0] is J_all[1]}")

            # Check if they point to the same underlying data
            J0_arr = J_all[0].to_array(dense=True)
            J1_arr = J_all[1].to_array(dense=True)
            debug_print(f"    J0[0,0] = {J0_arr[0,0]:.12f}")
            debug_print(f"    J1[0,0] = {J1_arr[0,0]:.12f}")
            debug_print(f"    J0 and J1 share memory? {np.shares_memory(J0_arr, J1_arr)}")

        # DEBUG: Check J/K results AFTER JK computation
        if iteration == 2 and DEBUG:
            debug_print(f"\n  === J/K MATRICES AFTER JK.compute() ===")
            debug_print(f"  len(J_all) = {len(J_all)}, len(K_all) = {len(K_all)}")
            for idx in range(len(J_all)):
                J_trace = sum([J_all[idx].get(i, i) for i in range(min(5, J_all[idx].shape[0]))])
                K_trace = sum([K_all[idx].get(i, i) for i in range(min(5, K_all[idx].shape[0]))])
                debug_print(f"  J[{idx}] trace(first 5) = {J_trace:.6f}, K[{idx}] trace(first 5) = {K_trace:.6f}")

        # Distribute J/K back to wavefunctions
        jk_index = 0
        for wfn_idx, (wfn, n_states) in enumerate(zip([wfn1, wfn2], wfn_state_counts)):
            J_list = [J_all[jk_index + i] for i in range(n_states)]
            K_list = [K_all[jk_index + i] for i in range(n_states)]

            if iteration == 2 and DEBUG:
                debug_print(f"\n  === DISTRIBUTING J/K TO WFN{wfn_idx+1} ===")
                debug_print(f"  jk_index={jk_index}, n_states={n_states}")
                debug_print(f"  Assigning J_all[{jk_index}], K_all[{jk_index}] to wfn{wfn_idx+1}")

            wfn.set_jk_matrices(J_list, K_list)
            jk_index += n_states

        # DEBUG: Before iteration - check if J/K were set correctly
        if iteration == 2 and DEBUG:
            debug_print(f"\n  === ITERATION 2 DEBUG (BEFORE _scf_iteration) ===")
            # Check if wfn1 has received J/K
            debug_print(f"  wfn1 use_precomputed_jk_: {wfn1.use_precomputed_jk_ if hasattr(wfn1, 'use_precomputed_jk_') else 'N/A'}")
            debug_print(f"  wfn2 use_precomputed_jk_: {wfn2.use_precomputed_jk_ if hasattr(wfn2, 'use_precomputed_jk_') else 'N/A'}")

        # Each wavefunction completes its iteration
        # WFN1
        continue_flag_1, reason_1 = wfn1._scf_iteration()
        E1 = wfn1.compute_E()
        converged_1 = (reason_1 == 'converged')
        dE1 = abs(E1 - E_prev_1) if iteration > 1 else 0.0
        E_prev_1 = E1

        # WFN2
        continue_flag_2, reason_2 = wfn2._scf_iteration()
        E2 = wfn2.compute_E()
        converged_2 = (reason_2 == 'converged')
        dE2 = abs(E2 - E_prev_2) if iteration > 1 else 0.0
        E_prev_2 = E2

        # DEBUG: After iteration - check Ca matrices
        if iteration == 2 and DEBUG:
            Ca1_after = wfn1.Ca()
            Ca2_after = wfn2.Ca()
            debug_print(f"\n  === ITERATION 2 DEBUG (AFTER _scf_iteration) ===")
            debug_print(f"  Ca1[0,0] = {Ca1_after.get(0,0):.12f}")
            debug_print(f"  Ca2[0,0] = {Ca2_after.get(0,0):.12f}")
            debug_print(f"  E1 = {E1:.12f}")
            debug_print(f"  E2 = {E2:.12f}")
            debug_print(f"  Are they locked? {abs(E1 - E2) < 1e-6}")

        # Get convergence information
        dD1 = wfn1._scf_Dnorm if hasattr(wfn1, '_scf_Dnorm') else 0.0
        dD2 = wfn2._scf_Dnorm if hasattr(wfn2, '_scf_Dnorm') else 0.0

        # Status
        status1 = "CONV" if converged_1 else "RUN"
        status2 = "CONV" if converged_2 else "RUN"
        status = f"{status1}/{status2}"

        # Get matrix norms for diagnostics
        C_occ1 = wfn1.get_orbital_matrices()[0]
        D1 = wfn1.Da()
        F1 = wfn1.Fa()

        C_norm1 = np.linalg.norm(C_occ1.to_array())
        D_norm1 = D1.vector_dot(D1)
        F_norm1 = F1.vector_dot(F1)
        J_norm1 = J_all[0].vector_dot(J_all[0])
        K_norm1 = K_all[0].vector_dot(K_all[0])

        C_occ2 = wfn2.get_orbital_matrices()[0]
        D2 = wfn2.Da()
        F2 = wfn2.Fa()

        C_norm2 = np.linalg.norm(C_occ2.to_array())
        D_norm2 = D2.vector_dot(D2)
        F_norm2 = F2.vector_dot(F2)
        J_norm2 = J_all[1].vector_dot(J_all[1])
        K_norm2 = K_all[1].vector_dot(K_all[1])

        # DEBUG: Check if C matrices and J/K matrices are identical (iteration 1-3 only)
        if iteration <= 3 and DEBUG:
            C1_arr = C_occ1.to_array()
            C2_arr = C_occ2.to_array()
            C_diff = np.linalg.norm(C1_arr - C2_arr)

            J1_arr = J_all[0].to_array(dense=True)
            J2_arr = J_all[1].to_array(dense=True)
            J_diff = np.linalg.norm(J1_arr - J2_arr)

            K1_arr = K_all[0].to_array(dense=True)
            K2_arr = K_all[1].to_array(dense=True)
            K_diff = np.linalg.norm(K1_arr - K2_arr)

            debug_print(f"  Iteration {iteration}: ||C1 - C2|| = {C_diff:.6e}, ||J1 - J2|| = {J_diff:.6e}, ||K1 - K2|| = {K_diff:.6e}")

        # Print iteration summary
        print(f"{iteration:4d} | {E1:16.10f} | {dE1:10.2e} | {dD1:10.2e} | {E2:16.10f} | {dE2:10.2e} | {dD2:10.2e} | {status:>15}")

        # WFN1 details with DIIS info
        diis_str1 = f"err={dD1:.2e}"
        print(f"       WFN1: C={C_norm1:.4f} D={D_norm1:.4f} J={J_norm1:.2f} K={K_norm1:.2f} F={F_norm1:.2f} | Conv: ΔE={dE1:.2e} | DIIS: {diis_str1}")

        # WFN2 details with DIIS info
        diis_str2 = f"err={dD2:.2e}"
        print(f"       WFN2: C={C_norm2:.4f} D={D_norm2:.4f} J={J_norm2:.2f} K={K_norm2:.2f} F={F_norm2:.2f} | Conv: ΔE={dE2:.2e} | DIIS: {diis_str2}")

        # Check convergence
        if converged_1 and converged_2:
            print()
            print(f"✓ Both converged in {iteration} iterations")
            energies_multi = [E1, E2]
            break
    else:
        # Did not converge in max_iter iterations
        print()
        print(f"⚠ Did not fully converge: WFN1={'CONV' if converged_1 else 'NO'}, WFN2={'CONV' if converged_2 else 'NO'}")
        energies_multi = [E1, E2]
except Exception as e:
    print("\n" + "=" * 70)
    print("ERROR IN multi_scf():")
    print("=" * 70)
    print(f"Exception type: {type(e).__name__}")
    print(f"Exception message: {str(e)}")
    print()
    print("Test aborted due to error.")
    sys.exit(1)

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
pass_overall = pass1 and pass2 and converged_1 and converged_2

print("\n" + "="*70)
print("  Summary")
print("="*70)
print(f"  WFN1: {'CONV ✓' if converged_1 else 'NO'}, E_diff={diff1:.2e}")
print(f"  WFN2: {'CONV ✓' if converged_2 else 'NO'}, E_diff={diff2:.2e}")
print(f"  Speedup: {speedup:.2f}x")
if pass_overall:
    print(f"  Status: PASS ✓")
else:
    print(f"  Status: FAIL (energies don't match or not converged)")
print("="*70 + "\n")
