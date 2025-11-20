#!/usr/bin/env python3
"""
Universal SCF Tester for Psi4
Tests different reference types (RHF, UHF, ROHF) with various methods and molecules
Logs each calculation to separate file and saves results to SQLite database
"""

import sys
sys.path.insert(0, '/Users/stan/Documents/Workspaces/nanoreactor/codebase/psi4_dir/psi4_install/lib/')
import psi4
from datetime import datetime
import sqlite3
import os
import json
import time

# Configuration
psi4.set_memory('2 GB')
psi4.set_num_threads(8)

# Paths for logs and database
LOG_DIR = 'test_logs'
DB_FILE = 'scf_test_results.db'
ENERGY_TOLERANCE = 1.0e-6  # Hartree

# Create log directory if it doesn't exist
os.makedirs(LOG_DIR, exist_ok=True)

# Define test molecules
MOLECULES = {
    'H2O': {
        'geometry': """
            O       0.00000000     0.00000000     0.00000000
            H       1.00000000     0.50000000     0.00000000
            H      -1.00000000     0.50000000     0.00000000
        """,
        'charge': 0,
        'multiplicity': 1,
        'type': 'closed-shell'
    },
    'OH': {
        'geometry': """
            O       0.00000000     0.00000000     0.00000000
            H       0.00000000     0.00000000     0.96000000
        """,
        'charge': 0,
        'multiplicity': 2,
        'type': 'open-shell'
    },
    'CH3': {
        'geometry': """
            C       0.00000000     0.00000000     0.00000000
            H       1.07000000     0.00000000     0.00000000
            H      -0.53500000     0.92661000     0.00000000
            H      -0.53500000    -0.92661000     0.00000000
        """,
        'charge': 0,
        'multiplicity': 2,
        'type': 'open-shell'
    }
}

# Define test configurations
TEST_CONFIGS = [
    # Basic reference tests
    {'reference': 'RHF', 'methods': ['hf', 'pbe', 'b3lyp'], 'molecules': ['H2O'], 'options': {}},
    {'reference': 'UHF', 'methods': ['hf', 'pbe'], 'molecules': ['H2O', 'OH', 'CH3'], 'options': {}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH', 'CH3'], 'options': {}},

    # === RHF: SCF_TYPE tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'PK'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'MEM_DF'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'DISK_DF'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'OUT_OF_CORE'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'CD'}},

    # === UHF: SCF_TYPE tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'PK'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'DF'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'MEM_DF'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'DISK_DF'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'OUT_OF_CORE'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'CD'}},

    # === ROHF: SCF_TYPE tests ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'PK'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'DF'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'MEM_DF'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'DISK_DF'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'OUT_OF_CORE'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'CD'}},

    # === RHF: SCREENING tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'CSAM'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'DENSITY'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'NONE'}},

    # === UHF: SCREENING tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'CSAM'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'DENSITY'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'NONE'}},

    # === ROHF: SCREENING tests ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'CSAM'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'DENSITY'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'screening': 'NONE'}},

    # === RHF: GUESS tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'CORE'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'SAD'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'AUTO'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'MODHUCKEL'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'GWH'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'SADNO'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'SAP'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'SAPGAU'}},

    # === UHF: GUESS tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'CORE'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAD'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'AUTO'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'MODHUCKEL'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'GWH'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SADNO'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAP'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAPGAU'}},

    # === ROHF: GUESS tests ===
    # Note: SAD/AUTO require df_scf_guess=true for ROHF open-shell systems to enable DF fallback
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'CORE'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAD', 'df_scf_guess': True}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'AUTO', 'df_scf_guess': True}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'MODHUCKEL'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'GWH'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SADNO'}},  # Natural orbitals
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAP'}},
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'guess': 'SAPGAU'}},

    # === RHF: INCFOCK tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'incfock': 'true'}},

    # === UHF: INCFOCK tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'incfock': 'true'}},

    # === ROHF: INCFOCK tests ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'incfock': 'true'}},

    # === RHF: DIIS tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'diis': 'false'}},

    # === UHF: DIIS tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'diis': 'false'}},

    # === ROHF: DIIS tests ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'diis': 'false'}},

    # === RHF: SCF_INITIAL_ACCELERATOR tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_initial_accelerator': 'EDIIS'}},
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_initial_accelerator': 'NONE'}},

    # === UHF: SCF_INITIAL_ACCELERATOR tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_initial_accelerator': 'EDIIS'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_initial_accelerator': 'NONE'}},

    # === ROHF: SCF_INITIAL_ACCELERATOR tests (note: may not work for ROHF) ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_initial_accelerator': 'NONE'}},

    # === RHF: DAMPING tests ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'damping_percentage': 20.0}},

    # === UHF: DAMPING tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'damping_percentage': 20.0}},

    # === ROHF: DAMPING tests ===
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'damping_percentage': 20.0}},

    # === UHF-specific: MOM tests (Maximum Overlap Method for excited states) ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'mom_start': 2, 'mom_occ': [5], 'mom_vir': [6]}},

    # === UHF-specific: Fractional occupation tests ===
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'frac_start': 2, 'frac_occ': [5], 'frac_val': [0.5], 'mom_start': 0}},

    # ========================================================================
    # MULTI-SCF VALIDATION EQUIVALENTS (for comparison with test_ultimate_validation.py)
    # These tests provide reference data for multi-SCF validation
    # Each config below corresponds to configs in test_ultimate_validation.py POSITIVE_TESTS
    # ========================================================================

    # === MULTI-SCF EQUIV: UHF variants for H2O (Basic test coverage) ===
    # Corresponds to: "Basic: RHF HF+PBE" - UHF variant
    {'reference': 'UHF', 'methods': ['hf', 'pbe'], 'molecules': ['H2O'], 'options': {}},

    # === MULTI-SCF EQUIV: UHF with SCF_TYPE variants ===
    # Corresponds to: "SCF_TYPE: PK", "SCF_TYPE: DF"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'PK'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF'}},

    # === MULTI-SCF EQUIV: UHF with SCREENING variants ===
    # Corresponds to: "SCREENING: CSAM", "SCREENING: DENSITY"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'CSAM'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'DENSITY'}},

    # === MULTI-SCF EQUIV: UHF with GUESS variants ===
    # Corresponds to: "GUESS: SAD", "GUESS: CORE"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'SAD'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'guess': 'CORE'}},

    # === MULTI-SCF EQUIV: UHF with DIIS disabled ===
    # Corresponds to: "DIIS: false"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'diis': 'false'}},

    # === MULTI-SCF EQUIV: UHF with DAMPING ===
    # Corresponds to: "DAMPING: 20%"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'damping_percentage': 20.0}},

    # === MULTI-SCF EQUIV: UHF with INCFOCK ===
    # Corresponds to: "INCFOCK: true"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'incfock': 'true'}},

    # === MULTI-SCF EQUIV: UHF with SCF_INITIAL_ACCELERATOR ===
    # Corresponds to: "ACCELERATOR: EDIIS"
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_initial_accelerator': 'EDIIS'}},

    # === MULTI-SCF EQUIV: Combination tests ===
    # Corresponds to: "COMBO: DF+SAD+EDIIS" - RHF variant
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF', 'guess': 'SAD', 'scf_initial_accelerator': 'EDIIS'}},
    # Corresponds to: "COMBO: DF+SAD+EDIIS" - UHF variant
    {'reference': 'UHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF', 'guess': 'SAD', 'scf_initial_accelerator': 'EDIIS'}},

    # ========================================================================
    # GRADIENT-SPECIFIC TESTS (for multi-grad() validation)
    # These tests exercise gradient calculation paths with different options
    # ========================================================================

    # === DFT Grid density (affects DFT gradients) ===
    # Coarse grid
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'dft_spherical_points': 110, 'dft_radial_points': 20}},
    # Medium grid (default-like)
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'dft_spherical_points': 302, 'dft_radial_points': 75}},
    # Fine grid
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'dft_spherical_points': 590, 'dft_radial_points': 99}},

    # === Different DFT functional types (different gradient implementations) ===
    # GGA functionals
    {'reference': 'RHF', 'methods': ['pbe', 'bp86', 'blyp'], 'molecules': ['H2O'], 'options': {}},
    # Hybrid functionals
    {'reference': 'RHF', 'methods': ['b3lyp', 'pbe0'], 'molecules': ['H2O'], 'options': {}},
    # Meta-GGA (if available)
    {'reference': 'RHF', 'methods': ['tpss'], 'molecules': ['H2O'], 'options': {}},
    # Range-separated hybrids
    {'reference': 'RHF', 'methods': ['wb97x'], 'molecules': ['H2O'], 'options': {}},

    # UHF with different functionals (open-shell gradients)
    {'reference': 'UHF', 'methods': ['pbe', 'b3lyp', 'pbe0'], 'molecules': ['OH'], 'options': {}},
    {'reference': 'UHF', 'methods': ['tpss'], 'molecules': ['OH'], 'options': {}},

    # === SCF_TYPE + Functional combinations (gradient code paths) ===
    # DF with various functionals
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF'}},
    {'reference': 'RHF', 'methods': ['b3lyp'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF'}},
    # CD with functionals
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'scf_type': 'CD'}},
    # PK with functionals
    {'reference': 'UHF', 'methods': ['pbe'], 'molecules': ['OH'], 'options': {'scf_type': 'PK'}},

    # === Gradient with different basis sets (basis set dependencies) ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'basis': 'cc-pvdz'}},
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'basis': 'cc-pvdz'}},
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'basis': 'cc-pvdz'}},

    # === DF gradient with different auxiliary basis ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF', 'df_basis_scf': 'cc-pvdz-jkfit'}},
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'scf_type': 'DF', 'df_basis_scf': 'cc-pvdz-jkfit'}},

    # === Gradient + convergence options ===
    # Tighter convergence (affects gradient accuracy)
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'e_convergence': 1.0e-10, 'd_convergence': 1.0e-10}},
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'e_convergence': 1.0e-10, 'd_convergence': 1.0e-10}},

    # === Gradient + screening (affects JK build, thus gradient) ===
    {'reference': 'RHF', 'methods': ['hf'], 'molecules': ['H2O'], 'options': {'screening': 'CSAM'}},
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'], 'options': {'screening': 'CSAM'}},
    {'reference': 'UHF', 'methods': ['pbe'], 'molecules': ['OH'], 'options': {'screening': 'DENSITY'}},

    # === Open-shell gradient tests (ROHF and UHF) ===
    # ROHF gradients with different functionals (if available)
    {'reference': 'ROHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'scf_type': 'DF'}},
    # UHF gradients with tight convergence
    {'reference': 'UHF', 'methods': ['hf'], 'molecules': ['OH'], 'options': {'e_convergence': 1.0e-10}},
    {'reference': 'UHF', 'methods': ['pbe'], 'molecules': ['CH3'], 'options': {'scf_type': 'DF'}},

    # === Combination: DF + DFT + Grid ===
    {'reference': 'RHF', 'methods': ['pbe'], 'molecules': ['H2O'],
     'options': {'scf_type': 'DF', 'dft_spherical_points': 302, 'dft_radial_points': 75}},
    {'reference': 'UHF', 'methods': ['b3lyp'], 'molecules': ['OH'],
     'options': {'scf_type': 'DF', 'dft_spherical_points': 302, 'dft_radial_points': 75}},
]

# Common SCF options
BASE_OPTIONS = {
    'debug': 3,
    'basis': '6-31g(d)',
    'puream': 'false',
    'guess': 'huckel',
    'e_convergence': 1.0e-8,
    'd_convergence': 1.0e-8,
    'scf_type': 'DIRECT',
    'screening': 'schwarz',
    'incfock': 'false',
    'diis': 'true',
    'damping_percentage': 0.0,
    'dft_basis_tolerance': 3.59e-10,
    'df_scf_guess': 'false',
    'pk_max_buckets': 1200,
    'maxiter': 100,
    'dft_spherical_points': 302,
    'dft_radial_points': 96
}


class TestDatabase:
    """Database manager for test results and reference values"""

    def __init__(self, db_file):
        self.db_file = db_file
        self.conn = None
        self.init_database()

    def init_database(self):
        """Initialize database with required tables"""
        self.conn = sqlite3.connect(self.db_file)
        cursor = self.conn.cursor()

        # Table for test results
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS test_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp TEXT NOT NULL,
                molecule TEXT NOT NULL,
                reference TEXT NOT NULL,
                method TEXT NOT NULL,
                basis TEXT NOT NULL,
                energy REAL,
                iterations INTEGER,
                computation_time REAL,
                status TEXT NOT NULL,
                error_message TEXT,
                log_file TEXT,
                deviation_from_ref REAL
            )
        ''')

        # Table for reference values (ground truth)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS reference_values (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule TEXT NOT NULL,
                reference TEXT NOT NULL,
                method TEXT NOT NULL,
                basis TEXT NOT NULL,
                energy REAL NOT NULL,
                iterations INTEGER,
                computation_time REAL,
                notes TEXT,
                created_at TEXT NOT NULL,
                UNIQUE(molecule, reference, method, basis)
            )
        ''')

        self.conn.commit()

    def save_result(self, result_data):
        """Save test result to database"""
        cursor = self.conn.cursor()
        cursor.execute('''
            INSERT INTO test_results
            (timestamp, molecule, reference, method, basis, energy, iterations,
             computation_time, status, error_message, log_file,
             deviation_from_ref)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            result_data['timestamp'],
            result_data['molecule'],
            result_data['reference'],
            result_data['method'],
            result_data['basis'],
            result_data.get('energy'),
            result_data.get('iterations'),
            result_data.get('computation_time'),
            result_data['status'],
            result_data.get('error_message'),
            result_data.get('log_file'),
            result_data.get('deviation_from_ref')
        ))
        self.conn.commit()

    def get_reference(self, molecule, reference, method, basis):
        """Get reference value for comparison"""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT energy, iterations, computation_time FROM reference_values
            WHERE molecule=? AND reference=? AND method=? AND basis=?
        ''', (molecule, reference, method, basis))
        result = cursor.fetchone()
        return result if result else None

    def set_reference(self, molecule, reference, method, basis, energy, iterations, computation_time=None, notes=''):
        """Set a reference value (ground truth)"""
        cursor = self.conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO reference_values
            (molecule, reference, method, basis, energy, iterations, computation_time, notes, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (molecule, reference, method, basis, energy, iterations, computation_time, notes,
              datetime.now().isoformat()))
        self.conn.commit()
        time_str = f", {computation_time*1000:.1f}ms" if computation_time else ""
        print(f"Reference set: {molecule}/{reference}/{method}/{basis} = {energy:.15f} Ha{time_str}")

    def get_test_history(self, molecule=None, reference=None, method=None, limit=10):
        """Get recent test history"""
        cursor = self.conn.cursor()
        query = "SELECT * FROM test_results WHERE 1=1"
        params = []

        if molecule:
            query += " AND molecule=?"
            params.append(molecule)
        if reference:
            query += " AND reference=?"
            params.append(reference)
        if method:
            query += " AND method=?"
            params.append(method)

        query += " ORDER BY timestamp DESC LIMIT ?"
        params.append(limit)

        cursor.execute(query, params)
        return cursor.fetchall()

    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


def run_test(molecule_name, reference, method, db, extra_options=None):
    """Run a single SCF test with logging and database storage"""
    mol_data = MOLECULES[molecule_name]
    basis = BASE_OPTIONS['basis']

    if extra_options is None:
        extra_options = {}

    # Generate log file name with timestamp (include option signature if present)
    timestamp = datetime.now()
    option_suffix = ""
    if extra_options:
        # Format option values properly for filenames (no brackets, spaces, etc)
        formatted_opts = []
        for k, v in sorted(extra_options.items()):
            if isinstance(v, list):
                v_str = "_".join(str(x).replace(".", "p") for x in v)
            else:
                v_str = str(v).replace(" ", "_").replace(".", "p")
            formatted_opts.append(f"{k}_{v_str}")
        option_suffix = "_" + "_".join(formatted_opts)
    log_filename = f"{molecule_name}_{reference}_{method.upper()}{option_suffix}_{timestamp.strftime('%Y%m%d_%H%M%S')}.log"
    log_path = os.path.join(LOG_DIR, log_filename)

    # Prepare result data
    result_data = {
        'timestamp': timestamp.isoformat(),
        'molecule': molecule_name,
        'reference': reference,
        'method': method.upper(),
        'basis': basis,
        'log_file': log_path,
        'status': 'FAILED',
        'error_message': None,
        'deviation_from_ref': None
    }

    # Set up geometry with charge and multiplicity
    geom_string = f"{mol_data['charge']} {mol_data['multiplicity']}\n{mol_data['geometry']}"

    try:
        # Clean previous calculation and options
        psi4.core.clean()
        psi4.core.clean_options()  # Clear options from previous test to ensure isolation

        # Set psi4 output to log file
        psi4.core.set_output_file(log_path, False)

        # Create molecule
        mol = psi4.geometry(geom_string)

        # Set options (merge base options with extra options)
        options = BASE_OPTIONS.copy()
        options['reference'] = reference
        options.update(extra_options)
        psi4.set_options(options)

        # Run calculation and time it
        start_time = time.perf_counter()
        energy = psi4.energy(method)
        computation_time = time.perf_counter() - start_time

        # Get iterations
        iterations = int(psi4.core.variable('SCF ITERATIONS'))

        # Update result data
        result_data['status'] = 'SUCCESS'
        result_data['energy'] = energy
        result_data['iterations'] = iterations
        result_data['computation_time'] = computation_time

        # Check against reference value if exists
        ref_value = db.get_reference(molecule_name, reference, method.upper(), basis)
        if ref_value:
            ref_energy, ref_iterations, ref_time = ref_value

            # Energy deviation
            deviation = abs(energy - ref_energy)
            result_data['deviation_from_ref'] = deviation

            # Check if deviation is within tolerance
            if deviation > ENERGY_TOLERANCE:
                result_data['error_message'] = f"Energy deviation {deviation:.2e} exceeds tolerance {ENERGY_TOLERANCE:.2e}"

            # Performance comparison
            if ref_time and ref_time > 0:
                time_ratio = computation_time / ref_time
                result_data['time_ratio'] = time_ratio
                result_data['ref_time'] = ref_time

        # Save to database
        db.save_result(result_data)

        return result_data

    except Exception as e:
        result_data['error_message'] = str(e)
        result_data['status'] = 'FAILED'

        # Save failed result to database
        db.save_result(result_data)

        return result_data

    finally:
        # Reset psi4 output to stdout (for next iteration)
        psi4.core.set_output_file('stdout', False)

def print_header():
    """Print test header"""
    print("=" * 100)
    print(f"PSI4 Universal SCF Tester")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 100)
    print()

def print_results(results):
    """Print test results in table format with validation status"""
    print("\n" + "=" * 100)
    print("VALIDATION RESULTS")
    print("=" * 100)

    # Header
    print(f"{'Molecule':<10} {'Ref':<6} {'Method':<8} {'Status':<3} {'Energy (Hartree)':<22} {'Iter':<5} {'Time(s)':<6} {'Perf':<8}")
    print("-" * 100)

    # Results
    failed_tests = []
    perf_degradations = []
    perf_improvements = []

    for result in results:
        mol = result['molecule']
        ref = result['reference']
        method = result['method']
        status = result['status']

        if status == 'SUCCESS':
            energy = f"{result['energy']:18.12f}"
            iterations = f"{result['iterations']}"
            comp_time = f"{result.get('computation_time', 0):.2f}"

            # Performance indicator
            perf_str = ""
            if 'time_ratio' in result and result['time_ratio']:
                ratio = result['time_ratio']
                if ratio < 0.95:  # >5% faster
                    perf_str = f"🚀 {ratio:.2f}x"
                    perf_improvements.append((result, ratio))
                elif ratio > 1.10:  # >10% slower
                    perf_str = f"⚠️  {ratio:.2f}x"
                    perf_degradations.append((result, ratio))
                else:
                    perf_str = f"≈ {ratio:.2f}x"

            # Green checkmark for success
            status_symbol = "\033[92m✓\033[0m"  # Green
        else:
            energy = "N/A"
            iterations = "N/A"
            comp_time = "N/A"
            perf_str = ""
            # Red cross for failed
            status_symbol = "\033[91m✗\033[0m"  # Red
            failed_tests.append(result)

        print(f"{mol:<10} {ref:<6} {method:<8} {status_symbol:<3} {energy:<22} {iterations:<5} {comp_time:<6} {perf_str:<8}")

    # Summary
    print("-" * 100)
    total = len(results)
    success = sum(1 for r in results if r['status'] == 'SUCCESS')
    failed = total - success

    print(f"Total: {total} | Success: {success} | Failed: {failed}")
    print(f"Database: {DB_FILE}")
    print(f"Logs: {LOG_DIR}/")

    # Performance summary
    if perf_improvements or perf_degradations:
        print("\n" + "=" * 100)
        print("PERFORMANCE ANALYSIS")
        print("=" * 100)

        if perf_improvements:
            print(f"\n🚀 IMPROVEMENTS ({len(perf_improvements)} tests faster):")
            for result, ratio in sorted(perf_improvements, key=lambda x: x[1])[:5]:
                speedup = 1/ratio
                print(f"  {result['molecule']}/{result['reference']}/{result['method']}: "
                      f"{speedup:.2f}x faster ({result['computation_time']*1000:.1f}ms vs {result['ref_time']*1000:.1f}ms)")

        if perf_degradations:
            print(f"\n⚠️  SLOWDOWNS ({len(perf_degradations)} tests slower):")
            for result, ratio in sorted(perf_degradations, key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {result['molecule']}/{result['reference']}/{result['method']}: "
                      f"{ratio:.2f}x slower ({result['computation_time']*1000:.1f}ms vs {result['ref_time']*1000:.1f}ms)")

    print("=" * 100)

    # Print failed test logs
    if failed_tests:
        print("\n" + "=" * 100)
        print("FAILED TESTS - LOG FILES")
        print("=" * 100)
        for result in failed_tests:
            print(f"{result['molecule']}/{result['reference']}/{result['method']}: {result.get('log_file', 'N/A')}")
            if result.get('error_message'):
                print(f"  Error: {result['error_message']}")
        print("=" * 100)

def set_references_from_results(db, results):
    """Interactively set reference values from successful test results"""
    print("\n" + "=" * 100)
    print("SET REFERENCE VALUES")
    print("=" * 100)

    # Check if running in interactive mode (terminal available)
    if not sys.stdin.isatty():
        print("Non-interactive mode detected - skipping reference value prompt")
        print("To set references, run: python test_simple.py (in interactive terminal)")
        print("=" * 100)
        return

    print("Do you want to set current successful results as reference values? (y/n): ", end='', flush=True)

    try:
        choice = input().strip().lower()
        if choice == 'y':
            count = 0
            for result in results:
                if result['status'] == 'SUCCESS':
                    db.set_reference(
                        result['molecule'],
                        result['reference'],
                        result['method'],
                        result['basis'],
                        result['energy'],
                        result['iterations'],
                        computation_time=result.get('computation_time'),
                        notes=f"Set from test run on {result['timestamp']}"
                    )
                    count += 1
            print(f"\nSet {count} reference values")
        else:
            print("Skipped setting reference values")
    except (EOFError, KeyboardInterrupt):
        print("\nSkipped setting reference values (interrupted)")

    print("=" * 100)


def main():
    """Main test runner"""
    print_header()

    # Initialize database
    db = TestDatabase(DB_FILE)

    results = []
    test_num = 0
    exit_code = 0

    try:
        # Iterate through all test configurations
        for config in TEST_CONFIGS:
            reference = config['reference']
            methods = config['methods']
            molecules = config['molecules']
            extra_options = config.get('options', {})

            # Build description string for options
            option_desc = ""
            if extra_options:
                option_desc = " [" + ", ".join([f"{k}={v}" for k, v in extra_options.items()]) + "]"

            print(f"\nTesting {reference} with molecules: {', '.join(molecules)}{option_desc}")
            print("-" * 100)

            for molecule in molecules:
                for method in methods:
                    test_num += 1
                    print(f"[{test_num}] Running: {molecule} | {reference} | {method.upper()}{option_desc}... ", end='', flush=True)

                    result = run_test(molecule, reference, method, db, extra_options)
                    results.append(result)

                    if result['status'] == 'SUCCESS':
                        deviation_str = ""
                        if result.get('deviation_from_ref') is not None:
                            dev = result['deviation_from_ref']
                            if dev > ENERGY_TOLERANCE:
                                deviation_str = f" WARNING: deviation {dev:.2e}"
                            else:
                                deviation_str = f" deviation {dev:.2e}"

                        print(f"OK E = {result['energy']:.10f} ({result['iterations']} iter, {result['computation_time']:.2f}s){deviation_str}")
                    else:
                        print(f"FAILED: {result.get('error_message', 'Unknown error')}")

        # Print final results
        print_results(results)

        # Offer to set references
        set_references_from_results(db, results)

        # Check if there were failures
        failed_count = sum(1 for r in results if r['status'] != 'SUCCESS')
        if failed_count > 0:
            exit_code = 1

    finally:
        # Close database
        db.close()

        # Clean up Psi4 resources (threads, memory, file handles)
        psi4.core.clean()
        psi4.core.clean_options()

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return exit_code


if __name__ == '__main__':
    import gc
    import os

    try:
        exit_code = main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        exit_code = 130
    except Exception as e:
        print(f"\n\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        exit_code = 1
    finally:
        # Final cleanup: release all Psi4 resources
        try:
            psi4.core.clean()
            psi4.core.clean_options()
        except:
            pass

    # Flush output buffers before os._exit()
    # os._exit() bypasses normal cleanup, including buffer flushing
    sys.stdout.flush()
    sys.stderr.flush()

    # Force immediate exit without waiting for cleanup
    # This is necessary because Psi4 may leave background threads/resources
    os._exit(exit_code)
