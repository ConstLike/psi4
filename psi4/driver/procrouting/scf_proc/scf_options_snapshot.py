"""
Options Snapshot Pattern for Multi-SCF

Solves non-determinism bug where global options are modified between wfn creations,
causing different SCF convergence behavior.

Problem:
    wfn1 created → reads global DIIS_START=0
    baseline test runs → changes global DIIS_START=14
    wfn2 created → reads global DIIS_START=14  # Different!
    multi_scf([wfn1, wfn2]) → non-deterministic behavior

Solution:
    Freeze global options once BEFORE creating all wfn, then apply frozen
    snapshot to each wfn independently.
"""

from psi4 import core


def snapshot_scf_options():
    """
    Freeze current global SCF options into a dictionary.

    This captures the state of all SCF options that affect convergence behavior,
    preventing pollution from global state changes between wfn creation.

    Returns
    -------
    dict
        Dictionary of option_name -> option_value for all critical SCF options

    Example
    -------
    >>> snapshot = snapshot_scf_options()
    >>> snapshot['DIIS_START']
    0
    """
    snapshot = {}

    # ===== DIIS Options (CRITICAL - main source of non-determinism) =====
    snapshot['DIIS_START'] = core.get_option('SCF', 'DIIS_START')
    snapshot['DIIS_MIN_VECS'] = core.get_option('SCF', 'DIIS_MIN_VECS')
    snapshot['DIIS_MAX_VECS'] = core.get_option('SCF', 'DIIS_MAX_VECS')
    snapshot['DIIS_RMS_ERROR'] = core.get_option('SCF', 'DIIS_RMS_ERROR')

    # ===== Damping Options =====
    snapshot['DAMPING_PERCENTAGE'] = core.get_option('SCF', 'DAMPING_PERCENTAGE')
    snapshot['DAMPING_CONVERGENCE'] = core.get_option('SCF', 'DAMPING_CONVERGENCE')

    # ===== SOSCF Options =====
    snapshot['SOSCF_START_CONVERGENCE'] = core.get_option('SCF', 'SOSCF_START_CONVERGENCE')
    snapshot['SOSCF_CONV'] = core.get_option('SCF', 'SOSCF_CONV')
    snapshot['SOSCF_MIN_ITER'] = core.get_option('SCF', 'SOSCF_MIN_ITER')
    snapshot['SOSCF_MAX_ITER'] = core.get_option('SCF', 'SOSCF_MAX_ITER')

    # ===== MOM Options =====
    snapshot['MOM_START'] = core.get_option('SCF', 'MOM_START')
    snapshot['MOM_OCC'] = core.get_option('SCF', 'MOM_OCC')

    # ===== FRAC Options =====
    snapshot['FRAC_START'] = core.get_option('SCF', 'FRAC_START')
    snapshot['FRAC_OCC'] = core.get_option('SCF', 'FRAC_OCC')
    snapshot['FRAC_VAL'] = core.get_option('SCF', 'FRAC_VAL')
    snapshot['FRAC_RENORMALIZE'] = core.get_option('SCF', 'FRAC_RENORMALIZE')

    # ===== Convergence Options =====
    snapshot['MAXITER'] = core.get_option('SCF', 'MAXITER')
    snapshot['E_CONVERGENCE'] = core.get_option('SCF', 'E_CONVERGENCE')
    snapshot['D_CONVERGENCE'] = core.get_option('SCF', 'D_CONVERGENCE')
    snapshot['FAIL_ON_MAXITER'] = core.get_option('SCF', 'FAIL_ON_MAXITER')

    # ===== Other Options =====
    snapshot['LEVEL_SHIFT'] = core.get_option('SCF', 'LEVEL_SHIFT')
    snapshot['PRINT'] = core.get_option('SCF', 'PRINT')

    # ===== COSX Options =====
    snapshot['COSX_MAXITER_FINAL'] = core.get_option('SCF', 'COSX_MAXITER_FINAL')

    return snapshot


def apply_options_snapshot(wfn, snapshot):
    """
    Apply frozen options snapshot to a wavefunction.

    This stores the snapshot in wfn._options_snapshot so that
    _scf_initialize_iteration_state() can read from it instead of
    reading from global state.

    Parameters
    ----------
    wfn : psi4.core.HF
        Wavefunction object (RHF, UHF, ROHF, etc.)
    snapshot : dict
        Dictionary from snapshot_scf_options()

    Example
    -------
    >>> snapshot = snapshot_scf_options()
    >>> wfn1 = core.RHF(ref_wfn, functional)
    >>> apply_options_snapshot(wfn1, snapshot)
    >>> # wfn1 will now use frozen options, immune to global changes
    """
    # Store snapshot in wfn for _scf_initialize_iteration_state() to use
    wfn._options_snapshot = snapshot


def get_option_from_snapshot(wfn, option_name, fallback_scope='SCF'):
    """
    Helper to read option from snapshot or fall back to global.

    This provides backward compatibility: if wfn doesn't have snapshot
    (single-cycle SCF before migration), read from global as before.

    Parameters
    ----------
    wfn : psi4.core.HF
        Wavefunction object
    option_name : str
        Name of the option (e.g., 'DIIS_START')
    fallback_scope : str
        Scope to use for core.get_option() fallback (default: 'SCF')

    Returns
    -------
    value
        Option value from snapshot or global state

    Example
    -------
    >>> diis_start = get_option_from_snapshot(wfn, 'DIIS_START')
    """
    if hasattr(wfn, '_options_snapshot') and wfn._options_snapshot is not None:
        return wfn._options_snapshot[option_name]
    else:
        # Fallback to global (backward compatibility for single-cycle SCF)
        return core.get_option(fallback_scope, option_name)
