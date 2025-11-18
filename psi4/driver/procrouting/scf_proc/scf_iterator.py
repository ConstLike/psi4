
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
The SCF iteration functions
"""
import numpy as np

from psi4 import core

from ... import p4util
from ...constants import constants
from ...p4util.exceptions import SCFConvergenceError, ValidationError
from ..solvent.efp import get_qm_atoms_opts, modify_Fock_induced, modify_Fock_permanent
from .scf_options_snapshot import get_option_from_snapshot, snapshot_scf_options, apply_options_snapshot

#import logging
#logger = logging.getLogger("scf.scf_iterator")
#logger.setLevel(logging.DEBUG)

# Q: I expect more local settings of options for part of SCF.
#    For convcrit, do we want:
#   (A) easy to grep
#    with p4util.OptionsStateCM(['SCF', 'E_CONVERGENCE'], ['SCF', 'D_CONVERGENCE']):
#        core.set_local_option('SCF', 'E_CONVERGENCE', 1.e-5)
#        core.set_local_option('SCF', 'D_CONVERGENCE', 1.e-4)
#        self.iterations()
#
#   or (B) functional. options never touched
#    self.iterations(e_conv=1.e-5, d_conv=1.e-4)


def scf_compute_energy(self):
    """
    Compute SCF energy using multi_scf() framework.

    This is a wrapper that calls multi_scf([self]) to ensure single SCF
    uses the same code path as multi-cycle SCF. This ensures all features
    (DIIS, damping, SOSCF, MOM, FRAC, DF_SCF_GUESS, etc.) are automatically
    available in both modes.

    Single source of truth: All SCF calculations flow through multi_scf().

    Returns
    -------
    float
        SCF energy computed by finalize_energy()
    """
    # Initialize iteration_energies for backward compatibility
    self.iteration_energies = []

    # Single SCF is just multi_scf with one wavefunction
    # This ensures identical behavior and eliminates code duplication
    try:
        energies = multi_scf([self], verbose=True)
        scf_energy = energies[0]
    except SCFConvergenceError as e:
        if core.get_option("SCF", "FAIL_ON_MAXITER"):
            core.print_out("  Failed to converge.\n")
            raise e
        else:
            core.print_out("  Energy and/or wave function did not converge, but proceeding anyway.\n\n")
            scf_energy = self.get_energies("Total Energy")
    else:
        core.print_out("  Energy and wave function converged.\n\n")

    # Finalize and return
    scf_energy = self.finalize_energy()
    return scf_energy


def _build_jk(wfn, memory):
    jk = core.JK.build(wfn.get_basisset("ORBITAL"),
                       aux=wfn.get_basisset("DF_BASIS_SCF"),
                       do_wK=wfn.functional().is_x_lrc(),
                       memory=memory)
    return jk


def initialize_jk(self, memory, jk=None):

    functional = self.functional()
    if jk is None:
        jk = _build_jk(self, memory)

    self.set_jk(jk)

    jk.set_print(self.get_print())
    jk.set_memory(memory)
    jk.set_do_K(functional.is_x_hybrid())
    jk.set_do_wK(functional.is_x_lrc())
    jk.set_omega(functional.x_omega())

    jk.set_omega_alpha(functional.x_alpha())
    jk.set_omega_beta(functional.x_beta())

    jk.initialize()
    jk.print_header()


def scf_initialize(self):
    """Specialized initialization, compute integrals and does everything to prepare for iterations"""

    # Figure out memory distributions

    # Get memory in terms of doubles
    total_memory = (core.get_memory() / 8) * core.get_global_option("SCF_MEM_SAFETY_FACTOR")

    # Figure out how large the DFT collocation matrices are
    vbase = self.V_potential()
    if vbase:
        collocation_size = vbase.grid().collocation_size()
        if vbase.functional().ansatz() == 1:
            collocation_size *= 4  # First derivs
        elif vbase.functional().ansatz() == 2:
            collocation_size *= 10  # Second derivs
    else:
        collocation_size = 0

    # Change allocation for collocation matrices based on DFT type
    initialize_jk_obj = False
    if isinstance(self.jk(), core.JK):
        core.print_out("\nRe-using passed JK object instead of rebuilding\n")
        jk = self.jk()
    else:
        initialize_jk_obj = True
        jk = _build_jk(self, total_memory)
    jk_size = jk.memory_estimate()

    # Give remaining to collocation
    if total_memory > jk_size:
        collocation_memory = total_memory - jk_size
    # Give up to 10% to collocation
    elif (total_memory * 0.1) > collocation_size:
        collocation_memory = collocation_size
    else:
        collocation_memory = total_memory * 0.1

    if collocation_memory > collocation_size:
        collocation_memory = collocation_size

    # Set constants
    self.iteration_ = 0
    self.memory_jk_ = int(total_memory - collocation_memory)
    self.memory_collocation_ = int(collocation_memory)

    if self.get_print():
        core.print_out("  ==> Integral Setup <==\n\n")

    # Initialize EFP
    efp_enabled = hasattr(self.molecule(), 'EFP')
    if efp_enabled:
        # EFP: Set QM system, options, and callback. Display efp geom in [A]
        efpobj = self.molecule().EFP
        core.print_out(efpobj.banner())
        core.print_out(efpobj.geometry_summary(units_to_bohr=constants.bohr2angstroms))

        efpptc, efpcoords, efpopts = get_qm_atoms_opts(self.molecule())
        efpobj.set_point_charges(efpptc, efpcoords)
        efpobj.set_opts(efpopts, label='psi', append='psi')

        efpobj.set_electron_density_field_fn(efp_field_fn)

    # Initialize all integrals and perform the first guess
    if self.attempt_number_ == 1:
        mints = core.MintsHelper(self.basisset())

        if initialize_jk_obj:
            self.initialize_jk(self.memory_jk_, jk=jk)
        if self.V_potential():
            self.V_potential().build_collocation_cache(self.memory_collocation_)
        core.timer_on("HF: Form core H")
        self.form_H()
        core.timer_off("HF: Form core H")

        if efp_enabled:
            # EFP: Add in permanent moment contribution and cache
            core.timer_on("HF: Form Vefp")
            verbose = core.get_option('SCF', "PRINT")
            Vefp = modify_Fock_permanent(self.molecule(), mints, verbose=verbose - 1)
            Vefp = core.Matrix.from_array(Vefp)
            self.H().add(Vefp)
            Horig = self.H().clone()
            self.Horig = Horig
            core.print_out("  QM/EFP: iterating Total Energy including QM/EFP Induction\n")
            core.timer_off("HF: Form Vefp")

        core.timer_on("HF: Form S/X")
        self.form_Shalf()
        core.timer_off("HF: Form S/X")

        core.print_out("\n  ==> Pre-Iterations <==\n\n")

        # force SCF_SUBTYPE to AUTO during SCF guess
        optstash = p4util.OptionsState(["SCF", "SCF_SUBTYPE"])
        core.set_local_option("SCF", "SCF_SUBTYPE", "AUTO")

        core.timer_on("HF: Guess")
        self.guess()
        core.timer_off("HF: Guess")

        optstash.restore()

        # Print out initial docc/socc/etc data
        if self.get_print():
            lack_occupancy = core.get_local_option('SCF', 'GUESS') in ['SAD']
            if core.get_global_option('GUESS') in ['SAD']:
                lack_occupancy = core.get_local_option('SCF', 'GUESS') in ['AUTO']
                self.print_preiterations(small=lack_occupancy)
            else:
                self.print_preiterations(small=lack_occupancy)

    else:
        # We're reading the orbitals from the previous set of iterations.
        self.form_D()
        self.set_energies("Total Energy", self.compute_initial_E())

    # turn off VV10 for iterations
    if core.get_option('SCF', "DFT_VV10_POSTSCF") and self.functional().vv10_b() > 0.0:
        core.print_out("  VV10: post-SCF option active \n \n")
        self.functional().set_lock(False)
        self.functional().set_do_vv10(False)
        self.functional().set_lock(True)

    # Print iteration header
    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    diis_rms = get_option_from_snapshot(self, 'DIIS_RMS_ERROR')
    core.print_out("  ==> Iterations <==\n\n")
    core.print_out("%s                        Total Energy        Delta E     %s |[F,P]|\n\n" %
                   ("   " if is_dfjk else "", "RMS" if diis_rms else "MAX"))


def _scf_initialize_iteration_state(self, e_conv, d_conv):
    """
    Initialize state for SCF iterations.

    This method sets up all parameters and state variables needed for the
    SCF iteration loop, storing them as self._scf_* members to enable
    multi-cycle SCF coordination.

    Parameters
    ----------
    e_conv : float or None
        Energy convergence threshold
    d_conv : float or None
        Density convergence threshold
    """
    # Store convergence criteria
    self._scf_e_conv = e_conv
    self._scf_d_conv = d_conv

    # ========================================================================
    # Cache ALL options from snapshot (ONCE) to avoid scattered lookups
    # ========================================================================
    # This consolidates all option reads in ONE place, making code cleaner
    # and easier to maintain. Previously had 34 scattered get_option_from_snapshot()
    # calls throughout iteration code.

    # DIIS options
    self._scf_diis = get_option_from_snapshot(self, 'DIIS')
    self._scf_diis_start = get_option_from_snapshot(self, 'DIIS_START')
    self._scf_diis_min_vecs = get_option_from_snapshot(self, 'DIIS_MIN_VECS')
    self._scf_diis_max_vecs = get_option_from_snapshot(self, 'DIIS_MAX_VECS')
    self._scf_diis_rms_error = get_option_from_snapshot(self, 'DIIS_RMS_ERROR')

    # AEDIIS options
    self._scf_initial_accelerator = get_option_from_snapshot(self, 'SCF_INITIAL_ACCELERATOR')
    self._scf_initial_start_diis_transition = get_option_from_snapshot(self, 'SCF_INITIAL_START_DIIS_TRANSITION')
    self._scf_initial_finish_diis_transition = get_option_from_snapshot(self, 'SCF_INITIAL_FINISH_DIIS_TRANSITION')

    # Damping options
    self._scf_damping_percentage = get_option_from_snapshot(self, 'DAMPING_PERCENTAGE')
    self._scf_damping_convergence = get_option_from_snapshot(self, 'DAMPING_CONVERGENCE')

    # SOSCF options
    self._scf_soscf = get_option_from_snapshot(self, 'SOSCF')
    self._scf_soscf_start_convergence = get_option_from_snapshot(self, 'SOSCF_START_CONVERGENCE')
    self._scf_soscf_conv = get_option_from_snapshot(self, 'SOSCF_CONV')
    self._scf_soscf_min_iter = get_option_from_snapshot(self, 'SOSCF_MIN_ITER')
    self._scf_soscf_max_iter = get_option_from_snapshot(self, 'SOSCF_MAX_ITER')
    self._scf_soscf_print = get_option_from_snapshot(self, 'SOSCF_PRINT')

    # MOM options
    self._scf_mom_start = get_option_from_snapshot(self, 'MOM_START')
    self._scf_mom_occ = get_option_from_snapshot(self, 'MOM_OCC')

    # FRAC options
    self._scf_frac_start = get_option_from_snapshot(self, 'FRAC_START')
    self._scf_frac_occ = get_option_from_snapshot(self, 'FRAC_OCC')
    self._scf_frac_val = get_option_from_snapshot(self, 'FRAC_VAL')
    self._scf_frac_renormalize = get_option_from_snapshot(self, 'FRAC_RENORMALIZE')

    # Convergence options
    self._scf_maxiter = get_option_from_snapshot(self, 'MAXITER')
    self._scf_e_convergence = get_option_from_snapshot(self, 'E_CONVERGENCE')
    self._scf_d_convergence = get_option_from_snapshot(self, 'D_CONVERGENCE')
    self._scf_fail_on_maxiter = get_option_from_snapshot(self, 'FAIL_ON_MAXITER')

    # Other options
    self._scf_level_shift = get_option_from_snapshot(self, 'LEVEL_SHIFT')
    self._scf_level_shift_cutoff = get_option_from_snapshot(self, 'LEVEL_SHIFT_CUTOFF')
    self._scf_print = get_option_from_snapshot(self, 'PRINT')

    # COSX options
    self._scf_cosx_maxiter_final = get_option_from_snapshot(self, 'COSX_MAXITER_FINAL')

    # ========================================================================
    # Derived configuration flags (computed from cached options)
    # ========================================================================

    self._scf_is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    self._scf_verbose = self._scf_print
    self._scf_reference = core.get_option('SCF', "REFERENCE")  # OK from global (doesn't change)
    self._scf_damping_enabled = _validate_damping(self)
    self._scf_soscf_enabled = _validate_soscf(self)
    self._scf_frac_enabled = _validate_frac(self)
    self._scf_efp_enabled = hasattr(self.molecule(), 'EFP')
    self._scf_cosx_enabled = "COSX" in core.get_option('SCF', 'SCF_TYPE')  # OK from global

    # Set DIIS/MOM members used by iteration loop
    # Required by both scf_iterate() and multi_scf()
    self.diis_enabled_ = self.validate_diis()
    self.MOM_excited_ = _validate_MOM(self)
    self.diis_start_ = self._scf_diis_start

    # COSX early screening parameters
    self._scf_early_screening = False
    if self._scf_cosx_enabled:
        self._scf_early_screening = True
        self.jk().set_COSX_grid("Initial")

    self._scf_maxiter_post_screening = self._scf_cosx_maxiter_final
    if self._scf_maxiter_post_screening < -1:
        raise ValidationError('COSX_MAXITER_FINAL ({}) must be -1 or above. If you wish to attempt full SCF converge on the final COSX grid, set COSX_MAXITER_FINAL to -1.'.format(self._scf_maxiter_post_screening))

    self._scf_early_screening_disabled = False

    # Initialize iteration state variables
    self._scf_SCFE_old = 0.0
    self._scf_Dnorm = 0.0
    self._scf_Ediff = 0.0
    self._scf_iter_post_screening = 0


def scf_iterate(self, e_conv=None, d_conv=None):

    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    verbose = core.get_option('SCF', "PRINT")
    reference = core.get_option('SCF', "REFERENCE")

    # self.member_data_ signals are non-local, used internally by c-side fns
    self.diis_enabled_ = self.validate_diis()
    self.MOM_excited_ = _validate_MOM(self)
    self.diis_start_ = get_option_from_snapshot(self, 'DIIS_START')
    damping_enabled = _validate_damping(self)
    soscf_enabled = _validate_soscf(self)
    frac_enabled = _validate_frac(self)
    efp_enabled = hasattr(self.molecule(), 'EFP')
    cosx_enabled = "COSX" in core.get_option('SCF', 'SCF_TYPE')
    ooo_scf = core.get_option("SCF", "ORBITAL_OPTIMIZER_PACKAGE") in ["OOO", "OPENORBITALOPTIMIZER"]
    if ooo_scf:
        pcm_enabled = core.get_option('SCF', 'PCM')
        ddx_enabled = core.get_option('SCF', 'DDX')
        pe_enabled = core.get_option('SCF', 'PE')
        level_shift_enabled = core.get_option("SCF", "LEVEL_SHIFT") != 0.0
        autograc_enabled = core.get_option("SAPT", "SAPT_DFT_GRAC_COMPUTE") != "NONE"
        guessmix_enabled = core.get_option("SCF", "GUESS_MIX")
        if (reference in ["ROHF", "CUHF"] or soscf_enabled or self.MOM_excited_ or frac_enabled or
            efp_enabled or pcm_enabled or ddx_enabled or pe_enabled or autograc_enabled or
            level_shift_enabled or guessmix_enabled):
            core.print_out(f"    Note: OpenOrbitalOptimizer not compatible with at least one of the following. Falling back to orbital_optimizer_package=internal\n")
            core.print_out(f"          {reference=}, soscf={soscf_enabled}, mom={self.MOM_excited_}, frac={frac_enabled}, efp={efp_enabled},\n")
            core.print_out(f"          pcm={pcm_enabled}, ddx={ddx_enabled}, pe={pe_enabled}, autograc={autograc_enabled}, level_shift={level_shift_enabled},\n")
            core.print_out(f"          guess_mix={guessmix_enabled}\n")
        else:
            # SAD needs some special work since the guess doesn't actually make the orbitals in Psi4
            if self.sad_ and self.iteration_ <= 0:
                self.iteration_ += 1
                self.form_G()
                self.form_initial_F()
                self.form_initial_C()
                self.reset_occupation()
                self.find_occupation()
                ene_sad = self.compute_E()
                core.print_out(
                    "   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s\n" %
                    ("DF-" if is_dfjk else "", reference, "SAD", ene_sad, ene_sad, 0.0, ""))
            if core.get_option("SCF", "GUESS") == "READ" and self.iteration_ <= 0:
                self.form_G()
                self.form_initial_F()
                self.form_initial_C()
                self.reset_occupation()
                self.find_occupation()
                ene_sad = self.compute_E()

            try:
                self.openorbital_scf()
            except RuntimeError as ex:
                if "openorbital_scf is virtual; it has not been implemented for your class" in str(ex):
                    core.print_out(f"    Note: OpenOrbitalOptimizer NYI for {reference}. Falling back to Internal.\n")
                else:
                    raise ex
            else:
                SCFE = self.compute_E()
                self.set_energies("Total Energy", SCFE)
                self.set_variable("SCF ITERATION ENERGY", SCFE)
                self.iteration_energies.append(SCFE)  # note 1-len array, not niter-len array like INTERNAL

                self.form_G()
                self.form_F()
                self.form_C()
                self.form_D()
                return

    # Initialize iteration state and start iterations
    self._scf_initialize_iteration_state(e_conv, d_conv)

    # Main SCF iteration loop
    while True:
        should_continue, reason = self._scf_iteration()

        if not should_continue:
            break

        if self.iteration_ >= self._scf_maxiter:
            raise SCFConvergenceError("""SCF iterations""", self.iteration_, self, self._scf_Ediff, self._scf_Dnorm)


def _scf_iteration(self):
    """
    Performs ONE SCF iteration.

    This method is designed to be called externally by a multi-cycle SCF
    coordinator. It uses state stored in self._scf_* members initialized
    by _scf_initialize_iteration_state().

    Returns
    -------
    tuple: (continue_flag, reason)
        continue_flag: True to continue iterations, False to stop
        reason: 'converged' | 'max_iter' | 'early_screening_maxiter' | 'mom_not_started' | 'frac_not_started'
    """
    self.iteration_ += 1

    diis_performed = False
    soscf_performed = False
    self.frac_performed_ = False
    #self.MOM_performed_ = False  # redundant from common_init()

    self.save_density_and_energy()

    if self._scf_efp_enabled:
        # EFP: Add efp contribution to Fock matrix
        self.H().copy(self.Horig)
        global mints_psi4_yo
        mints_psi4_yo = core.MintsHelper(self.basisset())
        Vefp = modify_Fock_induced(self.molecule().EFP, mints_psi4_yo, verbose=self._scf_verbose - 1)
        Vefp = core.Matrix.from_array(Vefp)
        self.H().add(Vefp)

    SCFE = 0.0
    self.clear_external_potentials()

    # Two-electron contribution to Fock matrix from self.jk()
    core.timer_on("HF: Form G")
    self.form_G()
    core.timer_off("HF: Form G")

    # Check if special J/K construction algorithms were used
    incfock_performed = hasattr(self.jk(), "do_incfock_iter") and self.jk().do_incfock_iter()
    upcm = 0.0
    if core.get_option('SCF', 'PCM'):
        calc_type = core.PCM.CalcType.Total
        if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
            calc_type = core.PCM.CalcType.NucAndEle
        Dt = self.Da().clone()
        Dt.add(self.Db())
        upcm, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
        SCFE += upcm
        self.push_back_external_potential(Vpcm)
    self.set_variable("PCM POLARIZATION ENERGY", upcm)  # P::e PCM
    self.set_energies("PCM Polarization", upcm)

    uddx = 0.0
    if core.get_option('SCF', 'DDX'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        uddx, Vddx, self.ddx_state = self.ddx.get_solvation_contributions(Dt, self.ddx_state)
        SCFE += uddx
        self.push_back_external_potential(Vddx)
    self.set_variable("DD SOLVATION ENERGY", uddx)  # P::e DDX
    self.set_energies("DD Solvation Energy", uddx)

    upe = 0.0
    if core.get_option('SCF', 'PE'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        upe, Vpe = self.pe_state.get_pe_contribution(
            Dt, elec_only=False
        )
        SCFE += upe
        self.push_back_external_potential(Vpe)
    self.set_variable("PE ENERGY", upe)  # P::e PE
    self.set_energies("PE Energy", upe)

    core.timer_on("HF: Form F")
    # SAD: since we don't have orbitals yet, we might not be able
    # to form the real Fock matrix. Instead, build an initial one
    if (self.iteration_ == 0) and self.sad_:
        self.form_initial_F()
    else:
        self.form_F()
    core.timer_off("HF: Form F")

    if self._scf_verbose > 3:
        self.Fa().print_out()
        self.Fb().print_out()

    SCFE += self.compute_E()
    if self._scf_efp_enabled:
        global efp_Dt_psi4_yo

        # EFP: Add efp contribution to energy
        efp_Dt_psi4_yo = self.Da().clone()
        efp_Dt_psi4_yo.add(self.Db())
        SCFE += self.molecule().EFP.get_wavefunction_dependent_energy()

    self.set_energies("Total Energy", SCFE)
    core.set_variable("SCF ITERATION ENERGY", SCFE)
    self.iteration_energies.append(SCFE)

    self._scf_Ediff = SCFE - self._scf_SCFE_old
    self._scf_SCFE_old = SCFE

    status = []

    # Check if we are doing SOSCF
    if (self._scf_soscf_enabled and (self.iteration_ >= 3) and (self._scf_Dnorm < self._scf_soscf_start_convergence)):
        self._scf_Dnorm = self.compute_orbital_gradient(False, self._scf_diis_max_vecs)
        diis_performed = False
        if self.functional().needs_xc():
            base_name = "SOKS, nmicro="
        else:
            base_name = "SOSCF, nmicro="

        if not _converged(self._scf_Ediff, self._scf_Dnorm, e_conv=self._scf_e_conv, d_conv=self._scf_d_conv):
            nmicro = self.soscf_update(self._scf_soscf_conv,
                                       self._scf_soscf_min_iter,
                                       self._scf_soscf_max_iter,
                                       self._scf_soscf_print)
            # if zero, the soscf call bounced for some reason
            soscf_performed = (nmicro > 0)

            if soscf_performed:
                self.find_occupation()
                status.append(base_name + str(nmicro))
            else:
                if self._scf_verbose > 0:
                    core.print_out("Did not take a SOSCF step, using normal convergence methods\n")

        else:
            # need to ensure orthogonal orbitals and set epsilon
            status.append(base_name + "conv")
            core.timer_on("HF: Form C")
            self.form_C()
            core.timer_off("HF: Form C")
            soscf_performed = True  # Stops DIIS

    if not soscf_performed:
        # Normal convergence procedures if we do not do SOSCF

        # SAD: form initial orbitals from the initial Fock matrix, and
        # reset the occupations. The reset is necessary because SAD
        # nalpha_ and nbeta_ are not guaranteed physical.
        # From here on, the density matrices are correct.
        if (self.iteration_ == 0) and self.sad_:
            self.form_initial_C()
            self.reset_occupation()
            self.find_occupation()

        else:
            # Run DIIS
            core.timer_on("HF: DIIS")
            try:
                diis_performed = False
                add_to_diis_subspace = self.diis_enabled_ and self.iteration_ >= self.diis_start_

                self._scf_Dnorm = self.compute_orbital_gradient(add_to_diis_subspace, self._scf_diis_max_vecs)

                if add_to_diis_subspace:
                    for engine_used in self.diis(self._scf_Dnorm):
                        status.append(engine_used)
            finally:
                # Turn off timer even if exception occurs to prevent timer stack errors
                core.timer_off("HF: DIIS")

            if self._scf_verbose > 4 and diis_performed:
                core.print_out("  After DIIS:\n")
                self.Fa().print_out()
                self.Fb().print_out()

            # frac, MOM invoked here from Wfn::HF::find_occupation
            core.timer_on("HF: Form C")
            level_shift = self._scf_level_shift
            if level_shift > 0 and self._scf_Dnorm > self._scf_level_shift_cutoff:
                status.append("SHIFT")
                self.form_C(level_shift)
            else:
                self.form_C()
            core.timer_off("HF: Form C")

            if self.MOM_performed_:
                status.append("MOM")

            if self.frac_performed_:
                status.append("FRAC")

            if incfock_performed:
                status.append("INCFOCK")

            # Reset occupations if necessary
            if (self.iteration_ == 0) and self.reset_occ_:
                self.reset_occupation()
                self.find_occupation()

    # Form new density matrix
    core.timer_on("HF: Form D")
    self.form_D()
    core.timer_off("HF: Form D")

    self.set_variable("SCF ITERATION ENERGY", SCFE)
    core.set_variable("SCF D NORM", self._scf_Dnorm)

    # After we've built the new D, damp the update
    if (self._scf_damping_enabled and self.iteration_ > 1 and self._scf_Dnorm > self._scf_damping_convergence):
        damping_percentage = self._scf_damping_percentage
        self.damping_update(damping_percentage * 0.01)
        status.append("DAMP={}%".format(round(damping_percentage)))

    if core.has_option_changed("SCF", "ORBITALS_WRITE"):
        filename = core.get_option("SCF", "ORBITALS_WRITE")
        # Use wfn_name for unique filename in multi-SCF (empty wfn_name preserves original behavior)
        unique_filename = self.get_orbitals_filename(filename)
        self.to_file(unique_filename)

    if self._scf_verbose > 3:
        self.Ca().print_out()
        self.Cb().print_out()
        self.Da().print_out()
        self.Db().print_out()

    # Print out the iteration
    core.print_out(
        "   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s\n" %
        ("DF-" if self._scf_is_dfjk else "", self._scf_reference, "SAD" if
         ((self.iteration_ == 0) and self.sad_) else self.iteration_, SCFE, self._scf_Ediff, self._scf_Dnorm, '/'.join(status)))

    # if a an excited MOM is requested but not started, don't stop yet
    # Note that MOM_performed_ just checks initialization, and our convergence measures used the pre-MOM orbitals
    if self.MOM_excited_ and ((not self.MOM_performed_) or self.iteration_ == self._scf_mom_start):
        return (True, 'mom_not_started')

    # if a fractional occupation is requested but not started, don't stop yet
    if self._scf_frac_enabled and not self.frac_performed_:
        return (True, 'frac_not_started')

    # have we completed our post-early screening SCF iterations?
    if self._scf_early_screening_disabled:
        self._scf_iter_post_screening += 1
        if self._scf_iter_post_screening >= self._scf_maxiter_post_screening and self._scf_maxiter_post_screening > 0:
            return (False, 'early_screening_maxiter')

    # Call any postiteration callbacks
    if not ((self.iteration_ == 0) and self.sad_) and _converged(self._scf_Ediff, self._scf_Dnorm, e_conv=self._scf_e_conv, d_conv=self._scf_d_conv):

        if self._scf_early_screening:

            # we've reached convergence with early screning enabled; disable it
            self._scf_early_screening = False

            # make note of the change to early screening; next SCF iteration(s) will be the last
            self._scf_early_screening_disabled = True

            # cosx uses the largest grid for its final SCF iteration(s)
            if self._scf_cosx_enabled:
                self.jk().set_COSX_grid("Final")

            # clear any cached matrices associated with incremental fock construction
            # the change in the screening spoils the linearity in the density matrix
            if hasattr(self.jk(), 'clear_D_prev'):
                self.jk().clear_D_prev()

            if self._scf_maxiter_post_screening == 0:
                return (False, 'converged')
            else:
                core.print_out("  Energy and wave function converged with early screening.\n")
                core.print_out("  Continuing SCF iterations with tighter screening.\n\n")
        else:
            return (False, 'converged')

    # Continue iterating
    return (True, 'continue')


def scf_finalize_energy(self):
    """Performs stability analysis and calls back SCF with new guess
    if needed, Returns the SCF energy. This function should be called
    once orbitals are ready for energy/property computations, usually
    after iterations() is called.

    """

    # post-scf vv10 correlation
    if core.get_option('SCF', "DFT_VV10_POSTSCF") and self.functional().vv10_b() > 0.0:
        self.functional().set_lock(False)
        self.functional().set_do_vv10(True)
        self.functional().set_lock(True)
        core.print_out("  ==> Computing Non-Self-Consistent VV10 Energy Correction <==\n\n")
        SCFE = 0.0
        self.form_V()
        SCFE += self.compute_E()
        self.set_energies("Total Energy", SCFE)

    # Perform wavefunction stability analysis before doing
    # anything on a wavefunction that may not be truly converged.
    if core.get_option('SCF', 'STABILITY_ANALYSIS') != "NONE":

        # We need the integral file, make sure it is written and
        # compute it if needed
        if core.get_option('SCF', 'REFERENCE') not in {"UHF", "UKS"}:
            # Don't bother computing needed integrals if we can't do anything with them.
            if self.functional().needs_xc():
                raise ValidationError("Stability analysis not yet supported for XC functionals.")

            #psio = core.IO.shared_object()
            #psio.open(constants.PSIF_SO_TEI, 1)  # PSIO_OPEN_OLD
            #try:
            #    psio.tocscan(constants.PSIF_SO_TEI, "IWL Buffers")
            #except TypeError:
            #    # "IWL Buffers" actually found but psio_tocentry can't be returned to Py
            #    psio.close(constants.PSIF_SO_TEI, 1)
            #else:
            #    # tocscan returned None
            #    psio.close(constants.PSIF_SO_TEI, 1)

            # logic above foiled by psio_tocentry not returning None<--nullptr in pb11 2.2.1
            #   so forcibly recomputing for now until stability revamp
            core.print_out("    SO Integrals not on disk. Computing...")
            mints = core.MintsHelper(self.basisset())

            mints.integrals()
            core.print_out("done.\n")

            # Q: Not worth exporting all the layers of psio, right?

        follow = self.stability_analysis()

        while follow and self.attempt_number_ <= core.get_option('SCF', 'MAX_ATTEMPTS'):
            self.attempt_number_ += 1
            core.print_out("    Running SCF again with the rotated orbitals.\n")

            if self.initialized_diis_manager_:
                self.diis_manager_.reset_subspace()
            # reading the rotated orbitals in before starting iterations
            self.form_D()
            self.set_energies("Total Energy", self.compute_initial_E())
            self.iterations()
            follow = self.stability_analysis()

        if follow and self.attempt_number_ > core.get_option('SCF', 'MAX_ATTEMPTS'):
            core.print_out("    There's still a negative eigenvalue. Try modifying FOLLOW_STEP_SCALE\n")
            core.print_out("    or increasing MAX_ATTEMPTS (not available for PK integrals).\n")

    # At this point, we are not doing any more SCF cycles
    #   and we can compute and print final quantities.

    if hasattr(self.molecule(), 'EFP'):
        efpobj = self.molecule().EFP

        efpobj.compute()  # do_gradient=do_gradient)
        efpene = efpobj.get_energy(label='psi')
        efp_wfn_independent_energy = efpene['total'] - efpene['ind']
        self.set_energies("EFP", efpene['total'])

        SCFE = self.get_energies("Total Energy")
        SCFE += efp_wfn_independent_energy
        self.set_energies("Total Energy", SCFE)
        core.print_out(efpobj.energy_summary(scfefp=SCFE, label='psi'))

        self.set_variable("EFP ELST ENERGY", efpene['electrostatic'] + efpene['charge_penetration'] + efpene['electrostatic_point_charges'])  # P::e EFP
        self.set_variable("EFP IND ENERGY", efpene['polarization'])  # P::e EFP
        self.set_variable("EFP DISP ENERGY", efpene['dispersion'])  # P::e EFP
        self.set_variable("EFP EXCH ENERGY", efpene['exchange_repulsion'])  # P::e EFP
        self.set_variable("EFP TOTAL ENERGY", efpene['total'])  # P::e EFP
        self.set_variable("CURRENT ENERGY", efpene['total'])  # P::e EFP

    core.print_out("\n  ==> Post-Iterations <==\n\n")

    if self.V_potential():
        quad = self.V_potential().quadrature_values()
        rho_a = quad['RHO_A']/2 if self.same_a_b_dens() else quad['RHO_A']
        rho_b = quad['RHO_B']/2 if self.same_a_b_dens() else quad['RHO_B']
        rho_ab = (rho_a + rho_b)
        self.set_variable("GRID ELECTRONS TOTAL",rho_ab)  # P::e SCF
        self.set_variable("GRID ELECTRONS ALPHA",rho_a)  # P::e SCF
        self.set_variable("GRID ELECTRONS BETA",rho_b)  # P::e SCF
        dev_a = rho_a - self.nalpha()
        dev_b = rho_b - self.nbeta()
        core.print_out(f"   Electrons on quadrature grid:\n")
        if self.same_a_b_dens():
            core.print_out(f"      Ntotal   = {rho_ab:15.10f} ; deviation = {dev_b+dev_a:.3e} \n\n")
        else:
            core.print_out(f"      Nalpha   = {rho_a:15.10f} ; deviation = {dev_a:.3e}\n")
            core.print_out(f"      Nbeta    = {rho_b:15.10f} ; deviation = {dev_b:.3e}\n")
            core.print_out(f"      Ntotal   = {rho_ab:15.10f} ; deviation = {dev_b+dev_a:.3e} \n\n")
        if ((dev_b+dev_a) > 0.1):
            core.print_out("   WARNING: large deviation in the electron count on grid detected. Check grid size!")
    self.check_phases()
    self.compute_spin_contamination()
    self.frac_renormalize()
    reference = core.get_option("SCF", "REFERENCE")

    energy = self.get_energies("Total Energy")

    #    fail_on_maxiter = core.get_option("SCF", "FAIL_ON_MAXITER")
    #    if converged or not fail_on_maxiter:
    #
    #        if print_lvl > 0:
    #            self.print_orbitals()
    #
    #        if converged:
    #            core.print_out("  Energy converged.\n\n")
    #        else:
    #            core.print_out("  Energy did not converge, but proceeding anyway.\n\n")

    if core.get_option('SCF', 'PRINT') > 0:
        self.print_orbitals()

    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    core.print_out("  @%s%s Final Energy: %20.14f" % ('DF-' if is_dfjk else '', reference, energy))
    # if (perturb_h_) {
    #     core.print_out(" with %f %f %f perturbation" %
    #                    (dipole_field_strength_[0], dipole_field_strength_[1], dipole_field_strength_[2]))
    # }
    core.print_out("\n\n")
    self.print_energies()

    # force list into Matrix for storage
    iteration_energies = np.array(self.iteration_energies).reshape(-1, 1)
    iteration_energies = core.Matrix.from_array(iteration_energies)
    core.set_variable("SCF TOTAL ENERGIES", core.Matrix.from_array(iteration_energies))
    self.set_variable("SCF TOTAL ENERGIES", core.Matrix.from_array(iteration_energies))

    self.clear_external_potentials()
    if core.get_option('SCF', 'PCM'):
        calc_type = core.PCM.CalcType.Total
        if core.get_option("PCM", "PCM_SCF_TYPE") == "SEPARATE":
            calc_type = core.PCM.CalcType.NucAndEle
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpcm = self.get_PCM().compute_PCM_terms(Dt, calc_type)
        self.push_back_external_potential(Vpcm)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation("PCM", lambda pert_dm : self.get_PCM().compute_V(pert_dm))

    if core.get_option('SCF', 'PE'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        _, Vpe = self.pe_state.get_pe_contribution(
            Dt, elec_only=False
        )
        self.push_back_external_potential(Vpe)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation("PE", lambda pert_dm : self.pe_state.get_pe_contribution(pert_dm, elec_only=True)[1])

    if core.get_option('SCF', 'DDX'):
        Dt = self.Da().clone()
        Dt.add(self.Db())
        Vddx = self.ddx.get_solvation_contributions(Dt)[1]
        self.push_back_external_potential(Vddx)
        # Set callback function for CPSCF
        self.set_external_cpscf_perturbation(
            "DDX", lambda pert_dm : self.ddx.get_solvation_contributions(pert_dm, elec_only=True, nonequilibrium=True)[1])

    # Orbitals are always saved, in case an MO guess is requested later
    # save_orbitals()

    # Shove variables into global space
    for k, v in self.variables().items():
        core.set_variable(k, v)

    # TODO re-enable
    self.finalize()
    if self.V_potential():
        self.V_potential().clear_collocation_cache()

    core.print_out("\nComputation Completed\n")
    core.del_variable("SCF D NORM")

    return energy


def scf_print_energies(self):
    enuc = self.get_energies('Nuclear')
    e1 = self.get_energies('One-Electron')
    e2 = self.get_energies('Two-Electron')
    exc = self.get_energies('XC')
    ed = self.get_energies('-D')
    self.del_variable('-D Energy')
    evv10 = self.get_energies('VV10')
    eefp = self.get_energies('EFP')
    epcm = self.get_energies('PCM Polarization')
    edd = self.get_energies('DD Solvation Energy')
    epe = self.get_energies('PE Energy')
    ke = self.get_energies('Kinetic')

    hf_energy = enuc + e1 + e2
    dft_energy = hf_energy + exc + ed + evv10
    total_energy = dft_energy + eefp + epcm + edd + epe
    full_qm = (not core.get_option('SCF', 'PCM') and not core.get_option('SCF', 'DDX') and not core.get_option('SCF', 'PE')
               and not hasattr(self.molecule(), 'EFP'))

    core.print_out("   => Energetics <=\n\n")
    core.print_out("    Nuclear Repulsion Energy =        {:24.16f}\n".format(enuc))
    core.print_out("    One-Electron Energy =             {:24.16f}\n".format(e1))
    core.print_out("    Two-Electron Energy =             {:24.16f}\n".format(e2))
    if self.functional().needs_xc():
        core.print_out("    DFT Exchange-Correlation Energy = {:24.16f}\n".format(exc))
        core.print_out("    Empirical Dispersion Energy =     {:24.16f}\n".format(ed))
        core.print_out("    VV10 Nonlocal Energy =            {:24.16f}\n".format(evv10))
    if core.get_option('SCF', 'PCM'):
        core.print_out("    PCM Polarization Energy =         {:24.16f}\n".format(epcm))
    if core.get_option('SCF', 'DDX'):
        core.print_out("    DD Solvation Energy =            {:24.16f}\n".format(edd))
    if core.get_option('SCF', 'PE'):
        core.print_out("    PE Energy =                       {:24.16f}\n".format(epe))
    if hasattr(self.molecule(), 'EFP'):
        core.print_out("    EFP Energy =                      {:24.16f}\n".format(eefp))
    core.print_out("    Total Energy =                    {:24.16f}\n".format(total_energy))

    if core.get_option('SCF', 'PE'):
        core.print_out(self.pe_state.cppe_state.summary_string)

    self.set_variable("NUCLEAR REPULSION ENERGY", enuc)  # P::e SCF
    self.set_variable("ONE-ELECTRON ENERGY", e1)  # P::e SCF
    self.set_variable("TWO-ELECTRON ENERGY", e2)  # P::e SCF
    if self.functional().needs_xc():
        self.set_variable("DFT XC ENERGY", exc)  # P::e SCF
        self.set_variable("DFT VV10 ENERGY", evv10)  # P::e SCF
        self.set_variable("DFT FUNCTIONAL TOTAL ENERGY", hf_energy + exc + evv10)  # P::e SCF
        #self.set_variable(self.functional().name() + ' FUNCTIONAL TOTAL ENERGY', hf_energy + exc + evv10)
        self.set_variable("DFT TOTAL ENERGY", dft_energy)  # overwritten later for DH  # P::e SCF
    else:
        potential = total_energy - ke
        self.set_variable("HF KINETIC ENERGY", ke)  # P::e SCF
        self.set_variable("HF POTENTIAL ENERGY", potential)  # P::e SCF
        if full_qm:
            self.set_variable("HF VIRIAL RATIO", - potential / ke)  # P::e SCF
        self.set_variable("HF TOTAL ENERGY", hf_energy)  # P::e SCF
    if hasattr(self, "_disp_functor"):
        self.set_variable("DISPERSION CORRECTION ENERGY", ed)  # P::e SCF
    #if abs(ed) > 1.0e-14:
    #    for pv, pvv in self.variables().items():
    #        if abs(pvv - ed) < 1.0e-14:
    #            if pv.endswith('DISPERSION CORRECTION ENERGY') and pv.startswith(self.functional().name()):
    #                fctl_plus_disp_name = pv.split()[0]
    #                self.set_variable(fctl_plus_disp_name + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH
    #else:
    #    self.set_variable(self.functional().name() + ' TOTAL ENERGY', dft_energy)  # overwritten later for DH

    self.set_variable("SCF ITERATIONS", self.iteration_)  # P::e SCF


def scf_print_preiterations(self,small=False):
    # small version does not print Nalpha,Nbeta,Ndocc,Nsocc, e.g. for SAD guess where they are not
    # available
    ct = self.molecule().point_group().char_table()

    if not small:
        core.print_out("   -------------------------------------------------------\n")
        core.print_out("    Irrep   Nso     Nmo     Nalpha   Nbeta   Ndocc  Nsocc\n")
        core.print_out("   -------------------------------------------------------\n")

        for h in range(self.nirrep()):
            core.print_out(
                f"     {ct.gamma(h).symbol():<3s}   {self.nsopi()[h]:6d}  {self.nmopi()[h]:6d}  {self.nalphapi()[h]:6d}  {self.nbetapi()[h]:6d}  {self.doccpi()[h]:6d}  {self.soccpi()[h]:6d}\n"
            )

        core.print_out("   -------------------------------------------------------\n")
        core.print_out(
            f"    Total  {self.nso():6d}  {self.nmo():6d}  {self.nalpha():6d}  {self.nbeta():6d}  {self.nbeta():6d}  {self.nalpha() - self.nbeta():6d}\n"
        )
        core.print_out("   -------------------------------------------------------\n\n")
    else:
        core.print_out("   -------------------------\n")
        core.print_out("    Irrep   Nso     Nmo    \n")
        core.print_out("   -------------------------\n")

        for h in range(self.nirrep()):
            core.print_out(
                f"     {ct.gamma(h).symbol():<3s}   {self.nsopi()[h]:6d}  {self.nmopi()[h]:6d} \n"
            )

        core.print_out("   -------------------------\n")
        core.print_out(
            f"    Total  {self.nso():6d}  {self.nmo():6d}\n"
        )
        core.print_out("   -------------------------\n\n")


# Bind functions to core.HF class
core.HF.initialize = scf_initialize
core.HF.initialize_jk = initialize_jk
core.HF.iterations = scf_iterate
core.HF._scf_initialize_iteration_state = _scf_initialize_iteration_state
core.HF._scf_iteration = _scf_iteration
core.HF.compute_energy = scf_compute_energy
core.HF.finalize_energy = scf_finalize_energy
core.HF.print_energies = scf_print_energies
core.HF.print_preiterations = scf_print_preiterations
core.HF.iteration_energies = []


def _converged(e_delta, d_rms, e_conv=None, d_conv=None):
    if e_conv is None:
        e_conv = core.get_option("SCF", "E_CONVERGENCE")
    if d_conv is None:
        d_conv = core.get_option("SCF", "D_CONVERGENCE")

    return (abs(e_delta) < e_conv and d_rms < d_conv)


def _validate_damping(wfn=None):
    """Sanity-checks DAMPING control options

    Parameters
    ----------
    wfn : HF wavefunction, optional
        If provided, reads options from wfn._options_snapshot instead of global.
        This ensures deterministic behavior in multi_scf().

    Raises
    ------
    ValidationError
        If any of |scf__damping_percentage|, |scf__damping_convergence|
        don't play well together.

    Returns
    -------
    bool
        Whether DAMPING is enabled during scf.

    """
    # Read from snapshot if available, otherwise from global
    if wfn is not None:
        damping_pct = get_option_from_snapshot(wfn, 'DAMPING_PERCENTAGE')
        damping_conv = get_option_from_snapshot(wfn, 'DAMPING_CONVERGENCE')
    else:
        damping_pct = core.get_option('SCF', 'DAMPING_PERCENTAGE')
        damping_conv = core.get_option('SCF', 'DAMPING_CONVERGENCE')

    enabled = (damping_pct > 0.0)
    if enabled:
        if damping_pct < 0.0 or damping_pct > 100.0:
            raise ValidationError('SCF DAMPING_PERCENTAGE ({}) must be between 0 and 100'.format(damping_pct))

        if damping_conv < 0.0:
            raise ValidationError('SCF DAMPING_CONVERGENCE ({}) must be > 0'.format(damping_conv))

    return enabled


def _validate_diis(self):
    """Sanity-checks DIIS control options

    Reads options from wfn._options_snapshot if available, otherwise falls back to global.

    Raises
    ------
    psi4.driver.p4util.exceptions.ValidationError
        If any of DIIS options don't play well together.

    Returns
    -------
    bool
        Whether some form of DIIS is enabled during SCF.

    """

    restricted_open = self.same_a_b_orbs() and not self.same_a_b_dens()
    aediis_accelerator = get_option_from_snapshot(self, 'SCF_INITIAL_ACCELERATOR')
    aediis_active = aediis_accelerator != "NONE" and not restricted_open

    if aediis_active:
        start = get_option_from_snapshot(self, 'SCF_INITIAL_START_DIIS_TRANSITION')
        stop = get_option_from_snapshot(self, 'SCF_INITIAL_FINISH_DIIS_TRANSITION')
        if start < stop:
            raise ValidationError('SCF_INITIAL_START_DIIS_TRANSITION error magnitude cannot be less than SCF_INITIAL_FINISH_DIIS_TRANSITION.')
        elif start < 0:
            raise ValidationError('SCF_INITIAL_START_DIIS_TRANSITION cannot be negative.')
        elif stop < 0:
            raise ValidationError('SCF_INITIAL_FINISH_DIIS_TRANSITION cannot be negative.')

    diis_enabled = get_option_from_snapshot(self, 'DIIS')
    enabled = bool(diis_enabled) or aediis_active
    if enabled:
        start = get_option_from_snapshot(self, 'DIIS_START')
        if start < 1:
            raise ValidationError('SCF DIIS_START ({}) must be at least 1'.format(start))

    return enabled


def _validate_frac(wfn=None):
    """Sanity-checks FRAC control options

    Parameters
    ----------
    wfn : HF wavefunction, optional
        If provided, reads options from wfn._options_snapshot instead of global.

    Raises
    ------
    ValidationError
        If any of |scf__frac_start| don't play well together.

    Returns
    -------
    bool
        Whether FRAC is enabled during scf.

    """
    if wfn is not None:
        frac_start = get_option_from_snapshot(wfn, 'FRAC_START')
    else:
        frac_start = core.get_option('SCF', 'FRAC_START')

    enabled = (frac_start != 0)
    if enabled:
        if frac_start < 0:
            raise ValidationError('SCF FRAC_START ({}) must be at least 1'.format(frac_start))

    return enabled


def _validate_MOM(wfn=None):
    """Sanity-checks MOM control options

    Parameters
    ----------
    wfn : HF wavefunction, optional
        If provided, reads options from wfn._options_snapshot instead of global.

    Raises
    ------
    ValidationError
        If any of |scf__mom_start|, |scf__mom_occ| don't play well together.

    Returns
    -------
    bool
        Whether excited-state MOM (not just the plain stabilizing MOM) is enabled during scf.

    """
    if wfn is not None:
        mom_start = get_option_from_snapshot(wfn, 'MOM_START')
        mom_occ = get_option_from_snapshot(wfn, 'MOM_OCC')
    else:
        mom_start = core.get_option('SCF', "MOM_START")
        mom_occ = core.get_option('SCF', "MOM_OCC")

    enabled = (mom_start != 0 and len(mom_occ) > 0)
    if enabled:
        if mom_start < 0:
            raise ValidationError('SCF MOM_START ({}) must be at least 1'.format(mom_start))

    return enabled


def _validate_soscf(wfn=None):
    """Sanity-checks SOSCF control options

    Parameters
    ----------
    wfn : HF wavefunction, optional
        If provided, reads options from wfn._options_snapshot instead of global.

    Raises
    ------
    ValidationError
        If any of |scf__soscf|, |scf__soscf_start_convergence|,
        |scf__soscf_min_iter|, |scf__soscf_max_iter| don't play well together.

    Returns
    -------
    bool
        Whether SOSCF is enabled during scf.

    """
    # Read from snapshot if available
    if wfn is not None:
        enabled = get_option_from_snapshot(wfn, 'SOSCF')
        start = get_option_from_snapshot(wfn, 'SOSCF_START_CONVERGENCE')
        miniter = get_option_from_snapshot(wfn, 'SOSCF_MIN_ITER')
    else:
        enabled = core.get_option('SCF', 'SOSCF')
        start = core.get_option('SCF', 'SOSCF_START_CONVERGENCE')
        miniter = core.get_option('SCF', 'SOSCF_MIN_ITER')

    if enabled:
        if start < 0.0:
            raise ValidationError('SCF SOSCF_START_CONVERGENCE ({}) must be positive'.format(start))

        if miniter < 1:
            raise ValidationError('SCF SOSCF_MIN_ITER ({}) must be at least 1'.format(miniter))

        if wfn is not None:
            maxiter = get_option_from_snapshot(wfn, 'SOSCF_MAX_ITER')
            conv = get_option_from_snapshot(wfn, 'SOSCF_CONV')
        else:
            maxiter = core.get_option('SCF', 'SOSCF_MAX_ITER')
            conv = core.get_option('SCF', 'SOSCF_CONV')

        if maxiter < miniter:
            raise ValidationError('SCF SOSCF_MAX_ITER ({}) must be at least SOSCF_MIN_ITER ({})'.format(
                maxiter, miniter))

        if conv < 1.e-10:
            raise ValidationError('SCF SOSCF_CONV ({}) must be achievable'.format(conv))

    return enabled

core.HF.validate_diis = _validate_diis


def validate_multi_scf_compatibility(wfn_list):
    """
    Validate that wavefunctions can share a single JK object.

    For shared JK computation, the following must match across all wfn:
    - Primary basis set (JK built with specific basis dimensions)
    - Geometry (3-index integrals depend on atomic coordinates)
    - Auxiliary basis (if SCF_TYPE='DF')
    - LRC capability (do_wK flag)
    - LRC omega parameter (if LRC functional)
    - RSH alpha/beta parameters (if LRC and WCOMBINE=TRUE)

    Safe to differ: multiplicity, reference type, charge, convergence settings.

    Parameters
    ----------
    wfn_list : list of HF
        Wavefunctions to validate for shared JK compatibility

    Raises
    ------
    ValidationError
        If wavefunctions have incompatible parameters for shared JK
    """
    if len(wfn_list) < 2:
        return  # Single wfn, no compatibility issues

    ref_wfn = wfn_list[0]

    # ========================================================================
    # Category 1: JK Structure (baked into JK.build(), cannot change)
    # ========================================================================

    # Check 1.1: Primary Basis Set
    # Code: _build_jk() line 97: wfn.get_basisset("ORBITAL")
    # Why: JK object built with specific basis dimensions, 3-index integrals
    #      computed for this basis. Different basis requires different integrals.
    ref_basis = ref_wfn.basisset()
    ref_basis_name = ref_basis.name()
    ref_nbf = ref_basis.nbf()

    for i, wfn in enumerate(wfn_list[1:], 1):
        wfn_basis = wfn.basisset()
        wfn_basis_name = wfn_basis.name()
        wfn_nbf = wfn_basis.nbf()

        if wfn_basis_name != ref_basis_name or wfn_nbf != ref_nbf:
            raise ValidationError(
                f"\n"
                f"Shared JK Compatibility Error: Primary basis set mismatch\n"
                f"{'=' * 70}\n"
                f"Wavefunction {i} has different primary basis set:\n"
                f"  Reference (wfn 0): {ref_basis_name} ({ref_nbf} basis functions)\n"
                f"  Wavefunction {i}:  {wfn_basis_name} ({wfn_nbf} basis functions)\n"
                f"\n"
                f"Why this matters:\n"
                f"  All wavefunctions share a single JK object built with specific\n"
                f"  basis set dimensions. The 3-index integrals (Q|μν) are computed\n"
                f"  for this basis. Different basis sets require different integrals.\n"
                f"\n"
                f"Code location: _build_jk() line 97\n"
                f"  jk = core.JK.build(wfn.get_basisset('ORBITAL'), ...)\n"
                f"\n"
                f"Solution:\n"
                f"  Use the same primary basis set for all wavefunctions in multi_scf().\n"
                f"  Example: set basis cc-pVDZ for all molecules before creating wfn.\n"
                f"{'=' * 70}\n"
            )

    # Check 1.2: Geometry (Molecule)
    # Code: jk.initialize() line 121 - computes 3-index integrals
    # Why: Integrals (Q|μν) = ∫∫ χ_Q(r1) χ_μ(r1) r12^-1 χ_ν(r2) dr1 dr2
    #      Basis functions χ centered on atoms, integrals depend on geometry
    ref_mol = ref_wfn.molecule()

    for i, wfn in enumerate(wfn_list[1:], 1):
        wfn_mol = wfn.molecule()

        # Check if same molecule object (identity check is fastest)
        if wfn_mol is not ref_mol:
            # Different objects - check if geometries are identical
            # This can happen if molecules were cloned
            ref_geom = ref_mol.geometry().np.flatten()
            wfn_geom = wfn_mol.geometry().np.flatten()

            if ref_geom.shape != wfn_geom.shape or \
               not core.np.allclose(ref_geom, wfn_geom, atol=1e-10):
                raise ValidationError(
                    f"\n"
                    f"Shared JK Compatibility Error: Geometry mismatch\n"
                    f"{'=' * 70}\n"
                    f"Wavefunction {i} has different molecular geometry:\n"
                    f"  Reference and wfn {i} have different atomic coordinates.\n"
                    f"\n"
                    f"Why this matters:\n"
                    f"  The 3-index integrals (Q|μν) depend on atomic positions because\n"
                    f"  basis functions are centered on atoms. Different geometries\n"
                    f"  require completely different integrals.\n"
                    f"\n"
                    f"Code location: jk.initialize() called in initialize_jk() line 121\n"
                    f"  Computes integrals for current molecular geometry.\n"
                    f"\n"
                    f"Solution:\n"
                    f"  Use the same molecule object for all wavefunctions.\n"
                    f"  If you need different multiplicities/charges, use:\n"
                    f"    mol.set_multiplicity(...)  # Changes occupation\n"
                    f"    mol.set_molecular_charge(...)  # Changes charge\n"
                    f"  but keep the SAME molecule object (same geometry).\n"
                    f"{'=' * 70}\n"
                )

    # Check 1.3: Auxiliary Basis (for DF only)
    # Code: _build_jk() line 98: aux=wfn.get_basisset("DF_BASIS_SCF")
    # Why: For density fitting, different auxiliary basis yields different approximation
    #      (μν|ρσ) ≈ (μν|P) (P|Q)^-1 (Q|ρσ)
    scf_type = core.get_global_option('SCF_TYPE')
    if scf_type == 'DF':
        ref_aux = ref_wfn.get_basisset("DF_BASIS_SCF")
        ref_aux_name = ref_aux.name()

        for i, wfn in enumerate(wfn_list[1:], 1):
            wfn_aux = wfn.get_basisset("DF_BASIS_SCF")
            wfn_aux_name = wfn_aux.name()

            if wfn_aux_name != ref_aux_name:
                raise ValidationError(
                    f"\n"
                    f"Shared JK Compatibility Error: DF auxiliary basis mismatch\n"
                    f"{'=' * 70}\n"
                    f"Wavefunction {i} has different DF auxiliary basis:\n"
                    f"  Reference (wfn 0): {ref_aux_name}\n"
                    f"  Wavefunction {i}:  {wfn_aux_name}\n"
                    f"\n"
                    f"Why this matters:\n"
                    f"  For SCF_TYPE='DF', the auxiliary basis is used to approximate\n"
                    f"  4-index integrals: (μν|ρσ) ≈ (μν|P) (P|Q)^-1 (Q|ρσ)\n"
                    f"  Different auxiliary basis yields different J/K integrals.\n"
                    f"\n"
                    f"Code location: _build_jk() line 98\n"
                    f"  jk = core.JK.build(..., aux=wfn.get_basisset('DF_BASIS_SCF'))\n"
                    f"\n"
                    f"Solution:\n"
                    f"  Use the same DF auxiliary basis for all wavefunctions.\n"
                    f"  Example: set df_basis_scf cc-pVDZ-RI for all molecules.\n"
                    f"{'=' * 70}\n"
                )

    # Check 1.4: LRC Capability
    # Code: _build_jk() line 99: do_wK=wfn.functional().is_x_lrc()
    # Why: JK built with or without long-range K capability. If built without,
    #      cannot compute wK later!
    ref_is_lrc = ref_wfn.functional().is_x_lrc()

    for i, wfn in enumerate(wfn_list[1:], 1):
        wfn_is_lrc = wfn.functional().is_x_lrc()

        if wfn_is_lrc != ref_is_lrc:
            ref_func_name = ref_wfn.functional().name()
            wfn_func_name = wfn.functional().name()

            raise ValidationError(
                f"\n"
                f"Shared JK Compatibility Error: LRC capability mismatch\n"
                f"{'=' * 70}\n"
                f"Wavefunction {i} has incompatible functional:\n"
                f"  Reference (wfn 0): {ref_func_name} (is_x_lrc={ref_is_lrc})\n"
                f"  Wavefunction {i}:  {wfn_func_name} (is_x_lrc={wfn_is_lrc})\n"
                f"\n"
                f"Why this matters:\n"
                f"  The shared JK object is built with or without long-range K (wK)\n"
                f"  capability. If built without wK support (non-LRC functional), it\n"
                f"  cannot compute wK for LRC functionals later!\n"
                f"\n"
                f"Code location: _build_jk() line 99\n"
                f"  jk = core.JK.build(..., do_wK=wfn.functional().is_x_lrc())\n"
                f"\n"
                f"Solution:\n"
                f"  All wavefunctions must be either:\n"
                f"    - All LRC functionals (ωB97X, ωB97X-D, LC-ωPBE, etc.), OR\n"
                f"    - All non-LRC (HF, B3LYP, PBE0, etc.)\n"
                f"  Mixing LRC and non-LRC is not supported for shared JK.\n"
                f"{'=' * 70}\n"
            )

    # ========================================================================
    # Category 2: JK Configuration (set in initialize_jk(), not updated later)
    # ========================================================================

    # Check 2.1: LRC omega parameter
    # Code: initialize_jk() line 116
    # Why: omega is the range-separation parameter used in erf(ωr)/r operator
    #      for long-range wK integrals. Different omega yields different integrals.
    if ref_is_lrc:
        ref_omega = ref_wfn.functional().x_omega()

        for i, wfn in enumerate(wfn_list[1:], 1):
            wfn_omega = wfn.functional().x_omega()

            # Check omega MUST match
            if abs(wfn_omega - ref_omega) > 1e-10:
                ref_func_name = ref_wfn.functional().name()
                wfn_func_name = wfn.functional().name()

                raise ValidationError(
                    f"\n"
                    f"Shared JK Compatibility Error: LRC omega parameter mismatch\n"
                    f"{'=' * 70}\n"
                    f"Wavefunction {i} has different omega parameter:\n"
                    f"  Reference (wfn 0): {ref_func_name}, omega = {ref_omega}\n"
                    f"  Wavefunction {i}:  {wfn_func_name}, omega = {wfn_omega}\n"
                    f"\n"
                    f"Why this matters:\n"
                    f"  For LRC functionals, omega is the range-separation parameter:\n"
                    f"    wK_μν = Σ_ρσ P_ρσ (μρ|erf(ω r)/r|νσ)\n"
                    f"  Different omega yields different wK integrals.\n"
                    f"  The shared JK is configured with omega from wfn[0] only.\n"
                    f"  Other wavefunctions skip JK configuration due to idempotency,\n"
                    f"  so they inherit omega from wfn[0].\n"
                    f"\n"
                    f"Code location: initialize_jk() line 116\n"
                    f"  jk.set_omega(functional.x_omega())\n"
                    f"  Only called for wfn[0], skipped for wfn[1:] (scf_initialize\n"
                    f"  idempotency check at lines 146-148).\n"
                    f"\n"
                    f"Solution:\n"
                    f"  All wavefunctions must use LRC functionals with the SAME omega.\n"
                    f"  Different LRC functionals (ωB97X vs ωB97X-D) typically have\n"
                    f"  different omega values, so they cannot be mixed.\n"
                    f"{'=' * 70}\n"
                )

    # Checks 2.2-2.3: LRC alpha/beta parameters (CONDITIONAL on wcombine)
    # Code: initialize_jk() lines 118-119
    # Why: alpha/beta affect integrals ONLY if wcombine=TRUE (combined K mode)
    #
    # Physics:
    #   wcombine=FALSE (default): K and wK computed separately, alpha/beta used
    #                             during Fock assembly, can differ between wfn
    #   wcombine=TRUE: Combined matrix = alpha*K + beta*wK computed directly
    #                  with alpha/beta baked into integrals, must match
    #
    # See: INVESTIGATION_VALIDATION_CHECKS.md for C++ code analysis (DFHelper.cc)

    # Check if wcombine mode is enabled (global option)
    wcombine_enabled = core.get_option('SCF', 'WCOMBINE') if core.has_option_changed('SCF', 'WCOMBINE') else False

    if ref_is_lrc and wcombine_enabled:
        # wcombine mode: alpha/beta affect integrals, MUST validate
        ref_alpha = ref_wfn.functional().x_alpha()
        ref_beta = ref_wfn.functional().x_beta()

        for i, wfn in enumerate(wfn_list[1:], 1):
            wfn_alpha = wfn.functional().x_alpha()
            wfn_beta = wfn.functional().x_beta()

            # Check alpha
            if abs(wfn_alpha - ref_alpha) > 1e-10:
                ref_func_name = ref_wfn.functional().name()
                wfn_func_name = wfn.functional().name()

                raise ValidationError(
                    f"\n"
                    f"Shared JK Compatibility Error: RSH alpha parameter mismatch (wcombine=TRUE)\n"
                    f"{'=' * 70}\n"
                    f"Wavefunction {i} has different alpha parameter:\n"
                    f"  Reference (wfn 0): {ref_func_name}, alpha = {ref_alpha}\n"
                    f"  Wavefunction {i}:  {wfn_func_name}, alpha = {wfn_alpha}\n"
                    f"\n"
                    f"Why this matters:\n"
                    f"  You have enabled WCOMBINE=TRUE (combined K matrix mode).\n"
                    f"  In this mode, the combined matrix K_combined = alpha*K + beta*wK\n"
                    f"  is computed directly with alpha/beta baked into integrals.\n"
                    f"  Different alpha yields incorrect results.\n"
                    f"\n"
                    f"Code location: DFHelper.cc compute_sparse_pQq_blocking_p_symm_abw()\n"
                    f"  param_Mp = omega_alpha * buffer + omega_beta * wbuffer\n"
                    f"  (only called if wcombine=TRUE)\n"
                    f"\n"
                    f"Solution:\n"
                    f"  All wavefunctions must use functionals with the SAME alpha.\n"
                    f"  Or disable wcombine mode (set WCOMBINE=FALSE, which is default).\n"
                    f"{'=' * 70}\n"
                )

            # Check beta
            if abs(wfn_beta - ref_beta) > 1e-10:
                ref_func_name = ref_wfn.functional().name()
                wfn_func_name = wfn.functional().name()

                raise ValidationError(
                    f"\n"
                    f"Shared JK Compatibility Error: RSH beta parameter mismatch (wcombine=TRUE)\n"
                    f"{'=' * 70}\n"
                    f"Wavefunction {i} has different beta parameter:\n"
                    f"  Reference (wfn 0): {ref_func_name}, beta = {ref_beta}\n"
                    f"  Wavefunction {i}:  {wfn_func_name}, beta = {wfn_beta}\n"
                    f"\n"
                    f"Why this matters:\n"
                    f"  You have enabled WCOMBINE=TRUE (combined K matrix mode).\n"
                    f"  In this mode, the combined matrix K_combined = alpha*K + beta*wK\n"
                    f"  is computed directly with alpha/beta baked into integrals.\n"
                    f"  Different beta yields incorrect results.\n"
                    f"\n"
                    f"Code location: DFHelper.cc compute_sparse_pQq_blocking_p_symm_abw()\n"
                    f"  param_Mp = omega_alpha * buffer + omega_beta * wbuffer\n"
                    f"  (only called if wcombine=TRUE)\n"
                    f"\n"
                    f"Solution:\n"
                    f"  All wavefunctions must use functionals with the SAME beta.\n"
                    f"  Or disable wcombine mode (set WCOMBINE=FALSE, which is default).\n"
                    f"{'=' * 70}\n"
                )
    # else: wcombine=FALSE (default) - alpha/beta used in Fock assembly, not in integrals, OK to differ

    # All checks passed!
    # Parameters that CAN differ (safe) are NOT checked:
    # - Multiplicity, Reference type, Charge, Non-LRC XC, Hybrid fraction,
    #   Convergence settings - these don't affect shared JK
    # Parameters already protected:
    # - SCF_TYPE (options snapshot), do_J/do_K (overwritten in _multi_scf_inner)


def multi_scf(wfn_list, e_conv=None, d_conv=None, max_iter=None, verbose=True):
    """
    Run multiple SCF calculations with shared JK computation.

    This multi-cycle SCF coordinator uses the refactored _scf_iteration()
    method from each wavefunction. This version supports ALL SCF features:
    DIIS, damping, SOSCF, MOM, FRAC, convergence acceleration, etc.

    Algorithm
    ---------
    1. Initialize iteration state for each wfn via _scf_initialize_iteration_state()
    2. Main loop:
       a. Collect all C matrices from all wfn (via get_orbital_matrices())
       b. Single shared JK computation (jk.compute())
       c. Distribute J/K to each wfn (via set_jk_matrices())
       d. Each wfn._scf_iteration() uses precomputed J/K
       e. Check convergence for all wfn
    3. Return final energies

    Parameters
    ----------
    wfn_list : list of HF wavefunction objects
        List of RHF, UHF, or ROHF wavefunctions to converge simultaneously
    e_conv : float, optional
        Energy convergence threshold. Defaults to SCF E_CONVERGENCE option
    d_conv : float, optional
        Density convergence threshold. Defaults to SCF D_CONVERGENCE option
    max_iter : int, optional
        Maximum number of iterations. Defaults to SCF MAXITER option
    verbose : bool, optional
        Print iteration information. Default True

    Returns
    -------
    list of float
        Final energies for each wavefunction

    Raises
    ------
    ValidationError
        If wfn_list is empty or wavefunctions are incompatible
    SCFConvergenceError
        If any wavefunction fails to converge within max_iter

    Examples
    --------
    Converge singlet and triplet states simultaneously::

        import psi4
        mol = psi4.geometry('''
        0 1
        H 0 0 0
        H 0 0 0.74
        ''')
        psi4.set_options({'basis': 'cc-pVDZ'})

        # Create singlet wavefunction
        mol.set_multiplicity(1)
        wfn_s = psi4.core.RHF(psi4.core.Wavefunction.build(mol, 'cc-pVDZ'))

        # Create triplet wavefunction
        mol.set_multiplicity(3)
        wfn_t = psi4.core.UHF(psi4.core.Wavefunction.build(mol, 'cc-pVDZ'))

        # Converge both simultaneously with shared JK
        energies = multi_scf([wfn_s, wfn_t])

        print(f"Singlet energy: {energies[0]}")
        print(f"Triplet energy: {energies[1]}")

    Notes
    -----
    - All wavefunctions must use the same basis set
    - All wavefunctions must use the same JK algorithm (SCF_TYPE)
    - Each wfn maintains its own DIIS, damping, convergence state
    - ~1.8-2x faster than running SCF cycles independently
    - Uses existing C++ infrastructure (use_precomputed_jk_ flag)

    DF_SCF_GUESS Support (Andy Trick 2.0)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    When DF_SCF_GUESS=True (default) and SCF_TYPE='DIRECT':
    - Phase 1: All wfn converge using fast DF iterations
    - Phase 2: DIIS reset, JK reinitialized for DIRECT
    - Phase 3: Final convergence with DIRECT integrals

    This matches single SCF behavior and produces identical iteration counts.
    Total iterations shown = DF iterations + DIRECT iterations.
    """
    if len(wfn_list) == 0:
        raise ValidationError("multi_scf requires at least one wavefunction")

    # Validate wavefunction compatibility for shared JK
    # Checks: basis set, geometry, LRC parameters, etc.
    # Only checks what affects SHARED components (no half-measures!)
    validate_multi_scf_compatibility(wfn_list)

    # Auto-assign unique wavefunction names for file isolation
    # This enables proper DIIS, stability analysis, and orbital file separation
    wfn_names_used = set()
    for i, wfn in enumerate(wfn_list):
        current_name = wfn.get_wfn_name()

        # Auto-assign if empty (default for single-cycle SCF)
        if not current_name:
            wfn_name = f"wfn_{i}"
            wfn.set_wfn_name(wfn_name)
        else:
            wfn_name = current_name

        # Validate no duplicates
        if wfn_name in wfn_names_used:
            raise ValidationError(
                f"Duplicate wavefunction name '{wfn_name}' detected. "
                f"Each wavefunction in multi_scf must have a unique name for file isolation. "
                f"Either clear all names (for auto-naming) or ensure all custom names are unique."
            )
        wfn_names_used.add(wfn_name)

    # Get convergence criteria
    if e_conv is None:
        e_conv = core.get_option('SCF', 'E_CONVERGENCE')
    if d_conv is None:
        d_conv = core.get_option('SCF', 'D_CONVERGENCE')
    if max_iter is None:
        max_iter = core.get_option('SCF', 'MAXITER')

    # Snapshot global options once before creating/initializing any wfn
    # Prevents non-determinism from global state pollution between wfn creation
    options_snapshot = snapshot_scf_options()

    # Apply options snapshot to all wfn before any initialization
    # Ensures wfn.initialize() reads from frozen snapshot, not from global state
    for wfn in wfn_list:
        apply_options_snapshot(wfn, options_snapshot)

    # Andy Trick 2.0 for multi-SCF: DF guess for DIRECT algorithm
    # This speeds up DIRECT (which recomputes full integrals each iteration)
    # by first converging via fast DF iterations, then fully converging in fewer DIRECT iterations
    use_df_guess = (core.get_option('SCF', 'DF_SCF_GUESS') and
                    core.get_global_option('SCF_TYPE') == 'DIRECT')

    if use_df_guess:
        if verbose:
            core.print_out("  Starting with a DF guess for DIRECT algorithm...\n\n")

        # Phase 1: DF pre-iterations (fast convergence)
        with p4util.OptionsStateCM(['SCF_TYPE']):
            core.set_global_option('SCF_TYPE', 'DF')

            # Initialize all wfn with DF
            for wfn in wfn_list:
                if wfn.jk() is None:
                    wfn.initialize()

            # Run multi-SCF iterations with DF
            try:
                _multi_scf_inner(wfn_list, e_conv, d_conv, max_iter, verbose)
            except SCFConvergenceError:
                # DF guess failed to converge - not critical, continue to DIRECT
                pass

        if verbose:
            core.print_out("\n  DF guess converged. Switching to DIRECT...\n\n")

        # Phase 2: Reset for DIRECT
        for wfn in wfn_list:
            # Reset DIIS subspace (DF gradients not compatible with DIRECT)
            if wfn.initialized_diis_manager_:
                wfn.diis_manager_.reset_subspace()

            # Re-initialize JK with DIRECT (recomputes full integrals)
            wfn.initialize_jk(wfn.memory_jk_)
    else:
        # Normal initialization (no DF guess)
        # Shared JK pre-initialization reduces memory and computation time.
        #
        # Problem: Each wfn.initialize() creates its own JK via _build_jk()
        #          causing N× redundant 3-index integrals computation.
        #
        # Solution: Create single shared JK, then share it with all wfn.
        #          scf_initialize() is idempotent - reuses JK if already set (line 154-156)
        #
        # NOTE: 3-index integrals (Q|μν) depend only on basis + geometry,
        #       not on reference type (RHF/UHF/ROHF). JK object is generic
        #       and works with any number of density matrices.

        needs_jk_init = any(wfn.jk() is None for wfn in wfn_list)

        if needs_jk_init:
            # Build single shared JK for all wavefunctions
            ref_wfn = wfn_list[0]

            # Use integer arithmetic to avoid float-to-int type error
            # Python 3: `/` returns float, but C++ expects int
            # Use `//` (integer division) + explicit int() conversion
            total_memory_bytes = core.get_memory()
            safety_factor = core.get_global_option("SCF_MEM_SAFETY_FACTOR")
            total_memory = int((total_memory_bytes // 8) * safety_factor)

            # Create JK for reference wfn
            shared_jk = _build_jk(ref_wfn, total_memory)

            # Initialize JK (computes 3-index integrals - expensive!)
            # This happens ONCE for all wavefunctions!
            ref_wfn.initialize_jk(total_memory, jk=shared_jk)

            # Share JK with ALL other wavefunctions
            # Now they skip expensive _build_jk() and reuse shared JK!
            for wfn in wfn_list[1:]:
                wfn.set_jk(shared_jk)

            if verbose:
                core.print_out(f"  Shared JK object created for {len(wfn_list)} wavefunctions.\n")
                core.print_out("  Memory reduction: ~{}× for 3-index integrals!\n\n".format(len(wfn_list)))

        # Now initialize all wavefunctions
        # scf_initialize() sees JK already set → reuses it (idempotent!)
        # Only initializes per-wfn components: H, S^-1/2, guess, DIIS, PSIO
        for wfn in wfn_list:
            wfn.initialize()  # Reuses shared JK!

    # Phase 3: Main iterations (DIRECT if using DF guess, otherwise normal)
    return _multi_scf_inner(wfn_list, e_conv, d_conv, max_iter, verbose)


def _multi_scf_inner(wfn_list, e_conv, d_conv, max_iter, verbose):
    """
    Inner multi-SCF iteration loop.

    This function is called by multi_scf() and performs the actual SCF iterations.
    It is separated to allow DF_SCF_GUESS to run DF pre-iterations followed by
    DIRECT final iterations (Andy Trick 2.0).

    Parameters
    ----------
    wfn_list : list of Wavefunction
        List of wavefunctions to converge simultaneously
    e_conv : float
        Energy convergence threshold
    d_conv : float
        Density convergence threshold
    max_iter : int
        Maximum number of SCF iterations
    verbose : bool
        Verbosity level

    Returns
    -------
    list of float
        Final energies for each wavefunction
    """
    if verbose:
        core.print_out("\n  ==> Multi-Cycle SCF (NEW Architecture) <==\n\n")
        core.print_out("  Number of wavefunctions: {}\n".format(len(wfn_list)))
        core.print_out("  E convergence: {:.2e}\n".format(e_conv))
        core.print_out("  D convergence: {:.2e}\n".format(d_conv))
        core.print_out("  Maximum iterations: {}\n\n".format(max_iter))

        # Print wavefunction info
        for i, wfn in enumerate(wfn_list):
            wfn_type = wfn.__class__.__name__
            n_states = wfn.n_states()
            core.print_out("  Wavefunction {}: {} (n_states={})\n".format(i, wfn_type, n_states))
        core.print_out("\n")

    # Get JK object from first wavefunction (all must use same basis)
    jk = wfn_list[0].jk()
    if jk is None:
        raise ValidationError("JK object not initialized after wfn.initialize(). This should not happen.")

    # Ensure JK is configured to compute J and K matrices
    jk.set_do_J(True)
    jk.set_do_K(True)  # Needed for hybrid functionals (HF exchange)

    # Initialize iteration state for each wavefunction
    # Snapshot already applied above, so reads from frozen options
    for wfn in wfn_list:
        wfn._scf_initialize_iteration_state(e_conv, d_conv)

    # Track convergence status for each wfn
    converged_flags = [False] * len(wfn_list)
    final_energies = [0.0] * len(wfn_list)

    if verbose:
        header = "  " + "Iter".rjust(4)
        for i in range(len(wfn_list)):
            header += "  Energy_{}".format(i).rjust(20)
        header += "  " + "Status".rjust(25)
        core.print_out(header + "\n")
        core.print_out("  " + "-" * (len(header) - 2) + "\n")

    # Main multi-cycle SCF iteration loop
    for iteration in range(1, max_iter + 1):

        # Step 1: Collect occupied orbital matrices from all wavefunctions
        # All wfn participate in JK (including converged ones) to maintain
        # consistent indexing and prevent discontinuous Fock changes.
        # Converged wfn continue providing C matrices but don't iterate.
        all_C_occ_matrices = []
        wfn_state_counts = []  # Track how many states each wfn has
        active_wfn_indices = []  # Track which wfn need iteration (non-converged)

        for i, wfn in enumerate(wfn_list):
            # Get current C matrices for this wfn
            C_matrices = wfn.get_orbital_matrices()

            if verbose >= 2:
                status = "CONVERGED" if converged_flags[i] else "ACTIVE"
                core.print_out(f"  [DEBUG iter={iteration}] wfn {i}: {status}, collected {len(C_matrices)} C matrices\n")

            all_C_occ_matrices.extend(C_matrices)
            wfn_state_counts.append(len(C_matrices))

            # Track which wfn still need iteration (not converged)
            if not converged_flags[i]:
                active_wfn_indices.append(i)

        # Early exit if all converged (no active wfn left)
        if not active_wfn_indices:
            if verbose:
                core.print_out("\n  All wavefunctions converged!\n\n")
            break

        # Step 2: Shared JK computation for ALL wavefunctions (coupled convergence)
        # Use exported wrapper methods instead of direct vector manipulation
        # C_clear() clears both C_left and C_right
        # C_add() appends to both C_left and C_right (symmetric JK)
        jk.C_clear()
        for C_occ in all_C_occ_matrices:
            jk.C_add(C_occ)

        # Single JK call for ALL wavefunctions (including converged)
        jk.compute()

        # Step 3: Distribute J/K/wK results back to ALL wavefunctions
        # This maintains consistent Fock operators even after some wfn converge
        jk_index = 0
        J_all = jk.J()
        K_all = jk.K()
        wK_all = jk.wK()  # Long-range K for LRC functionals (empty if not LRC)

        # Diagnostic check (can be removed after testing)
        expected_matrices = len(all_C_occ_matrices)
        if len(J_all) != expected_matrices or len(K_all) != expected_matrices:
            raise ValidationError(
                f"JK compute failed: expected {expected_matrices} J/K matrices, "
                f"got {len(J_all)} J and {len(K_all)} K matrices. "
                f"C_left has {len(jk.C_left())} matrices, C_right has {len(jk.C_right())} matrices."
            )

        # Note: wK_all may be empty if not using LRC functional, that's OK
        # C++ set_jk_matrices() has default empty vector for wK_list

        # Distribute to ALL wfn (maintains consistent indexing)
        for i, wfn in enumerate(wfn_list):
            n_states = wfn_state_counts[i]
            J_list = [J_all[jk_index + j] for j in range(n_states)]
            K_list = [K_all[jk_index + j] for j in range(n_states)]
            wK_list = [wK_all[jk_index + j] for j in range(n_states)] if wK_all else []
            wfn.set_jk_matrices(J_list, K_list, wK_list)
            jk_index += n_states

        # Step 4: Each wavefunction completes its SCF iteration
        # Iterate through ALL wfn (to track status), but only active ones call _scf_iteration()
        # This is where ALL features work: DIIS, damping, SOSCF, MOM, FRAC, etc.
        all_converged = True
        status_strs = []

        for i, wfn in enumerate(wfn_list):
            if converged_flags[i]:
                # Converged, skip iteration
                status_strs.append("CONV")
                continue

            # Call the refactored _scf_iteration() method
            # This uses precomputed J/K automatically (use_precomputed_jk_ flag)
            should_continue, reason = wfn._scf_iteration()

            if not should_continue:
                if reason == 'converged':
                    # Mark as converged
                    converged_flags[i] = True
                    final_energies[i] = wfn.get_energies("Total Energy")
                    status_strs.append("CONV")
                    if verbose >= 2:
                        core.print_out(f"  [DEBUG iter={iteration}] wfn {i} CONVERGED\n")
                elif reason in ['mom_not_started', 'frac_not_started']:
                    # Special cases - keep iterating
                    all_converged = False
                    status_strs.append(reason[:4].upper())
                elif reason == 'early_screening_maxiter':
                    # COSX finished final grid iterations
                    converged_flags[i] = True
                    final_energies[i] = wfn.get_energies("Total Energy")
                    status_strs.append("COSX")
                else:
                    # Unknown reason
                    status_strs.append(reason[:4].upper())
            else:
                # Continue iterating
                all_converged = False
                status_strs.append("----")
                final_energies[i] = wfn.get_energies("Total Energy")

        # Print iteration info
        if verbose:
            line = "  " + str(iteration).rjust(4)
            for E in final_energies:
                line += "  {:20.14f}".format(E)
            line += "  " + "/".join(status_strs).rjust(25)
            core.print_out(line + "\n")

        # Check if all converged
        if all(converged_flags):
            if verbose:
                core.print_out("\n  All wavefunctions converged!\n\n")
            break

    else:
        # Max iterations reached without full convergence
        not_converged = [i for i, flag in enumerate(converged_flags) if not flag]
        raise SCFConvergenceError(
            "Multi-cycle SCF: {} wavefunction(s) did not converge in {} iterations (indices: {})".format(
                len(not_converged), max_iter, not_converged),
            iteration, wfn_list[0], 0.0, 0.0)

    return final_energies


def efp_field_fn(xyz):
    """Callback function for PylibEFP to compute electric field from electrons
    in ab initio part for libefp polarization calculation.

    Parameters
    ----------
    xyz : list
        (3 * npt, ) flat array of points at which to compute electric field

    Returns
    -------
    list
        (3 * npt, ) flat array of electric field at points in `xyz`.

    Notes
    -----
    Function signature defined by libefp, so function uses number of
    basis functions and integrals factory `mints_psi4_yo` and total density
    matrix `efp_Dt_psi4_yo` from global namespace.

    """
    points = core.Matrix.from_array(np.array(xyz).reshape(-1, 3))
    field = mints_psi4_yo.electric_field_value(points, efp_Dt_psi4_yo).np.flatten()
    return field
