
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
    """Base class Wavefunction requires this function. Here it is
    simply a wrapper around initialize(), iterations(), finalize_energy(). It
    returns the SCF energy computed by finalize_energy().

    """
    if core.get_option('SCF', 'DF_SCF_GUESS') and (core.get_global_option('SCF_TYPE') == 'DIRECT'):
        # speed up DIRECT algorithm (recomputes full (non-DF) integrals
        #   each iter) by first converging via fast DF iterations, then
        #   fully converging in fewer slow DIRECT iterations. aka Andy trick 2.0
        core.print_out("  Starting with a DF guess...\n\n")
        with p4util.OptionsStateCM(['SCF_TYPE']):
            core.set_global_option('SCF_TYPE', 'DF')
            self.initialize()
            try:
                self.iterations()
            except SCFConvergenceError:
                self.finalize()
                raise SCFConvergenceError("""SCF DF preiterations""", self.iteration_, self, 0, 0)
        core.print_out("\n  DF guess converged.\n\n")

        # reset the DIIS & JK objects in prep for DIRECT
        if self.initialized_diis_manager_:
            self.diis_manager_.reset_subspace()
        self.initialize_jk(self.memory_jk_)
    else:
        self.initialize()
    self.iteration_energies = []

    try:
        self.iterations()
    except SCFConvergenceError as e:
        if core.get_option("SCF", "FAIL_ON_MAXITER"):
            core.print_out("  Failed to converge.\n")
            # energy = 0.0
            # A P::e fn to either throw or protest upon nonconvergence
            # die_if_not_converged()
            raise e
        else:
            core.print_out("  Energy and/or wave function did not converge, but proceeding anyway.\n\n")
    else:
        core.print_out("  Energy and wave function converged.\n\n")

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

    # Store configuration flags
    self._scf_is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')
    self._scf_verbose = get_option_from_snapshot(self, 'PRINT')
    self._scf_reference = core.get_option('SCF', "REFERENCE")  # OK from global (doesn't change)
    self._scf_damping_enabled = _validate_damping(self)
    self._scf_soscf_enabled = _validate_soscf(self)
    self._scf_frac_enabled = _validate_frac(self)
    self._scf_efp_enabled = hasattr(self.molecule(), 'EFP')
    self._scf_cosx_enabled = "COSX" in core.get_option('SCF', 'SCF_TYPE')  # OK from global

    # CRITICAL: Set DIIS/MOM members that are used by _scf_iteration()
    # These were originally in scf_iterate() but needed for multi_scf() too
    self.diis_enabled_ = self.validate_diis()
    self.MOM_excited_ = _validate_MOM(self)
    self.diis_start_ = get_option_from_snapshot(self, 'DIIS_START')

    # COSX early screening parameters
    self._scf_early_screening = False
    if self._scf_cosx_enabled:
        self._scf_early_screening = True
        self.jk().set_COSX_grid("Initial")

    self._scf_maxiter_post_screening = get_option_from_snapshot(self, 'COSX_MAXITER_FINAL')
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

        if self.iteration_ >= get_option_from_snapshot(self, 'MAXITER'):
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
    if (self._scf_soscf_enabled and (self.iteration_ >= 3) and (self._scf_Dnorm < get_option_from_snapshot(self, 'SOSCF_START_CONVERGENCE'))):
        self._scf_Dnorm = self.compute_orbital_gradient(False, get_option_from_snapshot(self, 'DIIS_MAX_VECS'))
        diis_performed = False
        if self.functional().needs_xc():
            base_name = "SOKS, nmicro="
        else:
            base_name = "SOSCF, nmicro="

        if not _converged(self._scf_Ediff, self._scf_Dnorm, e_conv=self._scf_e_conv, d_conv=self._scf_d_conv):
            nmicro = self.soscf_update(get_option_from_snapshot(self, 'SOSCF_CONV'),
                                       get_option_from_snapshot(self, 'SOSCF_MIN_ITER'),
                                       get_option_from_snapshot(self, 'SOSCF_MAX_ITER'),
                                       get_option_from_snapshot(self, 'SOSCF_PRINT'))
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
            diis_performed = False
            add_to_diis_subspace = self.diis_enabled_ and self.iteration_ >= self.diis_start_

            self._scf_Dnorm = self.compute_orbital_gradient(add_to_diis_subspace, get_option_from_snapshot(self, 'DIIS_MAX_VECS'))

            if add_to_diis_subspace:
                for engine_used in self.diis(self._scf_Dnorm):
                    status.append(engine_used)

            core.timer_off("HF: DIIS")

            if self._scf_verbose > 4 and diis_performed:
                core.print_out("  After DIIS:\n")
                self.Fa().print_out()
                self.Fb().print_out()

            # frac, MOM invoked here from Wfn::HF::find_occupation
            core.timer_on("HF: Form C")
            level_shift = get_option_from_snapshot(self, 'LEVEL_SHIFT')
            if level_shift > 0 and self._scf_Dnorm > get_option_from_snapshot(self, 'LEVEL_SHIFT_CUTOFF'):
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
    if (self._scf_damping_enabled and self.iteration_ > 1 and self._scf_Dnorm > get_option_from_snapshot(self, 'DAMPING_CONVERGENCE')):
        damping_percentage = get_option_from_snapshot(self, 'DAMPING_PERCENTAGE')
        self.damping_update(damping_percentage * 0.01)
        status.append("DAMP={}%".format(round(damping_percentage)))

    if core.has_option_changed("SCF", "ORBITALS_WRITE"):
        filename = core.get_option("SCF", "ORBITALS_WRITE")
        self.to_file(filename)

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
    if self.MOM_excited_ and ((not self.MOM_performed_) or self.iteration_ == get_option_from_snapshot(self, 'MOM_START')):
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
    >>> # Converge singlet and triplet RHF simultaneously
    >>> rhf_singlet = psi4.core.RHF(wfn_singlet, ...)
    >>> rhf_triplet = psi4.core.RHF(wfn_triplet, ...)
    >>> energies = multi_scf([rhf_singlet, rhf_triplet])

    Notes
    -----
    - All wavefunctions must use the same basis set
    - All wavefunctions must use the same JK algorithm (SCF_TYPE)
    - Each wfn maintains its own DIIS, damping, convergence state
    - ~1.8-2x faster than running SCF cycles independently
    - Uses existing C++ infrastructure (use_precomputed_jk_ flag)
    """
    if len(wfn_list) == 0:
        raise ValidationError("multi_scf requires at least one wavefunction")

    # Get convergence criteria
    if e_conv is None:
        e_conv = core.get_option('SCF', 'E_CONVERGENCE')
    if d_conv is None:
        d_conv = core.get_option('SCF', 'D_CONVERGENCE')
    if max_iter is None:
        max_iter = core.get_option('SCF', 'MAXITER')

    # CRITICAL: Snapshot global options ONCE before creating/initializing any wfn
    # This prevents non-determinism from global state pollution between wfn creation
    options_snapshot = snapshot_scf_options()

    # CRITICAL: Apply options snapshot to ALL wfn BEFORE any initialization
    # This ensures wfn.initialize() reads from frozen snapshot, not from global state
    for wfn in wfn_list:
        apply_options_snapshot(wfn, options_snapshot)

    # Initialize wavefunctions if not already done
    # After snapshot applied, wfn.initialize() reads from frozen options
    for wfn in wfn_list:
        if wfn.jk() is None:
            wfn.initialize()

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

        # Step 1: Collect occupied orbital matrices from ALL wavefunctions
        # CRITICAL: ALL wfn must participate in JK (including converged ones)
        # to maintain consistent indexing and prevent discontinuous Fock changes.
        #
        # FREEZE PATTERN FIX (Step 1.5.1):
        # When a wfn converges, form_C() updates Ca_ one last time, causing
        # other wfn to see different J/K on the next iteration (discontinuity).
        # Solution: Freeze C matrices at the state BEFORE form_C() modified them.
        # Since get_orbital_matrices() already deep-copies via get_block(),
        # we just need to save references when wfn converges.
        all_C_occ_matrices = []
        wfn_state_counts = []  # Track how many states each wfn has
        active_wfn_indices = []  # Track which wfn need iteration (non-converged)

        for i, wfn in enumerate(wfn_list):
            if converged_flags[i]:
                # Use frozen C matrices from convergence iteration
                # This prevents other wfn from seeing discontinuous J/K
                C_matrices = wfn._frozen_C_for_jk
                if verbose >= 2:
                    core.print_out(f"  [DEBUG iter={iteration}] wfn {i}: FROZEN, using {len(C_matrices)} frozen C matrices\n")
            else:
                # Get current C matrices (get_block() already deep-copies data)
                C_matrices = wfn.get_orbital_matrices()

                # Save snapshot for potential freezing if wfn converges this iteration
                # No clone() needed - get_orbital_matrices() already returned new matrices
                wfn._C_snapshot_for_jk = C_matrices

                if verbose >= 2:
                    core.print_out(f"  [DEBUG iter={iteration}] wfn {i}: ACTIVE, collected {len(C_matrices)} fresh C matrices\n")

            all_C_occ_matrices.extend(C_matrices)
            wfn_state_counts.append(len(C_matrices))

            # Track which wfn still need iteration
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
                # Already converged, skip iteration
                status_strs.append("CONV")
                continue

            # Call the refactored _scf_iteration() method
            # This uses precomputed J/K automatically (use_precomputed_jk_ flag)
            should_continue, reason = wfn._scf_iteration()

            if not should_continue:
                if reason == 'converged':
                    # FREEZE PATTERN: Save C matrices from BEFORE form_C() modified them
                    # This prevents discontinuity in J/K for other wfn on next iteration
                    wfn._frozen_C_for_jk = wfn._C_snapshot_for_jk
                    converged_flags[i] = True
                    final_energies[i] = wfn.get_energies("Total Energy")
                    status_strs.append("CONV")
                    if verbose >= 2:
                        core.print_out(f"  [DEBUG iter={iteration}] wfn {i} CONVERGED, freezing {len(wfn._frozen_C_for_jk)} C matrices\n")
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
