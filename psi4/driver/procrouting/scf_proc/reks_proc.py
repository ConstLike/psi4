#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
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
REKS reference: SCF option setup and GUESS READ orbital/FON transport
"""

import re
from typing import List, Tuple

import numpy as np

from psi4 import core

from ...p4util.exceptions import ValidationError
from .subclass_methods import _RHF_orbital_gradient

core.REKS.compute_orbital_gradient = _RHF_orbital_gradient


_MO_TOKEN_RE = re.compile(r'^(HOMO|LUMO)([+-]\d+)?$')


def _parse_mo_token(token: str, n_alpha: int) -> int:
    """Resolve a single GUESS_ACTIVE_WINDOW MO token to a 0-based MO index."""
    tok = token.strip().upper()
    if tok.lstrip('-').isdigit():
        return int(tok)
    m = _MO_TOKEN_RE.match(tok)
    if not m:
        raise ValidationError(f"GUESS_ACTIVE_WINDOW bad token: {token!r}")
    base = (n_alpha - 1) if m.group(1) == 'HOMO' else n_alpha
    return base + (int(m.group(2)) if m.group(2) else 0)


def _window_to_swaps(window: List, dest: List[int], n_alpha: int, nmo: int) -> List[Tuple[int, int]]:
    """Resolve a GUESS_ACTIVE_WINDOW spec into the swap pairs that move the listed
    MOs into columns `dest`, in the listed order."""
    if len(window) != len(dest):
        raise ValidationError(
            f"GUESS_ACTIVE_WINDOW lists {len(window)} orbitals but the active space holds "
            f"{len(dest)}: {list(window)}")
    src = [_parse_mo_token(str(tok), n_alpha) for tok in window]
    for tok, i in zip(window, src):
        if not (0 <= i < nmo):
            raise ValidationError(
                f"GUESS_ACTIVE_WINDOW entry {tok!r} resolves to MO {i}, outside the {nmo} "
                "orbitals of the guess wavefunction")
    repeated = sorted({i for i in src if src.count(i) > 1})
    if repeated:
        raise ValidationError(
            f"GUESS_ACTIVE_WINDOW names MO {repeated} more than once: {list(window)}")
    if dest and dest[-1] >= nmo:
        raise ValidationError(
            f"GUESS_ACTIVE_WINDOW: active columns reach {dest[-1]}, past the {nmo} orbitals "
            "of the guess wavefunction")

    cur = list(range(nmo))  # cur[k] = guess MO currently in column k
    swaps: List[Tuple[int, int]] = []
    for p, i in zip(dest, src):
        k = cur.index(i)
        if k == p:
            continue
        swaps.append((p, k))
        cur[p], cur[k] = cur[k], cur[p]
    return swaps


def _apply_mo_swaps(C: core.Matrix, swaps: List[Tuple[int, int]]) -> None:
    """Swap MO columns in-place (C1)."""
    A = C.np
    for i, j in swaps:
        A[:, [i, j]] = A[:, [j, i]]


def _sym_lowdin(block, S):
    """Symmetric Lowdin orthonormalization in metric S: C -> C (C^T S C)^{-1/2}."""
    T = block.T @ S @ block
    w, V = np.linalg.eigh(T)
    w = np.maximum(w, 1e-10)
    return block @ ((V * (1.0 / np.sqrt(w))) @ V.T)


def _reks_block_basis_projection(Ca_in, old_wfn, scf_wfn):
    """Werner-Knowles projection to a new basis/geometry per [core | active | virtual] block, keeping the
    GVB pairs (active_k, active_{M-1-k}); C1 symmetry."""
    reks_spec = core.get_option("SCF", "REKS")
    N_act_e   = int(round(reks_spec[0]))
    n_act     = int(round(reks_spec[1]))
    G         = 2  # GVB(2): orbitals per geminal pair (structural invariant)
    nelectron = old_wfn.nalpha() + old_wfn.nbeta()
    Ncore     = (nelectron - N_act_e) // 2
    nmo       = old_wfn.nmopi()[0]
    n_virt    = nmo - Ncore - n_act
    n_pairs   = n_act // G

    old_basis = old_wfn.basisset()
    new_basis = scf_wfn.basisset()
    nso_new   = new_basis.nbf()

    Ca_np = np.asarray(Ca_in)
    active_cols = list(range(Ncore, Ncore + n_act))

    if core.get_option("SCF", "REKS_REPORT_LEVEL") >= 4:
        core.print_out(
            f"  REKS block projection: Ncore={Ncore} n_active={n_act} n_pairs={n_pairs} "
            f"n_virtual={n_virt}\n")

    # Projection per block (orthonormal within each block only).
    # Polly, R. et al. Mol. Phys. 2004, 102, 2311.
    blocks = []
    if Ncore > 0:
        blocks.append(list(range(0, Ncore)))
    blocks.append(active_cols)
    if n_virt > 0:
        blocks.append(list(range(Ncore + n_act, nmo)))

    C_proj = np.zeros((nso_new, nmo))
    for cols in blocks:
        ncol = len(cols)
        block_mat = core.Matrix.from_array(Ca_np[:, cols].copy())
        dim_col = core.Dimension([ncol])
        pblock = scf_wfn.basis_projection(block_mat, dim_col, old_basis, new_basis)
        C_proj[:, cols] = np.asarray(pblock)

    # Gram-Schmidt core > active > virtual with in-block Lowdin; the order makes the virtual step exact.
    # Werner, H.-J.; Knowles, P. J. J. Chem. Phys. 1985, 82, 5053.
    S_new = np.asarray(scf_wfn.S())

    # Active <- project out core, then in-block Lowdin.
    if Ncore > 0:
        core_cols = list(range(0, Ncore))
        C_core = C_proj[:, core_cols]
        overlap = C_core.T @ S_new @ C_proj[:, active_cols]
        C_proj[:, active_cols] -= C_core @ overlap
        C_proj[:, active_cols] = _sym_lowdin(C_proj[:, active_cols], S_new)

    # Virtual <- project out (core + active), then in-block Lowdin.
    if n_virt > 0:
        occ_act_cols = list(range(0, Ncore + n_act))
        virt_cols = list(range(Ncore + n_act, nmo))
        C_oa = C_proj[:, occ_act_cols]
        overlap = C_oa.T @ S_new @ C_proj[:, virt_cols]
        C_proj[:, virt_cols] -= C_oa @ overlap
        C_proj[:, virt_cols] = _sym_lowdin(C_proj[:, virt_cols], S_new)

    # Procrustes rotation of the active block onto the per-pair projections (keeps pair character and
    # orthogonality to core and virtual).
    C_ref_act = np.zeros((nso_new, n_act))
    for p in range(n_pairs):
        i_bond = Ncore + p
        i_anti = Ncore + (n_act - 1 - p)
        cols = [i_bond] if i_bond == i_anti else [i_bond, i_anti]
        block_mat = core.Matrix.from_array(Ca_np[:, cols].copy())
        dim_col = core.Dimension([len(cols)])
        pblock = np.asarray(scf_wfn.basis_projection(block_mat, dim_col, old_basis, new_basis))
        for k, c in enumerate(cols):
            C_ref_act[:, c - Ncore] = pblock[:, k]

    C_g_act = C_proj[:, active_cols]
    M = C_g_act.T @ S_new @ C_ref_act
    U_svd, sigma, Vt_svd = np.linalg.svd(M)
    R_procrustes = Vt_svd.T @ U_svd.T
    C_proj[:, active_cols] = C_g_act @ R_procrustes

    if core.get_option("SCF", "REKS_REPORT_LEVEL") >= 4:
        max_offdiag = np.max(np.abs(R_procrustes - np.diag(np.diag(R_procrustes))))
        core.print_out(
            f"  REKS Procrustes: max|R_offdiag| = {max_offdiag:.6e} "
            f"sigma = [{' '.join(f'{s:.8f}' for s in sigma)}]\n")

    return core.Matrix.from_array(C_proj)


def _geometry_follow_projection(Ca_in, target_wfn, n_low=None):
    """Same coefficients in target_wfn's metric: C_low (C_low^T S C_low)^{-1/2} over the first n_low columns
    (all if None), virtuals as their orthonormal complement."""
    S_new = np.asarray(target_wfn.S())
    C = np.asarray(Ca_in)
    nmo = C.shape[1]
    if n_low is None or n_low <= 0 or n_low >= nmo:
        n_low = nmo

    C_low = _sym_lowdin(C[:, :n_low], S_new)
    if n_low >= nmo:
        return core.Matrix.from_array(C_low)
    C_virt = C[:, n_low:]
    C_virt = C_virt - C_low @ (C_low.T @ S_new @ C_virt)
    C_virt = _sym_lowdin(C_virt, S_new)
    return core.Matrix.from_array(np.hstack([C_low, C_virt]))


def force_scf_options(post_scf):
    """Set the SCF options a REKS reference requires; refuse REKS as a post-SCF reference."""
    if post_scf:
        raise ValidationError("REFERENCE REKS is a complete calculation and is not a reference for "
                              "post-SCF methods.")
    # REKS runs its own GVB-DIIS/TRAH; the single-determinant accelerators are off.
    core.set_local_option("SCF", "DIIS", False)
    core.set_local_option("SCF", "SOSCF", False)
    core.set_local_option("SCF", "SCF_INITIAL_ACCELERATOR", "NONE")
    # Default level shift 0.6 held to convergence (no cutoff).
    if not core.has_option_changed("SCF", "LEVEL_SHIFT"):
        core.set_local_option("SCF", "LEVEL_SHIFT", 0.6)
    if not core.has_option_changed("SCF", "LEVEL_SHIFT_CUTOFF"):
        core.set_local_option("SCF", "LEVEL_SHIFT_CUTOFF", 0.0)


def build_wavefunction(ref_wfn, superfunc):
    """Build the REKS wavefunction, refusing core features REKS does not support."""
    unsupported = [name for name in ("PCM", "DDX", "PE") if core.get_option("SCF", name)]
    if hasattr(ref_wfn.molecule(), "EFP"):
        unsupported.append("EFP")
    if superfunc.is_c_hybrid():
        unsupported.append(f"double-hybrid functional {superfunc.name()}")
    if core.get_option("SCF", "MOM_START") != 0:
        unsupported.append("MOM_START")
    if core.get_option("SCF", "FRAC_START") != 0:
        unsupported.append("FRAC_START")
    if core.get_option("SCF", "DAMPING_PERCENTAGE") > 0.0:
        unsupported.append("DAMPING_PERCENTAGE")
    if core.get_option("SCF", "STABILITY_ANALYSIS") != "NONE":
        unsupported.append("STABILITY_ANALYSIS")
    if sum(core.get_option("SCF", "TDSCF_STATES")) > 0:
        unsupported.append("TDSCF_STATES")
    if superfunc.vv10_b() > 0.0:
        unsupported.append(f"VV10 nonlocal correlation ({superfunc.name()})")
    if core.get_option("SCF", "DFT_GRAC_SHIFT") != 0.0:
        unsupported.append("DFT_GRAC_SHIFT")
    if core.get_option("SCF", "ORBITAL_OPTIMIZER_PACKAGE") != "INTERNAL":
        unsupported.append("ORBITAL_OPTIMIZER_PACKAGE")
    if unsupported:
        raise ValidationError("Not supported for REFERENCE REKS: " + ", ".join(unsupported) + ".")
    return core.REKS(ref_wfn, superfunc)


def validate_active_window(scf_wfn):
    """Check GUESS_ACTIVE_WINDOW against the reference and the guess of scf_wfn."""
    if core.get_option('SCF', 'GUESS_ACTIVE_WINDOW'):
        if not isinstance(scf_wfn, core.REKS):
            core.print_out("  GUESS_ACTIVE_WINDOW ignored: it fills a REKS active space and "
                           f"this run is {core.get_option('SCF', 'REFERENCE')}.\n")
        elif core.get_option('SCF', 'GUESS') != 'READ':
            raise ValidationError(
                "GUESS_ACTIVE_WINDOW reorders the orbitals taken from the READ guess and needs "
                f"GUESS READ; this REKS run guesses with {core.get_option('SCF', 'GUESS')}, "
                "which builds its orbitals inside the SCF and leaves nothing to reorder.")


def read_guess(scf_wfn, old_wfn, read_filename, scf_molecule):
    """Set the REKS guess orbitals and FONs of scf_wfn from old_wfn."""
    # REKS needs the full MO set: active virtuals carry occupation.
    old_is_reks = (old_wfn.name() == "REKS")

    Ca_in = old_wfn.Ca_subset("SO", "ALL")
    Cb_in = old_wfn.Cb_subset("SO", "ALL")
    # Rows of the saved orbitals must match the source basis as rebuilt under this run's PUREAM.
    if Ca_in.rowdim().sum() != old_wfn.basisset().nbf():
        raise ValidationError(
            f"GUESS READ: the saved orbitals have {Ca_in.rowdim().sum()} basis-function rows, but the "
            f"source basis rebuilt under this run's PUREAM has {old_wfn.basisset().nbf()}; the file "
            "was written with the other PUREAM setting. Use the PUREAM of the source calculation.")

    window = core.get_option('SCF', 'GUESS_ACTIVE_WINDOW')
    if window:
        dest = scf_wfn.get_active_mo_indices()
        swaps = _window_to_swaps(window, dest, old_wfn.nalpha(), Ca_in.coldim()[0])
        _apply_mo_swaps(Ca_in, swaps)
        if Cb_in is not Ca_in:
            _apply_mo_swaps(Cb_in, swaps)
        core.print_out(f"  GUESS READ: active window {list(window)} -> MO columns {dest} "
                       f"(swaps {swaps})\n")

    # Direct copy only when basis name and geometry match.
    same_basis = (old_wfn.basisset().name() == scf_wfn.basisset().name())
    old_geom = np.asarray(old_wfn.molecule().geometry())
    new_geom = np.asarray(scf_molecule.geometry())
    same_geometry = (old_geom.shape == new_geom.shape
                     and np.allclose(old_geom, new_geom, atol=1e-12, rtol=0.0))

    # Dispatch on (same basis?, same geometry?, old wfn is REKS?): copy, transport, project, or both.
    if same_basis and same_geometry:
        core.print_out(f"  Reading orbitals from file {read_filename}, no projection.\n\n")
        pCa, pCb = Ca_in, Cb_in

    elif old_is_reks and same_basis and not same_geometry:
        core.print_out(
            f"  Reading orbitals from file {read_filename}, following orbitals to the new\n"
            "  geometry (diabatic transport: fixed MO coefficients + block Lowdin).\n\n")
        core.timer_on("REKS: guess geometry follow")
        n_low = scf_wfn.get_Ncore() + scf_wfn.n_active_orbitals()
        pCa = _geometry_follow_projection(Ca_in, scf_wfn, n_low)
        core.timer_off("REKS: guess geometry follow")
        pCb = pCa

    elif old_is_reks and same_geometry:
        # Pair-preserving Werner projection.
        core.print_out(
            f"  Reading orbitals from file {read_filename}, projecting to new basis.\n")
        core.print_out("  Computing basis projection from %s to %s\n\n" % (
            old_wfn.basisset().name(), scf_wfn.basisset().name()))
        core.timer_on("REKS: guess basis projection")
        pCa = _reks_block_basis_projection(Ca_in, old_wfn, scf_wfn)
        core.timer_off("REKS: guess basis projection")
        pCb = pCa

    elif old_is_reks:
        # Project at old geometry, then transport.
        core.print_out(
            f"  Reading orbitals from file {read_filename}, projecting to new basis then\n"
            "  following to new geometry.\n\n")
        core.timer_on("REKS: guess basis+geometry projection")
        old_geom_new_basis = core.Wavefunction.build(
            old_wfn.molecule(), core.get_global_option('BASIS'))
        pCa_old_geom = _reks_block_basis_projection(Ca_in, old_wfn, old_geom_new_basis)
        n_low = scf_wfn.get_Ncore() + scf_wfn.n_active_orbitals()
        pCa = _geometry_follow_projection(pCa_old_geom, scf_wfn, n_low)
        core.timer_off("REKS: guess basis+geometry projection")
        pCb = pCa

    else:
        if not same_basis:
            core.print_out(
                f"  Reading orbitals from file {read_filename}, projecting to new basis.\n")
            core.print_out("  Computing basis projection from %s to %s\n\n" % (
                old_wfn.basisset().name(), scf_wfn.basisset().name()))
        else:
            core.print_out(
                f"  Reading orbitals from file {read_filename}, projecting to new geometry.\n")
            core.print_out(
                "  Same basis name but atoms moved; cross-geometry basis_projection\n"
                "  (Polly et al., Mol. Phys. 102, 2311 (2004)).\n\n")
        noccpi_a = old_wfn.nmopi()
        noccpi_b = old_wfn.nmopi()
        pCa = scf_wfn.basis_projection(Ca_in, noccpi_a, old_wfn.basisset(), scf_wfn.basisset())
        if old_wfn.same_a_b_orbs():
            pCb = pCa
        else:
            pCb = scf_wfn.basis_projection(Cb_in, noccpi_b, old_wfn.basisset(), scf_wfn.basisset())

    scf_wfn.guess_Ca(pCa)
    scf_wfn.guess_Cb(pCb)

    # REKS -> REKS hot restart: transport the FON of every generation and run sector.
    if old_is_reks:
        if old_wfn.has_array_variable("REKS FON"):
            old_n_e = int(round(old_wfn.scalar_variable("REKS N ACTIVE ELECTRONS")))
            old_n_o = int(round(old_wfn.scalar_variable("REKS N ACTIVE ORBITALS")))
            new_reks_spec = core.get_option("SCF", "REKS")
            new_n_e = int(round(new_reks_spec[0]))
            new_n_o = int(round(new_reks_spec[1]))
            if (old_n_e, old_n_o) == (new_n_e, new_n_o):
                n_sectors = 1
                if old_wfn.has_scalar_variable("REKS N SECTORS"):
                    n_sectors = max(1, int(round(old_wfn.scalar_variable("REKS N SECTORS"))))

                def _sec_psivar(s, reks_name):
                    # Sector 0 keeps reks_name; s > 0 inserts "SECTOR s " after the "REKS " prefix.
                    if s == 0:
                        return reks_name
                    tail = reks_name[5:] if reks_name.startswith("REKS ") else reks_name
                    return f"REKS SECTOR {s} {tail}"

                for s in range(n_sectors):
                    sec_tag = "" if s == 0 else f"sector {s} "
                    fon_name = _sec_psivar(s, "REKS FON")
                    if old_wfn.has_array_variable(fon_name):
                        sec_fon = old_wfn.array_variable(fon_name)
                        scf_wfn.guess_fon(sec_fon, s)
                        if core.get_option("SCF", "REKS_REPORT_LEVEL") >= 4:
                            fon_flat = sec_fon.to_array().flatten().tolist()
                            core.print_out(
                                f"  REKS READ: transporting {sec_tag}FON {fon_flat} "
                                f"from saved wavefunction\n")
                    # Higher FON generations (m, u, v, w), restored where present.
                    for gen, name in ((1, "REKS M_FON"), (2, "REKS U_FON"),
                                      (3, "REKS V_FON"), (4, "REKS W_FON")):
                        upper_name = _sec_psivar(s, name)
                        if old_wfn.has_array_variable(upper_name):
                            upper_mat = old_wfn.array_variable(upper_name)
                            scf_wfn.guess_fon_upper(gen, upper_mat, s)
                            if core.get_option("SCF", "REKS_REPORT_LEVEL") >= 4:
                                upper_flat = upper_mat.to_array().flatten().tolist()
                                core.print_out(
                                    f"  REKS READ: transporting {sec_tag}{name} {upper_flat} "
                                    f"from saved wavefunction\n")
            else:
                core.print_out(
                    f"  REKS READ: active space mismatch "
                    f"({old_n_e},{old_n_o}) -> ({new_n_e},{new_n_o}); "
                    f"FON not restored, falling back to default\n")
        else:
            core.print_out(
                "  REKS READ: saved wavefunction has no FON slot; "
                "falling back to default (cold-start convergence expected)\n")


def missing_restart_error(name, read_filename):
    """Build the ValidationError for a REKS GUESS READ without a saved wavefunction."""
    return ValidationError(
        f"REKS GUESS=READ: no saved orbitals to read.\n"
        f"  looked for: {read_filename}\n"
        "Pick one:\n"
        "\n"
        "# A. two calls in one run: 1st saves, 2nd reads\n"
        "    psi4.set_options({'guess': 'sad'})\n"
        f"    psi4.energy('{name}')\n"
        "    psi4.set_options({'guess': 'read'})\n"
        f"    psi4.energy('{name}')\n"
        "\n"
        "# B. save to a file, restart in a later run\n"
        "    psi4.set_options({'guess': 'sad'})\n"
        f"    e, wfn = psi4.energy('{name}', return_wfn=True)\n"
        "    wfn.to_file('orbitals.npy')\n"
        f"    psi4.energy('{name}', restart_file='orbitals.npy')\n"
        "\n"
        "# C. Psithon input deck\n"
        "    set guess sad\n"
        f"    energy('{name}')\n"
        "    set guess read\n"
        f"    energy('{name}')\n"
        "\n"
        "# D. optimization from the file saved in B\n"
        f"    psi4.optimize('{name}', restart_file='orbitals.npy')\n"
        "\n"
        "# Tip. Choose the guess MOs that fill the active space, in column order\n"
        "    psi4.set_options({'guess_active_window': ['HOMO', 'LUMO+1']})")
