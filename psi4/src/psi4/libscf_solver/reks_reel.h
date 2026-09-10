/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#pragma once

#include "reks_arch.h"

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"

#include <memory>
#include <vector>

namespace psi {
namespace reks {

class ActiveEriTiles;

namespace studio {

struct Reel {
    /// Owned SCF-MO snapshot (clone/copy of Ca_).
    SharedMatrix Ca_fock_transform;

    /// One FONSnapshot per run sector, index = sector ordinal; each entry's layer
    /// count follows that sector's n_generations.
    std::vector<FONSnapshot> fon_state;

    /// Per-microstate energy E_L[L] and weight C_L[L]; E_SA = sum_L C_L[L] E_L[L].
    std::vector<double>       E_L;
    std::vector<double>       C_L;

    /// Per-microstate MO Fock, stored as diagonal + active rows over the MO column
    /// window the microstate's pool needs. Sized lazily per-L at first Fock build.
    ///   F_MO_diag_*_L[L]  : length mo_width (diagonal over the window).
    ///   F_MO_arows_*_L[L] : n_active x mo_width (active rows over the window).
    /// SA pool: mo_lo = 0, mo_width = nmo (full row). Every other microstate:
    /// mo_lo = Ncore, mo_width = n_active (active block only).
    std::vector<std::vector<double>> F_MO_diag_a_L;
    std::vector<std::vector<double>> F_MO_diag_b_L;
    std::vector<SharedMatrix>        F_MO_arows_a_L;
    std::vector<SharedMatrix>        F_MO_arows_b_L;

    /// Active-space block of the per-microstate MO Fock spin difference,
    /// (F_MO_a - F_MO_b)(p, Ncore + q) over p, q < n_active, row-major per microstate.
    ///   fock_aa_row[L] : row of fock_aa, -1 until L is first built
    ///   fock_aa        : rows of n_active^2, appended in first-build order
    std::vector<int>    fock_aa_row;
    std::vector<double> fock_aa;

    /// MO diagonal of each Fock slot, [slot * slot_diag_stride + q] over q < nmo,
    /// from the last full-diagonal AO->MO transform. Slot J_act[k] holds
    /// diag_q = (qq|kk), slot K_act[k] holds (qk|qk).
    ///   slot_diag_stride : nmo, or 0 when the last transform was windowed.
    ///   slot_diag_J/_K   : first J / first K slot; slot_diag_K < 0 when the
    ///                      functional applies no exact exchange.
    std::vector<double> slot_diag;
    int                 slot_diag_stride = 0;
    int                 slot_diag_J = -1;
    int                 slot_diag_K = -1;

    /// C_L-weighted accumulated AO Fock, Sum_L C_L (Fa_AO_L + Fb_AO_L).
    SharedMatrix F_acc_AO;   ///< AO accumulator (nso x nso)
    SharedMatrix F_acc_MO;   ///< C^T F_acc_AO C (nmo x nmo); sized lazily to nmo

    /// AO Fock buffer, reused each L.
    SharedMatrix F_alpha_AO_buffer;

    /// Cached per-L scalars.
    ///   E_Fa_L/E_Fb_L    : tr(D_a.F_a) / tr(D_b.F_b), incl. tr(D.Vxc) when XC is active.
    ///   tr_D_Vxc_a_L/b_L : tr(D_a.Vxc_a) / tr(D_b.Vxc_b).
    std::vector<double> E_Fa_L;
    std::vector<double> E_Fb_L;
    std::vector<double> tr_D_Vxc_a_L;
    std::vector<double> tr_D_Vxc_b_L;

    /// Per-microstate XC (grid-direct): MO quantities of C^T V_xc^L C, active rows +
    /// diagonal only; the full AO V_xc^L is never formed. Vxc_arows shares
    /// F_MO_arows' MO column window, sized lazily per-L at first Fock build;
    /// Vxc_diag stays empty when the full MO diagonal is not wanted.
    std::vector<SharedMatrix>        Vxc_arows_a_L;  ///< [L] n_active x mo_width
    std::vector<SharedMatrix>        Vxc_arows_b_L;
    std::vector<std::vector<double>> Vxc_diag_a_L;   ///< [L] length nmo, or empty
    std::vector<std::vector<double>> Vxc_diag_b_L;
    std::vector<double>       E_xc_L;

    /// C_L-weighted AO back-projection Sum_L C_L V_xc^L_{a,b}.
    SharedMatrix V_acc_xc_a;  ///< nso x nso
    SharedMatrix V_acc_xc_b;

    /// Coupling Fock.
    SharedMatrix F_reks;      ///< AO (nso x nso)
    SharedMatrix F_reks_MO;   ///< MO (nmo x nmo); lazy: null until first use

    /// Core density D_core = sum_{c < Ncore} c_c c_c^T (nso x nso).
    SharedMatrix D_core;

    /// tr(D_p . H) for the occupation pattern p, with
    /// D_p = D_core + sum_{i: p >> i & 1} c_act_i c_act_i^T. Bitmask-indexed,
    /// size = 2^M (M = active orbital count); only the referenced patterns are filled.
    std::vector<double> base_density_e1;

    /// JK cache (lazy: null until first build)
    SharedMatrix              J_core;
    SharedMatrix              K_core;
    SharedMatrix              wK_core;
    std::vector<SharedMatrix> J_act;
    std::vector<SharedMatrix> K_act;
    std::vector<SharedMatrix> wK_act;

    /// Lagrangians (size catalog.data->n_lagrangian_pairs)
    std::vector<double> lagrangians;

    /// True after SCF; assemble_F_reks_MO must not run again.
    bool lagrangians_frozen = false;

    /// Scratch
    SharedMatrix temp_buffer;
    SharedMatrix C_core_buffer;              ///< Ncore cols of Ca (nso x Ncore)
    std::vector<SharedMatrix> c_act_buffer;  ///< one MO column each (nso x 1)

    /// TRAH active-block rotation: run-constant shapes.
    SharedMatrix trah_K;           ///< antisymmetric K then exp(K_act) in place (n_act x n_act)
    SharedMatrix trah_Ca_act;      ///< gathered active columns of Ca (nso x n_act)
    SharedMatrix trah_Ca_act_rot;  ///< rotated active columns Ca_act * U (nso x n_act)
};

/// n_extra extends every per-L array to n_microstates + n_extra; L >= n_microstates
/// are the runtime SA extra determinants. Catalog-derived sizes (base_density_e1,
/// lagrangians) are unchanged. MO-sized objects (F_MO_diag/arows, Vxc_arows/diag,
/// F_acc_MO, F_reks_MO) are sized lazily at first Fock build, since nmo is set
/// only after form_Shalf.
void allocate_reel(Reel& reel,
                   Catalog catalog,
                   const Dimension& nsopi,
                   int nso,
                   int n_extra = 0);

}  // namespace studio
}  // namespace reks
}  // namespace psi
