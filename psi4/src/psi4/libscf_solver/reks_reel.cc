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

#include "reks_reel.h"

namespace psi {
namespace reks {
namespace studio {

void allocate_reel(Reel& reel,
                   Catalog catalog,
                   const Dimension& nsopi,
                   int nso,
                   int n_extra) {
    const int n_micro = catalog.data->n_microstates;
    const int n_micro_total = n_micro + n_extra;
    const int n_act   = catalog.data->n_active_orbitals;
    const int n_base  = 1 << n_act;  // 2^M occupation patterns over the active space

    reel.E_L.assign(n_micro_total, 0.0);
    reel.C_L.assign(n_micro_total, 0.0);

    // Outer global-L vectors only; per-L contents stay empty here.
    reel.F_MO_diag_a_L .assign(n_micro_total, std::vector<double>{});
    reel.F_MO_diag_b_L .assign(n_micro_total, std::vector<double>{});
    reel.F_MO_arows_a_L.assign(n_micro_total, SharedMatrix{});
    reel.F_MO_arows_b_L.assign(n_micro_total, SharedMatrix{});

    reel.fock_aa_row.assign(n_micro_total, -1);
    reel.fock_aa.clear();

    reel.F_alpha_AO_buffer = std::make_shared<Matrix>(
        "F_alpha_AO_buffer", nsopi, nsopi);

    reel.F_acc_AO = std::make_shared<Matrix>("F_acc_AO", nsopi, nsopi);
    // F_acc_MO is MO-sized (nmo x nmo); sized lazily at first Fock build.

    reel.E_Fa_L      .assign(n_micro_total, 0.0);
    reel.E_Fb_L      .assign(n_micro_total, 0.0);
    reel.tr_D_Vxc_a_L.assign(n_micro_total, 0.0);
    reel.tr_D_Vxc_b_L.assign(n_micro_total, 0.0);

    // Per-microstate XC: outer global-L vectors only; per-L MO arows/diag stay
    // null here.
    reel.Vxc_arows_a_L.assign(n_micro_total, SharedMatrix{});
    reel.Vxc_arows_b_L.assign(n_micro_total, SharedMatrix{});
    reel.Vxc_diag_a_L .assign(n_micro_total, std::vector<double>{});
    reel.Vxc_diag_b_L .assign(n_micro_total, std::vector<double>{});
    reel.E_xc_L       .assign(n_micro_total, 0.0);
    // Weighted-AO accumulators: single matrices, not per-L.
    reel.V_acc_xc_a = std::make_shared<Matrix>("V_acc_xc_a", nsopi, nsopi);
    reel.V_acc_xc_b = std::make_shared<Matrix>("V_acc_xc_b", nsopi, nsopi);

    reel.F_reks = std::make_shared<Matrix>("F_REKS (coupling)", nsopi, nsopi);

    // Only the core density is materialized; every other occupation pattern is
    // carried by its one-electron trace, indexed by the same bitmask.
    reel.D_core = std::make_shared<Matrix>("D_core", nsopi, nsopi);
    reel.base_density_e1.assign(n_base, 0.0);

    reel.lagrangians.assign(catalog.data->n_lagrangian_pairs, 0.0);
    reel.lagrangians_frozen = false;

    reel.temp_buffer = std::make_shared<Matrix>("temp_buffer", nso, nso);
}

}  // namespace studio
}  // namespace reks
}  // namespace psi
