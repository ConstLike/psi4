/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
#include "reks_cassette.h"

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/exception.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace psi {
namespace reks {
namespace studio {

/// Per-run REKS state (FONs, per-L Fock/energy caches, coupling Fock, JK caches, Lagrangians, TRAH
/// scratch); MO-sized members are sized at the first Fock build.
struct Reel {
    /// Owned copy of the SCF MO coefficients.
    SharedMatrix Ca_fock_transform;

    /// One FONSnapshot per run sector, index = sector ordinal; each entry's layer
    /// count follows that sector's n_generations.
    std::vector<FONSnapshot> fon_state;

    /// Per-microstate energy E_L[L] and weight C_L[L]; E_SA = sum_L C_L[L] E_L[L].
    std::vector<double>       E_L;
    std::vector<double>       C_L;

    /// Per-L MO Fock diagonal (mo_width) and active rows (n_active x mo_width); window = all nmo for the SA
    /// pool, the active block otherwise.
    std::vector<std::vector<double>> F_MO_diag_a_L;
    std::vector<std::vector<double>> F_MO_diag_b_L;
    std::vector<SharedMatrix>        F_MO_arows_a_L;
    std::vector<SharedMatrix>        F_MO_arows_b_L;

    /// Active block (F_MO_a - F_MO_b)(p, Ncore + q) per L: fock_aa rows of n_active^2 in first-build order,
    /// fock_aa_row[L] its row (-1 before).
    std::vector<int>    fock_aa_row;
    std::vector<double> fock_aa;

    /// MO diagonal per Fock slot [slot * stride + q] (stride nmo, 0 if windowed): J_act[k] -> (qq|kk),
    /// K_act[k] -> (qk|qk); slot_diag_K < 0 without exact exchange.
    std::vector<double> slot_diag;
    int                 slot_diag_stride = 0;
    int                 slot_diag_J = -1;
    int                 slot_diag_K = -1;

    /// C_L-weighted accumulated AO Fock, Sum_L C_L (Fa_AO_L + Fb_AO_L).
    SharedMatrix F_acc_AO;   ///< nso x nso
    SharedMatrix F_acc_MO;   ///< C^T F_acc_AO C (nmo x nmo); sized lazily to nmo

    /// AO Fock buffer.
    SharedMatrix F_alpha_AO_buffer;

    /// Per-L tr(D_s F_s) (incl. tr(D V_xc) with XC) and tr(D_s V_xc_s).
    std::vector<double> E_Fa_L;
    std::vector<double> E_Fb_L;
    std::vector<double> tr_D_Vxc_a_L;
    std::vector<double> tr_D_Vxc_b_L;

    /// Per-L C^T V_xc^L C: active rows (F_MO_arows window) and diagonal (empty unless wanted).
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

    /// tr(D_p H), D_p = D_core + sum_{i in p} c_i c_i^T, per bitmask p < 2^M (referenced ones filled).
    std::vector<double> base_density_e1;

    /// JK cache; null until first build.
    SharedMatrix              J_core;
    SharedMatrix              K_core;
    SharedMatrix              wK_core;
    std::vector<SharedMatrix> J_act;
    std::vector<SharedMatrix> K_act;
    std::vector<SharedMatrix> wK_act;

    /// Lagrangians, size catalog.data->n_lagrangian_pairs.
    std::vector<double> lagrangians;

    SharedMatrix temp_buffer;
    SharedMatrix C_core_buffer;              ///< Ncore cols of Ca (nso x Ncore)
    std::vector<SharedMatrix> c_act_buffer;  ///< one MO column each (nso x 1)

    /// TRAH active-block rotation buffers.
    SharedMatrix trah_K;           ///< antisymmetric K then exp(K_act) in place (n_act x n_act)
    SharedMatrix trah_Ca_act;      ///< gathered active columns of Ca (nso x n_act)
    SharedMatrix trah_Ca_act_rot;  ///< rotated active columns Ca_act * U (nso x n_act)
};

/// FONs of block (s, gen): position i = geminal sa_cassette.geminals_active(s, gen)[i] of reel.fon_state[s].
inline std::vector<double> fon_vector(const Reel& reel, const Cassette& sa_cassette, int s, int gen) {
    const auto& act = sa_cassette.geminals_active(s, gen);
    std::vector<double> v;
    v.reserve(act.size());
    for (int g : act) v.push_back(reel.fon_state[s].layers[gen][g].p);
    return v;
}
inline void set_fon_vector(Reel& reel, const Cassette& sa_cassette, int s, int gen,
                           const std::vector<double>& fons) {
    const auto& act = sa_cassette.geminals_active(s, gen);
    if (fons.size() != act.size())
        throw PSIEXCEPTION("set_fon_vector: size does not match active geminal pool");
    for (size_t i = 0; i < act.size(); ++i)
        reel.fon_state[s].layers[gen][act[i]] = {fons[i], 2.0 - fons[i]};
}

/// FON snapshot indexed by FON block (parallel to fon_blocks); entry b holds
/// fon_vector(fon_blocks[b].s, fon_blocks[b].gen).
inline std::vector<std::vector<double>> capture_fons(const Reel& reel, const Cassette& sa_cassette,
                                                     const std::vector<FonBlock>& fon_blocks) {
    std::vector<std::vector<double>> out(fon_blocks.size());
    for (size_t b = 0; b < fon_blocks.size(); ++b)
        out[b] = fon_vector(reel, sa_cassette, fon_blocks[b].s, fon_blocks[b].gen);
    return out;
}
inline void restore_fons(Reel& reel, const Cassette& sa_cassette, const std::vector<FonBlock>& fon_blocks,
                         const std::vector<std::vector<double>>& fons) {
    for (size_t b = 0; b < fons.size() && b < fon_blocks.size(); ++b)
        set_fon_vector(reel, sa_cassette, fon_blocks[b].s, fon_blocks[b].gen, fons[b]);
}

/// Per-L arrays sized n_microstates + n_extra (extras at the tail); MO-sized objects wait for the first
/// Fock build (nmo is known after form_Shalf).
void allocate_reel(Reel& reel,
                   Catalog catalog,
                   const Dimension& nsopi,
                   int nso,
                   int n_extra = 0);

}  // namespace studio
}  // namespace reks
}  // namespace psi
