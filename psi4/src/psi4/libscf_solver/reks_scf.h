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
#include "reks_cassette.h"
#include "reks_reel.h"

#include "psi4/libmints/dimension.h"

#include <memory>
#include <vector>

namespace psi {
class JK;
class SuperFunctional;
class RV;
class UV;

namespace reks {
namespace scf {

/// Aggregate per-K weights into the per-microstate C_L vector (sized
/// n_total_microstates: catalog L's from compute_C_L, extra determinants at
/// their fixed ensemble weight). Pre: reel.fon_state is current.
void compute_weights(const studio::Reel&       reel,
                     const studio::Cassette&   sa_cassette,
                     std::vector<double>&    C_L);

/// State-averaged energy scalar, summed over the SCF pool:
///   E_SA = sum_{L in sa_cassette.microstates()} C_L[L] * E_L[L]
double compute_E_SA(const studio::Cassette&      sa_cassette,
                    const std::vector<double>& E_L,
                    const std::vector<double>& C_L);

struct MicrostateEnergy {
    double total;          ///< E_1e + E_2e + E_nuc (+ xc_correction)
    double E_1e;           ///< tr(Da H) + tr(Db H)
    double E_2e;           ///< 1/2 (E_Fa + E_Fb - E_1e)
    double xc_correction;  ///< E_xc[L] - 1/2 tr(D V_xc)[L]; 0 when apply_xc is false
};

/// Energy of microstate L from the cached Fock scalars (reel.E_Fa_L/E_Fb_L) and
/// one-electron traces (reel.base_density_e1[alpha_idx/beta_idx]).
MicrostateEnergy microstate_energy(const studio::Reel& reel, int L,
                                   int alpha_idx, int beta_idx, double E_nuc,
                                   bool apply_xc);

/// Refresh the core density and the one-electron trace of every referenced
/// occupation pattern:
///   reel.D_core            = sum_{c<Ncore} c_c c_c^T
///   D_p                    = D_core + sum_{i: p>>i & 1} c_act_i c_act_i^T
///   reel.base_density_e1[p] = tr(D_p . H) = tr(D_core . H) + sum_{i in p} c_i^T H c_i
/// Only D_core is materialized; every other pattern is carried by its scalar trace.
///
/// Pre:
///   reel.D_core allocated and reel.base_density_e1 sized to 2^n_active_orbitals.
///   referenced_bitmasks: alpha/beta occupation bitmasks to populate (caller-selected).
///   active_mo: MO indices of the active orbitals (size = n_active).
///   Ca: SO-basis MO coefficients (only the core columns and active_mo read).
///   H: one-electron Hamiltonian.
/// Post:
///   reel.D_core rebuilt from Ca; reel.base_density_e1[p] populated for p in
///   referenced_bitmasks (and 0).
void build_base_densities(const studio::Cassette&   cassette,
                          const std::vector<int>& referenced_bitmasks,
                          studio::Reel&             reel,
                          int                     Ncore,
                          const std::vector<int>& active_mo,
                          const SharedMatrix&     Ca,
                          const SharedMatrix&     H);

/// JK build + per-orbital cache + closed-shell HF G assembly (HF part only;
/// XC not included in G_out):
///   G_out = 2 J_at_c_occ - alpha K_at_c_occ - beta wK_at_c_occ
/// Clones JK output into reel.J_core/.K_core/.wK_core/.J_act/.K_act/.wK_act;
/// aliases (not clones) the closed-shell HF slot J_at_c_occ/K_at_c_occ/wK_at_c_occ.
///
/// `pair_k_demand` requests reel.K_act irrespective of the functional: a build
/// without exact exchange then pays one extra J-off pass over the n_active columns.
///
/// Pre:
///   C_occ_cache: SO-basis occupied orbital block.
/// Post:
///   reel.J_core/.K_core/.wK_core/.J_act/.K_act/.wK_act populated.
///   reel.K_act holds plain exchange K(c_i c_i^T) or is empty; a wcombine JK folds K
///     into wK and leaves it empty.
///   reel.K_core populated only when the Fock itself applies exact exchange.
///   J_at_c_occ/K_at_c_occ/wK_at_c_occ and G_out populated only when
///     functional.needs_xc() && functional.is_gga() && functional.x_alpha() < 0.01;
///     left untouched otherwise.
void build_jk_cache(
    const studio::Cassette&         cassette,
    studio::Reel&                   reel,
    int                           Ncore,
    const std::vector<int>&       active_mo,
    const SharedMatrix&           Ca,
    const SharedMatrix&           C_occ_cache,
    psi::JK&                      jk,
    psi::SuperFunctional&         functional,
    SharedMatrix&                 J_at_c_occ,
    SharedMatrix&                 K_at_c_occ,
    SharedMatrix&                 wK_at_c_occ,
    SharedMatrix                  G_out,
    bool                          pair_k_demand);

/// Per-L grid-direct microstate XC, without forming a per-L AO V_xc:
///   reel.Vxc_arows_a/b_L = active rows of C^T V_xc^L C
///   reel.Vxc_diag_a/b_L  = MO diagonal (need_diag only)
///   reel.tr_D_Vxc_a/b_L  = tr(rho^L_s V_xc^L_s)
///   reel.E_xc_L          = grid quadrature XC energy
/// Accumulates the C_L-weighted AO back-projection Sum_L C_L V_xc^L into
/// reel.V_acc_xc_a/b.
///
/// Allocations of reel.Vxc_arows_*_L/Vxc_diag_*_L/E_xc_L/tr_D_Vxc_*_L stay sized to
/// n_total_microstates; only the cassette-listed slots are filled per call.
///
/// Pre:
///   Ca_full: Fock-basis MOs (reel.Ca_fock_transform); active_mo contiguous (Ncore+i).
///   microstate_Ls: work set in the run-wide (global) L numbering (catalog
///     records, then extras tail).
/// Post:
///   no-op (nothing written) when !needs_xc, uv_potential null, or microstate_Ls empty.
///   uv_potential->quadrature_values() RHO_A/RHO_B refreshed with the C_L-weighted
///     density integral (grid-electron diagnostic).
/// `mo_lo` / `mo_width` select the MO column window the active rows carry; the SA pool
/// takes the full row (0, nmo), every other pool the active block (Ncore, n_active).
/// `want_v_acc` writes reel.V_acc_xc_{a,b}; without it the grid sweep drops the whole
/// C_L-weighted AO back-projection.
void build_microstate_xc(
    studio::Reel&                   reel,
    const studio::Cassette&         sa_cassette,
    const std::vector<int>&         microstate_Ls,
    int                           Ncore,
    psi::SuperFunctional&         functional,
    const std::shared_ptr<psi::UV>&  uv_potential,
    const SharedMatrix&           Ca_full,
    const std::vector<int>&       active_mo,
    const std::vector<double>&    C_L,
    bool                          need_diag,
    int                           mo_lo,
    int                           mo_width,
    bool                          want_v_acc);

/// Which exchange terms enter the Fock, and at what scale.
struct FockExchangePolicy {
    bool   apply_K  = false;   ///< is_x_hybrid && !reel.K_act.empty() && !wcombine
    double k_scale  = 0.0;     ///< functional.x_alpha()
    bool   apply_wK = false;   ///< is_x_lrc && !reel.wK_act.empty()
    double wk_scale = 0.0;     ///< wcombine ? 1.0 : functional.x_beta()
};
FockExchangePolicy fock_exchange_policy(const studio::Reel& reel,
                                        psi::SuperFunctional& functional,
                                        psi::JK& jk);

/// Add the Coulomb pattern to out: out += wJcore*J_core + sum_i wJ[i]*J_act[i].
/// J_core added only when Ncore>0 and present; wJ has length n_active.
void add_coulomb_pattern(SharedMatrix out, const studio::Reel& reel,
                         int Ncore, int n_active, double wJcore, const double* wJ);

/// Subtract the exchange pattern from out, per policy:
///   out -= k_scale *(wKcore*K_core  + sum_i wK[i]*K_act[i])   when pol.apply_K
///   out -= wk_scale*(wKcore*wK_core + sum_i wK[i]*wK_act[i])  when pol.apply_wK.
/// wK has length n_active; core terms added only when Ncore>0 and present.
void add_exchange_pattern(SharedMatrix out, const studio::Reel& reel,
                          int Ncore, int n_active,
                          const FockExchangePolicy& pol, double wKcore, const double* wK);

/// Per-microstate AO Fock assembly + AO->MO transform.
/// For each L in microstate_Ls:
///   - assembles F_AO = H + J_pattern + (hybrid)*K_pattern + (lrc)*wK_pattern
///     (V_xc is NOT in the AO Fock; it enters the MO Fock as the grid-direct
///     diagonal + active rows)
///   - caches reel.E_Fa/Fb_L[L] = D_a/b . F_AO (+ tr(D.V_xc) when needs_xc)
///   - extracts the MO Fock (diagonal + active rows) into
///     reel.F_MO_diag_*_L[L] / reel.F_MO_arows_*_L[L].
/// Only the diagonal and active rows of F_MO are formed; the off-diagonal
/// (uniform-weight) bulk is never materialized here.
/// `microstate_Ls`: work set in the run-wide (global) L numbering.
/// `need_full_diag` selects the diagonal extent: true builds it over all nmo,
/// false over the active window [Ncore, Ncore + n_active) only. The energy trace
/// is exact either way. Pass the same value as build_microstate_xc's `need_diag`.
///
/// Pre:
///   reel.Ca_fock_transform = clone of the Fock-basis Ca (set by caller).
///   reel.J_core/.K_core/.wK_core/.J_act/.K_act/.wK_act populated by
///     scf::build_jk_cache.
///   reel.Vxc_arows_a/b_L / Vxc_diag_a/b_L populated for those L's by
///     build_microstate_xc when functional.needs_xc().
///   D_core = C_core C_core^T over the FIRST Ncore columns of
///     reel.Ca_fock_transform. Read only when !need_full_diag && Ncore > 0; may
///     be null otherwise.
/// `mo_lo` / `mo_width` select the MO column window every per-L row carries, and must
/// match the window build_microstate_xc was given for the same L's. The window has to
/// cover the active block, so mo_lo <= Ncore and mo_lo + mo_width >= Ncore + n_active.
/// Post:
///   reel.E_Fa_L[L] / reel.E_Fb_L[L] populated for the L's.
///   reel.F_MO_arows_*_L[L] populated for the L's, n_active x mo_width.
///   reel.F_MO_diag_*_L[L] sized mo_width for the L's, valid over the extent
///     `need_full_diag` selects and zero outside it.
///   reel.fock_aa_row[L] / reel.fock_aa carry the active block of F_MO_a - F_MO_b.
///   reel.temp_buffer clobbered when need_full_diag.
void build_sa_focks_MO(
    studio::Reel&                        reel,
    const studio::Cassette&              sa_cassette,
    const std::vector<int>&              microstate_Ls,
    int                                Ncore,
    psi::SuperFunctional&              functional,
    psi::JK&                           jk,
    const SharedMatrix&                H,
    const SharedMatrix&                D_core,
    bool                               need_full_diag,
    int                                mo_lo,
    int                                mo_width);

/// Post-SCF: fills the per-L Fock cache for si_cassettes' microstates not in
/// the SA cassette's active set, occupation read via sa_cassette. Reuses the
/// JK cache from the last SCF iteration (Ca/Da unchanged); no build_jk_cache
/// call. Built over the active MO window with no MO diagonal; the active
/// block is stored in reel.fock_aa and the per-L rows released before return.
///
/// Pre:
///   reel.Ca_fock_transform: SCF-snapshot Ca (set per SCF iter by caller);
///     shared MO basis for SA and SI L's.
///   reel JK cache (J_core/J_act/...) valid from the last SCF iteration.
/// Post:
///   reel.fock_aa_row[L] / reel.fock_aa carry the active-block Fock difference
///   for each si_cassette microstate not in the SA cassette's active set;
///   reel.E_Fa/Fb_L, reel.E_xc_L, reel.tr_D_Vxc_*_L filled for them.
///   reel.F_MO_*/Vxc_* rows of those L's released again.
///   reel.temp_buffer clobbered.
void fill_missing_si_focks(
    studio::Reel&                          reel,
    const studio::Cassette&                sa_cassette,
    const std::vector<studio::Cassette>&   si_cassettes,
    int                                Ncore,
    psi::SuperFunctional&              functional,
    const std::shared_ptr<psi::UV>&    uv_potential,
    const SharedMatrix&                Ca_full,
    const std::vector<int>&            active_mo,
    psi::JK&                           jk,
    const SharedMatrix&                H);

/// Populate E_L[Ls] for catalog L's listed in cassette.microstates_active
/// but NOT in the SA cassette's active set, from cached per-L Fock and XC
/// scalars. No JK or XC compute calls.
///
/// Pre:
///   reel.E_Fa_L / .E_Fb_L / .tr_D_Vxc_a/b_L sized n_total_microstates and
///   populated for SI L's by scf::build_sa_focks_MO + build_microstate_xc.
///   reel.E_xc_L sized n_total_microstates when needs_xc.
///   reel.base_density_e1 indexed by occupation bitmask (size 2^n_active_orbitals),
///   populated by scf::build_base_densities.
/// Post:
///   E_L[Ls] written for Ls in cassette \ SA cassette's active set.
void build_si_microstate_focks(const studio::Reel&                  reel,
                               const studio::Cassette&              cassette,
                               const studio::Cassette&              sa_cassette,
                               double                             nuclear_repulsion,
                               bool                               needs_xc,
                               std::vector<double>&               E_L);

/// Build coupling Fock F_reks_MO and accumulate Lagrangian multipliers in
/// reel.lagrangians. Driven over the SCF pool sa_cassette.microstates().
///
/// Pre:
///   reel.fon_state is current.
///   reel.F_MO_diag_*_L[L] / F_MO_arows_*_L[L] populated by scf::build_sa_focks_MO,
///   reel.F_acc_MO populated by REKS::build_F_acc_MO (supplies the uniform bulk).
///   C_L[L] populated by scf::compute_weights; size = n_total_microstates.
/// Post:
///   F_reks_MO is symmetric, MO-basis coupling Fock.
///   reel.lagrangians resized to sa_cassette.n_lagrangian_pairs() and filled.
void assemble_F_reks_MO(studio::Reel&                        reel,
                        const studio::Cassette&              sa_cassette,
                        const std::vector<double>&         C_L,
                        const std::vector<int>&            active_mo,
                        int                                Ncore,
                        int                                N,
                        SharedMatrix                       F_reks_MO);

}  // namespace scf
}  // namespace reks
}  // namespace psi
