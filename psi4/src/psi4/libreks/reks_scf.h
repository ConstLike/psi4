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

#include "reks_cassette.h"
#include "reks_reel.h"

#include <memory>
#include <vector>

namespace psi {
class JK;
class SuperFunctional;

namespace reks {
namespace xc {
class MicrostateXC;
}
namespace scf {

/// C_L over all microstates: catalog L from the config weights at reel.fon_state, extras at their fixed weight.
void compute_weights(const studio::Reel&       reel,
                     const studio::Cassette&   sa_cassette,
                     std::vector<double>&    C_L);

/// State-averaged energy scalar, summed over the SCF pool:
///   E_SA = sum_{L in sa_cassette.microstates()} C_L[L] * E_L[L]
double compute_E_SA(const studio::Cassette&      sa_cassette,
                    const std::vector<double>& E_L,
                    const std::vector<double>& C_L);

/// True once compute_sa_energies has filled E_L: a computed microstate energy carries
/// E_const, so |E_L| of the first SA microstate is nonzero.
bool sa_energies_ready(const studio::Cassette& sa_cassette, const std::vector<double>& E_L);

struct MicrostateEnergy {
    double total;          ///< E_1e + E_2e + E_const (+ xc_correction)
    double E_1e;           ///< tr(Da H) + tr(Db H)
    double E_2e;           ///< 1/2 (E_Fa + E_Fb - E_1e)
    double tr_D_Vxc;       ///< tr(D V_xc)[L]; 0 when apply_xc is false
    double xc_correction;  ///< E_xc[L] - 1/2 tr(D V_xc)[L]; 0 when apply_xc is false
};

/// Energy of microstate L from the cached Fock scalars (reel.E_Fa_L/E_Fb_L) and
/// one-electron traces (reel.base_density_e1[alpha_idx/beta_idx]); E_const = E_nuc + E_-D.
MicrostateEnergy microstate_energy(const studio::Reel& reel, int L,
                                   int alpha_idx, int beta_idx, double E_const,
                                   bool apply_xc);

/// D_core = sum_{c<Ncore} c_c c_c^T and, per referenced bitmask p (and 0),
/// base_density_e1[p] = tr(D_core H) + sum_{i in p} c_i^T H c_i.
void build_base_densities(const studio::Cassette&   cassette,
                          const std::vector<int>& referenced_bitmasks,
                          studio::Reel&             reel,
                          int                     Ncore,
                          const std::vector<int>& active_mo,
                          const SharedMatrix&     Ca,
                          const SharedMatrix&     H);

/// One JK pass into reel J/K/wK core and per-active-orbital caches (K_act also on pair_k_demand); for a pure
/// GGA also G_out = 2 J - alpha K - beta wK of the closed-shell occupied block.
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

/// Grid-direct XC per L, no AO V_xc^L: active rows over [mo_lo, mo_lo + mo_width) and diagonal (need_diag) of
/// C^T V_xc^L C, tr(rho^L V_xc^L), E_xc_L; V_acc_xc = sum_L C_L V_xc^L if want_v_acc.
void build_microstate_xc(
    studio::Reel&                   reel,
    const studio::Cassette&         sa_cassette,
    const std::vector<int>&         microstate_Ls,
    int                           Ncore,
    psi::SuperFunctional&         functional,
    reks::xc::MicrostateXC*       xc,
    const SharedMatrix&           Ca_full,
    const std::vector<double>&    C_L,
    bool                          need_diag,
    int                           mo_lo,
    int                           mo_width,
    bool                          want_v_acc);

/// Which exchange terms enter the Fock, and at what scale.
struct FockExchangePolicy {
    bool   apply_K  = false;   ///< is_x_hybrid && !reel.K_act.empty()
    double k_scale  = 0.0;     ///< functional.x_alpha()
    bool   apply_wK = false;   ///< is_x_lrc && !reel.wK_act.empty()
    double wk_scale = 0.0;     ///< functional.x_beta()
};
FockExchangePolicy fock_exchange_policy(const studio::Reel& reel,
                                        psi::SuperFunctional& functional);

/// Add the Coulomb pattern to out: out += wJcore*J_core + sum_i wJ[i]*J_act[i].
/// J_core added only when Ncore>0 and present; wJ has length n_active.
void add_coulomb_pattern(SharedMatrix out, const studio::Reel& reel,
                         int Ncore, int n_active, double wJcore, const double* wJ);

/// out -= k_scale (wKcore K_core + sum_i wK[i] K_act[i]) if apply_K, and the same with wk_scale and the wK
/// matrices if apply_wK.
void add_exchange_pattern(SharedMatrix out, const studio::Reel& reel,
                          int Ncore, int n_active,
                          const FockExchangePolicy& pol, double wKcore, const double* wK);

/// Per L: E_Fa/Fb_L = tr(D F_L) (+ tr(D V_xc)), MO diagonal (all MOs if need_full_diag, else active) and active
/// rows over the build_microstate_xc window of F_L + V_xc^L, active block of F_a - F_b; D_core used windowed.
void build_sa_focks_MO(
    studio::Reel&                        reel,
    const studio::Cassette&              sa_cassette,
    const std::vector<int>&              microstate_Ls,
    int                                Ncore,
    psi::SuperFunctional&              functional,
    const SharedMatrix&                H,
    const SharedMatrix&                D_core,
    bool                               need_full_diag,
    int                                mo_lo,
    int                                mo_width);

/// Per-L Fock scalars and active-block fock_aa of SI microstates outside the SA set, from the last SCF JK
/// cache and Ca_fock_transform; their per-L rows are released again.
void fill_missing_si_focks(
    studio::Reel&                          reel,
    const studio::Cassette&                sa_cassette,
    const std::vector<studio::Cassette>&   si_cassettes,
    int                                Ncore,
    psi::SuperFunctional&              functional,
    reks::xc::MicrostateXC*            xc,
    const SharedMatrix&                Ca_full,
    const SharedMatrix&                H);

/// E_L of the cassette's microstates outside the SA set, from the cached per-L Fock/XC scalars and
/// base_density_e1.
void build_si_microstate_focks(const studio::Reel&                  reel,
                               const studio::Cassette&              cassette,
                               const studio::Cassette&              sa_cassette,
                               double                             E_const,
                               bool                               needs_xc,
                               std::vector<double>&               E_L);

/// Symmetric MO coupling Fock F_reks_MO over the SA microstates (bulk from reel.F_acc_MO) and the
/// Lagrangian multipliers reel.lagrangians.
void assemble_F_reks_MO(studio::Reel&                        reel,
                        const studio::Cassette&              sa_cassette,
                        const std::vector<double>&         C_L,
                        const std::vector<int>&            active_mo,
                        int                                Ncore,
                        int                                N,
                        SharedMatrix                       F_reks_MO);

/// 2e Hessian diagonals gamma_aa, gamma_ca = 2 sum_L C_L A2e(L,pq) from the pair ERIs (pp|kk), (pk|pk) in
/// reel.slot_diag (full MO extent).
void compute_gamma_vectors(int                        Ncore,
                           int                        n_act,
                           const studio::Reel&        reel,
                           const studio::Cassette&    sa_cassette,
                           double                     x_alpha,
                           int                        iteration,
                           std::vector<double>&       gamma_aa,
                           std::vector<double>&       gamma_ca);

}  // namespace scf
}  // namespace reks
}  // namespace psi
