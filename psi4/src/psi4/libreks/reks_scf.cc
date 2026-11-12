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

#include "reks_scf.h"

#include "reks_math.h"
#include "reks_report_level.h"
#include "reks_xc.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <set>
#include <string>
#include <utility>

namespace psi {
namespace reks {
namespace scf {

void compute_weights(const studio::Reel&       reel,
                     const studio::Cassette&   sa_cassette,
                     std::vector<double>&    C_L) {
    const int n_catalog     = sa_cassette.n_catalog_microstates();
    const int n_micro_total = sa_cassette.n_total_microstates();
    C_L.assign(static_cast<size_t>(n_micro_total), 0.0);

    // C_K scratch holds only the catalog-L weights (compute_C_L's write range).
    thread_local std::vector<double> C_K;
    if (static_cast<int>(C_K.size()) < n_catalog) C_K.resize(n_catalog);
    for (int s = 0; s < sa_cassette.n_sectors(); ++s) {
        for (int i : sa_cassette.sector_config_positions(s)) {
            const int    K   = sa_cassette.K_indices[i];
            const double w_K = sa_cassette.K_weights[i];
            sa_cassette.compute_C_L(K, reel.fon_state[s], C_K);
            for (int L : sa_cassette.config_microstate_writes(K))
                C_L[L] += w_K * C_K[L];
        }
    }

    // Extra determinants carry a fixed, FON-independent ensemble weight.
    const int n_extra = sa_cassette.n_extra_microstates();
    for (int e = 0; e < n_extra; ++e)
        C_L[n_catalog + e] = sa_cassette.extra_weights[e];
}

double compute_E_SA(const studio::Cassette&      sa_cassette,
                    const std::vector<double>& E_L,
                    const std::vector<double>& C_L) {
    double E_SA = 0.0;
    for (int L : sa_cassette.microstates()) E_SA += C_L[L] * E_L[L];
    return E_SA;
}

bool sa_energies_ready(const studio::Cassette& sa_cassette, const std::vector<double>& E_L) {
    const auto& sa_micro = sa_cassette.microstates();
    return !sa_micro.empty() && std::abs(E_L[sa_micro.front()]) >= constants::ENERGY_THRESHOLD;
}

/// f_i = 1/2 [sum_K w_K sum_L C_L^(K)(n) (alpha_i + beta_i) + sum_e w_e^extra (alpha_i + beta_i)].
static void compute_all_f_eff(const studio::Cassette&                 sa_cassette,
                              const std::vector<studio::FONSnapshot>& fon_state,
                              std::vector<double>&                    f_eff) {
    const int n_micro  = sa_cassette.n_catalog_microstates();
    const int n_active = sa_cassette.n_active_orbitals();

    thread_local std::vector<double> C_K;
    if (static_cast<int>(C_K.size()) < n_micro) C_K.resize(n_micro);

    f_eff.assign(n_active, 0.0);
    for (int s = 0; s < sa_cassette.n_sectors(); ++s) {
        for (int i : sa_cassette.sector_config_positions(s)) {
            const int    K   = sa_cassette.K_indices[i];
            const double w_K = sa_cassette.K_weights[i];
            sa_cassette.compute_C_L(K, fon_state[s], C_K);
            for (int orbital_idx = 0; orbital_idx < n_active; ++orbital_idx) {
                double occ = 0.0;
                for (int L : sa_cassette.config_microstate_writes(K)) {
                    const auto& m = sa_cassette.microstate(L);
                    occ += C_K[L] * static_cast<double>(
                               m.alpha[orbital_idx] + m.beta[orbital_idx]);
                }
                f_eff[orbital_idx] += w_K * occ;
            }
        }
    }

    // Extra determinants enter the ensemble occupation at their fixed weight.
    const int n_extra = sa_cassette.n_extra_microstates();
    for (int e = 0; e < n_extra; ++e) {
        const auto& m = sa_cassette.microstate(n_micro + e);
        for (int orbital_idx = 0; orbital_idx < n_active; ++orbital_idx)
            f_eff[orbital_idx] += sa_cassette.extra_weights[e] * static_cast<double>(
                                      m.alpha[orbital_idx] + m.beta[orbital_idx]);
    }

    for (int orbital_idx = 0; orbital_idx < n_active; ++orbital_idx)
        f_eff[orbital_idx] *= 0.5;
}

MicrostateEnergy microstate_energy(const studio::Reel& reel, int L,
                                   int alpha_idx, int beta_idx, double E_const,
                                   bool apply_xc) {
    MicrostateEnergy e;
    e.E_1e = reel.base_density_e1[alpha_idx] + reel.base_density_e1[beta_idx];
    e.E_2e = 0.5 * (reel.E_Fa_L[L] + reel.E_Fb_L[L] - e.E_1e);
    e.tr_D_Vxc = apply_xc ? reel.tr_D_Vxc_a_L[L] + reel.tr_D_Vxc_b_L[L] : 0.0;
    e.xc_correction = apply_xc ? reel.E_xc_L[L] - 0.5 * e.tr_D_Vxc : 0.0;
    e.total = e.E_1e + e.E_2e + E_const + e.xc_correction;
    return e;
}

void build_base_densities(const studio::Cassette&   cassette,
                          const std::vector<int>& referenced_bitmasks,
                          studio::Reel&             reel,
                          int                     Ncore,
                          const std::vector<int>& active_mo,
                          const SharedMatrix&     Ca,
                          const SharedMatrix&     H) {
    const int nso      = Ca->rowspi()[0];
    const int nmo      = Ca->colspi()[0];
    const int n_active = cassette.n_active_orbitals();
    double**  Cp       = Ca->pointer(0);

    // The GEMM writes every element at beta = 0.0; only the coreless case needs
    // the explicit zero.
    if (Ncore > 0) {
        double** Dp = reel.D_core->pointer(0);
        C_DGEMM('N', 'T', nso, nso, Ncore, 1.0, Cp[0], nmo, Cp[0], nmo, 0.0,
                Dp[0], nso);
    } else {
        reel.D_core->zero();
    }
    const double e1_core = reel.D_core->vector_dot(H);
    reel.base_density_e1[0] = e1_core;

    // Active MO columns gathered as a contiguous nso x n_active block, then
    // T = H C_act, so h[i] = c_i^T H c_i is the dot of matching columns.
    std::vector<double> C_act(static_cast<size_t>(nso) * n_active);
    std::vector<double> T(static_cast<size_t>(nso) * n_active);
    for (int i = 0; i < n_active; ++i)
        C_DCOPY(nso, &Cp[0][active_mo[i]], nmo, C_act.data() + i, n_active);
    C_DGEMM('N', 'N', nso, n_active, nso, 1.0, H->pointer(0)[0], nso,
            C_act.data(), n_active, 0.0, T.data(), n_active);
    std::vector<double> h(n_active);
    for (int i = 0; i < n_active; ++i)
        h[i] = C_DDOT(nso, C_act.data() + i, n_active, T.data() + i, n_active);

    for (int p : referenced_bitmasks) {
        if (p == 0) continue;  // core-only pattern handled above
        double e1 = e1_core;
        for (int i = 0; i < n_active; ++i)
            if (p & (1 << i)) e1 += h[i];
        reel.base_density_e1[p] = e1;
    }
}

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
    bool                          pair_k_demand) {
    timer_on("REKS: build_jk_cache");
    const int nso       = Ca->rowspi()[0];
    const int nmo       = Ca->colspi()[0];
    const int n_active  = cassette.n_active_orbitals();

    double** Cp = Ca->pointer(0);
    // Persistent Reel buffer (run-constant shape), refreshed from Ca on entry.
    SharedMatrix C_core;
    if (Ncore > 0) {
        if (!reel.C_core_buffer)
            reel.C_core_buffer = std::make_shared<Matrix>("C_core", nso, Ncore);
        C_core = reel.C_core_buffer;
        double** Ccp = C_core->pointer(0);
        // Core MOs are Ca columns [0, Ncore).
        for (int mu = 0; mu < nso; ++mu)
            std::memcpy(Ccp[mu], Cp[mu], static_cast<size_t>(Ncore) * sizeof(double));
    }
    if (static_cast<int>(reel.c_act_buffer.size()) != n_active)
        reel.c_act_buffer.resize(n_active);
    std::vector<SharedMatrix>& c_act = reel.c_act_buffer;
    for (int i = 0; i < n_active; ++i) {
        if (!c_act[i])
            c_act[i] = std::make_shared<Matrix>("c_act_" + std::to_string(i), nso, 1);
        // Column active_mo[i] of Ca: row-major stride nmo into a unit-stride column.
        C_DCOPY(nso, &Cp[0][active_mo[i]], nmo, c_act[i]->pointer(0)[0], 1);
    }

    // One closed-shell density carries G whole only where no exact exchange enters it;
    // x_alpha < 0.01 is that test, a numerical zero on the exchange fraction.
    const bool build_c_occ = functional.needs_xc() && functional.is_gga() &&
                             functional.x_alpha() < 0.01;

    const int core_present = (Ncore > 0) ? 1 : 0;
    const int c_occ_slot   = core_present + n_active;
    const int n_pushed     = c_occ_slot + (build_c_occ ? 1 : 0);

    // clone() on first use, copy() thereafter: copy() deep-copies so the next
    // jk.compute() cannot alias the cache.
    auto refresh_one = [](SharedMatrix& dst, const SharedMatrix& src) {
        if (!dst) dst = src->clone();
        else      dst->copy(src);
    };
    auto refresh_vec = [](std::vector<SharedMatrix>& dst,
                          const std::vector<SharedMatrix>& src, int base, int n) {
        if (static_cast<int>(dst.size()) != n) dst.resize(n);
        for (int i = 0; i < n; ++i) {
            if (!dst[i]) dst[i] = src[base + i]->clone();
            else         dst[i]->copy(src[base + i]);
        }
    };

    auto& C_left  = jk.C_left();
    auto& C_right = jk.C_right();

    const bool have_K    = jk.get_do_K();
    // K(c_i c_i^T) = (mu i|nu i) for the pair ERIs; without exact exchange in the Fock it is a K-only pass
    // over the active columns.
    const bool supplement_K = pair_k_demand && !have_K;
    if (supplement_K) {
        C_left.clear();
        C_right.clear();
        for (int i = 0; i < n_active; ++i) C_left.push_back(c_act[i]);
        jk.set_do_J(false);
        jk.set_do_K(true);
        jk.compute();
        refresh_vec(reel.K_act, jk.K(), 0, n_active);
        jk.set_do_K(false);
        jk.set_do_J(true);
    }

    C_left.clear();
    C_right.clear();
    if (Ncore > 0) C_left.push_back(C_core);
    for (int i = 0; i < n_active; ++i) C_left.push_back(c_act[i]);
    if (build_c_occ) C_left.push_back(C_occ_cache);
    jk.compute();

    // JK::compute allocates K_/wK_ only for tasks requested via set_do_K/set_do_wK;
    // indexing an unset one reads past the end of an empty vector.
    const bool have_wK = jk.get_do_wK();
    if (functional.is_x_hybrid() && !jk.get_do_K())
        throw PSIEXCEPTION(
            "scf::build_jk_cache: the functional carries exact exchange but the JK object "
            "was not asked for K. set_do_K must follow SuperFunctional::is_x_hybrid().");
    if (functional.is_x_lrc() && !have_wK)
        throw PSIEXCEPTION(
            "scf::build_jk_cache: the functional is range-separated but the JK object was "
            "not asked for wK. set_do_wK must follow SuperFunctional::is_x_lrc().");

    const auto& J_ref = jk.J();
    if (static_cast<int>(J_ref.size()) != n_pushed) {
        timer_off("REKS: build_jk_cache");
        throw PSIEXCEPTION(
            "scf::build_jk_cache: jk J size " +
            std::to_string(J_ref.size()) + " != pushed density count " +
            std::to_string(n_pushed));
    }
    if (Ncore > 0) refresh_one(reel.J_core, J_ref[0]);
    else           reel.J_core.reset();
    refresh_vec(reel.J_act, J_ref, core_present, n_active);
    if (have_K) {
        const auto& K_ref = jk.K();
        if (Ncore > 0) refresh_one(reel.K_core, K_ref[0]);
        else           reel.K_core.reset();
        refresh_vec(reel.K_act, K_ref, core_present, n_active);
    } else {
        // K_core enters the Fock alone, so the supplement leaves it absent.
        reel.K_core.reset();
        if (!supplement_K) reel.K_act.clear();
    }
    if (have_wK) {
        const auto& wK_ref = jk.wK();
        if (Ncore > 0) refresh_one(reel.wK_core, wK_ref[0]);
        else           reel.wK_core.reset();
        refresh_vec(reel.wK_act, wK_ref, core_present, n_active);
    } else {
        reel.wK_core.reset();
        reel.wK_act.clear();
    }

    // Closed-shell C_occ slot: alias jk's J/K/wK outputs (not cloned).
    if (build_c_occ) {
        const double alpha = functional.x_alpha();
        J_at_c_occ = J_ref[c_occ_slot];
        if (have_K) K_at_c_occ = jk.K()[c_occ_slot];
        else        K_at_c_occ.reset();
        if (have_wK) wK_at_c_occ = jk.wK()[c_occ_slot];
        else         wK_at_c_occ.reset();

        // G_out = 2 J - alpha K - beta wK.
        G_out->zero();
        G_out->axpy(2.0, J_at_c_occ);
        if (functional.is_x_hybrid()) {
            G_out->axpy(-alpha, K_at_c_occ);
        }
        if (have_wK) {
            G_out->axpy(-functional.x_beta(), wK_at_c_occ);
        }
    }
    timer_off("REKS: build_jk_cache");
}

FockExchangePolicy fock_exchange_policy(const studio::Reel& reel,
                                        psi::SuperFunctional& functional) {
    const bool is_x_lrc = functional.is_x_lrc();
    FockExchangePolicy pol;
    pol.apply_K  = functional.is_x_hybrid() && !reel.K_act.empty();
    pol.k_scale  = functional.x_alpha();
    pol.apply_wK = is_x_lrc && !reel.wK_act.empty();
    pol.wk_scale = is_x_lrc ? functional.x_beta() : 0.0;
    return pol;
}

void add_coulomb_pattern(SharedMatrix out, const studio::Reel& reel,
                         int Ncore, int n_active, double wJcore, const double* wJ) {
    if (Ncore > 0 && reel.J_core) out->axpy(wJcore, reel.J_core);
    for (int i = 0; i < n_active; ++i)
        if (wJ[i] != 0.0) out->axpy(wJ[i], reel.J_act[i]);
}

void add_exchange_pattern(SharedMatrix out, const studio::Reel& reel,
                          int Ncore, int n_active,
                          const FockExchangePolicy& pol, double wKcore, const double* wK) {
    if (pol.apply_K) {
        if (Ncore > 0 && reel.K_core) out->axpy(-pol.k_scale * wKcore, reel.K_core);
        for (int i = 0; i < n_active; ++i)
            if (wK[i] != 0.0) out->axpy(-pol.k_scale * wK[i], reel.K_act[i]);
    }
    if (pol.apply_wK) {
        if (Ncore > 0 && reel.wK_core) out->axpy(-pol.wk_scale * wKcore, reel.wK_core);
        for (int i = 0; i < n_active; ++i)
            if (wK[i] != 0.0) out->axpy(-pol.wk_scale * wK[i], reel.wK_act[i]);
    }
}

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
    bool                          want_v_acc) {
    timer_on("REKS: build_microstate_xc");
    const int n_active = sa_cassette.n_active_orbitals();
    const int n_L      = static_cast<int>(microstate_Ls.size());
    const bool needs_xc = functional.needs_xc();
    if (!needs_xc || !xc || n_L == 0) {
        timer_off("REKS: build_microstate_xc");
        return;
    }

    // Per-L occupations (Ncore doubly + active) and weights, in L-set order.
    const int n_p = Ncore + n_active;
    std::vector<std::vector<double>> occ_a(n_L, std::vector<double>(n_p, 0.0));
    std::vector<std::vector<double>> occ_b(n_L, std::vector<double>(n_p, 0.0));
    std::vector<double> cl(n_L, 0.0);
    std::vector<SharedMatrix> ar_a(n_L), ar_b(n_L);
    std::vector<std::vector<double>> di_a(n_L), di_b(n_L);
    std::vector<double> tra(n_L, 0.0), trb(n_L, 0.0), exc(n_L, 0.0);
    for (int idx = 0; idx < n_L; ++idx) {
        const int L = microstate_Ls[idx];
        // Sized at first touch and whenever mo_width changes.
        if (!reel.Vxc_arows_a_L[L] || reel.Vxc_arows_a_L[L]->coldim(0) != mo_width) {
            reel.Vxc_arows_a_L[L] = std::make_shared<Matrix>(
                "Vxc_arows_a L=" + std::to_string(L), n_active, mo_width);
            reel.Vxc_arows_b_L[L] = std::make_shared<Matrix>(
                "Vxc_arows_b L=" + std::to_string(L), n_active, mo_width);
        }
        const auto& m = sa_cassette.microstate(L);
        for (int i = 0; i < Ncore; ++i) { occ_a[idx][i] = 1.0; occ_b[idx][i] = 1.0; }
        for (int i = 0; i < n_active; ++i) {
            occ_a[idx][Ncore + i] = static_cast<double>(m.alpha[i]);
            occ_b[idx][Ncore + i] = static_cast<double>(m.beta[i]);
        }
        cl[idx]   = C_L[L];
        ar_a[idx] = reel.Vxc_arows_a_L[L];
        ar_b[idx] = reel.Vxc_arows_b_L[L];
    }

    // Closed-shell microstates (alpha==beta) go through the same unrestricted call.
    xc->compute(
        Ca_full, Ncore, n_active, occ_a, occ_b, cl,
        ar_a, ar_b, di_a, di_b, tra, trb, exc,
        want_v_acc ? reel.V_acc_xc_a : nullptr, want_v_acc ? reel.V_acc_xc_b : nullptr, need_diag,
        mo_lo);

    for (int idx = 0; idx < n_L; ++idx) {
        const int L = microstate_Ls[idx];
        reel.Vxc_diag_a_L[L] = std::move(di_a[idx]);
        reel.Vxc_diag_b_L[L] = std::move(di_b[idx]);
        reel.E_xc_L[L]       = exc[idx];
        reel.tr_D_Vxc_a_L[L] = tra[idx];
        reel.tr_D_Vxc_b_L[L] = trb[idx];
    }

    timer_off("REKS: build_microstate_xc");
}

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
    int                                mo_width) {
    timer_on("REKS: SA Fock AO->MO transform");
    const int n_active = sa_cassette.n_active_orbitals();
    // act_col + i: column of active orbital i inside the window [mo_lo, mo_lo + mo_width).
    const int act_col = Ncore - mo_lo;

    // Fixed Fock-basis MOs snapshot for the AO->MO transform, not the live Ca_.
    const SharedMatrix& Ca = reel.Ca_fock_transform;
    const int N            = Ca->rowspi()[0];
    const int Mn           = Ca->colspi()[0];

    double** Cp = Ca->pointer(0);
    double** Tp = reel.temp_buffer->pointer(0);

    const bool needs_xc = functional.needs_xc();
    const FockExchangePolicy pol = fock_exchange_policy(reel, functional);

    // F_L^s = 1 (H + 2 J_core - x K_core - b wK_core) + sum_i [(na_i + nb_i) J_act[i] - x n^s_i K_act[i]
    // - b n^s_i wK_act[i]], so C^T F_L^s C sums the slot transforms; K slots also give (pk|pk).
    const bool have_K_slots = !reel.K_act.empty();
    const int slot_J  = 1;
    const int slot_K  = slot_J + n_active;
    const int slot_wK = slot_K + (have_K_slots ? n_active : 0);
    const int n_slot  = slot_wK + (pol.apply_wK ? n_active : 0);

    // transform_slot writes every entry it declares.
    const size_t n_diag = need_full_diag ? static_cast<size_t>(n_slot) * Mn : 0;
    if (reel.slot_diag.size() < n_diag) reel.slot_diag.resize(n_diag);
    reel.slot_diag_stride = need_full_diag ? Mn : 0;
    reel.slot_diag_J      = need_full_diag ? slot_J : -1;
    reel.slot_diag_K      = (need_full_diag && have_K_slots) ? slot_K : -1;
    double* const slot_diag = reel.slot_diag.data();

    std::vector<double> slot_arows(static_cast<size_t>(n_slot) * n_active * mo_width);
    std::vector<double> slot_core_trace(need_full_diag ? 0 : static_cast<size_t>(n_slot), 0.0);
    std::vector<double> act_scratch(
        need_full_diag ? 0 : static_cast<size_t>(n_active) * N);

    // Windowed transform: core block of the MO diagonal as tr(D_core M) (zero at Ncore == 0).
    const double* Dcp =
        (!need_full_diag && Ncore > 0) ? D_core->pointer(0)[0] : nullptr;

    // diag_q = sum_mu C[mu,q] (M C)[mu,q] over diag_cols, active rows = C_act^T M C over [mo_lo, mo_lo + mo_width);
    // a zero-weight slot gives only the core/active diagonal (pair ERIs).
    auto transform_slot = [&](int slot, const SharedMatrix& M, int diag_cols, bool want_arows) {
        double** Mp = M->pointer(0);
        double* arows = slot_arows.data() + static_cast<size_t>(slot) * n_active * mo_width;
        if (need_full_diag) {
            C_DGEMM('N', 'N', N, diag_cols, N, 1.0, Mp[0], N, Cp[0], Mn, 0.0, Tp[0], N);
            // mu-outer accumulation keeps both operands and the accumulator at unit stride.
            double* d = slot_diag + static_cast<size_t>(slot) * Mn;
            std::fill(d, d + Mn, 0.0);
            for (int mu = 0; mu < N; ++mu) {
                const double* cmu = Cp[mu];
                const double* tmu = Tp[0] + static_cast<size_t>(mu) * N;
                for (int q = 0; q < diag_cols; ++q) d[q] += cmu[q] * tmu[q];
            }
            if (want_arows)
                C_DGEMM('T', 'N', n_active, mo_width, N, 1.0, &Cp[0][Ncore], Mn,
                        &Tp[0][mo_lo], N, 0.0, arows, mo_width);
            return;
        }
        C_DGEMM('T', 'N', n_active, N, N, 1.0, &Cp[0][Ncore], Mn, Mp[0], N, 0.0,
                act_scratch.data(), N);
        C_DGEMM('N', 'N', n_active, mo_width, N, 1.0, act_scratch.data(), N,
                &Cp[0][mo_lo], Mn, 0.0, arows, mo_width);
        // D_core and M are both symmetric AO matrices, so their elementwise dot is
        // tr(D_core M).
        slot_core_trace[slot] =
            Dcp ? C_DDOT(static_cast<size_t>(N) * N, Mp[0], 1, Dcp, 1) : 0.0;
    };

    // Slot 0 is the only one that has to be accumulated in AO; the rest are the
    // cached JK outputs themselves.
    SharedMatrix F_base = reel.F_alpha_AO_buffer;
    F_base->copy(H);
    if (Ncore > 0 && reel.J_core) F_base->axpy(2.0, reel.J_core);
    if (pol.apply_K  && Ncore > 0 && reel.K_core)  F_base->axpy(-pol.k_scale,  reel.K_core);
    if (pol.apply_wK && Ncore > 0 && reel.wK_core) F_base->axpy(-pol.wk_scale, reel.wK_core);
    transform_slot(0, F_base, Mn, true);
    for (int i = 0; i < n_active; ++i) transform_slot(slot_J + i, reel.J_act[i], Mn, true);
    if (have_K_slots) {
        // pol.apply_K=false: weight zero in fill_coef, so only the diagonal over
        // [0, Ncore+n_active) is needed, with no active rows.
        const int  k_cols  = pol.apply_K ? Mn : Ncore + n_active;
        for (int i = 0; i < n_active; ++i)
            transform_slot(slot_K + i, reel.K_act[i], k_cols, pol.apply_K);
    }
    if (pol.apply_wK)
        for (int i = 0; i < n_active; ++i) transform_slot(slot_wK + i, reel.wK_act[i], Mn, true);

    std::vector<double> w_ab(n_active);

    // Out(2 nc x W) = Coef(2 nc x n_slot) Slot(n_slot x W); row 2j = alpha, 2j + 1 = beta.
    const int    nL     = static_cast<int>(microstate_Ls.size());
    const size_t arow_w = static_cast<size_t>(n_active) * mo_width;
    const int    n_act2 = n_active * n_active;

    // Staging buffers are sized to a chunk of nc microstates, bounded by kStageBytes.
    constexpr size_t kStageBytes = 8u << 20;
    const size_t row_w = arow_w + (need_full_diag ? static_cast<size_t>(Mn) : 0);
    int chunk = static_cast<int>(kStageBytes / (2 * row_w * sizeof(double)));
    chunk = std::max(1, std::min(chunk, nL));

    std::vector<double> coef_all(static_cast<size_t>(2 * chunk) * n_slot, 0.0);
    std::vector<double> arow_stage(static_cast<size_t>(2 * chunk) * arow_w);
    std::vector<double> diag_stage(need_full_diag ? static_cast<size_t>(2 * chunk) * Mn : 0);
    std::vector<std::vector<double>> w_a_L(chunk, std::vector<double>(n_active));
    std::vector<std::vector<double>> w_b_L(chunk, std::vector<double>(n_active));

    reel.fock_aa.reserve(reel.fock_aa.size() + static_cast<size_t>(nL) * n_act2);

    auto fill_coef = [&](double* c, const double* w_s) {
        c[0] = 1.0;
        for (int i = 0; i < n_active; ++i) c[slot_J + i] = w_ab[i];
        if (pol.apply_K)
            for (int i = 0; i < n_active; ++i) c[slot_K + i] = -pol.k_scale * w_s[i];
        if (pol.apply_wK)
            for (int i = 0; i < n_active; ++i) c[slot_wK + i] = -pol.wk_scale * w_s[i];
    };

    // tr(D_L^s F_L^s) = sum_{k<Ncore} diag[k] + sum_i n_i diag[Ncore+i]; windowed, the core sum is the
    // slot-weighted core trace and Fdiag is zero outside the active block.
    auto emit_spin = [&](int row, const double* w_s, double* Fdiag, double* Farows) {
        std::copy(arow_stage.begin() + static_cast<size_t>(row) * arow_w,
                  arow_stage.begin() + static_cast<size_t>(row + 1) * arow_w, Farows);
        double E_core = 0.0;
        if (need_full_diag) {
            std::copy(diag_stage.begin() + static_cast<size_t>(row) * Mn,
                      diag_stage.begin() + static_cast<size_t>(row + 1) * Mn, Fdiag);
            for (int k = 0; k < Ncore; ++k) E_core += Fdiag[k];
        } else {
            std::fill(Fdiag, Fdiag + mo_width, 0.0);
            const double* c = coef_all.data() + static_cast<size_t>(row) * n_slot;
            for (int s = 0; s < n_slot; ++s) E_core += c[s] * slot_core_trace[s];
            for (int i = 0; i < n_active; ++i)
                Fdiag[act_col + i] = Farows[static_cast<size_t>(i) * mo_width + act_col + i];
        }
        double E_F = E_core;
        for (int i = 0; i < n_active; ++i) E_F += w_s[i] * Fdiag[act_col + i];
        return E_F;
    };

    for (int base = 0; base < nL; base += chunk) {
        const int nc = std::min(chunk, nL - base);

        for (int j = 0; j < nc; ++j) {
            const int L = microstate_Ls[base + j];
            // Sized at first touch and whenever mo_width changes.
            if (!reel.F_MO_arows_a_L[L] || reel.F_MO_arows_a_L[L]->coldim(0) != mo_width) {
                reel.F_MO_arows_a_L[L] = std::make_shared<Matrix>(
                    "F_MO_arows_a L=" + std::to_string(L), n_active, mo_width);
                reel.F_MO_arows_b_L[L] = std::make_shared<Matrix>(
                    "F_MO_arows_b L=" + std::to_string(L), n_active, mo_width);
                reel.F_MO_diag_a_L[L].assign(mo_width, 0.0);
                reel.F_MO_diag_b_L[L].assign(mo_width, 0.0);
            }
            // Spin weights (0/1): Coulomb is spin-shared (total density), exchange
            // is spin-resolved (alpha reads na, beta reads nb).
            const auto& m = sa_cassette.microstate(L);
            for (int i = 0; i < n_active; ++i) {
                w_a_L[j][i] = static_cast<double>(m.alpha[i]);
                w_b_L[j][i] = static_cast<double>(m.beta[i]);
                w_ab[i] = w_a_L[j][i] + w_b_L[j][i];
            }
            fill_coef(coef_all.data() + static_cast<size_t>(2 * j) * n_slot, w_a_L[j].data());
            fill_coef(coef_all.data() + static_cast<size_t>(2 * j + 1) * n_slot, w_b_L[j].data());
        }

        C_DGEMM('N', 'N', 2 * nc, static_cast<int>(arow_w), n_slot, 1.0, coef_all.data(),
                n_slot, slot_arows.data(), static_cast<int>(arow_w), 0.0, arow_stage.data(),
                static_cast<int>(arow_w));
        if (need_full_diag)
            C_DGEMM('N', 'N', 2 * nc, Mn, n_slot, 1.0, coef_all.data(), n_slot, slot_diag, Mn,
                    0.0, diag_stage.data(), Mn);

        for (int j = 0; j < nc; ++j) {
            const int L = microstate_Ls[base + j];

            reel.E_Fa_L[L] = emit_spin(2 * j, w_a_L[j].data(), reel.F_MO_diag_a_L[L].data(),
                                       reel.F_MO_arows_a_L[L]->pointer(0)[0]);
            reel.E_Fb_L[L] = emit_spin(2 * j + 1, w_b_L[j].data(), reel.F_MO_diag_b_L[L].data(),
                                       reel.F_MO_arows_b_L[L]->pointer(0)[0]);
            if (needs_xc) {
                reel.E_Fa_L[L] += reel.tr_D_Vxc_a_L[L];
                reel.E_Fb_L[L] += reel.tr_D_Vxc_b_L[L];

                // Grid-direct XC: add the per-microstate MO diagonal + active rows of
                // C^T V_xc^L C into the MO Fock, which already holds the H+J-K-wK part.
                const double* Vxa = reel.Vxc_arows_a_L[L]->pointer(0)[0];
                const double* Vxb = reel.Vxc_arows_b_L[L]->pointer(0)[0];
                if (need_full_diag) {
                    C_DAXPY(Mn, 1.0, reel.Vxc_diag_a_L[L].data(), 1,
                            reel.F_MO_diag_a_L[L].data(), 1);
                    C_DAXPY(Mn, 1.0, reel.Vxc_diag_b_L[L].data(), 1,
                            reel.F_MO_diag_b_L[L].data(), 1);
                } else {
                    // The MO diagonal was not formed on the grid; its active entries are
                    // the diagonal of the active rows.
                    for (int i = 0; i < n_active; ++i) {
                        const size_t d = static_cast<size_t>(i) * mo_width + act_col + i;
                        reel.F_MO_diag_a_L[L][act_col + i] += Vxa[d];
                        reel.F_MO_diag_b_L[L][act_col + i] += Vxb[d];
                    }
                }
                C_DAXPY(arow_w, 1.0, Vxa, 1, reel.F_MO_arows_a_L[L]->pointer(0)[0], 1);
                C_DAXPY(arow_w, 1.0, Vxb, 1, reel.F_MO_arows_b_L[L]->pointer(0)[0], 1);
            }

            // Active block of the spin difference. Appended on first build,
            // overwritten on every later one.
            int aa_row = reel.fock_aa_row[L];
            if (aa_row < 0) {
                aa_row = static_cast<int>(reel.fock_aa.size() / n_act2);
                reel.fock_aa_row[L] = aa_row;
                reel.fock_aa.resize(reel.fock_aa.size() + n_act2);
            }
            double*       aa = reel.fock_aa.data() + static_cast<size_t>(aa_row) * n_act2;
            const double* Fa = reel.F_MO_arows_a_L[L]->pointer(0)[0] + act_col;
            const double* Fb = reel.F_MO_arows_b_L[L]->pointer(0)[0] + act_col;
            for (int p = 0; p < n_active; ++p) {
                const size_t off = static_cast<size_t>(p) * mo_width;
                for (int q = 0; q < n_active; ++q)
                    aa[p * n_active + q] = Fa[off + q] - Fb[off + q];
            }
        }
    }
    timer_off("REKS: SA Fock AO->MO transform");
}

void fill_missing_si_focks(
    studio::Reel&                          reel,
    const studio::Cassette&                sa_cassette,
    const std::vector<studio::Cassette>&   si_cassettes,
    int                                Ncore,
    psi::SuperFunctional&              functional,
    reks::xc::MicrostateXC*            xc,
    const SharedMatrix&                Ca_full,
    const SharedMatrix&                H) {
    // L_to_build: union of the SI cassettes' missing L; per-L outputs are independent.
    std::set<int> built;
    std::vector<int> L_to_build;
    for (const studio::Cassette& si_cassette : si_cassettes) {
        for (int L : si_cassette.microstates()) {
            if (sa_cassette.active(L)) continue;       // built during SCF
            if (!built.insert(L).second) continue;     // built for an earlier SI cassette
            L_to_build.push_back(L);
        }
    }
    if (L_to_build.empty()) return;

    const int n_active = sa_cassette.n_active_orbitals();

    // C_L does not enter the per-L MO outputs; no AO back-projection.
    build_microstate_xc(reel, sa_cassette, L_to_build, Ncore, functional,
                              xc, Ca_full, reel.C_L, /*need_diag=*/false,
                              /*mo_lo=*/Ncore, /*mo_width=*/n_active, /*want_v_acc=*/false);

    // D_core over Ca_full in temp_buffer (reel.D_core belongs to the live Ca_).
    SharedMatrix D_core;
    if (Ncore > 0) {
        const int N  = Ca_full->rowspi()[0];
        const int Mn = Ca_full->colspi()[0];
        double** Cp = Ca_full->pointer(0);
        double** Tp = reel.temp_buffer->pointer(0);
        C_DGEMM('N', 'T', N, N, Ncore, 1.0, Cp[0], Mn, Cp[0], Mn, 0.0, Tp[0], N);
        D_core = reel.temp_buffer;
    }
    build_sa_focks_MO(reel, sa_cassette, L_to_build, Ncore, functional, H,
                      D_core, /*need_full_diag=*/false,
                      /*mo_lo=*/Ncore, /*mo_width=*/n_active);

    // Release per-L rows; the active block stays in reel.fock_aa.
    for (int L : L_to_build) {
        reel.F_MO_arows_a_L[L].reset();
        reel.F_MO_arows_b_L[L].reset();
        reel.Vxc_arows_a_L[L].reset();
        reel.Vxc_arows_b_L[L].reset();
        reel.F_MO_diag_a_L[L] = std::vector<double>{};
        reel.F_MO_diag_b_L[L] = std::vector<double>{};
        reel.Vxc_diag_a_L[L]  = std::vector<double>{};
        reel.Vxc_diag_b_L[L]  = std::vector<double>{};
    }
}

void build_si_microstate_focks(const studio::Reel&                  reel,
                               const studio::Cassette&              cassette,
                               const studio::Cassette&              sa_cassette,
                               double                             E_const,
                               bool                               needs_xc,
                               std::vector<double>&               E_L) {
    const int n_active = cassette.n_active_orbitals();

    for (int Ls : cassette.microstates()) {
        if (sa_cassette.active(Ls)) continue;  // already built in the SCF live set

        const auto& m = cassette.microstate(Ls);
        const int alpha_idx = studio::base_density_index(m, n_active, /*beta=*/false);
        const int beta_idx  = studio::base_density_index(m, n_active, /*beta=*/true);

        E_L[Ls] = microstate_energy(reel, Ls, alpha_idx, beta_idx,
                                    E_const, needs_xc).total;
    }
}

void assemble_F_reks_MO(studio::Reel&                        reel,
                        const studio::Cassette&              sa_cassette,
                        const std::vector<double>&         C_L,
                        const std::vector<int>&            active_mo,
                        int                                Ncore,
                        int                                N,
                        SharedMatrix                       F_reks_MO) {
    F_reks_MO->zero();
    reel.lagrangians.assign(
        static_cast<size_t>(sa_cassette.n_lagrangian_pairs()), 0.0);

    const int n_active       = sa_cassette.n_active_orbitals();
    const int n_pairs        = sa_cassette.n_lagrangian_pairs();
    double**  F              = F_reks_MO->pointer(0);

    // f_eff is independent of L (depends only on reel.fon_state).
    std::vector<double> f_eff;
    compute_all_f_eff(sa_cassette, reel.fon_state, f_eff);
    std::vector<int>    act(n_active);
    int last_active = 0;
    for (int i = 0; i < n_active; ++i) {
        act[i] = active_mo[i];
        if (act[i] > last_active) last_active = act[i];
    }

    std::vector<int>    n_alpha(n_active), n_beta(n_active);
    std::vector<double> W_diag_a(n_active), W_diag_b(n_active);

    // Coupling Fock, sum over L, s: ii 0.5 C_L n_s/f F; ci 0.5 C_L (1 - n_s)/(1 - f) F[i][c]; iv 0.5 C_L n_s/f F[i][v];
    // i<j C_L (n_s[i] - n_s[j]) sign(f_i - f_j) F[i][aj]; 1/f and 1/(1 - f) gated at FON_THRESHOLD.
    for (int L : sa_cassette.microstates()) {
        const auto&   micro = sa_cassette.microstate(L);
        const double  Cl    = C_L[L];
        // Per-L store: full diagonal + active rows (row i = MO Ncore+i).
        // F_MO is symmetric, so core-active F_MO(c,ai) = Farows_a[i][c].
        const double* Fdiag_a  = reel.F_MO_diag_a_L[L].data();
        const double* Fdiag_b  = reel.F_MO_diag_b_L[L].data();
        double**      Farows_a = reel.F_MO_arows_a_L[L]->pointer(0);
        double**      Farows_b = reel.F_MO_arows_b_L[L]->pointer(0);

        for (int i = 0; i < n_active; ++i) {
            n_alpha[i] = micro.alpha[i];
            n_beta[i]  = micro.beta[i];
            W_diag_a[i] = (f_eff[i] > constants::FON_THRESHOLD)
                              ? 0.5 * Cl * n_alpha[i] / f_eff[i] : 0.0;
            W_diag_b[i] = (f_eff[i] > constants::FON_THRESHOLD)
                              ? 0.5 * Cl * n_beta[i]  / f_eff[i] : 0.0;
        }
        for (int i = 0; i < n_active; ++i) {
            int ai = act[i];
            F[ai][ai] += W_diag_a[i] * Fdiag_a[ai] + W_diag_b[i] * Fdiag_b[ai];
        }
        for (int i = 0; i < n_active; ++i) {
            int ai = act[i];
            double omf = 1.0 - f_eff[i];
            double Wca_a = (omf > constants::FON_THRESHOLD)
                              ? 0.5 * Cl * (1 - n_alpha[i]) / omf : 0.0;
            double Wca_b = (omf > constants::FON_THRESHOLD)
                              ? 0.5 * Cl * (1 - n_beta[i])  / omf : 0.0;
            for (int c = 0; c < Ncore; ++c) {
                F[c][ai] += Wca_a * Farows_a[i][c] + Wca_b * Farows_b[i][c];
            }
        }
        for (int i = 0; i < n_active; ++i) {
            int ai = act[i];
            for (int v = last_active + 1; v < N; ++v) {
                F[ai][v] += W_diag_a[i] * Farows_a[i][v] + W_diag_b[i] * Farows_b[i][v];
            }
        }
        for (int i = 0; i < n_active; ++i) {
            for (int j = i + 1; j < n_active; ++j) {
                int ai = act[i];
                int aj = act[j];
                int sign_ij = (f_eff[i] > f_eff[j]) ? 1
                              : ((f_eff[i] < f_eff[j]) ? -1 : 0);
                double Wij_a = Cl * (n_alpha[i] - n_alpha[j]) * sign_ij;
                double Wij_b = Cl * (n_beta[i]  - n_beta[j])  * sign_ij;
                F[ai][aj] += Wij_a * Farows_a[i][aj] + Wij_b * Farows_b[i][aj];
            }
        }

        // Lagrange multiplier per pair (orb_from -> orb_to):
        //   lambda_p += sum_s C_L n_s[orb_from] * F_arows_s[orb_from][act[orb_to]]
        for (int p = 0; p < n_pairs; ++p) {
            const auto& pr = sa_cassette.lagrangian_pair(p);
            int aj = act[pr.orb_to];
            reel.lagrangians[p] +=
                Cl * (n_alpha[pr.orb_from] * Farows_a[pr.orb_from][aj]
                    + n_beta [pr.orb_from] * Farows_b[pr.orb_from][aj]);
        }
    }

    // Uniform-weight bulk (core-core, virt-virt, core-virt) = 0.5 * F_acc_MO.
    {
        double** Facc = reel.F_acc_MO->pointer(0);
        for (int i = 0; i < Ncore; ++i)
            for (int j = i; j < Ncore; ++j)
                F[i][j] += 0.5 * Facc[i][j];
        for (int i = last_active + 1; i < N; ++i)
            for (int j = i; j < N; ++j)
                F[i][j] += 0.5 * Facc[i][j];
        for (int i = 0; i < Ncore; ++i)
            for (int j = last_active + 1; j < N; ++j)
                F[i][j] += 0.5 * Facc[i][j];
    }

    // Symmetrize F (lower triangle from upper) in kTile x kTile blocks.
    constexpr int kTile = 64;
    for (int i0 = 0; i0 < N; i0 += kTile) {
        const int i1 = std::min(i0 + kTile, N);
        for (int j0 = i0; j0 < N; j0 += kTile) {
            const int j1 = std::min(j0 + kTile, N);
            for (int i = i0; i < i1; ++i)
                for (int j = std::max(j0, i + 1); j < j1; ++j) F[j][i] = F[i][j];
        }
    }
}

void compute_gamma_vectors(int                        Ncore,
                           int                        n_act,
                           const studio::Reel&        reel,
                           const studio::Cassette&    sa_cassette,
                           double                     x_alpha,
                           int                        iteration,
                           std::vector<double>&       gamma_aa,
                           std::vector<double>&       gamma_ca) {
    const int Mn = reel.slot_diag_stride;

    // J(p,k) = (pp|kk), K(p,k) = (pk|pk): MO-diagonal slices of the J_act/K_act
    // Fock slots (reks_reel.h), over MO index p and active index k.
    auto J_pk = [&](int p, int k) {
        return reel.slot_diag[static_cast<size_t>(reel.slot_diag_J + k) * Mn + p];
    };
    auto K_pk = [&](int p, int k) {
        return reel.slot_diag[static_cast<size_t>(reel.slot_diag_K + k) * Mn + p];
    };

    // loc < Ncore -> core (1,1); else active occupations.
    auto get_occ = [&](int L, int loc) -> std::pair<int, int> {
        if (loc < Ncore) return {1, 1};
        int act = loc - Ncore;
        return {sa_cassette.microstate(L).alpha[act], sa_cassette.microstate(L).beta[act]};
    };

    // 2e Hessian diagonal (A2e at r2 = r1): gamma_ij = 2 sum_L C_L [(dna^2 + dnb^2)(K_ij - x (J_ij + K_ij)/2)
    // + 2 dna dnb K_ij], x = x_alpha; no range-separated exchange, no XC kernel.
    auto compute_gamma = [&](int i_loc, int j_loc, double J_ij, double K_ij) -> double {
        const double coul = K_ij;
        const double exch = 0.5 * x_alpha * (J_ij + K_ij);
        double h2e = 0.0;

        for (int L : sa_cassette.microstates()) {
            const double CL = reel.C_L[L];
            if (std::abs(CL) < 1e-14) continue;

            auto [na_i, nb_i] = get_occ(L, i_loc);
            auto [na_j, nb_j] = get_occ(L, j_loc);
            const double dna = na_i - na_j;
            const double dnb = nb_i - nb_j;

            h2e += CL * ((dna * dna + dnb * dnb) * (coul - exch) + 2.0 * dna * dnb * coul);
        }

        return 2.0 * h2e;
    };

    int n_rot_aa = n_act * (n_act - 1) / 2;
    gamma_aa.resize(n_rot_aa);
    {
        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                gamma_aa[idx++] = compute_gamma(Ncore + i, Ncore + j,
                                                J_pk(Ncore + i, j), K_pk(Ncore + i, j));
            }
        }
    }

    // Core orbitals are (1, 1) in every microstate: core-active coefficients depend on the active index only.
    gamma_ca.resize(Ncore * n_act);
    if (Ncore > 0) {
        for (int k = 0; k < n_act; ++k) {
            const double g_J = compute_gamma(0, Ncore + k, 1.0, 0.0);
            const double g_K = compute_gamma(0, Ncore + k, 0.0, 1.0);
            for (int c = 0; c < Ncore; ++c)
                gamma_ca[c * n_act + k] = g_J * J_pk(c, k) + g_K * K_pk(c, k);
        }
    }

    if (report::reports(5)) {
        outfile->Printf("  [GVB-DIIS] iter=%d gamma_aa (%d pairs):", iteration, n_rot_aa);
        for (int i = 0; i < n_rot_aa; ++i) outfile->Printf(" %.6f", gamma_aa[i]);
        outfile->Printf("\n");
        if (Ncore > 0) {
            outfile->Printf("  [GVB-DIIS] iter=%d gamma_ca (%d pairs): min=%.6f max=%.6f\n", iteration, Ncore * n_act,
                            *std::min_element(gamma_ca.begin(), gamma_ca.end()),
                            *std::max_element(gamma_ca.begin(), gamma_ca.end()));
        }
    }
}

}  // namespace scf
}  // namespace reks
}  // namespace psi
