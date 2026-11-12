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

#ifndef REKS_DIIS_ADAPTERS_H
#define REKS_DIIS_ADAPTERS_H

#include "reks_report_level.h"
#include "reks_diis_step.h"
#include "reks_gvb_diis.h"
#include "reks_orbital_diis.h"
#include "reks_cassette.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace psi {
namespace reks {

/// Inputs and targets for one collapsed DIIS step.
/// Handles into caller storage; valid only for the duration of the call.
struct StepContext {
    // Frame handles: ORB only (CFM is frame-free).
    SharedMatrix*        Ca = nullptr;          ///< live orbitals
    SharedMatrix*        Ca_ref = nullptr;      ///< ORB DIIS frame anchor
    std::vector<double>* cum_kappa = nullptr;   ///< ORB cumulative kappa
    double               orb_base_map_damp = 0.6;  ///< ORB damped-base-map alpha

    // Shared build inputs (F^g / Newton step).
    SharedMatrix F_gen;                                         ///< generalized Fock (nmo x nmo, MO)
    const std::vector<std::vector<double>>* F_MO_diag_a_L = nullptr;
    const std::vector<std::vector<double>>* F_MO_diag_b_L = nullptr;
    const std::vector<double>*              C_L = nullptr;
    const studio::Cassette*                 sa_cassette = nullptr;
    const std::vector<int>*                 active_mo_indices = nullptr;
    int Ncore = 0;
    int nmo = 0;
    // Optional Hessian corrections; null = off.
    const std::vector<double>* gamma_aa = nullptr;
    const std::vector<double>* gamma_ca = nullptr;

    // Apply targets / AO transform.
    SharedMatrix  SC;            ///< AO<-MO map S*Ca
    SharedMatrix  F_reks_MO;     ///< ORB placeholder Fock (MO)
    SharedMatrix  reel_F_reks;   ///< ORB placeholder Fock (AO)
    SharedMatrix* Fa = nullptr;  ///< SCF alpha Fock target
    SharedMatrix* Fb = nullptr;  ///< SCF beta Fock target

    // Diagnostics.
    int    iteration = -1;
    double gorb = 0.0;                            ///< ||g_orb|| at this iterate
    const std::vector<double>* n_fon_dbg = nullptr;  ///< n-FON snapshot; null = none
};

/// Outcome of one collapsed DIIS step.
struct StepResult {
    bool extrap_ok = false;  ///< an extrapolated step was taken (vs the base map)
    bool cap_fired = false;  ///< ORB step cap fired this step
};

/// Frame captured for a best-iterate rewind. Empty for CFM, whose AO history is frame-free.
struct EngineMemento {
    SharedMatrix        Ca_ref;      ///< ORB frame anchor; null for CFM
    std::vector<double> cum_kappa;   ///< ORB cumulative kappa
    /// Anchor-frame id at capture (new on a new anchor, restored by restore()); history valid iff it equals
    /// the live epoch. CFM keeps 0.
    unsigned long       epoch = 0;
};

/// Port isolating the CFM/ORBITAL engine asymmetry: one implementation per formulation.
class IEngineAdapter {
   public:
    virtual ~IEngineAdapter() = default;

    /// One collapsed DIIS step: build the payload + error, store, solve (extrapolate else base
    /// map), and apply the result to the SCF Fock / cumulative rotation.
    virtual StepResult run_step(const StepContext& ctx) = 0;

    /// Abandon extrapolation and re-establish the base map.
    virtual void restart(const StepContext& ctx) = 0;
    /// Best-iterate rewind: reconcile history/frame after Ca_/FONs are rewound to best.
    virtual void on_rewind(const StepContext& ctx) = 0;
    /// Retain-mode rewind: keep the extrapolation subspace across the rewind. Caller
    /// guarantees the frame matches the retained history.
    virtual void on_rewind_retain(const StepContext& ctx) { on_rewind(ctx); }
    /// Capture / restore the frame memento for a best-iterate rewind.
    virtual EngineMemento snapshot(const StepContext& ctx) const = 0;
    virtual void restore(const EngineMemento& m, const StepContext& ctx) = 0;

    /// Identity of the live anchor frame. Retained history is valid iff a memento's
    /// epoch equals this.
    virtual unsigned long epoch() const { return 0; }
};

/// *p, or an empty vector when p is null.
inline const std::vector<double>& deref(const std::vector<double>* p) {
    static const std::vector<double> empty;
    return p ? *p : empty;
}

/// CFM (Muller composite-Fock) engine wrapper. Payload = F^g in AO basis.
class CfmAdapter final : public IEngineAdapter {
   public:
    explicit CfmAdapter(GVBDIISEngine* engine) : engine_(engine) {}

    void reset() { engine_->reset(); }

    /// Build F^g + occupied-row error, store, extrapolate (else the raw F^g base map), and copy
    /// F^opt into the SCF Fock.
    StepResult run_step(const StepContext& ctx) override {
        const int nmo = ctx.nmo;
        auto F_g_MO = engine_->build_multishell_fock(
            ctx.F_gen, *ctx.F_MO_diag_a_L, *ctx.F_MO_diag_b_L, *ctx.C_L, *ctx.sa_cassette,
            *ctx.active_mo_indices, ctx.Ncore, nmo, deref(ctx.gamma_aa), deref(ctx.gamma_ca),
            ctx.iteration);

        // E_MO[o,j] = E_MO[j,o] = F^g[o,j], o occupied, j > o: E_AO = U Yt^T + Yt U^T (DSYR2K, U = SC[:,O],
        // Yt = SC R'^T, R_OO halved) is flattened, then F^g_AO = E_AO + SC diag(F^g) SC^T (DSYRK) in place.
        const int nso   = ctx.SC->rowdim(0);
        const int n_occ = ctx.active_mo_indices->empty() ? ctx.Ncore
                                                         : ctx.active_mo_indices->back() + 1;
        double** Fg  = F_g_MO->pointer(0);
        double** SCp = ctx.SC->pointer(0);

        if (report::reports(5)) {
            // Checks the congruence below: occupied rows symmetric off the diagonal, virt-virt zero, diagonal of
            // one sign (asym, vv, mixed_diag).
            double asym = 0.0, vv = 0.0;
            for (int i = 0; i < n_occ; ++i)
                for (int j = i + 1; j < nmo; ++j)
                    asym = std::max(asym, std::abs(Fg[i][j] - Fg[j][i]));
            for (int i = n_occ; i < nmo; ++i)
                for (int j = i + 1; j < nmo; ++j) vv = std::max(vv, std::abs(Fg[i][j]));
            double dmin = Fg[0][0], dmax = Fg[0][0];
            for (int i = 1; i < nmo; ++i) {
                dmin = std::min(dmin, Fg[i][i]);
                dmax = std::max(dmax, Fg[i][i]);
            }
            outfile->Printf("  [GVB-DIIS] iter=%d Fg asym=%.3e vv=%.3e mixed_diag=%d\n",
                            ctx.iteration, asym, vv, (dmin * dmax < 0.0) ? 1 : 0);
        }

        auto Rp = std::make_shared<Matrix>("F^g occupied rows", n_occ, nmo);
        double** Rp_p = Rp->pointer(0);
        for (int o = 0; o < n_occ; ++o) {
            for (int q = 0; q < n_occ; ++q) Rp_p[o][q] = 0.5 * Fg[o][q];
            for (int q = n_occ; q < nmo; ++q) Rp_p[o][q] = Fg[o][q];
            Rp_p[o][o] = 0.0;
        }
        auto Yt = linalg::doublet(ctx.SC, Rp, false, true);   // nso x n_occ
        auto F_g_AO = std::make_shared<Matrix>("F^g AO", nso, nso);
        double** Fao = F_g_AO->pointer(0);
        C_DSYR2K('U', 'N', nso, n_occ, 1.0, SCp[0], nmo, Yt->pointer(0)[0], n_occ, 0.0,
                 Fao[0], nso);
        auto error = GVBDIISEngine::flatten_upper_triangle(F_g_AO, nso);

        // SC diag(F^g) SC^T as one-signed DSYRK on columns scaled by sqrt|F^g[i,i]|, kb columns per chunk.
        const double sgn = (Fg[nmo - 1][nmo - 1] < 0.0) ? -1.0 : 1.0;
        std::vector<double> sd(nmo);
        for (int q = 0; q < nmo; ++q) sd[q] = std::sqrt(sgn * Fg[q][q]);
        const int kb = std::min(nmo, 256);
        if (static_cast<int>(scaled_.size()) < nso * kb) scaled_.resize(static_cast<size_t>(nso) * kb);
        double* sp = scaled_.data();
        for (int q0 = 0; q0 < nmo; q0 += kb) {
            const int k = std::min(kb, nmo - q0);
            for (int mu = 0; mu < nso; ++mu)
                for (int q = 0; q < k; ++q)
                    sp[static_cast<size_t>(mu) * kb + q] = SCp[mu][q0 + q] * sd[q0 + q];
            C_DSYRK('U', 'N', nso, k, sgn, sp, kb, 1.0, Fao[0], nso);
        }

        // Mirror to full storage in one tiled pass: the strided lower-triangle writes
        // stay inside a tile.
        {
            constexpr int T = 64;
            for (int i0 = 0; i0 < nso; i0 += T)
                for (int j0 = i0; j0 < nso; j0 += T) {
                    const int i1 = std::min(i0 + T, nso);
                    const int j1 = std::min(j0 + T, nso);
                    for (int i = i0; i < i1; ++i)
                        for (int j = std::max(j0, i + 1); j < j1; ++j) Fao[j][i] = Fao[i][j];
                }
        }

        if (report::reports(4)) print_error_debug(error, nso, ctx);
        engine_->store(F_g_AO, std::move(error));
        if (report::reports(5)) print_gradient_breakdown(ctx);
        if (report::reports(4))
            outfile->Printf("  [GVB-DIIS] iter=%d count=%d max_err=%.3e gorb=%.3e\n", ctx.iteration,
                            engine_->count(), engine_->last_max_error(), ctx.gorb);

        SharedMatrix F_extrap = engine_->extrapolate();
        const bool extrap_ok = (F_extrap != nullptr);

        apply(extrap_ok ? F_extrap : F_g_AO, ctx);

        if (report::reports(4)) {
            outfile->Printf("  [GVB-DIIS] iter=%d extrapolated=%s |c|^2=%.3e\n", ctx.iteration,
                            extrap_ok ? "yes" : "no", engine_->last_c_norm_sq());
            if (extrap_ok) {
                auto q = engine_->last_coefficients();
                outfile->Printf("  [GVB-DIIS] iter=%d DIIS coefficients:", ctx.iteration);
                for (size_t k = 0; k < q.size(); ++k) outfile->Printf(" %.4f", q[k]);
                outfile->Printf("\n");
            }
        }
        return {extrap_ok, false};
    }

    /// AO error is frame-independent: a restart is a plain history clear.
    void restart(const StepContext&) override { engine_->reset(); }
    /// AO payloads are frame-independent: only the diverged (newest) entry needs dropping.
    void on_rewind(const StepContext&) override { engine_->drop_newest(); }
    /// Frame-free: no memento.
    EngineMemento snapshot(const StepContext&) const override { return {}; }
    void restore(const EngineMemento&, const StepContext&) override {}

   private:
    void apply(SharedMatrix F_opt, const StepContext& ctx) {
        (*ctx.Fa)->copy(F_opt);
        (*ctx.Fb)->copy(*ctx.Fa);
    }

    /// Max/RMS error location (decoded with the packing dimension dim, AO/SO) and n-FON echo.
    void print_error_debug(const std::vector<double>& error, int dim, const StepContext& ctx) const {
        double err_max = 0.0, err_rms = 0.0;
        int err_max_idx = 0;
        for (int k = 0; k < static_cast<int>(error.size()); ++k) {
            double ae = std::abs(error[k]);
            err_rms += error[k] * error[k];
            if (ae > err_max) { err_max = ae; err_max_idx = k; }
        }
        err_rms = std::sqrt(err_rms / error.size());
        // Inverse of the ii-outer strict-upper-triangle packing: rows before ii hold
        // ii*(dim-1) - ii*(ii-1)/2 pairs; ii is the smaller root of ii^2 - (2*dim-1)*ii + 2*k.
        const double half_span = 2.0 * dim - 1.0;
        const int ei = static_cast<int>(
            0.5 * (half_span - std::sqrt(half_span * half_span - 8.0 * err_max_idx)));
        const int ej = ei + 1 + (err_max_idx - (ei * (dim - 1) - ei * (ei - 1) / 2));
        outfile->Printf("  [GVB-DIIS] iter=%d error: max=%.3e at (%d,%d), rms=%.3e, size=%d\n",
                        ctx.iteration, err_max, ei, ej, err_rms, static_cast<int>(error.size()));
        if (ctx.n_fon_dbg) {
            const auto& fv = *ctx.n_fon_dbg;
            outfile->Printf("  [GVB-DIIS] iter=%d FON=(", ctx.iteration);
            for (size_t p = 0; p < fv.size(); ++p)
                outfile->Printf("%.6f%s", fv[p], p + 1 < fv.size() ? ", " : "");
            outfile->Printf(")\n");
        }
    }

    /// ||g_orb|| per rotation block (aa, ca, cv, av), with no aa pair frozen.
    void print_gradient_breakdown(const StepContext& ctx) const {
        GradientBlocks blocks;
        orbital_gradient_norm(ctx.F_gen->pointer(0), *ctx.active_mo_indices, ctx.Ncore, ctx.nmo, {},
                              &blocks);
        outfile->Printf("  [GVB-DIIS] iter=%d gorb_components: aa=%.3e ca=%.3e cv=%.3e av=%.3e\n",
                        ctx.iteration, blocks.aa, blocks.ca, blocks.cv, blocks.av);
    }

    GVBDIISEngine* engine_;
    /// sqrt-scaled SC column panel, persisted across run_step() calls to avoid reallocation.
    std::vector<double> scaled_;
};

/// ORBITAL (Ionova-Carter rotation DIIS + Sethio r-GDIIS) engine wrapper. Payload = cumulative
/// kappa.
class OrbAdapter final : public IEngineAdapter {
   public:
    explicit OrbAdapter(OrbitalDIISEngine* engine) : engine_(engine) {}

    /// Non-redundant rotation-parameter count for the packed kappa vector.
    static int n_nonred(int Ncore, int n_act, int n_virt) {
        return OrbitalDIISEngine::n_nonred(Ncore, n_act, n_virt);
    }

    void apply_rotation(const SharedMatrix& Ca_ref, const SharedMatrix& Ca_out,
                        const std::vector<double>& kappa,
                        const std::vector<int>& active_mo_indices,
                        int Ncore, int nmo) const {
        engine_->apply_rotation(Ca_ref, Ca_out, kappa, active_mo_indices, Ncore, nmo);
    }

    void reset() { engine_->reset(); }

    /// Newton step + r-GDIIS cumulative trial, store, extrapolate (else the base map), cap, and
    /// commit the cumulative kappa. Sets the SCF placeholder Fock; rotates nothing.
    StepResult run_step(const StepContext& ctx) override {
        auto error = engine_->compute_newton_rotation(
            ctx.F_gen, *ctx.F_MO_diag_a_L, *ctx.F_MO_diag_b_L, *ctx.C_L, *ctx.sa_cassette,
            *ctx.active_mo_indices, ctx.Ncore, ctx.nmo, deref(ctx.gamma_aa), deref(ctx.gamma_ca),
            ctx.iteration);

        // r-GDIIS: cum_trial = previous cumulative + Newton (zeroed if the frame resized).
        std::vector<double> cum_trial = *ctx.cum_kappa;
        const int n_kappa = static_cast<int>(error.size());
        if (static_cast<int>(cum_trial.size()) != n_kappa) cum_trial.assign(n_kappa, 0.0);
        for (int k = 0; k < n_kappa; ++k) cum_trial[k] += error[k];

        engine_->store(std::move(cum_trial), error);
        std::vector<double> kappa_opt = engine_->extrapolate();
        const bool extrap_ok = !kappa_opt.empty();
        if (!extrap_ok) kappa_opt = base_map(error, ctx);

        const bool capped = engine_->apply_step_cap(kappa_opt);
        if (capped && report::reports(4)) {
            double max_abs = 0.0;
            for (double v : kappa_opt) max_abs = std::max(max_abs, std::abs(v));
            outfile->Printf("  [ORB-DIIS] iter=%d step cap fired (max=%.4f)\n", ctx.iteration, max_abs);
        }

        *ctx.cum_kappa = std::move(kappa_opt);

        if (report::reports(4)) {
            outfile->Printf(
                "  [ORB-DIIS] iter=%d count=%d err_max=%.3e gorb=%.3e extrap=%s |c|^2=%.3e\n",
                ctx.iteration, engine_->count(), engine_->last_max_error(), ctx.gorb,
                extrap_ok ? "yes" : "no", engine_->last_c_norm_sq());
            if (extrap_ok) {
                auto q = engine_->last_coefficients();
                outfile->Printf("  [ORB-DIIS] iter=%d coeffs:", ctx.iteration);
                for (size_t k = 0; k < q.size(); ++k) outfile->Printf(" %.4f", q[k]);
                outfile->Printf("\n");
            }
        }

        auto F_reks_AO = linalg::triplet(ctx.SC, ctx.F_reks_MO, ctx.SC, false, false, true);
        ctx.reel_F_reks->copy(F_reks_AO);
        (*ctx.Fa)->copy(ctx.reel_F_reks);
        (*ctx.Fb)->copy(*ctx.Fa);
        return {extrap_ok, capped};
    }

    /// cum_kappa is frame-relative: restart clears history AND re-anchors the frame to the current
    /// orbitals with zero cumulative rotation. The next base map is then damped.
    void restart(const StepContext& ctx) override {
        engine_->reset();
        reanchor(ctx);
        damp_next_ = true;
    }

    /// The frame-relative history is stale once Ca_ is rewound. History only; restore()
    /// sets the frame.
    void on_rewind(const StepContext&) override { engine_->reset(); }
    /// Drops only the newest vector; the next base map is damped.
    void on_rewind_retain(const StepContext&) override {
        engine_->drop_newest();
        damp_next_ = true;
    }

    /// Frame memento = the Ca_ref anchor + cumulative kappa + the frame identity.
    EngineMemento snapshot(const StepContext& ctx) const override {
        EngineMemento m;
        if (ctx.Ca_ref && *ctx.Ca_ref) m.Ca_ref = (*ctx.Ca_ref)->clone();
        if (ctx.cum_kappa) m.cum_kappa = *ctx.cum_kappa;
        m.epoch = epoch_;
        return m;
    }
    /// Restores the frame the memento was taken in, epoch included: returning to a frame
    /// is not the establishment of a new one.
    void restore(const EngineMemento& m, const StepContext& ctx) override {
        if (ctx.Ca_ref && m.Ca_ref) {
            *ctx.Ca_ref = m.Ca_ref->clone();
            (*ctx.Ca_ref)->set_name("Ca_ref (orbital DIIS)");
        }
        if (ctx.cum_kappa) *ctx.cum_kappa = m.cum_kappa;
        epoch_ = m.epoch;
    }

    /// Establishes a new anchor frame at the current orbitals.
    void anchor_here(const StepContext& ctx) { reanchor(ctx); }

    unsigned long epoch() const override { return epoch_; }

   private:
    /// orb_base_map_damp * Newton after a restart or retain-rewind, else the r-GDIIS trial cum_kappa + Newton.
    std::vector<double> base_map(const std::vector<double>& nw, const StepContext& ctx) {
        std::vector<double> kappa(nw.size(), 0.0);
        if (damp_next_) {
            for (size_t k = 0; k < nw.size(); ++k) kappa[k] = ctx.orb_base_map_damp * nw[k];
            damp_next_ = false;
        } else {
            const bool size_ok = ctx.cum_kappa && ctx.cum_kappa->size() == nw.size();
            for (size_t k = 0; k < nw.size(); ++k)
                kappa[k] = (size_ok ? (*ctx.cum_kappa)[k] : 0.0) + nw[k];
        }
        return kappa;
    }

    /// Sole mint site for anchor frames: clones the live orbitals into the anchor, zeroes
    /// the cumulative rotation and issues a fresh epoch.
    void reanchor(const StepContext& ctx) {
        if (ctx.Ca && *ctx.Ca && ctx.Ca_ref) {
            *ctx.Ca_ref = (*ctx.Ca)->clone();
            (*ctx.Ca_ref)->set_name("Ca_ref (orbital DIIS)");
        }
        if (ctx.cum_kappa) std::fill(ctx.cum_kappa->begin(), ctx.cum_kappa->end(), 0.0);
        ++epoch_;
    }

    OrbitalDIISEngine* engine_;
    bool damp_next_ = false;  ///< next base map is damped (set by restart or retain-rewind)
    unsigned long epoch_ = 0; ///< live anchor-frame identity; owned here, minted in reanchor
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_ADAPTERS_H
