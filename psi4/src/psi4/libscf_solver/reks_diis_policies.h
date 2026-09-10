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

#ifndef REKS_DIIS_POLICIES_H
#define REKS_DIIS_POLICIES_H

#include "reks_math.h"
#include "reks_diis_config.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <limits>

namespace psi {
namespace reks {

/// Read-only view of the DIIS history algebra. e_i = the Pulay error vector of slot i.
struct DiisView {
    int n = 0;
    int err_size = 0;
    const double* const* E = nullptr;  ///< n row pointers; E[i] = e_i, err_size long
    const long*   seq = nullptr;       ///< store sequence per slot (larger = newer)
};

/// The conditioned linear system: F = survivor residual matrix, column a = e_{survivors[a]},
/// column-scaled by weights.
struct ConditionedSystem {
    std::vector<int>    survivors;          ///< kept slots, ASCENDING physical order
    std::vector<double> weights;            ///< per-survivor column scale, INDEX-PARALLEL to survivors
    double tikhonov_delta = 0.0;            ///< ridge on the engaged SVD (0 = none)
    bool   engaged = false;                 ///< false = bypass, true = SVD-on-F
    double cond_before = 0.0;               ///< decision cond(E): raw-E SVD at entry
};

/// A chainable conditioning stage. MUST be the identity (return `in` unchanged) when its own
/// criterion is not tripped.
class IConditioningPolicy {
   public:
    virtual ~IConditioningPolicy() = default;
    virtual ConditionedSystem refine(const DiisView& v, ConditionedSystem in) const = 0;
    virtual int priority() const = 0;       ///< fixed chain slot (lower runs first)
};

/// Pollock-Rebholz angle filter with folded Chupin stale-drop. Prunes near-collinear survivors
/// on the error subspace (MGS-on-E), co-pruning the matching weights.
class AngleFilterPolicy : public IConditioningPolicy {
   public:
    AngleFilterPolicy(double angle_tol, double stale_delta)
        : angle_tol_(angle_tol), stale_delta_(stale_delta) {}

    int priority() const override { return 20; }

    ConditionedSystem refine(const DiisView& v, ConditionedSystem in) const override {
        const int m = static_cast<int>(in.survivors.size());
        if (m <= 1) return in;

        // subE/subSeq index only the current survivors; rows are referenced from v, not copied.
        std::vector<const double*> subE(m);
        std::vector<long> subSeq(m);
        for (int a = 0; a < m; ++a) {
            const int s = in.survivors[a];
            subSeq[a] = v.seq[s];
            subE[a] = v.E[s];
        }

        std::vector<int> kept_local;
        std::vector<double> sigmas;
        diis_angle_filter(subE.data(), subSeq.data(), m, v.err_size, angle_tol_, stale_delta_,
                          kept_local, sigmas, q_scratch_);

        ConditionedSystem out = in;
        out.survivors.clear();
        out.weights.clear();
        for (int loc : kept_local) {  // ascending local order -> ascending physical order
            out.survivors.push_back(in.survivors[loc]);
            out.weights.push_back(in.weights[loc]);
        }
        return out;
    }

   private:
    double angle_tol_;
    double stale_delta_;
    /// MGS scratch, persisted across refine() calls to avoid reallocation.
    mutable std::vector<double> q_scratch_;
};

/// Sets the augmented-LS ridge delta = tikhonov_scale on the engaged SVD; scale <= 0 -> identity.
class TikhonovPolicy : public IConditioningPolicy {
   public:
    explicit TikhonovPolicy(double tikhonov_scale) : scale_(tikhonov_scale) {}

    int priority() const override { return 30; }

    ConditionedSystem refine(const DiisView&, ConditionedSystem in) const override {
        if (scale_ > 0.0) in.tikhonov_delta = scale_;
        return in;
    }

   private:
    double scale_;
};

/// Composes the conditioning stages in ascending priority order behind the bypass short-circuit:
/// decision cond(E) <= kappa_bypass returns the identity system with no stage run; above it every
/// stage refines and the system comes back engaged.
class ConditioningChain {
   public:
    ConditioningChain() = default;
    ConditioningChain(double kappa_bypass, double svd_rcond, double coeff_norm_max)
        : kappa_bypass_(kappa_bypass), svd_rcond_(svd_rcond), coeff_norm_max_(coeff_norm_max) {}

    /// Insert a policy, keeping stages_ sorted ascending by priority().
    void add(std::shared_ptr<IConditioningPolicy> p) {
        auto it = std::upper_bound(stages_.begin(), stages_.end(), p,
                                   [](const std::shared_ptr<IConditioningPolicy>& a,
                                      const std::shared_ptr<IConditioningPolicy>& b) {
                                       return a->priority() < b->priority();
                                   });
        stages_.insert(it, std::move(p));
    }

    double svd_rcond() const { return svd_rcond_; }
    double coeff_norm_max() const { return coeff_norm_max_; }

    /// Decision cond(E): sigma_max/sigma_min from a values-only thin SVD (C_DGESDD) of the RAW,
    /// unweighted error matrix (es x n). The row-major n x es view IS the column-major es x n
    /// matrix, copied because DGESDD destroys its input. sigma_min <= 0 or an SVD failure ->
    /// +inf (engaged is the safe side), never NaN.
    double cond_raw(const DiisView& v) const {
        if (v.n <= 0 || v.err_size <= 0) return 0.0;
        const int m = v.err_size;
        const int n = v.n;
        if (cond_scratch_.size() < static_cast<size_t>(n) * m)
            cond_scratch_.resize(static_cast<size_t>(n) * m);
        double* A = cond_scratch_.data();
        for (int i = 0; i < n; ++i)
            std::copy(v.E[i], v.E[i] + m, A + static_cast<size_t>(i) * m);
        const int k = std::min(m, n);
        std::vector<double> s(k, 0.0);
        std::vector<int> iwork(8 * k);
        double wkopt = 0.0;
        int info = C_DGESDD('N', m, n, A, m, s.data(), nullptr, 1, nullptr, 1, &wkopt, -1,
                            iwork.data());
        if (info == 0) {
            const int lwork = std::max(1, static_cast<int>(wkopt));
            std::vector<double> work(lwork);
            info = C_DGESDD('N', m, n, A, m, s.data(), nullptr, 1, nullptr, 1, work.data(),
                            lwork, iwork.data());
        }
        if (info != 0) return std::numeric_limits<double>::infinity();
        const double smax = s.front();
        const double smin = s.back();
        if (!(smin > 0.0)) return std::numeric_limits<double>::infinity();
        return smax / smin;
    }

    ConditionedSystem apply(const DiisView& v) const {
        ConditionedSystem s = identity_system(v);
        s.cond_before = cond_raw(v);
        if (s.cond_before <= kappa_bypass_) {
            s.engaged = false;
            return s;
        }
        return run_stages(v, s.cond_before);
    }

    /// Run the stages irrespective of the bypass short-circuit; the result is always engaged.
    /// cond_before: precomputed decision cond(E), not recomputed here.
    ConditionedSystem apply_engaged(const DiisView& v, double cond_before) const {
        return run_stages(v, cond_before);
    }

   private:
    ConditionedSystem run_stages(const DiisView& v, double cond_before) const {
        ConditionedSystem s = identity_system(v);
        s.cond_before = cond_before;
        for (const auto& st : stages_) s = st->refine(v, std::move(s));
        s.engaged = true;
        return s;
    }

    static ConditionedSystem identity_system(const DiisView& v) {
        ConditionedSystem s;
        s.survivors.resize(v.n);
        s.weights.assign(v.n, 1.0);
        for (int i = 0; i < v.n; ++i) s.survivors[i] = i;
        return s;
    }

    std::vector<std::shared_ptr<IConditioningPolicy>> stages_;
    double kappa_bypass_ = 1e9;
    double svd_rcond_ = 1e-8;
    double coeff_norm_max_ = 1e3;
    /// DGESDD input copy, persisted across cond_raw() calls to avoid reallocation.
    mutable std::vector<double> cond_scratch_;
};

inline ConditioningChain assemble_conditioning_chain(const DiisConfig& cfg) {
    ConditioningChain chain(cfg.kappa_bypass, cfg.svd_rcond, cfg.coeff_norm_max);
    if (cfg.cond_angle_filter)
        chain.add(std::make_shared<AngleFilterPolicy>(cfg.angle_tol, cfg.stale_delta));
    if (cfg.cond_tikhonov)
        chain.add(std::make_shared<TikhonovPolicy>(cfg.tikhonov_scale));
    return chain;
}

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_POLICIES_H
