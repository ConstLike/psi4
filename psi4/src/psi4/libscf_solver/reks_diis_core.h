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

#ifndef REKS_DIIS_CORE_H
#define REKS_DIIS_CORE_H

#include "reks_math.h"
#include "reks_diis_policies.h"

#include <vector>
#include <utility>
#include <cmath>
#include <cstddef>

namespace psi {
namespace reks {

/// @class DiisCore
/// @brief Payload-agnostic DIIS history and solve core.
///
/// Owns only the error vector, store sequence, and coefficient per slot. The payload is
/// owner-held, keyed by the same physical slot.
class DiisCore {
   public:
    /// Physical (not survivor-local) index written by the last store().
    struct StoreResult {
        int wrote = -1;
    };

    DiisCore() = default;
    explicit DiisCore(int max_vectors) : max_vectors_(max_vectors) {}

    /// Appends error to the ring, setting last_max_error_ = max_k |error[k]|. Once the ring holds
    /// max_vectors_ entries, overwrites the smallest-|q| slot. Precondition: every stored error
    /// has the same length; err_size_ is taken from the last store and indexes every slot.
    StoreResult store(std::vector<double> error) {
        const int n_err = static_cast<int>(error.size());
        last_max_error_ = n_err > 0 ? std::abs(error[C_IDAMAX(n_err, error.data(), 1)]) : 0.0;
        err_size_ = n_err;

        Slot slot;
        slot.error = std::move(error);
        slot.q = 0.0;
        slot.seq = seq_counter_++;

        StoreResult sr;
        if (static_cast<int>(hist_.size()) < max_vectors_) {
            hist_.push_back(std::move(slot));
            sr.wrote = static_cast<int>(hist_.size()) - 1;
        } else {
            int min_idx = 0;
            double min_q = std::abs(hist_[0].q);
            for (int i = 1; i < static_cast<int>(hist_.size()); ++i) {
                double qi = std::abs(hist_[i].q);
                if (qi < min_q) {
                    min_q = qi;
                    min_idx = i;
                }
            }
            hist_[min_idx] = std::move(slot);
            sr.wrote = min_idx;
        }
        last_store_idx_ = sr.wrote;
        if (gram_err_size_ == err_size_)
            refresh_gram_slot(sr.wrote);
        else
            gram_err_size_ = -1;  // err_size_ moved: the whole cache is stale
        return sr;
    }

    /// Erases the slot written by the last store(), shifting the slots above it down. Returns the
    /// erased index, or -1 if no store is tracked.
    int drop_newest() {
        int dropped = -1;
        if (last_store_idx_ >= 0 && last_store_idx_ < static_cast<int>(hist_.size())) {
            const int n_before = static_cast<int>(hist_.size());
            hist_.erase(hist_.begin() + last_store_idx_);
            dropped = last_store_idx_;
            if (gram_err_size_ == err_size_) erase_gram_slot(dropped, n_before);
        }
        last_store_idx_ = -1;
        return dropped;
    }

    void reset() {
        hist_.clear();
        gram_err_size_ = -1;
        last_max_error_ = 0.0;
        last_c_norm_sq_ = 0.0;
        last_store_idx_ = -1;
    }

    int    count() const { return static_cast<int>(hist_.size()); }
    double last_max_error() const { return last_max_error_; }
    double last_c_norm_sq() const { return last_c_norm_sq_; }

    std::vector<double> last_coefficients() const {
        std::vector<double> q;
        q.reserve(hist_.size());
        for (const auto& s : hist_) q.push_back(s.q);
        return q;
    }

    /// Full error-Gram G_ij = <e_i, e_j>.
    /// Served from gram_, which store()/drop_newest() keep current: one store touches a single
    /// slot, so only that slot's row and column are recomputed, and drop_newest() shifts the
    /// cache the same way it shifts hist_. A change of err_size_ invalidates the whole cache.
    std::vector<double> build_gram() const {
        const int n = count();
        if (gram_err_size_ != err_size_) rebuild_gram();
        std::vector<double> G(static_cast<size_t>(n) * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                G[static_cast<size_t>(i) * n + j] = gram_[gram_at(i, j)];
        return G;
    }

    /// Bordered-B coefficients over all live slots, unpruned. False on a singular system.
    ///   min_c c^T G c  s.t.  1^T c = 1,  G = build_gram(); solved as the bordered
    ///   linear system in diis_bordered_coefficients (reks_math.h).
    bool bordered_B(std::vector<double>& c, double& c_norm_sq) const {
        const std::vector<double> G = build_gram();
        return diis_bordered_coefficients(G.data(), count(), c, c_norm_sq);
    }

    /// Constrained least squares over the survivor slots: their error vectors form the columns of
    /// the residual matrix F, anchored on the newest-seq survivor. Returns c indexed by survivor
    /// (size survivors.size()) and |c|^2.
    std::pair<std::vector<double>, double> constrained_lstsq_svd(
        const std::vector<int>& survivors, const std::vector<double>& weights,
        double tikhonov_delta, double svd_rcond) const {
        const int m = static_cast<int>(survivors.size());
        const int es = err_size_;
        fcol_scratch_.resize(static_cast<size_t>(m) * es);  // column-major, col a = survivor a
        for (int a = 0; a < m; ++a) {
            const auto& e = hist_[survivors[a]].error;
            std::copy(e.begin(), e.end(), fcol_scratch_.begin() + static_cast<size_t>(a) * es);
        }
        return solve_packed(m, survivors, weights, tikhonov_delta, svd_rcond);
    }

    /// Re-solve after removing survivor drop_local from the set the last solve packed. Columns are
    /// contiguous, so the removal is a shift of the tail -- the surviving columns keep the bytes
    /// they were packed with, which a repack would have reproduced exactly.
    std::pair<std::vector<double>, double> constrained_lstsq_svd_drop(
        int drop_local, const std::vector<int>& survivors, const std::vector<double>& weights,
        double tikhonov_delta, double svd_rcond) const {
        const int m = static_cast<int>(survivors.size());
        const int es = err_size_;
        const size_t at = static_cast<size_t>(drop_local) * es;
        std::copy(fcol_scratch_.begin() + at + es,
                  fcol_scratch_.begin() + static_cast<size_t>(m + 1) * es,
                  fcol_scratch_.begin() + at);
        fcol_scratch_.resize(static_cast<size_t>(m) * es);
        return solve_packed(m, survivors, weights, tikhonov_delta, svd_rcond);
    }

    /// Solves for the coefficients through the ConditioningChain, then returns combine(c).
    ///   decision cond(E) <= kappa_bypass -> BYPASS: bordered_B over all slots, bit-identical to
    ///                                       the bare bordered solve.
    ///   otherwise                        -> ENGAGED: constrained_lstsq_svd on the survivors.
    /// The |c|^2 <= coeff_norm_max guard applies to EVERY solve. On a trip the ladder escalates,
    /// each rung strictly reducing the survivor set:
    ///   (a) engaged re-solve (angle-filter survivor pruning + ridge when configured) --
    ///       entered from a bypass trip; a solve that was already engaged skips to (b), a
    ///       repeat would be identical;
    ///   (b) transient removal of the oldest-seq SURVIVOR + re-solve, repeating down to m = 2
    ///       (hist_ is never modified);
    ///   (c) raw base map for this iteration (no extrapolation).
    /// A singular bypassed solve escalates through the same ladder from (a).
    /// c is size-n and slot-indexed on both paths, 0 on the slots the chain pruned. Returns a
    /// default-constructed payload (no extrapolation) when n < 2, the chain leaves < 2
    /// survivors, or ladder rung (c) fires.
    template <class CombineFn>
    auto extrapolate(const ConditioningChain& chain, CombineFn&& combine)
        -> decltype(combine(std::declval<const std::vector<double>&>())) {
        using Ret = decltype(combine(std::declval<const std::vector<double>&>()));

        const int n = count();
        if (n < 2) return Ret{};

        const int es = err_size_;
        std::vector<long> seq(n);
        std::vector<double*> E(n);
        for (int i = 0; i < n; ++i) {
            seq[i] = hist_[i].seq;
            E[i] = hist_[i].error.data();
        }

        DiisView view{n, es, E.data(), seq.data()};
        ConditionedSystem sys = chain.apply(view);

        if (static_cast<int>(sys.survivors.size()) < 2) return Ret{};

        std::vector<double> c(n, 0.0);
        double c_norm_sq = 0.0;
        bool over_guard = false;
        if (!sys.engaged) {
            const bool solved = bordered_B(c, c_norm_sq);

            // Bypass is trusted only if DGESV solves the system and |c|^2 <= coeff_norm_max;
            // cond(E) <= kappa_bypass alone does not bound |c|.
            if (!solved || c_norm_sq > chain.coeff_norm_max()) {
                ConditionedSystem esys = chain.apply_engaged(view, sys.cond_before);
                if (static_cast<int>(esys.survivors.size()) < 2) return Ret{};
                auto re = constrained_lstsq_svd(esys.survivors, esys.weights, esys.tikhonov_delta,
                                                chain.svd_rcond());
                std::fill(c.begin(), c.end(), 0.0);
                for (size_t a = 0; a < esys.survivors.size(); ++a)
                    c[esys.survivors[a]] = re.first[a];
                c_norm_sq = re.second;
                sys = std::move(esys);
                over_guard = (c_norm_sq > chain.coeff_norm_max());
            }
        } else {
            // This engaged solve already is rung (a); a guard trip below enters the ladder at (b).
            auto re = constrained_lstsq_svd(sys.survivors, sys.weights, sys.tikhonov_delta,
                                            chain.svd_rcond());
            for (size_t a = 0; a < sys.survivors.size(); ++a) c[sys.survivors[a]] = re.first[a];
            c_norm_sq = re.second;
            over_guard = (c_norm_sq > chain.coeff_norm_max());
        }

        // Rung (b): transiently drop the oldest-seq survivor and re-solve, down to m = 2. The
        // drops are per-solve only -- hist_ keeps every slot.
        if (over_guard) {
            std::vector<int> surv = sys.survivors;
            std::vector<double> wts = sys.weights;
            while (c_norm_sq > chain.coeff_norm_max() && static_cast<int>(surv.size()) > 2) {
                int drop = 0;
                for (int a = 1; a < static_cast<int>(surv.size()); ++a)
                    if (seq[surv[a]] < seq[surv[drop]]) drop = a;
                surv.erase(surv.begin() + drop);
                wts.erase(wts.begin() + drop);
                auto re = constrained_lstsq_svd_drop(drop, surv, wts, sys.tikhonov_delta,
                                                     chain.svd_rcond());
                std::fill(c.begin(), c.end(), 0.0);
                for (size_t a = 0; a < surv.size(); ++a) c[surv[a]] = re.first[a];
                c_norm_sq = re.second;
            }
            if (c_norm_sq > chain.coeff_norm_max()) {
                // Rung (c): no over-long solution is ever accepted -- raw base map this iteration.
                return Ret{};
            }
        }

        last_c_norm_sq_ = c_norm_sq;
        for (int i = 0; i < n; ++i) hist_[i].q = c[i];
        return combine(c);
    }

   private:
    /// Solve over the m columns already packed in fcol_scratch_, anchored on the newest-seq
    /// survivor.
    std::pair<std::vector<double>, double> solve_packed(
        int m, const std::vector<int>& survivors, const std::vector<double>& weights,
        double tikhonov_delta, double svd_rcond) const {
        int anchor = 0;
        long best_seq = (m > 0) ? hist_[survivors[0]].seq : 0;
        for (int a = 0; a < m; ++a)
            if (hist_[survivors[a]].seq > best_seq) {
                best_seq = hist_[survivors[a]].seq;
                anchor = a;
            }
        std::vector<double> c;
        double c_norm_sq = 0.0;
        diis_constrained_lstsq_svd(fcol_scratch_.data(), err_size_, m, weights.data(), anchor,
                                   tikhonov_delta, svd_rcond, c, c_norm_sq);
        return {std::move(c), c_norm_sq};
    }

    /// Row-major index into gram_, whose leading dimension is max_vectors_ so that a slot's
    /// position never moves when count() changes.
    size_t gram_at(int i, int j) const {
        return static_cast<size_t>(i) * max_vectors_ + j;
    }

    /// <e_i, e_j>.
    double gram_dot(int i, int j) const {
        return C_DDOT(err_size_, const_cast<double*>(hist_[i].error.data()), 1,
                      const_cast<double*>(hist_[j].error.data()), 1);
    }

    void rebuild_gram() const {
        const int n = count();
        gram_.assign(static_cast<size_t>(max_vectors_) * max_vectors_, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) {
                const double dot = gram_dot(i, j);
                gram_[gram_at(i, j)] = dot;
                gram_[gram_at(j, i)] = dot;
            }
        gram_err_size_ = err_size_;
    }

    /// Slot w was appended or overwritten: its row and column are the only stale entries.
    void refresh_gram_slot(int w) {
        const int n = count();
        if (static_cast<int>(gram_.size()) != max_vectors_ * max_vectors_) {
            rebuild_gram();
            return;
        }
        for (int i = 0; i < n; ++i) {
            const double dot = (i <= w) ? gram_dot(i, w) : gram_dot(w, i);
            gram_[gram_at(i, w)] = dot;
            gram_[gram_at(w, i)] = dot;
        }
    }

    /// hist_ erased slot w out of n_before entries, shifting the tail down; the cache follows,
    /// carrying its values across verbatim.
    void erase_gram_slot(int w, int n_before) {
        for (int i = 0; i < n_before; ++i) {
            if (i == w) continue;
            const int di = (i > w) ? i - 1 : i;
            for (int j = 0; j < n_before; ++j) {
                if (j == w) continue;
                const int dj = (j > w) ? j - 1 : j;
                gram_[gram_at(di, dj)] = gram_[gram_at(i, j)];
            }
        }
    }

    struct Slot {
        std::vector<double> error;
        double q = 0.0;  ///< coefficient from the last extrapolate(), 0 until then
        long   seq = 0;  ///< store order (larger = newer); not the slot order once eviction starts
    };

    int    max_vectors_ = 10;
    std::vector<Slot> hist_;
    long   seq_counter_ = 0;
    int    last_store_idx_ = -1;
    int    err_size_ = 0;
    double last_max_error_ = 0.0;
    double last_c_norm_sq_ = 0.0;

    /// Cached error-Gram, leading dimension max_vectors_. gram_err_size_ is the err_size_ it
    /// was built at; -1 means invalid. mutable so build_gram() stays const.
    mutable std::vector<double> gram_;
    mutable int gram_err_size_ = -1;

    /// Column-major residual matrix of the last solve; the ladder shortens it in place.
    mutable std::vector<double> fcol_scratch_;
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_DIIS_CORE_H
