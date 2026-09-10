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

#ifndef REKS_MINRES_H
#define REKS_MINRES_H

/// @file reks_minres.h
/// @brief MINRES (Paige-Saunders) with full Lanczos reorthogonalization, SPD diagonal
/// preconditioner and a Krylov trust region. Pure math: the operator enters only through
/// its action, no REKS types. Dense subproblems go through the LAPACK wrappers of libqt.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include "psi4/libqt/qt.h"

namespace psi {
namespace reks {

/// Outcome of one MINRES trust-region solve.
struct MinresResult {
    std::vector<double> k;     ///< solution of H k = -g restricted to the Krylov subspace
    double resid_2norm;        ///< ||g + H k||_2 achieved
    int n_matvec;              ///< operator applications spent
    bool hit_radius;           ///< the constraint ||k|| <= Delta is active
    double step_norm;          ///< ||k||_2
};

namespace minres_detail {

/// Rounds an LAPACK lwork = -1 query probe (work[0]) up to an integer workspace size.
inline int query_lwork(double probe) { return static_cast<int>(probe) + 1; }

/// LIWORK for DGELSD: 3*MINMN*NLVL + 11*MINMN with NLVL the divide-and-conquer depth for
/// SMLSIZ = 25 (dgelsd.f); the LAPACK workspace query returns LWORK only.
inline int gelsd_liwork(int m) {
    const int smlsiz = 25;
    const int nlvl = std::max(0, static_cast<int>(std::log2(static_cast<double>(m) /
                                                            static_cast<double>(smlsiz + 1))) + 1);
    return 3 * m * nlvl + 11 * m + 1;
}

/// Symmetric eigendecomposition of a column-major m x m matrix; on exit a holds the
/// eigenvectors as columns and w the eigenvalues in ascending order.
inline void sym_eig(int m, double* a, double* w, char jobz, double* work, int lwork) {
    int info = C_DSYEV(jobz, 'U', m, a, m, w, work, lwork);
    if (info != 0) throw std::runtime_error("reks::minres_tr: C_DSYEV failed");
}

/// Minimum-norm least-squares solution of A y = b for column-major A (n x m), rhs of
/// length n. A and b are overwritten; y lands in b[0..m-1]. s takes the m singular values.
inline void min_norm_lstsq(int n, int m, double* a, double* b, double* s, double* work, int lwork,
                           int* iwork) {
    int rank = 0;
    int info = C_DGELSD(n, m, 1, a, n, b, n, s, -1.0, &rank, work, lwork, iwork);
    if (info != 0) throw std::runtime_error("reks::minres_tr: C_DGELSD failed");
}

/// Thin SVD A = U diag(s) Vt of a column-major A (n x m, n >= m). A is overwritten.
inline void thin_svd(int n, int m, double* a, double* u, double* s, double* vt, double* work,
                     int lwork, int* iwork) {
    int info = C_DGESDD('S', n, m, a, n, s, u, n, vt, m, work, lwork, iwork);
    if (info != 0) throw std::runtime_error("reks::minres_tr: C_DGESDD failed");
}

/// ||z(mu)||_2 for z(mu) = V diag(s/(s^2+mu)) c, the ridge family of the thin SVD.
inline double ridge_norm(int m, const double* s, const double* c, double mu) {
    double acc = 0.0;
    for (int i = 0; i < m; ++i) {
        const double t = s[i] * c[i] / (s[i] * s[i] + mu);
        acc += t * t;
    }
    return std::sqrt(acc);
}

}  // namespace minres_detail

/// MINRES (Paige, C. C.; Saunders, M. A. SIAM J. Numer. Anal. 1975, 12, 617) with full
/// reorthogonalization of the Lanczos basis, SPD preconditioner and a trust region on the
/// step norm.
///
/// Solves H k = -g with H symmetric indefinite, available only as an action.
///
/// Scheme. With M = diag(1/Minv_diag) SPD and S = M^{-1/2} = diag(sqrt(Minv_diag)),
///
///     H_hat = S H S,   g_hat = S g,   k = S x,
///
/// Lanczos runs on H_hat and produces a 2-orthonormal basis V; the step basis
/// V_k = S V is M-orthonormal, so ||k|| = ||V_k y|| = sqrt(y^T G y) with G = V_k^T V_k
/// and NOT ||y||. Full reorthogonalization (twice-modified Gram-Schmidt against the whole
/// stored basis) replaces the three-term recursion, which loses the basis under matvec
/// noise. The preconditioner must be SPD: with an indefinite M the M-inner product is not
/// an inner product and no equivalent symmetric system exists (Wathen, A. J. Acta Numer.
/// 2015, 24, 329, Sec. 5).
///
/// Iterate and exit criterion live in the 2-norm, computed from the stored R = H V_k:
///
///     min_y ||g + R y||_2   subject to   y^T G y <= Delta^2.
///
/// The radius is imposed by solving that constrained least-squares problem inside the
/// subspace already built, not by truncating the Krylov iteration.
///
/// With G = W diag(gamma) W^T and X = W diag(gamma^{-1/2}) the substitution y = X z makes
/// the constraint isotropic, y^T G y = z^T z; the thin SVD A = R X = U diag(s) V^T then
/// gives the whole regularization path in closed form,
///
///     z(mu) = V diag(s/(s^2+mu)) U^T (-g),   ||z(mu)|| strictly decreasing in mu,
///
/// so mu is found by bisection on ||z(mu)|| = Delta.
///
/// @param n         problem dimension
/// @param g         gradient, length n
/// @param Hv        operator action y = H x
/// @param Minv_diag diagonal of M^-1, strictly positive, length n
/// @param eta       target ||g + H k||_2 <= eta * ||g||_2
/// @param delta     trust radius; +infinity requests the unconstrained solve
/// @param max_inner cap on Lanczos steps, one operator application each
inline MinresResult minres_tr(int n, const double* g, const std::function<void(const double*, double*)>& Hv,
                              const double* Minv_diag, double eta, double delta, int max_inner) {
    if (n <= 0) throw std::invalid_argument("reks::minres_tr: n must be positive");
    if (g == nullptr) throw std::invalid_argument("reks::minres_tr: g is null");
    if (Minv_diag == nullptr) throw std::invalid_argument("reks::minres_tr: Minv_diag is null");
    if (!Hv) throw std::invalid_argument("reks::minres_tr: Hv is empty");
    if (max_inner <= 0) throw std::invalid_argument("reks::minres_tr: max_inner must be positive");
    if (!(eta > 0.0)) throw std::invalid_argument("reks::minres_tr: eta must be positive");
    if (!(delta > 0.0)) throw std::invalid_argument("reks::minres_tr: delta must be positive");
    for (int i = 0; i < n; ++i) {
        if (!(Minv_diag[i] > 0.0) || !std::isfinite(Minv_diag[i]))
            throw std::invalid_argument("reks::minres_tr: Minv_diag must be finite and strictly positive");
    }

    const std::size_t ns = static_cast<std::size_t>(n);
    const bool bounded = std::isfinite(delta);

    MinresResult res;
    res.k.assign(ns, 0.0);
    res.resid_2norm = C_DNRM2(ns, const_cast<double*>(g), 1);
    res.n_matvec = 0;
    res.hit_radius = false;
    res.step_norm = 0.0;

    const double gnorm = res.resid_2norm;
    if (gnorm == 0.0) return res;

    const int mmax = std::min(max_inner, n);
    const std::size_t mm = static_cast<std::size_t>(mmax);

    // Half preconditioner S = M^{-1/2}; Lanczos runs on S H S.
    std::vector<double> sq(ns);
    for (int i = 0; i < n; ++i) sq[i] = std::sqrt(Minv_diag[i]);

    std::vector<double> V(ns * mm, 0.0);   // Lanczos basis of S H S, 2-orthonormal
    std::vector<double> HP(ns * mm, 0.0);  // R = H V_k, columns H (S v_j)
    std::vector<double> G(mm * mm, 0.0);   // Gram matrix V_k^T V_k, the step metric
    std::vector<double> alpha(mm, 0.0);
    std::vector<double> beta(mm, 0.0);

    std::vector<double> w(ns), pj(ns), dv(ns), rvec(ns), kvec(ns);
    std::vector<double> Amat(ns * mm), Umat(ns * mm), svals(mm), Vt(mm * mm);
    std::vector<double> Gw(mm * mm), Xmat(mm * mm), gamma(mm);
    std::vector<double> yvec(mm), zvec(mm), tvec(mm), cvec(mm);

    // LAPACK scratch, sized once by a workspace query at the largest subproblem m = mmax. The
    // required workspace is nondecreasing in m, so this covers every Lanczos step.
    std::vector<int> iwork(
        std::max(static_cast<std::size_t>(minres_detail::gelsd_liwork(mmax)), 8u * mm));
    int lwork = 1;
    {
        double probe = 0.0;
        int rank = 0;
        C_DSYEV('V', 'U', mmax, Gw.data(), mmax, gamma.data(), &probe, -1);
        lwork = std::max(lwork, minres_detail::query_lwork(probe));
        C_DGESDD('S', n, mmax, Amat.data(), n, svals.data(), Umat.data(), n, Vt.data(), mmax,
                 &probe, -1, iwork.data());
        lwork = std::max(lwork, minres_detail::query_lwork(probe));
        C_DGELSD(n, mmax, 1, Amat.data(), n, rvec.data(), n, svals.data(), -1.0, &rank, &probe, -1,
                 iwork.data());
        lwork = std::max(lwork, minres_detail::query_lwork(probe));
    }
    std::vector<double> work(static_cast<std::size_t>(lwork));

    for (int i = 0; i < n; ++i) w[i] = sq[i] * g[i];
    const double beta1 = C_DNRM2(ns, w.data(), 1);
    C_DCOPY(ns, w.data(), 1, V.data(), 1);
    C_DSCAL(ns, 1.0 / beta1, V.data(), 1);

    double anorm_est = 0.0;
    const double eps = std::numeric_limits<double>::epsilon();

    for (int j = 0; j < mmax; ++j) {
        double* vj = &V[static_cast<std::size_t>(j) * ns];
        double* hpj = &HP[static_cast<std::size_t>(j) * ns];

        // Lanczos step on H_hat = S H S: p_j = S v_j, w = H_hat v_j = S (H p_j) (H p_j stored
        // in HP for the 2-norm residual), alpha_j = v_j^T w, w <- w - alpha_j v_j - beta_{j-1} v_{j-1}.
        for (int i = 0; i < n; ++i) pj[i] = sq[i] * vj[i];
        Hv(pj.data(), hpj);
        ++res.n_matvec;
        for (int i = 0; i < n; ++i) w[i] = sq[i] * hpj[i];

        alpha[j] = C_DDOT(ns, vj, 1, w.data(), 1);
        C_DAXPY(ns, -alpha[j], vj, 1, w.data(), 1);
        if (j > 0) C_DAXPY(ns, -beta[j - 1], &V[static_cast<std::size_t>(j - 1) * ns], 1, w.data(), 1);

        // Full reorthogonalization, two passes: one sweep leaves O(sqrt(eps)) components
        // when the recursion coefficients are inaccurate.
        for (int pass = 0; pass < 2; ++pass) {
            for (int i = 0; i <= j; ++i) {
                const double* vi = &V[static_cast<std::size_t>(i) * ns];
                const double c = C_DDOT(ns, vi, 1, w.data(), 1);
                C_DAXPY(ns, -c, vi, 1, w.data(), 1);
            }
        }
        const double bnew = C_DNRM2(ns, w.data(), 1);

        // Gershgorin bound on ||H_hat||: max row sum of the tridiagonal Lanczos matrix
        // (diagonal alpha_j, off-diagonals beta_{j-1}, beta_j).
        const double row = std::fabs(alpha[j]) + (j > 0 ? beta[j - 1] : 0.0) + bnew;
        anorm_est = std::max(anorm_est, row);

        // Column j of the step metric G_ij = (S v_i)^T (S v_j) = v_i^T M^-1 v_j. Column-major V
        // (n x mm, ld = ns) is row-major (mm x n, ld = ns) for the wrapper, and column j of G is
        // contiguous, so the first j+1 rows of V hit it in one pass.
        for (int i = 0; i < n; ++i) dv[i] = Minv_diag[i] * vj[i];
        C_DGEMV('N', j + 1, n, 1.0, V.data(), n, dv.data(), 1, 0.0,
                &G[static_cast<std::size_t>(j) * mm], 1);
        for (int i = 0; i < j; ++i)
            G[static_cast<std::size_t>(i) * mm + static_cast<std::size_t>(j)] =
                G[static_cast<std::size_t>(j) * mm + static_cast<std::size_t>(i)];

        const int m = j + 1;
        const std::size_t msz = static_cast<std::size_t>(m);

        // Constraint metric: y = X z with X^T G X = I turns y^T G y <= Delta^2 into
        // ||z|| <= Delta. Without a radius the metric is never needed.
        if (bounded) {
            for (std::size_t col = 0; col < msz; ++col)
                for (std::size_t row_i = 0; row_i < msz; ++row_i)
                    Gw[col * msz + row_i] = G[col * mm + row_i];
            minres_detail::sym_eig(m, Gw.data(), gamma.data(), 'V', work.data(), lwork);
            for (std::size_t col = 0; col < msz; ++col) {
                const double scale = 1.0 / std::sqrt(gamma[col]);
                for (std::size_t row_i = 0; row_i < msz; ++row_i)
                    Xmat[col * msz + row_i] = Gw[col * msz + row_i] * scale;
            }
            // A = HP X in column-major is A_rm = X_rm HP_rm for the row-major wrapper, with
            // A_rm and HP_rm (m x n, ld = ns) and X_rm (m x m, ld = m).
            C_DGEMM('N', 'N', m, n, m, 1.0, Xmat.data(), m, HP.data(), n, 0.0, Amat.data(), n);
        } else {
            C_DCOPY(msz * ns, HP.data(), 1, Amat.data(), 1);
        }

        double mu = 0.0;
        bool active = false;
        if (bounded) {
            minres_detail::thin_svd(n, m, Amat.data(), Umat.data(), svals.data(), Vt.data(),
                                    work.data(), lwork, iwork.data());
            // c = -U^T g: column-major U (n x m, ld = ns) is row-major (m x n, ld = ns).
            C_DGEMV('N', m, n, -1.0, Umat.data(), n, const_cast<double*>(g), 1, 0.0, cvec.data(), 1);

            const double stol = static_cast<double>(std::max(n, m)) * eps * svals[0];
            double nz = 0.0;
            for (std::size_t i = 0; i < msz; ++i) {
                const double t = (svals[i] > stol) ? cvec[i] / svals[i] : 0.0;
                nz += t * t;
            }
            nz = std::sqrt(nz);

            if (nz > delta) {
                // ||z(mu)|| decreases strictly from nz at mu = 0; ||z(mu)|| <= ||s.*c||/mu
                // brackets the root from above.
                active = true;
                double hi = 0.0;
                for (std::size_t i = 0; i < msz; ++i) hi += (svals[i] * cvec[i]) * (svals[i] * cvec[i]);
                hi = std::sqrt(hi) / delta;
                double lo = 0.0;
                while ((hi - lo) > 4.0 * eps * (1.0 + hi)) {
                    const double mid = 0.5 * (lo + hi);
                    if (minres_detail::ridge_norm(m, svals.data(), cvec.data(), mid) > delta)
                        lo = mid;
                    else
                        hi = mid;
                }
                mu = hi;
            }

            // t = diag(s/(s^2+mu)) c, with the zero-mu limit truncated at stol; z = V t.
            for (std::size_t i = 0; i < msz; ++i) {
                if (active)
                    tvec[i] = svals[i] * cvec[i] / (svals[i] * svals[i] + mu);
                else
                    tvec[i] = (svals[i] > stol) ? cvec[i] / svals[i] : 0.0;
            }
            // z = Vt t, Vt row-major m x m with lda = m.
            C_DGEMV('N', m, m, 1.0, Vt.data(), m, tvec.data(), 1, 0.0, zvec.data(), 1);
            // Back to Lanczos coordinates: y = X z. Xmat is column-major m x m, i.e.
            // row-major X^T, so transa='T' recovers X z.
            C_DGEMV('T', m, m, 1.0, Xmat.data(), m, zvec.data(), 1, 0.0, yvec.data(), 1);
        } else {
            for (int i = 0; i < n; ++i) rvec[i] = -g[i];
            minres_detail::min_norm_lstsq(n, m, Amat.data(), rvec.data(), svals.data(), work.data(),
                                          lwork, iwork.data());
            for (std::size_t i = 0; i < msz; ++i) yvec[i] = rvec[i];
        }

        // k = V_k y = S (V y), the diagonal applied once: column-major V (n x mm, ld = ns) is
        // row-major (mm x n, ld = ns), so V y is the transposed action.
        C_DGEMV('T', m, n, 1.0, V.data(), n, yvec.data(), 1, 0.0, kvec.data(), 1);
        for (int l = 0; l < n; ++l) kvec[l] *= sq[l];
        // r = g + HP y: column-major HP (ns x msz, ld = ns) is row-major (msz x ns,
        // ld = ns), so HP y is the transposed gemv, accumulated onto the copied g.
        C_DCOPY(ns, const_cast<double*>(g), 1, rvec.data(), 1);
        C_DGEMV('T', static_cast<int>(msz), ns, 1.0, HP.data(), ns, yvec.data(), 1, 1.0,
                rvec.data(), 1);
        const double rn = C_DNRM2(ns, rvec.data(), 1);

        // The subspace grows monotonically, so the previous iterate stays feasible; keep
        // the smallest residual seen, which under matvec noise need not be the last.
        if (rn <= res.resid_2norm) {
            res.resid_2norm = rn;
            res.k = kvec;
            res.step_norm = C_DNRM2(ns, kvec.data(), 1);
            res.hit_radius = active;
        }

        if (rn <= eta * gnorm) break;
        // Lanczos breakdown: bnew negligible relative to the operator norm estimate.
        if (bnew <= static_cast<double>(n) * eps * anorm_est) break;
        if (j + 1 < mmax) {
            beta[j] = bnew;
            double* vnext = &V[static_cast<std::size_t>(j + 1) * ns];
            C_DCOPY(ns, w.data(), 1, vnext, 1);
            C_DSCAL(ns, 1.0 / bnew, vnext, 1);
        }
    }

    return res;
}

}  // namespace reks
}  // namespace psi

#endif  // REKS_MINRES_H
