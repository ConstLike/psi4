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

#ifndef REKS_MATH_H
#define REKS_MATH_H

/// Math primitives of REKS: f(x) interpolant and derivatives, symmetric and generalized eigensolvers, DIIS.

#include <vector>
#include <cmath>
#include <string>
#include "reks_si_types.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include <algorithm>
#include <climits>
#include <utility>
#include <limits>

namespace psi {
namespace reks {

/// Interpolation parameter delta for the FON interpolant f(x). Process-wide.
inline double delta_interp = 0.4;

/// Numerical thresholds for FON optimization and convergence checks.
namespace constants {
constexpr double ENERGY_THRESHOLD = 1e-10;
constexpr double FON_THRESHOLD = 1e-10;

/// FON box-constraint slack: FON clamped to [eps, 2-eps]. 0.0 = full [0,2] box.
constexpr double FON_BOUNDARY_EPS = 0.0;
}  // namespace constants

/// Interpolant f(x) = x^(1 - (x+delta)/(2(1+delta))), x = n_p*n_q in [0,1].
inline double f_interp(double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    const double c = 0.5 / (1.0 + delta_interp);
    const double exponent = 1.0 - c * (x + delta_interp);

    return std::pow(x, exponent);
}

/// df/dx = f * (exponent/x - c*ln(x)) for x in (0,1); zero at boundaries.
inline double df_interp(double x) {
    if (x <= 1e-14) return 0.0;
    if (x >= 1.0) return 0.0;

    const double c = 0.5 / (1.0 + delta_interp);
    const double exponent = 1.0 - c * (x + delta_interp);
    const double f = std::pow(x, exponent);
    const double log_x = std::log(x);

    return f * (exponent / x - c * log_x);
}

/// d2f/dx2 = f * ((exp/x - c*ln x)^2 - exp/x^2 - 2c/x); zero at boundaries.
inline double d2f_interp(double x) {
    if (x <= 1e-14) return 0.0;
    if (x >= 1.0) return 0.0;

    const double c = 0.5 / (1.0 + delta_interp);
    const double exponent = 1.0 - c * (x + delta_interp);
    const double f = std::pow(x, exponent);
    const double log_x = std::log(x);

    const double term1 = exponent / x - c * log_x;
    const double term2 = -exponent / (x * x) - 2.0 * c / x;

    return f * (term1 * term1 + term2);
}

inline void lapack_diagonalize(const double* A, int n,
                               std::vector<double>& eigenvalues,
                               std::vector<double>& eigenvectors);

/// Eigenvectors of a block-diagonal matrix: block c = data[off[c], off[c+1]) (b_c x b_c row-major),
/// eigenvector j = its column local_col[j] on rows comps[block_of[j]].
struct BlockedEigenvectors {
    std::vector<double> data;
    std::vector<size_t> off;        ///< nb + 1 offsets into data
    std::vector<int>    block_of;   ///< n entries: eigenvector -> component
    std::vector<int>    local_col;  ///< n entries: eigenvector -> its column in that component

    /// Row `row` of component c, i.e. the coefficients of comps[c][row] over c's eigenvectors.
    const double* block_row(int c, int row, int b) const {
        return data.data() + off[c] + static_cast<size_t>(row) * b;
    }
};

/// Symmetric eigensolver per diagonal block (cost sum b^3); eigenvalues ascending as lapack_diagonalize,
/// eigenvectors in block form.
inline void lapack_diagonalize_blocked(const BlockedMatrix& A,
                                       std::vector<double>& eigenvalues,
                                       BlockedEigenvectors& eigenvectors) {
    const std::vector<std::vector<int>>& comps = A.blocks;
    const int n  = A.n;
    const int nb = static_cast<int>(comps.size());
    eigenvalues.assign(n, 0.0);
    eigenvectors.block_of.assign(n, 0);
    eigenvectors.local_col.assign(n, 0);
    eigenvectors.off.assign(std::max(nb, 1) + 1, 0);
    for (int c = 0; c < nb; ++c)
        eigenvectors.off[c + 1] = eigenvectors.off[c] + comps[c].size() * comps[c].size();
    eigenvectors.data.assign(eigenvectors.off[std::max(nb, 0)], 0.0);
    if (n == 0) return;

    // Independent per-block eigenproblems writing disjoint slices.
    std::vector<int> eval_off(nb + 1, 0);
    for (int c = 0; c < nb; ++c) eval_off[c + 1] = eval_off[c] + static_cast<int>(comps[c].size());
    std::vector<double> blk_evals(n, 0.0);
    const int nthread = Process::environment.get_n_threads();
#pragma omp parallel num_threads(nthread)
    {
        std::vector<double> c_evals, c_evecs;
#pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < nb; ++c) {
            const int b = static_cast<int>(comps[c].size());
            double* const dst = eigenvectors.data.data() + eigenvectors.off[c];
            if (b == 1) {
                // 1x1 block: the entry is the eigenvalue, the eigenvector is 1.
                dst[0] = 1.0;
                blk_evals[eval_off[c]] = A.tile(c)[0];
                continue;
            }
            lapack_diagonalize(A.tile(c), b, c_evals, c_evecs);
            std::copy(c_evecs.begin(), c_evecs.end(), dst);
            std::copy(c_evals.begin(), c_evals.end(), blk_evals.begin() + eval_off[c]);
        }
    }

    // Order all eigenpairs by eigenvalue as LAPACK would.
    std::vector<std::pair<double, std::pair<int, int>>> order;  // (eigenvalue, (comp, local col))
    order.reserve(n);
    for (int c = 0; c < nb; ++c)
        for (int q = 0; q < static_cast<int>(comps[c].size()); ++q)
            order.emplace_back(blk_evals[eval_off[c] + q], std::make_pair(c, q));
    // stable: ties keep deterministic (component, column) order.
    std::stable_sort(order.begin(), order.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });

    for (int col = 0; col < n; ++col) {
        eigenvalues[col] = order[col].first;
        eigenvectors.block_of[col]  = order[col].second.first;
        eigenvectors.local_col[col] = order[col].second.second;
    }
}

/// Symmetric eigensolver (C_DSYEVD): eigenvalues ascending, eigenvectors[i*n + j] = component i of
/// vector j; A is kept, throws on LAPACK failure.
inline void lapack_diagonalize(const double* A, int n,
                               std::vector<double>& eigenvalues,
                               std::vector<double>& eigenvectors) {
    eigenvalues.assign(n, 0.0);
    eigenvectors.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    if (n == 0) return;

    // C_DSYEVD overwrites its input; A is symmetric, so the row-major copy needs no transpose.
    std::vector<double> a(A, A + static_cast<size_t>(n) * static_cast<size_t>(n));

    // DSYEVD needs lwork = 2n^2 + 6n + 1, beyond 32-bit int past n = 32768; DSYEV (3n - 1) covers that.
    const long long lwork_dc = 1 + 6LL * n + 2LL * n * n;
    const bool use_divide_conquer = lwork_dc <= static_cast<long long>(INT_MAX);

    int info = 0;
    if (use_divide_conquer) {
        // lwork/liwork = -1 asks DSYEVD for its optimal workspace sizes.
        double wkopt = 0.0;
        int iwkopt = 0;
        C_DSYEVD('V', 'U', n, a.data(), n, eigenvalues.data(), &wkopt, -1, &iwkopt, -1);
        int lwork  = std::max(static_cast<int>(lwork_dc), static_cast<int>(wkopt));
        int liwork = std::max(3 + 5 * n, iwkopt);
        std::vector<double> work(lwork);
        std::vector<int>    iwork(liwork);
        info = C_DSYEVD('V', 'U', n, a.data(), n, eigenvalues.data(), work.data(), lwork,
                        iwork.data(), liwork);
        if (info != 0)
            throw PSIEXCEPTION("reks::lapack_diagonalize: C_DSYEVD failed, info=" +
                               std::to_string(info));
    } else {
        double wkopt = 0.0;
        C_DSYEV('V', 'U', n, a.data(), n, eigenvalues.data(), &wkopt, -1);
        int lwork = std::max(3 * n - 1, static_cast<int>(wkopt));
        std::vector<double> work(lwork);
        info = C_DSYEV('V', 'U', n, a.data(), n, eigenvalues.data(), work.data(), lwork);
        if (info != 0)
            throw PSIEXCEPTION("reks::lapack_diagonalize: C_DSYEV failed, info=" +
                               std::to_string(info));
    }

    // On exit a holds eigenvectors as column-major columns: a[j*n + i] = i-th comp of eigvec j.
    const int nthread = Process::environment.get_n_threads();
#pragma omp parallel for schedule(static) num_threads(nthread)
    for (int i = 0; i < n; ++i) {
        double* const row_i = eigenvectors.data() + static_cast<size_t>(i) * n;
        const double* const col_i = a.data() + i;
        for (int j = 0; j < n; ++j) row_i[j] = col_i[static_cast<size_t>(j) * n];
    }
}

/// DIIS min c^T G c, 1^T c = 1, from the bordered system [[0, -1^T], [-1, G]] [l; c] = [-1; 0] (C_DGESV);
/// m >= 2; false, c = 0 if singular (Pulay, P. Chem. Phys. Lett. 1980, 73, 393; J. Comput. Chem. 1982, 3, 556).
inline bool diis_bordered_coefficients(const double* G, int m,
                                       std::vector<double>& c, double& c_norm_sq) {
    c.assign(m, 0.0);
    c_norm_sq = 0.0;
    const int dim = m + 1;
    std::vector<double> A(static_cast<size_t>(dim) * dim, 0.0);
    std::vector<double> rhs(dim, 0.0);
    std::vector<int> ipiv(dim);

    for (int i = 1; i <= m; ++i) {
        A[i] = -1.0;                                  // row 0
        A[static_cast<size_t>(i) * dim] = -1.0;       // col 0
    }
    rhs[0] = -1.0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            A[static_cast<size_t>(i + 1) * dim + (j + 1)] = G[static_cast<size_t>(i) * m + j];

    const int info = C_DGESV(dim, 1, A.data(), dim, ipiv.data(), rhs.data(), dim);
    if (info != 0) return false;

    for (int i = 0; i < m; ++i) {
        c[i] = rhs[i + 1];
        c_norm_sq += c[i] * c[i];
    }
    return true;
}

/// Newest-first, keeping the anchor: drop |e_i| > |e_new| / sqrt(stale_delta) (Chupin, M. et al. ESAIM: M2AN 2021,
/// 55, 2785), then sin(e_i, newer) < angle_tol (Pollock, S.; Rebholz, L. G. SIAM J. Sci. Comput. 2023, 45, A1571).
inline void diis_angle_filter(const double* const* E, const long* seq, int n, int err_size,
                              double angle_tol, double stale_delta,
                              std::vector<int>& survivors, std::vector<double>& sigmas,
                              std::vector<double>& q_scratch) {
    survivors.clear();
    sigmas.clear();
    if (n <= 0) return;

    std::vector<int> ord(n);
    for (int i = 0; i < n; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [seq](int a, int b) { return seq[a] > seq[b]; });
    const int anchor = ord[0];

    // norm2[i] = |e_i|^2; enew2 = |e_anchor|^2.
    std::vector<double> norm2(n);
    for (int i = 0; i < n; ++i) {
        const double* ei = E[i];
        norm2[i] = C_DDOT(err_size, const_cast<double*>(ei), 1, const_cast<double*>(ei), 1);
    }
    const double enew2 = norm2[anchor];

    // Orthonormal residuals q of the kept vectors; the first n_q of n rows are written.
    q_scratch.resize(static_cast<size_t>(n) * err_size);
    std::vector<double>& Q = q_scratch;
    int n_q = 0;
    std::vector<double> r(err_size);
    for (int idx : ord) {
        const double enorm = std::sqrt(norm2[idx]);
        if (idx == anchor && !(enorm > 0.0)) {
            survivors.clear();
            sigmas.clear();
            return;
        }
        if (idx != anchor && !(enorm > 0.0)) continue;
        if (!survivors.empty() && stale_delta * enorm * enorm >= enew2) continue;  // Chupin stale-drop

        const double* ei = E[idx];
        C_DCOPY(err_size, const_cast<double*>(ei), 1, r.data(), 1);
        // Modified Gram-Schmidt: r <- r - (r.q) q, sequentially over kept q.
        for (int q = 0; q < n_q; ++q) {
            double* qk = Q.data() + static_cast<size_t>(q) * err_size;
            const double proj = C_DDOT(err_size, r.data(), 1, qk, 1);
            C_DAXPY(err_size, -proj, qk, 1, r.data(), 1);
        }
        // sqrt of the plain sum of squares, not C_DNRM2: the scaled recurrence returns
        // a slightly different value.
        const double rnorm = std::sqrt(C_DDOT(err_size, r.data(), 1, r.data(), 1));
        const double sigma = rnorm / enorm;

        const bool keep = survivors.empty()  // anchor: kept unconditionally (rnorm == enorm > 0)
                          || sigma >= angle_tol;
        if (keep) {
            // Element-wise divide, not a C_DSCAL by the reciprocal: that would add a rounding.
            double* qnew = Q.data() + static_cast<size_t>(n_q) * err_size;
            for (int k = 0; k < err_size; ++k) qnew[k] = r[k] / rnorm;
            ++n_q;
            sigmas.push_back(survivors.empty() ? std::numeric_limits<double>::infinity() : sigma);
            survivors.push_back(idx);
        }
    }

    // Co-sort (survivors, sigmas) to ascending physical index.
    std::vector<int> perm(survivors.size());
    for (size_t i = 0; i < perm.size(); ++i) perm[i] = static_cast<int>(i);
    std::sort(perm.begin(), perm.end(),
              [&](int a, int b) { return survivors[a] < survivors[b]; });
    std::vector<int> s2(survivors.size());
    std::vector<double> g2(sigmas.size());
    for (size_t i = 0; i < perm.size(); ++i) {
        s2[i] = survivors[perm[i]];
        g2[i] = sigmas[perm[i]];
    }
    survivors.swap(s2);
    sigmas.swap(g2);
}

/// DIIS min |F diag(w) c|, 1^T c = 1: c = e_anchor + Z y, Z = Householder basis of 1^T z = 0, y by C_DGELSS
/// (rcond svd_rcond) with ridge [F_w Z; delta I]; Fcol column-major es x m, m >= 2 (Pollock, Rebholz 2023, Thm. 3.1).
inline void diis_constrained_lstsq_svd(const double* Fcol, int es, int m,
                                       const double* weights, int anchor,
                                       double tikhonov_delta, double svd_rcond,
                                       std::vector<double>& c_out, double& c_norm_sq) {
    c_out.assign(m, 0.0);
    c_norm_sq = 0.0;

    // Householder reflector mapping ones -> multiple of e1; its columns 1..m-1 form Z.
    const double alpha = std::sqrt(static_cast<double>(m));  // ||ones||
    std::vector<double> hv(m, 1.0);
    hv[0] = 1.0 + alpha;  // ones + sign(ones[0])*alpha*e1, sign = +1
    double vtv = 0.0;
    for (int i = 0; i < m; ++i) vtv += hv[i] * hv[i];
    auto Z = [&](int i, int t) {
        double e = (i == t + 1) ? 1.0 : 0.0;
        return e - 2.0 * hv[i] * hv[t + 1] / vtv;
    };

    const int ncol = m - 1;
    // Zt[t*m + a] = Z(a,t)
    std::vector<double> Zt(static_cast<size_t>(m) * ncol, 0.0);
    for (int t = 0; t < ncol; ++t)
        for (int a = 0; a < m; ++a) Zt[static_cast<size_t>(t) * m + a] = Z(a, t);

    const bool ridge = (tikhonov_delta > 0.0);
    const int la_m = ridge ? (es + ncol) : es;
    const int la_n = ncol;
    const int lda = la_m;
    const int ldb = std::max(la_m, la_n);

    // Zw[t*m + a] = weights[a] * Z(a,t)
    std::vector<double> Zw(static_cast<size_t>(m) * ncol, 0.0);
    for (int t = 0; t < ncol; ++t)
        for (int a = 0; a < m; ++a)
            Zw[static_cast<size_t>(t) * m + a] = weights[a] * Zt[static_cast<size_t>(t) * m + a];

    // A = F_w Z column-major la_m x ncol = row-major Zw (ncol x m) * Fcol (m x es); the ridge rows
    // delta I sit past the es entries the GEMM writes.
    std::vector<double> A(static_cast<size_t>(lda) * la_n, 0.0);
    if (es > 0)
        C_DGEMM('N', 'N', ncol, es, m, 1.0, Zw.data(), m, const_cast<double*>(Fcol), es, 0.0,
                A.data(), lda);
    if (ridge)
        for (int t = 0; t < ncol; ++t) A[static_cast<size_t>(t) * lda + es + t] = tikhonov_delta;

    // rhs = -F_w c0 = -weights[anchor] * F[:,anchor]; ridge tail = 0.
    std::vector<double> b(ldb, 0.0);
    for (int k = 0; k < es; ++k)
        b[k] = -weights[anchor] * Fcol[static_cast<size_t>(anchor) * es + k];

    std::vector<double> s(std::min(la_m, la_n), 0.0);
    int rank = 0;

    double wkopt = 0.0;
    C_DGELSS(la_m, la_n, 1, A.data(), lda, b.data(), ldb, s.data(), svd_rcond, &rank, &wkopt, -1);
    int lwork = std::max(1, static_cast<int>(wkopt));
    std::vector<double> work(lwork);
    C_DGELSS(la_m, la_n, 1, A.data(), lda, b.data(), ldb, s.data(), svd_rcond, &rank, work.data(),
             lwork);

    // c = c0 + Z y,  y = b[0..ncol-1] (DGELSS solution, in place).
    for (int a = 0; a < m; ++a) {
        double val = (a == anchor) ? 1.0 : 0.0;
        for (int t = 0; t < ncol; ++t) val += Zt[static_cast<size_t>(t) * m + a] * b[t];
        c_out[a] = val;
        c_norm_sq += val * val;
    }
}

/// H c = E S c by canonical orthogonalization: X = U_s diag(s^-1/2) over s > sigma (blocked S from
/// lapack_diagonalize_blocked), Hp = X^T H X; roots na..n-1 = S-null sentinels; n_out leading columns.
inline void generalized_diagonalize(double* H_in, int n,
                                     const std::vector<double>& s_evals,
                                     const BlockedEigenvectors& s_evecs,
                                     const std::vector<std::vector<int>>& comps,
                                     double sigma, double mu, int n_out,
                                     std::vector<double>& eigenvalues,
                                     std::vector<double>& eigenvectors) {
    const int nb = static_cast<int>(comps.size());
    const int nthread = Process::environment.get_n_threads();

    // Active S eigenvectors grouped by block (c = X cp needs no permutation back); the rest is null.
    std::vector<int> null_idx;
    std::vector<std::vector<int>> cols_of(nb);
    for (int k = 0; k < n; ++k) {
        if (s_evals[k] > sigma) cols_of[s_evecs.block_of[k]].push_back(k);
        else                    null_idx.push_back(k);
    }
    std::vector<int> active_idx;
    active_idx.reserve(n);
    std::vector<int> col_off(nb + 1, 0);
    for (int c = 0; c < nb; ++c) {
        col_off[c] = static_cast<int>(active_idx.size());
        active_idx.insert(active_idx.end(), cols_of[c].begin(), cols_of[c].end());
    }
    const int na = static_cast<int>(active_idx.size());
    col_off[nb] = na;
    const int nn = static_cast<int>(null_idx.size());

    eigenvalues.assign(n, 0.0);
    eigenvectors.assign(static_cast<size_t>(n) * std::max(n_out, 0), 0.0);

    if (na > 0) {
        timer_on("REKS: si_diag_lowdin");
        // X_c[r, q] = U[comps[c][r], k] / sqrt(s[k]), k = active_idx[col_off[c] + q] (b_c x p_c row-major).
        std::vector<size_t> x_off(nb + 1, 0);
        for (int c = 0; c < nb; ++c)
            x_off[c + 1] = x_off[c] +
                comps[c].size() * static_cast<size_t>(col_off[c + 1] - col_off[c]);
        std::vector<double> X(x_off[nb], 0.0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread)
        for (int c = 0; c < nb; ++c) {
            const int b = static_cast<int>(comps[c].size());
            const int p = col_off[c + 1] - col_off[c];
            double* const Xc = X.data() + x_off[c];
            for (int q = 0; q < p; ++q) {
                const int k = active_idx[col_off[c] + q];
                const double inv = 1.0 / std::sqrt(s_evals[k]);
                const int lc = s_evecs.local_col[k];
                for (int r = 0; r < b; ++r)
                    Xc[static_cast<size_t>(r) * p + q] = s_evecs.block_row(c, r, b)[lc] * inv;
            }
        }

        // Hp[rows_a, cols_b] = X_a^T H[a,b] X_b per block pair; rows disjoint across a (parallel loop).
        std::vector<double> Hp(static_cast<size_t>(na) * na, 0.0);
#pragma omp parallel num_threads(nthread)
        {
            std::vector<double> H_tile, W;
#pragma omp for schedule(dynamic, 1)
            for (int a = 0; a < nb; ++a) {
                const std::vector<int>& idx_a = comps[a];
                const int ba = static_cast<int>(idx_a.size());
                const int pa = col_off[a + 1] - col_off[a];
                if (pa == 0) continue;
                for (int b = 0; b < nb; ++b) {
                    const std::vector<int>& idx_b = comps[b];
                    const int bb = static_cast<int>(idx_b.size());
                    const int pb = col_off[b + 1] - col_off[b];
                    if (pb == 0) continue;

                    // A single component covers every row and column in order, so H is
                    // already the tile.
                    double* Hab = H_in;
                    int ldh = n;
                    if (nb > 1) {
                        H_tile.assign(static_cast<size_t>(ba) * bb, 0.0);
                        for (int r = 0; r < ba; ++r) {
                            const double* const h_row = H_in + static_cast<size_t>(idx_a[r]) * n;
                            double* const dst = H_tile.data() + static_cast<size_t>(r) * bb;
                            for (int s = 0; s < bb; ++s) dst[s] = h_row[idx_b[s]];
                        }
                        Hab = H_tile.data();
                        ldh = bb;
                    }
                    W.assign(static_cast<size_t>(pa) * bb, 0.0);
                    C_DGEMM('T', 'N', pa, bb, ba, 1.0, X.data() + x_off[a], pa, Hab, ldh,
                            0.0, W.data(), bb);
                    C_DGEMM('N', 'N', pa, pb, bb, 1.0, W.data(), bb, X.data() + x_off[b], pb,
                            0.0, Hp.data() + static_cast<size_t>(col_off[a]) * na + col_off[b],
                            na);
                }
            }
        }

        // Symmetrize: column-major C_DSYEV with uplo 'U' reads the lower triangle of this buffer.
#pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int i = 0; i < na; ++i) {
            double* const hp_row_i = Hp.data() + static_cast<size_t>(i) * na;
            for (int j = i + 1; j < na; ++j) {
                double& upper = hp_row_i[j];
                double& lower = Hp[static_cast<size_t>(j) * na + i];
                const double avg = 0.5 * (upper + lower);
                upper = lower = avg;
            }
        }

        timer_off("REKS: si_diag_lowdin");

        std::vector<double> hp_evals, hp_evecs;
        timer_on("REKS: si_diag_eigen");
        lapack_diagonalize(Hp.data(), na, hp_evals, hp_evecs);
        timer_off("REKS: si_diag_eigen");

        // c = X cp per block (rows comps[c]), only for the requested columns.
        for (int j = 0; j < na; ++j) eigenvalues[j] = hp_evals[j];
        const int nv = std::min(na, std::max(n_out, 0));
        timer_on("REKS: si_diag_backtransform");
#pragma omp parallel num_threads(nthread)
        {
            std::vector<double> C_block;
#pragma omp for schedule(dynamic, 1)
            for (int c = 0; c < nb; ++c) {
                const std::vector<int>& idx = comps[c];
                const int b = static_cast<int>(idx.size());
                const int p = col_off[c + 1] - col_off[c];
                if (p == 0 || nv == 0) continue;
                double* const cp = hp_evecs.data() + static_cast<size_t>(col_off[c]) * na;
                if (b == n) {
                    C_DGEMM('N', 'N', n, nv, p, 1.0, X.data() + x_off[c], p, cp, na, 0.0,
                            eigenvectors.data(), n_out);
                    continue;
                }
                C_block.assign(static_cast<size_t>(b) * nv, 0.0);
                C_DGEMM('N', 'N', b, nv, p, 1.0, X.data() + x_off[c], p, cp, na, 0.0,
                        C_block.data(), nv);
                for (int r = 0; r < b; ++r)
                    std::copy_n(C_block.data() + static_cast<size_t>(r) * nv, nv,
                                eigenvectors.data() + static_cast<size_t>(idx[r]) * n_out);
            }
        }
        timer_off("REKS: si_diag_backtransform");
    }

    // Null states: column = original S-null eigvec; eigval = sentinel mu/sigma. A null
    // eigenvector is nonzero only on its own block's rows, so it is scattered there.
    const double sentinel = mu / sigma;
    for (int kn = 0; kn < nn; ++kn) eigenvalues[na + kn] = sentinel;
    for (int kn = 0; kn < nn && na + kn < n_out; ++kn) {
        const int k = null_idx[kn];
        const int c = s_evecs.block_of[k];
        const std::vector<int>& idx = comps[c];
        const int b = static_cast<int>(idx.size());
        const int lc = s_evecs.local_col[k];
        for (int r = 0; r < b; ++r)
            eigenvectors[static_cast<size_t>(idx[r]) * n_out + na + kn] =
                s_evecs.block_row(c, r, b)[lc];
    }
}

}  // namespace reks
}  // namespace psi

#endif  // REKS_MATH_H
