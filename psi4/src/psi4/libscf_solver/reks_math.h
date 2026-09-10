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

#ifndef REKS_MATH_H
#define REKS_MATH_H

/// @file reks_math.h
/// @brief Pure math primitives for REKS(N,M): f(x) interpolant + derivatives,
/// symmetric and generalized (Lowdin) eigensolvers. All functions inline.

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include "reks_si_types.h"
#include "psi4/libqt/qt.h"
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

/// Eigenvectors of a block-diagonal matrix, kept one block at a time: component c owns
/// `data[off[c] .. off[c+1])` as a b_c x b_c row-major matrix, and global eigenvector j is
/// its column `local_col[j]`, nonzero only on the rows comps[block_of[j]]. Dense storage
/// would be n x n; this is sum_c b_c^2, far smaller whenever the b_c stay well below n.
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

/// Diagonalize a real symmetric matrix stored by its blocks. Each block's eigenvectors,
/// zero outside the block, are eigenvectors of the whole matrix, so blocks solve
/// independently at cost sum(b^3) instead of n^3. A partition coarser than the nonzeros
/// only merges blocks and stays correct. Eigenvalues come out ascending as from
/// lapack_diagonalize; the eigenvectors keep their block form.
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

    // Components are independent eigenproblems writing disjoint slices; block sizes span
    // orders of magnitude, hence dynamic scheduling.
    std::vector<int> eval_off(nb + 1, 0);
    for (int c = 0; c < nb; ++c) eval_off[c + 1] = eval_off[c] + static_cast<int>(comps[c].size());
    std::vector<double> blk_evals(n, 0.0);
#pragma omp parallel
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

/// Diagonalize a real symmetric matrix via LAPACK C_DSYEVD (divide and conquer).
/// Eigenvalues ascending; eigenvectors[i*n + j] = i-th component of the j-th eigenvector.
/// Input A is not modified. Throws on LAPACK failure.
inline void lapack_diagonalize(const double* A, int n,
                               std::vector<double>& eigenvalues,
                               std::vector<double>& eigenvectors) {
    eigenvalues.assign(n, 0.0);
    eigenvectors.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    if (n == 0) return;

    // C_DSYEVD overwrites its input; A is symmetric, so the row-major copy needs no transpose.
    std::vector<double> a(A, A + static_cast<size_t>(n) * static_cast<size_t>(n));

    // Divide and conquer needs 2n^2 + 6n + 1 doubles of workspace, and LAPACK takes lwork
    // as a 32-bit int: past n = 32768 the request cannot be expressed at all. QR iteration
    // asks for 3n - 1, so it is the driver that remains reachable there.
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
            throw std::runtime_error("reks::lapack_diagonalize: C_DSYEVD failed, info=" +
                                     std::to_string(info));
    } else {
        double wkopt = 0.0;
        C_DSYEV('V', 'U', n, a.data(), n, eigenvalues.data(), &wkopt, -1);
        int lwork = std::max(3 * n - 1, static_cast<int>(wkopt));
        std::vector<double> work(lwork);
        info = C_DSYEV('V', 'U', n, a.data(), n, eigenvalues.data(), work.data(), lwork);
        if (info != 0)
            throw std::runtime_error("reks::lapack_diagonalize: C_DSYEV failed, info=" +
                                     std::to_string(info));
    }

    // On exit a holds eigenvectors as column-major columns: a[j*n + i] = i-th comp of eigvec j.
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        double* const row_i = eigenvectors.data() + static_cast<size_t>(i) * n;
        const double* const col_i = a.data() + i;
        for (int j = 0; j < n; ++j) row_i[j] = col_i[static_cast<size_t>(j) * n];
    }
}

/// Bordered constrained-least-squares DIIS coefficients (Pulay equations).
///
/// Solves  min_c c^T G c  s.t.  1^T c = 1,  G the error-Gram (G_ij = <e_i,e_j>),
/// dense m x m as G[i*m + j], via the (m+1) x (m+1) bordered system solved by C_DGESV:
///     [  0  -1^T ] [ l ]   [ -1 ]
///     [ -1   G   ] [ c ] = [  0 ]
/// l the Lagrange multiplier. c is sized m and c_norm_sq set on every path; false (c all-zero)
/// if the system is singular or m <= 0. Throws on an illegal DGESV argument.
///
/// Ref: Pulay, P. Chem. Phys. Lett. 1980, 73, 393; J. Comput. Chem. 1982, 3, 556.
inline bool diis_bordered_coefficients(const double* G, int m,
                                       std::vector<double>& c, double& c_norm_sq) {
    c.assign(std::max(m, 0), 0.0);
    c_norm_sq = 0.0;
    if (m <= 0) return false;
    if (m == 1) {  // constraint forces c = [1]
        c[0] = 1.0;
        c_norm_sq = 1.0;
        return true;
    }
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

    int info = C_DGESV(dim, 1, A.data(), dim, ipiv.data(), rhs.data(), dim);
    if (info < 0)
        throw std::runtime_error("reks::diis_bordered_coefficients: C_DGESV illegal argument, "
                                 "info=" + std::to_string(info));
    if (info > 0) return false;

    for (int i = 0; i < m; ++i) {
        c[i] = rhs[i + 1];
        c_norm_sq += c[i] * c[i];
    }
    return true;
}

/// Pollock-Rebholz angle filter over a set of Pulay error vectors.
///
///   sigma_i = |orthogonal component of e_i vs span of newer survivors| / |e_i|
///           = sin(angle to that span); drop e_i when sigma_i < angle_tol (near-collinear).
///
/// Swept newest-first (largest seq); the newest vector is the anchor, always kept
/// (Pollock-Rebholz Thm 3.1 -> survivor Gram condition bounded). A vector whose norm
/// exceeds 1/sqrt(stale_delta) times the newest is dropped before the angle test
/// (Chupin stale-drop). survivors are returned in ascending physical-index order with
/// co-sorted sigmas (anchor sigma = +inf). Empty survivors signals a zero-norm anchor.
///
/// E: n row pointers, E[i] = e_i and err_size long; seq[i]: store sequence of slot i.
///
/// Ref: Pollock, S.; Rebholz, L. G. SIAM J. Sci. Comput. 2023, 45, A1571, Alg. 2.4
///      (angle filter); Chupin, M. et al. ESAIM: M2AN 2021, 55, 2785, Alg. 4 (stale-drop).
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

    // norm2[i] = |e_i|^2; enew2 is the anchor's, used by the stale-drop test below.
    std::vector<double> norm2(n);
    for (int i = 0; i < n; ++i) {
        const double* ei = E[i];
        norm2[i] = C_DDOT(err_size, const_cast<double*>(ei), 1, const_cast<double*>(ei), 1);
    }
    const double enew2 = norm2[anchor];

    // Orthonormal residuals of the kept vectors, n x err_size, sized once so the row pointers
    // taken below stay valid as n_q grows; only the first n_q rows are ever written.
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
        // MODIFIED Gram-Schmidt: each projection is subtracted from the running r right
        // after the preceding one, keeping orthogonality loss O(kappa).
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

/// Constrained-least-squares DIIS coefficients by NULL-SPACE elimination on the residual
/// matrix F (columns = error vectors), solved via SVD (C_DGELSS): the solve conditions on
/// cond(F).
///
/// Solves  min_c || F_w c ||_2  s.t.  1^T c = 1,  F_w = F diag(w),  by eliminating the
/// constraint:
///   c0 = e_anchor (unit at column anchor);  Z = orthonormal basis of {z : 1^T z = 0}
///   (a fixed Householder from the all-ones vector);  c = c0 + Z y with y from the unconstrained
///   LS  min_y || (F_w Z) y - (-F_w c0) ||_2.  delta > 0 ridges via [F_w Z; delta I] over
///   [-F_w c0; 0] (augmented-LS Tikhonov). Truncation drops singular values with
///   sigma_k/sigma_max < svd_rcond.
///   Z orthonormal => |y| = |c - c0|: ridge and rank truncation pull c toward c0.
///
///   Fcol: column-major es x m (Fcol[a*es + k] = component k of survivor column a).
///   weights: size m COLUMN scale of F.
///   anchor: local column index setting c0 = e_anchor.
///   c_out: size-m coefficients (survivor-indexed); c_norm_sq = |c|^2.
///
/// Ref: Pollock, S.; Rebholz, L. G. SIAM J. Sci. Comput. 2023, 45, A1571, Thm. 3.1;
///      Golub, G. H.; Van Loan, C. F. Matrix Computations, 4th ed.; Johns Hopkins University
///      Press: Baltimore, 2013; constrained LS via null-space basis.
inline void diis_constrained_lstsq_svd(const double* Fcol, int es, int m,
                                       const double* weights, int anchor,
                                       double tikhonov_delta, double svd_rcond,
                                       std::vector<double>& c_out, double& c_norm_sq) {
    c_out.assign(m, 0.0);
    c_norm_sq = 0.0;
    if (m <= 0) return;
    if (m == 1) {  // single survivor: constraint forces c = [1].
        c_out[0] = 1.0;
        c_norm_sq = 1.0;
        return;
    }

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
    // Zt[t*m + a] = Z(a,t), independent of k; read in the A-build and c-reconstruction below.
    std::vector<double> Zt(static_cast<size_t>(m) * std::max(ncol, 0), 0.0);
    for (int t = 0; t < ncol; ++t)
        for (int a = 0; a < m; ++a) Zt[static_cast<size_t>(t) * m + a] = Z(a, t);

    const bool ridge = (tikhonov_delta > 0.0);
    const int la_m = ridge ? (es + ncol) : es;
    const int la_n = ncol;
    const int lda = la_m;
    const int ldb = std::max(la_m, la_n);

    // Zw[t*m + a] = weights[a] * Z(a,t); Zt itself stays unweighted for the
    // c reconstruction below.
    std::vector<double> Zw(static_cast<size_t>(m) * std::max(ncol, 0), 0.0);
    for (int t = 0; t < ncol; ++t)
        for (int a = 0; a < m; ++a)
            Zw[static_cast<size_t>(t) * m + a] = weights[a] * Zt[static_cast<size_t>(t) * m + a];

    // A = F_w Z, column-major la_m x ncol (bottom ncol rows = delta*I when ridged).
    // Row-major view: A(ncol x es, ld = lda) = Zw(ncol x m, ld = m) * Fcol(m x es, ld = es);
    // the column-major es x m Fcol buffer is that row-major m x es block byte for byte.
    // beta = 0 writes only the leading es entries of each row, so the ridge tail below
    // lands in storage the GEMM never touches.
    std::vector<double> A(static_cast<size_t>(lda) * la_n, 0.0);
    if (ncol > 0 && es > 0)
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

/// Solve H*c = E*S*c via rank-revealing Lowdin (canonical orthogonalization
/// with truncation of S-null modes).
///
/// Diagonalize S = U diag(s) U^T. Partition s_i into active (s_i > sigma) and
/// null (s_i <= sigma). Build X = U_active diag(1/sqrt(s_active)). Project
/// Hp = X^T H X (size na x na) and diagonalize Hp.
///
/// Output eigvecs: column j of eigenvectors[i*n+j] is the j-th eigvec.
/// First na entries: physical eigvals (sorted ascending); last n-na entries:
/// sentinel value (mu/sigma) with eigvec set to the S-null direction, giving
/// c^T S c ~ 0.
///
/// H_in: row-major n*n symmetric (not destroyed). Eigenvalues are sized n; eigenvectors
/// are n x n_out, holding only the `n_out` leading columns: physical roots up to n_out,
/// then null directions from column na on. The full n x n form is n_out = n.
/// s_evals/s_evecs are the S eigendecomposition and `comps` the block structure it came
/// out of -- all three are what lapack_diagonalize_blocked returns.
/// An S eigenvector is zero outside its own block, so X is block-structured and every
/// product below runs on b_c x p_c factors: the transform costs sum_c b_c*p_c per column
/// instead of the dense n*na, and only the na x na eigenproblem stays dense.
inline void generalized_diagonalize(double* H_in, int n,
                                     const std::vector<double>& s_evals,
                                     const BlockedEigenvectors& s_evecs,
                                     const std::vector<std::vector<int>>& comps,
                                     double sigma, double mu, int n_out,
                                     std::vector<double>& eigenvalues,
                                     std::vector<double>& eigenvectors) {
    const int nb = static_cast<int>(comps.size());

    // Active S-eigenvectors, grouped by block; the rest is the null space. Grouping
    // permutes the columns of X and with them the basis of Hp, and c = X cp comes out in
    // that same basis, so nothing has to be permuted back.
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
        // Block c holds X_c[r, q] = U[comps[c][r], k] / sqrt(s[k]), b_c x p_c row-major,
        // k = active_idx[col_off[c] + q]. The reciprocal is formed once per column:
        // x * (1/y) and x / y differ in the last ulp.
        std::vector<size_t> x_off(nb + 1, 0);
        for (int c = 0; c < nb; ++c)
            x_off[c + 1] = x_off[c] +
                comps[c].size() * static_cast<size_t>(col_off[c + 1] - col_off[c]);
        std::vector<double> X(x_off[nb], 0.0);
#pragma omp parallel for schedule(dynamic, 1)
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

        // Hp = X^T H X, one block pair at a time: Hp[rows_a, cols_b] = X_a^T H[a,b] X_b.
        // Nothing of size n x na is ever formed -- the largest live buffer is one b_a x b_b
        // tile of H. Rows of Hp are disjoint across a, so the a-loop is the parallel one.
        std::vector<double> Hp(static_cast<size_t>(na) * na, 0.0);
#pragma omp parallel
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

        // Symmetrize: FP noise can break exact symmetry. Load-bearing, not hygiene --
        // lapack_diagonalize hands this row-major buffer to a column-major C_DSYEV with
        // uplo='U', which physically reads the lower triangle.
#pragma omp parallel for schedule(guided)
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

        // c = X cp; column j of eigenvectors is the j-th physical eigvec in n-dim. Block c
        // owns rows comps[c] and reads the cp rows its own columns occupy. Only the
        // columns the caller asked for are built; the rest of [0, na) is never written.
        for (int j = 0; j < na; ++j) eigenvalues[j] = hp_evals[j];
        const int nv = std::min(na, std::max(n_out, 0));
        timer_on("REKS: si_diag_backtransform");
#pragma omp parallel
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
