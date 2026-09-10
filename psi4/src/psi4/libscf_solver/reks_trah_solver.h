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

#ifndef REKS_TRAH_SOLVER_H
#define REKS_TRAH_SOLVER_H

/// @file reks_trah_solver.h
/// @brief Trust-Region solver via level-shift on H eigendecomposition -- pure math
///
/// Solves the box-constrained trust-region subproblem:
///
///     min_s   g^T s + 0.5 s^T H s
///     s.t.    ||s|| <= trust_radius
///             lower <= x + s <= upper
///
/// Algorithm:
///   1. Eigendecompose H = V diag(lambda) V^T
///   2. If H is positive definite and Newton step fits: use Newton (mu = 0)
///   3. Otherwise: bisect on level shift mu < min(lambda) such that
///      ||s(mu)|| = trust_radius, where s(mu) = -(H - mu*I)^{-1} g
///
/// The step s(mu) is ALWAYS a descent direction (g^T s < 0) because
/// (H - mu*I) is positive definite for mu < min(lambda).
/// ||s(mu)|| is monotonically increasing in mu (as mu -> min(lambda)),
/// guaranteeing bisection convergence.
///
/// Box constraints handled by single-pass active-set: fix violated variables
/// at their bounds, adjust effective gradient for coupling, re-solve reduced TRS.
///
/// References:
///   Nocedal, J.; Wright, S. J. Numerical Optimization, 2nd ed.; Springer: New York, 2006;
///     Ch. 4 (trust region)
///   More, J. J.; Sorensen, D. C. SIAM J. Sci. Stat. Comput. 1983, 4, 553 (level shift for TRS)
///   Helmich-Paris, B. J. Chem. Phys. 2021, 154, 164104
///   Helmich-Paris, B. J. Chem. Phys. 2022, 156, 204104

#include <vector>
#include <cmath>
#include <algorithm>

#include "reks_math.h"

namespace psi {
namespace reks {

/// Result of a single TRAH step computation.
///
/// predicted_decrease = pred_orb + sum_b pred_fon[b] + sum_b pred_cross[b].
/// Block-diagonal Hessian (FON blocks decouple) gives zero FON x FON cross
/// terms. pred_fon/pred_cross indexed by FON block b, parallel to
/// CombinedGradient::fon_blocks.
struct TRAHResult {
    std::vector<double> step;       ///< Step vector (size N)
    double predicted_decrease;      ///< Predicted energy decrease (positive = good)
    double pred_orb = 0.0;          ///< Orbital subspace contribution
    std::vector<double> pred_fon;     ///< FON subspace contribution, per FON block
    std::vector<double> pred_cross;   ///< Orbital x FON cross contribution, per FON block
    double mu;                      ///< AH eigenvalue (level shift)
    int n_active_bounds;            ///< Number of variables fixed at bounds
    bool step_on_boundary;          ///< Step reached trust region boundary
};

/// Persistent state across macro-iterations for trust-radius management
struct TRAHState {
    double trust_radius = 0.4;     ///< Current trust radius
    double prev_energy = 0.0;      ///< E_SA from previous macro-iteration
    double prev_predicted = 0.0;   ///< Predicted decrease from previous step
    double prev_pred_orb = 0.0;        ///< Previous pred: orbital contribution
    std::vector<double> prev_pred_fon;     ///< Previous pred: FON per FON block
    std::vector<double> prev_pred_cross;   ///< Previous pred: orbital x FON cross per FON block
    bool initialized = false;      ///< Has prev_energy been set?
    double prev_rho = 1.0;         ///< rho from previous step (actual/predicted)
    bool prev_step_on_boundary = false;  ///< Whether previous step hit TR boundary

    static constexpr double TR_MIN = 0.01;  ///< Minimum trust radius
    static constexpr double TR_MAX = 2.0;   ///< Maximum trust radius
};

/// Pure-math trust-region augmented Hessian solver.
/// All methods are static. Cost: O(N^3) for the eigendecomposition.
class TrustRegionSolver {
   public:
    /// Compute trust-region step with box constraints.
    ///
    /// @param g       Gradient vector (size N)
    /// @param H       Hessian matrix (N*N, row-major)
    /// @param lower   Lower bounds on x (size N; use -1e30 for unconstrained)
    /// @param upper   Upper bounds on x (size N; use +1e30 for unconstrained)
    /// @param x       Current point (size N)
    /// @param state   Trust region state (read-only: uses trust_radius)
    /// @return        TRAHResult with step, predicted decrease, diagnostics
    static TRAHResult compute_step(const std::vector<double>& g, const std::vector<double>& H,
                                   const std::vector<double>& lower, const std::vector<double>& upper,
                                   const std::vector<double>& x, const TRAHState& state) {
        int N = static_cast<int>(g.size());

        TRAHResult result;
        result.step.assign(N, 0.0);
        result.predicted_decrease = 0.0;
        result.mu = 0.0;
        result.n_active_bounds = 0;
        result.step_on_boundary = false;

        if (N == 0 || state.trust_radius < 1e-14) return result;

        // Unconstrained TRS solve over all N variables.
        double mu_val = 0.0;
        solve_ah(g.data(), H.data(), N, state.trust_radius, result.step.data(), mu_val);
        result.mu = mu_val;

        // Box-constraint check: clamp steps that exit [lower, upper].
        std::vector<bool> fixed(N, false);
        for (int i = 0; i < N; ++i) {
            double x_new = x[i] + result.step[i];
            if (x_new < lower[i]) {
                result.step[i] = lower[i] - x[i];
                fixed[i] = true;
            } else if (x_new > upper[i]) {
                result.step[i] = upper[i] - x[i];
                fixed[i] = true;
            }
        }

        int n_fixed = 0;
        for (int i = 0; i < N; ++i) {
            if (fixed[i]) n_fixed++;
        }
        result.n_active_bounds = n_fixed;

        // Re-solve on the free subspace, folding fixed-variable coupling into the gradient.
        if (n_fixed > 0 && n_fixed < N) {
            std::vector<int> free_idx;
            for (int i = 0; i < N; ++i) {
                if (!fixed[i]) free_idx.push_back(i);
            }
            int M = static_cast<int>(free_idx.size());

            // Effective gradient: g_eff[a] = g[ia] + sum_{fixed j} H[ia,j] * step[j]
            std::vector<double> g_eff(M);
            for (int a = 0; a < M; ++a) {
                int ia = free_idx[a];
                g_eff[a] = g[ia];
                for (int j = 0; j < N; ++j) {
                    if (fixed[j]) {
                        g_eff[a] += H[ia * N + j] * result.step[j];
                    }
                }
            }

            std::vector<double> H_red(M * M);
            for (int a = 0; a < M; ++a) {
                for (int b = 0; b < M; ++b) {
                    H_red[a * M + b] = H[free_idx[a] * N + free_idx[b]];
                }
            }

            // h_eff = sqrt(trust_radius^2 - fixed_norm_sq): L2 step budget left for the free block
            double fixed_norm_sq = 0.0;
            for (int i = 0; i < N; ++i) {
                if (fixed[i]) fixed_norm_sq += result.step[i] * result.step[i];
            }
            double h_sq = state.trust_radius * state.trust_radius - fixed_norm_sq;
            if (h_sq > 1e-14) {
                double h_eff = std::sqrt(h_sq);
                std::vector<double> s_red(M, 0.0);
                double mu2 = 0.0;
                solve_ah(g_eff.data(), H_red.data(), M, h_eff, s_red.data(), mu2);
                for (int a = 0; a < M; ++a) {
                    result.step[free_idx[a]] = s_red[a];
                }
            }
        }

        double step_norm_sq = 0.0;
        for (int i = 0; i < N; ++i) step_norm_sq += result.step[i] * result.step[i];
        result.step_on_boundary = (std::sqrt(step_norm_sq) >= 0.99 * state.trust_radius);

        // Predicted decrease on the actual (post-clamp) step, so rho stays meaningful.
        result.predicted_decrease = predicted_decrease(g.data(), H.data(), result.step.data(), N);

        return result;
    }

    /// Update trust radius after observing actual energy.
    ///
    /// Schedule (Nocedal-Wright with boundary-aware expansion):
    ///   rho > 0.75 AND step on boundary :  TR *= 2.0  (model accurate, TR was binding)
    ///   rho > 0.75 AND step interior     :  TR *= 1.2  (model accurate, TR not binding)
    ///   0.25 < rho <= 0.75               :  TR *= 1.05 (gentle growth)
    ///   0 <= rho <= 0.25                 :  TR *= 0.7  (poor model, shrink)
    ///   rho < 0                          :  TR *= 0.5, step always accepted (no-reject policy)
    ///
    /// No-reject: a rejected step wastes an SCF iteration; TR shrinkage alone is
    /// sufficient safety.
    ///
    /// @param actual_energy  Current E_SA = sum(C_L * E_L)
    /// @param state          Modified: trust_radius updated, prev_energy stored
    static void update_trust_radius(double actual_energy, TRAHState& state) {
        if (!state.initialized) {
            state.prev_energy = actual_energy;
            state.initialized = true;
            return;
        }

        double actual_decrease = state.prev_energy - actual_energy;

        if (state.prev_predicted > 1e-14) {
            double rho = actual_decrease / state.prev_predicted;

            if (rho > 0.75) {
                state.trust_radius *= state.prev_step_on_boundary ? 2.0 : 1.2;
            } else if (rho > 0.25) {
                state.trust_radius *= 1.05;
            } else if (rho >= 0.0) {
                state.trust_radius *= 0.7;
            } else {
                state.trust_radius *= 0.5;
            }
        }
        // If prev_predicted ~ 0: gradient was negligible (converged), no TR change

        state.trust_radius = std::clamp(state.trust_radius, TRAHState::TR_MIN, TRAHState::TR_MAX);
        state.prev_energy = actual_energy;
    }

    /// Compute predicted energy decrease for a given step.
    /// pred = -(g^T s + 0.5 s^T H s).  Positive if step decreases energy.
    static double predicted_decrease(const double* g, const double* H, const double* step, int N) {
        if (N == 0) return 0.0;
        // C_DGEMV, not C_DSYMV: H is only symmetric up to rounding and the full
        // product averages the two triangles, which DSYMV would discard.
        std::vector<double> Hs(N, 0.0);
        C_DGEMV('N', N, N, 1.0, const_cast<double*>(H), N, const_cast<double*>(step), 1, 0.0,
                Hs.data(), 1);
        const double gs  = C_DDOT(N, const_cast<double*>(g), 1, const_cast<double*>(step), 1);
        const double sHs = C_DDOT(N, const_cast<double*>(step), 1, Hs.data(), 1);
        return -(gs + 0.5 * sHs);
    }

   private:
    /// Level-shift solve of the unconstrained TRS (file-level algorithm above).
    ///
    /// @param g     Gradient (size N)
    /// @param H     Hessian (NxN row-major)
    /// @param N     Problem dimension
    /// @param h     Trust radius
    /// @param step  Output: step vector (size N, caller-allocated)
    /// @param mu    Output: level shift (mu=0 for Newton, else mu < lambda_min(H))
    static void solve_ah(const double* g, const double* H, int N, double h, double* step, double& mu) {
        double g_norm_sq = 0.0;
        for (int i = 0; i < N; ++i) g_norm_sq += g[i] * g[i];
        if (g_norm_sq < 1e-28) {
            std::fill(step, step + N, 0.0);
            mu = 0.0;
            return;
        }

        std::vector<double> lambda, V;
        lapack_diagonalize(H, N, lambda, V);
        // lambda: eigenvalues sorted ascending
        // V[k*N+i] = k-th component of i-th eigenvector

        // gt = V^T g, gradient rotated into eigenvector basis; transa='T' since
        // eigenvectors are columns of row-major V.
        std::vector<double> gt(N, 0.0);
        C_DGEMV('T', N, N, 1.0, V.data(), N, const_cast<double*>(g), 1, 0.0, gt.data(), 1);

        double lambda_min = lambda[0];

        // Newton step (mu=0) if H is PD and fits: ||s||^2 = sum_i gt[i]^2/lambda[i]^2 <= h^2
        if (lambda_min > 1e-10) {
            double newton_norm_sq = 0.0;
            for (int i = 0; i < N; ++i) {
                newton_norm_sq += gt[i] * gt[i] / (lambda[i] * lambda[i]);
            }
            if (std::sqrt(newton_norm_sq) <= h) {
                mu = 0.0;
                build_step_from_eigdecomp(gt, lambda, V, N, 0.0, step);
                return;
            }
        }

        // Level-shift bisection to find ||s(mu)|| = h
        //
        // ||s(mu)||^2 = sum_i gt[i]^2 / (lambda[i] - mu)^2
        // is monotonically increasing as mu -> lambda_min from below.
        // At mu = -inf: ||s|| -> 0.  At mu -> lambda_min^-: ||s|| -> inf.

        auto step_norm_at = [&](double trial_mu) -> double {
            double norm_sq = 0.0;
            for (int i = 0; i < N; ++i) {
                double d = lambda[i] - trial_mu;
                norm_sq += gt[i] * gt[i] / (d * d);
            }
            return std::sqrt(norm_sq);
        };

        double mu_hi = lambda_min - 1e-8;  // Just below lambda_min: large step
        double mu_lo = lambda_min - 10.0;  // Far below: small step

        // Ensure mu_lo gives step norm < h
        for (int i = 0; i < 30 && step_norm_at(mu_lo) > h; ++i) {
            mu_lo = 2.0 * mu_lo - lambda_min;
        }

        // Handle hard case: even near lambda_min, step is small
        // (gradient orthogonal to eigenvector of lambda_min)
        if (step_norm_at(mu_hi) <= h) {
            mu = mu_hi;
            build_step_from_eigdecomp(gt, lambda, V, N, mu, step);
            return;
        }

        mu = mu_lo;
        static constexpr int MAX_BISECT = 50;

        for (int iter = 0; iter < MAX_BISECT; ++iter) {
            double mu_mid = 0.5 * (mu_lo + mu_hi);
            double sn = step_norm_at(mu_mid);

            if (sn > h) {
                mu_hi = mu_mid;  // Step too large: shift further from lambda_min
            } else {
                mu_lo = mu_mid;  // Step too small: shift closer to lambda_min
            }

            mu = mu_mid;

            // Converged: step norm within 0.1% of trust radius
            if (std::abs(sn / h - 1.0) < 1e-3) break;
        }

        build_step_from_eigdecomp(gt, lambda, V, N, mu, step);
    }

    /// Compute step from H eigendecomposition and level shift:
    ///   s = -sum_i (gt[i] / (lambda[i] - mu)) * V[:,i]
    static void build_step_from_eigdecomp(const std::vector<double>& gt, const std::vector<double>& lambda,
                                          const std::vector<double>& V, int N, double mu, double* step) {
        std::vector<double> coeff(N);
        for (int i = 0; i < N; ++i) coeff[i] = -gt[i] / (lambda[i] - mu);
        // s = V coeff; transa='N' since eigenvector i is column i of row-major V.
        // beta = 0 overwrites every entry of step.
        C_DGEMV('N', N, N, 1.0, const_cast<double*>(V.data()), N, coeff.data(), 1, 0.0, step, 1);
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_TRAH_SOLVER_H
