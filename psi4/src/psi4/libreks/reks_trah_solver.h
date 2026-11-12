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

#ifndef REKS_TRAH_SOLVER_H
#define REKS_TRAH_SOLVER_H

/// min g^T s + s^T H s / 2, |s| <= h, lower <= x + s <= upper: Newton if H > 0 and it fits, else s = -(H - mu)^-1 g
/// with mu < lambda_min bisected to |s| = h; violated bounds are fixed once and the rest re-solved.
/// Nocedal, J.; Wright, S. J. Numerical Optimization, 2nd ed.; Springer: New York, 2006; Ch. 4.
/// More, J. J.; Sorensen, D. C. SIAM J. Sci. Stat. Comput. 1983, 4, 553.
/// Helmich-Paris, B. J. Chem. Phys. 2021, 154, 164104; J. Chem. Phys. 2022, 156, 204104.

#include <vector>
#include <cmath>
#include <algorithm>

#include "reks_math.h"

namespace psi {
namespace reks {

struct TRAHResult {
    std::vector<double> step;
    double predicted_decrease;      ///< -(g^T s + s^T H s / 2); positive = decrease
    double mu;                      ///< level shift (0 = Newton)
    int n_active_bounds;            ///< variables fixed at bounds
    bool step_on_boundary;          ///< |s| >= 0.99 trust_radius
};

/// Trust-region state carried between steps; predicted decrease = orb + sum_b (fon[b] + cross[b]) over FON blocks b.
struct TRAHState {
    double trust_radius = 0.4;
    double prev_energy = 0.0;      ///< E_SA at the previous step
    double prev_predicted = 0.0;
    double prev_pred_orb = 0.0;
    std::vector<double> prev_pred_fon;     ///< per FON block
    std::vector<double> prev_pred_cross;   ///< orbital x FON cross term per FON block
    bool initialized = false;      ///< prev_energy set
    double prev_rho = 1.0;         ///< actual / predicted decrease of the previous step
    bool prev_step_on_boundary = false;

    static constexpr double TR_MIN = 0.01;
    static constexpr double TR_MAX = 2.0;
};

/// Trust-region step solver, O(N^3) in the Hessian eigendecomposition; all methods static.
class TrustRegionSolver {
   public:
    /// Box-constrained trust-region step at x (H row-major N x N; +-1e30 bounds = free), radius from state.
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

        double mu_val = 0.0;
        solve_ah(g.data(), H.data(), N, state.trust_radius, result.step.data(), mu_val);
        result.mu = mu_val;

        // Clamp steps that leave [lower, upper].
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

        // Re-solve on the free subspace with g_eff[a] = g[ia] + sum_{fixed j} H[ia,j] step[j]
        // and L2 budget h_eff = sqrt(trust_radius^2 - sum_{fixed} step^2).
        if (n_fixed > 0 && n_fixed < N) {
            std::vector<int> free_idx;
            for (int i = 0; i < N; ++i) {
                if (!fixed[i]) free_idx.push_back(i);
            }
            int M = static_cast<int>(free_idx.size());

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

        result.predicted_decrease = predicted_decrease(g.data(), H.data(), result.step.data(), N);

        return result;
    }

    /// TR *= 2 (rho > 0.75 on the boundary), 1.2 (rho > 0.75), 1.05 (rho > 0.25), 0.7 (rho >= 0), 0.5 (rho < 0,
    /// step kept); stores actual_energy as prev_energy.
    static void update_trust_radius(double actual_energy, TRAHState& state) {
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

        state.trust_radius = std::clamp(state.trust_radius, TRAHState::TR_MIN, TRAHState::TR_MAX);
        state.prev_energy = actual_energy;
    }

    /// pred = -(g^T s + 0.5 s^T H s).
    static double predicted_decrease(const double* g, const double* H, const double* step, int N) {
        if (N == 0) return 0.0;
        // Full H s product: H is symmetric only up to rounding.
        std::vector<double> Hs(N, 0.0);
        C_DGEMV('N', N, N, 1.0, const_cast<double*>(H), N, const_cast<double*>(step), 1, 0.0,
                Hs.data(), 1);
        const double gs  = C_DDOT(N, const_cast<double*>(g), 1, const_cast<double*>(step), 1);
        const double sHs = C_DDOT(N, const_cast<double*>(step), 1, Hs.data(), 1);
        return -(gs + 0.5 * sHs);
    }

   private:
    /// Unconstrained TRS by level shift: step (caller-allocated) and mu (0 = Newton, else < lambda_min).
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
        // lambda ascending; V[k*N+i] = k-th component of eigenvector i.
        // gt = V^T g (transa='T': eigenvectors are columns of row-major V).
        std::vector<double> gt(N, 0.0);
        C_DGEMV('T', N, N, 1.0, V.data(), N, const_cast<double*>(g), 1, 0.0, gt.data(), 1);

        double lambda_min = lambda[0];

        // Newton step (mu = 0) if H > 0 and |s|^2 = sum_i gt_i^2 / lambda_i^2 <= h^2.
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

        // Bisection for |s(mu)| = h, |s(mu)|^2 = sum_i gt_i^2 / (lambda_i - mu)^2, increasing as mu -> lambda_min.
        auto step_norm_at = [&](double trial_mu) -> double {
            double norm_sq = 0.0;
            for (int i = 0; i < N; ++i) {
                double d = lambda[i] - trial_mu;
                norm_sq += gt[i] * gt[i] / (d * d);
            }
            return std::sqrt(norm_sq);
        };

        double mu_hi = lambda_min - 1e-8;
        double mu_lo = lambda_min - 10.0;

        // Widen mu_lo until |s(mu_lo)| <= h.
        for (int i = 0; i < 30 && step_norm_at(mu_lo) > h; ++i) {
            mu_lo = 2.0 * mu_lo - lambda_min;
        }

        // Hard case: g orthogonal to the lambda_min eigenvector, |s| <= h even at mu_hi.
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
                mu_hi = mu_mid;
            } else {
                mu_lo = mu_mid;
            }

            mu = mu_mid;

            // |s| within 0.1% of h.
            if (std::abs(sn / h - 1.0) < 1e-3) break;
        }

        build_step_from_eigdecomp(gt, lambda, V, N, mu, step);
    }

    /// s = -sum_i gt[i] / (lambda[i] - mu) V[:,i]
    static void build_step_from_eigdecomp(const std::vector<double>& gt, const std::vector<double>& lambda,
                                          const std::vector<double>& V, int N, double mu, double* step) {
        std::vector<double> coeff(N);
        for (int i = 0; i < N; ++i) coeff[i] = -gt[i] / (lambda[i] - mu);
        // s = V coeff (transa='N': eigenvector i is column i of row-major V).
        C_DGEMV('N', N, N, 1.0, const_cast<double*>(V.data()), N, coeff.data(), 1, 0.0, step, 1);
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_TRAH_SOLVER_H
