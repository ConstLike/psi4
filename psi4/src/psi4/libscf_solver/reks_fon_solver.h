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

#ifndef REKS_FON_SOLVER_H
#define REKS_FON_SOLVER_H

#include "reks_report_level.h"
#include "reks_math.h"
#include "psi4/libqt/qt.h"

#include <functional>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include <numeric>

namespace psi {
namespace reks {

/// FON optimization objective. Each callback sets FON state from its trial
/// point internally before evaluating; the solver holds no FON state.
struct FONObjective {
    int ndim;
    std::function<double(const std::vector<double>&)> energy;
    std::function<std::vector<double>(const std::vector<double>&)> gradient;
    std::function<std::vector<double>(const std::vector<double>&)> hessian;
};


struct Config1DNR {
    double lower_bound = 0.0;
    double grad_tol = 1e-10;
    double energy_tol_factor = 0.0;
    double hess_threshold = 1e-10;
    bool reject_positive_hess = false;
    bool boundary_masking = false;
    double max_step = 0.5;
    int max_iter = 50;
    double step_tol = 0.0;
    bool check_step_convergence = false;
    bool proportional_fallback = false;
    double fallback_fraction = 0.5;

    static Config1DNR post_scf_preset() {
        Config1DNR c;
        c.lower_bound = 1.0;
        c.grad_tol = 1e-10;
        c.energy_tol_factor = 10.0;
        c.hess_threshold = 1e-12;
        c.reject_positive_hess = false;
        c.boundary_masking = true;
        c.max_step = 0.5;
        c.max_iter = 50;
        c.check_step_convergence = false;
        c.proportional_fallback = true;
        c.fallback_fraction = 0.1;
        return c;
    }
};

struct MicroSolverIterInfo {
    int iter;
    double energy;
    double grad_norm;
    double step_norm;
    std::vector<double> fon;
};

struct MicroSolverConfig {
    int max_nr_iter = 50;
    double lower_bound = 0.0;
    double upper_bound = 2.0;
    bool collect_history = false;    ///< fills MicroSolverResult::history

    double gradient_tol = 0.0;       ///< ||grad|| (0 = disabled)
    double convergence_tol = 1e-12;  ///< ||dx|| (0 = disabled)
    double energy_tol_factor = 0.0;  ///< predicted |dE| < factor*eps*max(1,|E|) (0 = disabled)

    double hess_threshold = 1e-10;
    bool reject_positive_hess = false;
    bool proportional_fallback = false;  ///< true = -fraction*grad; false = fixed-size SD
    double fallback_fraction = 0.5;
    double max_step = 0.0;               ///< 0 = unlimited
    bool boundary_masking = false;

    bool use_line_search = true;         ///< false = raw NR + clamp
    int max_ls_iter = 20;
    double min_eigenvalue = 0.01;        ///< LS-mode eigenvalue floor
    bool two_pass_ls = true;             ///< Try negative direction on failure
};

struct MicroSolverResult {
    std::vector<double> fon;
    double energy = 0.0;
    bool converged = false;
    int nr_iterations = 0;
    std::vector<MicroSolverIterInfo> history;  ///< populated if cfg.collect_history
};

/// N-dim Newton-Raphson + golden-ratio line search FON solver: drives the
/// FON state to local-min within a single call.
class FONMicroSolver {
public:
    explicit FONMicroSolver(int ndim) : ndim_(ndim) {}

    // Single-start Newton descent on the box [LB, UB]^ndim. Per iter:
    //   g = grad(x),  dx = -H_reg^{-1} g  (eigenvalue-regularized, or raw NR),
    //   x <- clamp(x + alpha dx)          (alpha from line search, or alpha = 1).
    // Converges on masked ||g|| < gradient_tol, ||dx|| < convergence_tol, or
    // predicted |dE| < energy_tol_factor*eps*max(1,|E|); an exact zero step
    // exits non-converged.
    MicroSolverResult solve(const FONObjective& obj,
                            const std::vector<double>& start,
                            const MicroSolverConfig& cfg) {
        MicroSolverResult result;
        result.fon = start;
        clamp(result.fon, cfg.lower_bound, cfg.upper_bound);

        double E_current = obj.energy(result.fon);
        result.energy = E_current;

        for (int iter = 0; iter < cfg.max_nr_iter; ++iter) {
            auto grad = obj.gradient(result.fon);

            // Boundary masking: convergence check only; Newton step uses unmasked grad.
            double gnorm;
            if (cfg.boundary_masking) {
                double masked_sq = 0.0;
                for (int i = 0; i < ndim_; ++i) {
                    double g = grad[i];
                    if (result.fon[i] >= cfg.upper_bound - 1e-10 && g < 0.0) g = 0.0;
                    if (result.fon[i] <= cfg.lower_bound + 1e-10 && g > 0.0) g = 0.0;
                    masked_sq += g * g;
                }
                gnorm = std::sqrt(masked_sq);
            } else {
                gnorm = vec_norm(grad);
            }

            if (cfg.gradient_tol > 0.0 && gnorm < cfg.gradient_tol) {
                result.energy = E_current;
                result.converged = true;
                result.nr_iterations = iter;
                return result;
            }

            // Hessian only needed for the Newton step; skip it on the converged iter.
            auto hess = obj.hessian(result.fon);

            std::vector<double> dx;
            if (cfg.use_line_search) {
                dx = compute_newton_step(grad, hess, cfg.min_eigenvalue);
            } else {
                dx = compute_raw_nr_step(grad, hess, cfg);
            }

            if (cfg.max_step > 0.0) {
                double dx_norm = vec_norm(dx);
                if (dx_norm > cfg.max_step) {
                    double scale = cfg.max_step / dx_norm;
                    for (int i = 0; i < ndim_; ++i) dx[i] *= scale;
                }
            }

            // dE_lin = sum_i grad[i] * (clamp(x_i+dx_i, LB, UB) - x_i): first-order
            // energy change of the box-clamped trial step.
            if (cfg.energy_tol_factor > 0.0) {
                double dE_lin = 0.0;
                for (int i = 0; i < ndim_; ++i) {
                    const double xi = std::max(cfg.lower_bound,
                                               std::min(cfg.upper_bound,
                                                        result.fon[i] + dx[i]));
                    dE_lin += grad[i] * (xi - result.fon[i]);
                }
                const double dE_tol = cfg.energy_tol_factor *
                                      std::numeric_limits<double>::epsilon() *
                                      std::max(1.0, std::abs(E_current));
                if (std::abs(dE_lin) < dE_tol) {
                    result.energy = E_current;
                    result.converged = true;
                    result.nr_iterations = iter;
                    return result;
                }
            }

            double actual_step_norm;
            if (cfg.use_line_search) {
                auto ls = line_search(obj, result.fon, dx, E_current, cfg);
                result.fon = ls.fon;
                E_current = ls.energy;
                actual_step_norm = ls.actual_step_norm;
            } else {
                std::vector<double> x_prev = result.fon;
                for (int i = 0; i < ndim_; ++i) result.fon[i] += dx[i];
                clamp(result.fon, cfg.lower_bound, cfg.upper_bound);
                double sq = 0.0;
                for (int i = 0; i < ndim_; ++i) {
                    double d = result.fon[i] - x_prev[i];
                    sq += d * d;
                }
                actual_step_norm = std::sqrt(sq);
                E_current = obj.energy(result.fon);
            }

            if (cfg.collect_history) {
                MicroSolverIterInfo info;
                info.iter = iter;
                info.energy = E_current;
                info.grad_norm = gnorm;
                info.step_norm = actual_step_norm;
                info.fon = result.fon;
                result.history.push_back(info);
            }

            if (cfg.convergence_tol > 0.0 && actual_step_norm > 0.0 &&
                actual_step_norm <= cfg.convergence_tol) {
                result.energy = E_current;
                result.converged = true;
                result.nr_iterations = iter + 1;
                return result;
            }

            if (actual_step_norm == 0.0) {
                result.energy = E_current;
                result.converged = false;
                result.nr_iterations = iter + 1;
                return result;
            }
        }

        result.energy = E_current;
        result.converged = false;
        result.nr_iterations = cfg.max_nr_iter;
        return result;
    }

private:
    int ndim_;

    static constexpr double GOLDEN = 1.6180339887498948482;
    static constexpr double GOLDEN_SQ = GOLDEN * GOLDEN;

    struct LSResult {
        std::vector<double> fon;
        double energy;
        double actual_step_norm;
    };

    /// Regularized Newton dx = -H_reg^{-1} g via DSYEV; eigenvalues clamped to
    /// min_eigenvalue. Scalar fast path for N=1.
    std::vector<double> compute_newton_step(
        const std::vector<double>& grad,
        const std::vector<double>& hess,
        double min_eigenvalue) const {

        if (ndim_ == 1) {
            double h = hess[0];
            if (h < min_eigenvalue) h = min_eigenvalue;
            return {-grad[0] / h};
        }

        std::vector<double> V(ndim_ * ndim_);
        for (int i = 0; i < ndim_ * ndim_; ++i) V[i] = hess[i];

        std::vector<double> w(ndim_);
        int lwork = std::max(3 * ndim_ - 1, 1);
        std::vector<double> work(lwork);
        int info = C_DSYEV('V', 'L', ndim_, V.data(), ndim_, w.data(), work.data(), lwork);

        // DSYEV failure: fixed-length steepest descent.
        if (info != 0) {
            double gnorm = vec_norm(grad);
            if (gnorm < 1e-14) return std::vector<double>(ndim_, 0.0);
            std::vector<double> step(ndim_);
            double fallback_step = 0.01;
            for (int i = 0; i < ndim_; ++i) {
                step[i] = -fallback_step * grad[i] / gnorm;
            }
            return step;
        }

        for (int i = 0; i < ndim_; ++i) {
            if (w[i] < min_eigenvalue) w[i] = min_eigenvalue;
        }

        // dx = V diag(-1/w) V^T g in the eigenbasis: g_e = V^T g, dx_e = -g_e/w, dx = V dx_e.
        std::vector<double> g_e(ndim_, 0.0);
        for (int i = 0; i < ndim_; ++i) {
            for (int j = 0; j < ndim_; ++j) {
                g_e[i] += V[j + i * ndim_] * grad[j];
            }
        }

        std::vector<double> dx_e(ndim_);
        for (int i = 0; i < ndim_; ++i) {
            dx_e[i] = -g_e[i] / w[i];
        }

        std::vector<double> dx(ndim_, 0.0);
        for (int i = 0; i < ndim_; ++i) {
            for (int j = 0; j < ndim_; ++j) {
                dx[i] += V[i + j * ndim_] * dx_e[j];
            }
        }

        return dx;
    }

    /// Raw Newton dx = -H^{-1} g; falls back to gradient step on singular H or
    /// (with reject_positive_hess) convex H.
    std::vector<double> compute_raw_nr_step(
        const std::vector<double>& grad,
        const std::vector<double>& hess,
        const MicroSolverConfig& cfg) const {

        std::vector<double> V(ndim_ * ndim_);
        for (int i = 0; i < ndim_ * ndim_; ++i) V[i] = hess[i];

        std::vector<double> w(ndim_);
        int lwork = std::max(3 * ndim_ - 1, 1);
        std::vector<double> work(lwork);
        int info = C_DSYEV('V', 'L', ndim_, V.data(), ndim_, w.data(), work.data(), lwork);

        if (info != 0) {
            return compute_fallback_step(grad, cfg);
        }

        bool use_fallback = false;
        bool all_positive = true;
        for (int i = 0; i < ndim_; ++i) {
            if (std::abs(w[i]) < cfg.hess_threshold) use_fallback = true;
            if (w[i] <= 0.0) all_positive = false;
        }
        if (cfg.reject_positive_hess && all_positive) use_fallback = true;

        if (use_fallback) {
            return compute_fallback_step(grad, cfg);
        }

        // dx = V diag(-1/w) V^T g in the eigenbasis: g_e = V^T g, dx_e = -g_e/w, dx = V dx_e.
        std::vector<double> g_e(ndim_, 0.0);
        for (int i = 0; i < ndim_; ++i)
            for (int j = 0; j < ndim_; ++j)
                g_e[i] += V[j + i * ndim_] * grad[j];

        std::vector<double> dx_e(ndim_);
        for (int i = 0; i < ndim_; ++i)
            dx_e[i] = -g_e[i] / w[i];

        std::vector<double> dx(ndim_, 0.0);
        for (int i = 0; i < ndim_; ++i)
            for (int j = 0; j < ndim_; ++j)
                dx[i] += V[i + j * ndim_] * dx_e[j];

        return dx;
    }

    /// Gradient fallback: proportional (-fraction*grad) or fixed-size SD.
    std::vector<double> compute_fallback_step(
        const std::vector<double>& grad,
        const MicroSolverConfig& cfg) const {

        if (cfg.proportional_fallback) {
            std::vector<double> dx(ndim_);
            for (int i = 0; i < ndim_; ++i) dx[i] = -cfg.fallback_fraction * grad[i];
            return dx;
        }
        double gnorm = vec_norm(grad);
        if (gnorm < 1e-10) return std::vector<double>(ndim_, 0.0);
        std::vector<double> dx(ndim_);
        double step = cfg.max_step * cfg.fallback_fraction;
        for (int i = 0; i < ndim_; ++i) dx[i] = -step * grad[i] / gnorm;
        return dx;
    }

    /// Golden-ratio backtracking LS, factor /= GOLDEN^2 (accelerated past midpoint).
    /// Two-pass tries the negative direction if positive fails.
    LSResult line_search(const FONObjective& obj,
                         const std::vector<double>& x_current,
                         const std::vector<double>& dx,
                         double E_current,
                         const MicroSolverConfig& cfg) const {

        double dx_norm = vec_norm(dx);

        int n_passes = cfg.two_pass_ls ? 2 : 1;

        for (int pass = 0; pass < n_passes; ++pass) {
            double sign = (pass == 0) ? 1.0 : -1.0;
            double factor = sign;

            for (int ls = 0; ls < cfg.max_ls_iter; ++ls) {
                std::vector<double> x_trial(ndim_);
                for (int i = 0; i < ndim_; ++i) {
                    x_trial[i] = x_current[i] + factor * dx[i];
                }
                clamp(x_trial, cfg.lower_bound, cfg.upper_bound);

                double E_trial = obj.energy(x_trial);

                if (E_trial < E_current) {
                    double actual_norm = 0.0;
                    for (int i = 0; i < ndim_; ++i) {
                        double d = x_trial[i] - x_current[i];
                        actual_norm += d * d;
                    }
                    return {x_trial, E_trial, std::sqrt(actual_norm)};
                }

                double ls_bailout = (cfg.convergence_tol > 0.0) ? cfg.convergence_tol : 1e-14;
                if (std::abs(factor) * dx_norm < ls_bailout) {
                    break;
                }

                factor /= GOLDEN_SQ;
                if (ls > cfg.max_ls_iter / 2) {
                    factor /= GOLDEN_SQ;
                }
            }
        }

        return {x_current, E_current, 0.0};
    }

    static void clamp(std::vector<double>& x, double lb, double ub) {
        for (double& v : x) {
            v = std::max(lb, std::min(ub, v));
        }
    }

    static double vec_norm(const std::vector<double>& v) {
        double s = 0.0;
        for (double x : v) s += x * x;
        return std::sqrt(s);
    }
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_FON_SOLVER_H
