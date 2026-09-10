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

#include "reks_fon_relax.h"

#include "reks_fon_solver.h"

#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace psi {
namespace reks {
namespace fon_relax {

namespace {

constexpr double OUTER_TOL = 1e-7;

// Stays inside 2.0: scheme>0 C_K weight has a df_interp singularity at
// occupancy = 2.0 (dC/d2C undefined there). Gap (1e-8) is < OUTER_TOL.
constexpr double POSTSCF_UPPER_BOUND = 2.0 - 1e-8;

// Accessors into the generation-indexed snapshot/cassette/DiagDeps storage,
// keyed by `gen`. `name` is this side's FON channel display name.
struct SideOps {
    int  gen;
    std::string name;

    std::vector<studio::GeminalFON>& fons(studio::FONSnapshot& s) const {
        return s.layers[gen];
    }
    const std::vector<int>& active(const studio::Cassette& c) const {
        return c.geminals_active_gen(gen);
    }
    studio::IntSpan reads(const studio::DiagDepsView& d) const {
        return d.geminal_reads[gen];
    }
};

MicroSolverConfig make_microsolver_config(const Config1DNR& preset,
                                          bool use_line_search) {
    MicroSolverConfig mcfg;
    mcfg.max_nr_iter           = preset.max_iter;
    mcfg.gradient_tol          = preset.grad_tol;
    mcfg.energy_tol_factor     = preset.energy_tol_factor;
    mcfg.convergence_tol       = preset.check_step_convergence ? preset.step_tol : 0.0;
    mcfg.max_step              = preset.max_step;
    mcfg.lower_bound           = preset.lower_bound;
    mcfg.upper_bound           = POSTSCF_UPPER_BOUND;
    mcfg.boundary_masking      = preset.boundary_masking;
    mcfg.hess_threshold        = preset.hess_threshold;
    mcfg.reject_positive_hess  = preset.reject_positive_hess;
    mcfg.proportional_fallback = preset.proportional_fallback;
    mcfg.fallback_fraction     = preset.fallback_fraction;
    mcfg.use_line_search       = use_line_search;
    mcfg.collect_history       = report::reports(4);
    return mcfg;
}

// Primary-K: minimum global config index in si_cassette.K_indices that reads
// geminal g on the selected side. The catalog orders a geminal's reference
// config ahead of its excited mirrors; the minimum is invariant under pool
// reordering. Returns -1 if no config reads g (geminal inactive in this
// cassette).
int find_primary_K(const studio::Cassette&  si_cassette,
                   const SideOps&         side,
                   int                    g) {
    int primary = -1;
    for (int K : si_cassette.K_indices) {
        const auto& reads = side.reads(si_cassette.diag_deps(K));
        if (std::find(reads.begin(), reads.end(), g) != reads.end()) {
            if (primary < 0 || K < primary) primary = K;
        }
    }
    return primary;
}

// N-dimensional NR for one primary-K group (K, gs) on generation side.gen.
// Minimises E_K(x) = sum_L C_L^{K}(x) * E_L[L] over x in [LB, UB]^N where
// x[i] = side_fons[gs[i]].p.
// Multistart from 2^N corner-points + start point; keeps the lowest-energy
// result. Returns max_i |x_new[i] - x_old[i]| as the outer-GS termination delta.
//
// Hessian diag d2E/dx_i^2 via compute_weight_derivs, off-diag d2E/dx_i dx_j
// via compute_mixed_weight_derivs.
double run_nd_group(const studio::Cassette&      si_cassette,
                    int                        K,
                    const std::vector<int>&    gs,
                    const SideOps&             side,
                    const std::vector<double>& E_L,
                    const Config1DNR&          preset,
                    studio::FONSnapshot&         snap) {
    const int N       = static_cast<int>(gs.size());
    const int n_micro = static_cast<int>(E_L.size());
    // compute_C_L / weight_derivs touch only diag_deps[K].microstate_writes;
    // elsewhere outputs are 0. Sum only over that set.
    const studio::IntSpan Ls = si_cassette.config_microstate_writes(K);
    const int gen = side.gen;
    auto weight_derivs = [&si_cassette, gen](int Kc, int g,
                                             const studio::FONSnapshot& f,
                                             std::vector<double>& dC_,
                                             std::vector<double>& d2C_) {
        si_cassette.compute_weight_derivs(gen, Kc, g, f, dC_, d2C_);
    };
    auto mixed_derivs = [&si_cassette, gen](int Kc, int gi, int gj,
                                            const studio::FONSnapshot& f,
                                            std::vector<double>& d2Cm) {
        si_cassette.compute_mixed_weight_derivs(gen, Kc, gi, gj, f, d2Cm);
    };

    // Scratch (n_micro) for compute_C_L / mixed_derivs.
    std::vector<double> C_K(n_micro);
    std::vector<double> d2C_mixed(n_micro);

    auto write_x = [&side, &gs, &snap, N](const std::vector<double>& x) {
        auto& fons = side.fons(snap);
        for (int i = 0; i < N; ++i) fons[gs[i]] = { x[i], 2.0 - x[i] };
    };

    // dC[i]/d2C[i]: first/second derivative of C_K[L] wrt geminal gs[i].
    std::vector<std::vector<double>> dC(N, std::vector<double>(n_micro));
    std::vector<std::vector<double>> d2C(N, std::vector<double>(n_micro));
    std::vector<double> x_swept;
    // Recomputes dC/d2C only if x != x_swept; assumes write_x(x) already ran.
    auto sweep = [&](const std::vector<double>& x) {
        if (x_swept == x) return;
        for (int i = 0; i < N; ++i) weight_derivs(K, gs[i], snap, dC[i], d2C[i]);
        x_swept = x;
    };

    FONObjective obj;
    obj.ndim = N;
    obj.energy = [&si_cassette, K, &E_L, &snap, &Ls, write_x, &C_K]
                 (const std::vector<double>& x) {
        write_x(x);
        si_cassette.compute_C_L(K, snap, C_K);
        double E = 0.0;
        for (int L : Ls) E += C_K[L] * E_L[L];
        return E;
    };
    obj.gradient = [&, N]
                   (const std::vector<double>& x) {
        write_x(x);
        sweep(x);
        std::vector<double> grad(N, 0.0);
        for (int i = 0; i < N; ++i) {
            double dE = 0.0;
            for (int L : Ls) dE += dC[i][L] * E_L[L];
            grad[i] = dE;
        }
        return grad;
    };
    obj.hessian = [&, K, mixed_derivs, N]
                  (const std::vector<double>& x) {
        write_x(x);
        sweep(x);
        std::vector<double> H(N * N, 0.0);
        for (int i = 0; i < N; ++i) {
            double d2E = 0.0;
            for (int L : Ls) d2E += d2C[i][L] * E_L[L];
            H[i * N + i] = d2E;
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                mixed_derivs(K, gs[i], gs[j], snap, d2C_mixed);
                double d2E_ij = 0.0;
                for (int L : Ls) d2E_ij += d2C_mixed[L] * E_L[L];
                H[i * N + j] = d2E_ij;
                H[j * N + i] = d2E_ij;
            }
        }
        return H;
    };

    std::vector<double> x_old(N);
    {
        auto& fons = side.fons(snap);
        for (int i = 0; i < N; ++i) x_old[i] = fons[gs[i]].p;
    }

    const double lb = preset.lower_bound;
    const double ub = POSTSCF_UPPER_BOUND;

    // Multistart seeds: 2^N corners of [LB, UB]^N plus x_old; no corner sits
    // at the C_K singularity (see POSTSCF_UPPER_BOUND).
    std::vector<std::vector<double>> seeds;
    const int n_corners = 1 << N;
    seeds.reserve(static_cast<size_t>(n_corners) + 1);
    for (int mask = 0; mask < n_corners; ++mask) {
        std::vector<double> s(N);
        for (int i = 0; i < N; ++i) s[i] = ((mask >> i) & 1) ? ub : lb;
        seeds.push_back(std::move(s));
    }
    seeds.push_back(x_old);

    // 1D uses raw NR + clamp; 2D+ uses line-search.
    auto mcfg = make_microsolver_config(preset, /*use_line_search=*/(N >= 2));

    auto fmt_g_list = [&gs, N]() {
        std::string s;
        for (int i = 0; i < N; ++i) {
            if (i) s += ",";
            s += std::to_string(gs[i]);
        }
        return s;
    };
    // Comma-joined FON vector at .8f -- resolves the near-boundary regime (UB = 2 - 1e-8).
    auto fmt_vec = [](const std::vector<double>& v) {
        std::string s;
        char buf[32];
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s += ",";
            std::snprintf(buf, sizeof(buf), "%.8f", v[i]);
            s += buf;
        }
        return s;
    };

    double best_E = std::numeric_limits<double>::infinity();
    std::vector<double> best_x = x_old;
    for (const auto& seed : seeds) {
        FONMicroSolver solver(N);
        auto result = solver.solve(obj, seed, mcfg);
        if (result.energy < best_E) {
            best_E = result.energy;
            best_x = result.fon;
        }
        if (report::reports(4)) {
            const std::string seed_str = fmt_vec(seed);
            const std::string g_str = fmt_g_list();
            for (const auto& it : result.history) {
                outfile->Printf("  fon_relax %s K=%d g=[%s] seed=[%s] "
                                "iter=%d fon=[%s] energy=%.10f grad_norm=%.2e step_norm=%.2e\n",
                                side.name.c_str(), K, g_str.c_str(), seed_str.c_str(),
                                it.iter, fmt_vec(it.fon).c_str(), it.energy,
                                it.grad_norm, it.step_norm);
            }
            outfile->Printf("  fon_relax %s K=%d g=[%s] seed=[%s] "
                            "FINAL iter=%d fon=[%s] energy=%.10f conv=%d\n",
                            side.name.c_str(), K, g_str.c_str(), seed_str.c_str(),
                            result.nr_iterations, fmt_vec(result.fon).c_str(),
                            result.energy, result.converged ? 1 : 0);
        }
    }
    {
        std::vector<double> tmp = best_x;
        write_x(tmp);
    }
    double max_delta = 0.0;
    for (int i = 0; i < N; ++i) {
        max_delta = std::max(max_delta, std::abs(best_x[i] - x_old[i]));
    }
    return max_delta;
}

// One side (generation): groups its SI-only geminals by primary K, then outer
// Gauss-Seidel sweeps the groups until max_delta < OUTER_TOL. No iteration
// cap: block-coordinate descent with a descent guarantee converges.
void optimize_side(const studio::Cassette&      si_cassette,
                   const studio::Cassette&      sa_cassette,
                   int                        sa_sector,
                   const SideOps&             side,
                   const std::vector<double>& E_L,
                   studio::FONSnapshot&         snap) {
    std::map<int, std::vector<int>> groups;
    // Frozen set = this sector's SA active geminals only; a geminal SA-active in
    // another sector must still relax here if it is SI-only for this sector.
    const std::vector<int>& sa_active = sa_cassette.geminals_active(sa_sector, side.gen);
    for (int g : side.active(si_cassette)) {
        if (std::find(sa_active.begin(), sa_active.end(), g) != sa_active.end()) continue;
        int K = find_primary_K(si_cassette, side, g);
        if (K < 0) continue;
        groups[K].push_back(g);
    }
    if (groups.empty()) return;

    auto preset = Config1DNR::post_scf_preset();
    for (int outer = 0; ; ++outer) {
        double max_delta = 0.0;
        for (const auto& kv : groups) {
            max_delta = std::max(max_delta,
                run_nd_group(si_cassette, kv.first, kv.second, side,
                             E_L, preset, snap));
        }
        if (report::reports(4)) {
            outfile->Printf("  fon_relax %s outer=%d max_delta=%.2e\n",
                            side.name.c_str(), outer, max_delta);
        }
        if (max_delta < OUTER_TOL) break;
    }
}

}  // namespace

Result optimize(const studio::Cassette&      si_cassette,
                           const studio::Cassette&      sa_cassette,
                           int                        sa_sector,
                           const studio::FONSnapshot&   fons_initial,
                           const std::vector<double>& E_L) {
    studio::FONSnapshot snap = fons_initial;
    // Relax every FON generation the SI cassette carries (n, m, u, ...): each is
    // an independent block-coordinate side, all driven by the same N-dim NR.
    for (int gen = 0; gen < si_cassette.n_generations(); ++gen) {
        SideOps side{ gen, studio::fon_channel_name(gen, si_cassette.sector()) };
        optimize_side(si_cassette, sa_cassette, sa_sector, side, E_L, snap);
    }
    return Result{ std::move(snap.layers) };
}

}  // namespace fon_relax
}  // namespace reks
}  // namespace psi
