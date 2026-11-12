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

#include "reks_fon_step.h"

#include "reks_math.h"
#include "reks_objectives.h"
#include "reks_report_level.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace psi {
namespace reks {

using report::reports;

void FonStep::configure(const MicroSolverConfig& cfg, double boundary_skin) {
    cfg_ = cfg;
    boundary_skin_ = boundary_skin;
}

void FonStep::relax(FonStepContext& ctx) {
    timer_on("REKS: gvb_fon_step");

    const int n_sec = ctx.sa_cassette.n_sectors();
    const int n_gen = ctx.sa_cassette.n_generations();
    const int P = static_cast<int>(ctx.sa_cassette.geminals_active_gen(0).size());
    const auto* templ = ctx.sa_cassette.geminal_templates();
    const auto& gen0_pool = ctx.sa_cassette.geminals_active_gen(0);

    // The swap dedup is per orbital pair, so each gen-0 geminal must own a distinct
    // pair (union of every sector's gen-0 geminals; templates are M-scoped shared).
#ifndef NDEBUG
    {
        std::set<std::pair<int, int>> seen_pairs;
        for (int p = 0; p < P; ++p) {
            const auto& o = templ[gen0_pool[p]].orbitals;
            auto key = std::minmax(o[0], o[1]);
            assert(seen_pairs.emplace(key.first, key.second).second &&
                   "gen-0 geminal <-> orbital-pair is not 1:1");
        }
    }
#endif

    // Relabel of a swapped pair: FON -> 2 - FON for every (sector, generation) geminal on it (columns are
    // shared); generations >= 1 report the new value.
    auto reflect_pair = [&](int orb0, int orb1, int pair_idx) -> int {
        int count = 0;
        for (int sec = 0; sec < n_sec; ++sec) {
            for (int gen = 0; gen < n_gen; ++gen) {
                for (int gg : ctx.sa_cassette.geminals_active(sec, gen)) {
                    const auto& go = templ[gg].orbitals;
                    if (!((go[0] == orb0 && go[1] == orb1) ||
                          (go[0] == orb1 && go[1] == orb0))) continue;
                    const double f_new = 2.0 - ctx.reel.fon_state[sec].layers[gen][gg].p;
                    ctx.reel.fon_state[sec].layers[gen][gg] = {f_new, 2.0 - f_new};
                    ++count;
                    if (gen >= 1 && reports(4))
                        outfile->Printf("  [GVB-FON] iter=%d %s-geminal %d synced (-> %.6f) with pair %d swap\n",
                                        ctx.iteration, reks::studio::fon_channel_name(gen, sec).c_str(),
                                        gg, f_new, pair_idx);
                }
            }
        }
        return count;
    };

    if (!initialized_) {
#ifndef NDEBUG
        // Snapshot to check that a single-sector fgap swap leaves every FON unchanged.
        std::vector<std::vector<double>> dbg_fon_pre;
        if (n_sec == 1) dbg_fon_pre = studio::capture_fons(ctx.reel, ctx.sa_cassette, ctx.fon_blocks);
#endif
        // Pre-NR swap of bonding/antibonding columns where eps[mo_q] < eps[mo_p]; multi-sector runs reflect
        // every sector's FONs with it.
        bool any_swap = false;
        if (ctx.epsilon) {
            double* eps = ctx.epsilon->pointer(0);
            for (int p = 0; p < P; ++p) {
                int g = gen0_pool[p];
                const auto& orbs = templ[g].orbitals;
                int i_idx = ctx.active[orbs[0]];
                int j_idx = ctx.active[orbs[1]];
                double gap = eps[j_idx] - eps[i_idx];
                if (gap < 0.0) {
                    if (reports(4)) {
                        outfile->Printf(
                            "  [GVB-FON] iter=%d INV-SWAP pair %d (fgap=%+.6f -> swap MO %d <-> %d)\n",
                            ctx.iteration, p, gap, i_idx, j_idx);
                    }
                    ctx.Ca->swap_columns(0, i_idx, j_idx);
                    double e_tmp = eps[i_idx];
                    eps[i_idx] = eps[j_idx];
                    eps[j_idx] = e_tmp;
                    ctx.swapped = true;
                    any_swap = true;
                    if (n_sec > 1) reflect_pair(orbs[0], orbs[1], p);
                }
            }
        }
#ifndef NDEBUG
        if (n_sec == 1)
            assert(studio::capture_fons(ctx.reel, ctx.sa_cassette, ctx.fon_blocks) == dbg_fon_pre &&
                   "fgap swap reflected a FON in a single-sector run");
#endif

        // No active n-geminal (e.g. OSS-only pool): no n-FON to solve, and no zero-size DSYEV.
        if (P == 0) {
            initialized_ = true;
            timer_off("REKS: gvb_fon_step");
            return;
        }

        // Per-sector gen-0 seed: multi-start (all-1.0 vs GVB-PP CI), preferring the
        // interior minimum over a boundary basin.
        for (int sec = 0; sec < n_sec; ++sec) {
            const int P_s = static_cast<int>(ctx.sa_cassette.geminals_active(sec, 0).size());
            if (P_s == 0) continue;

            if (reports(4) && n_sec > 1)
                outfile->Printf("  [GVB-FON] iter=%d sector %d init\n", ctx.iteration, sec);

            auto snap = std::make_shared<reks::studio::FONSnapshot>(ctx.reel.fon_state[sec]);
            auto obj  = reks::build_geminal_objective(sec, 0, ctx.sa_cassette, snap, ctx.reel.E_L);
            reks::FONMicroSolver init_solver(P_s);

            auto is_boundary = [&](const std::vector<double>& f) {
                for (int p = 0; p < P_s; ++p)
                    if (std::min(f[p], 2.0 - f[p]) <= boundary_skin_) return true;
                return false;
            };

            std::vector<std::vector<double>> starts;
            starts.push_back(std::vector<double>(P_s, 1.0));
            auto fon_ci = reks::fon_init_from_ci(sec, ctx.sa_cassette, ctx.reel.E_L, ctx.iteration,
                                                 /*narrate=*/false);
            if (static_cast<int>(fon_ci.size()) == P_s) {
                bool different = false;
                for (int p = 0; p < P_s; ++p)
                    if (std::abs(fon_ci[p] - 1.0) > 1e-6) { different = true; break; }
                if (different) starts.push_back(fon_ci);
            }

            std::vector<double> best_interior_fon;
            double best_interior_E = std::numeric_limits<double>::infinity();
            std::vector<double> best_any_fon = starts[0];
            double best_any_E = std::numeric_limits<double>::infinity();
            int best_idx_any = -1, best_idx_int = -1;
            for (size_t st = 0; st < starts.size(); ++st) {
                auto res = init_solver.solve(obj, starts[st], cfg_);
                double E_try = obj.energy(res.fon);
                bool boundary = is_boundary(res.fon);
                if (reports(4)) {
                    outfile->Printf("  [GVB-FON] iter=%d multi-start %zu: start=(", ctx.iteration, st);
                    for (int p = 0; p < P_s; ++p)
                        outfile->Printf("%.4f%s", starts[st][p], p + 1 < P_s ? "," : "");
                    outfile->Printf(") -> FON=(");
                    for (int p = 0; p < P_s; ++p)
                        outfile->Printf("%.6f%s", res.fon[p], p + 1 < P_s ? "," : "");
                    outfile->Printf(") E=%.10f%s\n", E_try, boundary ? " [boundary]" : "");
                }
                if (E_try < best_any_E) {
                    best_any_E = E_try;
                    best_any_fon = res.fon;
                    best_idx_any = static_cast<int>(st);
                }
                if (!boundary && E_try < best_interior_E) {
                    best_interior_E = E_try;
                    best_interior_fon = res.fon;
                    best_idx_int = static_cast<int>(st);
                }
            }
            std::vector<double> best_fon = best_interior_fon.empty() ? best_any_fon : best_interior_fon;
            double best_E = best_interior_fon.empty() ? best_any_E : best_interior_E;
            int best_idx = best_interior_fon.empty() ? best_idx_any : best_idx_int;
            if (reports(4) && starts.size() > 1) {
                outfile->Printf("  [GVB-FON] iter=%d multi-start picked branch %d (E=%.10f, %s)\n",
                                ctx.iteration, best_idx, best_E,
                                best_interior_fon.empty() ? "no interior found, using best E" : "interior preferred");
            }

            studio::set_fon_vector(ctx.reel, ctx.sa_cassette, sec, 0, best_fon);
            if (reports(4)) {
                outfile->Printf("  [GVB-FON] iter=%d init: FON=(", ctx.iteration);
                for (int p = 0; p < P_s; ++p)
                    outfile->Printf("%.6f%s", best_fon[p], p + 1 < P_s ? "," : "");
                outfile->Printf(")%s\n", any_swap ? " (after fgap inv-swap)" : "");
            }
        }

        initialized_ = true;
        timer_off("REKS: gvb_fon_step");
        return;
    }

    constexpr double fon_bnd_margin = 1e-8;
    // Per-generation FON floor: REKS_<G>_FON if set (>= 0), else open lower 1e-8.
    auto fon_lower = [&](int gen) {
        return ctx.fon_lower[gen] >= 0.0 ? ctx.fon_lower[gen] : 1e-8;
    };

    // gen-0 crossing record, keyed by the aggregate pair position (union of all
    // sectors' gen-0 pairs); ORs across sectors that share the pair.
    std::vector<int> g0_to_pos(ctx.sa_cassette.n_geminals(), -1);
    for (int p = 0; p < P; ++p) g0_to_pos[gen0_pool[p]] = p;
    std::vector<char> crossed(P, 0);
    std::vector<double> swap_before(P, 0.0), swap_after(P, 0.0);

    // Newton-Raphson per (sector, generation) block against its sector's snapshot; the FON Hessian is
    // block-diagonal, so the order is free.
    for (int sec = 0; sec < n_sec; ++sec) {
        if (reports(4) && n_sec > 1)
            outfile->Printf("  [GVB-FON] iter=%d sector %d\n", ctx.iteration, sec);
        for (int gen = 0; gen < n_gen; ++gen) {
            const auto& gem = ctx.sa_cassette.geminals_active(sec, gen);
            const int Pg = static_cast<int>(gem.size());
            if (Pg == 0) continue;  // no active geminals: skip the zero-dim DSYEV.

            const double lo = fon_lower(gen);
            std::vector<double> x0 = studio::fon_vector(ctx.reel, ctx.sa_cassette, sec, gen);

            auto snap = std::make_shared<reks::studio::FONSnapshot>(ctx.reel.fon_state[sec]);
            auto obj  = reks::build_geminal_objective(sec, gen, ctx.sa_cassette, snap, ctx.reel.E_L);

            reks::FONMicroSolver solver(Pg);
            reks::MicroSolverConfig cfg = cfg_;
            cfg.lower_bound = lo;
            cfg.collect_history = reports(5);

            auto res = solver.solve(obj, x0, cfg);

            // Snap to the boundary only on a fresh crossing.
            for (int k = 0; k < Pg; ++k) {
                if (res.fon[k] >= 2.0 - fon_bnd_margin && x0[k] < 2.0 - fon_bnd_margin)
                    res.fon[k] = 2.0 - fon_bnd_margin;
                if (res.fon[k] <= lo && x0[k] > lo)
                    res.fon[k] = lo;
            }

            studio::set_fon_vector(ctx.reel, ctx.sa_cassette, sec, gen, res.fon);

            // gen-0 FON crossing 1.0 flags the shared pair for an orbital swap.
            if (gen == 0) {
                for (int k = 0; k < Pg; ++k) {
                    const bool cr =
                        (x0[k] > 1.0 && res.fon[k] < 1.0) ||
                        (x0[k] < 1.0 && res.fon[k] > 1.0);
                    if (!cr) continue;
                    const int p = g0_to_pos[gem[k]];
                    crossed[p] = 1;
                    swap_before[p] = x0[k];
                    swap_after[p] = res.fon[k];
                }
            }

            if (reports(4)) {
                const std::string chan = reks::studio::fon_channel_name(gen, sec);
                outfile->Printf("  [GVB-FON] iter=%d %s-NR_iters=%d %s-FON=(",
                                ctx.iteration, chan.c_str(), res.nr_iterations, chan.c_str());
                for (int k = 0; k < Pg; ++k)
                    outfile->Printf("%.8f%s", res.fon[k], k + 1 < Pg ? ", " : "");
                outfile->Printf(") %s-prev=(", chan.c_str());
                for (int k = 0; k < Pg; ++k)
                    outfile->Printf("%.8f%s", x0[k], k + 1 < Pg ? ", " : "");
                outfile->Printf(")\n");
            }
        }
    }

    // A gen-0 FON crossing 1.0 swaps its pair's columns once and reflects every FON on the pair.
    for (int p = 0; p < P; ++p) {
        if (!crossed[p]) continue;

        if (reports(4)) {
            outfile->Printf("  [GVB-FON] iter=%d swap pair %d (FON %.6f -> %.6f)\n",
                            ctx.iteration, p, swap_before[p], swap_after[p]);
        }
        const int g0 = gen0_pool[p];
        const auto& orbs = templ[g0].orbitals;
        const int i_idx = ctx.active[orbs[0]];
        const int j_idx = ctx.active[orbs[1]];
        ctx.Ca->swap_columns(0, i_idx, j_idx);
        ctx.swapped = true;
        reflect_pair(orbs[0], orbs[1], p);
    }

    timer_off("REKS: gvb_fon_step");
}

}  // namespace reks
}  // namespace psi
