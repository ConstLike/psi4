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

#include "reks_objectives.h"

#include "reks_gradient_engine.h"
#include "reks_math.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace psi {
namespace reks {

// File-local FON evaluators. They read the FONs the snapshot holds; the trial
// values must be written into it before any call.
namespace {

// Run sector s's partial SA energy at a trial FON snapshot, over all generations
// the snapshot carries:
//   E = sum_L C_L E_L,  C_L = sum_{K in sector s} w_K C_L^(K)(fons)
double geminal_target_energy(int s,
                             const studio::Cassette& sa_cassette,
                             const studio::FONSnapshot& fons,
                             const std::vector<double>& E_L,
                             std::vector<double>& C_K_buf,
                             std::vector<double>& C_L_buf) {
    const int n_micro = static_cast<int>(E_L.size());
    const auto& positions = sa_cassette.sector_config_positions(s);

    if (static_cast<int>(C_K_buf.size()) < n_micro) C_K_buf.resize(n_micro);
    C_L_buf.assign(n_micro, 0.0);
    for (int pos : positions) {
        const int    K   = sa_cassette.K_indices[pos];
        const double w_K = sa_cassette.K_weights[pos];
        sa_cassette.compute_C_L(K, fons, C_K_buf);
        for (int L : sa_cassette.diag_deps(K).microstate_writes)
            C_L_buf[L] += w_K * C_K_buf[L];
    }
    const int n_catalog = sa_cassette.n_catalog_microstates();
    const int n_extra   = sa_cassette.n_extra_microstates();
    for (int e = 0; e < n_extra; ++e)
        C_L_buf[n_catalog + e] = sa_cassette.extra_weights[e];
    double E = 0.0;
    for (int L : sa_cassette.microstates_active) E += C_L_buf[L] * E_L[L];
    return E;
}

// Weight-derivative sweep cached at its trial point x.
struct WeightDerivsMemo {
    std::vector<double>              x;
    REKSGradientEngine::WeightDerivs wd;
    bool                             valid = false;
};

}  // anonymous namespace

FONObjective
build_geminal_objective(int s,
                        int gen,
                        const studio::Cassette& sa_cassette,
                        std::shared_ptr<studio::FONSnapshot> snap,
                        const std::vector<double>& E_L) {
    FONObjective obj;
    obj.ndim = static_cast<int>(sa_cassette.geminals_active(s, gen).size());

    // Energy-evaluator scratch; shared_ptr keeps it alive in the returned
    // lambdas. Gradient/hessian use the engine's own scratch.
    const int n_micro = sa_cassette.n_catalog_microstates();
    auto C_K_buf = std::make_shared<std::vector<double>>(n_micro);
    auto C_L_buf = std::make_shared<std::vector<double>>(n_micro);

    auto write = [snap, &sa_cassette, s, gen](const std::vector<double>& x) {
        const auto& act = sa_cassette.geminals_active(s, gen);
        for (size_t i = 0; i < act.size(); ++i)
            snap->layers[gen][act[i]] = { x[i], 2.0 - x[i] };
    };

    auto memo  = std::make_shared<WeightDerivsMemo>();
    auto sweep = [snap, &sa_cassette, s, gen, memo](const std::vector<double>& x) {
        if (!memo->valid || memo->x != x ||
            memo->wd.gen_ != REKSGradientEngine::sweep_generation()) {
            memo->wd    = REKSGradientEngine::accumulate_weight_derivs(s, gen, sa_cassette, *snap);
            memo->x     = x;
            memo->valid = true;
        }
        return memo->wd;
    };

    // energy, gradient, hessian evaluate the same sector-s partial SA energy;
    // extras are a FON-independent offset, omitted from gradient/hessian.
    obj.energy = [snap, &sa_cassette, &E_L, write, s, C_K_buf, C_L_buf]
                 (const std::vector<double>& x) {
        write(x);
        return geminal_target_energy(s, sa_cassette, *snap, E_L, *C_K_buf, *C_L_buf);
    };
    obj.gradient = [&sa_cassette, &E_L, write, sweep]
                   (const std::vector<double>& x) {
        write(x);
        return REKSGradientEngine::fon_gradient_from(sweep(x), sa_cassette, E_L);
    };
    obj.hessian = [snap, &sa_cassette, &E_L, write, sweep, s, gen]
                  (const std::vector<double>& x) {
        write(x);
        return REKSGradientEngine::fon_hessian_from(sweep(x), s, gen, sa_cassette, *snap, E_L);
    };
    return obj;
}

std::vector<double> fon_init_from_ci(int s,
                                     const studio::Cassette&      sa_cassette,
                                     const std::vector<double>& E_L,
                                     int iteration, bool narrate) {
    // GVB-PP CI Hamiltonian H (N x N, N = 2^U): CS microstate energies on the
    // diagonal, OS energy differences as off-diagonal couplings between CS
    // states that differ in one geminal. FON is read from the ground-state
    // eigenvector.
    const int U = static_cast<int>(sa_cassette.geminals_active(s, 0).size());  // active geminal count
    const int M = sa_cassette.n_active_orbitals();                              // active orbital count

    if (U == 0) return std::vector<double>(0, 1.0);

    const auto* n_geminals  = sa_cassette.geminal_templates();
    // Runtime extra determinants (index >= n_catalog) carry no FON; the loops
    // below run over the catalog GVB-PP structure only.
    const int n_catalog = sa_cassette.n_catalog_microstates();

    // in_sector[L] = 1 iff L is in sector s's catalog microstate set: union of
    // diag_deps(K).microstate_writes over sector s's configs.
    std::vector<char> in_sector(sa_cassette.n_total_microstates(), 0);
    for (int pos : sa_cassette.sector_config_positions(s)) {
        const int K = sa_cassette.K_indices[pos];
        for (int L : sa_cassette.diag_deps(K).microstate_writes) in_sector[L] = 1;
    }

    // POD-microstate helpers over active orbitals [0, M).
    auto occ = [&](int L, int o) {
        const auto& ms = sa_cassette.microstate(L);
        return static_cast<int>(ms.alpha[o]) + static_cast<int>(ms.beta[o]);
    };
    auto is_cs = [&](int L) {
        const auto& ms = sa_cassette.microstate(L);
        for (int i = 0; i < M; ++i)
            if (ms.alpha[i] != ms.beta[i]) return false;
        return true;
    };
    auto sz = [&](int L) {
        const auto& ms = sa_cassette.microstate(L);
        int na = 0, nb = 0;
        for (int i = 0; i < M; ++i) { na += ms.alpha[i]; nb += ms.beta[i]; }
        return 0.5 * (na - nb);
    };
    // Absolute orbital index j (0 or 1) of active n-geminal gi (pool position).
    auto gem_orb = [&](int gi, int j) {
        return n_geminals[sa_cassette.geminals_active(s, 0)[gi]].orbitals[j];
    };

    const int N = 1 << U;

    // GVB-CS microstates: each geminal has exactly one orbital doubly occupied.
    // cs_bits[i][gi] = which orbital (0/1) of geminal gi is doubly occupied in CS i.
    std::vector<int> cs_indices;
    std::vector<std::vector<int>> cs_bits;
    cs_indices.reserve(N);
    cs_bits.reserve(N);
    for (int L : sa_cassette.microstates_active) {
        if (L >= n_catalog) continue;  // skip runtime extra determinants
        if (!in_sector[L]) continue;   // sector s's determinants only
        bool gvb_cs = true;
        std::vector<int> bits(U);
        for (int gi = 0; gi < U; ++gi) {
            int n_doubly = 0;
            bits[gi] = 0;
            for (int j = 0; j < 2; ++j) {
                int o = gem_orb(gi, j);
                if (occ(L, o) == 2) { bits[gi] = j; n_doubly++; }
                else if (occ(L, o) != 0) { gvb_cs = false; break; }
            }
            if (!gvb_cs || n_doubly != 1) { gvb_cs = false; break; }
        }
        if (!gvb_cs) continue;
        cs_indices.push_back(L);
        cs_bits.push_back(std::move(bits));
    }

    if (static_cast<int>(cs_indices.size()) != N) {
        if (narrate && report::reports(4))
            outfile->Printf("  [CI-INIT] iter=%d Expected %d CS microstates (U=%d), found %zu; returning FON=1.0\n",
                            iteration, N, U, cs_indices.size());
        return std::vector<double>(U, 1.0);
    }

    // Per-pair coupling = E[Sz=0] - E[Sz=1] over OS microstates that single-occupy
    // orbitals (j1, j2) of geminal gi while spectator geminals stay bonding
    // (orbital 0 doubly occupied).
    auto matches_coupling_pattern = [&](int L, int gi, int j1, int j2) -> bool {
        const auto& ms = sa_cassette.microstate(L);
        int o_j1 = gem_orb(gi, j1), o_j2 = gem_orb(gi, j2);
        if (ms.alpha[o_j1] != 1) return false;
        if (ms.alpha[o_j2] + ms.beta[o_j2] != 1) return false;
        for (int m = 0; m < 2; ++m) {
            if (m == j1 || m == j2) continue;
            if (occ(L, gem_orb(gi, m)) != 0) return false;
        }
        for (int gj = 0; gj < U; ++gj) {
            if (gj == gi) continue;
            if (occ(L, gem_orb(gj, 0)) != 2) return false;
            for (int m = 1; m < 2; ++m)
                if (occ(L, gem_orb(gj, m)) != 0) return false;
        }
        return true;
    };

    // coupling_mat[gi][j1*2 + j2] symmetric, 2x2 per geminal.
    std::vector<std::vector<double>> coupling_mat(U);
    for (int gi = 0; gi < U; ++gi) {
        coupling_mat[gi].assign(2 * 2, 0.0);
        for (int j1 = 0; j1 < 2; ++j1) {
            for (int j2 = j1 + 1; j2 < 2; ++j2) {
                int Sz0_L = -1, Sz1_L = -1;
                for (int L : sa_cassette.microstates_active) {
                    if (L >= n_catalog) continue;
                    if (!in_sector[L]) continue;
                    if (is_cs(L)) continue;
                    if (!matches_coupling_pattern(L, gi, j1, j2)) continue;
                    double s = sz(L);
                    if (std::fabs(s) < 0.01) Sz0_L = L;
                    if (std::fabs(s - 1.0) < 0.01) Sz1_L = L;
                }
                if (Sz0_L >= 0 && Sz1_L >= 0) {
                    double c = E_L[Sz0_L] - E_L[Sz1_L];
                    coupling_mat[gi][j1 * 2 + j2] = c;
                    coupling_mat[gi][j2 * 2 + j1] = c;
                } else if (narrate && report::reports(4)) {
                    outfile->Printf("  [CI-INIT] iter=%d Missing coupling for geminal[%d] pair (%d,%d)\n",
                                    iteration, gi, j1, j2);
                }
            }
        }
    }

    // Weak-coupling check: max coupling per geminal vs threshold.
    const double coupling_threshold = 1.0e-8;
    std::vector<bool> fix_unit(U, false);
    bool all_weak = true;
    for (int gi = 0; gi < U; ++gi) {
        double max_coupling = 0.0;
        for (int j1 = 0; j1 < 2; ++j1)
            for (int j2 = j1 + 1; j2 < 2; ++j2)
                max_coupling = std::max(max_coupling, std::fabs(coupling_mat[gi][j1 * 2 + j2]));
        fix_unit[gi] = (max_coupling <= coupling_threshold);
        if (!fix_unit[gi]) all_weak = false;
    }
    if (all_weak) {
        if (narrate && report::reports(4))
            outfile->Printf("  [CI-INIT] iter=%d All couplings below threshold; defaulting to FON=1.0 (diradical)\n",
                            iteration);
        return std::vector<double>(U, 1.0);
    }

    // H stored column-major for C_DSYEV; both triangles filled explicitly.
    std::vector<double> H(N * N, 0.0);
    for (int i = 0; i < N; i++) {
        H[i + N * i] = E_L[cs_indices[i]];
        for (int j = i + 1; j < N; j++) {
            int diff_unit = -1, n_diff = 0;
            for (int gi = 0; gi < U; gi++) {
                if (cs_bits[i][gi] != cs_bits[j][gi]) { diff_unit = gi; n_diff++; }
            }
            if (n_diff == 1) {
                double c = coupling_mat[diff_unit][cs_bits[i][diff_unit] * 2 + cs_bits[j][diff_unit]];
                H[j + N * i] = c;
                H[i + N * j] = c;
            }
        }
    }

    std::vector<double> w(N);
    int lwork = std::max(3 * N - 1, 1);
    std::vector<double> work(lwork);
    int info = C_DSYEV('V', 'L', N, H.data(), N, w.data(), work.data(), lwork);
    if (info != 0) {
        if (narrate && report::reports(4))
            outfile->Printf("  [CI-INIT] iter=%d dsyev failed (info=%d), returning FON=1.0\n", iteration, info);
        return std::vector<double>(U, 1.0);
    }

    // Per-geminal FON = 2 * sum_{i: cs_bits[i][gi]==0} c_i^2 from the ground-state
    // eigenvector (column 0 of H after C_DSYEV).
    std::vector<double> fon(U, 0.0);
    for (int gi = 0; gi < U; gi++) {
        double sum_sq = 0.0;
        for (int i = 0; i < N; i++)
            if (cs_bits[i][gi] == 0) sum_sq += H[i] * H[i];
        fon[gi] = 2.0 * sum_sq;
    }

    // Weak-coupling override + interior clamping.
    for (int gi = 0; gi < U; gi++) {
        if (fix_unit[gi]) fon[gi] = 2.0 - constants::FON_BOUNDARY_EPS;
        fon[gi] = std::clamp(fon[gi], constants::FON_BOUNDARY_EPS,
                            2.0 - constants::FON_BOUNDARY_EPS);
    }

    if (narrate && report::reports(4)) {
        outfile->Printf("  [CI-INIT] iter=%d result: FON =", iteration);
        for (int gi = 0; gi < U; gi++) outfile->Printf(" %.10f", fon[gi]);
        outfile->Printf("\n");
    }

    return fon;
}

}  // namespace reks
}  // namespace psi
