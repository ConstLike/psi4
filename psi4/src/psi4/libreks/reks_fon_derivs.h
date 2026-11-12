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

#ifndef REKS_FON_DERIVS_H
#define REKS_FON_DERIVS_H

/// FON gradient and Hessian of one (sector, generation) block from sum_K w_K dC_L^K/dx.

#include <cstdint>
#include <vector>

#include "reks_arch.h"
#include "reks_cassette.h"

namespace psi {
namespace reks {
namespace fon {

/// Counter bumped by every accumulate_weight_derivs call on this thread.
inline std::uint64_t& sweep_generation() {
    thread_local std::uint64_t g = 0;
    return g;
}

/// Thread-local rows dC[fi][L] = sum_K w_K dC_L^K/dx_fi, d2C[fi][L] = sum_K w_K d2C_L^K/dx_fi^2, valid
/// until the next sweep on this thread (gen_ = its sweep_generation()).
struct WeightDerivs {
    const double* dC  = nullptr;
    const double* d2C = nullptr;
    int n_fon  = 0;
    int stride = 0;
    std::uint64_t gen_ = 0;

    const double* dC_row(int fi)  const { return dC  + static_cast<size_t>(fi) * stride; }
    const double* d2C_row(int fi) const { return d2C + static_cast<size_t>(fi) * stride; }
};

/// dC and d2C over sector s's configs, one evaluation per active geminal of generation gen.
inline WeightDerivs accumulate_weight_derivs(int s, int gen,
                                             const studio::Cassette& sa_cassette,
                                             const studio::FONSnapshot& fons) {
    const auto& gem_active = sa_cassette.geminals_active(s, gen);
    const auto& positions  = sa_cassette.sector_config_positions(s);
    const int n_fon   = static_cast<int>(gem_active.size());
    const int stride  = sa_cassette.n_total_microstates();
    const int n_micro = sa_cassette.n_catalog_microstates();

    thread_local std::vector<double> dC_all, d2C_all, dC_K, d2C_K;
    const size_t need = static_cast<size_t>(n_fon) * stride;
    if (dC_all.size()  < need) dC_all.resize(need);
    if (d2C_all.size() < need) d2C_all.resize(need);
    if (static_cast<int>(dC_K.size())  < n_micro) dC_K.resize(n_micro);
    if (static_cast<int>(d2C_K.size()) < n_micro) d2C_K.resize(n_micro);

    for (int fi = 0; fi < n_fon; ++fi) {
        double* dC  = dC_all.data()  + static_cast<size_t>(fi) * stride;
        double* d2C = d2C_all.data() + static_cast<size_t>(fi) * stride;
        for (int L : sa_cassette.microstates_active) { dC[L] = 0.0; d2C[L] = 0.0; }
        const int g = gem_active[fi];
        for (int pos : positions) {
            const int    K   = sa_cassette.K_indices[pos];
            const double w_K = sa_cassette.K_weights[pos];
            sa_cassette.compute_weight_derivs(gen, K, g, fons, dC_K, d2C_K);
            for (int L : sa_cassette.diag_deps(K).microstate_writes) {
                dC[L]  += w_K * dC_K[L];
                d2C[L] += w_K * d2C_K[L];
            }
        }
    }
    return WeightDerivs{dC_all.data(), d2C_all.data(), n_fon, stride, ++sweep_generation()};
}

/// g_fon[i] = sum_L wd.dC_row(i)[L] * E_L[L]
inline std::vector<double> fon_gradient_from(const WeightDerivs& wd,
                                             const studio::Cassette& sa_cassette,
                                             const std::vector<double>& E_L) {
    std::vector<double> g_fon(wd.n_fon, 0.0);
    for (int fi = 0; fi < wd.n_fon; ++fi) {
        const double* dC_acc = wd.dC_row(fi);
        double acc = 0.0;
        for (int L : sa_cassette.microstates_active)
            acc += dC_acc[L] * E_L[L];
        g_fon[fi] = acc;
    }
    return g_fon;
}

inline std::vector<double> compute_fon_gradient(int s, int gen,
                                                const studio::Cassette& sa_cassette,
                                                const studio::FONSnapshot& fons,
                                                const std::vector<double>& E_L) {
    return fon_gradient_from(accumulate_weight_derivs(s, gen, sa_cassette, fons),
                             sa_cassette, E_L);
}

/// H_ff of block (s, gen): H_ii = sum_L d2C_i[L] E_L, H_ij = sum_K w_K sum_L d2C_L^K/dx_i dx_j E_L over sector
/// s's configs (row-major n_fon x n_fon).
inline std::vector<double> fon_hessian_from(const WeightDerivs& wd, int s, int gen,
                                            const studio::Cassette& sa_cassette,
                                            const studio::FONSnapshot& fons,
                                            const std::vector<double>& E_L) {
    const auto& gem_active = sa_cassette.geminals_active(s, gen);
    const auto& positions  = sa_cassette.sector_config_positions(s);
    const int n_fon   = wd.n_fon;
    const int n_micro = static_cast<int>(E_L.size());

    thread_local std::vector<double> d2C_mix;
    if (static_cast<int>(d2C_mix.size()) < n_micro) d2C_mix.resize(n_micro);

    std::vector<double> H_ff(n_fon * n_fon, 0.0);
    for (int fi = 0; fi < n_fon; ++fi) {
        const double* d2C_acc = wd.d2C_row(fi);
        double acc = 0.0;
        for (int L : sa_cassette.microstates_active)
            acc += d2C_acc[L] * E_L[L];
        H_ff[fi * n_fon + fi] = acc;
    }

    for (int fi = 0; fi < n_fon; ++fi) {
        const int g_i = gem_active[fi];
        for (int fj = fi + 1; fj < n_fon; ++fj) {
            const int g_j = gem_active[fj];
            double acc = 0.0;
            for (int pos : positions) {
                const int    K   = sa_cassette.K_indices[pos];
                const double w_K = sa_cassette.K_weights[pos];
                sa_cassette.compute_mixed_weight_derivs(gen, K, g_i, g_j, fons, d2C_mix);
                for (int L : sa_cassette.diag_deps(K).microstate_writes)
                    acc += w_K * d2C_mix[L] * E_L[L];
            }
            H_ff[fi * n_fon + fj] = acc;
            H_ff[fj * n_fon + fi] = acc;
        }
    }
    return H_ff;
}

inline std::vector<double> compute_fon_hessian(int s, int gen,
                                               const studio::Cassette& sa_cassette,
                                               const studio::FONSnapshot& fons,
                                               const std::vector<double>& E_L) {
    return fon_hessian_from(accumulate_weight_derivs(s, gen, sa_cassette, fons),
                            s, gen, sa_cassette, fons, E_L);
}

}  // namespace fon
}  // namespace reks
}  // namespace psi

#endif  // REKS_FON_DERIVS_H
