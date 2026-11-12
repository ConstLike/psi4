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

#pragma once

#include "psi4/libmints/typedefs.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace psi {
class VBase;
class SuperFunctional;

namespace reks {
namespace xc {

/// Grid-direct spin-polarized XC of REKS microstates over a host RV potential.
class MicrostateXC {
   public:
    /// polarized_ = host->functional()->build_polarized().
    explicit MicrostateXC(std::shared_ptr<VBase> host);

    /// Per L without an AO V_xc^L: E_xc, tr(rho_s V_xc_s), active rows over [mo_lo, mo_lo + mo_width) and diagonal
    /// (need_diag) of C^T V_xc^L C; V_acc = sum_L C_L V_xc^L if given; grid_electrons() of the weighted density.
    /// Active MOs are the Ca_full columns [Ncore, Ncore + n_active).
    void compute(SharedMatrix Ca_full, int Ncore, int n_active,
                 const std::vector<std::vector<double>>& occ_alpha,
                 const std::vector<std::vector<double>>& occ_beta, const std::vector<double>& C_L,
                 std::vector<SharedMatrix>& arows_a, std::vector<SharedMatrix>& arows_b,
                 std::vector<std::vector<double>>& diag_a, std::vector<std::vector<double>>& diag_b,
                 std::vector<double>& tr_a, std::vector<double>& tr_b, std::vector<double>& E_xc,
                 SharedMatrix V_acc_a, SharedMatrix V_acc_b, bool need_diag, int mo_lo = 0);

    /// Quadrature integral of the C_L-weighted total density; empty until the first compute.
    std::optional<double> grid_electrons() const;

   private:
    /// Per-thread scratch, retained across calls and rebuilt only when Shape changes.
    /// Buffers are zeroed before use; retention alone does not clear them.
    struct Scratch {
        /// Every extent the buffers depend on. A mismatch rebuilds the whole set.
        struct Shape {
            int nbf = -1;
            int nmo = -1;
            int n_p = -1;
            int nact = -1;
            int n_micro = -1;
            int n_strings = -1;
            int max_points = -1;
            int max_functions = -1;
            int chunk = -1;
            int ansatz = -1;
            int num_threads = -1;
            int need_diag = -1;
            int acc_w = -1;             ///< second axis of the weighted active-row accumulator
            int project_in_block = -1;  ///< 1 when MO window == active block: TapWide is one slice wide
            bool operator==(const Shape& o) const {
                return nbf == o.nbf && nmo == o.nmo && n_p == o.n_p && nact == o.nact &&
                       n_micro == o.n_micro && n_strings == o.n_strings && max_points == o.max_points &&
                       max_functions == o.max_functions && chunk == o.chunk && ansatz == o.ansatz &&
                       num_threads == o.num_threads && need_diag == o.need_diag &&
                       acc_w == o.acc_w && project_in_block == o.project_in_block;
            }
        };
        /// One thread's buffers.
        struct Thread {
            SharedMatrix Cloc, chi, chix, chiy, chiz;
            SharedMatrix Sd, TapAcc_a, TapAcc_b, Vloc, MbpT;
            SharedMatrix SDrho, SDgx, SDgy, SDgz, SDgam, SDtau;
            SharedMatrix TapWide, TAwide, M1w, M2S, M3, Vcof;
            std::map<std::string, SharedVector> bin;
            std::shared_ptr<SuperFunctional> bworker;
            std::vector<double> diA, diB, trA, trB, exc, stauA, stauB;
            /// MO-diagonal staging: point-local products, and the per-slice diagonal.
            std::vector<double> Zd, Dd;
        };
        Shape shape;
        std::vector<Thread> threads;
        /// Call-level scratch shared by every thread: weighted-AO active rows per spin,
        /// written under one lock per microstate, and the MO-projection buffer.
        std::vector<double> wa_a_all, wa_b_all, proj;
    };

    std::shared_ptr<VBase> host_;                 ///< RV: grid, basis, point workers with the collocation cache
    std::shared_ptr<SuperFunctional> polarized_;  ///< host functional's spin-polarized twin
    Scratch scratch_;
    std::optional<double> grid_electrons_;
};

}  // namespace xc
}  // namespace reks
}  // namespace psi
