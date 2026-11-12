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

#include "reks_xc.h"

#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace reks {
namespace xc {

MicrostateXC::MicrostateXC(std::shared_ptr<VBase> host) : host_(std::move(host)) {
    polarized_ = host_->functional()->build_polarized();
}

std::optional<double> MicrostateXC::grid_electrons() const { return grid_electrons_; }

namespace {

// Project C_subset onto the block's local AO rows and form chi = phi @ Cloc (value, and for
// ansatz>=1 the three gradient tensors).
void build_microstate_chi(int npoints, int nlocal, int n_p, int ansatz, int coll_funcs,
                          double** C_subset_p, const std::vector<int>& function_map,
                          double** phi, double** phix, double** phiy, double** phiz,
                          double** Cloc, double** chi, double** chix, double** chiy, double** chiz) {
    for (int ml = 0; ml < nlocal; ++ml) {
        const int mg = function_map[ml];
        for (int p = 0; p < n_p; ++p) Cloc[ml][p] = C_subset_p[mg][p];
    }
    C_DGEMM('N', 'N', npoints, n_p, nlocal, 1.0, phi[0], coll_funcs, Cloc[0], n_p, 0.0, chi[0], n_p);
    if (ansatz >= 1) {
        C_DGEMM('N', 'N', npoints, n_p, nlocal, 1.0, phix[0], coll_funcs, Cloc[0], n_p, 0.0, chix[0], n_p);
        C_DGEMM('N', 'N', npoints, n_p, nlocal, 1.0, phiy[0], coll_funcs, Cloc[0], n_p, 0.0, chiy[0], n_p);
        C_DGEMM('N', 'N', npoints, n_p, nlocal, 1.0, phiz[0], coll_funcs, Cloc[0], n_p, 0.0, chiz[0], n_p);
    }
}

// Per-point sums of one channel over an MO column range: rho, grad rho, gamma_ss = |grad|^2 and the
// square-gradient sum (tau before its factor 0.5).
struct MicrostateChannelSums {
    double* rho = nullptr;
    double* gx = nullptr;
    double* gy = nullptr;
    double* gz = nullptr;
    double* gamma = nullptr;
    double* tau = nullptr;
};

// rho (and grad, gamma_ss, tau by ansatz) of one spin channel over columns [p0, n_p), continuing the seed
// sums over [0, p0) (null = 0); occ_scale 1 = UKS, 2 = RKS; tau_factor 0.5 = tau, 1 = raw seed.
void accumulate_microstate_rho_channel(int npoints, int p0, int n_p, int ansatz, double occ_scale,
                                       double tau_factor, const double* occ,
                                       double** chi, double** chix, double** chiy, double** chiz,
                                       const MicrostateChannelSums& seed,
                                       const MicrostateChannelSums& out) {
    for (int P = 0; P < npoints; ++P) {
        double r = seed.rho ? seed.rho[P] : 0.0;
        const double* chi_row = chi[P];
        for (int p = p0; p < n_p; ++p) {
            const double c2 = chi_row[p] * chi_row[p];
            r += occ_scale * occ[p] * c2;
        }
        out.rho[P] = r;
    }
    if (ansatz >= 1) {
        for (int P = 0; P < npoints; ++P) {
            double rx = seed.gx ? seed.gx[P] : 0.0;
            double ry = seed.gy ? seed.gy[P] : 0.0;
            double rz = seed.gz ? seed.gz[P] : 0.0;
            const double* c = chi[P];
            const double* cx = chix[P];
            const double* cy = chiy[P];
            const double* cz = chiz[P];
            for (int p = p0; p < n_p; ++p) {
                const double tw = 2.0 * occ_scale * occ[p] * c[p];
                rx += tw * cx[p];
                ry += tw * cy[p];
                rz += tw * cz[p];
            }
            out.gx[P] = rx;
            out.gy[P] = ry;
            out.gz[P] = rz;
            out.gamma[P] = rx * rx + ry * ry + rz * rz;
        }
    }
    if (ansatz >= 2) {
        for (int P = 0; P < npoints; ++P) {
            double t = seed.tau ? seed.tau[P] : 0.0;
            const double* cx = chix[P];
            const double* cy = chiy[P];
            const double* cz = chiz[P];
            for (int p = p0; p < n_p; ++p)
                t += occ_scale * occ[p] * (cx[p] * cx[p] + cy[p] * cy[p] + cz[p] * cz[p]);
            out.tau[P] = tau_factor * t;
        }
    }
}

}  // namespace

void MicrostateXC::compute(SharedMatrix Ca_full,
                                  int Ncore, int n_active,
                                  const std::vector<std::vector<double>>& occ_alpha,
                                  const std::vector<std::vector<double>>& occ_beta,
                                  const std::vector<double>& C_L,
                                  std::vector<SharedMatrix>& arows_a,
                                  std::vector<SharedMatrix>& arows_b,
                                  std::vector<std::vector<double>>& diag_a,
                                  std::vector<std::vector<double>>& diag_b,
                                  std::vector<double>& tr_a,
                                  std::vector<double>& tr_b,
                                  std::vector<double>& E_xc,
                                  SharedMatrix V_acc_a,
                                  SharedMatrix V_acc_b,
                                  bool need_diag,
                                  int mo_lo) {
    timer_on("REKS: microstate xc");

    const int n_micro = static_cast<int>(occ_alpha.size());
    if (n_micro == 0) { timer_off("REKS: microstate xc"); return; }

    const int nbf = host_->basis()->nbf();
    const auto grid = host_->grid();
    const auto workers = host_->properties();
    const int nthread = static_cast<int>(workers.size());
    const int nmo   = Ca_full->colspi()[0];
    const int n_p   = Ncore + n_active;        // occupied+active columns participating in rho
    const int nact  = n_active;
    const int ansatz = polarized_->ansatz();
    const int max_functions = grid->max_functions();
    const int max_points = grid->max_points();
    // V_acc_a/V_acc_b are shared nbf x nbf AO back-projection accumulators; a null pair
    // drops the whole E1 chain.
    const bool want_e1 = (V_acc_a != nullptr);
    // Projection window: the arows carry MO columns [mo_lo, mo_lo + mo_width).
    const int mo_width = arows_a[0]->colspi()[0];
    // project_in_block: the MO window equals the active block, so the weighted-AO rows
    // accumulate directly in the nact-wide active basis instead of the nbf-wide AO buffer.
    const bool project_in_block = (mo_lo == Ncore && mo_width == nact);
    const int  acc_w = project_in_block ? nact : nbf;

    // Distinct spin strings: each density is accumulated once per block and shared by its channels.
    std::vector<std::vector<double>> string_occ;
    std::vector<int> str_a(n_micro), str_b(n_micro);
    {
        std::map<std::vector<double>, int> string_id;
        auto intern = [&](const std::vector<double>& o) {
            auto it = string_id.find(o);
            if (it != string_id.end()) return it->second;
            const int id = static_cast<int>(string_occ.size());
            string_id.emplace(o, id);
            string_occ.push_back(o);
            return id;
        };
        for (int L = 0; L < n_micro; ++L) {
            str_a[L] = intern(occ_alpha[L]);
            str_b[L] = intern(occ_beta[L]);
        }
    }
    const int n_strings = static_cast<int>(string_occ.size());

    // Unit core occupations: the core partial sum is accumulated once per block and seeded into each string.
    const std::vector<double> core_occ(Ncore, 1.0);

    // Closed-shell microstate (occ_alpha==occ_beta): beta XC channel equals alpha. Compute alpha
    // only and mirror beta; skips the redundant beta Tap build, arows GEMMs, and diagonal.
    std::vector<char> closed_ms(n_micro);
    for (int L = 0; L < n_micro; ++L) closed_ms[L] = (str_a[L] == str_b[L]) ? 1 : 0;

    // Microstate chunk: one batched functional call and one batched active-row GEMM set per
    // chunk (up to 2 channel slices per microstate). Cap keeps per-thread Tap scratch bounded.
    int chunk = static_cast<int>((8ull << 20) / ((size_t)max_points * max_functions * 8 * 2));
    chunk = std::max(4, std::min(chunk, 32));
    chunk = std::min(chunk, n_micro);
    const int max_slices = 2 * chunk;
    // Slices staged in TapWide/TAwide per GEMM, within a 2 MiB budget for TapWide (AO row per point per slice).
    int stage = static_cast<int>((2ull << 20) / ((size_t)max_points * max_functions * 8));
    stage = std::max(2, std::min(stage, max_slices));

    for (int L = 0; L < n_micro; ++L) {
        diag_a[L].assign(need_diag ? nmo : 0, 0.0);
        diag_b[L].assign(need_diag ? nmo : 0, 0.0);
        tr_a[L] = tr_b[L] = E_xc[L] = 0.0;
    }
    if (want_e1) { V_acc_a->zero(); V_acc_b->zero(); }

    // Weighted-AO active rows (already the active MO block when project_in_block), reduced
    // over threads and moved into arows once after the sweep.
    scratch_.wa_a_all.assign((size_t)n_micro * nact * acc_w, 0.0);
    scratch_.wa_b_all.assign((size_t)n_micro * nact * acc_w, 0.0);
    std::vector<double>& wa_a_all = scratch_.wa_a_all;
    std::vector<double>& wa_b_all = scratch_.wa_b_all;

    // MO columns the chi tensor carries: the full MO range only when the MO diagonal is wanted,
    // otherwise the core+active window [0, n_p) that feeds rho and the active rows.
    const int nchi = need_diag ? nmo : n_p;
    const int bpts = chunk * max_points;
    // Tap coefficients per point: the phi multiplier and the three gradient multipliers.
    const int ncof = 4;
    // Points staged at a time when forming the MO diagonal.
    constexpr int kZChunk = 64;

    Scratch::Shape shape;
    shape.nbf = nbf;
    shape.nmo = nmo;
    shape.n_p = n_p;
    shape.nact = nact;
    shape.n_micro = n_micro;
    shape.n_strings = n_strings;
    shape.max_points = max_points;
    shape.max_functions = max_functions;
    shape.chunk = chunk;
    shape.ansatz = ansatz;
    shape.num_threads = nthread;
    shape.need_diag = need_diag ? 1 : 0;
    shape.acc_w = acc_w;
    shape.project_in_block = project_in_block ? 1 : 0;

    // Thread scratch, rebuilt only when the shape changes. The per-call zeroing of the
    // accumulators happens inside the parallel region, each thread on its own buffers.
    Scratch& sc = scratch_;
    if (!(sc.shape == shape)) {
        sc.threads.clear();
        sc.threads.resize(nthread);
        sc.shape = shape;
        for (int t = 0; t < nthread; ++t) {
            auto& ts = sc.threads[t];
            ts.Cloc = std::make_shared<Matrix>("Cloc", max_functions, nchi);
            ts.chi = std::make_shared<Matrix>("chi", max_points, nchi);
            ts.TapAcc_a = std::make_shared<Matrix>("TapAcc_a", max_points, max_functions);
            ts.TapAcc_b = std::make_shared<Matrix>("TapAcc_b", max_points, max_functions);
            ts.Vloc = std::make_shared<Matrix>("Vloc", max_functions, max_functions);  // E1 back-proj
            if (ansatz >= 2) ts.MbpT = std::make_shared<Matrix>("MbpT", max_points, max_functions);
            if (ansatz >= 1) {
                ts.chix = std::make_shared<Matrix>("chix", max_points, nchi);
                ts.chiy = std::make_shared<Matrix>("chiy", max_points, nchi);
                ts.chiz = std::make_shared<Matrix>("chiz", max_points, nchi);
            }
            if (ansatz >= 2) {
                ts.Sd = std::make_shared<Matrix>("Sd", max_points, nact);
                // Holds either an AO-local row (nlocal) or the projected active block (nact).
                ts.M3 = std::make_shared<Matrix>("M3", n_active, std::max(max_functions, nact));
            }
            // Distinct-string channel densities (row = string): rho, grad, gamma_ss, tau. Row
            // n_strings holds the core partial sums the strings are seeded from.
            const int nsd = n_strings + 1;
            ts.SDrho = std::make_shared<Matrix>("SDrho", nsd, max_points);
            if (ansatz >= 1) {
                ts.SDgx = std::make_shared<Matrix>("SDgx", nsd, max_points);
                ts.SDgy = std::make_shared<Matrix>("SDgy", nsd, max_points);
                ts.SDgz = std::make_shared<Matrix>("SDgz", nsd, max_points);
                ts.SDgam = std::make_shared<Matrix>("SDgam", nsd, max_points);
            }
            if (ansatz >= 2) ts.SDtau = std::make_shared<Matrix>("SDtau", nsd, max_points);
            // Staged slices: TAwide in active MO (offset j*nact), TapWide in local AO (j*nlocal), M1w/M2S hold
            // wa_s(L)[i][nu] = sum_P chi_act[P][i] TapWide[P][nu] + TAwide[P][i] phi[P][nu]; Cloc = local rows of Ca.
            ts.TapWide = std::make_shared<Matrix>("TapWide", max_points,
                                                  max_functions * (project_in_block ? 1 : stage));
            // Vcof row per slice: Tap coefficients (a0, ax, ay, az) per point; 2*n_micro rows for the MO
            // diagonal, max_slices otherwise.
            ts.Vcof = std::make_shared<Matrix>("Vcof", need_diag ? 2 * n_micro : max_slices,
                                               ncof * max_points);
            ts.TAwide = std::make_shared<Matrix>("TAwide", max_points, n_active * stage);
            // Per slice, M1w holds either an AO-local row (nlocal) or the projected active
            // block (nact).
            ts.M1w = std::make_shared<Matrix>("M1w", n_active,
                                              std::max(max_functions, nact) * stage);
            ts.M2S = std::make_shared<Matrix>("M2S", n_active * stage, max_functions);
            // Batched functional evaluation: one compute_functional call per chunk over stacked
            // per-microstate point slices (slot j at offset j*npoints).
            auto& bin = ts.bin;
            bin["RHO_A"] = std::make_shared<Vector>("RHO_A", bpts);
            bin["RHO_B"] = std::make_shared<Vector>("RHO_B", bpts);
            if (ansatz >= 1) {
                bin["GAMMA_AA"] = std::make_shared<Vector>("GAMMA_AA", bpts);
                bin["GAMMA_AB"] = std::make_shared<Vector>("GAMMA_AB", bpts);
                bin["GAMMA_BB"] = std::make_shared<Vector>("GAMMA_BB", bpts);
            }
            if (ansatz >= 2) {
                bin["TAU_A"] = std::make_shared<Vector>("TAU_A", bpts);
                bin["TAU_B"] = std::make_shared<Vector>("TAU_B", bpts);
            }
            ts.bworker = polarized_->build_worker();
            ts.bworker->set_max_points(bpts);
            ts.bworker->allocate();
            // Per-L MO diagonal (flat, n_micro*nmo), unallocated when no diagonal is wanted.
            if (need_diag) {
                ts.diA.resize((size_t)n_micro * nmo);
                ts.diB.resize((size_t)n_micro * nmo);
            }
            ts.trA.resize(n_micro);
            ts.trB.resize(n_micro);
            ts.exc.resize(n_micro);
            if (need_diag) {
                ts.Zd.assign((size_t)ncof * kZChunk * nmo, 0.0);
                ts.Dd.resize((size_t)2 * n_micro * nmo);
            }
            if (ansatz >= 2) { ts.stauA.resize(max_points); ts.stauB.resize(max_points); }  // meta E1 per-point
        }
    }
    std::vector<double> nqT(nthread, 0.0);  // C_L-weighted quadrature electron count

    double** Ca_p = Ca_full->pointer();
    double** V_acc_a_p = want_e1 ? V_acc_a->pointer() : nullptr;
    double** V_acc_b_p = want_e1 ? V_acc_b->pointer() : nullptr;

    // One lock per microstate for the shared weighted-AO rows; std::mutex, since omp_lock_t's size
    // depends on which omp.h the include path finds.
    std::vector<std::mutex> wa_lock(n_micro);
    auto wa_set = [&wa_lock](int L) { wa_lock[L].lock(); };
    auto wa_unset = [&wa_lock](int L) { wa_lock[L].unlock(); };

    const char* t_grid = need_diag ? "REKS: xc_grid[full]" : "REKS: xc_grid[probe]";
    const char* t_proj = need_diag ? "REKS: xc_proj[full]" : "REKS: xc_proj[probe]";

    timer_on(t_grid);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto pworker = workers[rank];
        auto& ts = sc.threads[rank];

        // Per-call reset of this thread's accumulators. Every other buffer is written in full
        // before it is read, either by a beta=0 GEMM or by an explicit per-block fill.
        std::fill(ts.diA.begin(), ts.diA.end(), 0.0);
        std::fill(ts.diB.begin(), ts.diB.end(), 0.0);
        std::fill(ts.trA.begin(), ts.trA.end(), 0.0);
        std::fill(ts.trB.begin(), ts.trB.end(), 0.0);
        std::fill(ts.exc.begin(), ts.exc.end(), 0.0);

        auto bworker = ts.bworker;
        double** Cloc = ts.Cloc->pointer();
        double** chi  = ts.chi->pointer();
        double** chix = (ansatz >= 1) ? ts.chix->pointer() : nullptr;
        double** chiy = (ansatz >= 1) ? ts.chiy->pointer() : nullptr;
        double** chiz = (ansatz >= 1) ? ts.chiz->pointer() : nullptr;
        double** Sd   = (ansatz >= 2) ? ts.Sd->pointer() : nullptr;
        double** TapAcc_a = ts.TapAcc_a->pointer();
        double** TapAcc_b = ts.TapAcc_b->pointer();
        double** Vloc_p   = ts.Vloc->pointer();
        double** MbpT     = (ansatz >= 2) ? ts.MbpT->pointer() : nullptr;
        double** SDrho = ts.SDrho->pointer();
        double** SDgx  = (ansatz >= 1) ? ts.SDgx->pointer() : nullptr;
        double** SDgy  = (ansatz >= 1) ? ts.SDgy->pointer() : nullptr;
        double** SDgz  = (ansatz >= 1) ? ts.SDgz->pointer() : nullptr;
        double** SDgam = (ansatz >= 1) ? ts.SDgam->pointer() : nullptr;
        double** SDtau = (ansatz >= 2) ? ts.SDtau->pointer() : nullptr;
        double** TapWide = ts.TapWide->pointer();
        double** Vcof = ts.Vcof->pointer();
        double** TAwide  = ts.TAwide->pointer();
        double** M1w  = ts.M1w->pointer();
        double** M2S  = ts.M2S->pointer();
        double** M3   = (ansatz >= 2) ? ts.M3->pointer() : nullptr;
        // Tap/TapAcc below feed only the E1 back-projection and, outside project_in_block,
        // the active-row GEMM; skip that per-point AO pass when neither is wanted.
        const bool need_ao_pass = want_e1 || !project_in_block;
        const int ld_tw  = max_functions * (project_in_block ? 1 : stage);
        const int ld_taw = nact * stage;
        const int ld_m1w = std::max(max_functions, nact) * stage;
        const int ld_m3  = std::max(max_functions, nact);
        double* waA = wa_a_all.data();
        double* waB = wa_b_all.data();
        double* diA = need_diag ? ts.diA.data() : nullptr;
        double* diB = need_diag ? ts.diB.data() : nullptr;
        // MO diagonal: Zd = point-local products of kZChunk points, Dd = per-slice diagonal (rows 1..3
        // zero for LSDA).
        std::vector<double>& Zd = ts.Zd;
        std::vector<double>& Dd = ts.Dd;
        double* trAr = ts.trA.data();
        double* trBr = ts.trB.data();
        double* excr = ts.exc.data();
        double* stauA = (ansatz >= 2) ? ts.stauA.data() : nullptr;
        double* stauB = (ansatz >= 2) ? ts.stauB.data() : nullptr;
        auto& bin = ts.bin;
        double* brho_a = bin.at("RHO_A")->pointer();
        double* brho_b = bin.at("RHO_B")->pointer();
        double* bgaa = (ansatz >= 1) ? bin.at("GAMMA_AA")->pointer() : nullptr;
        double* bgab = (ansatz >= 1) ? bin.at("GAMMA_AB")->pointer() : nullptr;
        double* bgbb = (ansatz >= 1) ? bin.at("GAMMA_BB")->pointer() : nullptr;
        double* btau_a = (ansatz >= 2) ? bin.at("TAU_A")->pointer() : nullptr;
        double* btau_b = (ansatz >= 2) ? bin.at("TAU_B")->pointer() : nullptr;
        std::vector<int> slice_L(2 * n_micro), slice_beta(2 * n_micro);
        // Channel rows of the per-string density tables; row n_strings carries the core prefix.
        auto channel_rows = [&](int row) {
            MicrostateChannelSums c;
            c.rho = SDrho[row];
            if (ansatz >= 1) { c.gx = SDgx[row]; c.gy = SDgy[row]; c.gz = SDgz[row]; c.gamma = SDgam[row]; }
            if (ansatz >= 2) c.tau = SDtau[row];
            return c;
        };
        const MicrostateChannelSums core_sums = channel_rows(n_strings);

#pragma omp for schedule(guided)
        for (size_t Q = 0; Q < grid->blocks().size(); ++Q) {
            auto block = grid->blocks()[Q];
            const int npoints = block->npoints();
            const double* w = block->w();
            const auto& function_map = block->functions_local_to_global();
            const int nlocal = static_cast<int>(function_map.size());

            pworker->prepare_basis_only(block, false);
            auto phi_mat = pworker->basis_value("PHI");
            const int coll_funcs = phi_mat->ncol();
            double** phi = phi_mat->pointer();
            double** phix = nullptr; double** phiy = nullptr; double** phiz = nullptr;
            if (ansatz >= 1) {
                phix = pworker->basis_value("PHI_X")->pointer();
                phiy = pworker->basis_value("PHI_Y")->pointer();
                phiz = pworker->basis_value("PHI_Z")->pointer();
            }
            build_microstate_chi(npoints, nlocal, nchi, ansatz, coll_funcs, Ca_p, function_map,
                                 phi, phix, phiy, phiz, Cloc, chi, chix, chiy, chiz);

            // Per-block E1 weighted-AO accumulator (reset).
            if (need_ao_pass) for (int P = 0; P < npoints; ++P) {
                std::fill(TapAcc_a[P], TapAcc_a[P] + nlocal, 0.0);
                std::fill(TapAcc_b[P], TapAcc_b[P] + nlocal, 0.0);
            }
            if (want_e1 && ansatz >= 2) {
                std::fill(stauA, stauA + npoints, 0.0);
                std::fill(stauB, stauB + npoints, 0.0);
            }

            // Core sums over [0, Ncore) once per block, then each distinct string over [Ncore, n_p) from them.
            accumulate_microstate_rho_channel(npoints, 0, Ncore, ansatz, 1.0, 1.0, core_occ.data(),
                                              chi, chix, chiy, chiz,
                                              MicrostateChannelSums(), core_sums);
            for (int s = 0; s < n_strings; ++s) {
                accumulate_microstate_rho_channel(npoints, Ncore, n_p, ansatz, 1.0, 0.5,
                                                  string_occ[s].data(), chi, chix, chiy, chiz,
                                                  core_sums, channel_rows(s));
            }

            int sbase = 0;  // slices staged so far in this block
            for (int c0 = 0; c0 < n_micro; c0 += chunk) {
                const int nc = std::min(chunk, n_micro - c0);
                // vbase: Vcof rows are addressed block-wide (offset sbase) only for the MO
                // diagonal; otherwise rows are chunk-local (offset 0).
                const int vbase = need_diag ? sbase : 0;

                // Pack stacked functional inputs (slot j = microstate c0+j).
                for (int j = 0; j < nc; ++j) {
                    const int L = c0 + j;
                    const size_t off = (size_t)j * npoints;
                    const double* ra = SDrho[str_a[L]];
                    const double* rb = SDrho[str_b[L]];
                    std::copy(ra, ra + npoints, brho_a + off);
                    std::copy(rb, rb + npoints, brho_b + off);
                    if (ansatz >= 1) {
                        std::copy(SDgam[str_a[L]], SDgam[str_a[L]] + npoints, bgaa + off);
                        std::copy(SDgam[str_b[L]], SDgam[str_b[L]] + npoints, bgbb + off);
                        const double* gax = SDgx[str_a[L]];
                        const double* gay = SDgy[str_a[L]];
                        const double* gaz = SDgz[str_a[L]];
                        const double* gbx = SDgx[str_b[L]];
                        const double* gby = SDgy[str_b[L]];
                        const double* gbz = SDgz[str_b[L]];
                        double* gab = bgab + off;
                        for (int P = 0; P < npoints; ++P)
                            gab[P] = gax[P]*gbx[P] + gay[P]*gby[P] + gaz[P]*gbz[P];
                    }
                    if (ansatz >= 2) {
                        std::copy(SDtau[str_a[L]], SDtau[str_a[L]] + npoints, btau_a + off);
                        std::copy(SDtau[str_b[L]], SDtau[str_b[L]] + npoints, btau_b + off);
                    }
                }

                // One functional evaluation for the whole chunk: the functional is pointwise,
                // stacked slices give the per-microstate values.
                auto& bvals = bworker->compute_functional(bin, nc * npoints);
                double* zk = bvals["V"]->pointer();
                double* v_rho_a = bvals["V_RHO_A"]->pointer();
                double* v_rho_b = bvals["V_RHO_B"]->pointer();
                double* v_gamma_aa=nullptr;double* v_gamma_ab=nullptr;double* v_gamma_bb=nullptr;
                double* v_tau_a=nullptr;double* v_tau_b=nullptr;
                if (ansatz >= 1) {
                    v_gamma_aa=bvals["V_GAMMA_AA"]->pointer();
                    v_gamma_ab=bvals["V_GAMMA_AB"]->pointer();
                    v_gamma_bb=bvals["V_GAMMA_BB"]->pointer();
                }
                if (ansatz >= 2) { v_tau_a=bvals["V_TAU_A"]->pointer(); v_tau_b=bvals["V_TAU_B"]->pointer(); }

                int ns = 0;  // channel slices opened this chunk; indexes Vcof and slice_L
                int nw = 0;  // slices currently held in TapWide/TAwide
                int w0 = 0;  // chunk-relative index of the first slice held

                // wa_s(L)[i][function_map[nu]] += (chi_act^T TapWide + TAwide^T phi)[i][nu]; in MO the active
                // block is A + A^T, A(i,j) = sum_P chi_act[P][i] TAwide[P][sj,j], with no AO row formed.
                auto flush_stage = [&]() {
                    if (nw == 0) return;
                    // Point-major pass over the slices held, q = Ncore+i:
                    //   TAwide[P][sj,i] = a0 chi[P][q] + a.grad chi[P][q]
                    for (int P = 0; P < npoints; ++P) {
                        const double* ch = chi[P];
                        const double* cx = (ansatz >= 1) ? chix[P] : nullptr;
                        const double* cy = (ansatz >= 1) ? chiy[P] : nullptr;
                        const double* cz = (ansatz >= 1) ? chiz[P] : nullptr;
                        double* taw = TAwide[P];
                        for (int sj = 0; sj < nw; ++sj) {
                            const double* cf = Vcof[vbase + w0 + sj] + (size_t)ncof*P;
                            double* t = taw + (size_t)sj*nact;
                            const double a0 = cf[0];
                            if (ansatz >= 1) {
                                const double ax = cf[1], ay = cf[2], az = cf[3];
                                for (int i = 0; i < nact; ++i) {
                                    const int q = Ncore + i;
                                    t[i] = a0*ch[q] + ax*cx[q] + ay*cy[q] + az*cz[q];
                                }
                            } else {
                                for (int i = 0; i < nact; ++i) t[i] = a0*ch[Ncore + i];
                            }
                        }
                    }
                    if (nlocal > 0) {
                        if (project_in_block) {
                            C_DGEMM('T','N', nact, nw*nact, npoints, 1.0, &chi[0][Ncore], nchi,
                                    TAwide[0], ld_taw, 0.0, M1w[0], ld_m1w);
                            for (int sj = 0; sj < nw; ++sj) {
                                const int Ls = slice_L[sbase + w0 + sj];
                                double* tgt = (slice_beta[sbase + w0 + sj] ? waB : waA)
                                            + (size_t)Ls*nact*nact;
                                wa_set(Ls);
                                for (int i = 0; i < nact; ++i) {
                                    const double* Ai = M1w[i] + (size_t)sj*nact;
                                    double* trow = tgt + (size_t)i*nact;
                                    for (int j = 0; j < nact; ++j)
                                        trow[j] += Ai[j] + M1w[j][(size_t)sj*nact + i];
                                }
                                wa_unset(Ls);
                            }
                        } else {
                            C_DGEMM('T','N', nact, nw*nlocal, npoints, 1.0, &chi[0][Ncore], nchi,
                                    TapWide[0], ld_tw, 0.0, M1w[0], ld_m1w);
                            C_DGEMM('T','N', nw*nact, nlocal, npoints, 1.0, TAwide[0], ld_taw,
                                    phi[0], coll_funcs, 0.0, M2S[0], max_functions);
                            for (int sj = 0; sj < nw; ++sj) {
                                const int Ls = slice_L[sbase + w0 + sj];
                                double* tgt = (slice_beta[sbase + w0 + sj] ? waB : waA)
                                            + (size_t)Ls*nact*nbf;
                                wa_set(Ls);
                                for (int i = 0; i < nact; ++i) {
                                    const double* r1 = M1w[i] + (size_t)sj*nlocal;
                                    const double* r2 = M2S[(size_t)sj*nact + i];
                                    double* trow = tgt + (size_t)i*nbf;
                                    for (int k = 0; k < nlocal; ++k) trow[function_map[k]] += r1[k] + r2[k];
                                }
                                wa_unset(Ls);
                            }
                        }
                    }
                    w0 += nw;
                    nw = 0;
                };

                for (int j = 0; j < nc; ++j) {
                    const int L = c0 + j;
                    const double cl = C_L[L];
                    const bool closed = closed_ms[L];
                    const size_t off = (size_t)j * npoints;
                    const double* ra = SDrho[str_a[L]];
                    const double* rb = SDrho[str_b[L]];
                    const double* vra = v_rho_a + off;
                    const double* vrb = v_rho_b + off;
                    const double* vgaa = (ansatz >= 1) ? v_gamma_aa + off : nullptr;
                    const double* vgab = (ansatz >= 1) ? v_gamma_ab + off : nullptr;
                    const double* vgbb = (ansatz >= 1) ? v_gamma_bb + off : nullptr;
                    const double* gab = (ansatz >= 1) ? bgab + off : nullptr;

                    excr[L] += C_DDOT(npoints, w, 1, zk + off, 1);

                    // tr(D^s V_xc^s) = sum_P w [v_rho rho + 2 v_gamma_ss |grad|^2 + v_gamma_os grad.grad + 2 v_tau tau]
                    // (tau = 0.5 sum_d D phid phid).
                    double ta = 0.0, tb = 0.0;
                    for (int P = 0; P < npoints; ++P) { ta += w[P]*vra[P]*ra[P]; tb += w[P]*vrb[P]*rb[P]; }
                    if (ansatz >= 1) {
                        const double* gaa = SDgam[str_a[L]];
                        const double* gbb = SDgam[str_b[L]];
                        for (int P = 0; P < npoints; ++P) {
                            ta += w[P]*(2.0*vgaa[P]*gaa[P] + vgab[P]*gab[P]);
                            tb += w[P]*(2.0*vgbb[P]*gbb[P] + vgab[P]*gab[P]);
                        }
                    }
                    if (ansatz >= 2) {
                        const double* tua = SDtau[str_a[L]];
                        const double* tub = SDtau[str_b[L]];
                        const double* vta = v_tau_a + off;
                        const double* vtb = v_tau_b + off;
                        for (int P = 0; P < npoints; ++P) {
                            ta += 2.0*w[P]*vta[P]*tua[P];
                            tb += 2.0*w[P]*vtb[P]*tub[P];
                        }
                    }
                    trAr[L] += ta; trBr[L] += tb;

                    double nq = 0.0;
                    for (int P = 0; P < npoints; ++P) nq += w[P]*(ra[P] + rb[P]);
                    nqT[rank] += cl * nq;

                    // Tap[P][k] = 0.5 v_rho w phi_k + w (2 v_gamma_ss grad rho_s + v_gamma_os grad rho_os) . grad phi_k
                    // into TapWide and Vcof, TapAcc += C_L Tap (beta mirrored for a closed microstate).
                    const double* gax = (ansatz >= 1) ? SDgx[str_a[L]] : nullptr;
                    const double* gay = (ansatz >= 1) ? SDgy[str_a[L]] : nullptr;
                    const double* gaz = (ansatz >= 1) ? SDgz[str_a[L]] : nullptr;
                    const double* gbx = (ansatz >= 1) ? SDgx[str_b[L]] : nullptr;
                    const double* gby = (ansatz >= 1) ? SDgy[str_b[L]] : nullptr;
                    const double* gbz = (ansatz >= 1) ? SDgz[str_b[L]] : nullptr;

                    {  // alpha slice
                        if (nw == stage) flush_stage();
                        const size_t scol = project_in_block ? 0 : (size_t)nw * nlocal;
                        double* cof = Vcof[vbase + ns];
                        slice_L[sbase + ns] = L; slice_beta[sbase + ns] = 0; ++ns; ++nw;
                        for (int P = 0; P < npoints; ++P) {
                            const double a0 = 0.5 * vra[P] * w[P];
                            double ax = 0.0, ay = 0.0, az = 0.0;
                            if (ansatz >= 1) {
                                ax = w[P]*(2.0*vgaa[P]*gax[P] + vgab[P]*gbx[P]);
                                ay = w[P]*(2.0*vgaa[P]*gay[P] + vgab[P]*gby[P]);
                                az = w[P]*(2.0*vgaa[P]*gaz[P] + vgab[P]*gbz[P]);
                            }
                            double* cf = cof + (size_t)ncof*P;
                            cf[0] = a0; cf[1] = ax; cf[2] = ay; cf[3] = az;
                            if (!need_ao_pass) continue;
                            double* tap = TapWide[P] + scol;
                            double* acc = TapAcc_a[P];
                            double* accb = TapAcc_b[P];
                            const double* ph = phi[P];
                            if (ansatz >= 1) {
                                const double* px = phix[P]; const double* py = phiy[P]; const double* pz = phiz[P];
                                if (closed) {
                                    for (int k = 0; k < nlocal; ++k) {
                                        const double tv = a0*ph[k] + ax*px[k] + ay*py[k] + az*pz[k];
                                        tap[k] = tv;
                                        acc[k] += cl*tv;
                                        accb[k] += cl*tv;
                                    }
                                } else {
                                    for (int k = 0; k < nlocal; ++k) {
                                        const double tv = a0*ph[k] + ax*px[k] + ay*py[k] + az*pz[k];
                                        tap[k] = tv;
                                        acc[k] += cl*tv;
                                    }
                                }
                            } else if (closed) {
                                for (int k = 0; k < nlocal; ++k) {
                                    const double tv = a0*ph[k];
                                    tap[k] = tv;
                                    acc[k] += cl*tv;
                                    accb[k] += cl*tv;
                                }
                            } else {
                                for (int k = 0; k < nlocal; ++k) {
                                    const double tv = a0*ph[k];
                                    tap[k] = tv;
                                    acc[k] += cl*tv;
                                }
                            }
                        }
                    }
                    if (!closed) {  // beta slice
                        if (nw == stage) flush_stage();
                        const size_t scol = project_in_block ? 0 : (size_t)nw * nlocal;
                        double* cof = Vcof[vbase + ns];
                        slice_L[sbase + ns] = L; slice_beta[sbase + ns] = 1; ++ns; ++nw;
                        for (int P = 0; P < npoints; ++P) {
                            const double b0 = 0.5 * vrb[P] * w[P];
                            double bx = 0.0, by = 0.0, bz = 0.0;
                            if (ansatz >= 1) {
                                bx = w[P]*(2.0*vgbb[P]*gbx[P] + vgab[P]*gax[P]);
                                by = w[P]*(2.0*vgbb[P]*gby[P] + vgab[P]*gay[P]);
                                bz = w[P]*(2.0*vgbb[P]*gbz[P] + vgab[P]*gaz[P]);
                            }
                            double* cf = cof + (size_t)ncof*P;
                            cf[0] = b0; cf[1] = bx; cf[2] = by; cf[3] = bz;
                            if (!need_ao_pass) continue;
                            double* tap = TapWide[P] + scol;
                            double* acc = TapAcc_b[P];
                            const double* ph = phi[P];
                            if (ansatz >= 1) {
                                const double* px = phix[P]; const double* py = phiy[P]; const double* pz = phiz[P];
                                for (int k = 0; k < nlocal; ++k) {
                                    const double tv = b0*ph[k] + bx*px[k] + by*py[k] + bz*pz[k];
                                    tap[k] = tv;
                                    acc[k] += cl*tv;
                                }
                            } else {
                                for (int k = 0; k < nlocal; ++k) {
                                    const double tv = b0*ph[k];
                                    tap[k] = tv;
                                    acc[k] += cl*tv;
                                }
                            }
                        }
                    }
                }

                // meta: sum_d (v_tau w) chi_d_i chi_d_q via Sd[P][i] = v_tau w chi_d[P][Ncore+i] against phi_d;
                // the diagonal adds (v_tau w) chi_d[P][q]^2 directly.
                if (ansatz >= 2) {
                    for (int j = 0; j < nc; ++j) {
                        const int L = c0 + j;
                        const double cl = C_L[L];
                        const bool closed = closed_ms[L];
                        const size_t off = (size_t)j * npoints;
                        const double* vta = v_tau_a + off;
                        const double* vtb = v_tau_b + off;
                        double* waAL = waA + (size_t)L*nact*acc_w;
                        double* waBL = waB + (size_t)L*nact*acc_w;
                        double** chid_arr[3] = {chix, chiy, chiz};
                        double** phid_arr[3] = {phix, phiy, phiz};
                        double* diAL = need_diag ? (diA + (size_t)L*nmo) : nullptr;
                        // The three directions accumulate into M3 before the single locked
                        // scatter below.
                        const int m3w = project_in_block ? nact : nlocal;
                        auto meta_dir = [&](const double* vt, double* waL, double* diL) {
                            for (int i = 0; i < nact; ++i) std::fill(M3[i], M3[i] + m3w, 0.0);
                            for (int d = 0; d < 3; ++d) {
                                double** chid = chid_arr[d];
                                double** phid = phid_arr[d];
                                for (int P = 0; P < npoints; ++P) {
                                    const double s = vt[P]*w[P];
                                    const double* cd = chid[P];
                                    double* sd = Sd[P];
                                    for (int i = 0; i < nact; ++i) sd[i] = s*cd[Ncore + i];
                                    if (diL) {
                                        for (int q = 0; q < nmo; ++q) diL[q] += cd[q]*(s*cd[q]);
                                    }
                                }
                                if (project_in_block) {
                                    // sum_nu phid[P][nu] Ca[nu][Ncore+j] is chid[P][Ncore+j],
                                    // so the active block is Sd^T chid_act directly.
                                    C_DGEMM('T','N', nact, nact, npoints, 1.0, Sd[0], nact,
                                            &chid[0][Ncore], nchi, 1.0, M3[0], ld_m3);
                                } else {
                                    C_DGEMM('T','N', nact, nlocal, npoints, 1.0, Sd[0], nact,
                                            phid[0], coll_funcs, 1.0, M3[0], ld_m3);
                                }
                            }
                            wa_set(L);
                            if (project_in_block) {
                                for (int i = 0; i < nact; ++i) {
                                    const double* r = M3[i];
                                    double* trow = waL + (size_t)i*nact;
                                    for (int j = 0; j < nact; ++j) trow[j] += r[j];
                                }
                            } else {
                                for (int i = 0; i < nact; ++i) {
                                    const double* r = M3[i];
                                    double* trow = waL + (size_t)i*nbf;
                                    for (int k = 0; k < nlocal; ++k) trow[function_map[k]] += r[k];
                                }
                            }
                            wa_unset(L);
                        };
                        meta_dir(vta, waAL, diAL);
                        if (!closed) {
                            double* diBL = need_diag ? (diB + (size_t)L*nmo) : nullptr;
                            meta_dir(vtb, waBL, diBL);
                        }
                        // E1 meta: STau += C_L v_tau w (beta == alpha for closed).
                        if (want_e1) for (int P = 0; P < npoints; ++P) {
                            stauA[P] += cl*vta[P]*w[P];
                            stauB[P] += cl*(closed ? vta[P] : vtb[P])*w[P];
                        }
                    }
                }

                flush_stage();
                sbase += ns;
            }

            // Full-MO diagonal: V_qq(slice) = 2 sum_{P,d} Vcof[slice][ncof P + d] Z[ncof P + d][q],
            // Z[ncof P + d][q] = chi[P][q] {chi, chix, chiy, chiz}_d[P][q], one GEMM per point chunk.
            if (need_diag && sbase > 0) {
                std::fill(Dd.begin(), Dd.begin() + (size_t)sbase*nmo, 0.0);
                for (int P0 = 0; P0 < npoints; P0 += kZChunk) {
                    const int np = std::min(kZChunk, npoints - P0);
                    for (int p = 0; p < np; ++p) {
                        const double* ch = chi[P0 + p];
                        double* z0 = Zd.data() + (size_t)(ncof*p + 0)*nmo;
                        for (int q = 0; q < nmo; ++q) z0[q] = ch[q]*ch[q];
                        if (ansatz >= 1) {
                            const double* cx = chix[P0 + p];
                            const double* cy = chiy[P0 + p];
                            const double* cz = chiz[P0 + p];
                            double* z1 = Zd.data() + (size_t)(ncof*p + 1)*nmo;
                            double* z2 = Zd.data() + (size_t)(ncof*p + 2)*nmo;
                            double* z3 = Zd.data() + (size_t)(ncof*p + 3)*nmo;
                            for (int q = 0; q < nmo; ++q) z1[q] = ch[q]*cx[q];
                            for (int q = 0; q < nmo; ++q) z2[q] = ch[q]*cy[q];
                            for (int q = 0; q < nmo; ++q) z3[q] = ch[q]*cz[q];
                        }
                    }
                    C_DGEMM('N','N', sbase, nmo, ncof*np, 2.0, &Vcof[0][(size_t)ncof*P0],
                            ncof*max_points, Zd.data(), nmo, 1.0, Dd.data(), nmo);
                }
                for (int sj = 0; sj < sbase; ++sj) {
                    double* diL = (slice_beta[sj] ? diB : diA) + (size_t)slice_L[sj]*nmo;
                    C_DAXPY(nmo, 1.0, Dd.data() + (size_t)sj*nmo, 1, diL, 1);
                }
            }

            // E1 back-projection of the C_L-weighted AO potential, once per block.
            for (int sp = 0; want_e1 && sp < 2; ++sp) {
                double** TapAcc = (sp == 0) ? TapAcc_a : TapAcc_b;
                double** Vp = (sp == 0) ? V_acc_a_p : V_acc_b_p;
                double* stau = (sp == 0) ? stauA : stauB;
                // V_loc = sym(phi^T TapAcc) (LSDA+GGA) + sum_d phid^T diag(STau) phid (meta).
                C_DGEMM('T','N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, TapAcc[0], max_functions,
                        0.0, Vloc_p[0], max_functions);
                for (int m = 0; m < nlocal; ++m)
                    for (int n = 0; n <= m; ++n) { Vloc_p[m][n] = Vloc_p[n][m] = Vloc_p[m][n] + Vloc_p[n][m]; }
                if (ansatz >= 2) {
                    double** phid_arr[3] = {phix, phiy, phiz};
                    for (int d = 0; d < 3; ++d) {
                        double** phid = phid_arr[d];
                        for (int P = 0; P < npoints; ++P)
                            for (int ml = 0; ml < nlocal; ++ml) MbpT[P][ml] = stau[P]*phid[P][ml];
                        C_DGEMM('T','N', nlocal, nlocal, npoints, 1.0, phid[0], coll_funcs, MbpT[0], max_functions,
                                1.0, Vloc_p[0], max_functions);
                    }
                }
                // The AO accumulator is shared; only the scatter of this block is serialized.
#pragma omp critical(reks_vacc)
                {
                    for (int ml = 0; ml < nlocal; ++ml) {
                        const int mg = function_map[ml];
                        for (int nl = 0; nl < nlocal; ++nl) Vp[mg][function_map[nl]] += Vloc_p[ml][nl];
                    }
                }
            }
        }

#pragma omp critical
        {
            // Per-L scalars and, when it is wanted, the MO diagonal; the weighted-AO rows and
            // the AO accumulator were already scattered into their shared buffers.
            for (int L = 0; L < n_micro; ++L) {
                E_xc[L] += excr[L]; tr_a[L] += trAr[L]; tr_b[L] += trBr[L];
                if (need_diag) {
                    C_DAXPY(nmo, 1.0, diA + (size_t)L*nmo, 1, diag_a[L].data(), 1);
                    C_DAXPY(nmo, 1.0, diB + (size_t)L*nmo, 1, diag_b[L].data(), 1);
                }
            }
        }
    }
    timer_off(t_grid);

    timer_on(t_proj);
    // arows = wa Ca over columns [mo_lo, mo_lo + mo_width), one GEMM per contiguous microstate run;
    // project_in_block copies the active block as is.
    sc.proj.resize(project_in_block ? 0 : (size_t)n_micro * nact * mo_width);
    std::vector<double>& proj = sc.proj;
    auto project_run = [&](double* wa, std::vector<SharedMatrix>& arows, int L0, int nL) {
        if (nL <= 0 || mo_width <= 0) return;
        if (project_in_block) {
            for (int L = 0; L < nL; ++L)
                C_DCOPY((size_t)nact*nact, wa + (size_t)(L0 + L)*nact*nact, 1,
                        arows[L0 + L]->pointer()[0], 1);
            return;
        }
        C_DGEMM('N','N', nL*nact, mo_width, nbf, 1.0, wa + (size_t)L0*nact*nbf, nbf,
                &Ca_p[0][mo_lo], nmo, 0.0, proj.data(), mo_width);
        for (int L = 0; L < nL; ++L)
            C_DCOPY((size_t)nact*mo_width, proj.data() + (size_t)L*nact*mo_width, 1,
                    arows[L0 + L]->pointer()[0], 1);
    };
    project_run(wa_a_all.data(), arows_a, 0, n_micro);
    // Closed-shell beta is mirrored from alpha below (its weighted-AO block stays zero);
    // project_run here runs only over the stretches of consecutive open microstates.
    for (int L0 = 0; L0 < n_micro;) {
        if (closed_ms[L0]) { ++L0; continue; }
        int L1 = L0;
        while (L1 < n_micro && !closed_ms[L1]) ++L1;
        project_run(wa_b_all.data(), arows_b, L0, L1 - L0);
        L0 = L1;
    }

    const double nq_total = std::accumulate(nqT.begin(), nqT.end(), 0.0);
    grid_electrons_ = nq_total;

    // Closed-shell microstates: mirror beta arows = alpha (and the full MO diagonal when need_diag).
    for (int L = 0; L < n_micro; ++L) {
        if (!closed_ms[L]) continue;
        arows_b[L]->copy(arows_a[L]);
        if (need_diag) diag_b[L] = diag_a[L];
    }

    for (int L = 0; L < n_micro; ++L)
        if (std::isnan(E_xc[L]))
            throw PSIEXCEPTION("MicrostateXC::compute: XC energy NaN.");

    timer_off(t_proj);
    timer_off("REKS: microstate xc");
}

}  // namespace xc
}  // namespace reks
}  // namespace psi
