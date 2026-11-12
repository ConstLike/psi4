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

#include "reks_properties.h"

#include "reks_arch.h"
#include "reks_naturals.h"
#include "reks_report_level.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libqt/qt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace psi {
namespace reks {
namespace properties {

using report::active_orbital_labels;
using report::reports;
using report::sector_psivar_name;

DipoleIntegralCache build_dipole_integral_cache(const report::ReportFrame& frame,
                                                std::shared_ptr<Molecule> molecule,
                                                std::shared_ptr<BasisSet> basisset) {
    DipoleIntegralCache dip;
    const int N_act  = frame.sa_cassette.n_active_orbitals();
    const size_t N2  = static_cast<size_t>(N_act) * N_act;

    auto mints = std::make_shared<MintsHelper>(basisset);
    std::vector<SharedMatrix> ao_dip = mints->ao_dipole();  // [0]=x, [1]=y, [2]=z

    for (int A = 0; A < molecule->natom(); ++A) {
        double Z = molecule->Z(A);
        dip.mu_nuc[0] += Z * molecule->x(A);
        dip.mu_nuc[1] += Z * molecule->y(A);
        dip.mu_nuc[2] += Z * molecule->z(A);
    }

    // AO -> MO dipole over Ca_fock_transform, the orbitals frame.active indexes.
    auto C = frame.reel.Ca_fock_transform ? frame.reel.Ca_fock_transform : frame.Ca;
    int nso = C->rowspi()[0];
    int nmo = C->colspi()[0];

    // Dipole per axis over the occupied MO columns [0, n_occ) only.
    const int n_occ = frame.active.empty() ? frame.Ncore : frame.active.back() + 1;
    double** Cp_d = C->pointer(0);
    std::vector<double> T(static_cast<size_t>(nso) * n_occ);
    for (int axis = 0; axis < 3; ++axis) {
        // T = D_ao C_occ   (nso x n_occ)
        C_DGEMM('N', 'N', nso, n_occ, nso, 1.0, ao_dip[axis]->pointer(0)[0], nso, Cp_d[0], nmo,
                0.0, T.data(), n_occ);
        // Core contribution: 2 * sum_{p<frame.Ncore} d_pp.
        for (int p = 0; p < frame.Ncore; ++p)
            dip.mu_core[axis] +=
                2.0 * C_DDOT(nso, &Cp_d[0][p], nmo, T.data() + p, n_occ);
        // d_act[axis][p*N_act+q] = c_{a_p} . (D_ao c_{a_q}).
        dip.d_act[axis].resize(N2);
        for (int p = 0; p < N_act; ++p)
            for (int q = 0; q < N_act; ++q)
                dip.d_act[axis][p * N_act + q] =
                    C_DDOT(nso, &Cp_d[0][frame.active[p]], nmo,
                           T.data() + frame.active[q], n_occ);
    }
    return dip;
}

SIProperties compute_state_properties(const report::ReportFrame& frame,
                                      const studio::Cassette& cassette,
                                      const SIResult& sir,
                                      const rdm::DiabaticRdm& rho_diab,
                                      const DipoleIntegralCache& dip) {
    timer_on("REKS: compute_state_properties");
    SIProperties out;

    // Without a diabatic RDM there are no state dipoles.
    if (!reks::studio::has_diabatic_rdm(frame.catalog)) {
        timer_off("REKS: compute_state_properties");
        return out;
    }

    const int n_rep   = sir.n_report();
    const int N_act   = frame.sa_cassette.n_active_orbitals();
    const int n_pairs = n_rep * (n_rep - 1) / 2;
    const size_t N2   = static_cast<size_t>(N_act) * N_act;

    // Adiabatic 1-RDM rho^(IJ)_adiab = sum_{K,L} c_I^K c_J^L rho^(KL)_diab.
    std::vector<double> rho_adiab =
        reks::rdm::adiabatic(rho_diab, cassette, sir);

    out.perm_dipole.assign(n_rep, {0.0, 0.0, 0.0});
    out.trans_dipole.assign(n_pairs, {0.0, 0.0, 0.0});
    out.osc_str.assign(n_pairs, 0.0);

    // mu_IJ[axis] = Tr[rho^(IJ)_adiab d_act[axis]] = sum_pq rho^(IJ)_pq d_pq.
    auto contract_axis = [&](int I, int J, int axis) -> double {
        const double* blk = &rho_adiab[(static_cast<size_t>(I) * n_rep + J) * N2];
        double s = 0.0;
        for (size_t pq = 0; pq < N2; ++pq) s += blk[pq] * dip.d_act[axis][pq];
        return s;
    };

    // Permanent dipoles: mu_nuc + mu_core + Tr[rho^(II) d_active].
    for (int I = 0; I < n_rep; ++I)
        for (int axis = 0; axis < 3; ++axis)
            out.perm_dipole[I][axis] =
                dip.mu_nuc[axis] + dip.mu_core[axis] + contract_axis(I, I, axis);

    // Transition dipoles mu_IJ; oscillator strength f = (2/3) dE |mu_IJ|^2.
    int idx = 0;
    for (int I = 0; I < n_rep; ++I) {
        for (int J = I + 1; J < n_rep; ++J) {
            for (int axis = 0; axis < 3; ++axis)
                out.trans_dipole[idx][axis] = contract_axis(I, J, axis);
            double dE  = sir.energies[J] - sir.energies[I];
            double mu2 = 0.0;
            for (int a = 0; a < 3; ++a)
                mu2 += out.trans_dipole[idx][a] * out.trans_dipole[idx][a];
            out.osc_str[idx] = (2.0 / 3.0) * dE * mu2;
            ++idx;
        }
    }

    out.rho_adiab = std::move(rho_adiab);

    timer_off("REKS: compute_state_properties");
    return out;
}

void print_state_dipoles(const SIProperties& props, const SIResult& sir) {
    if (!reports(2)) return;
    const int n_rep = static_cast<int>(props.perm_dipole.size());
    if (n_rep < 1) return;
    const int n_pairs = n_rep * (n_rep - 1) / 2;
    if (n_pairs <= 0) return;

    // PsiOutStream::Printf flushes the stream on every call, so lines are batched
    // into one call per kFlushBytes of text.
    constexpr size_t kFlushBytes = 64u << 10;
    std::string buf;
    buf.reserve(kFlushBytes + 128);
    char cell[128];
    auto emit = [&](bool force) {
        if (buf.size() >= kFlushBytes || (force && !buf.empty())) {
            outfile->Printf(buf);
            buf.clear();
        }
    };

    outfile->Printf("\n    %s Transition Dipole Moments (a.u.)\n",
                    reks::report::next_table_tag().c_str());
    outfile->Printf("    --------------------------------\n");
    outfile->Printf("      %-12s %12s %12s %12s %12s\n", "Pair", "mu_x", "mu_y", "mu_z", "|mu|");
    for (int I = 0; I < n_rep; I++) {
        for (int J = I + 1; J < n_rep; J++) {
            const int idx = reks::si_pair_index(I, J, n_rep);
            double mx = props.trans_dipole[idx][0];
            double my = props.trans_dipole[idx][1];
            double mz = props.trans_dipole[idx][2];
            double mag = std::sqrt(mx * mx + my * my + mz * mz);
            snprintf(cell, sizeof(cell), "      S%d -> S%-6d %12.6f %12.6f %12.6f %12.6f\n",
                     I, J, mx, my, mz, mag);
            buf.append(cell);
            emit(false);
        }
    }
    emit(true);
    outfile->Printf("    --------------------------------\n");

    outfile->Printf("\n    %s Oscillator Strengths\n", reks::report::next_table_tag().c_str());
    outfile->Printf("    --------------------\n");
    outfile->Printf("      %-12s %12s %12s\n", "Pair", "f", "dE (eV)");
    for (int I = 0; I < n_rep; I++) {
        for (int J = I + 1; J < n_rep; J++) {
            double dE = sir.energies[J] - sir.energies[I];
            snprintf(cell, sizeof(cell), "      S%d -> S%-6d %12.5f %12.4f\n",
                     I, J, props.osc_str[reks::si_pair_index(I, J, n_rep)], dE * pc_hartree2ev);
            buf.append(cell);
            emit(false);
        }
    }
    emit(true);
    outfile->Printf("    --------------------\n");
    outfile->Printf("\n");
}

void register_si_wfn_variables(Wavefunction& wfn, const report::ReportFrame& frame,
                               int K, int s,
                               const SIResult& sir,
                               const SIProperties& props,
                               bool primary_alias,
                               const rdm::DiabaticRdm& rho_diab) {
    // Array names must not end in "DIPOLE" (Psi4 reshapes to (3,)); use "DIPOLES".
    const std::string suf = " K=" + std::to_string(K);
    const int n     = sir.n_states;    // config-space dimension
    const int n_rep = sir.n_report();  // reported states
    if (n < 1) return;

    // set_array_variable clones, so the bare alias costs a second full copy of what it names.
    auto record = [&](const std::string& base, SharedMatrix M) {
        wfn.set_array_variable(sector_psivar_name(s, base + suf), M);
        if (primary_alias) wfn.set_array_variable(sector_psivar_name(s, base), M);
    };

    // Bulk arrays carry the K-suffixed name alone; the alias would clone them.
    auto record_k_only = [&](const std::string& base, SharedMatrix M) {
        wfn.set_array_variable(sector_psivar_name(s, base + suf), M);
    };

    // Arrays carrying the config axis (coefficients, Hamiltonian, overlap, diabatic RDM)
    // scale as O(n^2); every other array here is linear in n_rep.
    const bool export_config_space = reports(3);

    auto E = std::make_shared<Matrix>("SSR ENERGIES" + suf, 1, n_rep);
    for (int i = 0; i < n_rep; ++i) E->set(0, 0, i, sir.energies[i]);
    record("SSR ENERGIES", E);

    if (export_config_space) {
        // flat is row-major nrow x n; each of the nrow destination rows is one contiguous run.
        auto make_sq = [&](const std::string& name, const std::vector<double>& flat, int nrow) {
            auto M = std::make_shared<Matrix>(name + suf, nrow, n);
            for (int i = 0; i < nrow; ++i)
                C_DCOPY(n, const_cast<double*>(&flat[static_cast<size_t>(i) * n]), 1,
                        M->pointer(0)[i], 1);
            return M;
        };
        record_k_only("SSR COEFFICIENTS", make_sq("SSR COEFFICIENTS", sir.coeffs, sir.n_show(n)));
        // pre-diagonalization
        record_k_only("SSR HAMILTONIAN", make_sq("SSR HAMILTONIAN", sir.hamiltonian, n));

        // S is nonzero only where a catalog pair with a nonempty overlap channel joins two
        // configs; exported as (flat index, value) pairs rather than a dense n x n matrix.
        const reks::BlockedMatrix& S = sir.overlap;
        size_t nnz = 0;
        for (double v : S.data)
            if (v != 0.0) ++nnz;
        auto S_sp = std::make_shared<Matrix>("SSR OVERLAP SPARSE" + suf,
                                             static_cast<int>(nnz), 2);
        double** dst = S_sp->pointer(0);
        size_t j = 0;
        for (int c = 0; c < S.n_blocks(); ++c) {
            const std::vector<int>& idx = S.blocks[c];
            const int b = static_cast<int>(idx.size());
            const double* const tile = S.tile(c);
            for (int r = 0; r < b; ++r)
                for (int s = 0; s < b; ++s) {
                    const double v = tile[static_cast<size_t>(r) * b + s];
                    if (v == 0.0) continue;
                    // row-major i*n + j, exact below 2^53
                    dst[j][0] = static_cast<double>(static_cast<size_t>(idx[r]) * n + idx[s]);
                    dst[j][1] = v;
                    ++j;
                }
        }
        record_k_only("SSR OVERLAP SPARSE", S_sp);  // pre-diagonalization
    }

    // State properties (physical states only); empty for catalogs without a diabatic RDM.
    const bool have_props = !props.perm_dipole.empty();
    if (have_props) {
        auto pdm = std::make_shared<Matrix>("SSR PERMANENT DIPOLES" + suf, n_rep, 3);
        for (int I = 0; I < n_rep; ++I)
            for (int a = 0; a < 3; ++a)
                pdm->set(0, I, a, props.perm_dipole[I][a]);
        record("SSR PERMANENT DIPOLES", pdm);
    }

    // props.trans_dipole / props.osc_str are packed upper-triangular over n_rep.
    const int n_pairs = n_rep * (n_rep - 1) / 2;
    if (have_props && n_pairs > 0) {
        auto tdm = std::make_shared<Matrix>("SSR TRANSITION DIPOLES" + suf, n_pairs, 3);
        for (int I = 0; I < n_rep; ++I)
            for (int J = I + 1; J < n_rep; ++J) {
                const int idx = reks::si_pair_index(I, J, n_rep);
                for (int a = 0; a < 3; ++a) tdm->set(0, idx, a, props.trans_dipole[idx][a]);
            }
        record("SSR TRANSITION DIPOLES", tdm);

        auto osc = std::make_shared<Matrix>("SSR OSCILLATOR STRENGTHS" + suf, n_rep, n_rep);
        for (int I = 0; I < n_rep; ++I)
            for (int J = I + 1; J < n_rep; ++J) {
                const double f = props.osc_str[reks::si_pair_index(I, J, n_rep)];
                osc->set(0, I, J, f);
                osc->set(0, J, I, f);
            }
        record("SSR OSCILLATOR STRENGTHS", osc);
    }

    if (have_props && n_rep > 1) {
        auto ddm = std::make_shared<Matrix>("SSR DIFFERENCE DIPOLES" + suf, n_rep - 1, 3);
        for (int I = 1; I < n_rep; ++I)
            for (int a = 0; a < 3; ++a)
                ddm->set(0, I - 1, a, props.perm_dipole[I][a] - props.perm_dipole[0][a]);
        record("SSR DIFFERENCE DIPOLES", ddm);
    }

    // Diabatic 1-RDM rho^(KL) nonzeros as rows [(K n + L) N_act^2 + p N_act + q, value].
    if (export_config_space && reks::studio::has_diabatic_rdm(frame.catalog)) {
        const size_t N2 = static_cast<size_t>(rho_diab.n_active) * rho_diab.n_active;
        size_t nnz = 0;
        for (double v : rho_diab.val)
            if (v != 0.0) ++nnz;
        if (nnz > 0) {
            auto rho_mat = std::make_shared<Matrix>(
                "SSR 1-RDM DIABATIC SPARSE" + suf, static_cast<int>(nnz), 2);
            double** dst = rho_mat->pointer(0);
            size_t j = 0;
            for (int K = 0; K < rho_diab.n_si; ++K)
                for (size_t z = rho_diab.row_ptr[K]; z < rho_diab.row_ptr[K + 1]; ++z) {
                    if (rho_diab.val[z] == 0.0) continue;
                    const size_t flat =
                        (static_cast<size_t>(K) * rho_diab.n_si + rho_diab.col_L[z]) * N2 +
                        rho_diab.col_pq[z];
                    dst[j][0] = static_cast<double>(flat);
                    dst[j][1] = rho_diab.val[z];
                    ++j;
                }
            record_k_only("SSR 1-RDM DIABATIC SPARSE", rho_mat);
        }
    }

    // Adiabatic 1-RDM rho^(IJ): rows I*n_rep+J, cols p*N_act+q.
    if (reks::studio::has_diabatic_rdm(frame.catalog) && !props.rho_adiab.empty()) {
        const int N_act = frame.sa_cassette.n_active_orbitals();
        const size_t N2 = static_cast<size_t>(N_act) * N_act;
        auto rho_mat = std::make_shared<Matrix>(
            "SSR 1-RDM ADIABATIC" + suf, n_rep * n_rep, N_act * N_act);
        for (int I = 0; I < n_rep; ++I)
            for (int J = 0; J < n_rep; ++J)
                C_DCOPY(static_cast<int>(N2),
                        const_cast<double*>(
                            &props.rho_adiab[static_cast<size_t>(I * n_rep + J) * N2]), 1,
                        rho_mat->pointer(0)[I * n_rep + J], 1);
        record_k_only("SSR 1-RDM ADIABATIC", rho_mat);
    }

    // MOs (Ca_fock_transform, else frame.Ca) + Ncore.
    {
        auto C_path_b = frame.reel.Ca_fock_transform ? frame.reel.Ca_fock_transform : frame.Ca;
        record("REKS C PATH B", C_path_b);
        auto ncore_mat = std::make_shared<Matrix>("REKS NCORE" + suf, 1, 1);
        ncore_mat->set(0, 0, 0, static_cast<double>(frame.Ncore));
        record("REKS NCORE", ncore_mat);
    }

    // Per generation: SI-relaxed FONs of this sector and META rows (scheme, orb_p, orb_q).
    {
        const int ng = frame.sa_cassette.n_geminals();
        const reks::studio::FONSnapshot& si_fons = frame.reel.fon_state[s];
        for (int gen = 0; ng > 0 && gen < frame.sa_cassette.n_generations(); ++gen) {
            // Sector s populates only its own generations; skip layers it leaves empty.
            if (static_cast<int>(si_fons.layers[gen].size()) < ng) continue;
            const std::string base = reks::studio::fon_array_name(gen);
            auto fon_si = std::make_shared<Matrix>(base + suf, 1, ng);
            auto fon_si_meta = std::make_shared<Matrix>(base + " META" + suf, 3, ng);
            for (int g = 0; g < ng; ++g) {
                const auto& tpl = frame.sa_cassette.geminal_templates()[g];
                const auto& fon = si_fons.layers[gen][g];
                fon_si->set(0, 0, g, fon.p);
                fon_si_meta->set(0, 0, g, static_cast<double>(tpl.scheme));
                fon_si_meta->set(0, 1, g, static_cast<double>(tpl.orbitals[0]));
                fon_si_meta->set(0, 2, g, static_cast<double>(tpl.orbitals[1]));
            }
            wfn.set_array_variable(sector_psivar_name(s, base + suf), fon_si);
            wfn.set_array_variable(sector_psivar_name(s, base + " META" + suf), fon_si_meta);
        }
    }
}

StateNaturalOrbitals compute_state_natural(const studio::Catalog& catalog,
                                           const studio::Cassette& cassette,
                                           const SIResult& sir,
                                           const std::vector<double>& rho_adiab) {
    timer_on("REKS: compute_state_natural");
    if (!reks::studio::has_diabatic_rdm(catalog) || rho_adiab.empty()) {
        timer_off("REKS: compute_state_natural");
        return {};
    }
    auto result = reks::naturals::compute(rho_adiab, cassette, sir);
    timer_off("REKS: compute_state_natural");
    return result;
}

void print_state_natural(const report::ReportFrame& frame,
                         const StateNaturalOrbitals& nos, int K,
                         const studio::Cassette& cassette,
                         const SIResult& sir,
                         const rdm::DiabaticRdm& rho_diab, int spin2) {
    if (!reports(2)) return;
    // n_rep (nos.n_states) bounds the state axis; n_cfg (sir.n_states) is the
    // independent config axis.
    const int n_rep = nos.n_states;
    const int n_cfg = sir.n_states;
    const int N    = nos.n_active;
    if (n_rep < 1 || N < 1) return;
    // active_orbital_labels(frame) dereferences frame.epsilon.
    if (!frame.epsilon) return;

    std::string ml = reks::report::method_label(
        frame.sa_cassette,
        static_cast<int>(cassette.K_indices.size()));
    if (frame.sa_cassette.n_sectors() > 1) {
        ml += "  " + reks::report::spin_tag(spin2);
    }

    const std::vector<int>& slice_K = cassette.K_indices;
    auto slice_name = [&](int J) -> std::string {
        const int idx = (static_cast<int>(slice_K.size()) == n_cfg) ? slice_K[J] : J;
        return reks::report::si_config_name(frame.sa_cassette, idx);
    };

    // Frame width = 31 char label + 14 per coeff cell, over the columns shown
    // per block (n_cfg capped at kSIBlockColumns; wider n_cfg wraps to more blocks).
    const int ncol_si = std::min(n_cfg, reks::kSIBlockColumns);
    const int frame_w = std::max(60, 31 + ncol_si * 14);
    const std::string heq(static_cast<size_t>(frame_w), '=');
    const std::string hdh(static_cast<size_t>(frame_w), '-');

    outfile->Printf("\n  %s\n", heq.c_str());
    outfile->Printf("    SSR Natural Occupations / Orbitals (unrelaxed)  K=%d  %s\n",
                    K, ml.c_str());
    outfile->Printf("  %s\n", hdh.c_str());

    // Table A coeffs[k*n_cfg+J]: rows = presented states, cols = configs.
    if (!sir.coeffs.empty() &&
        sir.coeffs.size() == static_cast<size_t>(sir.n_show(n_cfg)) * n_cfg) {
        const std::string tag_a = reks::report::next_table_tag();
        if (sir.report_cap == 0)
            outfile->Printf("\n    %s Table A: SI eigenvector composition (full)\n",
                            tag_a.c_str());
        else
            outfile->Printf("\n    %s Table A: SI eigenvector composition (first %d state(s))\n",
                            tag_a.c_str(), n_rep);
        outfile->Printf("    ------------------------------------------\n");

        const int ncol = std::max(1, reks::kSIBlockColumns);
        const int nblocks = (n_cfg + ncol - 1) / ncol;

        std::vector<std::string> col_names(n_cfg);
        for (int J = 0; J < n_cfg; ++J) col_names[J] = slice_name(J);

        // One Printf per row: PsiOutStream::Printf flushes the stream on every call.
        std::string line;
        line.reserve(static_cast<size_t>(ncol + 2) * 16);
        char cell[64];

        for (int b = 0; b < nblocks; ++b) {
            const int j0 = b * ncol;
            const int j1 = std::min(j0 + ncol, n_cfg);

            snprintf(cell, sizeof(cell), "      %-3s  %-20s", "k", "E_k (Ha)");
            line.assign(cell);
            for (int J = j0; J < j1; ++J) {
                snprintf(cell, sizeof(cell), "  %12s", col_names[J].c_str());
                line.append(cell);
            }
            line.push_back('\n');
            outfile->Printf(line);

            for (int k = 0; k < n_rep; ++k) {
                const double E_k =
                    (k < static_cast<int>(sir.energies.size())) ? sir.energies[k] : 0.0;
                snprintf(cell, sizeof(cell), "      %-3d  %+20.12f", k, E_k);
                line.assign(cell);
                const double* row = &sir.coeffs[static_cast<size_t>(k) * n_cfg];
                for (int J = j0; J < j1; ++J) {
                    snprintf(cell, sizeof(cell), "  %+12.6f", row[J]);
                    line.append(cell);
                }
                line.push_back('\n');
                outfile->Printf(line);
            }
            if (b + 1 < nblocks) outfile->Printf("\n");
        }
        outfile->Printf("    ------------------------------------------\n");
    }

    outfile->Printf("\n    %s Table B: per-state FONs and natural orbitals\n",
                    reks::report::next_table_tag().c_str());
    outfile->Printf("    --------------------------------------------\n");
    outfile->Printf("      columns = natural orbitals alpha with occupation n_alpha\n");
    outfile->Printf("      rows = coefficient on each active SCF MO by energy-rank label\n");
    // Slot-letter -> energy-rank key; state-independent, printed once.
    const std::vector<std::string> rk_lab = active_orbital_labels(frame);
    outfile->Printf("      active SCF MO rows: ");
    for (int p = 0; p < N; ++p)
        outfile->Printf(" %s=%s", reks::report::orbital_label(p).c_str(), rk_lab[p].c_str());
    outfile->Printf("\n");
    for (int k = 0; k < n_rep; ++k) {
        outfile->Printf("\n      SSR state %d:\n", k);

        // Columns ordered by greedy decreasing |c| (dominant coefficients on the diagonal, ties row-major).
        std::vector<int> col_order(N);
        {
            const double* blk = &nos.orbitals[static_cast<size_t>(k) * N * N];
            std::vector<std::pair<int, int>> cells;
            cells.reserve(static_cast<size_t>(N) * N);
            for (int p = 0; p < N; ++p)
                for (int alpha = 0; alpha < N; ++alpha) cells.emplace_back(p, alpha);
            std::stable_sort(cells.begin(), cells.end(),
                             [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                                 return std::abs(blk[a.first * N + a.second]) >
                                        std::abs(blk[b.first * N + b.second]);
                             });
            std::vector<char> row_taken(N, 0), col_taken(N, 0);
            for (const auto& cell : cells) {
                if (row_taken[cell.first] || col_taken[cell.second]) continue;
                col_order[cell.first] = cell.second;
                row_taken[cell.first] = 1;
                col_taken[cell.second] = 1;
            }
        }

        outfile->Printf("        %-12s", "n_alpha:");
        for (int alpha = 0; alpha < N; ++alpha)
            outfile->Printf("  %10.6f", nos.occupations[k * N + col_order[alpha]]);
        outfile->Printf("\n");

        // Row p = active SCF MO; column alpha = natural orbital (occupation n_alpha above).
        for (int p = 0; p < N; ++p) {
            outfile->Printf("        MO %-9s", rk_lab[p].c_str());
            for (int alpha = 0; alpha < N; ++alpha) {
                const double cpa =
                    nos.orbitals[k * N * N + p * N + col_order[alpha]];
                outfile->Printf("  %+10.4f", cpa);
            }
            outfile->Printf("\n");
        }
    }
    outfile->Printf("    --------------------------------------------\n");

    if (reports(5)) {
        outfile->Printf("\n  [DEBUG rho^(k) before diag, K=%d]\n", K);
        for (int k = 0; k < n_rep; ++k) {
            outfile->Printf("    state %2d:\n", k);
            for (int p = 0; p < N; ++p) {
                outfile->Printf("     ");
                for (int q = 0; q < N; ++q) {
                    outfile->Printf(" %+12.8f",
                                    nos.rdm[k * N * N + p * N + q]);
                }
                outfile->Printf("\n");
            }
        }

        // CSR row entries are sorted ascending in (L, pq): each block is the contiguous
        // run of equal L, in J order. Unbuilt for catalogs without a diabatic RDM.
        outfile->Printf("\n  [DEBUG diabatic rho^(IJ), nonzero blocks, K=%d]\n", K);
        const size_t N2 = static_cast<size_t>(N) * N;
        std::vector<double> blk(N2);
        for (int I = 0; I < rho_diab.n_si; ++I) {
            size_t z = rho_diab.row_ptr[I];
            while (z < rho_diab.row_ptr[I + 1]) {
                const int J = rho_diab.col_L[z];
                std::fill(blk.begin(), blk.end(), 0.0);
                double bmax = 0.0;
                for (; z < rho_diab.row_ptr[I + 1] && rho_diab.col_L[z] == J; ++z) {
                    blk[rho_diab.col_pq[z]] = rho_diab.val[z];
                    bmax = std::max(bmax, std::abs(rho_diab.val[z]));
                }
                if (bmax < 1.0e-14) continue;
                outfile->Printf("    rho^[%s,%s]:\n",
                                slice_name(I).c_str(), slice_name(J).c_str());
                for (int p = 0; p < N; ++p) {
                    outfile->Printf("     ");
                    for (int q = 0; q < N; ++q) {
                        outfile->Printf(" %+12.8f", blk[p * N + q]);
                    }
                    outfile->Printf("\n");
                }
            }
        }
    }

    // Integrity: trace = n_act_e, asym ~ 1e-14, occ in [0,2].
    outfile->Printf("\n    %s Integrity check\n", reks::report::next_table_tag().c_str());
    outfile->Printf("    ---------------\n");
    const double n_act_e = static_cast<double>(frame.sa_cassette.n_electrons());
    for (int k = 0; k < n_rep; ++k) {
        double trace = 0.0;
        for (int p = 0; p < N; ++p) trace += nos.rdm[k * N * N + p * N + p];

        double max_asym = 0.0;
        for (int p = 0; p < N; ++p) {
            for (int q = p + 1; q < N; ++q) {
                const double a = nos.rdm[k * N * N + p * N + q];
                const double b = nos.rdm[k * N * N + q * N + p];
                max_asym = std::max(max_asym, std::abs(a - b));
            }
        }

        double n_min = +2.0, n_max = -1.0;
        for (int alpha = 0; alpha < N; ++alpha) {
            n_min = std::min(n_min, nos.occupations[k * N + alpha]);
            n_max = std::max(n_max, nos.occupations[k * N + alpha]);
        }

        outfile->Printf("      state %2d: trace=%.10f (expect %.1f)  asym=%.2e  "
                        "n in [%+0.4f, %+0.4f]\n",
                        k, trace, n_act_e, max_asym, n_min, n_max);
    }
    outfile->Printf("    ---------------\n");

    outfile->Printf("  %s\n", heq.c_str());
    outfile->Printf("  End of %s SSR output.\n", ml.c_str());
    outfile->Printf("  %s\n\n", heq.c_str());
}

void register_natural_wfn_variables(Wavefunction& wfn, const report::ReportFrame& frame,
                                    int K, int s,
                                    const StateNaturalOrbitals& nos,
                                    bool primary_alias) {
    const int n_rep = nos.n_states;
    const int N     = nos.n_active;
    if (n_rep < 1 || N < 1) return;

    const std::string suf_K = " K=" + std::to_string(K);

    auto record = [&](const std::string& base_with_state, SharedMatrix M) {
        wfn.set_array_variable(sector_psivar_name(s, base_with_state + suf_K), M);
        if (primary_alias) wfn.set_array_variable(sector_psivar_name(s, base_with_state), M);
    };

    {
        auto occ = std::make_shared<Matrix>("SSR NATURAL OCCUPATIONS" + suf_K,
                                            n_rep, N);
        for (int k = 0; k < n_rep; ++k)
            for (int alpha = 0; alpha < N; ++alpha)
                occ->set(0, k, alpha, nos.occupations[k * N + alpha]);
        record("SSR NATURAL OCCUPATIONS", occ);
    }

    // Per-state matrices: rho (active), U (active), C * U (AO).
    auto C_full = frame.reel.Ca_fock_transform ? frame.reel.Ca_fock_transform : frame.Ca;
    const int nso = C_full->rowspi()[0];

    auto C_active = std::make_shared<Matrix>("C_active", nso, N);
    {
        double** Cp = C_full->pointer(0);
        for (int p = 0; p < N; ++p) {
            const int mo = frame.active[p];
            for (int mu = 0; mu < nso; ++mu) {
                C_active->set(0, mu, p, Cp[mu][mo]);
            }
        }
    }

    for (int k = 0; k < n_rep; ++k) {
        const std::string suf_S = " STATE=" + std::to_string(k);

        // rho^(k): active-basis pre-diag, symmetric.
        auto rho_k = std::make_shared<Matrix>(
            "SSR 1-RDM ACTIVE" + suf_K + suf_S, N, N);
        for (int p = 0; p < N; ++p)
            for (int q = 0; q < N; ++q)
                rho_k->set(0, p, q, nos.rdm[k * N * N + p * N + q]);
        record("SSR 1-RDM ACTIVE" + suf_S, rho_k);

        // U^(k): natural orbitals in active MO basis.
        auto Uk = std::make_shared<Matrix>(
            "SSR NATURAL ORBITALS ACTIVE" + suf_K + suf_S, N, N);
        for (int p = 0; p < N; ++p)
            for (int alpha = 0; alpha < N; ++alpha)
                Uk->set(0, p, alpha, nos.orbitals[k * N * N + p * N + alpha]);
        record("SSR NATURAL ORBITALS ACTIVE" + suf_S, Uk);

        // C * U^(k): natural orbitals expanded in AO basis.
        auto C_NO_k = std::make_shared<Matrix>(
            "SSR NATURAL ORBITALS AO" + suf_K + suf_S, nso, N);
        C_NO_k->gemm(false, false, 1.0, C_active, Uk, 0.0);
        record("SSR NATURAL ORBITALS AO" + suf_S, C_NO_k);
    }
}

}  // namespace properties
}  // namespace reks
}  // namespace psi
