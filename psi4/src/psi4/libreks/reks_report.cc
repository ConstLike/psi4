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

#include "reks_report.h"
#include "reks_scf.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/wavefunction.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace psi {
namespace reks {
namespace report {

std::string input_usage_examples() {
    return
"\n"
"  How to set up a REKS calculation ([N, M] = active electrons,\n"
"  active orbitals):\n"
"\n"
"  Python API:\n"
"      psi4.set_options({\n"
"          'reference':       'reks',\n"
"          'reks':            [4, 4],\n"
"          'sa_reks_configs': ['PPS1', 'OSS1', 'OSS2'],\n"
"          'sa_reks_weights': [0.5, 0.25, 0.25],\n"
"          'si_reks_configs': [['PPS1', 'OSS1'], ['PPS2', 'OSS3']],\n"
"      })\n"
"      psi4.energy('scf')\n"
"\n"
"  PSithon (.dat) input:\n"
"      set {\n"
"          reference        reks\n"
"          reks             [ 4, 4 ]\n"
"          sa_reks_configs  [ PPS1, OSS1, OSS2 ]\n"
"          sa_reks_weights  [ 0.5, 0.25, 0.25 ]\n"
"          si_reks_configs  [ [PPS1, OSS1], [PPS2, OSS3] ]\n"
"      }\n"
"      energy('scf')\n"
"\n"
"  Optional: sa_reks_configs and si_reks_configs default to the variant's\n"
"  built-in pool; sa_reks_weights defaults to a uniform average; sa_reks_extra\n"
"  adds raw determinants [alpha(M)|beta(M)] to the SA ensemble only.\n"
"\n"
"  Bulk tokens (in place of indices/names): 'full' selects the whole manifold\n"
"  config block; a config type ('PPS', 'OSS', ...) selects every config of that\n"
"  type, e.g. si_reks_configs [['full']] or [['PPS', 'OSS']].\n"
"  'single-scheme:N' selects the configs of pairing scheme N and 'mixed-scheme:N'\n"
"  the SSR state set on scheme N (open-shell singlets from every scheme); the\n"
"  spin-projected configs are the base and ride along with either. Both default to 0\n"
"  and take several indices ('single-scheme-1-2-3' unions the three) or scheme\n"
"  classes relative to scheme 0: 'T' (one transposition of antibonding partners),\n"
"  'D<k>' (k disjoint transpositions), 'C<l>' (an l-cycle), e.g.\n"
"  si_reks_configs [[['mixed-scheme', 'single-scheme-0-T-D2']]]. Quote them, and\n"
"  in a PSithon input write the indices behind '-' or '_', since the PSithon array\n"
"  parser splits an element on ':'.\n"
"  'sps' selects the spin-projected base and a pattern with '*' or '?' the names\n"
"  it matches. In an SI cassette, 'exclude:X' subtracts what X names from the rest\n"
"  of the cassette, e.g. si_reks_configs [[['full', 'exclude:sps']]].\n";
}

void input_error(const std::string& detail, bool show_usage) {
    const std::string usage = show_usage ? input_usage_examples() : std::string();
    outfile->Printf("\n");
    outfile->Printf(
        "  ========================= REKS INPUT ERROR =========================\n");
    outfile->Printf("  %s\n", detail.c_str());
    if (show_usage) outfile->Printf("%s", usage.c_str());
    outfile->Printf(
        "  ====================================================================\n\n");
    throw PSIEXCEPTION("REKS input error: " + detail + usage);
}

namespace {

// Report printing is serial; no synchronization needed.
int g_table_tag = 1;

}  // namespace

std::string next_table_tag() {
    return "[tag" + std::to_string(g_table_tag++) + "]";
}

void reset_table_tags() { g_table_tag = 1; }

std::string orbital_label(int p) {
    if (p < 0) return "?";
    return std::string(1, static_cast<char>('a' + p));
}

std::string microstate_label(const studio::Cassette& sa_cassette, int L) {
    if (L < 0 || L >= sa_cassette.n_total_microstates()) return "";
    const auto& ms = sa_cassette.microstate(L);
    const int N = sa_cassette.n_active_orbitals();
    std::string s = "|";
    bool first = true;
    for (int i = 0; i < N; ++i) {
        if (ms.alpha[i]) {
            if (!first) s += ",";
            s += orbital_label(i);
            first = false;
        }
    }
    for (int i = 0; i < N; ++i) {
        if (ms.beta[i]) {
            if (!first) s += ",";
            s += orbital_label(i);
            s += "_";
            first = false;
        }
    }
    s += "|";
    return s;
}

std::string si_config_name(const studio::Cassette& cassette, int K) {
    return cassette.si_config_name(K);
}

std::string method_label(const studio::Cassette& sa_cassette, int n_si) {
    const int n_sa = static_cast<int>(sa_cassette.K_indices.size())
                   + sa_cassette.n_extra_microstates();
    // N == M assumed (electrons not tracked separately): both REKS(.,.) slots use M.
    const int M = sa_cassette.n_active_orbitals();
    std::string s;
    if (n_si > 0 && n_si != n_sa) s = std::to_string(n_si) + "SI-";
    s += std::to_string(n_sa) + "SA-REKS(" + std::to_string(M) + "," +
         std::to_string(M) + ")";
    return s;
}

std::string spin_name(int spin2) {
    static const char* const names[] = {"singlet", "doublet", "triplet",
                                         "quartet", "quintet", "sextet", "septet"};
    const int mult = spin2 + 1;
    if (mult >= 1 && mult <= 7) return names[mult - 1];
    return "2S+1=" + std::to_string(mult);
}

std::string spin_tag(int spin2) {
    std::string name = spin_name(spin2);
    if (!name.empty() && name[0] >= 'a' && name[0] <= 'z')
        name[0] = static_cast<char>(name[0] - 'a' + 'A');
    return "2S=" + std::to_string(spin2) + " (" + name + ")";
}

namespace {

std::string config_def(const studio::Cassette& c, int K) {
    const char* d = c.si_config_def(K);
    return (d != nullptr) ? std::string(d) : std::string();
}

}  // namespace

void print_studio_banner() {
    if (!reports(2)) return;
    outfile->Printf(
        "   =========================================================================\n\n");
    outfile->Printf("                                REKS Studio\n\n");
    outfile->Printf("           Spin-Restricted Ensemble-referenced Kohn-Sham DFT\n");
    outfile->Printf("             State-Averaged (SA) and State-Interaction (SI)\n\n");
    outfile->Printf("            Konstantin Komarov, Michael Filatov, Seung Kyu Min\n\n");
    outfile->Printf(
        "    Ulsan National Institute of Science and Technology (UNIST), South Korea\n");
    outfile->Printf(
        "   =========================================================================\n\n");
}

void print_studio_setup(const studio::Cassette& sa, int spin, int n_core,
                        const std::vector<int>& active_mo_indices) {
    if (!reports(2)) return;
    const int N = sa.n_electrons();
    const int M = sa.n_active_orbitals();
    outfile->Printf("    REKS Studio Setup\n");
    outfile->Printf("    -----------------\n");
    outfile->Printf("      Variant            REKS(%d,%d)  scheme %d  spin %d  (%s)\n",
                    N, M, sa.scheme(), spin, spin_name(spin).c_str());
    outfile->Printf("      Active space       %d electrons in %d orbitals,  core %d\n",
                    N, M, n_core);
    // MO indices print 1-based; irrep is always A (REKS runs in C1).
    outfile->Printf("      Active MO indices  ");
    for (int idx : active_mo_indices) outfile->Printf("%dA ", idx + 1);
    outfile->Printf("\n\n");
}

void print_sector_map(const std::vector<int>& sector_2S,
                      const std::vector<int>& sa_counts,
                      const std::vector<int>& si_counts,
                      const std::vector<int>& n_generations) {
    if (!reports(2)) return;
    outfile->Printf("    REKS Sector Map\n");
    outfile->Printf("    ---------------\n");
    for (size_t s = 0; s < sector_2S.size(); ++s) {
        const std::string sa = sa_counts[s] > 0
            ? "SA configs " + std::to_string(sa_counts[s]) : std::string("SA-empty");
        const std::string si = si_counts[s] > 0 ? std::to_string(si_counts[s])
                                                : std::string("-");
        outfile->Printf("      sector %d : 2S=%d  %-13s SI cassettes %-3s generations %d\n",
                        static_cast<int>(s), sector_2S[s], sa.c_str(), si.c_str(),
                        n_generations[s]);
    }
    outfile->Printf("\n");
}

void print_pairing_schemes(const studio::Cassette& sa, bool show_triplet) {
    if (!reports(2)) return;
    outfile->Printf("    Pairing schemes (GVB(2,2) geminal pairs)\n\n");

    const studio::GeminalTemplate* tpl = sa.geminal_templates();
    const int ng = sa.n_geminals();
    std::vector<int> scheme_order;
    for (int g = 0; g < ng; ++g) {
        const int s = tpl[g].scheme;
        if (std::find(scheme_order.begin(), scheme_order.end(), s) == scheme_order.end())
            scheme_order.push_back(s);
    }
    for (int s : scheme_order) {
        std::string pairs;
        bool first = true;
        for (int g = 0; g < ng; ++g) {
            if (tpl[g].scheme != s) continue;
            if (!first) pairs += " & ";
            pairs += "(" + orbital_label(tpl[g].orbitals[0]) + "," +
                     orbital_label(tpl[g].orbitals[1]) + ")";
            first = false;
        }
        outfile->Printf("      scheme %d : %s\n", s, pairs.c_str());
    }
    outfile->Printf("\n");

    outfile->Printf("        GVB(2,2) wavefunctions for pair (p,q)\n");
    outfile->Printf("          Singlet closed-shell\n");
    outfile->Printf("            Psi_0(p,q) = +1/sqrt(2.0) * sqrt(n0_p) * |p,p_|\n");
    outfile->Printf("                         -1/sqrt(2.0) * sqrt(n0_q) * |q,q_|\n");
    outfile->Printf("            Psi_2(p,q) = +1/sqrt(2.0) * sqrt(n0_q) * |p,p_|\n");
    outfile->Printf("                         +1/sqrt(2.0) * sqrt(n0_p) * |q,q_|\n");
    outfile->Printf("          Singlet open-shell\n");
    outfile->Printf("            Psi_1(p,q) = +1/sqrt(2.0) * |p,q_|\n");
    outfile->Printf("                         +1/sqrt(2.0) * |q,p_|\n");
    if (show_triplet) {
        outfile->Printf("          Triplet open-shell\n");
        outfile->Printf("            Psi_To(p,q) = +1/sqrt(2.0) * |p,q_|\n");
        outfile->Printf("                          -1/sqrt(2.0) * |q,p_|  [ms=0]\n");
    }
    outfile->Printf("          FON constraint: n_p + n_q = 2.0\n\n");
}

void print_sa_pool(const studio::Cassette& sa_cassette) {
    if (!reports(2)) return;
    outfile->Printf("    SA pool\n");
    outfile->Printf("    -------\n\n");

    const studio::Cassette& sa = sa_cassette;

    struct Row { std::string k, name, def, ename; double w; };
    std::vector<Row> rows;
    double wsum = 0.0;
    for (int K : sa.K_indices) {
        const double w = sa.w_configs[K];
        wsum += w;
        const std::string name = si_config_name(sa, K);
        rows.push_back({std::to_string(K), name, config_def(sa, K), name, w});
    }
    for (int e = 0; e < sa_cassette.n_extra_microstates(); ++e) {
        const double w = sa_cassette.extra_weights[e];
        wsum += w;
        const std::string label = microstate_label(sa_cassette, sa_cassette.n_catalog_microstates() + e);
        rows.push_back({"", "Extra" + std::to_string(e + 1),
                        "A[(core) " + label + "]", label, w});
    }

    int wname = 6;  // at least as wide as the "config" header
    for (const Row& r : rows) wname = std::max(wname, static_cast<int>(r.name.size()));

    outfile->Printf("      %-4s %-*s %-11s %s\n", "K", wname, "config", "weight", "definition");
    for (const Row& r : rows)
        outfile->Printf("      %-4s %-*s %-11.6f %s\n", r.k.c_str(), wname, r.name.c_str(),
                        r.w, r.def.c_str());
    outfile->Printf("      %-4s %-*s ------------\n", "", wname, "");
    outfile->Printf("      %-4s %-*s sum = %.6f\n\n", "", wname, "", wsum);

    for (size_t i = 0; i < rows.size(); ++i)
        outfile->Printf(i == 0 ? "      E[SA-REKS] = %.6f * E(%s)\n"
                               : "                 + %.6f * E(%s)\n",
                        rows[i].w, rows[i].ename.c_str());
    outfile->Printf("\n");
}

void print_sa_energy_decomposition(const studio::Cassette&    sa_cassette,
                                   const std::vector<double>& e_config,
                                   const std::vector<double>& e_extra,
                                   double e_sa) {
    if (!reports(2)) return;
    const studio::Cassette& sa = sa_cassette;
    struct Row { std::string name; double w, e; };
    std::vector<Row> rows;
    for (size_t i = 0; i < sa.K_indices.size(); ++i) {
        const int K = sa.K_indices[i];
        rows.push_back({si_config_name(sa, K), sa.w_configs[K],
                        i < e_config.size() ? e_config[i] : 0.0});
    }
    for (int e = 0; e < sa_cassette.n_extra_microstates(); ++e) {
        const std::string label = microstate_label(sa_cassette, sa_cassette.n_catalog_microstates() + e);
        rows.push_back({label, sa_cassette.extra_weights[e],
                        e < static_cast<int>(e_extra.size()) ? e_extra[e] : 0.0});
    }
    if (rows.empty()) return;

    int wname = 6;  // at least as wide as the "config" header
    for (const Row& r : rows) wname = std::max(wname, static_cast<int>(r.name.size()));

    // Field widths: name(wname) sp weight(11) sp E(20) sp w*E(20).
    const int prefix_w = wname + 1 + 11 + 1 + 20 + 1;
    const int rule_w   = prefix_w + 20;
    const std::string rule(rule_w, '-');

    outfile->Printf("\n   => SA-REKS Energy Decomposition <=\n\n");
    outfile->Printf("    %-*s %11s %20s %20s\n", wname, "config", "weight",
                    "E(config) [Eh]", "w * E(config)");
    outfile->Printf("    %s\n", rule.c_str());
    for (const Row& r : rows)
        outfile->Printf("    %-*s %11.6f %20.12f %20.12f\n", wname, r.name.c_str(),
                        r.w, r.e, r.w * r.e);
    outfile->Printf("    %s\n", rule.c_str());
    outfile->Printf("    %-*s%20.12f\n", prefix_w, "E[SA-REKS]", e_sa);
    outfile->Printf("\n");
}

void print_si_pool(const std::vector<studio::Cassette>& si_cassettes,
                   const studio::Cassette& sa_cassette,
                   const std::vector<int>& sector_2S) {
    if (!reports(2)) return;
    std::vector<int> sa_set(sa_cassette.K_indices);
    std::sort(sa_set.begin(), sa_set.end());
    const bool multi_sector = sa_cassette.n_sectors() > 1;

    std::vector<size_t> show;
    for (size_t e = 0; e < si_cassettes.size(); ++e) {
        std::vector<int> ks(si_cassettes[e].K_indices);
        std::sort(ks.begin(), ks.end());
        if (ks != sa_set) show.push_back(e);
    }
    if (show.empty()) return;

    outfile->Printf("    SI pool\n");
    outfile->Printf("    -------\n\n");
    for (size_t e : show) {
        const studio::Cassette& c = si_cassettes[e];
        const int dim = static_cast<int>(c.K_indices.size());
        std::string label = method_label(sa_cassette, dim);
        if (multi_sector) {
            const int s = c.sector();
            const int spin2 = (s >= 0 && s < static_cast<int>(sector_2S.size()))
                                  ? sector_2S[s] : 0;
            label += "  " + spin_tag(spin2);
        }
        std::string cap_note;
        if (c.report_states > 0)
            cap_note = "  report: first " + std::to_string(c.report_states) + " state(s)";
        outfile->Printf("      cassette  dim  label\n");
        outfile->Printf("      %-8zu  %-4d %s%s\n", e, dim, label.c_str(), cap_note.c_str());

        int wname = 6;  // at least as wide as the "config" header
        for (int K : c.K_indices)
            wname = std::max(wname, static_cast<int>(si_config_name(c, K).size()));
        outfile->Printf("\n      %-4s %-*s %s\n", "K", wname, "config", "definition");
        for (int K : c.K_indices)
            outfile->Printf("      %-4d %-*s %s\n", K, wname, si_config_name(c, K).c_str(),
                            config_def(c, K).c_str());
        outfile->Printf("\n");
    }
}

void print_adiabatic_energies(const reks::SIResult& sir) {
    if (!reports(2)) return;
    const int n_print = sir.n_report();
    if (n_print < 1) return;

    outfile->Printf("\n    %s Adiabatic State Energies\n", next_table_tag().c_str());
    outfile->Printf("    ------------------------\n");
    for (int i = 0; i < n_print; ++i) {
        outfile->Printf("      SSR State S%-2d  %20.12f Ha\n", i, sir.energies[i]);
    }
    const int n_null = sir.n_states - sir.n_physical;
    if (n_null > 0) {
        outfile->Printf("\n      Note: %d of %d SI configurations linearly dependent "
                        "(overlap null-space); %d root(s) removed, %d physical state(s).\n",
                        n_null, sir.n_states, n_null, sir.n_physical);
    }
    if (n_print < sir.n_physical && sir.report_cap > 0) {
        outfile->Printf("\n      Note: reporting the first %d of %d physical state(s) "
                        "(SI_REKS_REPORT_STATES).\n", n_print, sir.n_physical);
    }
    outfile->Printf("    ------------------------\n");
}

void print_excitation_energies(const reks::SIResult& sir) {
    if (!reports(2)) return;
    const int n_print = sir.n_report();
    if (n_print < 2) return;

    outfile->Printf("\n    %s Excitation Energies (from S0)\n", next_table_tag().c_str());
    outfile->Printf("    -----------------------------\n");
    for (int i = 1; i < n_print; ++i) {
        const double dE = sir.energies[i] - sir.energies[0];
        outfile->Printf("      S0 -> S%-2d: Delta E = %12.6f Ha = %12.6f eV\n",
                        i, dE, dE * 27.21138624598);
    }
    outfile->Printf("    -----------------------------\n");
}

void print_microstate_energies(const studio::Cassette&      sa_cassette,
                               const std::vector<int>&    studio_L_indices,
                               const std::vector<double>& E_L) {
    if (!reports(2)) return;
    const int n = static_cast<int>(studio_L_indices.size());
    if (n == 0) return;

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return E_L[studio_L_indices[a]] < E_L[studio_L_indices[b]]; });

    const int n_rows = reports(3) ? n : std::min(n, kMicrostateReportRows);

    std::vector<std::string> labels(n_rows);
    int max_w = 0;
    for (int k = 0; k < n_rows; ++k) {
        labels[k] = microstate_label(sa_cassette, studio_L_indices[order[k]]);
        max_w = std::max(max_w, static_cast<int>(labels[k].size()));
    }

    outfile->Printf("    %s Microstate Energies (sorted by energy)\n", next_table_tag().c_str());
    outfile->Printf("    --------------------------------------\n");
    for (int k = 0; k < n_rows; ++k) {
        outfile->Printf("      E[%-*s] = %20.12f Ha\n",
                        max_w, labels[k].c_str(), E_L[studio_L_indices[order[k]]]);
    }
    if (n_rows < n)
        outfile->Printf("      ... %d more not shown; needs REKS_REPORT_LEVEL 3\n", n - n_rows);
    outfile->Printf("    --------------------------------------\n\n");
}

void print_lagrangian_eri_table(const studio::Cassette&               cassette,
                                const std::vector<double>&          lagrangians,
                                const ActiveEriTiles&               pair_eris) {
    if (!reports(2)) return;
    const int N = cassette.n_active_orbitals();

    const std::vector<studio::LagrangianPair> lagrangian_pairs =
        cassette.lagrangian_pairs_view();

    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(static_cast<size_t>(N) * (N - 1) / 2);
    for (int p = 0; p < N; ++p)
        for (int q = p + 1; q < N; ++q) pairs.emplace_back(p, q);
    const int M = static_cast<int>(pairs.size());

    std::vector<std::string> pair_label(M);
    for (int i = 0; i < M; ++i)
        pair_label[i] = orbital_label(pairs[i].first) + orbital_label(pairs[i].second);

    std::vector<double> eps_orb(static_cast<size_t>(N) * N, 0.0);
    for (size_t k = 0; k < lagrangian_pairs.size(); ++k) {
        const int i = lagrangian_pairs[k].orb_from;
        const int j = lagrangian_pairs[k].orb_to;
        if (i < 0 || j < 0 || i >= N || j >= N) continue;
        eps_orb[i * N + j] = lagrangians[k];
        eps_orb[j * N + i] = lagrangians[k];
    }

    const int eri_cell_w = 11;
    const int header_w = 2 + 4 + 4 + 10 + 3 + eri_cell_w * M + 1;
    const int rule_w = std::max(85, header_w);
    std::string rule(rule_w, '-');

    outfile->Printf("\n    %s Lagrangian eps_pq and ERI (pq|st) for building SSR(%d,%d) Hamiltonian matrix\n",
                    next_table_tag().c_str(), N, N);
    outfile->Printf("    %s\n", rule.c_str());
    outfile->Printf("      pq      eps_pq    st:");
    for (int j = 0; j < M; ++j) outfile->Printf(" %10s", pair_label[j].c_str());
    outfile->Printf("\n");
    outfile->Printf("    %s\n", rule.c_str());
    for (int i = 0; i < M; ++i) {
        const int p = pairs[i].first;
        const int q = pairs[i].second;
        const double eps_pq = eps_orb[p * N + q];
        outfile->Printf("      %-2s  %10.6f      ", pair_label[i].c_str(), eps_pq);
        for (int j = 0; j < M; ++j) {
            if (i == j) {
                outfile->Printf("     --    ");
            } else {
                outfile->Printf(" %10.6f", pair_eris.at(p, q, pairs[j].first, pairs[j].second));
            }
        }
        outfile->Printf("\n");
    }
    outfile->Printf("    %s\n", rule.c_str());
}

namespace {

/// True if the report level allows printing now; otherwise emits the "needs level N"
/// line and returns false. nsi is the full config-axis dimension, uncapped by report_cap.
bool si_matrix_reports(const char* title, int nsi) {
    const int min_level = nsi > kSIMatrixReportDim ? 4 : 2;
    if (reports(min_level)) return true;
    if (reports(2))
        outfile->Printf("\n    %s (%dx%d) not shown; needs REKS_REPORT_LEVEL %d\n", title,
                        nsi, nsi, min_level);
    return false;
}

/// Shared nsi x nsi table printer; elem(i, j) supplies each entry (dense array or
/// blocked storage).
template <class Elem>
void print_si_matrix_body(const studio::Cassette&   cassette,
                          const char*             title,
                          Elem                    elem,
                          int                     nsi,
                          const std::vector<int>& K_indices) {
    const int ncol = std::max(1, kSIBlockColumns);
    const int nblocks = (nsi + ncol - 1) / ncol;
    const bool use_slice = static_cast<int>(K_indices.size()) == nsi;

    std::vector<std::string> labels(nsi);
    for (int i = 0; i < nsi; ++i)
        labels[i] = si_config_name(cassette, use_slice ? K_indices[i] : i);

    // PsiOutStream::Printf flushes the stream on every call; each row is assembled
    // into one string, then Printf'd once.
    std::string line;
    line.reserve(static_cast<size_t>(ncol + 1) * 16 + 8);
    char cell[64];
    auto pad_label = [&](const std::string& lbl) {
        line.append(lbl);
        if (lbl.size() < 12) line.append(12 - lbl.size(), ' ');
    };

    outfile->Printf("\n    %s %s (%dx%d)\n", next_table_tag().c_str(), title, nsi, nsi);
    outfile->Printf("    -----------------------------------------\n");
    for (int b = 0; b < nblocks; ++b) {
        const int j0 = b * ncol;
        const int j1 = std::min(j0 + ncol, nsi);

        line.assign("      ");
        line.append(12, ' ');
        for (int j = j0; j < j1; ++j) {
            snprintf(cell, sizeof(cell), "  %12s", labels[j].c_str());
            line.append(cell);
        }
        line.push_back('\n');
        outfile->Printf(line);

        for (int i = 0; i < nsi; ++i) {
            line.assign("      ");
            pad_label(labels[i]);
            for (int j = j0; j < j1; ++j) {
                snprintf(cell, sizeof(cell), "  %+12.6f", elem(i, j));
                line.append(cell);
            }
            line.push_back('\n');
            outfile->Printf(line);
        }
        if (b + 1 < nblocks) outfile->Printf("\n");
    }
    outfile->Printf("    -----------------------------------------\n");
}

}  // namespace

void print_si_matrix(const studio::Cassette&      cassette,
                     const char*                title,
                     const std::vector<double>& M,
                     int                        nsi,
                     const std::vector<int>&    K_indices) {
    if (nsi <= 0 || M.size() < static_cast<size_t>(nsi) * nsi) return;
    if (!si_matrix_reports(title, nsi)) return;
    print_si_matrix_body(cassette, title, [&](int i, int j) {
        return M[static_cast<size_t>(i) * nsi + j];
    }, nsi, K_indices);
}

void print_si_matrix(const studio::Cassette&   cassette,
                     const char*             title,
                     const BlockedMatrix&    M,
                     int                     nsi,
                     const std::vector<int>& K_indices) {
    if (nsi <= 0 || M.n != nsi) return;
    if (!si_matrix_reports(title, nsi)) return;
    print_si_matrix_body(cassette, title, [&](int i, int j) { return M.at(i, j); },
                         nsi, K_indices);
}

std::string sector_psivar_name(int s, const std::string& name) {
    if (s == 0) return name;
    const std::string reks_prefix = "REKS ";
    const std::string tail =
        (name.rfind(reks_prefix, 0) == 0) ? name.substr(reks_prefix.size()) : name;
    return "REKS SECTOR " + std::to_string(s) + " " + tail;
}

// rank[mo] = position of mo when all MOs are sorted by (epsilon, mo_index)
// ascending, on irrep 0.
std::vector<int> mo_energy_rank(const SharedVector& epsilon, int nmo) {
    const int h = 0;
    std::vector<std::pair<double, int>> orb_e(nmo);
    for (int a = 0; a < nmo; ++a) orb_e[a] = std::make_pair(epsilon->get(h, a), a);
    std::sort(orb_e.begin(), orb_e.end());
    std::vector<int> rank(nmo);
    for (int a = 0; a < nmo; ++a) rank[orb_e[a].second] = a;
    return rank;
}

// Irrep is always "A": REKS forces C1.
std::vector<std::string> active_orbital_labels(const ReportFrame& frame) {
    const auto rank = mo_energy_rank(frame.epsilon, frame.nmo);
    const std::string irrep = frame.irrep;
    std::vector<std::string> out(frame.active.size());
    for (size_t k = 0; k < frame.active.size(); ++k)
        out[k] = std::to_string(rank[frame.active[k]] + 1) + irrep;
    return out;
}

void print_active_orbitals(const ReportFrame& frame, const std::vector<int>& rank,
                           const std::vector<std::tuple<double, int, int>>& active) {
    const int n_active = static_cast<int>(frame.active.size());
    char hdr[80];

    if (!active.empty()) {
        std::snprintf(hdr, sizeof(hdr), "Active Fractionally Occupied (%d):",
                      static_cast<int>(active.size()));
        outfile->Printf("    %-70s\n\n", hdr);

        // One column per (sector, generation, scheme) with an active geminal, labelled n, m, ... (+ scheme and
        // sector tags); NaN outside its partition.
        struct FonColumn {
            int sector;
            int gen;
            int scheme;
            std::vector<double> fon;
        };
        std::vector<FonColumn> cols;
        for (int sec = 0; sec < frame.sa_cassette.n_sectors(); ++sec) {
            for (int gen = 0; gen < frame.sa_cassette.n_generations(); ++gen) {
                const auto& act = frame.sa_cassette.geminals_active(sec, gen);
                if (act.empty()) continue;
                std::vector<int> gen_schemes;
                for (int g : act) {
                    const int s = frame.sa_cassette.geminal_templates()[g].scheme;
                    if (std::find(gen_schemes.begin(), gen_schemes.end(), s) == gen_schemes.end())
                        gen_schemes.push_back(s);
                }
                std::sort(gen_schemes.begin(), gen_schemes.end());
                for (int s : gen_schemes) {
                    FonColumn col;
                    col.sector = sec;
                    col.gen = gen;
                    col.scheme = s;
                    col.fon.assign(n_active, std::numeric_limits<double>::quiet_NaN());
                    for (int g : act) {
                        const auto& tpl = frame.sa_cassette.geminal_templates()[g];
                        if (tpl.scheme != s) continue;
                        col.fon[tpl.orbitals[0]] = frame.reel.fon_state[sec].layers[gen][g].p;
                        col.fon[tpl.orbitals[1]] = frame.reel.fon_state[sec].layers[gen][g].q;
                    }
                    cols.push_back(std::move(col));
                }
            }
        }

        std::vector<int> all_schemes;
        for (const auto& c : cols)
            if (std::find(all_schemes.begin(), all_schemes.end(), c.scheme) == all_schemes.end())
                all_schemes.push_back(c.scheme);
        std::sort(all_schemes.begin(), all_schemes.end());
        const bool multi_scheme = all_schemes.size() > 1;
        const int ncol = static_cast<int>(cols.size());

        // Slot-letter partition "(a,f) & ..." for scheme s and its energy-rank
        // twin, pooled over all generations and de-duplicated. Returns pair count.
        auto partition_strings = [&](int s, std::string& slots, std::string& ranks) {
            std::vector<std::pair<int, int>> pairs;
            for (int gen = 0; gen < frame.sa_cassette.n_generations(); ++gen)
                for (int g : frame.sa_cassette.geminals_active_gen(gen)) {
                    const auto& tpl = frame.sa_cassette.geminal_templates()[g];
                    if (tpl.scheme != s) continue;
                    std::pair<int, int> pr{tpl.orbitals[0], tpl.orbitals[1]};
                    if (std::find(pairs.begin(), pairs.end(), pr) == pairs.end())
                        pairs.push_back(pr);
                }
            std::sort(pairs.begin(), pairs.end());
            slots.clear();
            ranks.clear();
            for (size_t i = 0; i < pairs.size(); ++i) {
                const int p = pairs[i].first;
                const int q = pairs[i].second;
                if (i) { slots += " & "; ranks += " & "; }
                slots += "(" + reks::report::orbital_label(p) + ","
                             + reks::report::orbital_label(q) + ")";
                const std::string r_p =
                    std::to_string(rank[frame.active[p]] + 1) + frame.irrep;
                const std::string r_q =
                    std::to_string(rank[frame.active[q]] + 1) + frame.irrep;
                ranks += "(" + r_p + "," + r_q + ")";
            }
            return pairs.size();
        };

        if (ncol > 1) {
            // Header: one shared scheme line when a single pairing scheme spans
            // all generations, else a per-scheme legend.
            if (!multi_scheme) {
                std::string slots, ranks;
                if (partition_strings(all_schemes[0], slots, ranks) > 1)
                    outfile->Printf("      Geminal Pairs Scheme: %s  =  %s\n\n",
                                    slots.c_str(), ranks.c_str());
            } else {
                for (size_t i = 0; i < all_schemes.size(); ++i) {
                    std::string slots, ranks;
                    partition_strings(all_schemes[i], slots, ranks);
                    outfile->Printf("      %-14s Scheme %d: %s\n",
                                    i == 0 ? "Geminal Pairs" : "", all_schemes[i], slots.c_str());
                }
                outfile->Printf("\n");
            }

            // Cell body "<label> = <fon>"; label carries a scheme index only when
            // more than one pairing scheme is present.
            auto make_body = [&](const FonColumn& c, const std::string& orb, double v) {
                std::string lab = reks::studio::fon_channel_name(
                    c.gen, c.sector, multi_scheme ? c.scheme : -1);
                lab += "_" + orb;
                char buf[64];
                std::snprintf(buf, sizeof(buf), "%s = %9.6f", lab.c_str(), v);
                return std::string(buf);
            };

            // Uniform value-cell width so the columns line up.
            size_t body_w = 0;
            for (const auto& it : active) {
                const int apos = std::get<2>(it);
                const std::string orb = reks::report::orbital_label(apos);
                for (const auto& c : cols)
                    if (!std::isnan(c.fon[apos]))
                        body_w = std::max(body_w, make_body(c, orb, c.fon[apos]).size());
            }

            constexpr int col_gap = 4;
            for (const auto& it : active) {
                const double e = std::get<0>(it);
                const int rk = std::get<1>(it);
                const int apos = std::get<2>(it);
                const std::string orb = reks::report::orbital_label(apos);
                char prefix[64];
                std::snprintf(prefix, sizeof(prefix), "    %4d%-4s(%s) %9.6f",
                              rk, frame.irrep.c_str(), orb.c_str(), e);
                std::string row = prefix;
                for (const auto& c : cols) {
                    row += std::string(col_gap, ' ');
                    std::string body = std::isnan(c.fon[apos])
                                           ? std::string()
                                           : make_body(c, orb, c.fon[apos]);
                    body.resize(body_w, ' ');
                    row += body;
                }
                while (!row.empty() && row.back() == ' ') row.pop_back();
                outfile->Printf("%s\n", row.c_str());
            }
        } else {
            std::vector<double> fon_of_active_pos(
                n_active, std::numeric_limits<double>::quiet_NaN());
            for (int g : frame.sa_cassette.geminals_active_gen(0)) {
                const auto& orbs = frame.sa_cassette.geminal_templates()[g].orbitals;
                fon_of_active_pos[orbs[0]] = frame.reel.fon_state[frame.sa_cassette.sector()].layers[0][g].p;
                fon_of_active_pos[orbs[1]] = frame.reel.fon_state[frame.sa_cassette.sector()].layers[0][g].q;
            }
            const int n_pairs = static_cast<int>(frame.sa_cassette.geminals_active_gen(0).size());
            if (n_pairs > 1) {
                std::vector<std::pair<int, int>> pairs;
                pairs.reserve(n_pairs);
                for (int g : frame.sa_cassette.geminals_active_gen(0)) {
                    const auto& orbs = frame.sa_cassette.geminal_templates()[g].orbitals;
                    pairs.emplace_back(orbs[0], orbs[1]);
                }
                std::sort(pairs.begin(), pairs.end());

                std::string scheme;
                std::string scheme_rank;
                for (size_t i = 0; i < pairs.size(); ++i) {
                    const int p = pairs[i].first;
                    const int q = pairs[i].second;
                    const std::string lbl_p = reks::report::orbital_label(p);
                    const std::string lbl_q = reks::report::orbital_label(q);
                    if (i > 0) scheme += " & ";
                    scheme += "(" + lbl_p + "," + lbl_q + ")";
                    const std::string r_p = std::to_string(rank[frame.active[p]] + 1) + frame.irrep;
                    const std::string r_q = std::to_string(rank[frame.active[q]] + 1) + frame.irrep;
                    if (i > 0) scheme_rank += " & ";
                    scheme_rank += "(" + r_p + "," + r_q + ")";
                }
                outfile->Printf("      Geminal Pairs Scheme: %s  =  %s\n\n",
                                scheme.c_str(), scheme_rank.c_str());
            }

            for (const auto& it : active) {
                const double e = std::get<0>(it);
                const int rk = std::get<1>(it);
                const int apos = std::get<2>(it);
                const std::string lbl = reks::report::orbital_label(apos);
                outfile->Printf("    %4d%-4s(%s) %9.6f    n_%s = %9.6f\n",
                                rk, frame.irrep.c_str(), lbl.c_str(), e, lbl.c_str(),
                                fon_of_active_pos[apos]);
            }
        }
        outfile->Printf("\n");
    }
}

void print_fon_table(const ReportFrame& frame, const studio::Cassette& cassette) {
    if (!reports(2)) return;
    const int n_active = static_cast<int>(frame.active.size());
    if (n_active <= 0 || !frame.epsilon) return;

    const int n_gen = cassette.n_generations();
    const auto* tpl_tab = cassette.geminal_templates();

    auto add_unique = [](std::vector<int>& v, int x) {
        if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
    };

    // Distinct scheme tags active in any generation of the cassette.
    std::vector<int> scheme_tags;
    for (int gen = 0; gen < n_gen; ++gen)
        for (int g : cassette.geminals_active_gen(gen))
            add_unique(scheme_tags, tpl_tab[g].scheme);
    std::sort(scheme_tags.begin(), scheme_tags.end());

    auto gen_has_tag = [&](int gen, int tag) {
        for (int g : cassette.geminals_active_gen(gen))
            if (tpl_tab[g].scheme == tag) return true;
        return false;
    };
    auto orb_letter = [&](int p) { return reks::report::orbital_label(p); };

    // Pair partition is generation-independent; collect the scheme's distinct orbital pairs.
    auto pair_partition_label = [&](int tag) {
        std::string out;
        std::vector<std::pair<int, int>> seen;
        for (int gen = 0; gen < n_gen; ++gen) {
            for (int g : cassette.geminals_active_gen(gen)) {
                const auto& gem = tpl_tab[g];
                if (gem.scheme != tag) continue;
                std::pair<int, int> pr{gem.orbitals[0], gem.orbitals[1]};
                if (std::find(seen.begin(), seen.end(), pr) != seen.end()) continue;
                seen.push_back(pr);
                if (!out.empty()) out += " & ";
                out += "(" + orb_letter(pr.first) + "," + orb_letter(pr.second) + ")";
            }
        }
        return out;
    };

    // Per-orbital FON for one (generation, scheme tag); filled[k] flags orbitals used by the scheme.
    auto build_fon = [&](int gen, int tag) {
        std::vector<double> fon(n_active, 0.0);
        std::vector<int> filled(n_active, 0);
        const auto& cass_fons = frame.reel.fon_state[cassette.sector()];
        for (int g : cassette.geminals_active_gen(gen)) {
            const auto& tmpl = tpl_tab[g];
            if (tmpl.scheme != tag) continue;
            const int p = tmpl.orbitals[0];
            const int q = tmpl.orbitals[1];
            if (p < 0 || p >= n_active || q < 0 || q >= n_active) continue;
            const auto& fon_g = cass_fons.layers[gen][g];
            fon[p] = fon_g.p;
            fon[q] = fon_g.q;
            filled[p] = 1;
            filled[q] = 1;
        }
        return std::make_pair(std::move(fon), std::move(filled));
    };

    const int label_w = 22;
    const int cell_w = 12;
    const int rule_w = label_w + cell_w * n_active + 4;
    const std::string rule(rule_w, '-');

    outfile->Printf("\n    %s Fractional Occupation Numbers (active orbitals)\n",
                    reks::report::next_table_tag().c_str());
    outfile->Printf("    %s\n", rule.c_str());

    for (size_t i = 0; i < scheme_tags.size(); ++i) {
        outfile->Printf("      %-14s Scheme %d: %s\n",
                        i == 0 ? "Geminal Pairs" : "",
                        scheme_tags[i],
                        pair_partition_label(scheme_tags[i]).c_str());
    }
    outfile->Printf("    %s\n", rule.c_str());

    // Energy-rank labels in slot order.
    const std::vector<std::string> rk_lab = active_orbital_labels(frame);
    outfile->Printf("      %-14s ", "Orbital");
    for (int k = 0; k < n_active; ++k)
        outfile->Printf(" %10s ", rk_lab[k].c_str());
    outfile->Printf("\n");

    // frame.epsilon irrep 0: REKS is C1.
    outfile->Printf("      %-14s ", "Energy (Ha)");
    for (int k = 0; k < n_active; ++k) {
        outfile->Printf(" %10.6f ", frame.epsilon->get(0, frame.active[k]));
    }
    outfile->Printf("\n");

    outfile->Printf("      %-14s ", "REKS label");
    for (int k = 0; k < n_active; ++k) {
        const std::string lab = reks::report::orbital_label(k);
        outfile->Printf(" %10s ", lab.c_str());
    }
    outfile->Printf("\n");

    for (int tag : scheme_tags) {
        for (int gen = 0; gen < n_gen; ++gen) {
            if (!gen_has_tag(gen, tag)) continue;
            auto [fon, filled] = build_fon(gen, tag);
            const std::string chan = reks::studio::fon_channel_name(gen, cassette.sector());
            char header[40];
            std::snprintf(header, sizeof(header), "%s-FON (Scheme%d)", chan.c_str(), tag);
            outfile->Printf("      %-13s", header);
            for (int k = 0; k < n_active; ++k) {
                if (filled[k]) outfile->Printf(" %10.6f ", fon[k]);
                else           outfile->Printf(" %10s ", "--");
            }
            outfile->Printf("\n");
        }
    }

    outfile->Printf("    %s\n", rule.c_str());
}

void print_memory_footprint(const std::string& tag, const ReportFrame& frame,
                            const Wavefunction& wfn) {
    auto mat_elems = [](const SharedMatrix& M) -> std::size_t {
        if (!M) return 0;
        std::size_t n = 0;
        for (int h = 0; h < M->nirrep(); ++h) {
            n += static_cast<std::size_t>(M->rowspi()[h]) *
                 static_cast<std::size_t>(M->colspi()[h]);
        }
        return n;
    };
    auto vec_elems = [&](const std::vector<SharedMatrix>& v) -> std::size_t {
        std::size_t n = 0;
        for (const auto& M : v) n += mat_elems(M);
        return n;
    };
    auto elems_to_mib = [](std::size_t n) -> double {
        return static_cast<double>(n) * 8.0 / 1048576.0;
    };

    const int nso = wfn.nsopi()[0];
    const int n_sa = static_cast<int>(frame.sa_cassette.microstates().size());
    const int n_micro = frame.sa_cassette.n_total_microstates();  // per-L reel arrays are this long
    const int n_active = frame.sa_cassette.n_active_orbitals();
    const int n_base = (1 << frame.sa_cassette.n_active_orbitals());

    auto diag_elems = [](const std::vector<std::vector<double>>& v) -> std::size_t {
        std::size_t n = 0;
        for (const auto& d : v) n += d.size();
        return n;
    };
    const std::size_t e_F_AO_a = mat_elems(frame.reel.F_alpha_AO_buffer);
    const std::size_t e_F_MO_a = vec_elems(frame.reel.F_MO_arows_a_L)
                               + diag_elems(frame.reel.F_MO_diag_a_L);
    const std::size_t e_F_MO_b = vec_elems(frame.reel.F_MO_arows_b_L)
                               + diag_elems(frame.reel.F_MO_diag_b_L);
    const std::size_t e_F_acc  = mat_elems(frame.reel.F_acc_AO) + mat_elems(frame.reel.F_acc_MO);
    const std::size_t e_Vxc_a  = vec_elems(frame.reel.Vxc_arows_a_L) + diag_elems(frame.reel.Vxc_diag_a_L)
                               + mat_elems(frame.reel.V_acc_xc_a);
    const std::size_t e_Vxc_b  = vec_elems(frame.reel.Vxc_arows_b_L) + diag_elems(frame.reel.Vxc_diag_b_L)
                               + mat_elems(frame.reel.V_acc_xc_b);
    const std::size_t e_Dcore  = mat_elems(frame.reel.D_core);
    const std::size_t e_J_act  = vec_elems(frame.reel.J_act);
    const std::size_t e_K_act  = vec_elems(frame.reel.K_act);
    const std::size_t e_wK_act = vec_elems(frame.reel.wK_act);
    const std::size_t e_J_core = mat_elems(frame.reel.J_core);
    const std::size_t e_K_core = mat_elems(frame.reel.K_core);
    const std::size_t e_wK_core= mat_elems(frame.reel.wK_core);
    const std::size_t e_Freks  = mat_elems(frame.reel.F_reks);
    const std::size_t e_Twork  = mat_elems(frame.reel.temp_buffer);
    const std::size_t e_Ca     = mat_elems(wfn.Ca());
    const std::size_t e_Cb     = mat_elems(wfn.Cb());
    const std::size_t e_Da     = mat_elems(wfn.Da());
    const std::size_t e_Db     = mat_elems(wfn.Db());
    const std::size_t e_H      = mat_elems(wfn.H());
    const std::size_t e_S      = mat_elems(wfn.S());

    const std::size_t e_total =
        e_F_AO_a + e_F_MO_a + e_F_MO_b + e_F_acc +
        e_Vxc_a  + e_Vxc_b  +
        e_Dcore  +
        e_J_act  + e_K_act  + e_wK_act +
        e_J_core + e_K_core + e_wK_core +
        e_Freks  + e_Twork  +
        e_Ca + e_Cb + e_Da + e_Db + e_H + e_S;

    outfile->Printf("\n  ==> REKS Memory Footprint (%s) <==\n\n", tag.c_str());
    outfile->Printf("    Dimensions: nso=%d  Ncore=%d  n_active=%d  n_sa=%d  n_micro=%d  n_base=%d\n",
                    nso, frame.Ncore, n_active, n_sa, n_micro, n_base);
    outfile->Printf("    %-22s  %6s  %20s  %12s\n", "Pool", "Count", "Dim", "MiB");
    outfile->Printf("    ---------------------------------------------------------------------\n");

    auto row_vec = [&](const char* name, const std::vector<SharedMatrix>& v, std::size_t e) {
        const int cnt = static_cast<int>(v.size());
        int r = 0, c = 0;
        for (const auto& M : v) { if (M) { r = M->rowspi()[0]; c = M->colspi()[0]; break; } }
        char dim[24];
        std::snprintf(dim, sizeof(dim), "%d x %d", r, c);
        outfile->Printf("    %-22s  %6d  %20s  %12.2f\n", name, cnt, dim, elems_to_mib(e));
    };
    auto row_mat = [&](const char* name, const SharedMatrix& M, std::size_t e) {
        int r = M ? M->rowspi()[0] : 0;
        int c = M ? M->colspi()[0] : 0;
        char dim[24];
        std::snprintf(dim, sizeof(dim), "%d x %d", r, c);
        outfile->Printf("    %-22s  %6d  %20s  %12.2f\n", name, M ? 1 : 0, dim, elems_to_mib(e));
    };

    row_mat("reel_.F_alpha_AO_buffer",  frame.reel.F_alpha_AO_buffer, e_F_AO_a);
    row_vec("F_MO_arows_a (MO)",     frame.reel.F_MO_arows_a_L, e_F_MO_a);
    row_vec("F_MO_arows_b (MO)",     frame.reel.F_MO_arows_b_L, e_F_MO_b);
    row_mat("reel_.F_acc_MO",        frame.reel.F_acc_MO,       e_F_acc);
    row_vec("Vxc_arows_a (MO)",   frame.reel.Vxc_arows_a_L, e_Vxc_a);
    row_vec("Vxc_arows_b (MO)",   frame.reel.Vxc_arows_b_L, e_Vxc_b);
    row_mat("reel_.D_core",            frame.reel.D_core, e_Dcore);
    row_vec("reel_.J_act",              frame.reel.J_act,  e_J_act);
    row_vec("reel_.K_act",              frame.reel.K_act,  e_K_act);
    row_vec("reel_.wK_act",             frame.reel.wK_act, e_wK_act);
    row_mat("reel_.J_core",             frame.reel.J_core,  e_J_core);
    row_mat("reel_.K_core",             frame.reel.K_core,  e_K_core);
    row_mat("reel_.wK_core",            frame.reel.wK_core, e_wK_core);
    row_mat("reel_.F_reks",             frame.reel.F_reks,  e_Freks);
    row_mat("reel_.temp_buffer",          frame.reel.temp_buffer, e_Twork);
    row_mat("Ca_",                 wfn.Ca(), e_Ca);
    row_mat("Cb_",                 wfn.Cb(), e_Cb);
    row_mat("Da_",                 wfn.Da(), e_Da);
    row_mat("Db_",                 wfn.Db(), e_Db);
    row_mat("H_",                  wfn.H(),  e_H);
    row_mat("S_",                  wfn.S(),  e_S);
    outfile->Printf("    ---------------------------------------------------------------------\n");
    outfile->Printf("    %-22s  %6s  %20s  %12.2f\n", "TOTAL", "", "", elems_to_mib(e_total));
    outfile->Printf("\n");
}

void print_mo_consistency(const ReportFrame& frame) {
    if (!frame.reel.Ca_fock_transform) return;

    const int nao = frame.Ca->nrow();
    const int nmo = frame.Ca->ncol();
    double** Cp = frame.Ca->pointer(0);
    double** Cp_fock = frame.reel.Ca_fock_transform->pointer(0);

    int act_min = frame.active[0], act_max = frame.active[0];
    for (int idx : frame.active) {
        act_min = std::min(act_min, idx);
        act_max = std::max(act_max, idx);
    }
    const int mo_lo = std::max(0, act_min - 4);
    const int mo_hi = std::min(nmo - 1, act_max + 4);
    const int nshow = mo_hi - mo_lo + 1;

    outfile->Printf("\n  [MO-CONSISTENCY] Ca_ (post form_C) vs reel_.Ca_fock_transform (used for F_MO / SI):\n");

    outfile->Printf("%21s", "");
    for (int j = 0; j < nshow; ++j) {
        int mo_j = mo_lo + j;
        const char* tag = (mo_j < act_min) ? "core" : (mo_j > act_max ? "virt" : "act ");
        outfile->Printf("    %s[%4d]      ", tag, mo_j);
    }
    outfile->Printf("\n");

    outfile->Printf("%21s", "");
    for (int j = 0; j < nshow; ++j) outfile->Printf("        Ca     Cfock");
    outfile->Printf("\n");

    for (int mu = 0; mu < nao; ++mu) {
        outfile->Printf("           AO  [%4d]", mu);
        for (int j = 0; j < nshow; ++j) {
            int mo_j = mo_lo + j;
            outfile->Printf("   %+7.4f   %+7.4f", Cp[mu][mo_j], Cp_fock[mu][mo_j]);
        }
        outfile->Printf("\n");
    }
}

void print_f_reks_mo_debug(const ReportFrame& frame) {
    if (!frame.reel.F_reks_MO) return;

    const int nmo = frame.reel.F_reks_MO->ncol();
    double** Fp = frame.reel.F_reks_MO->pointer(0);

    int act_min = frame.active[0], act_max = frame.active[0];
    for (int idx : frame.active) {
        act_min = std::min(act_min, idx);
        act_max = std::max(act_max, idx);
    }
    const int mo_lo = std::max(0, act_min - 4);
    const int mo_hi = std::min(nmo - 1, act_max + 4);
    const int nshow = mo_hi - mo_lo + 1;

    outfile->Printf("\n  [F_REKS_MO] SA coupling Fock in MO basis\n");

    outfile->Printf("%21s", "");
    for (int j = 0; j < nshow; ++j) {
        int mo_j = mo_lo + j;
        const char* tag = (mo_j < act_min) ? "core" : (mo_j > act_max ? "virt" : "act ");
        outfile->Printf("  %s[%4d]", tag, mo_j);
    }
    outfile->Printf("\n");

    for (int i = 0; i < nshow; ++i) {
        int mo_i = mo_lo + i;
        const char* tag = (mo_i < act_min) ? "core" : (mo_i > act_max ? "virt" : "act ");
        outfile->Printf("           %s[%4d]", tag, mo_i);
        for (int j = 0; j < nshow; ++j) {
            int mo_j = mo_lo + j;
            outfile->Printf("   %+9.5f", Fp[mo_i][mo_j]);
        }
        outfile->Printf("\n");
    }
}

void print_si_hamiltonian_overlap(const ReportFrame& frame, const studio::Cassette& cassette,
                                  const SIResult& sir) {
    if (sir.hamiltonian.empty()) return;
    const int nsi = sir.n_states;
    if (nsi <= 0) return;

    const std::vector<int>& K_indices = cassette.K_indices;

    reks::report::print_si_matrix(frame.sa_cassette, "Pre-diagonalization SI Hamiltonian",
                                    sir.hamiltonian, nsi, K_indices);
    if (!sir.overlap.empty()) {
        reks::report::print_si_matrix(frame.sa_cassette, "SI Overlap S",
                                        sir.overlap, nsi, K_indices);
    }
}

void print_state_properties_header(const ReportFrame& frame, const studio::Cassette& cassette,
                                   int K, int spin2) {
    if (!reports(2)) return;
    std::string ml = reks::report::method_label(
        frame.sa_cassette,
        static_cast<int>(cassette.K_indices.size()));
    // Multi-sector runs tag the title with the SI section's spin manifold.
    if (frame.sa_cassette.n_sectors() > 1) {
        ml += " " + reks::report::spin_tag(spin2);
    }
    outfile->Printf("\n  ==> %s State Properties (K=%d) <==\n\n", ml.c_str(), K);
}

void print_sa_diagnostics(const ReportFrame& frame, int iteration, ConvergenceTrace& trace) {
    int n_sa = static_cast<int>(frame.sa_cassette.microstates().size());
    if (n_sa < 2) return;
    if (!frame.catalog.data || frame.reel.E_L.empty()) return;

    double E_SA = reks::scf::compute_E_SA(frame.sa_cassette, frame.reel.E_L, frame.reel.C_L);

    // trace.last_computed_dE must be set before trace.prev_E_total is overwritten below.
    trace.last_computed_dE =
        (iteration > 1 && std::abs(trace.prev_E_total) > 1e-10) ? (E_SA - trace.prev_E_total) : 0.0;
    if (iteration > 1 && std::abs(trace.prev_E_total) > 1e-10) {
        if (iteration > 2 && std::abs(trace.prev_dE) > 1e-10 && trace.last_computed_dE * trace.prev_dE < 0) {
            trace.energy_oscillation_count++;
            outfile->Printf("  *** E_SA OSCILLATION (total: %d) ***\n", trace.energy_oscillation_count);
        }
        trace.prev_dE = trace.last_computed_dE;
    }
    trace.prev_E_total = E_SA;
}

void print_iteration_diagnostics(const ReportFrame& frame, int iteration, ConvergenceTrace& trace) {
    double E_SA = reks::scf::compute_E_SA(frame.sa_cassette, frame.reel.E_L, frame.reel.C_L);
    double dE = trace.last_computed_dE;

    outfile->Printf("  [REKS] iter=%d E_SA=%.14f dE=%+.6e", iteration, E_SA, dE);

    // Gen-0 FONs per sector and their change since the last iteration (trace.prev_fon, sector order).
    const int n_sec = frame.sa_cassette.n_sectors();
    std::vector<double> cur_fon;
    int flat = 0;
    for (int sec = 0; sec < n_sec; ++sec) {
        const auto& gem0 = frame.sa_cassette.geminals_active(sec, 0);
        const int nps = static_cast<int>(gem0.size());
        if (nps == 0) continue;
        if (n_sec > 1) outfile->Printf(" sector%d", sec);
        for (int p = 0; p < nps; ++p) {
            double fp = frame.reel.fon_state[sec].layers[0][gem0[p]].p;
            double prev_fp = (flat < static_cast<int>(trace.prev_fon.size())) ? trace.prev_fon[flat] : fp;
            outfile->Printf(" f%d=%.8f df%d=%+.4e", p, fp, p, fp - prev_fp);
            cur_fon.push_back(fp);
            ++flat;
        }
    }

    // Fock eigenvalue gap eps[q] - eps[p] over the shared active spectrum; the gen-0
    // pair union is sector-agnostic (one orbital set).
    {
        double* eps = frame.epsilon->pointer(0);
        const auto& gen0_pool = frame.sa_cassette.geminals_active_gen(0);
        const int npair = static_cast<int>(gen0_pool.size());
        for (int p = 0; p < npair; ++p) {
            const auto& orbs = frame.sa_cassette.geminal_templates()[gen0_pool[p]].orbitals;
            int mo_p = frame.active[orbs[0]];
            int mo_q = frame.active[orbs[1]];
            outfile->Printf(" fgap%d=%+.6f", p, eps[mo_q] - eps[mo_p]);
        }
    }

    if (trace.energy_oscillation_count > 0) outfile->Printf(" osc=%d", trace.energy_oscillation_count);
    outfile->Printf("\n");

    trace.prev_fon = std::move(cur_fon);
}

}  // namespace report
}  // namespace reks
}  // namespace psi
