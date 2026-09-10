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

#include "reks_report.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <algorithm>
#include <cstdio>
#include <numeric>
#include <string>
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
"  config block, a config type ('PPS', 'OSS', ...) selects every config of that\n"
"  type, e.g. si_reks_configs [['full']] or [['PPS', 'OSS']].\n"
"  'single-scheme:N' selects the configs of pairing scheme N and 'mixed-scheme:N'\n"
"  the SSR state set on scheme N (open-shell singlets from every scheme); the\n"
"  spin-projected configs are the base and ride along with either. Both default to 0\n"
"  and take several indices ('single-scheme-1-2-3' unions the three). Quote them, and\n"
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
    const char* name = cassette.si_config_name(K);
    if (name != nullptr) return name;
    return "cfg" + std::to_string(K);
}

std::string method_label(const studio::Cassette& sa_cassette, int n_si) {
    const int n_sa = static_cast<int>(sa_cassette.K_indices.size())
                   + sa_cassette.n_extra_microstates();
    // Assumes REKS(N,N): the (M,M) label uses n_active_orbitals as electron count.
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
    // 1-based MO index + irrep ("A": REKS forces C1).
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
                                   double e_sa, double e_pen, double e_vv10) {
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

    if (e_pen != 0.0 || e_vv10 != 0.0) {
        if (e_pen != 0.0)
            outfile->Printf("    %-*s%20.12f\n", prefix_w, "IPR Penalty Energy", e_pen);
        if (e_vv10 != 0.0)
            outfile->Printf("    %-*s%20.12f\n", prefix_w, "VV10 Nonlocal Energy", e_vv10);
        outfile->Printf("    %s\n", rule.c_str());
        outfile->Printf("    %-*s%20.12f\n", prefix_w, "Total Energy",
                        e_sa + e_pen + e_vv10);
    }
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
        // Display-only report-cap note; dim above remains the true cassette dimension.
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
    if (N < 2) return;

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
    if (lagrangian_pairs.size() != lagrangians.size())
        throw PSIEXCEPTION("reks::report: lagrangian_pairs / lagrangians size mismatch");
    for (size_t k = 0; k < lagrangian_pairs.size(); ++k) {
        const int i = lagrangian_pairs[k].orb_from;
        const int j = lagrangian_pairs[k].orb_to;
        if (i < 0 || j < 0 || i >= N || j >= N) continue;
        eps_orb[i * N + j] = lagrangians[k];
        eps_orb[j * N + i] = lagrangians[k];
    }

    const bool have_eris = !pair_eris.empty() && pair_eris.n_active() == N;

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
            } else if (have_eris) {
                const int s = pairs[j].first;
                const int t = pairs[j].second;
                const double v = pair_eris.at(p, q, s, t);
                outfile->Printf(" %10.6f", v);
            } else {
                outfile->Printf("        - ");
            }
        }
        outfile->Printf("\n");
    }
    outfile->Printf("    %s\n", rule.c_str());
    if (!have_eris) {
        outfile->Printf("    (ERI (pq|st) not provided; showing eps_pq column only)\n");
    }
}

namespace {

/// True when the caller may print now; otherwise prints the "needs level N" line and
/// returns false. nsi is the full config-axis dimension, uncapped by report_cap.
bool si_matrix_reports(const char* title, int nsi) {
    const int min_level = nsi > kSIMatrixReportDim ? 4 : 2;
    if (reports(min_level)) return true;
    if (reports(2))
        outfile->Printf("\n    %s (%dx%d) not shown; needs REKS_REPORT_LEVEL %d\n", title,
                        nsi, nsi, min_level);
    return false;
}

/// Body of print_si_matrix over any element accessor elem(i, j).
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

}  // namespace report
}  // namespace reks
}  // namespace psi
