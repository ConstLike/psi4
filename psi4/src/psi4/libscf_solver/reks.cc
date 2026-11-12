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

#include "reks.h"
#include "psi4/libreks/reks_arch.h"
#include "psi4/libreks/reks_diis_step.h"
#include "psi4/libreks/reks_objectives.h"
#include "psi4/libreks/reks_si.h"
#include "psi4/libreks/reks_rdm.h"
#include "psi4/libreks/reks_report.h"
#include "psi4/libreks/reks_properties.h"
#include "psi4/libreks/reks_fon_relax.h"
#include "psi4/libreks/reks_fon_derivs.h"
#include "psi4/libreks/reks_fon_step.h"
#include "psi4/libreks/reks_gvb_driver.h"
#include "psi4/libreks/reks_math.h"
#include "psi4/libreks/reks_scf.h"
#include "psi4/libreks/reks_trah.h"
#include "psi4/libreks/reks_xc.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/factory.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <array>
#include <cctype>
#include <set>
#include <tuple>
#include <cstdio>

namespace psi {
namespace scf {

using reks::report::reports;

namespace {

// Returns true iff user assigned a value to either form of an alias pair.
bool alias_changed(Options& opts,
                   const std::string& primary,
                   const std::string& alias) {
    return opts[primary].has_changed() || opts[alias].has_changed();
}

// Data entry of whichever alias the user set, throws if both are set; without
// alias_changed() it is the primary's default.
Data& read_alias(Options& opts,
                 const std::string& primary,
                 const std::string& alias) {
    const bool p = opts[primary].has_changed();
    const bool a = opts[alias].has_changed();
    if (p && a) {
        reks::report::input_error("Specify only one of '" + primary + "' and '" + alias +
                                  "' (both were set).");
    }
    return a ? opts[alias] : opts[primary];
}

// One config-pool entry -> its global K indices: a bulk token ("full", a config type, "sps", a scheme
// token, a '*'/'?' pattern) names many, a config name one.
std::vector<int> resolve_config_entry(const reks::studio::Catalog& cat,
                                      const std::string& tok) {
    std::vector<int> bulk = reks::studio::expand_config_token(cat, tok);
    if (!bulk.empty()) return bulk;
    const int K = reks::studio::find_config_index(cat, tok);
    if (K < 0)
        reks::report::input_error(
            "Unknown REKS config token '" + tok +
            "'. Use the keyword 'full', 'sps', 'single-scheme[:N]', 'mixed-scheme[:N]', "
            "a name pattern with '*' or '?', a config type (" +
            reks::studio::config_types_csv(cat) +
            "), or a config name (" +
            reks::studio::config_names_csv(cat) + ").");
    return {K};
}

// "exclude:X" (also 'exclude-X' / 'exclude_X', since the PSithon array parser
// splits an element on ':') -> X, else empty.
std::string exclude_payload(const std::string& tok) {
    const std::string key = "exclude";
    if (tok.size() <= key.size() + 1) return {};
    for (size_t i = 0; i < key.size(); ++i)
        if (std::tolower(static_cast<unsigned char>(tok[i])) != key[i]) return {};
    const char sep = tok[key.size()];
    if (sep != ':' && sep != '-' && sep != '_') return {};
    return tok.substr(key.size() + 1);
}

// Config-pool array (SA list or SI row) -> global K: strings are entries or "exclude:<entry>" (subtracted,
// order-free), numbers are sector-local indices (+ config_start). Every kept K is range-checked against the
// sector; no dup check.
std::vector<int> parse_config_indices(Data& arr, const reks::studio::Catalog& cat,
                                      bool allow_exclude, const std::string& optname,
                                      const std::string& where) {
    const int config_start = cat.sector_config_start();
    const int config_count = cat.sector_config_count();
    std::vector<int> out;
    std::set<int> dropped;
    const int n = static_cast<int>(arr.size());
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        const std::string t = arr[i].type();
        if (t != "string" && t != "istring") {
            out.push_back(config_start + static_cast<int>(arr[i].to_double()));
            continue;
        }
        const std::string tok = arr[i].to_string();
        const std::string drop = exclude_payload(tok);
        if (drop.empty()) {
            const std::vector<int> add = resolve_config_entry(cat, tok);
            out.insert(out.end(), add.begin(), add.end());
            continue;
        }
        if (!allow_exclude)
            reks::report::input_error(
                "SA_REKS_CONFIGS: '" + tok +
                "' is not allowed in an SA block, whose weights follow the typed "
                "order; list the configs the block keeps instead.");
        const std::vector<int> gone = resolve_config_entry(cat, drop);
        dropped.insert(gone.begin(), gone.end());
    }
    if (!dropped.empty()) {
        std::vector<int> kept;
        kept.reserve(out.size());
        for (int K : out)
            if (dropped.count(K) == 0) kept.push_back(K);
        out.swap(kept);
    }
    for (int K : out)
        if (K < config_start || K >= config_start + config_count)
            reks::report::input_error(
                optname + ": config index (local) " + std::to_string(K - config_start) +
                " is out of range [0, " + std::to_string(config_count) + ") for " + where + ".");
    return out;
}

// Nesting depth of one option entry: scalar = 0, list of scalars (or empty
// list) = 1, list of lists = 2, and so on.
int entry_depth(Data& d) {
    if (!d.is_array()) return 0;
    int mx = 0;
    for (int i = 0; i < static_cast<int>(d.size()); ++i) mx = std::max(mx, entry_depth(d[i]));
    return 1 + mx;
}

}  // anonymous namespace

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    reks_common_init();
}

REKS::~REKS() {
    if (rv_potential_) {
        rv_potential_->finalize();
    }
}

void REKS::setup_potential() {
    if (functional_->needs_xc()) {
        rv_potential_ = std::make_shared<psi::RV>(functional_, basisset_, options_);
        rv_potential_->initialize();

        ms_xc_ = std::make_unique<reks::xc::MicrostateXC>(rv_potential_);

        if (reports(4)) {
            outfile->Printf("\n  REKS: Initialized microstate XC on the RV grid, functional %s\n",
                            functional_->name().c_str());
            outfile->Printf("  REKS: Initialized RV potential for closed-shell XC on Da_\n");
        }
    }
}

void REKS::save_density_and_energy() {}

double REKS::compute_initial_E() { return nuclearrep_ + Da_->vector_dot(H_); }

int REKS::complete_orbital_basis(SharedMatrix C, SharedMatrix S, SharedMatrix X,
                                  int n_loaded, int h) {
    const int nso = C->rowspi()[h];
    const int nmo = C->colspi()[h];

    double** Cp = C->pointer(h);
    double** Sp = S->pointer(h);
    double** Xp = X->pointer(h);

    // S-metric MGS: c <- c - sum_k (c_k^T S c) c_k, c <- c / sqrt(c^T S c). Loaded columns seed the set (kept
    // iff c^T S c > DEP_TOL, untouched if orthonormal to SKIP_TOL); X fills the rest at FILL_TOL.
    constexpr double FILL_TOL = 1.0e-8;
    // Squared S-norm below which a loaded column is spanned by its predecessors. Looser
    // than FILL_TOL: normalizing such a column amplifies noise by 1/sqrt(DEP_TOL).
    constexpr double DEP_TOL = 1.0e-6;
    // max|Cl^T S Cl - I|, Cl = C[:, 0..n_loaded-1], under which the loaded block needs
    // no repair.
    constexpr double SKIP_TOL = 1.0e-10;

    // G = Cl^T S Cl, worst = max|G - I| (lda = nmo selects Cl in C); T = S Cl is reused by loaded_ok.
    bool loaded_ok = false;
    std::vector<double> T;
    if (n_loaded > 0) {
        T.assign((size_t)nso * n_loaded, 0.0);
        std::vector<double> G((size_t)n_loaded * n_loaded);
        C_DGEMM('N', 'N', nso, n_loaded, nso, 1.0, Sp[0], nso, Cp[0], nmo, 0.0,
                T.data(), n_loaded);
        C_DGEMM('T', 'N', n_loaded, n_loaded, nso, 1.0, Cp[0], nmo, T.data(), n_loaded,
                0.0, G.data(), n_loaded);
        double worst = 0.0;
        for (int i = 0; i < n_loaded; ++i)
            for (int j = 0; j < n_loaded; ++j)
                worst = std::max(worst,
                                 std::abs(G[(size_t)i * n_loaded + j] - (i == j ? 1.0 : 0.0)));
        loaded_ok = (worst <= SKIP_TOL);
    }
    if (loaded_ok && n_loaded == nmo) return n_loaded;

    // Column-major scratch: C columns are strided by nmo, MGS inner loops need them
    // contiguous.
    std::vector<double> Ck((size_t)nmo * nso);
    std::vector<double> SCk((size_t)nmo * nso);
    std::vector<double> c(nso);

    std::vector<double> Sc(nso);

    // Sc = S c (S symmetric).
    auto form_Sc = [&]() {
        C_DGEMV('N', nso, nso, 1.0, Sp[0], nso, c.data(), 1, 0.0, Sc.data(), 1);
    };
    // Accept column k as inv*c, taking S c from the Sc formed for the same c.
    // Callers must have run form_Sc() (directly or through s_norm_sq) on this c.
    auto accept = [&](int k, double inv) {
        double* ck  = &Ck[(size_t)k * nso];
        double* sck = &SCk[(size_t)k * nso];
        C_DCOPY(nso, c.data(), 1, ck, 1);
        C_DCOPY(nso, Sc.data(), 1, sck, 1);
        if (inv != 1.0) {
            C_DSCAL(nso, inv, ck, 1);
            C_DSCAL(nso, inv, sck, 1);
        }
    };
    // Modified Gram-Schmidt, two passes, against accepted cols 0..upto-1. Sequential in
    // k, c updated in place:  c <- c - (sck . c) ck,  sck = S ck.
    auto orthogonalize = [&](int upto) {
        for (int pass = 0; pass < 2; ++pass) {
            for (int k = 0; k < upto; ++k) {
                const double* sck = &SCk[(size_t)k * nso];
                const double ov = C_DDOT(nso, sck, 1, c.data(), 1);
                if (ov == 0.0) continue;
                C_DAXPY(nso, -ov, &Ck[(size_t)k * nso], 1, c.data(), 1);
            }
        }
    };
    auto s_norm_sq = [&]() {
        form_Sc();
        return C_DDOT(nso, c.data(), 1, Sc.data(), 1);
    };

    // loaded_ok trusts the loaded columns as orthonormal; else each is re-projected via MGS.
    int kept = 0;
    if (loaded_ok) {
        // S c_j = column j of T.
        for (int j = 0; j < n_loaded; ++j) {
            for (int mu = 0; mu < nso; ++mu) c[mu] = Cp[mu][j];
            C_DCOPY(nso, &T[j], n_loaded, Sc.data(), 1);
            accept(j, 1.0);
        }
        kept = n_loaded;
    } else {
        for (int j = 0; j < n_loaded; ++j) {
            for (int mu = 0; mu < nso; ++mu) c[mu] = Cp[mu][j];
            orthogonalize(kept);
            const double sq = s_norm_sq();
            if (sq > DEP_TOL) {
                accept(kept, 1.0 / std::sqrt(sq));
                ++kept;
            }
        }
    }

    int x_cursor = 0;
    for (int j = kept; j < nmo; ++j) {
        bool placed = false;

        while (!placed && x_cursor < X->colspi()[h]) {
            for (int mu = 0; mu < nso; ++mu) c[mu] = Xp[mu][x_cursor];
            ++x_cursor;

            orthogonalize(j);

            const double sq = s_norm_sq();
            if (sq > FILL_TOL) {
                accept(j, 1.0 / std::sqrt(sq));
                placed = true;
            }
        }

        if (!placed) {
            throw PSIEXCEPTION(
                "REKS::complete_orbital_basis: ran out of canonical-orthogonalizer "
                "columns before filling Ca_. Possible degenerate AO basis.");
        }
    }

    // Every column 0..nmo-1 is written, so dropped columns leave no stale data.
    for (int mu = 0; mu < nso; ++mu)
        for (int j = 0; j < nmo; ++j) Cp[mu][j] = Ck[(size_t)j * nso + mu];
    return kept;
}

void REKS::guess() {
    std::string guess_type = options_.get_str("GUESS");

    if (guess_type != "READ" || !guess_Ca_) {
        HF::guess();
        return;
    }

    if (print_ && reports(2)) {
        outfile->Printf(
            "  SCF Guess: REKS orbitals restored from previous computation "
            "(REKS::guess override).\n\n");
    }

    // Cb_ aliases Ca_ (reks_common_init, same_a_b_orbs_ holds): only alpha is read.
    // A restricted source loses nothing; an unrestricted one drops its beta.
    if (guess_Cb_ && guess_Cb_ != guess_Ca_ && !guess_Cb_->equal(guess_Ca_) && reports(2)) {
        outfile->Printf(
            "  REKS READ: WARNING -- the guess wavefunction is unrestricted (Cb != Ca).\n"
            "    REKS optimizes a single orbital set, so only the ALPHA coefficients are\n"
            "    read; the beta set is discarded.\n");
    }

    int n_loaded_a[8] = {0};
    for (int h = 0; h < nirrep_; ++h) {
        const int ncols = std::min(guess_Ca_->colspi()[h], nmopi_[h]);
        n_loaded_a[h] = ncols;
        // Both sides are row-major: a row segment of the leading ncols columns is contiguous.
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            C_DCOPY(ncols, guess_Ca_->pointer(h)[mu], 1, Ca_->pointer(h)[mu], 1);
        }
    }

    auto report_dropped = [&](const char* label, int h, int kept, int n_loaded) {
        if (kept >= n_loaded || !reports(2)) return;
        outfile->Printf(
            "  REKS READ: WARNING -- %s irrep %d: %d of %d loaded MOs were linearly "
            "dependent in the current AO metric and were dropped; refilled from the\n"
            "    canonical orthogonalizer. Columns after the first drop shifted, so the\n"
            "    active-window layout may no longer match guess_active_window.\n",
            label, h, n_loaded - kept, n_loaded);
    };
    for (int h = 0; h < nirrep_; ++h) {
        int kept = complete_orbital_basis(Ca_, S_, X_, n_loaded_a[h], h);
        report_dropped("Ca", h, kept, n_loaded_a[h]);
    }

    // READ guess -> reel_.fon_state[s]: gen 0 from guess_fon_[s], gen >= 1 from guess_fon_upper_[s][gen];
    // FONs clamp into (0, 2), a size mismatch drops that generation.
    if (catalog_.data) {
        constexpr double fon_bnd_margin = 1e-8;
        auto restore_fon_layer = [&](int s, int gen, SharedMatrix& src) {
            if (!src) return;
            const auto& act = sa_cassette_.geminals_active(s, gen);
            const int n    = static_cast<int>(act.size());
            const int rows = src->nrow();
            const int got  = rows * src->ncol();
            const char letter = reks::studio::generation_letter(gen);
            const std::string chan = reks::studio::fon_channel_name(gen, s);
            if (got != n) {
                if (reports(2))
                    outfile->Printf(
                        "  REKS READ: WARNING -- %s-FON guess size %d does not match active "
                        "count %d; falling back to default FON=1.0.\n", chan.c_str(), got, n);
                src.reset();
                return;
            }
            double** Gp = src->pointer(0);
            const bool narrate = reports(2);
            if (narrate)
                outfile->Printf("  REKS READ: restoring %s-FON from saved wavefunction:\n",
                                chan.c_str());
            for (int k = 0; k < n; ++k) {
                double f_raw = (rows == 1) ? Gp[0][k] : Gp[k][0];
                double f = f_raw;
                if (f > 2.0 - fon_bnd_margin) f = 2.0 - fon_bnd_margin;
                if (f < 0.0 + fon_bnd_margin) f = 0.0 + fon_bnd_margin;
                // A FON saved at the boundary reloads equal to it; report only real moves.
                const bool clamped = (f != f_raw);
                reel_.fon_state[s].layers[gen][act[k]] = {f, 2.0 - f};
                if (!narrate) continue;
                if (clamped)
                    outfile->Printf("    %c[%d] = %.12f  (clamped from %.12f)\n", letter, k, f, f_raw);
                else
                    outfile->Printf("    %c[%d] = %.12f\n", letter, k, f);
            }
        };
        for (int s = 0; s < sa_cassette_.n_sectors(); ++s) {
            if (s < static_cast<int>(guess_fon_.size()))
                restore_fon_layer(s, 0, guess_fon_[s]);
            for (int gen = 1; gen < sa_cassette_.n_generations(); ++gen)
                if (s < static_cast<int>(guess_fon_upper_.size()))
                    restore_fon_layer(s, gen, guess_fon_upper_[s][gen]);
        }
    }

    format_guess();
    form_D();
    iteration_ = -1;
}

void REKS::finalize() {
    // SA base densities on Ca_fock_transform, which differs from the exported Ca_ on a swap exit.
    if (catalog_.data && reel_.Ca_fock_transform &&
        !call_sheet_.sa_referenced_bitmask.empty()) {
        reks::scf::build_base_densities(
            sa_cassette_, call_sheet_.sa_referenced_bitmask, reel_, Ncore_,
            active_mo_indices_, reel_.Ca_fock_transform, H_);
    }

    // Print the converged E[SA-REKS] energy decomposition (per-config + extra-microstate terms).
    if (catalog_.data && !sa_cassette_.K_indices.empty()) {
        std::vector<double> e_config;
        e_config.reserve(sa_cassette_.K_indices.size());
        // Each config's energy is evaluated against its owning sector's FON snapshot,
        // resolved through the init-time run-sector map (identity at n_sectors==1).
        for (int K : sa_cassette_.K_indices)
            e_config.push_back(get_config_energy(K, run_sector_of_config(K)));
        std::vector<double> e_extra;
        const int n_extra = sa_cassette_.n_extra_microstates();
        e_extra.reserve(n_extra);
        for (int e = 0; e < n_extra; ++e)
            e_extra.push_back(reel_.E_L[sa_cassette_.n_catalog_microstates() + e]);
        reks::report::print_sa_energy_decomposition(
            sa_cassette_, e_config, e_extra, energies_["REKS SA"]);
    }

    if (gvb_ && gvb_->verdict_trips() > 0 && reports(4))
        outfile->Printf("  [GVB-VERDICT] iter=%d rollbacks=%d best_gorb=%.3e\n", iteration_,
                        gvb_->verdict_trips(), gvb_->best_gorb());

    if (gvb_) report_convergence_quality();

    // FONs per generation, active-space metadata and "REKS N SECTORS"; sector 0 keeps unqualified names,
    // s > 0 prefixes "SECTOR s ".
    if (catalog_.data) {
        set_scalar_variable("REKS N SECTORS", static_cast<double>(sa_cassette_.n_sectors()));
        for (int s = 0; s < sa_cassette_.n_sectors(); ++s) {
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
                const auto& act = sa_cassette_.geminals_active(s, gen);
                const int n = static_cast<int>(act.size());
                if (n == 0) continue;
                const std::string base =
                    reks::report::sector_psivar_name(s, reks::studio::fon_array_name(gen));
                auto fon_mat  = std::make_shared<Matrix>(base, 1, n);
                auto fon_meta = std::make_shared<Matrix>(base + " META", 3, n);
                for (int k = 0; k < n; ++k) {
                    const int g = act[k];
                    const auto& tmpl = sa_cassette_.geminal_templates()[g];
                    fon_mat->set(0, k, reel_.fon_state[s].layers[gen][g].p);
                    fon_meta->set(0, 0, k, static_cast<double>(tmpl.scheme));
                    fon_meta->set(0, 1, k, static_cast<double>(tmpl.orbitals[0]));
                    fon_meta->set(0, 2, k, static_cast<double>(tmpl.orbitals[1]));
                }
                set_array_variable(base, fon_mat);
                set_array_variable(base + " META", fon_meta);
            }
        }
        set_scalar_variable("REKS N ACTIVE ELECTRONS",
                            static_cast<double>(sa_cassette_.n_electrons()));
        set_scalar_variable("REKS N ACTIVE ORBITALS",
                            static_cast<double>(sa_cassette_.n_active_orbitals()));
        set_scalar_variable("REKS SCHEME",
                            static_cast<double>(sa_cassette_.scheme()));
        set_scalar_variable("REKS G",
                            static_cast<double>(sa_cassette_.n_orbitals_per_geminal()));
        // W_PPS only exists for catalogs that carry a PPS1 config (S=0 sectors).
        const int k_pps = reks::studio::find_config_index(catalog_, "PPS1");
        if (k_pps >= 0)
            set_scalar_variable("REKS W_PPS", sa_cassette_.w_configs[k_pps]);
        if (!sa_cassette_.K_weights.empty()) {
            const auto& saw = sa_cassette_.K_weights;
            auto saw_mat = std::make_shared<Matrix>("REKS SA_WEIGHTS", 1,
                                                    static_cast<int>(saw.size()));
            for (size_t i = 0; i < saw.size(); ++i) {
                saw_mat->set(0, static_cast<int>(i), saw[i]);
            }
            set_array_variable("REKS SA_WEIGHTS", saw_mat);
        }
        if (reports(4)) {
            outfile->Printf("  REKS SAVE: persisted FON vectors and active-space metadata (");
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen)
                outfile->Printf("%s%s=%zu", gen ? ", " : "",
                                reks::studio::fon_channel_name(gen, sa_sector()).c_str(),
                                sa_cassette_.geminals_active_gen(gen).size());
            outfile->Printf(")\n");
        }
    }

    if (!si_cassettes_.empty()) {
        {
            // Entry i = the cassette's position in SI_REKS_CONFIGS (0..N-1).
            const int N = static_cast<int>(si_cassettes_.size());
            auto list = std::make_shared<Matrix>("SSR SI LIST", 1, N);
            for (int i = 0; i < N; ++i)
                list->set(0, 0, i, static_cast<double>(i));
            set_array_variable("SSR SI LIST", list);
        }

        si_computed_ = false;
        compute_si();

        const reks::report::ReportFrame frame = report_frame();

        // Cassette-independent dipole integrals, built once (frozen orbitals).
        reks::properties::DipoleIntegralCache dip_cache;
        if (reks::studio::has_diabatic_rdm(catalog_)) {
            timer_on("REKS: dipole_mo_cache");
            dip_cache = reks::properties::build_dipole_integral_cache(frame, molecule_, basisset_);
            timer_off("REKS: dipole_mo_cache");
        }

        bool mo_consistency_printed = false;
        // Cassettes are grouped by sector; table tags restart on each sector change.
        int tagged_sector = -1;
        for (size_t i = 0; i < si_cassettes_.size(); ++i) {
            const reks::studio::Cassette& cassette = si_cassettes_[i];
            const reks::SIResult&         sir      = si_results_[i];
            const int si_cassette_pos = static_cast<int>(i);
            const int spin2 = catalog_.sector_entry(run_to_blob_[cassette.sector()]).spin2;

            if (cassette.sector() != tagged_sector) {
                reks::report::reset_table_tags();
                tagged_sector = cassette.sector();
            }

            reks::rdm::DiabaticRdm rho_diab;
            if (reks::studio::has_diabatic_rdm(catalog_)) {
                timer_on("REKS: build_diabatic_rdm");
                rho_diab = reks::rdm::diabatic(reel_, cassette);
                timer_off("REKS: build_diabatic_rdm");
            }

            reks::properties::SIProperties props =
                reks::properties::compute_state_properties(frame, cassette, sir, rho_diab, dip_cache);

            if (reks::scf::sa_energies_ready(sa_cassette_, reel_.E_L)) {
                if (reports(5) && !mo_consistency_printed) {
                    reks::report::print_mo_consistency(frame);
                    mo_consistency_printed = true;
                }
                if (reports(5)) reks::report::print_f_reks_mo_debug(frame);

                reks::report::print_state_properties_header(frame, cassette, si_cassette_pos, spin2);

                reks::report::print_microstate_energies(
                    sa_cassette_, cassette.microstates_active, reel_.E_L);
                reks::report::print_fon_table(frame, cassette);

                reks::report::print_lagrangian_eri_table(
                    sa_cassette_, reel_.lagrangians, active_eri_pool_);

                reks::report::print_si_hamiltonian_overlap(frame, cassette, sir);
                reks::report::print_adiabatic_energies(sir);
                reks::report::print_excitation_energies(sir);

                reks::properties::print_state_dipoles(props, sir);
            }

            reks::properties::register_si_wfn_variables(*this, frame, si_cassette_pos,
                                                        cassette.sector(), sir, props,
                                                        /*primary_alias=*/ i == 0, rho_diab);

            auto nos = reks::properties::compute_state_natural(catalog_, cassette, sir,
                                                               props.rho_adiab);
            reks::properties::print_state_natural(frame, nos, si_cassette_pos, cassette, sir,
                                                  rho_diab, spin2);
            reks::properties::register_natural_wfn_variables(*this, frame, si_cassette_pos,
                                                             cassette.sector(), nos,
                                                             /*primary_alias=*/ i == 0);
        }
    }

    ms_xc_.reset();
    if (rv_potential_) {
        rv_potential_->finalize();
    }

    // L_mn = sum_{i < nalphapi_[h]} eps_i C_mi C_ni over the leading columns of Ca_ (stride nmopi_[h]);
    // zeroed first, C_DGEMM skips nocc == 0.
    Lagrangian_->zero();
    for (int h = 0; h < nirrep_; ++h) {
        const int nrow = Lagrangian_->rowdim(h);
        const int ncol = Lagrangian_->coldim(h);
        const int nocc = nalphapi_[h];
        if (nrow == 0 || ncol == 0 || nocc == 0) continue;

        double**      Cp  = Ca_->pointer(h);
        const double* eps = epsilon_a_->pointer(h);
        std::vector<double> Ceps(static_cast<size_t>(nrow) * nocc);
        for (int m = 0; m < nrow; ++m)
            for (int i = 0; i < nocc; ++i)
                Ceps[static_cast<size_t>(m) * nocc + i] = eps[i] * Cp[m][i];

        C_DGEMM('N', 'T', nrow, ncol, nocc, 1.0, Ceps.data(), nocc, Cp[0],
                Ca_->coldim(h), 0.0, Lagrangian_->pointer(h)[0], ncol);
    }

    stage_Ca_.reset();
    stage_eps_.reset();
    stage_fons_.clear();
    staged_ = false;
    G_.reset();
    J_.reset();
    K_.reset();
    wK_.reset();

    HF::finalize();
}

SharedMatrix REKS::build_generalized_fock() const {
    timer_on("REKS: build_generalized_fock");

    int nmo = Ca_->colspi()[0];
    int n_active = sa_cassette_.n_active_orbitals();

    if (!F_gen_workspace_ || F_gen_workspace_->rowdim(0) != nmo ||
        F_gen_workspace_->coldim(0) != nmo) {
        F_gen_workspace_ = std::make_shared<Matrix>("Generalized Fock", nmo, nmo);
    }
    SharedMatrix F_gen = F_gen_workspace_;
    F_gen->zero();
    double** Fgp = F_gen->pointer(0);

    // L is the global microstate index.
    for (int L : sa_cassette_.microstates()) {
        double C_L = reel_.C_L[L];

        if (std::abs(C_L) < 1e-14) continue;

        // row a_idx = MO active_mo_indices_[a_idx].
        double** Farows_a = reel_.F_MO_arows_a_L[L]->pointer(0);
        double** Farows_b = reel_.F_MO_arows_b_L[L]->pointer(0);
        const auto& m = sa_cassette_.microstate(L);

        for (int a_idx = 0; a_idx < n_active; ++a_idx) {
            const int n_alpha = m.alpha[a_idx];
            const int n_beta  = m.beta[a_idx];
            if (!n_alpha && !n_beta) continue;

            // n_alpha / n_beta are 0/1 occupations: an unoccupied spin contributes nothing.
            const int a = active_mo_indices_[a_idx];
            if (n_alpha) C_DAXPY(nmo, C_L, Farows_a[a_idx], 1, Fgp[a], 1);
            if (n_beta) C_DAXPY(nmo, C_L, Farows_b[a_idx], 1, Fgp[a], 1);
        }
    }

    // F_gen[i,q] = F_acc_MO[i,q] for i < Ncore_: n_alpha = n_beta = 1 for core in
    // every microstate, so sum_L C_L (Fa_L[i,q]+Fb_L[i,q]) is exactly F_acc_MO[i,q].
    {
        double** Facc = reel_.F_acc_MO->pointer(0);
        for (int i = 0; i < Ncore_; ++i) C_DAXPY(nmo, 1.0, Facc[i], 1, Fgp[i], 1);
    }

    if (reports(5)) {
        outfile->Printf("\n  === Generalized Fock Matrix Asymmetry ===\n");

        double max_asym_ca = 0.0;
        for (int i = 0; i < Ncore_; ++i) {
            for (int a_idx = 0; a_idx < n_active; ++a_idx) {
                int a = active_mo_indices_[a_idx];
                double asym = std::abs(Fgp[i][a] - Fgp[a][i]);
                if (asym > max_asym_ca) max_asym_ca = asym;
            }
        }
        outfile->Printf("  Max core-active asymmetry: %.6e\n", max_asym_ca);

        double max_asym_aa = 0.0;
        for (int i = 0; i < n_active; ++i) {
            for (int j = i + 1; j < n_active; ++j) {
                int ai = active_mo_indices_[i];
                int aj = active_mo_indices_[j];
                double asym = std::abs(Fgp[ai][aj] - Fgp[aj][ai]);
                if (asym > max_asym_aa) max_asym_aa = asym;
            }
        }
        outfile->Printf("  Max active-active asymmetry: %.6e\n", max_asym_aa);

        int last_active = active_mo_indices_.back();
        double max_asym_av = 0.0;
        for (int a_idx = 0; a_idx < n_active; ++a_idx) {
            int a = active_mo_indices_[a_idx];
            for (int v = last_active + 1; v < nmo; ++v) {
                double asym = std::abs(Fgp[a][v] - Fgp[v][a]);
                if (asym > max_asym_av) max_asym_av = asym;
            }
        }
        outfile->Printf("  Max active-virtual asymmetry: %.6e\n", max_asym_av);
    }

    timer_off("REKS: build_generalized_fock");
    return F_gen;
}

void REKS::build_F_acc_MO() {
    timer_on("REKS: build_F_acc_MO");

    const int  n_active = sa_cassette_.n_active_orbitals();
    const bool needs_xc = functional_->needs_xc();
    const reks::scf::FockExchangePolicy pol =
        reks::scf::fock_exchange_policy(reel_, *functional_);

    // sumC = Sum_L C_L; wJ[i] = Sum_L C_L*(na_i+nb_i).
    double sumC = 0.0;
    std::vector<double> wJ(n_active, 0.0);
    for (int L : sa_cassette_.microstates()) {
        const double C_L = reel_.C_L[L];
        sumC += C_L;
        const auto& m = sa_cassette_.microstate(L);
        for (int i = 0; i < n_active; ++i)
            wJ[i] += C_L * (m.alpha[i] + m.beta[i]);
    }
    std::vector<double> wJ2(n_active);
    for (int i = 0; i < n_active; ++i) wJ2[i] = 2.0 * wJ[i];

    // F_acc_AO = Sum_L C_L (Fa_AO_L + Fb_AO_L). Coulomb weight 2*wJ (spin-shared,
    // doubled by the both-spin add), exchange weight wJ; core weights 4*sumC : 2*sumC.
    SharedMatrix F_acc_AO = reel_.F_acc_AO;
    F_acc_AO->copy(H_);
    F_acc_AO->scale(2.0 * sumC);
    reks::scf::add_coulomb_pattern(F_acc_AO, reel_, Ncore_, n_active, 4.0 * sumC, wJ2.data());
    reks::scf::add_exchange_pattern(F_acc_AO, reel_, Ncore_, n_active, pol, 2.0 * sumC, wJ.data());
    if (needs_xc) {
        // V_acc_xc_{a,b} = Sum_L C_L V_xc_{a,b}^L (AO).
        F_acc_AO->add(reel_.V_acc_xc_a);
        F_acc_AO->add(reel_.V_acc_xc_b);
    }

    // AO -> MO: F_acc_MO = C^T F_acc_AO C (both spins folded in).
    const SharedMatrix& Ca = reel_.Ca_fock_transform;
    const int N  = Ca->rowspi()[0];
    const int Mn = Ca->colspi()[0];
    if (!reel_.F_acc_MO || reel_.F_acc_MO->coldim(0) != Mn) {
        reel_.F_acc_MO = std::make_shared<Matrix>("F_acc_MO", Mn, Mn);
    }
    double** Cp  = Ca->pointer(0);
    double** Tp  = reel_.temp_buffer->pointer(0);
    double** Aap = F_acc_AO->pointer(0);
    double** Amp = reel_.F_acc_MO->pointer(0);
    // T(nso x nmo) = F_AO(nso x nso) C(nso x nmo); F_acc_MO(nmo x nmo) = C^T T.
    C_DGEMM('N', 'N', N, Mn, N, 1.0, Aap[0], N, Cp[0], Mn, 0.0, Tp[0], N);
    C_DGEMM('T', 'N', Mn, Mn, N, 1.0, Cp[0], Mn, Tp[0], N, 0.0, Amp[0], Mn);

    timer_off("REKS: build_F_acc_MO");
}

double REKS::compute_full_gradient_norm(SharedMatrix F_gen, const std::vector<char>& frozen) const {
    return reks::orbital_gradient_norm(F_gen->pointer(0), active_mo_indices_, Ncore_, nmopi_[0],
                                       frozen, nullptr);
}

std::vector<char> REKS::same_role_aa_frozen() const {
    const int n_act = static_cast<int>(active_mo_indices_.size());
    const int n_rot = n_act * (n_act - 1) / 2;
    std::vector<char> frozen(n_rot, 0);
    if (n_rot == 0 || !catalog_.data) return frozen;

    // unit = geminal of each active orbital, role = its slot in that geminal (0 = high occupancy).
    const auto& gems = sa_cassette_.geminals_active_gen(0);
    if (gems.size() < 2) return frozen;
    std::vector<int> unit(n_act, -1), role(n_act, -1);
    for (int u = 0; u < static_cast<int>(gems.size()); ++u) {
        const auto& orbs = sa_cassette_.geminal_templates()[gems[u]].orbitals;
        for (int r = 0; r < static_cast<int>(orbs.size()); ++r) {
            unit[orbs[r]] = u;
            role[orbs[r]] = r;
        }
    }

    int idx = 0;
    for (int i = 0; i < n_act; ++i)
        for (int j = i + 1; j < n_act; ++j, ++idx)
            if (unit[i] >= 0 && unit[j] >= 0 && unit[i] != unit[j] && role[i] == role[j])
                frozen[idx] = 1;
    return frozen;
}

SharedMatrix REKS::build_uniform_active_fock() const {
    const int n_act = static_cast<int>(active_mo_indices_.size());
    if (n_act == 0 || static_cast<int>(reel_.J_act.size()) != n_act) return nullptr;

    // Uniform active occupation q = N_act / n_act, split evenly over the two spins.
    const double q = static_cast<double>(nelectron_ - 2 * Ncore_) / static_cast<double>(n_act);
    std::vector<double> w_ab(n_act, q), w_a(n_act, 0.5 * q);

    reks::scf::FockExchangePolicy pol = reks::scf::fock_exchange_policy(reel_, *functional_);

    auto F_avg = std::make_shared<Matrix>("F_avg AO", nsopi_, nsopi_);
    F_avg->copy(H_);
    reks::scf::add_coulomb_pattern(F_avg, reel_, Ncore_, n_act, 2.0, w_ab.data());
    reks::scf::add_exchange_pattern(F_avg, reel_, Ncore_, n_act, pol, 1.0, w_a.data());
    return F_avg;
}

SharedMatrix REKS::same_role_canonicalizer(std::vector<double>* eigs, double* offdiag_norm) const {
    if (eigs) eigs->clear();
    if (offdiag_norm) *offdiag_norm = 0.0;
    const int n_act = static_cast<int>(active_mo_indices_.size());
    if (n_act < 2 || !catalog_.data) return nullptr;
    const auto& gems = sa_cassette_.geminals_active_gen(0);
    if (gems.size() < 2) return nullptr;

    auto F_avg = build_uniform_active_fock();
    if (!F_avg) return nullptr;

    // Role blocks, one orbital of every geminal each, collected by ascending active slot.
    constexpr int kRoles = 2;
    std::array<std::vector<int>, kRoles> block;
    {
        std::vector<int> role(n_act, -1);
        for (int u = 0; u < static_cast<int>(gems.size()); ++u) {
            const auto& orbs = sa_cassette_.geminal_templates()[gems[u]].orbitals;
            for (int r = 0; r < static_cast<int>(orbs.size()) && r < kRoles; ++r)
                role[orbs[r]] = r;
        }
        for (int k = 0; k < n_act; ++k)
            if (role[k] >= 0) block[role[k]].push_back(k);
    }

    const int nso = Ca_->rowspi()[0];
    double** Cp = Ca_->pointer(0);
    double** Fp = F_avg->pointer(0);

    // T = F_avg C_act, F_act = C_act^T T.
    auto C_act = std::make_shared<Matrix>("C_act", nso, n_act);
    double** Cap = C_act->pointer(0);
    for (int k = 0; k < n_act; ++k)
        for (int mu = 0; mu < nso; ++mu) Cap[mu][k] = Cp[mu][active_mo_indices_[k]];
    auto T = std::make_shared<Matrix>("F_avg C_act", nso, n_act);
    double** Tp = T->pointer(0);
    C_DGEMM('N', 'N', nso, n_act, nso, 1.0, Fp[0], nso, Cap[0], n_act, 0.0, Tp[0], n_act);
    auto F_act = std::make_shared<Matrix>("F_avg act", n_act, n_act);
    double** Fap = F_act->pointer(0);
    C_DGEMM('T', 'N', n_act, n_act, nso, 1.0, Cap[0], n_act, Tp[0], n_act, 0.0, Fap[0], n_act);

    // V = identity outside the role blocks, role-block eigenvectors inside.
    auto V = std::make_shared<Matrix>("role canonicalizer", n_act, n_act);
    double** Vp = V->pointer(0);
    for (int k = 0; k < n_act; ++k) Vp[k][k] = 1.0;
    double offdiag_sq = 0.0;
    for (int r = 0; r < kRoles; ++r) {
        const auto& sl = block[r];
        const int m = static_cast<int>(sl.size());
        if (m < 2) continue;
        auto B = std::make_shared<Matrix>("F_avg role block", m, m);
        double** Bp = B->pointer(0);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j) {
                Bp[i][j] = 0.5 * (Fap[sl[i]][sl[j]] + Fap[sl[j]][sl[i]]);
                if (i != j) offdiag_sq += Bp[i][j] * Bp[i][j];
            }
        auto Bv = std::make_shared<Matrix>("role evecs", m, m);
        auto Be = std::make_shared<Vector>("role evals", m);
        B->diagonalize(Bv, Be, ascending);
        double** Bvp = Bv->pointer(0);
        // Eigenvector phase fixed by a positive diagonal.
        for (int c = 0; c < m; ++c)
            if (Bvp[c][c] < 0.0)
                for (int i = 0; i < m; ++i) Bvp[i][c] = -Bvp[i][c];
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j) Vp[sl[i]][sl[j]] = Bvp[i][j];
        if (eigs)
            for (int i = 0; i < m; ++i) eigs->push_back(Be->get(i));
    }
    if (offdiag_norm) *offdiag_norm = std::sqrt(offdiag_sq);
    return V;
}

void REKS::canonicalize_same_role_frame() {
    std::vector<double> eig_report;
    double offdiag = 0.0;
    auto V = same_role_canonicalizer(&eig_report, &offdiag);
    if (!V) return;

    const int n_act = static_cast<int>(active_mo_indices_.size());
    const int nso = Ca_->rowspi()[0];
    const int nmo = Ca_->colspi()[0];
    double** Cp = Ca_->pointer(0);
    double** Vp = V->pointer(0);

    double dev_max = 0.0;
    for (int i = 0; i < n_act; ++i)
        for (int j = 0; j < n_act; ++j)
            dev_max = std::max(dev_max, std::abs(Vp[i][j] - (i == j ? 1.0 : 0.0)));

    // Active columns of Ca_ <- active columns * V.
    {
        std::vector<double> col(n_act, 0.0);
        for (int mu = 0; mu < nso; ++mu) {
            for (int j = 0; j < n_act; ++j) {
                double s = 0.0;
                for (int i = 0; i < n_act; ++i) s += Cp[mu][active_mo_indices_[i]] * Vp[i][j];
                col[j] = s;
            }
            for (int j = 0; j < n_act; ++j) Cp[mu][active_mo_indices_[j]] = col[j];
        }
    }

    // eps of an active slot = Fock diagonal c_a^T Fa c_a in the rotated frame.
    {
        double** Fa_p = Fa_->pointer(0);
        double* eps_a = epsilon_a_->pointer(0);
        std::vector<double> t(nso, 0.0);
        for (int k = 0; k < n_act; ++k) {
            const int a = active_mo_indices_[k];
            // t = Fa c_a; incx = nmo walks column a of the row-major Ca_.
            C_DGEMV('N', nso, nso, 1.0, Fa_p[0], nso, &Cp[0][a], nmo, 0.0, t.data(), 1);
            eps_a[a] = C_DDOT(nso, t.data(), 1, &Cp[0][a], nmo);
        }
        epsilon_b_->copy(*epsilon_a_);
    }

    if (reports(4)) {
        outfile->Printf("  [CANON] iter=%d F_avg role offdiag=%.3e max|V-I|=%.3e eps_avg=",
                        iteration_, offdiag, dev_max);
        for (double e : eig_report) outfile->Printf(" %+.6f", e);
        outfile->Printf("\n");
    }
}

bool REKS::trah_live() const { return use_trah_ && !(gvb_ && gvb_->active()); }

bool REKS::convergence_held() const {
    return gvb_ && (gvb_->gate_confirm_pending() || gvb_->landed_iter() == iteration_);
}

void REKS::print_orbitals() {
    if (!reports(2)) return;

    const int h = 0;
    const int nmo = nmopi_[h];
    const auto labels_irrep = molecule_->irrep_labels();
    const int n_active = static_cast<int>(active_mo_indices_.size());

    std::vector<int> mo_to_active(nmo, -1);
    for (int k = 0; k < n_active; ++k) {
        const int mo = active_mo_indices_[k];
        if (mo >= 0 && mo < nmo) mo_to_active[mo] = k;
    }

    const std::vector<int> rank = reks::report::mo_energy_rank(epsilon_a_, nmo);

    using OrbEntry = std::pair<double, std::pair<std::string, int>>;
    std::vector<OrbEntry> core, vir;
    std::vector<std::tuple<double, int, int>> active;
    for (int a = 0; a < nmo; ++a) {
        const double e = epsilon_a_->get(h, a);
        const int rk = rank[a] + 1;
        if (mo_to_active[a] >= 0) {
            active.emplace_back(e, rk, mo_to_active[a]);
        } else if (a < Ncore_) {
            core.emplace_back(e, std::make_pair(labels_irrep[h], rk));
        } else {
            vir.emplace_back(e, std::make_pair(labels_irrep[h], rk));
        }
    }
    std::sort(core.begin(), core.end());
    std::sort(vir.begin(), vir.end());
    // Sorted by active-window slot (apos), not energy: labels a,b,c,... stay
    // fixed to their geminal slots.
    std::sort(active.begin(), active.end(),
              [](const std::tuple<double, int, int>& x, const std::tuple<double, int, int>& y) {
                  return std::get<2>(x) < std::get<2>(y);
              });

    outfile->Printf("    Orbital Energies [Eh]\n    ---------------------\n\n");

    char hdr[80];
    std::snprintf(hdr, sizeof(hdr), "Core Doubly Occupied (%d):", static_cast<int>(core.size()));
    print_orbital_pairs(hdr, core);

    reks::report::print_active_orbitals(report_frame(), rank, active);

    std::snprintf(hdr, sizeof(hdr), "Virtual (%d):", static_cast<int>(vir.size()));
    print_orbital_pairs(hdr, vir);
}

reks::report::ReportFrame REKS::report_frame() const {
    return {catalog_, sa_cassette_, reel_, active_mo_indices_, Ncore_, nmopi_[0],
            Ca_, epsilon_a_, molecule_->irrep_labels()[0]};
}

void REKS::build_fon_blocks() {
    // fon_blocks_ = {(s, gen, 0, count) : count = |geminals_active(s, gen)| > 0}, ordered by (s, gen);
    // offsets are assigned later.
    fon_blocks_.clear();
    for (int s = 0; s < sa_cassette_.n_sectors(); ++s)
        for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
            const int count = static_cast<int>(sa_cassette_.geminals_active(s, gen).size());
            if (count == 0) continue;
            fon_blocks_.push_back(reks::studio::FonBlock{s, gen, 0, count});
        }

    if (trah_) trah_->resize_fon_blocks(fon_blocks_);
}

// One-shot setup: reads REKS options, builds the (N,M) catalog/cassette pool and
// run-sector partition, allocates REKS matrices, and configures the FON/DIIS/TRAH solvers.
void REKS::reks_common_init() {
    Fa_ = SharedMatrix(factory_->create_matrix("F"));
    Ca_ = SharedMatrix(factory_->create_matrix("MO coefficients (C)"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF density"));
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("orbital energies");
    Lagrangian_ = SharedMatrix(factory_->create_matrix("X"));
    Va_ = SharedMatrix(factory_->create_matrix("V"));

    // REKS is restricted: beta aliases alpha for the whole run; no path gives the beta
    // handles their own storage.
    Fb_ = Fa_;
    Cb_ = Ca_;
    Db_ = Da_;
    epsilon_b_ = epsilon_a_;
    Vb_ = Va_;

    G_ = SharedMatrix(factory_->create_matrix("G"));
    J_ = SharedMatrix(factory_->create_matrix("J"));
    K_ = SharedMatrix(factory_->create_matrix("K"));
    wK_ = SharedMatrix(factory_->create_matrix("wK"));

    same_a_b_orbs_ = true;
    same_a_b_dens_ = true;

    // Set before subclass_init() so setup_potential()'s diagnostics honor it.
    const int report_level = options_.get_int("REKS_REPORT_LEVEL");
    if (report_level < reks::report::kMinLevel || report_level > reks::report::kMaxLevel)
        reks::report::input_error(
            "REKS_REPORT_LEVEL must be " + std::to_string(reks::report::kMinLevel) + ".." +
            std::to_string(reks::report::kMaxLevel) + ", got " + std::to_string(report_level),
            false);
    reks::report::set_level(report_level);

    subclass_init();

    name_ = "REKS";

    si_computed_ = false;

    // FON interpolant delta (process-wide global).
    reks::delta_interp = options_.get_double("REKS_FON_INTERP_DELTA");
    if (reks::delta_interp < 0.0)
        reks::report::input_error("REKS_FON_INTERP_DELTA must be >= 0.", false);

    // Per-generation FON floor from REKS_<G>_FON (G = N/M/U/V/W/X/Y/Z for gen 0..7).
    // floor >= 0 pins generation g FONs to [floor, fon_floor_max]; -1 disables.
    constexpr double fon_floor_max = 2.0 - 1e-8;
    for (int gen = 0; gen < reks::studio::kMaxGen; ++gen) {
        char up = static_cast<char>(reks::studio::generation_letter(gen) - 'a' + 'A');
        std::string opt = std::string("REKS_") + up + "_FON";
        double floor = options_.get_double(opt);
        if (floor > fon_floor_max) {
            if (reports(2))
                outfile->Printf("  [REKS] WARNING: %s=%.8f exceeds maximum valid floor "
                                "%.8f; clamped (FON frozen at upper boundary).\n",
                                opt.c_str(), floor, fon_floor_max);
            floor = fon_floor_max;
        }
        reks_fon_lower_[gen] = floor;
        if (floor >= 0.0 && reports(4))
            outfile->Printf("  [REKS] %s lower bound = %.8f\n", opt.c_str(), floor);
    }
    // Micro-Newton FON solver defaults; lower_bound is replaced per-generation before each solve.
    reks::MicroSolverConfig fon_micro_cfg;
    fon_micro_cfg.max_nr_iter     = options_.get_int   ("REKS_FON_MICRO_NR_MAXITER");
    fon_micro_cfg.max_ls_iter     = options_.get_int   ("REKS_FON_MICRO_LS_MAXITER");
    fon_micro_cfg.convergence_tol = options_.get_double("REKS_FON_MICRO_CONVERGENCE");
    fon_micro_cfg.min_eigenvalue  = options_.get_double("REKS_FON_MICRO_MIN_EIGENVALUE");
    fon_micro_cfg.lower_bound     = options_.get_double("REKS_FON_MICRO_BOUND_MARGIN");
    fon_micro_cfg.upper_bound     = 2.0 - fon_micro_cfg.lower_bound;
    fon_step_ = std::make_unique<reks::FonStep>();
    fon_step_->configure(fon_micro_cfg, options_.get_double("REKS_FON_MICRO_BOUNDARY_SKIN"));

    orbital_guard_.set_adaptive(options_.get_bool("LEVEL_SHIFT_ADAPT"));

    if (!options_["REKS"].has_changed()) {
        reks::report::input_error(
            "REKS active space not specified: set the 'reks' option to [N, M].");
    }
    // Active space [N, M]; spin manifolds come from SA/SI_REKS_2SPIN or the ladder order of nested blocks.
    if (options_["REKS"].size() != 2) {
        reks::report::input_error(
            "REKS: the 'reks' option must have exactly 2 elements [N, M], got " +
            std::to_string(options_["REKS"].size()) + ".");
    }
    const int reks_N = static_cast<int>(options_["REKS"][0].to_double());
    const int reks_M = static_cast<int>(options_["REKS"][1].to_double());
    // Supported: REKS(N,N), N even >= 2, singlet molecule, C1.
    if (reks_N != reks_M || reks_N < 2 || reks_N % 2 != 0 || multiplicity_ != 1 || nirrep_ != 1) {
        reks::report::input_error(
            "REKS supports only REKS(N,N) with N even and N >= 2, molecule multiplicity 1 and C1 "
            "symmetry; got reks [" + std::to_string(reks_N) + ", " + std::to_string(reks_M) +
            "], multiplicity " + std::to_string(multiplicity_) + ", " + std::to_string(nirrep_) +
            " irreps. Keep multiplicity 1 and select spin manifolds with SA_REKS_2SPIN / "
            "SI_REKS_2SPIN; add 'symmetry c1' to the molecule.");
    }

    if (!reks::studio::variant_supported(reks_N, reks_M)) {
        reks::report::input_error(
            "REKS active space [N=" + std::to_string(reks_N) + ", M=" + std::to_string(reks_M) +
            "] is not installed. Installed [N,M] active spaces: " +
            reks::studio::supported_variants_str() + ".");
    }
    timer_on("REKS: make_catalog");
    // Merged (N,M) blob bound to its lowest manifold; the run binding is installed once the run set is known.
    catalog_ = reks::studio::make_catalog(reks_N, reks_M);
    timer_off("REKS: make_catalog");

    // Core = Ncore doubly occupied orbitals holding nelectron - N electrons (even for a singlet):
    // require Ncore >= 0 and Ncore + M <= nmo.
    Ncore_ = (nelectron_ - reks_N) / 2;
    if (Ncore_ < 0) {
        reks::report::input_error(
            "REKS active space mismatch: active electrons N=" + std::to_string(reks_N) +
            " exceed the molecule's " + std::to_string(nelectron_) +
            " total electrons (implied core orbitals = " + std::to_string(Ncore_) + " < 0).");
    }
    if (Ncore_ + reks_M > nmo_) {
        reks::report::input_error(
            "REKS active space does not fit the basis: core (" + std::to_string(Ncore_) +
            ") + active orbitals M=" + std::to_string(reks_M) + " = " +
            std::to_string(Ncore_ + reks_M) + " exceed the " + std::to_string(nmo_) +
            " molecular orbitals available. Enlarge the basis or shrink the active space.");
    }

    // Catalog view bound to blob sector b: sector-local config indices and names
    // resolve within that manifold's config block.
    auto blob_view = [this](int b) { return catalog_.bound_to(b); };

    const int n_avail = catalog_.n_sectors_available();
    // Available manifolds as "0, 2, 4" for error text.
    auto avail_2s_csv = [&]() {
        std::string avail;
        for (int b = 0; b < n_avail; ++b) {
            if (b) avail += ", ";
            avail += std::to_string(catalog_.sector_entry(b).spin2);
        }
        return avail;
    };
    // Ladder assignment as "block 0 -> 2S=0, ..." for error text.
    auto ladder_map_csv = [&](const char* unit) {
        std::string s;
        for (int b = 0; b < n_avail; ++b) {
            if (b) s += ", ";
            s += std::string(unit) + " " + std::to_string(b) + " -> 2S=" +
                 std::to_string(catalog_.sector_entry(b).spin2);
        }
        return s;
    };

    // Label/ladder resolution trace appended to resolution-failure error text.
    std::string label_trace;
    auto trace_label = [&](const char* unit, int i, int spin2, bool ladder) {
        label_trace += std::string("  ") + unit + " " + std::to_string(i) + " -> 2S=" +
                       std::to_string(spin2) + (ladder ? " (ladder)" : " (explicit)") + "\n";
    };
    auto with_trace = [&](const std::string& msg) {
        return label_trace.empty() ? msg : msg + "\n  Resolved so far:\n" + label_trace;
    };

    // Nesting depth of the first non-empty element, required of every other; -1 if all are empty.
    auto probe_form = [&](const char* optname, Data& src, int max_depth) {
        int form = -1, first_i = -1;
        for (int i = 0; i < static_cast<int>(src.size()); ++i) {
            Data& e = src[i];
            if (e.is_array() && e.size() == 0) continue;
            const int d = entry_depth(e);
            if (form < 0) {
                form = d;
                first_i = i;
            } else if (d != form) {
                reks::report::input_error(
                    std::string(optname) + ": mixed nesting depth: element " +
                    std::to_string(first_i) + " has depth " + std::to_string(form) +
                    " but element " + std::to_string(i) + " has depth " +
                    std::to_string(d) + ".");
            }
        }
        if (form > max_depth) {
            reks::report::input_error(
                std::string(optname) + ": element " + std::to_string(first_i) +
                " is nested too deeply (depth " + std::to_string(form) + ", maximum " +
                std::to_string(max_depth) + " for this option).");
        }
        return form;
    };

    // Flat list of 2S labels: >= 0, pairwise distinct; manifold existence is
    // checked at blob resolution.
    auto parse_2spin = [&](const char* optname, Data& src) {
        const int n = static_cast<int>(src.size());
        if (n == 0) {
            reks::report::input_error(std::string(optname) +
                ": set but empty; omit it or list one 2S value per block/group.");
        }
        std::vector<int> labels;
        std::set<int> seen;
        for (int i = 0; i < n; ++i) {
            const int t = static_cast<int>(src[i].to_double());
            if (t < 0) {
                reks::report::input_error(std::string(optname) +
                    ": entries are 2S >= 0; got " + std::to_string(t) + ".");
            }
            if (!seen.insert(t).second) {
                reks::report::input_error(std::string(optname) +
                    ": entries must be pairwise distinct; 2S=" + std::to_string(t) +
                    " is repeated.");
            }
            labels.push_back(t);
        }
        return labels;
    };
    // Blob sector whose manifold is 2S=t, fail-closed with the available listing.
    auto blob_of_2s = [&](const char* optname, int t) {
        for (int b = 0; b < n_avail; ++b)
            if (catalog_.sector_entry(b).spin2 == t) return b;
        reks::report::input_error(with_trace(
            std::string(optname) + ": active space [" + std::to_string(reks_N) + ", " +
            std::to_string(reks_M) + "] has no 2S=" + std::to_string(t) +
            " spin manifold (available 2S: " + avail_2s_csv() + ")."));
    };

    // Parse+validate one SA block's local-index/name list into global K indices.
    auto parse_block_sa = [&](Data& list, int b, const std::string& where) {
        const reks::studio::Catalog cv = blob_view(b);
        const int start = cv.sector_config_start();
        std::vector<int> Ks =
            parse_config_indices(list, cv, /*allow_exclude=*/false, "SA_REKS_CONFIGS", where);
        std::set<int> seen;
        for (int K : Ks) {
            if (!seen.insert(K).second) {
                reks::report::input_error(
                    "SA_REKS_CONFIGS: duplicate config (local) " + std::to_string(K - start) +
                    " in " + where + " (each config at most once per block).");
            }
        }
        return Ks;
    };
    // One SI cassette's list -> global K as a SET: a repeat keeps its first position (in the weighted SA
    // pool a repeat is an error).
    auto parse_cassette_si = [&](Data& list, int b, const std::string& where) {
        const reks::studio::Catalog cv = blob_view(b);
        const std::vector<int> parsed =
            parse_config_indices(list, cv, /*allow_exclude=*/true, "SI_REKS_CONFIGS", where);
        std::set<int> seen;
        std::vector<int> Ks;
        Ks.reserve(parsed.size());
        for (int K : parsed) {
            if (seen.insert(K).second) Ks.push_back(K);
        }
        return Ks;
    };

    // SA blocks bound to blob sectors in typed order; an empty ladder block skips its rung, an empty
    // labeled block still declares its manifold.
    struct SaBlockBind {
        int blob;
        std::vector<int> Ks;
    };
    std::vector<SaBlockBind> sa_blocks;
    bool sa_default_deferred = false;  // default pool lands on run sector 0 (resolved below)

    const bool have_sa       = alias_changed(options_, "SA_REKS_CONFIGS", "REKS_SA_CONFIGS");
    const bool have_sa_2spin = alias_changed(options_, "SA_REKS_2SPIN", "REKS_SA_2SPIN");
    std::vector<int> sa_labels;
    if (have_sa_2spin)
        sa_labels = parse_2spin("SA_REKS_2SPIN",
                                read_alias(options_, "SA_REKS_2SPIN", "REKS_SA_2SPIN"));

    if (have_sa) {
        Data& src = read_alias(options_, "SA_REKS_CONFIGS", "REKS_SA_CONFIGS");
        const int form = probe_form("SA_REKS_CONFIGS", src, 1);
        const bool flat = (form == 0);  // flat list of scalars/names = one block
        const int n_blocks = flat ? 1 : static_cast<int>(src.size());
        if (have_sa_2spin) {
            if (static_cast<int>(sa_labels.size()) != n_blocks) {
                std::string msg =
                    "SA_REKS_2SPIN: " + std::to_string(sa_labels.size()) + " label(s) for " +
                    std::to_string(n_blocks) + " typed SA block(s) in SA_REKS_CONFIGS; "
                    "lengths must match.";
                if (n_blocks == 0)
                    msg += " To request the default pool of manifold x, OMIT sa_reks_configs "
                           "and set sa_reks_2spin [x].";
                reks::report::input_error(msg);
            }
            for (int i = 0; i < n_blocks; ++i) {
                const int b = blob_of_2s("SA_REKS_2SPIN", sa_labels[i]);
                trace_label("SA block", i, sa_labels[i], false);
                Data& blk = flat ? src : src[i];
                sa_blocks.push_back({b, parse_block_sa(
                    blk, b, "SA block " + std::to_string(i) + " (2S=" +
                    std::to_string(sa_labels[i]) + ")")});
            }
        } else {
            if (n_blocks > n_avail) {
                reks::report::input_error(
                    "SA_REKS_CONFIGS: active space [" + std::to_string(reks_N) + ", " +
                    std::to_string(reks_M) + "] has " + std::to_string(n_avail) +
                    " spin manifolds, got " + std::to_string(n_blocks) + " SA blocks. "
                    "Ladder mode maps block i to the i-th available manifold (" +
                    ladder_map_csv("block") + "); use SA_REKS_2SPIN for explicit labels.");
            }
            for (int i = 0; i < n_blocks; ++i) {
                Data& blk = flat ? src : src[i];
                if (static_cast<int>(blk.size()) == 0) continue;  // ladder []: skip rung
                const int spin2 = catalog_.sector_entry(i).spin2;
                trace_label("SA block", i, spin2, true);
                sa_blocks.push_back({i, parse_block_sa(
                    blk, i, "SA block " + std::to_string(i) + " (2S=" +
                    std::to_string(spin2) + ")")});
            }
        }
    } else if (have_sa_2spin) {
        if (sa_labels.size() != 1) {
            reks::report::input_error(
                "SA_REKS_2SPIN: " + std::to_string(sa_labels.size()) + " labels but "
                "SA_REKS_CONFIGS is absent (one synthesized default block); a single "
                "label [x] puts the default SA pool on manifold x.");
        }
        const int b = blob_of_2s("SA_REKS_2SPIN", sa_labels[0]);
        trace_label("SA block", 0, sa_labels[0], false);
        sa_blocks.push_back(SaBlockBind{b, catalog_.sector_sa_default(b)});
    } else {
        sa_default_deferred = true;
    }

    // SI groups bound to blob sectors: 2-level = one implicit group, 3-level = one per outer entry; an
    // empty labeled group declares its manifold, ladder mode only through content.
    struct SiGroupBind {
        int blob;
        std::vector<std::vector<int>> cassettes;
    };
    std::vector<SiGroupBind> si_groups;
    bool si_default_deferred = false;  // default cassette on run sector 0 (resolved below)
    Data* si_implicit_src = nullptr;   // 2-level form; the group's sector is run sector 0

    const bool have_si       = alias_changed(options_, "SI_REKS_CONFIGS", "REKS_SI_CONFIGS");
    const bool have_si_2spin = alias_changed(options_, "SI_REKS_2SPIN", "REKS_SI_2SPIN");
    std::vector<int> si_labels;
    if (have_si_2spin)
        si_labels = parse_2spin("SI_REKS_2SPIN",
                                read_alias(options_, "SI_REKS_2SPIN", "REKS_SI_2SPIN"));

    // Per-cassette cap on reported SI states; entry i governs the cassette printed
    // as K=i. Validated once the final cassette list is known.
    const bool have_si_report =
        alias_changed(options_, "SI_REKS_REPORT_STATES", "REKS_SI_REPORT_STATES");
    std::vector<int> si_report;
    if (have_si_report) {
        Data& src = read_alias(options_, "SI_REKS_REPORT_STATES", "REKS_SI_REPORT_STATES");
        const int n_entries = static_cast<int>(src.size());
        si_report.reserve(n_entries);
        for (int i = 0; i < n_entries; ++i)
            si_report.push_back(static_cast<int>(src[i].to_double()));
    }

    // One 3-level group entry -> its non-empty cassettes.
    auto parse_si_group = [&](Data& grp, int b, int gi, int spin2) {
        std::vector<std::vector<int>> cassettes;
        for (int j = 0; j < static_cast<int>(grp.size()); ++j) {
            Data& cass = grp[j];
            if (!cass.is_array()) {
                reks::report::input_error(
                    "SI_REKS_CONFIGS: in the 3-level form, group " + std::to_string(gi) +
                    "'s entry must be a list of SI cassettes ([[...], ...]); each cassette "
                    "is itself a list of config indices or names.");
            }
            std::vector<int> Ks = parse_cassette_si(
                cass, b, "SI group " + std::to_string(gi) + " (2S=" +
                std::to_string(spin2) + ")");
            if (Ks.empty()) continue;
            cassettes.push_back(std::move(Ks));
        }
        return cassettes;
    };

    if (have_si) {
        Data& src = read_alias(options_, "SI_REKS_CONFIGS", "REKS_SI_CONFIGS");
        const int form = probe_form("SI_REKS_CONFIGS", src, 2);
        const int n_elems = static_cast<int>(src.size());
        if (form == 0) {
            reks::report::input_error(
                "SI_REKS_CONFIGS: scalar elements are not a cassette form; use a list of "
                "cassettes [[...], ...] (2-level) or cassette groups [[[...]], ...] "
                "(3-level).");
        }
        if (static_cast<int>(si_labels.size()) > 1 && form == 1) {
            reks::report::input_error(
                "SI_REKS_2SPIN has " + std::to_string(si_labels.size()) + " labels (one per "
                "SI group), which requires the 3-level SI_REKS_CONFIGS form [[[...]], ...] "
                "(one group of cassettes per label); got a 2-level list of cassettes.");
        }
        if (have_si_2spin) {
            if (form == 1) {
                // 2-level with one label: the implicit group lands on that manifold.
                const int b = blob_of_2s("SI_REKS_2SPIN", si_labels[0]);
                trace_label("SI group", 0, si_labels[0], false);
                SiGroupBind g{b, {}};
                for (int j = 0; j < n_elems; ++j) {
                    std::vector<int> Ks = parse_cassette_si(
                        src[j], b, "SI group 0 (2S=" + std::to_string(si_labels[0]) + ")");
                    if (Ks.empty()) continue;
                    g.cassettes.push_back(std::move(Ks));
                }
                si_groups.push_back(std::move(g));
            } else {
                // 3-level (or all-empty): one group per outer entry, one label per group.
                if (static_cast<int>(si_labels.size()) != n_elems) {
                    std::string msg =
                        "SI_REKS_2SPIN: " + std::to_string(si_labels.size()) +
                        " label(s) for " + std::to_string(n_elems) + " SI group(s) in "
                        "SI_REKS_CONFIGS; lengths must match.";
                    if (n_elems == 0)
                        msg += " To request the default SI cassette of manifold x, OMIT "
                               "si_reks_configs and set si_reks_2spin [x].";
                    reks::report::input_error(msg);
                }
                for (int i = 0; i < n_elems; ++i) {
                    const int b = blob_of_2s("SI_REKS_2SPIN", si_labels[i]);
                    trace_label("SI group", i, si_labels[i], false);
                    si_groups.push_back({b, parse_si_group(src[i], b, i, si_labels[i])});
                }
            }
        } else {
            if (form < 0) {
                // Present but all-empty: explicitly no SI anywhere; no default.
            } else if (form == 1) {
                si_implicit_src = &src;  // sector known once the run set is built
            } else {
                if (n_elems > n_avail) {
                    reks::report::input_error(
                        "SI_REKS_CONFIGS: active space [" + std::to_string(reks_N) + ", " +
                        std::to_string(reks_M) + "] has " + std::to_string(n_avail) +
                        " spin manifolds, got " + std::to_string(n_elems) + " SI groups. "
                        "Ladder mode maps group i to the i-th available manifold (" +
                        ladder_map_csv("group") + "); use SI_REKS_2SPIN for explicit labels.");
                }
                for (int i = 0; i < n_elems; ++i) {
                    Data& grp = src[i];
                    if (static_cast<int>(grp.size()) == 0) continue;  // ladder []: skip rung
                    const int spin2 = catalog_.sector_entry(i).spin2;
                    trace_label("SI group", i, spin2, true);
                    SiGroupBind g{i, parse_si_group(grp, i, i, spin2)};
                    if (g.cassettes.empty()) continue;  // content-free: declares nothing
                    si_groups.push_back(std::move(g));
                }
            }
        }
    } else if (have_si_2spin) {
        if (si_labels.size() != 1) {
            reks::report::input_error(
                "SI_REKS_2SPIN: " + std::to_string(si_labels.size()) + " labels but "
                "SI_REKS_CONFIGS is absent (one default cassette); a single label [x] "
                "puts the default SI cassette on manifold x.");
        }
        const int b = blob_of_2s("SI_REKS_2SPIN", si_labels[0]);
        trace_label("SI group", 0, si_labels[0], false);
        SiGroupBind g{b, {}};
        for (auto& sl : reks::studio::default_si_indices(blob_view(b))) {
            if (sl.empty()) continue;
            g.cassettes.push_back(std::move(sl));
        }
        si_groups.push_back(std::move(g));
    } else {
        si_default_deferred = true;
    }

    // Run set: manifolds declared by SA blocks and SI groups (empty -> the lowest one), ascending in 2S.
    std::set<int> run_blobs;
    for (const auto& blk : sa_blocks) run_blobs.insert(blk.blob);
    for (const auto& g : si_groups) run_blobs.insert(g.blob);
    if (run_blobs.empty()) run_blobs.insert(0);

    std::vector<int> declared_2S;  // run sector s -> 2S of its manifold
    const int n_run_sectors = static_cast<int>(run_blobs.size());
    run_to_blob_.assign(n_run_sectors, -1);
    blob_to_run_sector_.assign(n_avail, -1);
    {
        int s = 0;
        for (int b : run_blobs) {
            run_to_blob_[s] = b;
            blob_to_run_sector_[b] = s;
            declared_2S.push_back(catalog_.sector_entry(b).spin2);
            ++s;
        }
    }
    catalog_ = catalog_.bound_to(run_to_blob_[0]);

    // Deferred defaults land on run sector 0 (the run's lowest manifold).
    if (sa_default_deferred)
        sa_blocks.push_back(
            SaBlockBind{run_to_blob_[0], catalog_.sector_sa_default(run_to_blob_[0])});
    if (si_implicit_src) {
        const int b = run_to_blob_[0];
        SiGroupBind g{b, {}};
        Data& src = *si_implicit_src;
        for (int j = 0; j < static_cast<int>(src.size()); ++j) {
            std::vector<int> Ks = parse_cassette_si(
                src[j], b, "SI group 0 (2S=" + std::to_string(declared_2S[0]) + ")");
            if (Ks.empty()) continue;
            g.cassettes.push_back(std::move(Ks));
        }
        if (!g.cassettes.empty()) si_groups.push_back(std::move(g));
    }
    if (si_default_deferred) {
        SiGroupBind g{run_to_blob_[0], {}};
        for (auto& sl : reks::studio::default_si_indices(catalog_sector_view(0))) {
            if (sl.empty()) continue;
            g.cassettes.push_back(std::move(sl));
        }
        si_groups.push_back(std::move(g));
    }

    // SA pool in canonical order (blocks by ascending manifold); weight_perm maps canonical -> typed
    // positions for explicit SA_REKS_WEIGHTS.
    std::vector<int> sa_Ks;
    std::vector<int> sa_sector_counts(n_run_sectors, 0);
    std::vector<int> weight_perm;
    {
        const int nb = static_cast<int>(sa_blocks.size());
        std::vector<int> typed_off(nb, 0);
        for (int j = 1; j < nb; ++j)
            typed_off[j] = typed_off[j - 1] + static_cast<int>(sa_blocks[j - 1].Ks.size());
        std::vector<int> order(nb);
        for (int j = 0; j < nb; ++j) order[j] = j;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return sa_blocks[a].blob < sa_blocks[b].blob; });
        for (int j : order) {
            const auto& blk = sa_blocks[j];
            sa_sector_counts[blob_to_run_sector_[blk.blob]] = static_cast<int>(blk.Ks.size());
            for (int k = 0; k < static_cast<int>(blk.Ks.size()); ++k) {
                sa_Ks.push_back(blk.Ks[k]);
                weight_perm.push_back(typed_off[j] + k);
            }
        }
    }

    // SA extra determinants: one microstate [alpha(M) | beta(M)] each, in the SA ensemble only (not SI).
    const int n_active_orb = catalog_.data->n_active_orbitals;
    const int n_active_elec = catalog_.data->n_electrons;
    std::vector<reks::studio::Microstate> extra_microstates;
    if (alias_changed(options_, "SA_REKS_EXTRA", "REKS_SA_EXTRA")) {
        Data& src = read_alias(options_, "SA_REKS_EXTRA", "REKS_SA_EXTRA");
        const int n = static_cast<int>(src.size());
        extra_microstates.reserve(n);
        for (int i = 0; i < n; ++i) {
            Data& occ = src[i];
            if (static_cast<int>(occ.size()) != 2 * n_active_orb) {
                reks::report::input_error(
                    "SA_REKS_EXTRA[" + std::to_string(i) + "]: expected " +
                    std::to_string(2 * n_active_orb) +
                    " occupation entries [alpha(M)|beta(M)] for the M=" +
                    std::to_string(n_active_orb) + " active orbitals, got " +
                    std::to_string(occ.size()) + ".");
            }
            reks::studio::Microstate m{};
            int n_e = 0;
            for (int j = 0; j < 2 * n_active_orb; ++j) {
                const int v = static_cast<int>(occ[j].to_double());
                if (v != 0 && v != 1) {
                    reks::report::input_error(
                        "SA_REKS_EXTRA[" + std::to_string(i) + "][" + std::to_string(j) +
                        "] = " + std::to_string(v) + ": occupations must be 0 or 1.");
                }
                if (j < n_active_orb) m.alpha[j] = static_cast<int8_t>(v);
                else                  m.beta[j - n_active_orb] = static_cast<int8_t>(v);
                n_e += v;
            }
            if (n_e != n_active_elec) {
                reks::report::input_error(
                    "SA_REKS_EXTRA[" + std::to_string(i) + "]: determinant has " +
                    std::to_string(n_e) + " electrons but the active space holds N=" +
                    std::to_string(n_active_elec) + ".");
            }
            extra_microstates.push_back(m);
        }
    }
    const int n_extra = static_cast<int>(extra_microstates.size());

    // SA weights over the full pool: config weights first, then extra-determinant
    // weights. Length = n_configs_selected + n_extra; non-negative, sum to 1.
    const int n_pool = static_cast<int>(sa_Ks.size()) + n_extra;
    if (n_pool == 0) {
        reks::report::input_error(
            "The REKS SA ensemble is empty: SA_REKS_CONFIGS selected no configurations "
            "and SA_REKS_EXTRA supplied no determinants. The SA ensemble needs at least one.");
    }
    std::vector<double> sa_weights;
    std::vector<double> extra_weights;
    if (alias_changed(options_, "SA_REKS_WEIGHTS", "REKS_SA_WEIGHTS")) {
        Data& src = read_alias(options_, "SA_REKS_WEIGHTS", "REKS_SA_WEIGHTS");
        if (static_cast<int>(src.size()) != n_pool) {
            reks::report::input_error(
                "SA_REKS_WEIGHTS: expected " + std::to_string(n_pool) +
                " weights (SA_REKS_CONFIGS=" + std::to_string(sa_Ks.size()) +
                " + SA_REKS_EXTRA=" + std::to_string(n_extra) +
                "), got " + std::to_string(src.size()) + ".");
        }
        std::vector<double> all;
        all.reserve(n_pool);
        for (int i = 0; i < n_pool; ++i) all.push_back(src[i].to_double());
        double sum = 0.0;
        for (int i = 0; i < n_pool; ++i) {
            if (all[i] < 0.0) {
                reks::report::input_error(
                    "SA_REKS_WEIGHTS[" + std::to_string(i) +
                    "] = " + std::to_string(all[i]) + ": weights must be non-negative.");
            }
            sum += all[i];
        }
        if (std::abs(sum - 1.0) > 1.0e-8) {
            reks::report::input_error(
                "SA_REKS_WEIGHTS: weights sum to " + std::to_string(sum) +
                ", which differs from 1.0 by more than 1e-8.");
        }
        // Config weights arrive in typed block order; place each at its canonical
        // position. Extras stay at the tail.
        sa_weights.resize(sa_Ks.size());
        for (size_t i = 0; i < sa_Ks.size(); ++i) sa_weights[i] = all[weight_perm[i]];
        extra_weights.assign(all.begin() + sa_Ks.size(), all.end());
    } else {
        std::vector<double> all = reks::studio::uniform_weights(n_pool);
        sa_weights.assign(all.begin(), all.begin() + sa_Ks.size());
        extra_weights.assign(all.begin() + sa_Ks.size(), all.end());
    }

    // SI cassettes in canonical group order; si_lists[i] = global K of cassette i, si_sectors[i] = its run
    // sector.
    std::vector<std::vector<int>> si_lists;
    std::vector<int>              si_sectors;
    {
        const int ng = static_cast<int>(si_groups.size());
        std::vector<int> order(ng);
        for (int j = 0; j < ng; ++j) order[j] = j;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return si_groups[a].blob < si_groups[b].blob; });
        for (int j : order) {
            const int s = blob_to_run_sector_[si_groups[j].blob];
            for (auto& Ks : si_groups[j].cassettes) {
                si_lists.push_back(std::move(Ks));
                si_sectors.push_back(s);
            }
        }
    }

    // One report cap per SI cassette, each within its cassette's dimension.
    if (have_si_report) {
        if (si_report.empty()) {
            reks::report::input_error(
                "SI_REKS_REPORT_STATES: set but empty; omit it or list one value per SI cassette.");
        }
        if (si_lists.empty()) {
            reks::report::input_error(
                "SI_REKS_REPORT_STATES: set but the run has no SI cassettes.");
        }
        if (si_report.size() != si_lists.size()) {
            reks::report::input_error(
                "SI_REKS_REPORT_STATES: " + std::to_string(si_report.size()) + " value(s) for " +
                std::to_string(si_lists.size()) +
                " SI cassette(s); list exactly one per cassette.");
        }
        for (size_t i = 0; i < si_report.size(); ++i) {
            if (si_report[i] < 1) {
                reks::report::input_error(
                    "SI_REKS_REPORT_STATES: entry " + std::to_string(i) + " is " +
                    std::to_string(si_report[i]) + "; each value must be at least 1.");
            }
            const int dim = static_cast<int>(si_lists[i].size());
            if (si_report[i] > dim) {
                reks::report::input_error(
                    "SI_REKS_REPORT_STATES: entry " + std::to_string(i) + " is " +
                    std::to_string(si_report[i]) + " but SI cassette " + std::to_string(i) +
                    " has dimension " + std::to_string(dim) + ".");
            }
        }
    }

    if (!(options_.get_double("SI_REKS_OVERLAP_CUTOFF") > 0.0)) {
        reks::report::input_error(
            "SI_REKS_OVERLAP_CUTOFF is " + std::to_string(options_.get_double("SI_REKS_OVERLAP_CUTOFF")) +
            "; it must be positive (relative cutoff on the SI overlap eigenvalues, default 1e-10).", false);
    }

    use_trah_ = options_.get_bool("REKS_USE_TRAH");
    if (use_trah_) {
        trah_ = std::make_unique<reks::TrahDriver>();
        trah_->configure();
    }

    use_gvb_diis_ = options_.get_bool("REKS_GVB_DIIS");
    cosx_ = options_.get_str("SCF_TYPE").find("COSX") != std::string::npos;

    // Plain solver: no accelerator, plain F_reks diag + FON Newton.
    use_plain_ = !use_trah_ && !use_gvb_diis_;
    if (use_plain_ && reports(2)) {
        outfile->Printf("  REKS: plain solver (no GVB-DIIS/TRAH).\n");
    }
    if (use_gvb_diis_) {
        if (use_trah_) {
            if (reports(4))
                outfile->Printf("  GVB-DIIS hybrid: TRAH warmup, switching to GVB-DIIS at iteration %d.\n",
                                options_.get_int("REKS_GVB_DIIS_START"));
        } else {
            if (reports(4))
                outfile->Printf("  GVB-DIIS: pure mode (no TRAH warmup).\n");
        }
        reks::GvbDriverConfig gcfg;
        reks::DiisConfig& dcfg = gcfg.diis;
        dcfg.max_vectors       = options_.get_int("DIIS_MAX_VECS");
        dcfg.level_shift       = options_.get_double("REKS_GVB_LEVEL_SHIFT");
        dcfg.cond_angle_filter = options_.get_bool("REKS_DIIS_COND_ANGLE_FILTER");
        dcfg.angle_tol         = options_.get_double("REKS_DIIS_COND_ANGLE_CUTOFF");
        dcfg.stale_delta       = options_.get_double("REKS_DIIS_COND_STALE_DELTA");
        dcfg.kappa_bypass      = options_.get_double("REKS_DIIS_COND_KAPPA_BYPASS");
        dcfg.svd_rcond         = options_.get_double("REKS_DIIS_COND_SVD_RCOND");
        dcfg.coeff_norm_max    = options_.get_double("REKS_DIIS_COND_COEFF_NORM_MAX");
        dcfg.cond_tikhonov     = options_.get_bool("REKS_DIIS_COND_TIKHONOV");
        dcfg.tikhonov_scale    = options_.get_double("REKS_DIIS_COND_TIKHONOV_SCALE");
        dcfg.mon_verdict       = options_.get_bool("REKS_GVB_ROBUST_VERDICT");
        dcfg.verdict_rho       = options_.get_double("REKS_DIIS_MON_VERDICT_RHO");
        dcfg.verdict_rewind_streak = options_.get_int("REKS_DIIS_MON_VERDICT_STREAK");
        dcfg.mon_suppress_window = options_.get_int("REKS_DIIS_MON_SUPPRESS_WINDOW");
        dcfg.mon_nonmonotone_m = options_.get_int("REKS_DIIS_MON_NONMONOTONE_M");
        dcfg.guard_latch_grace = options_.get_int("REKS_DIIS_GUARD_LATCH_GRACE");
        dcfg.guard_band_B      = options_.get_double("REKS_DIIS_GUARD_BAND_B");
        dcfg.guard_landing_rtol = options_.get_double("REKS_DIIS_GUARD_LANDING_TOLERANCE");
        dcfg.mon_progress_eps  = options_.get_double("REKS_DIIS_MON_PROGRESS_EPS");
        dcfg.fon_branch_tol    = options_.get_double("REKS_DIIS_FON_BRANCH_CUTOFF");
        dcfg.restart_budget    = options_.get_int("REKS_DIIS_RESTART_BUDGET");
        dcfg.mon_cycle         = options_.get_bool("REKS_DIIS_MON_CYCLE");
        dcfg.cycle_window      = options_.get_int("REKS_DIIS_MON_CYCLE_WINDOW");
        dcfg.cycle_cap_streak  = options_.get_int("REKS_DIIS_MON_CYCLE_CAP_STREAK");
        dcfg.orb_base_map_damp = options_.get_double("REKS_DIIS_ORB_BASE_MAP_DAMP");
        gcfg.start_iteration      = options_.get_int("REKS_GVB_DIIS_START");
        gcfg.d_convergence        = options_.get_double("D_CONVERGENCE");
        gcfg.e_convergence        = options_.get_double("E_CONVERGENCE");
        gcfg.maxiter              = options_.get_int("MAXITER");
        gcfg.gate_streak          = options_.get_int("REKS_GVB_GATE_STREAK");
        if (gcfg.gate_streak < 1)
            reks::report::input_error("REKS_GVB_GATE_STREAK must be >= 1.", false);
        gcfg.orbital_formulation  = (options_.get_str("REKS_DIIS_FORMULATION") == "ORBITAL");
        gcfg.cfm_activation_gorb  = options_.get_double("REKS_DIIS_CFM_ACTIVATION_GORB");
        gvb_ = std::make_unique<reks::GvbDiisDriver>();
        gvb_->configure(gcfg);
    }

    active_mo_indices_.resize(catalog_.data->n_active_orbitals);
    for (int i = 0; i < catalog_.data->n_active_orbitals; ++i)
        active_mo_indices_[i] = Ncore_ + i;

    reks::studio::allocate_reel(reel_, catalog_, nsopi_, nso_, n_extra);

    si_cassettes_.clear();

    // SA cassette: the SCF ensemble pool with the sector-less extras at the tail.
    sa_cassette_ = reks::studio::make_cassette(
        catalog_, std::move(sa_Ks), std::move(sa_weights),
        reks::studio::CassetteRole::SA,
        std::move(extra_microstates), std::move(extra_weights));

    // Multi-sector partition: positions[s] = contiguous K_indices positions of sector s (empty if SA-empty);
    // single-sector keeps make_cassette's default.
    if (n_run_sectors > 1) {
        std::vector<std::vector<int>> positions(n_run_sectors);
        std::vector<int> ngen(n_run_sectors);
        int off = 0;
        for (int s = 0; s < n_run_sectors; ++s) {
            positions[s].reserve(sa_sector_counts[s]);
            for (int i = 0; i < sa_sector_counts[s]; ++i) positions[s].push_back(off + i);
            off += sa_sector_counts[s];
            ngen[s] = catalog_.sector_entry(run_to_blob_[s]).n_generations;
        }
        sa_cassette_.install_sector_partition(std::move(positions), ngen);
    }

    // One FON snapshot per run sector (layers span its generations), reset to closed shell (1,1).
    reel_.fon_state.assign(n_run_sectors, reks::studio::FONSnapshot{});
    for (int s = 0; s < n_run_sectors; ++s) {
        const int ngen_s = catalog_.sector_entry(run_to_blob_[s]).n_generations;
        for (int gen = 0; gen < ngen_s; ++gen)
            reel_.fon_state[s].layers[gen].assign(
                catalog_.data->n_geminals, reks::studio::GeminalFON{1.0, 1.0});
    }

    build_fon_blocks();

    // Per-config geminal weight dependencies (DiagDeps) and the unioned active
    // geminal set per FON generation.
    if (reports(4)) {
        const auto* gt = sa_cassette_.geminal_templates();
        auto print_geminal_set = [&](const std::vector<int>& gset, const auto* tpl) {
            outfile->Printf("{");
            for (size_t i = 0; i < gset.size(); ++i) {
                const int g = gset[i];
                outfile->Printf("g%d(%c,%c)%s", g, 'a' + tpl[g].orbitals[0],
                                'a' + tpl[g].orbitals[1], i + 1 < gset.size() ? ", " : "");
            }
            outfile->Printf("}");
        };
        outfile->Printf("  [POOL] SA pool: %zu configs -> FON dependency union\n",
                        sa_cassette_.K_indices.size());
        for (int K : sa_cassette_.K_indices) {
            const auto& d = sa_cassette_.diag_deps(K);
            const char* nm = sa_cassette_.si_config_name(K);
            outfile->Printf("  [POOL]   K=%-2d %-12s reads", K, nm ? nm : "?");
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
                const auto& reads = d.geminal_reads[gen];
                outfile->Printf(" %s={", reks::studio::fon_channel_name(gen, sa_sector()).c_str());
                for (size_t j = 0; j < reads.size(); ++j)
                    outfile->Printf("%d%s", reads[j], j + 1 < reads.size() ? "," : "");
                outfile->Printf("}");
            }
            outfile->Printf("\n");
        }
        for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
            const auto& act = sa_cassette_.geminals_active_gen(gen);
            outfile->Printf("  [POOL]   => %s-FON (%zu) active = ",
                            reks::studio::fon_channel_name(gen, sa_sector()).c_str(), act.size());
            print_geminal_set(act, gt);
            outfile->Printf("\n");
        }
    }

    si_cassettes_.reserve(si_lists.size());
    for (size_t i = 0; i < si_lists.size(); ++i) {
        const int s = si_sectors[i];
        std::vector<int> sl_copy(si_lists[i]);
        reks::studio::Cassette c = reks::studio::make_cassette(
            catalog_sector_view(s), std::move(sl_copy), {},
            reks::studio::CassetteRole::SI);
        c.run_sector = s;
        if (!si_report.empty()) c.report_states = si_report[i];
        si_cassettes_.push_back(std::move(c));
    }

    call_sheet_ = reks::studio::build_call_sheet(sa_cassette_, si_cassettes_);

    if (reports(4)) {
        if (n_run_sectors == 1) {
            outfile->Printf("  [REKS_INIT] N=%d M=%d S=%d | "
                            "sa_Ks.size=%zu si_cassettes=%zu G=%d\n",
                            reks_N, reks_M, declared_2S[0],
                            sa_cassette_.K_indices.size(), si_cassettes_.size(),
                            sa_cassette_.n_orbitals_per_geminal());
        } else {
            std::string sec2S;
            for (int s = 0; s < n_run_sectors; ++s)
                sec2S += (s ? "," : "") + std::to_string(declared_2S[s]);
            outfile->Printf("  [REKS_INIT] N=%d M=%d sectors=%d 2S=[%s] | "
                            "sa_Ks.size=%zu si_cassettes=%zu G=%d\n",
                            reks_N, reks_M, n_run_sectors, sec2S.c_str(),
                            sa_cassette_.K_indices.size(), si_cassettes_.size(),
                            sa_cassette_.n_orbitals_per_geminal());
        }
    }

    if (reports(4))
        reks::report::print_memory_footprint("after reel allocation", report_frame(), *this);

    if (functional_->needs_xc()) {
        Vb_ = SharedMatrix(factory_->create_matrix("V beta"));
    }

    reks::report::print_studio_banner();
    reks::report::print_studio_setup(sa_cassette_, declared_2S[0], Ncore_, active_mo_indices_);
    if (n_run_sectors > 1) {
        std::vector<int> si_counts(n_run_sectors, 0);
        for (int s : si_sectors) ++si_counts[s];
        std::vector<int> ngen(n_run_sectors);
        for (int s = 0; s < n_run_sectors; ++s)
            ngen[s] = catalog_.sector_entry(run_to_blob_[s]).n_generations;
        reks::report::print_sector_map(declared_2S, sa_sector_counts, si_counts, ngen);
    }
    const bool show_triplet =
        std::any_of(declared_2S.begin(), declared_2S.end(), [](int s2) { return s2 > 0; });
    reks::report::print_pairing_schemes(sa_cassette_, show_triplet);
    reks::report::print_sa_pool(sa_cassette_);
    reks::report::print_si_pool(si_cassettes_, sa_cassette_, declared_2S);

    if (reports(4)) {
        outfile->Printf("    Run diagnostics\n");
        outfile->Printf("    ---------------\n");
        if (use_plain_) {
            outfile->Printf("      FON solver:          plain F_reks + FON Newton (no accelerator)\n");
        } else if (use_gvb_diis_ && !use_trah_) {
            outfile->Printf("      FON solver:          GVB-DIIS (no TRAH warmup)\n");
        } else if (use_gvb_diis_) {
            outfile->Printf("      FON solver:          TRAH warmup + GVB-DIIS\n");
        } else {
            outfile->Printf("      FON solver:          TRAH (combined orbital+FON)\n");
        }
        outfile->Printf("      Interpolation delta: %.4f\n", reks::delta_interp);
        outfile->Printf("      Initial FONs:       ");
        for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
            int g = sa_cassette_.geminals_active_gen(0)[u];
            double f0 = sa_fons().layers[0][g].p;
            outfile->Printf("  unit %d: (%.6f, %.6f)", u, f0, 2.0 - f0);
        }
        outfile->Printf("\n");
        outfile->Printf("      Report level:        %d\n",
                        static_cast<int>(reks::report::level()));
        outfile->Printf("      f_interp sanity:     f(0)=%.4f f(0.25)=%.4f f(0.5)=%.4f f(1)=%.4f\n",
                        reks::f_interp(0.0), reks::f_interp(0.25), reks::f_interp(0.5), reks::f_interp(1.0));
        outfile->Printf("\n");
    }

}

// Occupied block [0, nalphapi_[0]); its first Ncore_ columns are in D_core, the rest are the active
// pairs (none -> no update).
void REKS::form_Da_from_core() {
    const int nso  = Ca_->rowspi()[0];
    const int nmo  = Ca_->colspi()[0];
    const int nocc = nalphapi_[0];
    double**  Cp   = Ca_->pointer(0);

    Da_->copy(reel_.D_core);
    double** Dp = Da_->pointer(0);
    // D_core + C_occ C_occ^T over the active-occupied columns, one rank-k update.
    if (nocc > Ncore_)
        C_DGEMM('N', 'T', nso, nso, nocc - Ncore_, 1.0, &Cp[0][Ncore_], nmo, &Cp[0][Ncore_], nmo,
                1.0, Dp[0], nso);
}

void REKS::build_base_densities() {
    timer_on("REKS: build_base_densities");
    reks::scf::build_base_densities(sa_cassette_, call_sheet_.sa_referenced_bitmask,
                                         reel_, Ncore_, active_mo_indices_, Ca_, H_);
    timer_off("REKS: build_base_densities");
}

int REKS::base_idx(int L, bool alpha) const {
    return reks::studio::base_density_index(
        sa_cassette_.microstate(L), sa_cassette_.n_active_orbitals(), /*beta=*/!alpha);
}

void REKS::assemble_F_reks_MO() {
    const int nmo = Ca_->colspi()[0];
    // Sized at first use: nmo is set only after form_Shalf.
    if (!reel_.F_reks_MO || reel_.F_reks_MO->coldim(0) != nmo) {
        reel_.F_reks_MO = std::make_shared<Matrix>("F_REKS_MO", nmo, nmo);
    }

    reks::scf::assemble_F_reks_MO(
        reel_, sa_cassette_,
        reel_.C_L, active_mo_indices_, Ncore_, nmo, reel_.F_reks_MO);
}

void REKS::assemble_F_reks() {
    assemble_F_reks_MO();

    // Pure GGA: virt-virt block of F_reks <- closed-shell RKS Fock (negative microstate weights give
    // unphysical virtual energies); only the virtual columns are transformed.
    if (functional_ && functional_->needs_xc() && functional_->is_gga() && functional_->x_alpha() < 0.01) {
        const int nmo_c = Ca_->colspi()[0];
        const int first_virt = active_mo_indices_.back() + 1;
        const int n_virt = nmo_c - first_virt;
        if (n_virt > 0) {
            const int nso = Ca_->rowspi()[0];
            // G_ is the closed-shell Fock contribution (2J - x_alpha K - .. wK + Vxc).
            auto F_rhf_AO = std::make_shared<Matrix>("F_rhf_AO", nsopi_, nsopi_);
            F_rhf_AO->copy(H_);
            F_rhf_AO->add(G_);

            // Cv is the trailing column block of Ca_, used in place with lda = nmo; the
            // result (Cv^T F_rhf_AO) Cv is written into F_reks_MO with ldc = nmo.
            double** Cap = Ca_->pointer(0);
            double* Cv0 = &Cap[0][first_virt];
            std::vector<double> Tvv(static_cast<size_t>(n_virt) * nso);
            double** Freks = reel_.F_reks_MO->pointer(0);
            C_DGEMM('T', 'N', n_virt, nso, nso, 1.0, Cv0, nmo_c, F_rhf_AO->pointer(0)[0], nso, 0.0,
                    Tvv.data(), nso);
            C_DGEMM('N', 'N', n_virt, n_virt, nso, 1.0, Tvv.data(), nso, Cv0, nmo_c, 0.0,
                    &Freks[first_virt][first_virt], nmo_c);
        }
    }
}

void REKS::compute_weighting_factors() {
    timer_on("REKS: compute_weighting_factors");
    reks::scf::compute_weights(reel_, sa_cassette_, reel_.C_L);

    if (reports(5)) {
        outfile->Printf("\n  === REKS Weighting Factors C_L ===\n");
        for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
            int g = sa_cassette_.geminals_active_gen(0)[u];
            double f0 = sa_fons().layers[0][g].p;
            double f1 = 2.0 - f0;
            outfile->Printf("  unit %d: fon_p = %.6f, fon_q = %.6f, f = %.6f\n", u, f0, f1, reks::f_interp(f0 * f1));
        }
        outfile->Printf("  SA config weights:");
        for (int K : sa_cassette_.K_indices) {
            outfile->Printf(" w[%d]=%.4f", K, sa_cassette_.w_configs[K]);
        }
        outfile->Printf("\n");

        int n_micro = sa_cassette_.n_total_microstates();
        double sum = 0.0;
        for (int L = 0; L < n_micro; ++L) {
            outfile->Printf("  C_L[%d] = %12.8f\n", L, reel_.C_L[L]);
            sum += reel_.C_L[L];
        }
        outfile->Printf("  Sum C_L = %12.8f (should be ~1.0)\n", sum);
    }
    timer_off("REKS: compute_weighting_factors");
}

void REKS::build_sa_focks() {
    timer_on("REKS: build_sa_focks");

    // Refreshes C_L from the current FON.
    reks::scf::compute_weights(reel_, sa_cassette_, reel_.C_L);

    reks::scf::build_jk_cache(
        sa_cassette_, reel_, Ncore_, active_mo_indices_, Ca_, C_occ_cache_,
        *jk_, *functional_, J_, K_, wK_, G_, /*need_full_diag=*/use_gvb_diis_);

    // Closed-shell RV grid sweep for pure-GGA virt-virt XC (x_alpha < 0.01):
    // G_ += Va_[Da_].
    const bool rv_G_gga = functional_->is_gga() && functional_->x_alpha() < 0.01;
    if (rv_G_gga && functional_->needs_xc() && rv_potential_) {
        timer_on("REKS: xc_compute (sa rv)");
        rv_potential_->set_D({Da_});
        rv_potential_->compute_V({Va_});
        G_->add(Va_);
        timer_off("REKS: xc_compute (sa rv)");
    }

    // Refresh reel_.Ca_fock_transform to Ca_, reallocating if the shape changed.
    if (!reel_.Ca_fock_transform ||
        reel_.Ca_fock_transform->rowspi() != Ca_->rowspi() ||
        reel_.Ca_fock_transform->colspi() != Ca_->colspi()) {
        reel_.Ca_fock_transform = Ca_->clone();
    } else {
        reel_.Ca_fock_transform->copy(Ca_);
    }

    // need_diag gates the per-L V_xc MO diagonal and its full-MO projection.
    reks::scf::build_microstate_xc(
        reel_, sa_cassette_, sa_cassette_.microstates(), Ncore_, *functional_,
        ms_xc_.get(), reel_.Ca_fock_transform, reel_.C_L,
        /*need_diag=*/use_gvb_diis_,
        /*mo_lo=*/0, /*mo_width=*/reel_.Ca_fock_transform->colspi()[0],
        /*want_v_acc=*/true);

    // Copies the microstate sweep's grid-electron count into RV's quadrature_values,
    // since the RV branch above may have skipped its own compute_V pass.
    if (ms_xc_ && rv_potential_) {
        if (const auto nq = ms_xc_->grid_electrons()) {
            auto& rvq = rv_potential_->quadrature_values();
            rvq["RHO_A"] = *nq;
            rvq["RHO_B"] = *nq;
        }
    }

    // D_core = C_core C_core^T of the Ca_ snapshot held in Ca_fock_transform.
    reks::scf::build_sa_focks_MO(
        reel_, sa_cassette_, sa_cassette_.microstates(), Ncore_, *functional_, H_,
        reel_.D_core, /*need_full_diag=*/use_gvb_diis_,
        /*mo_lo=*/0, /*mo_width=*/reel_.Ca_fock_transform->colspi()[0]);

    timer_off("REKS: build_sa_focks");
}

// Post-SCF SI: SI-only densities/Focks, per-sector SI FON, then H_si c = E S_si c per cassette.
void REKS::compute_si() {
    if (si_computed_) return;
    timer_on("REKS: compute_si");

    si_results_.assign(si_cassettes_.size(), reks::SIResult{});
    const bool needs_xc = functional_->needs_xc();

    // SI-only patterns, in the Fock-snapshot frame of the frozen E_L and Focks.
    if (!call_sheet_.si_only_bitmask.empty()) {
        reks::scf::build_base_densities(sa_cassette_, call_sheet_.si_only_bitmask,
                                        reel_, Ncore_, active_mo_indices_,
                                        reel_.Ca_fock_transform ? reel_.Ca_fock_transform : Ca_,
                                        H_);
    }

    // Per-L Focks outside the SA cassette's L set, from the SCF JK cache and Ca_fock_transform.
    reks::scf::fill_missing_si_focks(
        reel_, sa_cassette_, si_cassettes_, Ncore_, *functional_,
        ms_xc_.get(), reel_.Ca_fock_transform, H_);

    // E_L is FON-independent and L-keyed; filled from real cassettes only (union pairs have no per-L Fock).
    for (const reks::studio::Cassette& cassette : si_cassettes_) {
        reks::scf::build_si_microstate_focks(
            reel_, cassette, sa_cassette_,
            nuclearrep_ + energies_["-D"], needs_xc, reel_.E_L);
    }

    // Per-sector SI FON: a sector's SI cassettes union into one cassette relaxed once (primary K = sector
    // min, order-invariant); sectors decouple exactly at frozen orbitals.
    {
        timer_on("REKS: fon_relax");
        const int n_run_sectors = static_cast<int>(reel_.fon_state.size());
        for (int s = 0; s < n_run_sectors; ++s) {
            std::vector<int> union_K;
            for (const reks::studio::Cassette& cassette : si_cassettes_)
                if (cassette.sector() == s)
                    union_K.insert(union_K.end(),
                                   cassette.K_indices.begin(), cassette.K_indices.end());
            if (union_K.empty()) continue;
            std::sort(union_K.begin(), union_K.end());
            union_K.erase(std::unique(union_K.begin(), union_K.end()), union_K.end());

            reks::studio::Cassette union_cassette = reks::studio::make_cassette(
                catalog_sector_view(s), std::move(union_K), {},
                reks::studio::CassetteRole::SI);
            union_cassette.run_sector = s;

            // Seed from the sector's SA-converged FON, relax the sector's SI-only
            // geminals, write back into its snapshot; SA geminals keep their value.
            auto post = reks::fon_relax::optimize(
                union_cassette, sa_cassette_, s, reel_.fon_state[s], reel_.E_L);
            reel_.fon_state[s].layers = std::move(post.layers);
        }
        timer_off("REKS: fon_relax");
    }

    // Sorted union of call_sheet_.eri_pairs and, if sa_ran and n_act >= 2, every active pair p < q, in the
    // k <= l order compute_active_mo_eri needs; empty -> no JK pass.
    {
        std::set<std::pair<int,int>> pair_set(call_sheet_.eri_pairs.begin(),
                                              call_sheet_.eri_pairs.end());
        const bool sa_ran = reks::scf::sa_energies_ready(sa_cassette_, reel_.E_L);
        const int n_act = static_cast<int>(active_mo_indices_.size());
        if (sa_ran && n_act >= 2)
            for (int p = 0; p < n_act; ++p)
                for (int q = p + 1; q < n_act; ++q) pair_set.emplace(p, q);
        active_eri_pool_ = compute_active_mo_eri(
            active_mo_indices_,
            std::vector<std::pair<int,int>>(pair_set.begin(), pair_set.end()));
    }

    // SI H/S built from the catalog's coupling tables via the cassette's pair rows.
    for (size_t k = 0; k < si_cassettes_.size(); ++k) {
        const reks::studio::Cassette& cassette = si_cassettes_[k];

        std::vector<double> H_si;
        reks::BlockedMatrix S_si;
        int n_out = 0;
        {
            timer_on("REKS: compute_SI_energies");
            {
                // Per-cassette selector_packs, released with this scope; serial, timed apart from
                // build_hamiltonian.
                timer_on("REKS: si_selector_demand");
                const reks::studio::SelectorDemand demand =
                    reks::studio::selector_demand(cassette);
                timer_off("REKS: si_selector_demand");
                timer_on("REKS: si_selector_packs");
                const reks::studio::SelectorPacks selector_packs =
                    reks::studio::decode_selector_packs(
                        catalog_sector_view(cassette.sector()), demand);
                timer_off("REKS: si_selector_packs");
                timer_on("REKS: si_build_hamiltonian");
                reks::si::build_hamiltonian(
                    reel_, cassette, selector_packs, reel_.E_L,
                    reel_.lagrangians, functional_->x_alpha(), H_si, n_out,
                    &active_eri_pool_);
                timer_off("REKS: si_build_hamiltonian");
            }
            timer_on("REKS: si_build_overlap");
            reks::si::build_overlap(reel_, cassette, S_si);
            timer_off("REKS: si_build_overlap");
            si_results_[k] = reks::si::diagonalize(
                cassette, std::move(H_si), std::move(S_si), n_out);
            timer_off("REKS: compute_SI_energies");
        }
    }

    si_computed_ = true;
    timer_off("REKS: compute_si");
}

reks::ActiveEriTiles REKS::compute_active_mo_eri(
    const std::vector<int>& active_mo_indices,
    const std::vector<std::pair<int,int>>& rank1_pairs) {
    const int N = static_cast<int>(active_mo_indices.size());
    if (N == 0 || rank1_pairs.empty()) {
        return reks::ActiveEriTiles(N, {});
    }

    reks::ActiveEriTiles result(N, rank1_pairs);

    // Active MO ERI per pair k <= l, one batched jk_->compute(): R(u) = C_act^T J(u) C_act,
    // tile(k,k) = R(c_k), tile(k,l) = (R(u_+) - R(u_-))/2 with u_+- = (c_l +- c_k)/sqrt(2).
    timer_on("REKS: compute_active_mo_eri");

    const int nso = Ca_->rowspi()[0];

    // Same MOs used to build F_MO and Lagrangians.
    auto C_source = reel_.Ca_fock_transform ? reel_.Ca_fock_transform : Ca_;
    auto C_active = std::make_shared<Matrix>("C_active", nso, N);
    for (int i = 0; i < N; ++i) {
        const int mo = active_mo_indices[i];
        for (int mu = 0; mu < nso; ++mu) {
            C_active->set(mu, i, C_source->get(mu, mo));
        }
    }
    double** Ca = C_active->pointer();

    // slot_meta tags each column's role (diag/plus/minus) for the tile-fill step below.
    enum SlotKind { Diag, Plus, Minus };
    struct SlotMeta { SlotKind kind; int k, l; };

    std::vector<SharedMatrix> col_alpha;
    col_alpha.reserve(rank1_pairs.size() * 2);
    std::vector<SlotMeta> slot_meta;
    slot_meta.reserve(rank1_pairs.size() * 2);

    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    for (const auto& kv : rank1_pairs) {
        const int k = kv.first;
        const int l = kv.second;
        if (k == l) {
            auto col = std::make_shared<Matrix>("u_diag", nso, 1);
            double* up = col->pointer()[0];
            for (int mu = 0; mu < nso; ++mu) up[mu] = Ca[mu][k];
            col_alpha.push_back(col);
            slot_meta.push_back({Diag, k, l});
        } else {
            auto cp = std::make_shared<Matrix>("u_plus",  nso, 1);
            auto cm = std::make_shared<Matrix>("u_minus", nso, 1);
            double* p_plus  = cp->pointer()[0];
            double* p_minus = cm->pointer()[0];
            for (int mu = 0; mu < nso; ++mu) {
                p_plus [mu] = inv_sqrt2 * (Ca[mu][k] + Ca[mu][l]);
                p_minus[mu] = inv_sqrt2 * (Ca[mu][l] - Ca[mu][k]);
            }
            col_alpha.push_back(cp);
            slot_meta.push_back({Plus,  k, l});
            col_alpha.push_back(cm);
            slot_meta.push_back({Minus, k, l});
        }
    }

    // Push to C_left only -> JK runs lr_symmetric_=true.
    bool prev_do_K = functional_->is_x_hybrid() || functional_->is_x_lrc();

    auto& C_left  = jk_->C_left();
    auto& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();
    for (const auto& cm : col_alpha) C_left.push_back(cm);

    jk_->set_do_J(true);
    jk_->set_do_K(false);
    timer_on("REKS: jk_compute (active_eri)");
    jk_->compute();
    timer_off("REKS: jk_compute (active_eri)");

    const std::vector<SharedMatrix>& Jvec = jk_->J();

    auto Tmp = std::make_shared<Matrix>("Tmp", nso, N);
    auto Tij = std::make_shared<Matrix>("Tij", N,   N);

    auto compute_T_alpha = [&](size_t alpha, double** Rij) {
        Tmp->gemm(false, false, 1.0, Jvec[alpha], C_active, 0.0);
        Tij->gemm(true,  false, 1.0, C_active,    Tmp,      0.0);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                Rij[i][j] = Tij->get(i, j);
    };

    auto Rplus_buf  = std::make_shared<Matrix>("Rij+", N, N);
    auto Rminus_buf = std::make_shared<Matrix>("Rij-", N, N);
    double** Rp = Rplus_buf->pointer();
    double** Rm = Rminus_buf->pointer();

    // slot_meta order: each off-diag pair emits Plus then Minus contiguously.
    for (size_t alpha = 0; alpha < slot_meta.size(); ++alpha) {
        const SlotMeta& sm = slot_meta[alpha];
        if (sm.kind == Diag) {
            compute_T_alpha(alpha, Rp);
            double* tile = result.tile_data(sm.k, sm.l);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    tile[static_cast<size_t>(i) * N + j] = Rp[i][j];
        } else if (sm.kind == Plus) {
            compute_T_alpha(alpha, Rp);
            // Defer write until paired Minus is available (next slot).
        } else {  // Minus
            compute_T_alpha(alpha, Rm);
            double* tile = result.tile_data(sm.k, sm.l);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    tile[static_cast<size_t>(i) * N + j] = 0.5 * (Rp[i][j] - Rm[i][j]);
        }
    }

    C_left.clear();
    C_right.clear();
    jk_->set_do_K(prev_do_K);
    timer_off("REKS: compute_active_mo_eri");

    return result;
}

void REKS::compute_sa_energies() {
    timer_on("REKS: compute_sa_energies");

    energies_["-D"] = scalar_variable("-D Energy");
    double E_const = nuclearrep_ + energies_["-D"];
    bool needs_xc = functional_->needs_xc();
    double alpha = functional_->x_alpha();

    if (reports(5)) {
        outfile->Printf("\n  Microstate energy computation:\n");
        outfile->Printf("  E_const = %.10f, needs_xc = %d, alpha = %.4f\n", E_const, needs_xc, alpha);
    }

    // E_L = tr(Da H) + tr(Db H) + 1/2 [tr(Da Fa) + tr(Db Fb) - E_1e] + E_const + (E_xc - 1/2 tr(D V_xc));
    // base-density scalars are keyed by alpha_idx/beta_idx, the per-L terms by L.
    for (int L : sa_cassette_.microstates()) {
        const int alpha_idx = base_idx(L, /*alpha=*/true);
        const int beta_idx  = base_idx(L, /*alpha=*/false);
        const auto e = reks::scf::microstate_energy(
            reel_, L, alpha_idx, beta_idx, E_const, needs_xc);
        reel_.E_L[L] = e.total;

        if (needs_xc && reports(5)) {
            const double tr_D_Vxc = reel_.tr_D_Vxc_a_L[L] + reel_.tr_D_Vxc_b_L[L];
            outfile->Printf("  L=%d: E_1e=%.10f, E_2e=%.10f, E_const=%.10f\n", L, e.E_1e, e.E_2e, E_const);
            outfile->Printf("        E_L(before XC)=%.10f\n", e.E_1e + e.E_2e + E_const);
            outfile->Printf("        E_xc=%.10f, tr(D*Vxc)=%.10f, 0.5*tr=%.10f\n", reel_.E_xc_L[L], tr_D_Vxc,
                            0.5 * tr_D_Vxc);
            outfile->Printf("        xc_correction=%.10f\n", e.xc_correction);
            outfile->Printf("        E_L(after XC)=%.10f\n", reel_.E_L[L]);
        }
    }

    if (reports(5)) {
        outfile->Printf("  Microstate energies:");
        for (int L : sa_cassette_.microstates()) {
            outfile->Printf(" E[%s]=%.10f",
                            reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.E_L[L]);
        }
        outfile->Printf("\n");
        outfile->Printf("  Microstate weights: ");
        for (int L : sa_cassette_.microstates()) {
            outfile->Printf(" C[%s]=%.10f",
                            reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.C_L[L]);
        }
        outfile->Printf("\n");
    }

    if (reports(4)) {
        outfile->Printf("  [E_MICRO] iter=%d", iteration_);
        for (int L : sa_cassette_.microstates()) {
            outfile->Printf(" E[%s]=%.10f",
                            reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.E_L[L]);
        }
        outfile->Printf("\n");
        outfile->Printf("  [C_MICRO] iter=%d", iteration_);
        for (int L : sa_cassette_.microstates()) {
            outfile->Printf(" C[%s]=%.10f",
                            reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.C_L[L]);
        }
        outfile->Printf("\n");
        // (fon_p, fon_q) per generation.
        outfile->Printf("  [FON] iter=%d", iteration_);
        for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
            const std::string chan = reks::studio::fon_channel_name(gen, sa_sector());
            const auto& gem = sa_cassette_.geminals_active_gen(gen);
            for (size_t p = 0; p < gem.size(); ++p) {
                const auto& fon = sa_fons().layers[gen][gem[p]];
                outfile->Printf(" %s-FON%zu=(%.8f, %.8f)", chan.c_str(), p, fon.p, fon.q);
            }
        }
        outfile->Printf("\n");
    }
    timer_off("REKS: compute_sa_energies");
}

void REKS::form_D() {
    timer_on("REKS: form_D");
    C_occ_cache_ = Ca_subset("SO", "OCC");
    build_base_densities();
    form_Da_from_core();

    if (reports(5)) {
        outfile->Printf("\n  === Active Orbitals (after form_D) ===\n");

        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];

        int nso = nsopi_[0];
        double** Cp = Ca_->pointer(0);
        double** Sp = S_->pointer(0);
        double overlap_rs = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                overlap_rs += Cp[mu][active_r] * Sp[mu][nu] * Cp[nu][active_s];
            }
        }
        outfile->Printf("  <r|s> = %12.8f (should be ~0)\n", overlap_rs);

        double norm_r = 0.0, norm_s = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                norm_r += Cp[mu][active_r] * Sp[mu][nu] * Cp[nu][active_r];
                norm_s += Cp[mu][active_s] * Sp[mu][nu] * Cp[nu][active_s];
            }
        }
        outfile->Printf("  <r|r> = %12.8f, <s|s> = %12.8f (should be 1.0)\n", norm_r, norm_s);
    }
    timer_off("REKS: form_D");
}

void REKS::evaluate_frame() {
    form_D();

    // SA microstate Focks in a single JK batch.
    build_sa_focks();

    compute_sa_energies();
}

void REKS::form_G() {
    timer_on("REKS: form_G");
    // One-shot: a stand-on-best reported pair lives exactly one iteration.
    if (gvb_) gvb_->disarm_report_override();
    // SAD guess iteration: inline JK build.
    if (sad_ && iteration_ <= 0) {
        auto C_occ_sad = Ca_subset("SO", "OCC");
        jk_->C_left().clear();
        jk_->C_left().push_back(C_occ_sad);
        timer_on("REKS: jk_compute (sad)");
        jk_->compute();
        timer_off("REKS: jk_compute (sad)");
        J_->copy(jk_->J()[0]);
        G_->copy(J_);
        G_->scale(2.0);
        if (functional_->is_x_hybrid()) {
            K_->copy(jk_->K()[0]);
            G_->axpy(-functional_->x_alpha(), K_);
        }
        if (functional_->is_x_lrc()) {
            wK_->copy(jk_->wK()[0]);
            G_->axpy(-functional_->x_beta(), wK_);
        }
        if (functional_->needs_xc() && rv_potential_) {
            rv_potential_->set_D({Da_});
            timer_on("REKS: xc_compute (sad)");
            rv_potential_->compute_V({Va_});
            timer_off("REKS: xc_compute (sad)");
            G_->add(Va_);
        }
        timer_off("REKS: form_G");
        return;
    }

    // Stepped state of the previous form_C becomes the live state.
    if (staged_) {
        Ca_->copy(stage_Ca_);
        epsilon_a_->copy(*stage_eps_);
        std::swap(reel_.fon_state, stage_fons_);
        staged_ = false;
    }

    evaluate_frame();

    if (reports(4)) {
        reks::report::print_sa_diagnostics(report_frame(), iteration_, trace_);
    }

    // FON relaxation: TRAH already moved it in combined_step -> refresh weights only;
    // otherwise the FON Newton step relaxes it.
    const bool trah_owns_fon = use_trah_ && !(gvb_ && gvb_->active());
    if (trah_owns_fon) {
        compute_weighting_factors();
    } else {
        reks::FonStepContext ctx{iteration_, active_mo_indices_, sa_cassette_, reel_, fon_blocks_,
                                 reks_fon_lower_, Ca_, epsilon_a_};
        fon_step_->relax(ctx);
        compute_weighting_factors();
        if (ctx.swapped && gvb_) gvb_->note_swaps();
    }

    if (use_gvb_diis_ && reports(4)) {
        const int n_sec = sa_cassette_.n_sectors();
        // Gen-0 (n) FON, gradient, and Hessian diagonal per sector; the sector token
        // is empty in a single-sector run.
        for (int sec = 0; sec < n_sec; ++sec) {
            const int np = static_cast<int>(sa_cassette_.geminals_active(sec, 0).size());
            if (np == 0) continue;
            auto fon_vec = reks::studio::fon_vector(reel_, sa_cassette_, sec, 0);
            auto fon_grad = reks::fon::compute_fon_gradient(
                sec, 0, sa_cassette_, reel_.fon_state[sec], reel_.E_L);
            auto fon_hess = reks::fon::compute_fon_hessian(
                sec, 0, sa_cassette_, reel_.fon_state[sec], reel_.E_L);
            const std::string sec_tok = (n_sec > 1) ? ("sector " + std::to_string(sec) + " ") : "";
            outfile->Printf("  [GVB-FON] iter=%d %sFON=(", iteration_, sec_tok.c_str());
            for (int p = 0; p < np; ++p)
                outfile->Printf("%.8f%s", fon_vec[p], p + 1 < np ? ", " : "");
            outfile->Printf(") grad=(");
            for (int p = 0; p < np; ++p)
                outfile->Printf("%+.6e%s", fon_grad[p], p + 1 < np ? ", " : "");
            outfile->Printf(") hess_diag=(");
            for (int p = 0; p < np; ++p)
                outfile->Printf("%+.6e%s", fon_hess[p * np + p], p + 1 < np ? ", " : "");
            outfile->Printf(")\n");

            for (int p = 0; p < np; ++p) {
                if (fon_vec[p] > 1.99 || fon_vec[p] < 1.01) {
                    double newton_step = (std::abs(fon_hess[p * np + p]) > 1e-10)
                        ? -fon_grad[p] / fon_hess[p * np + p] : 0.0;
                    outfile->Printf("  [GVB-FON] iter=%d %sBOUNDARY pair %d: f=%.8f g=%+.6e h=%+.6e "
                                    "newton_step=%+.6e\n", iteration_, sec_tok.c_str(), p, fon_vec[p],
                                    fon_grad[p], fon_hess[p * np + p], newton_step);
                }
            }
        }
    }

    if (reports(4)) {
        reks::report::print_iteration_diagnostics(report_frame(), iteration_, trace_);
    }

    timer_off("REKS: form_G");
}

// GVB-DIIS convergence quality as wfn variables and a trace line; crossing CLEAN = both gates cleared
// with margin, GRAZING = within a decade of one.
void REKS::report_convergence_quality() {
    const reks::GvbDiisDriver::Quality q = gvb_->quality(nmopi_[0]);

    set_scalar_variable("REKS QUALITY GORB FRO", q.gorb_fro);
    set_scalar_variable("REKS QUALITY GORB RMS", q.gorb_rms);
    set_scalar_variable("REKS QUALITY LANDINGS", static_cast<double>(q.landings));
    set_scalar_variable("REKS QUALITY VERDICT TRIPS", static_cast<double>(q.verdict_trips));
    set_scalar_variable("REKS QUALITY MAX RUNG", static_cast<double>(q.rung));
    set_scalar_variable("REKS QUALITY BASIN FIRES",
                        static_cast<double>(q.basin_fires));
    set_scalar_variable("REKS QUALITY GRAZING", q.grazing ? 1.0 : 0.0);

    if (reports(4)) {
        static const char* kRung[] = {"ARMED", "RETAIN", "FALLBACK"};
        outfile->Printf(
            "  [QUALITY] iter=%d gorb_fro=%.3e gorb_rms=%.3e landings=%d trips=%d rung=%s "
            "basin_fires=%d crossing=%s\n",
            iteration_, q.gorb_fro, q.gorb_rms, q.landings, q.verdict_trips, kRung[q.rung],
            q.basin_fires, q.grazing ? "GRAZING" : "CLEAN");
    }
}

// SAD iteration: Fa_ = H_ + G_ from the SAD density.
void REKS::form_initial_F() {
    timer_on("REKS: form_F");
    if (reports(4)) {
        outfile->Printf("\n  form_F: SAD iteration\n");
    }
    Fa_->copy(H_);
    Fa_->add(G_);
    timer_off("REKS: form_F");
}

// Builds the REKS coupling Fock (reel_.F_reks_MO/F_reks); once GVB-DIIS is active
// also runs the orbital-update step, else leaves Fa_/Fb_ = F_reks.
void REKS::form_F() {
    timer_on("REKS: form_F");
    if (reports(4) && iteration_ <= 2) {
        outfile->Printf("  [E_MICRO] iter=%d", iteration_);
        { int shown = 0;
          for (int L : sa_cassette_.microstates()) {
              if (shown++ >= 8) { outfile->Printf(" ..."); break; }
              outfile->Printf(" E[%s]=%.6f",
                              reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.E_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [C_MICRO] iter=%d", iteration_);
        { int shown = 0;
          for (int L : sa_cassette_.microstates()) {
              if (shown++ >= 8) { outfile->Printf(" ..."); break; }
              outfile->Printf(" C[%s]=%.6f",
                              reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.C_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [FON] iter=%d", iteration_);
        for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
            const std::string chan = reks::studio::fon_channel_name(gen, sa_sector());
            const auto& gem = sa_cassette_.geminals_active_gen(gen);
            for (size_t k = 0; k < gem.size(); ++k) {
                const auto& fon = sa_fons().layers[gen][gem[k]];
                outfile->Printf(" %s-FON%zu=(%.8f, %.8f)", chan.c_str(), k, fon.p, fon.q);
            }
        }
        outfile->Printf("\n");
    }

    if (reports(4)) {
        outfile->Printf("\n  form_F: REKS iteration %d\n", iteration_);
    }

    int N = nmopi_[0];

    // C_L-weighted accumulated MO Fock for the uniform off-diagonal blocks.
    build_F_acc_MO();

    if (reports(4) && iteration_ > 0) {
        auto F_gen = build_generalized_fock();
        double** Fg = F_gen->pointer(0);
        int n_act = active_mo_indices_.size();

        double g_orb_sq = 0.0;
        std::vector<std::tuple<int, int, double>> g_orb_elems;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                int ai = active_mo_indices_[i];
                int aj = active_mo_indices_[j];
                double g_ij = -2.0 * (Fg[ai][aj] - Fg[aj][ai]);
                g_orb_elems.emplace_back(ai, aj, g_ij);
                g_orb_sq += g_ij * g_ij;
            }
        }

        auto g_fon = reks::fon::compute_fon_gradient(
            0, 0, sa_cassette_, sa_fons(), reel_.E_L);

        outfile->Printf("  [TRAH] iter=%d grad_check ||g_orb||=%.5e", iteration_, std::sqrt(g_orb_sq));
        for (size_t k = 0; k < g_fon.size(); ++k) {
            outfile->Printf(" g_fon%zu=%+.5e", k, g_fon[k]);
        }
        outfile->Printf(" g_orb:");
        for (auto& [ai, aj, g] : g_orb_elems) {
            outfile->Printf(" (%d,%d)=%+.5e", ai, aj, g);
        }
        outfile->Printf("\n");
    }

    // TRAH gate: orbital gradient of the Fock being judged, over TRAH's variational rotations.
    if (trah_live() && iteration_ > 0)
        trah_gorb_norm_ = compute_full_gradient_norm(build_generalized_fock(), same_role_aa_frozen());

    assemble_F_reks();

    // SC = S*Ca, the AO<-MO map: F_AO = SC F_MO SC^T.
    auto SC = linalg::doublet(S_, Ca_, false, false);

    // GVB-DIIS activation gate + accelerator step on this Fock; F_gen_ aliases the
    // build_generalized_fock workspace for the duration of the iteration.
    bool gvb_stepped = false;
    if (gvb_) {
        reks::GvbFockContext ctx{iteration_,
                                 use_trah_,
                                 active_mo_indices_,
                                 Ncore_,
                                 Ca_->colspi()[0],
                                 sa_cassette_,
                                 reel_,
                                 fon_blocks_,
                                 Ca_,
                                 Cb_,
                                 epsilon_a_,
                                 Fa_,
                                 Fb_,
                                 SC,
                                 functional_ ? functional_->x_alpha() : 1.0,
                                 [this]() { F_gen_ = build_generalized_fock(); return F_gen_; },
                                 same_role_aa_frozen()};
        // COSX grid switch changes the energy functional: the driver restarts as at activation.
        if (cosx_) {
            const std::string grid = std::dynamic_pointer_cast<CompositeJK>(jk_)->get_COSX_grid();
            if (!cosx_grid_.empty() && grid != cosx_grid_) gvb_->restart(ctx);
            cosx_grid_ = grid;
        }
        gvb_stepped = gvb_->on_fock(ctx);

        // Landing at MAXITER: re-evaluate the landed frame; its gorb becomes the reported gorb.
        if (gvb_->terminal_landing(iteration_)) {
            evaluate_frame();
            build_F_acc_MO();
            assemble_F_reks();
            gvb_->set_terminal_report(compute_full_gradient_norm(ctx.build_F_gen(), {}));
        }
    }

    if (!gvb_stepped) {
        auto F_reks_AO = linalg::triplet(SC, reel_.F_reks_MO, SC, false, false, true);
        reel_.F_reks->copy(F_reks_AO);
        Fa_->copy(reel_.F_reks);
        Fb_->copy(Fa_);
    }

    if (reports(4)) {
        double** Fmo = reel_.F_reks_MO->pointer(0);
        int n_act = active_mo_indices_.size();

        outfile->Printf("  [F_REKS] iter=%d        ", iteration_);
        for (int j = 0; j < n_act; ++j) {
            outfile->Printf("    %4d        ", active_mo_indices_[j]);
        }
        outfile->Printf("\n");
        for (int i = 0; i < n_act; ++i) {
            int mi = active_mo_indices_[i];
            outfile->Printf("  [F_REKS] iter=%d %4d", iteration_, mi);
            for (int j = 0; j < n_act; ++j) {
                int mj = active_mo_indices_[j];
                outfile->Printf("  %12.6f  ", Fmo[mi][mj]);
            }
            outfile->Printf("\n");
        }

        // Pair gaps + max off-diag (active block) + max core-active/active-virt couplings.
        int np_local = static_cast<int>(sa_cassette_.geminals_active_gen(0).size());
        outfile->Printf("  [F_REKS] iter=%d", iteration_);
        for (int p = 0; p < np_local; ++p) {
            int g = sa_cassette_.geminals_active_gen(0)[p];
            const auto& orbs = sa_cassette_.geminal_templates()[g].orbitals;
            int pp = orbs[0];
            int qq = orbs[1];
            int ip = active_mo_indices_[pp];
            int iq = active_mo_indices_[qq];
            outfile->Printf(" GVB(%d,%d)[%d,%d]=%.6f", pp, qq, ip, iq, Fmo[ip][ip] - Fmo[iq][iq]);
        }

        double maxOff = 0.0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                maxOff = std::max(maxOff, std::abs(Fmo[active_mo_indices_[i]][active_mo_indices_[j]]));
            }
        }
        outfile->Printf(" maxoff=%.6f", maxOff);

        double max_ca = 0.0, max_av = 0.0;
        int first_active = active_mo_indices_.front();
        int last_active = active_mo_indices_.back();
        for (int a : active_mo_indices_) {
            for (int c = 0; c < first_active; ++c) {
                max_ca = std::max(max_ca, std::abs(Fmo[c][a]));
            }
            for (int v = last_active + 1; v < N; ++v) {
                max_av = std::max(max_av, std::abs(Fmo[a][v]));
            }
        }
        outfile->Printf(" max_ca=%.6f max_av=%.6f\n", max_ca, max_av);
    }

    if (reports(5)) {
        outfile->Printf("\n  === REKS Coupling Fock Matrix ===\n");
        outfile->Printf("  tr(F_reks_AO) = %15.10f\n", reel_.F_reks->trace());
        outfile->Printf("  tr(F_reks_MO) = %15.10f\n", reel_.F_reks_MO->trace());
        const auto& lagr = reel_.lagrangians;
        for (size_t i = 0; i < lagr.size(); ++i) {
            outfile->Printf("  lagrangian[%zu] = %15.10f\n", i, lagr[i]);
        }

        double** Fmo = reel_.F_reks_MO->pointer(0);
        outfile->Printf("  REKS Fock diagonal: ");
        for (int i = 0; i < std::min(8, N); ++i) {
            outfile->Printf("F(%d)=%10.6f ", i, Fmo[i][i]);
        }
        outfile->Printf("\n");
        outfile->Printf("  REKS Fock active block:\n");
        for (size_t i = 0; i < active_mo_indices_.size(); ++i) {
            int mi = active_mo_indices_[i];
            for (size_t j = 0; j < active_mo_indices_.size(); ++j) {
                int mj = active_mo_indices_[j];
                outfile->Printf("    F(%d,%d)=%10.6f", mi, mj, Fmo[mi][mj]);
            }
            outfile->Printf("\n");
        }

        double max_asym = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double asym = std::abs(Fmo[i][j] - Fmo[j][i]);
                if (asym > max_asym) max_asym = asym;
            }
        }
        outfile->Printf("  Max MO asymmetry: %.3e\n", max_asym);
    }

    timer_off("REKS: form_F");
}

// [F,D]_S via HF::form_FDSmSDF, rescaled so rms(result) == ||g_orb||_F / nmo under TRAH or GVB-DIIS,
// with the landed frame's ||g_orb||_F on a GVB-DIIS landing at MAXITER.
SharedMatrix REKS::form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso) {
    auto result = HF::form_FDSmSDF(Fso, Dso);

    if (trah_live() && iteration_ > 0) {
        const double true_rms = gorb_rms(trah_gorb_norm_);
        set_scalar_variable("REKS GORB FRO", trah_gorb_norm_);
        set_scalar_variable("REKS GATE RMS", true_rms);
        const double current_rms = result->rms();
        if (current_rms > 1e-15)
            result->scale(true_rms / current_rms);
        else
            result->set(true_rms);
    }

    // TRAH: rms floored at ||V - I||_F / nmo of the same-role canonicalizer the orbitals still owe.
    if (trah_live()) {
        double canon_residual = 0.0;
        if (auto V = same_role_canonicalizer()) {
            const int n_act = static_cast<int>(active_mo_indices_.size());
            double** Vp = V->pointer(0);
            double dev_sq = 0.0;
            for (int i = 0; i < n_act; ++i)
                for (int j = 0; j < n_act; ++j) {
                    const double dev = Vp[i][j] - (i == j ? 1.0 : 0.0);
                    dev_sq += dev * dev;
                }
            canon_residual = std::sqrt(dev_sq);
        }
        const double canon_rms = gorb_rms(canon_residual);
        const double current_rms = result->rms();
        if (canon_rms > current_rms) {
            if (current_rms > 1e-15)
                result->scale(canon_rms / current_rms);
            else
                result->set(canon_rms);
            if (reports(4)) {
                outfile->Printf("  [CANON] iter=%d Dnorm %.2e -> %.2e (frame still rotating)\n",
                                iteration_, current_rms, canon_rms);
            }
        }
    }

    // Gate: current gorb. Exported: witnessed gorb; the returned rms takes it on a terminal landing.
    if (gvb_ && gvb_->active() && gvb_->gorb() > 0.0) {
        gvb_->advance_gate(iteration_, gorb_rms(gvb_->gorb()));
        const double report_rms = gorb_rms(gvb_->report_gorb());
        set_scalar_variable("REKS GORB FRO", gvb_->report_gorb());
        set_scalar_variable("REKS GATE RMS", report_rms);
        const double true_rms =
            gvb_->terminal_landing(iteration_) ? report_rms : gorb_rms(gvb_->gorb());
        double current_rms = result->rms();
        if (current_rms > 1e-15) {
            result->scale(true_rms / current_rms);
        }
    }

    return result;
}

// iter<=0: guess Ca_ from Fa_. Afterwards: the step goes to the stage buffers, the evaluated state stays live.
void REKS::form_C(double shift) {
    timer_on("REKS: form_C");
    // iter<=0 dispatch: READ-guess path if guess_Ca_ provided (keeps user MO layout),
    // else cold-start blind-diag of Fa_.
    if (iteration_ <= 0) {
        if (guess_Ca_) {
            if (reports(4)) {
                outfile->Printf(
                    "  [FORM_C] iter=%d READ-guess path: frozen active block (guess_fon_=%s)\n",
                    iteration_, has_read_fon_guess() ? "SET" : "nullptr");
            }
            form_C_read_guess_update();
            timer_off("REKS: form_C");
            return;
        }

        if (reports(4)) {
            outfile->Printf("  [FORM_C] iter=%d cold path: blind-diag of Fa_\n", iteration_);
        }
        diagonalize_F(Fa_, Ca_, epsilon_a_);
        find_occupation();

        if (reports(4)) {
            int N_guess = nmopi_[0];
            int first_virt_g = active_mo_indices_.empty() ? Ncore_ : (active_mo_indices_.back() + 1);
            double* eps_g = epsilon_a_->pointer(0);
            outfile->Printf("  [LEVEL_SHIFT] iter=%d GUESS eps:", iteration_);
            for (int c = std::max(0, Ncore_ - 4); c < Ncore_; ++c)
                outfile->Printf(" core[%d]=%+.6f", c, eps_g[c]);
            for (int a : active_mo_indices_) outfile->Printf(" act[%d]=%+.6f", a, eps_g[a]);
            for (int v = first_virt_g; v < N_guess && v < first_virt_g + 4; ++v)
                outfile->Printf(" virt[%d]=%+.6f", v, eps_g[v]);
            outfile->Printf("\n");
        }

        timer_off("REKS: form_C");
        return;
    }

    // Null step: form_F already restored the state -- hold Ca_/Cb_/epsilons untouched.
    if (gvb_ && gvb_->take_null_step()) {
        if (reports(4))
            outfile->Printf("  [GVB-DIIS] iter=%d form_C: null step, orbitals held\n", iteration_);
        timer_off("REKS: form_C");
        return;
    }

    stage_in();
    form_C_step(shift);
    stage_out();
    timer_off("REKS: form_C");
}

// Stage buffers := live Ca_, eps, FONs.
void REKS::stage_in() {
    if (!stage_Ca_) {
        stage_Ca_ = Ca_->clone();
    } else {
        stage_Ca_->copy(Ca_);
    }
    if (!stage_eps_) {
        stage_eps_ = std::make_shared<Vector>(*epsilon_a_);
    } else {
        stage_eps_->copy(*epsilon_a_);
    }
    stage_fons_ = reel_.fon_state;
}

// Swaps live and stage buffers in place: live = state at stage_in, buffers = stepped state.
void REKS::stage_out() {
    const size_t n_C = static_cast<size_t>(Ca_->rowspi()[0]) * Ca_->colspi()[0];
    double* C_live = Ca_->pointer(0)[0];
    std::swap_ranges(C_live, C_live + n_C, stage_Ca_->pointer(0)[0]);
    double* eps_live = epsilon_a_->pointer(0);
    std::swap_ranges(eps_live, eps_live + epsilon_a_->dimpi()[0], stage_eps_->pointer(0));
    std::swap(reel_.fon_state, stage_fons_);
    staged_ = true;
}

// Orbital step on the live state: ORB-DIIS rotation, or TRAH step and F_MO diag with staircase shift and swap tracking.
void REKS::form_C_step(double shift) {
    // ORBITAL-DIIS path (Ionova-Carter): cumulative kappa rotates Ca_, bypassing Fock-diag.
    if (gvb_ && gvb_->active() && gvb_->orbital_runtime_active()) {
        const int N_orb = Ca_->colspi()[0];
        gvb_->apply_orbital_rotation(Ca_, active_mo_indices_, Ncore_, N_orb);

        // Eigenvalue placeholder: eps_i = (Ca_^T Fa_ Ca_)_ii via Y = Fa_ Ca_ once and a
        // column dot per i; no full nmo x nmo congruence product is built.
        {
            const int nso_orb = Ca_->rowspi()[0];
            std::vector<double> Y(static_cast<size_t>(nso_orb) * N_orb);
            double** Cap_orb = Ca_->pointer(0);
            C_DGEMM('N', 'N', nso_orb, N_orb, nso_orb, 1.0, Fa_->pointer(0)[0], nso_orb,
                    Cap_orb[0], N_orb, 0.0, Y.data(), N_orb);
            double* eps_a_orb = epsilon_a_->pointer(0);
            for (int i = 0; i < N_orb; ++i)
                eps_a_orb[i] = C_DDOT(nso_orb, &Cap_orb[0][i], N_orb, Y.data() + i, N_orb);
        }
        epsilon_b_->copy(*epsilon_a_);

        find_occupation();
        if (reports(4)) {
            double max_kappa = 0.0;
            for (double v : gvb_->cum_kappa()) max_kappa = std::max(max_kappa, std::abs(v));
            outfile->Printf("  [ORB-DIIS] iter=%d form_C: applied expm(K), max|kappa_cum|=%.4f\n",
                            iteration_, max_kappa);
        }
        return;
    }

    // TRAH joint kappa+FON step before Fock diag, so F_MO sees rotated Ca_.
    if (trah_live() && iteration_ > 0) {
        const int n_act = static_cast<int>(active_mo_indices_.size());
        reks::ActiveEriTiles hess_eri;
        reks::TrahStepContext ctx{iteration_,
                                  active_mo_indices_,
                                  sa_cassette_,
                                  reel_,
                                  fon_blocks_,
                                  reks_fon_lower_,
                                  nullptr,
                                  &hess_eri,
                                  functional_->x_alpha(),
                                  Ca_,
                                  epsilon_a_,
                                  guess_Ca_ && has_read_fon_guess(),
                                  guess_Ca_ != nullptr,
                                  gvb_ && gvb_->active(),
                                  n_sa_microstates()};
        ctx.F_gen = build_generalized_fock();
        ctx.aa_frozen = same_role_aa_frozen();

        // Active-MO ERI tiles of all pairs p <= q: a full JK pass per iteration, report level 5 only.
        if (reports(5) && n_act >= 2) {
            std::vector<std::pair<int, int>> canonical_pairs;
            canonical_pairs.reserve(static_cast<size_t>(n_act) * (n_act + 1) / 2);
            for (int p = 0; p < n_act; ++p)
                for (int q = p; q < n_act; ++q) canonical_pairs.emplace_back(p, q);
            hess_eri = compute_active_mo_eri(active_mo_indices_, canonical_pairs);
        }

        trah_->step(ctx);
        if (ctx.clear_read_guess) clear_read_fon_guess();
    }

    int N = nmopi_[0];

    // F_MO = Ca^T Fa Ca, diagonalized in the current MO basis.
    auto F_MO_diis = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);

    // Freeze active-active off-diag so Fock diag does not undo TRAH rotations.
    if (use_trah_ && iteration_ > 0 && !(gvb_ && gvb_->active())) {
        int n_act = static_cast<int>(active_mo_indices_.size());
        double** Fmo = F_MO_diis->pointer(0);
        for (int i = 0; i < n_act; ++i)
            for (int j = 0; j < n_act; ++j)
                if (i != j) Fmo[active_mo_indices_[i]][active_mo_indices_[j]] = 0.0;
    }

    auto U = std::make_shared<Matrix>("Eigenvectors", N, N);
    auto eps = std::make_shared<Vector>("Eigenvalues", N);

    double** Fmo_guard = F_MO_diis->pointer(0);

    std::vector<int> ls_window;
    for (int c = std::max(0, Ncore_ - 4); c < Ncore_; ++c) ls_window.push_back(c);
    for (int a : active_mo_indices_) ls_window.push_back(a);
    int ls_first_virt = active_mo_indices_.empty() ? Ncore_ : (active_mo_indices_.back() + 1);
    for (int v = ls_first_virt; v < N && v < ls_first_virt + 4; ++v) ls_window.push_back(v);

    auto ls_print_diag = [&](const char* stage, double* eps_or_diag, bool is_eps) {
        outfile->Printf("  [LEVEL_SHIFT] iter=%d %s", iteration_, stage);
        for (int idx : ls_window) {
            const char* tag = (idx < Ncore_) ? "core" : (idx >= ls_first_virt) ? "virt" : "act";
            double val = is_eps ? eps_or_diag[idx] : Fmo_guard[idx][idx];
            outfile->Printf(" %s[%d]=%+.6f", tag, idx, val);
        }
        outfile->Printf("\n");
    };

    if (reports(4)) {
        ls_print_diag("before_shift F_diag:", nullptr, false);
    }

    // Staircase shift: F_MO[i][i] += s_i (s = (k+1) ls on active slot k, swapped in a GVB pair if F_p > F_q;
    // (n_act+1) ls on virtuals); diagonalize, eps_j -= sum_i U_ij^2 s_i, MOM orders U by max_i |U_ij|.
    std::vector<std::pair<int, int>> gvb_pairs;
    for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
        int g = sa_cassette_.geminals_active_gen(0)[u];
        const auto& orbs = sa_cassette_.geminal_templates()[g].orbitals;
        gvb_pairs.emplace_back(orbs[0], orbs[1]);
    }
    double applied_shift = orbital_guard_.pre_diag(Fmo_guard, N, active_mo_indices_, gvb_pairs,
                                                   iteration_, shift);

    if (reports(4)) {
        ls_print_diag("after_shift  F_diag:", nullptr, false);
    }

    if (reports(4) && applied_shift > 0.0) {
        if (std::abs(applied_shift - shift) > 1e-6) {
            outfile->Printf("  [FORM_C] iter=%d adaptive_shift: %.4f (user=%.3f)\n", iteration_, applied_shift, shift);
        }
        for (const auto& inv : orbital_guard_.inverted_pairs()) {
            outfile->Printf("  [FORM_C] iter=%d staircase_swap: pair(%d,%d) F(%d)=%.6f > F(%d)=%.6f\n", iteration_,
                            inv.p, inv.q, active_mo_indices_[inv.p], inv.F_p, active_mo_indices_[inv.q], inv.F_q);
        }
    }

    F_MO_diis->diagonalize(U, eps);

    double** Up = U->pointer(0);
    double* eps_raw = eps->pointer();

    if (reports(4)) {
        outfile->Printf("  [LEVEL_SHIFT] iter=%d after_diag (shifted) eps:", iteration_);
        for (int idx : ls_window) {
            const char* tag = (idx < Ncore_) ? "core" : (idx >= ls_first_virt) ? "virt" : "act";
            outfile->Printf(" %s[%d]=%+.6f", tag, idx, eps_raw[idx]);
        }
        outfile->Printf("\n");

        // U eigenvector block: rows = old MO (basis index), cols = new MO (slot index)
        outfile->Printf("  [LEVEL_SHIFT] iter=%d U overlap (rows=old MO, cols=new MO):\n", iteration_);
        outfile->Printf("                       ");
        for (int j : ls_window) outfile->Printf(" %9d", j);
        outfile->Printf("\n");
        for (int i : ls_window) {
            const char* tag = (i < Ncore_) ? "core" : (i >= ls_first_virt) ? "virt" : "act ";
            outfile->Printf("           %s[%4d]", tag, i);
            for (int j : ls_window) outfile->Printf(" %+9.4f", Up[i][j]);
            outfile->Printf("\n");
        }
    }

    // Exact level-shift undo: eps[j] -= sum_i |U[i][j]|^2 * shift[i].
    if (applied_shift > 0.0) {
        const std::vector<double>& shifts = orbital_guard_.applied_shifts();
        // i-outer: U is row-major, so each row is one contiguous stream.
        std::vector<double> sub(N, 0.0);
        for (int i = 0; i < N; ++i) {
            if (shifts[i] == 0.0) continue;
            const double s = shifts[i];
            const double* ui = Up[i];
            for (int j = 0; j < N; ++j) sub[j] += s * ui[j] * ui[j];
        }
        for (int j = 0; j < N; ++j) eps_raw[j] -= sub[j];
        if (reports(4)) {
            outfile->Printf("  [LEVEL_SHIFT] iter=%d shift_undo  eps:", iteration_);
            for (int idx : ls_window) {
                const char* tag = (idx < Ncore_) ? "core" : (idx >= ls_first_virt) ? "virt" : "act";
                outfile->Printf(" %s[%d]=%+.6f", tag, idx, eps_raw[idx]);
            }
            outfile->Printf("\n");
        }
    }

    auto mixing = orbital_guard_.post_diag(Up, eps_raw, N, Ncore_, active_mo_indices_);

    if (reports(4)) {
        outfile->Printf("  [LEVEL_SHIFT] iter=%d after_reorder eps:", iteration_);
        for (int idx : ls_window) {
            const char* tag = (idx < Ncore_) ? "core" : (idx >= ls_first_virt) ? "virt" : "act";
            outfile->Printf(" %s[%d]=%+.6f", tag, idx, eps_raw[idx]);
        }
        outfile->Printf("\n");

        // Final U after phase fix + MOM reorder (rows=old MO, cols=new slot).
        outfile->Printf("  [LEVEL_SHIFT] iter=%d U overlap final (after phase+MOM):\n", iteration_);
        outfile->Printf("                       ");
        for (int j : ls_window) outfile->Printf(" %9d", j);
        outfile->Printf("\n");
        for (int i : ls_window) {
            const char* tag = (i < Ncore_) ? "core" : (i >= ls_first_virt) ? "virt" : "act ";
            outfile->Printf("           %s[%4d]", tag, i);
            for (int j : ls_window) outfile->Printf(" %+9.4f", Up[i][j]);
            outfile->Printf("\n");
        }
    }

    if (gvb_ && mixing.n_swaps > 0) {
        gvb_->note_swaps();
    }

    if (reports(4)) {
        int n_act = static_cast<int>(active_mo_indices_.size());
        outfile->Printf("  [FORM_C] iter=%d post_diag: eps=", iteration_);
        for (int k = 0; k < n_act; ++k) outfile->Printf(" %.6f", eps_raw[active_mo_indices_[k]]);
        outfile->Printf(" max_mix=%.4f mix_target=%d swaps=%d\n", mixing.max_active_mixing, mixing.max_mixing_target,
                        mixing.n_swaps);

        if (mixing.max_active_mixing > 0.1) {
            int first_virt = active_mo_indices_.back() + 1;
            for (int a : active_mo_indices_) {
                double max_mix = 0.0;
                int max_j = -1;
                for (int j = 0; j < N; ++j) {
                    if (j == a) continue;
                    if (std::abs(Up[a][j]) > max_mix) {
                        max_mix = std::abs(Up[a][j]);
                        max_j = j;
                    }
                }
                if (max_mix > 0.1) {
                    const char* block = (max_j < Ncore_) ? "core" : (max_j >= first_virt) ? "virt" : "active";
                    outfile->Printf("  [FORM_C] iter=%d mixing: MO %d diag=%.3f max_offdiag=%.3f -> MO %d (%s)\n",
                                    iteration_, a, Up[a][a], max_mix, max_j, block);
                }
            }
        }
    }

    auto C_new = linalg::doublet(Ca_, U, false, false);
    Ca_->copy(C_new);

    // eps already shift-clean: level-shift undo + reorder applied.
    double* eps_p = eps->pointer();
    double* eps_a = epsilon_a_->pointer(0);
    if (gvb_ && gvb_->active() && !gvb_->orbital_runtime_active() && reel_.F_reks_MO) {
        // CFM runtime: Fa_ = F^g carries an artificial diagonal, so eps_i = u_i . (F_reks_MO U)[:,i], the
        // diagonal of U^T F_reks_MO U.
        double** Up_c = U->pointer(0);
        std::vector<double> Y(static_cast<size_t>(N) * N);
        C_DGEMM('N', 'N', N, N, N, 1.0, reel_.F_reks_MO->pointer(0)[0], N, Up_c[0], N, 0.0,
                Y.data(), N);
        for (int i = 0; i < N; ++i) eps_a[i] = C_DDOT(N, &Up_c[0][i], N, Y.data() + i, N);
    } else {
        for (int i = 0; i < N; ++i) {
            eps_a[i] = eps_p[i];
        }
    }
    epsilon_b_->copy(*epsilon_a_);

    // Same-role frame TRAH excludes from its variational set.
    if (trah_live()) canonicalize_same_role_frame();

    if (reports(4)) {
        int n_act = active_mo_indices_.size();
        // U diagonal: 1.0 = no rotation, cos(theta) for rotation angle theta.
        char buf[48];
        std::string u_str;
        snprintf(buf, sizeof(buf), "  [FORM_C] iter=%d U_diag:", iteration_);
        u_str = buf;
        for (int i = 0; i < n_act; ++i) {
            int mi = active_mo_indices_[i];
            snprintf(buf, sizeof(buf), " U(%d,%d)=%.6f", mi, mi, Up[mi][mi]);
            u_str += buf;
        }
        u_str += " | offdiag:";
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                int mi = active_mo_indices_[i];
                int mj = active_mo_indices_[j];
                snprintf(buf, sizeof(buf), " U(%d,%d)=%.6f", mi, mj, Up[mi][mj]);
                u_str += buf;
            }
        }
        outfile->Printf("%s\n", u_str.c_str());
    }

    if (reports(5)) {
        outfile->Printf("\n  === MO-basis Orbital Update ===\n");
        int first_active = active_mo_indices_.front();
        int last_active = active_mo_indices_.back();
        outfile->Printf("  Eigenvalues around active space:\n");
        for (int i = std::max(0, first_active - 2); i < std::min(N, last_active + 3); ++i) {
            outfile->Printf("    eps(%d) = %10.6f\n", i, eps_p[i]);
        }
    }
}

void REKS::form_C_read_guess_update() {
    timer_on("REKS: form_C_read_guess_update");
    const int N = nmopi_[0];
    const int n_act = static_cast<int>(active_mo_indices_.size());

    // F_MO is built only for diagnostics + epsilon_a_ estimate; Ca_ untouched.
    auto F_MO = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
    double** Fmo = F_MO->pointer(0);

    double max_off_aa = 0.0;
    double max_off_av = 0.0;
    int last_active = (n_act > 0) ? active_mo_indices_.back() : -1;
    for (int i = 0; i < n_act; ++i) {
        int mi = active_mo_indices_[i];
        for (int j = 0; j < n_act; ++j) {
            if (i == j) continue;
            int mj = active_mo_indices_[j];
            max_off_aa = std::max(max_off_aa, std::abs(Fmo[mi][mj]));
        }
        for (int v = last_active + 1; v < N; ++v) {
            max_off_av = std::max(max_off_av, std::abs(Fmo[mi][v]));
        }
    }

    // F_MO diagonal -> epsilon_a_ (exact for selfread, 1st-order otherwise).
    double* eps_a = epsilon_a_->pointer(0);
    for (int i = 0; i < N; ++i) {
        eps_a[i] = Fmo[i][i];
    }
    epsilon_b_->copy(*epsilon_a_);

    find_occupation();

    if (reports(4)) {
        outfile->Printf(
            "  [FORM_C] iter=%d read_guess_update: NO-OP on Ca_ (joint kappa+FON "
            "delegated to combined_step iter>=1)\n",
            iteration_);
        outfile->Printf(
            "  [FORM_C] iter=%d read_guess_update: F_MO max_off_aa=%.3e max_off_av=%.3e\n",
            iteration_, max_off_aa, max_off_av);
        outfile->Printf("  [FORM_C] iter=%d read_guess_update: F_MO active diag =",
                        iteration_);
        for (int i = 0; i < n_act; ++i) {
            int mi = active_mo_indices_[i];
            outfile->Printf(" %+.6f", Fmo[mi][mi]);
        }
        outfile->Printf("\n");

        int first_virt_r = active_mo_indices_.empty() ? Ncore_ : (active_mo_indices_.back() + 1);
        outfile->Printf("  [LEVEL_SHIFT] iter=%d GUESS_read eps:", iteration_);
        for (int c = std::max(0, Ncore_ - 4); c < Ncore_; ++c)
            outfile->Printf(" core[%d]=%+.6f", c, eps_a[c]);
        for (int a : active_mo_indices_) outfile->Printf(" act[%d]=%+.6f", a, eps_a[a]);
        for (int v = first_virt_r; v < N && v < first_virt_r + 4; ++v)
            outfile->Printf(" virt[%d]=%+.6f", v, eps_a[v]);
        outfile->Printf("\n");
    }
    timer_off("REKS: form_C_read_guess_update");
}

// SAD guess: RHF-style energy of Da_; otherwise E_SA = sum_L C_L E_L or the armed stand-on-best override.
double REKS::compute_E() {
    if (sad_ && iteration_ <= 0) {
        double one_electron_E = 2.0 * Da_->vector_dot(H_);
        double kinetic_E = 2.0 * Da_->vector_dot(T_);
        double coulomb_E = 2.0 * Da_->vector_dot(J_);

        double XC_E = 0.0;
        if (functional_->needs_xc() && rv_potential_) {
            XC_E = rv_potential_->quadrature_values()["FUNCTIONAL"];
        }

        double exchange_E = 0.0;
        double alpha = functional_->x_alpha();
        double beta_lrc = functional_->x_beta();

        if (functional_->is_x_hybrid() && K_) {
            exchange_E -= alpha * Da_->vector_dot(K_);
        }
        if (functional_->is_x_lrc()) {
            exchange_E -= beta_lrc * Da_->vector_dot(wK_);
        }

        double dashD_E = scalar_variable("-D Energy");

        energies_["Nuclear"] = nuclearrep_;
        energies_["Kinetic"] = kinetic_E;
        energies_["One-Electron"] = one_electron_E;
        energies_["Two-Electron"] = coulomb_E + exchange_E;
        energies_["XC"] = XC_E;
        energies_["-D"] = dashD_E;

        return nuclearrep_ + one_electron_E + coulomb_E + exchange_E + XC_E + dashD_E;
    }

    double kinetic_E = 2.0 * Da_->vector_dot(T_);
    energies_["Kinetic"] = kinetic_E;

    // E_SA = sum over the SCF pool of C_L * E_L (SI-only microstates are not
    // in sa_cassette_.microstates() and so do not contribute).
    double E_SA = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);

    // C_L-weighted components: E_1e_L = tr(Da H) + tr(Db H),
    // E_2e_pure_L = 0.5 (Da Fa + Db Fb - E_1e_L) - 0.5 tr(D_L V_xc_L).
    const bool needs_xc = functional_->needs_xc();
    double one_electron_E = 0.0;
    double two_electron_E = 0.0;
    double xc_E = 0.0;
    for (int L : sa_cassette_.microstates()) {
        double C_L = reel_.C_L[L];
        if (C_L == 0.0) continue;
        const auto e = reks::scf::microstate_energy(reel_, L, base_idx(L, /*alpha=*/true),
                                                    base_idx(L, /*alpha=*/false), 0.0, needs_xc);
        one_electron_E += C_L * e.E_1e;
        two_electron_E += C_L * (e.E_2e - 0.5 * e.tr_D_Vxc);
        if (needs_xc) {
            xc_E += C_L * reel_.E_xc_L[L];
        }
    }

    energies_["Nuclear"]      = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"]           = xc_E;
    energies_["REKS SA"]      = E_SA;
    energies_["Total Energy"] = E_SA;

    if (reports(5)) {
        outfile->Printf("\n  === REKS State-Averaged Energy ===\n");
        outfile->Printf("  E_SA    = %20.12f\n", E_SA);
    }

    // GVB-DIIS reports the best iterate's energies under a stand-on-best override and latches dE.
    return gvb_ ? gvb_->report_energy(energies_, iteration_) : E_SA;
}

SharedVector REKS::occupation_a() const {
    SharedVector occ = HF::occupation_a();
    if (reel_.fon_state.empty()) return occ;
    for (int g : sa_cassette_.geminals_active_gen(0)) {    // gen 0 = n
        const auto& tpl = sa_cassette_.geminal_templates()[g];
        if (tpl.scheme != 0) continue;                 // scheme 0 = canonical SCF pairs
        const int mo_p = active_mo_indices_[tpl.orbitals[0]];
        const int mo_q = active_mo_indices_[tpl.orbitals[1]];
        occ->set(0, mo_p, 0.5 * reel_.fon_state[0].layers[0][g].p);
        occ->set(0, mo_q, 0.5 * reel_.fon_state[0].layers[0][g].q);
    }
    return occ;
}

const reks::SIResult& REKS::primary_si_result() const {
    static const reks::SIResult empty;
    return si_results_.empty() ? empty : si_results_.front();
}

std::string REKS::get_microstate_label(int L) const {
    return reks::report::microstate_label(sa_cassette_, L);
}

double REKS::get_f_value(int g, int s) const {
    if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
    const auto& layer = reel_.fon_state[s].layers[0];
    if (g < 0 || g >= static_cast<int>(layer.size())) return 0.0;
    const auto& fon = layer[g];
    return reks::f_interp(fon.p * fon.q);
}

SharedMatrix REKS::SI_hamiltonian() const {
    const reks::SIResult& sir = primary_si_result();
    int n = sir.n_states;
    auto H = std::make_shared<Matrix>("SI Hamiltonian", n, n);
    for (int i = 0; i < n; i++)
        C_DCOPY(n, const_cast<double*>(sir.hamiltonian.data() + static_cast<size_t>(i) * n),
                1, H->pointer(0)[i], 1);
    return H;
}

SharedVector REKS::SI_energies() const {
    const reks::SIResult& sir = primary_si_result();
    int n = sir.n_report();
    auto E = std::make_shared<Vector>("SI Energies", n);
    for (int i = 0; i < n; i++) E->set(0, i, sir.energies[i]);
    return E;
}

SharedMatrix REKS::SI_coefficients() const {
    // Rows are states and follow the report cap; columns stay at the config dimension.
    const reks::SIResult& sir = primary_si_result();
    int n = sir.n_states;
    int n_rep = sir.n_show(n);
    auto C = std::make_shared<Matrix>("SI Coefficients", n_rep, n);
    for (int i = 0; i < n_rep; i++)
        C_DCOPY(n, const_cast<double*>(sir.coeffs.data() + static_cast<size_t>(i) * n),
                1, C->pointer(0)[i], 1);
    return C;
}

SharedMatrix REKS::SI_overlap() const {
    // sir.overlap is stored by blocks; assembled here into a dense n x n matrix.
    const reks::SIResult& sir = primary_si_result();
    const reks::BlockedMatrix& blocked = sir.overlap;
    const int n = sir.n_states;
    auto S = std::make_shared<Matrix>("SI Overlap", n, n);
    for (int c = 0; c < blocked.n_blocks(); ++c) {
        const std::vector<int>& idx = blocked.blocks[c];
        const int b = static_cast<int>(idx.size());
        const double* const tile = blocked.tile(c);
        for (int r = 0; r < b; ++r)
            for (int s = 0; s < b; ++s)
                S->set(0, idx[r], idx[s], tile[static_cast<size_t>(r) * b + s]);
    }
    return S;
}

std::vector<std::string> REKS::iteration_labels() const {
    std::vector<std::string> labels;
    if (iteration_ <= 0) return labels;

    if (gvb_ && gvb_->active()) {
        if (!gvb_->orbital_formulation()) {
            labels.emplace_back("CFM-GVB-DIIS");
        } else if (gvb_->orbital_runtime_active()) {
            labels.emplace_back("ORB-GVB-DIIS");
        } else {
            labels.emplace_back("CFM-GVB-DIIS");
        }
    } else if (use_trah_) {
        labels.emplace_back("TRAH");
    } else if (use_plain_) {
        labels.emplace_back("PLAIN");
    }
    if (trah_ && trah_->kdiis_active()) {
        labels.emplace_back("kDIIS");
    }
    return labels;
}

}  // namespace scf
}  // namespace psi
