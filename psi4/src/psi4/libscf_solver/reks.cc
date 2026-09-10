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

#include "reks.h"
#include "reks_arch.h"
#include "reks_diis_step.h"
#include "reks_objectives.h"
#include "reks_minres.h"
#include "reks_si.h"
#include "reks_rdm.h"
#include "reks_naturals.h"
#include "reks_report.h"
#include "reks_fon_relax.h"

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
#include "psi4/libmints/mintshelper.h"
#include "psi4/physconst.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <array>
#include <limits>
#include <cctype>
#include <set>
#include <cassert>
#include <chrono>
#include <random>
#include <tuple>
#include <cstdio>
#include <cstdlib>

namespace psi {
namespace scf {

using reks::report::reports;

namespace {

// RAII wall-time accumulator into a named bucket. nullptr bucket -> no-op.
class ScopedStage {
   public:
    ScopedStage(std::map<std::string, double>* bucket, const char* key)
        : bucket_(bucket), key_(key) {
        if (bucket_) start_ = std::chrono::steady_clock::now();
    }
    ~ScopedStage() { release(); }
    void release() {
        if (bucket_) {
            (*bucket_)[key_] += std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - start_)
                                    .count();
            bucket_ = nullptr;
        }
    }
    ScopedStage(const ScopedStage&) = delete;
    ScopedStage& operator=(const ScopedStage&) = delete;

   private:
    std::map<std::string, double>* bucket_;
    const char* key_;
    std::chrono::steady_clock::time_point start_;
};

// Returns true iff user assigned a value to either form of an alias pair.
bool alias_changed(Options& opts,
                   const std::string& primary,
                   const std::string& alias) {
    return opts[primary].has_changed() || opts[alias].has_changed();
}

// Returns the Data entry for whichever alias the user set; throws if both
// were set. Precondition: alias_changed() holds, else returns primary's
// (default-valued) Data.
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

// One config-pool entry -> the global K indices it names. A bulk token ("full" =
// whole sector block; a config type = every config of that type; "sps" = the
// spin-projected base; "single-scheme[:N]" / "mixed-scheme[:N]" = the pairing-
// scheme selections; a pattern with '*' or '?' = the names it matches) expands to
// many, a config name to one.
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

// One config-pool array (flat SA list or one SI row) -> global K indices. A
// string element is a config-pool entry, or "exclude:<entry>", which subtracts
// what that entry names from everything the list adds, whatever their order. A
// numeric element is a sector-LOCAL int index (shifted by the selected sector's
// config_start). No range/dup validation here.
std::vector<int> parse_config_indices(Data& arr, const reks::studio::Catalog& cat,
                                      bool allow_exclude) {
    const int config_start = cat.sector_config_start();
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

// Qualify a wavefunction-variable name by run sector s. Sector 0 keeps the exact
// existing name (byte-stable single-sector output). Sector s > 0 becomes
// "REKS SECTOR s <tail>", where <tail> drops a leading "REKS " when present so
// REKS-prefixed names do not stutter ("REKS FON" -> "REKS SECTOR 1 FON"; SSR-
// prefixed names keep their prefix -> "REKS SECTOR 1 SSR ENERGIES K=0").
std::string sector_psivar_name(int s, const std::string& name) {
    if (s == 0) return name;
    const std::string reks_prefix = "REKS ";
    const std::string tail =
        (name.rfind(reks_prefix, 0) == 0) ? name.substr(reks_prefix.size()) : name;
    return "REKS SECTOR " + std::to_string(s) + " " + tail;
}

}  // anonymous namespace

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    reks_common_init();
}

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    reks_common_init();
}

REKS::~REKS() {
    if (uv_potential_) {
        uv_potential_->finalize();
    }
    if (rv_potential_) {
        rv_potential_->finalize();
    }
}

void REKS::setup_potential() {
    if (functional_->needs_xc()) {
        polarized_functional_ = functional_->build_polarized();
        uv_potential_ = std::make_shared<psi::UV>(polarized_functional_, basisset_, options_);
        uv_potential_->initialize();

        // Both potentials integrate the same molecule with the same options, so one grid serves both.
        rv_potential_ = std::make_shared<psi::RV>(functional_, basisset_, options_);
        rv_potential_->set_grid(uv_potential_->grid());
        rv_potential_->initialize();

        uv_potential_->share_collocation_cache_from(*rv_potential_);

        if (reports(4)) {
            outfile->Printf("\n  REKS: Initialized UV potential for spin-polarized XC\n");
            outfile->Printf("        Functional: %s -> %s (polarized)\n", functional_->name().c_str(),
                            polarized_functional_->name().c_str());
            outfile->Printf("  REKS: Initialized RV potential for closed-shell XC on Da_\n");
        }
    }
}

void REKS::save_density_and_energy() { Dold_->copy(Da_); }

double REKS::compute_initial_E() { return nuclearrep_ + Da_->vector_dot(H_); }

void REKS::damping_update(double damp) {
    Da_->scale(1.0 - damp);
    Da_->axpy(damp, Dold_);
}

int REKS::complete_orbital_basis(SharedMatrix C, SharedMatrix S, SharedMatrix X,
                                  int n_loaded, int h) {
    const int nso = C->rowspi()[h];
    const int nmo = C->colspi()[h];
    if (!S || !X) {
        throw PSIEXCEPTION(
            "REKS::complete_orbital_basis: S_ or X_ is null at guess time");
    }

    double** Cp = C->pointer(h);
    double** Sp = S->pointer(h);
    double** Xp = X->pointer(h);

    // S-metric MGS on irrep h: c <- c - sum_k (c_k^T S c) c_k over accepted cols c_k,
    // then c <- c / sqrt(c^T S c). Loaded cols 0..n_loaded-1 seed the set (kept iff
    // c^T S c > DEP_TOL after projection, skipped entirely if already orthonormal to
    // SKIP_TOL); cols kept..nmo-1 fill from X through the same projection, threshold
    // FILL_TOL.
    constexpr double FILL_TOL = 1.0e-8;
    // Squared S-norm below which a loaded column is spanned by its predecessors. Looser
    // than FILL_TOL: normalizing such a column amplifies noise by 1/sqrt(DEP_TOL).
    constexpr double DEP_TOL = 1.0e-6;
    // max|Cl^T S Cl - I|, Cl = C[:, 0..n_loaded-1], under which the loaded block needs
    // no repair.
    constexpr double SKIP_TOL = 1.0e-10;

    // G = Cl^T S Cl, worst = max|G - I|; the lda=nmo stride selects Cl's columns out of C.
    // T = S Cl survives the check: the loaded_ok branch below reads its columns instead
    // of re-forming S c per column.
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

    // Sc = S c. Row-major wrapper over S (nso x nso); S is symmetric so the
    // transpose flag is immaterial.
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
        // Already orthonormal; the fill below still needs S times each loaded column,
        // and column j of T is exactly that.
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

// The whole proxy occupation lands in irrep 0: REKS forces C1.
void REKS::find_occupation() {
    if (multiplicity_ > 1) {
        const int n_occ_proxy = Ncore_ + catalog_.data->n_active_orbitals;
        for (int h = 0; h < nirrep_; ++h) {
            nalphapi_[h] = (h == 0) ? n_occ_proxy : 0;
            nbetapi_[h] = nalphapi_[h];
        }
    } else {
        HF::find_occupation();
    }
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

    if (guess_Ca_->nirrep() != nirrep_) {
        throw PSIEXCEPTION("REKS::guess READ: irrep count of guess orbitals does not match wavefunction.");
    }
    if (guess_Ca_->rowspi() != nsopi_) {
        throw PSIEXCEPTION("REKS::guess READ: nsopi of guess orbitals does not match wavefunction.");
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
            "    active-window layout may no longer match guess_mo_swaps.\n",
            label, h, n_loaded - kept, n_loaded);
    };
    for (int h = 0; h < nirrep_; ++h) {
        int kept = complete_orbital_basis(Ca_, S_, X_, n_loaded_a[h], h);
        report_dropped("Ca", h, kept, n_loaded_a[h]);
    }

    // Write the READ guess into reel_.fon_state[s] per FON generation and run
    // sector: gen 0 from guess_fon_[s], gen >= 1 from guess_fon_upper_[s][gen].
    // FONs clamp into (0, 2). A size mismatch clears that generation's guess.
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

    if (reports(4) && catalog_.data) {
        for (int s = 0; s < sa_cassette_.n_sectors(); ++s)
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
                const auto& act = sa_cassette_.geminals_active(s, gen);
                if (act.empty()) continue;
                const char letter = reks::studio::generation_letter(gen);
                const std::string chan = reks::studio::fon_channel_name(gen, s);
                outfile->Printf("  [READ_DBG] REKS::guess() post-restore %s-FON:", chan.c_str());
                for (int k = 0; k < static_cast<int>(act.size()); ++k)
                    outfile->Printf(" %c%d=%.15f", letter, k,
                                    reel_.fon_state[s].layers[gen][act[k]].p);
                outfile->Printf("\n");
            }
    }

    // Seed Fa_ from the saved wavefunction; skip on shape mismatch.
    if (guess_Fa_) {
        if (guess_Fa_->nirrep() == nirrep_ && guess_Fa_->rowspi() == nsopi_
                && guess_Fa_->colspi() == nsopi_) {
            Fa_->copy(guess_Fa_);
            if (Fb_ && Fb_.get() != Fa_.get()) {
                Fb_->copy(guess_Fa_);
            }
            if (reports(4)) {
                // Fa_ is square (C1): rms() * nrow() = sqrt(sum A_ij^2) = ||Fa||_F.
                outfile->Printf("  [READ_DBG] REKS::guess() seeded Fa_ from guess_Fa_, "
                                "||Fa||_F=%.4e trace=%.6e\n",
                                Fa_->rms() * Fa_->nrow(), Fa_->trace());
            }
        } else if (reports(2)) {
            outfile->Printf(
                "  REKS READ: WARNING -- guess_Fa_ shape (%dx%d) does not "
                "match Fa_; Fa transport skipped.\n",
                guess_Fa_->nrow(), guess_Fa_->ncol());
        }
    }

    if (S_) {
        // Orthonormality in the AO metric: G = C^T S C must equal I.
        // Diagnostic only, no side effects on the wavefunction.
        auto G_ortho = linalg::triplet(Ca_, S_, Ca_, true, false, false);
        const int nmo = G_ortho->colspi()[0];
        double** Gp = G_ortho->pointer(0);
        double max_diag_err = 0.0, max_offdiag = 0.0;
        for (int i = 0; i < nmo; ++i) {
            for (int j = 0; j < nmo; ++j) {
                const double s = Gp[i][j];
                if (i == j) max_diag_err = std::max(max_diag_err, std::abs(s - 1.0));
                else max_offdiag = std::max(max_offdiag, std::abs(s));
            }
        }
        constexpr double ORTHO_TOL = 1e-8;
        if ((max_diag_err > ORTHO_TOL || max_offdiag > ORTHO_TOL) && reports(2)) {
            outfile->Printf("  REKS READ: WARNING -- Ca is NOT orthonormal in the current "
                            "AO metric after re-orthonormalization:\n");
            outfile->Printf("    max|diag(C^T S C) - 1| = %.3e\n", max_diag_err);
            outfile->Printf("    max|offdiag(C^T S C)|  = %.3e\n", max_offdiag);
            outfile->Printf("    complete_orbital_basis failed to repair the guess metric; "
                            "the SCF cannot recover it.\n");
        } else if (reports(4)) {
            outfile->Printf("  [READ_DBG] REKS::guess() orthonormality check: "
                            "max|diag-1|=%.2e  max|off|=%.2e (OK)\n",
                            max_diag_err, max_offdiag);
        }
    }

    format_guess();
    form_D();
    iteration_ = -1;
}

void REKS::form_V() {
    if (!rv_potential_) return;

    rv_potential_->set_D({Da_});
    rv_potential_->compute_V({Va_});  // Vb_ aliases Va_
}

void REKS::finalize() {
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
            sa_cassette_, e_config, e_extra,
            energies_["REKS SA"], energies_["IPR Penalty"], energies_["VV10"]);
    }

    print_ipr_final_diagnostic();

    const double pr_min = active_pr_min();
    if (pr_min >= 0.0 && reports(2))
        outfile->Printf("  [PR] final: PR_min=%.6f over %d active MOs (%s populations)\n",
                        pr_min, static_cast<int>(active_mo_indices_.size()),
                        ipr_method_.c_str());

    if (use_gvb_diis_ && gvb_verdict_trips_ > 0 && reports(4))
        outfile->Printf("  [GVB-VERDICT] iter=%d rollbacks=%d best_gorb=%.3e\n", iteration_,
                        gvb_verdict_trips_, gvb_best_gorb_);

    if (use_gvb_diis_) report_convergence_quality();

    reel_.lagrangians_frozen = true;

    // Persists each active FON generation + active-space metadata, and the sector
    // count as "REKS N SECTORS": sector 0 keeps unqualified names (byte-stable
    // single-sector output), s > 0 prefixes "SECTOR s ".
    if (catalog_.data) {
        set_scalar_variable("REKS N SECTORS", static_cast<double>(sa_cassette_.n_sectors()));
        for (int s = 0; s < sa_cassette_.n_sectors(); ++s) {
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
                const auto& act = sa_cassette_.geminals_active(s, gen);
                const int n = static_cast<int>(act.size());
                if (n == 0) continue;
                const std::string base =
                    sector_psivar_name(s, reks::studio::fon_array_name(gen));
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

        // Cassette-independent dipole integrals, built once (frozen orbitals).
        DipoleIntegralCache dip_cache;
        if (reks::studio::has_diabatic_rdm(catalog_)) {
            timer_on("REKS: dipole_mo_cache");
            dip_cache = build_dipole_integral_cache();
            timer_off("REKS: dipole_mo_cache");
        }

        bool mo_consistency_printed = false;
        // Cassettes are grouped by sector; table tags restart on each sector change.
        int tagged_sector = -1;
        for (size_t i = 0; i < si_cassettes_.size(); ++i) {
            const reks::studio::Cassette& cassette = si_cassettes_[i];
            const reks::SIResult&         sir      = si_results_[i];
            const int si_cassette_pos = static_cast<int>(i);

            if (cassette.sector() != tagged_sector) {
                reks::report::reset_table_tags();
                tagged_sector = cassette.sector();
            }

            reks::rdm::DiabaticRdm rho_diab;
            if (reks::studio::has_diabatic_rdm(catalog_)) {
                timer_on("REKS: build_diabatic_rdm");
                ScopedStage _ss_bdr(reports(4) ? &post_scf_times_ : nullptr,
                                    "build_diabatic_rdm");
                rho_diab = reks::rdm::diabatic(reel_, cassette);
                timer_off("REKS: build_diabatic_rdm");
            }

            SIProperties props = compute_state_properties(cassette, sir, rho_diab, dip_cache);

            const auto& sa_micro = sa_cassette_.microstates();
            if (!sa_micro.empty() &&
                std::abs(reel_.E_L[sa_micro.front()]) >= reks::constants::ENERGY_THRESHOLD) {
                if (reports(5) && !mo_consistency_printed) {
                    print_mo_consistency();
                    mo_consistency_printed = true;
                }
                if (reports(5)) print_f_reks_mo_debug();

                print_state_properties_header(cassette, si_cassette_pos);

                reks::report::print_microstate_energies(
                    sa_cassette_, cassette.microstates_active, reel_.E_L);
                print_fon_table(cassette);

                reks::report::print_lagrangian_eri_table(
                    sa_cassette_, reel_.lagrangians, active_eri_pool_);

                print_si_hamiltonian_overlap(cassette, sir);
                reks::report::print_adiabatic_energies(sir);
                reks::report::print_excitation_energies(sir);

                print_state_dipoles(props, si_cassette_pos, sir);
            }

            register_si_wfn_variables(si_cassette_pos, cassette.sector(), sir, props,
                                      /*primary_alias=*/ i == 0, rho_diab);

            auto nos = compute_state_natural(cassette, sir, props.rho_adiab);
            print_state_natural(nos, si_cassette_pos, cassette, sir, rho_diab);
            register_natural_wfn_variables(si_cassette_pos, cassette.sector(), nos,
                                           /*primary_alias=*/ i == 0);
        }
    }

    if (uv_potential_) {
        uv_potential_->finalize();
    }
    if (rv_potential_) {
        rv_potential_->finalize();
    }

    // L = (C_occ diag(eps_occ)) C_occ^T over the occupied orbitals, i.e.
    // L_mn = sum_{i < nalphapi_[h]} eps_i C_mi C_ni. C_occ is the leading
    // nalphapi_[h] columns of Ca_, read in place at row stride nmopi_[h].
    // C_DGEMM leaves its output untouched at nocc == 0, so the zero comes first.
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

    Dold_.reset();
    G_.reset();
    J_.reset();
    K_.reset();
    wK_.reset();

    HF::finalize();
}

SharedMatrix REKS::build_generalized_fock() const {
    ScopedStage _ss_bgf(reports(4) ? &scf_iter_times_ : nullptr, "build_generalized_fock");

    timer_on("REKS: build_generalized_fock");

    int nmo = Ca_->colspi()[0];
    int n_active = sa_cassette_.n_active_orbitals();

    if (!F_gen_workspace_ || F_gen_workspace_->rowdim(0) != nmo ||
        F_gen_workspace_->coldim(0) != nmo) {
        F_gen_workspace_ = std::make_shared<Matrix>("Generalized Fock", nmo, nmo);
    }
    // Single-writer workspace: only F_gen_workspace_ and, when F_gen_ aliases it,
    // F_gen_ may hold the SharedMatrix; any extra ref would be silently overwritten.
    const long expected_refs = 1 + (F_gen_.get() == F_gen_workspace_.get() ? 1 : 0);
    if (F_gen_workspace_.use_count() > expected_refs) {
        throw PSIEXCEPTION(
            "build_generalized_fock: workspace SharedMatrix is still held by an "
            "outside caller (use_count=" +
            std::to_string(F_gen_workspace_.use_count()) + ", expected<=" +
            std::to_string(expected_refs) +
            "). Do not retain the returned matrix across another call.");
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

        // Virtual rows: n_alpha = n_beta = 0, no contribution.
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

    // IPR delocalization penalty (no-op when lambda_ipr_ == 0).
    add_ipr_contribution_to_F_gen(F_gen);

    timer_off("REKS: build_generalized_fock");
    return F_gen;
}

void REKS::build_F_acc_MO(bool core_rows_only) {
    timer_on("REKS: build_F_acc_MO");

    const int  n_active = sa_cassette_.n_active_orbitals();
    const bool needs_xc = functional_->needs_xc();
    const reks::scf::FockExchangePolicy pol =
        reks::scf::fock_exchange_policy(reel_, *functional_, *jk_);

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
        // reel_.V_acc_xc_{a,b} already hold Sum_L C_L V_xc_{a,b}^L as a single AO
        // back-projection; no per-L AO V_xc to sum here.
        F_acc_AO->add(reel_.V_acc_xc_a);
        F_acc_AO->add(reel_.V_acc_xc_b);
    }

    // AO -> MO: F_acc_MO = C^T F_acc_AO C (single transform, both spins folded in).
    const SharedMatrix& Ca = reel_.Ca_fock_transform;
    const int N  = Ca->rowspi()[0];   // nso (AO extent)
    const int Mn = Ca->colspi()[0];   // nmo (MO extent, reduced after form_Shalf)
    if (!reel_.F_acc_MO || reel_.F_acc_MO->coldim(0) != Mn) {
        reel_.F_acc_MO = std::make_shared<Matrix>("F_acc_MO", Mn, Mn);
    }
    double** Cp  = Ca->pointer(0);
    double** Tp  = reel_.temp_buffer->pointer(0);
    double** Aap = F_acc_AO->pointer(0);
    double** Amp = reel_.F_acc_MO->pointer(0);
    // T(nso x nmo) = F_AO(nso x nso) C(nso x nmo); F_acc_MO(nmo x nmo) = C^T T.
    // core_rows_only limits the second product to the first Ncore_ rows of F_acc_MO;
    // rows past Ncore_ keep their previous values. Cost: 2*Ncore_*nso*(nso+nmo) vs
    // 2*nso*nmo*(nso+Ncore_) for the full-width product.
    if (core_rows_only) {
        if (Ncore_ > 0) {
            C_DGEMM('T', 'N', Ncore_, N, N, 1.0, Cp[0], Mn, Aap[0], N, 0.0, Tp[0], N);
            C_DGEMM('N', 'N', Ncore_, Mn, N, 1.0, Tp[0], N, Cp[0], Mn, 0.0, Amp[0], Mn);
        }
    } else {
        C_DGEMM('N', 'N', N, Mn, N, 1.0, Aap[0], N, Cp[0], Mn, 0.0, Tp[0], N);
        C_DGEMM('T', 'N', Mn, Mn, N, 1.0, Cp[0], Mn, Tp[0], N, 0.0, Amp[0], Mn);
    }

    timer_off("REKS: build_F_acc_MO");
}

double REKS::compute_full_gradient_norm(SharedMatrix F_gen, reks::GradientBlocks* blocks) const {
    return reks::orbital_gradient_norm(F_gen->pointer(0), active_mo_indices_, Ncore_, nmopi_[0],
                                       blocks);
}

int REKS::nk_n_rot() const {
    const int nmo = Ca_->colspi()[0];
    const int n_act = static_cast<int>(active_mo_indices_.size());
    const int first_virt = (n_act > 0) ? active_mo_indices_.back() + 1 : Ncore_;
    return reks::OrbitalDIISEngine::n_nonred(Ncore_, n_act, nmo - first_virt);
}

std::vector<int> REKS::nk_fon_active_set() const {
    // Same bound test the micro-solver masks its convergence check with
    // (reks_fon_solver.h boundary_masking), over the fon_blocks_ order of capture_all_fons().
    const double upper = fon_micro_cfg_.upper_bound;
    std::vector<int> out;
    for (const auto& blk : fon_blocks_) {
        const double lower = reks_fon_lower_[blk.gen] >= 0.0 ? reks_fon_lower_[blk.gen] : 1e-8;
        for (double f : fon_vector(blk.s, blk.gen)) {
            if (f >= upper - 1e-10) out.push_back(1);
            else if (f <= lower + 1e-10) out.push_back(-1);
            else out.push_back(0);
        }
    }
    return out;
}

bool REKS::reduced_gradient_at(const std::vector<double>& kappa, const NKCurvature& curv,
                               std::vector<double>& g_out, std::vector<double>* denom_out,
                               bool commit) {
    const int nmo = Ca_->colspi()[0];

    // Probe state: everything else in reel_ is rebuilt below.
    SharedMatrix Ca_save = Ca_->clone();
    const auto fon_save = capture_all_fons();
    const auto active_set_base = nk_fon_active_set();
    const bool had_swaps_save = gvb_diis_had_swaps_;
    const int fon_swaps_save = gvb_fon_step_swaps_;

    Ca_->copy(orbital_diis_.apply_rotation(nk_Ca_anchor_, kappa, active_mo_indices_, Ncore_, nmo));

    // form_D: the occupied-column cache and the base densities the Fock build reads.
    C_occ_cache_ = Ca_subset("SO", "OCC");
    build_base_densities();
    form_Da_from_core();

    // form_G: SA Focks at the displaced orbitals, then the FON re-minimization that makes
    // the gradient reduced. No SA-Fock rebuild follows the FON move.
    nk_probe_gradient_only_ = (denom_out == nullptr);
    build_sa_focks();
    nk_probe_gradient_only_ = false;

    compute_sa_energies();

    // FON re-minimization: up to nk_fon_micro_ passes of gvb_fon_step(), stopping early once
    // the FON move is below 1e-10. The FON left by pass k is the pass-k+1 reference.
    auto fon_before = capture_all_fons();
    for (int fon_micro = 0; fon_micro < nk_fon_micro_; ++fon_micro) {
        gvb_fon_step();
        auto fon_after = capture_all_fons();
        double dn = 0.0;
        for (size_t b = 0; b < fon_after.size() && b < fon_before.size(); ++b)
            for (size_t i = 0; i < fon_after[b].size() && i < fon_before[b].size(); ++i)
                dn = std::max(dn, std::abs(fon_after[b][i] - fon_before[b][i]));
        fon_before = std::move(fon_after);
        if (dn < 1e-10) break;
    }

    // Probe path (denom_out == nullptr) needs only the core rows of F_acc_MO.
    build_F_acc_MO(/*core_rows_only=*/denom_out == nullptr);

    {
        SharedMatrix F_gen = build_generalized_fock();
        // Null denom_out is gradient-only: build_sa_focks() ran with need_full_mo_diag()
        // false, so reel_.F_MO_diag_*_L is valid over [Ncore, Ncore + n_act) and zero elsewhere.
        orbital_diis_.compute_newton_rotation(
            F_gen, reel_.F_MO_diag_a_L, reel_.F_MO_diag_b_L, reel_.C_L, sa_cassette_,
            active_mo_indices_, Ncore_, nmo, curv.gamma_aa, curv.gamma_ca, iteration_,
            curv.ipr_diag_aa, &g_out, denom_out, /*narrate=*/false);
    }  // F_gen released here: build_generalized_fock owns a single workspace

    // Usability: a column swap or a FON active-set move puts g_plus in a different problem
    // than g0. gvb_fon_step() zeroes its swap counter on entry, so a positive count is this
    // probe's own.
    const bool usable = (gvb_fon_step_swaps_ == 0) && (nk_fon_active_set() == active_set_base);

    // commit=true leaves orbitals, FON and reel_ at the probe's point; commit=false
    // restores them below, recomputing the weights derived from FON.
    if (!commit) {
        Ca_->copy(Ca_save);
        restore_all_fons(fon_save);
        compute_weighting_factors();
        gvb_diis_had_swaps_ = had_swaps_save;
        gvb_fon_step_swaps_ = fon_swaps_save;
    }

    return usable;
}

bool REKS::nk_Hv(const std::vector<double>& v, const std::vector<double>& g0, const SharedMatrix& Ghat,
                 double h, const NKCurvature& curv, std::vector<double>& Hv_out) {
    const int n = static_cast<int>(v.size());
    const double vnorm = C_DNRM2(v.size(), const_cast<double*>(v.data()), 1);
    if (vnorm == 0.0) {
        Hv_out.assign(n, 0.0);
        return true;
    }

    std::vector<double> kappa(n);
    const double scale = h / vnorm;
    for (int i = 0; i < n; ++i) kappa[i] = scale * v[i];

    std::vector<double> g_plus;
    if (!reduced_gradient_at(kappa, curv, g_plus, nullptr)) return false;

    // Chart correction, without which the operator is not symmetric and Lanczos is invalid.
    // The difference above is taken in the local chart C_anchor exp(K(kappa)), whose exponential
    // differential adds a purely antisymmetric term to the fixed-chart Hessian:
    //   J_ij = H_ij + (1/2) <g0, [E_j, E_i]>,  E_pq = e_p e_q^T - e_q e_p^T
    //   (S v)_i = (1/2) [Ghat, V]_(p_i, q_i),  V = unpack(v),  Ghat = unpack(g0)
    const int nmo = Ca_->colspi()[0];
    const int n_act_nk = static_cast<int>(active_mo_indices_.size());
    const int fv = (n_act_nk > 0) ? active_mo_indices_.back() + 1 : Ncore_;  // first virtual
    const int nv = nmo - fv;
    auto Vmat = std::make_shared<Matrix>("nk chart direction", nmo, nmo);
    reks::OrbitalDIISEngine::unpack_kappa_to_K(v, Vmat, active_mo_indices_, Ncore_, nmo);
    // Ghat and Vmat are skew, so (Ghat Vmat)^T = Vmat Ghat and the commutator is P - P^T.
    //
    // Both are built by unpack_kappa_to_K: [[A, B], [-B^T, 0]] over O = [0, fv) and
    // V = [fv, nmo), (V,V) block identically zero. Blockwise product needs three GEMMs,
    //   W_OO = Ag Av - Bg Bv^T,   W_OV = Ag Bv,   W_VO = -Bg^T Av,
    // dropping W_VV = -Bg^T Bv (unused by pack_K_to_kappa).
    // The sub-blocks are used in place with lda = nmo; W is written into Vmat's own storage
    // only after Vmat has been read, so no aliasing.
    auto Wmat = std::make_shared<Matrix>("nk chart commutator", nmo, nmo);
    {
        double** Gp = Ghat->pointer(0);
        double** Vp = Vmat->pointer(0);
        double** W = Wmat->pointer(0);
        // W_OO = Ag Av - Bg Bv^T
        if (fv > 0) {
            C_DGEMM('N', 'N', fv, fv, fv, 1.0, Gp[0], nmo, Vp[0], nmo, 0.0, W[0], nmo);
            if (nv > 0)
                C_DGEMM('N', 'T', fv, fv, nv, -1.0, &Gp[0][fv], nmo, &Vp[0][fv], nmo, 1.0,
                        W[0], nmo);
        }
        if (fv > 0 && nv > 0) {
            // W_OV = Ag Bv
            C_DGEMM('N', 'N', fv, nv, fv, 1.0, Gp[0], nmo, &Vp[0][fv], nmo, 0.0, &W[0][fv], nmo);
            // W_VO = -Bg^T Av
            C_DGEMM('T', 'N', nv, fv, fv, -1.0, &Gp[0][fv], nmo, Vp[0], nmo, 0.0, &W[fv][0], nmo);
        }
    }
    double** Wp = Wmat->pointer(0);
    // Antisymmetrise over the blocks pack_K_to_kappa reads; (V,V) stays zero (W_VV unused).
    for (int p = 0; p < fv; ++p) {
        Wp[p][p] = 0.0;
        for (int q = p + 1; q < nmo; ++q) {
            const double w = Wp[p][q] - Wp[q][p];
            Wp[p][q] = +w;
            Wp[q][p] = -w;
        }
    }
    std::vector<double> chart;
    reks::OrbitalDIISEngine::pack_K_to_kappa(Wmat, chart, active_mo_indices_, Ncore_, nmo);

    Hv_out.resize(n);
    const double inv = vnorm / h;
    for (int i = 0; i < n; ++i) Hv_out[i] = inv * (g_plus[i] - g0[i]) - 0.5 * chart[i];
    return true;
}

double REKS::gvb_step_norm(SharedMatrix Ca_prev, const SharedMatrix& SC) const {
    const int N = Ca_->colspi()[0];
    if (!Ca_prev || Ca_prev->colspi()[0] != N || Ca_prev->rowspi()[0] != Ca_->rowspi()[0])
        return -1.0;
    // U = Ca_prev^T (S Ca_) is the MO-basis map between the two iterates; its skew part is the
    // first-order generator of the movement, exact while U is near identity. Assumes
    // SC = S*Ca_ for the current Ca_.
    // Every non-redundant rotation pair carries at least one occupied index, so only the
    // occupied rows and columns of U are formed: P = U[O, :] and Q = U[:, O] over
    // O = [0, n_occ). A pair with both indices in O appears in both orders inside P; a pair
    // with one virtual index appears once and counts twice.
    const int nso = Ca_->rowspi()[0];
    const int n_occ = active_mo_indices_.empty() ? Ncore_ : active_mo_indices_.back() + 1;
    double** Cprev = Ca_prev->pointer(0);
    double** SCp = SC->pointer(0);
    std::vector<double> P(static_cast<size_t>(n_occ) * N);
    std::vector<double> Q(static_cast<size_t>(N) * n_occ);
    C_DGEMM('T', 'N', n_occ, N, nso, 1.0, Cprev[0], N, SCp[0], N, 0.0, P.data(), N);
    C_DGEMM('T', 'N', N, n_occ, nso, 1.0, Cprev[0], N, SCp[0], N, 0.0, Q.data(), n_occ);

    double sum = 0.0;
    for (int p = 0; p < n_occ; ++p) {
        const double* Pp = P.data() + static_cast<size_t>(p) * N;
        for (int q = 0; q < n_occ; ++q) {
            if (p == q || (p < Ncore_ && q < Ncore_)) continue;
            const double d = 0.5 * (Pp[q] - Q[static_cast<size_t>(q) * n_occ + p]);
            sum += d * d;
        }
        for (int q = n_occ; q < N; ++q) {
            const double d = 0.5 * (Pp[q] - Q[static_cast<size_t>(q) * n_occ + p]);
            sum += 2.0 * d * d;
        }
    }
    // |K1|_F/sqrt(2) is the norm of the equivalent packed non-redundant step.
    return std::sqrt(0.5 * sum);
}

std::vector<double> REKS::nk_fon_kkt_multipliers() const {
    const auto flags = nk_fon_active_set();
    std::vector<double> out;
    for (const auto& blk : fon_blocks_) {
        const auto gv = reks::REKSGradientEngine::compute_fon_gradient(
            blk.s, blk.gen, sa_cassette_, reel_.fon_state[blk.s], reel_.E_L);
        out.insert(out.end(), gv.begin(), gv.end());
    }
    for (size_t i = 0; i < out.size(); ++i)
        if (i >= flags.size() || flags[i] == 0) out[i] = 0.0;
    return out;
}

bool REKS::nk_run_episode(const NKCurvature& curv, const std::vector<double>& denom_live) {
    const int n = nk_n_rot();
    const int nmo = Ca_->colspi()[0];
    // Frobenius-norm gate equivalent to the RMS gate gvb_d_conv_ (gorb_rms(x) = x/nmo).
    const double gate_fro = gvb_d_conv_ * static_cast<double>(nmo);

    auto norm2 = [](const std::vector<double>& a) {
        return C_DNRM2(a.size(), const_cast<double*>(a.data()), 1);
    };

    // Preconditioner M = 2|denom| over all four blocks, floored at 1e-3 of its largest entry:
    // the same shaped denominator the engine divides its own Newton step by, so H k = -g is
    // solved in the sign the step is taken in.
    std::vector<double> minv(n);
    double d_max = 0.0;
    for (double d : denom_live) d_max = std::max(d_max, std::abs(d));
    for (int i = 0; i < n; ++i)
        minv[i] = 1.0 / (2.0 * std::max(std::abs(denom_live[i]), 1e-3 * d_max));

    // Episode frame: its own anchor and its own kappa, applied by the episode itself through
    // the packing and the Pade rotation the matvec uses. Ca_ref_, orbital_diis_cum_kappa_ and
    // orbital_runtime_active_ are untouched, so the frame is the same on both DIIS
    // formulations. The anchor moves to every accepted point, so the frame origin kappa_ep
    // stays at zero and is the point the outer loop linearizes at.
    nk_Ca_anchor_ = Ca_->clone();
    SharedMatrix nk_Ca_entry = Ca_->clone();
    const std::vector<double> kappa_ep(n, 0.0);
    nk_grad_evals_ = 0;
    const char* abort_reason = "";

    auto grad_at = [&](const std::vector<double>& k, std::vector<double>& g) {
        ++nk_grad_evals_;
        return reduced_gradient_at(k, curv, g, nullptr);
    };

    // Every probe leaves Ca_ and the relaxed FON at the point it evaluated, so an episode that
    // does not land has to put both back: re-anchoring on the entry orbitals and taking one
    // gradient there restores the pair.
    auto discard = [&](const char* why) {
        nk_Ca_anchor_ = nk_Ca_entry;
        std::vector<double> g_back;
        grad_at(kappa_ep, g_back);
        if (reports(4))
            outfile->Printf("  [GVB-NK] iter=%d episode DISCARDED (%s, %d gradients)\n",
                            iteration_, why, nk_grad_evals_);
        return false;
    };

    // Entry gradient on the probe path: the residual, rho and the exit test then all read the
    // one map every finite difference below is taken of.
    std::vector<double> g0;
    const bool g0_ok = grad_at(kappa_ep, g0);
    double ng = norm2(g0);
    const double ng_entry = ng;

    // Finite-difference step h = 2 sqrt(sigma_g / T_scale), the one-sided optimum under
    // gradient noise (Gill, P. E. et al. Practical Optimization; Academic Press: London,
    // 1981; section 8.6). Measured here unless REKS_GVB_NK_FD_H > 0.
    double fd_h = nk_fd_h_;
    double sigma_g = 0.0, T_scale = 0.0, hv_norm = 0.0, eta_floor = 0.0;
    if (g0_ok && !(fd_h > 0.0)) {
        std::vector<double> g0b;
        if (grad_at(kappa_ep, g0b)) {
            std::vector<double> dg(n);
            for (int i = 0; i < n; ++i) dg[i] = g0b[i] - g0[i];
            sigma_g = norm2(dg);
        }
        std::mt19937_64 rng(20260727ULL);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) v[i] = uni(rng);
        const double vs = norm2(v);
        for (int i = 0; i < n; ++i) v[i] /= vs;
        const double hT = 1.0e-3;
        std::vector<double> k1(n), k2(n), g1, g2;
        for (int i = 0; i < n; ++i) {
            k1[i] = hT * v[i];
            k2[i] = 2.0 * hT * v[i];
        }
        if (grad_at(k1, g1) && grad_at(k2, g2)) {
            std::vector<double> d2(n), d1(n);
            for (int i = 0; i < n; ++i) d2[i] = (g2[i] - 2.0 * g1[i] + g0[i]) / (hT * hT);
            for (int i = 0; i < n; ++i) d1[i] = (g1[i] - g0[i]) / hT;
            T_scale = norm2(d2);
            hv_norm = norm2(d1);
        }
        const double sigma = (sigma_g > 0.0) ? sigma_g
                                             : std::numeric_limits<double>::epsilon() * ng;
        if (T_scale > 0.0) fd_h = 2.0 * std::sqrt(sigma / T_scale);
        // Relative error of one matvec at fd_h; floors the Krylov forcing term.
        if (hv_norm > 0.0) eta_floor = 2.0 * std::sqrt(sigma * T_scale) / hv_norm;
    }

    if (!g0_ok) abort_reason = "entry-probe-unusable";
    else if (!(fd_h > 0.0)) abort_reason = "fd-step-unmeasured";

    // Outer inexact-Newton loop on g(kappa) = 0 with a trust region on the step norm.
    double eta_prev = 0.5;
    double ng_prev = ng;
    double delta = 0.5;
    int macro = 0;
    int macros_run = 0;
    int n_reject = 0;
    bool model_rejected = false;
    bool accepted_any = false;

    for (; abort_reason[0] == '\0' && macro < nk_max_macro_; ++macro) {
        if (gorb_rms(ng) < nk_exit_factor_ * gvb_d_conv_) {
            abort_reason = "gate";
            break;
        }
        if (nk_grad_evals_ >= nk_max_grad_) {
            abort_reason = "max-grad";
            break;
        }

        ++macros_run;

        // Eisenstat-Walker choice 2 (Eisenstat, S. C.; Walker, H. F. SIAM J. Sci. Comput. 1996,
        // 17, 16-32, eq. 2.6 and section 3.1) with the anti-oversolve safeguard, a 0.5 ceiling
        // tightened from their 0.9, the last-step guard of Kelley, C. T. (Solving Nonlinear
        // Equations with Newton's Method; SIAM: Philadelphia, 2003) and the finite-difference
        // noise floor.
        double eta = 0.5;
        if (macro > 0) {
            const double r = ng / ng_prev;
            eta = 0.9 * r * r;
            if (0.9 * eta_prev * eta_prev > 0.1) eta = std::max(eta, 0.9 * eta_prev * eta_prev);
            eta = std::min(eta, 0.5);
        }
        eta = std::min(0.5, std::max(eta, 0.5 * gate_fro / ng));
        eta = std::max(eta, eta_floor);

        int n_unusable = 0;
        int n_starved = 0;

        // Base point of the chart correction in nk_Hv. g0 is fixed until the macro-iteration
        // accepts a step, so the unpacking is done once here instead of once per matvec.
        auto Ghat = std::make_shared<Matrix>("nk chart gradient", nmo, nmo);
        reks::OrbitalDIISEngine::unpack_kappa_to_K(g0, Ghat, active_mo_indices_, Ncore_, nmo);

        // An unusable direction or an exhausted gradient budget returns a zero Hv: minres_tr
        // reads it as no curvature and the loop breaks right after (n_starved/n_unusable below).
        auto Hv = [&](const double* v, double* out) {
            if (nk_grad_evals_ >= nk_max_grad_) {
                ++n_starved;
                for (int i = 0; i < n; ++i) out[i] = 0.0;
                return;
            }
            std::vector<double> vv(v, v + n), hv;
            const bool ok = nk_Hv(vv, g0, Ghat, fd_h, curv, hv);
            ++nk_grad_evals_;
            if (ok) {
                for (int i = 0; i < n; ++i) out[i] = hv[i];
            } else {
                ++n_unusable;
                for (int i = 0; i < n; ++i) out[i] = 0.0;
            }
        };

        reks::MinresResult res =
            reks::minres_tr(n, g0.data(), Hv, minv.data(), eta, delta, nk_max_inner_);

        std::vector<double> g_new;
        const bool try_ok = grad_at(res.k, g_new);
        const double ng_new = try_ok ? norm2(g_new) : ng;

        const double pred = ng - res.resid_2norm;
        const double act = ng - ng_new;
        const double rho = (try_ok && pred > 0.0) ? act / pred : -1.0;
        const bool accepted = (rho > 0.0);

        if (reports(4))
            outfile->Printf(
                "  [GVB-NK] iter=%d macro=%d eta=%.2e resid=%.3e nmv=%d delta=%.3e "
                "rho=%+.4f |g| %.6e -> %.6e %s\n",
                iteration_, macro, eta, res.resid_2norm, res.n_matvec, delta, rho, ng, ng_new,
                accepted ? "ACCEPT" : "reject");

        // A rejected macro-iteration means the quadratic model disagreed with the measured
        // reduction. Rejections are counted over the whole episode, not consecutively: an
        // episode whose model is wrong interleaves its rejections with weak accepts, and a
        // consecutive counter is reset by each of them (Conn, A. R. et al. Trust-Region
        // Methods; SIAM: Philadelphia, 2000; section 6.4).
        if (!accepted) ++n_reject;
        if (n_reject >= nk_max_reject_) {
            abort_reason = "model-reject";
            model_rejected = true;
            break;
        }

        ng_prev = ng;
        if (accepted) {
            // Re-anchor on the accepted point: the finite difference is taken about the frame
            // origin, so the next macro-iteration starts from a frame whose origin is the point
            // it linearizes at. Composing the rotations keeps this exact.
            nk_Ca_anchor_ = orbital_diis_.apply_rotation(nk_Ca_anchor_, res.k,
                                                         active_mo_indices_, Ncore_, nmo);
            g0 = g_new;
            ng = ng_new;
            eta_prev = eta;
            accepted_any = true;
        }

        // Two-tier expansion on a boundary step: near-exact model agreement means the radius,
        // not the model, is what limits the step, so the doubling of the standard rule costs a
        // whole macro-iteration per factor of two (Conn, A. R. et al. Trust-Region Methods;
        // SIAM: Philadelphia, 2000; section 17.1 allows [2, 4]).
        if (rho < 0.25)
            delta *= 0.5;
        else if (rho > 0.9 && res.hit_radius)
            delta *= 4.0;
        else if (rho > 0.75 && res.hit_radius)
            delta *= 2.0;
        // Delta_0 = min(||k_newton||, 0.5): the first solve runs at the standing cap, so its
        // step is the unconstrained Newton step whenever that is shorter (Conn, A. R. et al.
        // Trust-Region Methods; SIAM: Philadelphia, 2000).
        if (macro == 0 && accepted) delta = std::min(res.step_norm, 0.5);

        if (n_starved > 0) {
            abort_reason = "max-grad";
            break;
        }
        if (n_unusable > 0) {
            abort_reason = "probe-unusable";
            break;
        }
    }
    if (abort_reason[0] == '\0') abort_reason = "max-macro";

    // Closing probe at the frame origin, which re-anchoring has left standing on the last
    // accepted point (or on the entry point if nothing was accepted). Unconditional: every
    // probe overwrites reel_, G_ and the single F_gen_ workspace. Full build, not gradient-only:
    // the latter fills the per-microstate MO Fock diagonal over the active block alone, not all
    // of nmo.
    std::vector<double> g_final, denom_final;
    reduced_gradient_at(kappa_ep, curv, g_final, &denom_final, true);
    ++nk_grad_evals_;
    const double ng_final = norm2(g_final);

    if (!accepted_any) return discard(abort_reason);

    // A model the trust region kept rejecting is no basis for keeping the steps it accepted
    // before it broke: the episode is discarded whole, not landed on its last accepted point.
    if (model_rejected) return discard(abort_reason);

    // Ca_ and the relaxed FON already stand at the accepted kappa (closing probe);
    // epsilon_a_ takes the diagonal of the freshly assembled coupling Fock.
    assemble_F_reks_MO();
    {
        double** Fmo = reel_.F_reks_MO->pointer(0);
        double* eps_a = epsilon_a_->pointer(0);
        for (int i = 0; i < nmo; ++i) eps_a[i] = Fmo[i][i];
        epsilon_b_->copy(*epsilon_a_);
    }
    find_occupation();

    if (reports(4))
        outfile->Printf(
            "  [GVB-NK] iter=%d episode LANDED (%s): macros=%d gradients=%d "
            "|g|_F %.6e -> %.6e (rms %.3e -> %.3e)\n",
            iteration_, abort_reason, macros_run, nk_grad_evals_, ng_entry, ng_final,
            gorb_rms(ng_entry), gorb_rms(ng_final));
    return true;
}

void REKS::print_orbitals() {
    if (!reports(2)) return;
    if (nirrep_ > 1) {
        throw PSIEXCEPTION("REKS::print_orbitals requires c1 symmetry (nirrep_=1)");
    }
    if (!catalog_.data || Ncore_ < 0) {
        HF::print_orbitals();
        return;
    }

    const int h = 0;
    const int nmo = nmopi_[h];
    const auto labels_irrep = molecule_->irrep_labels();
    const int n_active = static_cast<int>(active_mo_indices_.size());

    std::vector<int> mo_to_active(nmo, -1);
    for (int k = 0; k < n_active; ++k) {
        const int mo = active_mo_indices_[k];
        if (mo >= 0 && mo < nmo) mo_to_active[mo] = k;
    }

    const std::vector<int> rank = mo_energy_rank();

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

    if (!active.empty()) {
        std::snprintf(hdr, sizeof(hdr), "Active Fractionally Occupied (%d):",
                      static_cast<int>(active.size()));
        outfile->Printf("    %-70s\n\n", hdr);

        // One column per (sector, generation, scheme) with an active geminal; a
        // (sector, generation) with none is skipped. A single global scheme labels
        // columns n, m, u, ...; extra schemes or sectors add tags (nt1, mt1, ...).
        // Entries outside a column's partition read NaN.
        struct FonColumn {
            int sector;
            int gen;
            int scheme;
            std::vector<double> fon;
        };
        std::vector<FonColumn> cols;
        for (int sec = 0; sec < sa_cassette_.n_sectors(); ++sec) {
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
                const auto& act = sa_cassette_.geminals_active(sec, gen);
                if (act.empty()) continue;
                std::vector<int> gen_schemes;
                for (int g : act) {
                    const int s = sa_cassette_.geminal_templates()[g].scheme;
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
                        const auto& tpl = sa_cassette_.geminal_templates()[g];
                        if (tpl.scheme != s) continue;
                        col.fon[tpl.orbitals[0]] = reel_.fon_state[sec].layers[gen][g].p;
                        col.fon[tpl.orbitals[1]] = reel_.fon_state[sec].layers[gen][g].q;
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
            for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen)
                for (int g : sa_cassette_.geminals_active_gen(gen)) {
                    const auto& tpl = sa_cassette_.geminal_templates()[g];
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
                    std::to_string(rank[active_mo_indices_[p]] + 1) + labels_irrep[h];
                const std::string r_q =
                    std::to_string(rank[active_mo_indices_[q]] + 1) + labels_irrep[h];
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
                              rk, labels_irrep[h].c_str(), orb.c_str(), e);
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
            for (int g : sa_cassette_.geminals_active_gen(0)) {
                const auto& orbs = sa_cassette_.geminal_templates()[g].orbitals;
                fon_of_active_pos[orbs[0]] = sa_fons().layers[0][g].p;
                fon_of_active_pos[orbs[1]] = sa_fons().layers[0][g].q;
            }
            const int n_pairs = static_cast<int>(sa_cassette_.geminals_active_gen(0).size());
            if (n_pairs > 1) {
                std::vector<std::pair<int, int>> pairs;
                pairs.reserve(n_pairs);
                for (int g : sa_cassette_.geminals_active_gen(0)) {
                    const auto& orbs = sa_cassette_.geminal_templates()[g].orbitals;
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
                    // Energy-rank twin of the slot-letter pair.
                    const std::string r_p = std::to_string(rank[active_mo_indices_[p]] + 1) + labels_irrep[h];
                    const std::string r_q = std::to_string(rank[active_mo_indices_[q]] + 1) + labels_irrep[h];
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
                                rk, labels_irrep[h].c_str(), lbl.c_str(), e, lbl.c_str(),
                                fon_of_active_pos[apos]);
            }
        }
        outfile->Printf("\n");
    }

    std::snprintf(hdr, sizeof(hdr), "Virtual (%d):", static_cast<int>(vir.size()));
    print_orbital_pairs(hdr, vir);
}

// rank[mo] = position of mo when all MOs are sorted by (epsilon_a_, mo_index)
// ascending, on irrep 0.
std::vector<int> REKS::mo_energy_rank() const {
    const int h = 0;
    const int nmo = nmopi_[h];
    std::vector<std::pair<double, int>> orb_e(nmo);
    for (int a = 0; a < nmo; ++a) orb_e[a] = std::make_pair(epsilon_a_->get(h, a), a);
    std::sort(orb_e.begin(), orb_e.end());
    std::vector<int> rank(nmo);
    for (int a = 0; a < nmo; ++a) rank[orb_e[a].second] = a;
    return rank;
}

// Irrep is always "A": REKS forces C1.
std::vector<std::string> REKS::active_orbital_labels() const {
    const auto rank = mo_energy_rank();
    const std::string irrep = molecule_->irrep_labels()[0];
    std::vector<std::string> out(active_mo_indices_.size());
    for (size_t k = 0; k < active_mo_indices_.size(); ++k)
        out[k] = std::to_string(rank[active_mo_indices_[k]] + 1) + irrep;
    return out;
}

void REKS::build_fon_blocks() {
    // fon_blocks_ = { (s, gen, 0, count) : count = |geminals_active(s, gen)| > 0 },
    // ordered lexicographically by (s, gen); offset is assigned later and left 0
    // here.
    fon_blocks_.clear();
    for (int s = 0; s < sa_cassette_.n_sectors(); ++s)
        for (int gen = 0; gen < sa_cassette_.n_generations(); ++gen) {
            const int count = static_cast<int>(sa_cassette_.geminals_active(s, gen).size());
            if (count == 0) continue;
            fon_blocks_.push_back(reks::studio::FonBlock{s, gen, 0, count});
        }

    // ahc_pred_fon_frac_smooth_[b] = 0.5 if fon_blocks_[b].gen == 0 else 0.0.
    // Seeded once; a re-init with the same block layout keeps the accumulated EMA.
    if (ahc_pred_fon_frac_smooth_.size() != fon_blocks_.size()) {
        ahc_pred_fon_frac_smooth_.assign(fon_blocks_.size(), 0.0);
        for (size_t b = 0; b < fon_blocks_.size(); ++b)
            ahc_pred_fon_frac_smooth_[b] = (fon_blocks_[b].gen == 0) ? 0.5 : 0.0;
    }
}

void REKS::combined_step() {
    timer_on("REKS: combined_step");
    ScopedStage _ss_cs(reports(4) ? &scf_iter_times_ : nullptr, "combined_step");

    int n_act = static_cast<int>(active_mo_indices_.size());
    int n_rot = n_act * (n_act - 1) / 2;
    // Joint vector [n_rot | FON blocks in lexicographic (s, gen) order]; each block
    // spans [blk.offset, blk.offset + blk.count). N == gh.N.
    std::vector<reks::studio::FonBlock> blocks = fon_blocks_;
    int N = n_rot;
    for (auto& blk : blocks) { blk.offset = N; N += blk.count; }
    // Per-generation FON floor: REKS_<G>_FON if set (>= 0), else open lower 0.
    // Floors are uniform across sectors, so they key on generation only.
    auto fon_lower = [&](int gen) {
        return reks_fon_lower_[gen] >= 0.0 ? reks_fon_lower_[gen] : 0.0;
    };
    // SA-sector and gen 0 (n, PPS) FON count.
    const int sa_s = sa_sector();
    const int n_fon = static_cast<int>(sa_cassette_.geminals_active_gen(0).size());

    if (reports(4) && iteration_ <= 2) {
        int n_sa = n_sa_microstates();
        std::string fon_counts;
        for (const auto& blk : blocks)
            fon_counts += " " + reks::studio::fon_channel_name(blk.gen, blk.s) +
                          "fon=" + std::to_string(blk.count);
        outfile->Printf("  [COMBINED_STEP] entry iter=%d  n_act=%d n_rot=%d%s n_sa=%d N=%d\n",
                        iteration_, n_act, n_rot, fon_counts.c_str(), n_sa, N);
        outfile->Printf("  [COMBINED_STEP]   trah_state_.initialized=%d trust_radius=%.6f bfgs_update_count=%d\n",
                        trah_state_.initialized ? 1 : 0, trah_state_.trust_radius, bfgs_update_count_);
        outfile->Printf("  [COMBINED_STEP]   guess_Ca_=%s  kappa_disabled=%d gvb_diis_active=%d\n",
                        guess_Ca_ ? "SET" : "nullptr",
                        trah_kappa_disabled_ ? 1 : 0, ctrl_state_.is_active() ? 1 : 0);
        outfile->Printf("  [E_MICRO] iter=%d", iteration_);
        { int shown = 0;
          for (int L : sa_cassette_.microstates()) {
              if (shown++ >= 8) break;
              outfile->Printf(" E[%s]=%.6f",
                              reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.E_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [C_MICRO] iter=%d", iteration_);
        { int shown = 0;
          for (int L : sa_cassette_.microstates()) {
              if (shown++ >= 8) break;
              outfile->Printf(" C[%s]=%.6f",
                              reks::report::microstate_label(sa_cassette_, L).c_str(), reel_.C_L[L]);
          } }
        outfile->Printf("\n");
        outfile->Printf("  [FON] iter=%d", iteration_);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            const auto& gem = sa_cassette_.geminals_active(blk.s, blk.gen);
            for (int k = 0; k < blk.count; ++k) {
                const auto& fon = reel_.fon_state[blk.s].layers[blk.gen][gem[k]];
                outfile->Printf(" %s-FON%d=(%.8f, %.8f)", chan.c_str(), k, fon.p, fon.q);
            }
        }
        outfile->Printf("\n");
    }

    // TR update from PREVIOUS step's prediction. Skip when predicted
    // decrease is negligible (DIIS dominates energy change -> bogus rho).
    if (trah_state_.initialized) {
        double E_SA_actual = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);

        constexpr double pred_threshold = 1e-10;
        if (trah_state_.prev_predicted > pred_threshold) {
            double actual_decrease = trah_state_.prev_energy - E_SA_actual;
            double rho = actual_decrease / trah_state_.prev_predicted;
            trah_state_.prev_rho = rho;
            if (reports(4)) {
                outfile->Printf("  [TRAH] iter=%d E_prev=%14.8f E_now=%14.8f actual=%12.8f pred=%12.8f rho=%8.4f %s\n",
                                iteration_, trah_state_.prev_energy, E_SA_actual, actual_decrease,
                                trah_state_.prev_predicted, rho, (rho < 0.0 ? "E_UP" : "accept"));
                outfile->Printf("  [TRAH] iter=%d pred_decomp: orb=%+.6f", iteration_,
                                trah_state_.prev_pred_orb);
                for (size_t bi = 0; bi < blocks.size() && bi < trah_state_.prev_pred_fon.size(); ++bi) {
                    const std::string lt = reks::studio::fon_channel_name(blocks[bi].gen, blocks[bi].s);
                    outfile->Printf(" %sfon=%+.6f cross_%s=%+.6f", lt.c_str(),
                                    trah_state_.prev_pred_fon[bi], lt.c_str(),
                                    trah_state_.prev_pred_cross[bi]);
                }
                outfile->Printf("\n");
            }
            reks::TrustRegionSolver::update_trust_radius(E_SA_actual, trah_state_);

            if (actual_decrease >= 0) {
                // AHC PI-controller state update from accepted-step feedback:
                //   EMA s <- (1-AHC_BETA_)*s + AHC_BETA_*x  for rho and the per-block
                //   prediction fractions frac = |pred_block| / sum_blocks |pred|;
                //   integral winds up by AHC_I_RATE_*(1-rho_smooth) (capped AHC_I_MAX_)
                //   when rho_smooth low, decays when rho_smooth high.
                ahc_rho_smooth_ = (1.0 - AHC_BETA_) * ahc_rho_smooth_ + AHC_BETA_ * rho;
                double pred_total = std::abs(trah_state_.prev_pred_orb) + 1e-15;
                for (size_t bi = 0; bi < blocks.size() && bi < trah_state_.prev_pred_fon.size(); ++bi)
                    pred_total += std::abs(trah_state_.prev_pred_fon[bi]);
                double orb_frac  = std::abs(trah_state_.prev_pred_orb)  / pred_total;
                ahc_pred_orb_frac_smooth_  = (1.0 - AHC_BETA_) * ahc_pred_orb_frac_smooth_  + AHC_BETA_ * orb_frac;
                for (size_t bi = 0; bi < blocks.size(); ++bi) {
                    double fon_frac = (bi < trah_state_.prev_pred_fon.size())
                        ? std::abs(trah_state_.prev_pred_fon[bi]) / pred_total : 0.0;
                    ahc_pred_fon_frac_smooth_[bi] =
                        (1.0 - AHC_BETA_) * ahc_pred_fon_frac_smooth_[bi] + AHC_BETA_ * fon_frac;
                }
                if (ahc_rho_smooth_ < 0.4) {
                    ahc_integral_ += AHC_I_RATE_ * (1.0 - ahc_rho_smooth_);
                    ahc_integral_ = std::min(ahc_integral_, AHC_I_MAX_);
                } else if (ahc_rho_smooth_ > 0.8) {
                    ahc_integral_ *= 0.7;
                }
                consecutive_bad_rho_ = 0;
            } else {
                // rho < 0: do not corrupt AHC. 2 consecutive bad -> BFGS reset.
                consecutive_bad_rho_++;
                if (consecutive_bad_rho_ >= 2 && bfgs_update_count_ > 0) {
                    bfgs_initialized_ = false;
                    bfgs_update_count_ = 0;
                    bfgs_B_.clear();
                    bfgs_prev_step_.clear();
                    consecutive_bad_rho_ = 0;
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d RESET: 2 consecutive energy increases -> diagonal H_oo\n",
                                        iteration_);
                    }
                }
            }
            // 2 consecutive kDIIS energy increases -> reset kDIIS state.
            if (kdiis_used_prev_) {
                if (actual_decrease < 0) {
                    kdiis_energy_fail_++;
                    if (kdiis_energy_fail_ >= 2) {
                        kdiis_active_ = false;
                        kdiis_count_ = 0;
                        kdiis_oldest_ = 0;
                        kdiis_stag_count_ = 0;
                        kdiis_energy_fail_ = 0;
                        if (reports(4)) {
                            outfile->Printf("  [kDIIS] iter=%d RESET: 2 consecutive energy increases\n", iteration_);
                        }
                    }
                } else {
                    kdiis_energy_fail_ = 0;
                }
            }
        } else {
            if (reports(4)) {
                outfile->Printf("  [TRAH] iter=%d pred=%.2e < threshold -- skipping TR update (DIIS-dominated)\n",
                                iteration_, trah_state_.prev_predicted);
            }
        }

        trah_state_.prev_energy = E_SA_actual;
    }

    // First-call FON seed (gen 0): READ guess kept as-is, else set to 1.0 (equal
    // occupancy); the CI-Hamiltonian-diag seed is computed for debug print only.
    if (!trah_state_.initialized) {
        const auto& sa_gem0 = sa_cassette_.geminals_active_gen(0);
        if (guess_Ca_ && has_read_fon_guess() && !read_boundary_trap_) {
            // FON already holds the restored READ guess; keep as-is.
            if (reports(4)) {
                outfile->Printf(
                    "  [TRAH] iter=%d CI-init SKIPPED: FON restored from READ guess, keeping",
                    iteration_);
                for (int k = 0; k < n_fon; ++k) {
                    outfile->Printf(" %.10f", reel_.fon_state[sa_s].layers[0][sa_gem0[k]].p);
                }
                outfile->Printf("\n");
            }
        } else if (guess_Ca_ && has_read_fon_guess() && read_boundary_trap_) {
            // READ boundary trap escape: prev step was zero (KKT + g_orb=0).
            // Reset FON to 1.0 so TRAH explores the full landscape.
            for (int k = 0; k < n_fon; ++k) {
                reel_.fon_state[sa_s].layers[0][sa_gem0[k]] = {1.0, 2.0 - 1.0};
            }
            read_boundary_trap_ = false;
            clear_read_fon_guess();
            if (reports(4)) {
                auto ci_fon = reks::fon_init_from_ci(sa_sector(), sa_cassette_, reel_.E_L, iteration_);
                outfile->Printf("  [TRAH] iter=%d BOUNDARY TRAP ESCAPE: CI-init FON (raw):", iteration_);
                for (int k = 0; k < n_fon; ++k) outfile->Printf(" %.6f", ci_fon[k]);
                outfile->Printf(" -> all set to 1.0 (reset to explore full landscape)\n");
            }
        } else {
            // Cold path: CI-init from noisy early E_L is unreliable;
            // start FON at 1.0 (equal occupancy) and let gradients find it.
            for (int k = 0; k < n_fon; ++k) {
                reel_.fon_state[sa_s].layers[0][sa_gem0[k]] = {1.0, 2.0 - 1.0};
            }
            if (reports(4)) {
                auto ci_fon = reks::fon_init_from_ci(sa_sector(), sa_cassette_, reel_.E_L, iteration_);
                outfile->Printf("  [TRAH] iter=%d CI-init FON (raw):", iteration_);
                for (int k = 0; k < n_fon; ++k) outfile->Printf(" %.6f", ci_fon[k]);
                outfile->Printf(" -> all set to 1.0 (cold path)\n");
            }
        }
    }

    auto F_gen = build_generalized_fock();

    // Active-MO ERI tiles for all canonical pairs p <= q.
    // Full JK pass over n_act*(n_act+1)/2 densities every SCF iteration;
    // gated behind report level 5.
    reks::ActiveEriTiles hess_eri;
    if (reports(5) && n_act >= 2) {
        std::vector<std::pair<int, int>> canonical_pairs;
        canonical_pairs.reserve(static_cast<size_t>(n_act) * (n_act + 1) / 2);
        for (int p = 0; p < n_act; ++p)
            for (int q = p; q < n_act; ++q) canonical_pairs.emplace_back(p, q);
        hess_eri = compute_active_mo_eri(active_mo_indices_, canonical_pairs,
                                         &scf_iter_times_);
    }

    // Joint gradient + Hessian, one FON block per active (sector, generation).
    auto gh = reks::REKSGradientEngine::compute(F_gen, reel_, sa_cassette_,
                                                reel_.F_MO_arows_a_L, reel_.F_MO_arows_b_L,
                                                reel_.E_L, reel_.C_L, active_mo_indices_,
                                                reports(4), &hess_eri,
                                                functional_->x_alpha());

    // Engine and host must agree on the joint-vector length.
    if (gh.N != N)
        throw PSIEXCEPTION(
            "combined_step: gradient engine joint size " + std::to_string(gh.N) +
            " != combined_step joint size " + std::to_string(N) +
            "; engine/host FON layout mismatch.");

    // IPR penalty Hessian (gradient already routed through F_gen).
    add_ipr_contribution_to_hessian(gh.H, N);

    if (reports(4) && lambda_ipr_ > 0.0) {
        double g_orb_max = 0.0;
        for (int i = 0; i < n_rot; ++i)
            g_orb_max = std::max(g_orb_max, std::abs(gh.g[i]));
        print_ipr_iteration_diagnostic(iteration_, g_orb_max);
    }

    // x_current = [kappa=0..., FON blocks in (s, gen) order...].
    std::vector<double> x_current(N, 0.0);
    for (const auto& blk : blocks) {
        const auto& gem = sa_cassette_.geminals_active(blk.s, blk.gen);
        for (int k = 0; k < blk.count; ++k)
            x_current[blk.offset + k] = reel_.fon_state[blk.s].layers[blk.gen][gem[k]].p;
    }

    // Bounds: kappa in [-max_kappa, max_kappa]; FON in [fon_lower(gen), 2].
    double max_kappa = 0.3;
    std::vector<double> lower(N), upper(N);
    for (int i = 0; i < n_rot; ++i) {
        lower[i] = -max_kappa;
        upper[i] = max_kappa;
    }
    for (const auto& blk : blocks) {
        const double lo = fon_lower(blk.gen);
        for (int k = 0; k < blk.count; ++k) {
            lower[blk.offset + k] = lo;
            upper[blk.offset + k] = 2.0;
        }
    }

    // Boundary stagnation detection. Requires ALL FON (every (s, gen) block)
    // pinned at boundary with gradient pointing further out.
    bool all_at_boundary = (N > n_rot);
    bool do_odb = false;
    {
        for (size_t bi = 0; bi < blocks.size() && all_at_boundary; ++bi) {
            const auto& blk = blocks[bi];
            const auto& gem = sa_cassette_.geminals_active(blk.s, blk.gen);
            for (int k = 0; k < blk.count && all_at_boundary; ++k) {
                int idx = blk.offset + k;
                double fon = reel_.fon_state[blk.s].layers[blk.gen][gem[k]].p;
                double g_fon_k = gh.g[idx];
                bool at_upper = (fon >= upper[idx] - 1e-8) && (g_fon_k < 0);
                bool at_lower = (fon <= lower[idx] + 1e-8) && (g_fon_k > 0);
                if (!at_upper && !at_lower) all_at_boundary = false;
            }
        }

        if (all_at_boundary) {
            double g_orb_sq = 0.0;
            for (int i = 0; i < n_rot; ++i) g_orb_sq += gh.g[i] * gh.g[i];
            double g_orb_norm = std::sqrt(g_orb_sq);

            if (prev_gorb_at_boundary_ > 0.0) {
                // Stagnation: orbital gradient norm dropped < 5% since last boundary iter.
                do_odb = (g_orb_norm >= 0.95 * prev_gorb_at_boundary_);
            }
            prev_gorb_at_boundary_ = g_orb_norm;
        } else {
            prev_gorb_at_boundary_ = -1.0;
        }
    }

    if (do_odb && all_at_boundary) {
        kdiis_stag_count_++;
        if (reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag=%d/%d buf=%d/%d\n", iteration_, kdiis_stag_count_, KDIIS_TRIGGER,
                            kdiis_count_, KDIIS_MAX_VEC);
        }
    } else if (!all_at_boundary) {
        if (kdiis_stag_count_ > 0 && reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag RESET: left boundary (was %d)\n", iteration_, kdiis_stag_count_);
        }
        kdiis_stag_count_ = 0;
        if (kdiis_active_) {
            kdiis_active_ = false;
            kdiis_count_ = 0;
            kdiis_oldest_ = 0;
            kdiis_energy_fail_ = 0;
            if (reports(4)) {
                outfile->Printf("  [kDIIS] iter=%d deactivated: left boundary\n", iteration_);
            }
        }
    } else {
        if (kdiis_stag_count_ > 0 && reports(4)) {
            outfile->Printf("  [kDIIS] iter=%d stag RESET: g_orb improved (was %d)\n", iteration_, kdiis_stag_count_);
        }
        kdiis_stag_count_ = 0;
    }

    // AHC: PI-controller diagonal Hessian shift from rho feedback, one lambda
    // per block (orbital, then per (sector, generation) FON block):
    //   P = max(0, 1/rho_smooth - 1),  correction = P + I (= ahc_integral_)
    //   lambda_block = mean(block's H diagonal) * correction * frac_block
    //   orbital block: signed mean (curvature stays > 0); FON block: abs mean
    //   (H_ff diagonal can go negative at a saddle)
    //   H_diag[block] += lambda_block      (applied only when correction > 1e-4)
    {
        double ahc_p = std::max(0.0, 1.0 / std::max(ahc_rho_smooth_, 0.01) - 1.0);
        double ahc_correction = ahc_p + ahc_integral_;
        if (ahc_correction > 1e-4) {
            double h_oo_mean = 0.0;
            for (int i = 0; i < n_rot; ++i) h_oo_mean += gh.H[i * N + i];
            if (n_rot > 0) h_oo_mean /= n_rot;

            double lambda_orb = h_oo_mean * ahc_correction * ahc_pred_orb_frac_smooth_;
            for (int i = 0; i < n_rot; ++i) gh.H[i * N + i] += lambda_orb;

            // One lambda per FON block: H_ff diagonal scale differs by ~O(10)
            // across generations.
            std::vector<double> lambda_fon(blocks.size(), 0.0);
            for (size_t bi = 0; bi < blocks.size(); ++bi) {
                const int nf = blocks[bi].count;
                if (nf == 0) continue;
                const int off = blocks[bi].offset;
                double h_mean = 0.0;
                for (int k = 0; k < nf; ++k)
                    h_mean += std::abs(gh.H[(off + k) * N + (off + k)]);
                h_mean /= nf;
                lambda_fon[bi] = h_mean * ahc_correction * ahc_pred_fon_frac_smooth_[bi];
                for (int k = 0; k < nf; ++k)
                    gh.H[(off + k) * N + (off + k)] += lambda_fon[bi];
            }

            if (reports(4)) {
                outfile->Printf(
                    "  [AHC] iter=%d rho_s=%.4f P=%.4f I=%.4f corr=%.4f lam_orb=%.6f",
                    iteration_, ahc_rho_smooth_, ahc_p, ahc_integral_, ahc_correction, lambda_orb);
                for (size_t bi = 0; bi < blocks.size(); ++bi)
                    outfile->Printf(" lam_%sfon=%.6f",
                                    reks::studio::fon_channel_name(blocks[bi].gen, blocks[bi].s).c_str(),
                                    lambda_fon[bi]);
                outfile->Printf("\n");
            }
        }
    }

    // BFGS quasi-Newton update for H_oo only. Activation max|g_orb| < 0.25
    // (Chaban, G. et al. Theor. Chem. Acc. 1997, 97, 88.)
    // Powell-damped rank-2; H_ff / H_of recomputed exactly.
    {
        double g_orb_max = 0.0;
        for (int i = 0; i < n_rot; ++i) g_orb_max = std::max(g_orb_max, std::abs(gh.g[i]));
        constexpr double BFGS_ACTIVATION_THRESHOLD = 0.25;
        if (g_orb_max >= BFGS_ACTIVATION_THRESHOLD) {
            bfgs_initialized_ = false;
            bfgs_update_count_ = 0;
            bfgs_B_.clear();
        }
    }
    if (n_rot > 0) {
        if (bfgs_initialized_ && !bfgs_prev_step_.empty()) {
            // s = prev kappa step, y = orbital gradient difference.
            std::vector<double> s_orb(n_rot), y_orb(n_rot);
            for (int i = 0; i < n_rot; ++i) {
                s_orb[i] = bfgs_prev_step_[i];
                y_orb[i] = gh.g[i] - bfgs_prev_g_[i];
            }

            double s_norm2 = 0.0;
            for (int i = 0; i < n_rot; ++i) s_norm2 += s_orb[i] * s_orb[i];

            // READ curvature correction: scale B_oo by rho before first update.
            if (bfgs_update_count_ == 0 && guess_Ca_ && has_read_fon_guess()) {
                double scale = std::max(trah_state_.prev_rho, 1.0);
                if (scale > 1.5) {
                    for (int i = 0; i < n_rot; ++i)
                        bfgs_B_[i * n_rot + i] *= scale;
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d pre-update H_oo scaled by rho=%.2f (READ curvature correction)\n",
                                        iteration_, scale);
                    }
                }
            }

            if (std::sqrt(s_norm2) > 1e-8) {
                std::vector<double> Bs(n_rot, 0.0);
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) Bs[i] += bfgs_B_[i * n_rot + j] * s_orb[j];

                double sTBs = 0.0, yTs = 0.0;
                for (int i = 0; i < n_rot; ++i) {
                    sTBs += s_orb[i] * Bs[i];
                    yTs += y_orb[i] * s_orb[i];
                }

                // Powell damping (preserves positive-definiteness).
                double theta = 1.0;
                if (sTBs > 1e-15 && yTs < 0.2 * sTBs) {
                    theta = 0.8 * sTBs / (sTBs - yTs);
                    for (int i = 0; i < n_rot; ++i) y_orb[i] = theta * y_orb[i] + (1.0 - theta) * Bs[i];
                    yTs = 0.0;
                    for (int i = 0; i < n_rot; ++i) yTs += y_orb[i] * s_orb[i];
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d Powell damping: theta=%.4f yTs_new=%.6f\n", iteration_, theta,
                                        yTs);
                    }
                }

                // B <- B - (Bs)(Bs)^T/(s^TBs) + (y)(y)^T/(y^Ts)
                if (sTBs > 1e-15 && yTs > 1e-15) {
                    for (int i = 0; i < n_rot; ++i)
                        for (int j = 0; j < n_rot; ++j)
                            bfgs_B_[i * n_rot + j] += -Bs[i] * Bs[j] / sTBs + y_orb[i] * y_orb[j] / yTs;
                    bfgs_update_count_++;

                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d update #%d: sTBs=%.6f yTs=%.6f theta=%.4f B_oo_diag=[",
                                        iteration_, bfgs_update_count_, sTBs, yTs, theta);
                        for (int i = 0; i < n_rot; ++i) outfile->Printf(" %.4f", bfgs_B_[i * n_rot + i]);
                        outfile->Printf(" ]\n");
                    }
                } else if (reports(4)) {
                    outfile->Printf("  [BFGS] iter=%d update SKIPPED: sTBs=%.2e yTs=%.2e\n", iteration_, sTBs, yTs);
                }
            } else if (reports(4)) {
                outfile->Printf("  [BFGS] iter=%d update SKIPPED: ||s_orb||=%.2e < 1e-8\n", iteration_,
                                std::sqrt(s_norm2));
            }

            // Use BFGS H_oo only after at least one successful update.
            if (bfgs_update_count_ > 0) {
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) gh.H[i * N + j] = bfgs_B_[i * n_rot + j];
            } else {
                for (int i = 0; i < n_rot; ++i)
                    for (int j = 0; j < n_rot; ++j) bfgs_B_[i * n_rot + j] = gh.H[i * N + j];
                constexpr double RHO_SCALE_THRESHOLD_SYNC = 2.0;
                if (guess_Ca_ && has_read_fon_guess() && trah_state_.prev_rho > RHO_SCALE_THRESHOLD_SYNC) {
                    double scale = trah_state_.prev_rho;
                    for (int i = 0; i < n_rot; ++i) {
                        bfgs_B_[i * n_rot + i] *= scale;
                        gh.H[i * N + i] *= scale;
                    }
                    if (reports(4)) {
                        outfile->Printf("  [BFGS] iter=%d re-sync H_oo scaled by rho=%.2f (READ curvature correction)\n",
                                        iteration_, scale);
                    }
                }
            }
        } else if (!bfgs_initialized_) {
            // Initialize B_oo from current H_oo. READ guess: scale by prev_rho
            // to compensate frozen-Fock H^{1e} curvature underestimate.
            bfgs_B_.resize(n_rot * n_rot);
            for (int i = 0; i < n_rot; ++i)
                for (int j = 0; j < n_rot; ++j) bfgs_B_[i * n_rot + j] = gh.H[i * N + j];

            constexpr double RHO_SCALE_THRESHOLD = 2.0;
            if (guess_Ca_ && has_read_fon_guess() && trah_state_.prev_rho > RHO_SCALE_THRESHOLD) {
                double scale = trah_state_.prev_rho;
                for (int i = 0; i < n_rot; ++i) {
                    bfgs_B_[i * n_rot + i] *= scale;
                    gh.H[i * N + i] *= scale;
                }
                if (reports(4)) {
                    outfile->Printf("  [BFGS] iter=%d initialized B_oo(%dx%d), H_oo scaled by rho=%.2f (READ curvature correction)\n",
                                    iteration_, n_rot, n_rot, scale);
                }
            } else if (reports(4)) {
                outfile->Printf("  [BFGS] iter=%d initialized B_oo(%dx%d) from diagonal Hessian\n", iteration_, n_rot,
                                n_rot);
            }
            bfgs_initialized_ = true;
        }
    }

    // Snapshot gradient BEFORE KKT elimination zeroes locked entries.
    bfgs_prev_g_.assign(gh.g.begin(), gh.g.end());

    // kDIIS ring-buffer for orbital gradient (unaffected by KKT on FON).
    if (all_at_boundary && n_rot > 0) {
        if (static_cast<int>(kdiis_g_history_.size()) < KDIIS_MAX_VEC) {
            kdiis_g_history_.resize(KDIIS_MAX_VEC);
        }
        int write_idx;
        if (kdiis_count_ < KDIIS_MAX_VEC) {
            write_idx = kdiis_count_;
            kdiis_count_++;
        } else {
            write_idx = kdiis_oldest_;
            kdiis_oldest_ = (kdiis_oldest_ + 1) % KDIIS_MAX_VEC;
        }
        kdiis_g_history_[write_idx].assign(gh.g.begin(), gh.g.begin() + n_rot);
        if (reports(5)) {
            double gn = 0;
            for (int r = 0; r < n_rot; ++r) gn += gh.g[r] * gh.g[r];
            outfile->Printf("  [kDIIS] iter=%d stored g_orb ||g||=%.6e buf=%d/%d\n", iteration_, std::sqrt(gn),
                            kdiis_count_, KDIIS_MAX_VEC);
        }
    }

    // Boundary ODB: zero all orbital x FON cross-coupling on the BFGS Hessian
    // (every FON column block, all generations).
    const int n_fon_total = N - n_rot;
    if (do_odb && all_at_boundary && n_fon_total > 0) {
        for (int i = 0; i < n_rot; ++i) {
            for (int c = n_rot; c < N; ++c) {
                gh.H[i * N + c] = 0.0;
                gh.H[c * N + i] = 0.0;
            }
        }
        if (reports(4)) {
            double g_orb_sq = 0.0;
            for (int i = 0; i < n_rot; ++i) g_orb_sq += gh.g[i] * gh.g[i];
            double g_fon_sq = 0.0;
            for (int c = n_rot; c < N; ++c) g_fon_sq += gh.g[c] * gh.g[c];
            outfile->Printf("  [AUTOPILOT] iter=%d BOUNDARY ODB: ||g_orb||=%.4e ||g_fon||=%.4e"
                            " (zeroed orbital-FON cross-coupling, all generations)\n",
                            iteration_, std::sqrt(g_orb_sq), std::sqrt(g_fon_sq));
        }
    }

    // KKT elimination: zero gradient/row/col for boundary-locked FON;
    // prevents phantom leakage through H_of into other variables. Each
    // generation honours its own lower bound (e.g. the REKS_M_FON floor on m).
    for (const auto& blk : blocks) {
        const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
        for (int k = 0; k < blk.count; ++k) {
            int idx = blk.offset + k;
            double fon_k = x_current[idx];
            double g_fon_k = gh.g[idx];
            bool locked_upper = (fon_k >= upper[idx] - 1e-8) && (g_fon_k < 0);
            bool locked_lower = (fon_k <= lower[idx] + 1e-8) && (g_fon_k > 0);
            if (locked_upper || locked_lower) {
                gh.g[idx] = 0.0;
                for (int j = 0; j < N; ++j) {
                    gh.H[idx * N + j] = 0.0;
                    gh.H[j * N + idx] = 0.0;
                }
                gh.H[idx * N + idx] = 1.0;
                if (reports(4)) {
                    outfile->Printf("  [TRAH] iter=%d KKT-eliminate %sfon%d (f=%.6f, g=%+.6f)\n",
                                    iteration_, chan.c_str(), k, fon_k, g_fon_k);
                }
            }
        }
    }

    auto step_result = reks::TrustRegionSolver::compute_step(gh.g, gh.H, lower, upper, x_current, trah_state_);

    if (reports(4) && step_result.n_active_bounds > 0) {
        outfile->Printf("  [AUTOPILOT] iter=%d BOUNDS: %d vars clamped.", iteration_, step_result.n_active_bounds);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            for (int k = 0; k < blk.count; ++k) {
                int idx = blk.offset + k;
                double fon = x_current[idx];
                double s_fon = step_result.step[idx];
                double fon_new = fon + s_fon;
                if (std::abs(s_fon) < 1e-12 && (fon >= 2.0 - 1e-8 || fon <= 1e-8)) {
                    outfile->Printf(" %sfon%d=%.4f(clamped,delta=0)", chan.c_str(), k, fon);
                } else {
                    outfile->Printf(" %sfon%d=%.4f->%.4f(delta=%+.6f)", chan.c_str(), k, fon, fon_new, s_fon);
                }
            }
        }
        double s_orb_sq = 0.0, s_fon_sq = 0.0;
        for (int i = 0; i < n_rot; ++i) s_orb_sq += step_result.step[i] * step_result.step[i];
        for (int c = n_rot; c < N; ++c) s_fon_sq += step_result.step[c] * step_result.step[c];
        double s_total = std::sqrt(s_orb_sq + s_fon_sq);
        outfile->Printf(" orb_frac=%.4f\n", s_total > 1e-15 ? std::sqrt(s_orb_sq) / s_total : 0.0);
    }

    // kDIIS: orbital-gradient DIIS on kappa only at boundary stagnation
    // (Fischer, T. H.; Almlof, J. J. Phys. Chem. 1992, 96, 9768. SO/DIIS, simplified for
    // kappa-reset basis.)
    bool kdiis_used = false;
    bool aitken_fired_step = false;
    if (kdiis_stag_count_ >= KDIIS_TRIGGER && kdiis_count_ >= KDIIS_MIN_VEC && n_rot > 0 && !trah_kappa_disabled_) {
        kdiis_active_ = true;
        int m = kdiis_count_;

        // Error-Gram B_ij = g_i^T g_j (m x m, symmetric).
        std::vector<double> B_diis(m * m, 0.0);
        for (int i = 0; i < m; ++i) {
            const int ii = (kdiis_oldest_ + i) % KDIIS_MAX_VEC;
            for (int j = i; j < m; ++j) {
                const int jj = (kdiis_oldest_ + j) % KDIIS_MAX_VEC;
                const double dot = C_DDOT(n_rot, kdiis_g_history_[ii].data(), 1,
                                          kdiis_g_history_[jj].data(), 1);
                B_diis[i * m + j] = dot;
                B_diis[j * m + i] = dot;
            }
        }

        // Bordered Pulay system, solved by C_DGESV. The repo convention carries the
        // opposite sign on the Lagrange multiplier; the coefficients are the same.
        std::vector<double> diis_c;
        double c_norm_sq = 0.0;
        bool solve_ok = reks::diis_bordered_coefficients(B_diis.data(), m, diis_c, c_norm_sq);

        // Coefficient guard local to kDIIS: DiisCore's |c|^2 bound is a different rule.
        if (solve_ok) {
            for (int i = 0; i < m; ++i) {
                if (std::abs(diis_c[i]) > 10.0) {
                    solve_ok = false;
                    break;
                }
            }
        }

        if (solve_ok) {
            // Extrapolated orbital gradient g_interp = sum_i c_i g_i.
            std::vector<double> g_interp(n_rot, 0.0);
            for (int i = 0; i < m; ++i) {
                int ii = (kdiis_oldest_ + i) % KDIIS_MAX_VEC;
                for (int r = 0; r < n_rot; ++r) {
                    g_interp[r] += diis_c[i] * kdiis_g_history_[ii][r];
                }
            }

            // H_oo kappa = -g_interp by Cholesky on the BFGS H_oo, diagonal fallback
            // when it is not positive definite.
            std::vector<double> kappa_diis(n_rot, 0.0);
            if (bfgs_update_count_ > 0) {
                // C_DPOSV overwrites both arguments and gh.H is reused below, so
                // the leading block is copied out first. The copy is row-major and
                // C_DPOSV is a column-major passthrough, so uplo='U' addresses the
                // lower triangle of gh.H.
                std::vector<double> A_oo(static_cast<size_t>(n_rot) * n_rot);
                for (int i = 0; i < n_rot; ++i)
                    C_DCOPY(n_rot, &gh.H[static_cast<size_t>(i) * N], 1,
                            &A_oo[static_cast<size_t>(i) * n_rot], 1);
                for (int i = 0; i < n_rot; ++i) kappa_diis[i] = -g_interp[i];

                const int info = C_DPOSV('U', n_rot, 1, A_oo.data(), n_rot,
                                         kappa_diis.data(), n_rot);
                if (info < 0)
                    throw PSIEXCEPTION("combined_step: C_DPOSV illegal argument, info=" +
                                       std::to_string(info));
                if (info > 0) {
                    // Not positive definite: same verdict as the leading-minor test.
                    for (int i = 0; i < n_rot; ++i) {
                        double h_ii = gh.H[i * N + i];
                        kappa_diis[i] = (std::abs(h_ii) > 1e-10) ? -g_interp[i] / h_ii : 0.0;
                    }
                }
            } else {
                for (int i = 0; i < n_rot; ++i) {
                    double h_ii = gh.H[i * N + i];
                    kappa_diis[i] = (std::abs(h_ii) > 1e-10) ? -g_interp[i] / h_ii : 0.0;
                }
            }

            for (int i = 0; i < n_rot; ++i) {
                kappa_diis[i] = std::max(-max_kappa, std::min(max_kappa, kappa_diis[i]));
            }

            double trs_kappa_norm = 0.0;
            for (int i = 0; i < n_rot; ++i) trs_kappa_norm += step_result.step[i] * step_result.step[i];
            trs_kappa_norm = std::sqrt(trs_kappa_norm);
            double trs_pred = step_result.predicted_decrease;

            for (int i = 0; i < n_rot; ++i) {
                step_result.step[i] = kappa_diis[i];
            }

            step_result.predicted_decrease =
                reks::TrustRegionSolver::predicted_decrease(gh.g.data(), gh.H.data(), step_result.step.data(), N);

            kdiis_used = true;

            if (reports(4)) {
                double kdiis_norm = 0.0;
                for (int i = 0; i < n_rot; ++i) kdiis_norm += kappa_diis[i] * kappa_diis[i];
                kdiis_norm = std::sqrt(kdiis_norm);
                double g_interp_norm = 0.0, g_orb_norm = 0.0;
                for (int i = 0; i < n_rot; ++i) {
                    g_interp_norm += g_interp[i] * g_interp[i];
                    g_orb_norm += gh.g[i] * gh.g[i];
                }
                g_interp_norm = std::sqrt(g_interp_norm);
                g_orb_norm = std::sqrt(g_orb_norm);
                outfile->Printf("  [kDIIS] iter=%d active: m=%d c=[", iteration_, m);
                for (int i = 0; i < m; ++i) outfile->Printf(" %.3f", diis_c[i]);
                outfile->Printf(" ]\n");
                outfile->Printf("  [kDIIS] iter=%d ||g_interp||=%.6e ||g_orb||=%.6e ratio=%.4f\n", iteration_,
                                g_interp_norm, g_orb_norm, g_orb_norm > 1e-15 ? g_interp_norm / g_orb_norm : 0.0);
                outfile->Printf(
                    "  [kDIIS] iter=%d ||kappa_diis||=%.6f ||kappa_trs||=%.6f "
                    "pred_diis=%.8f pred_trs=%.8f\n",
                    iteration_, kdiis_norm, trs_kappa_norm, step_result.predicted_decrease, trs_pred);
            }
        } else {
            if (reports(4)) {
                outfile->Printf("  [kDIIS] iter=%d solve FAILED (m=%d) -- using TRS step\n", iteration_, m);
            }
        }
    }

    // Aitken/secant extrapolation of the orbital step along a single soft mode.
    // When the TRS step under-shoots one near-degenerate rotation the raw step
    // sequence is ~collinear and geometric (kappa_n ~ rho*kappa_{n-1}); sum the
    // remaining geometric series: kappa_inf = kappa_n / (1 - rho).
    // rho = dot(s_n, s_{n-1}) / |s_{n-1}|^2; gated on tight collinearity
    // (cosang > 0.99) so it only completes a step to the same stationary point.
    std::vector<double> cur_raw_orb;
    if (n_rot > 0) cur_raw_orb.assign(step_result.step.begin(), step_result.step.begin() + n_rot);
    if (!kdiis_used && !trah_kappa_disabled_ && n_rot > 0 &&
        static_cast<int>(prev_orb_step_.size()) == n_rot) {
        double dot = 0.0, np2 = 0.0, nn2 = 0.0;
        for (int i = 0; i < n_rot; ++i) {
            dot += cur_raw_orb[i] * prev_orb_step_[i];
            np2 += prev_orb_step_[i] * prev_orb_step_[i];
            nn2 += cur_raw_orb[i] * cur_raw_orb[i];
        }
        if (np2 > 1e-20 && nn2 > 1e-20) {
            double rho = dot / np2;
            double cosang = dot / std::sqrt(nn2 * np2);
            if (cosang > 0.99 && rho > 0.6 && rho < 0.999) {
                double factor = std::min(1.0 / (1.0 - rho), 50.0);
                for (int i = 0; i < n_rot; ++i) step_result.step[i] *= factor;
                aitken_fired_step = true;
                if (reports(4)) {
                    outfile->Printf("  [AITKEN] iter=%d cos=%.4f rho=%.4f factor=%.2f soft-mode step extrapolation\n",
                                    iteration_, cosang, rho, factor);
                }
            }
        }
    }
    // Store the RAW (pre-extrapolation) orbital step. Reset when
    // kDIIS/DIIS-takeover broke the natural step sequence.
    if (kdiis_used || trah_kappa_disabled_) {
        prev_orb_step_.clear();
    } else if (n_rot > 0) {
        prev_orb_step_ = cur_raw_orb;
    }

    // trah_kappa_disabled (DIIS takeover): kappa zeroed, FON-only step.
    if (trah_kappa_disabled_) {
        for (int k = 0; k < n_rot; ++k) {
            step_result.step[k] = 0.0;
        }
        step_result.predicted_decrease =
            reks::TrustRegionSolver::predicted_decrease(gh.g.data(), gh.H.data(), step_result.step.data(), N);
    }

    // Damp intra-pair orbital rotations near degeneracy: damp = gap/THR.
    // Suppresses spontaneous localization at flat-surface directions.
    {
        constexpr double GAP_THRESHOLD = 0.1;
        double* eps = epsilon_a_->pointer(0);
        bool any_damped = false;

        std::vector<int> mo_to_unit(n_act, -1);
        for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
            int g = sa_cassette_.geminals_active_gen(0)[u];
            for (int orb : sa_cassette_.geminal_templates()[g].orbitals) mo_to_unit[orb] = u;
        }

        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                bool is_intra = (mo_to_unit[i] >= 0 && mo_to_unit[i] == mo_to_unit[j]);
                if (is_intra) {
                    int mo_i = active_mo_indices_[i];
                    int mo_j = active_mo_indices_[j];
                    double gap = std::abs(eps[mo_i] - eps[mo_j]);

                    if (gap < GAP_THRESHOLD) {
                        double damp = gap / GAP_THRESHOLD;
                        double kappa_before = step_result.step[idx];
                        step_result.step[idx] *= damp;
                        any_damped = true;

                        if (reports(4)) {
                            outfile->Printf(
                                "  [DELOC] iter=%d damp pair(%d,%d) gap=%.6f damp=%.4f "
                                "kappa=%.6f->%.6f\n",
                                iteration_, mo_i, mo_j, gap, damp,
                                kappa_before, step_result.step[idx]);
                        }
                    }
                }
                ++idx;
            }
        }

        if (any_damped) {
            step_result.predicted_decrease =
                reks::TrustRegionSolver::predicted_decrease(gh.g.data(), gh.H.data(), step_result.step.data(), N);
        }
    }

    // READ boundary trap: zero step + all-FON KKT-locked = wrong basin.
    if (guess_Ca_ && has_read_fon_guess() && !trah_state_.initialized) {
        double step_norm_trap = 0.0;
        for (int i = 0; i < N; ++i) step_norm_trap += step_result.step[i] * step_result.step[i];
        step_norm_trap = std::sqrt(step_norm_trap);
        if (step_norm_trap < 1e-10) {
            read_boundary_trap_ = true;
            trah_state_.prev_energy = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);
            // Do NOT mark initialized; CI-init must re-run.
            if (reports(4)) {
                outfile->Printf(
                    "  [TRAH] iter=%d BOUNDARY TRAP DETECTED: ||step||=%.2e, "
                    "will restart FON via CI-init on next iteration\n",
                    iteration_, step_norm_trap);
            }
            timer_off("REKS: combined_step");
            return;
        }
    }

    // First call: seed trah_state_ with the current SA energy baseline.
    if (!trah_state_.initialized) {
        trah_state_.prev_energy = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);
        trah_state_.initialized = true;
    }
    // Aitken-extrapolated step is off the quadratic model; zero predicted decrease
    // so next iter skips the TR/AHC update (rho from an off-model step is bogus).
    if (aitken_fired_step) step_result.predicted_decrease = 0.0;
    trah_state_.prev_predicted = step_result.predicted_decrease;

    // Decompose the quadratic model -dE = -(g^T s + 0.5 s^T H s) by subspace block:
    // orb (rotation block), one FON term per (sector, generation) block (its own
    // diagonal block), and one orbital x FON cross term per block (s_orb^T H_of
    // s_fon). Cross-block FON x FON entries are identically zero (block-diagonal
    // Hessian).
    {
        const auto& s = step_result.step;
        const auto& g = gh.g;
        const auto& H = gh.H;

        double pred_orb = 0.0;
        for (int i = 0; i < n_rot; ++i) {
            pred_orb += g[i] * s[i];
            for (int j = 0; j < n_rot; ++j) {
                pred_orb += 0.5 * s[i] * H[i * N + j] * s[j];
            }
        }

        std::vector<double> pred_fon(blocks.size(), 0.0), pred_cross(blocks.size(), 0.0);
        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            const int b = blocks[bi].offset;
            const int e = b + blocks[bi].count;
            for (int i = b; i < e; ++i) {
                pred_fon[bi] += g[i] * s[i];
                for (int j = b; j < e; ++j) {
                    pred_fon[bi] += 0.5 * s[i] * H[i * N + j] * s[j];
                }
            }
            // Cross orbital-FON predicted-decrease term: debug print only.
            if (reports(4)) {
                for (int i = 0; i < n_rot; ++i) {
                    for (int j = b; j < e; ++j) {
                        pred_cross[bi] += s[i] * H[i * N + j] * s[j];
                    }
                }
            }
        }

        // Sign flip: positive = decrease.
        step_result.pred_orb = -pred_orb;
        trah_state_.prev_pred_orb = step_result.pred_orb;
        step_result.pred_fon.assign(blocks.size(), 0.0);
        step_result.pred_cross.assign(blocks.size(), 0.0);
        trah_state_.prev_pred_fon.assign(blocks.size(), 0.0);
        trah_state_.prev_pred_cross.assign(blocks.size(), 0.0);
        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            step_result.pred_fon[bi]   = -pred_fon[bi];
            step_result.pred_cross[bi] = -pred_cross[bi];
            trah_state_.prev_pred_fon[bi]   = step_result.pred_fon[bi];
            trah_state_.prev_pred_cross[bi] = step_result.pred_cross[bi];
        }
        trah_state_.prev_step_on_boundary = step_result.step_on_boundary;
    }

    // Apply FON step, one block per (sector, generation): fon_post = clamp(fon_pre
    // + step, [0,2]). Margin guard caps the first entry into the margin zone at
    // +/- fon_bnd_margin from a bound; once inside, small steps reach it naturally.
    constexpr double fon_bnd_margin = 1e-8;
    std::vector<std::vector<double>> fon_pre(blocks.size()), fon_post(blocks.size());
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        const auto& blk = blocks[bi];
        const auto& gem = sa_cassette_.geminals_active(blk.s, blk.gen);
        fon_pre[bi].resize(blk.count);
        fon_post[bi].resize(blk.count);
        for (int k = 0; k < blk.count; ++k) {
            fon_pre[bi][k] = x_current[blk.offset + k];
            double new_fon = fon_pre[bi][k] + step_result.step[blk.offset + k];
            new_fon = std::max(0.0, std::min(2.0, new_fon));
            if (new_fon >= 2.0 - fon_bnd_margin && fon_pre[bi][k] < 2.0 - fon_bnd_margin) {
                new_fon = 2.0 - fon_bnd_margin;
            }
            if (new_fon <= fon_bnd_margin && fon_pre[bi][k] > fon_bnd_margin) {
                new_fon = fon_bnd_margin;
            }
            fon_post[bi][k] = new_fon;
            reel_.fon_state[blk.s].layers[blk.gen][gem[k]] = {new_fon, 2.0 - new_fon};
        }
    }

    // Apply orbital rotation Ca_ <- Ca_ * exp(K) (antisymmetric K).
    // K is nonzero only in the active x active block, so exp(K) is identity
    // outside it and the rotation only mixes the n_act active columns of Ca_.
    // Work the n_act x n_act block directly: O(nso*n_act^2), not O(nso^3).
    if (n_rot > 0) {
        int nso = Ca_->rowspi()[0];

        // Persistent scratch (run-constant shapes).
        if (!reel_.trah_K)
            reel_.trah_K = std::make_shared<Matrix>("trah K", n_act, n_act);
        if (!reel_.trah_Ca_act)
            reel_.trah_Ca_act = std::make_shared<Matrix>("trah Ca_act", nso, n_act);
        if (!reel_.trah_Ca_act_rot)
            reel_.trah_Ca_act_rot = std::make_shared<Matrix>("trah Ca_act_rot", nso, n_act);

        double** Kp = reel_.trah_K->pointer(0);
        reel_.trah_K->zero();

        int idx = 0;
        double max_kappa = 0.0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                double kappa_ij = step_result.step[idx];
                Kp[i][j] = kappa_ij;
                Kp[j][i] = -kappa_ij;
                if (std::abs(kappa_ij) > max_kappa) max_kappa = std::abs(kappa_ij);
                ++idx;
            }
        }

        reel_.trah_K->expm(4, true);  // K_act -> U = exp(K_act)

        // Gather active columns, rotate (Ca_act * U), scatter back.
        double** Cap = Ca_->pointer(0);
        double** Aa = reel_.trah_Ca_act->pointer(0);
        for (int i = 0; i < n_act; ++i) {
            int ai = active_mo_indices_[i];
            for (int mu = 0; mu < nso; ++mu) Aa[mu][i] = Cap[mu][ai];
        }

        double** Rp = reel_.trah_Ca_act_rot->pointer(0);
        C_DGEMM('N', 'N', nso, n_act, n_act, 1.0, Aa[0], n_act, Kp[0], n_act, 0.0, Rp[0], n_act);

        for (int j = 0; j < n_act; ++j) {
            int aj = active_mo_indices_[j];
            for (int mu = 0; mu < nso; ++mu) Cap[mu][aj] = Rp[mu][j];
        }

        if (reports(4)) {
            outfile->Printf("  [TRAH] iter=%d Applied orbital rotation: max|kappa|=%.6f\n", iteration_, max_kappa);
        }
    }

    {
        bfgs_prev_step_.resize(N);
        for (int i = 0; i < n_rot; ++i) bfgs_prev_step_[i] = step_result.step[i];
        for (size_t bi = 0; bi < blocks.size(); ++bi)
            for (int k = 0; k < blocks[bi].count; ++k)
                bfgs_prev_step_[blocks[bi].offset + k] = fon_post[bi][k] - fon_pre[bi][k];
    }

    kdiis_used_prev_ = kdiis_used;

    if (reports(4)) {
        int it = iteration_;

        outfile->Printf("  [TRAH] iter=%d", it);
        if (trah_kappa_disabled_) outfile->Printf(" (kappa disabled)");
        if (kdiis_used)
            outfile->Printf(" kDIIS=ON(stag=%d,buf=%d)", kdiis_stag_count_, kdiis_count_);
        else if (kdiis_stag_count_ > 0)
            outfile->Printf(" kDIIS=arming(%d/%d)", kdiis_stag_count_, KDIIS_TRIGGER);
        outfile->Printf("\n");

        double g_norm = 0.0;
        for (int i = 0; i < N; ++i) g_norm += gh.g[i] * gh.g[i];
        g_norm = std::sqrt(g_norm);

        double step_norm = 0.0;
        for (int i = 0; i < N; ++i) step_norm += step_result.step[i] * step_result.step[i];
        step_norm = std::sqrt(step_norm);

        double step_over_tr = (trah_state_.trust_radius > 1e-15) ? step_norm / trah_state_.trust_radius : 0.0;
        outfile->Printf("  [TRAH] iter=%d ||g||=%10.6f ||step||=%10.6f TR=%6.3f step/TR=%.4f mu=%10.6f pred=%12.8f", it,
                        g_norm, step_norm, trah_state_.trust_radius, step_over_tr, step_result.mu,
                        step_result.predicted_decrease);
        if (step_result.n_active_bounds > 0) outfile->Printf(" bounds=%d", step_result.n_active_bounds);
        outfile->Printf("\n");

        // Step-gradient alignment per subspace (cos=-1 steepest descent, ~0 TR fighting Newton).
        {
            double dot_orb = 0.0, sn_orb = 0.0, gn_orb = 0.0;
            for (int i = 0; i < n_rot; ++i) {
                dot_orb += step_result.step[i] * gh.g[i];
                sn_orb  += step_result.step[i] * step_result.step[i];
                gn_orb  += gh.g[i] * gh.g[i];
            }
            double dot_total = 0.0;
            for (int i = 0; i < N; ++i) dot_total += step_result.step[i] * gh.g[i];
            double denom_orb   = std::sqrt(sn_orb) * std::sqrt(gn_orb);
            double denom_total = step_norm * g_norm;
            outfile->Printf("  [TRAH] iter=%d cos(step,g): total=%+.4f orb=%+.4f", it,
                            denom_total > 1e-15 ? dot_total / denom_total : 0.0,
                            denom_orb   > 1e-15 ? dot_orb   / denom_orb   : 0.0);
            for (const auto& blk : blocks) {
                if (blk.count == 0) continue;
                double dot = 0.0, sn = 0.0, gn = 0.0;
                for (int i = blk.offset; i < blk.offset + blk.count; ++i) {
                    dot += step_result.step[i] * gh.g[i];
                    sn  += step_result.step[i] * step_result.step[i];
                    gn  += gh.g[i] * gh.g[i];
                }
                double denom = std::sqrt(sn) * std::sqrt(gn);
                outfile->Printf(" %sfon=%+.4f",
                                reks::studio::fon_channel_name(blk.gen, blk.s).c_str(),
                                denom > 1e-15 ? dot / denom : 0.0);
            }
            outfile->Printf("\n");
        }

        std::vector<std::string> labels(N);
        {
            int idx = 0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    labels[idx] = "k(" + std::to_string(active_mo_indices_[i]) + "," +
                                  std::to_string(active_mo_indices_[j]) + ")";
                    ++idx;
                }
            }
            for (const auto& blk : blocks) {
                const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
                const auto& gem = sa_cassette_.geminals_active(blk.s, blk.gen);
                for (int k = 0; k < blk.count; ++k) {
                    const auto& orbs = sa_cassette_.geminal_templates()[gem[k]].orbitals;
                    labels[blk.offset + k] = chan + "(" +
                        std::to_string(active_mo_indices_[orbs[0]]) + "," +
                        std::to_string(active_mo_indices_[orbs[1]]) + ")";
                }
            }
        }

        int lw = 0;
        for (const auto& lb : labels) lw = std::max(lw, static_cast<int>(lb.size()));
        lw = std::max(lw, 10);

        outfile->Printf("  [TRAH] iter=%d g = [", it);
        for (int i = 0; i < N; ++i) {
            outfile->Printf(" %10.6f", gh.g[i]);
        }
        outfile->Printf(" ]\n");

        outfile->Printf("  [TRAH] iter=%d H_diag = [", it);
        for (int i = 0; i < N; ++i) {
            outfile->Printf(" %10.6f", gh.H[i * N + i]);
        }
        outfile->Printf(" ]\n");

        outfile->Printf("  [TRAH] iter=%d h_oo_intra:", it);
        for (int k = 0; k < n_fon; ++k) {
            const auto& orbs = sa_cassette_.geminal_templates()[sa_cassette_.geminals_active_gen(0)[k]].orbitals;
            int p = orbs[0], q = orbs[1];
            int idx = p * (2 * n_act - p - 1) / 2 + (q - p - 1);
            outfile->Printf(" unit%d=%.6f", k, gh.H[idx * N + idx]);
        }
        outfile->Printf("\n");

        if (bfgs_initialized_ && n_rot > 1) {
            outfile->Printf("  [BFGS_DIAG] iter=%d B_oo off-diag per row:", it);
            for (int i = 0; i < n_rot; ++i) {
                double max_off = 0.0;
                for (int j = 0; j < n_rot; ++j)
                    if (i != j) max_off = std::max(max_off, std::abs(bfgs_B_[i * n_rot + j]));
                outfile->Printf(" r%d=%.4f", i, max_off);
            }
            outfile->Printf("\n");

            outfile->Printf("  [BFGS_DIAG] iter=%d step_diag_only = [", it);
            for (int i = 0; i < n_rot; ++i) {
                double hii = gh.H[i * N + i];
                double s_diag = (std::abs(hii) > 1e-15) ? -gh.g[i] / hii : 0.0;
                outfile->Printf(" %+.6f", s_diag);
            }
            outfile->Printf(" ]\n");
            outfile->Printf("  [BFGS_DIAG] iter=%d step_actual   = [", it);
            for (int i = 0; i < n_rot; ++i) {
                outfile->Printf(" %+.6f", step_result.step[i]);
            }
            outfile->Printf(" ]\n");
        }

        outfile->Printf("  [TRAH] iter=%d h_of_cross:", it);
        for (const auto& blk : blocks) {
            const std::string chan = reks::studio::fon_channel_name(blk.gen, blk.s);
            for (int k = 0; k < blk.count; ++k) {
                int fon_idx = blk.offset + k;
                double max_cross = 0.0;
                for (int r = 0; r < n_rot; ++r) {
                    double v = std::abs(gh.H[r * N + fon_idx]);
                    if (v > max_cross) max_cross = v;
                }
                outfile->Printf(" %sfon%d_max=%.6f", chan.c_str(), k, max_cross);
            }
        }
        outfile->Printf("\n");

        if (!gh.hoo_evals_raw.empty()) {
            outfile->Printf("  [TRAH] iter=%d H_oo evals:", it);
            for (size_t i = 0; i < gh.hoo_evals_raw.size(); ++i) {
                outfile->Printf(" %.6f", gh.hoo_evals_raw[i]);
            }
            outfile->Printf("\n");
        }

        // Per (sector, generation) block: FON-Hessian eigenvalue spectrum, then the block matrix.
        for (size_t bi = 0; bi < gh.fon_blocks.size(); ++bi) {
            const auto& fblk = gh.fon_blocks[bi];
            const std::string chan = reks::studio::fon_channel_name(fblk.gen, fblk.s);
            if (!gh.fon_hess_evals[bi].empty()) {
                outfile->Printf("  [TRAH] iter=%d H_%sf evals:", it, chan.c_str());
                for (double ev : gh.fon_hess_evals[bi]) outfile->Printf(" %.6f", ev);
                outfile->Printf("\n");
            }
            if (gh.fon_hess_block[bi].empty()) continue;
            const int nf  = fblk.count;
            const int off = fblk.offset;
            const std::vector<double>& hblk = gh.fon_hess_block[bi];
            outfile->Printf("  [TRAH] iter=%d H_%sf: %*s", it, chan.c_str(), lw, "");
            for (int j = 0; j < nf; ++j) {
                outfile->Printf(" %*s", lw, labels[off + j].c_str());
            }
            outfile->Printf("\n");
            for (int i = 0; i < nf; ++i) {
                outfile->Printf("  [TRAH] iter=%d H_%sf: %*s", it, chan.c_str(), lw, labels[off + i].c_str());
                for (int j = 0; j < nf; ++j) {
                    outfile->Printf(" %*.6f", lw, hblk[i * nf + j]);
                }
                outfile->Printf("\n");
            }
        }

        if (N <= 25) {
            outfile->Printf("  [TRAH] iter=%d H: %*s", it, lw, "");
            for (int j = 0; j < N; ++j) {
                outfile->Printf(" %*s", lw, labels[j].c_str());
            }
            outfile->Printf("\n");
            for (int i = 0; i < N; ++i) {
                outfile->Printf("  [TRAH] iter=%d H: %*s", it, lw, labels[i].c_str());
                for (int j = 0; j < N; ++j) {
                    outfile->Printf(" %*.6f", lw, gh.H[i * N + j]);
                }
                outfile->Printf("\n");
            }
        }

        outfile->Printf("  [TRAH] iter=%d kappa:", it);
        for (int i = 0; i < n_rot; ++i) outfile->Printf(" %s=%+.6f", labels[i].c_str(), step_result.step[i]);
        outfile->Printf("\n");

        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            const auto& blk = blocks[bi];
            if (blk.count == 0) continue;
            outfile->Printf("  [TRAH] iter=%d %s-FON:", it,
                            reks::studio::fon_channel_name(blk.gen, blk.s).c_str());
            for (int k = 0; k < blk.count; ++k) {
                int idx = blk.offset + k;
                outfile->Printf(" %s=%.6f%+.6f->%.6f(delta=%+.6f)",
                                labels[idx].c_str(), fon_pre[bi][k], step_result.step[idx],
                                fon_post[bi][k], fon_post[bi][k] - fon_pre[bi][k]);
            }
            outfile->Printf("\n");
        }
    }

    timer_off("REKS: combined_step");
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

    Dold_ = SharedMatrix(factory_->create_matrix("D old"));
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
            outfile->Printf("  [REKS] %s lower bound = %.8f (experimental limiter active)\n",
                            opt.c_str(), floor);
    }
    // Micro-Newton FON solver defaults; lower_bound and collect_history are placeholders,
    // replaced per-generation and per-call before each solve.
    fon_micro_cfg_ = reks::MicroSolverConfig();
    fon_micro_cfg_.max_nr_iter     = options_.get_int   ("REKS_FON_MICRO_MAX_NR_ITER");
    fon_micro_cfg_.max_ls_iter     = options_.get_int   ("REKS_FON_MICRO_MAX_LS_ITER");
    fon_micro_cfg_.convergence_tol = options_.get_double("REKS_FON_MICRO_CONV_TOL");
    fon_micro_cfg_.min_eigenvalue  = options_.get_double("REKS_FON_MICRO_MIN_EIGENVALUE");
    fon_micro_cfg_.lower_bound     = options_.get_double("REKS_FON_MICRO_BOUND_MARGIN");
    fon_micro_cfg_.upper_bound     = 2.0 - fon_micro_cfg_.lower_bound;
    fon_micro_cfg_.two_pass_ls     = true;
    fon_micro_cfg_.collect_history = false;
    fon_micro_boundary_skin_ = options_.get_double("REKS_FON_MICRO_BOUNDARY_SKIN");

    if (nirrep_ != 1)
        reks::report::input_error(
            "REKS requires C1 symmetry, but the molecule has " + std::to_string(nirrep_) +
            " irreps. Add 'symmetry c1' to the molecule.", false);

    lambda_ipr_ = options_.get_double("REKS_DELOC_IPR_PENALTY");
    ipr_method_ = options_.get_str   ("REKS_DELOC_IPR_METHOD");
    if (lambda_ipr_ < 0.0)
        reks::report::input_error("REKS_DELOC_IPR_PENALTY must be >= 0.", false);

    ipr_prev_ = -1.0;
    ipr_prevprev_ = -1.0;
    ipr_guard_streak_ = 0;
    ipr_guard_streak_start_ipr_ = -1.0;
    ipr_guard_disabled_ = false;
    ipr_guard_streak_limit_ = options_.get_int("REKS_DELOC_IPR_GUARD_STREAK_LIMIT");

    // S^{1/2} + AO->atom maps; built unconditionally since PR diagnostics print
    // regardless of lambda_ipr_.
    build_ipr_static_data();
    if (lambda_ipr_ > 0.0 && reports(4)) {
        outfile->Printf("\n  === REKS IPR Penalty ===\n");
        outfile->Printf("  [IPR] init: method=%s lambda=%.3e nso=%d n_atoms=%d\n",
                        ipr_method_.c_str(), lambda_ipr_,
                        static_cast<int>(ao2atom_.size()), n_atoms_ipr_);
    }

    orbital_guard_.set_adaptive(options_.get_bool("LEVEL_SHIFT_ADAPT"));

    if (!options_["REKS"].has_changed()) {
        reks::report::input_error(
            "REKS active space not specified: set the 'reks' option to [N, M].");
    }
    // Active space [N, M]. Spin manifolds are declared through SA_REKS_2SPIN /
    // SI_REKS_2SPIN or the ladder order of nested config blocks; the molecule
    // multiplicity is always canonical (even N -> 1, odd N -> 2).
    if (options_["REKS"].size() == 3) {
        reks::report::input_error(
            "REKS: the spin element was removed; 'reks' is [N, M] (2 elements). "
            "Select spin manifolds with SA_REKS_2SPIN / SI_REKS_2SPIN "
            "(e.g. a triplet-manifold run: sa_reks_2spin [2]).");
    }
    if (options_["REKS"].size() != 2) {
        reks::report::input_error(
            "REKS: the 'reks' option must have exactly 2 elements [N, M], got " +
            std::to_string(options_["REKS"].size()) + ".");
    }
    const int reks_N = static_cast<int>(options_["REKS"][0].to_double());
    const int reks_M = static_cast<int>(options_["REKS"][1].to_double());
    const int expected_mult = (reks_N % 2 == 0) ? 1 : 2;
    if (multiplicity_ != expected_mult) {
        reks::report::input_error(
            "REKS requires the canonical symmetric-Ms reference: N=" + std::to_string(reks_N) +
            " implies molecule multiplicity " + std::to_string(expected_mult) +
            " (even N -> 1, odd N -> 2), but the molecule multiplicity is " +
            std::to_string(multiplicity_) + ". Do not set the multiplicity to 2S+1; keep "
            "the canonical value and pick spin manifolds via SA_REKS_2SPIN / SI_REKS_2SPIN.");
    }
    // Build the SCF reference at symmetric Ms for odd-electron active spaces.
    const bool use_symmetric_proxy = (reks_N % 2 != 0);

    if (!reks::studio::variant_supported(reks_N, reks_M)) {
        reks::report::input_error(
            "REKS active space [N=" + std::to_string(reks_N) + ", M=" + std::to_string(reks_M) +
            "] is not installed. Installed [N,M] active spaces: " +
            reks::studio::supported_variants_str() + ".");
    }
    timer_on("REKS: make_catalog");
    // Load the merged (N,M) blob bound to its lowest manifold. Declared 2S labels
    // resolve against the sector table; the run binding is installed once the run
    // set is known.
    catalog_ = reks::studio::make_catalog(reks_N, reks_M);
    timer_off("REKS: make_catalog");

    // Core = Ncore doubly occupied orbitals holding nelectron - N electrons:
    // require (nelectron - N) even and >= 0, and Ncore + M <= nmo.
    if ((nelectron_ - reks_N) % 2 != 0) {
        reks::report::input_error(
            "REKS active space mismatch: total electrons (" + std::to_string(nelectron_) +
            ") minus active electrons N=" + std::to_string(reks_N) + " = " +
            std::to_string(nelectron_ - reks_N) +
            " is odd, so the closed core cannot be doubly occupied. Check N in 'reks'.");
    }
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
    auto blob_view = [this](int b) {
        reks::studio::Catalog c = catalog_;
        c.active_sector = b;
        return c;
    };

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

    // Form probe: the first non-empty element fixes the nesting depth; every other
    // non-empty element must match it. Returns -1 when all elements are empty lists
    // (or the option itself is an empty list).
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
        const int count = cv.sector_config_count();
        std::vector<int> Ks = parse_config_indices(list, cv, /*allow_exclude=*/false);
        std::set<int> seen;
        for (int K : Ks) {
            if (K < start || K >= start + count) {
                reks::report::input_error(
                    "SA_REKS_CONFIGS: config index (local) " + std::to_string(K - start) +
                    " is out of range [0, " + std::to_string(count) + ") for " + where + ".");
            }
            if (!seen.insert(K).second) {
                reks::report::input_error(
                    "SA_REKS_CONFIGS: duplicate config (local) " + std::to_string(K - start) +
                    " in " + where + " (each config at most once per block).");
            }
        }
        return Ks;
    };
    // Parse+validate one SI cassette's local-index/name list into global K indices.
    // A cassette is the basis of one SI subspace, so it is a SET: an entry repeated
    // by overlapping bulk tokens or by hand keeps its first position and drops the
    // rest. The SA pool stays a weighted list, where a repeat is an error.
    auto parse_cassette_si = [&](Data& list, int b, const std::string& where) {
        const reks::studio::Catalog cv = blob_view(b);
        const int start = cv.sector_config_start();
        const int count = cv.sector_config_count();
        const std::vector<int> parsed = parse_config_indices(list, cv, /*allow_exclude=*/true);
        std::set<int> seen;
        std::vector<int> Ks;
        Ks.reserve(parsed.size());
        for (int K : parsed) {
            if (K < start || K >= start + count) {
                reks::report::input_error(
                    "SI_REKS_CONFIGS: config index (local) " + std::to_string(K - start) +
                    " is out of range [0, " + std::to_string(count) + ") for " + where + ".");
            }
            if (seen.insert(K).second) Ks.push_back(K);
        }
        return Ks;
    };

    // Typed SA blocks bound to blob sectors, in typed order. A ladder-mode empty
    // block skips its rung; an explicitly labeled block declares its manifold even
    // when empty.
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
        const reks::studio::Catalog cv = blob_view(b);
        const int* sa_def = cv.sector_sa_default(cv.base, cv.active_sector);
        SaBlockBind bind{b, {}};
        for (int i = 0; i < cv.sector_n_sa_default(); ++i) bind.Ks.push_back(sa_def[i]);
        sa_blocks.push_back(std::move(bind));
    } else {
        sa_default_deferred = true;
    }

    // SI groups bound to blob sectors. 2-level = one implicit group of cassettes;
    // 3-level = one group per outer entry. An explicitly labeled empty group
    // declares its manifold; ladder mode declares only via content.
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

    // RUN SET: manifolds declared by SA blocks and SI groups; empty -> the lowest
    // available manifold. Run ordinals ascend in 2S (std::set iterates ascending,
    // and the blob sector table is ascending in 2S).
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
    catalog_.active_sector = run_to_blob_[0];

    // Deferred defaults land on run sector 0 (the run's lowest manifold).
    if (sa_default_deferred) {
        const reks::studio::Catalog cv = catalog_sector_view(0);
        const int* sa_def = cv.sector_sa_default(cv.base, cv.active_sector);
        SaBlockBind bind{run_to_blob_[0], {}};
        for (int i = 0; i < cv.sector_n_sa_default(); ++i) bind.Ks.push_back(sa_def[i]);
        sa_blocks.push_back(std::move(bind));
    }
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

    // SA pool in canonical order: blocks sorted by ascending manifold, config order
    // preserved within each block. weight_perm maps each canonical config position
    // to its typed position so explicit SA_REKS_WEIGHTS (given as typed) follow.
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

    // SA extra determinants: runtime Slater determinants added to the SA ensemble
    // (orbital optimization only; absent from SI). One microstate each, occupation
    // vector [alpha(M) | beta(M)].
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

    // SI cassettes in canonical group order (ascending manifold), cassette order
    // preserved within each group. si_lists[i] holds cassette i's global-K configs;
    // si_sectors[i] its run sector for FON dispatch.
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

    // Symmetric-Ms proxy: the reference occupation would silently ignore user DOCC/SOCC.
    if ((input_docc_ || input_socc_) && use_symmetric_proxy) {
        reks::report::input_error(
            "DOCC/SOCC cannot be combined with a symmetric-Ms REKS reference: "
            "occupations are fixed by the symmetric proxy. Remove DOCC/SOCC.");
    }
    if (use_symmetric_proxy) find_occupation();

    use_trah_ = options_.get_bool("REKS_USE_TRAH");
    if (use_trah_) {
        trah_state_ = reks::TRAHState();
        trah_state_.trust_radius = 0.1;
    }

    use_gvb_diis_ = options_.get_bool("REKS_GVB_DIIS");

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
        reks::DiisConfig dcfg;
        dcfg.max_vectors       = options_.get_int("DIIS_MAX_VECS");
        dcfg.level_shift       = options_.get_double("REKS_GVB_LEVEL_SHIFT");
        dcfg.cond_angle_filter = options_.get_bool("REKS_DIIS_COND_ANGLE_FILTER");
        dcfg.angle_tol         = options_.get_double("REKS_DIIS_COND_ANGLE_TOL");
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
        dcfg.guard_landing_rtol = options_.get_double("REKS_DIIS_GUARD_LANDING_RTOL");
        dcfg.mon_progress_eps  = options_.get_double("REKS_DIIS_MON_PROGRESS_EPS");
        dcfg.fon_branch_tol    = options_.get_double("REKS_DIIS_FON_BRANCH_TOL");
        dcfg.restart_budget    = options_.get_int("REKS_DIIS_RESTART_BUDGET");
        dcfg.mon_cycle         = options_.get_bool("REKS_DIIS_MON_CYCLE");
        dcfg.cycle_window      = options_.get_int("REKS_DIIS_MON_CYCLE_WINDOW");
        dcfg.cycle_cap_streak  = options_.get_int("REKS_DIIS_MON_CYCLE_CAP_STREAK");
        dcfg.orb_base_map_damp = options_.get_double("REKS_DIIS_ORB_BASE_MAP_DAMP");
        gvb_diis_ = reks::GVBDIISEngine(dcfg);
        gvb_diis_.reset();
        gvb_diis_start_ = options_.get_int("REKS_GVB_DIIS_START");
        diis_orb_base_map_damp_ = dcfg.orb_base_map_damp;
        diis_controller_.configure(dcfg);
        diis_controller_.bind(&cfm_adapter_, &orb_adapter_);
        gvb_d_conv_ = options_.get_double("D_CONVERGENCE");
        gvb_e_conv_ = options_.get_double("E_CONVERGENCE");
        gvb_scf_maxiter_ = options_.get_int("MAXITER");
        gvb_gate_streak_req_ = options_.get_int("REKS_GVB_GATE_STREAK");
        nk_enabled_      = options_.get_bool("REKS_GVB_NK");
        nk_trigger_      = options_.get_int("REKS_GVB_NK_TRIGGER");
        nk_min_landings_ = options_.get_int("REKS_GVB_NK_MIN_LANDINGS");
        nk_max_reject_   = options_.get_int("REKS_GVB_NK_MAX_REJECT");
        nk_fon_micro_    = options_.get_int("REKS_GVB_NK_FON_MICRO");
        nk_require_active_fon_ = options_.get_bool("REKS_GVB_NK_REQUIRE_ACTIVE_FON");
        nk_max_grad_     = options_.get_int("REKS_GVB_NK_MAX_GRAD");
        nk_max_macro_    = options_.get_int("REKS_GVB_NK_MAX_MACRO");
        nk_max_inner_    = options_.get_int("REKS_GVB_NK_MAX_INNER");
        nk_fd_h_         = options_.get_double("REKS_GVB_NK_FD_H");
        nk_exit_factor_  = options_.get_double("REKS_GVB_NK_EXIT_FACTOR");
        nk_trigger_streak_ = 0;
        nk_trigger_streak_max_ = 0;
        nk_episode_done_ = false;
        nk_grad_evals_ = 0;
        gvb_prev_Ca_.reset();
        gvb_step_norm_ = -1.0;
        gvb_prev_step_norm_ = -1.0;
        gvb_prev_gorb_ = -1.0;
        // The adjudication band is the quadratic-model energy of one Newton step at
        // curvature guard_band_B, the guard's own threshold.
        outcome_guard_.configure(options_.get_double("E_CONVERGENCE"), dcfg.fon_branch_tol,
                                 dcfg.guard_band_B, dcfg.guard_latch_grace, dcfg.restart_budget,
                                 dcfg.mon_progress_eps, dcfg.guard_landing_rtol,
                                 dcfg.mon_nonmonotone_m);

        // ORBITAL (default; Ionova, I. V.; Carter, E. A. J. Chem. Phys. 1995, 102, 1251;
        // Sethio, D. et al. J. Phys. Chem. A 2024, 128, 2472) vs CFM (Muller, R. P. et al.
        // J. Chem. Phys. 1994, 100, 1226).
        std::string diis_form = options_.get_str("REKS_DIIS_FORMULATION");
        diis_formulation_ = (diis_form == "ORBITAL") ? DIISFormulation::ORBITAL
                                                     : DIISFormulation::CFM;
        orbital_runtime_active_ = false;
        orbital_cap_fired_this_iter_ = false;
        if (diis_formulation_ == DIISFormulation::ORBITAL) {
            if (reports(4)) {
                outfile->Printf("  GVB-DIIS formulation: ORBITAL (Ionova & Carter 1995; Sethio et al. 2024) "
                                "with limit-cycle re-anchored restart\n");
            }
            orbital_diis_ = reks::OrbitalDIISEngine(dcfg);
            orbital_diis_.reset();
            orbital_diis_cum_kappa_.clear();
            Ca_ref_.reset();
        }
    }

    active_mo_indices_.resize(catalog_.data->n_active_orbitals);
    for (int i = 0; i < catalog_.data->n_active_orbitals; ++i)
        active_mo_indices_[i] = Ncore_ + i;

    reks::studio::allocate_reel(reel_, catalog_, nsopi_, nso_, n_extra);

    si_cassettes_.clear();

    // SA cassette: the SCF ensemble pool with the sector-less extras at the tail.
    // Source vectors are not read again, so move them in.
    sa_cassette_ = reks::studio::make_cassette(
        catalog_, std::move(sa_Ks), std::move(sa_weights),
        reks::studio::CassetteRole::SA,
        std::move(extra_microstates), std::move(extra_weights));

    // Multi-sector: install the host's declared-sector partition. positions[s] are
    // the contiguous positions into K_indices of sector s's config block (in the
    // sector-block-first concatenation order); an SA-empty sector has an empty list.
    // Single-sector keeps make_cassette's construction default (bit-identical).
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

    // One FON snapshot per run sector; each entry's populated layers span its
    // sector's generations (shared n_geminals). Reset to closed-shell (1,1) here,
    // including SA-empty sectors.
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

    allocate_reks_matrices();
    if (reports(4)) print_memory_footprint("after allocate_reks_matrices");

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

void REKS::allocate_reks_matrices() {
    int n_sa = n_sa_microstates();
    int n_base = (1 << sa_cassette_.n_active_orbitals());

    if (static_cast<int>(reel_.base_density_e1.size()) != n_base) {
        throw PSIEXCEPTION(
            "REKS aliasing: reel_.base_density_e1 size " +
            std::to_string(reel_.base_density_e1.size()) +
            " != n_base " + std::to_string(n_base));
    }
    if (static_cast<int>(reel_.F_MO_arows_a_L.size()) < n_sa) {
        throw PSIEXCEPTION(
            "REKS aliasing: reel_.F_MO_arows_a_L size " +
            std::to_string(reel_.F_MO_arows_a_L.size()) +
            " < n_sa_microstates " + std::to_string(n_sa));
    }
}

void REKS::print_memory_footprint(const std::string& tag) const {
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

    const int nso = nsopi_[0];
    const int n_sa = n_sa_microstates();
    const int n_micro = sa_cassette_.n_total_microstates();  // per-L reel arrays are this long
    const int n_active = sa_cassette_.n_active_orbitals();
    const int n_base = (1 << sa_cassette_.n_active_orbitals());

    auto diag_elems = [](const std::vector<std::vector<double>>& v) -> std::size_t {
        std::size_t n = 0;
        for (const auto& d : v) n += d.size();
        return n;
    };
    const std::size_t e_F_AO_a = mat_elems(reel_.F_alpha_AO_buffer);
    const std::size_t e_F_MO_a = vec_elems(reel_.F_MO_arows_a_L)
                               + diag_elems(reel_.F_MO_diag_a_L);
    const std::size_t e_F_MO_b = vec_elems(reel_.F_MO_arows_b_L)
                               + diag_elems(reel_.F_MO_diag_b_L);
    const std::size_t e_F_acc  = mat_elems(reel_.F_acc_AO) + mat_elems(reel_.F_acc_MO);
    const std::size_t e_Vxc_a  = vec_elems(reel_.Vxc_arows_a_L) + diag_elems(reel_.Vxc_diag_a_L)
                               + mat_elems(reel_.V_acc_xc_a);
    const std::size_t e_Vxc_b  = vec_elems(reel_.Vxc_arows_b_L) + diag_elems(reel_.Vxc_diag_b_L)
                               + mat_elems(reel_.V_acc_xc_b);
    const std::size_t e_Dcore  = mat_elems(reel_.D_core);
    const std::size_t e_J_act  = vec_elems(reel_.J_act);
    const std::size_t e_K_act  = vec_elems(reel_.K_act);
    const std::size_t e_wK_act = vec_elems(reel_.wK_act);
    const std::size_t e_J_core = mat_elems(reel_.J_core);
    const std::size_t e_K_core = mat_elems(reel_.K_core);
    const std::size_t e_wK_core= mat_elems(reel_.wK_core);
    const std::size_t e_Freks  = mat_elems(reel_.F_reks);
    const std::size_t e_Twork  = mat_elems(reel_.temp_buffer);
    const std::size_t e_Ca     = mat_elems(Ca_);
    const std::size_t e_Cb     = mat_elems(Cb_);
    const std::size_t e_Da     = mat_elems(Da_);
    const std::size_t e_Db     = mat_elems(Db_);
    const std::size_t e_H      = mat_elems(H_);
    const std::size_t e_S      = mat_elems(S_);

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
                    nso, Ncore_, n_active, n_sa, n_micro, n_base);
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

    row_mat("reel_.F_alpha_AO_buffer",  reel_.F_alpha_AO_buffer, e_F_AO_a);
    row_vec("F_MO_arows_a (MO)",     reel_.F_MO_arows_a_L, e_F_MO_a);
    row_vec("F_MO_arows_b (MO)",     reel_.F_MO_arows_b_L, e_F_MO_b);
    row_mat("reel_.F_acc_MO",        reel_.F_acc_MO,       e_F_acc);
    row_vec("Vxc_arows_a (MO)",   reel_.Vxc_arows_a_L, e_Vxc_a);
    row_vec("Vxc_arows_b (MO)",   reel_.Vxc_arows_b_L, e_Vxc_b);
    row_mat("reel_.D_core",            reel_.D_core, e_Dcore);
    row_vec("reel_.J_act",              reel_.J_act,  e_J_act);
    row_vec("reel_.K_act",              reel_.K_act,  e_K_act);
    row_vec("reel_.wK_act",             reel_.wK_act, e_wK_act);
    row_mat("reel_.J_core",             reel_.J_core,  e_J_core);
    row_mat("reel_.K_core",             reel_.K_core,  e_K_core);
    row_mat("reel_.wK_core",            reel_.wK_core, e_wK_core);
    row_mat("reel_.F_reks",             reel_.F_reks,  e_Freks);
    row_mat("reel_.temp_buffer",          reel_.temp_buffer, e_Twork);
    row_mat("Ca_",                 Ca_, e_Ca);
    row_mat("Cb_",                 Cb_, e_Cb);
    row_mat("Da_",                 Da_, e_Da);
    row_mat("Db_",                 Db_, e_Db);
    row_mat("H_",                  H_,  e_H);
    row_mat("S_",                  S_,  e_S);
    outfile->Printf("    ---------------------------------------------------------------------\n");
    outfile->Printf("    %-22s  %6s  %20s  %12.2f\n", "TOTAL", "", "", elems_to_mib(e_total));
    outfile->Printf("\n");
}

// The occupied MO block is [0, nalphapi_[0]); its first Ncore_ columns already sit
// in D_core. nalphapi_[0] - Ncore_ is the active electron pair count; zero pairs
// skips the update below.
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
    ScopedStage _ss(reports(4) ? &scf_iter_times_ : nullptr,
                    "build_base_densities");
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
    // F_reks_MO is the MO-space coupling Fock (nmo x nmo); sized lazily here at
    // its first touch, since nmo is set only after form_Shalf.
    if (!reel_.F_reks_MO || reel_.F_reks_MO->coldim(0) != nmo) {
        reel_.F_reks_MO = std::make_shared<Matrix>("F_REKS_MO", nmo, nmo);
    }

    reks::scf::assemble_F_reks_MO(
        reel_, sa_cassette_,
        reel_.C_L, active_mo_indices_, Ncore_, nmo, reel_.F_reks_MO);
}

void REKS::compute_weighting_factors() {
    timer_on("REKS: compute_weighting_factors");
    ScopedStage _ss_cwf(reports(4) ? &scf_iter_times_ : nullptr, "compute_weighting_factors");
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
    ScopedStage _ss(reports(4) ? &scf_iter_times_ : nullptr,
                    "build_sa_focks");

    if (!C_occ_cache_) {
        throw PSIEXCEPTION(
            "build_sa_focks: C_occ_cache_ not populated; form_D() "
            "must run first");
    }

    // Refreshes C_L from the current FON.
    reks::scf::compute_weights(reel_, sa_cassette_, reel_.C_L);

    {
        ScopedStage _ss_jk(reports(4) ? &scf_iter_times_ : nullptr, "jk_compute_sa");
        reks::scf::build_jk_cache(
            sa_cassette_, reel_, Ncore_, active_mo_indices_, Ca_, C_occ_cache_,
            *jk_, *functional_, J_, K_, wK_, G_, need_full_mo_diag());
    }

    // Closed-shell RV grid sweep: needed for pure-GGA virt-virt XC (folds Va_
    // into G_, x_alpha < 0.01) or the VV10 nonlocal energy vv10_E_rv_. A
    // hybrid without VV10 reads neither, so the sweep is skipped.
    const bool rv_G_gga = functional_->is_gga() && functional_->x_alpha() < 0.01;
    const bool rv_xc_live = functional_->needs_xc() && rv_potential_ &&
        (rv_G_gga || functional_->needs_vv10());
    if (rv_xc_live) {
        timer_on("REKS: xc_compute (sa rv)");
        ScopedStage _ss_rv(reports(4) ? &scf_iter_times_ : nullptr, "xc_compute_sa_rv");
        rv_potential_->set_D({Da_});
        rv_potential_->compute_V({Va_});
        vv10_E_rv_ = functional_->needs_vv10()
                         ? rv_potential_->quadrature_values()["VV10"]
                         : 0.0;
        if (rv_G_gga) G_->add(Va_);
        timer_off("REKS: xc_compute (sa rv)");
    } else {
        vv10_E_rv_ = 0.0;
        // This functional never sweeps RV; its point workers rebuild lazily if one is asked for.
        if (rv_potential_) rv_potential_->release_point_workers();
    }

    // Snapshot Ca once per SCF iter so this iteration's XC build and the
    // AO->MO transform below share one MO basis.
    if (!reel_.Ca_fock_transform ||
        reel_.Ca_fock_transform->rowspi() != Ca_->rowspi() ||
        reel_.Ca_fock_transform->colspi() != Ca_->colspi()) {
        reel_.Ca_fock_transform = Ca_->clone();
    } else {
        reel_.Ca_fock_transform->copy(Ca_);
    }

    // need_diag gates the per-L V_xc MO diagonal and its full-MO projection;
    // only GVB-DIIS needs it, the default path reads only the active rows.
    // The Newton-Krylov gradient probe reads neither (F_gen comes from the
    // active rows and F_acc_MO).
    {
        ScopedStage _ss_uv(reports(4) ? &scf_iter_times_ : nullptr, "xc_compute_sa_uv");
        reks::scf::build_microstate_xc(
            reel_, sa_cassette_, sa_cassette_.microstates(), Ncore_, *functional_,
            uv_potential_, reel_.Ca_fock_transform, active_mo_indices_, reel_.C_L,
            /*need_diag=*/need_full_mo_diag(),
            /*mo_lo=*/0, /*mo_width=*/reel_.Ca_fock_transform->colspi()[0],
            /*want_v_acc=*/true);
    }

    // Copies UV's converged RHO_A/RHO_B into RV's quadrature_values, since
    // the RV branch above may have skipped its own compute_V pass.
    if (functional_->needs_xc() && uv_potential_ && rv_potential_) {
        const auto& uvq = uv_potential_->quadrature_values();
        const auto it = uvq.find("RHO_A");
        if (it != uvq.end()) {
            auto& rvq = rv_potential_->quadrature_values();
            rvq["RHO_A"] = it->second;
            rvq["RHO_B"] = uvq.at("RHO_B");
        }
    }

    {
        ScopedStage _ss_mo(reports(4) ? &scf_iter_times_ : nullptr, "build_sa_focks_MO");
        // Ca_fock_transform is the copy of Ca_ taken above, and D_core is
        // C_core C_core^T over that same Ca_, so it is the core density of the
        // transform basis.
        reks::scf::build_sa_focks_MO(
            reel_, sa_cassette_, sa_cassette_.microstates(), Ncore_, *functional_, *jk_, H_,
            reel_.D_core, /*need_full_diag=*/need_full_mo_diag(),
            /*mo_lo=*/0, /*mo_width=*/reel_.Ca_fock_transform->colspi()[0]);
    }

    timer_off("REKS: build_sa_focks");
}

// Post-SCF state interaction: builds SI-only base densities/Focks, relaxes
// per-sector SI FON, then for each SI cassette solves the generalized
// eigenproblem H_si c = E S_si c over its configurations K.
void REKS::compute_si() {
    if (si_computed_) return;
    timer_on("REKS: compute_si");
    // SCF is finished here; flush the SCF [TIME] summary once.
    if (reports(4) && !scf_summary_printed_) {
        log_scf_summary_();
        scf_summary_printed_ = true;
        post_scf_times_.clear();
        post_scf_start_time_ = std::chrono::steady_clock::now();
    }
    ScopedStage _ss_csi(reports(4) ? &post_scf_times_ : nullptr, "compute_si");

    si_results_.assign(si_cassettes_.size(), reks::SIResult{});
    const bool needs_xc = functional_->needs_xc();

    // form_D built only the SA-referenced patterns; build the SI-only ones
    // here from Ca_ (unchanged since convergence) -- bit-identical to a
    // per-iteration form_D build.
    if (!call_sheet_.si_only_bitmask.empty()) {
        reks::scf::build_base_densities(sa_cassette_, call_sheet_.si_only_bitmask,
                                        reel_, Ncore_, active_mo_indices_, Ca_, H_);
    }

    // SCF cached per-L Focks only for the SA cassette's L set; SI Hamiltonians
    // also need L's outside it. Reuses the SCF JK cache and the Ca snapshot
    // in Ca_fock_transform.
    {
        ScopedStage _ss_xsi(reports(4) ? &post_scf_times_ : nullptr, "xc_compute_si");
        reks::scf::fill_missing_si_focks(
            reel_, sa_cassette_, si_cassettes_, Ncore_, *functional_,
            uv_potential_, reel_.Ca_fock_transform, active_mo_indices_, *jk_, H_);
    }

    // E_L[L] is FON-independent and L-keyed, so filling it per real cassette
    // is dedup-safe. Driven by the real cassettes only -- the union cassette's
    // cross-cassette pairs reference microstates with no per-L Fock.
    for (const reks::studio::Cassette& cassette : si_cassettes_) {
        reks::scf::build_si_microstate_focks(
            reel_, cassette, sa_cassette_,
            nuclearrep_, needs_xc, reel_.E_L);
    }

    // Per-sector SI FON: each sector's SI cassettes union into one cassette,
    // relaxed once so every real cassette in the sector shares the same FON.
    // Order-invariant (primary-K is the sector's global min, find_primary_K,
    // independent of pool partitioning). Sectors decouple exactly at frozen
    // orbitals, so per-sector relaxation is the exact block decomposition;
    // it reads only allocated E_L scalars.
    {
        timer_on("REKS: fon_relax");
        ScopedStage _ss_fr(reports(4) ? &post_scf_times_ : nullptr, "fon_relax");
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

    // Sorted, deduplicated union of call_sheet_.eri_pairs and, when sa_ran
    // and n_act >= 2, every canonical (p < q) active-MO pair; std::set
    // enforces the (k <= l) order compute_active_mo_eri requires. Empty
    // union -> empty pool, no JK pass.
    {
        std::set<std::pair<int,int>> pair_set(call_sheet_.eri_pairs.begin(),
                                              call_sheet_.eri_pairs.end());
        const auto& sa_micro = sa_cassette_.microstates();
        const bool sa_ran = !sa_micro.empty() &&
            std::abs(reel_.E_L[sa_micro.front()]) >= reks::constants::ENERGY_THRESHOLD;
        const int n_act = static_cast<int>(active_mo_indices_.size());
        if (sa_ran && n_act >= 2)
            for (int p = 0; p < n_act; ++p)
                for (int q = p + 1; q < n_act; ++q) pair_set.emplace(p, q);
        active_eri_pool_ = compute_active_mo_eri(
            active_mo_indices_,
            std::vector<std::pair<int,int>>(pair_set.begin(), pair_set.end()),
            &post_scf_times_);
    }

    // SI H/S built from the catalog's coupling tables via the cassette's pair rows.
    for (size_t k = 0; k < si_cassettes_.size(); ++k) {
        const reks::studio::Cassette& cassette = si_cassettes_[k];

        std::vector<double> H_si;
        reks::BlockedMatrix S_si;
        int n_out = 0, n_overlap = 0;
        {
            timer_on("REKS: compute_SI_energies");
            ScopedStage _ss_cse(reports(4) ? &post_scf_times_ : nullptr, "compute_SI_energies");
            {
                // selector_packs is sized to this cassette, decoded here and released
                // with this scope; serial, timed apart from the parallel build_hamiltonian
                // pass below.
                ScopedStage _ss_sd(reports(4) ? &post_scf_times_ : nullptr,
                                   "si_selector_demand");
                timer_on("REKS: si_selector_demand");
                const reks::studio::SelectorDemand demand =
                    reks::studio::selector_demand(cassette);
                timer_off("REKS: si_selector_demand");
                _ss_sd.release();
                ScopedStage _ss_sp(reports(4) ? &post_scf_times_ : nullptr,
                                   "si_selector_packs");
                timer_on("REKS: si_selector_packs");
                const reks::studio::SelectorPacks selector_packs =
                    reks::studio::decode_selector_packs(
                        catalog_sector_view(cassette.sector()), demand);
                timer_off("REKS: si_selector_packs");
                _ss_sp.release();
                ScopedStage _ss_bh(reports(4) ? &post_scf_times_ : nullptr,
                                   "si_build_hamiltonian");
                timer_on("REKS: si_build_hamiltonian");
                reks::si::build_hamiltonian(
                    reel_, cassette, selector_packs, reel_.E_L,
                    reel_.lagrangians, functional_->x_alpha(), H_si, n_out,
                    &active_eri_pool_);
                timer_off("REKS: si_build_hamiltonian");
            }
            {
                ScopedStage _ss_bo(reports(4) ? &post_scf_times_ : nullptr,
                                   "si_build_overlap");
                timer_on("REKS: si_build_overlap");
                reks::si::build_overlap(
                    reel_, cassette, S_si, n_overlap);
                timer_off("REKS: si_build_overlap");
            }
            {
                ScopedStage _ss_dg(reports(4) ? &post_scf_times_ : nullptr,
                                   "si_diagonalize");
                si_results_[k] = reks::si::diagonalize(
                    cassette, std::move(H_si), std::move(S_si), n_out);
            }
            timer_off("REKS: compute_SI_energies");
        }
    }

    si_computed_ = true;
    timer_off("REKS: compute_si");
    _ss_csi.release();
    if (reports(4)) log_post_scf_summary_("compute_si");
}

reks::ActiveEriTiles REKS::compute_active_mo_eri(
    const std::vector<int>& active_mo_indices,
    const std::vector<std::pair<int,int>>& rank1_pairs,
    std::map<std::string, double>* stage_bucket) {
    const int N = static_cast<int>(active_mo_indices.size());
    if (N == 0 || rank1_pairs.empty()) {
        return reks::ActiveEriTiles(N, {});
    }

    // Validate canonical sorted unique input (codegen invariant).
    for (size_t p = 0; p < rank1_pairs.size(); ++p) {
        const auto& kv = rank1_pairs[p];
        if (kv.first < 0 || kv.second < kv.first || kv.second >= N) {
            throw PSIEXCEPTION(
                "compute_active_mo_eri: pair (" + std::to_string(kv.first) + "," +
                std::to_string(kv.second) + ") out of range for N=" + std::to_string(N));
        }
        if (p > 0 && !(rank1_pairs[p - 1] < kv)) {
            throw PSIEXCEPTION(
                "compute_active_mo_eri: pairs must be canonical (k<=l) and "
                "lexicographically strictly sorted; got duplicate or unsorted entry");
        }
    }

    // Pair-driven active-space MO ERI build. For each canonical pair (k<=l):
    //   k == l : 1 density (c_k)
    //   k <  l : 2 densities, u_+ = (c_k+c_l)/sqrt(2), u_- = (c_l-c_k)/sqrt(2)
    // All densities for all pairs batch into one jk_->compute() call.
    // R(u) = C_active^T J(u) C_active; tile(k,l) = R(u) if k==l, else (R(u_+) - R(u_-))/2.
    timer_on("REKS: compute_active_mo_eri");
    reks::ActiveEriTiles result(N, rank1_pairs);

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
    {
        timer_on("REKS: jk_compute (active_eri)");
        ScopedStage _ss_jk(reports(4) ? stage_bucket : nullptr, "jk_compute_active_eri");
        jk_->compute();
        timer_off("REKS: jk_compute (active_eri)");
    }

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

void REKS::compute_gamma_vectors(std::vector<double>& gamma_aa,
                                 std::vector<double>& gamma_ca) const {
    const int n_act = static_cast<int>(active_mo_indices_.size());
    const int Mn = reel_.slot_diag_stride;

    // J(p,k) = (pp|kk), K(p,k) = (pk|pk): MO-diagonal slices of the J_act/K_act
    // Fock slots (reks_reel.h), over MO index p and active index k.
    if (reel_.slot_diag_K < 0)
        throw PSIEXCEPTION(
            "compute_gamma_vectors: no plain (pk|pk) available. A wcombine JK folds "
            "exact exchange into wK, leaving no K to read; set WCOMBINE false.");
    if (Mn != Ca_->colspi()[0])
        throw PSIEXCEPTION(
            "compute_gamma_vectors: the SA Fock MO transform ran windowed this iteration, "
            "so the slot diagonals cover no full MO extent.");
    auto J_pk = [&](int p, int k) {
        return reel_.slot_diag[static_cast<size_t>(reel_.slot_diag_J + k) * Mn + p];
    };
    auto K_pk = [&](int p, int k) {
        return reel_.slot_diag[static_cast<size_t>(reel_.slot_diag_K + k) * Mn + p];
    };

    // loc < Ncore -> core (1,1); else active occupations.
    auto get_occ = [&](int L, int loc) -> std::pair<int, int> {
        if (loc < Ncore_) return {1, 1};
        int act = loc - Ncore_;
        return {sa_cassette_.microstate(L).alpha[act], sa_cassette_.microstate(L).beta[act]};
    };

    // Diagonal of the 2e orbital Hessian on the pair (i,j), SA-averaged over the
    // microstates with weight C_L. Setting r2 = r1 in A2e (reks_gradient_engine.h)
    // collapses its ERIs onto the pair quantities held here, (ij|ij) = K_ij and
    // (ii|jj) = J_ij, so with dna = na_i - na_j and dnb = nb_i - nb_j:
    //   coul = K_ij,  exch = x_alpha (J_ij + K_ij) / 2
    //   gamma_ij = 2 sum_L C_L [ (dna^2 + dnb^2)(coul - exch) + 2 dna dnb coul ]
    // coul carries full weight for every functional; exact exchange enters the energy
    // at x_alpha and leaves the curvature scaled the same way. Not represented: the
    // erf-attenuated exchange of a range-separated functional (no (pp|kk)_omega slot),
    // and the XC kernel -- the same limits A2e itself carries.
    const double x_alpha = functional_ ? functional_->x_alpha() : 1.0;
    auto compute_gamma = [&](int i_loc, int j_loc, double J_ij, double K_ij) -> double {
        const double coul = K_ij;
        const double exch = 0.5 * x_alpha * (J_ij + K_ij);
        double h2e = 0.0;

        for (int L : sa_cassette_.microstates()) {
            const double CL = reel_.C_L[L];
            if (std::abs(CL) < 1e-14) continue;

            auto [na_i, nb_i] = get_occ(L, i_loc);
            auto [na_j, nb_j] = get_occ(L, j_loc);
            const double dna = na_i - na_j;
            const double dnb = nb_i - nb_j;

            h2e += CL * ((dna * dna + dnb * dnb) * (coul - exch) + 2.0 * dna * dnb * coul);
        }

        return 2.0 * h2e;
    };

    int n_rot_aa = n_act * (n_act - 1) / 2;
    gamma_aa.resize(n_rot_aa);
    {
        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                gamma_aa[idx++] = compute_gamma(Ncore_ + i, Ncore_ + j,
                                                J_pk(Ncore_ + i, j), K_pk(Ncore_ + i, j));
            }
        }
    }

    // Core orbitals carry (na, nb) = (1, 1) in every microstate, so the SA-averaged
    // coefficients of a core-active pair depend on the active index alone: one core row
    // fixes them for all of them, and only the pair ERIs vary with c.
    gamma_ca.resize(Ncore_ * n_act);
    if (Ncore_ > 0) {
        for (int k = 0; k < n_act; ++k) {
            const double g_J = compute_gamma(0, Ncore_ + k, 1.0, 0.0);
            const double g_K = compute_gamma(0, Ncore_ + k, 0.0, 1.0);
            for (int c = 0; c < Ncore_; ++c)
                gamma_ca[c * n_act + k] = g_J * J_pk(c, k) + g_K * K_pk(c, k);
        }
    }

    if (reports(5)) {
        outfile->Printf("  [GVB-DIIS] iter=%d gamma_aa (%d pairs):", iteration_, n_rot_aa);
        for (int i = 0; i < n_rot_aa; ++i) outfile->Printf(" %.6f", gamma_aa[i]);
        outfile->Printf("\n");
        if (Ncore_ > 0) {
            outfile->Printf("  [GVB-DIIS] iter=%d gamma_ca (%d pairs): min=%.6f max=%.6f\n", iteration_, Ncore_ * n_act,
                            *std::min_element(gamma_ca.begin(), gamma_ca.end()),
                            *std::max_element(gamma_ca.begin(), gamma_ca.end()));
        }
    }
}

void REKS::compute_sa_energies() {
    timer_on("REKS: compute_sa_energies");
    ScopedStage _ss_csae(reports(4) ? &scf_iter_times_ : nullptr, "compute_sa_energies");
    // E_L = tr(Da H) + tr(Db H) + 1/2 [tr(Da Fa) + tr(Db Fb) - E_1e] + E_nuc
    //       + (E_xc - 1/2 tr(D V_xc))   (DFT correction)

    double E_nuc = nuclearrep_;
    bool needs_xc = functional_->needs_xc();
    double alpha = functional_->x_alpha();

    if (reports(5)) {
        outfile->Printf("\n  Microstate energy computation:\n");
        outfile->Printf("  E_nuc = %.10f, needs_xc = %d, alpha = %.4f\n", E_nuc, needs_xc, alpha);
    }

    // alpha_idx/beta_idx (from base_idx) key the base-density scalar table;
    // E_Fa_L/E_Fb_L, reel_.E_L and the XC scalars are indexed directly by L.
    const bool apply_xc = needs_xc && uv_potential_;
    for (int L : sa_cassette_.microstates()) {
        const int alpha_idx = base_idx(L, /*alpha=*/true);
        const int beta_idx  = base_idx(L, /*alpha=*/false);
        const auto e = reks::scf::microstate_energy(
            reel_, L, alpha_idx, beta_idx, E_nuc, apply_xc);
        reel_.E_L[L] = e.total;

        if (apply_xc && reports(5)) {
            const double tr_D_Vxc = reel_.tr_D_Vxc_a_L[L] + reel_.tr_D_Vxc_b_L[L];
            outfile->Printf("  L=%d: E_1e=%.10f, E_2e=%.10f, E_nuc=%.10f\n", L, e.E_1e, e.E_2e, E_nuc);
            outfile->Printf("        E_L(before XC)=%.10f\n", e.E_1e + e.E_2e + E_nuc);
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
        // FON pairs (fon_p, fon_q) per generation: n (PPS), m (OSS), u (DOSS), ...
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

// GVB-DIIS FON relaxation. First entry: per-sector multi-start NR over
// {(1,...,1), CI-init}, keep the lowest-E_SA interior candidate (boundary basins
// rejected). Thereafter: per-(sector, generation) joint NR from the current FON
// snapshot.
void REKS::gvb_fon_step() {
    timer_on("REKS: gvb_fon_step");
    ScopedStage _ss_gfs(reports(4) ? &scf_iter_times_ : nullptr, "gvb_fon_step");
    gvb_fon_step_swaps_ = 0;

    // A computed microstate energy carries E_nuc: a near-zero first-active-microstate
    // energy means compute_sa_energies has not run yet.
    const auto& sa_micro = sa_cassette_.microstates();
    if (sa_micro.empty() ||
        std::abs(reel_.E_L[sa_micro.front()]) < reks::constants::ENERGY_THRESHOLD) {
        timer_off("REKS: gvb_fon_step");
        return;
    }

    const int n_sec = sa_cassette_.n_sectors();
    const int n_gen = sa_cassette_.n_generations();
    const int P = static_cast<int>(sa_cassette_.geminals_active_gen(0).size());
    const auto* templ = sa_cassette_.geminal_templates();
    const auto& gen0_pool = sa_cassette_.geminals_active_gen(0);

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

    // Energy-preserving relabel of a shared orbital pair: reflect FON -> 2 - FON
    // for every (sector, generation) geminal on that pair. All sectors reflect,
    // since the swapped Ca_ columns are shared by all sectors' configs; gen-0
    // (the pair itself) reflects silently, higher generations report their
    // synced value.
    auto reflect_pair = [&](int orb0, int orb1, int pair_idx) -> int {
        int count = 0;
        for (int sec = 0; sec < n_sec; ++sec) {
            for (int gen = 0; gen < n_gen; ++gen) {
                for (int gg : sa_cassette_.geminals_active(sec, gen)) {
                    const auto& go = templ[gg].orbitals;
                    if (!((go[0] == orb0 && go[1] == orb1) ||
                          (go[0] == orb1 && go[1] == orb0))) continue;
                    const double f_new = 2.0 - reel_.fon_state[sec].layers[gen][gg].p;
                    reel_.fon_state[sec].layers[gen][gg] = {f_new, 2.0 - f_new};
                    ++count;
                    if (gen >= 1 && reports(4))
                        outfile->Printf("  [GVB-FON] iter=%d %s-geminal %d synced (-> %.6f) with pair %d swap\n",
                                        iteration_, reks::studio::fon_channel_name(gen, sec).c_str(),
                                        gg, f_new, pair_idx);
                }
            }
        }
        return count;
    };

    if (!gvb_fon_initialized_) {
#ifndef NDEBUG
        // Site-1 must reflect no FON in a single-sector run: an unconditional
        // reflect would mutate a pre-SCF-seeded higher-generation FON. Snapshot
        // before the fgap loop and re-check after.
        std::vector<std::vector<double>> dbg_fon_pre;
        if (n_sec == 1) dbg_fon_pre = capture_all_fons();
#endif
        // Pre-NR fgap inv-swap: if eps[mo_q] < eps[mo_p] for any pair, swap
        // bonding/antibonding columns in Ca_ (READ guess across state boundary).
        // Multi-sector reflects all sectors' FONs so restart-loaded higher-sector
        // FONs stay consistent with the shared column swap.
        bool any_swap = false;
        if (epsilon_a_) {
            double* eps = epsilon_a_->pointer(0);
            for (int p = 0; p < P; ++p) {
                int g = gen0_pool[p];
                const auto& orbs = templ[g].orbitals;
                int i_idx = active_mo_indices_[orbs[0]];
                int j_idx = active_mo_indices_[orbs[1]];
                double gap = eps[j_idx] - eps[i_idx];
                if (gap < 0.0) {
                    if (reports(4)) {
                        outfile->Printf(
                            "  [GVB-FON] iter=%d INV-SWAP pair %d (fgap=%+.6f -> swap MO %d <-> %d)\n",
                            iteration_, p, gap, i_idx, j_idx);
                    }
                    for (int row = 0; row < nso_; ++row) {
                        double tmp = Ca_->get(row, i_idx);
                        Ca_->set(row, i_idx, Ca_->get(row, j_idx));
                        Ca_->set(row, j_idx, tmp);
                    }
                    double e_tmp = eps[i_idx];
                    eps[i_idx] = eps[j_idx];
                    eps[j_idx] = e_tmp;
                    gvb_diis_had_swaps_ = true;
                    ++gvb_fon_step_swaps_;
                    any_swap = true;
                    if (n_sec > 1) reflect_pair(orbs[0], orbs[1], p);
                }
            }
        }
#ifndef NDEBUG
        if (n_sec == 1)
            assert(capture_all_fons() == dbg_fon_pre &&
                   "site-1 fgap swap reflected a FON in a single-sector run");
#endif

        // SA pool with no active n-geminals (OSS-only, extra-determinant-only,
        // DOSS_SPS-only): no n-FON to seed. A zero-dimension FONMicroSolver would
        // call LAPACK DSYEV with N=0.
        if (P == 0) {
            gvb_fon_initialized_ = true;
            compute_weighting_factors();
            timer_off("REKS: gvb_fon_step");
            return;
        }

        const reks::MicroSolverConfig init_cfg = fon_micro_cfg_;

        const double bnd_skin = fon_micro_boundary_skin_;

        // Per-sector gen-0 seed: multi-start (all-1.0 vs GVB-PP CI), preferring the
        // interior minimum over a boundary basin.
        for (int sec = 0; sec < n_sec; ++sec) {
            const int P_s = static_cast<int>(sa_cassette_.geminals_active(sec, 0).size());
            if (P_s == 0) continue;

            if (reports(4) && n_sec > 1)
                outfile->Printf("  [GVB-FON] iter=%d sector %d init\n", iteration_, sec);

            auto snap = std::make_shared<reks::studio::FONSnapshot>(reel_.fon_state[sec]);
            auto obj  = reks::build_geminal_objective(sec, 0, sa_cassette_, snap, reel_.E_L);
            reks::FONMicroSolver init_solver(P_s);

            auto is_boundary = [&](const std::vector<double>& f) {
                for (int p = 0; p < P_s; ++p)
                    if (std::min(f[p], 2.0 - f[p]) <= bnd_skin) return true;
                return false;
            };

            std::vector<std::vector<double>> starts;
            starts.push_back(std::vector<double>(P_s, 1.0));
            auto fon_ci = reks::fon_init_from_ci(sec, sa_cassette_, reel_.E_L, iteration_,
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
                auto res = init_solver.solve(obj, starts[st], init_cfg);
                double E_try = obj.energy(res.fon);
                bool boundary = is_boundary(res.fon);
                if (reports(4)) {
                    outfile->Printf("  [GVB-FON] iter=%d multi-start %zu: start=(", iteration_, st);
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
                                iteration_, best_idx, best_E,
                                best_interior_fon.empty() ? "no interior found, using best E" : "interior preferred");
            }

            set_fon_vector(sec, 0, best_fon);
            if (reports(4)) {
                outfile->Printf("  [GVB-FON] iter=%d init: FON=(", iteration_);
                for (int p = 0; p < P_s; ++p)
                    outfile->Printf("%.6f%s", best_fon[p], p + 1 < P_s ? "," : "");
                outfile->Printf(")%s\n", any_swap ? " (after fgap inv-swap)" : "");
            }
        }

        gvb_fon_initialized_ = true;
        compute_weighting_factors();
        timer_off("REKS: gvb_fon_step");
        return;
    }

    constexpr double fon_bnd_margin = 1e-8;
    // Per-generation FON floor: REKS_<G>_FON if set (>= 0), else open lower 1e-8.
    auto fon_lower = [&](int gen) {
        return reks_fon_lower_[gen] >= 0.0 ? reks_fon_lower_[gen] : 1e-8;
    };

    // gen-0 crossing record, keyed by the aggregate pair position (union of all
    // sectors' gen-0 pairs); ORs across sectors that share the pair.
    std::vector<int> g0_to_pos(sa_cassette_.n_geminals(), -1);
    for (int p = 0; p < P; ++p) g0_to_pos[gen0_pool[p]] = p;
    std::vector<char> crossed(P, 0);
    std::vector<double> swap_before(P, 0.0), swap_after(P, 0.0);

    // Per-(sector, generation) joint multidim Newton-Raphson, mirroring the TRAH
    // fon_blocks_ walk. The FON-FON Hessian is block-diagonal across (sector,
    // generation), so the loop order is irrelevant; each solve holds the other
    // blocks fixed at their current snapshot value. A config's weight depends only
    // on its own sector's FONs, so each block solves against its sector's snapshot.
    for (int sec = 0; sec < n_sec; ++sec) {
        if (reports(4) && n_sec > 1)
            outfile->Printf("  [GVB-FON] iter=%d sector %d\n", iteration_, sec);
        for (int gen = 0; gen < n_gen; ++gen) {
            const auto& gem = sa_cassette_.geminals_active(sec, gen);
            const int Pg = static_cast<int>(gem.size());
            if (Pg == 0) continue;  // no active geminals: skip the zero-dim DSYEV.

            const double lo = fon_lower(gen);
            std::vector<double> x0 = fon_vector(sec, gen);

            auto snap = std::make_shared<reks::studio::FONSnapshot>(reel_.fon_state[sec]);
            auto obj  = reks::build_geminal_objective(sec, gen, sa_cassette_, snap, reel_.E_L);

            reks::FONMicroSolver solver(Pg);
            reks::MicroSolverConfig cfg = fon_micro_cfg_;
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

            set_fon_vector(sec, gen, res.fon);

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
                                iteration_, chan.c_str(), res.nr_iterations, chan.c_str());
                for (int k = 0; k < Pg; ++k)
                    outfile->Printf("%.8f%s", res.fon[k], k + 1 < Pg ? ", " : "");
                outfile->Printf(") %s-prev=(", chan.c_str());
                for (int k = 0; k < Pg; ++k)
                    outfile->Printf("%.8f%s", x0[k], k + 1 < Pg ? ", " : "");
                outfile->Printf(")\n");
            }
        }
    }

    // Orbital swap when a gen-0 (n) FON crosses 1.0: swap the shared bonding/
    // antibonding columns once per pair and reflect every (sector, generation) FON
    // on that pair.
    for (int p = 0; p < P; ++p) {
        if (!crossed[p]) continue;

        if (reports(4)) {
            outfile->Printf("  [GVB-FON] iter=%d swap pair %d (FON %.6f -> %.6f)\n",
                            iteration_, p, swap_before[p], swap_after[p]);
        }
        const int g0 = gen0_pool[p];
        const auto& orbs = templ[g0].orbitals;
        const int i_idx = active_mo_indices_[orbs[0]];
        const int j_idx = active_mo_indices_[orbs[1]];
        for (int row = 0; row < nso_; ++row) {
            double tmp = Ca_->get(row, i_idx);
            Ca_->set(row, i_idx, Ca_->get(row, j_idx));
            Ca_->set(row, j_idx, tmp);
        }
        gvb_diis_had_swaps_ = true;
        ++gvb_fon_step_swaps_;
        reflect_pair(orbs[0], orbs[1], p);
    }

    compute_weighting_factors();
    timer_off("REKS: gvb_fon_step");
}

void REKS::form_D() {
    timer_on("REKS: form_D");
    ScopedStage _ss_form_D(reports(4) ? &scf_iter_times_ : nullptr, "form_D");
    // C_occ_cache_ stays valid for the rest of this SCF iter: Ca_ does not
    // change between form_D and build_sa_focks.
    C_occ_cache_ = Ca_subset("SO", "OCC");
    if (reports(4) && iteration_ <= 1) {
        outfile->Printf("  [READ_DBG] form_D iter=%d: C_occ cols=%d  nalphapi_[0]=%d  Ncore+n_active=%d\n",
                        iteration_, C_occ_cache_->colspi()[0], nalphapi_[0],
                        Ncore_ + sa_cassette_.n_active_orbitals());
        int nso = Ca_->rowspi()[0];
        int ncol = std::min(nmopi_[0], 8);
        outfile->Printf("  [READ_DBG]   Ca_ col norms at form_D iter=%d:", iteration_);
        double** Cp = Ca_->pointer(0);
        for (int i = 0; i < ncol; ++i) {
            double nn = 0.0;
            for (int mu = 0; mu < nso; ++mu) nn += Cp[mu][i] * Cp[mu][i];
            outfile->Printf(" c%d=%.4f", i, std::sqrt(nn));
        }
        outfile->Printf("\n");
    }
    build_base_densities();
    form_Da_from_core();

    if (reports(5)) {
        outfile->Printf("\n  === Active Orbitals (after form_D) ===\n");

        int active_r = active_mo_indices_[0];
        int active_s = active_mo_indices_[1];

        // Overlap <r|s> = C_r^T * S * C_s
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

        // Norms <r|r>, <s|s>
        double norm_r = 0.0, norm_s = 0.0;
        for (int mu = 0; mu < nso; ++mu) {
            for (int nu = 0; nu < nso; ++nu) {
                norm_r += Cp[mu][active_r] * Sp[mu][nu] * Cp[nu][active_r];
                norm_s += Cp[mu][active_s] * Sp[mu][nu] * Cp[nu][active_s];
            }
        }
        outfile->Printf("  <r|r> = %12.8f, <s|s> = %12.8f (should be 1.0)\n", norm_r, norm_s);
    }
    _ss_form_D.release();
    if (reports(4) && iteration_ > 0 && iteration_ != last_logged_iter_) {
        double iter_total = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - iter_start_time_).count();
        log_iter_timings_(iteration_, iter_total);
        for (const auto& kv : scf_iter_times_) scf_total_times_[kv.first] += kv.second;
        ++n_iters_logged_;
        last_logged_iter_ = iteration_;
    }
    timer_off("REKS: form_D");
}

void REKS::form_G() {
    timer_on("REKS: form_G");
    // One-shot: a stand-on-best reported pair lives exactly one iteration.
    gvb_report_E_override_armed_ = false;
    // SCF-iter boundary: reset per-iter timer.
    if (reports(4) && iteration_ > 0 && iteration_ != last_logged_iter_) {
        if (n_iters_logged_ == 0) {
            scf_start_time_ = std::chrono::steady_clock::now();
        }
        scf_iter_times_.clear();
        iter_start_time_ = std::chrono::steady_clock::now();
    }
    ScopedStage _ss_form_G(reports(4) ? &scf_iter_times_ : nullptr, "form_G");
    if (reports(4) && iteration_ <= 2) {
        outfile->Printf("  [READ_DBG] form_G entry: iter=%d sad_=%d  Ca_subset(OCC).colspi=%d  nalphapi_[0]=%d\n",
                        iteration_, sad_ ? 1 : 0,
                        Ca_subset("SO", "OCC")->colspi()[0], nalphapi_[0]);
    }
    // SAD guess iteration: inline JK build.
    if (sad_ && iteration_ <= 0) {
        auto C_occ_sad = Ca_subset("SO", "OCC");
        jk_->C_left().clear();
        jk_->C_left().push_back(C_occ_sad);
        {
            timer_on("REKS: jk_compute (sad)");
            ScopedStage _ss_jk(reports(4) ? &scf_iter_times_ : nullptr, "jk_compute_sad");
            jk_->compute();
            timer_off("REKS: jk_compute (sad)");
        }
        J_->copy(jk_->J()[0]);
        G_->copy(J_);
        G_->scale(2.0);
        if (functional_->is_x_hybrid()) {
            K_->copy(jk_->K()[0]);
            G_->axpy(-functional_->x_alpha(), K_);
        }
        if (functional_->is_x_lrc()) {
            wK_->copy(jk_->wK()[0]);
            G_->axpy(jk_->get_wcombine() ? -1.0 : -functional_->x_beta(), wK_);
        }
        if (functional_->needs_xc() && rv_potential_) {
            rv_potential_->set_D({Da_});
            {
                timer_on("REKS: xc_compute (sad)");
                ScopedStage _ss_xc(reports(4) ? &scf_iter_times_ : nullptr, "xc_compute_sad");
                rv_potential_->compute_V({Va_});
                timer_off("REKS: xc_compute (sad)");
            }
            G_->add(Va_);
        }
        timer_off("REKS: form_G");
        return;
    }

    // SA microstate Focks in a single JK batch.
    build_sa_focks();

    compute_sa_energies();

    if (reports(4)) {
        print_sa_diagnostics();
    }

    // FON relaxation: TRAH already moved it in combined_step -> refresh weights only;
    // otherwise gvb_fon_step (a Newton solver despite the GVB name) relaxes it.
    const bool trah_owns_fon = use_trah_ && !ctrl_state_.is_active();
    if (trah_owns_fon) {
        compute_weighting_factors();
    } else {
        gvb_fon_step();
    }

    if (use_gvb_diis_ && reports(4)) {
        const int n_sec = sa_cassette_.n_sectors();
        // Gen-0 (n) FON, gradient, and Hessian diagonal per sector; the sector token
        // is empty in a single-sector run.
        for (int sec = 0; sec < n_sec; ++sec) {
            const int np = static_cast<int>(sa_cassette_.geminals_active(sec, 0).size());
            if (np == 0) continue;
            auto fon_vec = fon_vector(sec, 0);
            auto fon_grad = reks::REKSGradientEngine::compute_fon_gradient(
                sec, 0, sa_cassette_, reel_.fon_state[sec], reel_.E_L);
            auto fon_hess = reks::REKSGradientEngine::compute_fon_hessian(
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
                    outfile->Printf("  [GVB-FON] iter=%d %sBOUNDARY pair %d: f=%.8f g=%+.6e h=%+.6e newton_step=%+.6e\n",
                                    iteration_, sec_tok.c_str(), p, fon_vec[p], fon_grad[p], fon_hess[p * np + p], newton_step);
                }
            }
        }
    }

    if (reports(4)) {
        print_iteration_diagnostics();
    }

    timer_off("REKS: form_G");
}

// Convergence quality of the finished GVB-DIIS run: writes wfn variables and one
// trace line.
//
//   crossing  CLEAN   the final iterate cleared both gates with margin
//             GRAZING it landed within a decade of one of them
void REKS::report_convergence_quality() {
    const double gorb_fro = gvb_report_gorb_;
    const double rms = gorb_rms(gorb_fro);
    const int rung = static_cast<int>(ctrl_state_.rung());

    const bool grazing = (gvb_gate_dE_valid_ && std::abs(gvb_gate_dE_) > 0.1 * gvb_e_conv_) ||
                         (rms > 0.1 * gvb_d_conv_);

    set_scalar_variable("REKS QUALITY GORB FRO", gorb_fro);
    set_scalar_variable("REKS QUALITY GORB RMS", rms);
    set_scalar_variable("REKS QUALITY LANDINGS", static_cast<double>(gvb_landing_count_));
    set_scalar_variable("REKS QUALITY VERDICT TRIPS", static_cast<double>(gvb_verdict_trips_));
    set_scalar_variable("REKS QUALITY MAX RUNG", static_cast<double>(rung));
    set_scalar_variable("REKS QUALITY BASIN FIRES",
                        static_cast<double>(outcome_guard_.b_fire_count()));
    set_scalar_variable("REKS QUALITY GRAZING", grazing ? 1.0 : 0.0);
    set_scalar_variable("REKS QUALITY NK TRIGGER STREAK",
                        static_cast<double>(nk_trigger_streak_max_));

    if (reports(4)) {
        static const char* kRung[] = {"ARMED", "RETAIN", "FALLBACK"};
        outfile->Printf(
            "  [QUALITY] iter=%d gorb_fro=%.3e gorb_rms=%.3e landings=%d trips=%d rung=%s "
            "basin_fires=%d crossing=%s nk_streak_max=%d nk_grad_evals=%d\n",
            iteration_, gorb_fro, rms, gvb_landing_count_, gvb_verdict_trips_, kRung[rung],
            outcome_guard_.b_fire_count(), grazing ? "GRAZING" : "CLEAN", nk_trigger_streak_max_,
            nk_grad_evals_);
    }
}

// Rejected step = null step (Conn, A. R.; Gould, N. I. M.; Toint, P. L. Trust-Region
// Methods; SIAM: Philadelphia, 2000).
void REKS::execute_landed_rewind(GvbReportWitness witness, bool retain_history) {
    Ca_->copy(gvb_best_Ca_);
    Cb_->copy(gvb_best_Ca_);
    restore_all_fons(gvb_best_fons_);
    reks::StepContext ctx;
    ctx.Ca = &Ca_;
    ctx.Ca_ref = &Ca_ref_;
    ctx.cum_kappa = &orbital_diis_cum_kappa_;
    ctx.Ncore = Ncore_;
    // Frame guard: a frame-relative (ORB) history survives only in its own anchor frame,
    // compared by frame identity (epoch) and evaluated before restore().
    bool retain = retain_history;
    if (retain && gvb_best_engine_memento_.Ca_ref &&
        gvb_best_engine_memento_.epoch != diis_controller_.active()->epoch())
        retain = false;
    diis_controller_.active()->restore(gvb_best_engine_memento_, ctx);
    if (retain)
        diis_controller_.active()->on_rewind_retain(ctx);
    else
        diis_controller_.active()->on_rewind(ctx);
    ctrl_state_.on_landing_or_restart(iteration_);
    ++gvb_landing_count_;
    gvb_null_step_ = true;
    gvb_report_witness_ = witness;
}

void REKS::gvb_report_pair(GvbReportWitness witness) {
    switch (witness) {
        case GvbReportWitness::CURRENT:
        case GvbReportWitness::REJECTED:
            gvb_report_gorb_ = gvb_diis_gorb_norm_;
            gvb_report_E_override_armed_ = false;
            break;
        case GvbReportWitness::BEST:
            gvb_report_gorb_ = gvb_best_gorb_;
            gvb_report_E_override_armed_ = true;
            gvb_report_E_override_ = gvb_best_E_;
            break;
    }
}

// Builds the REKS coupling Fock (reel_.F_reks_MO/F_reks); once GVB-DIIS is active
// (ctrl_state_) also runs the orbital-update step, else leaves Fa_/Fb_ = F_reks.
void REKS::form_F() {
    timer_on("REKS: form_F");
    ScopedStage _ss_form_F(reports(4) ? &scf_iter_times_ : nullptr, "form_F");
    // Snapshot INCOMING Fa_ on iter<=0 before form_F overwrites it.
    if (reports(4) && iteration_ <= 0 && Fa_) {
        int nmo = Ca_->colspi()[0];
        int n_act = static_cast<int>(active_mo_indices_.size());
        if (n_act > 0 && nmo >= n_act) {
            auto F_MO_in = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
            double** Fmo = F_MO_in->pointer(0);
            double max_off_aa = 0.0;
            for (int i = 0; i < n_act; ++i) {
                for (int j = i + 1; j < n_act; ++j) {
                    max_off_aa = std::max(max_off_aa,
                                          std::abs(Fmo[active_mo_indices_[i]][active_mo_indices_[j]]));
                }
            }
            double max_av_in = 0.0;
            int last_active = active_mo_indices_.back();
            for (int ai : active_mo_indices_) {
                for (int v = last_active + 1; v < nmo; ++v) {
                    max_av_in = std::max(max_av_in, std::abs(Fmo[ai][v]));
                }
            }
            outfile->Printf("  [READ_DBG] form_F iter=%d INCOMING Fa_ in MO: "
                            "max_off_aa=%.3e max_av=%.3e (sad=%d guess_Ca=%s)\n",
                            iteration_, max_off_aa, max_av_in,
                            sad_ ? 1 : 0, guess_Ca_ ? "SET" : "null");
        }
    }

    if (reports(4) && iteration_ <= 2) {
        int n_act = static_cast<int>(active_mo_indices_.size());
        outfile->Printf("  [READ_DBG] form_F entry: iter=%d sad_=%d guess_Ca_=%s\n",
                        iteration_, sad_ ? 1 : 0, guess_Ca_ ? "SET" : "nullptr");
        outfile->Printf("  [READ_DBG]   nalphapi_[0]=%d nbetapi_[0]=%d\n", nalphapi_[0], nbetapi_[0]);
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
        int nso = Ca_->rowspi()[0];
        int nmo = Ca_->colspi()[0];
        double** Cp = Ca_->pointer(0);
        auto col_norm = [&](int a) {
            double nn = 0.0;
            for (int mu = 0; mu < nso; ++mu) nn += Cp[mu][a] * Cp[mu][a];
            return std::sqrt(nn);
        };
        outfile->Printf("  [READ_DBG] iter=%d Ca_ col norms:", iteration_);
        for (int c = 0; c < Ncore_; ++c)
            outfile->Printf(" core[%d]=%.4f", c, col_norm(c));
        for (int k = 0; k < n_act; ++k) {
            int a = active_mo_indices_[k];
            outfile->Printf(" act[%d]=%.4f", a, col_norm(a));
        }
        int first_virt = n_act > 0 ? active_mo_indices_.back() + 1 : Ncore_;
        for (int v = first_virt; v < nmo; ++v)
            outfile->Printf(" virt[%d]=%.4f", v, col_norm(v));
        outfile->Printf("\n");
    }
    if (sad_ && iteration_ <= 0) {
        if (reports(4)) {
            outfile->Printf("\n  form_F: SAD iteration\n");
        }
        Fa_->copy(H_);
        Fa_->add(G_);
        timer_off("REKS: form_F");
        return;
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

        auto g_fon = reks::REKSGradientEngine::compute_fon_gradient(
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

    // PGC: KKT boundary stationarity check. all_boundary holds when every pool
    // FON is pinned at a bound with the gradient pointing into the wall:
    //   FON >= 2 - eps and g < 0 (upper), or FON <= eps and g > 0 (lower).
    // Suppressed once GVB-DIIS is active: its form_FDSmSDF branch owns the reported rms.
    if (use_trah_ && !ctrl_state_.is_active() && iteration_ > 0) {
        int n_fon_pgc = static_cast<int>(sa_cassette_.geminals_active_gen(0).size());
        // Pool-sized engine gradient: g_fon_pgc[k] is geminal geminals_active_gen(0)[k].
        auto g_fon_pgc = reks::REKSGradientEngine::compute_fon_gradient(
            0, 0, sa_cassette_, sa_fons(), reel_.E_L);
        bool all_boundary = (n_fon_pgc > 0);
        for (int k = 0; k < n_fon_pgc && all_boundary; ++k) {
            double fon = sa_fons().layers[0][sa_cassette_.geminals_active_gen(0)[k]].p;
            bool at_upper = (fon >= 2.0 - 1e-8) && (g_fon_pgc[k] < 0);
            bool at_lower = (fon <= 1e-8) && (g_fon_pgc[k] > 0);
            if (!at_upper && !at_lower) all_boundary = false;
        }
        double g_orb_sq_pgc = 0.0;
        {
            auto F_gen_pgc = build_generalized_fock();
            double** Fg = F_gen_pgc->pointer(0);
            int n_act = static_cast<int>(active_mo_indices_.size());
            for (int i = 0; i < n_act; ++i)
                for (int j = i + 1; j < n_act; ++j) {
                    int ai = active_mo_indices_[i], aj = active_mo_indices_[j];
                    double g_ij = -2.0 * (Fg[ai][aj] - Fg[aj][ai]);
                    g_orb_sq_pgc += g_ij * g_ij;
                }
        }
        pgc_gorb_norm_ = std::sqrt(g_orb_sq_pgc);
        // Suppress PGC boundary convergence on un-validated READ guess:
        // converging here would short-circuit the boundary-trap escape.
        if (all_boundary && guess_Ca_ && has_read_fon_guess() && !trah_state_.initialized) {
            pgc_boundary_stationary_ = false;
            if (reports(4)) {
                outfile->Printf("  [AUTOPILOT] iter=%d PGC: boundary detected but SUPPRESSED "
                                "(READ guess, TRAH not yet initialized, ||g_orb||=%.6e)\n",
                                iteration_, pgc_gorb_norm_);
            }
        } else {
            pgc_boundary_stationary_ = all_boundary;
            if (all_boundary && reports(4)) {
                outfile->Printf("  [AUTOPILOT] iter=%d PGC: boundary stationary, ||g_orb||=%.6e\n", iteration_,
                                pgc_gorb_norm_);
            }
        }
    } else {
        pgc_boundary_stationary_ = false;
    }

    assemble_F_reks_MO();

    // Pure-GGA virt-virt: replace F_reks block with the closed-shell RKS Fock to avoid
    // pathological virtual energies from negative microstate weights. Transform only the
    // virtual columns (O(nso^2 n_virt)), not the full MO Fock.
    if (functional_ && functional_->needs_xc() && functional_->is_gga() && functional_->x_alpha() < 0.01) {
        const int first_virt = active_mo_indices_.back() + 1;
        const int n_virt = N - first_virt;
        if (n_virt > 0) {
            const int nso = Ca_->rowspi()[0];
            // G_ is the closed-shell Fock contribution (2J - x_alpha K - .. wK + Vxc).
            auto F_rhf_AO = std::make_shared<Matrix>("F_rhf_AO", nsopi_, nsopi_);
            F_rhf_AO->copy(H_);
            F_rhf_AO->add(G_);

            // Cv is the trailing column block of Ca_, used in place with lda = nmo; the
            // result (Cv^T F_rhf_AO) Cv is written into F_reks_MO with ldc = nmo.
            const int nmo_c = Ca_->colspi()[0];
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

    // SC = S*Ca, the AO<-MO map: F_AO = SC F_MO SC^T.
    auto SC = linalg::doublet(S_, Ca_, false, false);

    // Activation gate: ORBITAL activates immediately; CFM waits for gorb < gorb_thresh
    // unless TRAH is off.
    // True once the CFM activation gate already built this iteration's generalized Fock
    // into F_gen_ (build_generalized_fock's single-writer workspace).
    bool f_gen_current = false;
    if (use_gvb_diis_ && !ctrl_state_.is_active() && iteration_ >= gvb_diis_start_) {
        bool gate_passes = false;
        double gorb_check = 0.0;
        if (diis_formulation_ == DIISFormulation::ORBITAL) {
            gate_passes = true;
        } else {
            F_gen_ = build_generalized_fock();
            f_gen_current = true;
            gorb_check = compute_full_gradient_norm(F_gen_);
            const double gorb_thresh = options_.get_double("REKS_DIIS_CFM_ACTIVATION_GORB");
            gate_passes = (gorb_check < gorb_thresh || !use_trah_);
            if (!gate_passes && reports(4)) {
                outfile->Printf("  [GVB-DIIS] iter=%d WAITING: gorb=%.3e > thresh=%.3e\n",
                                iteration_, gorb_check, gorb_thresh);
            }
        }

        if (gate_passes) {
            // Activation resets both axes and the landing accounting in one place.
            ctrl_state_.activate();
            cfm_adapter_.reset();
            gvb_best_Ca_.reset();
            gvb_best_fons_.clear();
            gvb_best_engine_memento_ = reks::EngineMemento{};
            gvb_best_gorb_ = -1.0;
            gvb_best_E_ = 0.0;
            gvb_prev_did_extrap_ = false;
            gvb_verdict_trips_ = 0;
            gvb_landing_count_ = 0;
            gvb_null_step_ = false;
            gvb_report_witness_ = GvbReportWitness::CURRENT;
            gvb_report_E_override_armed_ = false;
            diis_controller_.arm();
            outcome_guard_.arm();

            if (diis_formulation_ == DIISFormulation::ORBITAL) {
                int nmo_act = Ca_->colspi()[0];
                int n_act = static_cast<int>(active_mo_indices_.size());
                int first_virt = (n_act > 0) ? active_mo_indices_.back() + 1 : Ncore_;
                int n_virt = nmo_act - first_virt;
                orb_adapter_.reset();
                orbital_diis_cum_kappa_.assign(
                    reks::OrbAdapter::n_nonred(Ncore_, n_act, n_virt), 0.0);
                // Mints the activation anchor through the adapter: same epoch owner
                // and projection as every other anchor event.
                {
                    reks::StepContext anchor_ctx;
                    anchor_ctx.Ca = &Ca_;
                    anchor_ctx.Ca_ref = &Ca_ref_;
                    anchor_ctx.cum_kappa = &orbital_diis_cum_kappa_;
                    anchor_ctx.Ncore = Ncore_;
                    orb_adapter_.anchor_here(anchor_ctx);
                }
                orbital_runtime_active_ = true;
                orbital_cap_fired_this_iter_ = false;
                if (reports(4)) {
                    outfile->Printf(
                        "  [GVB-DIIS] iter=%d SWITCH: %s -> ORBITAL-DIIS (n_kappa=%d)\n",
                        iteration_, use_trah_ ? "TRAH" : "burn-in",
                        static_cast<int>(orbital_diis_cum_kappa_.size()));
                }
            } else if (reports(4)) {
                outfile->Printf("  [GVB-DIIS] iter=%d SWITCH: %s -> GVB-DIIS (gorb=%.3e)\n",
                                iteration_, use_trah_ ? "TRAH" : "burn-in", gorb_check);
            }
        }
    }

    if (ctrl_state_.is_active()) {
        timer_on("REKS: gvb_diis_step");
        // Orbital swaps change MO ordering, invalidating the stored error vectors and the
        // cumulative kappa.
        const bool swap_detected = gvb_diis_had_swaps_;
        gvb_diis_had_swaps_ = false;
        if (swap_detected && reports(4))
            outfile->Printf("  [GVB-DIIS] iter=%d RESET: orbital swap detected\n", iteration_);

        diis_controller_.select(orbital_runtime_active_);

        if (!f_gen_current) F_gen_ = build_generalized_fock();
        int nmo = Ca_->colspi()[0];

        // Orbital gradient at the current point (= the result of the previous step).
        reks::GradientBlocks gblocks;
        gvb_diis_gorb_norm_ = compute_full_gradient_norm(F_gen_, &gblocks);
        gvb_report_witness_ = GvbReportWitness::CURRENT;

        // Stall signature: nk_signature holds when gvb_step_norm_ < kStepCollapse *
        // gvb_prev_step_norm_ and gvb_diis_gorb_norm_ > kStepCollapse * gvb_prev_gorb_
        // while gorb_rms >= gvb_d_conv_. Sampled every active iteration -- a rewind
        // moves the orbitals back a full step, so it breaks the streak on its own.
        constexpr double kStepCollapse = 0.6;
        gvb_step_norm_ = gvb_step_norm(gvb_prev_Ca_, SC);
        if (!gvb_prev_Ca_ || gvb_prev_Ca_->rowspi() != Ca_->rowspi() ||
            gvb_prev_Ca_->colspi() != Ca_->colspi())
            gvb_prev_Ca_ = Ca_->clone();
        else
            gvb_prev_Ca_->copy(Ca_);

        const bool nk_sample = (gvb_step_norm_ >= 0.0) && (gvb_prev_step_norm_ > 0.0) &&
                               (gvb_prev_gorb_ > 0.0);
        const bool nk_signature =
            nk_sample && gorb_rms(gvb_diis_gorb_norm_) >= gvb_d_conv_ &&
            gvb_step_norm_ < kStepCollapse * gvb_prev_step_norm_ &&
            gvb_diis_gorb_norm_ > kStepCollapse * gvb_prev_gorb_;
        nk_trigger_streak_ = nk_signature ? nk_trigger_streak_ + 1 : 0;
        nk_trigger_streak_max_ = std::max(nk_trigger_streak_max_, nk_trigger_streak_);
        if (gvb_step_norm_ >= 0.0) {
            gvb_prev_step_norm_ = gvb_step_norm_;
            gvb_prev_gorb_ = gvb_diis_gorb_norm_;
        }

        // Frame handles only: the verdict path passes ctx to restart(). The build inputs are
        // filled in the solve branch, which a null step skips.
        reks::StepContext ctx;
        ctx.Ca = &Ca_;
        ctx.Ca_ref = &Ca_ref_;
        ctx.cum_kappa = &orbital_diis_cum_kappa_;
        ctx.orb_base_map_damp = diis_orb_base_map_damp_;

        // e_sa_cur = SA energy of the current frame;
        // fon_disp = max_{b,i} |fon_cur[b][i] - fon_best[b][i]|, 0 until a best exists.
        const double e_sa_cur = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);
        double fon_disp = 0.0;
        if (gvb_best_gorb_ >= 0.0) {
            const auto fons_cur = capture_all_fons();
            for (size_t b = 0; b < fons_cur.size() && b < gvb_best_fons_.size(); ++b)
                for (size_t i = 0; i < fons_cur[b].size() && i < gvb_best_fons_[b].size(); ++i)
                    fon_disp =
                        std::max(fon_disp, std::abs(fons_cur[b][i] - gvb_best_fons_[b][i]));
        }
        outcome_guard_.observe(e_sa_cur, fon_disp);

        // Convergence refusal: the D-gate would pass at a localization signature (a FON-silent
        // lower well against a stationary-grade best) -- land on best (null step) and restart the
        // accelerator instead. Supersedes the iterate, so the monitor review below is skipped.
        bool s_adjudicated = false;
        if (gorb_rms(gvb_diis_gorb_norm_) < gvb_d_conv_ && outcome_guard_.localization_signature()) {
            const double gorb_refused = gvb_diis_gorb_norm_;
            execute_landed_rewind(GvbReportWitness::BEST);
            // Past the restart budget the refusal keeps the rewind but drops the accelerator
            // restart and the budget count.
            if (!outcome_guard_.budget_exhausted()) {
                diis_controller_.active()->restart(ctx);
                outcome_guard_.on_monitor_restart();
            }
            s_adjudicated = true;
            if (reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d CONVERGENCE REFUSED: orbital symmetry breaking detected "
                    "(gorb %.3e, E_sa %.10f below stationary best %.10f, fon_disp %.3e)\n",
                    iteration_, gorb_refused, e_sa_cur, gvb_best_E_, fon_disp);
        }

        // Basin transition: the FON branch moved off the tracked one -- invalidate best (capture
        // re-seeds on the current iterate this iteration) and re-arm detection. Structural
        // directive: the base map restarts from the current frame, with no rewind.
        bool basin_fire = false;
        if (!s_adjudicated && outcome_guard_.basin_transition_fire()) {
            basin_fire = true;
            gvb_best_gorb_ = -1.0;
            gvb_best_E_ = 0.0;
            gvb_best_Ca_.reset();
            gvb_best_fons_.clear();
            gvb_best_engine_memento_ = reks::EngineMemento{};
            outcome_guard_.on_best_invalidated();
            diis_controller_.arm();
            ctrl_state_.set_rung(QuiesceLadder::ARMED);
            if (reports(4))
                outfile->Printf(
                    "  [DIIS-CTRL] iter=%d BASIN TRANSITION #%d: FON branch moved (fon_disp %.3e), "
                    "re-baseline from current frame\n",
                    iteration_, outcome_guard_.b_fire_count(), fon_disp);
        }

        // Trajectory verdict on the step just applied: VerdictMonitor rewinds to the best iterate
        // on a ||g_orb|| increase beyond rho and escalates a rewind streak to a restart;
        // CycleMonitor restarts the base map on a no-progress stall or an ORB cap-streak.
        // Pre-step directives (basin transition, orbital swap, IPR basin guard) are evaluated
        // before the review and compose rewind-before-restart.
        if (!s_adjudicated) {
            reks::TrajectoryView traj;
            traj.iteration = iteration_;
            traj.gorb = gvb_diis_gorb_norm_;
            traj.best_gorb = gvb_best_gorb_;
            traj.prev_did_extrap = gvb_prev_did_extrap_;
            traj.best_available = (gvb_best_Ca_ != nullptr);
            traj.cap_fired = orbital_runtime_active_ && orbital_cap_fired_this_iter_;
            traj.last_rollback_iter = ctrl_state_.last_rollback_iter();
            traj.energy_descent = outcome_guard_.energy_descent();

            reks::MonitorAction pre = diis_controller_.pre_step_directive(
                swap_detected, ipr_basin_guard_prestep(), basin_fire);
            const reks::ResolvedAction act = diis_controller_.review(traj);

            // Quiescence ladder: past the restart budget, monitor verdicts escalate
            // deterministically. RETAIN: rewind to best keeping the DIIS subspace minus its
            // newest vector (damped base map) -- the surviving history breaks the rewind
            // replay. FALLBACK: once a landing recurs on an unchanged best in RETAIN, DIIS
            // keeps running with monitor verdicts discarded until convergence or MAXITER.
            // Structural directives (basin transition, swap, IPR) stay live throughout.
            reks::QuiesceInputs qin;
            qin.rung = ctrl_state_.rung();
            qin.quiesced = outcome_guard_.budget_exhausted();
            qin.landing_repeat = outcome_guard_.landing_repeat(gvb_best_E_, gvb_best_gorb_);
            qin.best_below_current = (gvb_best_Ca_ != nullptr) && gvb_best_gorb_ >= 0.0 &&
                                     gvb_best_gorb_ < gvb_diis_gorb_norm_;

            const reks::ResolvedAction dec = diis_controller_.apply_quiescence(act, qin);

            const double gorb_bad = gvb_diis_gorb_norm_;

            if (dec.set_rung) ctrl_state_.set_rung(dec.rung_target);
            if (dec.rewind)
                execute_landed_rewind(
                    dec.landed_on_best ? GvbReportWitness::BEST : GvbReportWitness::REJECTED,
                    dec.retain_history);
            if (dec.record_landing) outcome_guard_.record_landing(gvb_best_E_, gvb_best_gorb_);
            if (dec.count_toward_budget) outcome_guard_.on_monitor_restart();
            if (dec.verdict_trip) ++gvb_verdict_trips_;

            if (dec.set_rung && dec.rung_target == QuiesceLadder::FALLBACK) {
                if (reports(4))
                    outfile->Printf(
                        "  [DIIS-CTRL] iter=%d QUIESCE FALLBACK: retain restart landed on an "
                        "unchanged best; monitor verdicts discarded until convergence or "
                        "MAXITER (%s)\n",
                        iteration_, dec.reason);
            } else if (dec.set_rung && dec.rung_target == QuiesceLadder::RETAIN) {
                if (reports(4))
                    outfile->Printf(
                        "  [DIIS-CTRL] iter=%d QUIESCE RETAIN: restart budget exhausted (%d); "
                        "rewind to best keeping DIIS history (%s)\n",
                        iteration_, outcome_guard_.restarts_since_improvement(), dec.reason);
            } else if (dec.rewind) {
                if (dec.repeat_landing && reports(4))
                    outfile->Printf(
                        "  [DIIS-CTRL] iter=%d REWIND_REPEAT: landing on unchanged best "
                        "(gorb %.3e) counted toward restart budget (%d)\n",
                        iteration_, gvb_best_gorb_, outcome_guard_.restarts_since_improvement());
                if (reports(4))
                    outfile->Printf(
                        "  [DIIS-CTRL] iter=%d REWIND_TO_BEST (%s): gorb %.3e "
                        "-> best %.3e (null step)\n",
                        iteration_, dec.reason, gorb_bad, gvb_best_gorb_);
            }

            if (dec.restart) {
                diis_controller_.active()->restart(ctx);
                ctrl_state_.on_landing_or_restart(iteration_);
                if (reports(4))
                    outfile->Printf("  [DIIS-CTRL] iter=%d RESTART_BASE_MAP (%s)\n", iteration_,
                                    dec.reason);
            }

            // Structural directives run regardless of the monitor verdict: idempotent history
            // flush + re-anchor from the current frame.
            if (pre.verdict == reks::MonitorVerdict::RESTART_BASE_MAP) {
                diis_controller_.active()->restart(ctx);
                ctrl_state_.on_landing_or_restart(iteration_);
                ctrl_state_.set_gate_streak(0);
                if (reports(4))
                    outfile->Printf("  [DIIS-CTRL] iter=%d RESTART_BASE_MAP (%s)\n", iteration_,
                                    pre.reason);
            }
        }

        // Terminal MAXITER landing: exhaustion without convergence ends on the best iterate in a
        // consistent frame.
        if (iteration_ >= gvb_scf_maxiter_ && gorb_rms(gvb_diis_gorb_norm_) >= gvb_d_conv_ &&
            gvb_best_Ca_ && gvb_best_gorb_ >= 0.0) {
            if (gvb_null_step_) {
                gvb_report_witness_ = GvbReportWitness::BEST;
            } else if (gvb_best_gorb_ < gvb_diis_gorb_norm_) {
                const double gorb_last = gvb_diis_gorb_norm_;
                execute_landed_rewind(GvbReportWitness::BEST);
                if (reports(4))
                    outfile->Printf("  [DIIS-CTRL] iter=%d MAXITER: landing on best (gorb %.3e -> %.3e)\n",
                                    iteration_, gorb_last, gvb_best_gorb_);
            }
        }

        // Best-iterate capture (lowest gorb seen this SCF): the rewind target. Skipped on a null
        // step -- the reported gorb belongs to the rejected iterate, not the restored frame.
        if (!gvb_null_step_ && !outcome_guard_.refuse_capture() &&
            (gvb_best_gorb_ < 0.0 || gvb_diis_gorb_norm_ < gvb_best_gorb_)) {
            gvb_best_gorb_ = gvb_diis_gorb_norm_;
            gvb_best_E_ = e_sa_cur;
            gvb_best_Ca_ = Ca_->clone();
            gvb_best_fons_ = capture_all_fons();
            reks::StepContext snap_ctx;
            snap_ctx.Ca_ref = &Ca_ref_;
            snap_ctx.cum_kappa = &orbital_diis_cum_kappa_;
            gvb_best_engine_memento_ = diis_controller_.active()->snapshot(snap_ctx);
            outcome_guard_.on_capture(gvb_best_E_, gvb_best_gorb_);
        }

        if (reports(4) && lambda_ipr_ > 0.0) {
            double** Fg_ptr = F_gen_->pointer(0);
            int n_act_loc = static_cast<int>(active_mo_indices_.size());
            double g_orb_max = 0.0;
            for (int i = 0; i < n_act_loc; ++i)
                for (int j = i + 1; j < n_act_loc; ++j) {
                    int ai = active_mo_indices_[i];
                    int aj = active_mo_indices_[j];
                    double g = -2.0 * (Fg_ptr[ai][aj] - Fg_ptr[aj][ai]);
                    g_orb_max = std::max(g_orb_max, std::abs(g));
                }
            print_ipr_iteration_diagnostic(iteration_, g_orb_max);
        }

        if (!gvb_null_step_) {
            // Collapsed solve through the active adapter: build payload + error, store,
            // extrapolate (else base map), apply -- CFM diagonalizes F^opt, ORB accumulates
            // kappa for the deferred expm rotation.
            std::vector<double> gamma_aa, gamma_ca;
            compute_gamma_vectors(gamma_aa, gamma_ca);
            std::vector<double> ipr_diag_aa;
            ipr_hessian_diag_aa(ipr_diag_aa);

            // Two independent stall signatures open the episode: nk_trigger_streak_ counts
            // consecutive iterations whose step collapses faster than the gradient (flat
            // stall); gvb_landing_count_ counts rollbacks (a limit cycle). A rollback
            // resets nk_trigger_streak_.
            const bool nk_stall = nk_trigger_streak_ >= nk_trigger_;
            const bool nk_cycle = gvb_landing_count_ >= nk_min_landings_;

            // miniTRAH episode, at most one per SCF, on the same side of ctx.F_gen as the
            // probe. A landing sets gvb_null_step_ (the accepted point persists) and
            // re-anchors the engine, whose history the episode invalidated.
            if (nk_enabled_ && !nk_episode_done_ && (nk_stall || nk_cycle)) {
                const NKCurvature curv{gamma_aa, gamma_ca, ipr_diag_aa};
                std::vector<double> denom_live;
                orbital_diis_.compute_newton_rotation(
                    F_gen_, reel_.F_MO_diag_a_L, reel_.F_MO_diag_b_L, reel_.C_L, sa_cassette_,
                    active_mo_indices_, Ncore_, nmo, gamma_aa, gamma_ca, iteration_,
                    ipr_diag_aa, nullptr, &denom_live, /*narrate=*/false);

                // The probe relaxes the FON at every displaced geometry, so the map the finite
                // difference samples is g_orb(kappa, n*(kappa)) and its Jacobian is the Schur
                // complement H_oo - H_of H_ff^-1 H_fo. An active FON bound pins n*, killing
                // dn*/dkappa and leaving H_oo alone; with none active the correction is live.
                // A weakly active bound is the worst of both: it breaks strict complementary
                // slackness, and with it the linearity the Krylov method assumes.
                double kkt_min = std::numeric_limits<double>::max();
                for (double m : nk_fon_kkt_multipliers())
                    if (m != 0.0) kkt_min = std::min(kkt_min, std::abs(m));
                const bool has_active_bound = kkt_min != std::numeric_limits<double>::max();
                const bool kkt_ok = nk_require_active_fon_
                                        ? (has_active_bound && kkt_min >= 1e-3)
                                        : (!has_active_bound || kkt_min >= 1e-3);

                if (reports(4))
                    outfile->Printf(
                        "  [GVB-NK] iter=%d TRIGGER(%s) streak=%d landings=%d step_norm=%.3e "
                        "gorb_rms=%.3e kkt_min=%.3e %s\n",
                        iteration_, nk_stall ? (nk_cycle ? "stall+cycle" : "stall") : "cycle",
                        nk_trigger_streak_, gvb_landing_count_, gvb_step_norm_,
                        gorb_rms(gvb_diis_gorb_norm_),
                        (kkt_min == std::numeric_limits<double>::max()) ? 0.0 : kkt_min,
                        kkt_ok ? "ENTER" : "BLOCKED (no strictly active FON bound)");

                // A blocked entry costs no Fock build and is not an episode; it only spends
                // the streak.
                if (!kkt_ok) nk_trigger_streak_ = 0;
                bool landed = false;
                if (kkt_ok) {
                    nk_episode_done_ = true;
                    landed = nk_run_episode(curv, denom_live);
                }

                if (landed) {
                    gvb_null_step_ = true;
                    ctx.Ncore = Ncore_;
                    diis_controller_.active()->restart(ctx);
                    ctrl_state_.on_landing_or_restart(iteration_);
                    gvb_prev_did_extrap_ = false;
                    orbital_cap_fired_this_iter_ = false;
                    nk_trigger_streak_ = 0;
                    gvb_report_pair(gvb_report_witness_);
                    timer_off("REKS: gvb_diis_step");
                    timer_off("REKS: form_F");
                    return;
                }
            }

            std::vector<double> fon_dbg;
            if (reports(4) && !orbital_runtime_active_) fon_dbg = n_fon_vector();

            ctx.F_gen             = F_gen_;
            ctx.F_MO_diag_a_L     = &reel_.F_MO_diag_a_L;
            ctx.F_MO_diag_b_L     = &reel_.F_MO_diag_b_L;
            ctx.C_L               = &reel_.C_L;
            ctx.sa_cassette       = &sa_cassette_;
            ctx.active_mo_indices = &active_mo_indices_;
            ctx.Ncore             = Ncore_;
            ctx.nmo               = nmo;
            ctx.gamma_aa          = &gamma_aa;
            ctx.gamma_ca          = &gamma_ca;
            ctx.ipr_diag_aa       = &ipr_diag_aa;
            ctx.SC                = SC;
            ctx.F_reks_MO         = reel_.F_reks_MO;
            ctx.reel_F_reks       = reel_.F_reks;
            ctx.Fa                = &Fa_;
            ctx.Fb                = &Fb_;
            ctx.iteration         = iteration_;
            ctx.gorb              = gvb_diis_gorb_norm_;
            ctx.n_fon_dbg         = fon_dbg.empty() ? nullptr : &fon_dbg;

            reks::StepResult sr = diis_controller_.step(ctx);
            gvb_prev_did_extrap_         = sr.extrap_ok;
            orbital_cap_fired_this_iter_ = sr.cap_fired;
        } else {
            // Null step: no solve.
            gvb_prev_did_extrap_ = false;
            orbital_cap_fired_this_iter_ = false;
        }
        gvb_report_pair(gvb_report_witness_);
        timer_off("REKS: gvb_diis_step");
    } else {
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

        const double pr_min = active_pr_min();
        if (pr_min >= 0.0)
            outfile->Printf("  [PR] iter=%d PR_min=%.6f\n", iteration_, pr_min);
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

// [F,D]_S = F*D*S - S*D*F via HF::form_FDSmSDF, rescaled here to gorb-consistent
// RMS units for GVB-DIIS/PGC/TRAH convergence reporting.
SharedMatrix REKS::form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso) {
    auto result = HF::form_FDSmSDF(Fso, Dso);

    // GVB-DIIS: rescale the commutator so rms(result) == ||g_orb||_F / nmo
    // (stock D_CONVERGENCE units).
    if (ctrl_state_.is_active() && gvb_diis_gorb_norm_ > 0.0) {
        set_scalar_variable("REKS GORB FRO", gvb_diis_gorb_norm_);
        const double true_rms = gorb_rms(gvb_diis_gorb_norm_);
        // Streak advances once per iteration_; repeat calls this iteration_ are no-ops.
        if (gvb_gate_streak_iter_ != iteration_) {
            gvb_gate_streak_iter_ = iteration_;
            const bool joint = gvb_gate_dE_valid_ && std::abs(gvb_gate_dE_) < gvb_e_conv_ &&
                               true_rms < gvb_d_conv_;
            ctrl_state_.set_gate_streak(joint ? ctrl_state_.gate_streak() + 1 : 0);
        }
        set_scalar_variable("REKS GATE RMS", true_rms);
        double current_rms = result->rms();
        if (current_rms > 1e-15) {
            result->scale(true_rms / current_rms);
        }
    }

    // KKT boundary: replace Dnorm by ||g_orb|| (FON gradient is irrelevant
    // when FON is pinned). Safe -- DIIS is disabled for REKS.
    if (pgc_boundary_stationary_ && use_trah_) {
        double current_rms = result->rms();
        if (current_rms > 1e-15 && pgc_gorb_norm_ < current_rms) {
            result->scale(pgc_gorb_norm_ / current_rms);
            if (reports(4)) {
                outfile->Printf("  [AUTOPILOT] PGC: Dnorm %.2e -> ||g_orb|| %.2e (KKT boundary)\n", current_rms,
                                pgc_gorb_norm_);
            }
        }
    }

    return result;
}

// Build Ca_ from Fa_: iter<=0 dispatch, optional TRAH/ORB-DIIS, F_MO diag with staircase shift, SwapOrbs.
void REKS::form_C(double shift) {
    timer_on("REKS: form_C");
    ScopedStage _ss_form_C(reports(4) ? &scf_iter_times_ : nullptr, "form_C");
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
    if (ctrl_state_.is_active() && gvb_null_step_) {
        gvb_null_step_ = false;
        if (reports(4))
            outfile->Printf("  [GVB-DIIS] iter=%d form_C: null step, orbitals held\n", iteration_);
        timer_off("REKS: form_C");
        return;
    }

    // ORBITAL-DIIS path (Ionova-Carter): cumulative kappa rotates Ca_, bypassing Fock-diag.
    if (ctrl_state_.is_active() && orbital_runtime_active_) {
        const int N_orb = Ca_->colspi()[0];
        auto Ca_new = orb_adapter_.apply_rotation(
            Ca_ref_, orbital_diis_cum_kappa_, active_mo_indices_, Ncore_, N_orb);
        Ca_->copy(Ca_new);

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
            for (double v : orbital_diis_cum_kappa_) max_kappa = std::max(max_kappa, std::abs(v));
            outfile->Printf("  [ORB-DIIS] iter=%d form_C: applied expm(K), max|kappa_cum|=%.4f\n",
                            iteration_, max_kappa);
        }
        timer_off("REKS: form_C");
        return;
    }

    // TRAH joint kappa+FON step before Fock diag, so F_MO sees rotated Ca_.
    if (trah_live() && iteration_ > 0) {
        combined_step();
    }

    int N = nmopi_[0];

    // F_MO = Ca^T Fa Ca: diagonalize in current MO basis (numerical stability vs raw Fa).
    auto F_MO_diis = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);

    // Freeze active-active off-diag so Fock diag does not undo TRAH rotations.
    if (use_trah_ && iteration_ > 0 && !trah_kappa_disabled_ && !ctrl_state_.is_active()) {
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

    // Staircase level-shift (orbital_guard_.pre_diag/post_diag) prevents near-degenerate
    // active/virtual F_MO eigenvalues from swapping orbital identity across iterations:
    //   F_MO[i][i] += shift[i]     shift[i] = (k+1)*user_shift for active slot k, swapped
    //                               within a GVB pair (p,q) if F_p > F_q; flat virtual
    //                               shift = (n_act+1)*user_shift
    //   diagonalize F_MO -> U, eps_raw
    //   eps_raw[j] -= sum_i U[i][j]^2 * shift[i]        (undo the shift)
    //   MOM: reorder U columns / eps_raw by max_i |U[i][j]|
    std::vector<std::pair<int, int>> gvb_pairs;
    for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
        int g = sa_cassette_.geminals_active_gen(0)[u];
        const auto& orbs = sa_cassette_.geminal_templates()[g].orbitals;
        gvb_pairs.emplace_back(orbs[0], orbs[1]);
    }
    double applied_shift = orbital_guard_.pre_diag(Fmo_guard, N, Ncore_, active_mo_indices_, gvb_pairs,
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

    if (use_gvb_diis_ && mixing.n_swaps > 0) {
        gvb_diis_had_swaps_ = true;
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
    if (ctrl_state_.is_active() && !orbital_runtime_active_ && reel_.F_reks_MO) {
        // CFM runtime (native CFM or ORBITAL fallen back to CFM): Fa_=F^g carries an
        // artificial diag i*ls, so take eps from F_reks projected to the new MO basis:
        // F_phys_new = U^T reel_.F_reks_MO U, eps_i = u_i . (F_reks_MO U)[:,i]. Only the
        // diagonal is read; no full nmo x nmo congruence product is built.
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

    timer_off("REKS: form_C");
}

void REKS::form_C_read_guess_update() {
    timer_on("REKS: form_C_read_guess_update");
    ScopedStage _ss_fcrgu(reports(4) ? &scf_iter_times_ : nullptr, "form_C_read_guess_update");
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

// SAD-guess branch returns the RHF-style component energy of Da_; the converged
// branch returns the SA total E_total = E_SA + E_pen + vv10 (E_SA = sum_L C_L E_L),
// or the stand-on-best override when armed.
double REKS::compute_E() {
    if (sad_ && iteration_ <= 0) {
        double one_electron_E = 2.0 * Da_->vector_dot(H_);
        double kinetic_E = 2.0 * Da_->vector_dot(T_);
        double coulomb_E = 2.0 * Da_->vector_dot(J_);

        double XC_E = 0.0;
        double VV10_E = 0.0;
        if (functional_->needs_xc() && rv_potential_) {
            XC_E = rv_potential_->quadrature_values()["FUNCTIONAL"];
        }
        if (functional_->needs_vv10() && rv_potential_) {
            VV10_E = rv_potential_->quadrature_values()["VV10"];
        }

        double exchange_E = 0.0;
        double alpha = functional_->x_alpha();
        double beta_lrc = functional_->x_beta();

        if (functional_->is_x_hybrid() && K_) {
            exchange_E -= alpha * Da_->vector_dot(K_);
        }
        if (functional_->is_x_lrc()) {
            if (jk_->get_do_wK() && jk_->get_wcombine()) {
                exchange_E -= Da_->vector_dot(wK_);
            } else {
                exchange_E -= beta_lrc * Da_->vector_dot(wK_);
            }
        }

        double dashD_E = scalar_variable("-D Energy");

        energies_["Nuclear"] = nuclearrep_;
        energies_["Kinetic"] = kinetic_E;
        energies_["One-Electron"] = one_electron_E;
        energies_["Two-Electron"] = coulomb_E + exchange_E;
        energies_["XC"] = XC_E;
        energies_["VV10"] = VV10_E;
        energies_["-D"] = dashD_E;

        return nuclearrep_ + one_electron_E + coulomb_E + exchange_E + XC_E + VV10_E + dashD_E;
    }

    double kinetic_E = 2.0 * Da_->vector_dot(T_);
    energies_["Kinetic"] = kinetic_E;

    // E_SA = sum over the SCF pool of C_L * E_L (SI-only microstates are not
    // in sa_cassette_.microstates() and so do not contribute).
    double E_SA = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);

    // SA-decomposed energy components, weighted by C_L:
    //   E_1e_L = tr(Da*H) + tr(Db*H)
    //   E_2e_pure_L = 0.5*(Da*Fa + Db*Fb - E_1e_L) - 0.5*tr(D_L*V_xc_L)   (pure 2e, RHF/UHF convention)
    const bool needs_xc = functional_->needs_xc();
    double one_electron_E = 0.0;
    double two_electron_E = 0.0;
    double xc_E = 0.0;
    for (int L : sa_cassette_.microstates()) {
        double C_L = reel_.C_L[L];
        if (C_L == 0.0) continue;
        int alpha_idx = base_idx(L, /*alpha=*/true);
        int beta_idx  = base_idx(L, /*alpha=*/false);

        double E_1e_L = reel_.base_density_e1[alpha_idx] + reel_.base_density_e1[beta_idx];
        double E_Fa   = reel_.E_Fa_L[L];
        double E_Fb   = reel_.E_Fb_L[L];
        double E_2e_L_code = 0.5 * (E_Fa + E_Fb - E_1e_L);

        double tr_D_Vxc = 0.0;
        if (needs_xc && uv_potential_) {
            tr_D_Vxc = reel_.tr_D_Vxc_a_L[L] + reel_.tr_D_Vxc_b_L[L];
        }
        double E_2e_pure_L = E_2e_L_code - 0.5 * tr_D_Vxc;

        one_electron_E += C_L * E_1e_L;
        two_electron_E += C_L * E_2e_pure_L;
        if (needs_xc) {
            xc_E += C_L * reel_.E_xc_L[L];
        }
    }

    double dashD_E = scalar_variable("-D Energy");
    double E_pen   = compute_ipr_penalty_energy();
    double E_total = E_SA + E_pen + vv10_E_rv_;

    energies_["Nuclear"]      = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"]           = xc_E;
    energies_["-D"]           = dashD_E;
    energies_["REKS SA"]      = E_SA;
    energies_["IPR Penalty"]  = E_pen;
    energies_["VV10"]         = vv10_E_rv_;
    energies_["Total Energy"] = E_total;

    if (reports(5)) {
        outfile->Printf("\n  === REKS State-Averaged Energy ===\n");
        outfile->Printf("  E_SA    = %20.12f\n", E_SA);
        outfile->Printf("  E_pen   = %20.12f\n", E_pen);
        outfile->Printf("  E_total = %20.12f\n", E_total);
    }

    // Stand-on-best: report the best iterate's energy so the SCF gate sees one coherent
    // (E, Dnorm) pair.
    const double E_report = gvb_report_E_override_armed_ ? gvb_report_E_override_ : E_total;

    // dE latch: shifts prev <- cur only on the first compute_E() call of a new iteration_.
    if (gvb_gate_E_iter_ != iteration_) {
        gvb_gate_dE_valid_ = (gvb_gate_E_iter_ >= 0);
        gvb_gate_E_prev_ = gvb_gate_E_cur_;
        gvb_gate_E_iter_ = iteration_;
    }
    gvb_gate_E_cur_ = E_report;
    if (gvb_gate_dE_valid_) gvb_gate_dE_ = E_report - gvb_gate_E_prev_;

    return E_report;
}

void REKS::print_sa_diagnostics() {
    int n_sa = n_sa_microstates();
    if (n_sa < 2) return;
    if (!catalog_.data || reel_.E_L.empty()) return;

    double E_SA = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);

    // last_computed_dE_ must be set before prev_E_total_ is overwritten below.
    last_computed_dE_ = (iteration_ > 1 && std::abs(prev_E_total_) > 1e-10) ? (E_SA - prev_E_total_) : 0.0;
    if (iteration_ > 1 && std::abs(prev_E_total_) > 1e-10) {
        if (iteration_ > 2 && std::abs(prev_dE_) > 1e-10 && last_computed_dE_ * prev_dE_ < 0) {
            energy_oscillation_count_++;
            outfile->Printf("  *** E_SA OSCILLATION (total: %d) ***\n", energy_oscillation_count_);
        }
        prev_dE_ = last_computed_dE_;
    }
    prev_E_total_ = E_SA;
}

void REKS::print_iteration_diagnostics() {
    double E_SA = reks::scf::compute_E_SA(sa_cassette_, reel_.E_L, reel_.C_L);
    double dE = last_computed_dE_;

    outfile->Printf("  [REKS] iter=%d E_SA=%.14f dE=%+.6e", iteration_, E_SA, dE);

    // Per-sector gen-0 FON value + delta vs previous iteration. prev_fon_ stores
    // every sector's gen-0 FON concatenated in sector order (df is a printed delta
    // only, no gate); the sector tag is emitted only in a multi-sector run.
    const int n_sec = sa_cassette_.n_sectors();
    std::vector<double> cur_fon;
    int flat = 0;
    for (int sec = 0; sec < n_sec; ++sec) {
        const auto& gem0 = sa_cassette_.geminals_active(sec, 0);
        const int nps = static_cast<int>(gem0.size());
        if (nps == 0) continue;
        if (n_sec > 1) outfile->Printf(" sector%d", sec);
        for (int p = 0; p < nps; ++p) {
            double fp = reel_.fon_state[sec].layers[0][gem0[p]].p;
            double prev_fp = (flat < static_cast<int>(prev_fon_.size())) ? prev_fon_[flat] : fp;
            outfile->Printf(" f%d=%.8f df%d=%+.4e", p, fp, p, fp - prev_fp);
            cur_fon.push_back(fp);
            ++flat;
        }
    }

    // Fock eigenvalue gap eps[q] - eps[p] over the shared active spectrum; the gen-0
    // pair union is sector-agnostic (one orbital set).
    {
        double* eps = epsilon_a_->pointer(0);
        const auto& gen0_pool = sa_cassette_.geminals_active_gen(0);
        const int npair = static_cast<int>(gen0_pool.size());
        for (int p = 0; p < npair; ++p) {
            const auto& orbs = sa_cassette_.geminal_templates()[gen0_pool[p]].orbitals;
            int mo_p = active_mo_indices_[orbs[0]];
            int mo_q = active_mo_indices_[orbs[1]];
            outfile->Printf(" fgap%d=%+.6f", p, eps[mo_q] - eps[mo_p]);
        }
    }

    if (energy_oscillation_count_ > 0) outfile->Printf(" osc=%d", energy_oscillation_count_);
    outfile->Printf("\n");

    prev_fon_ = std::move(cur_fon);
}

SharedVector REKS::fon_occupation(int s) const {
    SharedVector occ = occupation_a();
    if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return occ;
    for (int g : sa_cassette_.geminals_active_gen(0)) {    // gen 0 = n
        const auto& tpl = sa_cassette_.geminal_templates()[g];
        if (tpl.scheme != 0) continue;                 // scheme 0 = canonical SCF pairs
        const int mo_p = active_mo_indices_[tpl.orbitals[0]];
        const int mo_q = active_mo_indices_[tpl.orbitals[1]];
        occ->set(0, mo_p, 0.5 * reel_.fon_state[s].layers[0][g].p);
        occ->set(0, mo_q, 0.5 * reel_.fon_state[s].layers[0][g].q);
    }
    return occ;
}

const reks::SIResult& REKS::primary_si_result() const {
    static const reks::SIResult empty;
    return si_results_.empty() ? empty : si_results_.front();
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

REKS::DipoleIntegralCache REKS::build_dipole_integral_cache() const {
    DipoleIntegralCache dip;
    const int N_act  = sa_cassette_.n_active_orbitals();
    const size_t N2  = static_cast<size_t>(N_act) * N_act;

    // AO dipole integrals + nuclear dipole.
    auto mints = std::make_shared<MintsHelper>(basisset_);
    std::vector<SharedMatrix> ao_dip = mints->ao_dipole();  // [0]=x, [1]=y, [2]=z

    for (int A = 0; A < molecule_->natom(); ++A) {
        double Z = molecule_->Z(A);
        dip.mu_nuc[0] += Z * molecule_->x(A);
        dip.mu_nuc[1] += Z * molecule_->y(A);
        dip.mu_nuc[2] += Z * molecule_->z(A);
    }

    // AO -> MO dipole. Must use reel_.Ca_fock_transform (same MOs as F_MO / active_mo_indices_).
    auto C = reel_.Ca_fock_transform ? reel_.Ca_fock_transform : Ca_;
    int nso = C->rowspi()[0];
    int nmo = C->colspi()[0];

    // Each axis is half-transformed onto the core+active columns of C (contiguous
    // prefix [0, n_occ)) and contracted back over those columns only; the full MO
    // dipole is never formed.
    const int n_occ = active_mo_indices_.empty() ? Ncore_ : active_mo_indices_.back() + 1;
    double** Cp_d = C->pointer(0);
    std::vector<double> T(static_cast<size_t>(nso) * n_occ);
    for (int axis = 0; axis < 3; ++axis) {
        // T = D_ao C_occ   (nso x n_occ)
        C_DGEMM('N', 'N', nso, n_occ, nso, 1.0, ao_dip[axis]->pointer(0)[0], nso, Cp_d[0], nmo,
                0.0, T.data(), n_occ);
        // Core contribution: 2 * sum_{p<Ncore_} d_pp.
        for (int p = 0; p < Ncore_; ++p)
            dip.mu_core[axis] +=
                2.0 * C_DDOT(nso, &Cp_d[0][p], nmo, T.data() + p, n_occ);
        // d_act[axis][p*N_act+q] = c_{a_p} . (D_ao c_{a_q}).
        dip.d_act[axis].resize(N2);
        for (int p = 0; p < N_act; ++p)
            for (int q = 0; q < N_act; ++q)
                dip.d_act[axis][p * N_act + q] =
                    C_DDOT(nso, &Cp_d[0][active_mo_indices_[p]], nmo,
                           T.data() + active_mo_indices_[q], n_occ);
    }
    return dip;
}

REKS::SIProperties REKS::compute_state_properties(
    const reks::studio::Cassette& cassette, const reks::SIResult& sir,
    const reks::rdm::DiabaticRdm& rho_diab, const DipoleIntegralCache& dip) const {
    timer_on("REKS: compute_state_properties");
    ScopedStage _ss_csp(reports(4) ? &post_scf_times_ : nullptr, "state_properties");
    SIProperties out;

    if (sir.n_states < 1 || sir.coeffs.empty()) {
        timer_off("REKS: compute_state_properties");
        return out;
    }
    // Diabatic-RDM availability is variant-global; absent it there are no
    // cross-state dipoles.
    if (!reks::studio::has_diabatic_rdm(catalog_)) {
        timer_off("REKS: compute_state_properties");
        return out;
    }

    const int n_rep   = sir.n_report();
    const int N_act   = sa_cassette_.n_active_orbitals();
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
    _ss_csp.release();
    if (reports(4)) log_post_scf_summary_("state_properties");
    return out;
}

void REKS::print_state_properties_header(const reks::studio::Cassette& cassette, int K) const {
    if (!reports(2)) return;
    std::string ml = reks::report::method_label(
        sa_cassette_,
        static_cast<int>(cassette.K_indices.size()));
    // Multi-sector runs tag the title with the SI section's spin manifold.
    if (sa_cassette_.n_sectors() > 1) {
        const int spin2 = catalog_.sector_entry(run_to_blob_[cassette.sector()]).spin2;
        ml += " " + reks::report::spin_tag(spin2);
    }
    outfile->Printf("\n  ==> %s State Properties (K=%d) <==\n\n", ml.c_str(), K);
}

void REKS::print_state_dipoles(const SIProperties& props, int K,
                               const reks::SIResult& sir) const {
    if (!reports(2)) return;
    (void)K;
    // props is already sized to the reported states.
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

void REKS::register_si_wfn_variables(int K, int s,
                                     const reks::SIResult& sir,
                                     const SIProperties& props,
                                     bool primary_alias,
                                     const reks::rdm::DiabaticRdm& rho_diab) {
    // Array names must not end in "DIPOLE" (Psi4 reshapes to (3,)); use "DIPOLES".
    const std::string suf = " K=" + std::to_string(K);
    const int n    = sir.n_states;         // config-space dimension (C/H/S matrices, diabatic RDM)
    const int n_rep     = sir.n_report();  // reported states; also the extent of props
    if (n < 1) return;

    // set_array_variable clones, so the bare alias costs a second full copy of what it names.
    auto record = [&](const std::string& base, SharedMatrix M) {
        set_array_variable(sector_psivar_name(s, base + suf), M);
        if (primary_alias) set_array_variable(sector_psivar_name(s, base), M);
    };

    // Bulk arrays carry the K-suffixed name alone; the alias would clone them.
    auto record_k_only = [&](const std::string& base, SharedMatrix M) {
        set_array_variable(sector_psivar_name(s, base + suf), M);
    };

    // Arrays carrying the config axis (coefficients, Hamiltonian, overlap, diabatic RDM)
    // scale as O(n^2) where every other array here is linear in n_rep; at large active
    // spaces they dominate the property dump. They ride the extended report.
    const bool export_config_space = reports(3);

    // SI eigensystem.
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
                pdm->set(0, I, a, props.perm_dipole[I][a] * pc_dipmom_au2debye);
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
                ddm->set(0, I - 1, a,
                        (props.perm_dipole[I][a] - props.perm_dipole[0][a]) * pc_dipmom_au2debye);
        record("SSR DIFFERENCE DIPOLES", ddm);
    }

    // Cross-state diabatic 1-RDM rho^(KL), K/L over config space. Row j is
    // [flat index (K*n + L)*N_act^2 + p*N_act + q, value], exact in double up to 2^53
    // (far above n^2*N_act^2). Exported nonzeros only; the catalog's declared cells can
    // vanish exactly at a given FON snapshot.
    if (export_config_space && reks::studio::has_diabatic_rdm(catalog_)) {
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

    // Adiabatic 1-RDM rho^(IJ): rows I*n_rep+J, cols p*N_act+q. rdm::adiabatic stores
    // the blocks at the same stride.
    if (reks::studio::has_diabatic_rdm(catalog_) && !props.rho_adiab.empty()) {
        const int N_act = sa_cassette_.n_active_orbitals();
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

    // Path B MOs (reel_.Ca_fock_transform) + Ncore.
    {
        auto C_path_b = reel_.Ca_fock_transform ? reel_.Ca_fock_transform : Ca_;
        record("REKS C PATH B", C_path_b);
        auto ncore_mat = std::make_shared<Matrix>("REKS NCORE" + suf, 1, 1);
        ncore_mat->set(0, 0, 0, static_cast<double>(Ncore_));
        record("REKS NCORE", ncore_mat);
    }

    // Per generation: FON array (ng) + META array (3 rows: scheme, orb_p, orb_q),
    // from this run sector's SI-relaxed FON snapshot; geminal templates are
    // variant-global.
    {
        const int ng = sa_cassette_.n_geminals();
        const reks::studio::FONSnapshot& si_fons = reel_.fon_state[s];
        for (int gen = 0; ng > 0 && gen < sa_cassette_.n_generations(); ++gen) {
            // Sector s populates only its own generations; skip layers it leaves empty.
            if (static_cast<int>(si_fons.layers[gen].size()) < ng) continue;
            const std::string base = reks::studio::fon_array_name(gen);
            auto fon_si = std::make_shared<Matrix>(base + suf, 1, ng);
            auto fon_si_meta = std::make_shared<Matrix>(base + " META" + suf, 3, ng);
            for (int g = 0; g < ng; ++g) {
                const auto& tpl = sa_cassette_.geminal_templates()[g];
                const auto& fon = si_fons.layers[gen][g];
                fon_si->set(0, 0, g, fon.p);
                fon_si_meta->set(0, 0, g, static_cast<double>(tpl.scheme));
                fon_si_meta->set(0, 1, g, static_cast<double>(tpl.orbitals[0]));
                fon_si_meta->set(0, 2, g, static_cast<double>(tpl.orbitals[1]));
            }
            set_array_variable(sector_psivar_name(s, base + suf), fon_si);
            set_array_variable(sector_psivar_name(s, base + " META" + suf), fon_si_meta);
        }
    }
}

reks::StateNaturalOrbitals
REKS::compute_state_natural(const reks::studio::Cassette& cassette,
                            const reks::SIResult& sir,
                            const std::vector<double>& rho_adiab) const {
    timer_on("REKS: compute_state_natural");
    ScopedStage _ss_csn(reports(4) ? &post_scf_times_ : nullptr, "state_natural");
    if (sir.n_states < 1 || !catalog_.data) {
        timer_off("REKS: compute_state_natural");
        return {};
    }
    if (!reks::studio::has_diabatic_rdm(catalog_) || rho_adiab.empty()) {
        timer_off("REKS: compute_state_natural");
        return {};
    }
    auto result = reks::naturals::compute(rho_adiab, cassette, sir);
    timer_off("REKS: compute_state_natural");
    _ss_csn.release();
    if (reports(4)) log_post_scf_summary_("state_natural");
    return result;
}

void REKS::print_state_natural(
    const reks::StateNaturalOrbitals& nos, int K,
    const reks::studio::Cassette& cassette, const reks::SIResult& sir,
    const reks::rdm::DiabaticRdm& rho_diab) const {
    if (!reports(2)) return;
    // n_rep (nos.n_states) bounds the state axis; n_cfg (sir.n_states) is the
    // independent config axis.
    const int n_rep = nos.n_states;
    const int n_cfg = sir.n_states;   // config-space dimension; stride of coeffs
    const int N    = nos.n_active;
    if (n_rep < 1 || N < 1) return;
    // active_orbital_labels() dereferences epsilon_a_.
    if (!epsilon_a_) return;

    std::string ml = reks::report::method_label(
        sa_cassette_,
        static_cast<int>(cassette.K_indices.size()));
    if (sa_cassette_.n_sectors() > 1) {
        const int spin2 = catalog_.sector_entry(run_to_blob_[cassette.sector()]).spin2;
        ml += "  " + reks::report::spin_tag(spin2);
    }

    const std::vector<int>& slice_K = cassette.K_indices;
    auto slice_name = [&](int J) -> std::string {
        const int idx = (static_cast<int>(slice_K.size()) == n_cfg) ? slice_K[J] : J;
        return reks::report::si_config_name(sa_cassette_, idx);
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

        // kSIBlockColumns configs per block; n_cfg beyond that wraps to more blocks.
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
    const std::vector<std::string> rk_lab = active_orbital_labels();
    outfile->Printf("      active SCF MO rows: ");
    for (int p = 0; p < N; ++p)
        outfile->Printf(" %s=%s", reks::report::orbital_label(p).c_str(), rk_lab[p].c_str());
    outfile->Printf("\n");
    for (int k = 0; k < n_rep; ++k) {
        outfile->Printf("\n      SSR state %d:\n", k);

        // Column display order: (MO, natural orbital) pairs are claimed in order of
        // decreasing |c|, putting the dominant coefficients on the diagonal. Ties
        // keep row-major order (reproducible).
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

    // Integrity: trace = n_act_e, asym ~ 1e-14, occ in [0,2]. Deviation = codegen bug.
    outfile->Printf("\n    %s Integrity check\n", reks::report::next_table_tag().c_str());
    outfile->Printf("    ---------------\n");
    const double n_act_e = static_cast<double>(sa_cassette_.n_electrons());
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

void REKS::register_natural_wfn_variables(
    int K, int s,
    const reks::StateNaturalOrbitals& nos,
    bool primary_alias) {
    const int n_rep = nos.n_states;
    const int N     = nos.n_active;
    if (n_rep < 1 || N < 1) return;

    const std::string suf_K = " K=" + std::to_string(K);

    auto record = [&](const std::string& base_with_state, SharedMatrix M) {
        set_array_variable(sector_psivar_name(s, base_with_state + suf_K), M);
        if (primary_alias) set_array_variable(sector_psivar_name(s, base_with_state), M);
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
    auto C_full = reel_.Ca_fock_transform ? reel_.Ca_fock_transform : Ca_;
    const int nso = C_full->rowspi()[0];

    auto C_active = std::make_shared<Matrix>("C_active", nso, N);
    {
        double** Cp = C_full->pointer(0);
        for (int p = 0; p < N; ++p) {
            const int mo = active_mo_indices_[p];
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

void REKS::print_mo_consistency() const {
    if (!reel_.Ca_fock_transform) return;

    const int nao = Ca_->nrow();
    const int nmo = Ca_->ncol();
    double** Cp = Ca_->pointer(0);
    double** Cp_fock = reel_.Ca_fock_transform->pointer(0);

    int act_min = active_mo_indices_[0], act_max = active_mo_indices_[0];
    for (int idx : active_mo_indices_) {
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

void REKS::print_f_reks_mo_debug() const {
    if (!reel_.F_reks_MO) return;

    const int nmo = reel_.F_reks_MO->ncol();
    double** Fp = reel_.F_reks_MO->pointer(0);

    int act_min = active_mo_indices_[0], act_max = active_mo_indices_[0];
    for (int idx : active_mo_indices_) {
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

void REKS::print_si_hamiltonian_overlap(const reks::studio::Cassette& cassette,
                                        const reks::SIResult& sir) const {
    if (sir.hamiltonian.empty()) return;
    const int nsi = sir.n_states;
    if (nsi <= 0) return;

    const std::vector<int>& K_indices = cassette.K_indices;

    reks::report::print_si_matrix(sa_cassette_, "Pre-diagonalization SI Hamiltonian",
                                    sir.hamiltonian, nsi, K_indices);
    if (!sir.overlap.empty()) {
        reks::report::print_si_matrix(sa_cassette_, "SI Overlap S",
                                        sir.overlap, nsi, K_indices);
    }
}

void REKS::print_fon_table(const reks::studio::Cassette& cassette) const {
    if (!reports(2)) return;
    const int n_active = static_cast<int>(active_mo_indices_.size());
    if (n_active <= 0 || !epsilon_a_) return;

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
        const auto& cass_fons = reel_.fon_state[cassette.sector()];
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

    // Energy-rank labels (slot order), matching the orbital-energy block.
    const std::vector<std::string> rk_lab = active_orbital_labels();
    outfile->Printf("      %-14s ", "Orbital");
    for (int k = 0; k < n_active; ++k)
        outfile->Printf(" %10s ", rk_lab[k].c_str());
    outfile->Printf("\n");

    // epsilon_a_ irrep 0: REKS is C1.
    outfile->Printf("      %-14s ", "Energy (Ha)");
    for (int k = 0; k < n_active; ++k) {
        outfile->Printf(" %10.6f ", epsilon_a_->get(0, active_mo_indices_[k]));
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

std::shared_ptr<REKS> REKS::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
    auto reks_wfn = std::make_shared<REKS>(wfn, functional_, wfn->options(), wfn->psio());

    if (Ca_) {
        reks_wfn->Ca_ = Ca_subset("AO", "ALL");
        reks_wfn->Cb_ = reks_wfn->Ca_;
    }
    if (Da_) {
        reks_wfn->Da_ = Da_subset("AO");
        reks_wfn->Db_ = reks_wfn->Da_;
    }
    if (Fa_) {
        reks_wfn->Fa_ = Fa_subset("AO");
        reks_wfn->Fb_ = reks_wfn->Fa_;
    }
    if (epsilon_a_) {
        reks_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nalphapi_, "AO", "ALL");
        reks_wfn->epsilon_b_ = reks_wfn->epsilon_a_;
    }

    // H_ and X_ come back empty from base construction; restore them.
    auto SO2AO = aotoso()->transpose();
    if (H_) reks_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) reks_wfn->X_->remove_symmetry(X_, SO2AO);

    return reks_wfn;
}

// IPR delocalization penalty:
//   E_pen  = lambda * sum_i sum_A p_A(i)^2
//   p_A(i) = sum_{mu in A} Ctil[mu,i]^2,  Ctil = S^{1/2} * Ca[:, active]

void REKS::build_ipr_static_data() {
    if (!S_)
        throw PSIEXCEPTION("[IPR] S_ not allocated");
    if (nirrep_ != 1)
        throw PSIEXCEPTION("[IPR] requires c1 symmetry (nirrep_=1)");

    S_half_ = S_->clone();
    S_half_->set_name("S^{1/2} for IPR penalty");
    S_half_->power(0.5, 1.0e-12);

    int nso = basisset_->nbf();
    ao2atom_.resize(nso);
    for (int mu = 0; mu < nso; ++mu)
        ao2atom_[mu] = basisset_->function_to_center(mu);

    n_atoms_ipr_ = molecule_->natom();
    atom_to_ao_.assign(n_atoms_ipr_, std::vector<int>());
    for (int mu = 0; mu < nso; ++mu)
        atom_to_ao_[ao2atom_[mu]].push_back(mu);
}

void REKS::compute_loewdin_populations(
    std::vector<std::vector<double>>& p_A,
    std::vector<std::vector<double>>& a_A) const {
    timer_on("REKS: compute_loewdin_populations");

    int nso    = basisset_->nbf();
    int n_act  = static_cast<int>(active_mo_indices_.size());
    int npairs = n_act * (n_act - 1) / 2;

    auto Ca_act = std::make_shared<Matrix>("Ca_active", nso, n_act);
    double** Cap      = Ca_->pointer(0);
    double** Ca_act_p = Ca_act->pointer(0);
    for (int mu = 0; mu < nso; ++mu)
        for (int i = 0; i < n_act; ++i)
            Ca_act_p[mu][i] = Cap[mu][active_mo_indices_[i]];

    // Reuse the memo when the active columns are bit-identical to the last build.
    if (ipr_pop_Ca_act_ && ipr_pop_Ca_act_->rowdim(0) == nso &&
        ipr_pop_Ca_act_->coldim(0) == n_act) {
        double** Cached = ipr_pop_Ca_act_->pointer(0);
        bool same = true;
        for (int mu = 0; mu < nso && same; ++mu)
            for (int i = 0; i < n_act; ++i)
                if (Ca_act_p[mu][i] != Cached[mu][i]) { same = false; break; }
        if (same) {
            p_A = ipr_pop_p_A_;
            a_A = ipr_pop_a_A_;
            timer_off("REKS: compute_loewdin_populations");
            return;
        }
    }

    // W = S^{1/2}*Ca (LOWDIN -> Ctil) or S*Ca (MULLIKEN).
    auto W = std::make_shared<Matrix>("W (IPR)", nso, n_act);
    double** Wp = W->pointer(0);

    const bool use_lowdin = (ipr_method_ == "LOWDIN");
    if (use_lowdin) {
        double** Sh = S_half_->pointer(0);
        C_DGEMM('N', 'N', nso, n_act, nso,
                1.0, Sh[0], nso,
                     Ca_act_p[0], n_act,
                0.0, Wp[0], n_act);
    } else {
        double** Sp = S_->pointer(0);
        C_DGEMM('N', 'N', nso, n_act, nso,
                1.0, Sp[0], nso,
                     Ca_act_p[0], n_act,
                0.0, Wp[0], n_act);
    }

    p_A.assign(n_act, std::vector<double>(n_atoms_ipr_, 0.0));
    a_A.assign(npairs, std::vector<double>(n_atoms_ipr_, 0.0));

    for (int A = 0; A < n_atoms_ipr_; ++A) {
        for (int mu : atom_to_ao_[A]) {
            if (use_lowdin) {
                for (int i = 0; i < n_act; ++i) {
                    double wi = Wp[mu][i];
                    p_A[i][A] += wi * wi;
                }
                int idx = 0;
                for (int i = 0; i < n_act; ++i) {
                    double wi = Wp[mu][i];
                    for (int j = i + 1; j < n_act; ++j) {
                        a_A[idx][A] += wi * Wp[mu][j];
                        ++idx;
                    }
                }
            } else {
                for (int i = 0; i < n_act; ++i) {
                    p_A[i][A] += Ca_act_p[mu][i] * Wp[mu][i];
                }
                int idx = 0;
                for (int i = 0; i < n_act; ++i) {
                    for (int j = i + 1; j < n_act; ++j) {
                        a_A[idx][A] += 0.5 * (Ca_act_p[mu][i] * Wp[mu][j]
                                            + Ca_act_p[mu][j] * Wp[mu][i]);
                        ++idx;
                    }
                }
            }
        }
    }

    // Closure invariant: sum_A p_A(i) = 1 for every active MO.
    for (int i = 0; i < n_act; ++i) {
        double s = 0.0;
        for (int A = 0; A < n_atoms_ipr_; ++A) s += p_A[i][A];
        if (std::abs(s - 1.0) > 1.0e-6) {
            timer_off("REKS: compute_loewdin_populations");
            throw PSIEXCEPTION(
                "[IPR] population closure failed for active MO; "
                "basis or S_half_ is corrupt");
        }
    }

    ipr_pop_Ca_act_ = Ca_act;
    ipr_pop_p_A_    = p_A;
    ipr_pop_a_A_    = a_A;
    timer_off("REKS: compute_loewdin_populations");
}

// Participation ratio of an active MO over the atomic populations p_A(i):
//   PR_i = 1 / sum_A p_A(i)^2
// Runs from 1 (all population on one centre) to n_atoms (even spread).
double REKS::active_pr_min() const {
    const int n_act = static_cast<int>(active_mo_indices_.size());
    if (n_act == 0 || n_atoms_ipr_ == 0 || !S_half_) return -1.0;

    std::vector<std::vector<double>> p_A, a_A;
    compute_loewdin_populations(p_A, a_A);

    double pr_min = -1.0;
    for (int i = 0; i < n_act; ++i) {
        double s2 = 0.0;
        for (int A = 0; A < n_atoms_ipr_; ++A) s2 += p_A[i][A] * p_A[i][A];
        if (s2 <= 0.0) return -1.0;
        const double pr = 1.0 / s2;
        if (pr_min < 0.0 || pr < pr_min) pr_min = pr;
    }
    return pr_min;
}

double REKS::compute_ipr_penalty_energy() const {
    if (lambda_ipr_ == 0.0) {
        last_ipr_total_ = 0.0;
        last_ipr_epen_  = 0.0;
        return 0.0;
    }
    std::vector<std::vector<double>> p_A, a_A;
    compute_loewdin_populations(p_A, a_A);

    int n_act = static_cast<int>(active_mo_indices_.size());
    double ipr = 0.0;
    for (int i = 0; i < n_act; ++i)
        for (int A = 0; A < n_atoms_ipr_; ++A)
            ipr += p_A[i][A] * p_A[i][A];

    last_ipr_total_ = ipr;
    last_ipr_epen_  = lambda_ipr_ * ipr;
    return last_ipr_epen_;
}

void REKS::add_ipr_contribution_to_F_gen(SharedMatrix F_gen) const {
    if (lambda_ipr_ == 0.0) return;

    std::vector<std::vector<double>> p_A, a_A;
    compute_loewdin_populations(p_A, a_A);

    int n_act = static_cast<int>(active_mo_indices_.size());
    double** Fg = F_gen->pointer(0);

    // Transpose entry [q,p] uses p_A(q), not p_A(p): the fill is asymmetric.
    // Diagonal stays zero (no shift of orbital energies).
    double gmax = 0.0;
    int idx = 0;
    for (int i = 0; i < n_act; ++i) {
        for (int j = i + 1; j < n_act; ++j) {
            int mi = active_mo_indices_[i];
            int mj = active_mo_indices_[j];
            double fij = 0.0, fji = 0.0;
            for (int A = 0; A < n_atoms_ipr_; ++A) {
                double a = a_A[idx][A];
                fij += p_A[i][A] * a;
                fji += p_A[j][A] * a;
            }
            fij *= 2.0 * lambda_ipr_;
            fji *= 2.0 * lambda_ipr_;
            Fg[mi][mj] += fij;
            Fg[mj][mi] += fji;
            double g = -2.0 * (fij - fji);
            gmax = std::max(gmax, std::abs(g));
            ++idx;
        }
    }
    last_ipr_gmax_ = gmax;
}

void REKS::add_ipr_contribution_to_hessian(
    std::vector<double>& hess, int N) const {

    if (lambda_ipr_ == 0.0) return;

    std::vector<std::vector<double>> p_A, a_A;
    compute_loewdin_populations(p_A, a_A);

    int n_act = static_cast<int>(active_mo_indices_.size());
    int n_rot = n_act * (n_act - 1) / 2;
    if (n_rot == 0) return;

    // IPR penalty Hessian in the orbital-rotation basis (rows/cols = active pairs
    // (p<q)):
    //   Hess[(pq),(rs)] = lambda * sum_A (H1 + H2)
    //   H1 = first-order-squared term (8 * a_pq * a_rs * Kronecker pattern)
    //   H2 = exp(K) second-order curvature term (double commutator)
    // M_A(i,j) below unifies the two population kinds: p_A on the diagonal, a_A
    // off-diagonal.

    // Pair index for i<j on a flat upper-triangle row-major layout.
    auto pidx = [n_act](int i, int j) {
        return i * n_act - i * (i + 1) / 2 + (j - i - 1);
    };
    auto Mget = [&](int A, int i, int j) -> double {
        if (i == j) return p_A[i][A];
        return (i < j) ? a_A[pidx(i, j)][A] : a_A[pidx(j, i)][A];
    };

    for (int p = 0; p < n_act; ++p) {
        for (int q = p + 1; q < n_act; ++q) {
            int row = pidx(p, q);
            for (int r = 0; r < n_act; ++r) {
                for (int s = r + 1; s < n_act; ++s) {
                    int col = pidx(r, s);
                    if (col < row) continue;

                    double H1 = 0.0, H2 = 0.0;
                    for (int A = 0; A < n_atoms_ipr_; ++A) {
                        double apq = a_A[pidx(p, q)][A];
                        double ars = a_A[pidx(r, s)][A];

                        // H^{(1)} = 8 * a_pq * a_rs * (d_pr + d_qs - d_ps - d_qr)
                        double sign = static_cast<double>(
                            (p == r) - (p == s) + (q == s) - (q == r));
                        H1 += 8.0 * apq * ars * sign;

                        // H^{(2)}, the exp(K) second-order term
                        // sum_i P^{ii} ([[P,E_pq],E_rs] + [[P,E_rs],E_pq])^{ii}: the
                        // double commutator symmetrized over both orderings of the
                        // rotation pair. G1 is the (r,s) ordering, G2 the (p,q)
                        // ordering; G1 == G2 when (p,q) == (r,s).
                        //   G1 = d_ps M(r,q) - d_pr M(s,q) + d_qs M(p,r) - d_qr M(p,s)
                        //   G2 = d_sq M(r,p) + d_rq M(p,s) - d_sp M(r,q) - d_rp M(q,s)
                        double G1 = 0.0;
                        if (p == s) G1 += Mget(A, r, q);
                        if (p == r) G1 -= Mget(A, s, q);
                        if (q == s) G1 += Mget(A, p, r);
                        if (q == r) G1 -= Mget(A, p, s);

                        double G2 = 0.0;
                        if (s == q) G2 += Mget(A, r, p);
                        if (r == q) G2 += Mget(A, p, s);
                        if (s == p) G2 -= Mget(A, r, q);
                        if (r == p) G2 -= Mget(A, q, s);

                        H2 += 2.0 * G1 * (p_A[q][A] - p_A[p][A])
                            + 2.0 * G2 * (p_A[s][A] - p_A[r][A]);
                    }

                    double Hpq_rs = lambda_ipr_ * (H1 + H2);
                    hess[row * N + col] += Hpq_rs;
                    if (row != col)
                        hess[col * N + row] += Hpq_rs;
                }
            }
        }
    }
}

void REKS::ipr_hessian_diag_aa(std::vector<double>& diag_out) const {
    int n_act = static_cast<int>(active_mo_indices_.size());
    int n_rot = n_act * (n_act - 1) / 2;
    diag_out.assign(n_rot, 0.0);
    if (lambda_ipr_ == 0.0 || n_rot == 0) return;

    // N = n_rot packs the active-active Hessian block densely; [idx][idx] is
    // then the per-pair IPR curvature.
    std::vector<double> H(static_cast<size_t>(n_rot) * n_rot, 0.0);
    add_ipr_contribution_to_hessian(H, n_rot);
    for (int idx = 0; idx < n_rot; ++idx)
        diag_out[idx] = H[static_cast<size_t>(idx) * n_rot + idx];
}

void REKS::print_ipr_iteration_diagnostic(int iter, double g_orb_max) const {
    int n_act = static_cast<int>(active_mo_indices_.size());

    std::vector<std::vector<double>> p_A, a_A;
    compute_loewdin_populations(p_A, a_A);
    double ipr_tot = 0.0;
    for (int i = 0; i < n_act; ++i)
        for (int A = 0; A < n_atoms_ipr_; ++A)
            ipr_tot += p_A[i][A] * p_A[i][A];
    last_ipr_total_ = ipr_tot;

    double ipr_can_ref = (n_atoms_ipr_ > 0)
                         ? static_cast<double>(n_act) / n_atoms_ipr_
                         : 0.0;
    double ipr_loc_ref = 0.5 * static_cast<double>(n_act);
    double threshold   = 0.5 * (ipr_can_ref + ipr_loc_ref);
    const char* basin  = (ipr_tot < threshold) ? "CAN" : "LOC";

    double lam   = lambda_ipr_;
    double E_pen = lam * ipr_tot;
    double gmax  = last_ipr_gmax_;
    double ratio = (g_orb_max > 1e-15) ? (gmax / g_orb_max) : 0.0;

    outfile->Printf(
        "  [IPR] iter=%d lambda=%.3e IPR=%.4f Epen=%.3e "
        "|g|max=%.3e g/gSA=%.2f basin=%s(can=%.3f loc=%.3f)\n",
        iter, lam, ipr_tot, E_pen, gmax, ratio,
        basin, ipr_can_ref, ipr_loc_ref);

    if (catalog_.data) {
        auto n_fon_vec = n_fon_vector();
        auto m_fon_vec = m_fon_vector();
        outfile->Printf("  [IPR] iter=%d n-FON=(", iter);
        for (size_t k = 0; k < n_fon_vec.size(); ++k)
            outfile->Printf("%.6f%s", n_fon_vec[k],
                            k + 1 < n_fon_vec.size() ? ", " : "");
        outfile->Printf(") m-FON=(");
        for (size_t k = 0; k < m_fon_vec.size(); ++k)
            outfile->Printf("%.6f%s", m_fon_vec[k],
                            k + 1 < m_fon_vec.size() ? ", " : "");
        outfile->Printf(")\n");
    }

    if (reports(5) && catalog_.data) {
        std::vector<int> mo_to_unit(n_act, -1);
        for (int u = 0; u < static_cast<int>(sa_cassette_.geminals_active_gen(0).size()); ++u) {
            int g = sa_cassette_.geminals_active_gen(0)[u];
            for (int orb : sa_cassette_.geminal_templates()[g].orbitals) mo_to_unit[orb] = u;
        }
        int idx = 0;
        for (int i = 0; i < n_act; ++i) {
            for (int j = i + 1; j < n_act; ++j) {
                double g = 0.0;
                for (int A = 0; A < n_atoms_ipr_; ++A)
                    g += (p_A[j][A] - p_A[i][A]) * a_A[idx][A];
                g *= 4.0 * lam;
                bool intra = (mo_to_unit[i] >= 0 &&
                              mo_to_unit[i] == mo_to_unit[j]);
                outfile->Printf(
                    "    [IPR] iter=%d pair(%d,%d) %s g=%.3e\n",
                    iter, active_mo_indices_[i], active_mo_indices_[j],
                    intra ? "intra" : "cross", g);
                ++idx;
            }
        }
    }

    // Warn only if penalty dominates SA pull (ratio>10) AND |g_IPR| is
    // SCF-relevant (>1e-3 Ha) -- the absolute gate suppresses false alarms
    // near convergence where g_SA -> 0.
    if (ratio > 10.0 && gmax > 1.0e-3) {
        outfile->Printf(
            "  [IPR] WARNING iter=%d: g_IPR=%.3e > 10x g_SA -- "
            "penalty may over-bias SCF.\n",
            iter, gmax);
    }
}

// IPR basin guard: IPR rising across three consecutive iterations, with the IPR past drift_lim,
// means the accelerator is crossing the basin divide.
//   can_ref   = n_act / n_atoms   (IPR of fully delocalized active MOs)
//   loc_ref   = n_act / 2         (IPR of active MOs each localized on 2 atoms)
//   drift_lim = can_ref + 0.6*(midpoint - can_ref), midpoint = (can_ref + loc_ref)/2
bool REKS::ipr_basin_guard_prestep() {
    bool ipr_guard_fires = false;
    if (lambda_ipr_ > 0.0 && n_atoms_ipr_ > 0 && ipr_prevprev_ > 0.0 && !ipr_guard_disabled_) {
        int    n_act_g   = static_cast<int>(active_mo_indices_.size());
        double can_ref   = static_cast<double>(n_act_g) / n_atoms_ipr_;
        double loc_ref   = 0.5 * static_cast<double>(n_act_g);
        double midpoint  = 0.5 * (can_ref + loc_ref);
        double drift_lim = can_ref + 0.6 * (midpoint - can_ref);

        bool rising_3 = (ipr_prevprev_ < ipr_prev_) && (ipr_prev_ < last_ipr_total_);
        bool in_zone  = (last_ipr_total_ > drift_lim);
        ipr_guard_fires = rising_3 && in_zone;
    }

    // Streak fail-safe: ipr_guard_streak_limit_ fires with no net IPR drop => basin unreachable,
    // disable.
    if (ipr_guard_fires) {
        if (ipr_guard_streak_ == 0) ipr_guard_streak_start_ipr_ = last_ipr_total_;
        ipr_guard_streak_++;
        if (ipr_guard_streak_ >= ipr_guard_streak_limit_ &&
            last_ipr_total_ >= ipr_guard_streak_start_ipr_) {
            ipr_guard_disabled_ = true;
            ipr_guard_fires     = false;
            if (reports(4))
                outfile->Printf(
                    "  [IPR-GUARD] iter=%d DISABLED: fired %d iters, IPR %.4f -> %.4f no net "
                    "decrease. Target basin unreachable from this trajectory; letting SCF converge "
                    "on its natural basin.\n",
                    iteration_, ipr_guard_streak_, ipr_guard_streak_start_ipr_, last_ipr_total_);
        }
    } else {
        ipr_guard_streak_           = 0;
        ipr_guard_streak_start_ipr_ = -1.0;
    }

    if (ipr_guard_fires && reports(4))
        outfile->Printf(
            "  [IPR-GUARD] iter=%d trigger: IPR rising %.4f -> %.4f -> %.4f, in transition zone; "
            "flushing DIIS history and using the raw base map\n",
            iteration_, ipr_prevprev_, ipr_prev_, last_ipr_total_);

    ipr_prevprev_ = ipr_prev_;
    ipr_prev_     = last_ipr_total_;
    return ipr_guard_fires;
}

void REKS::print_ipr_final_diagnostic() const {
    if (!reports(2)) return;
    if (lambda_ipr_ == 0.0) return;

    int n_act = static_cast<int>(active_mo_indices_.size());
    double ipr_can_ref = (n_atoms_ipr_ > 0)
                         ? static_cast<double>(n_act) / n_atoms_ipr_
                         : 0.0;
    double ipr_loc_ref = 0.5 * static_cast<double>(n_act);

    outfile->Printf("\n  === REKS IPR Penalty Final ===\n");
    outfile->Printf("  [IPR] final: IPR_total=%.4f  (can_ref=%.3f, "
                    "loc_ref=%.3f)\n",
                    last_ipr_total_, ipr_can_ref, ipr_loc_ref);
    outfile->Printf("  [IPR] final: E_pen=%.6e Ha (lambda=%.3e)\n",
                    lambda_ipr_ * last_ipr_total_, lambda_ipr_);
    const char* verdict =
        (last_ipr_total_ < 0.5 * (ipr_loc_ref + ipr_can_ref))
        ? "CANONICAL" : "LOCALIZED (WARNING)";
    outfile->Printf("  [IPR] final: basin verdict = %s\n", verdict);
}

// [TIME] logging helpers. Sub-stages nest inside top-level stages -- the
// breakdown is NOT a disjoint partition.

namespace {

const char* const kScfStageOrder[] = {
    "form_G",       "build_sa_focks",        "build_sa_focks_MO",
    "compute_sa_energies",
    "compute_weighting_factors",
    "form_F",       "build_generalized_fock",
    "form_C",       "combined_step",         "form_C_read_guess_update",
    "form_D",       "build_base_densities",
    "gvb_fon_step",
    "jk_compute_sa", "jk_compute_sad", "jk_compute_active_eri",
    "xc_compute_sa_rv", "xc_compute_sa_uv", "xc_compute_sad",
};

const char* const kPostScfStageOrder[] = {
    "compute_si",        "fon_relax",         "compute_SI_energies",
    "si_selector_demand", "si_selector_packs", "si_build_hamiltonian", "si_build_overlap",
    "si_diagonalize",
    "state_properties",  "state_natural",     "build_diabatic_rdm",
    "jk_compute_active_eri", "xc_compute_si",
};

const char* const kJkSubstages[] = {
    "jk_compute_sa", "jk_compute_sad", "jk_compute_active_eri",
};
const char* const kXcSubstages[] = {
    "xc_compute_sa_rv", "xc_compute_sa_uv", "xc_compute_sad", "xc_compute_si",
};

double sum_keys(const std::map<std::string, double>& bucket,
                const char* const* keys, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        auto it = bucket.find(keys[i]);
        if (it != bucket.end()) s += it->second;
    }
    return s;
}

void print_breakdown(const std::map<std::string, double>& bucket,
                     const char* const* keys, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        auto it = bucket.find(keys[i]);
        if (it == bucket.end() || it->second <= 0.0) continue;
        outfile->Printf("  %s=%.3fms", keys[i], it->second * 1000.0);
    }
}

}  // anonymous namespace

void REKS::log_iter_timings_(int iter, double iter_total_seconds) const {
    outfile->Printf("[TIME] iter=%d  total=%.3fms", iter, iter_total_seconds * 1000.0);
    print_breakdown(scf_iter_times_,
                    kScfStageOrder,
                    sizeof(kScfStageOrder) / sizeof(kScfStageOrder[0]));
    outfile->Printf("\n");
}

void REKS::log_scf_summary_() const {
    double scf_total = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - scf_start_time_)
                           .count();
    outfile->Printf("[TIME] SCF SUMMARY  iters=%d  total=%.3fs", n_iters_logged_, scf_total);
    for (auto key : kScfStageOrder) {
        auto it = scf_total_times_.find(key);
        if (it == scf_total_times_.end() || it->second <= 0.0) continue;
        outfile->Printf("  %s=%.3fs", key, it->second);
    }
    outfile->Printf("\n");
    log_jk_xc_summary_("SCF SUMMARY", scf_total_times_, scf_total);
}

void REKS::log_post_scf_summary_(const char* phase_label) const {
    // total = wall span since post-SCF start; stages below are nested and do not sum to it.
    double total = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - post_scf_start_time_)
                       .count();
    outfile->Printf("[TIME] POST-SCF (%s done)  total=%.3fms", phase_label, total * 1000.0);
    print_breakdown(post_scf_times_,
                    kPostScfStageOrder,
                    sizeof(kPostScfStageOrder) / sizeof(kPostScfStageOrder[0]));
    outfile->Printf("\n");
    log_jk_xc_summary_(phase_label, post_scf_times_, total);
}

void REKS::log_jk_xc_summary_(const char* phase_label,
                              const std::map<std::string, double>& bucket,
                              double phase_total_seconds) const {
    const std::size_t n_jk = sizeof(kJkSubstages) / sizeof(kJkSubstages[0]);
    const std::size_t n_xc = sizeof(kXcSubstages) / sizeof(kXcSubstages[0]);
    const double t_jk = sum_keys(bucket, kJkSubstages, n_jk);
    const double t_xc = sum_keys(bucket, kXcSubstages, n_xc);
    const double t_total = (phase_total_seconds > 0.0) ? phase_total_seconds : (t_jk + t_xc);
    const double t_other = std::max(0.0, t_total - t_jk - t_xc);
    auto pct = [t_total](double v) -> double {
        return (t_total > 0.0) ? (100.0 * v / t_total) : 0.0;
    };
    outfile->Printf("[TIME] %s f_jk/f_xc  total=%.3fs  jk=%.3fs (%.1f%%)  xc=%.3fs (%.1f%%)  other=%.3fs (%.1f%%)\n",
                    phase_label, t_total, t_jk, pct(t_jk), t_xc, pct(t_xc), t_other, pct(t_other));
    if (t_jk > 0.0) {
        outfile->Printf("[TIME] %s   jk: ", phase_label);
        for (std::size_t i = 0; i < n_jk; ++i) {
            auto it = bucket.find(kJkSubstages[i]);
            if (it == bucket.end() || it->second <= 0.0) continue;
            outfile->Printf(" %s=%.3fs", kJkSubstages[i], it->second);
        }
        outfile->Printf("\n");
    }
    if (t_xc > 0.0) {
        outfile->Printf("[TIME] %s   xc: ", phase_label);
        for (std::size_t i = 0; i < n_xc; ++i) {
            auto it = bucket.find(kXcSubstages[i]);
            if (it == bucket.end() || it->second <= 0.0) continue;
            outfile->Printf(" %s=%.3fs", kXcSubstages[i], it->second);
        }
        outfile->Printf("\n");
    }
}

std::vector<std::string> REKS::get_iter_accel_labels() const {
    std::vector<std::string> labels;
    if (iteration_ <= 0) return labels;

    if (ctrl_state_.is_active()) {
        if (diis_formulation_ != DIISFormulation::ORBITAL) {
            labels.emplace_back("CFM-GVB-DIIS");
        } else if (orbital_runtime_active_) {
            labels.emplace_back("ORB-GVB-DIIS");
        } else {
            labels.emplace_back("CFM-GVB-DIIS");
        }
    } else if (use_trah_) {
        labels.emplace_back("TRAH");
    } else if (use_plain_) {
        labels.emplace_back("PLAIN");
    }
    if (kdiis_active_) {
        labels.emplace_back("kDIIS");
    }
    return labels;
}

}  // namespace scf
}  // namespace psi
