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

#ifndef REKS_H
#define REKS_H

#include "hf.h"
#include "psi4/libreks/reks_si_types.h"
#include "psi4/libreks/reks_arch.h"
#include "psi4/libreks/reks_cassette.h"
#include "psi4/libreks/reks_reel.h"
#include "psi4/libreks/reks_convergence.h"
#include "psi4/libfock/v.h"
#include <vector>
#include <memory>
#include <utility>
#include <array>
#include <string>

namespace psi {
namespace reks {
class FonStep;
class GvbDiisDriver;
class TrahDriver;
namespace report {
struct ReportFrame;
}
namespace xc {
class MicrostateXC;
}
}
namespace scf {

/// Spin-restricted ensemble-referenced Kohn-Sham (REKS) SCF over the Catalog/Cassette/Reel active space,
/// one flow for every REKS(N,M).
class REKS : public HF {
   protected:
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
    SharedMatrix wK_;

    /// Offset-ABI view of the loaded catalog blob (base pointer + header).
    reks::studio::Catalog catalog_{};

    /// Runtime state (FONs, per-L Focks and energies, Lagrangians, base densities), sized from catalog_.
    reks::studio::Reel reel_;

    /// SCF ensemble pool over the run's sectors, extra determinants at the tail; dispatches over global L.
    reks::studio::Cassette sa_cassette_;

    /// SI cassettes: one per user SI_REKS_CONFIGS entry, in order.
    std::vector<reks::studio::Cassette> si_cassettes_;

    /// Call sheet over the SA cassette and si_cassettes_.
    reks::studio::CallSheet call_sheet_;

    /// FON blocks of the joint optimizer, one per (run sector, generation); offsets set by each pass.
    std::vector<reks::studio::FonBlock> fon_blocks_;

    /// Blob sector -> run sector ordinal; -1 when the blob sector is not in this run.
    std::vector<int> blob_to_run_sector_;

    /// Run sector ordinal -> blob sector index.
    std::vector<int> run_to_blob_;

    /// Builds fon_blocks_ from the SA cassette and (re)sizes the per-(s, gen) AHC FON smoothing state. Idempotent.
    void build_fon_blocks();

    /// Run sector ordinal owning global config index K; -1 if K is outside the run's sectors.
    [[nodiscard]] int run_sector_of_config(int K) const {
        const int b = catalog_.sector_of_config(K);
        return (b >= 0 && b < static_cast<int>(blob_to_run_sector_.size()))
            ? blob_to_run_sector_[b] : -1;
    }

    /// Catalog view bound to run sector s (local indices and names resolve in its block); a cheap copy.
    [[nodiscard]] reks::studio::Catalog catalog_sector_view(int s) const {
        return catalog_.bound_to(run_to_blob_[s]);
    }

    /// reel_.F_reks_MO = symmetrized sum_L C_L Fock_L (MO basis).
    void assemble_F_reks_MO();

    /// assemble_F_reks_MO(), then the pure-GGA virt-virt block of F_reks_MO from H + G.
    void assemble_F_reks();

    /// reel_.D_core and the one-electron trace of every SA occupation pattern into reel_.base_density_e1.
    void build_base_densities();

    /// Da_ = reel_.D_core + sum_{i in [Ncore_, nalphapi_[0])} c_i c_i^T; D_core must match Ca_.
    void form_Da_from_core();

    /// SA per-L Focks (AO + MO) for the SA cassette's microstate set; other L stay unfilled until post-SCF.
    void build_sa_focks();

    /// Alpha/beta base-density bitmask of microstate L: sum_i occ[i] << i over the active orbitals.
    [[nodiscard]] int base_idx(int L, bool alpha) const;

    int Ncore_ = -1;  ///< doubly occupied core orbitals
    std::vector<int> active_mo_indices_;

    bool si_computed_ = false;

    std::vector<reks::SIResult> si_results_;

    /// Primary (first) cassette result; empty while compute_si() has not run.
    [[nodiscard]] const reks::SIResult& primary_si_result() const;

    /// (ij|kl) tiles over the union of the SI coupling pairs and all off-diagonal active pairs.
    reks::ActiveEriTiles active_eri_pool_;

    /// Workspace aliased by the return value of build_generalized_fock().
    mutable SharedMatrix F_gen_workspace_;

    /// Ca_subset("SO","OCC") of the current iteration.
    SharedMatrix C_occ_cache_;

    /// Stepped state of the last form_C step (Ca, eps, FONs).
    SharedMatrix stage_Ca_;
    SharedVector stage_eps_;
    std::vector<reks::studio::FONSnapshot> stage_fons_;
    /// True between a form_C step and the head of the next form_G.
    bool staged_ = false;

    /// Spin-polarized XC for microstates whose D_alpha may differ from D_beta.
    std::unique_ptr<reks::xc::MicrostateXC> ms_xc_;

    /// RKS XC for the closed-shell SCF density Da_ (Db_ aliases Da_).
    std::shared_ptr<psi::RV> rv_potential_;

    /// FON lower bound per generation (REKS_<G>_FON): >= 0 confines generation g to [bound, 2], < 0 is off.
    std::array<double, reks::studio::kMaxGen> reks_fon_lower_{};

    /// E_SA / FON iteration trace.
    reks::ConvergenceTrace trace_;

    /// Adaptive level-shift + MOM-style reorder applied around the F_MO diagonalization.
    reks::OrbitalGuard orbital_guard_;

    /// TRAH combined solver (joint orbital+FON) toggle.
    bool use_trah_ = false;
    /// True only while GVB-DIIS has not activated; use_trah_ itself never resets after a
    /// TRAH warmup hands off.
    bool trah_live() const;
    /// Joint kappa+FON step driver with its cross-iteration state; present only when use_trah_.
    std::unique_ptr<reks::TrahDriver> trah_;

    double trah_gorb_norm_ = 0.0;              ///< ||g_orb||_F at the current Fock (TRAH convergence gate)

    /// Base-sector ordinal of the SA cassette.
    int sa_sector() const { return sa_cassette_.sector(); }
    /// SA-context FON snapshot (reel_.fon_state[sa_sector()]).
    const reks::studio::FONSnapshot& sa_fons() const { return reel_.fon_state[sa_sector()]; }
    reks::studio::FONSnapshot&       sa_fons()       { return reel_.fon_state[sa_sector()]; }

    /// GVB-DIIS (Muller, R. P. et al. J. Chem. Phys. 1994, 100, 1226).
    bool use_gvb_diis_ = false;                ///< Use GVB-DIIS convergence accelerator
    bool use_plain_ = false;                   ///< No accelerator: plain F_reks diag + FON Newton
    /// GVB-DIIS driver with its cross-iteration state; present only when use_gvb_diis_.
    std::unique_ptr<reks::GvbDiisDriver> gvb_;
    bool cosx_ = false;                        ///< SCF_TYPE contains COSX
    std::string cosx_grid_;                    ///< COSX grid of the previous Fock build; empty = none yet
    SharedMatrix F_gen_;                       ///< Cached generalized Fock for reuse within the iter

    /// FON Newton relaxation with its seed state; relaxes the FON whenever TRAH does not.
    std::unique_ptr<reks::FonStep> fon_step_;

    /// Frobenius orbital gradient -> rms per matrix element (D_CONVERGENCE units).
    double gorb_rms(double gorb_fro) const { return gorb_fro / nmopi_[0]; }

    /// Per-run convergence quality: wfn variables + the [QUALITY] trace line.
    void report_convergence_quality();

    /// Orbital-gradient norm over all rotation pairs (aa, ca, cv, av). frozen flags aa pairs
    /// to exclude, over the i<j packing; empty measures every rotation.
    double compute_full_gradient_norm(SharedMatrix F_gen, const std::vector<char>& frozen) const;

    /// Non-variational active rotations over the i<j packing: 1 for same-role orbitals of two geminals;
    /// all zero with fewer than two open geminals.
    std::vector<char> same_role_aa_frozen() const;

    /// Active-rotation-invariant AO Fock h + J[2 D_core + q sum_k D_k] - x_alpha K[D_core + q/2 sum_k D_k],
    /// q = N_act/n_act, from the reel J/K cache (null without it).
    SharedMatrix build_uniform_active_fock() const;

    /// Active rotation to the eigenvectors of build_uniform_active_fock() per role block (ascending, positive
    /// diagonal); null without J/K cache or two open geminals. Optional: eigenvalues, off-diagonal norm.
    SharedMatrix same_role_canonicalizer(std::vector<double>* eigs = nullptr,
                                         double* offdiag_norm = nullptr) const;

    /// Rotates the active orbitals of Ca_ by same_role_canonicalizer() and sets their eps to
    /// the Fock diagonal. No-op when that returns null.
    void canonicalize_same_role_frame();

    /// READ guess of the gen-0 FONs per run sector; with guess_Ca_ it marks the loaded basis as trusted.
    std::vector<SharedMatrix> guess_fon_;

    /// READ guess of the FONs of generation gen >= 1 per run sector, null when absent.
    std::vector<std::array<SharedMatrix, reks::studio::kMaxGen>> guess_fon_upper_;

    /// True when any run sector carries a gen-0 READ FON guess (trusted-basis sentinel).
    [[nodiscard]] bool has_read_fon_guess() const {
        for (const auto& m : guess_fon_) if (m) return true;
        return false;
    }
    /// Clear every sector's READ FON guess (boundary-trap escape).
    void clear_read_fon_guess() {
        for (auto& m : guess_fon_) m.reset();
    }

    /// Allocates SCF matrices with beta aliased to alpha (REKS is restricted);
    /// sets the report level from REKS_REPORT_LEVEL, then dispatches to subclass_init().
    void reks_common_init();

    /// iter<=0 form_C on a trusted READ guess: no-op on Ca_, copies F_MO diagonal
    /// into epsilon_a_.
    void form_C_read_guess_update();

    /// Copies the live Ca_, epsilon_a_ and reel_.fon_state into the stage buffers.
    void stage_in();
    /// Swaps live state and stage buffers in place and sets staged_.
    void stage_out();
    /// Orbital step of form_C for iteration > 0: ORB-DIIS rotation, or TRAH step and F_MO diag.
    void form_C_step(double shift);

    /// form_D, then the SA microstate Focks, C_L and E_L of the live state.
    void evaluate_frame();

    /// Microstate energies E_L.
    void compute_sa_energies();

    /// Microstate weighting factors C_L.
    void compute_weighting_factors();

    /// Builds each SI cassette's Hamiltonian from the converged orbitals; no-op once si_computed_.
    void compute_si();

    void setup_potential() override;  ///< RV potential plus the microstate XC evaluator
    void save_density_and_energy() override;
    double compute_initial_E() override;
    /// READ propagates the full Ca (incl. active virtuals) without touching nalphapi_/nalpha_;
    /// other guess types delegate to HF::guess().
    void guess() override;
    /// Complete C to an S-orthonormal basis: MGS keeps the independent loaded columns in front, X fills the
    /// rest; returns the kept count.
    static int complete_orbital_basis(SharedMatrix C, SharedMatrix S, SharedMatrix X,
                                      int n_loaded, int h);
    void form_D() override;
    void form_G() override;
    void form_F() override;
    void form_initial_F() override;  ///< SAD iteration: H + G of the SAD density
    /// TRAH/GVB-DIIS: rms == ||g_orb||_F / nmo; on a GVB-DIIS landing at MAXITER, of the landed frame.
    SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso) override;
    void form_C(double shift = 0.0) override;
    double compute_E() override;
    void finalize() override;         ///< Release the microstate XC scratch
    void print_orbitals() override;   ///< 3-section print: Core / Active (with FON) / Virtual

    /// Const snapshot of the SCF state for the libreks report/property functions.
    [[nodiscard]] reks::report::ReportFrame report_frame() const;

    /// (ij|kl) MO ERI tiles for the canonical pairs (k<=l) in rank1_pairs.
    /// Precondition: rank1_pairs canonical and strictly sorted.
    [[nodiscard]] reks::ActiveEriTiles compute_active_mo_eri(
        const std::vector<int>& active_mo_indices,
        const std::vector<std::pair<int,int>>& rank1_pairs);

    /// F_gen[p,q] = sum_L C_L (n^a_L[p] Fa_L[p,q] + n^b_L[p] Fb_L[p,q]); gradient G[i,a] = 2 (F_gen[i,a] - F_gen[a,i]).
    SharedMatrix build_generalized_fock() const;

    /// reel_.F_acc_MO = MO transform of sum_L C_L (Fa_AO_L + Fb_AO_L).
    void build_F_acc_MO();

   public:
    REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    ~REKS() override;

    std::shared_ptr<VBase> V_potential() const override { return rv_potential_; }

    /// READ guess of the gen-0 FONs of sector s (one per active geminal), clamped 1e-8 inside {0, 2}.
    void guess_fon(SharedMatrix fon, int s = 0) {
        if (s >= static_cast<int>(guess_fon_.size())) guess_fon_.resize(s + 1);
        guess_fon_[s] = fon;
    }

    /// READ guess of the generation gen >= 1 FONs of sector s; size must equal
    /// sa_cassette_.geminals_active_gen(gen).size().
    void guess_fon_upper(int gen, SharedMatrix fon, int s = 0) {
        if (s >= static_cast<int>(guess_fon_upper_.size())) guess_fon_upper_.resize(s + 1);
        guess_fon_upper_[s][gen] = fon;
    }

    /// Total microstate count (Catalog basis).
    [[nodiscard]] int n_microstates() const {
        return sa_cassette_.n_catalog_microstates();
    }

    /// SA microstate count = cardinality of sa_cassette_'s F_MO build set.
    [[nodiscard]] int n_sa_microstates() const {
        return static_cast<int>(sa_cassette_.microstates().size());
    }

    [[nodiscard]] double get_microstate_energy(int L) const {
        return reel_.E_L[L];
    }

    /// Occupation label of microstate L, e.g. "|a,b,c,a_,b_,d_|" (active MOs a.., beta barred).
    [[nodiscard]] std::string get_microstate_label(int L) const;

    /// f_interp of n-geminal g of run sector s.
    [[nodiscard]] double get_f_value(int g = 0, int s = 0) const;

    /// SI Hamiltonian before diagonalization.
    [[nodiscard]] SharedMatrix SI_hamiltonian() const;
    /// SI eigenvalues (adiabatic state energies).
    [[nodiscard]] SharedVector SI_energies() const;
    /// SI eigenvectors, rows = states.
    [[nodiscard]] SharedMatrix SI_coefficients() const;
    [[nodiscard]] SharedMatrix SI_overlap() const;

    /// Bonding-orbital FON of geminal g, generation gen (0=n, 1=m, 2=u, ...), run sector s.
    [[nodiscard]] double get_fon_p(int g, int gen = 0, int s = 0) const {
        if (gen < 0 || gen >= reks::studio::kMaxGen) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        const auto& layer = reel_.fon_state[s].layers[gen];
        return (g >= 0 && g < static_cast<int>(layer.size())) ? layer[g].p : 0.0;
    }

    /// Antibonding-orbital FON of geminal g, generation gen (0=n, 1=m, 2=u, ...), run sector s.
    [[nodiscard]] double get_fon_q(int g, int gen = 0, int s = 0) const {
        if (gen < 0 || gen >= reks::studio::kMaxGen) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        const auto& layer = reel_.fon_state[s].layers[gen];
        return (g >= 0 && g < static_cast<int>(layer.size())) ? layer[g].q : 0.0;
    }

    /// Per-spin MO occupations with sector 0's n-FON written onto its active n-geminal orbitals.
    [[nodiscard]] SharedVector occupation_a() const override;
    [[nodiscard]] SharedVector occupation_b() const override { return occupation_a(); }

    [[nodiscard]] double get_lagrangian(int idx) const {
        const auto& l = reel_.lagrangians;
        return (idx >= 0 && idx < static_cast<int>(l.size())) ? l[idx] : 0.0;
    }

    /// E_K = sum_L C_L(K, n) E_L at run sector s's FONs.
    [[nodiscard]] double get_config_energy(int K, int s = 0) const {
        if (K < 0 || K >= sa_cassette_.n_configs()) return 0.0;
        if (s < 0 || s >= static_cast<int>(reel_.fon_state.size())) return 0.0;
        std::vector<double> C_K(sa_cassette_.n_catalog_microstates(), 0.0);
        sa_cassette_.compute_C_L(K, reel_.fon_state[s], C_K);
        double E = 0.0;
        for (int L : sa_cassette_.config_microstate_writes(K))
            E += C_K[L] * reel_.E_L[L];
        return E;
    }

    /// Public per-sector FON vector accessors (active bonding p values of run sector s,
    /// generation gen). set mutates the live snapshot without re-running SCF.
    [[nodiscard]] std::vector<double> get_fon_vector(int s, int gen) const {
        return reks::studio::fon_vector(reel_, sa_cassette_, s, gen);
    }
    void put_fon_vector(int s, int gen, const std::vector<double>& fons) {
        reks::studio::set_fon_vector(reel_, sa_cassette_, s, gen, fons);
    }

    /// Total geminal slots in the catalog (generation-independent).
    [[nodiscard]] int n_pairs() const { return sa_cassette_.n_geminals(); }

    /// Scheme tag of geminal slot g (0 = SCF-optimized, >0 = post-SCF).
    [[nodiscard]] int get_n_scheme(int g) const {
        return (g >= 0 && g < sa_cassette_.n_geminals())
            ? sa_cassette_.geminal_templates()[g].scheme : -1;
    }

    /// Total configuration count.
    [[nodiscard]] int n_configs() const { return sa_cassette_.n_configs(); }

    /// Number of Lagrangian elements.
    [[nodiscard]] int n_lagrangians() const {
        return sa_cassette_.n_lagrangian_pairs();
    }

    /// Dimension of the SI Hamiltonian (variant config count, HS excluded).
    [[nodiscard]] int n_si_states() const { return sa_cassette_.n_configs(); }

    /// Number of doubly-occupied core orbitals (MOs below the active space).
    [[nodiscard]] int get_Ncore() const { return Ncore_; }

    [[nodiscard]] int n_active_orbitals() const {
        return sa_cassette_.n_active_orbitals();
    }

    /// MO indices of active orbitals (length n_active_orbitals).
    [[nodiscard]] std::vector<int> get_active_mo_indices() const {
        return active_mo_indices_;
    }

    /// Active accelerators of this SCF iteration, in {"CFM-GVB-DIIS","ORB-GVB-DIIS","TRAH","PLAIN","kDIIS"};
    /// empty at iter<=0.
    [[nodiscard]] std::vector<std::string> iteration_labels() const override;

    /// True while the GVB-DIIS (dE, rms) gate streak is short of REKS_GVB_GATE_STREAK or on a landing
    /// iteration; false with the controller inactive.
    [[nodiscard]] bool convergence_held() const override;
};

}  // namespace scf
}  // namespace psi

#endif  // REKS_H
