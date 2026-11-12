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

#ifndef REKS_TRAH_H
#define REKS_TRAH_H

/// @file reks_trah.h
/// @brief TRAH joint orbital+FON step driver and its cross-iteration state.

#include <array>
#include <vector>

#include "reks_arch.h"
#include "reks_cassette.h"
#include "reks_reel.h"
#include "reks_trah_solver.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace psi {
namespace reks {

class ActiveEriTiles;

/// Inputs of one TrahDriver::step call; the references are valid for that call only.
struct TrahStepContext {
    int iteration;
    const std::vector<int>& active;                          ///< active MO indices (absolute, 0-based)
    const studio::Cassette& sa_cassette;
    /// E_L, C_L, F_MO_arows_*_L read; fon_state and TRAH scratch written.
    studio::Reel& reel;
    const std::vector<studio::FonBlock>& fon_blocks;         ///< (s, gen) blocks; offsets are reassigned by the driver
    const std::array<double, studio::kMaxGen>& fon_lower;    ///< REKS_<G>_FON per generation, < 0 = unset
    SharedMatrix F_gen;                                      ///< generalized Fock at the current MOs
    /// Active-MO ERI tiles for the exact H_oo diagnostic; empty skips it.
    const ActiveEriTiles* hess_eri;
    double x_alpha;                                          ///< hybrid exact-exchange coefficient
    SharedMatrix Ca;                                         ///< MO coefficients, rotated in place
    SharedVector epsilon;                                    ///< orbital energies (intra-pair damping)
    bool read_guess;                                         ///< a READ FON guess is restored on entry
    bool has_guess_Ca;                                       ///< READ MO guess present
    bool ctrl_active;                                        ///< GVB-DIIS controller active
    int n_sa_microstates;
    std::vector<char> aa_frozen;                             ///< determined aa rotations (i<j packing), not stepped
    bool clear_read_guess = false;                           ///< OUT: READ FON guess consumed by a boundary-trap escape
};

/// Joint kappa+FON trust-region step: AHC shift, BFGS H_oo, kDIIS, Aitken, intra-pair damping,
/// READ boundary-trap escape.
class TrahDriver {
   public:
    /// Fresh trust-region state with the starting radius.
    void configure();

    /// Size the per-block AHC FON fractions to fon_blocks.
    void resize_fon_blocks(const std::vector<studio::FonBlock>& fon_blocks);

    /// One step: rotates ctx.Ca in place and writes ctx.reel.fon_state.
    void step(TrahStepContext& ctx);

    /// kDIIS extrapolation engaged at the last step.
    bool kdiis_active() const { return kdiis_active_; }

   private:
    TRAHState state_;

    /// AHC: rho-driven PI diagonal shift of H_oo and of each FON block.
    double ahc_rho_smooth_ = 1.0;                 ///< EMA of rho
    double ahc_pred_orb_frac_smooth_ = 0.5;       ///< EMA orbital fraction of the predicted decrease
    /// EMA FON fraction per block, seeded 0.5 (gen 0) or 0 (gen >= 1).
    std::vector<double> ahc_pred_fon_frac_smooth_;
    double ahc_integral_ = 0.0;                   ///< I-term
    static constexpr double AHC_BETA_ = 0.3;      ///< EMA weight
    static constexpr double AHC_I_RATE_ = 0.15;
    static constexpr double AHC_I_MAX_ = 20.0;    ///< I-term cap

    double prev_gorb_at_boundary_ = -1.0;      ///< ||g_orb|| at the previous boundary iteration (-1 = unset)
    std::vector<double> prev_orb_step_;        ///< previous raw orbital step (Aitken)

    /// READ guess with all FON pinned at a bound gave a zero step; the next call resets FON.
    bool read_boundary_trap_ = false;

    /// BFGS H_oo: rank-2 update with Powell damping, started from the diagonal H_oo.
    std::vector<double> bfgs_B_;           ///< n_rot x n_rot, row-major
    std::vector<double> bfgs_prev_g_;      ///< previous gradient (N)
    std::vector<double> bfgs_prev_step_;   ///< previous applied step (N)
    bool bfgs_initialized_ = false;
    int bfgs_update_count_ = 0;            ///< successful updates
    int consecutive_bad_rho_ = 0;          ///< consecutive rho < 0; BFGS reset at 2

    /// kDIIS on kappa at boundary stagnation (FON locked).
    /// Fischer, T. H.; Almlof, J. J. Phys. Chem. 1992, 96, 9768.
    std::vector<std::vector<double>> kdiis_g_history_;  ///< ring buffer of orbital gradients (n_rot each)
    int kdiis_count_ = 0;              ///< stored gradients
    int kdiis_oldest_ = 0;             ///< ring-buffer index of the oldest
    bool kdiis_active_ = false;
    int kdiis_stag_count_ = 0;         ///< consecutive boundary-stagnation iterations
    int kdiis_energy_fail_ = 0;        ///< consecutive energy increases while active
    bool kdiis_used_prev_ = false;
    static constexpr int KDIIS_MAX_VEC = 6;
    static constexpr int KDIIS_MIN_VEC = 3;    ///< minimum stored gradients to extrapolate
    static constexpr int KDIIS_TRIGGER = 3;    ///< stagnation iterations that activate kDIIS
};

}  // namespace reks
}  // namespace psi

#endif  // REKS_TRAH_H
