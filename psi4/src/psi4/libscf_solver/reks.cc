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

/**
 * @file reks.cc
 * @brief Implementation of Restricted Ensemble Kohn-Sham (REKS) method.
 *
 * Current Implementation Status:
 * - Phase 1: REKS class structure inheriting from RHF
 * - All calculations delegate to RHF (functionally equivalent to RHF)
 * - Foundation for REKS-specific features in future phases
 *
 * Future Development Phases:
 * - Phase 2: Active space definition and fractional occupations
 * - Phase 3: Ensemble energy functional
 * - Phase 4: SI-SA-REKS (state-interaction state-averaged REKS)
 * - Phase 5: Analytical gradients
 */

#include "reks.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"

namespace psi {
namespace scf {

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : RHF(ref_wfn, func) {
    reks_common_init();
}

REKS::REKS(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
           Options& options, std::shared_ptr<PSIO> psio)
    : RHF(ref_wfn, func, options, psio) {
    reks_common_init();
}

REKS::~REKS() {}

void REKS::reks_common_init() {
    // Update wavefunction name to REKS
    name_ = "REKS";

    // REKS-specific initialization will go here in future phases:
    // - Active space setup
    // - Fractional occupation initialization
    // - State-averaging weights
    // - Ensemble density construction parameters

    // For now, REKS behaves identically to RHF
    // All computational methods are inherited from RHF
}

std::shared_ptr<REKS> REKS::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    // Create base wavefunction copy
    auto wfn = Wavefunction::c1_deep_copy(basis);

    // Create REKS wavefunction from the copy
    auto reks_wfn = std::make_shared<REKS>(wfn, functional_, wfn->options(), wfn->psio());

    // Copy matrices that REKS/RHF initializes
    // Include only persistent matrices (some are deleted in finalize())
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

    // H_ and X_ are reset in HF constructor, copy them over
    auto SO2AO = aotoso()->transpose();
    if (H_) reks_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) reks_wfn->X_->remove_symmetry(X_, SO2AO);

    return reks_wfn;
}

}  // namespace scf
}  // namespace psi
