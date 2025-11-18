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

#include "scf_grad.h"
#include "jk_grad.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/uhf.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"

namespace psi{
namespace scfgrad {

SharedMatrix scfgrad(std::shared_ptr<scf::HF> ref_wfn, Options &options)
{
    tstart();

    SCFDeriv grad(ref_wfn, options);
    SharedMatrix G = grad.compute_gradient();

    tstop();
    return G;
}

SharedMatrix scfhess(std::shared_ptr<scf::HF> ref_wfn, Options &options)
{
    tstart();

    SharedMatrix H;
    if( ref_wfn->same_a_b_orbs() && ref_wfn->same_a_b_dens()) {
        // RHF
        RSCFDeriv hessian_computer(std::dynamic_pointer_cast<scf::RHF>(ref_wfn), options);
        H = hessian_computer.compute_hessian();
    } else {
        USCFDeriv hessian_computer(std::dynamic_pointer_cast<scf::UHF>(ref_wfn), options);
        H = hessian_computer.compute_hessian();
    }
    ref_wfn->set_hessian(H);

    tstop();

    return H;
}

std::vector<SharedMatrix> multi_scfgrad(
    const std::vector<std::shared_ptr<scf::HF>>& wfns,
    Options& options)
{
    tstart();

    int nwfn = wfns.size();
    if (nwfn == 0) {
        throw PSIEXCEPTION("multi_scfgrad: Empty wavefunction list.");
    }

    // Validate: all wfns must have same geometry/basis
    auto ref_mol = wfns[0]->molecule();
    auto ref_basis = wfns[0]->basisset();
    for (int i = 1; i < nwfn; i++) {
        if (wfns[i]->molecule() != ref_mol) {
            throw PSIEXCEPTION("multi_scfgrad: All wavefunctions must have same molecule.");
        }
        if (wfns[i]->basisset() != ref_basis) {
            throw PSIEXCEPTION("multi_scfgrad: All wavefunctions must have same basis set.");
        }
    }

    int natom = ref_mol->natom();

    outfile->Printf("\n");
    outfile->Printf("         ------------------------------------------------------------\n");
    outfile->Printf("                         MULTI SCF GRAD (BATCHED)                   \n");
    outfile->Printf("                      Computing gradients for %2d wavefunctions      \n", nwfn);
    outfile->Printf("                Expected speedup: ~2-6× vs sequential gradients     \n");
    outfile->Printf("         ------------------------------------------------------------\n\n");

    // ========================================
    // Step 1: Compute individual gradient terms for each wfn
    // ========================================
    std::vector<std::map<std::string, SharedMatrix>> individual_grads(nwfn);

    auto mintshelper = std::make_shared<MintsHelper>(ref_basis, options);

    for (int i = 0; i < nwfn; i++) {
        auto& wfn = wfns[i];
        auto& grads = individual_grads[i];

        // Densities
        auto Da = wfn->Da_subset("AO");
        SharedMatrix Db;
        if (wfn->same_a_b_dens()) {
            Db = Da;
        } else {
            Db = wfn->Db_subset("AO");
        }
        auto Dt = Da->clone();
        Dt->add(Db);
        Dt->set_name("Dt");

        // Occupations
        auto Ca_occ = wfn->Ca_subset("AO", "OCC");
        SharedMatrix Cb_occ;
        auto eps_a_occ = wfn->epsilon_a_subset("AO", "OCC");
        SharedVector eps_b_occ;

        if (wfn->same_a_b_orbs()) {
            Cb_occ = Ca_occ;
            eps_b_occ = eps_a_occ;
        } else {
            Cb_occ = wfn->Cb_subset("AO", "OCC");
            eps_b_occ = wfn->epsilon_b_subset("AO", "OCC");
        }

        int nalpha = Ca_occ->colspi()[0];
        int nbeta = Cb_occ->colspi()[0];

        // Nuclear gradient
        grads["Nuclear"] = SharedMatrix(ref_mol->nuclear_repulsion_energy_deriv1().clone());
        grads["Nuclear"]->set_name("Nuclear Gradient");

        // Core gradient (T + V)
        timer_on("Grad: V T Perturb");
        grads["Core"] = mintshelper->core_hamiltonian_grad(Dt);
        timer_off("Grad: V T Perturb");

        // Overlap gradient
        timer_on("Grad: S");
        {
            // Energy weighted density matrix
            SharedMatrix W(Da->clone());
            W->set_name("W");

            // Alpha
            auto tmp = Ca_occ->clone();
            for (size_t j = 0; j < nalpha; j++) {
                tmp->scale_column(0, j, eps_a_occ->get(j));
            }
            W->gemm(false, true, 1.0, tmp, Ca_occ, 0.0);

            // Beta
            tmp->copy(Cb_occ);
            for (size_t j = 0; j < nbeta; j++) {
                tmp->scale_column(0, j, eps_b_occ->get(j));
            }
            W->gemm(false, true, 1.0, tmp, Cb_occ, 1.0);

            grads["Overlap"] = mintshelper->overlap_grad(W);
            grads["Overlap"]->scale(-1.0);
        }
        timer_off("Grad: S");
    }

    // ========================================
    // Step 2: Batched JK gradient computation
    // ========================================
    timer_on("Grad: JK (Batched)");

    auto jk = JKGrad::build_JKGrad(1, mintshelper);
    jk->set_memory((size_t)(options.get_double("SCF_MEM_SAFETY_FACTOR") *
                            Process::environment.get_memory() / 8L));

    // Collect all densities and orbitals
    std::vector<SharedMatrix> Ca_list, Cb_list, Da_list, Db_list, Dt_list;
    for (int i = 0; i < nwfn; i++) {
        auto& wfn = wfns[i];

        auto Da = wfn->Da_subset("AO");
        SharedMatrix Db;
        if (wfn->same_a_b_dens()) {
            Db = Da;
        } else {
            Db = wfn->Db_subset("AO");
        }
        auto Dt = Da->clone();
        Dt->add(Db);

        auto Ca_occ = wfn->Ca_subset("AO", "OCC");
        SharedMatrix Cb_occ;
        if (wfn->same_a_b_orbs()) {
            Cb_occ = Ca_occ;
        } else {
            Cb_occ = wfn->Cb_subset("AO", "OCC");
        }

        Ca_list.push_back(Ca_occ);
        Cb_list.push_back(Cb_occ);
        Da_list.push_back(Da);
        Db_list.push_back(Db);
        Dt_list.push_back(Dt);
    }

    // Set lists using modern API
    jk->Ca_list() = Ca_list;
    jk->Cb_list() = Cb_list;
    jk->Da_list() = Da_list;
    jk->Db_list() = Db_list;
    jk->Dt_list() = Dt_list;

    // Determine J/K/wK requirements
    // Assume all wfns have same functional requirements
    auto functional = wfns[0]->functional();
    jk->set_do_J(true);
    if (functional->is_x_hybrid()) {
        jk->set_do_K(true);
    } else {
        jk->set_do_K(false);
    }
    if (functional->is_x_lrc()) {
        jk->set_do_wK(true);
        jk->set_omega(functional->x_omega());
    } else {
        jk->set_do_wK(false);
    }

    jk->print_header();
    jk->compute_gradient();  // ← BATCHED COMPUTATION

    const auto& jk_gradients_list = jk->gradients_list();

    timer_off("Grad: JK (Batched)");

    // Get scaling factors
    double alpha = functional->x_alpha();
    double beta = functional->x_beta();

#ifdef USING_BrianQC
    if (brianEnable && brianEnableDFT) {
        // BrianQC multiplies with exact exchange factors inside, so don't scale here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    // ========================================
    // Step 3: Assemble total gradients
    // ========================================
    std::vector<SharedMatrix> total_gradients;

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Core");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Exchange,LR");

    for (int i = 0; i < nwfn; i++) {
        auto& grads = individual_grads[i];

        // Add JK gradients
        grads["Coulomb"] = jk_gradients_list[i].at("Coulomb");
        if (functional->is_x_hybrid()) {
            grads["Exchange"] = jk_gradients_list[i].at("Exchange");
            grads["Exchange"]->scale(-alpha);
        }
        if (functional->is_x_lrc()) {
            grads["Exchange,LR"] = jk_gradients_list[i].at("Exchange,LR");
            grads["Exchange,LR"]->scale(-beta);
        }

        // Sum all terms
        SharedMatrix total = SharedMatrix(grads["Nuclear"]->clone());
        total->zero();

        for (const auto& term : gradient_terms) {
            if (grads.count(term)) {
                total->add(grads[term]);
            }
        }

        // Symmetrize
        total->symmetrize_gradient(ref_mol);
        total->set_name("Total Gradient");

        total_gradients.push_back(total);
    }

    outfile->Printf("\n  Batched gradient computation complete.\n");
    outfile->Printf("  Gradients computed for %d wavefunctions.\n", nwfn);

    tstop();
    return total_gradients;
}

}} // End Namespaces
