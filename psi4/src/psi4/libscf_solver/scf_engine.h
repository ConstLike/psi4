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

#ifndef SCF_ENGINE_H
#define SCF_ENGINE_H

#include <vector>
#include <memory>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace scf {

// Forward declaration
class HF;

/**
 * @brief Base class for SCF iteration engines
 *
 * SCFEngine coordinates the SCF iteration process by calling virtual methods
 * on a theory implementation (HF subclass). This is Layer 1 of the multi-cycle
 * architecture, providing a foundation for more advanced iteration schemes.
 *
 * Design goals:
 * - Non-breaking: opt-in, doesn't modify existing RHF/UHF code
 * - Reusable: provides common iteration logic for any HF theory
 * - Extensible: can be subclassed for multi-state, multi-cycle, etc.
 *
 * Usage:
 *   auto rhf = std::make_shared<RHF>(...);
 *   SCFEngine engine(rhf.get(), options);
 *   int iterations = engine.iterate();
 */
class SCFEngine {
   protected:
    // Reference to theory implementation (not owned)
    HF* theory_;

    // Iteration control parameters
    int max_iterations_;
    double e_convergence_;
    double d_convergence_;

    // State tracking
    int iteration_;
    bool converged_;
    std::vector<double> energy_history_;

    // Density RMS for convergence
    double rms_density_;

   public:
    /**
     * @brief Construct an SCF engine
     * @param theory Pointer to HF theory (RHF, UHF, etc.)
     * @param options Options object containing convergence thresholds
     */
    SCFEngine(HF* theory, Options& options);

    virtual ~SCFEngine() = default;

    /**
     * @brief Main iteration loop
     * @return Number of iterations performed
     *
     * Iterates until convergence or max iterations reached.
     * Calls iterate_step() and check_convergence() each iteration.
     */
    virtual int iterate();

   protected:
    /**
     * @brief Perform a single SCF iteration
     *
     * Default implementation calls theory virtual methods:
     *   1. form_G() - Build two-electron part (J, K, XC)
     *   2. form_F() - Assemble Fock matrix
     *   3. form_C() - Diagonalize Fock
     *   4. form_D() - Build density from orbitals
     *
     * Can be overridden for alternative iteration schemes.
     */
    virtual void iterate_step();

    /**
     * @brief Check convergence criteria
     *
     * Checks:
     *   - Energy change < e_convergence
     *   - Density RMS change < d_convergence
     *
     * Sets converged_ flag if both criteria met.
     */
    virtual void check_convergence();

    /**
     * @brief Compute density RMS change
     * @return RMS change in density matrix
     *
     * Computes ||D_new - D_old||_RMS for convergence checking.
     */
    virtual double compute_density_rms();

   public:
    // Accessors
    bool is_converged() const { return converged_; }
    int get_iteration() const { return iteration_; }
    double get_energy() const { return energy_history_.empty() ? 0.0 : energy_history_.back(); }
    const std::vector<double>& get_energy_history() const { return energy_history_; }
    double get_rms_density() const { return rms_density_; }
};

}  // namespace scf
}  // namespace psi

#endif  // SCF_ENGINE_H
