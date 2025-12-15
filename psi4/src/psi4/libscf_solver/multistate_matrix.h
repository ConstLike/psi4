/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * @END LICENSE
 */

#ifndef MULTISTATE_MATRIX_H
#define MULTISTATE_MATRIX_H

#include <vector>
#include <string>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/dimension.h"

namespace psi {

/**
 * MultiStateMatrix: Zero-overhead wrapper for N Matrix objects
 *
 * Provides unified interface for managing multiple Matrix objects,
 * typically used for multi-state SCF calculations (RHF/UHF/ROHF).
 *
 * Use cases:
 * - RHF: n_states=1 (single density/Fock)
 * - UHF: n_states=2 (Da/Db matrices)
 * - ROHF: n_states=2 (Da/Db matrices)
 * - Multi-state: n_states=N (multiple ensemble states)
 *
 * Design:
 * - Simple wrapper around std::vector<SharedMatrix>
 * - No memory duplication - each Matrix manages its own storage
 * - Zero overhead - direct access to underlying Matrix objects
 * - Compatible with JK::compute() batch operations via matrices()
 */
class PSI_API MultiStateMatrix {
   protected:
    /// State matrices - each Matrix manages its own memory
    std::vector<SharedMatrix> states_;

    /// Number of irreps (point group symmetry)
    int nirrep_;

    /// Symmetry (0 for symmetric, 1 for antisymmetric)
    int symmetry_;

    /// Rows per irrep
    Dimension rowspi_;

    /// Columns per irrep
    Dimension colspi_;

    /// Name for this multi-state container
    std::string name_;

   public:
    /**
     * Constructor
     * @param name Name for this container
     * @param n_states Number of states (1, 2, N)
     * @param nirrep Number of irreps
     * @param rowspi Rows per irrep (Dimension)
     * @param colspi Columns per irrep (Dimension)
     * @param symmetry Symmetry (0 or 1)
     */
    MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                     const Dimension& rowspi, const Dimension& colspi, int symmetry = 0);

    /**
     * Constructor (square matrices, rowspi == colspi)
     * @param name Name for this container
     * @param n_states Number of states
     * @param nirrep Number of irreps
     * @param dims Dimensions per irrep
     * @param symmetry Symmetry (0 or 1)
     */
    MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                     const Dimension& dims, int symmetry = 0);

    /// Default destructor
    ~MultiStateMatrix() = default;

    // Allow move operations
    MultiStateMatrix(MultiStateMatrix&&) = default;
    MultiStateMatrix& operator=(MultiStateMatrix&&) = default;

    // Disable copy (use copy_all_from for explicit copying)
    MultiStateMatrix(const MultiStateMatrix&) = delete;
    MultiStateMatrix& operator=(const MultiStateMatrix&) = delete;

    /// Get Matrix for a specific state
    SharedMatrix get(int state) const;

    /// Array-style access to state matrices
    SharedMatrix operator[](int state) { return states_[state]; }
    const SharedMatrix operator[](int state) const { return states_[state]; }

    /// Number of states
    int n_states() const { return static_cast<int>(states_.size()); }

    /// Number of irreps
    int nirrep() const { return nirrep_; }

    /// Symmetry
    int symmetry() const { return symmetry_; }

    /// Rows per irrep
    const Dimension& rowspi() const { return rowspi_; }

    /// Columns per irrep
    const Dimension& colspi() const { return colspi_; }

    /// Name
    const std::string& name() const { return name_; }

    /// Direct access to underlying vector (for JK batch operations)
    std::vector<SharedMatrix>& matrices() { return states_; }
    const std::vector<SharedMatrix>& matrices() const { return states_; }

    /// Zero all states
    void zero_all();

    /// Copy all states from another MultiStateMatrix
    void copy_all_from(const MultiStateMatrix& other);

    /// Print information
    void print() const;
};

}  // namespace psi

#endif  // MULTISTATE_MATRIX_H
