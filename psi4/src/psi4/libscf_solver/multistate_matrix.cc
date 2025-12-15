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

#include "multistate_matrix.h"
#include "psi4/libmints/matrix.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

MultiStateMatrix::MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                                   const Dimension& rowspi, const Dimension& colspi, int symmetry)
    : nirrep_(nirrep),
      symmetry_(symmetry),
      rowspi_(rowspi),
      colspi_(colspi),
      name_(name) {

    if (n_states <= 0) {
        throw PSIEXCEPTION("MultiStateMatrix: n_states must be positive");
    }

    // Reserve space and create Matrix objects
    states_.reserve(n_states);
    for (int s = 0; s < n_states; ++s) {
        std::string state_name = name_ + "[" + std::to_string(s) + "]";
        auto mat = std::make_shared<Matrix>(state_name, nirrep_, rowspi_, colspi_, symmetry_);
        mat->zero();
        states_.push_back(std::move(mat));
    }
}

MultiStateMatrix::MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                                   const Dimension& dims, int symmetry)
    : MultiStateMatrix(name, n_states, nirrep, dims, dims, symmetry) {}

SharedMatrix MultiStateMatrix::get(int state) const {
    if (state < 0 || state >= static_cast<int>(states_.size())) {
        throw PSIEXCEPTION("MultiStateMatrix::get: state index out of bounds");
    }
    return states_[state];
}

void MultiStateMatrix::zero_all() {
    for (auto& mat : states_) {
        mat->zero();
    }
}

void MultiStateMatrix::copy_all_from(const MultiStateMatrix& other) {
    if (static_cast<int>(states_.size()) != other.n_states() ||
        nirrep_ != other.nirrep_ ||
        rowspi_ != other.rowspi_ ||
        colspi_ != other.colspi_) {
        throw PSIEXCEPTION("MultiStateMatrix::copy_all_from: dimension mismatch");
    }

    for (size_t s = 0; s < states_.size(); ++s) {
        states_[s]->copy(other.states_[s]);
    }
}

void MultiStateMatrix::print() const {
    outfile->Printf("\n  MultiStateMatrix: %s\n", name_.c_str());
    outfile->Printf("  Number of states: %zu\n", states_.size());
    outfile->Printf("  Number of irreps: %d\n", nirrep_);
    outfile->Printf("  Symmetry: %d\n\n", symmetry_);

    for (size_t s = 0; s < states_.size(); ++s) {
        outfile->Printf("  State %zu:\n", s);
        if (states_[s]) {
            states_[s]->print();
        }
    }
}

}  // namespace psi
