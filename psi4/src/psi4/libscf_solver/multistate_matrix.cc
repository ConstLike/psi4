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
#include <cstring>
#include <cstdlib>

namespace psi {

MultiStateMatrix::MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                                   const Dimension& rowspi, const Dimension& colspi, int symmetry)
    : n_states_(n_states),
      nirrep_(nirrep),
      symmetry_(symmetry),
      rowspi_(rowspi),
      colspi_(colspi),
      name_(name),
      total_elements_(0),
      data_contiguous_(nullptr) {

    if (n_states_ <= 0) {
        throw PSIEXCEPTION("MultiStateMatrix: n_states must be positive");
    }

    state_views_.resize(n_states_);
    allocate();
}

MultiStateMatrix::MultiStateMatrix(const std::string& name, int n_states, int nirrep,
                                   const Dimension& dims, int symmetry)
    : MultiStateMatrix(name, n_states, nirrep, dims, dims, symmetry) {}

MultiStateMatrix::~MultiStateMatrix() {
    deallocate();
}

void MultiStateMatrix::allocate() {
    // Calculate total size needed for all states
    total_elements_ = 0;
    for (int s = 0; s < n_states_; ++s) {
        for (int h = 0; h < nirrep_; ++h) {
            size_t block_size = rowspi_[h] * colspi_[h ^ symmetry_];
            total_elements_ += block_size;
        }
    }

    // Allocate single contiguous aligned block (64-byte alignment)
    if (total_elements_ > 0) {
        constexpr size_t CACHE_LINE_SIZE = 64;
        void* aligned_ptr = nullptr;
        int ret = posix_memalign(&aligned_ptr, CACHE_LINE_SIZE,
                                total_elements_ * sizeof(double));
        if (ret != 0) {
            throw PSIEXCEPTION("MultiStateMatrix::allocate: posix_memalign failed");
        }
        data_contiguous_ = static_cast<double*>(aligned_ptr);

        // Zero initialize
        std::memset(data_contiguous_, 0, total_elements_ * sizeof(double));
    }

    // Create Matrix views pointing into contiguous block
    size_t offset = 0;
    for (int s = 0; s < n_states_; ++s) {
        std::string state_name = name_ + "[" + std::to_string(s) + "]";
        auto mat = std::make_shared<Matrix>(state_name, nirrep_, rowspi_, colspi_, symmetry_);

        // Copy data from contiguous block to Matrix using public API
        for (int h = 0; h < nirrep_; ++h) {
            int nrow = rowspi_[h];
            int ncol = colspi_[h ^ symmetry_];
            if (nrow > 0 && ncol > 0) {
                // Use Matrix::pointer() to get data pointer (public method)
                double** mat_data = mat->pointer(h);
                for (int i = 0; i < nrow; ++i) {
                    std::memcpy(mat_data[i], &data_contiguous_[offset], ncol * sizeof(double));
                    offset += ncol;
                }
            }
        }

        state_views_[s] = mat;
    }
}

void MultiStateMatrix::deallocate() {
    if (data_contiguous_) {
        ::free(data_contiguous_);
        data_contiguous_ = nullptr;
    }
    state_views_.clear();
}

SharedMatrix MultiStateMatrix::get(int state) const {
    if (state < 0 || state >= n_states_) {
        throw PSIEXCEPTION("MultiStateMatrix::get: state index out of bounds");
    }
    return state_views_[state];
}

void MultiStateMatrix::zero_all() {
    if (data_contiguous_) {
        std::memset(data_contiguous_, 0, total_elements_ * sizeof(double));
    }

    // Also zero the Matrix views
    for (int s = 0; s < n_states_; ++s) {
        if (state_views_[s]) {
            state_views_[s]->zero();
        }
    }
}

void MultiStateMatrix::copy_all_from(const MultiStateMatrix& other) {
    if (n_states_ != other.n_states_ || nirrep_ != other.nirrep_ ||
        rowspi_ != other.rowspi_ || colspi_ != other.colspi_) {
        throw PSIEXCEPTION("MultiStateMatrix::copy_all_from: dimension mismatch");
    }

    // Copy contiguous data
    if (data_contiguous_ && other.data_contiguous_) {
        std::memcpy(data_contiguous_, other.data_contiguous_,
                   total_elements_ * sizeof(double));
    }

    // Copy to Matrix views
    for (int s = 0; s < n_states_; ++s) {
        if (state_views_[s] && other.state_views_[s]) {
            state_views_[s]->copy(other.state_views_[s]);
        }
    }
}

void MultiStateMatrix::print() const {
    outfile->Printf("\n  MultiStateMatrix: %s\n", name_.c_str());
    outfile->Printf("  Number of states: %d\n", n_states_);
    outfile->Printf("  Number of irreps: %d\n", nirrep_);
    outfile->Printf("  Total elements: %zu\n", total_elements_);
    outfile->Printf("  Contiguous storage: %s\n\n",
                   data_contiguous_ ? "allocated" : "not allocated");

    for (int s = 0; s < n_states_; ++s) {
        outfile->Printf("  State %d:\n", s);
        if (state_views_[s]) {
            state_views_[s]->print();
        }
    }
}

}  // namespace psi
