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

#ifndef JK_GRAD_H
#define JK_GRAD_H

#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"
#include <map>
#include <vector>
#include <string>

namespace psi {

class BasisSet;
class PSIO;
class TwoBodyAOInt;
class MintsHelper;

namespace scfgrad {

class JKGrad {

protected:
    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    size_t memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, Process::environment.get_n_threads() otherwise)
    int omp_num_threads_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;
    /// Maximum derivative level
    int deriv_;

    std::shared_ptr<BasisSet> primary_;

    // === Phase B: Multi-wfn Architecture (lists-only, like JK) === //
    // Design Philosophy:
    // - JKGrad works with lists ONLY (following JK::C_left_ pattern)
    // - Single-wfn mode = list of size 1
    // - No batch_mode_ flag needed - always iterate over list
    // - Legacy set_Da() API provided as wrappers (clear + push_back)

    /// Alpha orbital matrices for each wavefunction
    std::vector<SharedMatrix> Ca_list_;
    /// Beta orbital matrices for each wavefunction
    std::vector<SharedMatrix> Cb_list_;
    /// Alpha density matrices for each wavefunction
    std::vector<SharedMatrix> Da_list_;
    /// Beta density matrices for each wavefunction
    std::vector<SharedMatrix> Db_list_;
    /// Total density matrices (Da + Db) for each wavefunction
    std::vector<SharedMatrix> Dt_list_;

    bool do_J_;
    bool do_K_;
    bool do_wK_;

    double omega_;

    /// Gradient results for each wavefunction: gradients_list_[wfn_idx]["J"], ["K"], ["Total"]
    /// This is the primary output - contains gradients for ALL wavefunctions
    std::vector<std::map<std::string, SharedMatrix>> gradients_list_;

    /// Legacy single-wfn output (for backward compatibility)
    /// This is populated by gradients() getter as a copy of gradients_list_[0]
    std::map<std::string, SharedMatrix> gradients_;

    std::map<std::string, SharedMatrix> hessians_;

    void common_init();

public:
    JKGrad(int deriv, std::shared_ptr<BasisSet> primary);
    virtual ~JKGrad();

    /**
    * Static instance constructor, used to get prebuilt DFJK/DirectJK objects
    * using knobs in options.
    * @param options Options reference, with preset parameters
    * @return abstract JK object, tuned in with preset options
    */
    static std::shared_ptr<JKGrad> build_JKGrad(int deriv, std::shared_ptr<MintsHelper> mints);

    // === Legacy Single-Wfn API (wrappers for backward compatibility) === //
    // Design: These wrappers clear the list and add single matrix (list of size 1)
    // This allows existing code like jk_grad->set_Da(Da) to continue working

    /// Legacy: Set single Ca (internally stored as list of size 1)
    void set_Ca(SharedMatrix Ca) {
        Ca_list_.clear();
        Ca_list_.push_back(Ca);
    }

    /// Legacy: Set single Cb
    void set_Cb(SharedMatrix Cb) {
        Cb_list_.clear();
        Cb_list_.push_back(Cb);
    }

    /// Legacy: Set single Da
    void set_Da(SharedMatrix Da) {
        Da_list_.clear();
        Da_list_.push_back(Da);
    }

    /// Legacy: Set single Db
    void set_Db(SharedMatrix Db) {
        Db_list_.clear();
        Db_list_.push_back(Db);
    }

    /// Legacy: Set single Dt
    void set_Dt(SharedMatrix Dt) {
        Dt_list_.clear();
        Dt_list_.push_back(Dt);
    }

    // === Modern List-Based API (like JK::C_left()) === //

    /// Direct access to Ca list (for multi-wfn usage)
    std::vector<SharedMatrix>& Ca_list() { return Ca_list_; }

    /// Direct access to Cb list
    std::vector<SharedMatrix>& Cb_list() { return Cb_list_; }

    /// Direct access to Da list
    std::vector<SharedMatrix>& Da_list() { return Da_list_; }

    /// Direct access to Db list
    std::vector<SharedMatrix>& Db_list() { return Db_list_; }

    /// Direct access to Dt list
    std::vector<SharedMatrix>& Dt_list() { return Dt_list_; }

    /**
     * Cutoff for individual contributions to the J/K matrices
     * Eventually we hope to use Schwarz/MBIE/Density cutoffs,
     * for now just Schwarz
     * @param cutoff ceiling of magnitude of elements to be
     *        ignored if possible
     */
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    /**
     * Maximum memory to use, in doubles (for tensor-based methods,
     * integral generation objects typically ignore this)
     * @param memory maximum number of doubles to allocate
     */
    void set_memory(size_t memory) { memory_ = memory; }
    /**
     * Maximum number of OpenMP threads to use. It may be necessary
     * to clamp this to some value smaller than the total number of
     * cores for machines with a high core-to-memory ratio to avoid
     * running out of memory due to integral generation objects
     * @param omp_nthread Maximum number of threads to use in
     *        integral generation objects (BLAS/LAPACK can still
     *        run with their original maximum number)
     */
    void set_omp_num_threads(int omp_nthread) { omp_num_threads_ = omp_nthread; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
    /**
    * Set to do J tasks
    * @param do_J do J matrices or not,
    *        defaults to true
    */
    void set_do_J(bool do_J) { do_J_ = do_J; }
    /**
    * Set to do K tasks
    * @param do_K do K matrices or not,
    *        defaults to true
    */
    void set_do_K(bool do_K) { do_K_ = do_K; }
    /**
    * Set to do wK tasks
    * @param do_wK do wK matrices or not,
    *        defaults to false
    */
    void set_do_wK(bool do_wK) { do_wK_ = do_wK; }
    /**
    * Set the omega value for wK
    * @param omega range-separation parameter
    */
    void set_omega(double omega) { omega_ = omega; }

    // === Gradient Results Accessors === //

    /// Modern: Get gradients for all wavefunctions (primary interface)
    /// Returns: gradients_list_[wfn_idx]["J"], ["K"], ["Total"]
    const std::vector<std::map<std::string, SharedMatrix>>& gradients_list() const {
        return gradients_list_;
    }

    /// Legacy: Get gradient for first wavefunction only (backward compatibility)
    /// This copies gradients_list_[0] into gradients_ and returns it
    /// For multi-wfn usage, prefer gradients_list() instead
    std::map<std::string, SharedMatrix>& gradients() {
        if (gradients_list_.empty()) {
            throw PSIEXCEPTION("JKGrad::gradients(): No gradients computed. Call compute_gradient() first.");
        }
        gradients_ = gradients_list_[0];  // Copy first wfn gradients
        return gradients_;
    }

    std::map<std::string, SharedMatrix>& hessians() { return hessians_; }

    virtual void compute_gradient() = 0;
    virtual void compute_hessian() = 0;

    virtual void print_header() const = 0;
};

class DFJKGrad : public JKGrad {

protected:
    std::shared_ptr<BasisSet> auxiliary_;
    std::shared_ptr<MintsHelper> mints_;

    std::shared_ptr<PSIO> psio_;

    /// Number of threads for DF integrals
    int df_ints_num_threads_;
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;

    void common_init();

    void build_Amn_terms();
    void build_Amn_lr_terms();
    void build_AB_inv_terms();
    void build_UV_terms();
    void build_AB_x_terms();
    void build_Amn_x_terms();

    /// File number for Alpha (Q|mn) tensor (single-wfn legacy)
    size_t unit_a_;
    /// File number for Beta (Q|mn) tensor (single-wfn legacy)
    size_t unit_b_;
    /// File number for J tensors (single-wfn legacy)
    size_t unit_c_;

    /// Multi-wfn PSIO units: [wfn_idx] -> {unit_a, unit_b, unit_c}
    /// unit_a: Alpha (Q|ij) for this wfn
    /// unit_b: Beta (Q|ij) for this wfn
    /// unit_c: c, W, V vectors/matrices for this wfn
    std::vector<std::tuple<size_t, size_t, size_t>> wfn_units_;

    /// Shared metric inverse unit (same for all wfn)
    size_t unit_metric_inv_;

public:
    DFJKGrad(int deriv, std::shared_ptr<MintsHelper> mints);
    ~DFJKGrad() override;

    void compute_gradient() override;
    void compute_hessian() override;

    void print_header() const override;

    /**
     * Minimum relative eigenvalue to retain in fitting inverse
     * All eigenvectors with \epsilon_i < condition * \epsilon_max
     * will be discarded
     * @param condition, minimum relative eigenvalue allowed,
     *        defaults to 1.0E-12
     */
    void set_condition(double condition) { condition_ = condition; }
    /**
     * Which file number should the Alpha (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit_a(size_t unit) { unit_a_ = unit; }
    /**
     * Which file number should the Beta (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit_b(size_t unit) { unit_b_ = unit; }
    /**
     * Which file number should the J tensors go in
     * @param unit Unit number
     */
    void set_unit_c(size_t unit) { unit_c_ = unit; }

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }
};

class DirectJKGrad : public JKGrad {

protected:
    // Number of threads to use
    int ints_num_threads_;

    void common_init();

    std::map<std::string, std::shared_ptr<Matrix> > compute1(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints);
    std::map<std::string, std::shared_ptr<Matrix> > compute2(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints);
public:
    DirectJKGrad(int deriv, std::shared_ptr<BasisSet> primary);
    ~DirectJKGrad() override;

    void compute_gradient() override;
    void compute_hessian() override;

    void print_header() const override;

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_ints_num_threads(int val) { ints_num_threads_ = val; }


};

}} // Namespaces
#endif
