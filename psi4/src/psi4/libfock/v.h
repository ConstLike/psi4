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

#ifndef LIBFOCK_DFT_H
#define LIBFOCK_DFT_H
#include "psi4/libmints/typedefs.h"
#include "psi4/pragma.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace psi {
class BasisSet;
class Options;
class DFTGrid;
class PointFunctions;
class SuperFunctional;
class BlockOPoints;

// => BASE CLASS <= //

/**
 * Class VBase
 *
 * Class to compute KS-V matrices and
 * K-matrix-vector products
 **/

class PSI_API VBase {
   protected:
    /// Debug flag
    int debug_;
    /// Print flag
    int print_;
    /// Number of threads
    int num_threads_;
    /// Number of basis functions;
    int nbf_;
    /// Rho threshold for the second derivative;
    double v2_rho_cutoff_;
    /// VV10 interior kernel threshold
    double vv10_rho_cutoff_;
    /// Options object, used to build grid
    Options& options_;
    /// Basis set used in the integration
    std::shared_ptr<BasisSet> primary_;
    /// Desired superfunctional kernel
    std::shared_ptr<SuperFunctional> functional_;
    /// Desired superfunctional kernel
    std::vector<std::shared_ptr<SuperFunctional>> functional_workers_;
    /// Point function computer (densities, gammas, basis values)
    std::vector<std::shared_ptr<PointFunctions>> point_workers_;
    /// Integration grid, built by KSPotential
    std::shared_ptr<DFTGrid> grid_;
    /// Quadrature values obtained during integration
    std::map<std::string, double> quad_values_;
    // Caches collocation grids
    std::unordered_map<size_t, std::map<std::string, SharedMatrix>> cache_map_;
    int cache_map_deriv_;

    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;
    SharedMatrix USO2AO_;

    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;

    // GRAC data
    bool grac_initialized_;

    // VV10 dispersion, return vv10_nlc energy
    void prepare_vv10_cache(DFTGrid& nlgrid, SharedMatrix D,
                            std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                            std::vector<std::shared_ptr<PointFunctions>>& nl_point_workers, int ansatz = 1);
    double vv10_nlc(SharedMatrix D, SharedMatrix ret);
    SharedMatrix vv10_nlc_gradient(SharedMatrix D);

    /// Rebuild the per-thread point workers if a subclass drops them; a no-op by default.
    virtual void ensure_point_workers() {}
    /// Set things up
    void common_init();

   public:
    VBase(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    virtual ~VBase();

    static std::shared_ptr<VBase> build_V(std::shared_ptr<BasisSet> primary,
                                          std::shared_ptr<SuperFunctional> functional, Options& options,
                                          const std::string& type = "RV");

    std::shared_ptr<BasisSet> basis() const { return primary_; }
    std::shared_ptr<SuperFunctional> functional() const { return functional_; }
    std::vector<std::shared_ptr<PointFunctions>> properties() const { return point_workers_; }
    std::shared_ptr<DFTGrid> grid() const { return grid_; }
    /// Adopt a grid built for the same molecule, basis and options instead of building one.
    /// Must be called before initialize(), which only builds a grid when none is set.
    void set_grid(std::shared_ptr<DFTGrid> grid) { grid_ = grid; }
    std::shared_ptr<BlockOPoints> get_block(int block);
    size_t nblocks();
    std::map<std::string, double>& quadrature_values() { return quad_values_; }

    // Creates a collocation cache map based on stride
    void build_collocation_cache(size_t memory);
    void clear_collocation_cache() { cache_map_.clear(); }

    void share_collocation_cache_from(VBase& source);

    // Set the D matrix, get it back if needed
    void set_D(std::vector<SharedMatrix> Dvec);
    const std::vector<SharedMatrix>& Dao() const { return D_AO_; }

    // Set the site of the grac shift
    void set_grac_shift(double value);

    /// Throws by default
    virtual void compute_V(std::vector<SharedMatrix> ret);
    /// Throws by default. Compute the orbital derivative of the KS potential for each spin,
    /// contract against Dx, and putting the result in ret.
    virtual void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    virtual std::vector<SharedMatrix> compute_fock_derivatives();
    virtual SharedMatrix compute_gradient();
    virtual SharedMatrix compute_hessian();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    virtual void initialize();
    virtual void finalize();

    virtual void print_header() const;
};

// => Derived Classes <= //
class SAP : public VBase {
   protected:
   public:
    SAP(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~SAP() override;

    void initialize() override;
    void finalize() override;

    void compute_V(std::vector<SharedMatrix> ret) override;
    void print_header() const override;
};

class RV : public VBase {
   protected:
    /// Build one RKSFunctions per thread over the current grid extents.
    void build_point_workers();
    void ensure_point_workers() override { if (point_workers_.empty()) build_point_workers(); }

   public:
    /// Drop the per-thread point workers; every RV sweep rebuilds them if they are gone.
    void release_point_workers() { point_workers_.clear(); }
    RV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~RV() override;

    void initialize() override;
    void finalize() override;

    // compute_V assuming same orbitals for different spin. Computes V_alpha, not spin-summed V.
    void compute_V(std::vector<SharedMatrix> ret) override;
    /// Compute the orbital derivative of the KS potential, contract against Dx, and
    /// putting the result in ret. ret[i] is Vx where x = Dx[i]. The "true" vector has
    /// 2^-0.5 Dx[i] for each input spin case and returns **half** the α component of the output.
    /// The singlet flag controls whether to assume singlet spin-integration (β components
    /// are the α components) or triplet (β components are -α components)
    void compute_Vx_full(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret, bool singlet);
    /// A convenience function to call compute_Vx_full for singlets.
    /// And no, we can't just make singlet a default argument. Then compute_Vx has different signatures for
    /// different VBase subclasses, so we can't call compute_Vx from VBase, which breaks the hessian code.
    void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) override { compute_Vx_full(Dx, ret, true); };
    std::vector<SharedMatrix> compute_fock_derivatives() override;
    SharedMatrix compute_gradient() override;
    SharedMatrix compute_hessian() override;

    void print_header() const override;
};

class UV : public VBase {
   protected:
    /// Per-thread scratch, retained across calls and rebuilt only when Shape changes.
    /// Buffers are zeroed before use; retention alone does not clear them.
    struct MicrostateScratch {
        /// Every extent the buffers depend on. A mismatch rebuilds the whole set.
        struct Shape {
            int nbf = -1;
            int nmo = -1;
            int n_p = -1;
            int nact = -1;
            int n_micro = -1;
            int n_strings = -1;
            int max_points = -1;
            int max_functions = -1;
            int chunk = -1;
            int ansatz = -1;
            int num_threads = -1;
            int need_diag = -1;
            int acc_w = -1;      ///< second axis of the weighted active-row accumulator
            int project_in_block = -1;  ///< 1 when MO window == active block; TapWide is then one slice wide, not several
            bool operator==(const Shape& o) const {
                return nbf == o.nbf && nmo == o.nmo && n_p == o.n_p && nact == o.nact &&
                       n_micro == o.n_micro && n_strings == o.n_strings && max_points == o.max_points &&
                       max_functions == o.max_functions && chunk == o.chunk && ansatz == o.ansatz &&
                       num_threads == o.num_threads && need_diag == o.need_diag &&
                       acc_w == o.acc_w && project_in_block == o.project_in_block;
            }
        };
        /// One thread's buffers.
        struct Thread {
            SharedMatrix Cloc, chi, chix, chiy, chiz;
            SharedMatrix Sd, TapAcc_a, TapAcc_b, Vloc, MbpT;
            SharedMatrix SDrho, SDgx, SDgy, SDgz, SDgam, SDtau;
            SharedMatrix TapWide, TAwide, M1w, M2S, M3, Vcof;
            std::map<std::string, SharedVector> bin;
            std::shared_ptr<SuperFunctional> bworker;
            std::vector<double> diA, diB, trA, trB, exc, stauA, stauB;
            /// MO-diagonal staging: point-local products, and the per-slice diagonal.
            std::vector<double> Zd, Dd;
        };
        Shape shape;
        std::vector<Thread> threads;
        /// Call-level scratch shared by every thread: weighted-AO active rows per spin,
        /// written under one lock per microstate, and the MO-projection buffer.
        std::vector<double> wa_a_all, wa_b_all, proj;
    };
    MicrostateScratch ms_scratch_;

   public:
    UV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~UV() override;

    void initialize() override;
    void finalize() override;

    void compute_V(std::vector<SharedMatrix> ret) override;
    /// Orbital derivative of the KS potential, contracted against Dx.
    /// ret[i] is Vx where x = Dx[i].
    /// ret[2n], ret[2n+1] are alpha and beta Vx where x concatenates Dx[2n] (alpha) and Dx[2n+1] (beta).
    void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) override;
    std::vector<SharedMatrix> compute_fock_derivatives() override;
    SharedMatrix compute_gradient() override;
    SharedMatrix compute_hessian() override;

    /// Grid-direct microstate XC. For K microstates sharing the MO basis Ca_full
    /// (all nmo columns), computes WITHOUT materializing a per-microstate AO V_xc:
    ///   - E_xc[L], tr_a[L]=tr(rho^L_a V_xc^L_a), tr_b[L]   (self-trace, grid form);
    ///   - arows_{a,b}[L] [n_active x mo_width]: active MO rows of C^T V_xc^L C over the
    ///     MO column window [mo_lo, mo_lo + mo_width), mo_width taken from the arows;
    ///   - diag_{a,b}[L]  [nmo]: diagonal of C^T V_xc^L C, only when need_diag; left
    ///     empty otherwise;
    ///   - V_acc_{a,b} [nbf x nbf]: ONE C_L-weighted AO back-projection Sum_L C_L V_xc^L_{a,b},
    ///     skipped entirely when the pair is null;
    ///   - quad_values_["RHO_A"/"RHO_B"] = quadrature integral of the C_L-weighted
    ///     total density (grid-electron diagnostic, closed-shell RHO convention).
    /// occ over the first n_p = Ncore+n_active columns (active contiguous: active_mo[i]=Ncore+i).
    void compute_V_microstates_mo(SharedMatrix Ca_full,
                                  int Ncore, int n_active,
                                  const std::vector<int>& active_mo,
                                  const std::vector<std::vector<double>>& occ_alpha,
                                  const std::vector<std::vector<double>>& occ_beta,
                                  const std::vector<double>& C_L,
                                  std::vector<SharedMatrix>& arows_a,
                                  std::vector<SharedMatrix>& arows_b,
                                  std::vector<std::vector<double>>& diag_a,
                                  std::vector<std::vector<double>>& diag_b,
                                  std::vector<double>& tr_a,
                                  std::vector<double>& tr_b,
                                  std::vector<double>& E_xc,
                                  SharedMatrix V_acc_a,
                                  SharedMatrix V_acc_b,
                                  bool need_diag,
                                  int mo_lo = 0);

    void print_header() const override;
};
}
#endif
