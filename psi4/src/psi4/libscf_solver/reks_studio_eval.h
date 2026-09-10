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

#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "reks_arch.h"
#include "psi4/libmints/typedefs.h"   // SharedMatrix

namespace psi {
namespace reks {

class ActiveEriTiles;

namespace studio {

template <typename FonsT>
inline double fon_slot_value(const FonSlot& s, const FonsT& fons) {
    const auto& cell = fons.layers[s.gen][s.g];
    return s.use_q ? cell.q : cell.p;
}

// Value of one FON factor at this snapshot; sqrt-group accumulates its radicand under
// one sqrt (rad_n==1 gives sqrt(p) exactly).
template <typename FonsT>
inline double fon_factor_value(const FonFactor& f, const FonSlot* sp, const int* slot_idx,
                               const FonsT& fons, double (*f_interp)(double)) {
    switch (f.kind) {
        case kFonLinear:
            return fon_slot_value(sp[slot_idx[f.slot]], fons);
        case kFonSqrtGroup: {
            double r = 1.0;
            for (int j = 0; j < f.rad_n; ++j)
                r *= fon_slot_value(sp[slot_idx[f.rad_off + j]], fons);
            return std::sqrt(r);
        }
        default: {                                   // kFonFinterp
            const FonSlot& s = sp[slot_idx[f.slot]];
            const auto& c = fons.layers[s.gen][s.g];
            return f_interp(c.p * c.q);
        }
    }
}

/// True when a slot addresses a geminal this snapshot carries.
template <typename FonsT>
inline bool fon_slot_in_snapshot(const FonSlot& s, const FonsT& fons) {
    return s.gen >= 0 && s.gen < kMaxGen && s.g >= 0 &&
           s.g < static_cast<int>(fons.layers[s.gen].size());
}

/// Every pool factor's value at this snapshot, indexed by pool id; a factor for another
/// spin sector is NaN.
template <typename FonsT>
inline void build_fon_cache(const FonFactor* fp, int n_fon_pool, const FonSlot* sp,
                            const int* slot_idx, const FonsT& fons,
                            double (*f_interp)(double), std::vector<double>& out) {
    out.assign(static_cast<size_t>(n_fon_pool < 0 ? 0 : n_fon_pool),
               std::numeric_limits<double>::quiet_NaN());
    for (int k = 0; k < n_fon_pool; ++k) {
        const FonFactor& f = fp[k];
        bool ok = true;
        if (f.kind == kFonSqrtGroup) {
            for (int j = 0; j < f.rad_n && ok; ++j)
                ok = fon_slot_in_snapshot(sp[slot_idx[f.rad_off + j]], fons);
        } else {
            ok = fon_slot_in_snapshot(sp[slot_idx[f.slot]], fons);
        }
        if (ok) out[k] = fon_factor_value(f, sp, slot_idx, fons, f_interp);
    }
}

// Product of n FON factors indexed by fon_idx[off, off+n).
template <typename FonsT>
inline double fon_product(const FonFactor* fp, const FonSlot* sp,
                          const int* fon_idx, const int* slot_idx,
                          int64_t off, int n,
                          const FonsT& fons, double (*f_interp)(double),
                          const double* fac = nullptr) {
    double v = 1.0;
    if (fac != nullptr) {
        for (int i = 0; i < n; ++i) v *= fac[fon_idx[off + i]];
        return v;
    }
    for (int i = 0; i < n; ++i)
        v *= fon_factor_value(fp[fon_idx[off + i]], sp, slot_idx, fons, f_interp);
    return v;
}

// fac is an optional build_fon_cache result; null evaluates each factor directly.
struct EvalCtx {
    double (*f_interp)(double);
    const double* fac;
};

template <typename FonsT>
inline double eval_term(const Term& tm, const FonFactor* fp, const FonSlot* sp,
                        const int* fon_idx, const int* slot_idx,
                        const FonsT& fons, const EvalCtx& ctx) {
    return tm.coeff * fon_product(fp, sp, fon_idx, slot_idx, tm.fon_off, tm.fon_n, fons,
                                  ctx.f_interp, ctx.fac);
}

// Cells of one RDM block, in cell-pool order; out must hold 2*cell_n entries
// (a p != q cell also emits its mirror at q*N + p). Returns the count written.
template <typename FonsT>
inline int eval_rdm_cells(const BlockPools& bp, int cell_off, int cell_n,
                          const FonsT& fons, const EvalCtx& ctx, int N, RdmCell* out) {
    int n = 0;
    for (int c = 0; c < cell_n; ++c) {
        const Cell cell = load_cell(bp, cell_off + c);
        double acc = 0.0;
        for (int t = 0; t < cell.term_n; ++t)
            acc += eval_term(load_term(bp, cell.term_off + t), bp.fon_pool, bp.fon_slot_pool,
                             bp.fon_idx, bp.fon_slot_idx, fons, ctx);
        out[n++] = RdmCell{cell.p * N + cell.q, acc};
        if (cell.p != cell.q) out[n++] = RdmCell{cell.q * N + cell.p, acc};
    }
    return n;
}

template <typename FonsT>
inline void eval_weight(const BlockPools& bp, int cell_off, int cell_n, const FonsT& fons,
                        const EvalCtx& ctx, std::vector<double>& C_L) {
    for (int c = 0; c < cell_n; ++c) {
        const Cell cell = load_cell(bp, cell_off + c);
        double acc = 0.0;
        for (int t = 0; t < cell.term_n; ++t)
            acc += eval_term(load_term(bp, cell.term_off + t), bp.fon_pool, bp.fon_slot_pool,
                             bp.fon_idx, bp.fon_slot_idx, fons, ctx);
        C_L[cell.p] = acc;               // cell.p == microstate L
    }
}

// FON sign convention: x = p-occupation of geminal (gen,g), q = 2 - x;
// dp/dx = +1, dq/dx = -1, d2/dx2 = 0 for linear slots. Leibniz product
// rule over Term factors below.
struct DerivCtx {
    double (*f_interp)(double);
    double (*df_interp)(double);
    double (*d2f_interp)(double);
};

struct FacJet { double v, d1, d2; };

// d/dx of one linear slot w.r.t. the FON DOF of geminal (gt,ggt); d2/dx2 == 0.
inline double fon_slot_lin_deriv(const FonSlot& s, int gt, int ggt) {
    return (s.gen == gt && s.g == ggt) ? (s.use_q ? -1.0 : 1.0) : 0.0;
}

// Value and 1st/2nd derivative of one factor w.r.t. the FON DOF of geminal (gt,ggt).
template <typename FonsT>
inline FacJet fon_factor_jet(const FonFactor& f, const FonSlot* sp, const int* slot_idx,
                             const FonsT& fons, int gt, int ggt, const DerivCtx& ctx) {
    switch (f.kind) {
        case kFonLinear: {
            const FonSlot& s = sp[slot_idx[f.slot]];
            double v  = fon_slot_value(s, fons);
            double d1 = fon_slot_lin_deriv(s, gt, ggt);
            return { v, d1, 0.0 };
        }
        case kFonSqrtGroup: {                       // grouped sqrt: one sqrt of the group; chain rule on R
            double R = 1.0, dR = 0.0, d2R = 0.0;    // R, dR/dx, d2R/dx2 by Horner fold (sub-factors linear)
            for (int j = 0; j < f.rad_n; ++j) {
                const FonSlot& s = sp[slot_idx[f.rad_off + j]];
                double vs = fon_slot_value(s, fons);
                double ds = fon_slot_lin_deriv(s, gt, ggt);
                d2R = d2R * vs + 2.0 * dR * ds;     // (old dR, old R) before they are updated
                dR  = dR  * vs + R * ds;
                R   = R   * vs;
            }
            if (R <= 0.0) return { 0.0, 0.0, 0.0 };
            double sq = std::sqrt(R);
            double d1 = 0.5 * dR / sq;
            double d2 = 0.5 * d2R / sq - 0.25 * dR * dR / (R * sq);
            return { sq, d1, d2 };
        }
        default: {                                   // kFonFinterp: f_interp(p*q)
            const FonSlot& s = sp[slot_idx[f.slot]];
            const auto& c = fons.layers[s.gen][s.g];
            double prod = c.p * c.q;
            double v = ctx.f_interp(prod);
            if (s.gen != gt || s.g != ggt) return { v, 0.0, 0.0 };
            double dx = c.q - c.p;                   // d(p*q)/dx at q = 2 - p
            double df = ctx.df_interp(prod);
            // d2(p*q)/dx2 = -2
            return { v, df * dx, ctx.d2f_interp(prod) * dx * dx - 2.0 * df };
        }
    }
}

struct MixedFacJet { double v, di, dj, dij; };

// Value, both 1st derivatives and the mixed 2nd derivative of one factor w.r.t.
// the FON DOFs of geminals (gt,gi) and (gt,gj), gi != gj. dij is nonzero only for
// kFonSqrtGroup whose radicand spans both geminals.
template <typename FonsT>
inline MixedFacJet fon_factor_mixed_jet(const FonFactor& f, const FonSlot* sp, const int* slot_idx,
                                        const FonsT& fons, int gt, int gi, int gj,
                                        const DerivCtx& ctx) {
    switch (f.kind) {
        case kFonLinear: {
            const FonSlot& s = sp[slot_idx[f.slot]];
            return { fon_slot_value(s, fons),
                     fon_slot_lin_deriv(s, gt, gi),
                     fon_slot_lin_deriv(s, gt, gj), 0.0 };
        }
        case kFonSqrtGroup: {                       // one fold carries R and its three derivatives
            double R = 1.0, dRi = 0.0, dRj = 0.0, d2Rij = 0.0;
            for (int k = 0; k < f.rad_n; ++k) {
                const FonSlot& s = sp[slot_idx[f.rad_off + k]];
                double vs = fon_slot_value(s, fons);
                double di = fon_slot_lin_deriv(s, gt, gi);
                double dj = fon_slot_lin_deriv(s, gt, gj);
                d2Rij = d2Rij * vs + dRi * dj + dRj * di;    // old dRi,dRj,R before update
                dRi   = dRi   * vs + R * di;
                dRj   = dRj   * vs + R * dj;
                R     = R     * vs;
            }
            if (R <= 0.0) return { 0.0, 0.0, 0.0, 0.0 };
            double sq = std::sqrt(R);
            return { sq, 0.5 * dRi / sq, 0.5 * dRj / sq,
                     0.5 * d2Rij / sq - 0.25 * dRi * dRj / (R * sq) };
        }
        default: {                                   // kFonFinterp: f_interp(p*q)
            const FonSlot& s = sp[slot_idx[f.slot]];
            const auto& c = fons.layers[s.gen][s.g];
            double prod = c.p * c.q;
            double v = ctx.f_interp(prod);
            const bool hit_i = (s.gen == gt && s.g == gi);
            const bool hit_j = (s.gen == gt && s.g == gj);
            if (!hit_i && !hit_j) return { v, 0.0, 0.0, 0.0 };
            double d = ctx.df_interp(prod) * (c.q - c.p);
            return { v, hit_i ? d : 0.0, hit_j ? d : 0.0, 0.0 };
        }
    }
}

// dC_L/dx and d2C_L/dx2 over a weight row's cells; x = FON DOF of geminal (gen,g).
template <typename FonsT>
inline void eval_weight_deriv(const BlockPools& bp, int cell_off, int cell_n,
                              const FonsT& fons, const DerivCtx& ctx, int gen, int g,
                              std::vector<double>& dC, std::vector<double>& d2C) {
    for (int c = 0; c < cell_n; ++c) {
        const Cell cell = load_cell(bp, cell_off + c);
        double dacc = 0.0, d2acc = 0.0;
        for (int t = 0; t < cell.term_n; ++t) {
            const Term tm = load_term(bp, cell.term_off + t);
            double P = 1.0, dP = 0.0, d2P = 0.0;     // product, d/dx, d2/dx2 by Horner fold
            for (int k = 0; k < tm.fon_n; ++k) {
                FacJet j = fon_factor_jet(bp.fon_pool[bp.fon_idx[tm.fon_off + k]], bp.fon_slot_pool, bp.fon_slot_idx, fons, gen, g, ctx);
                d2P = d2P * j.v + 2.0 * dP * j.d1 + P * j.d2;
                dP  = dP  * j.v + P * j.d1;
                P   = P   * j.v;
            }
            dacc  += tm.coeff * dP;
            d2acc += tm.coeff * d2P;
        }
        dC[cell.p]  = dacc;                          // cell.p == microstate L
        d2C[cell.p] = d2acc;
    }
}

// d2C_L/(dx_i dx_j) over a weight row's cells; x_i,x_j = FON DOFs of geminals (gen,g_i),(gen,g_j).
template <typename FonsT>
inline void eval_weight_mixed_deriv(const BlockPools& bp, int cell_off, int cell_n,
                                    const FonsT& fons, const DerivCtx& ctx,
                                    int gen, int g_i, int g_j,
                                    std::vector<double>& d2C_mixed) {
    for (int c = 0; c < cell_n; ++c) {
        const Cell cell = load_cell(bp, cell_off + c);
        double acc = 0.0;
        for (int t = 0; t < cell.term_n; ++t) {
            const Term tm = load_term(bp, cell.term_off + t);
            double P = 1.0, dPi = 0.0, dPj = 0.0, d2Pij = 0.0;   // product, d/dxi, d/dxj, d2/dxidxj by Horner fold
            for (int k = 0; k < tm.fon_n; ++k) {
                const FonFactor& f = bp.fon_pool[bp.fon_idx[tm.fon_off + k]];
                MixedFacJet j = fon_factor_mixed_jet(f, bp.fon_slot_pool, bp.fon_slot_idx,
                                                     fons, gen, g_i, g_j, ctx);
                d2Pij = d2Pij * j.v + dPi * j.dj + dPj * j.di + P * j.dij;   // old dPi,dPj,P before update
                dPi   = dPi   * j.v + P * j.di;
                dPj   = dPj   * j.v + P * j.dj;
                P     = P     * j.v;
            }
            acc += tm.coeff * d2Pij;
        }
        d2C_mixed[cell.p] = acc;
    }
}

// Off-diagonal SI coupling evaluators (H_ij, S_ij).
namespace coupling {

// H_ij = <K_i|H|K_j> off-diagonal SI-Hamiltonian element. `win` carries the pair's
// four value channels (e/fock/lagr/eri). A null fon_cache evaluates each FON factor
// directly. The per-microstate Fock enters only as the active block of its spin
// difference, row fock_aa_row[L] of fock_aa.
double eval_pair_H(const PairValues& win,
                   const FonFactor* fon_pool,
                   const int* fon_idx,
                   const FonSlot* slot_pool,
                   const int* slot_idx,
                   const Coeff* coeff_pool,
                   const ValRow* e_val_pool,
                   const FockIdxPart* fock_idxpart,
                   const ValRow* fock_valpart,
                   const int* fock_pack,
                   const EriIdxPart* eri_idxpart,
                   const ValRow* eri_valpart,
                   const int* eri_pack,
                   const LagrTerm* lagr_row_pool,
                   const FONSnapshot& fons,
                   const std::vector<double>& E_L,
                   const std::vector<int>& fock_aa_row,
                   const std::vector<double>& fock_aa,
                   int n_active,
                   const reks::ActiveEriTiles& active_eri,
                   const std::vector<double>& lagrangians,
                   const double* fon_cache = nullptr);

// S_ij = <K_i|K_j> pure-FON SI overlap over the pair's s channel. A null fon_cache
// evaluates each FON factor directly.
double eval_pair_S(const IntSpan& s_idx, const FonFactor* fon_pool,
                   const int* fon_idx, const FonSlot* slot_pool, const int* slot_idx,
                   const Coeff* coeff_pool, const ValRow* s_row_pool,
                   const FONSnapshot& fons, const double* fon_cache = nullptr);

}  // namespace coupling

}  // namespace studio
}  // namespace reks
}  // namespace psi
