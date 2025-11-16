# CRITICAL FINDING - NO PHYSICAL COUPLING!

## JK Source Code Analysis

**File:** `psi4/src/psi4/libfock/jk.cc`

**Key Functions:**
- `compute_D()` (lines 310-340): Builds density D[N] from C_left[N] × C_right[N]^T
- `AO2USO()` (lines 567-577): Transforms J_ao[N], K_ao[N] for each N independently

**Formula:**
```
J[N]_mn = (mn|ls) C[N]_li C[N]_si
K[N]_mn = (ml|ns) C[N]_li C[N]_si
```

**Code evidence:**
```cpp
for (size_t N = 0; N < D_.size(); ++N) {
    D_[N]->zero();
    // Compute D[N] from C_left[N] and C_right[N] ONLY
    // ...
}

// Later:
if (do_J_) {
    double** JAOp = J_ao_[N]->pointer();  // J for THIS N only
    // ... transform J[N]
}
```

## Conclusion

**There is NO physical coupling between wavefunctions in multi_scf!**

Each wavefunction receives J/K computed from its OWN density only.
The batching is purely computational (efficiency), not physical.

## Implication

**multi_scf([UHF, ROHF]) should converge in SAME iterations as independent UHF+ROHF!**

But reality:
- Independent UHF: 7 iterations
- multi_scf([UHF, ROHF]): 16 iterations (+9)

**Therefore: The +9 iterations is NOT due to JK coupling!**

## Where is the problem?

Possible locations:
1. **Python layer** - Options/initialization differences
2. **Convergence checks** - Different thresholds or logic
3. **DIIS management** - Interference between wfn's DIIS
4. **Precomputed JK path** - Bug in our set_jk_matrices implementation
5. **Grace pattern side effects** - Unintended consequences

## Next Investigation

Since there's no physical coupling, the +9 iterations is a BUG in our implementation!

**Action:** Debug Python layer to find where UHF behavior diverges from independent SCF.

**Key questions:**
1. Why does multi_scf([UHF only]) add +1 iteration?
2. What changes between iter 4 and iter 5 that causes UHF divergence?
3. Are options/convergence checks identical?
4. Is DIIS working correctly?
5. Are J/K matrices set correctly via set_jk_matrices()?

