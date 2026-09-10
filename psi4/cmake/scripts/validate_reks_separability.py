#!/usr/bin/env python3
"""Validate REKS SA separability: enforce H_nm = 0 theorem at build time.

For every REKS SA SCF module, the theorem proven in _CRITIQUE_MFON_TRAH.md sec 7.1
requires three structural invariants:

  (1) compute_weights_PPS* / compute_weight_derivs_PPS* must NOT read m_geminals_.
      (Psi_0 branch is pure n-FON.)
  (2) compute_weight_derivs_OSS[0-9]+_oss_fon* must NOT read n_geminals_.
      (Psi_1 m-branch is pure m-FON.)
  (3) compute_weight_derivs_OSS[0-9]+ (non-oss_fon) bodies must not emit any
      indexed dC_dfon / d2C_dfon2 assignments.
      (Psi_1 n-branch is identically zero.)

If any check fails, H_nm is no longer guaranteed zero and the zero-block
assumption in REKSGradientEngine::compute is wrong. Remediation is to
explicitly compute H_nm; do NOT re-introduce a default-zero virtual hook.
"""

import re
import sys
from pathlib import Path

FILES = [
    "reks22_scheme0_spin0_2sa_scf.cc",
    "reks44_scheme0_spin0_3sa_scf.cc",
    "reks44_scheme1_spin0_3sa_scf.cc",
    "reks66_scheme0_spin0_4sa_scf.cc",
]

RE_PPS_SIG = re.compile(
    r"void\s+REKS[0-9]+_[A-Za-z0-9_]+::"
    r"(compute_weights_PPS|compute_weight_derivs_PPS)[A-Za-z0-9_]*\s*\("
)
RE_OSS_OSS_FON_SIG = re.compile(
    r"void\s+REKS[0-9]+_[A-Za-z0-9_]+::"
    r"compute_weight_derivs_OSS[0-9]+_oss_fon[A-Za-z0-9_]*\s*\("
)
# Match OSS[0-9]+ followed immediately by opening paren — excludes _oss_fon variants.
RE_OSS_N_SIG = re.compile(
    r"void\s+REKS[0-9]+_[A-Za-z0-9_]+::"
    r"compute_weight_derivs_OSS[0-9]+\s*\("
)

RE_M_GEMINALS = re.compile(r"\bm_geminals_\b")
RE_N_GEMINALS = re.compile(r"\bn_geminals_\b")
RE_DC_INDEX = re.compile(r"\b(dC_dfon|d2C_dfon2)\[")


def iter_fn_bodies(lines, sig_re):
    """Yield (start_line, [lines...]) for every function whose signature matches."""
    in_fn = False
    start = 0
    body = []
    for idx, line in enumerate(lines, 1):
        if not in_fn:
            if sig_re.search(line):
                in_fn = True
                start = idx
                body = [line]
        else:
            body.append(line)
            if line.startswith("}"):
                yield start, body
                in_fn = False
                body = []


def check_body_for(body, start, bad_re, description, filename):
    violations = []
    for off, line in enumerate(body):
        if bad_re.search(line):
            violations.append(f"  {filename}:{start + off}: {line.rstrip()}")
    if violations:
        return [f"[validator] {description}"] + violations
    return []


def validate_file(path: Path):
    lines = path.read_text().splitlines()
    name = path.name
    errors = []

    for start, body in iter_fn_bodies(lines, RE_PPS_SIG):
        errors.extend(check_body_for(
            body, start, RE_M_GEMINALS,
            f"{name}: PPS helper at line {start} reads m_geminals_ (breaks n-FON purity)",
            name,
        ))
    for start, body in iter_fn_bodies(lines, RE_OSS_OSS_FON_SIG):
        errors.extend(check_body_for(
            body, start, RE_N_GEMINALS,
            f"{name}: OSSx_oss_fon helper at line {start} reads n_geminals_ (breaks m-FON purity)",
            name,
        ))
    for start, body in iter_fn_bodies(lines, RE_OSS_N_SIG):
        errors.extend(check_body_for(
            body, start, RE_DC_INDEX,
            f"{name}: OSSx n-layer helper at line {start} emits indexed dC/d2C assignment (breaks separability)",
            name,
        ))
    return errors


def main():
    if len(sys.argv) >= 2:
        src_dir = Path(sys.argv[1])
    else:
        src_dir = Path.cwd() / "psi4" / "src" / "psi4" / "libscf_solver"

    all_errors = []
    checked = 0
    for rel in FILES:
        path = src_dir / rel
        if not path.is_file():
            print(f"[validator] missing source: {path}", file=sys.stderr)
            return 2
        checked += 1
        all_errors.extend(validate_file(path))

    if all_errors:
        for line in all_errors:
            print(line, file=sys.stderr)
        print("", file=sys.stderr)
        print("=" * 80, file=sys.stderr)
        print("REKS SA separability VIOLATED.", file=sys.stderr)
        print("", file=sys.stderr)
        print(
            "H_nm = 0 theorem no longer holds across all SA variants.\n"
            "H_nm is the n-FON / m-FON cross Hessian block in\n"
            "REKSGradientEngine::compute; it is assumed identically zero.\n"
            "A violation means that zero block is now wrong and must be\n"
            "replaced by an explicit computation.\n\n"
            "Remediation:\n"
            "  - Add compute_nm_cross_hessian(int g_n, int g_m) to\n"
            "    REKSGradientEngine, parallel to compute_cross_hessian.\n"
            "  - Add virtual compute_nm_mixed_weight_derivs(g_n, g_m) to\n"
            "    REKSActiveSpace returning the mixed derivative vector per\n"
            "    config for each new variant.\n"
            "  - Fill the (n_rot+g_n, n_rot+n_fon+g_m) symmetric entries\n"
            "    in the assembled H matrix.\n\n"
            "See _CRITIQUE_MFON_TRAH.md section 7.1 for the full theorem.",
            file=sys.stderr,
        )
        print("=" * 80, file=sys.stderr)
        return 1

    print(f"[validator] REKS SA separability: OK ({checked} modules checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
