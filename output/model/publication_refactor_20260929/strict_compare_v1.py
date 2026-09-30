"""Array-by-array comparison of a lab solution with the saved reference.

    python -m refactor_lab.compare --lab OUT/solution_arrays.npz \
        --reference EXPORT/verification/reference_solution.npz --out OUT/comparison.json

Also used for same-start lab vs old-engine GE (`--reference` = old solution).
Every array common to both is reported as exact, or with max abs/rel
difference. Missing arrays on either side are listed, never skipped
silently. `--strict-paths` additionally requires exact key equality, including
keys normally checked by another oracle. Exit status is nonzero unless every
common array is exact and no required key is missing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# Reference-export objects that are not attributes of a solution. Each is
# checked elsewhere; the name and its checker are recorded in the output.
NON_SOLUTION = {"stationary_g_pre": "acceptance_oracle.fixed_price: exact equality of the frozen "
                "stationary reconstruction from the lab policy with the checkpoint stationary_g_pre"}


def compare(lab: Path, reference: Path, *, strict_paths: bool = False) -> dict:
    rows, missing, extra = {}, [], []
    with np.load(lab, allow_pickle=False) as a, np.load(reference, allow_pickle=False) as b:
        excluded = set() if strict_paths else set(NON_SOLUTION)
        missing = sorted(set(b.files) - set(a.files) - excluded)
        extra = sorted(set(a.files) - set(b.files))
        for k in sorted(set(a.files) & set(b.files)):
            x, y = a[k], b[k]
            if x.shape != y.shape or x.dtype != y.dtype:
                rows[k] = dict(exact=False, reason=f"shape/dtype {x.shape}/{x.dtype} vs {y.shape}/{y.dtype}")
                continue
            exact = bool(np.array_equal(x, y, equal_nan=x.dtype.kind == "f"))
            row = dict(exact=exact, shape=list(x.shape))
            if not exact and x.dtype.kind in "fc":
                d = np.abs(x - y)
                row.update(max_abs=float(np.nanmax(d)),
                           max_rel=float(np.nanmax(d / np.maximum(np.abs(y), 1e-300))),
                           nan_mismatch=int(np.count_nonzero(np.isnan(x) != np.isnan(y))))
            rows[k] = row
    passed = not missing and all(r["exact"] for r in rows.values()) and (not strict_paths or not extra)
    return dict(passed=passed, strict_paths=strict_paths, compared=len(rows), missing_in_lab=missing, lab_only=extra,
                checked_elsewhere=NON_SOLUTION,
                not_exact=sorted(k for k, r in rows.items() if not r["exact"]), arrays=rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lab", type=Path, required=True)
    ap.add_argument("--reference", type=Path, help="flat reference arrays (checkpoint export or old GE)")
    ap.add_argument("--baseline-dir", type=Path, help="pinned same-machine original-engine baseline")
    ap.add_argument("--baseline-pin")
    ap.add_argument("--root", type=Path, help="reference ROOT for baseline identity")
    ap.add_argument("--record-only", action="store_true",
                    help="diagnostic: write the comparison and exit 0 (e.g. Mac lab vs Torch checkpoint)")
    ap.add_argument("--strict-paths", action="store_true",
                    help="require exact array-key equality, including NON_SOLUTION keys")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.baseline_dir:
        from refactor_lab import baseline_identity
        if not (a.baseline_pin and a.root) or a.reference:
            raise SystemExit("--baseline-dir needs --baseline-pin and --root, and excludes --reference")
        baseline_identity.verify(a.root, a.baseline_dir, a.baseline_pin)   # before any use
        a.reference = a.baseline_dir / "solution_arrays.npz"
    if a.reference is None:
        raise SystemExit("--reference or --baseline-dir required")
    result = compare(a.lab, a.reference, strict_paths=a.strict_paths)
    result["mode"] = "diagnostic_record_only" if a.record_only else "required_exact"
    a.out.write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("passed", "compared", "missing_in_lab", "not_exact", "mode")}, indent=1))
    if (not result["passed"] and not a.record_only) or (
            a.strict_paths and (result["missing_in_lab"] or result["lab_only"])):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
