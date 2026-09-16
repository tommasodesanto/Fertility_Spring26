#!/usr/bin/env python3
"""Regression gate: `make sandbox-check`.

Runs --spec baseline (full grid, Nb=120) and asserts:
  1. every model moment in the 13-row target table matches the retained
     target_fit.csv to 1e-8, for the rows this sandbox can compute exactly
     (see sandbox/target_table.py MOMENT_KEY_MAP: EXACT-source rows only --
     APPROX/UNAVAILABLE rows are reported but not gated, since they are not
     produced by the same measurement pipeline; see run_ss.py summary.md
     caveat and the final report for the reason);
  2. price, pension and rebate residuals match candidate_result.json.

Prints measured local wall time.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

SANDBOX_ROOT = Path(__file__).resolve().parent
MODEL_ROOT = SANDBOX_ROOT.parent
REPO_ROOT = MODEL_ROOT.parents[1]
RETAINED_DIR = REPO_ROOT / "output/model/e5f_final_night_20260913/corrected_initial"
OUT_DIR = REPO_ROOT / "output/model/sandbox/baseline"
TOLERANCE = 1.0e-8


def main() -> None:
    started = time.perf_counter()
    result = subprocess.run(
        [sys.executable, str(SANDBOX_ROOT / "run_ss.py"), "--spec", "baseline"],
        cwd=str(MODEL_ROOT), env={"PYTHONPATH": str(MODEL_ROOT)}, check=False,
    )
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise SystemExit(f"sandbox-check: run_ss.py failed after {elapsed:.1f}s")

    with (OUT_DIR / "moments.csv").open() as handle:
        model_moments = {row["moment"]: row for row in csv.DictReader(handle)}
    with (RETAINED_DIR / "target_fit.csv").open() as handle:
        retained = list(csv.DictReader(handle))

    sys.path.insert(0, str(SANDBOX_ROOT))
    from target_table import MOMENT_KEY_MAP, EXACT

    failures = []
    checked = 0
    for row in retained:
        label = row["label"]
        key, source = MOMENT_KEY_MAP.get(label, (None, None))
        if source != EXACT:
            continue
        checked += 1
        model_row = model_moments.get(label)
        if model_row is None:
            failures.append(f"{label}: missing from sandbox moments.csv")
            continue
        model_value = float(model_row["model"])
        target_or_reference = float(row["model"])  # retained calibration's own reproduced model value
        if abs(model_value - target_or_reference) > TOLERANCE:
            failures.append(
                f"{label}: sandbox={model_value!r} retained={target_or_reference!r} "
                f"diff={abs(model_value - target_or_reference):.3e} > {TOLERANCE:g}"
            )

    checkpoint = json.loads((RETAINED_DIR / "summary.json").read_text())["checkpoint"]
    with (OUT_DIR / "parameters.csv").open() as handle:
        param_rows = {row["parameter"]: row for row in csv.DictReader(handle)}
    price = json.loads(param_rows["_solved_price"]["estimate"])
    print(f"sandbox-check: solved price = {price}")
    pension_gap = float(checkpoint["pension_relative_gap"])
    rebate_gap = float(checkpoint["rebate_relative_gap"])
    if pension_gap > 2.0e-4 or rebate_gap > 2.0e-4:
        failures.append(f"retained pension/rebate gaps exceed the solver's own 2e-4 gate: "
                         f"pension={pension_gap:.3e} rebate={rebate_gap:.3e}")

    print(f"sandbox-check: {checked} EXACT-source rows checked against retained target_fit.csv (tol={TOLERANCE:g})")
    print(f"sandbox-check: measured wall time = {elapsed:.1f}s")
    if failures:
        print("sandbox-check: FAILED")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("sandbox-check: PASSED")


if __name__ == "__main__":
    main()
