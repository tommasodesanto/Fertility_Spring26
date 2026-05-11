"""Compare Python and MATLAB JSON benchmark outputs."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_benchmarks.py PYTHON_JSON MATLAB_JSON")
    py = json.loads(Path(sys.argv[1]).read_text())
    ml = json.loads(Path(sys.argv[2]).read_text())
    keys = [
        "tfr",
        "own_rate",
        "mean_age_first_birth",
        "migration_rate_2245",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "total_mass",
    ]
    rows = []
    max_abs = 0.0
    for key in keys:
        a = float(py[key])
        b = float(ml[key])
        diff = a - b
        max_abs = max(max_abs, abs(diff))
        rows.append((key, a, b, diff))
    for idx, (a, b) in enumerate(zip(py["p_eq"], ml["p_eq"])):
        diff = float(a) - float(b)
        max_abs = max(max_abs, abs(diff))
        rows.append((f"p_eq[{idx}]", float(a), float(b), diff))
    for idx, (a, b) in enumerate(zip(py["pop_share"], ml["pop_share"])):
        diff = float(a) - float(b)
        max_abs = max(max_abs, abs(diff))
        rows.append((f"pop_share[{idx}]", float(a), float(b), diff))

    print(f"{'moment':40s} {'python':>14s} {'matlab':>14s} {'diff':>14s}")
    print("-" * 86)
    for key, a, b, diff in rows:
        print(f"{key:40s} {a:14.8f} {b:14.8f} {diff:14.8g}")
    print("-" * 86)
    print(f"max_abs_diff = {max_abs:.8g}")


if __name__ == "__main__":
    main()

