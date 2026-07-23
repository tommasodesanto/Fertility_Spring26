#!/usr/bin/env python3
"""Certify the timing-repaired target readout at a supplied parameter vector."""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .calibration import extract_moments
from .calibration_search import (
    load_start_theta,
    parameter_rows,
    target_fit_rows,
    wealth_components,
)
from .new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
    new_moment_target_system,
)
from .solver import run_model_cp_dt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--jacobian-summary", type=Path, default=None)
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.repeats < 2:
        raise ValueError("--repeats must be at least two for certification")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    system = new_moment_target_system()
    theta = load_start_theta(args.source)

    records: list[dict[str, Any]] = []
    for repeat in range(args.repeats):
        start = time.perf_counter()
        solution, parameters, price = run_model_cp_dt(
            {**new_moment_overrides(tight=True, optimized=True), **theta},
            verbose=False,
        )
        moments_all = extract_moments(solution, parameters)
        moments = {
            name: float(moments_all[name])
            for name in system.moment_names
        }
        records.append(
            {
                "repeat": repeat + 1,
                "elapsed_seconds": time.perf_counter() - start,
                "loss": float(system.loss(moments)),
                "market_residual": float(solution.best_max_abs_rel_excess),
                "price": float(np.asarray(price, dtype=float).reshape(-1)[0]),
                "moments": moments,
                "wealth_components": wealth_components(solution),
                "wealth_timing": str(solution.wealth_moment_timing),
                "bequest_timing": str(solution.bequest_moment_timing),
            }
        )
        print(
            f"repeat {repeat + 1}/{args.repeats}: "
            f"loss={records[-1]['loss']:.9f}, "
            f"residual={records[-1]['market_residual']:.3e}",
            flush=True,
        )

    first = records[0]
    repeat_identical = all(
        record["loss"] == first["loss"]
        and record["market_residual"] == first["market_residual"]
        and record["price"] == first["price"]
        and record["moments"] == first["moments"]
        and record["wealth_components"] == first["wealth_components"]
        for record in records[1:]
    )
    fit = target_fit_rows(first["moments"], system)
    parameters_table = parameter_rows(theta)
    jacobian = None
    if args.jacobian_summary is not None:
        jacobian = json.loads(args.jacobian_summary.read_text())

    summary = {
        "status": (
            "certified_repaired_readout_not_recalibrated"
            if repeat_identical
            else "failed_repeat_gate"
        ),
        "profile": NEW_MOMENT_PROFILE_NAME,
        "target_fingerprint": system.fingerprint,
        "source_parameter_file": str(args.source.resolve()),
        "theta": theta,
        "strict_repeats_bit_identical": repeat_identical,
        "strict_repeats": records,
        "target_fit": fit,
        "parameter_table": parameters_table,
        "measurement_contract": {
            "cross_sectional_wealth": "beginning_of_period_b_plus_gross_housing",
            "living_old_wealth_dispersion": "beginning_of_period_living_psid_analogue",
            "annual_bequest_flow": "post_saving_b_prime_plus_gross_housing_at_death",
        },
        "jacobian": jacobian,
    }
    write_json(outdir / "summary.json", summary)
    write_json(
        outdir / "results.json",
        {
            "status": summary["status"],
            "profile": NEW_MOMENT_PROFILE_NAME,
            "target_fingerprint": system.fingerprint,
            "theta": theta,
            "loss": first["loss"],
            "market_residual": first["market_residual"],
            "moments": first["moments"],
        },
    )
    write_csv(outdir / "target_fit.csv", fit)
    write_csv(outdir / "parameters.csv", parameters_table)

    normalized_gaps = np.asarray(
        [(float(row["model"]) - float(row["target"])) / abs(float(row["target"])) for row in fit]
    )
    labels = [str(row["moment"]) for row in fit]
    figure, axis = plt.subplots(figsize=(9.0, 6.4))
    positions = np.arange(len(labels))
    colors = np.where(normalized_gaps >= 0.0, "#B24C3D", "#366E8A")
    axis.barh(positions, 100.0 * normalized_gaps, color=colors)
    axis.axvline(0.0, color="black", linewidth=0.8)
    axis.set_yticks(positions, labels)
    axis.invert_yaxis()
    axis.set_xlabel("Model minus target (percent of target)")
    axis.set_title("Supplemental: timing-repaired target gaps at the old parameter vector")
    axis.grid(axis="x", alpha=0.2)
    figure.tight_layout()
    figure.savefig(outdir / "target_gaps.png", dpi=180)
    plt.close(figure)

    top = sorted(fit, key=lambda row: float(row["loss_contribution"]), reverse=True)
    jacobian_text = ""
    if jacobian is not None:
        jacobian_text = (
            f"The repaired local Jacobian is numerically full rank "
            f"({jacobian['rank']}/{jacobian['required_rank']}) but has condition "
            f"number {float(jacobian['condition_number']):,.0f}. This is weak, "
            "joint local identification—not a clean one-to-one mapping.\n"
        )
    report = f"""# Timing-repaired calibration audit

This packet evaluates the previously reported parameter vector under the
repaired measurement contract. It is **not** a recalibration.

## Result

- Strict loss: `{float(first["loss"]):.9f}`.
- Strict repeats bit-identical: `{repeat_identical}`.
- Market residual: `{float(first["market_residual"]):.3e}`.
- Wealth timing: beginning-of-period living-household balance sheet.
- Bequest timing: post-saving resources at the death node.
- Living-old dispersion sample: living PSID households aged 76--84.

The three largest loss contributions are:

1. `{top[0]["moment"]}`: model `{float(top[0]["model"]):.6f}`, target `{float(top[0]["target"]):.6f}`, contribution `{float(top[0]["loss_contribution"]):.6f}`.
2. `{top[1]["moment"]}`: model `{float(top[1]["model"]):.6f}`, target `{float(top[1]["target"]):.6f}`, contribution `{float(top[1]["loss_contribution"]):.6f}`.
3. `{top[2]["moment"]}`: model `{float(top[2]["model"]):.6f}`, target `{float(top[2]["target"]):.6f}`, contribution `{float(top[2]["loss_contribution"]):.6f}`.

{jacobian_text}
The old loss `0.021954` cannot be compared as a valid fit because its three
wealth/bequest rows used a hybrid balance sheet. The complete repaired target
table is in `target_fit.csv`; all free parameters, bounds, and bound flags are
in `parameters.csv`.
"""
    (outdir / "REPORT.md").write_text(report)

    if not repeat_identical:
        raise RuntimeError("strict repaired readout failed the exact-repeat gate")


if __name__ == "__main__":
    main()
