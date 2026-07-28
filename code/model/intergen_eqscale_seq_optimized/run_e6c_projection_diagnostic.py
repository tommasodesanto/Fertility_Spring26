#!/usr/bin/env python3
"""Strictly project the certified E6ABC winner onto the E6AB architecture."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any


MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.e6c_profile import e6c_overrides  # noqa: E402
from intergen_eqscale_seq_optimized.local_panel import jsonable  # noqa: E402
from intergen_eqscale_seq_optimized.parameters import (  # noqa: E402
    readiness_cumulative_probability,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    payload = json.loads(args.source.read_text())
    winner = ((payload.get("winners") or {}).get("E1"))
    metadata = payload.get("winner_metadata") or {}
    repeat = payload.get("winner_repeat_check") or {}
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise ValueError(f"{args.source} does not contain winners.E1.theta")
    if metadata.get("arm", winner.get("arm")) != "E6ABC":
        raise ValueError("projection source must be a certified E6ABC winner")
    if not (
        repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", math.inf)) == 0.0
        and float(repeat.get("max_abs_moment_difference", math.inf)) == 0.0
    ):
        raise ValueError("projection source must preserve exact strict repeats")

    os.environ.update(
        {
            "E3_L4": "1",
            "E5": "1",
            "E6A": "1",
            "E6B": "1",
            "E3_TFR_TOP_BIN_WEIGHT": "3.602359422009",
        }
    )
    os.environ.pop("E6C", None)
    from intergen_eqscale_seq_optimized import run_e1_chain
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    runtime_args = argparse.Namespace(
        J=int(metadata.get("J", 17)),
        Nb=int(metadata.get("Nb", 120)),
        max_iter_eq=40,
        tol_eq=2.5e-5,
    )
    common = chain.common_overrides(runtime_args)
    theta = {name: float(value) for name, value in winner["theta"].items()}
    system = e5_target_system()
    targets, weights = system.targets_dict(), system.weights_dict()
    records: list[dict[str, Any]] = []
    fits: list[dict[str, Any]] = []

    for label, gate in (("E6ABC", True), ("E6AB_projection", False)):
        overrides = {**common, **theta}
        if gate:
            overrides.update(e6c_overrides())
        started = time.perf_counter()
        solution, parameters, prices = chain.run_model_cp_dt(
            overrides,
            verbose=False,
        )
        elapsed = time.perf_counter() - started
        moments = chain.extract_moments(solution, parameters)
        residual = float(
            getattr(solution, "best_max_abs_rel_excess", math.inf)
        )
        timings = dict(getattr(solution, "timings", {}))
        strict = bool(
            timings.get(
                "strict_converged",
                getattr(solution, "converged", False),
            )
            and residual <= 2.5e-5
        )
        if not strict:
            raise RuntimeError(
                f"{label} is not a strict solve: residual={residual:.6g}"
            )
        loss = float(
            chain.diagnostic_loss(
                moments,
                targets=targets,
                weights=weights,
            )
        )
        records.append(
            {
                "configuration": label,
                "readiness_gate_enabled": gate,
                "loss": loss,
                "market_residual": residual,
                "strict_converged": strict,
                "elapsed_sec": elapsed,
                "price": float(prices[0]),
                "entry_settled_probability": (
                    readiness_cumulative_probability(
                        parameters,
                        float(parameters.age_start),
                    )
                    if gate
                    else 1.0
                ),
            }
        )
        for name, target in targets.items():
            model = float(moments[name])
            gap = model - float(target)
            fits.append(
                {
                    "configuration": label,
                    "moment": name,
                    "target": float(target),
                    "model": model,
                    "gap": gap,
                    "weight": float(weights[name]),
                    "loss_contribution": float(weights[name]) * gap**2,
                }
            )

    certified = records[0]
    if not math.isclose(
        float(certified["loss"]),
        float(winner["rank_loss"]),
        rel_tol=0.0,
        abs_tol=1e-8,
    ):
        raise RuntimeError("E6ABC projection diagnostic fails loss reproduction")
    certified_fit = {
        str(row["moment"]): float(row["model"])
        for row in winner["target_fit"]
    }
    reproduced_fit = {
        str(row["moment"]): float(row["model"])
        for row in fits
        if row["configuration"] == "E6ABC"
    }
    maximum_difference = max(
        abs(reproduced_fit[name] - certified_fit[name])
        for name in certified_fit
    )
    if maximum_difference > 1e-6:
        raise RuntimeError(
            "E6ABC projection diagnostic fails moment reproduction: "
            f"{maximum_difference:.6g}"
        )

    projection = records[1]
    summary = {
        "source": str(args.source),
        "certified_e6abc_loss": certified["loss"],
        "projected_e6ab_loss": projection["loss"],
        "projected_minus_certified_loss": (
            float(projection["loss"]) - float(certified["loss"])
        ),
        "entry_settled_probability": certified["entry_settled_probability"],
        "maximum_absolute_certified_moment_difference": maximum_difference,
        "all_cases_strict": all(bool(record["strict_converged"]) for record in records),
        "cases": records,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "strict_cases.csv", records)
    write_csv(args.outdir / "target_fit_full.csv", fits)
    (args.outdir / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )


if __name__ == "__main__":
    main()
