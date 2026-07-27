#!/usr/bin/env python3
"""Strictly reproduce a certified E6 winner and write the standard graph packet."""

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

import numpy as np


MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.diagnostics import write_diagnostics  # noqa: E402
from intergen_eqscale_seq_optimized.local_panel import jsonable  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--arm", choices=("e6a", "e6b", "e6ab"), required=True)
    return parser.parse_args()


def configure_environment(arm: str) -> None:
    for name in ("E6A", "E6B"):
        os.environ.pop(name, None)
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    if arm in {"e6a", "e6ab"}:
        os.environ["E6A"] = "1"
    if arm in {"e6b", "e6ab"}:
        os.environ["E6B"] = "1"


def expected_arm(arm: str) -> str:
    return {"e6a": "E6A", "e6b": "E6B", "e6ab": "E6AB"}[arm]


def read_certified(source: Path, arm: str) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = json.loads(source.read_text())
    winner = ((payload.get("winners") or {}).get("E1"))
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise ValueError(f"{source} does not contain winners.E1.theta")
    metadata = payload.get("winner_metadata") or {}
    reported_arm = metadata.get("arm", winner.get("arm"))
    if reported_arm != expected_arm(arm):
        raise ValueError(
            f"{source} reports arm {reported_arm!r}, expected {expected_arm(arm)!r}"
        )
    repeat = payload.get("winner_repeat_check") or {}
    if not (
        repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", math.inf)) == 0.0
        and float(repeat.get("max_abs_moment_difference", math.inf)) == 0.0
    ):
        raise ValueError(f"{source} does not preserve an exact strict repeat")
    return winner, metadata


def certified_rows(source: Path) -> list[dict[str, Any]]:
    path = source.parent / "target_fit_full.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 12 or len({row["moment"] for row in rows}) != 12:
        raise ValueError(f"{path} must contain twelve unique target rows")
    return rows


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
    winner, metadata = read_certified(args.source, args.arm)
    fit = certified_rows(args.source)
    configure_environment(args.arm)
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    runtime_args = argparse.Namespace(
        J=int(metadata.get("J", 17)),
        Nb=int(metadata.get("Nb", 120)),
        max_iter_eq=40,
        tol_eq=2.5e-5,
    )
    overrides = chain.common_overrides(runtime_args)
    theta = {name: float(value) for name, value in winner["theta"].items()}
    started = time.perf_counter()
    sol, parameters, prices = chain.run_model_cp_dt(
        {**overrides, **theta},
        verbose=False,
    )
    elapsed = time.perf_counter() - started
    moments = chain.extract_moments(sol, parameters)
    residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
    timings = dict(getattr(sol, "timings", {}))
    strict = bool(
        timings.get("strict_converged", getattr(sol, "converged", False))
        and residual <= 2.5e-5
    )
    if not strict:
        raise RuntimeError(
            f"diagnostic reproduction is not strict: residual={residual:.6g}"
        )

    reproduced_rows = []
    differences = []
    targets: dict[str, float] = {}
    weights: dict[str, float] = {}
    for row in fit:
        name = row["moment"]
        if name not in moments:
            raise RuntimeError(f"diagnostic reproduction omits {name}")
        target, weight = float(row["target"]), float(row["weight"])
        model = float(moments[name])
        certified_model = float(row["model"])
        difference = abs(model - certified_model)
        gap = model - target
        differences.append(difference)
        targets[name], weights[name] = target, weight
        reproduced_rows.append(
            {
                "moment": name,
                "target": target,
                "model": model,
                "certified_model": certified_model,
                "abs_difference_from_certified": difference,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap**2,
            }
        )
    maximum_difference = max(differences)
    if maximum_difference > 1e-6:
        raise RuntimeError(
            "certified moment verification failed: "
            f"maximum absolute difference={maximum_difference:.6g}"
        )
    loss = float(chain.diagnostic_loss(moments, targets=targets, weights=weights))
    if not math.isclose(
        loss,
        float(winner["rank_loss"]),
        rel_tol=0.0,
        abs_tol=1e-8,
    ):
        raise RuntimeError(
            f"loss reproduction failed: {loss:.15g} versus {winner['rank_loss']}"
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "target_fit_reproduced_full.csv", reproduced_rows)
    write_diagnostics(sol, parameters, args.outdir / "standard")
    permanent_diagnostics = {
        name: moments[name]
        for name in (
            "permanent_income_level_values",
            "childless_rate_by_permanent_income_level",
            "completed_fertility_by_permanent_income_level",
            "own_rate_3055_by_permanent_income_level",
            "childless_rate_high_minus_low_permanent_income",
            "completed_fertility_high_minus_low_permanent_income",
            "own_rate_3055_high_minus_low_permanent_income",
        )
        if name in moments
    }
    summary = {
        "arm": expected_arm(args.arm),
        "source": str(args.source),
        "strict_converged": strict,
        "market_residual": residual,
        "loss": loss,
        "certified_loss": float(winner["rank_loss"]),
        "maximum_absolute_certified_moment_difference": maximum_difference,
        "elapsed_sec": elapsed,
        "prices": np.asarray(prices, dtype=float).reshape(-1),
        "theta": theta,
        "permanent_income_diagnostics": permanent_diagnostics,
    }
    (args.outdir / "verification_summary.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )
    (args.outdir / "README.md").write_text(
        "\n".join(
            [
                "# Certified E6 diagnostic packet",
                "",
                "The winner is solved once at the strict collector settings.",
                "All twelve calibrated moments must reproduce within `1e-6` before",
                "the existing standard graph set or verification tables are written.",
                "No policy experiment is run.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
