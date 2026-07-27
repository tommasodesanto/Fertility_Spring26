#!/usr/bin/env python3
"""Strict fixed-E5b diagnostics for E6b alone and combined with E6a."""

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

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
DEFAULT_SOURCE = (
    ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
)
DEFAULT_CERTIFIED = DEFAULT_SOURCE.parent / "target_fit_full.csv"
DEFAULT_OUTDIR = (
    ROOT / "output/model/eqscale_seq_e6b_income_level_diagnostic_20260727"
)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty output {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_chain() -> Any:
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    if str(MODEL_ROOT) not in sys.path:
        sys.path.insert(0, str(MODEL_ROOT))
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    return chain


def load_winner(source: Path) -> dict[str, float]:
    payload = json.loads(source.read_text(encoding="utf-8"))
    raw = payload["winners"]["E1"]["theta"]
    return {name: float(value) for name, value in raw.items()}


def variants() -> list[tuple[str, dict[str, Any]]]:
    from intergen_eqscale_seq_optimized.e6a_profile import e6a_overrides
    from intergen_eqscale_seq_optimized.e6b_profile import e6b_overrides

    return [
        ("e5b_current", {}),
        ("e6b_permanent_income", e6b_overrides()),
        ("e6ab_tail_and_permanent_income", {**e6a_overrides(), **e6b_overrides()}),
    ]


def verify_baseline(
    moments: dict[str, Any],
    certified_path: Path,
    tolerance: float = 1e-6,
) -> None:
    with certified_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 12:
        raise ValueError(f"Expected 12 certified rows, found {len(rows)}")
    failures = []
    for row in rows:
        name = row["moment"]
        difference = abs(float(moments[name]) - float(row["model"]))
        if difference > tolerance:
            failures.append((name, difference))
    if failures:
        raise RuntimeError(f"Certified E5b reproduction gate failed: {failures}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--certified", type=Path, default=DEFAULT_CERTIFIED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--Nb", type=int, default=120)
    parser.add_argument("--max-iter-eq", type=int, default=40)
    parser.add_argument("--tol-eq", type=float, default=2.5e-5)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25

    args.outdir.mkdir(parents=True, exist_ok=True)
    chain = load_chain()
    theta = load_winner(args.source)
    system = __import__(
        "intergen_eqscale_seq_optimized.e5_profile",
        fromlist=["e5_target_system"],
    ).e5_target_system()
    targets = system.targets_dict()
    weights = system.weights_dict()
    summary_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []

    for name, variant_overrides in variants():
        overrides = chain.common_overrides(
            argparse.Namespace(
                J=17,
                Nb=args.Nb,
                max_iter_eq=args.max_iter_eq,
                tol_eq=args.tol_eq,
            )
        )
        overrides.update(variant_overrides)
        overrides.update(theta)
        started = time.perf_counter()
        sol, parameters, price = chain.run_model_cp_dt(overrides, verbose=False)
        elapsed = time.perf_counter() - started
        moments = chain.extract_moments(sol, parameters)
        if name == "e5b_current" and not args.smoke:
            verify_baseline(moments, args.certified)
        fit = chain.target_fit(moments, targets, weights)
        loss = float(sum(float(row["loss_contribution"]) for row in fit))
        residual = float(sol.best_max_abs_rel_excess)
        strict = bool(
            getattr(sol, "converged", False)
            and residual <= float(args.tol_eq)
        )
        summary_rows.append(
            {
                "variant": name,
                "loss": loss,
                "market_residual": residual,
                "strict_converged": strict,
                "elapsed_seconds": elapsed,
                "income_states": int(parameters.Nz),
                "tfr": float(moments["tfr"]),
                "childless_rate": float(moments["childless_rate"]),
                "mean_age_first_birth": float(moments["mean_age_first_birth"]),
                "share_first_births_age30plus": float(
                    moments["share_first_births_age30plus"]
                ),
                "own_family_gap": float(moments["own_family_gap"]),
                "own_rate": float(moments["own_rate"]),
                "aggregate_wealth_to_annual_gross_labor_earnings": float(
                    moments["aggregate_wealth_to_annual_gross_labor_earnings"]
                ),
                "old_total_wealth_to_annual_income_p90_p50_7684": float(
                    moments["old_total_wealth_to_annual_income_p90_p50_7684"]
                ),
                "permanent_income_level_values": json.dumps(
                    np.asarray(
                        moments["permanent_income_level_values"], dtype=float
                    ).tolist()
                ),
                "permanent_income_childless_by_level": json.dumps(
                    np.asarray(
                        moments["permanent_income_childless_by_level"], dtype=float
                    ).tolist()
                ),
                "permanent_income_completed_fertility_by_level": json.dumps(
                    np.asarray(
                        moments["permanent_income_completed_fertility_by_level"],
                        dtype=float,
                    ).tolist()
                ),
                "permanent_income_own_rate_3055_by_level": json.dumps(
                    np.asarray(
                        moments["permanent_income_own_rate_3055_by_level"],
                        dtype=float,
                    ).tolist()
                ),
            }
        )
        for row in fit:
            target_rows.append({"variant": name, **row})
        print(
            f"{name}: loss={loss:.6f}, residual={residual:.3e}, "
            f"elapsed={elapsed:.1f}s, Nz={parameters.Nz}",
            flush=True,
        )

    write_csv(args.outdir / "summary.csv", summary_rows)
    write_csv(args.outdir / "target_fit_full.csv", target_rows)
    metadata = {
        "status": "smoke" if args.smoke else "strict_fixed_E5b_diagnostic",
        "target_set": system.name,
        "target_count": len(targets),
        "free_parameter_count": 10,
        "source": str(args.source),
        "e6b": __import__(
            "intergen_eqscale_seq_optimized.e6b_profile",
            fromlist=["e6b_metadata"],
        ).e6b_metadata(),
        "solver": {
            "Nb": args.Nb,
            "max_iter_eq": args.max_iter_eq,
            "tol_eq": args.tol_eq,
        },
    }
    (args.outdir / "metadata.json").write_text(
        json.dumps(chain.jsonable(metadata), indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
