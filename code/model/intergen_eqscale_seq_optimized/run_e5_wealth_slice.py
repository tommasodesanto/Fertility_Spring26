#!/usr/bin/env python3
"""Evaluate an E5 winner on a grid of annual discount factors and theta0."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
DEFAULT_SOURCE = ROOT / "output/model/eqscale_seq_e5_recalibration_20260724/report/results.json"


def comma_floats(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated float list")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--beta-annual", type=comma_floats, required=True)
    parser.add_argument("--theta0", type=comma_floats, required=True)
    parser.add_argument("--Nb", type=int, default=120)
    parser.add_argument("--max-iter-eq", type=int, default=40)
    parser.add_argument("--tol-eq", type=float, default=2.5e-5)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def load_e5_chain() -> Any:
    """Load the E5 chain after setting the same environment contract as E5."""
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    if str(MODEL_ROOT) not in sys.path:
        sys.path.insert(0, str(MODEL_ROOT))
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    return chain


def load_winner_theta(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    theta = ((payload.get("winners") or {}).get("E1") or {}).get("theta")
    if not isinstance(theta, dict):
        raise ValueError(f"{path} does not contain winners.E1.theta")
    return {name: float(value) for name, value in theta.items()}


def existing_nodes(path: Path) -> set[tuple[float, float]]:
    if not path.exists():
        return set()
    with path.open(newline="", encoding="utf-8") as handle:
        return {(float(row["beta_annual"]), float(row["theta0"])) for row in csv.DictReader(handle)}


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.Nb, args.max_iter_eq, args.tol_eq = 60, 2, 0.25
    chain = load_e5_chain()
    theta = load_winner_theta(args.source)
    e5_system = __import__("intergen_eqscale_seq_optimized.e5_profile", fromlist=["e5_target_system"]).e5_target_system()
    moment_names = list(e5_system.moment_names)
    args.outdir.mkdir(parents=True, exist_ok=True)
    outpath = args.outdir / "wealth_slice.csv"
    done = existing_nodes(outpath)
    fields = ["beta_annual", "theta0", "market_residual", *moment_names, "own_rate_65_75"]
    write_header = not outpath.exists() or outpath.stat().st_size == 0
    with outpath.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if write_header:
            writer.writeheader()
            handle.flush()
        for beta_annual in args.beta_annual:
            for theta0 in args.theta0:
                node = (float(beta_annual), float(theta0))
                if node in done:
                    print(f"skip beta_annual={beta_annual:.10g} theta0={theta0:.10g}", flush=True)
                    continue
                runtime_args = argparse.Namespace(J=17, Nb=args.Nb, max_iter_eq=args.max_iter_eq, tol_eq=args.tol_eq)
                overrides = chain.common_overrides(runtime_args)
                overrides.update(theta)
                overrides["beta"] = float(beta_annual) ** 4.0
                overrides["theta0"] = float(theta0)
                sol, P, _ = chain.run_model_cp_dt(overrides, verbose=False)
                moments = chain.extract_moments(sol, P)
                row = {"beta_annual": float(beta_annual), "theta0": float(theta0),
                       "market_residual": float(getattr(sol, "best_max_abs_rel_excess", math.nan)),
                       **{name: float(moments.get(name, math.nan)) for name in moment_names},
                       "own_rate_65_75": float(moments.get("own_rate_65_75", math.nan))}
                writer.writerow(row)
                handle.flush()
                print(f"beta_annual={beta_annual:.10g} theta0={theta0:.10g} residual={row['market_residual']:.3e}", flush=True)


if __name__ == "__main__":
    main()
