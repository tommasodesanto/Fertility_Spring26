#!/usr/bin/env python3
"""Strict fixed-E5b diagnostics for the E6a fecundity-tail alternatives."""

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
    ROOT
    / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
)
DEFAULT_CERTIFIED = DEFAULT_SOURCE.parent / "target_fit_full.csv"
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e6a_fecundity_tail_frontier_20260727"

LERIDON_OMEGA1 = 0.01331
LERIDON_OMEGA2 = 0.14960
TAIL_START_AGE = 40.0
TERMINAL_AGE = 45.0


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


def terminal_decay_for_probability(
    probability: float,
    *,
    omega1: float = LERIDON_OMEGA1,
    omega2: float = LERIDON_OMEGA2,
) -> float:
    base_terminal = 1.0 - omega1 * math.exp(
        omega2 * (TERMINAL_AGE - 18.0)
    )
    if not (0.0 < probability < base_terminal):
        raise ValueError("terminal probability must lie between zero and the base schedule")
    return math.log(base_terminal / probability) / (
        TERMINAL_AGE - TAIL_START_AGE
    )


def schedules() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = [
        {
            "case": "e5b_current",
            "description": "Certified E5b schedule",
            "fecundity_omega1": 0.02,
            "fecundity_omega2": 0.134,
            "fecundity_tail_start_age": TAIL_START_AGE,
            "fecundity_terminal_decay": 0.0,
            "terminal_probability_before_hard_close": math.nan,
        },
        {
            "case": "leridon_two_parameter",
            "description": "Two-parameter fit to four-year conception anchors",
            "fecundity_omega1": LERIDON_OMEGA1,
            "fecundity_omega2": LERIDON_OMEGA2,
            "fecundity_tail_start_age": TAIL_START_AGE,
            "fecundity_terminal_decay": 0.0,
            "terminal_probability_before_hard_close": math.nan,
        },
    ]
    for probability in (0.10, 0.05, 0.03):
        rows.append(
            {
                "case": f"terminal_p{int(round(100 * probability)):02d}",
                "description": (
                    "Leridon anchors plus exponential terminal decay to "
                    f"{probability:.0%} just before age 45"
                ),
                "fecundity_omega1": LERIDON_OMEGA1,
                "fecundity_omega2": LERIDON_OMEGA2,
                "fecundity_tail_start_age": TAIL_START_AGE,
                "fecundity_terminal_decay": terminal_decay_for_probability(
                    probability
                ),
                "terminal_probability_before_hard_close": probability,
            }
        )
    return rows


def target_fit(
    moments: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
) -> list[dict[str, float | str]]:
    return [
        {
            "moment": name,
            "target": float(target),
            "model": float(moments[name]),
            "gap": float(moments[name]) - float(target),
            "weight": float(weights[name]),
            "loss_contribution": float(weights[name])
            * (float(moments[name]) - float(target)) ** 2,
        }
        for name, target in targets.items()
    ]


def solve_case(
    chain: Any,
    theta: dict[str, float],
    schedule: dict[str, Any],
) -> tuple[Any, Any, np.ndarray, dict[str, float], float]:
    overrides = chain.common_overrides(
        argparse.Namespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5)
    )
    overrides.update(theta)
    overrides.update(
        {
            "fecundity_omega1": schedule["fecundity_omega1"],
            "fecundity_omega2": schedule["fecundity_omega2"],
            "fecundity_tail_start_age": schedule["fecundity_tail_start_age"],
            "fecundity_terminal_decay": schedule["fecundity_terminal_decay"],
        }
    )
    began = time.perf_counter()
    sol, P, price = chain.run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - began
    moments = chain.extract_moments(sol, P)
    return sol, P, np.asarray(price, dtype=float).reshape(-1), moments, elapsed


def verify_baseline(
    moments: dict[str, float],
    certified: Path,
    tolerance: float = 1e-6,
) -> None:
    with certified.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    failed = []
    for row in rows:
        name = row["moment"]
        expected = float(row["model"])
        actual = float(moments[name])
        difference = abs(actual - expected)
        if difference > tolerance:
            failed.append(
                {
                    "moment": name,
                    "actual": actual,
                    "expected": expected,
                    "abs_difference": difference,
                }
            )
    if failed:
        raise RuntimeError(f"Certified E5b reproduction gate failed: {failed}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--certified", type=Path, default=DEFAULT_CERTIFIED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

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
    fit_rows: list[dict[str, Any]] = []
    schedule_rows: list[dict[str, Any]] = []
    baseline_models: dict[str, float] | None = None

    for schedule in schedules():
        print(f"START {schedule['case']}", flush=True)
        sol, P, price, moments, elapsed = solve_case(chain, theta, schedule)
        if schedule["case"] == "e5b_current":
            verify_baseline(moments, args.certified)
            baseline_models = {name: float(moments[name]) for name in targets}
        if baseline_models is None:
            raise RuntimeError("The certified E5b case must be evaluated first")

        rows = target_fit(moments, targets, weights)
        strict = bool(
            getattr(sol, "converged", False)
            and float(getattr(sol, "best_max_abs_rel_excess", math.inf)) <= 2.5e-5
        )
        first_birth_distribution = np.asarray(
            moments["first_birth_age_distribution"], dtype=float
        )
        ages = float(P.age_start) + float(P.da) * np.arange(int(P.J))
        share_38plus = float(first_birth_distribution[ages >= 38.0].sum())
        share_42plus = float(first_birth_distribution[ages >= 42.0].sum())
        loss = float(sum(float(row["loss_contribution"]) for row in rows))
        summary_rows.append(
            {
                "case": schedule["case"],
                "description": schedule["description"],
                "strict_converged": strict,
                "market_residual": float(
                    getattr(sol, "best_max_abs_rel_excess", math.inf)
                ),
                "elapsed_sec": elapsed,
                "price": float(price[0]),
                "signed_loss_at_fixed_e5b_theta": loss,
                "tfr": float(moments["tfr"]),
                "childless_rate": float(moments["childless_rate"]),
                "mean_age_first_birth": float(moments["mean_age_first_birth"]),
                "share_first_births_age30plus": float(
                    moments["share_first_births_age30plus"]
                ),
                "share_first_births_age38plus": share_38plus,
                "share_first_births_age42plus": share_42plus,
                "fecundity_omega1": schedule["fecundity_omega1"],
                "fecundity_omega2": schedule["fecundity_omega2"],
                "fecundity_tail_start_age": schedule[
                    "fecundity_tail_start_age"
                ],
                "fecundity_terminal_decay": schedule[
                    "fecundity_terminal_decay"
                ],
                "terminal_probability_before_hard_close": schedule[
                    "terminal_probability_before_hard_close"
                ],
            }
        )
        for row in rows:
            fit_rows.append(
                {
                    "case": schedule["case"],
                    **row,
                    "delta_model_vs_e5b": float(row["model"])
                    - baseline_models[str(row["moment"])],
                }
            )
        from intergen_eqscale_seq_optimized.parameters import get_fecundity_by_age

        schedule_values = np.asarray(get_fecundity_by_age(P), dtype=float)
        for age, probability in zip(ages, schedule_values):
            schedule_rows.append(
                {
                    "case": schedule["case"],
                    "age": float(age),
                    "conception_probability": float(probability),
                }
            )
        print(
            f"DONE {schedule['case']} loss={loss:.6f} "
            f"share38plus={share_38plus:.6f} share42plus={share_42plus:.6f}",
            flush=True,
        )

    write_csv(args.outdir / "frontier_summary.csv", summary_rows)
    write_csv(args.outdir / "target_fit_deltas.csv", fit_rows)
    write_csv(args.outdir / "schedule_by_model_age.csv", schedule_rows)
    metadata = {
        "status": "strict_fixed_e5b_diagnostic_not_recalibration",
        "source": str(args.source),
        "certified_table": str(args.certified),
        "target_count": len(targets),
        "free_parameter_count_in_source": 10,
        "fixed_theta": theta,
        "tight_solver": {"J": 17, "Nb": 120, "max_iter_eq": 40, "tol_eq": 2.5e-5},
        "evidence_anchors": {
            "four_year_conception": {"age30": 0.91, "age35": 0.84, "age40": 0.64},
            "terminal_success_sensitivities": [0.10, 0.05, 0.03],
        },
    }
    (args.outdir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
