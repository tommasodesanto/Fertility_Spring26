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


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.diagnostics import write_diagnostics  # noqa: E402
from intergen_eqscale_seq_optimized.local_panel import jsonable  # noqa: E402

NCHS_COUNTS = ROOT / "code/data/nchs_natality_timing/first_birth_counts_year_age.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--arm",
        choices=("e6a", "e6b", "e6ab", "e6abc"),
        required=True,
    )
    return parser.parse_args()


def configure_environment(arm: str) -> None:
    for name in ("E6A", "E6B", "E6C"):
        os.environ.pop(name, None)
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    if arm in {"e6a", "e6ab", "e6abc"}:
        os.environ["E6A"] = "1"
    if arm in {"e6b", "e6ab", "e6abc"}:
        os.environ["E6B"] = "1"
    if arm == "e6abc":
        os.environ["E6C"] = "1"


def expected_arm(arm: str) -> str:
    return {
        "e6a": "E6A",
        "e6b": "E6B",
        "e6ab": "E6AB",
        "e6abc": "E6ABC",
    }[arm]


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


def first_birth_age_diagnostic(
    moments: dict[str, Any],
    parameters: Any,
    outdir: Path,
) -> dict[str, Any]:
    """Write the supplemental model-versus-NCHS first-birth age comparison."""
    model_distribution = np.asarray(
        moments.get("first_birth_age_distribution", np.array([], dtype=float)),
        dtype=float,
    ).reshape(-1)
    if model_distribution.size != int(parameters.J):
        raise RuntimeError(
            "first-birth age distribution must have one entry per model age"
        )
    ages = float(parameters.age_start) + float(parameters.da) * np.arange(
        int(parameters.J), dtype=float
    )
    model_by_start = {
        int(round(age)): float(share)
        for age, share in zip(ages, model_distribution)
        if 18.0 <= age <= 45.0
    }

    counts: dict[int, float] = {}
    with NCHS_COUNTS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            year, age = int(row["year"]), int(row["age"])
            if 1979 <= year - age <= 1984:
                counts[age] = counts.get(age, 0.0) + float(row["n_first_births"])
    denominator = float(sum(counts.values()))
    if denominator <= 0.0:
        raise RuntimeError("NCHS first-birth cohort counts are empty")

    rows = []
    nchs_by_start: dict[int, float] = {}
    for start in range(18, 43, 4):
        data_share = (
            sum(value for age, value in counts.items() if start <= age < start + 4)
            / denominator
        )
        model_share = model_by_start.get(start, 0.0)
        nchs_by_start[start] = data_share
        rows.append(
            {
                "age_bin_start": start,
                "age_bin_end": start + 3,
                "model_share": model_share,
                "nchs_share": data_share,
                "model_minus_nchs": model_share - data_share,
            }
        )
    write_csv(outdir / "first_birth_age_distribution_4year_bins.csv", rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    x = np.arange(len(rows), dtype=float)
    width = 0.38
    ax.bar(
        x - width / 2,
        [float(row["model_share"]) for row in rows],
        width,
        label="Model",
    )
    ax.bar(
        x + width / 2,
        [float(row["nchs_share"]) for row in rows],
        width,
        label="NCHS 1979--84 cohorts",
    )
    ax.set_xticks(
        x,
        [f"{row['age_bin_start']}--{row['age_bin_end']}" for row in rows],
    )
    ax.set_xlabel("Age at first birth")
    ax.set_ylabel("Share of first births")
    ax.set_title("Supplemental first-birth age distribution")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(
        outdir / "supplemental_first_birth_age_distribution.png",
        dpi=180,
    )
    plt.close(fig)

    def grouped_share(starts: tuple[int, ...], source: dict[int, float]) -> float:
        return float(sum(source.get(start, 0.0) for start in starts))

    early_model = grouped_share((18, 22), model_by_start)
    early_nchs = grouped_share((18, 22), nchs_by_start)
    middle_model = grouped_share((26, 30), model_by_start)
    middle_nchs = grouped_share((26, 30), nchs_by_start)
    late_model = grouped_share((38, 42), model_by_start)
    late_nchs = grouped_share((38, 42), nchs_by_start)
    summary = {
        "early_18_25": {
            "model": early_model,
            "nchs": early_nchs,
            "gap": early_model - early_nchs,
        },
        "middle_26_33": {
            "model": middle_model,
            "nchs": middle_nchs,
            "gap": middle_model - middle_nchs,
        },
        "late_38_45": {
            "model": late_model,
            "nchs": late_nchs,
            "gap": late_model - late_nchs,
        },
        "missing_25_30_shape_confirmed": bool(
            early_model > early_nchs and middle_model < middle_nchs
        ),
    }
    (outdir / "first_birth_age_shape_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True)
    )
    return summary


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
    first_birth_shape = first_birth_age_diagnostic(
        moments,
        parameters,
        args.outdir,
    )
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
        "first_birth_age_shape": first_birth_shape,
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
                "The age-bin table and plot are supplemental and do not alter the",
                "stable standard graph set.",
                "No policy experiment is run.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
