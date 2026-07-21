#!/usr/bin/env python3
"""Fixed-parameter M5 robustness to the owner-housing ladder density."""

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

import numpy as np
import matplotlib.pyplot as plt

from intergen_housing_fertility_optimized.calibration import diagnostic_loss, extract_moments
from intergen_housing_fertility_optimized.m5_profile import M5_THETA, m5_overrides, m5_target_system
from intergen_housing_fertility_optimized.solver import run_model_cp_dt


CASES = (
    ("baseline_2room", np.arange(2.0, 10.0 + 0.01, 2.0)),
    ("dense_1room", np.arange(2.0, 10.0 + 0.01, 1.0)),
    ("dense_halfroom", np.arange(2.0, 10.0 + 0.01, 0.5)),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("../../output/model/intergen_owner_ladder_robustness_20260721"),
    )
    parser.add_argument(
        "--fixed-price",
        action="store_true",
        help="Use the canonical M5 price instead of resolving general equilibrium.",
    )
    parser.add_argument(
        "--custom-points",
        type=int,
        default=None,
        help="Run one evenly spaced ladder on [2,10] with this many points.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def rung_rows(label: str, ladder: np.ndarray, distribution: np.ndarray) -> list[dict[str, Any]]:
    masses = np.sum(distribution, axis=tuple(i for i in range(distribution.ndim) if i != 1))
    owner_total = float(np.sum(masses[1:]))
    rows: list[dict[str, Any]] = []
    for index, rooms in enumerate(ladder, start=1):
        mass = float(masses[index])
        rows.append(
            {
                "case": label,
                "rung": index,
                "rooms": float(rooms),
                "mass": mass,
                "share_among_owners": mass / owner_total if owner_total > 0 else np.nan,
            }
        )
    return rows


def write_summary_artifacts(
    outdir: Path,
    summaries: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
) -> None:
    labels = [row["case"] for row in summaries]
    display_labels = {
        "baseline_2room": "2-room",
        "dense_1room": "1-room",
        "dense_halfroom": "0.5-room",
    }
    moments = {
        (row["case"], row["moment"]): float(row["model"])
        for row in target_rows
    }
    panels = (
        ("Equilibrium price", [float(row["price"]) for row in summaries]),
        ("M5 loss at fixed theta", [float(row["loss_at_fixed_M5_theta"]) for row in summaries]),
        ("Ownership rate", [moments[(label, "own_rate")] for label in labels]),
        ("Ownership, ages 25--34", [moments[(label, "own_rate_2534")] for label in labels]),
        (
            "Childless owners with 6+ rooms",
            [moments[(label, "prime30_55_childless_owner_share_rooms_ge6")] for label in labels],
        ),
        ("Completed fertility", [moments[(label, "tfr")] for label in labels]),
    )
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    colors = ("#304C89", "#648DE5", "#F2B134")
    for ax, (title, values) in zip(axes.flat, panels):
        ax.bar(range(len(labels)), values, color=colors[: len(labels)])
        ax.set_title(title)
        ax.set_xticks(
            range(len(labels)),
            [display_labels.get(label, label.replace("_", "-")) for label in labels],
        )
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Owner-ladder density robustness: M5 parameters held fixed")
    fig.tight_layout()
    fig.savefig(outdir / "owner_ladder_robustness_summary.png", dpi=180)
    plt.close(fig)

    readme = """# Owner-ladder density robustness

This is a fixed-parameter robustness exercise for the current M5 model using
the parity-verified optimized solver. It compares the canonical owner ladder
`[2, 4, 6, 8, 10]` with one-room and half-room ladders over the same support.
Each case resolves the one-market general equilibrium at the unchanged M5
parameter vector and active 15-moment target system.

The exercise is diagnostic, not a recalibration. A deterioration in loss shows
that the calibrated parameter vector is not invariant to the housing grid; it
does not show that a dense-grid model cannot be recalibrated.

Artifacts:

- `summary.csv`: price, loss, residual, convergence, and runtime by ladder.
- `target_fit_long.csv`: complete 15-row target fit for every ladder.
- `owner_rung_mass.csv`: equilibrium owner mass on every rung.
- `owner_ladder_robustness_summary.png`: supplemental visual comparison.
- `metadata.json`: exact M5 theta, solver path, and ladder definitions.
"""
    (outdir / "README.md").write_text(readme)


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    target_system = m5_target_system()
    cases = CASES
    if args.custom_points is not None:
        if int(args.custom_points) < 2:
            raise ValueError("--custom-points must be at least 2")
        cases = ((f"custom_{int(args.custom_points)}point", np.linspace(2.0, 10.0, int(args.custom_points))),)
    metadata = {
        "status": "running",
        "purpose": "fixed-M5-theta owner-ladder density robustness",
        "solver_package": "intergen_housing_fertility_optimized",
        "equilibrium": not args.fixed_price,
        "theta": M5_THETA,
        "cases": {label: ladder.tolist() for label, ladder in cases},
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    summaries: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    all_rung_rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for case_number, (label, ladder) in enumerate(cases, start=1):
        case_start = time.perf_counter()
        print(f"case {case_number}/{len(cases)}: {label}, {len(ladder)} owner rungs", flush=True)
        overrides = m5_overrides(tight=True, optimized=True)
        overrides["H_own"] = ladder.copy()
        if args.fixed_price:
            from intergen_housing_fertility_optimized.m5_profile import M5_PRICE

            overrides.update(
                solve_mode="pe",
                p_fixed=np.array([M5_PRICE]),
                w_fixed=np.array([1.0]),
                entry_shares_fixed=np.array([1.0]),
            )
        solution, parameters, price = run_model_cp_dt(overrides, verbose=False)
        moments = extract_moments(solution, parameters)
        loss = diagnostic_loss(
            moments,
            targets=target_system.targets_dict(),
            weights=target_system.weights_dict(),
        )
        elapsed = time.perf_counter() - case_start
        summary = {
            "case": label,
            "n_owner_rungs": len(ladder),
            "minimum_rooms": float(ladder[0]),
            "maximum_rooms": float(ladder[-1]),
            "rung_spacing": float(ladder[1] - ladder[0]),
            "price": float(np.asarray(price).reshape(-1)[0]),
            "market_residual": float(solution.best_max_abs_rel_excess),
            "loss_at_fixed_M5_theta": float(loss),
            "elapsed_seconds": elapsed,
            "strict_converged": bool(getattr(solution, "timings", {}).get("strict_converged", False)),
        }
        summaries.append(summary)
        for name, target, weight in zip(
            target_system.moment_names, target_system.target_values, target_system.weights
        ):
            model = float(moments[name])
            gap = model - float(target)
            target_rows.append(
                {
                    "case": label,
                    "moment": name,
                    "target": float(target),
                    "model": model,
                    "gap": gap,
                    "weight": float(weight),
                    "loss_contribution": float(weight) * gap * gap,
                }
            )
        all_rung_rows.extend(rung_rows(label, ladder, np.asarray(solution.g)))
        latest = {
            "completed_cases": case_number,
            "total_cases": len(cases),
            "latest": summary,
            "elapsed_seconds": time.perf_counter() - started,
        }
        (outdir / "latest_completed_case.json").write_text(
            json.dumps(jsonable(latest), indent=2, sort_keys=True) + "\n"
        )
        write_csv(outdir / "summary.csv", summaries)
        write_csv(outdir / "target_fit_long.csv", target_rows)
        write_csv(outdir / "owner_rung_mass.csv", all_rung_rows)
        print(
            f"completed {label}: price={summary['price']:.8f}, loss={loss:.6f}, "
            f"residual={summary['market_residual']:.3e}, elapsed={elapsed:.1f}s",
            flush=True,
        )

    metadata["status"] = "complete"
    metadata["elapsed_seconds"] = time.perf_counter() - started
    (outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n")
    write_summary_artifacts(outdir, summaries, target_rows)
    print(f"complete: {outdir}", flush=True)


if __name__ == "__main__":
    main()
