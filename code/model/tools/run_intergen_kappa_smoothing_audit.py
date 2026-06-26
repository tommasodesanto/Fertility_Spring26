#!/usr/bin/env python3
"""Run a fixed-theta tenure-smoothing audit for the intergen model.

The script varies only ``tenure_choice_kappa`` around a trusted saved cache,
re-solves the model, writes a solution cache for each case, and rebuilds the
standard housing-block audit packet for direct visual comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
import pickle
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = Path(__file__).resolve().parent
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from build_intergen_mechanics_packet import solve_candidate, write_solution_cache  # noqa: E402


DEFAULT_SOURCE_CACHE = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "housing_block_grid_audit_20260626/ge_dense_core_nb120/solution_cache.pkl"
)
DEFAULT_OUTDIR = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "housing_block_kappa_smoothing_audit_20260626"
)
DEFAULT_KAPPAS = "0,0.005,0.01,0.02,0.05"


def main() -> None:
    args = parse_args()
    source_cache = args.source_cache.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    payload = load_cache(source_cache)
    theta = dict(payload["theta"])
    grid = dict(payload["grid"])
    targets = dict(payload["targets"])
    weights = dict(payload["weights"])
    target_set = str(payload.get("target_set", ""))
    entry_b = float(args.entry_b) if args.entry_b is not None else float(theta.get("b_entry_fixed", math.nan))
    source = {
        "theta": theta,
        "source_meta": {
            **dict(payload.get("source_meta", {})),
            "source_cache": str(source_cache),
            "audit": "tenure_choice_kappa_smoothing",
        },
    }

    kappas = parse_float_list(args.kappas)
    summary_rows: list[dict[str, Any]] = []
    moment_rows: list[dict[str, Any]] = []
    curve_rows_by_case: dict[str, list[dict[str, float]]] = {}

    for kappa in kappas:
        label = f"kappa_{slug_float(kappa)}"
        case_dir = outdir / label
        case_dir.mkdir(parents=True, exist_ok=True)
        cache_path = case_dir / "solution_cache.pkl"
        if args.skip_existing and cache_path.exists():
            case_payload = load_cache(cache_path)
            baseline = dict(case_payload["baseline"])
        else:
            print(f"[run] {label}", flush=True)
            baseline = solve_candidate(
                theta=theta,
                grid=grid,
                extra_overrides={"tenure_choice_kappa": float(kappa)},
                targets=targets,
                weights=weights,
                label=label,
            )
            write_solution_cache(
                cache_path,
                baseline=baseline,
                source=source,
                target_set=target_set,
                targets=targets,
                weights=weights,
                grid={**grid, "tenure_choice_kappa": float(kappa)},
            )
        packet_dir = case_dir / "housing_block_audit"
        if not args.skip_packet or not (packet_dir / "README.md").exists():
            run_housing_packet(cache_path, packet_dir)

        curve_rows = read_float_csv(packet_dir / "aggregate_by_current_liquid_wealth.csv")
        curve_rows_by_case[label] = curve_rows
        patho = curve_pathology(curve_rows, mass_min=float(args.mass_min), entry_b=entry_b)
        moments = dict(baseline["moments"])
        summary_rows.append(
            {
                "case": label,
                "tenure_choice_kappa": float(kappa),
                "rank_loss": float(baseline["rank_loss"]),
                "market_residual": float(baseline["market_residual"]),
                "elapsed_sec": float(baseline["elapsed_sec"]),
                "p_eq": float(np.asarray(baseline["p_eq"], dtype=float).reshape(-1)[0]),
                "own_rate": moments.get("own_rate"),
                "own_rate_2534": moments.get("own_rate_2534"),
                "old_age_own_rate": moments.get("old_age_own_rate"),
                "tfr": moments.get("tfr"),
                "housing_increment_0to1": moments.get("housing_increment_0to1"),
                "housing_increment_1to2": moments.get("housing_increment_1to2"),
                "renter_share_rooms_ge6": moments.get("prime30_55_childless_renter_share_rooms_ge6"),
                "owner_share_rooms_ge6": moments.get("prime30_55_childless_owner_share_rooms_ge6"),
                **patho,
            }
        )
        for moment, target in sorted(targets.items()):
            model_value = moments.get(moment, math.nan)
            weight = float(weights.get(moment, 1.0))
            diff = float(model_value) - float(target) if math.isfinite(float(model_value)) else math.nan
            moment_rows.append(
                {
                    "case": label,
                    "tenure_choice_kappa": float(kappa),
                    "moment": moment,
                    "target": float(target),
                    "model": model_value,
                    "weight": weight,
                    "loss_contribution": weight * diff * diff if math.isfinite(diff) else math.nan,
                }
            )
        write_csv(outdir / "summary.csv", summary_rows)
        write_csv(outdir / "target_moments_by_kappa.csv", moment_rows)
        plot_consumption_comparison(curve_rows_by_case, outdir / "current_liquid_consumption_by_kappa.png", float(args.mass_min))
        write_readme(outdir, source_cache, summary_rows, entry_b=entry_b)

    print(outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cache", type=Path, default=DEFAULT_SOURCE_CACHE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--kappas", default=DEFAULT_KAPPAS)
    parser.add_argument("--mass-min", type=float, default=1.0e-8)
    parser.add_argument("--entry-b", type=float, default=None, help="Wealth point for the entry-kink drop diagnostic.")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--skip-packet", action="store_true")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} is not a solution cache")
    return payload


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def slug_float(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def run_housing_packet(cache_path: Path, outdir: Path) -> None:
    cmd = [
        sys.executable,
        str(TOOLS_DIR / "build_intergen_housing_block_packet.py"),
        "--cache",
        str(cache_path),
        "--outdir",
        str(outdir),
    ]
    subprocess.run(cmd, cwd=str(ROOT), check=True)


def read_float_csv(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append({k: float(v) if v not in ("", None) else math.nan for k, v in row.items()})
    return rows


def curve_pathology(rows: list[dict[str, float]], *, mass_min: float, entry_b: float) -> dict[str, float]:
    sel = [r for r in rows if float(r.get("mass", 0.0)) >= mass_min]
    largest_drop = 0.0
    largest_drop_at_b = math.nan
    largest_drop_mass = math.nan
    for prev, cur in zip(sel, sel[1:]):
        dy = float(cur["avg_consumption"] - prev["avg_consumption"])
        if dy < largest_drop:
            largest_drop = dy
            largest_drop_at_b = float(cur["current_liquid_wealth"])
            largest_drop_mass = float(cur["mass"])
    zero_row = min(rows, key=lambda r: abs(float(r["current_liquid_wealth"]))) if rows else {}
    entry_stats = entry_drop(rows, entry_b)
    return {
        "largest_consumption_drop": largest_drop,
        "largest_consumption_drop_at_b": largest_drop_at_b,
        "largest_consumption_drop_mass": largest_drop_mass,
        "mass_at_nearest_zero": float(zero_row.get("mass", math.nan)),
        "nearest_zero_b": float(zero_row.get("current_liquid_wealth", math.nan)),
        **entry_stats,
    }


def entry_drop(rows: list[dict[str, float]], entry_b: float) -> dict[str, float]:
    if not rows or not math.isfinite(entry_b):
        return {
            "entry_b_target": entry_b,
            "entry_b_nearest": math.nan,
            "entry_b_mass": math.nan,
            "entry_b_consumption_drop": math.nan,
            "entry_b_prev_consumption": math.nan,
            "entry_b_consumption": math.nan,
        }
    idx = min(range(len(rows)), key=lambda k: abs(float(rows[k]["current_liquid_wealth"]) - entry_b))
    cur = rows[idx]
    prev = rows[idx - 1] if idx > 0 else None
    prev_c = float(prev["avg_consumption"]) if prev is not None else math.nan
    cur_c = float(cur["avg_consumption"])
    return {
        "entry_b_target": entry_b,
        "entry_b_nearest": float(cur["current_liquid_wealth"]),
        "entry_b_mass": float(cur["mass"]),
        "entry_b_consumption_drop": cur_c - prev_c if math.isfinite(prev_c) else math.nan,
        "entry_b_prev_consumption": prev_c,
        "entry_b_consumption": cur_c,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_consumption_comparison(curves: dict[str, list[dict[str, float]]], path: Path, mass_min: float) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, rows in curves.items():
        sel = [r for r in rows if float(r.get("mass", 0.0)) >= mass_min]
        xs = np.asarray([r["current_liquid_wealth"] for r in sel], dtype=float)
        ys = np.asarray([r["avg_consumption"] for r in sel], dtype=float)
        ax.plot(xs, ys, linewidth=1.8, label=label.replace("kappa_", "k="))
    ax.axvline(0.0, color="0.65", linewidth=1)
    ax.set_title("Average consumption by current liquid wealth, tenure-kappa sweep")
    ax.set_xlabel("current liquid wealth b")
    ax.set_ylabel("average consumption")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_readme(outdir: Path, source_cache: Path, summary_rows: list[dict[str, Any]], *, entry_b: float) -> None:
    with (outdir / "README.md").open("w", encoding="utf-8") as fh:
        fh.write("# Tenure-Kappa Smoothing Audit\n\n")
        fh.write(f"- Source cache: `{source_cache}`\n")
        fh.write("- Grid/theta/targets are inherited from the source cache; only `tenure_choice_kappa` changes.\n")
        fh.write(f"- Entry-kink diagnostic wealth: `{entry_b:.8g}`\n")
        fh.write("- Each case writes `solution_cache.pkl` and a standard `housing_block_audit/` packet.\n\n")
        fh.write("| kappa | loss | residual | own | young own | old own | TFR | entry c drop | entry mass | largest c drop | drop at b |\n")
        fh.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in summary_rows:
            fh.write(
                f"| {float(row['tenure_choice_kappa']):.4g} | {float(row['rank_loss']):.6g} | "
                f"{float(row['market_residual']):.2e} | {float(row['own_rate']):.4f} | "
                f"{float(row['own_rate_2534']):.4f} | {float(row['old_age_own_rate']):.4f} | "
                f"{float(row['tfr']):.4f} | {float(row['entry_b_consumption_drop']):.6g} | "
                f"{float(row['entry_b_mass']):.6g} | {float(row['largest_consumption_drop']):.6g} | "
                f"{float(row['largest_consumption_drop_at_b']):.6g} |\n"
            )
        fh.write("\n## Files\n\n")
        fh.write("- `summary.csv`\n")
        fh.write("- `target_moments_by_kappa.csv`\n")
        fh.write("- `current_liquid_consumption_by_kappa.png`\n")


if __name__ == "__main__":
    main()
