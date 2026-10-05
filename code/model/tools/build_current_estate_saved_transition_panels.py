#!/usr/bin/env python3
"""Render the September 14 transition layouts from saved current-estate arrays only."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True, help="Saved candidate/latest_completed_full.json")
    parser.add_argument("--receipt", type=Path, required=True, help="Matching run/latest_completed.json")
    parser.add_argument("--initial-summary", type=Path, required=True, help="Fresh reference standard_diagnostics/summary.json")
    parser.add_argument("--initial-stationary", type=Path, required=True, help="Fresh reference stationary.json")
    parser.add_argument("--targets", type=Path, help="Optional JSON: birth-window end year to observed fertility")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-incomplete", action="store_true", help="Plot an explicitly labeled incomplete saved path")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    native, receipt = load(args.native), load(args.receipt)
    initial_summary, initial_stationary = load(args.initial_summary), load(args.initial_stationary)
    targets = {int(k): float(v) for k, v in load(args.targets).items()} if args.targets else {}
    rows, fertility = native["rows"], native["fertility"]
    if not len(rows) == len(fertility) == int(receipt["horizon"]):
        raise ValueError("Saved dated row count disagrees with receipt horizon")
    if rows != receipt["full_dated_rows"] or fertility != receipt["fertility_rows"]:
        raise ValueError("Native dated arrays do not match completion receipt")
    root_replay = all(bool(receipt[k]) for k in ("root_pass", "root_gate_passed", "replay_pass"))
    if not root_replay and not args.allow_incomplete:
        raise ValueError("Root/replay gates incomplete; --allow-incomplete required")
    years = np.asarray([int(r["calendar_year"]) for r in rows])
    if not np.array_equal(years, 2007 + 4*np.arange(len(rows))):
        raise ValueError("Expected saved four-year dates beginning in 2007")
    if any(int(r["calendar_year"]) != int(f["calendar_year"]) for r, f in zip(rows, fertility)):
        raise ValueError("Fertility and aggregate years disagree")
    get = lambda k: np.asarray([float(r[k]) for r in rows])
    if not receipt.get("two_shock_result", False) and not np.allclose(get("psi_child"), float(receipt["psi"]), rtol=0, atol=1e-15):
        raise ValueError("First-shock receipt has a varying dated preference")
    period = np.asarray([float(f["period_tfr_topcode_adjusted"]) for f in fertility])
    initial_fertility = float(initial_summary["tfr"])
    initial_population = float(initial_stationary["population_scale"])
    initial_housing = float(initial_stationary["absolute_housing_demand"])
    initial_price = float(initial_stationary["price"])
    summary_housing = float(initial_summary["aggregate_housing_demand"])
    summary_price = float(initial_summary["owner_asset_price"][0])
    if abs(summary_housing*initial_population-initial_housing) > 1e-10 or abs(summary_price-initial_price) > 1e-12:
        raise ValueError("Fresh reference summary and stationary scale disagree")
    housing = 100*get("housing_demand")/initial_housing
    heads = 100*get("adult_population")/initial_population
    prices = 100*get("asset_price")/initial_price
    for label, series in ("period fertility", period), ("household heads", heads), ("housing", housing), ("house price", prices):
        if not np.isfinite(series).all():
            raise ValueError(f"Nonfinite {label} series")
    if min(initial_fertility, initial_population, initial_housing, initial_price) <= 0:
        raise ValueError("Fresh reference normalization must be positive")
    stage = "two-shock" if receipt.get("two_shock_result", False) else "first-shock"
    status = "accepted" if receipt.get("production_ready", False) else "diagnostic"
    note = f"{stage.capitalize()} H{receipt['horizon']} {status}; horizon comparison pending" if receipt.get("pending_horizon_comparison", False) else f"{stage.capitalize()} H{receipt['horizon']} {status}; see receipt for horizon status"
    if not receipt.get("two_shock_result", False):
        note += "; no second shock"
    figure_note = "One permanent fertility-preference decline in 2007." if stage == "first-shock" else "Dated fertility-preference changes."
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 15, "axes.labelsize": 12})
    orange, blue, grey = "#d97815", "#245f99", "#777777"
    # The repeated 2007 date reproduces the vertical change from the fresh initial state.
    def anchored(values: np.ndarray, initial: float, keep: np.ndarray):
        return np.r_[2003, 2007, years[keep]], np.r_[initial, initial, values[keep]]
    def save(fig, name: str):
        for ext in ("png", "pdf"):
            fig.savefig(out / f"{name}.{ext}", dpi=180, facecolor="white")
        plt.close(fig)

    keep = years <= 2063
    if keep.sum() != 15:
        raise ValueError("Saved path must cover 2007–2063 for the September 14 fertility layout")
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    short_x, short_y = anchored(period, initial_fertility, keep)
    model_line, = ax.plot(short_x, short_y, color=orange, lw=2.6, label="Model")
    if targets:
        data_x = np.asarray(sorted(targets), dtype=int)-4
        data_y = np.asarray([targets[int(y+4)] for y in data_x])
        data_line, = ax.plot(data_x, data_y, "s--", color=blue, lw=1.6, ms=5, label="US data")
        np.testing.assert_array_equal(data_line.get_ydata(), data_y)
    ax.axhline(initial_fertility, color=grey, lw=.9, ls=":")
    ax.axvline(2007, color=".45", lw=1, ls=":")
    ax.text(2007, min(initial_fertility+.045, 2.145), "2007 initial steady state", ha="left", va="top", fontsize=10, color=".35")
    ax.set(xlim=(2003, 2065), ylim=(1.55, 2.16), xticks=[2007, 2023, 2043, 2063],
           xlabel="Start of four-year period", ylabel="Period fertility",
           title="Period fertility: first-shock path through 2063" if stage == "first-shock" else "Period fertility: path through 2063")
    ax.grid(axis="y", alpha=.16)
    ax.legend(frameon=False, loc="lower right")
    fig.text(.5, .02, figure_note, ha="center", fontsize=10, color=".3")
    fig.tight_layout(rect=(0, .05, 1, 1))
    np.testing.assert_array_equal(model_line.get_xdata(), short_x)
    np.testing.assert_array_equal(model_line.get_ydata(), short_y)
    save(fig, "fertility_fit_2007_2063")

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 7.8))
    axes = axes.ravel()
    plotted = ((period, initial_fertility, "Period fertility", "Period fertility"),
               (heads, 100.0, "Household heads", "Initial steady state = 100"),
               (housing, 100.0, "Total housing", "Initial steady state = 100"),
               (prices, 100.0, "House price", "Initial steady state = 100"))
    for ax, (values, initial, title, ylabel) in zip(axes, plotted):
        full_x, full_y = anchored(values, initial, np.ones(len(years), dtype=bool))
        artist, = ax.plot(full_x, full_y, color=orange, lw=2.3)
        np.testing.assert_array_equal(artist.get_xdata(), full_x)
        np.testing.assert_array_equal(artist.get_ydata(), full_y)
        ax.axhline(initial, color=grey, lw=.9, ls=":")
        ax.axvline(2007, color=".45", lw=1, ls=":")
        ax.set(title=title, ylabel=ylabel, xlabel="Start of four-year period",
               xlim=(2003, years[-1]+2), xticks=[2007, 2023, 2043, 2063, 2083])
        ax.grid(axis="y", alpha=.16)
        ax.tick_params(axis="x", labelsize=10.5)
    fig.suptitle("Transition after a fertility-preference change", fontsize=17, y=.97)
    fig.text(.5, .02, figure_note+" Terminal steady state not shown.", ha="center", fontsize=10, color=".3")
    fig.subplots_adjust(left=.075, right=.98, top=.88, bottom=.12, wspace=.25, hspace=.42)
    save(fig, "full_transition_four_panel")

    with (out / "standard_transition_series.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["calendar_year", "period_fertility", "household_heads_index", "total_housing_index", "house_price_index"])
        writer.writerows(zip(years, period, heads, housing, prices))
    sources = [args.native, args.receipt, args.initial_summary, args.initial_stationary]
    if args.targets:
        sources.append(args.targets)
    provenance = {
        "classification": note, "candidate": receipt.get("task_id"), "horizon": receipt["horizon"],
        "source_sha256": {str(p): digest(p) for p in sources},
        "initial_reference": {"period_fertility": initial_fertility, "adult_population": initial_population,
                              "absolute_housing_demand": initial_housing, "house_price": initial_price},
        "initial_scale_check": {"summary_housing_per_household": summary_housing,
                                "summary_housing_times_population": summary_housing*initial_population,
                                "summary_house_price": summary_price},
        "saved_path_last_year": int(years[-1]), "no_extrapolation": True,
        "root_replay_pass": root_replay, "production_ready": receipt.get("production_ready", False),
        "pending_horizon_comparison": receipt.get("pending_horizon_comparison", False),
        "two_shock_result": receipt.get("two_shock_result", False), "terminal_line_plotted": False,
        "model_solves": 0, "plots": ["fertility_fit_2007_2063", "full_transition_four_panel"],
    }
    (out / "standard_transition_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
