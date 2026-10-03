"""Plot population aggregates from the latest completed saved model run.

Open/run ``code/model/run_model.py`` first. This file only loads its validated
saved result; it never initializes reference parameters or solves the model.
Edit the settings below, then run this file to write plots and aggregate tables
inside the selected run's ``aggregate_plots`` directory.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


# Set a completed run directory to pin a particular run. None uses the loader's
# latest fully validated completed run.
RUN_DIRECTORY: str | Path | None = None
SHOW_PLOTS = False
WEALTH_RANGE = "central"  # "central" or "all"

_MODEL_DIR = Path(__file__).resolve().parent
_TOOLS_DIR = _MODEL_DIR / "tools"


def main() -> Path:
    """Load the completed saved run, write seven figures and aggregate data."""
    import sys

    import numpy as np
    import matplotlib.pyplot as plt

    if str(_TOOLS_DIR) not in sys.path:
        sys.path.insert(0, str(_TOOLS_DIR))
    from model_policy_tools import aggregate_solution
    from model_run_io import load_run

    try:
        result, run_directory = load_run(RUN_DIRECTORY)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "No completed saved model run is available. Open and run "
            "code/model/run_model.py, then rerun this plotting script."
        ) from exc

    sol = result.solution
    P = result.P
    houses = np.asarray(P.H_own, dtype=float).reshape(-1)
    age_start = int(P.age_start)
    period_years = int(getattr(P, "period_years", getattr(P, "da", 4)))
    aggregates = aggregate_solution(
        sol, houses=houses, age_start=age_start, period_years=period_years,
    )
    ages = np.asarray([row["age"] for row in aggregates["by_age"]], dtype=int)
    by_age = aggregates["by_age"]
    output_directory = Path(run_directory) / "aggregate_plots"
    output_directory.mkdir(parents=True, exist_ok=True)

    consumption = np.asarray([row["mean_consumption"] for row in by_age], dtype=float)
    next_assets = np.asarray([row["mean_next_assets"] for row in by_age], dtype=float)
    inherited_assets = np.asarray([row["mean_inherited_assets"] for row in by_age], dtype=float)
    asset_change = np.asarray([row["mean_asset_change"] for row in by_age], dtype=float)
    rooms = np.asarray([row["mean_rooms"] for row in by_age], dtype=float)
    ownership = np.asarray([row["ownership_rate"] for row in by_age], dtype=float)

    consumption_fig, ax = plt.subplots()
    ax.plot(ages, consumption, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel(
        f"Mean consumption per {period_years}-year period\n"
        "(mean annual gross-earnings units)"
    )
    ax.set_title("Mean consumption by age")
    consumption_fig.tight_layout()
    consumption_fig.savefig(output_directory / "consumption_by_age.png", dpi=160)

    next_assets_fig, ax = plt.subplots()
    ax.plot(ages, next_assets, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel("Mean next-period financial assets, b'\n(mean annual gross-earnings units)")
    ax.set_title("Mean next-period financial assets by age")
    next_assets_fig.tight_layout()
    next_assets_fig.savefig(output_directory / "next_assets_by_age.png", dpi=160)

    inherited_fig, ax = plt.subplots()
    ax.plot(ages, inherited_assets, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel("Mean inherited financial assets, b\n(mean annual gross-earnings units)")
    ax.set_title("Mean inherited financial assets by age")
    inherited_fig.tight_layout()
    inherited_fig.savefig(output_directory / "inherited_assets_by_age.png", dpi=160)

    change_fig, ax = plt.subplots()
    ax.plot(ages, asset_change, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel("Mean change: b' − inherited b\n(mean annual gross-earnings units)")
    ax.set_title("Mean financial asset change by age")
    ax.text(
        0.02, 0.98,
        "Includes housing transaction cash flows; not national-account saving.",
        transform=ax.transAxes, va="top", fontsize=8,
    )
    change_fig.tight_layout()
    change_fig.savefig(output_directory / "financial_asset_change_by_age.png", dpi=160)

    rooms_fig, ax = plt.subplots()
    ax.plot(ages, rooms, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel("Mean housing services (rooms)")
    ax.set_title("Mean housing services by age")
    rooms_fig.tight_layout()
    rooms_fig.savefig(output_directory / "rooms_by_age.png", dpi=160)

    ownership_fig, ax = plt.subplots()
    ax.plot(ages, ownership, marker=".")
    ax.set_xlabel("Age")
    ax.set_ylabel("Ownership rate")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Ownership by age")
    ownership_fig.tight_layout()
    ownership_fig.savefig(output_directory / "ownership_by_age.png", dpi=160)

    distribution = aggregates["inherited_asset_distribution"]
    asset_nodes = np.asarray(distribution["asset_nodes"], dtype=float)
    pooled_mass = np.asarray(distribution["mass"], dtype=float)
    if WEALTH_RANGE not in {"central", "all"}:
        raise ValueError("WEALTH_RANGE must be 'central' or 'all'")
    if WEALTH_RANGE == "central":
        cdf = np.cumsum(pooled_mass / pooled_mass.sum())
        first = max(0, int(np.searchsorted(cdf, 0.0001)) - 1)
        last = min(asset_nodes.size - 1, int(np.searchsorted(cdf, 0.995)) + 1)
        visible = (asset_nodes >= asset_nodes[first]) & (asset_nodes <= asset_nodes[last])
    else:
        visible = np.ones(asset_nodes.size, dtype=bool)
    wealth_fig, ax = plt.subplots()
    ax.plot(asset_nodes[visible], pooled_mass[visible], marker=".")
    ax.set_xlabel("Inherited financial assets b\n(mean annual gross-earnings units)")
    ax.set_ylabel("Population mass at node")
    ax.set_title(f"Pooled inherited-asset distribution ({WEALTH_RANGE} range)")
    if WEALTH_RANGE == "central":
        ax.set_xlim(asset_nodes[first], asset_nodes[last])
    wealth_fig.tight_layout()
    wealth_fig.savefig(output_directory / "inherited_asset_distribution.png", dpi=160)

    with (output_directory / "aggregates_by_age.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(by_age[0].keys()))
        writer.writeheader()
        writer.writerows(by_age)
    metadata = {
        "run_directory": str(Path(run_directory).resolve()),
        "run_label": str(getattr(result, "label", "")),
        "price": float(result.price),
        "wealth_range": WEALTH_RANGE,
        "units": aggregates["units"],
        "population_mass": aggregates["population_mass"],
        "overall": aggregates["overall"],
        "by_age": by_age,
    }
    (output_directory / "aggregates.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    )

    figures = [consumption_fig, next_assets_fig, inherited_fig, change_fig,
               rooms_fig, ownership_fig, wealth_fig]
    if SHOW_PLOTS:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)
    print(f"Wrote 7 plots and aggregate tables to {output_directory}")
    return output_directory


if __name__ == "__main__":
    main()
