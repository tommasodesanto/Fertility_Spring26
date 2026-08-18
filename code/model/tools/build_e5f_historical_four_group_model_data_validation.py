#!/usr/bin/env python3
"""Build historical model-data validation for the selected four-group candidate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_HISTORICAL_DATA = (
    ROOT
    / "output/model/e5f_no_policy_transition_report_"
    "jump11_polish_r2_20260818/historical_model_data.csv"
)

SERIES = (
    (
        "household_population_index",
        "Adult-household mass",
        "adult_household_index_2007",
        "Index (2007 = 1)",
    ),
    (
        "period_tfr_explicit_states",
        "Period fertility",
        "period_tfr_explicit_states",
        "Births per model exposure",
    ),
    (
        "period_first_birth_mean_age",
        "Mean age at first birth",
        "period_first_birth_mean_age",
        "Years",
    ),
    (
        "period_first_birth_share_age30plus",
        "First births at age 30+",
        "period_first_birth_share_age30plus",
        "Share",
    ),
    (
        "aggregate_ownership_rate_18_85",
        "Ownership",
        "owner_rate",
        "Share",
    ),
    (
        "mean_rooms_literal_18_85",
        "Mean occupied rooms",
        "mean_rooms",
        "Rooms",
    ),
    (
        "real_house_price_index",
        "Real house price",
        "asset_price_index_2007",
        "Index (2007 = 1)",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinate-report", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument(
        "--historical-data", type=Path, default=DEFAULT_HISTORICAL_DATA
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def comparison_role(series_id: str, year: int) -> str:
    if series_id == "household_population_index":
        if year == 2007:
            return "initial_condition"
        if year < 2023:
            return "demographic_closure_training_fit"
        return "demographic_closure_heldout_2023"
    if series_id == "period_tfr_explicit_states":
        return "directional_holdout_different_exposure_unit"
    if series_id == "mean_rooms_literal_18_85" and year == 2007:
        return "holdout_with_2007_rooms_topcode_break"
    return "untargeted_historical_holdout"


def fit_statistics(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for series_id in sorted({str(row["series_id"]) for row in rows}):
        subset = sorted(
            [row for row in rows if row["series_id"] == series_id],
            key=lambda row: int(row["calendar_year"]),
        )
        data = np.asarray([float(row["data_value"]) for row in subset])
        model = np.asarray([float(row["model_value"]) for row in subset])
        signs = []
        for data_change, model_change in zip(
            np.diff(data), np.diff(model), strict=True
        ):
            signs.append(
                int(
                    (data_change == 0.0 and model_change == 0.0)
                    or data_change * model_change > 0.0
                )
            )
        output.append(
            {
                "series_id": series_id,
                "series_label": subset[0]["series_label"],
                "rmse_levels": float(np.sqrt(np.mean((model - data) ** 2))),
                "gap_2007": float(model[0] - data[0]),
                "gap_2023": float(model[-1] - data[-1]),
                "change_gap_2007_2023": float(
                    (model[-1] - model[0]) - (data[-1] - data[0])
                ),
                "interval_direction_matches": int(sum(signs)),
                "interval_count": len(signs),
            }
        )
    return output


def make_figures(rows: list[dict[str, Any]], outdir: Path) -> None:
    plotted = [row for row in SERIES if row[0] != "mean_rooms_literal_18_85"]
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.1))
    for ax, (series_id, label, _, units) in zip(axes.flat, plotted, strict=True):
        subset = sorted(
            [row for row in rows if row["series_id"] == series_id],
            key=lambda row: int(row["calendar_year"]),
        )
        years = [int(row["calendar_year"]) for row in subset]
        ax.plot(
            years,
            [float(row["data_value"]) for row in subset],
            color="#222222",
            marker="o",
            lw=2.0,
            label="U.S. data",
        )
        ax.plot(
            years,
            [float(row["model_value"]) for row in subset],
            color="#26547c",
            marker="o",
            lw=2.0,
            label="Model",
        )
        ax.set(title=label, xlabel="Year", ylabel=units)
        ax.grid(alpha=0.18)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "historical_model_data_validation.png", dpi=220)
    fig.savefig(outdir / "historical_model_data_validation.pdf")
    plt.close(fig)

    subset = sorted(
        [row for row in rows if row["series_id"] == "mean_rooms_literal_18_85"],
        key=lambda row: int(row["calendar_year"]),
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    years = [int(row["calendar_year"]) for row in subset]
    ax.plot(
        years[1:],
        [float(row["data_value"]) for row in subset[1:]],
        color="#222222",
        marker="o",
        lw=2.0,
        label="U.S. data (2011--2023)",
    )
    ax.plot(
        years,
        [float(row["model_value"]) for row in subset],
        color="#26547c",
        marker="o",
        lw=2.0,
        label="Model",
    )
    ax.set(title="Mean occupied rooms", xlabel="Year", ylabel="Rooms")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "historical_rooms_validation.png", dpi=220)
    fig.savefig(outdir / "historical_rooms_validation.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    report_path = args.coordinate_report.resolve()
    results_dir = args.results_dir.resolve()
    historical_path = args.historical_data.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("status") != "complete_historical_four_group_coordinate_collection":
        raise RuntimeError("Coordinate report is incomplete")
    best_task = int(report["best"]["task_id"])
    model_path = results_dir / f"task_{best_task:03d}" / "transition_path.csv"
    if not model_path.is_file() or not historical_path.is_file():
        raise FileNotFoundError("Selected transition or historical data is missing")
    model_rows = read_csv(model_path)
    if [int(float(row["year"])) for row in model_rows] != [2007, 2011, 2015, 2019, 2023]:
        raise RuntimeError("Selected model path has the wrong dates")
    model_by_year = {int(float(row["year"])): row for row in model_rows}
    historical_rows = read_csv(historical_path)
    data_lookup = {
        (row["series_id"], int(float(row["calendar_year"]))): row
        for row in historical_rows
    }
    output: list[dict[str, Any]] = []
    for series_id, label, model_field, units in SERIES:
        for year in (2007, 2011, 2015, 2019, 2023):
            source = data_lookup.get((series_id, year))
            if source is None or source["data_value"] == "":
                raise RuntimeError(f"Historical data missing for {series_id}, {year}")
            model_value = float(model_by_year[year][model_field])
            if not math.isfinite(model_value):
                raise RuntimeError(f"Model value is nonfinite for {series_id}, {year}")
            output.append(
                {
                    "series_id": series_id,
                    "series_label": label,
                    "calendar_year": year,
                    "data_value": float(source["data_value"]),
                    "model_value": model_value,
                    "gap": model_value - float(source["data_value"]),
                    "units": units,
                    "comparison_role": comparison_role(series_id, year),
                    "data_source": source["data_source"],
                    "comparability": source["comparability"],
                }
            )
    write_csv(outdir / "historical_model_data.csv", output)
    statistics = fit_statistics(output)
    write_csv(outdir / "historical_fit_statistics.csv", statistics)
    make_figures(output, outdir)
    summary = {
        "status": "complete_historical_four_group_model_data_validation",
        "selected_task_id": best_task,
        "selected_candidate_id": report["best"]["candidate_id"],
        "selected_transition_loss": report["best"]["loss"],
        "historical_roles": {
            "population_2007": "initial condition",
            "population_2011_2019": "used to estimate four demographic rates",
            "population_2023": "held out from demographic-rate estimation",
            "other_series": "not used to estimate the four demographic rates",
        },
        "fit_statistics": statistics,
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (outdir / "README.md").write_text(
        "# Historical validation under the four-group closure\n\n"
        "The demographic closure has four time-invariant age-block rates. They are "
        "estimated on 2007--2019 household totals and age distributions; 2023 is "
        "held out. Household and preference parameters are selected from the "
        "separate 2023 twelve-moment objective. Ownership, rooms, period fertility, "
        "period first-birth timing, and real house prices are shown as historical "
        "validation series. The model period-fertility denominator is adult "
        "households, not national female exposure. The 2007 ACS rooms observation "
        "is omitted from the rooms graph because of its 9+ top-code break.\n",
        encoding="utf-8",
    )
    print(
        "HISTORICAL_FOUR_GROUP_MODEL_DATA_VALIDATION_COMPLETE "
        f"candidate={report['best']['candidate_id']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
