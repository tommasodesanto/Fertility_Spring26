#!/usr/bin/env python3
"""Build the saved-data-only 2007–2063 provisional fertility slide plot."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "output/model/transition_readiness_v1/current_baseline_20261003/local_continuation_v1/fit_job/run"
OUTPUT = ROOT / "output/model/transition_readiness_v1/current_baseline_20261003/monitor/slide_plot"
PATH_JSON = RUN / "candidate_0001/horizon_032/latest_completed.json"
COMPLETE_JSON = RUN / "candidate_0001/complete.json"
PLAN_JSON = RUN / "effective_plan.json"
REFERENCE_JSON = RUN / "native_reference/repeat_0/phase_b_ge/selected_root/standard_diagnostics/summary.json"
STATIONARY_JSON = RUN / "native_reference/repeat_0/stationary.json"
H24_JSON = RUN / "candidate_0001/horizon_024/latest_completed.json"
H24_ROOT = RUN / "candidate_0001/horizon_024/root.json"
H32_ROOT = RUN / "candidate_0001/horizon_032/root.json"
MAP009_JSON = RUN / "candidate_0001/horizon_032/map_009/mapping.json"
ANNUAL_CSV = ROOT / "output/model/e5f_matched_pf_20260909a/path_pilot_20260910/fertility_data/annual_fertility_2007_2023.csv"
REFERENCE_PNG = ROOT / "output/model/e5f_original_queue_20260913a/terminal_restart_v1/fertility_replay_iter3/output/fertility_fit_2007_2063.png"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    paths = (PATH_JSON, COMPLETE_JSON, PLAN_JSON, REFERENCE_JSON, STATIONARY_JSON, H24_JSON, H24_ROOT,
             H32_ROOT, MAP009_JSON, ANNUAL_CSV, REFERENCE_PNG)
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    saved = json.loads(PATH_JSON.read_text())
    completion = json.loads(COMPLETE_JSON.read_text())
    plan = json.loads(PLAN_JSON.read_text())
    reference = json.loads(REFERENCE_JSON.read_text())
    stationary = json.loads(STATIONARY_JSON.read_text())
    h24 = json.loads(H24_JSON.read_text())
    root24 = json.loads(H24_ROOT.read_text())
    root32 = json.loads(H32_ROOT.read_text())
    map009 = json.loads(MAP009_JSON.read_text())

    fert = saved["fertility"]
    years = np.asarray([int(row["calendar_year"]) for row in fert], dtype=int)
    model = np.asarray([float(row["period_tfr_topcode_adjusted"]) for row in fert])
    keep = years <= 2063
    years, model = years[keep], model[keep]
    if len(years) != 15 or years[0] != 2007 or years[-1] != 2063 or np.any(np.diff(years) != 4):
        raise ValueError("Saved fertility path must cover 2007–2063 in four-year steps")
    np.testing.assert_array_equal(model, [float(r["period_tfr_topcode_adjusted"]) for r in map009["fertility"][:15]])
    np.testing.assert_array_equal(
        completion["payload"]["models"],
        [float(r["period_tfr_topcode_adjusted"]) for r in fert[:4]],
    )
    for root, result in ((root24, h24), (root32, saved)):
        if not root["converged"] or not all(root["gates"].values()) or not all(result["gates"].values()):
            raise ValueError("24/32-horizon saved root/replay acceptance gates are incomplete")

    initial = float(reference["tfr"])
    if not np.isfinite(initial) or initial <= 0:
        raise ValueError("Saved native-reference TFR is invalid")
    psi = float(completion["psi"])
    original_psi = float(stationary["psi_child"])
    psi_decline = 100 * (1 - psi / original_psi)

    annual: dict[int, float] = {}
    with ANNUAL_CSV.open(newline="") as fh:
        for row in csv.DictReader(fh):
            annual[int(row["year"])] = float(row["period_tfr_births_per_woman"])
    target_rows = plan["target_contract"]["rows"]
    data_years, data_values = [], []
    for row in target_rows:
        y0, y1 = int(row["birth_year_start"]), int(row["birth_year_end"])
        vals = [annual[y] for y in range(y0, y1 + 1)]
        value = float(np.mean(vals))
        if abs(value - float(row["target"])) > 1e-12:
            raise ValueError(f"Annual source and target contract disagree for {y0}–{y1}")
        data_years.append(int(row["decision_year"]))
        data_values.append(value)
    data_years = np.asarray(data_years, dtype=int)
    data_values = np.asarray(data_values, dtype=float)

    # Repeated 2007 date draws the shock impact as a vertical change from initial SS.
    model_x = np.concatenate(([2003, 2007], years))
    model_y = np.concatenate(([initial, initial], model))
    if not (np.isfinite(model_x).all() and np.isfinite(model_y).all() and np.isfinite(data_values).all()):
        raise ValueError("Plot coordinates contain nonfinite values")
    if model_x.max() > 2063:
        raise ValueError("Refusing to extrapolate beyond the saved 2063 horizon")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 15, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    if ax.spines["top"].get_visible() or ax.spines["right"].get_visible():
        raise AssertionError("Expected top and right spines to be hidden")
    orange, blue, grey = "#d97815", "#245f99", ".5"
    model_line, = ax.plot(model_x, model_y, color=orange, lw=2.6, label="Model")
    data_line, = ax.plot(data_years, data_values, "s--", color=blue, lw=1.6, ms=5, label="US data")
    ax.axhline(initial, color=grey, lw=.9, ls=":")
    ax.axvline(2007, color=".45", lw=1, ls=":")
    ax.text(2007, initial + .015, "2007 initial steady state", ha="left", va="bottom", fontsize=10, color=".35")
    ax.set(xlim=(2003, 2065), ylim=(1.55, 2.16), xticks=[2007, 2023, 2043, 2063],
           xlabel="Start of four-year period", ylabel="Period fertility", title="Period fertility: 2007–2063")
    ax.grid(axis="y", alpha=.16)
    ax.legend(frameon=False)
    fig.text(.5, .02, f"Preference falls {psi_decline:.1f}% in 2007 and remains fixed.", ha="center", fontsize=10, color=".3")
    fig.tight_layout(rect=(0, .05, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(OUTPUT / f"fertility_2007_2063.{ext}", dpi=180, facecolor="white")

    # Verify the matplotlib artist coordinates against the exact plotted inputs.
    artist_model_x, artist_model_y = model_line.get_data()
    artist_data_x, artist_data_y = data_line.get_data()
    artist_check = bool(np.array_equal(artist_model_x, model_x) and np.array_equal(artist_model_y, model_y)
                        and np.array_equal(artist_data_x, data_years) and np.array_equal(artist_data_y, data_values))
    plt.close(fig)
    if not artist_check:
        raise AssertionError("Rendered line coordinates differ from inputs")

    with (OUTPUT / "fertility_2007_2063_plotted_values.csv").open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["series", "calendar_year", "period_fertility"])
        for x, y in zip(model_x, model_y):
            writer.writerow(["model", int(x), format(float(y), ".17g")])
        for x, y in zip(data_years, data_values):
            writer.writerow(["US data", int(x), format(float(y), ".17g")])

    target_indices = {int(y): i for i, y in enumerate(years)}
    h24_values = {int(r["calendar_year"]): float(r["period_tfr_topcode_adjusted"]) for r in h24["fertility"]}
    early_differences = [{"calendar_year": int(y), "horizon_024": h24_values[int(y)], "horizon_032": float(fert[i]["period_tfr_topcode_adjusted"]),
                         "absolute_difference": abs(h24_values[int(y)]-float(fert[i]["period_tfr_topcode_adjusted"]))}
                        for i, y in enumerate(years[:4])]
    horizon_max_2063 = max(abs(h24_values[int(row["calendar_year"])]-float(row["period_tfr_topcode_adjusted"]))
                           for row in fert if int(row["calendar_year"]) <= 2063)
    final_target_row = target_rows[-1]
    final_model = float(fert[target_indices[int(final_target_row["decision_year"])]]["period_tfr_topcode_adjusted"])
    final_target = float(final_target_row["target"])
    final_gap = final_model - final_target
    original_tolerance = float(plan["fit"]["fertility_tolerance"])
    source_paths = [PATH_JSON, COMPLETE_JSON, PLAN_JSON, REFERENCE_JSON, STATIONARY_JSON, H24_JSON,
                    H24_ROOT, H32_ROOT, MAP009_JSON, ANNUAL_CSV, REFERENCE_PNG]
    verification = {
        "classification": "provisional saved-data-only one-shock finite-horizon diagnostic projection",
        "production_ready": bool(completion.get("production_ready", False)),
        "full_path_certified": bool(completion.get("full_path_certified", False)),
        "terminal": False,
        "horizon_years_plotted": [int(years[0]), int(years[-1])],
        "no_extrapolation": bool(model_x.max() == 2063),
        "initial_steady_state_fertility": initial,
        "initial_fertility_source": str(REFERENCE_JSON.relative_to(ROOT)),
        "preference_decline_percent": 100*(1-float(completion["psi"])/float(stationary["psi_child"])),
        "plot_bounds": {"xlim": [2003, 2065], "ylim": [1.55, 2.16], "xticks": [2007, 2023, 2043, 2063], "figsize_inches": [11.5, 6.5]},
        "psi": {"value": psi, "original_value": original_psi, "lower_bound": original_psi * float(plan["psi_bound_ratios"][0]),
                "upper_bound": original_psi * float(plan["psi_bound_ratios"][1]), "decline_percent": psi_decline},
        "horizon_root_checks": {"horizon_024": {"converged": root24["converged"], "gates": root24["gates"], "latest_gates": h24["gates"]},
                                "horizon_032": {"converged": root32["converged"], "gates": root32["gates"], "latest_gates": saved["gates"]}},
        "horizon_032_matches_map_009": True,
        "candidate_first_four_match_completion_payload": True,
        "first_four_horizon_comparison": early_differences,
        "first_four_horizon_comparison_max_abs": max(r["absolute_difference"] for r in early_differences),
        "first_four_horizon_comparison_within_0_001": all(r["absolute_difference"] <= .001 for r in early_differences),
        "historical_horizon_comparison_passed": bool(completion["horizon_comparison"]["passed"]),
        "final_target_comparison": {"decision_year": int(final_target_row["decision_year"]), "target": final_target, "model": final_model,
                                    "gap": final_gap, "original_tolerance": original_tolerance,
                                    "within_original_tolerance": abs(final_gap) <= original_tolerance,
                                    "loss_contribution": float(completion["loss_contribution"])},
        "final_fit_accepted": bool(completion.get("shock_fit_complete", False)),
        "projection_max_24_32_difference_through_2063": horizon_max_2063,
        "projection_stability_check_0_001_passed": horizon_max_2063 <= .001,
        "projection_interpretation": "Projection beyond 2023 remains diagnostic; the historical acceptance gate passed.",
        "artist_coordinates_match_inputs": artist_check,
        "source_hashes_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in source_paths},
        "plotted_model": [{"calendar_year": int(x), "period_fertility": float(y)} for x, y in zip(years, model)],
        "plotted_data": [{"decision_year": int(x), "period_fertility": float(y)} for x, y in zip(data_years, data_values)],
        "reproduction_command": f"{ROOT}/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python {ROOT}/code/model/tools/build_current_estate_fertility_plot.py",
    }
    (OUTPUT / "verification.json").write_text(json.dumps(verification, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
