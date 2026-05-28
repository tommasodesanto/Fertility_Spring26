"""Plot model-vs-data ownership lifecycle diagnostics.

The model input is the saved diagnostic ownership-by-age table. The data input
is the ACS/MMS household-head ownership audit, which is the target-consistent
counterpart for model agents.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = Path(__file__).resolve().parent
for _path in (MODEL_ROOT, TOOLS_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.parameters import asdict  # noqa: E402
from dt_cp_model.solver import run_model_cp_dt  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402
from dt_cp_model.utils import make_grid  # noqa: E402
from export_distributional_discipline import current_children  # noqa: E402

DEFAULT_CASE = "renterhybrid_a93h52_m010"
DEFAULT_CASE_DIR = ROOT / "output/model/three_fit_diagnostics_20260528" / DEFAULT_CASE
DEFAULT_RECORD = ROOT / "output/model/three_fit_diagnostics_20260528/records/renterhybrid_a93h52_m010_best.json"
DEFAULT_MODEL_CSV = DEFAULT_CASE_DIR / "tables/ownership_by_age.csv"
DEFAULT_ACS_AGE_CSV = ROOT / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"
DEFAULT_ACS_LOC_CSV = ROOT / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_location_profiles.csv"
DEFAULT_MMS_AGE_CSV = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles_full.csv"
DEFAULT_OUTDIR = DEFAULT_CASE_DIR / "figures"
ACS_SAMPLE = "household_heads_hhwt_due_housing"

COLORS = {
    "model": "#1f5fa6",
    "data": "#a33a3a",
    "gap": "#3f6f43",
    "grid": "#c8ccd2",
    "band": "#b6bcc5",
}

LOC_LABEL = {"center": "Center", "periphery": "Periphery"}
LOC_INDEX = {"periphery": 0, "center": 1}


def read_model(path: Path) -> dict[str, dict[int, float]]:
    out: dict[str, dict[int, float]] = defaultdict(dict)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["current_child_bin"] != "all":
                continue
            loc = row["location"]
            if loc not in {"all", "center", "periphery"}:
                continue
            out[loc][int(float(row["age"]))] = float(row["own_rate"])
    return dict(out)


def read_acs_overall(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("source") != "ACS" or row.get("sample") != ACS_SAMPLE:
                continue
            if row.get("owner_rate", "") == "":
                continue
            out[int(float(row["age"]))] = float(row["owner_rate"])
    return out


def read_acs_location(path: Path) -> dict[str, dict[int, float]]:
    out: dict[str, dict[int, float]] = defaultdict(dict)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run: Rscript code/data/mms_center_periphery/audit_ownership_targets.R"
        )
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("source") != "ACS" or row.get("sample") != ACS_SAMPLE:
                continue
            loc = row.get("mms_location")
            if loc not in {"center", "periphery"}:
                continue
            if row.get("owner_rate", "") == "":
                continue
            out[loc][int(float(row["age"]))] = float(row["owner_rate"])
    return dict(out)


def read_mms_children_location(path: Path) -> dict[str, dict[int, float]]:
    out: dict[str, dict[int, float]] = defaultdict(dict)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            loc = row.get("mms_location")
            if loc not in {"center", "periphery"}:
                continue
            value = row.get("has_child_u18_rate") or row.get("has_children_rate")
            if value in (None, ""):
                continue
            out[loc][int(float(row["age"]))] = float(value)
    return dict(out)


def series_xy(values: dict[int, float], age_min: int, age_max: int) -> tuple[np.ndarray, np.ndarray]:
    ages = np.array(sorted(age for age in values if age_min <= age <= age_max), dtype=float)
    vals = np.array([values[int(age)] for age in ages], dtype=float)
    return ages, vals


def location_gap(values: dict[str, dict[int, float]], age_min: int, age_max: int) -> dict[int, float]:
    center = values.get("center", {})
    periphery = values.get("periphery", {})
    ages = sorted(set(center).intersection(periphery))
    return {age: periphery[age] - center[age] for age in ages if age_min <= age <= age_max}


def setup_axis(ax, title: str, age_min: int, age_max: int, ylabel: str = "Ownership rate (%)") -> None:
    ax.set_title(title)
    ax.set_xlim(age_min, age_max)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Age")
    ax.set_ylabel(ylabel)
    ax.grid(True, color=COLORS["grid"], alpha=0.55, linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def build_setup(record: dict):
    return build_direct_calibration_setup(
        "benchmark",
        geo_weight=100.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=record.get("hR_max"),
        alpha_cons=record.get("alpha_cons"),
        owner_h_bar_scale=record.get("owner_h_bar_scale"),
        owner_size_cost=record.get("owner_size_cost"),
        owner_size_cost_ref=record.get("owner_size_cost_ref"),
        owner_size_cost_power=record.get("owner_size_cost_power"),
        tenure_choice_kappa=record.get("tenure_choice_kappa"),
        weight_overrides=record.get("weight_overrides") or None,
        extra_targets=record.get("extra_targets") or None,
        H_own=record.get("H_own"),
    )


def solve_record(record: dict, setup):
    theta = np.asarray(record["theta"], dtype=float)
    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    P = apply_theta(P, [theta_dict[name] for name in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if OUTSIDE_VALUE_NAME in theta_dict:
        P.outside_value = theta_dict[OUTSIDE_VALUE_NAME]
        P.outside_value_is_calibrated = True
        P.allow_uncalibrated_outside_value = False
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    return run_model_cp_dt(P, verbose=False)


def compute_model_spatial_lifecycle(record_path: Path) -> dict[str, dict[str, dict[int, float]]]:
    record = json.loads(record_path.read_text())
    sol, P, _ = solve_record(record, build_setup(record))
    b_grid = make_grid(P)
    out: dict[str, dict[str, dict[int, float]]] = {
        "ownership": defaultdict(dict),
        "networth": defaultdict(dict),
        "children": defaultdict(dict),
    }

    for loc_name, loc in LOC_INDEX.items():
        for jj in range(P.J):
            age = int(P.age_start + jj)
            mass = float(np.sum(sol.g[:, :, loc, jj, :, :]))
            own_mass = float(np.sum(sol.g[:, 1:, loc, jj, :, :]))
            child_mass = 0.0
            wealth_sum = 0.0
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    if current_children(nn, cs, P) > 0:
                        child_mass += float(np.sum(sol.g[:, :, loc, jj, nn, cs]))
                    renter = sol.g[:, 0, loc, jj, nn, cs]
                    wealth_sum += float(np.sum(renter * b_grid))
                    for ten in range(1, 1 + P.n_house):
                        owner = sol.g[:, ten, loc, jj, nn, cs]
                        if float(np.sum(owner)) <= 1e-15:
                            continue
                        house_value = float(sol.p_eq[loc]) * float(P.H_own[ten - 1])
                        wealth_sum += float(np.sum(owner * (b_grid + house_value)))
            out["ownership"][loc_name][age] = own_mass / max(mass, 1e-12)
            out["children"][loc_name][age] = child_mass / max(mass, 1e-12)
            out["networth"][loc_name][age] = wealth_sum / max(mass, 1e-12)
    return {key: dict(value) for key, value in out.items()}


def plot_overall(model: dict[str, dict[int, float]], acs: dict[int, float], outbase: Path, age_min: int, age_max: int) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.axvspan(30, 55, color=COLORS["band"], alpha=0.16, linewidth=0, label="Target ages 30-55")

    ages, vals = series_xy(model["all"], age_min, age_max)
    ax.plot(ages, 100 * vals, color=COLORS["model"], linewidth=2.7, label="Model")

    ages, vals = series_xy(acs, age_min, age_max)
    ax.plot(ages, 100 * vals, color=COLORS["data"], linewidth=2.4, linestyle="--", label="ACS heads")

    setup_axis(ax, "Ownership Lifecycle: Model vs ACS Household Heads", age_min, age_max)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_spatial(
    model: dict[str, dict[int, float]],
    acs: dict[str, dict[int, float]],
    outbase: Path,
    age_min: int,
    age_max: int,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.5), sharex=False)
    panels = [("center", "Center"), ("periphery", "Periphery")]
    for ax, (loc, title) in zip(axes[:2], panels):
        ax.axvspan(30, 55, color=COLORS["band"], alpha=0.14, linewidth=0)
        ages, vals = series_xy(model[loc], age_min, age_max)
        ax.plot(ages, 100 * vals, color=COLORS["model"], linewidth=2.5, label="Model")
        ages, vals = series_xy(acs[loc], age_min, age_max)
        ax.plot(ages, 100 * vals, color=COLORS["data"], linewidth=2.2, linestyle="--", label="ACS heads")
        setup_axis(ax, title, age_min, age_max)

    ax = axes[2]
    m_gap = location_gap(model, age_min, age_max)
    d_gap = location_gap(acs, age_min, age_max)
    ages, vals = series_xy(m_gap, age_min, age_max)
    ax.plot(ages, 100 * vals, color=COLORS["model"], linewidth=2.5, label="Model")
    ages, vals = series_xy(d_gap, age_min, age_max)
    ax.plot(ages, 100 * vals, color=COLORS["data"], linewidth=2.2, linestyle="--", label="ACS heads")
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.axvspan(30, 55, color=COLORS["band"], alpha=0.14, linewidth=0)
    setup_axis(ax, "Periphery - Center Gap", age_min, age_max, "Ownership gap (pp)")
    ax.set_ylim(-20, 90)

    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("Spatial Ownership Lifecycle", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_spatial_lifecycle_panel(
    model: dict[str, dict[str, dict[int, float]]],
    acs_ownership: dict[str, dict[int, float]],
    data_children: dict[str, dict[int, float]],
    outbase: Path,
    age_min: int,
    age_max: int,
) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12.8, 10.2), sharex=True)
    row_specs = [
        ("ownership", "Ownership rate (%)"),
        ("networth", "Mean net worth\nmodel units: b or b+p_iH"),
        ("children", "Share with child under 18 (%)"),
    ]
    locs = ["center", "periphery"]

    for col, loc in enumerate(locs):
        axes[0, col].set_title(LOC_LABEL[loc], fontsize=14, fontweight="bold")
        for row, (metric, ylabel) in enumerate(row_specs):
            ax = axes[row, col]
            ax.axvspan(30, 55, color=COLORS["band"], alpha=0.12, linewidth=0)
            ages, vals = series_xy(model[metric][loc], age_min, age_max)
            scale = 100.0 if metric in {"ownership", "children"} else 1.0
            ax.plot(ages, scale * vals, color=COLORS["model"], linewidth=2.4, label="Model")

            if metric == "ownership":
                ages, vals = series_xy(acs_ownership[loc], age_min, age_max)
                ax.plot(ages, 100 * vals, color=COLORS["data"], linewidth=2.1, linestyle="--", label="ACS heads")
                ax.set_ylim(0, 100)
            elif metric == "children":
                ages, vals = series_xy(data_children[loc], age_min, age_max)
                ax.plot(ages, 100 * vals, color=COLORS["data"], linewidth=2.1, linestyle="--", label="ACS adults")
                ax.set_ylim(0, 75)
            else:
                ax.set_ylabel(ylabel)

            if row != 1:
                ax.set_ylabel(ylabel)
            ax.grid(True, color=COLORS["grid"], alpha=0.55, linewidth=0.7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if row == 2:
                ax.set_xlabel("Age")

    axes[0, 0].legend(frameon=False, loc="upper left")
    fig.suptitle("Spatial Lifecycle Diagnostics", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0.03, 0.02, 1, 0.97])
    fig.savefig(outbase.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_panel_csv(
    path: Path,
    case: str,
    model: dict[str, dict[str, dict[int, float]]],
    acs_ownership: dict[str, dict[int, float]],
    data_children: dict[str, dict[int, float]],
) -> None:
    rows = []
    for metric, by_loc in model.items():
        for loc, values in by_loc.items():
            for age, value in sorted(values.items()):
                rows.append({
                    "source": "model",
                    "case": case,
                    "metric": metric,
                    "sample": "",
                    "location": loc,
                    "age": age,
                    "value": value,
                })
    for loc, values in acs_ownership.items():
        for age, value in sorted(values.items()):
            rows.append({
                "source": "ACS",
                "case": "",
                "metric": "ownership",
                "sample": ACS_SAMPLE,
                "location": loc,
                "age": age,
                "value": value,
            })
    for loc, values in data_children.items():
        for age, value in sorted(values.items()):
            rows.append({
                "source": "ACS",
                "case": "",
                "metric": "children",
                "sample": "all_adults_mms_age_profiles_full",
                "location": loc,
                "age": age,
                "value": value,
            })

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "case", "metric", "sample", "location", "age", "value"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_combined_csv(
    path: Path,
    case: str,
    model: dict[str, dict[int, float]],
    acs_overall: dict[int, float],
    acs_loc: dict[str, dict[int, float]],
) -> None:
    rows = []
    for loc, values in sorted(model.items()):
        for age, rate in sorted(values.items()):
            rows.append({
                "source": "model",
                "case": case,
                "sample": "",
                "location": loc,
                "age": age,
                "owner_rate": rate,
            })
    for age, rate in sorted(acs_overall.items()):
        rows.append({
            "source": "ACS",
            "case": "",
            "sample": ACS_SAMPLE,
            "location": "all",
            "age": age,
            "owner_rate": rate,
        })
    for loc, values in sorted(acs_loc.items()):
        for age, rate in sorted(values.items()):
            rows.append({
                "source": "ACS",
                "case": "",
                "sample": ACS_SAMPLE,
                "location": loc,
                "age": age,
                "owner_rate": rate,
            })

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "case", "sample", "location", "age", "owner_rate"])
        writer.writeheader()
        writer.writerows(rows)


def print_key_numbers(model: dict[str, dict[int, float]], acs_overall: dict[int, float], acs_loc: dict[str, dict[int, float]]) -> None:
    ages = [22, 25, 27, 30, 35, 40, 45]
    print("Overall ownership, model vs ACS heads:")
    for age in ages:
        if age in model["all"] and age in acs_overall:
            print(f"  age {age}: model {model['all'][age]:.3f}, ACS {acs_overall[age]:.3f}")
    print("\nSpatial ownership at selected ages:")
    for age in [22, 27, 30, 35]:
        pieces = []
        for loc in ["center", "periphery"]:
            if age in model[loc] and age in acs_loc.get(loc, {}):
                pieces.append(f"{loc}: model {model[loc][age]:.3f}, ACS {acs_loc[loc][age]:.3f}")
        if pieces:
            print(f"  age {age}: " + "; ".join(pieces))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=DEFAULT_CASE)
    parser.add_argument("--record", type=Path, default=DEFAULT_RECORD)
    parser.add_argument("--model-csv", type=Path, default=DEFAULT_MODEL_CSV)
    parser.add_argument("--acs-age-csv", type=Path, default=DEFAULT_ACS_AGE_CSV)
    parser.add_argument("--acs-location-csv", type=Path, default=DEFAULT_ACS_LOC_CSV)
    parser.add_argument("--mms-age-csv", type=Path, default=DEFAULT_MMS_AGE_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--age-min", type=int, default=20)
    parser.add_argument("--age-max", type=int, default=55)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    model = read_model(args.model_csv)
    acs_overall = read_acs_overall(args.acs_age_csv)
    acs_loc = read_acs_location(args.acs_location_csv)
    data_children = read_mms_children_location(args.mms_age_csv)

    plot_overall(
        model,
        acs_overall,
        args.outdir / "ownership_lifecycle_heads_model_vs_data",
        args.age_min,
        args.age_max,
    )
    plot_spatial(
        model,
        acs_loc,
        args.outdir / "ownership_lifecycle_spatial_model_vs_data",
        args.age_min,
        args.age_max,
    )
    panel_model = compute_model_spatial_lifecycle(args.record)
    plot_spatial_lifecycle_panel(
        panel_model,
        acs_loc,
        data_children,
        args.outdir / "spatial_lifecycle_3x2_model_vs_data",
        args.age_min,
        args.age_max,
    )
    write_panel_csv(
        args.outdir / "spatial_lifecycle_3x2_model_vs_data.csv",
        args.case,
        panel_model,
        acs_loc,
        data_children,
    )
    write_combined_csv(
        args.outdir / "ownership_lifecycle_spatial_model_vs_data.csv",
        args.case,
        model,
        acs_overall,
        acs_loc,
    )
    print_key_numbers(model, acs_overall, acs_loc)
    print(f"\nWrote figures and CSV to {args.outdir}")


if __name__ == "__main__":
    main()
