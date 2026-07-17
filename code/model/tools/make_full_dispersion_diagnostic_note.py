#!/usr/bin/env python3
"""Build a full diagnostic packet for a dispersion-wave calibration record.

Unlike the quick snapshot note, this script re-solves the saved record and
exports the policy, distribution, and room-pattern objects from the realized
stationary distribution.  Outputs are kept together in one result subfolder.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = Path(__file__).resolve().parent
for path in (MODEL_ROOT, TOOLS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.objective import extract_moments  # noqa: E402
from dt_cp_model.parameters import asdict  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402
from dt_cp_model.utils import make_grid  # noqa: E402
from make_direct_fit_slide_note import (  # noqa: E402
    ROOM_BINS,
    build_room_compare_tables,
    compute_mean_parity_by_loc,
    compute_own_by_age,
    compute_own_by_loc,
    fertility_distribution_by_loc,
    latex_escape,
    mean_consumption_housing_by_loc,
    plot_benchmark_lifecycle_compare,
    plot_housing_choice_heatmap,
    plot_housing_size_fit,
    plot_market_equilibrium,
)
from export_room_pattern_diagnostics import (  # noqa: E402
    event_context_rows,
    owner_rung_rows,
    ownership_by_age_rows,
    ownership_profile_rows,
    renter_cap_rows,
    room_distribution_rows,
    write_csv,
)


TARGET_ROWS = [
    ("tfr", "TFR"),
    ("childless_rate", "Childless"),
    ("mean_age_first_birth", "Age first birth"),
    ("tfr_gradient", "TFR gradient"),
    ("own_rate", "Ownership"),
    ("own_gradient", "Ownership gradient"),
    ("own_family_gap", "Family ownership gap"),
    ("housing_increment_0to1", "Rooms 0 to 1 child"),
    ("housing_increment_1to2", "Rooms 1 to 2 children"),
    ("young_liquid_wealth_to_income", "Young liquid wealth/income"),
    ("center_share_nonparents", "Center share nonparents"),
    ("center_share_newparents", "Center share new parents"),
    ("migration_rate", "Migration rate"),
    ("old_age_own_rate", "Old-age ownership"),
    ("old_age_parent_childless_gap", "Old parent-childless gap"),
    ("inv_pop_share_C", "Center share"),
    ("inv_rent_ratio_C_over_P", "Rent ratio C/P"),
    ("housing_exp_share_usercost_total", "Housing expenditure share"),
    ("owner25_45_rooms_le6_share", "Prime owner rooms <=6"),
    ("owner25_45_rooms_7to8_share", "Prime owner rooms 7-8"),
    ("owner25_45_rooms_ge9_share", "Prime owner rooms >=9"),
]

PARAMETER_ORDER = [
    "beta",
    "b_entry_fixed",
    "psi_child",
    "h_bar_jump",
    "h_bar_n",
    "c_bar_n",
    "kappa_fert",
    "chi",
    "kappa_loc",
    "mu_move",
    "theta0",
    "theta_n",
    "h_bar_0",
    "E_C",
    "r_bar_C",
    "alpha_cons",
    "tenure_choice_kappa",
    "owner_h_bar_scale",
    "owner_size_cost",
    "owner_size_cost_ref",
    "owner_size_cost_power",
    "hR_max",
    "p_P",
    "p_C",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a full dispersion diagnostic note.")
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--hR-max", type=float, default=8.0)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--headline-moments",
        choices=("fresh", "record"),
        default="fresh",
        help="Moment source for headline target tables and fit bars.",
    )
    parser.add_argument(
        "--record-mismatch",
        choices=("fail", "warn", "ignore"),
        default="fail",
        help="What to do if saved record moments differ from the local re-solve.",
    )
    parser.add_argument(
        "--moment-check-tol",
        type=float,
        default=1e-4,
        help="Absolute tolerance for saved-vs-re-solve moment checks.",
    )
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args()

    outdir = args.outdir
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    records_dir = outdir / "records"
    for path in (tables_dir, figures_dir, records_dir):
        path.mkdir(parents=True, exist_ok=True)

    record = json.loads(args.record.read_text())
    (records_dir / "best.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    setup = build_setup(record, args.hR_max)
    sol, P, p_eq = solve_record(record, setup)
    b_grid = make_grid(P)
    fresh = extract_moments(sol, P, p_eq, 0.0, 0.0, 0.0, True)
    moments = dict(record.get("moments") or {})
    fresh_moments = {key: float(value) for key, value in vars(fresh).items() if is_number(value)}
    fresh_moments["inv_pop_share_C"] = first_finite(
        getattr(fresh, "inv_pop_share_C", np.nan),
        getattr(fresh, "pop_share_C", np.nan),
        moments.get("inv_pop_share_C", np.nan),
    )
    fresh_moments["inv_rent_ratio_C_over_P"] = first_finite(
        getattr(fresh, "inv_rent_ratio_C_over_P", np.nan),
        getattr(fresh, "rent_ratio_C_over_P", np.nan),
        moments.get("inv_rent_ratio_C_over_P", np.nan),
    )
    moment_check_rows = build_moment_check_rows(moments, fresh_moments)
    if moment_check_rows:
        write_dict_csv(tables_dir / "record_resolve_moment_check.csv", moment_check_rows)
        mismatches = [row for row in moment_check_rows if abs(row["gap"]) > args.moment_check_tol]
        if mismatches and args.record_mismatch != "ignore":
            max_gap = max(abs(row["gap"]) for row in mismatches)
            message = (
                f"Saved record moments differ from local re-solve for {len(mismatches)} moments "
                f"(max abs gap {max_gap:.6g}). See {tables_dir / 'record_resolve_moment_check.csv'}."
            )
            if args.record_mismatch == "fail":
                raise RuntimeError(message)
            print(f"WARNING: {message}", file=sys.stderr)
    if args.headline_moments == "fresh":
        moments.update(fresh_moments)
    else:
        for key, value in fresh_moments.items():
            moments.setdefault(key, value)

    targets = dict(setup.targets)
    targets.update(record.get("targets") or {})
    weights = dict(setup.weights)
    weights.update(record.get("weights") or {})
    target_rows = build_target_rows(targets, weights, moments)
    param_rows = build_parameter_rows(record, args.hR_max)
    write_dict_csv(tables_dir / "target_comparison.csv", target_rows)
    write_dict_csv(tables_dir / "parameters.csv", param_rows)

    room_rows = room_distribution_rows(args.label, sol, P)
    cap_rows = renter_cap_rows(args.label, sol, P)
    rung_rows = owner_rung_rows(args.label, sol, P)
    own_profile_rows = ownership_profile_rows(args.label, sol, P)
    own_age_rows = ownership_by_age_rows(args.label, sol, P)
    event_rows = event_context_rows(args.label, record, fresh, sol, P)
    write_csv(tables_dir / "room_distribution_full.csv", room_rows)
    write_csv(tables_dir / "renter_cap_shares.csv", cap_rows)
    write_csv(tables_dir / "owner_rung_shares.csv", rung_rows)
    write_csv(tables_dir / "ownership_profiles.csv", own_profile_rows)
    write_csv(tables_dir / "ownership_by_age.csv", own_age_rows)
    write_csv(tables_dir / "event_context.csv", event_rows)

    plot_target_table_bars(target_rows, figures_dir / "target_model_fit.png")
    plot_room_distribution_full(sol, P, b_grid, figures_dir / "full_room_distribution_current_children.png")
    plot_model_room_distribution_from_rows(room_rows, figures_dir / "model_room_distribution_current_vs_completed.png")
    plot_owner_rungs(rung_rows, figures_dir / "owner_rung_distribution.png")
    plot_ownership_age(own_age_rows, figures_dir / "ownership_by_age.png")
    plot_location_tenure_summary(sol, P, figures_dir / "location_tenure_distribution.png")

    plot_benchmark_lifecycle_compare(sol, P, b_grid, figures_dir / "benchmark_lifecycle_compare.png")
    plot_market_equilibrium(sol, P, p_eq, setup, figures_dir / "market_equilibrium.png", args.label)
    plot_housing_choice_heatmap(sol, P, b_grid, figures_dir / "housing_choice_heatmap.png", args.label)
    plot_housing_size_fit(sol, P, b_grid, setup, moments, figures_dir / "housing_size_fit.png", args.label)

    note_tex = outdir / "current_fit_dispersion_full_diagnostics_20260527.tex"
    write_note(note_tex, args.label, args.record, record, target_rows, param_rows, outdir, args.headline_moments)
    if args.compile:
        compile_latex(note_tex)
    print(f"Wrote {note_tex}")


def build_setup(record: dict, hR_max: float):
    alpha_cons_bounds = record.get("alpha_cons_bounds")
    if alpha_cons_bounds is not None:
        alpha_cons_bounds = tuple(float(x) for x in alpha_cons_bounds)
    return build_direct_calibration_setup(
        "benchmark",
        geo_weight=250.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=hR_max,
        alpha_cons=record.get("alpha_cons"),
        owner_h_bar_scale=record.get("owner_h_bar_scale"),
        owner_size_cost=record.get("owner_size_cost"),
        owner_size_cost_ref=record.get("owner_size_cost_ref"),
        owner_size_cost_power=record.get("owner_size_cost_power"),
        tenure_choice_kappa=record.get("tenure_choice_kappa"),
        weight_overrides=record.get("weight_overrides") or None,
        extra_targets=record.get("extra_targets") or None,
        calibrate_alpha_cons=bool(record.get("calibrate_alpha_cons", False)),
        alpha_cons_bounds=alpha_cons_bounds,
        H_own=record.get("H_own") or None,
    )


def solve_record(record: dict, setup):
    theta = np.asarray(record["theta"], dtype=float)
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
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
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    from dt_cp_model.solver import run_model_cp_dt

    return run_model_cp_dt(P, verbose=False)


def is_number(value) -> bool:
    return isinstance(value, (int, float, np.floating)) and np.isfinite(float(value))


def first_finite(*values) -> float:
    for value in values:
        if is_number(value):
            return float(value)
    return float("nan")


def build_target_rows(targets: dict, weights: dict, moments: dict) -> list[dict]:
    rows = []
    for key, label in TARGET_ROWS:
        if key not in targets:
            continue
        target = float(targets[key])
        model = float(moments.get(key, np.nan))
        rows.append(
            {
                "moment": key,
                "label": label,
                "target": target,
                "model": model,
                "gap": model - target,
                "weight": float(weights.get(key, np.nan)),
            }
        )
    return rows


def build_moment_check_rows(record_moments: dict, fresh_moments: dict) -> list[dict]:
    rows = []
    for key, _label in TARGET_ROWS:
        record_value = record_moments.get(key)
        fresh_value = fresh_moments.get(key)
        if not is_number(record_value) or not is_number(fresh_value):
            continue
        record_float = float(record_value)
        fresh_float = float(fresh_value)
        rows.append(
            {
                "moment": key,
                "record": record_float,
                "resolve": fresh_float,
                "gap": fresh_float - record_float,
                "abs_gap": abs(fresh_float - record_float),
            }
        )
    return rows


def build_parameter_rows(record: dict, hR_max: float) -> list[dict]:
    params = record.get("parameters") or {}
    rows = [{"parameter": key, "value": float(params[key])} for key in PARAMETER_ORDER if key in params]
    controls = {
        "alpha_cons": record.get("alpha_cons_effective")
        or (record.get("parameters") or {}).get("alpha_cons")
        or record.get("alpha_cons"),
        "tenure_choice_kappa": record.get("tenure_choice_kappa"),
        "owner_h_bar_scale": record.get("owner_h_bar_scale"),
        "owner_size_cost": record.get("owner_size_cost"),
        "owner_size_cost_ref": record.get("owner_size_cost_ref"),
        "owner_size_cost_power": record.get("owner_size_cost_power"),
        "hR_max": hR_max,
    }
    prices = record.get("p_eq") or []
    if len(prices) > 0:
        controls["p_P"] = prices[0]
    if len(prices) > 1:
        controls["p_C"] = prices[1]
    seen = {row["parameter"] for row in rows}
    for key, value in controls.items():
        if key not in seen and value is not None:
            rows.append({"parameter": key, "value": float(value)})
            seen.add(key)
    return rows


def write_dict_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_target_table_bars(rows: list[dict], path: Path) -> None:
    panel = [row for row in rows if row["moment"] not in {"mean_age_first_birth", "housing_increment_1to2"}]
    labels = [row["label"] for row in panel]
    x = np.arange(len(panel))
    width = 0.36
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.bar(x - width / 2, [row["target"] for row in panel], width, label="Target", color="#c7cbd1")
    ax.bar(x + width / 2, [row["model"] for row in panel], width, label="Model", color="#315f96")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_title("Target and Model Moments")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_room_distribution_full(sol, P, b_grid: np.ndarray, path: Path) -> None:
    _, bin_rows = build_room_compare_tables(sol, P, b_grid)
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), sharey=True)
    for r, tenure in enumerate(["Renter", "Owner"]):
        for c, child_bin in enumerate(["0", "1", "2+", "all"]):
            ax = axes[r, c]
            panel = [
                row
                for row in bin_rows
                if row["age_window"] == "25_45" and row["tenure"] == tenure and row["child_bin"] == child_bin
            ]
            by_bin = {row["bin"]: row for row in panel}
            support_bins = ROOM_BINS if tenure == "Owner" else ["<=4", "5", "6", "7-8"]
            labels = [b for b in support_bins if b in by_bin]
            x = np.arange(len(labels))
            width = 0.36
            ax.bar(x - width / 2, [by_bin[b]["share_acs"] for b in labels], width, color="#c7cbd1", label="ACS")
            ax.bar(x + width / 2, [by_bin[b]["share_model"] for b in labels], width, color="#7b3f8c", label="Model")
            ax.set_title(f"{tenure}, children {child_bin}")
            ax.set_xticks(x, labels, rotation=30)
            ax.grid(True, axis="y", alpha=0.25)
            if c == 0:
                ax.set_ylabel("Share")
            if r == 0 and c == 0:
                ax.legend(loc="upper right", fontsize=8)
    fig.suptitle("Full Prime-Age Room Distribution by Tenure and Current Children", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_model_room_distribution_from_rows(rows: list[dict], path: Path) -> None:
    subset = [
        row
        for row in rows
        if row["age_window"] == "25_45"
        and row["tenure"] in {"renter", "owner"}
        and row["child_bin"] in {"0", "1", "2plus", "all"}
        and row["child_concept"] in {"current_children", "completed_fertility"}
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.2), sharey=True)
    for r, tenure in enumerate(["renter", "owner"]):
        for c, concept in enumerate(["current_children", "completed_fertility"]):
            ax = axes[r, c]
            for child_bin, color in zip(["0", "1", "2plus", "all"], ["#315f96", "#5b8d45", "#9d5634", "#444444"]):
                matches = [
                    row
                    for row in subset
                    if row["tenure"] == tenure and row["child_concept"] == concept and row["child_bin"] == child_bin
                ]
                if not matches:
                    continue
                row = matches[0]
                if tenure == "renter":
                    labels = ["<=4", "5", "6", "7-8"]
                    vals = [
                        float(row["share_rooms_le4"]),
                        float(row["share_rooms_5"]),
                        float(row["share_rooms_6"]),
                        float(row["share_rooms_7_8"]),
                    ]
                else:
                    labels = ROOM_BINS
                    vals = [
                        float(row["share_rooms_le4"]),
                        float(row["share_rooms_5"]),
                        float(row["share_rooms_6"]),
                        float(row["share_rooms_7_8"]),
                        float(row["share_rooms_9_10"]),
                        float(row["share_rooms_11plus"]),
                    ]
                ax.plot(labels, vals, marker="o", label=child_bin, color=color)
            ax.set_title(f"{tenure.title()}, {concept.replace('_', ' ')}")
            ax.grid(True, axis="y", alpha=0.25)
            if c == 0:
                ax.set_ylabel("Model share")
            if r == 0 and c == 1:
                ax.legend(title="children", fontsize=8)
    fig.suptitle("Model Room Distribution: Current vs Completed Fertility State", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_owner_rungs(rows: list[dict], path: Path) -> None:
    panel = [
        row
        for row in rows
        if row["age_window"] == "25_45" and row["child_concept"] == "current_children" and row["child_bin"] == "all"
    ]
    labels = [f"{row['owner_rung']}: {float(row['rooms']):g}" for row in panel]
    shares = [float(row["rung_share"]) for row in panel]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.bar(labels, shares, color="#315f96")
    ax.set_title("Prime-Age Owner Rung Distribution")
    ax.set_xlabel("Owner rung")
    ax.set_ylabel("Share among owners")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_ownership_age(rows: list[dict], path: Path) -> None:
    panel = [row for row in rows if row["location"] == "all" and row["current_child_bin"] == "all"]
    ages = [float(row["age"]) for row in panel]
    rates = [float(row["own_rate"]) for row in panel]
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.plot(ages, rates, color="#315f96", linewidth=2.4)
    ax.set_title("Ownership by Age")
    ax.set_xlabel("Age")
    ax.set_ylabel("Ownership rate")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_location_tenure_summary(sol, P, path: Path) -> None:
    loc_labels = ["P", "C"][: P.I]
    own_loc = compute_own_by_loc(sol, P)
    parity_loc = compute_mean_parity_by_loc(sol, P)
    mean_c, mean_h = mean_consumption_housing_by_loc(sol, P)
    fert_dist = fertility_distribution_by_loc(sol, P)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes[0, 0].bar(loc_labels, own_loc, color="#315f96")
    axes[0, 0].set_title("Ownership by Location")
    axes[0, 0].set_ylim(0, 1)
    axes[0, 1].bar(loc_labels, parity_loc, color="#5b8d45")
    axes[0, 1].set_title("Mean Completed Fertility by Location")
    x = np.arange(P.I)
    width = 0.36
    axes[1, 0].bar(x - width / 2, mean_c, width, label="c", color="#315f96")
    axes[1, 0].bar(x + width / 2, mean_h, width, label="h", color="#9d5634")
    axes[1, 0].set_xticks(x, loc_labels)
    axes[1, 0].set_title("Mean Consumption and Housing")
    axes[1, 0].legend()
    bottom = np.zeros(P.I)
    for n in range(min(P.n_parity, fert_dist.shape[0])):
        axes[1, 1].bar(loc_labels, fert_dist[n, :], bottom=bottom, label=str(n))
        bottom += fert_dist[n, :]
    axes[1, 1].set_title("Completed Fertility Shares")
    axes[1, 1].legend(title="children", fontsize=8)
    for ax in axes.flat:
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def write_note(
    note_tex: Path,
    label: str,
    record_path: Path,
    record: dict,
    target_rows: list[dict],
    param_rows: list[dict],
    outdir: Path,
    headline_moments: str,
) -> None:
    fig_rel = Path("figures")
    table_rel = Path("tables")
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=0.65in]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{graphicx}",
        r"\usepackage{float}",
        r"\usepackage{longtable}",
        r"\title{Full Dispersion-Wave Fit Diagnostics}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Record}",
        r"\begin{itemize}",
        rf"\item Label: \textbf{{{latex_escape(label)}}}.",
        rf"\item Source record copied from: \texttt{{{latex_escape(str(record_path))}}}.",
        rf"\item Run tag: \texttt{{{latex_escape(record.get('run_tag', ''))}}}; job {record.get('job_id')}; eval {record.get('eval_id')}; loss {float(record.get('loss')):.3f}.",
        rf"\item Tenure/product logit $\kappa^T={float(record.get('tenure_choice_kappa') or math.nan):.3f}$; owner size cost {float(record.get('owner_size_cost') or 0.0):.4f}; owner $h$ scale {float(record.get('owner_h_bar_scale') or 1.0):.3f}.",
        rf"\item Owner ladder: [{latex_escape(', '.join(f'{float(x):.2g}' for x in (record.get('H_own') or [])))}].",
        rf"\item Headline target table uses \texttt{{{latex_escape(headline_moments)}}} moments. Policy and distribution figures are generated from a local re-solve at the saved parameters.",
        rf"\item Tables are in \texttt{{{latex_escape(str(table_rel))}}}; figures are in \texttt{{{latex_escape(str(fig_rel))}}}.",
        r"\end{itemize}",
        r"\section*{Target Fit}",
        target_table(target_rows),
        r"\section*{Parameters}",
        parameter_table(param_rows),
        r"\section*{Core Fit Plots}",
        figure_line(fig_rel / "target_model_fit.png", "Target and model moments, side by side."),
        figure_line(fig_rel / "benchmark_lifecycle_compare.png", "Lifecycle policy and distribution checks."),
        figure_line(fig_rel / "market_equilibrium.png", "Market equilibrium."),
        r"\section*{Housing Policy}",
        figure_line(fig_rel / "housing_choice_heatmap.png", "Housing tenure and size policy heatmap."),
        figure_line(fig_rel / "housing_size_fit.png", "Housing-size fit against ACS summaries."),
        r"\section*{Room Distribution}",
        figure_line(fig_rel / "full_room_distribution_current_children.png", "Full room distribution by tenure and current children, ACS vs model."),
        figure_line(fig_rel / "model_room_distribution_current_vs_completed.png", "Model full room distribution by current and completed-fertility states."),
        figure_line(fig_rel / "owner_rung_distribution.png", "Owner ladder/rung distribution."),
        r"\section*{Other Distributions}",
        figure_line(fig_rel / "ownership_by_age.png", "Ownership distribution by age."),
        figure_line(fig_rel / "location_tenure_distribution.png", "Location, tenure, housing, and completed-fertility distribution."),
        r"\end{document}",
    ]
    note_tex.write_text("\n".join(lines) + "\n")


def target_table(rows: list[dict]) -> str:
    out = [
        r"\begin{longtable}{lrrrr}",
        r"\toprule",
        r"Moment & Target & Model & Gap & Weight \\",
        r"\midrule",
    ]
    for row in rows:
        out.append(
            f"{latex_escape(row['label'])} & {row['target']:.3f} & {row['model']:.3f} & {row['gap']:.3f} & {row['weight']:.1f} \\\\"
        )
    out.extend([r"\bottomrule", r"\end{longtable}"])
    return "\n".join(out)


def parameter_table(rows: list[dict]) -> str:
    out = [
        r"\begin{longtable}{lr}",
        r"\toprule",
        r"Parameter & Value \\",
        r"\midrule",
    ]
    for row in rows:
        out.append(f"{latex_escape(row['parameter'])} & {float(row['value']):.5g} \\\\")
    out.extend([r"\bottomrule", r"\end{longtable}"])
    return "\n".join(out)


def figure_line(path: Path, caption: str) -> str:
    return "\n".join(
        [
            r"\begin{figure}[H]\centering",
            rf"\includegraphics[width=0.92\textwidth]{{{path.as_posix()}}}",
            rf"\caption{{{latex_escape(caption)}}}",
            r"\end{figure}",
        ]
    )


def compile_latex(note_tex: Path) -> None:
    for _ in range(2):
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", note_tex.name],
            cwd=note_tex.parent,
            check=True,
        )


if __name__ == "__main__":
    main()
