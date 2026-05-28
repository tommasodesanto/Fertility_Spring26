#!/usr/bin/env python3
"""Generate the slide-model figure bundle for direct-calibration best records."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]

import sys

if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.objective import extract_moments  # noqa: E402
from dt_cp_model.parameters import asdict  # noqa: E402
from dt_cp_model.solver import run_model_cp_dt  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402
from dt_cp_model.utils import make_grid  # noqa: E402
from export_distributional_discipline import age_index, cell_arrays, current_children, room_bin_masks  # noqa: E402


TARGET_ROWS = [
    ("tfr", "TFR"),
    ("childless_rate", "Childless"),
    ("mean_age_first_birth", "Age first birth"),
    ("own_rate", "Ownership"),
    ("own_gradient", "Ownership gradient"),
    ("own_family_gap", "Family ownership gap"),
    ("own_lifecycle_slope", "Lifecycle ownership slope"),
    ("prime_childless_renter_median_rooms", "Childless renter rooms"),
    ("prime_childless_owner_median_rooms", "Childless owner rooms"),
    ("housing_increment_0to1", "Rooms, 0 to 1 child"),
    ("housing_increment_1to2", "Rooms, 1 to 2 children"),
    ("old_age_own_rate", "Old-age ownership"),
    ("old_age_parent_childless_gap", "Old parent-childless gap"),
    ("inv_pop_share_C", "Center share"),
    ("inv_rent_ratio_C_over_P", "Rent ratio C/P"),
]

ROOM_BINS = ["<=4", "5", "6", "7-8", "9-10", "11+"]
STYLE = {
    "model_blue": (0.11, 0.34, 0.67),
    "model_green": (0.28, 0.56, 0.28),
    "model_purple": (0.62, 0.28, 0.55),
    "data_red": (0.70, 0.18, 0.12),
    "loc_periphery": (0.20, 0.40, 0.80),
    "loc_center": (0.85, 0.25, 0.10),
    "light_gray": (0.82, 0.84, 0.87),
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build slide-style fit figures and a comparison note.")
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help="Case spec LABEL=best.json,hR_max[,H1;H2;...]",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--note-tex", type=Path, required=True)
    parser.add_argument("--compile", action="store_true")
    parser.add_argument("--reuse-figures", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases = []
    for spec in args.case:
        label, record_path, hr_max, H_own = parse_case(spec)
        case = build_case(label, record_path, hr_max, H_own, args.outdir / label, reuse_figures=args.reuse_figures)
        cases.append(case)

    write_note(args.note_tex, args.outdir, cases)
    if args.compile:
        compile_latex(args.note_tex)
    print(f"Wrote note {args.note_tex}")


def parse_case(spec: str) -> tuple[str, Path, float, list[float] | None]:
    if "=" not in spec:
        raise ValueError(f"Case must be LABEL=path,hR_max, got {spec!r}")
    label, rhs = spec.split("=", 1)
    parts = rhs.split(",", 2)
    if len(parts) not in (2, 3):
        raise ValueError(f"Case must be LABEL=path,hR_max[,H1;H2;...], got {spec!r}")
    H_own = None
    if len(parts) == 3 and parts[2].strip():
        H_own = [float(x) for x in parts[2].split(";") if x.strip()]
    return label, Path(parts[0]), float(parts[1]), H_own


def build_case(
    label: str,
    record_path: Path,
    hr_max: float,
    H_own: list[float] | None,
    outdir: Path,
    *,
    reuse_figures: bool = False,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    record = json.loads(record_path.read_text())
    owner_h_bar_scale = record.get("owner_h_bar_scale")
    owner_size_cost = record.get("owner_size_cost")
    owner_size_cost_ref = record.get("owner_size_cost_ref")
    owner_size_cost_power = record.get("owner_size_cost_power")
    tenure_choice_kappa = record.get("tenure_choice_kappa")
    setup = build_direct_calibration_setup(
        "benchmark",
        geo_weight=100.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=hr_max,
        owner_h_bar_scale=owner_h_bar_scale,
        owner_size_cost=owner_size_cost,
        owner_size_cost_ref=owner_size_cost_ref,
        owner_size_cost_power=owner_size_cost_power,
        tenure_choice_kappa=tenure_choice_kappa,
        H_own=H_own,
    )
    moments = dict(record.get("moments", {}))
    figure_names = [
        "benchmark_lifecycle_compare.png",
        "market_equilibrium.png",
        "housing_choice_heatmap.png",
        "housing_size_fit.png",
    ]
    have_figures = all((outdir / name).exists() for name in figure_names)
    if not (reuse_figures and have_figures):
        sol, P, p_eq = solve_record(record, setup)
        b_grid = make_grid(P)
        fresh = extract_moments(sol, P, p_eq, 0.0, 0.0, 0.0, True)
        moments["fresh_own_rate"] = float(getattr(fresh, "own_rate", np.nan))
        moments["fresh_tfr"] = float(getattr(fresh, "tfr", np.nan))

        plot_benchmark_lifecycle_compare(sol, P, b_grid, outdir / "benchmark_lifecycle_compare.png")
        plot_market_equilibrium(sol, P, p_eq, setup, outdir / "market_equilibrium.png", label)
        plot_housing_choice_heatmap(sol, P, b_grid, outdir / "housing_choice_heatmap.png", label)
        plot_housing_size_fit(sol, P, b_grid, setup, moments, outdir / "housing_size_fit.png", label)
    write_case_summary(outdir / "target_comparison.csv", setup.targets, moments, label)

    return {
        "label": label,
        "outdir": outdir,
        "record_path": record_path,
        "hR_max": hr_max,
        "H_own": H_own,
        "record": record,
        "targets": setup.targets,
        "moments": moments,
        "prices": [float(x) for x in np.asarray(record.get("p_eq", [])).reshape(-1)],
    }


def solve_record(record: dict, setup) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
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


def plot_benchmark_lifecycle_compare(sol, P, b_grid: np.ndarray, out_path: Path) -> None:
    """Port of the slide alias `benchmark_lifecycle_compare.png`."""

    age_vec = np.arange(P.J) + P.age_start
    data = load_lifecycle_data()
    model_networth_age, _ = compute_model_age_profiles(sol, P, b_grid)
    model_children = compute_model_children_share(sol, P)
    model_wealth_norm = model_networth_age / max(float(np.nanmax(model_networth_age)), 1e-12)
    scf_ages = np.array([27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82], dtype=float)
    scf_wealth = np.array([31.5, 88.6, 138.6, 134.4, 213.6, 266.1, 321.1, 392.9, 393.5, 438.7, 338.2, 327.2])
    scf_wealth_norm = scf_wealth / max(float(np.max(scf_wealth)), 1e-12)
    own_by_age = np.asarray(getattr(sol, "own_by_age", compute_own_by_age(sol)), dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    fig.suptitle("Lifecycle: Model vs Data", fontweight="bold", fontsize=16)

    ax = axes[0]
    ax.plot(age_vec, 100 * own_by_age, "-", color=STYLE["model_blue"], linewidth=2.8, label="Model")
    ax.plot(data["own_age"], 100 * data["own"], "--", color=STYLE["data_red"], linewidth=2.4, label=data["own_legend"])
    ax.set_xlabel("Age")
    ax.set_ylabel("Ownership rate (%)")
    ax.set_title("Ownership")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    ax = axes[1]
    ax.plot(age_vec, model_wealth_norm, "-", color=STYLE["model_green"], linewidth=2.8, label="Model net worth")
    ax.plot(scf_ages, scf_wealth_norm, "--", color=STYLE["data_red"], linewidth=2.4, label="Data (SCF median)")
    ax.set_xlabel("Age")
    ax.set_ylabel("Normalized wealth")
    ax.set_title("Wealth Shape")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    family_age_cap = max(data["children_age"])
    m_mask = age_vec <= family_age_cap
    d_mask = data["children_age"] <= family_age_cap
    ax = axes[2]
    ax.plot(age_vec[m_mask], model_children[m_mask], "-", color=STYLE["model_purple"], linewidth=2.8, label="Model children-at-home share")
    ax.plot(
        data["children_age"][d_mask],
        data["children"][d_mask],
        "--",
        color=STYLE["data_red"],
        linewidth=2.4,
        label=data["children_legend"],
    )
    ax.set_xlabel("Age")
    ax.set_ylabel(data["children_label"])
    ax.set_title("Families With Children")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_market_equilibrium(sol, P, p_eq, setup, out_path: Path, label: str) -> None:
    """Port of the slide alias `market_equilibrium.png`."""

    locs = ["Periphery", "Center"]
    loc_colors = [STYLE["loc_periphery"], STYLE["loc_center"]]
    prices = np.asarray(p_eq, dtype=float)
    rents = (float(P.q) + float(P.delta) + float(P.tau_H)) * prices

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.5))
    fig.suptitle("Spatial Equilibrium", fontsize=14)

    ax = axes[0]
    x = np.arange(2)
    width = 0.35
    ax.bar(x - width / 2, [prices[0], rents[0]], width, color=loc_colors[0], label=locs[0])
    ax.bar(x + width / 2, [prices[1], rents[1]], width, color=loc_colors[1], label=locs[1])
    ax.set_xticks(x, ["Price p", "Rent r"])
    ax.set_ylabel("Level")
    ax.set_title("Equilibrium Prices")
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    own_by_loc = np.asarray(getattr(sol, "own_by_loc", compute_own_by_loc(sol, P)), dtype=float)
    mean_parity_by_loc = np.asarray(getattr(sol, "mean_parity_by_loc", compute_mean_parity_by_loc(sol, P)), dtype=float)
    mean_parity = max(float(getattr(sol, "mean_parity", np.nanmean(mean_parity_by_loc))), 0.01)
    bar_data = np.vstack([np.asarray(sol.pop_share), own_by_loc, mean_parity_by_loc / mean_parity])
    x = np.arange(3)
    ax.bar(x - width / 2, bar_data[:, 0], width, color=loc_colors[0], label=locs[0])
    ax.bar(x + width / 2, bar_data[:, 1], width, color=loc_colors[1], label=locs[1])
    ax.set_xticks(x, ["Pop share", "Own rate", "Rel. fertility"], rotation=12)
    ax.set_ylabel("Share / Rate")
    ax.set_title("Population & Outcomes")
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    mean_c, mean_h = mean_consumption_housing_by_loc(sol, P)
    ax = axes[2]
    x = np.arange(2)
    ax.bar(x - width / 2, [mean_c[0], mean_h[0]], width, color=loc_colors[0], label=locs[0])
    ax.bar(x + width / 2, [mean_c[1], mean_h[1]], width, color=loc_colors[1], label=locs[1])
    ax.set_xticks(x, ["Mean c", "Mean h"])
    ax.set_ylabel("Level")
    ax.set_title("Consumption & Housing")
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    fert_dist_loc = fertility_distribution_by_loc(sol, P)
    ax = axes[3]
    x = np.arange(P.n_parity)
    ax.bar(x - width / 2, fert_dist_loc[:, 0], width, color=loc_colors[0], label=locs[0])
    ax.bar(x + width / 2, fert_dist_loc[:, 1], width, color=loc_colors[1], label=locs[1])
    ax.set_xlabel("Number of children")
    ax.set_ylabel("Share")
    ax.set_title("Fertility Distribution")
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_housing_choice_heatmap(sol, P, b_grid: np.ndarray, out_path: Path, label: str) -> None:
    """Port of the slide alias `housing_choice_heatmap.png`."""

    age_axis = np.arange(P.J) + P.age_start
    ib = np.where((b_grid >= 0) & (b_grid <= 16))[0]
    jf = np.arange(P.A_f_start - 1, P.A_f_end)
    nt = P.n_house + 1
    labels = ["Rent"] + [f"H{k}" for k in range(1, P.n_house + 1)]
    cmap = plt.get_cmap("viridis", nt)
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2), sharey=True)
    for loc, ax in enumerate(axes):
        tc = np.asarray(sol.tenure_choice[np.ix_(ib, [0], [loc], jf, [0], [0])], dtype=float)
        tc = np.squeeze(tc).T + 1.0
        im = ax.imshow(
            tc,
            origin="lower",
            aspect="auto",
            extent=[b_grid[ib[0]], b_grid[ib[-1]], age_axis[jf[0]], age_axis[jf[-1]]],
            vmin=0.5,
            vmax=nt + 0.5,
            cmap=cmap,
        )
        ax.set_xlabel("Liquid wealth")
        if loc == 0:
            ax.set_ylabel("Liquid wealth")
            ax.set_ylabel("Age")
        ax.set_title(["Periphery", "Center"][loc])
        ax.set_xlim(b_grid[ib[0]], b_grid[ib[-1]])
        ax.set_ylim(age_axis[jf[0]], age_axis[jf[-1]])
        ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.88)
    cbar.set_ticks(np.arange(1, nt + 1))
    cbar.set_ticklabels(labels)
    fig.suptitle("Housing Choice Map: Childless Renters in the Fertile Window", fontweight="bold", fontsize=18)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_housing_size_fit(sol, P, b_grid: np.ndarray, setup, moments: dict, out_path: Path, label: str) -> None:
    """Port of the slide alias `housing_size_fit.png`."""

    cmp_summary, cmp_bins = build_room_compare_tables(sol, P, b_grid)
    fig, axes = plt.subplots(2, 2, figsize=(17, 9.5))
    fig.suptitle("Housing Size Fit: Model vs Data", fontweight="bold", fontsize=18)

    ax = axes[0, 0]
    rows = [pick_summary(cmp_summary, "all", tenure, "all") for tenure in ["Owner", "Renter"]]
    grouped_bar(
        ax,
        ["Owner", "Renter"],
        [row["p50_acs"] for row in rows],
        [row["p50_model"] for row in rows],
        "Median rooms",
        "All-Household Tenure Medians",
    )

    ax = axes[0, 1]
    labels = ["Owner, 0", "Owner, 1", "Owner, 2+", "Renter, 0", "Renter, 1", "Renter, 2+"]
    rows = [
        pick_summary(cmp_summary, "25_45", "Owner", "0"),
        pick_summary(cmp_summary, "25_45", "Owner", "1"),
        pick_summary(cmp_summary, "25_45", "Owner", "2+"),
        pick_summary(cmp_summary, "25_45", "Renter", "0"),
        pick_summary(cmp_summary, "25_45", "Renter", "1"),
        pick_summary(cmp_summary, "25_45", "Renter", "2+"),
    ]
    grouped_bar(
        ax,
        labels,
        [row["p50_acs"] for row in rows],
        [row["p50_model"] for row in rows],
        "Median rooms",
        "Prime-Age Medians by Tenure and Children",
        rotation=25,
    )

    plot_bin_panel(axes[1, 0], cmp_bins, "25_45", "Renter", "all", "Prime-Age Renters: Room Bin Shares")
    plot_bin_panel(axes[1, 1], cmp_bins, "25_45", "Owner", "all", "Prime-Age Owners: Room Bin Shares")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def load_acs_room_bins() -> dict[tuple[str, str, str], float]:
    path = ROOT / "output/model/acs_2023_rooms_bins.csv"
    out: dict[tuple[str, str, str], float] = {}
    if not path.exists():
        return out
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("age_window") != "25_45":
                continue
            out[(row["tenure"], row["child_bin"], row["bin"])] = float(row["share"])
    return out


def load_lifecycle_data() -> dict[str, np.ndarray | str]:
    mms_path = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles_full.csv"
    if not mms_path.exists():
        mms_path = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles.csv"
    fert_path = ROOT / "code/data/mms_center_periphery/output/mms_age_fertility_profiles.csv"
    own_path = ROOT / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"

    mms_rows = read_csv_rows(mms_path)
    fert_rows = read_csv_rows(fert_path)
    own_rows = [
        row
        for row in read_csv_rows(own_path)
        if row.get("source") == "ACS"
        and row.get("sample") == "household_heads_hhwt_due_housing"
        and row.get("owner_rate", "") != ""
    ]
    if not own_rows:
        raise ValueError(f"No ACS household-head DUE ownership lifecycle rows found in {own_path}")
    child_col = "has_child_u18_rate" if "has_child_u18_rate" in mms_rows[0] else "has_children_rate"
    child_label = "Share with child under 18" if child_col == "has_child_u18_rate" else "Share with children"
    child_legend = "Data child-under-18 share" if child_col == "has_child_u18_rate" else "Data has-children rate"

    own_age = np.asarray([float(row["age"]) for row in own_rows], dtype=float)
    own = np.asarray([float(row["owner_rate"]) for row in own_rows], dtype=float)
    own_order = np.argsort(own_age)
    own_age = own_age[own_order]
    own = own[own_order]
    rooms_age, rooms = collapse_age_profile(mms_rows, "age", "pop_weight", "mean_rooms")
    children_age, children = collapse_age_profile(mms_rows, "age", "pop_weight", child_col)
    fert_age, newparent = collapse_age_profile(fert_rows, "age", "weighted_n", "newparent_rate")
    return {
        "own_age": own_age,
        "own": own,
        "own_legend": "Data (ACS heads, DUE)",
        "rooms_age": rooms_age,
        "rooms": rooms,
        "children_age": children_age,
        "children": children,
        "fert_age": fert_age,
        "newparent": newparent,
        "children_label": child_label,
        "children_legend": child_legend,
    }


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def collapse_age_profile(rows: list[dict[str, str]], age_col: str, weight_col: str, value_col: str) -> tuple[np.ndarray, np.ndarray]:
    ages = sorted({int(float(row[age_col])) for row in rows if row.get(value_col, "") != ""})
    values = []
    for age in ages:
        wsum = 0.0
        vsum = 0.0
        for row in rows:
            if int(float(row[age_col])) != age or row.get(value_col, "") == "":
                continue
            w = float(row[weight_col])
            v = float(row[value_col])
            wsum += w
            vsum += w * v
        values.append(vsum / max(wsum, 1e-12))
    return np.asarray(ages, dtype=float), np.asarray(values, dtype=float)


def compute_model_age_profiles(sol, P, b_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    g = sol.g
    networth_age = np.zeros(P.J)
    rooms_age = np.zeros(P.J)
    for j in range(P.J):
        wealth_sum = 0.0
        room_sum = 0.0
        mass_sum = 0.0
        for loc in range(P.I):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    gr = g[:, 0, loc, j, nn, cs]
                    mr = float(np.sum(gr))
                    if mr > 1e-15:
                        wealth_sum += float(np.sum(gr * b_grid))
                        room_sum += float(np.sum(gr * sol.hR_pol[:, 0, loc, j, nn, cs]))
                        mass_sum += mr
                    for ten in range(1, 1 + P.n_house):
                        go = g[:, ten, loc, j, nn, cs]
                        mo = float(np.sum(go))
                        if mo <= 1e-15:
                            continue
                        hk = float(P.H_own[ten - 1])
                        wealth_sum += float(np.sum(go * (b_grid + sol.p_eq[loc] * hk)))
                        room_sum += mo * hk
                        mass_sum += mo
        networth_age[j] = wealth_sum / max(mass_sum, 1e-12)
        rooms_age[j] = room_sum / max(mass_sum, 1e-12)
    return networth_age, rooms_age


def compute_model_children_share(sol, P) -> np.ndarray:
    out = np.zeros(P.J)
    for j in range(P.J):
        mass = float(np.sum(sol.g[:, :, :, j, :, :]))
        child_mass = 0.0
        for nn in range(P.n_parity):
            for cs in range(P.n_child_states):
                if current_children(nn, cs, P) > 0:
                    child_mass += float(np.sum(sol.g[:, :, :, j, nn, cs]))
        out[j] = child_mass / max(mass, 1e-12)
    return out


def compute_own_by_age(sol) -> np.ndarray:
    g = sol.g
    out = np.zeros(g.shape[3])
    for j in range(g.shape[3]):
        mass = float(np.sum(g[:, :, :, j, :, :]))
        out[j] = float(np.sum(g[:, 1:, :, j, :, :])) / max(mass, 1e-12)
    return out


def compute_own_by_loc(sol, P) -> np.ndarray:
    out = np.zeros(P.I)
    for loc in range(P.I):
        mass = float(np.sum(sol.g[:, :, loc, :, :, :]))
        out[loc] = float(np.sum(sol.g[:, 1:, loc, :, :, :])) / max(mass, 1e-12)
    return out


def compute_mean_parity_by_loc(sol, P) -> np.ndarray:
    out = np.zeros(P.I)
    for loc in range(P.I):
        mass = float(np.sum(sol.g[:, :, loc, P.A_f_end :, :, :]))
        val = 0.0
        for nn in range(P.n_parity):
            val += nn * float(np.sum(sol.g[:, :, loc, P.A_f_end :, nn, :]))
        out[loc] = val / max(mass, 1e-12)
    return out


def mean_consumption_housing_by_loc(sol, P) -> tuple[np.ndarray, np.ndarray]:
    mean_c = np.zeros(P.I)
    mean_h = np.zeros(P.I)
    mass_loc = np.zeros(P.I)
    for loc in range(P.I):
        for j in range(P.J):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    gr = sol.g[:, 0, loc, j, nn, cs]
                    nz = gr > 1e-15
                    if np.any(nz):
                        mean_c[loc] += float(np.sum(gr[nz] * sol.c_pol[nz, 0, loc, j, nn, cs]))
                        mean_h[loc] += float(np.sum(gr[nz] * sol.hR_pol[nz, 0, loc, j, nn, cs]))
                        mass_loc[loc] += float(np.sum(gr[nz]))
                    for ten in range(1, 1 + P.n_house):
                        go = sol.g[:, ten, loc, j, nn, cs]
                        nz_o = go > 1e-15
                        if np.any(nz_o):
                            mean_c[loc] += float(np.sum(go[nz_o] * sol.c_pol[nz_o, ten, loc, j, nn, cs]))
                            mean_h[loc] += float(np.sum(go[nz_o])) * float(P.H_own[ten - 1])
                            mass_loc[loc] += float(np.sum(go[nz_o]))
    return mean_c / np.maximum(mass_loc, 1e-12), mean_h / np.maximum(mass_loc, 1e-12)


def fertility_distribution_by_loc(sol, P) -> np.ndarray:
    out = np.zeros((P.n_parity, P.I))
    for loc in range(P.I):
        mass = float(np.sum(sol.g[:, :, loc, P.A_f_end :, :, :]))
        for nn in range(P.n_parity):
            out[nn, loc] = float(np.sum(sol.g[:, :, loc, P.A_f_end :, nn, :])) / max(mass, 1e-12)
    return out


def model_room_bin_shares(sol, P, b_grid: np.ndarray, *, tenure: str, child_bin: str) -> list[float]:
    ten_filter = (lambda ten: ten == 0) if tenure == "Renter" else (lambda ten: ten > 0)
    if child_bin == "0":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 0
    elif child_bin == "1":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 1
    else:
        child_filter = lambda nn, cs: current_children(nn, cs, P) >= 2
    vals, wts = cell_arrays(
        sol,
        P,
        b_grid,
        ages=range(age_index(P, 25), age_index(P, 45) + 1),
        mask_fn=lambda ten, loc, nn, cs: ten_filter(ten) and child_filter(nn, cs),
    )
    if wts.size == 0 or float(np.sum(wts)) <= 0:
        return [float("nan")] * len(ROOM_BINS)
    masks = room_bin_masks(vals["rooms"])
    total = float(np.sum(wts))
    return [float(np.sum(wts[masks[name]]) / total) for name in ROOM_BINS]


def build_room_compare_tables(sol, P, b_grid: np.ndarray) -> tuple[list[dict], list[dict]]:
    acs_summary = load_acs_room_summary()
    acs_bins = load_acs_room_bins_full()
    summary = []
    bins = []
    for age_window, ages in {
        "all": range(P.J),
        "25_45": range(age_index(P, 25), age_index(P, 45) + 1),
    }.items():
        for tenure in ["Owner", "Renter"]:
            for child_bin in ["0", "1", "2+", "all"]:
                vals, wts = model_room_values(sol, P, b_grid, ages=ages, tenure=tenure, child_bin=child_bin)
                if wts.size == 0 or float(np.sum(wts)) <= 0:
                    continue
                key = (age_window, tenure, child_bin)
                acs = acs_summary.get(key, {})
                row = {
                    "age_window": age_window,
                    "tenure": tenure,
                    "child_bin": child_bin,
                    "p50_acs": acs.get("p50", np.nan),
                    "p50_model": weighted_quantile(vals, wts, 0.50),
                    "mean_acs": acs.get("mean", np.nan),
                    "mean_model": weighted_mean(vals, wts),
                }
                summary.append(row)
                masks = room_bin_masks(vals)
                total = float(np.sum(wts))
                for bin_name in ROOM_BINS:
                    acs_share = acs_bins.get((*key, bin_name), np.nan)
                    model_share = float(np.sum(wts[masks[bin_name]]) / max(total, 1e-12))
                    bins.append(
                        {
                            "age_window": age_window,
                            "tenure": tenure,
                            "child_bin": child_bin,
                            "bin": bin_name,
                            "share_acs": acs_share,
                            "share_model": model_share,
                        }
                    )
    return summary, bins


def load_acs_room_summary() -> dict[tuple[str, str, str], dict[str, float]]:
    path = ROOT / "output/model/acs_2023_rooms_summary.csv"
    out = {}
    if not path.exists():
        return out
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out[(row["age_window"], row["tenure"], row["child_bin"])] = {
                "mean": float(row["mean"]),
                "p50": float(row["p50"]),
            }
    return out


def load_acs_room_bins_full() -> dict[tuple[str, str, str, str], float]:
    path = ROOT / "output/model/acs_2023_rooms_bins.csv"
    out = {}
    if not path.exists():
        return out
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out[(row["age_window"], row["tenure"], row["child_bin"], row["bin"])] = float(row["share"])
    return out


def model_room_values(sol, P, b_grid: np.ndarray, *, ages, tenure: str, child_bin: str) -> tuple[np.ndarray, np.ndarray]:
    ten_filter = (lambda ten: ten > 0) if tenure == "Owner" else (lambda ten: ten == 0)
    if child_bin == "0":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 0
    elif child_bin == "1":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 1
    elif child_bin == "2+":
        child_filter = lambda nn, cs: current_children(nn, cs, P) >= 2
    else:
        child_filter = lambda nn, cs: True
    vals, wts = cell_arrays(
        sol,
        P,
        b_grid,
        ages=ages,
        mask_fn=lambda ten, loc, nn, cs: ten_filter(ten) and child_filter(nn, cs),
    )
    return vals["rooms"], wts


def pick_summary(rows: list[dict], age_window: str, tenure: str, child_bin: str) -> dict:
    matches = [row for row in rows if row["age_window"] == age_window and row["tenure"] == tenure and row["child_bin"] == child_bin]
    if len(matches) != 1:
        raise ValueError(f"Could not find unique summary row for {age_window} / {tenure} / {child_bin}")
    return matches[0]


def grouped_bar(ax, labels: list[str], acs_vals: list[float], model_vals: list[float], ylabel: str, title: str, rotation: float = 0) -> None:
    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, acs_vals, width, color=STYLE["light_gray"], label="ACS", linewidth=0.8)
    ax.bar(x + width / 2, model_vals, width, color=STYLE["model_blue"], label="Model", linewidth=0.8)
    ax.set_xticks(x, labels, rotation=rotation, ha="right" if rotation else "center")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)


def plot_bin_panel(ax, rows: list[dict], age_window: str, tenure: str, child_bin: str, title: str) -> None:
    panel = [row for row in rows if row["age_window"] == age_window and row["tenure"] == tenure and row["child_bin"] == child_bin]
    order = [b for b in ROOM_BINS if any(row["bin"] == b for row in panel)]
    panel_by_bin = {row["bin"]: row for row in panel}
    acs_vals = [panel_by_bin[b]["share_acs"] for b in order]
    model_vals = [panel_by_bin[b]["share_model"] for b in order]
    grouped_bar(ax, order, acs_vals, model_vals, "Share", title)
    ymax = max(0.75, 1.08 * max(max(acs_vals), max(model_vals)))
    ax.set_ylim(0, ymax)


def weighted_mean(x, w) -> float:
    return float(np.sum(np.asarray(x, dtype=float) * np.asarray(w, dtype=float)) / max(float(np.sum(w)), 1e-12))


def weighted_quantile(x, w, q: float) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    x = x[ok]
    w = w[ok]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cw = np.cumsum(w)
    return float(x[np.searchsorted(cw, q * cw[-1], side="left")])


def write_case_summary(path: Path, targets: dict, moments: dict, label: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "moment", "target", "model", "diff"])
        for key, _ in TARGET_ROWS:
            if key not in targets:
                continue
            target = float(targets[key])
            model = float(moments.get(key, np.nan))
            writer.writerow([label, key, target, model, model - target])


def write_note(note_tex: Path, outdir: Path, cases: list[dict]) -> None:
    note_tex.parent.mkdir(parents=True, exist_ok=True)
    rel_base = Path("../") / outdir.resolve().relative_to(ROOT)
    rows = []
    for key, name in TARGET_ROWS:
        vals = []
        target = None
        for case in cases:
            if key in case["targets"]:
                target = float(case["targets"][key])
                vals.append(float(case["moments"].get(key, np.nan)))
            else:
                vals.append(np.nan)
        if target is None:
            continue
        rows.append((name, target, vals))

    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=0.7in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{float}",
        r"\usepackage{hyperref}",
        r"\title{Current Best-Fit Slide Diagnostics}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Source Runs}",
        r"\begin{itemize}",
    ]
    for case in cases:
        rec = case["record"]
        lines.append(
            rf"\item \textbf{{{latex_escape(case['label'])}}}: loss {float(rec['loss']):.3f}, "
            rf"job {rec.get('job_id')}, eval {rec.get('eval_id')}, $h_R^{{\max}}={case['hR_max']:.1f}$."
        )
    lines.extend(
        [
            r"\end{itemize}",
            r"\section*{Moment Table}",
            r"\begin{table}[H]\centering\small",
            r"\begin{tabular}{l" + "r" * (2 + len(cases)) + r"}",
            r"\toprule",
            "Moment & Target & " + " & ".join(latex_escape(c["label"]) for c in cases) + r" \\",
            r"\midrule",
        ]
    )
    for name, target, vals in rows:
        lines.append(
            latex_escape(name)
            + f" & {target:.3f} & "
            + " & ".join("NA" if not np.isfinite(v) else f"{v:.3f}" for v in vals)
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])

    lines.append(r"\section*{Slide Figure Bundle}")
    for case in cases:
        case_rel = rel_base / case["label"]
        lines.extend(
            [
                rf"\subsection*{{{latex_escape(case['label'])}}}",
                figure_line(case_rel / "benchmark_lifecycle_compare.png", "Lifecycle fit"),
                figure_line(case_rel / "market_equilibrium.png", "Market equilibrium"),
                figure_line(case_rel / "housing_choice_heatmap.png", "Housing choice heatmap"),
                figure_line(case_rel / "housing_size_fit.png", "Appendix housing-size fit"),
            ]
        )
    lines.append(r"\end{document}")
    note_tex.write_text("\n".join(lines) + "\n")


def figure_line(path: Path, caption: str) -> str:
    return (
        r"\begin{figure}[H]\centering"
        + "\n"
        + rf"\includegraphics[width=0.92\textwidth]{{{path.as_posix()}}}"
        + "\n"
        + rf"\caption{{{latex_escape(caption)}}}"
        + "\n"
        + r"\end{figure}"
    )


def latex_escape(text: str) -> str:
    return (
        str(text)
        .replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("$", r"\$")
        .replace("#", r"\#")
        .replace("_", r"\_")
        .replace("{", r"\{")
        .replace("}", r"\}")
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
