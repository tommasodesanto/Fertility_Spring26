#!/usr/bin/env python3
"""Build the two established E4 eqscale note draft figures.

The baseline solve deliberately mirrors
``intergen_eqscale_seq_optimized/build_e2_packet.py --arm e4``.  The
resulting figure data are written beside a compact, deterministic verification
record so the rendered note figures remain auditable without a solution cache.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.build_e2_packet import (  # noqa: E402
    E4_SOURCE,
    arm_externals,
    load_winner_theta,
)
from intergen_eqscale_seq_optimized.calibration import extract_moments  # noqa: E402
from intergen_eqscale_seq_optimized.e1_profile import e1_overrides  # noqa: E402
from intergen_eqscale_seq_optimized.solver import run_model_cp_dt  # noqa: E402


OUTDIR = ROOT / "output/model/eqscale_note_draft_figures_20260724"
LIFECYCLE_FIGURE = ROOT / "latex/figures/eqscale_note_lifecycle_equilibrium.png"
DECISION_FIGURE = ROOT / "latex/figures/eqscale_note_decision_rules.png"
TARGET_FIT = ROOT / "output/model/eqscale_seq_e4_split_recalibration_20260723/report/target_fit_full.csv"
ACS = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles_full.csv"
ACS_OWNERSHIP = ROOT / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def certified_targets(path: Path) -> dict[str, float]:
    wanted = {"own_rate", "aggregate_mean_occupied_rooms_18_85", "tfr", "childless_rate"}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {row["moment"]: float(row["model"]) for row in csv.DictReader(handle) if row["moment"] in wanted}
    if set(rows) != wanted:
        raise ValueError(f"Missing certified E4 values in {path}: {sorted(wanted - set(rows))}")
    return rows


def model_age_profiles(sol: Any, P: Any) -> list[dict[str, float]]:
    """Age profiles from the stationary distribution and active housing policy."""
    g = np.asarray(sol.g, dtype=float)
    h_r = np.asarray(sol.hR_pol, dtype=float)
    if g.ndim != 7 or h_r.shape != g.shape:
        raise ValueError("Expected g and hR_pol with shape (b, tenure, market, age, z, parity, child-stage).")
    ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    owner_rooms = np.asarray(P.H_own, dtype=float)
    rows: list[dict[str, float]] = []
    for j, age in enumerate(ages[:-1]):
        gj = g[:, :, :, j, :, :, :]
        mass = float(np.sum(gj))
        owner_mass = float(np.sum(g[:, 1:, :, j, :, :, :]))
        child_mass = float(np.sum(g[:, :, :, j, :, 1:, 1 : int(P.n_child_stages) + 1]))
        housing = float(np.sum(g[:, 0, :, j, :, :, :] * h_r[:, 0, :, j, :, :, :]))
        for tenure, rooms in enumerate(owner_rooms, start=1):
            housing += float(np.sum(g[:, tenure, :, j, :, :, :])) * float(rooms)
        rows.append(
            {
                "age": float(age),
                "model_owner_rate": owner_mass / mass if mass > 1e-14 else math.nan,
                "model_current_child_rate": child_mass / mass if mass > 1e-14 else math.nan,
                "model_mean_rooms": housing / mass if mass > 1e-14 else math.nan,
                "model_population_mass": mass,
            }
        )
    return rows


def acs_age_profiles() -> list[dict[str, float]]:
    accum: dict[int, dict[str, float]] = {}
    with ACS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                age, weight = int(float(row["age"])), float(row["pop_weight"])
                child, rooms = float(row["has_child_u18_rate"]), float(row["mean_rooms"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (np.isfinite([weight, child, rooms]).all() and weight > 0.0):
                continue
            item = accum.setdefault(age, {"weight": 0.0, "child": 0.0, "rooms": 0.0})
            item["weight"] += weight
            item["child"] += weight * child
            item["rooms"] += weight * rooms
    ownership: dict[int, float] = {}
    with ACS_OWNERSHIP.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("source") == "ACS" and row.get("sample") == "household_heads_hhwt_due_housing" and row.get("owner_rate"):
                ownership[int(float(row["age"]))] = float(row["owner_rate"])
    result = []
    for age, item in sorted(accum.items()):
        result.append({"age": float(age), "acs_owner_rate": ownership.get(age, math.nan),
                       "acs_current_child_u18_rate": item["child"] / item["weight"],
                       "acs_mean_rooms": item["rooms"] / item["weight"],
                       "acs_population_weight": item["weight"]})
    if not result:
        raise ValueError("No usable ACS age-profile rows")
    return result


def decision_rule_rows(sol: Any, P: Any, requested_age: float = 42.0) -> list[dict[str, float | int]]:
    """Exact two series in diagnostics.plot_policy_childless_renter()."""
    ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    j = int(np.argmin(np.abs(ages - requested_age)))
    b_grid = np.asarray(sol.b_grid, dtype=float)
    z_grid = np.asarray(getattr(sol, "type_values", P.z_grid), dtype=float)
    h_r, tenure = np.asarray(sol.hR_pol), np.asarray(sol.tenure_choice)
    fert_probs, value = np.asarray(sol.fert_probs), np.asarray(sol.V)
    nvec = np.arange(int(P.n_parity), dtype=float)
    owner_rooms = np.asarray(P.H_own, dtype=float)
    rows: list[dict[str, float | int]] = []
    for zz, z in enumerate(z_grid):
        choice = tenure[:, 0, 0, j, zz, 0, 0]
        rooms = np.where(choice <= 0, h_r[:, 0, 0, j, zz, 0, 0], owner_rooms[np.maximum(choice - 1, 0)])
        family_size = fert_probs[:, 0, 0, j, zz, :] @ nvec
        valid = value[:, 0, 0, j, zz, 0, 0] > -1e9
        for bb, wealth in enumerate(b_grid):
            rows.append({"age": float(ages[j]), "income_state": float(z), "liquid_wealth": float(wealth),
                         "valid_state": int(valid[bb]),
                         "expected_physical_rooms_after_current_choices": float(rooms[bb]) if valid[bb] else math.nan,
                         "expected_family_size": float(family_size[bb]) if valid[bb] else math.nan})
    return rows


def plot_lifecycle(model: list[dict[str, float]], acs: list[dict[str, float]], path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "font.size": 13, "axes.titlesize": 16, "axes.labelsize": 13, "legend.fontsize": 12})
    navy, orange, green, black = "#17365D", "#D9791F", "#6E9D75", "#222222"
    ma = np.asarray([row["age"] for row in model])
    acs = [row for row in acs if float(row["age"]) <= float(np.max(ma))]
    da = np.asarray([row["age"] for row in acs])
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.25))
    for ax in axes:
        ax.axvspan(18.0, 45.0, color="#EFF1F4", zorder=0)
    axes[0].plot(da, [row["acs_owner_rate"] for row in acs], color=black, ls=(0, (2, 2)), lw=2.4, label="ACS data")
    axes[0].plot(ma, [row["model_owner_rate"] for row in model], color=navy, marker="o", ms=6.5, lw=2.8, label="Model")
    axes[0].set(title="Homeownership over the life cycle", xlabel="Age", ylabel="Homeownership rate", ylim=(0.0, 1.03))
    axes[0].yaxis.set_major_formatter(PercentFormatter(1.0)); axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    child_model = axes[1].plot(ma, [row["model_current_child_rate"] for row in model], color=green, marker="o", ms=5.5, lw=2.6, label="Children at home, model")[0]
    child_acs = axes[1].plot(da, [row["acs_current_child_u18_rate"] for row in acs], color=green, ls=(0, (2, 2)), lw=2.4, label="Own child under 18, ACS")[0]
    axes[1].set(title="Children at home and housing space", xlabel="Age", ylabel="Households with children at home", ylim=(0.0, 0.75))
    axes[1].yaxis.set_major_formatter(PercentFormatter(1.0)); axes[1].grid(axis="y", alpha=0.25)
    ax2 = axes[1].twinx()
    rooms_acs = ax2.plot(da, [row["acs_mean_rooms"] for row in acs], color=black, ls=(0, (2, 2)), lw=2.4, label="Rooms, ACS")[0]
    rooms_model = ax2.plot(ma, [row["model_mean_rooms"] for row in model], color=orange, marker="o", ms=5.5, lw=2.8, label="Rooms, model")[0]
    ax2.set_ylabel("Mean rooms", color=orange); ax2.set_ylim(2.5, 6.5)
    axes[1].axvline(45.0, color="#777777", ls="--", lw=1.5)
    axes[1].text(31.5, 0.70, "Fertile window", ha="center", va="top", color="#555555", fontsize=11)
    axes[1].legend([child_model, child_acs, rooms_acs, rooms_model], [line.get_label() for line in (child_model, child_acs, rooms_acs, rooms_model)], frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    axes[0].text(-0.17, 1.04, "(a)", transform=axes[0].transAxes, fontsize=16)
    axes[1].text(-0.14, 1.04, "(b)", transform=axes[1].transAxes, fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.95, top=0.90, bottom=0.22, wspace=0.20)
    path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=180, bbox_inches="tight"); plt.close(fig)


def plot_decision_rules(rows: list[dict[str, float | int]], path: Path, wealth_min: float = -2.2, wealth_max: float = 7.2) -> None:
    import matplotlib.pyplot as plt

    plt.rcdefaults()
    z_values = sorted({float(row["income_state"]) for row in rows})
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 5.9), sharex=True)
    for color, z in zip(["#355F82", "#D77A20", "#4C8C59", "#8061A7", "#A9413B"], z_values):
        group = sorted((row for row in rows if float(row["income_state"]) == z and int(row["valid_state"]) == 1 and wealth_min <= float(row["liquid_wealth"]) <= wealth_max), key=lambda row: float(row["liquid_wealth"]))
        x = [float(row["liquid_wealth"]) for row in group]
        axes[0].plot(x, [float(row["expected_physical_rooms_after_current_choices"]) for row in group], lw=2.0, color=color, label=rf"Income state $z={z:.6g}$")
        axes[1].plot(x, [float(row["expected_family_size"]) for row in group], lw=2.0, color=color)
    axes[0].set_title("(a) Expected physical housing after current choices"); axes[0].set_ylabel("Rooms"); axes[0].legend(frameon=False, loc="upper left", fontsize=10)
    axes[1].set_title("(b) Expected family size"); axes[1].set_ylabel("Expected children (3 denotes 3+)"); axes[1].set_ylim(-0.08, 3.05)
    for ax in axes:
        ax.set(xlabel="Liquid wealth (model units)", xlim=(wealth_min, wealth_max)); ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=180, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    # Exact baseline construction in build_e2_packet.solve_case(..., arm="e4"),
    # without its policy-loop and its unrelated packet writes.
    theta = load_winner_theta(E4_SOURCE, "e4")
    overrides = {**e1_overrides(tight=True, optimized=True), **arm_externals("e4"), **theta}
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    certified = certified_targets(TARGET_FIT)
    verification = {name: {"model": float(moments[name]), "certified": certified[name], "abs_diff": abs(float(moments[name]) - certified[name])} for name in certified}
    for name in ("own_rate", "aggregate_mean_occupied_rooms_18_85", "tfr", "childless_rate"):
        line = verification[name]
        print(f"{name}: model={line['model']:.15f} certified={line['certified']:.15f} abs_diff={line['abs_diff']:.3e}", flush=True)
    if any(item["abs_diff"] > 1e-6 for item in verification.values()):
        raise RuntimeError("Certified E4 moment verification failed; no figures written.")
    model, acs, decision = model_age_profiles(sol, P), acs_age_profiles(), decision_rule_rows(sol, P)
    write_rows(OUTDIR / "lifecycle_equilibrium.csv", [{"source": "model", **row} for row in model] + [{"source": "ACS_2012_2023", **row} for row in acs])
    write_rows(OUTDIR / "decision_rules_age42_childless_renter.csv", decision)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "summary.json").write_text(json.dumps({"source": str(E4_SOURCE), "theta": theta, "verification": verification, "figure_series": {"lifecycle": "stationary distribution and active housing policies", "decision_rules": "E packet childless-renter age-42 policy series"}}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plot_lifecycle(model, acs, LIFECYCLE_FIGURE)
    plot_decision_rules(decision, DECISION_FIGURE)


if __name__ == "__main__":
    main()
