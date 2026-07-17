"""Regenerate the standard intergenerational draft figures from a trusted cache.

The lifecycle figure compares model and ACS household-head profiles.  Its
family panel deliberately reports a *current stock*: the share of households
with a dependent child at home at each age.  In the model this is

    Pr[m(n,s) > 0 | age],

where ``m(n,s) > 0`` means positive completed parity and a dependent child
stage.  The ACS analogue is the weighted share of household heads with an own
child under age 18 (``has_child_u18_rate``).

``solution_cache.pkl`` files are Python pickles and must only be supplied when
they come from a trusted model run.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import math
import pickle
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ACS = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles_full.csv"
DEFAULT_ACS_OWNERSHIP = (
    ROOT / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"
)
DEFAULT_LIFECYCLE = Path("quant_lifecycle_equilibrium_repaired_nb120.png")
DEFAULT_DECISION = Path("quant_decision_rules_repaired_nb120.png")

# Paper-facing family-size labels used by the established decision-rule figure.
# The top category is displayed as 3+.
FAMILY_SIZE_VALUES = np.array([0.0, 2.0, 3.0])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solution-cache", type=Path, required=True)
    parser.add_argument("--acs-csv", type=Path, default=DEFAULT_ACS)
    parser.add_argument("--acs-ownership-csv", type=Path, default=DEFAULT_ACS_OWNERSHIP)
    parser.add_argument("--lifecycle-out", type=Path, default=DEFAULT_LIFECYCLE)
    parser.add_argument("--decision-rules-out", type=Path, default=DEFAULT_DECISION)
    parser.add_argument("--lifecycle-csv", type=Path, default=None)
    parser.add_argument("--decision-rules-csv", type=Path, default=None)
    parser.add_argument("--decision-age", type=float, default=42.0)
    parser.add_argument("--decision-wealth-min", type=float, default=-2.2)
    parser.add_argument("--decision-wealth-max", type=float, default=7.2)
    return parser.parse_args()


def install_numpy_pickle_compatibility_aliases() -> None:
    """Allow NumPy-2 pickles to load in the repository's NumPy-1.24 venv."""

    try:
        importlib.import_module("numpy._core")
        return
    except ImportError:
        pass
    sys.modules["numpy._core"] = np.core
    for module in ("multiarray", "numeric", "umath", "_multiarray_umath"):
        try:
            sys.modules[f"numpy._core.{module}"] = importlib.import_module(f"numpy.core.{module}")
        except ImportError:
            pass


def load_trusted_solution_cache(path: Path) -> tuple[Any, Any]:
    install_numpy_pickle_compatibility_aliases()
    with path.open("rb") as handle:
        payload = pickle.load(handle)  # noqa: S301 -- CLI explicitly requires a trusted cache.
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"Not an intergenerational solution cache: {path}")
    status = str(payload.get("status", ""))
    if "trusted" not in status:
        raise ValueError(f"Cache does not carry the trusted-cache status marker: {path}")
    baseline = payload["baseline"]
    if not isinstance(baseline, dict) or "sol" not in baseline or "P" not in baseline:
        raise ValueError(f"Cache baseline is missing sol/P objects: {path}")
    return baseline["sol"], baseline["P"]


def weighted_ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator > 1e-14 else math.nan


def model_age_profiles(sol: Any, P: Any) -> list[dict[str, float]]:
    g = np.asarray(sol.g, dtype=float)
    h_r = np.asarray(sol.hR_pol, dtype=float)
    if g.ndim != 7 or h_r.shape != g.shape:
        raise ValueError("Expected Markov-income g and hR_pol arrays with shape (b,ten,i,j,z,n,s).")
    ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    dep_last = int(P.n_child_stages)
    owner_rooms = np.asarray(P.H_own, dtype=float)
    rows: list[dict[str, float]] = []
    for j, age in enumerate(ages):
        gj = g[:, :, :, j, :, :, :]
        mass = float(np.sum(gj))
        owner_mass = float(np.sum(g[:, 1:, :, j, :, :, :]))
        child_mass = 0.0
        housing_sum = 0.0
        for nn in range(g.shape[5]):
            for cs in range(g.shape[6]):
                # m(n,s)>0 iff completed parity is positive and children are
                # still in a dependent stage.  Later completed-family states
                # therefore do not count as children currently at home.
                if nn > 0 and 1 <= cs <= dep_last:
                    child_mass += float(np.sum(g[:, :, :, j, :, nn, cs]))
        housing_sum += float(np.sum(g[:, 0, :, j, :, :, :] * h_r[:, 0, :, j, :, :, :]))
        for tenure in range(1, 1 + int(P.n_house)):
            housing_sum += float(np.sum(g[:, tenure, :, j, :, :, :])) * float(owner_rooms[tenure - 1])
        rows.append(
            {
                "age": float(age),
                "model_owner_rate": weighted_ratio(owner_mass, mass),
                "model_current_child_rate": weighted_ratio(child_mass, mass),
                "model_mean_rooms": weighted_ratio(housing_sum, mass),
                "model_population_mass": mass,
            }
        )
    # The final Bellman age is a terminal-value state without a current housing
    # decision.  Keep the standard figure on decision ages through J-1.
    return rows[:-1]


def acs_age_profiles(path: Path, ownership_path: Path) -> list[dict[str, float]]:
    accum: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                age = int(float(row["age"]))
                weight = float(row["pop_weight"])
                current_child = float(row["has_child_u18_rate"])
                rooms = float(row["mean_rooms"])
            except (KeyError, TypeError, ValueError):
                continue
            values = np.array([weight, current_child, rooms], dtype=float)
            if not np.all(np.isfinite(values)) or weight <= 0.0:
                continue
            item = accum.setdefault(age, {"weight": 0.0, "child": 0.0, "rooms": 0.0})
            item["weight"] += weight
            item["child"] += weight * current_child
            item["rooms"] += weight * rooms
    ownership: dict[int, float] = {}
    with ownership_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                row.get("source") != "ACS"
                or row.get("sample") != "household_heads_hhwt_due_housing"
                or not row.get("owner_rate")
            ):
                continue
            ownership[int(float(row["age"]))] = float(row["owner_rate"])
    rows: list[dict[str, float]] = []
    for age, item in sorted(accum.items()):
        weight = item["weight"]
        rows.append(
            {
                "age": float(age),
                "acs_owner_rate": ownership.get(age, math.nan),
                "acs_current_child_u18_rate": weighted_ratio(item["child"], weight),
                "acs_mean_rooms": weighted_ratio(item["rooms"], weight),
                "acs_population_weight": weight,
            }
        )
    if not rows:
        raise ValueError(f"No usable ACS age-profile rows found in {path}")
    return rows


def decision_rule_rows(sol: Any, P: Any, requested_age: float) -> list[dict[str, float | int]]:
    ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    j = int(np.argmin(np.abs(ages - requested_age)))
    b_grid = np.asarray(sol.b_grid, dtype=float)
    z_grid = np.asarray(getattr(sol, "type_values", P.z_grid), dtype=float)
    fert_probs = np.asarray(sol.fert_probs, dtype=float)
    tenure_choice = np.asarray(sol.tenure_choice, dtype=int)
    h_r = np.asarray(sol.hR_pol, dtype=float)
    owner_rooms = np.asarray(P.H_own, dtype=float)
    if int(P.n_parity) != FAMILY_SIZE_VALUES.size:
        raise ValueError("The established decision-rule labels require exactly three family-size states.")
    rows: list[dict[str, float | int]] = []
    for zz, z in enumerate(z_grid):
        for bb, wealth in enumerate(b_grid):
            probabilities = np.nan_to_num(fert_probs[bb, 0, 0, j, zz, :], nan=0.0, posinf=0.0, neginf=0.0)
            probability_sum = float(np.sum(probabilities))
            valid = probability_sum > 1e-14
            if valid:
                probabilities = probabilities / probability_sum
            else:
                probabilities = np.zeros_like(probabilities)
            branch_rooms: list[float] = []
            branch_tenures: list[int] = []
            for family_index in range(int(P.n_parity)):
                nn, cs = (0, 0) if family_index == 0 else (family_index, 1)
                target_tenure = int(tenure_choice[bb, 0, 0, j, zz, nn, cs])
                branch_tenures.append(target_tenure)
                rooms = (
                    float(h_r[bb, 0, 0, j, zz, nn, cs])
                    if target_tenure == 0
                    else float(owner_rooms[target_tenure - 1])
                )
                branch_rooms.append(rooms)
            expected_rooms = float(np.dot(probabilities, branch_rooms)) if valid else math.nan
            expected_family_size = float(np.dot(probabilities, FAMILY_SIZE_VALUES)) if valid else math.nan
            row: dict[str, float | int] = {
                "age": float(ages[j]),
                "income_state": float(z),
                "liquid_wealth": float(wealth),
                "valid_state": int(valid),
                "expected_physical_rooms_after_current_choices": expected_rooms,
                "expected_family_size": expected_family_size,
                "family_size_probability_sum": probability_sum,
            }
            for family_index, (family_size, probability, tenure, rooms) in enumerate(
                zip(FAMILY_SIZE_VALUES, probabilities, branch_tenures, branch_rooms)
            ):
                prefix = f"family_size_n{family_index}"
                row[f"{prefix}_children"] = float(family_size)
                row[f"{prefix}_probability"] = float(probability)
                row[f"{prefix}_tenure_choice"] = int(tenure)
                row[f"{prefix}_rooms"] = float(rooms)
            rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_lifecycle(model: list[dict[str, float]], acs: list[dict[str, float]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    plt.rcdefaults()
    navy = "#17365D"
    orange = "#D9791F"
    green = "#6E9D75"
    black = "#222222"
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 13,
            "axes.titlesize": 16,
            "axes.labelsize": 13,
            "legend.fontsize": 12,
        }
    )
    ma = np.asarray([row["age"] for row in model])
    acs = [row for row in acs if float(row["age"]) <= float(np.max(ma))]
    da = np.asarray([row["age"] for row in acs])
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.25))

    ax = axes[0]
    ax.axvspan(18.0, 42.0, color="#EFF1F4", zorder=0)
    ax.plot(da, [row["acs_owner_rate"] for row in acs], color=black, ls=(0, (2, 2)), lw=2.4, label="ACS data")
    ax.plot(ma, [row["model_owner_rate"] for row in model], color=navy, marker="o", ms=6.5, lw=2.8, label="Model")
    ax.set_title("Homeownership over the life cycle")
    ax.set_xlabel("Age")
    ax.set_ylabel("Homeownership rate")
    ax.set_ylim(0.0, 1.03)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)

    ax = axes[1]
    ax.axvspan(18.0, 42.0, color="#EFF1F4", zorder=0)
    child_model = ax.plot(
        ma,
        [row["model_current_child_rate"] for row in model],
        color=green,
        marker="o",
        ms=5.5,
        lw=2.6,
        label="Children at home, model",
    )[0]
    child_acs = ax.plot(
        da,
        [row["acs_current_child_u18_rate"] for row in acs],
        color=green,
        ls=(0, (2, 2)),
        lw=2.4,
        label="Own child under 18, ACS",
    )[0]
    ax.set_title("Children at home and housing space")
    ax.set_xlabel("Age")
    ax.set_ylabel("Households with children at home")
    ax.set_ylim(0.0, 0.75)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="y", alpha=0.25)
    ax2 = ax.twinx()
    rooms_acs = ax2.plot(
        da,
        [row["acs_mean_rooms"] for row in acs],
        color=black,
        ls=(0, (2, 2)),
        lw=2.4,
        label="Rooms, ACS",
    )[0]
    rooms_model = ax2.plot(
        ma,
        [row["model_mean_rooms"] for row in model],
        color=orange,
        marker="o",
        ms=5.5,
        lw=2.8,
        label="Rooms, model",
    )[0]
    ax2.set_ylabel("Mean rooms", color=orange)
    ax2.set_ylim(2.5, 6.5)
    ax.axvline(42.0, color="#777777", ls="--", lw=1.5)
    ax.legend(
        [child_model, child_acs, rooms_acs, rooms_model],
        [line.get_label() for line in (child_model, child_acs, rooms_acs, rooms_model)],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=2,
    )
    axes[0].text(-0.17, 1.04, "(a)", transform=axes[0].transAxes, fontsize=16)
    axes[1].text(-0.14, 1.04, "(b)", transform=axes[1].transAxes, fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.95, top=0.90, bottom=0.22, wspace=0.20)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_decision_rules(
    rows: list[dict[str, float | int]], path: Path, wealth_min: float, wealth_max: float
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcdefaults()
    z_values = sorted({float(row["income_state"]) for row in rows})
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 5.9), sharex=True)
    colors = ["#355F82", "#D77A20", "#4C8C59", "#8061A7", "#A9413B"]
    for color, z in zip(colors, z_values):
        group = sorted(
            (
                row
                for row in rows
                if float(row["income_state"]) == z
                and int(row["valid_state"]) == 1
                and wealth_min <= float(row["liquid_wealth"]) <= wealth_max
            ),
            key=lambda row: float(row["liquid_wealth"]),
        )
        x = [float(row["liquid_wealth"]) for row in group]
        axes[0].plot(
            x,
            [float(row["expected_physical_rooms_after_current_choices"]) for row in group],
            lw=2.0,
            color=color,
            label=rf"Income state $z={z:.6g}$",
        )
        axes[1].plot(
            x,
            [float(row["expected_family_size"]) for row in group],
            lw=2.0,
            color=color,
        )
    axes[0].set_title("(a) Expected physical housing after current choices")
    axes[0].set_ylabel("Rooms")
    axes[0].legend(frameon=False, loc="upper left", fontsize=10)
    axes[1].set_title("(b) Expected family size")
    axes[1].set_ylabel("Expected children (3 denotes 3+)")
    axes[1].set_ylim(-0.08, 3.05)
    for ax in axes:
        ax.set_xlabel("Liquid wealth (model units)")
        ax.set_xlim(wealth_min, wealth_max)
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    sol, P = load_trusted_solution_cache(args.solution_cache)
    model = model_age_profiles(sol, P)
    acs = acs_age_profiles(args.acs_csv, args.acs_ownership_csv)
    decision = decision_rule_rows(sol, P, args.decision_age)
    lifecycle_csv = args.lifecycle_csv or args.lifecycle_out.with_suffix(".csv")
    decision_csv = args.decision_rules_csv or args.decision_rules_out.with_suffix(".csv")
    lifecycle_rows = [{"source": "model", **row} for row in model] + [
        {"source": "ACS_2012_2023", **row} for row in acs
    ]
    write_rows(lifecycle_csv, lifecycle_rows)
    write_rows(decision_csv, decision)
    plot_lifecycle(model, acs, args.lifecycle_out)
    plot_decision_rules(decision, args.decision_rules_out, args.decision_wealth_min, args.decision_wealth_max)
    print(f"wrote {args.lifecycle_out}")
    print(f"wrote {lifecycle_csv}")
    print(f"wrote {args.decision_rules_out}")
    print(f"wrote {decision_csv}")


if __name__ == "__main__":
    main()
