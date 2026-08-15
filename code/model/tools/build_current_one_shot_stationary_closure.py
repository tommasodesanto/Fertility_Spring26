#!/usr/bin/env python3
"""Stationary population closure around the circulated one-shot household model.

The household problem is never changed.  At each candidate house price the
existing Markov-income solver returns normalized entrant flow E, mature-child
flow B, and housing demand D.  An external outer loop then solves either

    S = M / (E - rho B),       S D = Hs                         (open)

or

    B / E = 1,                S = Hs / D                       (closed).

This file is needed because the current Markov-income equilibrium routine does
not call the package's renewal-scale helper.  It is a stationary exercise, not
a calendar-time transition or a recalibration.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import numpy as np

from intergen_housing_fertility_optimized import solver as model
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.calibration import (
    diagnostic_loss,
    extract_moments,
)
from run_intergen_funded_policy_with_entry import solve_joint_case
from run_intergen_funded_property_tax_test import solve_balanced_case


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/current_one_shot_stationary_closure"
EXPECTED_SOURCE_SHA256 = "a2bba465b93c27cd358be93917a59af995ebdb2cd413520305d83a9014b1cfff"
EXPECTED_SOLVER_SHA256 = "4edf76ea08d07e4e1af7f98914afe1d7f7d1602bc805466b5847577d0cf92ce7"
FUNDED_HELPERS = {
    ROOT / "code/model/tools/run_intergen_funded_policy_with_entry.py":
        "53267a1f78d9361fe93f29bf9b551ccb8214aa671bb31d33310fa00518852c64",
    ROOT / "code/model/tools/run_intergen_funded_property_tax_test.py":
        "fe99765bfffd493845c77033e5f2406797f972b0aff2eeb1f271c59f4c6fa73b",
}

ROOT_TOL = 2.5e-5
CLOSED_TOL = 5e-7
BASELINE_NESTING_TOL = 5e-5


@dataclass(frozen=True)
class PolicyCase:
    name: str
    label: str
    description: str
    overrides: dict[str, Any]


@dataclass
class PriceReadout:
    price: float
    entry_E: float
    mature_B: float
    reproduction_B_over_E: float
    housing_demand_per_adult: float
    housing_supply: float
    tfr_equivalent: float
    births_per_adult: float
    owner_rate: float
    total_mass: float
    normalized_market_residual: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run the exact full grid for baseline and new-parent credit only.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_source(path: Path) -> tuple[dict[str, Any], dict[str, Any], str]:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            "Selected-result fingerprint changed: "
            f"expected {EXPECTED_SOURCE_SHA256}, observed {digest}."
        )
    payload = json.loads(raw)
    if payload.get("profile") != NEW_MOMENT_PROFILE_NAME:
        raise ValueError(
            f"Expected profile {NEW_MOMENT_PROFILE_NAME!r}; found {payload.get('profile')!r}."
        )
    selected = payload.get("selected")
    if not isinstance(selected, dict) or selected.get("status") != "ok":
        raise ValueError(f"No successful selected result in {path}.")
    if not isinstance(selected.get("theta"), dict):
        raise ValueError("Selected result has no parameter vector.")
    return payload, selected, digest


def policy_cases(base_parameters: Any, smoke: bool) -> list[PolicyCase]:
    cases = [
        PolicyCase(
            "baseline",
            "Baseline",
            "Saved circulated one-shot household estimate",
            {},
        ),
        PolicyCase(
            "new_parent_ltv95",
            "Dependent-child LTV 95",
            "Financed share rises from 80 to 95 percent while children are dependent",
            {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 0.95,
                "parent_dp_waiver_birth_state_only": True,
                "parent_dp_waiver_locations": np.array([], dtype=int),
                "parent_dp_waiver_owner_rungs": np.array([], dtype=int),
            },
        ),
        PolicyCase(
            "housing_supply_plus20",
            "Housing supply +20%",
            "Housing-supply intercept rises by 20 percent",
            {"H0": 1.20 * np.asarray(base_parameters.H0, dtype=float)},
        ),
        PolicyCase(
            "selling_cost_half",
            "Selling cost halved",
            "The proportional cost of selling an owned home is halved",
            {"psi": 0.50 * float(base_parameters.psi)},
        ),
    ]
    return cases[:2] if smoke else cases


def readout(solution: Any, price: float) -> PriceReadout:
    entry = float(solution.entry_rate)
    mature = float(solution.entrants_mature_total)
    demand = float(np.sum(np.asarray(solution.housing_demand, dtype=float)))
    supply = float(np.sum(np.asarray(solution.housing_supply, dtype=float)))
    return PriceReadout(
        price=float(price),
        entry_E=entry,
        mature_B=mature,
        reproduction_B_over_E=mature / max(entry, 1e-15),
        housing_demand_per_adult=demand,
        housing_supply=supply,
        tfr_equivalent=2.0 * float(solution.mean_parity),
        births_per_adult=float(solution.total_births_kfe),
        owner_rate=float(solution.own_rate),
        total_mass=float(solution.total_mass),
        normalized_market_residual=(demand - supply) / max(abs(supply), 1e-15),
    )


class CaseEvaluator:
    def __init__(
        self,
        case: PolicyCase,
        parameters: Any,
        b_grid: np.ndarray,
        normalized_solution: Any,
        normalized_price: float,
    ) -> None:
        self.case = case
        self.parameters = copy.deepcopy(parameters)
        self.b_grid = np.asarray(b_grid, dtype=float)
        self.shared = model.precompute_shared(self.parameters, self.b_grid)
        self.cache: dict[float, PriceReadout] = {}
        self.solve_count = 0
        self.cache[round(float(normalized_price), 12)] = readout(
            normalized_solution, float(normalized_price)
        )

    def evaluate(self, price: float) -> PriceReadout:
        price = float(price)
        key = round(price, 12)
        if key not in self.cache:
            solution = model.solve_markov_income_at_prices(
                np.array([price]),
                self.parameters,
                self.b_grid,
                verbose=False,
                fast_stats=False,
                SD=self.shared,
            )
            self.solve_count += 1
            self.cache[key] = readout(solution, price)
            print(
                f"  {self.case.name}: p={price:.7f} "
                f"B/E={self.cache[key].reproduction_B_over_E:.7f}",
                flush=True,
            )
        return self.cache[key]


def solve_decreasing_root(
    function: Callable[[float], float],
    start: float,
    *,
    tolerance: float,
    max_iterations: int = 18,
) -> tuple[float, float, int]:
    """Solve a locally decreasing scalar equation with a safeguarded bracket."""

    start = float(start)
    f_start = float(function(start))
    if not np.isfinite(f_start):
        f_start = 1e6
    if abs(f_start) <= tolerance:
        return start, f_start, 0

    lower = upper = start
    f_lower = f_upper = f_start
    expansions = 0
    if f_start > 0.0:
        for _ in range(20):
            upper *= 1.16
            f_upper = float(function(upper))
            expansions += 1
            if np.isfinite(f_upper) and f_upper <= 0.0:
                break
    else:
        for _ in range(20):
            lower /= 1.16
            f_lower = float(function(lower))
            expansions += 1
            if not np.isfinite(f_lower):
                f_lower = 1e6
            if f_lower >= 0.0:
                break
    if not (f_lower >= 0.0 and f_upper <= 0.0):
        raise RuntimeError(
            "Could not bracket decreasing root: "
            f"lower=({lower:.6g},{f_lower:.6g}), upper=({upper:.6g},{f_upper:.6g})."
        )

    midpoint = 0.5 * (lower + upper)
    f_mid = float(function(midpoint))
    for iteration in range(1, max_iterations + 1):
        midpoint = 0.5 * (lower + upper)
        f_mid = float(function(midpoint))
        if not np.isfinite(f_mid):
            f_mid = 1e6
        if abs(f_mid) <= tolerance:
            return midpoint, f_mid, expansions + iteration
        if f_mid > 0.0:
            lower, f_lower = midpoint, f_mid
        else:
            upper, f_upper = midpoint, f_mid
    return midpoint, f_mid, expansions + max_iterations


def open_row(
    evaluator: CaseEvaluator,
    case: PolicyCase,
    outside_share: float,
    retention: float,
    baseline_E: float,
    start_price: float,
) -> dict[str, Any]:
    outside_flow = float(outside_share * baseline_E)

    def residual(price: float) -> float:
        item = evaluator.evaluate(price)
        denominator = item.entry_E - retention * item.mature_B
        if denominator <= 1e-12:
            return 1e6
        scale = outside_flow / denominator
        return (
            scale * item.housing_demand_per_adult - item.housing_supply
        ) / max(item.housing_supply, 1e-15)

    price, market_residual, iterations = solve_decreasing_root(
        residual,
        start_price,
        tolerance=ROOT_TOL,
    )
    item = evaluator.evaluate(price)
    denominator = item.entry_E - retention * item.mature_B
    scale = outside_flow / denominator
    renewal_residual = scale * item.entry_E - (
        outside_flow + retention * scale * item.mature_B
    )
    return {
        "case": case.name,
        "label": case.label,
        "closure": "open_stationary",
        "outside_share_of_baseline_required_entrants": outside_share,
        "retention_rho": retention,
        "outside_flow_M": outside_flow,
        "price": price,
        "population_scale": scale,
        "tfr_equivalent": item.tfr_equivalent,
        "births_per_adult": item.births_per_adult,
        "total_births": scale * item.births_per_adult,
        "owner_rate": item.owner_rate,
        "entry_E_per_adult": item.entry_E,
        "mature_B_per_adult": item.mature_B,
        "reproduction_B_over_E": item.reproduction_B_over_E,
        "renewal_denominator_E_minus_rhoB": denominator,
        "housing_demand_per_adult": item.housing_demand_per_adult,
        "housing_supply": item.housing_supply,
        "market_residual": market_residual,
        "renewal_identity_residual": renewal_residual,
        "root_iterations_including_expansions": iterations,
    }


def closed_row(
    evaluator: CaseEvaluator,
    case: PolicyCase,
    start_price: float,
    current_housing_stock: float,
) -> dict[str, Any]:
    def residual(price: float) -> float:
        return evaluator.evaluate(price).reproduction_B_over_E - 1.0

    price, reproduction_residual, iterations = solve_decreasing_root(
        residual,
        start_price,
        tolerance=CLOSED_TOL,
        max_iterations=22,
    )
    item = evaluator.evaluate(price)
    static_scale = item.housing_supply / max(item.housing_demand_per_adult, 1e-15)
    fixed_stock_scale = current_housing_stock / max(
        item.housing_demand_per_adult, 1e-15
    )
    return {
        "case": case.name,
        "label": case.label,
        "closure": "closed_stationary_endpoint",
        "price": price,
        "price_ratio_to_current": price / start_price,
        "population_scale_static_supply": static_scale,
        "population_scale_fixed_current_stock": fixed_stock_scale,
        "tfr_equivalent": item.tfr_equivalent,
        "births_per_adult": item.births_per_adult,
        "total_births_static_supply": static_scale * item.births_per_adult,
        "total_births_fixed_current_stock": fixed_stock_scale * item.births_per_adult,
        "owner_rate": item.owner_rate,
        "entry_E_per_adult": item.entry_E,
        "mature_B_per_adult": item.mature_B,
        "reproduction_B_over_E": item.reproduction_B_over_E,
        "reproduction_root_residual": reproduction_residual,
        "housing_demand_per_adult": item.housing_demand_per_adult,
        "housing_supply": item.housing_supply,
        "root_iterations_including_expansions": iterations,
    }


def add_changes(
    rows: list[dict[str, Any]],
    group_keys: tuple[str, ...],
    value_keys: tuple[str, ...],
    *,
    baseline_case: str = "baseline",
) -> None:
    groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in rows:
        group = tuple(row[key] for key in group_keys)
        if row["case"] == baseline_case:
            groups[group] = row
    for row in rows:
        group = tuple(row[key] for key in group_keys)
        baseline = groups[group]
        for key in value_keys:
            base = float(baseline[key])
            row[f"{key}_change_percent"] = 100.0 * (float(row[key]) / base - 1.0)


def make_figures(
    evaluators: dict[str, CaseEvaluator],
    baseline_price: float,
    baseline_E: float,
    baseline_B: float,
    open_rows: list[dict[str, Any]],
    closed_rows: list[dict[str, Any]],
    outdir: Path,
    smoke: bool,
) -> list[dict[str, Any]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    baseline_eval = evaluators["baseline"]
    points = 7 if smoke else 15
    prices = np.linspace(0.38 * baseline_price, 1.35 * baseline_price, points)
    minimum_share = 1.0 - baseline_B / baseline_E
    outside_flow = minimum_share * baseline_E
    schedule_rows: list[dict[str, Any]] = []
    for price in prices:
        item = baseline_eval.evaluate(float(price))
        denominator = item.entry_E - item.mature_B
        schedule_rows.append(
            {
                "price": float(price),
                "price_index": float(price / baseline_price),
                "reproduction_B_over_E": item.reproduction_B_over_E,
                "housing_capacity_scale": item.housing_supply
                / max(item.housing_demand_per_adult, 1e-15),
                "renewal_scale_minimum_inflow": (
                    outside_flow / denominator if denominator > 1e-12 else np.nan
                ),
            }
        )

    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "font.size": 11})
    navy, orange, gray = "#17365D", "#C45A3B", "#777777"
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.25), constrained_layout=True)
    x = np.asarray([row["price_index"] for row in schedule_rows])
    axes[0].plot(
        x,
        [row["reproduction_B_over_E"] for row in schedule_rows],
        color=navy,
        linewidth=2.0,
    )
    axes[0].axhline(1.0, color=gray, linewidth=0.8, linestyle="--")
    axes[0].scatter([1.0], [baseline_B / baseline_E], color=orange, zorder=4)
    root = next(row for row in closed_rows if row["case"] == "baseline")
    axes[0].scatter(
        [root["price"] / baseline_price],
        [root["reproduction_B_over_E"]],
        color=navy,
        zorder=4,
    )
    axes[0].set(
        title="A. Fertility and replacement",
        xlabel="House price (current = 1)",
        ylabel="Mature households / required entrants",
    )

    capacity = np.asarray([row["housing_capacity_scale"] for row in schedule_rows])
    renewal = np.asarray([row["renewal_scale_minimum_inflow"] for row in schedule_rows])
    renewal[(renewal <= 0.0) | (renewal > 3.0)] = np.nan
    axes[1].plot(x, capacity, color=navy, linewidth=2.0, label="Housing capacity")
    axes[1].plot(x, renewal, color=orange, linewidth=2.0, label="Renewal scale")
    axes[1].scatter([1.0], [1.0], color="#222222", zorder=4)
    axes[1].set(
        title="B. Open steady state",
        xlabel="House price (current = 1)",
        ylabel="Population scale",
    )
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.grid(color="#E0E0E0", linewidth=0.5)
        axis.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"stationary_closure_schedules.{suffix}", dpi=220)
    plt.close(fig)

    policies = [name for name in evaluators if name != "baseline"]
    shares = sorted(
        {
            float(row["outside_share_of_baseline_required_entrants"])
            for row in open_rows
        }
    )
    selected_shares = [shares[0], shares[-1]]
    labels = [evaluators[name].case.label for name in policies]
    y = np.arange(len(policies))
    width = 0.36
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.7), constrained_layout=True)
    keys = (
        ("price_change_percent", "A. House price"),
        ("population_scale_change_percent", "B. Population"),
        ("total_births_change_percent", "C. Total births"),
    )
    colors = (navy, orange)
    for axis, (key, title) in zip(axes, keys):
        for index, share in enumerate(selected_shares):
            values = [
                next(
                    float(row[key])
                    for row in open_rows
                    if row["case"] == policy
                    and abs(
                        float(row["outside_share_of_baseline_required_entrants"])
                        - share
                    )
                    < 1e-10
                )
                for policy in policies
            ]
            axis.barh(
                y + (index - 0.5) * width,
                values,
                height=width,
                color=colors[index],
                label=(
                    "Full local retention"
                    if index == 0
                    else f"Outside share {100 * share:.0f}%"
                ),
            )
        axis.axvline(0.0, color="#444444", linewidth=0.7)
        axis.set_title(title, loc="left")
        axis.set_xlabel("Percent change")
        axis.grid(axis="x", color="#E0E0E0", linewidth=0.5)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
    axes[0].set_yticks(y, labels)
    axes[1].set_yticks(y, ["" for _ in y])
    axes[2].set_yticks(y, ["" for _ in y])
    axes[0].legend(frameon=False, fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"policy_closure_comparison.{suffix}", dpi=220)
    plt.close(fig)
    return schedule_rows


def write_readme(
    outdir: Path,
    summary: dict[str, Any],
    open_rows: list[dict[str, Any]],
    closed_rows: list[dict[str, Any]],
) -> None:
    baseline = summary["baseline"]
    closed = next(row for row in closed_rows if row["case"] == "baseline")
    shares = sorted(
        {
            float(row["outside_share_of_baseline_required_entrants"])
            for row in open_rows
        }
    )
    low_share = shares[0]
    funded = summary.get("funded_policy")
    lines = [
        "# Stationary population closure for the one-shot model",
        "",
        "This packet leaves the household model unchanged. An outer scalar root jointly "
        "clears stationary renewal and housing. It is not a transition or a recalibration.",
        "",
        "## Baseline",
        "",
        f"- Current price: `{baseline['price']:.8f}`.",
        f"- E: `{baseline['entry_E']:.12f}`; B: `{baseline['mature_B']:.12f}`; "
        f"B/E: `{baseline['reproduction_B_over_E']:.9f}`.",
        f"- Minimum outside share with full local retention: `{100 * low_share:.3f}%`.",
        f"- Closed root price: `{closed['price']:.8f}` "
        f"(`{100 * closed['price_ratio_to_current']:.2f}%` of current).",
        f"- Closed scale: `{closed['population_scale_static_supply']:.5f}` under the "
        "static supply curve and "
        f"`{closed['population_scale_fixed_current_stock']:.5f}` at fixed current stock.",
        "",
        "The closed scales are endpoint diagnostics, not forecasts. The saved "
        "one-shot estimate is held fixed for this mechanism calculation; its original "
        "first-child housing-response target has been withdrawn.",
        "",
        "## Policy readout",
        "",
        "`open_stationary.csv` reports every policy under the minimum-inflow closure, "
        "a 10% outside share, and a 20% outside share. The larger shares are sensitivity "
        "cases: local retention is chosen to make the same baseline scale stationary.",
        "",
    ]
    if isinstance(funded, dict):
        lines.extend(
            [
                "## Standing paper policies",
                "",
                "The purchase grant and the property-tax-plus-grant package are rerun with "
                "the original fiscal rule and the same stationary renewal closure. The tax "
                "base is all occupied housing, rentals and owners. Population is adult-household "
                "mass, and total births are a four-year stationary flow. Results "
                "are in `funded_policy_open_stationary.csv`; the corresponding fixed-population "
                "comparison is in `funded_policy_fixed_population.csv`.",
                f"The saved unrebated fixed-theta loss is "
                f"`{funded['baseline']['saved_unrebated_fixed_theta_loss']:.6f}`; "
                f"the rebated 1% tax reference has fixed-theta loss "
                f"`{funded['baseline']['rebated_tax_fixed_theta_loss']:.6f}`. "
                "It is a counterfactual reference, not a recalibration.",
                "",
            ]
        )
    lines.extend(
        [
        "## Files",
        "",
        "- `normalized_equilibrium.csv`: the old fixed-population comparison.",
        "- `open_stationary.csv`: price and population jointly determined.",
        "- `closed_endpoints.csv`: B/E=1 roots and implied scales.",
        "- `stationary_closure_schedules.pdf`: the two equilibrium schedules.",
        "- `policy_closure_comparison.pdf`: policy effects under two openness assumptions.",
        "- `funded_policy_closure_comparison.pdf`: the standing paper policies under the new closure.",
        "- `target_fit.csv` and `parameters.csv`: full calibration audit tables.",
        "",
        "## Regeneration",
        "",
        "```bash",
        "PYTHONPATH=code/model NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "code/model/.venv/bin/python "
        "code/model/tools/build_current_one_shot_stationary_closure.py",
        "```",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def renewal_scale_function(
    solution: Any,
    parameters: Any,
    outside: dict[str, float],
) -> tuple[SimpleNamespace, list[dict[str, Any]]]:
    """Return the fixed-inflow renewal scale expected by the funded-policy solver."""

    reference_population = max(float(solution.total_mass), 1e-14)
    entry_per_scale = float(solution.entry_rate) / reference_population
    mature_per_scale = float(solution.entrants_mature_total) / reference_population
    retention = float(outside["retention_rho"])
    outside_flow = float(outside["outside_flow_M"])
    denominator = entry_per_scale - retention * mature_per_scale
    finite = bool(denominator > 1e-12 and outside_flow >= 0.0)
    if finite:
        implied_population = outside_flow / denominator
        scale = implied_population / reference_population
        renewal_residual = (
            implied_population * entry_per_scale
            - outside_flow
            - retention * implied_population * mature_per_scale
        )
    else:
        implied_population = math.inf
        scale = math.inf
        renewal_residual = math.nan
    return (
        SimpleNamespace(
            finite=finite,
            scale_factor=float(scale),
            implied_total_population=float(implied_population),
            entry_per_scale=float(entry_per_scale),
            mature_per_scale=float(mature_per_scale),
            retention_rho=retention,
            outside_flow_M=outside_flow,
            denominator=float(denominator),
            renewal_residual=float(renewal_residual),
        ),
        [],
    )


def make_funded_policy_figure(
    outdir: Path,
    rows: list[dict[str, Any]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    policies = ["purchase_grant", "tax2_purchase_grant"]
    policy_labels = ["Purchase grant", "Tax + purchase grant"]
    shares = sorted(
        {
            float(row["outside_share_of_baseline_required_entrants"])
            for row in rows
        }
    )
    selected_shares = [shares[0], shares[-1]]
    y = np.arange(len(policies))
    width = 0.34
    navy, orange = "#17365D", "#C45A3B"
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 4.1), constrained_layout=True)
    fields = (
        ("price_change_percent", "A. House price"),
        ("population_scale_change_percent", "B. Population"),
        ("total_births_change_percent", "C. Total births"),
    )
    for axis, (field, title) in zip(axes, fields):
        for share_index, share in enumerate(selected_shares):
            values = [
                float(
                    next(
                        row[field]
                        for row in rows
                        if row["case"] == case
                        and abs(
                            float(row["outside_share_of_baseline_required_entrants"])
                            - share
                        )
                        < 1e-10
                    )
                )
                for case in policies
            ]
            offsets = y + (share_index - 0.5) * width
            bars = axis.barh(
                offsets,
                values,
                height=width,
                color=(navy, orange)[share_index],
                label=(
                    "Full local retention"
                    if share_index == 0
                    else f"Outside share {100 * share:.0f}%"
                ),
            )
            for bar, value in zip(bars, values):
                axis.annotate(
                    f"{value:+.2f}%",
                    (value, bar.get_y() + bar.get_height() / 2.0),
                    xytext=(4 if value >= 0 else -4, 0),
                    textcoords="offset points",
                    ha="left" if value >= 0 else "right",
                    va="center",
                    fontsize=9.0,
                )
        axis.axvline(0.0, color="#444444", linewidth=0.7)
        axis.set_title(title, loc="left")
        axis.set_xlabel("Percent change")
        axis.grid(axis="x", color="#E0E0E0", linewidth=0.5)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
    axes[0].set_yticks(y, policy_labels)
    axes[1].set_yticks(y, ["" for _ in y])
    axes[2].set_yticks(y, ["" for _ in y])
    axes[0].legend(frameon=False, fontsize=9.2, loc="best")
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"funded_policy_closure_comparison.{suffix}", dpi=220)
    plt.close(fig)


def run_funded_policy_portfolio(
    base_overrides: dict[str, Any],
    outdir: Path,
    *,
    smoke: bool,
    saved_unrebated_loss: float,
) -> dict[str, Any]:
    """Rerun the standing paper policies under the same renewal closure."""

    cases = [
        ("funded_baseline", "Rebated 1% property tax", 0.01, False),
        ("purchase_grant", "Purchase grant", 0.01, True),
        ("tax2_purchase_grant", "2% property tax + purchase grant", 0.02, True),
    ]
    target_system = new_moment_target_system()
    fixed_rows: list[dict[str, Any]] = []
    fixed_objects: dict[str, tuple[Any, Any]] = {}
    for name, label, annual_tax, grant in cases:
        print(f"funded fixed-population case: {label}", flush=True)
        fixed, solution, parameters = solve_balanced_case(
            name,
            annual_tax,
            grant,
            smoke,
            base_overrides,
            target_system,
        )
        reference_mass = max(float(solution.total_mass), 1e-14)
        fixed_payload = {"case": name, "label": label, **fixed}
        prime_age_owner_rate = float(fixed_payload.pop("own_rate"))
        property_tax_revenue = float(fixed_payload.pop("property_tax_revenue"))
        grant_recipient_mass = float(fixed_payload.pop("grant_recipient_mass"))
        grant_outlays = float(fixed_payload.pop("grant_outlays"))
        transfer_outlays = float(fixed_payload.pop("transfer_outlays"))
        budget_residual = float(fixed_payload.pop("government_budget_residual"))
        fixed_payload.update(
            {
                "prime_age_owner_rate_3055": prime_age_owner_rate,
                "aggregate_owner_rate": float(solution.own_rate),
                "births_per_adult": float(solution.total_births_kfe) / reference_mass,
                "entry_E_per_adult": float(solution.entry_rate) / reference_mass,
                "mature_B_per_adult": float(solution.entrants_mature_total)
                / reference_mass,
                "housing_demand_per_adult": float(
                    np.asarray(solution.housing_demand, dtype=float).reshape(-1)[0]
                )
                / reference_mass,
                "property_tax_revenue_per_adult": property_tax_revenue
                / reference_mass,
                "grant_recipient_mass_per_adult": grant_recipient_mass
                / reference_mass,
                "grant_outlays_per_adult": grant_outlays / reference_mass,
                "transfer_outlays_per_adult": transfer_outlays / reference_mass,
                "government_budget_residual_per_adult": budget_residual
                / reference_mass,
            }
        )
        fixed_payload["reproduction_B_over_E"] = (
            fixed_payload["mature_B_per_adult"]
            / max(fixed_payload["entry_E_per_adult"], 1e-15)
        )
        fixed_rows.append(fixed_payload)
        fixed_objects[name] = (solution, parameters)
        write_csv(outdir / "funded_policy_fixed_population.csv", fixed_rows)

    fixed_baseline = fixed_rows[0]
    baseline_solution, _ = fixed_objects["funded_baseline"]
    baseline_population = float(baseline_solution.total_mass)
    baseline_E = float(baseline_solution.entry_rate)
    baseline_B = float(baseline_solution.entrants_mature_total)
    baseline_ratio = baseline_B / baseline_E
    minimum_share = 1.0 - baseline_ratio
    if minimum_share <= 0.0:
        raise RuntimeError("Funded baseline has no positive fixed-inflow renewal closure.")
    shares = [minimum_share, 0.20] if smoke else [minimum_share, 0.10, 0.20]
    outside_objects = []
    for share in shares:
        retention = (1.0 - share) / baseline_ratio
        outside_flow = share * baseline_E
        if retention > 1.0 + 1e-10:
            raise RuntimeError("Funded-policy sensitivity requires retention above one.")
        outside_objects.append(
            {
                "outside_share": float(share),
                "retention_rho": float(retention),
                "outside_flow_M": float(outside_flow),
            }
        )

    open_rows: list[dict[str, Any]] = []
    for fixed, (name, label, annual_tax, grant) in zip(fixed_rows, cases):
        for outside in outside_objects:
            print(
                f"funded renewal case: {label}; outside share={100 * outside['outside_share']:.2f}%",
                flush=True,
            )
            joint = solve_joint_case(
                base_overrides,
                outside,
                label=f"{name}_renewal_{outside['outside_share']:.6f}",
                annual_tax=annual_tax,
                grant=grant,
                initial_price=float(fixed["price"]),
                initial_transfer=float(fixed["lump_sum_transfer_period_units"]),
                smoke=smoke,
                population_scale_function=renewal_scale_function,
            )
            item = joint.solution
            moments = extract_moments(item, joint.parameters)
            loss = diagnostic_loss(
                moments,
                targets=target_system.targets_dict(),
                weights=target_system.weights_dict(),
            )
            scale = float(joint.scale.scale_factor)
            reference_mass = max(float(item.total_mass), 1e-14)
            births_per_adult = float(item.total_births_kfe) / reference_mass
            entry_per_adult = float(item.entry_rate) / reference_mass
            mature_per_adult = float(item.entrants_mature_total) / reference_mass
            housing_demand_per_adult = float(
                np.asarray(item.housing_demand, dtype=float).reshape(-1)[0]
            ) / reference_mass
            tax_revenue_per_adult = float(item.property_tax_revenue) / reference_mass
            grant_mass_per_adult = (
                float(item.birth_entry_grant_recipient_mass) / reference_mass
            )
            grant_outlays_per_adult = (
                float(item.birth_entry_grant_outlays) / reference_mass
            )
            transfer_outlays_per_adult = (
                float(item.property_tax_transfer_outlays) / reference_mass
            )
            budget_residual_per_adult = (
                float(item.property_tax_budget_residual) / reference_mass
            )
            open_rows.append(
                {
                    "case": name,
                    "label": label,
                    "outside_share_of_baseline_required_entrants": outside[
                        "outside_share"
                    ],
                    "retention_rho": outside["retention_rho"],
                    "outside_flow_M": outside["outside_flow_M"],
                    "price": float(joint.price),
                    "population_scale": scale,
                    "tfr_equivalent": float(moments["tfr"]),
                    "childless_rate": float(moments["childless_rate"]),
                    "prime_age_owner_rate_3055": float(moments["own_rate"]),
                    "aggregate_owner_rate": float(item.own_rate),
                    "old_age_owner_rate": float(moments["old_age_own_rate"]),
                    "housing_increment_0to1": float(
                        moments["housing_increment_0to1"]
                    ),
                    "fixed_theta_loss": float(loss),
                    "births_per_adult": births_per_adult,
                    "total_births": scale * births_per_adult,
                    "lump_sum_transfer_period_units": float(joint.transfer),
                    "lump_sum_transfer_annual_income_units": float(joint.transfer)
                    / float(joint.parameters.period_years),
                    "entry_E_per_adult": entry_per_adult,
                    "mature_B_per_adult": mature_per_adult,
                    "reproduction_B_over_E": mature_per_adult
                    / max(entry_per_adult, 1e-15),
                    "housing_demand_per_adult": housing_demand_per_adult,
                    "property_tax_revenue_per_adult": tax_revenue_per_adult,
                    "grant_recipient_mass_per_adult": grant_mass_per_adult,
                    "grant_outlays_per_adult": grant_outlays_per_adult,
                    "transfer_outlays_per_adult": transfer_outlays_per_adult,
                    "government_budget_residual_per_adult": budget_residual_per_adult,
                    "property_tax_revenue_aggregate": scale
                    * tax_revenue_per_adult,
                    "grant_recipient_mass_aggregate": scale
                    * grant_mass_per_adult,
                    "grant_outlays_aggregate": scale * grant_outlays_per_adult,
                    "transfer_outlays_aggregate": scale
                    * transfer_outlays_per_adult,
                    "government_budget_residual_aggregate": scale
                    * budget_residual_per_adult,
                    "market_residual": float(joint.relative_excess),
                    "fiscal_residual": float(joint.scaled_fiscal_residual),
                    "renewal_residual": float(joint.scale.renewal_residual),
                    "model_evaluations": int(joint.joint_model_evaluations),
                }
            )
            write_csv(outdir / "funded_policy_open_stationary.csv", open_rows)

    for row in fixed_rows:
        for field in (
            "price",
            "tfr",
            "births_per_adult",
            "prime_age_owner_rate_3055",
            "aggregate_owner_rate",
        ):
            row[f"{field}_change_percent"] = 100.0 * (
                float(row[field]) / float(fixed_baseline[field]) - 1.0
            )
    add_changes(
        open_rows,
        ("outside_share_of_baseline_required_entrants",),
        (
            "price",
            "population_scale",
            "tfr_equivalent",
            "total_births",
            "prime_age_owner_rate_3055",
            "aggregate_owner_rate",
        ),
        baseline_case="funded_baseline",
    )
    write_csv(outdir / "funded_policy_fixed_population.csv", fixed_rows)
    write_csv(outdir / "funded_policy_open_stationary.csv", open_rows)
    make_funded_policy_figure(outdir, open_rows)

    baseline_rows = [row for row in open_rows if row["case"] == "funded_baseline"]
    nesting_pairs = (
        ("price", "price"),
        ("lump_sum_transfer_period_units", "lump_sum_transfer_period_units"),
        (
            "lump_sum_transfer_annual_income_units",
            "lump_sum_transfer_annual_income_units",
        ),
        ("tfr_equivalent", "tfr"),
        ("childless_rate", "childless_rate"),
        ("prime_age_owner_rate_3055", "prime_age_owner_rate_3055"),
        ("aggregate_owner_rate", "aggregate_owner_rate"),
        ("old_age_owner_rate", "old_age_own_rate"),
        ("housing_increment_0to1", "housing_increment_0to1"),
        ("fixed_theta_loss", "fixed_theta_loss"),
        ("births_per_adult", "births_per_adult"),
        ("entry_E_per_adult", "entry_E_per_adult"),
        ("mature_B_per_adult", "mature_B_per_adult"),
        ("reproduction_B_over_E", "reproduction_B_over_E"),
        ("housing_demand_per_adult", "housing_demand_per_adult"),
        ("property_tax_revenue_per_adult", "property_tax_revenue_per_adult"),
        ("grant_recipient_mass_per_adult", "grant_recipient_mass_per_adult"),
        ("grant_outlays_per_adult", "grant_outlays_per_adult"),
        ("transfer_outlays_per_adult", "transfer_outlays_per_adult"),
        (
            "government_budget_residual_per_adult",
            "government_budget_residual_per_adult",
        ),
        ("market_residual", "market_residual"),
    )

    def reported_object_nesting(row: dict[str, Any]) -> float:
        differences = [abs(float(row["population_scale"]) - 1.0)]
        for open_key, fixed_key in nesting_pairs:
            reference = float(fixed_baseline[fixed_key])
            differences.append(
                abs(float(row[open_key]) - reference) / max(1.0, abs(reference))
            )
        differences.append(
            abs(
                float(row["government_budget_residual_aggregate"])
                - float(fixed_baseline["government_budget_residual_per_adult"])
            )
        )
        return max(differences)

    nesting = max(reported_object_nesting(row) for row in baseline_rows)
    max_market = max(abs(float(row["market_residual"])) for row in open_rows)
    max_fiscal = max(abs(float(row["fiscal_residual"])) for row in open_rows)
    max_renewal = max(abs(float(row["renewal_residual"])) for row in open_rows)
    funded_tolerance = 2.0e-4 if smoke else ROOT_TOL
    gates = {
        "baseline_reported_objects_nesting_max_abs": nesting,
        "max_market_residual": max_market,
        "max_fiscal_residual": max_fiscal,
        "max_renewal_residual": max_renewal,
        "passed": bool(
            nesting <= BASELINE_NESTING_TOL
            and max_market <= funded_tolerance
            and max_fiscal <= funded_tolerance
            and max_renewal <= 1e-10
        ),
    }
    if not gates["passed"]:
        raise RuntimeError(f"Funded-policy acceptance gates failed: {gates}")
    return {
        "baseline": {
            "price": float(fixed_baseline["price"]),
            "population": baseline_population,
            "entry_E": baseline_E,
            "mature_B": baseline_B,
            "reproduction_B_over_E": baseline_ratio,
            "minimum_outside_share": minimum_share,
            "saved_unrebated_fixed_theta_loss": float(saved_unrebated_loss),
            "rebated_tax_fixed_theta_loss": float(
                fixed_baseline["fixed_theta_loss"]
            ),
        },
        "outside_sensitivity": outside_objects,
        "gates": gates,
        "cases": [name for name, _, _, _ in cases],
        "population_unit": "adult-household mass",
        "birth_flow_unit": "births over one four-year model period",
        "property_tax_base": "all occupied housing, rentals and owners",
    }


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    source = args.source.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    payload, selected, source_sha = load_source(source)
    solver_path = Path(model.__file__).resolve()
    solver_sha = hashlib.sha256(solver_path.read_bytes()).hexdigest()
    if solver_sha != EXPECTED_SOLVER_SHA256:
        raise ValueError(
            "Solver fingerprint changed: "
            f"expected {EXPECTED_SOLVER_SHA256}, observed {solver_sha}."
        )

    theta = {str(key): float(value) for key, value in selected["theta"].items()}
    base_overrides = new_moment_overrides(tight=True, optimized=True)
    base_overrides.update(theta)
    helper_hashes: dict[str, str] = {}
    for helper_path, expected_hash in FUNDED_HELPERS.items():
        observed_hash = hashlib.sha256(helper_path.read_bytes()).hexdigest()
        if observed_hash != expected_hash:
            raise ValueError(
                f"Funded-policy helper fingerprint changed for {helper_path}: "
                f"expected {expected_hash}, observed {observed_hash}."
            )
        helper_hashes[str(helper_path)] = observed_hash

    print("Solving saved one-shot normalized baseline...", flush=True)
    baseline_solution, baseline_parameters, baseline_price_array = model.run_model_cp_dt(
        base_overrides, verbose=False
    )
    baseline_price = float(np.asarray(baseline_price_array).reshape(-1)[0])
    if not bool(getattr(baseline_solution, "timings", {}).get("strict_converged", False)):
        raise RuntimeError("Retained normalized baseline did not strictly converge.")
    if abs(baseline_price - float(selected["price"])) > 5e-8:
        raise RuntimeError("Retained normalized baseline price no longer reproduces the source.")
    baseline_item = readout(baseline_solution, baseline_price)
    baseline_E = baseline_item.entry_E
    baseline_B = baseline_item.mature_B
    baseline_ratio = baseline_item.reproduction_B_over_E
    minimum_share = 1.0 - baseline_ratio
    if minimum_share <= 0.0:
        raise RuntimeError("Baseline has no positive minimum outside-flow closure.")
    outside_shares = [minimum_share, 0.10, 0.20]
    if args.smoke:
        outside_shares = [minimum_share, 0.20]
    retention_by_share = {
        share: (1.0 - share) / baseline_ratio for share in outside_shares
    }
    if max(retention_by_share.values()) > 1.0 + 1e-10:
        raise RuntimeError("An open sensitivity requires local retention above one.")

    cases = policy_cases(baseline_parameters, args.smoke)
    normalized_rows: list[dict[str, Any]] = []
    open_rows: list[dict[str, Any]] = []
    closed_rows: list[dict[str, Any]] = []
    evaluators: dict[str, CaseEvaluator] = {}
    current_housing_stock = baseline_item.housing_supply

    for index, case in enumerate(cases, start=1):
        case_started = time.perf_counter()
        print(f"case {index}/{len(cases)}: {case.label}", flush=True)
        if case.name == "baseline":
            solution = baseline_solution
            parameters = baseline_parameters
            price = baseline_price
        else:
            overrides = dict(base_overrides)
            overrides.update(case.overrides)
            solution, parameters, price_array = model.run_model_cp_dt(
                overrides, verbose=False
            )
            price = float(np.asarray(price_array).reshape(-1)[0])
            if not bool(getattr(solution, "timings", {}).get("strict_converged", False)):
                raise RuntimeError(f"{case.name}: normalized equilibrium did not converge.")
        item = readout(solution, price)
        normalized_rows.append(
            {
                "case": case.name,
                "label": case.label,
                "description": case.description,
                "price": price,
                "population_scale": 1.0,
                "tfr_equivalent": item.tfr_equivalent,
                "births_per_adult": item.births_per_adult,
                "total_births": item.births_per_adult,
                "owner_rate": item.owner_rate,
                "entry_E_per_adult": item.entry_E,
                "mature_B_per_adult": item.mature_B,
                "reproduction_B_over_E": item.reproduction_B_over_E,
                "housing_demand_per_adult": item.housing_demand_per_adult,
                "housing_supply": item.housing_supply,
                "market_residual": item.normalized_market_residual,
            }
        )
        evaluator = CaseEvaluator(
            case,
            parameters,
            np.asarray(solution.b_grid, dtype=float),
            solution,
            price,
        )
        evaluators[case.name] = evaluator
        for share in outside_shares:
            open_rows.append(
                open_row(
                    evaluator,
                    case,
                    share,
                    retention_by_share[share],
                    baseline_E,
                    price,
                )
            )
        closed_rows.append(
            closed_row(
                evaluator,
                case,
                baseline_price,
                current_housing_stock,
            )
        )
        write_csv(outdir / "normalized_equilibrium.csv", normalized_rows)
        write_csv(outdir / "open_stationary.csv", open_rows)
        write_csv(outdir / "closed_endpoints.csv", closed_rows)
        write_json_atomic(
            outdir / "latest_completed_case.json",
            {
                "status": "running",
                "completed_case": case.name,
                "completed_index": index,
                "cases_requested": len(cases),
                "fixed_price_solves": evaluator.solve_count,
                "elapsed_seconds": time.perf_counter() - case_started,
            },
        )

    add_changes(
        normalized_rows,
        (),
        ("price", "tfr_equivalent", "births_per_adult", "total_births", "owner_rate"),
    )
    add_changes(
        open_rows,
        ("outside_share_of_baseline_required_entrants",),
        ("price", "population_scale", "tfr_equivalent", "total_births", "owner_rate"),
    )
    add_changes(
        closed_rows,
        (),
        (
            "price",
            "population_scale_static_supply",
            "population_scale_fixed_current_stock",
            "tfr_equivalent",
            "total_births_static_supply",
        ),
    )
    write_csv(outdir / "normalized_equilibrium.csv", normalized_rows)
    write_csv(outdir / "open_stationary.csv", open_rows)
    write_csv(outdir / "closed_endpoints.csv", closed_rows)
    write_csv(outdir / "target_fit.csv", list(selected.get("target_fit", [])))
    write_csv(outdir / "parameters.csv", list(selected.get("parameters", [])))
    schedule_rows = make_figures(
        evaluators,
        baseline_price,
        baseline_E,
        baseline_B,
        open_rows,
        closed_rows,
        outdir,
        args.smoke,
    )
    write_csv(outdir / "baseline_schedules.csv", schedule_rows)
    saved_unrebated_loss = sum(
        float(row.get("loss_contribution", 0.0))
        for row in selected.get("target_fit", [])
    )
    funded_policy = run_funded_policy_portfolio(
        base_overrides,
        outdir,
        smoke=bool(args.smoke),
        saved_unrebated_loss=saved_unrebated_loss,
    )

    max_open_market = max(abs(float(row["market_residual"])) for row in open_rows)
    max_open_renewal = max(
        abs(float(row["renewal_identity_residual"])) for row in open_rows
    )
    max_closed_reproduction = max(
        abs(float(row["reproduction_root_residual"])) for row in closed_rows
    )
    baseline_open_rows = [row for row in open_rows if row["case"] == "baseline"]
    baseline_nesting = max(
        max(abs(float(row["price"]) / baseline_price - 1.0), abs(float(row["population_scale"]) - 1.0))
        for row in baseline_open_rows
    )
    gates = {
        "baseline_open_price_and_scale_nesting_max_abs": baseline_nesting,
        "max_open_market_residual": max_open_market,
        "max_open_renewal_identity_residual": max_open_renewal,
        "max_closed_reproduction_root_residual": max_closed_reproduction,
        "passed": bool(
            baseline_nesting <= BASELINE_NESTING_TOL
            and max_open_market <= ROOT_TOL
            and max_open_renewal <= 1e-10
            and max_closed_reproduction <= CLOSED_TOL
        ),
    }
    if not gates["passed"]:
        raise RuntimeError(f"Stationary-closure acceptance gates failed: {gates}")

    elapsed = time.perf_counter() - started
    driver_path = Path(__file__).resolve()
    summary = {
        "status": "complete",
        "purpose": "paper-facing stationary closure around the saved circulated one-shot estimate",
        "not_a_transition": True,
        "not_a_recalibration": True,
        "smoke": bool(args.smoke),
        "source": str(source),
        "source_sha256": source_sha,
        "solver": str(solver_path),
        "solver_sha256": solver_sha,
        "funded_policy_helper_sha256": helper_hashes,
        "driver": str(driver_path),
        "driver_sha256": hashlib.sha256(driver_path.read_bytes()).hexdigest(),
        "profile": payload.get("profile"),
        "theta": theta,
        "baseline": {
            "price": baseline_price,
            "entry_E": baseline_E,
            "mature_B": baseline_B,
            "reproduction_B_over_E": baseline_ratio,
            "minimum_outside_share": minimum_share,
            "minimum_outside_flow_per_four_years": baseline_E - baseline_B,
            "housing_stock": current_housing_stock,
            "tfr_equivalent": baseline_item.tfr_equivalent,
        },
        "open_sensitivity": [
            {
                "outside_share": share,
                "retention_rho": retention_by_share[share],
                "outside_flow": share * baseline_E,
            }
            for share in outside_shares
        ],
        "cases": [jsonable(case.__dict__) for case in cases],
        "gates": gates,
        "fixed_price_solve_counts": {
            name: evaluator.solve_count for name, evaluator in evaluators.items()
        },
        "funded_policy": funded_policy,
        "elapsed_seconds": elapsed,
        "limitations": [
            "The saved one-shot estimate is diagnostic: its original first-child housing-response target was withdrawn.",
            "Outside inflow and local retention are sensitivity parameters, not calibrated targets.",
            "The closed endpoint is a steady state, not a forecast or a historical transition.",
            "A calendar-time U.S. path still needs an age-specific household-formation bridge.",
        ],
    }
    write_json_atomic(outdir / "summary.json", summary)
    write_readme(outdir, summary, open_rows, closed_rows)
    write_json_atomic(
        outdir / "latest_completed_case.json",
        {
            "status": "complete",
            "cases": [case.name for case in cases],
            "elapsed_seconds": elapsed,
            "summary": str(outdir / "summary.json"),
        },
    )
    print(
        f"complete: cases={len(cases)}, elapsed={elapsed:.1f}s, output={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
