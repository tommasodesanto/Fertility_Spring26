#!/usr/bin/env python3
"""Build a generic age-structured population--housing transition packet.

This is a deliberately stylized numerical counterpart to the paper's
reproductive-equilibrium theory.  It is not calibrated and it never imports or
changes the maintained quantitative model.  The state contains adult household
cohorts and four five-year child stages.  An end-of-period birth enters the
newborn stage next date and becomes an adult household after 20 years.  A
contemporaneous housing-cost index clears a single market in every period.

Outputs
-------
transition_paths.csv
    One row per scenario and five-year period, including cohort accounting,
    market clearing, and the frozen-environment demographic eigenvalue.
cohort_distributions.csv
    Adult cohort masses and housing demand by age group.
momentum_condition_grid.csv
    Exact two-generation population/price momentum conditions.
summary.json, README.md
    Parameters, formulas, checks, and selected numerical results.
transition_paths.{pdf,png}, transition_accounting.{pdf,png},
momentum_conditions.{pdf,png}
    Three paper-ready figure concepts.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = ROOT / "output" / "model" / "dynamic_population_housing_paper"


@dataclass(frozen=True)
class Parameters:
    period_years: int = 5
    periods: int = 20
    adult_age_min: int = 20
    adult_age_max: int = 79
    maturation_periods: int = 4
    child_queue_survival: float = 0.992
    household_formation_per_birth: float = 0.5
    supply_elasticity: float = 0.65
    demand_price_elasticity: float = 1.0
    child_goods_cost: float = 0.72
    child_space: float = 0.22
    access_wedge_baseline: float = 0.24
    theta_initial: float = 1.16
    theta_long_run: float = 0.94
    theta_decay_periods: float = 2.75
    policy_start_period: float = 2.0
    policy_phasein_periods: float = 2.0
    young_access_reduction: float = 0.12
    old_retention_reduction: float = 0.18
    old_age_cutoff: int = 60


@dataclass(frozen=True)
class Scenario:
    name: str
    label: str
    access_reduction: float = 0.0
    old_retention_reduction: float = 0.0
    freeze_fertility_schedule: bool = False


SCENARIOS = (
    Scenario("baseline", "Baseline"),
    Scenario(
        "no_decline",
        "No fertility decline",
        freeze_fertility_schedule=True,
    ),
    Scenario("young_access", "Easier young access", access_reduction=1.0),
    Scenario("old_release", "Lower old retention", old_retention_reduction=1.0),
)


def age_grid(par: Parameters) -> np.ndarray:
    return np.arange(par.adult_age_min, par.adult_age_max + 1, par.period_years)


def adult_survival(ages: np.ndarray) -> np.ndarray:
    """Five-year survival/continuation probabilities for adult households."""
    rates = np.array(
        [0.995, 0.994, 0.992, 0.989, 0.985, 0.978,
         0.966, 0.946, 0.915, 0.870, 0.790, 0.0],
        dtype=float,
    )
    if rates.size != ages.size:
        raise ValueError("The survival schedule must match the adult age grid")
    return rates


def fertility_profile(ages: np.ndarray) -> np.ndarray:
    """Relative birth intensity by five-year age group."""
    return np.exp(-0.5 * ((ages - 30.0) / 7.2) ** 2)


def base_housing_profile(ages: np.ndarray) -> np.ndarray:
    """Housing-expenditure coefficients: rising through midlife, then persistent."""
    logistic = 0.27 + 0.43 / (1.0 + np.exp(-(ages - 42.0) / 8.0))
    late_retention = 0.06 / (1.0 + np.exp(-(ages - 63.0) / 4.0))
    return logistic + late_retention


def theta_path(period: float, par: Parameters) -> float:
    return par.theta_long_run + (par.theta_initial - par.theta_long_run) * np.exp(
        -period / par.theta_decay_periods
    )


def phase_in(period: float, par: Parameters) -> float:
    if period <= par.policy_start_period:
        return 0.0
    return 1.0 - np.exp(
        -(period - par.policy_start_period) / par.policy_phasein_periods
    )


def scenario_objects(
    period: float, scenario: Scenario, par: Parameters
) -> tuple[float, float, float]:
    progress = phase_in(period, par)
    access_wedge = par.access_wedge_baseline - (
        par.young_access_reduction * scenario.access_reduction * progress
    )
    old_retention = 1.0 - (
        par.old_retention_reduction
        * scenario.old_retention_reduction
        * progress
    )
    theta = (
        par.theta_initial
        if scenario.freeze_fertility_schedule
        else theta_path(period, par)
    )
    return theta, access_wedge, old_retention


def fertility_by_age(
    price: float,
    theta: float,
    access_wedge: float,
    ages: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> np.ndarray:
    common = theta / (
        par.child_goods_cost + par.child_space * (price + access_wedge)
    )
    return birth_scale * fertility_profile(ages) * common


def housing_by_age(
    price: float,
    theta: float,
    access_wedge: float,
    old_retention: float,
    ages: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> np.ndarray:
    base = base_housing_profile(ages).copy()
    base[ages >= par.old_age_cutoff] *= old_retention
    current_children = fertility_by_age(
        price, theta, access_wedge, ages, birth_scale, par
    ) / birth_scale
    return (
        base / price**par.demand_price_elasticity
        + par.child_space * current_children
    )


def demographic_operator(
    price: float,
    theta: float,
    access_wedge: float,
    ages: np.ndarray,
    survival: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> np.ndarray:
    """Linear transition for adult cohorts followed by the maturation queue."""
    j_count = ages.size
    # An end-of-period birth enters q_0 next date.  Four queued child stages
    # then represent ages 0--4, 5--9, 10--14, and 15--19 before adult entry.
    q_count = par.maturation_periods
    operator = np.zeros((j_count + q_count, j_count + q_count), dtype=float)

    # The oldest child cohort forms new adult households next period.
    operator[0, j_count + q_count - 1] = par.child_queue_survival
    for j in range(j_count - 1):
        operator[j + 1, j] = survival[j]

    births = fertility_by_age(
        price, theta, access_wedge, ages, birth_scale, par
    )
    # Queue mass is measured in future adult-household equivalents.
    operator[j_count, :j_count] = par.household_formation_per_birth * births
    for q in range(q_count - 1):
        operator[j_count + q + 1, j_count + q] = par.child_queue_survival
    return operator


def spectral_radius(operator: np.ndarray) -> float:
    values = np.linalg.eigvals(operator)
    if np.max(np.abs(np.imag(values))) > 1.0e-8:
        # A Leslie-type transition can have complex non-dominant roots.
        pass
    return float(np.max(np.abs(values)))


def replacement_birth_scale(
    ages: np.ndarray, survival: np.ndarray, par: Parameters
) -> float:
    """Normalize the pre-decline environment to demographic replacement at p=1."""
    lower, upper = 1.0e-6, 5.0
    for _ in range(100):
        trial = 0.5 * (lower + upper)
        matrix = demographic_operator(
            1.0,
            par.theta_initial,
            par.access_wedge_baseline,
            ages,
            survival,
            trial,
            par,
        )
        if spectral_radius(matrix) < 1.0:
            lower = trial
        else:
            upper = trial
    scale = 0.5 * (lower + upper)
    radius = spectral_radius(
        demographic_operator(
            1.0,
            par.theta_initial,
            par.access_wedge_baseline,
            ages,
            survival,
            scale,
            par,
        )
    )
    if abs(radius - 1.0) > 1.0e-10:
        raise RuntimeError("Could not normalize pre-decline reproduction to one")
    return scale


def normalized_stationary_state(
    operator: np.ndarray, adult_count: int
) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eig(operator)
    index = int(np.argmin(np.abs(values - 1.0)))
    vector = np.real(vectors[:, index])
    if vector.sum() < 0.0:
        vector *= -1.0
    if np.min(vector) <= 0.0:
        raise RuntimeError("Replacement operator has no positive stationary vector")
    vector /= vector[:adult_count].sum()
    return vector[:adult_count], vector[adult_count:]


def initial_state(
    stationary_adults: np.ndarray,
    stationary_queue: np.ndarray,
    ages: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Create an observed-style prime-age bulge without changing total adults."""
    tilt = (
        0.86
        + 0.34 * np.exp(-0.5 * ((ages - 39.0) / 10.0) ** 2)
        - 0.10 / (1.0 + np.exp(-(ages - 66.0) / 4.0))
    )
    adults = stationary_adults * tilt
    adults /= adults.sum()

    # Recent child cohorts are modestly smaller than the oldest queued cohort,
    # representing a fertility decline already under way before the sample date.
    queue_tilt = np.linspace(0.92, 1.00, stationary_queue.size)
    queue = stationary_queue * queue_tilt
    return adults, queue


def total_housing_demand(
    price: float,
    adults: np.ndarray,
    theta: float,
    access_wedge: float,
    old_retention: float,
    ages: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> float:
    demands = housing_by_age(
        price,
        theta,
        access_wedge,
        old_retention,
        ages,
        birth_scale,
        par,
    )
    return float(adults @ demands)


def clearing_price(
    adults: np.ndarray,
    theta: float,
    access_wedge: float,
    old_retention: float,
    supply_scale: float,
    ages: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> float:
    def residual(price: float) -> float:
        supply = supply_scale * price**par.supply_elasticity
        demand = total_housing_demand(
            price,
            adults,
            theta,
            access_wedge,
            old_retention,
            ages,
            birth_scale,
            par,
        )
        return supply - demand

    lower, upper = 1.0e-5, 1.0
    while residual(upper) < 0.0:
        upper *= 2.0
        if upper > 1.0e5:
            raise RuntimeError("Could not bracket the market-clearing price")
    for _ in range(160):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    price = 0.5 * (lower + upper)
    if abs(residual(price)) > 1.0e-10:
        raise RuntimeError("Housing market did not clear")
    return price


def advance_state(
    adults: np.ndarray,
    queue: np.ndarray,
    price: float,
    theta: float,
    access_wedge: float,
    ages: np.ndarray,
    survival: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    births = float(
        adults
        @ fertility_by_age(
            price, theta, access_wedge, ages, birth_scale, par
        )
    )
    entrants = par.child_queue_survival * queue[-1]
    next_adults = np.zeros_like(adults)
    next_adults[0] = entrants
    next_adults[1:] = survival[:-1] * adults[:-1]

    next_queue = np.zeros_like(queue)
    future_households = par.household_formation_per_birth * births
    next_queue[0] = future_households
    next_queue[1:] = par.child_queue_survival * queue[:-1]

    exits = float(np.sum((1.0 - survival) * adults))
    accounting_gap = float(next_adults.sum() - adults.sum() - (entrants - exits))
    if abs(accounting_gap) > 1.0e-12:
        raise RuntimeError("Adult cohort accounting failed")
    return next_adults, next_queue, {
        "birth_output": births,
        "future_household_output": float(future_households),
        "mature_entrants": float(entrants),
        "adult_exits": exits,
        "adult_accounting_gap": accounting_gap,
    }


def transition_decomposition(
    adults: np.ndarray,
    next_adults: np.ndarray,
    price: float,
    current_objects: tuple[float, float, float],
    next_objects: tuple[float, float, float],
    ages: np.ndarray,
    survival: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> dict[str, float]:
    theta, access, old_retention = current_objects
    next_theta, next_access, next_old_retention = next_objects
    current_d = housing_by_age(
        price, theta, access, old_retention, ages, birth_scale, par
    )
    next_schedule_d = housing_by_age(
        price,
        next_theta,
        next_access,
        next_old_retention,
        ages,
        birth_scale,
        par,
    )

    current_total = float(adults @ current_d)
    cohort_total = float(next_adults @ current_d)
    next_total_fixed_price = float(next_adults @ next_schedule_d)
    cohort_shift = cohort_total - current_total
    schedule_shift = next_total_fixed_price - cohort_total

    entrants = next_adults[0]
    entry_demand = float(entrants * current_d[0])
    aging_gain = float(
        np.sum(survival[:-1] * adults[:-1] * (current_d[1:] - current_d[:-1]))
    )
    exit_loss = float(
        -np.sum((1.0 - survival[:-1]) * adults[:-1] * current_d[:-1])
        - adults[-1] * current_d[-1]
    )
    component_gap = cohort_shift - (entry_demand + aging_gain + exit_loss)
    if abs(component_gap) > 1.0e-11:
        raise RuntimeError("Housing-demand decomposition failed")
    return {
        "fixed_price_demand_shift": next_total_fixed_price - current_total,
        "cohort_demand_shift": cohort_shift,
        "schedule_demand_shift": schedule_shift,
        "entry_demand_contribution": entry_demand,
        "aging_demand_contribution": aging_gain,
        "exit_demand_contribution": exit_loss,
        "demand_decomposition_gap": component_gap,
    }


def simulate(
    scenario: Scenario,
    initial_adults: np.ndarray,
    initial_queue: np.ndarray,
    supply_scale: float,
    ages: np.ndarray,
    survival: np.ndarray,
    birth_scale: float,
    par: Parameters,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    adults = initial_adults.copy()
    queue = initial_queue.copy()
    rows: list[dict[str, float | str]] = []
    cohort_rows: list[dict[str, float | str]] = []

    for period in range(par.periods + 1):
        objects = scenario_objects(period, scenario, par)
        theta, access, old_retention = objects
        price = clearing_price(
            adults,
            theta,
            access,
            old_retention,
            supply_scale,
            ages,
            birth_scale,
            par,
        )
        age_demand = housing_by_age(
            price,
            theta,
            access,
            old_retention,
            ages,
            birth_scale,
            par,
        )
        supply = supply_scale * price**par.supply_elasticity
        demand = float(adults @ age_demand)
        residual = supply - demand
        operator = demographic_operator(
            price, theta, access, ages, survival, birth_scale, par
        )
        growth_factor = spectral_radius(operator)
        young_mask = ages < 40
        prime_mask = (ages >= 40) & (ages < 60)
        old_mask = ages >= 60

        for age, mass, housing in zip(ages, adults, age_demand):
            cohort_rows.append(
                {
                    "scenario": scenario.name,
                    "scenario_label": scenario.label,
                    "period": period,
                    "year": period * par.period_years,
                    "age_group": f"{age}-{age + par.period_years - 1}",
                    "age_lower": int(age),
                    "adult_mass": float(mass),
                    "per_household_housing": float(housing),
                    "cohort_housing_demand": float(mass * housing),
                }
            )

        next_adults = adults.copy()
        next_queue = queue.copy()
        accounting = {
            "birth_output": np.nan,
            "mature_entrants": np.nan,
            "adult_exits": np.nan,
            "adult_accounting_gap": np.nan,
        }
        decomposition = {
            "fixed_price_demand_shift": np.nan,
            "cohort_demand_shift": np.nan,
            "schedule_demand_shift": np.nan,
            "entry_demand_contribution": np.nan,
            "aging_demand_contribution": np.nan,
            "exit_demand_contribution": np.nan,
            "demand_decomposition_gap": np.nan,
        }
        if period < par.periods:
            next_adults, next_queue, accounting = advance_state(
                adults,
                queue,
                price,
                theta,
                access,
                ages,
                survival,
                birth_scale,
                par,
            )
            next_objects = scenario_objects(period + 1, scenario, par)
            decomposition = transition_decomposition(
                adults,
                next_adults,
                price,
                objects,
                next_objects,
                ages,
                survival,
                birth_scale,
                par,
            )

        row: dict[str, float | str] = {
            "scenario": scenario.name,
            "scenario_label": scenario.label,
            "period": period,
            "year": period * par.period_years,
            "theta": theta,
            "access_wedge": access,
            "old_retention_factor": old_retention,
            "price": price,
            "adult_households": float(adults.sum()),
            "young_20_39_share": float(adults[young_mask].sum() / adults.sum()),
            "prime_40_59_share": float(adults[prime_mask].sum() / adults.sum()),
            "old_60_79_share": float(adults[old_mask].sum() / adults.sum()),
            "frozen_demographic_growth_factor": growth_factor,
            "housing_supply": supply,
            "housing_demand": demand,
            "market_residual": residual,
            **accounting,
            **decomposition,
        }
        rows.append(row)
        adults, queue = next_adults, next_queue

    return rows, cohort_rows


def build_condition_grid() -> pd.DataFrame:
    """Exact two-generation momentum conditions at below-replacement R=0.94."""
    reproduction = 0.94
    mu_values = np.linspace(0.65, 1.08, 174)
    kappa_values = np.linspace(0.25, 1.55, 196)
    rows = []
    for mu in mu_values:
        for kappa in kappa_values:
            population_rises = mu < reproduction
            price_rises = (1.0 - mu) * kappa > (1.0 - reproduction)
            rows.append(
                {
                    "reproduction": reproduction,
                    "old_to_young_mass_ratio_mu": float(mu),
                    "old_to_young_housing_ratio_kappa": float(kappa),
                    "population_rises": bool(population_rises),
                    "price_rises": bool(price_rises),
                    "joint_region": bool(population_rises and price_rises),
                    "region_code": int(population_rises) + 2 * int(price_rises),
                    "price_threshold_kappa": (
                        (1.0 - reproduction) / (1.0 - mu)
                        if mu < 1.0
                        else np.inf
                    ),
                }
            )
    return pd.DataFrame(rows)


def configure_plots() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.2,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.2,
            "legend.fontsize": 8.1,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
            "savefig.dpi": 320,
            "savefig.bbox": "tight",
        }
    )


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.pdf", pad_inches=0.15)
    fig.savefig(outdir / f"{stem}.png", pad_inches=0.15)
    plt.close(fig)


def plot_transition_paths(paths: pd.DataFrame, outdir: Path) -> None:
    colors = {
        "baseline": "#1f4e79",
        "no_decline": "#6b6b6b",
        "young_access": "#b54a36",
        "old_release": "#3d7a57",
    }
    styles = {
        "baseline": "-",
        "no_decline": ":",
        "young_access": "--",
        "old_release": "-.",
    }
    initial = float(
        paths.loc[
            (paths["scenario"] == "baseline") & (paths["period"] == 0),
            "adult_households",
        ].iloc[0]
    )
    initial_price = float(
        paths.loc[
            (paths["scenario"] == "baseline") & (paths["period"] == 0),
            "price",
        ].iloc[0]
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True)
    panels = (
        ("price", "Housing-cost index", initial_price),
        ("adult_households", "Adult-household index", initial),
        ("frozen_demographic_growth_factor", "Frozen demographic growth factor", 1.0),
        ("old_60_79_share", "Share of households age 60--79", None),
    )
    for axis, (column, label, normalizer) in zip(axes.flat, panels):
        for scenario in SCENARIOS:
            part = paths[paths["scenario"] == scenario.name]
            values = part[column].to_numpy()
            if normalizer is not None and column != "frozen_demographic_growth_factor":
                values = values / normalizer
            axis.plot(
                part["year"],
                values,
                label=scenario.label,
                color=colors[scenario.name],
                linestyle=styles[scenario.name],
                linewidth=1.8,
            )
        if column == "frozen_demographic_growth_factor":
            axis.axhline(1.0, color="0.35", linewidth=0.8, linestyle=":")
        elif normalizer is not None:
            axis.axhline(1.0, color="0.65", linewidth=0.7, linestyle=":")
        axis.set_ylabel(label)
        axis.grid(axis="y", color="0.90", linewidth=0.5)
    axes[1, 0].set_xlabel("Years from initial distribution")
    axes[1, 1].set_xlabel("Years from initial distribution")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("A fertility decline with inherited demographic momentum", y=0.985)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=4,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.875))
    save_figure(fig, outdir, "transition_paths")


def plot_transition_accounting(
    paths: pd.DataFrame,
    cohorts: pd.DataFrame,
    stationary_adults: np.ndarray,
    ages: np.ndarray,
    outdir: Path,
) -> None:
    baseline = paths[paths["scenario"] == "baseline"].copy()
    initial = cohorts[
        (cohorts["scenario"] == "baseline") & (cohorts["period"] == 0)
    ].sort_values("age_lower")
    stationary = stationary_adults / stationary_adults.sum()

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.25))
    width = 1.8
    axes[0].bar(
        ages - 1.1,
        stationary,
        width=width,
        color="#b8c6d1",
        label="Stationary age profile",
    )
    axes[0].bar(
        ages + 1.1,
        initial["adult_mass"].to_numpy() / initial["adult_mass"].sum(),
        width=width,
        color="#1f4e79",
        label="Inherited profile",
    )
    axes[0].set_xlabel("Age of household head")
    axes[0].set_ylabel("Share of adult households")
    axes[0].set_title("(a) Inherited cohort bulge")
    bar_peak = max(
        float(np.max(stationary)),
        float(np.max(initial["adult_mass"].to_numpy() / initial["adult_mass"].sum())),
    )
    axes[0].set_ylim(0.0, 1.32 * bar_peak)
    axes[0].legend(frameon=False, loc="upper center", ncol=2, fontsize=7.5)
    axes[0].grid(axis="y", color="0.90", linewidth=0.5)

    usable = baseline[baseline["period"] < baseline["period"].max()]
    x = usable["year"].to_numpy()
    components = (
        ("entry_demand_contribution", "Maturing entrants", "#3d7a57"),
        ("aging_demand_contribution", "Aging into larger housing", "#c48a2c"),
        ("exit_demand_contribution", "Deaths / household exit", "#9e4b4b"),
        ("schedule_demand_shift", "Fertility-preference schedule", "#7a68a6"),
    )
    positive_bottom = np.zeros_like(x, dtype=float)
    negative_bottom = np.zeros_like(x, dtype=float)
    for column, label, color in components:
        values = usable[column].to_numpy(dtype=float)
        bottom = np.where(values >= 0.0, positive_bottom, negative_bottom)
        axes[1].bar(
            x,
            values,
            bottom=bottom,
            width=3.6,
            color=color,
            label=label,
        )
        positive_bottom += np.where(values >= 0.0, values, 0.0)
        negative_bottom += np.where(values < 0.0, values, 0.0)
    axes[1].plot(
        x,
        usable["fixed_price_demand_shift"],
        color="black",
        linewidth=1.4,
        marker="o",
        markersize=2.8,
        label="Total at current price",
    )
    axes[1].axhline(0.0, color="0.35", linewidth=0.7)
    axes[1].set_xlabel("Years from initial distribution")
    axes[1].set_ylabel("Change in housing demand")
    axes[1].set_title("(b) Fixed-price demand decomposition")
    axes[1].legend(
        frameon=False,
        fontsize=6.8,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=2,
        columnspacing=0.9,
        handlelength=2.0,
    )
    axes[1].grid(axis="y", color="0.90", linewidth=0.5)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save_figure(fig, outdir, "transition_accounting")


def plot_momentum_conditions(grid: pd.DataFrame, outdir: Path) -> None:
    reproduction = float(grid["reproduction"].iloc[0])
    pivot = grid.pivot(
        index="old_to_young_housing_ratio_kappa",
        columns="old_to_young_mass_ratio_mu",
        values="region_code",
    )
    mus = pivot.columns.to_numpy(dtype=float)
    kappas = pivot.index.to_numpy(dtype=float)
    values = pivot.to_numpy(dtype=float)
    cmap = ListedColormap(["#ececec", "#b7d5e8", "#e9c6ad", "#4f8b6c"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    fig, axis = plt.subplots(figsize=(5.7, 3.6))
    axis.imshow(
        values,
        origin="lower",
        extent=(mus.min(), mus.max(), kappas.min(), kappas.max()),
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
        rasterized=True,
    )
    axis.axvline(reproduction, color="#1f4e79", linestyle="--", linewidth=1.2)
    valid_mu = mus[mus < 0.999]
    threshold = (1.0 - reproduction) / (1.0 - valid_mu)
    mask = threshold <= kappas.max()
    axis.plot(
        valid_mu[mask],
        threshold[mask],
        color="#8b3a2f",
        linewidth=1.5,
    )
    example_mu, example_kappa = 0.82, 0.85
    axis.scatter(
        [example_mu],
        [example_kappa],
        s=34,
        color="black",
        zorder=3,
    )
    axis.annotate(
        "Illustrative inherited bulge",
        (example_mu, example_kappa),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=8,
    )
    axis.text(
        0.665,
        1.49,
        "Population and price rise",
        fontsize=8.4,
        va="top",
        color="#244f3a",
    )
    axis.text(
        0.955,
        0.31,
        "Neither rises",
        fontsize=8.2,
        color="0.30",
    )
    axis.text(
        0.895,
        0.315,
        "Population\nonly",
        fontsize=7.8,
        ha="center",
        color="#315d78",
    )
    axis.text(
        0.947,
        1.40,
        "Price\nonly",
        fontsize=7.8,
        ha="center",
        color="#7b4936",
    )
    axis.set_xlim(mus.min(), mus.max())
    axis.set_ylim(kappas.min(), kappas.max())
    axis.set_xlabel(r"Old-to-young cohort mass, $\mu=O/Y$")
    axis.set_ylabel(r"Old-to-young housing demand, $\kappa=d_O/d_Y$")
    axis.set_title(
        r"Exact momentum regions with reproduction $R=0.94<1$"
    )
    axis.grid(False)
    fig.tight_layout()
    save_figure(fig, outdir, "momentum_conditions")


def selected_results(paths: pd.DataFrame) -> dict[str, object]:
    results: dict[str, object] = {}
    for scenario in SCENARIOS:
        part = paths[paths["scenario"] == scenario.name].copy()
        initial_price = float(part["price"].iloc[0])
        initial_pop = float(part["adult_households"].iloc[0])
        joint = part[
            (part["period"] < part["period"].max())
            & (part["frozen_demographic_growth_factor"] < 1.0 - 1.0e-8)
            & (part["adult_households"].shift(-1) > part["adult_households"])
            & (part["price"].shift(-1) > part["price"])
        ]
        results[scenario.name] = {
            "label": scenario.label,
            "minimum_frozen_growth_factor": float(
                part["frozen_demographic_growth_factor"].min()
            ),
            "peak_price_index": float(part["price"].max() / initial_price),
            "peak_adult_household_index": float(
                part["adult_households"].max() / initial_pop
            ),
            "year_of_peak_price": int(part.loc[part["price"].idxmax(), "year"]),
            "year_of_peak_adult_households": int(
                part.loc[part["adult_households"].idxmax(), "year"]
            ),
            "joint_momentum_interval_start_years": [
                int(v) for v in joint["year"].tolist()
            ],
            "year_25_price_index": float(
                part.loc[part["year"] == 25, "price"].iloc[0] / initial_price
            ),
            "year_25_adult_household_index": float(
                part.loc[part["year"] == 25, "adult_households"].iloc[0]
                / initial_pop
            ),
            "year_50_price_index": float(
                part.loc[part["year"] == 50, "price"].iloc[0] / initial_price
            ),
            "year_50_adult_household_index": float(
                part.loc[part["year"] == 50, "adult_households"].iloc[0]
                / initial_pop
            ),
            "terminal_price_index": float(part["price"].iloc[-1] / initial_price),
            "terminal_adult_household_index": float(
                part["adult_households"].iloc[-1] / initial_pop
            ),
        }
    return results


def write_readme(
    outdir: Path,
    par: Parameters,
    results: dict[str, object],
    checks: dict[str, float | int | bool],
) -> None:
    baseline = results["baseline"]
    text = f"""# Dynamic population--housing illustration

Status: **generic theory illustration; not calibrated and not a forecast.**

The object called `price` is the contemporaneous housing-market clearing cost.
It is not an asset-price transition with forward-looking capital gains. A full
asset-price implementation would additionally impose the intertemporal
no-arbitrage condition linking today's asset price, rent, and tomorrow's asset
price. This packet is therefore a demographic-transition diagnostic, not a
solved perfect-foresight equilibrium of the maintained lifecycle model.

## What is simulated

Adult households occupy twelve five-year age groups from 20--24 through
75--79. Births occur at the end of a period, enter the newborn stage next date,
and form new adult households after four five-year child stages. The child
queue is measured in future adult-household equivalents: each birth contributes
{par.household_formation_per_birth:.2f} equivalents before child survival. A smooth
decline in the value of children
lowers fertility. The initial distribution deliberately contains a prime-age
cohort bulge and smaller recent child cohorts, so it is not a steady-state age
distribution. There is no immigration after the initial date; earlier
immigration is one possible interpretation of the inherited bulge.

At each date the single housing-cost index solves

`H_s(p_t) = sum_j x[j,t] * d_j(p_t)`.

The frozen demographic growth factor is the spectral radius of the complete
adult-aging, fertility, and maturation operator evaluated at current prices
and preferences. A value below one means the current environment is
below replacement even when inherited cohort momentum keeps adult population
growing.

## Main numerical fact

In the baseline, the frozen growth factor falls below one, yet adult population
and housing costs both rise over the transition at years
`{baseline['joint_momentum_interval_start_years']}`. Each listed date starts a
five-year interval in which both objects rise. The peak adult-household index is
`{baseline['peak_adult_household_index']:.4f}` and the peak housing-cost index is
`{baseline['peak_price_index']:.4f}`. These magnitudes are illustrative.

The exact two-generation condition shown in `momentum_conditions.pdf` is:

- population rises with `R < 1` iff `mu = O/Y < R`;
- price rises iff `(1-mu) * kappa > 1-R`, where `kappa = d_O/d_Y`;
- both therefore require a sufficiently young inherited cohort mix and enough
  housing demand as that cohort ages.

## Policy comparisons

- `no_decline`: starts from the same inherited adult and child distributions
  but holds the fertility schedule at its initial value.
- `young_access`: the goods-equivalent access wedge for young families falls
  smoothly by {par.young_access_reduction:.2f}; this raises fertility and child-space demand.
- `old_release`: housing retention at ages 60+ falls smoothly by
  {100 * par.old_retention_reduction:.0f} percent; this releases housing and lowers prices.

The policy experiments are not welfare calculations. All four paths hold the
housing-supply schedule and initial adult and child distributions fixed.

## Files

- `transition_paths.csv`: aggregate paths and every accounting residual.
- `cohort_distributions.csv`: age-specific masses and housing demand.
- `momentum_condition_grid.csv`: exact two-generation sensitivity grid.
- `summary.json`: parameters, formulas, checks, and selected results.
- `transition_paths.pdf`: housing costs, adult population, reproduction, and aging.
- `transition_accounting.pdf`: initial cohort bulge and exact demand decomposition.
- `momentum_conditions.pdf`: conditions for joint population-price momentum.

## Verification

- Maximum absolute housing-market residual: `{checks['max_market_residual']:.3e}`.
- Maximum absolute adult-accounting gap: `{checks['max_adult_accounting_gap']:.3e}`.
- Maximum absolute demand-decomposition gap: `{checks['max_demand_decomposition_gap']:.3e}`.
- Fixed-price demand sign predicts the next price sign in all
  `{checks['price_sign_comparisons']}` nondegenerate transitions.

Rebuild from the repository root with:

```bash
code/model/.venv/bin/python code/model/tools/build_dynamic_population_housing_paper_illustrations.py
```
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def verify_paths(paths: pd.DataFrame) -> dict[str, float | int | bool]:
    max_market = float(paths["market_residual"].abs().max())
    max_accounting = float(paths["adult_accounting_gap"].abs().max(skipna=True))
    max_decomp = float(paths["demand_decomposition_gap"].abs().max(skipna=True))
    if max_market > 1.0e-9 or max_accounting > 1.0e-11 or max_decomp > 1.0e-10:
        raise RuntimeError("An accounting or market-clearing assertion failed")

    comparisons = 0
    for scenario in SCENARIOS:
        part = paths[paths["scenario"] == scenario.name].sort_values("period")
        current = part.iloc[:-1]
        next_price_change = part["price"].shift(-1).iloc[:-1].to_numpy() - current[
            "price"
        ].to_numpy()
        demand_shift = current["fixed_price_demand_shift"].to_numpy()
        nondegenerate = (np.abs(next_price_change) > 1.0e-9) | (
            np.abs(demand_shift) > 1.0e-9
        )
        if not np.all(np.sign(next_price_change[nondegenerate]) == np.sign(demand_shift[nondegenerate])):
            raise RuntimeError("Fixed-price demand failed to predict the next price sign")
        comparisons += int(np.sum(nondegenerate))

    baseline = paths[paths["scenario"] == "baseline"].copy()
    joint = (
        (baseline["period"] < baseline["period"].max())
        & (baseline["frozen_demographic_growth_factor"] < 1.0 - 1.0e-8)
        & (baseline["adult_households"].shift(-1) > baseline["adult_households"])
        & (baseline["price"].shift(-1) > baseline["price"])
    )
    if not bool(joint.any()):
        raise RuntimeError(
            "The baseline does not illustrate joint momentum below replacement"
        )
    return {
        "max_market_residual": max_market,
        "max_adult_accounting_gap": max_accounting,
        "max_demand_decomposition_gap": max_decomp,
        "price_sign_comparisons": comparisons,
        "baseline_has_joint_momentum_below_replacement": True,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUT,
        help="Output directory (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    par = Parameters()
    ages = age_grid(par)
    survival = adult_survival(ages)
    birth_scale = replacement_birth_scale(ages, survival, par)
    pre_operator = demographic_operator(
        1.0,
        par.theta_initial,
        par.access_wedge_baseline,
        ages,
        survival,
        birth_scale,
        par,
    )
    stationary_adults, stationary_queue = normalized_stationary_state(
        pre_operator, ages.size
    )
    adults0, queue0 = initial_state(
        stationary_adults, stationary_queue, ages
    )

    # Normalize the initial market-clearing price to one without modifying the
    # inherited state or any behavioral parameter.
    initial_objects = scenario_objects(0.0, SCENARIOS[0], par)
    supply_scale = total_housing_demand(
        1.0,
        adults0,
        *initial_objects,
        ages,
        birth_scale,
        par,
    )

    all_rows: list[dict[str, float | str]] = []
    all_cohorts: list[dict[str, float | str]] = []
    for scenario in SCENARIOS:
        rows, cohorts = simulate(
            scenario,
            adults0,
            queue0,
            supply_scale,
            ages,
            survival,
            birth_scale,
            par,
        )
        all_rows.extend(rows)
        all_cohorts.extend(cohorts)

    paths = pd.DataFrame(all_rows)
    cohorts = pd.DataFrame(all_cohorts)
    grid = build_condition_grid()
    paths.to_csv(outdir / "transition_paths.csv", index=False, float_format="%.12g")
    cohorts.to_csv(
        outdir / "cohort_distributions.csv", index=False, float_format="%.12g"
    )
    grid.to_csv(
        outdir / "momentum_condition_grid.csv", index=False, float_format="%.12g"
    )

    checks = verify_paths(paths)
    results = selected_results(paths)
    summary = {
        "status": "generic_theory_transition_diagnostic_not_calibrated_not_full_asset_price_equilibrium",
        "model": {
            "adult_transition": "x_0,t+1 = child_survival * queue_last,t; x_j+1,t+1 = survival_j * x_j,t",
            "maturation_queue": "q_0,t+1 = household_formation_per_birth * sum_j fertility_j(p_t,theta_t,access_t) * x_j,t; q_k+1,t+1 = child_survival * q_k,t",
            "birth_to_entry_delay_periods": par.maturation_periods,
            "housing_market": "supply_scale * p_t^supply_elasticity = sum_j x_j,t * d_j(p_t,theta_t,access_t,retention_t)",
            "price_interpretation": "contemporaneous housing-market clearing cost; not a forward-looking asset price",
            "reproduction": "spectral radius of the frozen adult-aging/fertility/maturation operator",
            "two_cohort_population_condition": "mu < R < 1",
            "two_cohort_price_condition": "(1-mu)*kappa > 1-R, kappa=d_O/d_Y",
        },
        "parameters": {**asdict(par), "birth_scale": birth_scale, "supply_scale": supply_scale},
        "age_grid": [int(v) for v in ages],
        "adult_survival": [float(v) for v in survival],
        "checks": checks,
        "results": results,
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_readme(outdir, par, results, checks)

    configure_plots()
    plot_transition_paths(paths, outdir)
    plot_transition_accounting(
        paths, cohorts, stationary_adults, ages, outdir
    )
    plot_momentum_conditions(grid, outdir)

    baseline = results["baseline"]
    print(f"Wrote dynamic transition packet to {outdir}")
    print(
        "Baseline joint momentum years: "
        f"{baseline['joint_momentum_interval_start_years']}; "
        f"peak price index={baseline['peak_price_index']:.4f}; "
        f"peak adult-household index={baseline['peak_adult_household_index']:.4f}."
    )
    print(
        "Checks: "
        f"market={checks['max_market_residual']:.2e}, "
        f"adult accounting={checks['max_adult_accounting_gap']:.2e}, "
        f"demand decomposition={checks['max_demand_decomposition_gap']:.2e}."
    )


if __name__ == "__main__":
    main()
