"""Build the linked housing-market and child-decision diagrams.

The parameters are chosen only to make the analytical transition easy to see.
They are not a calibration or a quantitative counterfactual.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = (
    REPO_ROOT / "output" / "model" / "fertility_population_housing_transition_note"
)


def child_demand(
    price: np.ndarray | float,
    child_value: float,
    child_cost: float,
    child_space: float,
) -> np.ndarray | float:
    """Return children per family at a given current housing cost."""

    return child_value / (child_cost + child_space * price)


def housing_per_family(
    price: np.ndarray | float,
    child_value: float,
    *,
    housing_weight: float,
    child_cost: float,
    child_space: float,
) -> np.ndarray | float:
    """Return adult plus child-related housing demand per family."""

    adult_housing = housing_weight / price
    child_housing = child_space * child_demand(
        price,
        child_value,
        child_cost,
        child_space,
    )
    return adult_housing + child_housing


def population_at_price(
    price: np.ndarray | float,
    child_value: float,
    *,
    supply_scale: float,
    supply_elasticity: float,
    housing_weight: float,
    old_housing_weight: float,
    child_cost: float,
    child_space: float,
) -> np.ndarray | float:
    """Return total young-plus-old population supported at a steady-state price."""

    housing_supply = supply_scale * np.power(price, supply_elasticity)
    young_housing = housing_per_family(
        price,
        child_value,
        housing_weight=housing_weight,
        child_cost=child_cost,
        child_space=child_space,
    )
    old_housing = old_housing_weight / price
    average_housing_per_adult = 0.5 * (young_housing + old_housing)
    return housing_supply / average_housing_per_adult


def price_at_population(
    population: float,
    child_value: float,
    *,
    supply_scale: float,
    supply_elasticity: float,
    housing_weight: float,
    old_housing_weight: float,
    child_cost: float,
    child_space: float,
) -> float:
    """Invert the monotone housing-clearing schedule."""

    price_grid = np.linspace(0.01, 4.0, 30_000)
    population_grid = population_at_price(
        price_grid,
        child_value,
        supply_scale=supply_scale,
        supply_elasticity=supply_elasticity,
        housing_weight=housing_weight,
        old_housing_weight=old_housing_weight,
        child_cost=child_cost,
        child_space=child_space,
    )
    if population < population_grid[0] or population > population_grid[-1]:
        raise ValueError("Population lies outside the inversion grid.")
    return float(np.interp(population, population_grid, price_grid))


def illustration() -> dict[str, float]:
    """Choose a normalization with visible impact and transition movements."""

    values = {
        "income": 3.0,
        "housing_weight": 0.6,
        "old_housing_weight": 0.4,
        "child_cost": 0.2,
        "child_space": 0.8,
        "replacement": 1.0,
        "old_theta": 1.0,
        "new_theta": 0.7,
        "supply_elasticity": 0.4,
        "old_population": 1.0,
    }
    values["old_price"] = (
        values["old_theta"] / values["replacement"] - values["child_cost"]
    ) / values["child_space"]
    values["new_price"] = (
        values["new_theta"] / values["replacement"] - values["child_cost"]
    ) / values["child_space"]

    young_housing = housing_per_family(
        values["old_price"],
        values["old_theta"],
        housing_weight=values["housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    old_housing = values["old_housing_weight"] / values["old_price"]
    average_housing_per_adult = 0.5 * (young_housing + old_housing)
    values["supply_scale"] = (
        values["old_population"]
        * average_housing_per_adult
        / values["old_price"] ** values["supply_elasticity"]
    )
    values["new_population"] = population_at_price(
        values["new_price"],
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    values["impact_price"] = price_at_population(
        values["old_population"],
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    values["impact_fertility"] = child_demand(
        values["impact_price"],
        values["new_theta"],
        values["child_cost"],
        values["child_space"],
    )
    return values


def verify(values: dict[str, float]) -> None:
    """Check every ordering displayed in the figure."""

    assert values["new_theta"] < values["old_theta"]
    assert values["new_theta"] > values["child_cost"] * values["replacement"]
    assert values["income"] > values["housing_weight"] + values["old_theta"]
    assert np.isclose(
        child_demand(
            values["old_price"],
            values["old_theta"],
            values["child_cost"],
            values["child_space"],
        ),
        values["replacement"],
    )
    assert np.isclose(
        child_demand(
            values["new_price"],
            values["new_theta"],
            values["child_cost"],
            values["child_space"],
        ),
        values["replacement"],
    )
    assert np.isclose(
        population_at_price(
            values["old_price"],
            values["old_theta"],
            supply_scale=values["supply_scale"],
            supply_elasticity=values["supply_elasticity"],
            housing_weight=values["housing_weight"],
            old_housing_weight=values["old_housing_weight"],
            child_cost=values["child_cost"],
            child_space=values["child_space"],
        ),
        values["old_population"],
    )
    assert values["new_population"] < values["old_population"]
    assert values["new_price"] < values["impact_price"] < values["old_price"]
    assert values["impact_fertility"] < values["replacement"]

    price_grid = np.linspace(0.2, 1.4, 500)
    old_population = population_at_price(
        price_grid,
        values["old_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    new_population = population_at_price(
        price_grid,
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    assert np.all(np.diff(old_population) > 0.0)
    assert np.all(np.diff(new_population) > 0.0)
    assert np.all(new_population > old_population)


def write_parameters(values: dict[str, float]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with (OUTDIR / "diagram_parameters.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["parameter", "value"])
        for name, value in values.items():
            writer.writerow([name, f"{value:.12g}"])


def add_arrow(
    axis: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = "#555555",
) -> None:
    """Add an unobtrusive directional arrow between two points."""

    axis.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={"arrowstyle": "->", "color": color, "linewidth": 1.25},
    )


def build_figure(values: dict[str, float]) -> None:
    """Draw population-to-price and price-to-children maps."""

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.2,
            "axes.labelsize": 10.5,
            "axes.titlesize": 11.0,
            "xtick.labelsize": 9.3,
            "ytick.labelsize": 9.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )

    old_color = "#244f74"
    new_color = "#a3473f"
    guide_color = "#a0a0a0"
    transition_color = "#555555"

    fig, (market_axis, child_axis) = plt.subplots(
        1,
        2,
        figsize=(9.2, 3.45),
        gridspec_kw={"wspace": 0.30},
    )

    price_grid = np.linspace(0.42, 1.16, 500)
    old_population_curve = population_at_price(
        price_grid,
        values["old_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    new_population_curve = population_at_price(
        price_grid,
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    market_axis.plot(
        old_population_curve,
        price_grid,
        color=old_color,
        linewidth=2.1,
    )
    market_axis.plot(
        new_population_curve,
        price_grid,
        color=new_color,
        linewidth=2.1,
        linestyle=(0, (5, 3)),
    )
    market_axis.text(1.08, 1.085, "before", color=old_color)
    market_axis.text(1.08, 0.80, "after", color=new_color)
    market_axis.set_title("(a) Population and the housing cost")
    market_axis.set_xlabel(r"Adult-household population, $S$")
    market_axis.set_ylabel(r"Housing cost, $p$")

    old_point = (values["old_population"], values["old_price"])
    impact_market_point = (values["old_population"], values["impact_price"])
    new_point = (values["new_population"], values["new_price"])
    market_axis.scatter(
        [old_point[0], impact_market_point[0], new_point[0]],
        [old_point[1], impact_market_point[1], new_point[1]],
        color=[old_color, new_color, new_color],
        s=[36, 34, 36],
        zorder=6,
    )
    market_axis.text(old_point[0] + 0.025, old_point[1] + 0.025, r"$A$", color=old_color)
    market_axis.text(
        impact_market_point[0] + 0.025,
        impact_market_point[1] - 0.055,
        r"$I$",
        color=new_color,
    )
    market_axis.text(new_point[0] - 0.065, new_point[1] + 0.025, r"$A'$", color=new_color)
    add_arrow(market_axis, old_point, impact_market_point)

    transition_start_price = values["impact_price"] - 0.055
    transition_end_price = values["new_price"] + 0.055
    transition_start_population = population_at_price(
        transition_start_price,
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    transition_end_population = population_at_price(
        transition_end_price,
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    add_arrow(
        market_axis,
        (transition_start_population, transition_start_price),
        (transition_end_population, transition_end_price),
    )
    market_axis.text(
        0.70,
        0.70,
        r"$S$ falls",
        rotation=34,
        color=transition_color,
        fontsize=9.2,
    )

    child_price_grid = np.linspace(0.46, 1.16, 500)
    old_child_curve = child_demand(
        child_price_grid,
        values["old_theta"],
        values["child_cost"],
        values["child_space"],
    )
    new_child_curve = child_demand(
        child_price_grid,
        values["new_theta"],
        values["child_cost"],
        values["child_space"],
    )
    child_axis.plot(child_price_grid, old_child_curve, color=old_color, linewidth=2.1)
    child_axis.plot(
        child_price_grid,
        new_child_curve,
        color=new_color,
        linewidth=2.1,
        linestyle=(0, (5, 3)),
    )
    child_axis.axhline(
        values["replacement"],
        color="#666666",
        linewidth=1.0,
        linestyle=(0, (2, 2)),
    )
    child_axis.text(1.045, 0.92, "before", color=old_color)
    child_axis.text(0.995, 0.62, "after", color=new_color)
    child_axis.set_title("(b) Housing cost and fertility")
    child_axis.set_xlabel(r"Housing cost, $p$")
    child_axis.set_ylabel(r"Reproduction rate, $n$")

    impact_child_point = (values["impact_price"], values["impact_fertility"])
    old_child_point = (values["old_price"], values["replacement"])
    new_child_point = (values["new_price"], values["replacement"])
    child_axis.scatter(
        [old_child_point[0], impact_child_point[0], new_child_point[0]],
        [old_child_point[1], impact_child_point[1], new_child_point[1]],
        color=[old_color, new_color, new_color],
        s=[36, 34, 36],
        zorder=6,
    )
    child_axis.text(old_child_point[0] + 0.014, old_child_point[1] + 0.035, r"$A$", color=old_color)
    child_axis.text(
        impact_child_point[0] + 0.015,
        impact_child_point[1] - 0.065,
        r"$I$",
        color=new_color,
    )
    child_axis.text(new_child_point[0] - 0.045, new_child_point[1] + 0.035, r"$A'$", color=new_color)
    add_arrow(child_axis, old_child_point, impact_child_point)

    child_transition_start_price = values["impact_price"] - 0.045
    child_transition_end_price = values["new_price"] + 0.045
    add_arrow(
        child_axis,
        (
            child_transition_start_price,
            child_demand(
                child_transition_start_price,
                values["new_theta"],
                values["child_cost"],
                values["child_space"],
            ),
        ),
        (
            child_transition_end_price,
            child_demand(
                child_transition_end_price,
                values["new_theta"],
                values["child_cost"],
                values["child_space"],
            ),
        ),
    )
    child_axis.text(
        0.72,
        0.88,
        r"$n$ rises",
        rotation=-31,
        color=transition_color,
        fontsize=9.2,
    )

    market_axis.set_xlim(0.35, 1.36)
    market_axis.set_ylim(0.42, 1.16)
    market_axis.set_xticks(
        [values["new_population"], values["old_population"]],
        labels=[r"$S_1$", r"$S_0$"],
    )
    market_axis.set_yticks(
        [values["new_price"], values["impact_price"], values["old_price"]],
        labels=[r"$p_1$", r"$\widetilde p$", r"$p_0$"],
    )

    child_axis.set_xlim(0.46, 1.16)
    child_axis.set_ylim(0.55, 1.42)
    child_axis.set_xticks(
        [values["new_price"], values["impact_price"], values["old_price"]],
        labels=[r"$p_1$", r"$\widetilde p$", r"$p_0$"],
    )
    child_axis.set_yticks(
        [values["impact_fertility"], values["replacement"]],
        labels=[r"$\widetilde n$", r"$1$"],
    )

    for axis in (market_axis, child_axis):
        axis.grid(False)
        axis.tick_params(direction="out", length=3.5)

    for price_value, population_value in (
        (values["old_price"], values["old_population"]),
        (values["impact_price"], values["old_population"]),
        (values["new_price"], values["new_population"]),
    ):
        market_axis.hlines(
            price_value,
            0.35,
            population_value,
            color=guide_color,
            linewidth=0.7,
            linestyle=(0, (3, 3)),
            zorder=0,
        )

    for price_value, fertility_value in (
        (values["old_price"], values["replacement"]),
        (values["impact_price"], values["impact_fertility"]),
        (values["new_price"], values["replacement"]),
    ):
        child_axis.vlines(
            price_value,
            0.55,
            fertility_value,
            color=guide_color,
            linewidth=0.7,
            linestyle=(0, (3, 3)),
            zorder=0,
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "steady_state_comparative_statics.pdf", pad_inches=0.05)
    fig.savefig(
        OUTDIR / "steady_state_comparative_statics.png",
        dpi=220,
        pad_inches=0.05,
    )
    plt.close(fig)


def build_equilibrium_figure(values: dict[str, float]) -> None:
    """Draw the two schedules that determine the initial steady state."""

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.2,
            "axes.labelsize": 10.5,
            "axes.titlesize": 11.0,
            "xtick.labelsize": 9.3,
            "ytick.labelsize": 9.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )

    old_color = "#244f74"
    guide_color = "#a0a0a0"
    fig, (market_axis, child_axis) = plt.subplots(
        1,
        2,
        figsize=(9.2, 3.45),
        gridspec_kw={"wspace": 0.30},
    )

    price_grid = np.linspace(0.42, 1.16, 500)
    population_curve = population_at_price(
        price_grid,
        values["old_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    market_axis.plot(population_curve, price_grid, color=old_color, linewidth=2.1)
    market_axis.set_title("(a) Housing-market clearing")
    market_axis.set_xlabel(r"Adult-household population, $S$")
    market_axis.set_ylabel(r"Housing cost, $p$")

    child_price_grid = np.linspace(0.46, 1.16, 500)
    child_curve = child_demand(
        child_price_grid,
        values["old_theta"],
        values["child_cost"],
        values["child_space"],
    )
    child_axis.plot(child_price_grid, child_curve, color=old_color, linewidth=2.1)
    child_axis.axhline(
        values["replacement"],
        color="#666666",
        linewidth=1.0,
        linestyle=(0, (2, 2)),
    )
    child_axis.set_title("(b) Replacement fertility")
    child_axis.set_xlabel(r"Housing cost, $p$")
    child_axis.set_ylabel(r"Reproduction rate, $n$")

    steady_market = (values["old_population"], values["old_price"])
    steady_child = (values["old_price"], values["replacement"])
    market_axis.scatter(*steady_market, color=old_color, s=38, zorder=6)
    child_axis.scatter(*steady_child, color=old_color, s=38, zorder=6)
    market_axis.text(
        steady_market[0] + 0.025,
        steady_market[1] + 0.025,
        r"$A$",
        color=old_color,
    )
    child_axis.text(
        steady_child[0] + 0.014,
        steady_child[1] + 0.035,
        r"$A$",
        color=old_color,
    )

    market_axis.set_xlim(0.35, 1.36)
    market_axis.set_ylim(0.42, 1.16)
    market_axis.set_xticks([values["old_population"]], labels=[r"$S_0$"])
    market_axis.set_yticks([values["old_price"]], labels=[r"$p_0$"])
    child_axis.set_xlim(0.46, 1.16)
    child_axis.set_ylim(0.55, 1.42)
    child_axis.set_xticks([values["old_price"]], labels=[r"$p_0$"])
    child_axis.set_yticks([values["replacement"]], labels=[r"$1$"])

    market_axis.hlines(
        values["old_price"],
        0.35,
        values["old_population"],
        color=guide_color,
        linewidth=0.7,
        linestyle=(0, (3, 3)),
        zorder=0,
    )
    child_axis.vlines(
        values["old_price"],
        0.55,
        values["replacement"],
        color=guide_color,
        linewidth=0.7,
        linestyle=(0, (3, 3)),
        zorder=0,
    )
    for axis in (market_axis, child_axis):
        axis.grid(False)
        axis.tick_params(direction="out", length=3.5)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "steady_state_equilibrium.pdf", pad_inches=0.05)
    fig.savefig(
        OUTDIR / "steady_state_equilibrium.png",
        dpi=220,
        pad_inches=0.05,
    )
    plt.close(fig)


def build_transition_stage_figure(
    values: dict[str, float],
    *,
    stage: str,
    output_stem: str,
) -> None:
    """Draw either the impact response or the later demographic adjustment."""

    if stage not in {"partial_equilibrium", "demographic_adjustment"}:
        raise ValueError(f"Unknown transition stage: {stage}")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.2,
            "axes.labelsize": 10.5,
            "axes.titlesize": 11.0,
            "xtick.labelsize": 9.3,
            "ytick.labelsize": 9.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )

    old_color = "#244f74"
    new_color = "#a3473f"
    guide_color = "#a0a0a0"
    transition_color = "#555555"
    fig, (market_axis, child_axis) = plt.subplots(
        1,
        2,
        figsize=(9.2, 3.45),
        gridspec_kw={"wspace": 0.30},
    )

    price_grid = np.linspace(0.42, 1.16, 500)
    old_population_curve = population_at_price(
        price_grid,
        values["old_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    new_population_curve = population_at_price(
        price_grid,
        values["new_theta"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        old_housing_weight=values["old_housing_weight"],
        child_cost=values["child_cost"],
        child_space=values["child_space"],
    )
    child_price_grid = np.linspace(0.46, 1.16, 500)
    old_child_curve = child_demand(
        child_price_grid,
        values["old_theta"],
        values["child_cost"],
        values["child_space"],
    )
    new_child_curve = child_demand(
        child_price_grid,
        values["new_theta"],
        values["child_cost"],
        values["child_space"],
    )

    old_market = (values["old_population"], values["old_price"])
    impact_market = (values["old_population"], values["impact_price"])
    new_market = (values["new_population"], values["new_price"])
    old_child = (values["old_price"], values["replacement"])
    impact_child = (values["impact_price"], values["impact_fertility"])
    new_child = (values["new_price"], values["replacement"])

    if stage == "partial_equilibrium":
        market_axis.plot(
            old_population_curve,
            price_grid,
            color=old_color,
            linewidth=2.1,
            label="initial",
        )
        market_axis.plot(
            new_population_curve,
            price_grid,
            color=new_color,
            linewidth=2.1,
            linestyle=(0, (5, 3)),
            label=r"after lower $v$",
        )
        child_axis.plot(
            child_price_grid,
            old_child_curve,
            color=old_color,
            linewidth=2.1,
            label="initial",
        )
        child_axis.plot(
            child_price_grid,
            new_child_curve,
            color=new_color,
            linewidth=2.1,
            linestyle=(0, (5, 3)),
            label=r"after lower $v$",
        )
        market_axis.scatter(
            [old_market[0], impact_market[0]],
            [old_market[1], impact_market[1]],
            color=[old_color, new_color],
            s=[38, 36],
            zorder=6,
        )
        child_axis.scatter(
            [old_child[0], impact_child[0]],
            [old_child[1], impact_child[1]],
            color=[old_color, new_color],
            s=[38, 36],
            zorder=6,
        )
        add_arrow(market_axis, old_market, impact_market)
        add_arrow(child_axis, old_child, impact_child)
        market_axis.text(old_market[0] + 0.025, old_market[1] + 0.025, r"$A$", color=old_color)
        market_axis.text(
            impact_market[0] + 0.025,
            impact_market[1] - 0.055,
            r"$I$",
            color=new_color,
        )
        child_axis.text(old_child[0] + 0.014, old_child[1] + 0.035, r"$A$", color=old_color)
        child_axis.text(
            impact_child[0] + 0.015,
            impact_child[1] - 0.065,
            r"$I$",
            color=new_color,
        )
        market_axis.set_xticks([values["old_population"]], labels=[r"$S_0$"])
        market_axis.set_yticks(
            [values["impact_price"], values["old_price"]],
            labels=[r"$\widetilde p$", r"$p_0$"],
        )
        child_axis.set_xticks(
            [values["impact_price"], values["old_price"]],
            labels=[r"$\widetilde p$", r"$p_0$"],
        )
        child_axis.set_yticks(
            [values["impact_fertility"], values["replacement"]],
            labels=[r"$\widetilde n$", r"$1$"],
        )
        market_axis.legend(frameon=False, fontsize=8, loc="upper left")
        child_axis.legend(frameon=False, fontsize=8, loc="upper right")
    else:
        market_axis.plot(
            new_population_curve,
            price_grid,
            color=new_color,
            linewidth=2.1,
        )
        child_axis.plot(
            child_price_grid,
            new_child_curve,
            color=new_color,
            linewidth=2.1,
        )
        market_axis.scatter(
            [impact_market[0], new_market[0]],
            [impact_market[1], new_market[1]],
            color=new_color,
            s=[36, 38],
            zorder=6,
        )
        child_axis.scatter(
            [impact_child[0], new_child[0]],
            [impact_child[1], new_child[1]],
            color=new_color,
            s=[36, 38],
            zorder=6,
        )
        add_arrow(market_axis, impact_market, new_market)
        add_arrow(child_axis, impact_child, new_child)
        market_axis.text(
            impact_market[0] + 0.025,
            impact_market[1] - 0.055,
            r"$I$",
            color=new_color,
        )
        market_axis.text(
            new_market[0] - 0.065,
            new_market[1] + 0.025,
            r"$A'$",
            color=new_color,
        )
        child_axis.text(
            impact_child[0] + 0.015,
            impact_child[1] - 0.065,
            r"$I$",
            color=new_color,
        )
        child_axis.text(
            new_child[0] - 0.045,
            new_child[1] + 0.035,
            r"$A'$",
            color=new_color,
        )
        market_axis.text(
            0.70,
            0.70,
            r"$S$ falls",
            rotation=34,
            color=transition_color,
            fontsize=9.2,
        )
        child_axis.text(
            0.72,
            0.88,
            r"$n$ rises",
            rotation=-31,
            color=transition_color,
            fontsize=9.2,
        )
        market_axis.set_xticks(
            [values["new_population"], values["old_population"]],
            labels=[r"$S_1$", r"$S_0$"],
        )
        market_axis.set_yticks(
            [values["new_price"], values["impact_price"]],
            labels=[r"$p_1$", r"$\widetilde p$"],
        )
        child_axis.set_xticks(
            [values["new_price"], values["impact_price"]],
            labels=[r"$p_1$", r"$\widetilde p$"],
        )
        child_axis.set_yticks(
            [values["impact_fertility"], values["replacement"]],
            labels=[r"$\widetilde n$", r"$1$"],
        )

    child_axis.axhline(
        values["replacement"],
        color="#666666",
        linewidth=1.0,
        linestyle=(0, (2, 2)),
    )
    market_axis.set_title("(a) Population and the housing cost")
    market_axis.set_xlabel(r"Adult-household population, $S$")
    market_axis.set_ylabel(r"Housing cost, $p$")
    child_axis.set_title("(b) Housing cost and fertility")
    child_axis.set_xlabel(r"Housing cost, $p$")
    child_axis.set_ylabel(r"Reproduction rate, $n$")
    market_axis.set_xlim(0.35, 1.36)
    market_axis.set_ylim(0.42, 1.16)
    child_axis.set_xlim(0.46, 1.16)
    child_axis.set_ylim(0.55, 1.42)

    market_guides = (
        (old_market, impact_market)
        if stage == "partial_equilibrium"
        else (impact_market, new_market)
    )
    child_guides = (
        (old_child, impact_child)
        if stage == "partial_equilibrium"
        else (impact_child, new_child)
    )
    for population_value, price_value in market_guides:
        market_axis.hlines(
            price_value,
            0.35,
            population_value,
            color=guide_color,
            linewidth=0.7,
            linestyle=(0, (3, 3)),
            zorder=0,
        )
    for price_value, fertility_value in child_guides:
        child_axis.vlines(
            price_value,
            0.55,
            fertility_value,
            color=guide_color,
            linewidth=0.7,
            linestyle=(0, (3, 3)),
            zorder=0,
        )
    for axis in (market_axis, child_axis):
        axis.grid(False)
        axis.tick_params(direction="out", length=3.5)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / f"{output_stem}.pdf", pad_inches=0.05)
    fig.savefig(OUTDIR / f"{output_stem}.png", dpi=220, pad_inches=0.05)
    plt.close(fig)


def main() -> None:
    values = illustration()
    verify(values)
    write_parameters(values)
    build_equilibrium_figure(values)
    build_transition_stage_figure(
        values,
        stage="partial_equilibrium",
        output_stem="partial_equilibrium_response",
    )
    build_transition_stage_figure(
        values,
        stage="demographic_adjustment",
        output_stem="demographic_adjustment",
    )
    build_figure(values)


if __name__ == "__main__":
    main()
