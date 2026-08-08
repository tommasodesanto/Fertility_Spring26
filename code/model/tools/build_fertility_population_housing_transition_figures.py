"""Build the linked steady-state diagrams for the fertility-housing note.

The figure is schematic.  Its two panels visualize the exact analytical maps
derived in ``latex/fertility_population_housing_transition_note.tex``:

1. household child demand maps replacement fertility into a house price; and
2. housing demand and supply map that house price into population.

The parameter values only determine the geometry of the drawing.  They are not
a calibration or a quantitative counterfactual.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = REPO_ROOT / "output" / "model" / "fertility_population_housing_transition_note"
LEGACY_TIME_PATH_OUTPUTS = (
    "population_transition.pdf",
    "housing_price_transition.pdf",
    "transition_paths.csv",
    "illustration_parameters.csv",
)


def preference_weight_for_price(
    price: float,
    *,
    alpha: float,
    income: float,
    child_cost: float,
    child_space: float,
    replacement_fertility: float,
) -> float:
    """Recover chi so that household fertility equals replacement at ``price``."""

    expenditure_share = (
        (child_cost + child_space * price) * replacement_fertility / income
    )
    if not 0.0 < expenditure_share < 1.0:
        raise ValueError("The illustrative replacement-fertility share must be in (0, 1).")
    return expenditure_share * (1.0 + alpha) / (1.0 - expenditure_share)


def inverse_fertility_demand(
    fertility: np.ndarray,
    *,
    chi: float,
    alpha: float,
    income: float,
    child_cost: float,
    child_space: float,
) -> np.ndarray:
    """Return the price consistent with a chosen fertility level."""

    utility_weight_sum = 1.0 + alpha + chi
    child_expenditure = chi * income / utility_weight_sum
    return (child_expenditure / fertility - child_cost) / child_space


def replacement_housing_demand(
    price: np.ndarray | float,
    *,
    alpha: float,
    income: float,
    child_cost: float,
    child_space: float,
    replacement_fertility: float,
) -> np.ndarray | float:
    """Per-household housing demand conditional on replacement fertility."""

    discretionary_component = (
        alpha
        * (income - child_cost * replacement_fertility)
        / ((1.0 + alpha) * price)
    )
    required_component = child_space * replacement_fertility / (1.0 + alpha)
    return discretionary_component + required_component


def population_supported_by_price(
    price: np.ndarray | float,
    *,
    supply_scale: float,
    supply_elasticity: float,
    alpha: float,
    income: float,
    child_cost: float,
    child_space: float,
    replacement_fertility: float,
) -> np.ndarray | float:
    """Return population from stationary housing-market clearing."""

    housing_supply = supply_scale * np.power(price, supply_elasticity)
    housing_per_household = replacement_housing_demand(
        price,
        alpha=alpha,
        income=income,
        child_cost=child_cost,
        child_space=child_space,
        replacement_fertility=replacement_fertility,
    )
    return housing_supply / housing_per_household


def illustration() -> dict[str, float]:
    """Construct a transparent normalization for the schematic diagram."""

    parameters = {
        "alpha": 1.0,
        "income": 4.0,
        "child_cost": 0.4,
        "child_space": 0.6,
        "replacement_fertility": 1.0,
        "old_price": 1.0,
        "new_price": 0.70,
        "supply_elasticity": 0.40,
    }
    parameters["old_chi"] = preference_weight_for_price(
        parameters["old_price"],
        alpha=parameters["alpha"],
        income=parameters["income"],
        child_cost=parameters["child_cost"],
        child_space=parameters["child_space"],
        replacement_fertility=parameters["replacement_fertility"],
    )
    parameters["new_chi"] = preference_weight_for_price(
        parameters["new_price"],
        alpha=parameters["alpha"],
        income=parameters["income"],
        child_cost=parameters["child_cost"],
        child_space=parameters["child_space"],
        replacement_fertility=parameters["replacement_fertility"],
    )
    old_housing_per_household = replacement_housing_demand(
        parameters["old_price"],
        alpha=parameters["alpha"],
        income=parameters["income"],
        child_cost=parameters["child_cost"],
        child_space=parameters["child_space"],
        replacement_fertility=parameters["replacement_fertility"],
    )
    parameters["supply_scale"] = old_housing_per_household / np.power(
        parameters["old_price"], parameters["supply_elasticity"]
    )
    parameters["old_population"] = population_supported_by_price(
        parameters["old_price"],
        supply_scale=parameters["supply_scale"],
        supply_elasticity=parameters["supply_elasticity"],
        alpha=parameters["alpha"],
        income=parameters["income"],
        child_cost=parameters["child_cost"],
        child_space=parameters["child_space"],
        replacement_fertility=parameters["replacement_fertility"],
    )
    parameters["new_population"] = population_supported_by_price(
        parameters["new_price"],
        supply_scale=parameters["supply_scale"],
        supply_elasticity=parameters["supply_elasticity"],
        alpha=parameters["alpha"],
        income=parameters["income"],
        child_cost=parameters["child_cost"],
        child_space=parameters["child_space"],
        replacement_fertility=parameters["replacement_fertility"],
    )
    return parameters


def write_parameters(parameters: dict[str, float]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with (OUTDIR / "diagram_parameters.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["parameter", "value"])
        for name, value in parameters.items():
            writer.writerow([name, f"{value:.12g}"])


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def build_figure(parameters: dict[str, float]) -> None:
    alpha = parameters["alpha"]
    income = parameters["income"]
    child_cost = parameters["child_cost"]
    child_space = parameters["child_space"]
    replacement = parameters["replacement_fertility"]
    old_price = parameters["old_price"]
    new_price = parameters["new_price"]
    old_population = parameters["old_population"]
    new_population = parameters["new_population"]

    old_color = "#245985"
    new_color = "#b4493f"
    neutral = "#555555"
    guide = "#a4a4a4"

    fig, (fertility_axis, population_axis) = plt.subplots(
        1,
        2,
        figsize=(8.65, 3.05),
        sharey=True,
        gridspec_kw={"wspace": 0.30},
    )

    fertility_grid = np.linspace(0.64, 1.48, 350)
    old_schedule = inverse_fertility_demand(
        fertility_grid,
        chi=parameters["old_chi"],
        alpha=alpha,
        income=income,
        child_cost=child_cost,
        child_space=child_space,
    )
    new_schedule = inverse_fertility_demand(
        fertility_grid,
        chi=parameters["new_chi"],
        alpha=alpha,
        income=income,
        child_cost=child_cost,
        child_space=child_space,
    )
    fertility_axis.plot(fertility_grid, old_schedule, color=old_color, linewidth=2.1)
    fertility_axis.plot(
        fertility_grid,
        new_schedule,
        color=new_color,
        linewidth=2.1,
        linestyle=(0, (5, 3)),
    )
    fertility_axis.axvline(
        replacement, color=neutral, linewidth=1.1, linestyle=(0, (2, 2))
    )
    fertility_axis.scatter(
        [replacement, replacement],
        [old_price, new_price],
        color=[old_color, new_color],
        s=31,
        zorder=5,
    )
    fertility_axis.annotate(
        "",
        xy=(replacement, new_price + 0.025),
        xytext=(replacement, old_price - 0.025),
        arrowprops={"arrowstyle": "->", "color": neutral, "lw": 1.1},
    )
    fertility_axis.text(
        0.69,
        1.49,
        r"old schedule $p_h^F(\ell;\chi_0)$",
        color=old_color,
        ha="left",
        va="center",
        fontsize=9.0,
    )
    fertility_axis.text(
        0.69,
        1.10,
        r"new schedule $p_h^F(\ell;\chi_1)$",
        color=new_color,
        ha="left",
        va="center",
        fontsize=9.0,
    )
    fertility_axis.text(
        replacement + 0.018,
        0.30,
        "replacement",
        color=neutral,
        rotation=90,
        ha="left",
        va="bottom",
        fontsize=8.7,
    )
    fertility_axis.set_xlim(0.62, 1.50)
    fertility_axis.set_xlabel(r"fertility $\ell$")
    fertility_axis.set_ylabel(r"house price $p_h$")
    fertility_axis.set_title("A. Household fertility schedule", loc="left", fontweight="semibold")
    fertility_axis.set_xticks([replacement], [r"$\bar\ell$"])

    price_grid = np.linspace(0.26, 1.56, 400)
    population_grid = population_supported_by_price(
        price_grid,
        supply_scale=parameters["supply_scale"],
        supply_elasticity=parameters["supply_elasticity"],
        alpha=alpha,
        income=income,
        child_cost=child_cost,
        child_space=child_space,
        replacement_fertility=replacement,
    )
    population_axis.plot(
        population_grid, price_grid, color=neutral, linewidth=2.2
    )
    population_axis.scatter(
        [old_population, new_population],
        [old_price, new_price],
        color=[old_color, new_color],
        s=31,
        zorder=5,
    )
    for population, price, color in (
        (old_population, old_price, old_color),
        (new_population, new_price, new_color),
    ):
        population_axis.hlines(
            price, 0.0, population, color=color, linewidth=1.0, linestyle=(0, (3, 3))
        )
        population_axis.vlines(
            population, 0.24, price, color=color, linewidth=1.0, linestyle=(0, (3, 3))
        )
    population_axis.annotate(
        "",
        xy=(new_population + 0.012, new_price + 0.012),
        xytext=(old_population - 0.012, old_price - 0.012),
        arrowprops={"arrowstyle": "->", "color": neutral, "lw": 1.1},
    )
    population_axis.text(
        old_population + 0.025,
        old_price + 0.035,
        r"$E_0$",
        color=old_color,
        ha="left",
        va="bottom",
    )
    population_axis.text(
        new_population - 0.025,
        new_price - 0.04,
        r"$E_1$",
        color=new_color,
        ha="right",
        va="top",
    )
    label_price = 1.30
    label_population = float(
        population_supported_by_price(
            label_price,
            supply_scale=parameters["supply_scale"],
            supply_elasticity=parameters["supply_elasticity"],
            alpha=alpha,
            income=income,
            child_cost=child_cost,
            child_space=child_space,
            replacement_fertility=replacement,
        )
    )
    population_axis.text(
        label_population + 0.035,
        label_price,
        r"housing-market locus $N^H(p_h)$",
        color=neutral,
        ha="left",
        va="center",
        fontsize=9.0,
    )
    population_axis.set_xlim(0.0, 1.65)
    population_axis.set_xlabel(r"adult households $N$")
    population_axis.set_title(
        "B. Housing clearing at replacement", loc="left", fontweight="semibold"
    )
    population_axis.set_xticks(
        [new_population, old_population], [r"$N_1^*$", r"$N_0^*$"]
    )

    for axis in (fertility_axis, population_axis):
        axis.set_ylim(0.24, 1.62)
        axis.set_yticks([new_price, old_price], [r"$p_1^*$", r"$p_0^*$"])
        axis.grid(False)

    fig.savefig(OUTDIR / "steady_state_comparative_statics.pdf", format="pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for filename in LEGACY_TIME_PATH_OUTPUTS:
        (OUTDIR / filename).unlink(missing_ok=True)
    parameters = illustration()
    write_parameters(parameters)
    apply_plot_style()
    build_figure(parameters)


if __name__ == "__main__":
    main()
