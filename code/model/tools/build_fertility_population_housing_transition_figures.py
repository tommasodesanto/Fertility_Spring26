"""Build the two linked diagrams for the simple fertility--housing note.

The parameters are chosen only to make the analytical relationships easy to
see.  They are not a calibration or a quantitative counterfactual.
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


def inverse_fertility(
    fertility: np.ndarray,
    beta: float,
    child_cost: float,
    child_space: float,
) -> np.ndarray:
    """House price at which a family chooses a given fertility level."""

    return (beta / fertility - child_cost) / child_space


def population_at_replacement(
    price: np.ndarray | float,
    *,
    supply_scale: float,
    supply_elasticity: float,
    housing_weight: float,
    child_space: float,
    replacement: float,
) -> np.ndarray | float:
    """Population supported by the housing market at replacement fertility."""

    supply = supply_scale * np.power(price, supply_elasticity)
    housing_per_family = housing_weight / price + child_space * replacement
    return supply / housing_per_family


def illustration() -> dict[str, float]:
    """Choose a transparent normalization for the schematic figure."""

    values = {
        "housing_weight": 1.0,
        "child_cost": 0.4,
        "child_space": 0.6,
        "replacement": 1.0,
        "old_price": 1.0,
        "new_price": 0.68,
        "supply_elasticity": 0.4,
    }
    values["old_beta"] = values["replacement"] * (
        values["child_cost"] + values["child_space"] * values["old_price"]
    )
    values["new_beta"] = values["replacement"] * (
        values["child_cost"] + values["child_space"] * values["new_price"]
    )

    old_housing = (
        values["housing_weight"] / values["old_price"]
        + values["child_space"] * values["replacement"]
    )
    values["supply_scale"] = old_housing
    values["old_population"] = population_at_replacement(
        values["old_price"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        child_space=values["child_space"],
        replacement=values["replacement"],
    )
    values["new_population"] = population_at_replacement(
        values["new_price"],
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        child_space=values["child_space"],
        replacement=values["replacement"],
    )
    return values


def verify(values: dict[str, float]) -> None:
    """Check that the drawing obeys the model's comparative statics."""

    assert values["new_beta"] < values["old_beta"]
    assert values["new_price"] < values["old_price"]
    assert values["new_population"] < values["old_population"]

    price_grid = np.linspace(0.2, 1.35, 300)
    population_grid = population_at_replacement(
        price_grid,
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        child_space=values["child_space"],
        replacement=values["replacement"],
    )
    assert np.all(np.diff(population_grid) > 0.0)


def write_parameters(values: dict[str, float]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with (OUTDIR / "diagram_parameters.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["parameter", "value"])
        for name, value in values.items():
            writer.writerow([name, f"{value:.12g}"])


def build_figure(values: dict[str, float]) -> None:
    """Draw the fertility--price and price--population maps."""

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )

    old_color = "#244f74"
    new_color = "#a3473f"
    guide_color = "#999999"

    fig, (left, right) = plt.subplots(
        1,
        2,
        figsize=(8.8, 2.75),
        sharey=True,
        gridspec_kw={"wspace": 0.26},
    )

    fertility = np.linspace(0.68, 1.42, 350)
    old_curve = inverse_fertility(
        fertility,
        values["old_beta"],
        values["child_cost"],
        values["child_space"],
    )
    new_curve = inverse_fertility(
        fertility,
        values["new_beta"],
        values["child_cost"],
        values["child_space"],
    )
    left.plot(fertility, old_curve, color=old_color, linewidth=2.0)
    left.plot(
        fertility,
        new_curve,
        color=new_color,
        linewidth=2.0,
        linestyle=(0, (5, 3)),
    )
    left.axvline(
        values["replacement"],
        color="#555555",
        linewidth=1.0,
        linestyle=(0, (2, 2)),
    )
    left.scatter(
        [values["replacement"], values["replacement"]],
        [values["old_price"], values["new_price"]],
        color=[old_color, new_color],
        s=28,
        zorder=5,
    )
    old_label_height = inverse_fertility(
        np.array([1.27]),
        values["old_beta"],
        values["child_cost"],
        values["child_space"],
    )[0]
    new_label_height = inverse_fertility(
        np.array([1.27]),
        values["new_beta"],
        values["child_cost"],
        values["child_space"],
    )[0]
    left.text(
        1.27,
        old_label_height + 0.04,
        r"before: $\beta_0$",
        color=old_color,
        ha="center",
    )
    left.text(
        1.27,
        new_label_height - 0.09,
        r"after: $\beta_1$",
        color=new_color,
        ha="center",
    )
    left.text(
        values["replacement"] + 0.012,
        1.28,
        "replacement",
        rotation=90,
        va="top",
        color="#555555",
    )
    left.set_title("(a) Fertility and the house price")
    left.set_xlabel(r"Children per family, $n$")
    left.set_ylabel(r"House price, $p$")
    left.set_xlim(0.68, 1.42)
    left.set_xticks([values["replacement"]], labels=[r"$\bar n$"])
    left.set_yticks(
        [values["new_price"], values["old_price"]],
        labels=[r"$p_1$", r"$p_0$"],
    )

    price = np.linspace(0.20, 1.35, 350)
    population = population_at_replacement(
        price,
        supply_scale=values["supply_scale"],
        supply_elasticity=values["supply_elasticity"],
        housing_weight=values["housing_weight"],
        child_space=values["child_space"],
        replacement=values["replacement"],
    )
    right.plot(population, price, color="#333333", linewidth=2.1)
    right.scatter(
        [values["old_population"], values["new_population"]],
        [values["old_price"], values["new_price"]],
        color=[old_color, new_color],
        s=28,
        zorder=5,
    )
    right.text(0.98, 1.04, "housing market", rotation=33, color="#333333")
    right.set_title("(b) House price and population")
    right.set_xlabel(r"Population (families), $N$")
    right.set_xlim(0.25, 1.35)
    right.set_xticks(
        [values["new_population"], values["old_population"]],
        labels=[r"$N_1$", r"$N_0$"],
    )

    for axis in (left, right):
        axis.set_ylim(0.20, 1.35)
        axis.grid(False)

    for price_value in (values["old_price"], values["new_price"]):
        left.hlines(
            price_value,
            0.68,
            values["replacement"],
            color=guide_color,
            linewidth=0.8,
            linestyle=(0, (3, 3)),
        )
        right.axhline(
            price_value,
            color=guide_color,
            linewidth=0.8,
            linestyle=(0, (3, 3)),
        )

    for population_value, price_value in (
        (values["old_population"], values["old_price"]),
        (values["new_population"], values["new_price"]),
    ):
        right.vlines(
            population_value,
            0.20,
            price_value,
            color=guide_color,
            linewidth=0.8,
            linestyle=(0, (3, 3)),
        )

    left.annotate(
        r"$\beta$ falls",
        xy=(
            0.84,
            inverse_fertility(
                np.array([0.84]),
                values["new_beta"],
                values["child_cost"],
                values["child_space"],
            )[0],
        ),
        xytext=(0.72, 1.23),
        arrowprops={"arrowstyle": "->", "color": "#555555", "linewidth": 1.0},
        ha="center",
        color="#555555",
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "steady_state_comparative_statics.pdf", pad_inches=0.04)
    fig.savefig(OUTDIR / "steady_state_comparative_statics.png", dpi=220, pad_inches=0.04)
    plt.close(fig)


def main() -> None:
    values = illustration()
    verify(values)
    write_parameters(values)
    build_figure(values)


if __name__ == "__main__":
    main()
