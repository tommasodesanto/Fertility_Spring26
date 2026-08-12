#!/usr/bin/env python3
"""Build the numerical illustration for the demographic-momentum theory note.

The script uses only the Python standard library.  It writes the transition
paths consumed directly by the note's PGFPlots figure and a small parameter
table.  The example is illustrative, not calibrated.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output" / "model" / "demographic_momentum_intergenerational_housing_note"


@dataclass(frozen=True)
class Parameters:
    alpha: float = 0.30
    beta_old: float = 0.80
    child_space: float = 0.50
    child_goods: float = 0.50
    access_wedge: float = 0.20
    supply_scale: float = 1.00
    supply_elasticity: float = 0.00
    replacement: float = 1.00
    theta_old: float = 1.10
    theta_new: float = 0.95
    young_initial: float = 0.72
    old_initial: float = 0.53


def fertility(price: float, theta: float, par: Parameters) -> float:
    return theta / (
        par.child_goods + par.child_space * (price + par.access_wedge)
    )


def young_housing(price: float, theta: float, par: Parameters) -> float:
    return par.alpha / price + par.child_space * fertility(price, theta, par)


def old_housing(price: float, par: Parameters) -> float:
    return par.beta_old / price


def excess_supply(
    price: float, young: float, old: float, theta: float, par: Parameters
) -> float:
    supply = par.supply_scale * price**par.supply_elasticity
    demand = young * young_housing(price, theta, par) + old * old_housing(price, par)
    return supply - demand


def clearing_price(young: float, old: float, theta: float, par: Parameters) -> float:
    lower = 1.0e-8
    upper = 1.0
    while excess_supply(upper, young, old, theta, par) < 0.0:
        upper *= 2.0
        if upper > 1.0e8:
            raise RuntimeError("Could not bracket the market-clearing price")
    for _ in range(200):
        midpoint = 0.5 * (lower + upper)
        if excess_supply(midpoint, young, old, theta, par) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def steady_state(theta: float, par: Parameters) -> tuple[float, float, float]:
    price = (
        theta / par.replacement - par.child_goods
    ) / par.child_space - par.access_wedge
    if price <= 0.0:
        raise ValueError("The illustrative parameters do not support a positive steady state")
    per_young = young_housing(price, theta, par)
    per_old = old_housing(price, par)
    young = (
        par.supply_scale * price**par.supply_elasticity / (per_young + per_old)
    )
    return price, young, 2.0 * young


def simulate(theta: float, periods: int, par: Parameters) -> list[dict[str, float]]:
    young = par.young_initial
    old = par.old_initial
    rows: list[dict[str, float]] = []
    for period in range(periods + 1):
        price = clearing_price(young, old, theta, par)
        births = fertility(price, theta, par)
        reproduction = births / par.replacement
        rows.append(
            {
                "period": period,
                "price": price,
                "population": young + old,
                "young": young,
                "old": old,
                "reproduction": reproduction,
            }
        )
        young, old = reproduction * young, young
    return rows


def main() -> None:
    par = Parameters()
    OUT.mkdir(parents=True, exist_ok=True)

    old_price, old_young, old_population = steady_state(par.theta_old, par)
    new_price, new_young, new_population = steady_state(par.theta_new, par)
    shock = simulate(par.theta_new, periods=20, par=par)
    no_shock = simulate(par.theta_old, periods=20, par=par)

    initial_old_price = clearing_price(
        par.young_initial, par.old_initial, par.theta_old, par
    )
    assert abs(par.young_initial + par.old_initial - old_population) < 1.0e-12
    assert abs(initial_old_price - old_price) < 1.0e-12
    assert shock[0]["reproduction"] < 1.0
    assert shock[1]["reproduction"] < 1.0
    assert shock[1]["population"] > shock[0]["population"]
    assert shock[1]["price"] > shock[0]["price"]
    assert shock[1]["price"] > old_price
    assert all(
        shocked["price"] < counterfactual["price"]
        for shocked, counterfactual in zip(shock[:11], no_shock[:11])
    )
    assert abs(shock[-1]["price"] - new_price) < 2.0e-3
    assert abs(shock[-1]["population"] - new_population) < 2.0e-3

    path_file = OUT / "transition_paths.csv"
    with path_file.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "period",
            "price_shock",
            "population_shock",
            "reproduction_shock",
            "price_no_shock",
            "population_no_shock",
            "reproduction_no_shock",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for shocked, counterfactual in zip(shock, no_shock):
            writer.writerow(
                {
                    "period": shocked["period"],
                    "price_shock": shocked["price"] / old_price,
                    "population_shock": shocked["population"] / old_population,
                    "reproduction_shock": shocked["reproduction"],
                    "price_no_shock": counterfactual["price"] / old_price,
                    "population_no_shock": counterfactual["population"] / old_population,
                    "reproduction_no_shock": counterfactual["reproduction"],
                }
            )

    parameter_rows = asdict(par)
    parameter_rows.update(
        {
            "old_ss_price": old_price,
            "old_ss_young": old_young,
            "old_ss_population": old_population,
            "new_ss_price": new_price,
            "new_ss_young": new_young,
            "new_ss_population": new_population,
            "shock_impact_price": shock[0]["price"],
            "shock_impact_reproduction": shock[0]["reproduction"],
            "shock_period1_price": shock[1]["price"],
            "shock_period1_population": shock[1]["population"],
        }
    )
    with (OUT / "illustration_parameters.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["object", "value"])
        writer.writerows(parameter_rows.items())

    print(f"Wrote {path_file}")
    print(
        "Illustration: "
        f"impact (p,R)=({shock[0]['price']:.4f},{shock[0]['reproduction']:.4f}); "
        f"next period (p,N)=({shock[1]['price']:.4f},{shock[1]['population']:.4f}); "
        f"new steady state (p,N)=({new_price:.4f},{new_population:.4f})."
    )


if __name__ == "__main__":
    main()
