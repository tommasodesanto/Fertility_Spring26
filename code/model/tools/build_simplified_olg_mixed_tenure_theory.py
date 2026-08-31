#!/usr/bin/env python3
"""Solve and plot the mixed-tenure version of the simplified OLG model.

The script is the constructive numerical companion to the paper-facing theory
section.  It keeps both renting and ownership active by adding an iid
type-I-extreme-value taste to the young tenure choice.  Conditional household
policies are the exact analytical policies in the theory appendix: young
renter and owner housing limits bind, young saving is interior, old renter
housing is capped, and old incumbent retention and estate constraints are
slack.  Every active-set inequality is verified after solving.  The dynamic
object is a terminally closed finite-horizon approximation: terminal prices,
transfers, and continuation values are fixed at the final steady state, while
the induced terminal state is recorded rather than imposed.

Outputs
-------
latex/figures/simplified_olg_reallocation_wedge.png
latex/figures/simplified_olg_mixed_tenure_transition.png
output/model/simplified_olg_mixed_tenure_theory/steady_states.csv
output/model/simplified_olg_mixed_tenure_theory/transition_path.csv
output/model/simplified_olg_mixed_tenure_theory/verification.json
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from math import log
from pathlib import Path
from typing import Callable

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
FIGURE_DIR = ROOT / "latex" / "figures"
OUTPUT_DIR = ROOT / "output" / "model" / "simplified_olg_mixed_tenure_theory"
WEDGE_FIGURE = FIGURE_DIR / "simplified_olg_reallocation_wedge.png"
TRANSITION_FIGURE = FIGURE_DIR / "simplified_olg_mixed_tenure_transition.png"


@dataclass(frozen=True)
class Parameters:
    q: float = 0.50
    tau_p: float = 0.05
    tau_g: float = 0.20
    phi: float = 0.80
    alpha: float = 0.40
    kappa: float = 0.50
    chi: float = 0.15
    beta: float = 0.40
    gamma: float = 0.30
    omega_b: float = 0.40
    income: float = 1.00
    liquid_wealth: float = 0.06
    nu: float = 2.00
    housing_stock: float = 2.00
    renter_cap: float = 0.45
    owner_cap: float = 2.00
    mean_ownership_taste: float = -0.15
    tenure_taste_scale: float = 0.12
    theta_initial: float = 0.35
    theta_final: float = 0.50

    @property
    def gross_return(self) -> float:
        return 1.0 / self.q

    @property
    def old_share_sum(self) -> float:
        return 1.0 + self.gamma + self.omega_b

    @property
    def k(self) -> float:
        return self.beta * self.old_share_sum


P = Parameters()


def logistic(value: float) -> float:
    """Numerically stable logistic tenure probability."""
    return float(1.0 / (1.0 + np.exp(-np.clip(value, -50.0, 50.0))))


def gain_tax_per_unit(
    price_change: float,
    appreciates: bool | None = None,
) -> float:
    """Evaluate the statutory max rule or a fixed smooth appreciation branch."""
    if appreciates is None:
        return P.tau_g * max(price_change, 0.0)
    return P.tau_g * price_change if appreciates else 0.0


def solve_fertility(theta: float, housing: float, wealth: float, user_cost: float) -> float:
    """Solve the unique binding-housing-limit fertility first-order condition."""
    lower = 1.0e-11
    upper = min(
        housing / P.kappa,
        (wealth - user_cost * housing) / P.chi,
    ) * (1.0 - 1.0e-10)
    if upper <= lower:
        raise ValueError("The binding-limit allocation has an empty log domain.")

    def residual(fertility: float) -> float:
        family_space = housing - P.kappa * fertility
        disposable = wealth - user_cost * housing - P.chi * fertility
        return (
            theta / fertility
            - P.alpha * P.kappa / family_space
            - (1.0 + P.k) * P.chi / disposable
        )

    if residual(lower) <= 0.0 or residual(upper) >= 0.0:
        raise ValueError("The fertility first-order condition is not bracketed.")
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def numerical_jacobian(
    function: Callable[[np.ndarray], np.ndarray],
    point: np.ndarray,
    step: float = 1.0e-5,
) -> np.ndarray:
    baseline = function(point)
    jacobian = np.empty((baseline.size, point.size))
    for column in range(point.size):
        shifted = point.copy()
        shifted[column] += step
        jacobian[:, column] = (function(shifted) - baseline) / step
    return jacobian


def damped_newton(
    function: Callable[[np.ndarray], np.ndarray],
    initial: np.ndarray,
    tolerance: float = 1.0e-10,
    max_iterations: int = 70,
) -> tuple[np.ndarray, int]:
    point = initial.copy()
    for iteration in range(max_iterations):
        residual = function(point)
        norm = float(np.max(np.abs(residual)))
        if norm < tolerance:
            return point, iteration
        jacobian = numerical_jacobian(function, point)
        try:
            direction = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError:
            direction = np.linalg.lstsq(jacobian, -residual, rcond=None)[0]
        damping = 1.0
        accepted = False
        for _ in range(35):
            candidate = point + damping * direction
            try:
                candidate_norm = float(np.max(np.abs(function(candidate))))
            except (FloatingPointError, OverflowError, ValueError):
                candidate_norm = np.inf
            if candidate_norm < norm:
                point = candidate
                accepted = True
                break
            damping *= 0.5
        if not accepted:
            raise RuntimeError(f"Newton line search failed at residual {norm:.3e}.")
    raise RuntimeError("Newton iteration limit reached.")


def old_renter_allocation(
    financial_wealth: float,
    price: float,
    next_price: float,
    transfer: float,
) -> dict[str, float | bool]:
    rent = (1.0 + P.tau_p) * price - P.q * next_price
    resources = financial_wealth + transfer
    if rent <= 0.0 or resources <= 0.0:
        raise ValueError("Old renter resources or rent are nonpositive.")
    unconstrained_housing = P.gamma * resources / (P.old_share_sum * rent)
    cap_binds = unconstrained_housing > P.renter_cap
    if cap_binds:
        housing = P.renter_cap
        residual_resources = resources - rent * housing
        consumption = residual_resources / (1.0 + P.omega_b)
        estate = P.omega_b * residual_resources / ((1.0 + P.omega_b) * P.q)
    else:
        housing = unconstrained_housing
        consumption = resources / P.old_share_sum
        estate = P.omega_b * resources / (P.old_share_sum * P.q)
    if min(consumption, housing, estate) <= 0.0:
        raise ValueError("Old renter allocation violates the log domain.")
    utility = (
        log(consumption)
        + P.gamma * log(housing)
        + P.omega_b * log(estate)
    )
    shadow_wedge = (
        P.gamma * consumption / housing - rent if cap_binds else 0.0
    )
    return {
        "consumption": consumption,
        "housing": housing,
        "estate": estate,
        "utility": utility,
        "rent": rent,
        "unconstrained_housing": unconstrained_housing,
        "cap_binds": cap_binds,
        "cap_shadow_wedge": shadow_wedge,
    }


def old_owner_allocation(
    financial_wealth: float,
    inherited_housing: float,
    basis_price: float,
    price: float,
    next_price: float,
    transfer: float,
    appreciates: bool | None = None,
) -> dict[str, float]:
    gain_tax = gain_tax_per_unit(price - basis_price, appreciates)
    rent = (1.0 + P.tau_p) * price - P.q * next_price
    incumbent_cost = rent - gain_tax
    resources = (
        financial_wealth
        + (price - gain_tax) * inherited_housing
        + transfer
    )
    if incumbent_cost <= 0.0 or resources <= 0.0:
        raise ValueError("Old incumbent resources or user cost are nonpositive.")
    consumption = resources / P.old_share_sum
    housing = P.gamma * resources / (P.old_share_sum * incumbent_cost)
    estate = P.omega_b * resources / (P.old_share_sum * P.q)
    retention_slack = inherited_housing - housing
    estate_slack = estate - next_price * housing
    if retention_slack <= 0.0 or estate_slack <= 0.0:
        raise ValueError("The old-incumbent interior active set fails.")
    utility = (
        log(consumption)
        + P.gamma * log(housing)
        + P.omega_b * log(estate)
    )
    return {
        "consumption": consumption,
        "housing": housing,
        "estate": estate,
        "utility": utility,
        "gain_tax": gain_tax,
        "rent": rent,
        "incumbent_cost": incumbent_cost,
        "sales": retention_slack,
        "retention_slack": retention_slack,
        "estate_slack": estate_slack,
    }


def unconstrained_young_housing(
    theta: float,
    wealth: float,
    user_cost: float,
) -> float:
    share_sum = 1.0 + P.k + P.alpha + theta
    family_space = P.alpha * wealth / (share_sum * user_cost)
    fertility = theta * wealth / (
        share_sum * (P.chi + P.kappa * user_cost)
    )
    return family_space + P.kappa * fertility


def conditional_young_policy(
    tenure: str,
    theta: float,
    price: float,
    next_price: float,
    price_after_next: float,
    transfer: float,
    next_transfer: float,
    next_appreciates: bool | None = None,
) -> dict[str, float | bool | dict[str, float | bool]]:
    rent = (1.0 + P.tau_p) * price - P.q * next_price
    next_gain_tax = gain_tax_per_unit(
        next_price - price, next_appreciates
    )
    owner = tenure == "owner"
    user_cost = rent + (P.q * next_gain_tax if owner else 0.0)
    lifetime_wealth = (
        P.income + P.liquid_wealth + transfer + P.q * next_transfer
    )
    if owner:
        housing_cap = min(
            P.owner_cap,
            P.liquid_wealth / ((1.0 - P.phi) * price),
        )
    else:
        housing_cap = P.renter_cap
    unconstrained_housing = unconstrained_young_housing(
        theta, lifetime_wealth, user_cost
    )
    cap_binds = unconstrained_housing > housing_cap
    if cap_binds:
        housing = housing_cap
        fertility = solve_fertility(
            theta, housing, lifetime_wealth, user_cost
        )
        consumption = (
            lifetime_wealth
            - user_cost * housing
            - P.chi * fertility
        ) / (1.0 + P.k)
    else:
        share_sum = 1.0 + P.k + P.alpha + theta
        consumption = lifetime_wealth / share_sum
        fertility = theta * lifetime_wealth / (
            share_sum * (P.chi + P.kappa * user_cost)
        )
        housing = unconstrained_housing

    if owner:
        saving = (
            P.k * consumption
            - P.q * next_transfer
            - (
                P.q * (next_price - next_gain_tax) - P.phi * price
            )
            * housing
        )
        next_financial_wealth = (
            P.gross_return * saving
            - P.gross_return * P.phi * price * housing
        )
        future_old = old_owner_allocation(
            next_financial_wealth,
            housing,
            price,
            next_price,
            price_after_next,
            next_transfer,
            appreciates=next_appreciates,
        )
    else:
        saving = P.k * consumption - P.q * next_transfer
        next_financial_wealth = P.gross_return * saving
        future_old = old_renter_allocation(
            next_financial_wealth,
            next_price,
            price_after_next,
            next_transfer,
        )
    if saving <= 0.0:
        raise ValueError("The young saving interior condition fails.")

    family_space = housing - P.kappa * fertility
    if min(consumption, fertility, family_space) <= 0.0:
        raise ValueError("The young allocation violates the log domain.")
    young_utility = (
        log(consumption)
        + P.alpha * log(family_space)
        + theta * log(fertility)
    )
    conditional_value = young_utility + P.beta * float(future_old["utility"])
    marginal_value_housing = P.alpha * consumption / family_space
    access_wedge = marginal_value_housing - user_cost
    return {
        "consumption": consumption,
        "fertility": fertility,
        "housing": housing,
        "saving": saving,
        "next_financial_wealth": next_financial_wealth,
        "conditional_value": conditional_value,
        "user_cost": user_cost,
        "marginal_value_housing": marginal_value_housing,
        "access_wedge": access_wedge,
        "housing_cap": housing_cap,
        "unconstrained_housing": unconstrained_housing,
        "cap_binds": cap_binds,
        "future_old": future_old,
    }


def young_mixture(
    theta: float,
    price: float,
    next_price: float,
    price_after_next: float,
    transfer: float,
    next_transfer: float,
    next_appreciates: bool | None = None,
) -> dict[str, float | dict[str, float | bool | dict[str, float | bool]]]:
    renter = conditional_young_policy(
        "renter",
        theta,
        price,
        next_price,
        price_after_next,
        transfer,
        next_transfer,
        next_appreciates,
    )
    owner = conditional_young_policy(
        "owner",
        theta,
        price,
        next_price,
        price_after_next,
        transfer,
        next_transfer,
        next_appreciates,
    )
    value_gap = (
        float(owner["conditional_value"])
        + P.mean_ownership_taste
        - float(renter["conditional_value"])
    )
    owner_share = logistic(value_gap / P.tenure_taste_scale)
    average_fertility = (
        (1.0 - owner_share) * float(renter["fertility"])
        + owner_share * float(owner["fertility"])
    )
    average_housing = (
        (1.0 - owner_share) * float(renter["housing"])
        + owner_share * float(owner["housing"])
    )
    return {
        "renter": renter,
        "owner": owner,
        "owner_share": owner_share,
        "average_fertility": average_fertility,
        "average_housing": average_housing,
        "value_gap": value_gap,
    }


def steady_state_values(theta: float, price: float, transfer: float) -> dict:
    young = young_mixture(theta, price, price, price, transfer, transfer)
    renter = young["renter"]
    owner = young["owner"]
    owner_share = float(young["owner_share"])
    old_renter = old_renter_allocation(
        float(renter["next_financial_wealth"]), price, price, transfer
    )
    old_owner = old_owner_allocation(
        float(owner["next_financial_wealth"]),
        float(owner["housing"]),
        price,
        price,
        price,
        transfer,
    )
    average_old_housing = (
        (1.0 - owner_share) * float(old_renter["housing"])
        + owner_share * float(old_owner["housing"])
    )
    housing_per_cohort_pair = (
        float(young["average_housing"]) + average_old_housing
    )
    cohort_mass = P.housing_stock / housing_per_cohort_pair
    return {
        "theta": theta,
        "price": price,
        "transfer": transfer,
        "young": young,
        "old_renter": old_renter,
        "old_owner": old_owner,
        "average_old_housing": average_old_housing,
        "housing_per_cohort_pair": housing_per_cohort_pair,
        "cohort_mass": cohort_mass,
    }


def solve_steady_state(theta: float, initial_price: float) -> dict:
    def residual(log_price_transfer: np.ndarray) -> np.ndarray:
        price, transfer = np.exp(log_price_transfer)
        values = steady_state_values(theta, price, transfer)
        young = values["young"]
        cohort = values["cohort_mass"]
        government_revenue = P.tau_p * price * P.housing_stock
        return np.array(
            [
                P.nu * float(young["average_fertility"]) - 1.0,
                transfer - government_revenue / (2.0 * cohort),
            ]
        )

    solution, iterations = damped_newton(
        residual,
        np.log(np.array([initial_price, 0.011])),
    )
    price, transfer = np.exp(solution)
    values = steady_state_values(theta, price, transfer)
    values["iterations"] = iterations
    values["max_residual"] = float(np.max(np.abs(residual(solution))))
    values["steady_state_jacobian_determinant"] = float(
        np.linalg.det(numerical_jacobian(residual, solution))
    )
    return values


def transition(
    initial_ss: dict,
    final_ss: dict,
    horizon: int,
) -> dict:
    dates = np.arange(horizon + 1)
    price_guess = initial_ss["price"] + (
        final_ss["price"] - initial_ss["price"]
    ) * (1.0 - np.exp(-0.32 * (dates + 1.0)))
    transfer_guess = initial_ss["transfer"] + (
        final_ss["transfer"] - initial_ss["transfer"]
    ) * (1.0 - np.exp(-0.32 * (dates + 1.0)))
    initial = np.concatenate([np.log(price_guess), np.log(transfer_guess)])
    fiscal_scale = P.tau_p * float(final_ss["price"]) * P.housing_stock

    def residual(
        log_paths: np.ndarray,
        return_details: bool = False,
        gain_branches: list[bool] | None = None,
    ):
        solved_prices = np.exp(log_paths[: horizon + 1])
        solved_transfers = np.exp(log_paths[horizon + 1 :])
        prices = np.concatenate(
            [
                solved_prices,
                np.array([final_ss["price"], final_ss["price"]]),
            ]
        )
        transfers = np.append(solved_transfers, final_ss["transfer"])

        young_mass = float(initial_ss["cohort_mass"])
        old_mass = float(initial_ss["cohort_mass"])
        previous_owner_share = float(initial_ss["young"]["owner_share"])
        previous_renter = initial_ss["young"]["renter"]
        previous_owner = initial_ss["young"]["owner"]
        lagged_price = float(initial_ss["price"])
        equations: list[float] = []
        rows: list[dict[str, float]] = []

        for date in range(horizon + 1):
            price = prices[date]
            next_price = prices[date + 1]
            price_after_next = prices[date + 2]
            transfer = transfers[date]
            next_transfer = transfers[date + 1]

            old_renter = old_renter_allocation(
                float(previous_renter["next_financial_wealth"]),
                price,
                next_price,
                transfer,
            )
            old_owner = old_owner_allocation(
                float(previous_owner["next_financial_wealth"]),
                float(previous_owner["housing"]),
                lagged_price,
                price,
                next_price,
                transfer,
                appreciates=(gain_branches[date] if gain_branches is not None else None),
            )
            average_old_housing = (
                (1.0 - previous_owner_share) * float(old_renter["housing"])
                + previous_owner_share * float(old_owner["housing"])
            )

            young = young_mixture(
                P.theta_final,
                price,
                next_price,
                price_after_next,
                transfer,
                next_transfer,
                next_appreciates=(
                    gain_branches[date + 1] if gain_branches is not None else None
                ),
            )
            owner_share = float(young["owner_share"])
            average_young_housing = float(young["average_housing"])
            average_fertility = float(young["average_fertility"])
            gain_tax = float(old_owner["gain_tax"])
            gain_revenue = (
                gain_tax
                * old_mass
                * previous_owner_share
                * float(old_owner["sales"])
            )
            property_revenue = P.tau_p * price * P.housing_stock
            housing_residual = (
                young_mass * average_young_housing
                + old_mass * average_old_housing
                - P.housing_stock
            ) / P.housing_stock
            government_residual = (
                (young_mass + old_mass) * transfer
                - property_revenue
                - gain_revenue
            ) / fiscal_scale
            equations.extend([housing_residual, government_residual])

            owner = young["owner"]
            renter = young["renter"]
            next_gain_tax = gain_tax_per_unit(
                next_price - price,
                gain_branches[date + 1] if gain_branches is not None else None,
            )
            gross_marginal_value_gap = (
                float(owner["access_wedge"]) + gain_tax + P.q * next_gain_tax
            )
            rows.append(
                {
                    "date": float(date),
                    "price": price,
                    "next_price": next_price,
                    "transfer": transfer,
                    "gain_tax_per_unit": gain_tax,
                    "next_gain_tax_per_unit": next_gain_tax,
                    "rent": float(old_owner["rent"]),
                    "incumbent_user_cost": float(old_owner["incumbent_cost"]),
                    "owner_user_cost": float(owner["user_cost"]),
                    "owner_access_wedge": float(owner["access_wedge"]),
                    "gross_marginal_value_gap": gross_marginal_value_gap,
                    "young_mass": young_mass,
                    "old_mass": old_mass,
                    "previous_owner_share": previous_owner_share,
                    "young_owner_share": owner_share,
                    "average_fertility": average_fertility,
                    "renter_fertility": float(renter["fertility"]),
                    "owner_fertility": float(owner["fertility"]),
                    "average_young_housing": average_young_housing,
                    "renter_housing": float(renter["housing"]),
                    "owner_housing": float(owner["housing"]),
                    "average_old_housing": average_old_housing,
                    "old_renter_housing": float(old_renter["housing"]),
                    "old_owner_housing": float(old_owner["housing"]),
                    "old_owner_sales": float(old_owner["sales"]),
                    "renter_saving": float(renter["saving"]),
                    "owner_saving": float(owner["saving"]),
                    "owner_unconstrained_housing": float(
                        owner["unconstrained_housing"]
                    ),
                    "owner_size_limit_slack": P.owner_cap
                    - float(owner["housing"]),
                    "renter_unconstrained_housing": float(
                        renter["unconstrained_housing"]
                    ),
                    "old_owner_retention_slack": float(old_owner["retention_slack"]),
                    "old_owner_estate_slack": float(old_owner["estate_slack"]),
                    "old_renter_cap_shadow_wedge": float(
                        old_renter["cap_shadow_wedge"]
                    ),
                    "housing_residual": housing_residual,
                    "government_residual": government_residual,
                }
            )

            next_young_mass = P.nu * average_fertility * young_mass
            old_mass = young_mass
            young_mass = next_young_mass
            previous_owner_share = owner_share
            previous_renter = renter
            previous_owner = owner
            lagged_price = price

        equation_array = np.asarray(equations)
        if return_details:
            terminal_state = {
                "young_mass": young_mass,
                "old_mass": old_mass,
                "previous_owner_share": previous_owner_share,
                "renter_financial_wealth": float(
                    previous_renter["next_financial_wealth"]
                ),
                "owner_financial_wealth": float(
                    previous_owner["next_financial_wealth"]
                ),
                "owner_housing": float(previous_owner["housing"]),
            }
            return equation_array, rows, terminal_state, log_paths
        return equation_array

    solution, iterations = damped_newton(
        residual,
        initial,
        tolerance=1.0e-9,
    )
    equations, rows, terminal_state, log_paths = residual(
        solution, return_details=True
    )
    solved_prices = np.exp(solution[: horizon + 1])
    price_path = np.concatenate(
        [solved_prices, np.array([final_ss["price"]])]
    )
    lagged_path = np.concatenate(
        [np.array([initial_ss["price"]]), solved_prices]
    )
    branch_price_changes = price_path - lagged_path
    gain_branches = [bool(change > 0.0) for change in branch_price_changes]
    branch_residual = lambda point: residual(
        point, gain_branches=gain_branches
    )
    jacobian = numerical_jacobian(branch_residual, solution)
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    return {
        "rows": rows,
        "terminal_state": terminal_state,
        "log_paths": log_paths,
        "iterations": iterations,
        "max_residual": float(np.max(np.abs(equations))),
        "jacobian_min_singular_value": float(np.min(singular_values)),
        "jacobian_condition_number": float(np.max(singular_values) / np.min(singular_values)),
        "minimum_absolute_branch_price_change": float(
            np.min(np.abs(branch_price_changes))
        ),
    }


def fertility_cap_derivative(
    theta: float,
    housing: float,
    fertility: float,
    wealth: float,
    user_cost: float,
) -> float:
    family_space = housing - P.kappa * fertility
    disposable = wealth - user_cost * housing - P.chi * fertility
    numerator = (
        P.alpha * P.kappa / family_space**2
        - (1.0 + P.k) * P.chi * user_cost / disposable**2
    )
    denominator = (
        theta / fertility**2
        + (1.0 + P.k) * P.chi**2 / disposable**2
        + P.alpha * P.kappa**2 / family_space**2
    )
    return numerator / denominator


def steady_state_row(label: str, values: dict) -> dict[str, float | str]:
    young = values["young"]
    renter = young["renter"]
    owner = young["owner"]
    return {
        "steady_state": label,
        "theta": values["theta"],
        "price": values["price"],
        "transfer": values["transfer"],
        "owner_share": young["owner_share"],
        "average_fertility": young["average_fertility"],
        "renter_fertility": renter["fertility"],
        "owner_fertility": owner["fertility"],
        "renter_housing": renter["housing"],
        "owner_housing": owner["housing"],
        "old_renter_housing": values["old_renter"]["housing"],
        "old_owner_housing": values["old_owner"]["housing"],
        "cohort_mass": values["cohort_mass"],
        "max_residual": values["max_residual"],
        "jacobian_determinant": values["steady_state_jacobian_determinant"],
    }


def write_figures(initial_ss: dict, final_ss: dict, transition_result: dict) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    blue = "#2B5C88"
    red = "#B1433E"
    green = "#3B755F"
    gray = "#666666"
    pale = "#E8EDF2"

    # Use date zero, when both current and anticipated realization taxes are positive.
    row = transition_result["rows"][0]
    incumbent_value = row["incumbent_user_cost"]
    market_cost = row["rent"]
    buyer_cost = row["owner_user_cost"]
    buyer_value = buyer_cost + row["owner_access_wedge"]
    fig, axis = plt.subplots(figsize=(6.7, 3.15), constrained_layout=True)
    values = [incumbent_value, market_cost, buyer_cost, buyer_value]
    labels = [
        r"Old incumbent: $r_t-\ell_t$",
        r"Market user cost: $r_t$",
        r"Young buyer's private cost: $r_t+q\ell_{t+1}$",
        r"Constrained young buyer: $r_t+q\ell_{t+1}+\zeta_t^{O,F}$",
    ]
    colors = [red, gray, green, blue]
    y_positions = np.arange(4.0)
    for y_position, value, color in zip(y_positions, values, colors):
        axis.hlines(y_position, min(values) - 0.012, value, color=color, linewidth=1.2)
        axis.scatter([value], [y_position], s=55, color=color, zorder=3)
    axis.annotate(
        "",
        xy=(buyer_value, 3.62),
        xytext=(incumbent_value, 3.62),
        arrowprops={"arrowstyle": "<->", "color": blue, "linewidth": 1.3},
    )
    axis.text(
        0.5 * (incumbent_value + buyer_value),
        3.73,
        r"gross marginal-value gap $=\zeta_t^{O,F}+\ell_t+q\ell_{t+1}$",
        ha="center",
        va="bottom",
        color=blue,
        fontsize=9,
    )
    axis.axvspan(incumbent_value, buyer_value, color=pale, alpha=0.50, zorder=0)
    axis.set_xlabel("Current-goods value of one housing unit")
    axis.set_yticks(y_positions, labels)
    axis.tick_params(axis="y", length=0)
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.5)
    axis.set_xlim(min(values) - 0.012, buyer_value + 0.018)
    axis.set_ylim(-0.35, 4.02)
    fig.savefig(WEDGE_FIGURE, dpi=260)
    plt.close(fig)

    rows = transition_result["rows"]
    dates = np.array([row_["date"] for row_ in rows])
    prices = np.array([row_["price"] for row_ in rows])
    gains = np.array([row_["gain_tax_per_unit"] for row_ in rows])
    fertility = np.array([row_["average_fertility"] for row_ in rows])
    young_mass = np.array([row_["young_mass"] for row_ in rows])
    old_mass = np.array([row_["old_mass"] for row_ in rows])
    owner_share = np.array([row_["young_owner_share"] for row_ in rows])

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 4.75), constrained_layout=True)
    axes[0, 0].plot(dates, prices, color=blue, linewidth=2.0)
    axes[0, 0].axhline(initial_ss["price"], color=gray, linestyle=":", linewidth=1)
    axes[0, 0].axhline(final_ss["price"], color=gray, linestyle="--", linewidth=1)
    axes[0, 0].set_title("House price")

    axes[0, 1].bar(dates, gains, color=red, width=0.72)
    axes[0, 1].set_title("Realized gains tax per sold unit")

    axes[1, 0].plot(dates, fertility, color=blue, linewidth=2.0)
    axes[1, 0].axhline(1.0 / P.nu, color=gray, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Average completed fertility")
    axes[1, 0].set_xlabel("Date after permanent shock")

    axes[1, 1].plot(dates, young_mass, color=blue, linewidth=2.0, label="Young mass")
    axes[1, 1].plot(
        dates, old_mass, color=red, linewidth=1.7, linestyle="--", label="Old mass"
    )
    share_axis = axes[1, 1].twinx()
    share_axis.plot(
        dates,
        owner_share,
        color=green,
        linewidth=1.5,
        linestyle=":",
        label="Owner share",
    )
    share_axis.set_ylabel("Owner share", color=green)
    share_axis.tick_params(axis="y", colors=green)
    share_axis.spines["top"].set_visible(False)
    axes[1, 1].set_title("Cohort masses and tenure")
    axes[1, 1].set_xlabel("Date after permanent shock")
    handles_1, labels_1 = axes[1, 1].get_legend_handles_labels()
    handles_2, labels_2 = share_axis.get_legend_handles_labels()
    axes[1, 1].legend(
        handles_1 + handles_2,
        labels_1 + labels_2,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.90,
        fontsize=7.5,
        loc="best",
    )

    for axis in axes.flat:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.5)
        axis.set_xlim(0, 12)
    fig.savefig(TRANSITION_FIGURE, dpi=260)
    plt.close(fig)


def write_outputs(
    initial_ss: dict,
    final_ss: dict,
    result_28: dict,
    horizon_stability_gap: float,
) -> dict[str, float | int | bool]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    steady_rows = [
        steady_state_row("initial", initial_ss),
        steady_state_row("final", final_ss),
    ]
    with (OUTPUT_DIR / "steady_states.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(steady_rows[0]))
        writer.writeheader()
        writer.writerows(steady_rows)
    rows = result_28["rows"]
    with (OUTPUT_DIR / "transition_path.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    final_young = final_ss["young"]
    terminal = result_28["terminal_state"]
    terminal_gaps = {
        "young_mass": abs(terminal["young_mass"] - final_ss["cohort_mass"]),
        "old_mass": abs(terminal["old_mass"] - final_ss["cohort_mass"]),
        "owner_share": abs(
            terminal["previous_owner_share"] - final_young["owner_share"]
        ),
        "renter_financial_wealth": abs(
            terminal["renter_financial_wealth"]
            - final_young["renter"]["next_financial_wealth"]
        ),
        "owner_financial_wealth": abs(
            terminal["owner_financial_wealth"]
            - final_young["owner"]["next_financial_wealth"]
        ),
        "owner_housing": abs(
            terminal["owner_housing"] - final_young["owner"]["housing"]
        ),
    }
    first_row = rows[0]
    owner = final_young["owner"]
    wealth = (
        P.income
        + P.liquid_wealth
        + final_ss["transfer"]
        + P.q * final_ss["transfer"]
    )
    analytic_derivative = fertility_cap_derivative(
        P.theta_final,
        float(owner["housing"]),
        float(owner["fertility"]),
        wealth,
        float(owner["user_cost"]),
    )
    step = 1.0e-5
    numerical_derivative = (
        solve_fertility(
            P.theta_final,
            float(owner["housing"]) + step,
            wealth,
            float(owner["user_cost"]),
        )
        - solve_fertility(
            P.theta_final,
            float(owner["housing"]) - step,
            wealth,
            float(owner["user_cost"]),
        )
    ) / (2.0 * step)
    verification = {
        "transition_horizon": 28,
        "transition_iterations": result_28["iterations"],
        "max_scaled_equilibrium_residual": result_28["max_residual"],
        "max_terminal_state_gap": max(terminal_gaps.values()),
        "terminal_state_gaps": terminal_gaps,
        "horizon_24_vs_28_first_nine_max_gap": horizon_stability_gap,
        "stacked_jacobian_min_singular_value": result_28[
            "jacobian_min_singular_value"
        ],
        "stacked_jacobian_condition_number": result_28[
            "jacobian_condition_number"
        ],
        "minimum_absolute_branch_price_change": result_28[
            "minimum_absolute_branch_price_change"
        ],
        "initial_ss_jacobian_determinant": initial_ss[
            "steady_state_jacobian_determinant"
        ],
        "final_ss_jacobian_determinant": final_ss[
            "steady_state_jacobian_determinant"
        ],
        "minimum_owner_share": min(row["young_owner_share"] for row in rows),
        "maximum_owner_share": max(row["young_owner_share"] for row in rows),
        "minimum_renter_share": min(1.0 - row["young_owner_share"] for row in rows),
        "minimum_owner_saving": min(row["owner_saving"] for row in rows),
        "minimum_renter_saving": min(row["renter_saving"] for row in rows),
        "minimum_old_owner_sales": min(row["old_owner_sales"] for row in rows),
        "minimum_incumbent_user_cost": min(
            row["incumbent_user_cost"] for row in rows
        ),
        "minimum_old_owner_estate_slack": min(
            row["old_owner_estate_slack"] for row in rows
        ),
        "minimum_owner_cap_bite": min(
            row["owner_unconstrained_housing"] - row["owner_housing"]
            for row in rows
        ),
        "minimum_owner_size_limit_slack": min(
            row["owner_size_limit_slack"] for row in rows
        ),
        "minimum_renter_cap_bite": min(
            row["renter_unconstrained_housing"] - row["renter_housing"]
            for row in rows
        ),
        "minimum_old_renter_cap_shadow_wedge": min(
            row["old_renter_cap_shadow_wedge"] for row in rows
        ),
        "minimum_gross_marginal_value_gap": min(
            row["gross_marginal_value_gap"] for row in rows
        ),
        "dates_with_positive_realized_gains_tax": sum(
            row["gain_tax_per_unit"] > 1.0e-8 for row in rows
        ),
        "fertility_cap_derivative_analytic": analytic_derivative,
        "fertility_cap_derivative_numerical": numerical_derivative,
        "fertility_derivative_absolute_error": abs(
            analytic_derivative - numerical_derivative
        ),
        "initial_owner_share": initial_ss["young"]["owner_share"],
        "final_owner_share": final_ss["young"]["owner_share"],
        "initial_price": initial_ss["price"],
        "final_price": final_ss["price"],
        "impact_price": first_row["price"],
    }
    with (OUTPUT_DIR / "verification.json").open("w") as stream:
        json.dump(
            {"parameters": asdict(P), "verification": verification},
            stream,
            indent=2,
        )
    write_figures(initial_ss, final_ss, result_28)
    return verification


def compare_horizons(result_24: dict, result_28: dict) -> float:
    keys = [
        "price",
        "transfer",
        "young_mass",
        "old_mass",
        "young_owner_share",
        "average_fertility",
        "average_young_housing",
        "average_old_housing",
    ]
    gaps = []
    for date in range(9):
        for key in keys:
            gaps.append(abs(result_24["rows"][date][key] - result_28["rows"][date][key]))
    return float(max(gaps))


def main() -> None:
    initial_ss = solve_steady_state(P.theta_initial, 0.34)
    final_ss = solve_steady_state(P.theta_final, 0.48)
    result_24 = transition(initial_ss, final_ss, 24)
    result_28 = transition(initial_ss, final_ss, 28)
    horizon_gap = compare_horizons(result_24, result_28)
    verification = write_outputs(initial_ss, final_ss, result_28, horizon_gap)
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
