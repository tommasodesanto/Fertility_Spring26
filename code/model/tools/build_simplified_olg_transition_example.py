#!/usr/bin/env python3
"""Solve and plot a reproducible special case of the simplified OLG model.

The calculation is deliberately small.  It conditions on ownership, a binding
down-payment limit, and interior old-age choices.  Its purpose is to verify the
transition equations in the technical appendix and to show that one permanent
preference shock can generate capital gains on several consecutive dates.

Outputs
-------
latex/figures/simplified_olg_transition_example.png
output/model/simplified_olg_transition_example/transition_path.csv
output/model/simplified_olg_transition_example/verification.json
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
FIGURE_PATH = ROOT / "latex" / "figures" / "simplified_olg_transition_example.png"
OUTPUT_DIR = ROOT / "output" / "model" / "simplified_olg_transition_example"


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


def solve_fertility(theta: float, housing: float, wealth: float, user_cost: float) -> float:
    """Solve the strictly decreasing binding-cap fertility first-order condition."""
    lower = 1.0e-10
    upper = min(
        housing / P.kappa,
        (wealth - user_cost * housing) / P.chi,
    ) * (1.0 - 1.0e-9)
    if upper <= lower:
        raise ValueError("The binding-cap allocation has an empty log domain.")

    def residual(fertility: float) -> float:
        family_space = housing - P.kappa * fertility
        disposable = wealth - user_cost * housing - P.chi * fertility
        return (
            theta / fertility
            - P.alpha * P.kappa / family_space
            - (1.0 + P.k) * P.chi / disposable
        )

    if not (residual(lower) > 0.0 and residual(upper) < 0.0):
        raise ValueError("The fertility first-order condition is not bracketed.")
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def numerical_jacobian(function, point: np.ndarray, step: float = 1.0e-5) -> np.ndarray:
    baseline = function(point)
    jacobian = np.empty((baseline.size, point.size))
    for column in range(point.size):
        shifted = point.copy()
        shifted[column] += step
        jacobian[:, column] = (function(shifted) - baseline) / step
    return jacobian


def damped_newton(function, initial: np.ndarray, tolerance: float = 1.0e-10) -> tuple[np.ndarray, int]:
    point = initial.copy()
    for iteration in range(60):
        residual = function(point)
        norm = float(np.max(np.abs(residual)))
        if norm < tolerance:
            return point, iteration
        step = np.linalg.solve(numerical_jacobian(function, point), -residual)
        damping = 1.0
        accepted = False
        for _ in range(30):
            candidate = point + damping * step
            try:
                candidate_norm = float(np.max(np.abs(function(candidate))))
            except (FloatingPointError, ValueError):
                candidate_norm = np.inf
            if candidate_norm < norm:
                point = candidate
                accepted = True
                break
            damping *= 0.5
        if not accepted:
            raise RuntimeError(f"Newton line search failed at residual {norm:.3e}.")
    raise RuntimeError("Newton iteration limit reached.")


def steady_state(theta: float, initial_price: float) -> dict[str, float]:
    def residual(log_price_transfer: np.ndarray) -> np.ndarray:
        price, transfer = np.exp(log_price_transfer)
        values = steady_state_values(theta, price, transfer)
        return np.array(
            [
                P.nu * values["fertility"] - 1.0,
                transfer
                - 0.5
                * P.tau_p
                * price
                * (values["young_housing"] + values["old_housing"]),
            ]
        )

    solution, iterations = damped_newton(
        residual,
        np.log(np.array([initial_price, 0.013])),
    )
    price, transfer = np.exp(solution)
    values = steady_state_values(theta, price, transfer)
    values["iterations"] = float(iterations)
    values["max_residual"] = float(np.max(np.abs(residual(solution))))
    return values


def steady_state_values(theta: float, price: float, transfer: float) -> dict[str, float]:
    rent = (1.0 + P.tau_p - P.q) * price
    young_housing = P.liquid_wealth / ((1.0 - P.phi) * price)
    lifetime_wealth = P.income + P.liquid_wealth + (1.0 + P.q) * transfer
    fertility = solve_fertility(theta, young_housing, lifetime_wealth, rent)
    disposable = lifetime_wealth - rent * young_housing - P.chi * fertility
    consumption = disposable / (1.0 + P.k)
    saving = (
        P.k * consumption
        - P.q * transfer
        - (P.q * price - P.phi * price) * young_housing
    )
    old_financial_wealth = (
        P.gross_return * saving
        - P.gross_return * P.phi * price * young_housing
    )
    old_resources = old_financial_wealth + price * young_housing + transfer
    old_housing = P.gamma * old_resources / (P.old_share_sum * rent)
    estate = P.omega_b * old_resources / (P.old_share_sum * P.q)
    cohort = P.housing_stock / (young_housing + old_housing)
    unconstrained_housing = unconstrained_young_housing(
        theta, lifetime_wealth, rent
    )
    return {
        "price": price,
        "transfer": transfer,
        "rent": rent,
        "fertility": fertility,
        "young_housing": young_housing,
        "old_housing": old_housing,
        "old_estate": estate,
        "old_financial_wealth": old_financial_wealth,
        "cohort": cohort,
        "unconstrained_young_housing": unconstrained_housing,
    }


def unconstrained_young_housing(theta: float, wealth: float, user_cost: float) -> float:
    share_sum = 1.0 + P.k + P.alpha + theta
    family_space = P.alpha * wealth / (share_sum * user_cost)
    fertility = theta * wealth / (
        share_sum * (P.chi + P.kappa * user_cost)
    )
    return family_space + P.kappa * fertility


def transition(initial_ss: dict[str, float], final_ss: dict[str, float], horizon: int):
    price_guess = initial_ss["price"] + (
        final_ss["price"] - initial_ss["price"]
    ) * (1.0 - np.exp(-0.30 * np.arange(horizon + 1)))
    transfer_guess = np.full(horizon + 1, final_ss["transfer"])
    initial = np.concatenate([np.log(price_guess), np.log(transfer_guess)])

    def residual(log_paths: np.ndarray, return_rows: bool = False):
        price = np.append(
            np.exp(log_paths[: horizon + 1]), final_ss["price"]
        )
        transfer = np.append(
            np.exp(log_paths[horizon + 1 :]), final_ss["transfer"]
        )
        young_mass = initial_ss["cohort"]
        old_mass = initial_ss["cohort"]
        old_financial_wealth = initial_ss["old_financial_wealth"]
        inherited_housing = initial_ss["young_housing"]
        lagged_price = initial_ss["price"]
        equations: list[float] = []
        rows: list[dict[str, float]] = []

        for date in range(horizon + 1):
            gain_tax = P.tau_g * max(price[date] - lagged_price, 0.0)
            next_gain_tax = P.tau_g * max(
                price[date + 1] - price[date], 0.0
            )
            rent = (1.0 + P.tau_p) * price[date] - P.q * price[date + 1]
            incumbent_cost = rent - gain_tax
            old_resources = (
                old_financial_wealth
                + (price[date] - gain_tax) * inherited_housing
                + transfer[date]
            )
            if incumbent_cost <= 0.0 or old_resources <= 0.0:
                raise ValueError("Old-age interior conditions fail.")
            old_housing = (
                P.gamma * old_resources / (P.old_share_sum * incumbent_cost)
            )
            estate = P.omega_b * old_resources / (P.old_share_sum * P.q)
            old_sales = inherited_housing - old_housing

            young_housing = P.liquid_wealth / (
                (1.0 - P.phi) * price[date]
            )
            owner_cost = rent + P.q * next_gain_tax
            lifetime_wealth = (
                P.income
                + P.liquid_wealth
                + transfer[date]
                + P.q * transfer[date + 1]
            )
            fertility = solve_fertility(
                P.theta_final, young_housing, lifetime_wealth, owner_cost
            )
            disposable = (
                lifetime_wealth
                - owner_cost * young_housing
                - P.chi * fertility
            )
            consumption = disposable / (1.0 + P.k)
            saving = (
                P.k * consumption
                - P.q * transfer[date + 1]
                - (
                    P.q * (price[date + 1] - next_gain_tax)
                    - P.phi * price[date]
                )
                * young_housing
            )
            next_old_financial_wealth = (
                P.gross_return * saving
                - P.gross_return * P.phi * price[date] * young_housing
            )
            unconstrained_housing = unconstrained_young_housing(
                P.theta_final, lifetime_wealth, owner_cost
            )

            housing_residual = (
                young_mass * young_housing
                + old_mass * old_housing
                - P.housing_stock
            )
            government_residual = (
                (young_mass + old_mass) * transfer[date]
                - P.tau_p * price[date] * P.housing_stock
                - gain_tax * old_mass * old_sales
            )
            equations.extend([housing_residual, government_residual])
            rows.append(
                {
                    "date": float(date),
                    "price": price[date],
                    "next_price": price[date + 1],
                    "transfer": transfer[date],
                    "gain_tax_per_unit": gain_tax,
                    "rent": rent,
                    "incumbent_user_cost": incumbent_cost,
                    "owner_user_cost": owner_cost,
                    "lifetime_wealth": lifetime_wealth,
                    "young_mass": young_mass,
                    "old_mass": old_mass,
                    "fertility": fertility,
                    "young_housing": young_housing,
                    "old_housing": old_housing,
                    "inherited_housing": inherited_housing,
                    "old_sales": old_sales,
                    "old_estate": estate,
                    "saving": saving,
                    "unconstrained_young_housing": unconstrained_housing,
                    "housing_residual": housing_residual,
                    "government_residual": government_residual,
                }
            )

            next_young_mass = P.nu * fertility * young_mass
            old_mass = young_mass
            young_mass = next_young_mass
            old_financial_wealth = next_old_financial_wealth
            inherited_housing = young_housing
            lagged_price = price[date]

        equations_array = np.asarray(equations)
        if return_rows:
            terminal_state = {
                "young_mass": young_mass,
                "old_mass": old_mass,
                "old_financial_wealth": old_financial_wealth,
                "inherited_housing": inherited_housing,
            }
            return equations_array, rows, terminal_state
        return equations_array

    solution, iterations = damped_newton(residual, initial, tolerance=1.0e-9)
    equations, rows, terminal_state = residual(solution, return_rows=True)
    return {
        "rows": rows,
        "terminal_state": terminal_state,
        "max_residual": float(np.max(np.abs(equations))),
        "iterations": iterations,
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


def write_outputs(
    initial_ss: dict[str, float],
    final_ss: dict[str, float],
    result,
    horizon_stability_gap: float,
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)

    rows = result["rows"]
    with (OUTPUT_DIR / "transition_path.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    dates = np.array([row["date"] for row in rows])
    prices = np.array([row["price"] for row in rows])
    gains = np.array([row["gain_tax_per_unit"] for row in rows])
    fertility = np.array([row["fertility"] for row in rows])
    young_mass = np.array([row["young_mass"] for row in rows])
    old_mass = np.array([row["old_mass"] for row in rows])

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.7), constrained_layout=True)
    blue = "#2B5C88"
    red = "#B1433E"
    gray = "#666666"

    axes[0, 0].plot(dates, prices, color=blue, linewidth=2.0)
    axes[0, 0].axhline(initial_ss["price"], color=gray, linestyle=":", linewidth=1)
    axes[0, 0].axhline(final_ss["price"], color=gray, linestyle="--", linewidth=1)
    axes[0, 0].set_title("House price")

    axes[0, 1].bar(dates, gains, color=red, width=0.72)
    axes[0, 1].set_title("Capital-gains tax per unit")

    axes[1, 0].plot(dates, fertility, color=blue, linewidth=2.0)
    axes[1, 0].axhline(1.0 / P.nu, color=gray, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Completed fertility")
    axes[1, 0].set_xlabel("Date after permanent shock")

    axes[1, 1].plot(dates, young_mass, color=blue, linewidth=2.0, label="Young")
    axes[1, 1].plot(dates, old_mass, color=red, linewidth=1.7, linestyle="--", label="Old")
    axes[1, 1].set_title("Cohort masses")
    axes[1, 1].set_xlabel("Date after permanent shock")
    axes[1, 1].legend(frameon=False, fontsize=8)

    for axis in axes.flat:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.5)
        axis.set_xlim(0, 12)
    fig.savefig(FIGURE_PATH, dpi=240)
    plt.close(fig)

    target_terminal = {
        "young_mass": final_ss["cohort"],
        "old_mass": final_ss["cohort"],
        "old_financial_wealth": final_ss["old_financial_wealth"],
        "inherited_housing": final_ss["young_housing"],
    }
    terminal_gap = max(
        abs(result["terminal_state"][key] - target_terminal[key])
        for key in target_terminal
    )
    derivative_errors = []
    for row in rows:
        housing_step = 1.0e-6 * row["young_housing"]
        numerical_derivative = (
            solve_fertility(
                P.theta_final,
                row["young_housing"] + housing_step,
                row["lifetime_wealth"],
                row["owner_user_cost"],
            )
            - solve_fertility(
                P.theta_final,
                row["young_housing"] - housing_step,
                row["lifetime_wealth"],
                row["owner_user_cost"],
            )
        ) / (2.0 * housing_step)
        analytical_derivative = fertility_cap_derivative(
            P.theta_final,
            row["young_housing"],
            row["fertility"],
            row["lifetime_wealth"],
            row["owner_user_cost"],
        )
        derivative_errors.append(abs(numerical_derivative - analytical_derivative))
    checks = {
        "max_equilibrium_residual": result["max_residual"],
        "max_terminal_state_gap": terminal_gap,
        "max_horizon_20_24_gap_first_9_dates": horizon_stability_gap,
        "max_fertility_derivative_fd_error": max(derivative_errors),
        "minimum_old_sales": min(row["old_sales"] for row in rows),
        "minimum_incumbent_user_cost": min(
            row["incumbent_user_cost"] for row in rows
        ),
        "minimum_estate_slack": min(
            row["old_estate"]
            - row["next_price"] * row["old_housing"]
            for row in rows
        ),
        "minimum_saving": min(row["saving"] for row in rows),
        "minimum_down_payment_cap_slack": min(
            row["unconstrained_young_housing"] - row["young_housing"]
            for row in rows
        ),
        "positive_gain_tax_dates": [
            int(row["date"])
            for row in rows
            if row["gain_tax_per_unit"] > 1.0e-4
        ],
    }
    payload = {
        "parameters": asdict(P),
        "initial_steady_state": initial_ss,
        "final_steady_state": final_ss,
        "transition_iterations": result["iterations"],
        "checks": checks,
    }
    (OUTPUT_DIR / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )

    if checks["max_equilibrium_residual"] >= 1.0e-8:
        raise RuntimeError("Transition equilibrium residual is too large.")
    if checks["max_terminal_state_gap"] >= 1.0e-6:
        raise RuntimeError("The terminal state has not converged.")
    if checks["max_horizon_20_24_gap_first_9_dates"] >= 1.0e-8:
        raise RuntimeError("The transition is not stable to the terminal horizon.")
    if checks["max_fertility_derivative_fd_error"] >= 1.0e-6:
        raise RuntimeError("The fertility derivative fails its finite-difference check.")
    if checks["minimum_old_sales"] <= 0.0:
        raise RuntimeError("The old retention constraint binds in the example.")
    if checks["minimum_estate_slack"] <= 0.0:
        raise RuntimeError("The old estate-composition constraint binds.")
    if checks["minimum_saving"] <= 0.0:
        raise RuntimeError("The young saving constraint binds.")
    if checks["minimum_down_payment_cap_slack"] <= 0.0:
        raise RuntimeError("The down-payment cap is not binding.")


def main() -> None:
    np.seterr(over="raise", invalid="raise", divide="raise")
    initial_ss = steady_state(P.theta_initial, initial_price=0.40)
    final_ss = steady_state(P.theta_final, initial_price=0.56)
    shorter_result = transition(initial_ss, final_ss, horizon=20)
    result = transition(initial_ss, final_ss, horizon=24)
    horizon_stability_gap = max(
        abs(shorter_result["rows"][date][key] - result["rows"][date][key])
        for date in range(9)
        for key in ("price", "transfer", "fertility", "young_mass")
    )
    write_outputs(
        initial_ss,
        final_ss,
        result,
        horizon_stability_gap,
    )
    print(f"Wrote {FIGURE_PATH}")
    print(f"Wrote {OUTPUT_DIR / 'transition_path.csv'}")
    print(f"Wrote {OUTPUT_DIR / 'verification.json'}")


if __name__ == "__main__":
    main()
