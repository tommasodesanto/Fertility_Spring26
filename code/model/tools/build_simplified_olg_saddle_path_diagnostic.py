#!/usr/bin/env python3
"""Build the local saddle-path diagnostic for the simplified OLG example.

The script imports the current construction driver, solves both steady states,
forms the implicit first-order descriptor system, and writes
``saddle_path_diagnostic.json`` beside the other theory outputs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig, ordqz


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MODEL_DIR = PROJECT_ROOT / "code" / "model" / "tools"
OUTPUT_DIR = (
    PROJECT_ROOT / "output" / "model" / "simplified_olg_mixed_tenure_theory"
)
OUTPUT_PATH = OUTPUT_DIR / "saddle_path_diagnostic.json"
sys.path.insert(0, str(MODEL_DIR))

import build_simplified_olg_mixed_tenure_theory as theory  # noqa: E402


STEP = 1.0e-5
STATE_NAMES = [
    "young_mass",
    "old_mass",
    "previous_owner_share",
    "old_renter_financial_wealth",
    "old_owner_financial_wealth",
    "old_owner_inherited_housing",
    "lagged_price_basis",
]
FORWARD_PRICE_NAMES = ["current_price", "next_price_auxiliary"]
STATIC_VARIABLE_NAMES = ["current_transfer"]
GOVERNMENT_ROW = 8
TRANSFER_COLUMN = 8


def log_stationary_vector(steady_state: dict) -> np.ndarray:
    young = steady_state["young"]
    return np.log(
        np.array(
            [
                steady_state["cohort_mass"],
                steady_state["cohort_mass"],
                young["owner_share"],
                young["renter"]["next_financial_wealth"],
                young["owner"]["next_financial_wealth"],
                young["owner"]["housing"],
                steady_state["price"],
                steady_state["price"],
                steady_state["transfer"],
                steady_state["price"],
            ],
            dtype=float,
        )
    )


def residual(
    current_logs: np.ndarray,
    next_logs: np.ndarray,
    theta: float,
    reference_price: float,
    current_gain_derivative: bool,
    next_gain_derivative: bool,
) -> np.ndarray:
    """One-date descriptor residual in log coordinates.

    Current and next vectors contain seven predetermined states followed by
    (P_t, T_t, Q_t), where Q_t is an auxiliary copy of P_{t+1}.  The equations
    are seven state laws, housing clearing, the government budget, and
    log(Q_t)-log(P_{t+1})=0.
    """

    current = np.exp(current_logs)
    nxt = np.exp(next_logs)
    (
        young_mass,
        old_mass,
        previous_owner_share,
        renter_financial_wealth,
        owner_financial_wealth,
        inherited_owner_housing,
        lagged_price,
        price,
        transfer,
        next_price_auxiliary,
    ) = current
    (
        next_young_mass,
        next_old_mass,
        next_previous_owner_share,
        next_renter_financial_wealth,
        next_owner_financial_wealth,
        next_inherited_owner_housing,
        next_lagged_price,
        next_price,
        next_transfer,
        price_after_next_auxiliary,
    ) = nxt

    old_renter = theory.old_renter_allocation(
        renter_financial_wealth,
        price,
        next_price_auxiliary,
        transfer,
    )
    old_owner = theory.old_owner_allocation(
        owner_financial_wealth,
        inherited_owner_housing,
        lagged_price,
        price,
        next_price_auxiliary,
        transfer,
        appreciates=current_gain_derivative,
    )
    young = theory.young_mixture(
        theta,
        price,
        next_price_auxiliary,
        price_after_next_auxiliary,
        transfer,
        next_transfer,
        next_appreciates=next_gain_derivative,
    )

    predicted_next_state = np.array(
        [
            theory.P.nu * young["average_fertility"] * young_mass,
            young_mass,
            young["owner_share"],
            young["renter"]["next_financial_wealth"],
            young["owner"]["next_financial_wealth"],
            young["owner"]["housing"],
            price,
        ],
        dtype=float,
    )
    supplied_next_state = np.array(
        [
            next_young_mass,
            next_old_mass,
            next_previous_owner_share,
            next_renter_financial_wealth,
            next_owner_financial_wealth,
            next_inherited_owner_housing,
            next_lagged_price,
        ],
        dtype=float,
    )

    average_old_housing = (
        (1.0 - previous_owner_share) * old_renter["housing"]
        + previous_owner_share * old_owner["housing"]
    )
    housing_residual = (
        young_mass * young["average_housing"]
        + old_mass * average_old_housing
        - theory.P.housing_stock
    ) / theory.P.housing_stock

    gains_revenue = (
        old_owner["gain_tax"]
        * old_mass
        * previous_owner_share
        * old_owner["sales"]
    )
    fiscal_scale = (
        theory.P.q
        * theory.P.tau_p
        * reference_price
        * theory.P.housing_stock
    )
    government_residual = (
        (young_mass + old_mass) * transfer
        - theory.P.q * theory.P.tau_p * price * theory.P.housing_stock
        - gains_revenue
    ) / fiscal_scale

    auxiliary_price_residual = np.log(next_price_auxiliary) - np.log(next_price)
    return np.concatenate(
        [
            np.log(supplied_next_state) - np.log(predicted_next_state),
            np.array(
                [housing_residual, government_residual, auxiliary_price_residual]
            ),
        ]
    )


def centered_jacobian(function, point: np.ndarray) -> np.ndarray:
    baseline = function(point)
    jacobian = np.empty((baseline.size, point.size))
    for column in range(point.size):
        higher = point.copy()
        lower = point.copy()
        higher[column] += STEP
        lower[column] -= STEP
        jacobian[:, column] = (function(higher) - function(lower)) / (2.0 * STEP)
    return jacobian


def pencil(
    steady_state: dict,
    current_gain_derivative: bool,
    next_gain_derivative: bool,
) -> tuple[np.ndarray, np.ndarray]:
    point = log_stationary_vector(steady_state)
    theta = float(steady_state["theta"])
    reference_price = float(steady_state["price"])
    current_jacobian = centered_jacobian(
        lambda current: residual(
            current,
            point,
            theta,
            reference_price,
            current_gain_derivative,
            next_gain_derivative,
        ),
        point,
    )
    next_jacobian = centered_jacobian(
        lambda nxt: residual(
            point,
            nxt,
            theta,
            reference_price,
            current_gain_derivative,
            next_gain_derivative,
        ),
        point,
    )
    return current_jacobian, next_jacobian


def analyze_pencil(current_jacobian: np.ndarray, next_jacobian: np.ndarray) -> dict:
    eigenvalues = eig(-current_jacobian, next_jacobian, right=False)
    finite_mask = np.isfinite(eigenvalues)
    finite_values = eigenvalues[finite_mask]
    stable_mask = np.abs(finite_values) < 1.0 - 1.0e-7
    unit_mask = np.abs(np.abs(finite_values) - 1.0) <= 1.0e-7
    unstable_mask = np.abs(finite_values) > 1.0 + 1.0e-7
    sorted_values = sorted(finite_values, key=abs)

    # The ordered generalized Schur vectors form an orthonormal basis for the
    # stable deflating subspace.  Its projection on the first seven coordinates
    # tests whether the stable manifold is locally a graph over predetermined
    # states.
    _, _, alpha, beta, _, right_schur_vectors = ordqz(
        -current_jacobian,
        next_jacobian,
        sort="iuc",
        output="complex",
    )
    ordered_values = np.full(alpha.shape, np.inf + 0.0j, dtype=complex)
    finite_beta = np.abs(beta) > 1.0e-14
    ordered_values[finite_beta] = alpha[finite_beta] / beta[finite_beta]
    stable_count = int(np.sum(np.isfinite(ordered_values) & (np.abs(ordered_values) < 1.0)))
    stable_basis = right_schur_vectors[:, :stable_count]
    projection_singular_values = np.linalg.svd(
        stable_basis[: len(STATE_NAMES), :], compute_uv=False
    )

    return {
        "finite_eigenvalues": [
            {
                "real": float(value.real),
                "imaginary": float(value.imag),
                "modulus": float(abs(value)),
            }
            for value in sorted_values
        ],
        "stable_count": int(np.sum(stable_mask)),
        "unit_count": int(np.sum(unit_mask)),
        "finite_unstable_count": int(np.sum(unstable_mask)),
        "infinite_count": int(np.sum(~finite_mask)),
        "largest_stable_modulus": float(np.max(np.abs(finite_values[stable_mask]))),
        "smallest_finite_unstable_modulus": float(
            np.min(np.abs(finite_values[unstable_mask]))
        ),
        "stable_projection_singular_values": [
            float(value) for value in projection_singular_values
        ],
        "stable_projection_min_singular_value": float(
            np.min(projection_singular_values)
        ),
        "stable_projection_condition_number": float(
            np.max(projection_singular_values)
            / np.min(projection_singular_values)
        ),
    }


def eliminate_static_transfer(
    current_jacobian: np.ndarray,
    next_jacobian: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Eliminate the date-t and date-(t+1) transfers using the fiscal row."""

    keep_rows = [row for row in range(current_jacobian.shape[0]) if row != GOVERNMENT_ROW]
    keep_columns = [
        column
        for column in range(current_jacobian.shape[1])
        if column != TRANSFER_COLUMN
    ]
    fiscal_current = current_jacobian[GOVERNMENT_ROW, :]
    fiscal_transfer_derivative = fiscal_current[TRANSFER_COLUMN]
    if abs(fiscal_transfer_derivative) <= 1.0e-10:
        raise ValueError("The government budget does not locally determine the transfer.")

    transfer_slope = (
        fiscal_current[keep_columns] / fiscal_transfer_derivative
    )
    current_reduced = (
        current_jacobian[np.ix_(keep_rows, keep_columns)]
        - np.outer(
            current_jacobian[keep_rows, TRANSFER_COLUMN],
            transfer_slope,
        )
    )
    next_reduced = (
        next_jacobian[np.ix_(keep_rows, keep_columns)]
        - np.outer(
            next_jacobian[keep_rows, TRANSFER_COLUMN],
            transfer_slope,
        )
    )
    diagnostics = {
        "full_next_jacobian_rank": int(np.linalg.matrix_rank(next_jacobian)),
        "government_next_jacobian_max_abs": float(
            np.max(np.abs(next_jacobian[GOVERNMENT_ROW, :]))
        ),
        "government_current_transfer_derivative": float(
            fiscal_transfer_derivative
        ),
        "reduced_next_jacobian_rank": int(np.linalg.matrix_rank(next_reduced)),
    }
    return current_reduced, next_reduced, diagnostics


def generalized_derivative_grid(corner_pencils: dict) -> dict:
    """Interpolate over Clarke derivatives of the two positive-part terms."""

    grid = np.linspace(0.0, 1.0, 21)
    count_tuples = set()
    largest_stable_modulus = 0.0
    smallest_unstable_modulus = np.inf
    smallest_distance_to_unit_circle = np.inf
    smallest_projection_singular_value = np.inf

    for current_weight in grid:
        for next_weight in grid:
            weights = {
                (0, 0): (1.0 - current_weight) * (1.0 - next_weight),
                (1, 0): current_weight * (1.0 - next_weight),
                (0, 1): (1.0 - current_weight) * next_weight,
                (1, 1): current_weight * next_weight,
            }
            current_jacobian = sum(
                weights[key] * corner_pencils[key][0] for key in weights
            )
            next_jacobian = sum(
                weights[key] * corner_pencils[key][1] for key in weights
            )
            result = analyze_pencil(current_jacobian, next_jacobian)
            count_tuples.add(
                (
                    result["stable_count"],
                    result["finite_unstable_count"],
                    result["infinite_count"],
                )
            )
            largest_stable_modulus = max(
                largest_stable_modulus, result["largest_stable_modulus"]
            )
            smallest_unstable_modulus = min(
                smallest_unstable_modulus,
                result["smallest_finite_unstable_modulus"],
            )
            smallest_projection_singular_value = min(
                smallest_projection_singular_value,
                result["stable_projection_min_singular_value"],
            )
            finite_moduli = np.array(
                [row["modulus"] for row in result["finite_eigenvalues"]]
            )
            smallest_distance_to_unit_circle = min(
                smallest_distance_to_unit_circle,
                float(np.min(np.abs(finite_moduli - 1.0))),
            )

    return {
        "grid_size_per_derivative": int(grid.size),
        "count_tuples_stable_finite_unstable_infinite": [
            list(values) for values in sorted(count_tuples)
        ],
        "largest_stable_modulus": largest_stable_modulus,
        "smallest_finite_unstable_modulus": smallest_unstable_modulus,
        "smallest_distance_to_unit_circle": smallest_distance_to_unit_circle,
        "smallest_stable_projection_singular_value": smallest_projection_singular_value,
    }


def main() -> None:
    endpoints = {
        "initial": theory.solve_steady_state(theory.P.theta_initial, 0.34),
        "final": theory.solve_steady_state(theory.P.theta_final, 0.48),
    }
    report = {
        "model_source": str(MODEL_DIR / "build_simplified_olg_mixed_tenure_theory.py"),
        "finite_difference_step_in_logs": STEP,
        "predetermined_states": STATE_NAMES,
        "forward_price_coordinates": FORWARD_PRICE_NAMES,
        "static_variables": STATIC_VARIABLE_NAMES,
        "descriptor_equations": [
            "seven state laws",
            "housing clearing",
            "government budget",
            "next-price auxiliary identity",
        ],
        "endpoints": {},
    }

    for endpoint_name, steady_state in endpoints.items():
        selections = {}
        corner_pencils = {}
        for current_derivative in (0, 1):
            for next_derivative in (0, 1):
                key = (current_derivative, next_derivative)
                current_jacobian, next_jacobian = pencil(
                    steady_state,
                    bool(current_derivative),
                    bool(next_derivative),
                )
                corner_pencils[key] = (current_jacobian, next_jacobian)
                full_analysis = analyze_pencil(
                    current_jacobian, next_jacobian
                )
                reduced_current, reduced_next, elimination = (
                    eliminate_static_transfer(current_jacobian, next_jacobian)
                )
                selections[f"current_{current_derivative}_next_{next_derivative}"] = {
                    "descriptor_pencil": full_analysis,
                    "after_static_transfer_elimination": analyze_pencil(
                        reduced_current, reduced_next
                    ),
                    "elimination_diagnostics": elimination,
                }
        current_additivity_error = float(
            np.linalg.norm(
                corner_pencils[(1, 1)][0]
                - corner_pencils[(1, 0)][0]
                - corner_pencils[(0, 1)][0]
                + corner_pencils[(0, 0)][0],
                ord=np.inf,
            )
        )
        next_additivity_error = float(
            np.linalg.norm(
                corner_pencils[(1, 1)][1]
                - corner_pencils[(1, 0)][1]
                - corner_pencils[(0, 1)][1]
                + corner_pencils[(0, 0)][1],
                ord=np.inf,
            )
        )
        generalized_grid = generalized_derivative_grid(corner_pencils)
        report["endpoints"][endpoint_name] = {
            "theta": float(steady_state["theta"]),
            "price": float(steady_state["price"]),
            "transfer": float(steady_state["transfer"]),
            "gain_derivative_selections": selections,
            "corner_additivity_error": {
                "current_jacobian_inf_norm": current_additivity_error,
                "next_jacobian_inf_norm": next_additivity_error,
            },
            "generalized_derivative_grid": generalized_grid,
        }

    for endpoint in report["endpoints"].values():
        for selection in endpoint["gain_derivative_selections"].values():
            descriptor = selection["descriptor_pencil"]
            reduced = selection["after_static_transfer_elimination"]
            elimination = selection["elimination_diagnostics"]
            assert descriptor["stable_count"] == len(STATE_NAMES)
            assert descriptor["unit_count"] == 0
            assert descriptor["finite_unstable_count"] == 2
            assert descriptor["infinite_count"] == 1
            assert descriptor["stable_projection_min_singular_value"] > 1.0e-4
            assert reduced["stable_count"] == len(STATE_NAMES)
            assert reduced["unit_count"] == 0
            assert reduced["finite_unstable_count"] == 2
            assert reduced["infinite_count"] == 0
            assert elimination["full_next_jacobian_rank"] == 9
            assert elimination["government_next_jacobian_max_abs"] < 1.0e-9
            assert abs(elimination["government_current_transfer_derivative"]) > 1.0e-4
            assert elimination["reduced_next_jacobian_rank"] == 9
        generalized_grid = endpoint["generalized_derivative_grid"]
        assert generalized_grid[
            "count_tuples_stable_finite_unstable_infinite"
        ] == [[7, 2, 1]]
        assert generalized_grid["smallest_stable_projection_singular_value"] > 1.0e-4
        assert endpoint["corner_additivity_error"][
            "current_jacobian_inf_norm"
        ] < 1.0e-8
        assert endpoint["corner_additivity_error"][
            "next_jacobian_inf_norm"
        ] < 1.0e-8

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
