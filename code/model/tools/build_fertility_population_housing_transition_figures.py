"""Build the two figures for the fertility-population-housing transition note.

The illustration solves the exact linear two-state system stated in
``latex/fertility_population_housing_transition_note.tex``.  It is a
transparent thought experiment, not a calibration or a model counterfactual.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = REPO_ROOT / "output" / "model" / "fertility_population_housing_transition_note"


def matrix_exponential_2x2(matrix: np.ndarray, time: float) -> np.ndarray:
    """Return exp(matrix * time) using the two distinct real eigenvalues."""

    roots = np.linalg.eigvals(matrix)
    roots = np.sort(np.real_if_close(roots).astype(float))
    slow, fast = roots[1], roots[0]
    if np.isclose(slow, fast):
        identity = np.eye(2)
        return np.exp(slow * time) * (identity + (matrix - slow * identity) * time)
    identity = np.eye(2)
    return (
        np.exp(slow * time) * (matrix - fast * identity)
        - np.exp(fast * time) * (matrix - slow * identity)
    ) / (slow - fast)


def solve_transition() -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    """Solve the normalized transition used in the note's illustration."""

    price_long_run = 0.90
    housing_supply_elasticity = 1.0
    fertility_price_semi_elasticity = 1.0
    population_adjustment = 0.02
    housing_depreciation = 0.025

    price_log_ss = np.log(price_long_run)
    fertility_shift = fertility_price_semi_elasticity * price_log_ss
    a = population_adjustment * fertility_price_semi_elasticity
    xi = housing_supply_elasticity
    delta = housing_depreciation

    transition_matrix = np.array(
        [[-a, a], [delta * xi, -delta * (1.0 + xi)]], dtype=float
    )
    forcing = np.array([population_adjustment * fertility_shift, 0.0])
    steady_state = -np.linalg.solve(transition_matrix, forcing)

    years = np.linspace(0.0, 500.0, 501)
    states = np.empty((years.size, 2))
    initial_state = np.zeros(2)
    for index, year in enumerate(years):
        states[index] = steady_state + matrix_exponential_2x2(
            transition_matrix, year
        ) @ (initial_state - steady_state)

    population_log = states[:, 0]
    housing_log = states[:, 1]
    price_log = population_log - housing_log
    fertility_gap = fertility_shift - fertility_price_semi_elasticity * price_log

    roots = np.sort(np.real_if_close(np.linalg.eigvals(transition_matrix)).astype(float))
    series = {
        "population_log": population_log,
        "housing_log": housing_log,
        "price_log": price_log,
        "population_ratio": np.exp(population_log),
        "housing_ratio": np.exp(housing_log),
        "price_ratio": np.exp(price_log),
        "net_reproduction_gap": fertility_gap,
    }
    parameters = {
        "price_long_run": price_long_run,
        "population_long_run": float(np.exp(steady_state[0])),
        "housing_long_run": float(np.exp(steady_state[1])),
        "housing_supply_elasticity": housing_supply_elasticity,
        "population_adjustment": population_adjustment,
        "housing_depreciation": housing_depreciation,
        "slow_root": float(roots[1]),
        "fast_root": float(roots[0]),
        "slow_half_life_years": float(np.log(2.0) / -roots[1]),
    }
    return years, series, parameters


def write_csv(
    years: np.ndarray, series: dict[str, np.ndarray], parameters: dict[str, float]
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with (OUTDIR / "transition_paths.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["year", *series.keys()])
        for index, year in enumerate(years):
            writer.writerow(
                [
                    f"{year:.8f}",
                    *(f"{values[index]:.12g}" for values in series.values()),
                ]
            )

    with (OUTDIR / "illustration_parameters.csv").open("w", newline="") as stream:
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
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def plot_path(
    years: np.ndarray,
    values: np.ndarray,
    long_run: float,
    ylabel: str,
    title: str,
    filename: str,
    ylim: tuple[float, float],
) -> None:
    color = "#1f5a85"
    fig, axis = plt.subplots(figsize=(4.25, 2.70))
    axis.plot(years, values, color=color, linewidth=2.2)
    axis.axhline(1.0, color="#777777", linestyle=(0, (2, 2)), linewidth=1.0)
    axis.axhline(long_run, color="#a33b33", linestyle=(0, (5, 3)), linewidth=1.2)
    axis.set_xlim(years[0], years[-1])
    axis.set_ylim(*ylim)
    axis.set_xlabel("Years after the fertility-schedule shift")
    axis.set_ylabel(ylabel)
    axis.set_title(title, loc="left", fontweight="semibold")
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.6)
    axis.text(
        years[-1] * 0.985,
        long_run + 0.006,
        "new steady state",
        color="#a33b33",
        fontsize=8.5,
        ha="right",
        va="bottom",
    )
    fig.savefig(OUTDIR / filename, format="pdf")
    plt.close(fig)


def main() -> None:
    years, series, parameters = solve_transition()
    write_csv(years, series, parameters)
    apply_plot_style()
    plot_path(
        years,
        series["population_ratio"],
        parameters["population_long_run"],
        r"$N_t/N_0$",
        "Population",
        "population_transition.pdf",
        (0.78, 1.02),
    )
    plot_path(
        years,
        series["price_ratio"],
        parameters["price_long_run"],
        r"$p_t/p_0$",
        "Housing-services price",
        "housing_price_transition.pdf",
        (0.88, 1.02),
    )


if __name__ == "__main__":
    main()
