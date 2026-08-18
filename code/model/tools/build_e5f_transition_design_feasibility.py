#!/usr/bin/env python3
"""Compare a no-bridge historical age law with a 2023-forward benchmark.

This is a demographic-accounting feasibility packet, not a household-model
recalibration.  It asks whether the active age law can reproduce the observed
2007--2023 household-head distribution without the date-by-date Census/ACS
reweighting.  It then estimates a time-invariant, five-age-group net household
formation/migration profile on the 2007--2019 transitions and evaluates 2023 as
an untouched holdout.  The packet also records the already-calibrated 2023
state and its no-policy forward continuation so the two proposed quantitative
designs can be compared without changing the core household solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from intergen_eqscale_seq_optimized.e1_profile import _survival_schedule
from run_e5f_open_population_transition import census_household_target_indices


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_AGE_PATH = (
    ROOT
    / "code/data/Spatial_aggregate_withmicrodata/output/"
    "national_householder_age_path/national_householder_age_path.csv"
)
DEFAULT_SELECTED_REPORT = (
    ROOT
    / "output/model/e5f_transition_ridge_refinement_"
    "jump11_polish_r2_20260818/report/summary.json"
)
DEFAULT_CONTINUATION_DIR = (
    ROOT
    / "output/model/e5f_no_policy_transition_report_"
    "jump11_polish_r2_20260818"
)
DEFAULT_OUTDIR = ROOT / "output/model/e5f_transition_design_feasibility"

YEARS = (2007, 2011, 2015, 2019, 2023)
AGE_LOWER = tuple(range(18, 83, 4))
FOUR_GROUPS: tuple[tuple[str, tuple[int, ...]], ...] = (
    ("18--21", (0,)),
    ("22--29", (1, 2)),
    ("30--45", (3, 4, 5, 6)),
    ("46--85", (7, 8, 9, 10, 11, 12, 13, 14, 15, 16)),
)
FIVE_GROUPS: tuple[tuple[str, tuple[int, ...]], ...] = (
    ("18--21", (0,)),
    ("22--29", (1, 2)),
    ("30--45", (3, 4, 5, 6)),
    ("46--65", (7, 8, 9, 10, 11)),
    ("66--85", (12, 13, 14, 15, 16)),
)
SEVEN_GROUPS: tuple[tuple[str, tuple[int, ...]], ...] = (
    ("18--21", (0,)),
    ("22--25", (1,)),
    ("26--29", (2,)),
    ("30--37", (3, 4)),
    ("38--49", (5, 6, 7)),
    ("50--65", (8, 9, 10, 11)),
    ("66--85", (12, 13, 14, 15, 16)),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--age-path", type=Path, default=DEFAULT_AGE_PATH)
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_SELECTED_REPORT)
    parser.add_argument("--continuation-dir", type=Path, default=DEFAULT_CONTINUATION_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def load_age_shares(path: Path) -> dict[int, np.ndarray]:
    rows = read_csv(path)
    by_year: dict[int, list[tuple[int, float]]] = {}
    for row in rows:
        year = int(row["year"])
        by_year.setdefault(year, []).append(
            (int(row["age_lower"]), float(row["share_conditional_age_range"]))
        )
    if tuple(sorted(by_year)) != YEARS:
        raise RuntimeError(f"Age-path years differ from {YEARS}: {sorted(by_year)}")
    result: dict[int, np.ndarray] = {}
    for year in YEARS:
        ordered = sorted(by_year[year])
        if tuple(age for age, _ in ordered) != AGE_LOWER:
            raise RuntimeError(f"Age cells are invalid in {year}")
        values = np.asarray([share for _, share in ordered], dtype=float)
        if not math.isclose(float(np.sum(values)), 1.0, rel_tol=0.0, abs_tol=2e-10):
            raise RuntimeError(f"Age shares do not sum to one in {year}")
        result[year] = values
    return result


def aging_operator() -> np.ndarray:
    survival = np.asarray(_survival_schedule(), dtype=float)
    if survival.shape != (16,) or np.any((survival < 0.0) | (survival > 1.0)):
        raise RuntimeError("The active survival schedule is invalid")
    operator = np.zeros((17, 17), dtype=float)
    operator[1:, :-1] = np.diag(survival)
    return operator


def age_metrics(model_mass: np.ndarray, target_mass: np.ndarray) -> dict[str, float]:
    model_total = float(np.sum(model_mass))
    target_total = float(np.sum(target_mass))
    model_share = model_mass / model_total
    target_share = target_mass / target_total
    return {
        "model_total": model_total,
        "target_total": target_total,
        "total_gap": model_total - target_total,
        "age_total_variation": 0.5 * float(np.sum(np.abs(model_share - target_share))),
        "maximum_age_share_gap": float(np.max(np.abs(model_share - target_share))),
        "age_mass_rmse": float(np.sqrt(np.mean((model_mass - target_mass) ** 2))),
    }


def simulate_no_bridge(
    initial: np.ndarray,
    operator: np.ndarray,
    inherited_entry_flow: float,
) -> list[np.ndarray]:
    states = [np.asarray(initial, dtype=float).copy()]
    for _ in range(len(YEARS) - 1):
        next_mass = operator @ states[-1]
        next_mass[0] = float(inherited_entry_flow)
        states.append(next_mass)
    return states


def fit_group_profile(
    targets: list[np.ndarray],
    operator: np.ndarray,
    groups: tuple[tuple[str, tuple[int, ...]], ...],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    # Estimate rates only on 2007->2011, 2011->2015, and 2015->2019.
    residual_rates = np.stack(
        [
            (targets[k + 1] - operator @ targets[k]) / float(np.sum(targets[k]))
            for k in range(3)
        ]
    )
    profile = np.zeros(17, dtype=float)
    rows: list[dict[str, Any]] = []
    for label, indices in groups:
        estimate = float(np.mean(residual_rates[:, list(indices)]))
        profile[list(indices)] = estimate
        rows.append(
            {
                "age_group": label,
                "age_cell_count": len(indices),
                "net_rate_per_age_cell_per_four_years": estimate,
                "net_group_rate_per_four_years": estimate * len(indices),
                "training_transitions": "2007-2011;2011-2015;2015-2019",
                "interpretation": "net household formation, migration, and dissolution",
            }
        )
    return profile, rows


def simulate_profile(
    initial: np.ndarray, operator: np.ndarray, profile: np.ndarray
) -> list[np.ndarray]:
    states = [np.asarray(initial, dtype=float).copy()]
    for _ in range(len(YEARS) - 1):
        current = states[-1]
        next_mass = operator @ current + profile * float(np.sum(current))
        if np.any(next_mass < -1e-12):
            raise RuntimeError("The fitted net-flow profile generates negative age mass")
        states.append(np.maximum(next_mass, 0.0))
    return states


def stable_growth_diagnostic(operator: np.ndarray, profile: np.ndarray) -> dict[str, Any]:
    law = operator + np.outer(profile, np.ones(17, dtype=float))
    values, vectors = np.linalg.eig(law)
    index = int(np.argmax(values.real))
    growth = float(values[index].real)
    vector = np.asarray(vectors[:, index].real, dtype=float)
    if float(np.sum(vector)) < 0.0:
        vector *= -1.0
    if np.any(vector < -1e-10) or growth <= 0.0:
        return {"usable": False, "dominant_eigenvalue": growth}
    shares = np.maximum(vector, 0.0)
    shares /= np.sum(shares)
    return {
        "usable": True,
        "dominant_eigenvalue_four_years": growth,
        "annual_growth_rate": growth ** 0.25 - 1.0,
        "stable_age_shares": shares,
    }


def selected_2023_and_forward(
    selected_report: dict[str, Any], continuation_dir: Path
) -> dict[str, Any]:
    best = dict(selected_report.get("best_candidate") or {})
    if selected_report.get("status") != "complete_refinement_with_two_independent_identity_repeats":
        raise RuntimeError("Selected transition report is not complete")
    continuation_summary_path = continuation_dir / "continuation_summary.json"
    path_file = continuation_dir / "continuation_paired_continuation_path.csv"
    if not continuation_summary_path.is_file() or not path_file.is_file():
        raise FileNotFoundError("The completed 2023-forward continuation packet is missing")
    continuation = json.loads(continuation_summary_path.read_text(encoding="utf-8"))
    if continuation.get("status") != "complete_no_policy_post2023_continuation_pair":
        raise RuntimeError("The 2023-forward continuation is not complete")
    rows = read_csv(path_file)
    by_scenario: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_scenario.setdefault(row["scenario"], []).append(row)
    result: dict[str, Any] = {
        "selected_transition_loss": float(best["transition_loss"]),
        "selected_2023_completed_fertility": float(best["terminal_tfr"]),
        "selected_2023_population_index_from_2007": float(best["terminal_population_index"]),
        "selected_2023_asset_price_index_from_2007": float(best["terminal_housing_cost_index"]),
        "selected_parameter_count": len(
            [
                name
                for name in best["theta"]
                if name
                not in {
                    "alpha_cons",
                    "delta_alpha",
                    "delta_alpha_jump",
                    "psi_child",
                    "tenure_choice_kappa",
                    "theta_n",
                }
            ]
        )
        + 1,
        "continuation_status": continuation["status"],
        "between_steady_states_statement": continuation["between_steady_states_statement"],
        "scenarios": {},
    }
    for scenario, scenario_rows in by_scenario.items():
        ordered = sorted(scenario_rows, key=lambda row: int(float(row["calendar_year"])))
        post = [row for row in ordered if int(float(row["calendar_year"])) >= 2023]
        if not post:
            continue
        first, last = post[0], post[-1]
        result["scenarios"][scenario] = {
            "initial_year": int(float(first["calendar_year"])),
            "terminal_year": int(float(last["calendar_year"])),
            "initial_population_index_2023": float(first["population_index_2023"]),
            "terminal_population_index_2023": float(last["population_index_2023"]),
            "initial_asset_price_index": float(first["asset_price_index"]),
            "terminal_asset_price_index": float(last["asset_price_index"]),
        }
    return result


def make_figures(
    comparison_rows: list[dict[str, Any]],
    age_rows: list[dict[str, Any]],
    outdir: Path,
) -> None:
    colors = {
        "data": "#222222",
        "active_no_bridge": "#b34a3c",
        "four_group_profile": "#26547c",
        "five_group_profile": "#26547c",
        "seven_group_profile": "#6a4c93",
    }
    labels = {
        "data": "U.S. data",
        "active_no_bridge": "Current age law, no bridge",
        "four_group_profile": "Four-group formation law",
        "five_group_profile": "Five-group sensitivity",
        "seven_group_profile": "Seven-group sensitivity",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1))
    for design in labels:
        rows = [row for row in comparison_rows if row["design"] == design]
        axes[0].plot(
            [row["year"] for row in rows],
            [row["model_total"] for row in rows],
            marker="o",
            lw=2.0,
            color=colors[design],
            label=labels[design],
        )
    axes[0].set(title="Adult-household mass", ylabel="Index (2007 = 1)", xlabel="Year")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.18)

    for design in ("active_no_bridge", "four_group_profile", "seven_group_profile"):
        rows = [row for row in comparison_rows if row["design"] == design]
        axes[1].plot(
            [row["year"] for row in rows],
            [row["age_total_variation"] for row in rows],
            marker="o",
            lw=2.0,
            color=colors[design],
            label=labels[design],
        )
    axes[1].set(
        title="Age-distribution error",
        ylabel="Total-variation distance",
        xlabel="Year",
    )
    axes[1].grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(outdir / "historical_demographic_feasibility.png", dpi=220)
    fig.savefig(outdir / "historical_demographic_feasibility.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    for design in ("data", "active_no_bridge", "four_group_profile", "seven_group_profile"):
        rows = [
            row for row in age_rows if row["design"] == design and row["year"] == 2023
        ]
        ax.plot(
            [row["age_lower"] for row in rows],
            [row["age_share"] for row in rows],
            marker="o",
            lw=2.0,
            color=colors[design],
            label=labels[design],
        )
    ax.set(
        title="2023 household-head age distribution",
        xlabel="Age-cell lower bound",
        ylabel="Share of adult households",
    )
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "historical_2023_age_holdout.png", dpi=220)
    fig.savefig(outdir / "historical_2023_age_holdout.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    age_path = args.age_path.resolve()
    selected_report_path = args.selected_report.resolve()
    continuation_dir = args.continuation_dir.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    for path in (age_path, selected_report_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    shares = load_age_shares(age_path)
    indices = census_household_target_indices()
    if tuple(sorted(indices)) != tuple(range(len(YEARS))):
        raise RuntimeError("The Census household index contract is incomplete")
    targets = [float(indices[k]) * shares[year] for k, year in enumerate(YEARS)]
    operator = aging_operator()

    selected_report = json.loads(selected_report_path.read_text(encoding="utf-8"))
    best = dict(selected_report.get("best_candidate") or {})
    old_state = dict(best.get("renewal_accounting_old_state") or {})
    inherited_entry_flow = float(old_state["old_entry_flow_E"])
    no_bridge = simulate_no_bridge(targets[0], operator, inherited_entry_flow)
    four_profile, four_parameters = fit_group_profile(targets, operator, FOUR_GROUPS)
    five_profile, five_parameters = fit_group_profile(targets, operator, FIVE_GROUPS)
    seven_profile, seven_parameters = fit_group_profile(targets, operator, SEVEN_GROUPS)
    four_path = simulate_profile(targets[0], operator, four_profile)
    five_path = simulate_profile(targets[0], operator, five_profile)
    seven_path = simulate_profile(targets[0], operator, seven_profile)
    designs = {
        "data": targets,
        "active_no_bridge": no_bridge,
        "four_group_profile": four_path,
        "five_group_profile": five_path,
        "seven_group_profile": seven_path,
    }

    comparison_rows: list[dict[str, Any]] = []
    age_rows: list[dict[str, Any]] = []
    for design, states in designs.items():
        for k, (year, state) in enumerate(zip(YEARS, states, strict=True)):
            metrics = age_metrics(state, targets[k])
            comparison_rows.append({"design": design, "year": year, **metrics})
            state_share = state / np.sum(state)
            for j, age in enumerate(AGE_LOWER):
                age_rows.append(
                    {
                        "design": design,
                        "year": year,
                        "age_lower": age,
                        "age_upper": age + 3,
                        "age_mass": float(state[j]),
                        "age_share": float(state_share[j]),
                    }
                )

    write_csv(outdir / "historical_path_comparison.csv", comparison_rows)
    write_csv(outdir / "historical_age_distributions.csv", age_rows)
    write_csv(outdir / "four_group_net_flow_parameters.csv", four_parameters)
    write_csv(outdir / "five_group_net_flow_parameters.csv", five_parameters)
    write_csv(outdir / "seven_group_net_flow_parameters.csv", seven_parameters)
    make_figures(comparison_rows, age_rows, outdir)

    four_stable = stable_growth_diagnostic(operator, four_profile)
    five_stable = stable_growth_diagnostic(operator, five_profile)
    seven_stable = stable_growth_diagnostic(operator, seven_profile)
    if four_stable.get("usable"):
        four_stable["distance_2007_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2007] - np.asarray(four_stable["stable_age_shares"])))
        )
        four_stable["distance_2023_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2023] - np.asarray(four_stable["stable_age_shares"])))
        )
    if five_stable.get("usable"):
        five_stable["distance_2007_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2007] - np.asarray(five_stable["stable_age_shares"])))
        )
        five_stable["distance_2023_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2023] - np.asarray(five_stable["stable_age_shares"])))
        )
    if seven_stable.get("usable"):
        seven_stable["distance_2007_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2007] - np.asarray(seven_stable["stable_age_shares"])))
        )
        seven_stable["distance_2023_to_stable_age_distribution"] = 0.5 * float(
            np.sum(np.abs(shares[2023] - np.asarray(seven_stable["stable_age_shares"])))
        )

    def endpoint(design: str) -> dict[str, Any]:
        return next(
            row
            for row in comparison_rows
            if row["design"] == design and int(row["year"]) == 2023
        )

    today_forward = selected_2023_and_forward(selected_report, continuation_dir)
    summary = {
        "status": "complete_transition_design_feasibility",
        "scope": (
            "demographic accounting and existing 2023-forward continuation; "
            "no household-model recalibration"
        ),
        "inputs": {
            "age_path": str(age_path),
            "age_path_sha256": file_sha256(age_path),
            "selected_report": str(selected_report_path),
            "selected_report_sha256": file_sha256(selected_report_path),
            "continuation_dir": str(continuation_dir),
        },
        "historical_no_bridge": {
            "inherited_entry_flow_per_2007_adult": inherited_entry_flow,
            "2023": endpoint("active_no_bridge"),
            "interpretation": (
                "The active law injects households only in the youngest age cell. "
                "Household preference parameters cannot repair this age-mass path "
                "before post-2007 births reach entry age."
            ),
        },
        "four_group_profile": {
            "estimated_parameter_count": len(FOUR_GROUPS),
            "training_dates": [2007, 2011, 2015, 2019],
            "holdout_date": 2023,
            "2023": endpoint("four_group_profile"),
            "stable_growth_diagnostic": four_stable,
        },
        "five_group_sensitivity": {
            "estimated_parameter_count": len(FIVE_GROUPS),
            "training_dates": [2007, 2011, 2015, 2019],
            "holdout_date": 2023,
            "2023": endpoint("five_group_profile"),
            "stable_growth_diagnostic": five_stable,
        },
        "seven_group_sensitivity": {
            "estimated_parameter_count": len(SEVEN_GROUPS),
            "training_dates": [2007, 2011, 2015, 2019],
            "holdout_date": 2023,
            "2023": endpoint("seven_group_profile"),
            "stable_growth_diagnostic": seven_stable,
        },
        "today_forward": today_forward,
        "next_gate": (
            "Embed the four-group proportional net-flow law in the full distribution "
            "operator, preserving within-age states, and run one exact no-reweighting "
            "GE smoke before any joint calibration."
        ),
    }
    write_json(outdir / "summary.json", summary)
    readme = f"""# Transition-design feasibility

This packet compares two quantitative designs without changing the household
solver.

The historical diagnostic starts from the observed 2007 household-head age
distribution and removes every later Census/ACS reweight.  Under the active age
law, the 2023 household index is `{endpoint('active_no_bridge')['model_total']:.6f}`
against `{endpoint('active_no_bridge')['target_total']:.6f}` in the data, and the
age-distribution total-variation distance is
`{endpoint('active_no_bridge')['age_total_variation']:.6f}`.

The four-group law estimates fixed proportional net household
formation/migration/dissolution rates from the first three transitions and
reserves 2023 as a holdout.  Its 2023 household index is
`{endpoint('four_group_profile')['model_total']:.6f}` and its age-distribution
distance is `{endpoint('four_group_profile')['age_total_variation']:.6f}`.  These
four quantities are demographic closure parameters, not household-preference
parameters.

The 2023-forward side records the completed calibrated 2023 state and the
existing closed/open no-policy continuations.  It does not claim that 2023 is a
steady state.

The next numerical gate is one full-model smoke with the four-group law and no
date-specific reweighting.  No joint SMM search should begin before that smoke
passes its mass, market, and holdout-age checks.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print(
        "TRANSITION_DESIGN_FEASIBILITY_COMPLETE "
        f"no_bridge_2023={endpoint('active_no_bridge')['model_total']:.6f} "
        f"four_group_2023={endpoint('four_group_profile')['model_total']:.6f} "
        f"target_2023={endpoint('four_group_profile')['target_total']:.6f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
