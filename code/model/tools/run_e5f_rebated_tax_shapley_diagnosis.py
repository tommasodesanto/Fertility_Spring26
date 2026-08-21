#!/usr/bin/env python3
"""Shapley decomposition of the 2023 rebated-property-tax response."""

from __future__ import annotations

import argparse
import copy
import csv
import itertools
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_coven_property_tax_smoke as impact
import run_e5f_post2023_no_policy_continuations as baseline
import run_e5f_post2023_policy_mechanisms as mechanism
import run_e5f_post2023_rebated_property_tax_smoke as tax


FACTORS = ("tax_rate", "asset_price", "equal_rebate")
METRICS = (
    "births_per_adult",
    "owner_rate",
    "dependent_child_owner_rate",
    "housing_demand_per_adult",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, required=True)
    parser.add_argument("--selected-case-transition", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--rebated-path-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--expected-case-transition-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256", required=True)
    parser.add_argument("--expected-selection-sha256", required=True)
    return parser.parse_args()


def read_policy_endpoints(path: Path) -> dict[str, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    result: dict[str, dict[str, float]] = {}
    for case in ("rebated-tax1-baseline", "rebated-tax2-reform"):
        found = [
            row
            for row in rows
            if row["policy_case"] == case and int(row["calendar_year"]) == 2023
        ]
        if len(found) != 1:
            raise RuntimeError(f"Expected one 2023 row for {case}; found {len(found)}")
        row = found[0]
        result[case] = {
            "asset_price": float(row["asset_price"]),
            "equal_rebate": float(row["equal_transfer_period_units"]),
            "annual_tax_rate": float(row["annual_property_tax_rate"]),
            "births_per_adult": float(row["topcode_adjusted_births_per_adult"]),
            "owner_rate": float(row["owner_rate"]),
            "dependent_child_owner_rate": float(row["dependent_child_owner_rate"]),
            "housing_demand_per_adult": float(row["housing_demand_per_adult"]),
        }
    return result


def shapley_decomposition(
    cell_values: dict[tuple[int, int, int], float]
) -> dict[str, float]:
    if set(cell_values) != set(itertools.product((0, 1), repeat=3)):
        raise ValueError("The Shapley decomposition requires all eight cells")
    contributions = {factor: 0.0 for factor in FACTORS}
    for order in itertools.permutations(range(3)):
        state = [0, 0, 0]
        previous = float(cell_values[tuple(state)])
        for index in order:
            state[index] = 1
            current = float(cell_values[tuple(state)])
            contributions[FACTORS[index]] += (current - previous) / 6.0
            previous = current
    total = sum(contributions.values())
    exact = float(cell_values[(1, 1, 1)] - cell_values[(0, 0, 0)])
    if not math.isclose(total, exact, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"Shapley add-up failed: {total} versus {exact}")
    return contributions


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    chain, model = transition.configure_sequential_model()
    contracts = baseline.validate_input_contracts(
        report_path=args.selected_report,
        task_summary_path=None,
        case_dir=None,
        case_transition_path=args.selected_case_transition,
        source_path=args.source,
        expected_report_sha256=args.expected_report_sha256,
        expected_task_summary_sha256=None,
        expected_case_transition_sha256=args.expected_case_transition_sha256,
        expected_source_sha256=args.expected_source_sha256,
        expected_target_fingerprint=args.expected_target_fingerprint,
        expected_code_bundle_sha256=args.expected_code_bundle_sha256,
        expected_renewal_contract_sha256=args.expected_renewal_contract_sha256,
        expected_scientific_contract_sha256=args.expected_scientific_contract_sha256,
        expected_selection_sha256=args.expected_selection_sha256,
        chain=chain,
        model=model,
    )
    prepared = baseline.prepare_model(contracts, args.source.resolve(), chain, model)
    history, state_2023, parameters_2023, mass_2023, reconstruction = (
        mechanism.reconstruct_matched_2023(
            prepared,
            contracts,
            market_tol=2.0e-4,
            market_max_iter=30,
            progress_dir=outdir,
        )
    )
    endpoints = read_policy_endpoints(args.rebated_path_csv)
    base = endpoints["rebated-tax1-baseline"]
    reform = endpoints["rebated-tax2-reform"]
    levels = {
        "tax_rate": (base["annual_tax_rate"], reform["annual_tax_rate"]),
        "asset_price": (base["asset_price"], reform["asset_price"]),
        "equal_rebate": (base["equal_rebate"], reform["equal_rebate"]),
    }

    counter = calendar.SolveCounter()
    cell_rows: list[dict[str, Any]] = []
    metric_cells = {metric: {} for metric in METRICS}
    for cell in itertools.product((0, 1), repeat=3):
        annual_tax = levels["tax_rate"][cell[0]]
        price = levels["asset_price"][cell[1]]
        transfer = levels["equal_rebate"][cell[2]]
        case = tax.TaxCase(
            name=f"tax{annual_tax:.2f}",
            label="fixed-component diagnostic",
            annual_tax_rate=annual_tax,
            rebate_revenue=True,
        )
        P = copy.deepcopy(parameters_2023)
        solved = impact.evaluate_fixed_price(
            state=state_2023,
            P=P,
            b_grid=prepared.b_grid,
            price=price,
            transfer=transfer,
            case=case,
            counter=counter,
        )
        evaluation = solved.evaluation
        adult_mass = float(np.sum(evaluation.g_post_fertility))
        current_mass = float(np.sum(evaluation.g_current))
        dependent_mass = float(np.sum(evaluation.g_current[..., 1:]))
        owner_mass = float(np.sum(evaluation.g_current[:, 1:, ...]))
        dependent_owner_mass = float(
            np.sum(evaluation.g_current[:, 1:, ..., 1:])
        )
        adjusted_births = transition.calendar_topcode_birth_accounting(
            evaluation.g_pre, evaluation.g_post_fertility, evaluation.births, P
        )["topcode_adjusted_birth_children"]
        values = {
            "births_per_adult": adjusted_births / adult_mass,
            "owner_rate": owner_mass / current_mass,
            "dependent_child_owner_rate": dependent_owner_mass / dependent_mass,
            "housing_demand_per_adult": float(np.sum(evaluation.demand_by_loc))
            / adult_mass,
        }
        for metric, value in values.items():
            metric_cells[metric][cell] = value
        cell_rows.append(
            {
                "tax_reform": cell[0],
                "price_reform": cell[1],
                "rebate_reform": cell[2],
                "annual_tax_rate": annual_tax,
                "asset_price": price,
                "equal_rebate": transfer,
                **values,
            }
        )
        print(f"SHAPLEY_CELL cell={cell} births={values['births_per_adult']:.8f}", flush=True)

    for cell, endpoint in (((0, 0, 0), base), ((1, 1, 1), reform)):
        row = next(
            item
            for item in cell_rows
            if tuple(item[key] for key in ("tax_reform", "price_reform", "rebate_reform"))
            == cell
        )
        for metric in METRICS:
            if not math.isclose(
                float(row[metric]), float(endpoint[metric]), rel_tol=0.0, abs_tol=5e-10
            ):
                raise RuntimeError(
                    f"Endpoint reproduction failed for {cell} {metric}: "
                    f"{row[metric]} versus {endpoint[metric]}"
                )

    shapley_rows: list[dict[str, Any]] = []
    for metric in METRICS:
        contributions = shapley_decomposition(metric_cells[metric])
        base_value = metric_cells[metric][(0, 0, 0)]
        for factor in FACTORS:
            contribution = contributions[factor]
            scale = 100.0 if "rate" in metric else 100.0 / base_value
            unit = "percentage_points" if "rate" in metric else "percent_of_baseline"
            shapley_rows.append(
                {
                    "metric": metric,
                    "factor": factor,
                    "level_contribution": contribution,
                    "reported_contribution": scale * contribution,
                    "reported_unit": unit,
                }
            )
    write_csv(outdir / "component_cells.csv", cell_rows)
    write_csv(outdir / "shapley_decomposition.csv", shapley_rows)
    baseline.write_csv(outdir / "shared_2007_2019_history.csv", history)
    baseline.write_csv(
        outdir / "history_reproduction_audit.csv",
        reconstruction.pop("history_reproduction_rows"),
    )

    plot_metrics = (
        ("births_per_adult", "Births per adult household"),
        ("dependent_child_owner_rate", "Dependent-child ownership"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.5), constrained_layout=True)
    colors = {"tax_rate": "#756bb1", "asset_price": "#2c7fb8", "equal_rebate": "#b0493d"}
    for axis, (metric, title) in zip(axes, plot_metrics, strict=True):
        rows = [row for row in shapley_rows if row["metric"] == metric]
        axis.barh(
            ["Tax rate", "House-price capitalization", "Equal rebate"],
            [row["reported_contribution"] for row in rows],
            color=[colors[row["factor"]] for row in rows],
        )
        axis.axvline(0.0, color="#777777", linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel(
            "Contribution (percentage points)"
            if "rate" in metric
            else "Contribution (% of baseline)"
        )
        axis.grid(axis="x", alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "rebated_tax_shapley.pdf", bbox_inches="tight")
    fig.savefig(outdir / "rebated_tax_shapley.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "status": "complete",
        "scope": (
            "2023 partial-equilibrium Shapley decomposition holding the fitted "
            "pre-fertility distribution fixed; all eight tax-price-rebate cells solved"
        ),
        "factors": list(FACTORS),
        "model_solve_count": counter.total,
        "endpoint_reproduction_tolerance": 5e-10,
        "reconstruction": reconstruction,
        "input_hashes": contracts["hashes"],
        "levels": levels,
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"REBATED_TAX_SHAPLEY_COMPLETE solves={counter.total} output={outdir}")


if __name__ == "__main__":
    main()
