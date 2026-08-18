#!/usr/bin/env python3
"""Evaluate two non-imposed 2007--2023 demographic closures in the full model.

The first design removes the Census/ACS bridge and otherwise leaves the active
calendar operator unchanged.  The second uses a fixed five-age-group law for
net household formation, migration, and dissolution.  Its five rates are
estimated on the 2007--2019 Census/ACS age-mass transitions, and 2023 remains
an untouched holdout.  At each age, net additions or exits are proportional to
the model's current within-age state distribution.

This is a diagnostic smoke, not a promoted calibration.  Household parameters
and the fertility-preference path are held at the selected transition estimate.
The purpose is to learn whether a no-date-specific-reweighting calibration is
numerically and economically viable before launching an SMM search.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import audit_closed_reproductive_closure as closure
import build_e5f_transition_design_feasibility as design
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as continuation
import run_e5f_transition_calibration as calibration

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REPORT = (
    ROOT
    / "output/model/e5f_transition_ridge_refinement_"
    "jump11_polish_r2_20260818/report/summary.json"
)
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_e5f_child_room_floor_"
    "psinneg_extended_20260806/report/results.json"
)
DEFAULT_AGE_PATH = design.DEFAULT_AGE_PATH
DEFAULT_OUTDIR = ROOT / "output/model/e5f_historical_demographic_closure_smoke"
PERIODS = 5


@dataclass
class SmokeState:
    g_pre: np.ndarray
    scheduled_entries: list[float]
    scheduled_raw_entries: list[float]
    price_guess: float
    initial_policy: Any | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--age-path", type=Path, default=DEFAULT_AGE_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_figures(
    outdir: Path,
    path_rows: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    selected_fit_rows: list[dict[str, Any]],
) -> None:
    colors = {
        "active_no_bridge": "#b34a3c",
        "four_group_net_flow": "#26547c",
        "selected_bridge": "#666666",
    }
    labels = {
        "active_no_bridge": "Current age law, no bridge",
        "four_group_net_flow": "Four-group formation law",
        "selected_bridge": "Selected fit with data bridge",
    }
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2))
    series = (
        ("adult_household_index_2007", "Adult-household mass", "Index (2007 = 1)"),
        ("asset_price_index_2007", "Model house price", "Index (2007 = 1)"),
        ("owner_rate", "Ownership", "Share"),
        ("mean_rooms", "Mean occupied rooms", "Rooms"),
    )
    for ax, (field, title, ylabel) in zip(axes.flat, series, strict=True):
        for design_name in ("active_no_bridge", "four_group_net_flow"):
            rows = [row for row in path_rows if row["design"] == design_name]
            ax.plot(
                [int(row["year"]) for row in rows],
                [float(row[field]) for row in rows],
                marker="o",
                lw=2.0,
                color=colors[design_name],
                label=labels[design_name],
            )
        if field == "adult_household_index_2007":
            targets = transition.census_household_target_indices()
            ax.plot(
                list(design.YEARS),
                [float(targets[index]) for index in range(len(design.YEARS))],
                marker="o",
                lw=2.0,
                color="#222222",
                label="U.S. household total",
            )
        ax.set(title=title, xlabel="Year", ylabel=ylabel)
        ax.grid(alpha=0.18)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "full_model_historical_closure_paths.png", dpi=220)
    fig.savefig(outdir / "full_model_historical_closure_paths.pdf")
    plt.close(fig)

    moments = [row["moment"] for row in selected_fit_rows]
    short = {
        "tfr": "Completed fertility",
        "childless_rate": "Childlessness",
        "mean_age_first_birth": "Mean age at first birth",
        "share_first_births_age30plus": "First births at 30+",
        "housing_increment_0to1": "First-birth rooms",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": "3+ vs 1--2 rooms",
        "own_family_gap": "Family ownership gap",
        "own_rate": "Ownership",
        "aggregate_mean_occupied_rooms_18_85": "Mean rooms",
        "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / earnings",
        "annual_bequest_flow_to_aggregate_wealth": "Bequests / wealth",
        "old_total_wealth_to_annual_income_p90_p50_7684": "Old wealth p90/p50",
    }
    values: dict[str, list[float]] = {
        "selected_bridge": [float(row["loss_contribution"]) for row in selected_fit_rows]
    }
    for name in ("active_no_bridge", "four_group_net_flow"):
        lookup = {
            row["moment"]: float(row["loss_contribution"])
            for row in fit_rows
            if row["design"] == name
        }
        values[name] = [lookup[moment] for moment in moments]
    y = np.arange(len(moments), dtype=float)
    width = 0.25
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    for offset, name in zip((-width, 0.0, width), values, strict=True):
        ax.barh(
            y + offset,
            values[name],
            height=width,
            color=colors[name],
            label=labels[name],
        )
    ax.set(
        yticks=y,
        yticklabels=[short.get(moment, moment) for moment in moments],
        xlabel="Weighted loss contribution",
        title="2023 target fit at the selected household parameters",
    )
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.18)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "terminal_target_loss_contributions.png", dpi=220)
    fig.savefig(outdir / "terminal_target_loss_contributions.pdf")
    plt.close(fig)


def age_masses(distribution: np.ndarray) -> np.ndarray:
    return np.sum(np.asarray(distribution, dtype=float), axis=(0, 1, 2, 4, 5, 6))


def apply_age_mass_target(
    distribution: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Rescale within-age states to a model-generated target age mass."""
    out = np.asarray(distribution, dtype=float).copy()
    target = np.asarray(target, dtype=float)
    before = age_masses(out)
    if target.shape != before.shape or np.any(target < -1e-12):
        raise RuntimeError("The demographic closure generated an invalid age target")
    factors = np.empty_like(target)
    for age_index, (mass, desired) in enumerate(zip(before, target, strict=True)):
        if mass <= 0.0:
            raise RuntimeError(f"Age cell {age_index} has no model mass to rescale")
        factors[age_index] = float(desired) / float(mass)
        out[:, :, :, age_index, :, :, :] *= factors[age_index]
    after = age_masses(out)
    maximum_gap = float(np.max(np.abs(after - target)))
    if maximum_gap > 2e-12:
        raise RuntimeError(f"Age-mass closure failed: max gap={maximum_gap:.3e}")
    return out, {
        "mass_before": float(np.sum(before)),
        "mass_after": float(np.sum(after)),
        "net_flow": float(np.sum(after) - np.sum(before)),
        "maximum_age_target_gap": maximum_gap,
        "minimum_rescaling_factor": float(np.min(factors)),
        "maximum_rescaling_factor": float(np.max(factors)),
    }


def load_contracts(
    report_path: Path, source_path: Path, chain: Any, model: Any
) -> dict[str, Any]:
    raw = json.loads(report_path.read_text(encoding="utf-8"))
    scientific = dict(raw["scientific_contract"])
    transition_path = report_path.parent / "selected_transition_path.csv"
    return continuation.validate_input_contracts(
        report_path=report_path,
        task_summary_path=None,
        case_dir=None,
        case_transition_path=transition_path,
        source_path=source_path,
        expected_report_sha256=continuation.file_sha256(report_path),
        expected_task_summary_sha256=None,
        expected_case_transition_sha256=continuation.file_sha256(transition_path),
        expected_source_sha256=str(scientific["source_sha256"]),
        expected_target_fingerprint=str(scientific["target_fingerprint"]),
        expected_code_bundle_sha256=str(
            scientific["code_fingerprints"]["bundle_sha256"]
        ),
        expected_renewal_contract_sha256=str(raw["renewal_contract_sha256"]),
        expected_scientific_contract_sha256=str(raw["scientific_contract_sha256"]),
        expected_selection_sha256=str(raw["selection_sha256"]),
        chain=chain,
        model=model,
    )


def advance_state(
    *,
    state: SmokeState,
    evaluation: calendar.PeriodEvaluation,
    parameters: Any,
    b_grid: np.ndarray,
    shared: Any,
    outside_flow: float,
    retention: float,
    demographic_profile: np.ndarray | None,
) -> tuple[SmokeState, dict[str, Any]]:
    empty_next, _, deaths, _ = transition.advance_sequential_calendar_distribution(
        evaluation, np.zeros(int(parameters.I)), parameters, b_grid, shared
    )
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        parameters,
    )
    conversion = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
    scheduled_B, next_queue = transition.advance_birth_vintage_queue(
        state.scheduled_entries,
        float(accounting["topcode_adjusted_birth_children"]),
        conversion,
    )
    _, next_raw_queue = transition.advance_birth_vintage_queue(
        state.scheduled_raw_entries,
        float(evaluation.births),
        conversion,
    )
    shares = np.asarray(parameters.entry_shares, dtype=float).reshape(-1)
    shares /= np.sum(shares)
    entrants = float(outside_flow) * shares
    entrants[0] += float(retention) * float(scheduled_B)
    empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        entrants, parameters, b_grid
    )
    raw_next = empty_next
    adjustment = {
        "mass_before": float(np.sum(raw_next)),
        "mass_after": float(np.sum(raw_next)),
        "net_flow": 0.0,
        "maximum_age_target_gap": 0.0,
        "minimum_rescaling_factor": 1.0,
        "maximum_rescaling_factor": 1.0,
    }
    if demographic_profile is None:
        next_pre = raw_next
    else:
        current_age_mass = age_masses(evaluation.g_post_fertility)
        target_age_mass = (
            design.aging_operator() @ current_age_mass
            + np.asarray(demographic_profile, dtype=float)
            * float(np.sum(current_age_mass))
        )
        next_pre, adjustment = apply_age_mass_target(raw_next, target_age_mass)

    health = calendar.distribution_health(
        {
            "pre_fertility": evaluation.g_pre,
            "post_fertility": evaluation.g_post_fertility,
            "current_choice": evaluation.g_current,
            "next_pre_fertility": next_pre,
        }
    )
    if (
        health["min_distribution_mass"] is None
        or int(health["nonfinite_distribution_count"]) != 0
        or float(health["min_distribution_mass"]) < -1e-14
    ):
        raise RuntimeError(f"Distribution-health gate failed: {health}")
    raw_expected = (
        float(np.sum(evaluation.g_post_fertility))
        - float(deaths)
        + float(np.sum(entrants))
    )
    raw_residual = float(np.sum(raw_next)) - raw_expected
    if abs(raw_residual) > 2e-8:
        raise RuntimeError(f"Raw mass accounting failed: {raw_residual:.3e}")
    return (
        SmokeState(
            g_pre=next_pre,
            scheduled_entries=list(next_queue),
            scheduled_raw_entries=list(next_raw_queue),
            price_guess=float(evaluation.policy.price[0]),
            initial_policy=evaluation.policy,
        ),
        {
            **adjustment,
            "raw_mass_accounting_residual": raw_residual,
            "scheduled_mature_entry_flow": float(scheduled_B),
            "entrant_flow_before_demographic_adjustment": float(np.sum(entrants)),
        },
    )


def run_design(
    *,
    label: str,
    prepared: continuation.PreparedModel,
    contracts: dict[str, Any],
    demographic_profile: np.ndarray | None,
    targets: list[np.ndarray],
    market_tol: float,
    market_max_iter: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    renewal = contracts["renewal_old_state"]
    raw_old_B = (
        transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        * prepared.births_old_raw
    )
    state = SmokeState(
        g_pre=prepared.initial_g_pre_2007.copy(),
        scheduled_entries=[prepared.B_old] * 4,
        scheduled_raw_entries=[raw_old_B] * 4,
        price_guess=prepared.old_price,
    )
    parameters = copy.deepcopy(prepared.parameters)
    counter = calendar.SolveCounter()
    path_rows: list[dict[str, Any]] = []
    age_rows: list[dict[str, Any]] = []
    period_records: dict[int, dict[str, Any]] = {}
    pending_branch: dict[str, Any] | None = None
    first_period_accounting: dict[str, np.ndarray] | None = None
    previous_psi: float | None = None
    max_market = 0.0
    max_raw_mass_residual = 0.0
    max_age_target_gap = 0.0

    for period, year in enumerate(design.YEARS):
        dated_psi = transition.preference_shifter_at_date(
            period, old_psi, new_psi, calibration.TRANSITION_PERIODS
        )
        parameters.psi_child = float(dated_psi)
        parameters.parent_dp_waiver = False
        parameters.parent_dp_waiver_phi = 1.0
        parameters.parent_dp_waiver_locations = np.array([], dtype=int)
        parameters.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
        parameters.parent_dp_waiver_birth_state_only = False
        if previous_psi is not None and not math.isclose(
            dated_psi, previous_psi, rel_tol=0.0, abs_tol=1e-14
        ):
            state.initial_policy = None
        previous_psi = float(dated_psi)
        dynamic_state = continuation.DynamicState(
            g_pre=state.g_pre,
            scheduled_entries=state.scheduled_entries,
            scheduled_raw_entries=state.scheduled_raw_entries,
            price_guess=state.price_guess,
            initial_policy=state.initial_policy,
        )
        evaluation, shared, fallback = continuation.evaluate_state(
            dynamic_state,
            parameters,
            prepared.b_grid,
            prepared.supply_rule,
            counter,
            market_tol,
            market_max_iter,
        )
        if float(evaluation.relative_market_residual) > 2e-4:
            raise RuntimeError("Housing-market residual exceeded the smoke tolerance")
        max_market = max(max_market, float(evaluation.relative_market_residual))

        if pending_branch is None:
            housing_response = calibration.first_birth_housing_response(
                evaluation, parameters, prepared.b_grid, shared
            )
            housing_detail = {
                "housing_response": float(housing_response),
                "status": "same-date diagnostic",
            }
        else:
            housing_detail = calibration.finish_dated_first_birth_housing_branch(
                pending_branch,
                evaluation,
                parameters,
                prepared.b_grid,
                shared,
                destination_period=period,
            )
            housing_response = float(housing_detail["housing_response"])
        moment_values = calibration.transition_cross_section_moments(
            evaluation,
            parameters,
            prepared.b_grid,
            shared,
            e5_target_system().moment_names,
            housing_increment_override=float(housing_response),
        )
        period_fertility = calibration.period_fertility_diagnostics(
            evaluation, parameters
        )
        accounting = calibration.first_birth_accounting_by_age(evaluation, parameters)
        if first_period_accounting is None:
            first_period_accounting = accounting
        period_records[period] = {
            "moments": moment_values,
            "first_birth_accounting_by_age": accounting,
            "housing": housing_detail,
        }
        pending_branch = (
            calibration.begin_dated_first_birth_housing_branch(
                evaluation,
                parameters,
                prepared.b_grid,
                shared,
                origin_period=period,
            )
            if period < PERIODS - 1
            else None
        )

        current_age = age_masses(evaluation.g_pre)
        current_metrics = design.age_metrics(current_age, targets[period])
        path_rows.append(
            {
                "design": label,
                "period": period,
                "year": year,
                "psi_child": float(dated_psi),
                "asset_price": float(evaluation.policy.price[0]),
                "asset_price_index_2007": float(evaluation.policy.price[0])
                / prepared.old_price,
                "adult_household_mass": float(np.sum(evaluation.g_pre)),
                "adult_household_index_2007": float(np.sum(evaluation.g_pre))
                / prepared.initial_mass_2007,
                "age_total_variation": current_metrics["age_total_variation"],
                "owner_rate": float(moment_values["own_rate"]),
                "mean_rooms": float(
                    moment_values["aggregate_mean_occupied_rooms_18_85"]
                ),
                "completed_fertility_age42": float(moment_values["tfr"]),
                "childless_rate_age42": float(moment_values["childless_rate"]),
                "period_tfr_explicit_states": float(
                    period_fertility["period_tfr_explicit_states"]
                ),
                "period_tfr_topcode_adjusted": float(
                    period_fertility["period_tfr_topcode_adjusted"]
                ),
                "period_first_birth_mean_age": float(
                    period_fertility["period_first_birth_mean_age"]
                ),
                "period_first_birth_share_age30plus": float(
                    period_fertility["period_first_birth_share_age30plus"]
                ),
                "first_birth_housing_response": float(housing_response),
                "relative_market_residual": float(evaluation.relative_market_residual),
                "grid_resolution_fallback": bool(fallback),
            }
        )
        shares = current_age / np.sum(current_age)
        for age_index, age_lower in enumerate(design.AGE_LOWER):
            age_rows.append(
                {
                    "design": label,
                    "year": year,
                    "age_lower": age_lower,
                    "model_share": float(shares[age_index]),
                    "data_share": float(targets[period][age_index] / np.sum(targets[period])),
                }
            )

        if period < PERIODS - 1:
            state, adjustment = advance_state(
                state=state,
                evaluation=evaluation,
                parameters=parameters,
                b_grid=prepared.b_grid,
                shared=shared,
                outside_flow=float(renewal["outside_flow_M"]),
                retention=float(renewal["retention_rho"]),
                demographic_profile=demographic_profile,
            )
            max_raw_mass_residual = max(
                max_raw_mass_residual,
                abs(float(adjustment["raw_mass_accounting_residual"])),
            )
            max_age_target_gap = max(
                max_age_target_gap,
                abs(float(adjustment["maximum_age_target_gap"])),
            )
            path_rows[-1].update(
                next_demographic_net_flow=float(adjustment["net_flow"]),
                next_minimum_age_rescaling_factor=float(
                    adjustment["minimum_rescaling_factor"]
                ),
                next_maximum_age_rescaling_factor=float(
                    adjustment["maximum_rescaling_factor"]
                ),
            )

    if first_period_accounting is None:
        raise RuntimeError("The first-period fertility accounting is missing")
    timing = calibration.cohort_timing_moments(
        first_period_accounting,
        period_records,
        parameters,
        terminal_period=calibration.TRANSITION_PERIODS,
    )
    terminal = dict(period_records[calibration.TRANSITION_PERIODS]["moments"])
    terminal["mean_age_first_birth"] = float(timing["mean_age_first_birth"])
    terminal["share_first_births_age30plus"] = float(
        timing["share_first_births_age30plus"]
    )
    target_system = e5_target_system()
    fit_rows: list[dict[str, Any]] = []
    total_loss = 0.0
    for name, target, weight in zip(
        target_system.moment_names,
        target_system.target_values,
        target_system.weights,
        strict=True,
    ):
        model_value = float(terminal[name])
        gap = model_value - float(target)
        contribution = float(weight) * gap * gap
        total_loss += contribution
        fit_rows.append(
            {
                "design": label,
                "moment": name,
                "target": float(target),
                "model": model_value,
                "gap": gap,
                "weight": float(weight),
                "loss_contribution": contribution,
            }
        )
    return path_rows, age_rows, {
        "design": label,
        "transition_loss_at_selected_parameters": total_loss,
        "max_market_residual": max_market,
        "max_raw_mass_accounting_residual": max_raw_mass_residual,
        "max_demographic_age_target_gap": max_age_target_gap,
        "model_solve_count": int(counter.total),
        "terminal_timing": timing,
        "fit_rows": fit_rows,
    }


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    report_path = args.selected_report.resolve()
    source_path = args.source.resolve()
    age_path = args.age_path.resolve()
    for path in (report_path, source_path, age_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    started = time.perf_counter()
    chain, model = transition.configure_sequential_model()
    contracts = load_contracts(report_path, source_path, chain, model)
    prepared = continuation.prepare_model(contracts, source_path, chain, model)
    shares = design.load_age_shares(age_path)
    population_indices = transition.census_household_target_indices()
    targets = [
        float(population_indices[index]) * shares[year]
        for index, year in enumerate(design.YEARS)
    ]
    four_profile, profile_rows = design.fit_group_profile(
        targets, design.aging_operator(), design.FOUR_GROUPS
    )

    all_paths: list[dict[str, Any]] = []
    all_ages: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    for label, profile in (
        ("active_no_bridge", None),
        ("four_group_net_flow", four_profile),
    ):
        print(f"HISTORICAL_DEMOGRAPHIC_CLOSURE_START design={label}", flush=True)
        paths, ages, summary = run_design(
            label=label,
            prepared=prepared,
            contracts=contracts,
            demographic_profile=profile,
            targets=targets,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
        )
        all_paths.extend(paths)
        all_ages.extend(ages)
        fit_rows.extend(summary.pop("fit_rows"))
        summaries.append(summary)
        print(
            "HISTORICAL_DEMOGRAPHIC_CLOSURE_DONE "
            f"design={label} loss={summary['transition_loss_at_selected_parameters']:.9g}",
            flush=True,
        )

    write_csv(outdir / "transition_paths.csv", all_paths)
    write_csv(outdir / "age_distribution_fit.csv", all_ages)
    write_csv(outdir / "terminal_target_fit.csv", fit_rows)
    write_csv(outdir / "four_group_net_flow_parameters.csv", profile_rows)
    raw_report = json.loads(report_path.read_text(encoding="utf-8"))
    make_figures(outdir, all_paths, fit_rows, list(raw_report["target_fit_table"]))
    payload = {
        "status": "complete_historical_demographic_closure_smoke",
        "scope": "selected household parameters and preference path; no refit",
        "selected_bridge_calibration_loss": float(raw_report["selected_transition_loss"]),
        "selected_report": str(report_path),
        "selected_report_sha256": continuation.file_sha256(report_path),
        "source": str(source_path),
        "source_sha256": continuation.file_sha256(source_path),
        "age_path": str(age_path),
        "age_path_sha256": continuation.file_sha256(age_path),
        "demographic_closure": {
            "training_dates": [2007, 2011, 2015, 2019],
            "holdout_date": 2023,
            "net_flow_interpretation": (
                "within-age proportional net household formation, migration, and dissolution"
            ),
            "date_specific_observed_reweighting": False,
        },
        "designs": summaries,
        "elapsed_seconds": time.perf_counter() - started,
    }
    continuation.write_json(outdir / "summary.json", payload)
    (outdir / "README.md").write_text(
        "# Full-model historical demographic-closure smoke\n\n"
        "This packet holds the selected household parameters and fertility-preference "
        "path fixed. It compares the raw calendar law with a five-parameter, "
        "four-parameter, time-invariant net household-formation/migration law "
        "estimated through 2019. "
        "No observed 2011--2023 population total or age share is inserted into either "
        "simulation. The full terminal target-fit table is in `terminal_target_fit.csv`.\n",
        encoding="utf-8",
    )
    print(
        "HISTORICAL_DEMOGRAPHIC_CLOSURE_SMOKE_COMPLETE "
        f"elapsed={payload['elapsed_seconds']:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
