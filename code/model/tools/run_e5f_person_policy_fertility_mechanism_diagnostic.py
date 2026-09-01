#!/usr/bin/env python3
"""Expose exact fertility-state and realized-housing diagnostics on one H128 path.

This default-off wrapper runs one immutable evaluation of the production
person-demography perfect-foresight policy path.  At selected dates it records
the full at-risk fertility state grid and compact realized cross-section
aggregates.  The production Bellman problem, transition law, and path update
are not changed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import run_dynamic_population_transition as calendar  # noqa: E402
import run_e5f_open_population_transition as transition  # noqa: E402
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402
import run_e5f_perfect_foresight_person_demography_policy as policy_driver  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402


MASS_FLOOR = 1.0e-15
DEFAULT_YEARS = "2023,2035,2051,2079,2103,2355,2531"


@dataclass
class CaptureContext:
    output_dir: Path
    target_years: set[int]
    horizon: int
    b_grid: np.ndarray
    active: bool = False
    evaluation_calls: int = 0
    current_period: int = 0
    captured_years: set[int] = field(default_factory=set)
    validations: list[dict[str, Any]] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        required=True,
        choices=("rebated-tax1-baseline", "rebated-tax2-reform"),
    )
    parser.add_argument("--terminal-summary", type=Path, required=True)
    parser.add_argument("--initial-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--horizon", type=int, default=128)
    parser.add_argument("--diagnostic-years", default=DEFAULT_YEARS)
    parser.add_argument("--post2023-tenure-choice-kappa", type=float, required=True)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR)
    return parser.parse_args()


def parse_years(text: str, horizon: int) -> set[int]:
    years = {int(item.strip()) for item in str(text).split(",") if item.strip()}
    allowed = {pf.CALENDAR_START_YEAR + 4 * period for period in range(int(horizon))}
    if not years or not years.issubset(allowed):
        raise ValueError(
            f"Diagnostic years must lie on the H{horizon} four-year grid: {sorted(years)}"
        )
    return years


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(pf.jsonable(payload), indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty diagnostic table: {path}")
    handle = (
        gzip.open(path, "wt", newline="", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open("w", newline="", encoding="utf-8")
    )
    with handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def age_group(age_years: float) -> str:
    if age_years < 25.0:
        return "under_25"
    if age_years < 35.0:
        return "25_34"
    if age_years < 45.0:
        return "35_44"
    if age_years < 65.0:
        return "45_64"
    return "65_plus"


def parent_status(parity: int, child_state: int) -> str:
    if parity == 0:
        return "childless"
    if child_state > 0:
        return "dependent_children"
    return "empty_nest_parent"


def binary_value_gap(probability: float, scale: float) -> float:
    probability = float(probability)
    if probability <= 0.0:
        return -math.inf
    if probability >= 1.0:
        return math.inf
    return float(scale) * math.log(probability / (1.0 - probability))


def fertility_attempt(
    policy: calendar.PolicyBundle,
    parameters: SimpleNamespace,
    *,
    wealth_index: int,
    origin_tenure: int,
    market: int,
    age_index: int,
    income_state: int,
    parity: int,
    child_state: int,
) -> tuple[float, float]:
    """Return attempt probability and its attempt-minus-wait value gap."""

    if not (int(parameters.A_f_start) <= age_index + 1 <= int(parameters.A_f_end)):
        return 0.0, math.nan
    if parity == 0:
        if child_state != int(calendar.model.readiness_settled_state(parameters)):
            return 0.0, math.nan
        probability = float(
            policy.fert_probs[
                wealth_index,
                origin_tenure,
                market,
                age_index,
                income_state,
                1,
            ]
        )
        scale = float(parameters.kappa_fert)
    elif 1 <= parity < int(parameters.n_parity) - 1 and child_state <= parity:
        continuation = getattr(parameters, "_fert2_probs", None)
        if continuation is None:
            raise RuntimeError("Continuation-birth probabilities were not retained")
        probability = float(
            continuation[
                wealth_index,
                origin_tenure,
                market,
                age_index,
                income_state,
                1,
                parity - 1,
                child_state,
            ]
        )
        raw_scale = getattr(parameters, "kappa_fert_continuation", None)
        scale = float(parameters.kappa_fert if raw_scale is None else raw_scale)
    else:
        return 0.0, math.nan
    return probability, binary_value_gap(probability, scale)


def fertility_state_rows(
    evaluation: calendar.PeriodEvaluation,
    parameters: SimpleNamespace,
    b_grid: np.ndarray,
    *,
    calendar_year: int,
    period: int,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    g_pre = np.asarray(evaluation.g_pre, dtype=float)
    fecundity = np.asarray(calendar.model.get_fecundity_by_age(parameters), dtype=float)
    z_grid, _, _ = calendar.model.income_transition_values(parameters)
    settled = int(calendar.model.readiness_settled_state(parameters))
    top_parity = int(parameters.n_parity) - 1
    top_weight = float(getattr(parameters, "tfr_top_bin_weight", top_parity))
    rows: list[dict[str, Any]] = []
    adjusted_births = 0.0
    explicit_births = 0.0
    at_risk_mass = 0.0
    for age in range(int(parameters.J)):
        if not (int(parameters.A_f_start) <= age + 1 <= int(parameters.A_f_end)):
            continue
        age_years = float(parameters.age_start) + float(age) * float(parameters.da)
        for parity in range(top_parity):
            valid_child_states = (settled,) if parity == 0 else range(parity + 1)
            for child_state in valid_child_states:
                child_units = 1.0 + (top_weight - top_parity if parity == top_parity - 1 else 0.0)
                for bb in range(len(b_grid)):
                    for origin in range(g_pre.shape[1]):
                        for market in range(g_pre.shape[2]):
                            for income_state in range(g_pre.shape[4]):
                                mass = float(
                                    g_pre[
                                        bb,
                                        origin,
                                        market,
                                        age,
                                        income_state,
                                        parity,
                                        child_state,
                                    ]
                                )
                                attempt, gap = fertility_attempt(
                                    evaluation.policy,
                                    parameters,
                                    wealth_index=bb,
                                    origin_tenure=origin,
                                    market=market,
                                    age_index=age,
                                    income_state=income_state,
                                    parity=parity,
                                    child_state=child_state,
                                )
                                realized = float(fecundity[age]) * attempt
                                contribution = mass * realized
                                explicit_births += contribution
                                adjusted_births += contribution * child_units
                                at_risk_mass += mass
                                rows.append(
                                    {
                                        "calendar_year": int(calendar_year),
                                        "period": int(period),
                                        "wealth_index": int(bb),
                                        "liquid_wealth": float(b_grid[bb]),
                                        "origin_tenure": int(origin),
                                        "market": int(market),
                                        "age_index": int(age),
                                        "age_years": age_years,
                                        "age_group": age_group(age_years),
                                        "income_state": int(income_state),
                                        "income_multiplier": float(z_grid[income_state]),
                                        "parity": int(parity),
                                        "child_state": int(child_state),
                                        "parent_status": parent_status(parity, child_state),
                                        "pre_fertility_mass": mass,
                                        "attempt_probability": attempt,
                                        "fecundity_probability": float(fecundity[age]),
                                        "realized_birth_probability": realized,
                                        "birth_child_units_if_realized": child_units,
                                        "expected_birth_child_units": realized * child_units,
                                        "birth_value_gap_attempt_minus_wait": gap,
                                    }
                                )
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        parameters,
    )
    target_adjusted = float(accounting["topcode_adjusted_birth_children"])
    if abs(explicit_births - float(evaluation.births)) > 2.0e-10:
        raise RuntimeError(
            f"Explicit birth reconstruction failed in {calendar_year}: "
            f"{explicit_births} versus {evaluation.births}"
        )
    if abs(adjusted_births - target_adjusted) > 2.0e-10:
        raise RuntimeError(
            f"Top-code birth reconstruction failed in {calendar_year}: "
            f"{adjusted_births} versus {target_adjusted}"
        )
    return rows, {
        "at_risk_mass": at_risk_mass,
        "explicit_births": explicit_births,
        "topcode_adjusted_births": adjusted_births,
        "birth_reconstruction_gap": adjusted_births - target_adjusted,
    }


def cross_section_rows(
    evaluation: calendar.PeriodEvaluation,
    parameters: SimpleNamespace,
    b_grid: np.ndarray,
    *,
    calendar_year: int,
    period: int,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    g_current = np.asarray(evaluation.g_current, dtype=float)
    shared = calendar.model.precompute_shared(parameters, b_grid)
    groups: dict[tuple[str, str, str], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    housing_total = 0.0
    for bb, tenure, market, age, income_state, parity, child_state in np.argwhere(
        g_current > MASS_FLOOR
    ):
        mass = float(
            g_current[bb, tenure, market, age, income_state, parity, child_state]
        )
        age_years = float(parameters.age_start) + float(age) * float(parameters.da)
        if tenure == 0:
            housing = float(
                evaluation.policy.hR_pol[
                    bb, tenure, market, age, income_state, parity, child_state
                ]
            )
        else:
            housing = float(parameters.H_own[tenure - 1])
        floor = float(shared.h_bar[parity, child_state])
        binding = float(floor > 0.0 and housing <= floor + 1.0e-8)
        tax = (
            float(parameters.tau_H)
            * float(evaluation.policy.price[market])
            * housing
            if tenure > 0
            else 0.0
        )
        key = (
            age_group(age_years),
            parent_status(int(parity), int(child_state)),
            "owner" if tenure > 0 else "renter",
        )
        item = groups[key]
        item["mass"] += mass
        item["housing_services"] += mass * housing
        item["child_floor_binding_mass"] += mass * binding
        item["property_tax"] += mass * tax
        item["liquid_wealth"] += mass * float(b_grid[bb])
        housing_total += mass * housing
    rows: list[dict[str, Any]] = []
    for (age_label, parent_label, tenure_label), item in sorted(groups.items()):
        mass = item["mass"]
        rows.append(
            {
                "calendar_year": int(calendar_year),
                "period": int(period),
                "age_group": age_label,
                "parent_status": parent_label,
                "tenure": tenure_label,
                "mass": mass,
                "mass_share": mass / max(float(np.sum(g_current)), 1.0e-15),
                "mean_housing_services": item["housing_services"] / mass,
                "child_floor_binding_share": item["child_floor_binding_mass"] / mass,
                "mean_property_tax_per_head": item["property_tax"] / mass,
                "mean_liquid_wealth": item["liquid_wealth"] / mass,
            }
        )
    target_housing = float(np.sum(evaluation.demand_by_loc))
    if abs(housing_total - target_housing) > 2.0e-9:
        raise RuntimeError(
            f"Housing reconstruction failed in {calendar_year}: "
            f"{housing_total} versus {target_housing}"
        )
    return rows, {
        "current_mass": float(np.sum(g_current)),
        "housing_demand": housing_total,
        "housing_reconstruction_gap": housing_total - target_housing,
    }


def write_year_diagnostic(
    context: CaptureContext,
    evaluation: calendar.PeriodEvaluation,
    parameters: SimpleNamespace,
) -> None:
    period = int(context.current_period)
    year = pf.CALENDAR_START_YEAR + 4 * period
    if year not in context.target_years:
        return
    b_grid = np.asarray(context.b_grid, dtype=float)
    fertility_rows, fertility_validation = fertility_state_rows(
        evaluation, parameters, b_grid, calendar_year=year, period=period
    )
    cross_rows, cross_validation = cross_section_rows(
        evaluation, parameters, b_grid, calendar_year=year, period=period
    )
    fertility_path = context.output_dir / f"fertility_states_{year}.csv.gz"
    cross_path = context.output_dir / f"cross_section_groups_{year}.csv"
    write_csv(fertility_path, fertility_rows)
    write_csv(cross_path, cross_rows)
    context.captured_years.add(year)
    context.validations.append(
        {
            "calendar_year": year,
            "fertility_state_rows": len(fertility_rows),
            "fertility_state_csv": str(fertility_path),
            "fertility_state_csv_sha256": sha256(fertility_path),
            "cross_section_rows": len(cross_rows),
            "cross_section_csv": str(cross_path),
            "cross_section_csv_sha256": sha256(cross_path),
            **fertility_validation,
            **cross_validation,
        }
    )


def install_capture(context: CaptureContext, b_grid: np.ndarray) -> Callable[[], None]:
    original_path = person_pf.evaluate_path_at_prices_person_demography
    original_period = calendar.evaluate_period

    def wrapped_period(*args: Any, **kwargs: Any) -> calendar.PeriodEvaluation:
        result = original_period(*args, **kwargs)
        if context.active:
            parameters = kwargs.get("P", args[2])
            write_year_diagnostic(context, result, parameters)
            context.current_period += 1
        return result

    def wrapped_path(*args: Any, **kwargs: Any) -> person_pf.PersonPathEvaluation:
        if context.evaluation_calls != 0:
            raise RuntimeError("Diagnostic expected exactly one immutable path evaluation")
        context.evaluation_calls += 1
        context.current_period = 0
        context.active = True
        try:
            result = original_path(*args, **kwargs)
        finally:
            context.active = False
        if context.current_period != context.horizon:
            raise RuntimeError(
                f"Captured {context.current_period} forward dates; expected {context.horizon}"
            )
        return result

    calendar.evaluate_period = wrapped_period
    person_pf.evaluate_path_at_prices_person_demography = wrapped_path

    def restore() -> None:
        calendar.evaluate_period = original_period
        person_pf.evaluate_path_at_prices_person_demography = original_path

    return restore


def main() -> None:
    args = parse_args()
    if int(args.horizon) < 2:
        raise ValueError("Horizon must be at least two")
    output_dir = args.output_dir.resolve()
    launcher_files = {"launch_contract.json", "driver.log", "heartbeat.json"}
    unexpected = (
        sorted(path.name for path in output_dir.iterdir() if path.name not in launcher_files)
        if output_dir.exists()
        else []
    )
    if unexpected:
        raise FileExistsError(f"Refusing to overwrite non-launcher output: {unexpected}")
    output_dir.mkdir(parents=True, exist_ok=True)
    years = parse_years(args.diagnostic_years, int(args.horizon))

    # Reconstruct only enough to learn the exact wealth grid before installing
    # the capture. The production driver independently repeats and audits this.
    chain, model = transition.configure_sequential_model()
    contracts, _ = pf.load_diagnostic_contracts(
        args.selected_report.resolve(),
        args.selected_transition.resolve(),
        args.source.resolve(),
    )
    import run_e5f_post2023_no_policy_continuations as continuation

    prepared = continuation.prepare_model(contracts, args.source.resolve(), chain, model)
    context = CaptureContext(
        output_dir=output_dir,
        target_years=years,
        horizon=int(args.horizon),
        b_grid=np.asarray(prepared.b_grid, dtype=float),
    )
    restore = install_capture(context, np.asarray(prepared.b_grid, dtype=float))
    original_argv = sys.argv[:]
    sys.argv = [
        str(Path(policy_driver.__file__).resolve()),
        "--case",
        args.case,
        "--terminal-summary",
        str(args.terminal_summary.resolve()),
        "--initial-path",
        str(args.initial_path.resolve()),
        "--output-dir",
        str(output_dir),
        "--horizon",
        str(int(args.horizon)),
        "--maximum-path-iterations",
        "1",
        "--post2023-tenure-choice-kappa",
        str(float(args.post2023_tenure_choice_kappa)),
        "--selected-report",
        str(args.selected_report.resolve()),
        "--selected-transition",
        str(args.selected_transition.resolve()),
        "--source",
        str(args.source.resolve()),
        "--source-dir",
        str(args.source_dir.resolve()),
        "--headship-dir",
        str(args.headship_dir.resolve()),
    ]
    try:
        policy_driver.main()
    finally:
        sys.argv = original_argv
        restore()
    if context.evaluation_calls != 1 or context.captured_years != years:
        raise RuntimeError(
            "Diagnostic capture incomplete: "
            f"calls={context.evaluation_calls}, requested={sorted(years)}, "
            f"captured={sorted(context.captured_years)}"
        )
    production_summary = output_dir / "summary.json"
    production = json.loads(production_summary.read_text(encoding="utf-8"))
    if not bool(production.get("path_converged")):
        raise RuntimeError("The immutable diagnostic seed did not reproduce a converged path")
    payload = {
        "status": "complete_unpromoted_fertility_mechanism_diagnostic",
        "interpretation": (
            "Default-off readout of exact fertility probabilities and realized housing "
            "along one fixed H128 path; it does not change the Bellman problem or law of motion."
        ),
        "case": args.case,
        "horizon": int(args.horizon),
        "diagnostic_years": sorted(years),
        "post2023_tenure_choice_kappa": float(args.post2023_tenure_choice_kappa),
        "initial_path": str(args.initial_path.resolve()),
        "initial_path_sha256": sha256(args.initial_path.resolve()),
        "terminal_summary": str(args.terminal_summary.resolve()),
        "terminal_summary_sha256": sha256(args.terminal_summary.resolve()),
        "production_summary_sha256": sha256(production_summary),
        "source_hashes": {
            "diagnostic_wrapper": sha256(Path(__file__).resolve()),
            "production_driver": sha256(Path(policy_driver.__file__).resolve()),
            "person_demography": sha256(Path(person_pf.__file__).resolve()),
            "perfect_foresight": sha256(Path(pf.__file__).resolve()),
            "solver": sha256(Path(calendar.model.__file__).resolve()),
        },
        "validations": sorted(context.validations, key=lambda row: row["calendar_year"]),
    }
    write_json(output_dir / "fertility_mechanism_summary.json", payload)


if __name__ == "__main__":
    main()
