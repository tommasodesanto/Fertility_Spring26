#!/usr/bin/env python3
"""Evaluate one fixed H128 path with an isolated tenure-mixing atom.

This default-off oracle leaves the Bellman problem and production transition
driver unchanged.  It replaces deterministic tenure realization only at one
predeclared calendar-time state by an owner share in ``[0,1]``.  The resulting
packet is a complementarity diagnostic, not an equilibrium solution.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

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
import run_e5f_person_policy_tenure_kink_diagnostic as kink  # noqa: E402


MIX_YEAR = 2351
WEALTH_INDEX = 49
ORIGIN_TENURE = 0
MARKET = 0
AGE_INDEX = 13
CHILD_STATE = 0
RENTER_TENURE = 0
OWNER_TENURE = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terminal-summary", type=Path, required=True)
    parser.add_argument("--initial-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--horizon", type=int, default=128)
    parser.add_argument("--mix-calendar-year", type=int, default=MIX_YEAR)
    parser.add_argument("--owner-share", type=float, required=True)
    parser.add_argument("--wealth-index", type=int, default=WEALTH_INDEX)
    parser.add_argument("--origin-tenure", type=int, default=ORIGIN_TENURE)
    parser.add_argument("--market", type=int, default=MARKET)
    parser.add_argument("--age-index", type=int, default=AGE_INDEX)
    parser.add_argument("--child-state", type=int, default=CHILD_STATE)
    parser.add_argument("--renter-tenure", type=int, default=RENTER_TENURE)
    parser.add_argument("--owner-tenure", type=int, default=OWNER_TENURE)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR)
    return parser.parse_args()


def array_sha256(array: np.ndarray) -> str:
    values = np.ascontiguousarray(array)
    return hashlib.sha256(values.view(np.uint8)).hexdigest()


def file_sha256(path: Path) -> str:
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


def deterministic_tenure_probabilities(choices: np.ndarray) -> np.ndarray:
    selected = np.asarray(choices)
    if selected.ndim < 1 or not np.issubdtype(selected.dtype, np.integer):
        raise ValueError("Tenure choices must be a nonempty integer array")
    number_of_tenures = int(selected.shape[1])
    if number_of_tenures < 2 or np.any(selected < 0) or np.any(selected >= number_of_tenures):
        raise ValueError("Tenure choices lie outside the tenure dimension")
    probabilities = np.zeros(selected.shape + (number_of_tenures,), dtype=float)
    np.put_along_axis(probabilities, selected[..., None], 1.0, axis=-1)
    return probabilities


def probabilities_with_atom_share(
    choices: np.ndarray,
    *,
    owner_share: float,
    wealth_index: int,
    origin_tenure: int,
    market: int,
    age_index: int,
    child_state: int,
    renter_tenure: int,
    owner_tenure: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    share = float(owner_share)
    if not math.isfinite(share) or not 0.0 <= share <= 1.0:
        raise ValueError("The owner share must lie in [0,1]")
    selected = np.asarray(choices)
    if selected.ndim != 7:
        raise ValueError(f"Expected seven-dimensional Markov-income choices, got {selected.shape}")
    probabilities = deterministic_tenure_probabilities(selected)
    Nb, nt, I, J, Nz, npar, ncs = selected.shape
    indices = {
        "wealth_index": (int(wealth_index), Nb),
        "origin_tenure": (int(origin_tenure), nt),
        "market": (int(market), I),
        "age_index": (int(age_index), J),
        "child_state": (int(child_state), ncs),
        "renter_tenure": (int(renter_tenure), nt),
        "owner_tenure": (int(owner_tenure), nt),
    }
    for name, (index, bound) in indices.items():
        if index < 0 or index >= bound:
            raise IndexError(f"{name}={index} lies outside [0,{bound})")
    if int(renter_tenure) == int(owner_tenure):
        raise ValueError("Renter and owner tenure indices must differ")

    atom = probabilities[
        int(wealth_index),
        int(origin_tenure),
        int(market),
        int(age_index),
        :,
        :,
        int(child_state),
        :,
    ]
    atom[...] = 0.0
    atom[..., int(renter_tenure)] = 1.0 - share
    atom[..., int(owner_tenure)] = share
    row_sums = np.sum(probabilities, axis=-1)
    if not np.array_equal(row_sums, np.ones_like(row_sums)):
        raise RuntimeError("Tenure probabilities do not sum exactly to one")
    mixed_rows = int(np.count_nonzero(np.sum(probabilities > 0.0, axis=-1) > 1))
    expected_mixed_rows = int(Nz * npar) if 0.0 < share < 1.0 else 0
    if mixed_rows != expected_mixed_rows:
        raise RuntimeError(
            f"Mixed-row count {mixed_rows} differs from {expected_mixed_rows}"
        )
    return probabilities, {
        "owner_share": share,
        "renter_share": 1.0 - share,
        "target_probability_rows": int(Nz * npar),
        "strictly_mixed_probability_rows": mixed_rows,
        "probability_row_sum_maximum_error": float(np.max(np.abs(row_sums - 1.0))),
    }


def complementarity_violation(owner_share: float, owner_minus_renter: np.ndarray) -> float:
    share = float(owner_share)
    gaps = np.asarray(owner_minus_renter, dtype=float)
    if gaps.size == 0 or np.any(~np.isfinite(gaps)):
        raise ValueError("Complementarity gaps must be finite and nonempty")
    if share <= 0.0:
        return float(max(np.max(gaps), 0.0))
    if share >= 1.0:
        return float(max(-np.min(gaps), 0.0))
    return float(np.max(np.abs(gaps)))


def validate_mixing_atom_options(
    option_values: np.ndarray,
    deterministic_choices: np.ndarray,
    *,
    renter_tenure: int,
    owner_tenure: int,
) -> np.ndarray:
    """Require a live renter/owner margin before assigning a fractional share."""

    options = np.asarray(option_values, dtype=float)
    choices = np.asarray(deterministic_choices)
    if options.ndim != 3 or options.shape[:2] != choices.shape:
        raise ValueError(
            "Mixing-atom options must be (income, parity, tenure) and align "
            "with the deterministic choices"
        )
    if not np.issubdtype(choices.dtype, np.integer):
        raise ValueError("Mixing-atom deterministic choices must be integers")
    number_of_tenures = int(options.shape[-1])
    for label, tenure in (
        ("renter_tenure", int(renter_tenure)),
        ("owner_tenure", int(owner_tenure)),
    ):
        if tenure < 0 or tenure >= number_of_tenures:
            raise IndexError(f"{label}={tenure} lies outside [0,{number_of_tenures})")
    if int(renter_tenure) == int(owner_tenure):
        raise ValueError("Renter and owner tenure indices must differ")

    renter_values = options[..., int(renter_tenure)]
    owner_values = options[..., int(owner_tenure)]
    both_live = (
        np.isfinite(renter_values)
        & np.isfinite(owner_values)
        & (renter_values > kink.NEGATIVE_INFINITY / 2.0)
        & (owner_values > kink.NEGATIVE_INFINITY / 2.0)
    )
    if not np.all(both_live):
        bad = np.argwhere(~both_live)
        first = tuple(int(index) for index in bad[0])
        raise RuntimeError(
            "The declared mixing atom contains an infeasible renter/owner branch: "
            f"income={first[0]}, parity={first[1]}, "
            f"renter_value={renter_values[first]:.12g}, "
            f"owner_value={owner_values[first]:.12g}, "
            f"invalid_rows={len(bad)}"
        )
    eligible_choices = (choices == int(renter_tenure)) | (
        choices == int(owner_tenure)
    )
    if not np.all(eligible_choices):
        bad = np.argwhere(~eligible_choices)
        first = tuple(int(index) for index in bad[0])
        raise RuntimeError(
            "The declared renter/owner atom has a third tenure as its deterministic "
            f"argmax: income={first[0]}, parity={first[1]}, "
            f"choice={int(choices[first])}, invalid_rows={len(bad)}"
        )
    return owner_values - renter_values


def load_path_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [int(row["period"]) for row in rows] != list(range(len(rows))):
        raise RuntimeError("The saved transition path is not contiguous")
    return rows


def main() -> None:
    args = parse_args()
    if int(args.horizon) < 2:
        raise ValueError("Horizon must be at least two")
    allowed_years = {
        pf.CALENDAR_START_YEAR + 4 * period for period in range(int(args.horizon))
    }
    if int(args.mix_calendar_year) not in allowed_years:
        raise ValueError("The mix year must lie on the requested model horizon")
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

    context = kink.CaptureContext(
        output_dir=output_dir, target_years={int(args.mix_calendar_year)}
    )
    restore_capture = kink.install_capture(context)
    captured_solve_date = pf.solve_date_policy
    captured_evaluate_period = calendar.evaluate_period
    transformation_records: list[dict[str, Any]] = []
    evaluation_records: list[dict[str, Any]] = []
    pending_mix = {"active": False, "parameters": None}

    def mixed_solve_date_policy(**kwargs: Any) -> calendar.PolicyBundle:
        policy = captured_solve_date(**kwargs)
        if context.current_year != int(args.mix_calendar_year):
            return policy
        if policy.tenure_probs is not None:
            raise RuntimeError("The complementarity oracle requires deterministic baseline tenure")
        choices_before = array_sha256(policy.tenure_choice)
        value_before = array_sha256(policy.V)
        option_rows: list[np.ndarray] = []
        choice_rows: list[np.ndarray] = []
        for income_state in range(int(context.current_nz)):
            captured = context.current_slices[(int(args.age_index), income_state)]
            options = np.asarray(captured["options"], dtype=float)
            option_rows.append(
                options[
                    int(args.wealth_index),
                    int(args.origin_tenure),
                    int(args.market),
                    :,
                    int(args.child_state),
                    :,
                ]
            )
            choice_rows.append(
                np.asarray(policy.tenure_choice)[
                    int(args.wealth_index),
                    int(args.origin_tenure),
                    int(args.market),
                    int(args.age_index),
                    income_state,
                    :,
                    int(args.child_state),
                ]
            )
        gap_array = validate_mixing_atom_options(
            np.stack(option_rows, axis=0),
            np.stack(choice_rows, axis=0),
            renter_tenure=int(args.renter_tenure),
            owner_tenure=int(args.owner_tenure),
        )
        probabilities, probability_record = probabilities_with_atom_share(
            policy.tenure_choice,
            owner_share=float(args.owner_share),
            wealth_index=int(args.wealth_index),
            origin_tenure=int(args.origin_tenure),
            market=int(args.market),
            age_index=int(args.age_index),
            child_state=int(args.child_state),
            renter_tenure=int(args.renter_tenure),
            owner_tenure=int(args.owner_tenure),
        )
        policy.tenure_probs = probabilities
        if choices_before != array_sha256(policy.tenure_choice) or value_before != array_sha256(policy.V):
            raise RuntimeError("The oracle changed the Bellman value or deterministic argmax")
        transformation_records.append(
            {
                **probability_record,
                "calendar_year": int(context.current_year),
                "wealth_index": int(args.wealth_index),
                "origin_tenure": int(args.origin_tenure),
                "market": int(args.market),
                "age_index": int(args.age_index),
                "child_state": int(args.child_state),
                "renter_tenure": int(args.renter_tenure),
                "owner_tenure": int(args.owner_tenure),
                "owner_minus_renter_value_gap_minimum": float(np.min(gap_array)),
                "owner_minus_renter_value_gap_maximum": float(np.max(gap_array)),
                "owner_minus_renter_value_gap_spread": float(np.ptp(gap_array)),
                "complementarity_violation": complementarity_violation(
                    float(args.owner_share), gap_array
                ),
                "deterministic_choice_sha256": choices_before,
                "bellman_value_sha256": value_before,
            }
        )
        pending_mix["active"] = True
        pending_mix["parameters"] = kwargs["P"]
        return policy

    def audited_evaluate_period(*positional: Any, **kwargs: Any) -> calendar.PeriodEvaluation:
        evaluation = captured_evaluate_period(*positional, **kwargs)
        if not pending_mix["active"]:
            return evaluation
        parameters = pending_mix["parameters"]
        if parameters is None or evaluation.policy.tenure_probs is None:
            raise RuntimeError("The mixed evaluation lacks parameters or probabilities")
        share = float(args.owner_share)
        endpoint_evaluations: dict[str, dict[str, float]] = {}
        next_distributions: dict[float, np.ndarray] = {}
        current_distributions: dict[float, np.ndarray] = {}
        for endpoint_share in (0.0, share, 1.0):
            probabilities, _ = probabilities_with_atom_share(
                evaluation.policy.tenure_choice,
                owner_share=endpoint_share,
                wealth_index=int(args.wealth_index),
                origin_tenure=int(args.origin_tenure),
                market=int(args.market),
                age_index=int(args.age_index),
                child_state=int(args.child_state),
                renter_tenure=int(args.renter_tenure),
                owner_tenure=int(args.owner_tenure),
            )
            policy = copy.copy(evaluation.policy)
            policy.tenure_probs = probabilities
            current = calendar.model.realize_current_cross_section(
                evaluation.g_post_fertility,
                policy.loc_probs,
                policy.tenure_choice,
                policy.tenure_probs,
                policy.maps.lmm_idx,
                policy.maps.lmm_wt,
                policy.maps.tmx_idx,
                policy.maps.tmx_wt,
                use_compiled_scatter=bool(getattr(parameters, "use_numba_scatter", False)),
            )
            demand = calendar.housing_demand_by_location(current, policy.hR_pol, parameters)
            tax = calendar.model.property_tax_revenue_from_distribution(
                current, policy.hR_pol, policy.price, parameters
            )
            endpoint_evaluation = replace(
                evaluation,
                policy=policy,
                g_current=current,
                demand_by_loc=demand,
            )
            raw_next, _, _, mass_residual = transition.advance_sequential_calendar_distribution(
                endpoint_evaluation,
                np.zeros(int(parameters.I)),
                parameters,
                kwargs.get("b_grid", positional[3]),
                kwargs.get("shared", positional[4]),
            )
            mass = float(np.sum(current))
            current_distributions[endpoint_share] = current
            next_distributions[endpoint_share] = raw_next
            endpoint_evaluations[str(endpoint_share)] = {
                "owner_rate": float(np.sum(current[:, 1:, ...])) / max(mass, 1e-15),
                "housing_demand": float(np.sum(demand)),
                "property_tax_revenue": float(tax),
                "raw_next_mass": float(np.sum(raw_next)),
                "raw_next_mass_residual": float(mass_residual),
            }
        expected_current = (
            (1.0 - share) * current_distributions[0.0]
            + share * current_distributions[1.0]
        )
        expected_next = (
            (1.0 - share) * next_distributions[0.0]
            + share * next_distributions[1.0]
        )
        current_gap = float(
            np.max(np.abs(current_distributions[share] - expected_current))
        )
        next_gap = float(np.max(np.abs(next_distributions[share] - expected_next)))
        actual_current_gap = float(
            np.max(np.abs(np.asarray(evaluation.g_current) - current_distributions[share]))
        )
        if current_gap > 2e-14 or next_gap > 2e-14 or actual_current_gap > 2e-14:
            raise RuntimeError(
                "Tenure mixing is inconsistent across current and forward laws: "
                f"current={current_gap:.6g}, next={next_gap:.6g}, "
                f"production_current={actual_current_gap:.6g}"
            )
        evaluation_records.append(
            {
                "calendar_year": int(args.mix_calendar_year),
                "endpoint_evaluations": endpoint_evaluations,
                "current_distribution_linearity_maximum_error": current_gap,
                "forward_distribution_linearity_maximum_error": next_gap,
                "production_current_distribution_maximum_error": actual_current_gap,
            }
        )
        pending_mix["active"] = False
        pending_mix["parameters"] = None
        return evaluation

    pf.solve_date_policy = mixed_solve_date_policy
    calendar.evaluate_period = audited_evaluate_period
    original_argv = sys.argv[:]
    sys.argv = [
        str(Path(policy_driver.__file__).resolve()),
        "--case",
        "rebated-tax2-reform",
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
        calendar.evaluate_period = captured_evaluate_period
        pf.solve_date_policy = captured_solve_date
        restore_capture()

    if len(transformation_records) != 1 or len(evaluation_records) != 1:
        raise RuntimeError(
            "The oracle must transform and audit exactly one calendar-time policy"
        )
    summary_path = output_dir / "summary.json"
    transition_path = output_dir / "transition_path.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = load_path_rows(transition_path)
    mix_rows = [row for row in rows if int(row["calendar_year"]) == int(args.mix_calendar_year)]
    if len(mix_rows) != 1:
        raise RuntimeError("The production path does not contain exactly one mix-year row")
    if summary.get("history_reproduction_status") != "passed" or not all(
        (summary.get("accounting_gates") or {}).values()
    ):
        raise RuntimeError("The mixed oracle failed a production accounting gate")
    if max(abs(float(row["feasibility_frontier_projection_mass"])) for row in rows) != 0.0:
        raise RuntimeError("The mixed oracle needed a feasibility projection")
    mix_row = mix_rows[0]
    payload = {
        "status": "complete_unaccepted_tenure_complementarity_oracle",
        "interpretation": (
            "One fixed-path evaluation with an atom-specific zero-dispersion tenure "
            "share. This is a complementarity oracle, not a solved transition."
        ),
        "case": "rebated-tax2-reform",
        "horizon": int(args.horizon),
        "initial_path": str(args.initial_path.resolve()),
        "initial_path_sha256": file_sha256(args.initial_path.resolve()),
        "terminal_summary": str(args.terminal_summary.resolve()),
        "terminal_summary_sha256": file_sha256(args.terminal_summary.resolve()),
        "production_summary_sha256": file_sha256(summary_path),
        "production_transition_path_sha256": file_sha256(transition_path),
        "source_hashes": {
            "oracle": file_sha256(Path(__file__).resolve()),
            "kink_capture": file_sha256(Path(kink.__file__).resolve()),
            "production_driver": file_sha256(Path(policy_driver.__file__).resolve()),
            "person_demography": file_sha256(Path(person_pf.__file__).resolve()),
            "perfect_foresight": file_sha256(Path(pf.__file__).resolve()),
            "solver": file_sha256(Path(calendar.model.__file__).resolve()),
        },
        "transformation": transformation_records[0],
        "evaluation_validation": evaluation_records[0],
        "mix_year_equilibrium_residuals": {
            "calendar_year": int(args.mix_calendar_year),
            "market": float(mix_row["relative_market_residual"]),
            "fiscal": float(mix_row["government_budget_residual"]),
            "owner_rate": float(mix_row["owner_rate"]),
            "housing_demand": float(mix_row["housing_demand"]),
            "property_tax_revenue": float(mix_row["property_tax_revenue"]),
        },
        "production_path_status": summary.get("path_status"),
        "production_accounting_gates": summary.get("accounting_gates"),
        "history_reproduction_status": summary.get("history_reproduction_status"),
    }
    write_json(output_dir / "tenure_complementarity_oracle_summary.json", payload)


if __name__ == "__main__":
    main()
