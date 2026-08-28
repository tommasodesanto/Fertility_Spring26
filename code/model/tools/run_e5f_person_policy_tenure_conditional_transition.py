#!/usr/bin/env python3
"""Solve a PF policy path conditional on one atom-specific tenure share.

This default-off wrapper leaves the household Bellman problem and the ordinary
production driver unchanged. During every path evaluation it replaces tenure
realization only at one predeclared live renter/owner atom. The resulting path
is conditional on the supplied share; it is a complementarity equilibrium only
if its final value gap also satisfies the relevant boundary/interior condition.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
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
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402
import run_e5f_perfect_foresight_person_demography_policy as policy_driver  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402
import run_e5f_person_policy_tenure_complementarity_oracle as oracle  # noqa: E402
import run_e5f_person_policy_tenure_kink_diagnostic as kink  # noqa: E402


@dataclass(frozen=True)
class AtomSpec:
    calendar_year: int = oracle.MIX_YEAR
    wealth_index: int = oracle.WEALTH_INDEX
    origin_tenure: int = oracle.ORIGIN_TENURE
    market: int = oracle.MARKET
    age_index: int = oracle.AGE_INDEX
    child_state: int = oracle.CHILD_STATE
    renter_tenure: int = oracle.RENTER_TENURE
    owner_tenure: int = oracle.OWNER_TENURE


@dataclass
class ConditionalCapture:
    output_dir: Path
    atom: AtomSpec
    owner_share: float
    horizon: int = 0
    active: bool = False
    solve_count: int = 0
    evaluation_count: int = 0
    capture_target: bool = False
    current_period: int | None = None
    current_parameters: SimpleNamespace | None = None
    current_nz: int = 0
    kernel_calls: int = 0
    target_slices: dict[int, np.ndarray] = field(default_factory=dict)
    live_value_errors: list[float] = field(default_factory=list)
    full_value_errors: list[float] = field(default_factory=list)
    transformations: list[dict[str, Any]] = field(default_factory=list)
    evaluation_records: list[dict[str, Any]] = field(default_factory=list)
    contract_hashes: dict[str, Any] = field(default_factory=dict)
    started: float = field(default_factory=time.perf_counter)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terminal-summary", type=Path, required=True)
    parser.add_argument("--initial-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--horizon", type=int, default=128)
    parser.add_argument("--owner-share", type=float, required=True)
    parser.add_argument("--mix-calendar-year", type=int, default=oracle.MIX_YEAR)
    parser.add_argument("--wealth-index", type=int, default=oracle.WEALTH_INDEX)
    parser.add_argument("--origin-tenure", type=int, default=oracle.ORIGIN_TENURE)
    parser.add_argument("--market", type=int, default=oracle.MARKET)
    parser.add_argument("--age-index", type=int, default=oracle.AGE_INDEX)
    parser.add_argument("--child-state", type=int, default=oracle.CHILD_STATE)
    parser.add_argument("--renter-tenure", type=int, default=oracle.RENTER_TENURE)
    parser.add_argument("--owner-tenure", type=int, default=oracle.OWNER_TENURE)
    parser.add_argument("--maximum-path-iterations", type=int, default=35)
    parser.add_argument("--price-damping", type=float, default=0.25)
    parser.add_argument("--transfer-damping", type=float, default=0.50)
    parser.add_argument("--maximum-log-price-step", type=float, default=0.12)
    parser.add_argument("--maximum-transfer-step", type=float, default=0.08)
    parser.add_argument("--psi-persistence", type=float, default=pf.DEFAULT_PSI_PERSISTENCE)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR)
    return parser.parse_args()


def write_json_atomic(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(pf.jsonable(payload), indent=2, sort_keys=True, default=str)
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def write_csv_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty checkpoint path: {path}")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_atom_options_and_choices(
    policy: calendar.PolicyBundle,
    slices: dict[int, np.ndarray],
    *,
    number_of_income_states: int,
    atom: AtomSpec,
) -> tuple[np.ndarray, np.ndarray]:
    if sorted(slices) != list(range(int(number_of_income_states))):
        raise RuntimeError("The conditional solver did not capture every income slice")
    option_rows: list[np.ndarray] = []
    choice_rows: list[np.ndarray] = []
    selected = np.asarray(policy.tenure_choice)
    for income_state in range(int(number_of_income_states)):
        options = np.asarray(slices[income_state], dtype=float)
        option_rows.append(
            options[
                atom.wealth_index,
                atom.origin_tenure,
                atom.market,
                :,
                atom.child_state,
                :,
            ]
        )
        choice_rows.append(
            selected[
                atom.wealth_index,
                atom.origin_tenure,
                atom.market,
                atom.age_index,
                income_state,
                :,
                atom.child_state,
            ]
        )
    return np.stack(option_rows, axis=0), np.stack(choice_rows, axis=0)


def path_array_hashes(rows: list[dict[str, Any]]) -> dict[str, str]:
    return {
        "price_path_sha256": oracle.array_sha256(
            np.asarray([float(row["asset_price"]) for row in rows], dtype=float)
        ),
        "transfer_path_sha256": oracle.array_sha256(
            np.asarray(
                [float(row["equal_transfer_period_units"]) for row in rows],
                dtype=float,
            )
        ),
    }


def source_hashes() -> dict[str, str]:
    return {
        "conditional_driver": file_sha256(Path(__file__).resolve()),
        "oracle_helper": file_sha256(Path(oracle.__file__).resolve()),
        "kink_helper": file_sha256(Path(kink.__file__).resolve()),
        "production_driver": file_sha256(Path(policy_driver.__file__).resolve()),
        "person_demography": file_sha256(Path(person_pf.__file__).resolve()),
        "perfect_foresight": file_sha256(Path(pf.__file__).resolve()),
        "funded_policy": file_sha256(Path(policy_driver.funded.__file__).resolve()),
        "continuation": file_sha256(Path(policy_driver.continuation.__file__).resolve()),
        "solver": file_sha256(
            MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
        ),
    }


def audit_production_packet(
    *,
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
    final_record: dict[str, Any],
    expected_horizon: int,
) -> dict[str, Any]:
    if len(rows) != int(expected_horizon):
        raise RuntimeError("Production path row count differs from the horizon")
    if [int(row["period"]) for row in rows] != list(range(int(expected_horizon))):
        raise RuntimeError("Production path periods are not contiguous")
    try:
        market = max(abs(float(row["relative_market_residual"])) for row in rows)
        fiscal = max(abs(float(row["government_budget_residual"])) for row in rows)
        feasibility = max(
            abs(float(row["feasibility_frontier_projection_mass"])) for row in rows
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("Production path lacks finite audit fields") from exc
    if not all(math.isfinite(value) for value in (market, fiscal, feasibility)):
        raise RuntimeError("Production path contains nonfinite audit fields")
    if not math.isclose(
        market,
        float(summary["maximum_market_residual"]),
        rel_tol=1e-12,
        abs_tol=1e-14,
    ) or not math.isclose(
        fiscal,
        float(summary["maximum_fiscal_residual"]),
        rel_tol=1e-12,
        abs_tol=1e-14,
    ):
        raise RuntimeError("Production summary residuals differ from its saved path")
    recomputed_converged = (
        market <= policy_driver.MARKET_TOLERANCE
        and fiscal <= policy_driver.FISCAL_ABSOLUTE_TOLERANCE
    )
    if bool(summary.get("path_converged")) != recomputed_converged:
        raise RuntimeError("Production convergence flag differs from saved residuals")
    if summary.get("history_reproduction_status") != "passed":
        raise RuntimeError("Historical-state reproduction did not pass")
    accounting = summary.get("accounting_gates") or {}
    if not accounting or not all(accounting.values()):
        raise RuntimeError("A production accounting gate failed")
    terminal = summary.get("terminal_convergence") or {}
    if not terminal or not bool(terminal.get("all_checks_pass")):
        raise RuntimeError("The H128 terminal-state gate did not pass")
    if feasibility != 0.0:
        raise RuntimeError("The conditional path needed a feasibility projection")
    if not math.isclose(
        market,
        float(final_record["maximum_market_residual"]),
        rel_tol=1e-12,
        abs_tol=1e-14,
    ) or not math.isclose(
        fiscal,
        float(final_record["maximum_fiscal_residual"]),
        rel_tol=1e-12,
        abs_tol=1e-14,
    ):
        raise RuntimeError("Final conditional checkpoint differs from production residuals")
    if float(final_record["maximum_feasibility_projection_mass"]) != 0.0:
        raise RuntimeError("Final conditional checkpoint used feasibility projection")
    for name, digest in path_array_hashes(rows).items():
        if final_record.get(name) != digest:
            raise RuntimeError(
                f"Final production path differs from final conditional record: {name}"
            )
    return {
        "maximum_market_residual": market,
        "maximum_fiscal_residual": fiscal,
        "maximum_feasibility_projection_mass": feasibility,
        "path_converged": recomputed_converged,
        "history_reproduction_passed": True,
        "accounting_gates_passed": True,
        "terminal_gates_passed": True,
    }


def install_conditional_share(context: ConditionalCapture) -> Callable[[], None]:
    original_evaluate_path = person_pf.evaluate_path_at_prices_person_demography
    original_solve_date = pf.solve_date_policy
    active_kernel: dict[str, Any] = {"model": None, "original": None}

    def wrapped_tenure_kernel(*args: Any) -> tuple[np.ndarray, np.ndarray]:
        original_kernel = active_kernel["original"]
        if original_kernel is None:
            raise RuntimeError("Conditional tenure kernel is not active")
        values, choices = original_kernel(*args)
        if context.capture_target:
            parameters = context.current_parameters
            if parameters is None:
                raise RuntimeError("Conditional capture lacks model parameters")
            call = context.kernel_calls
            context.kernel_calls += 1
            age, income_state = kink.bellman_slice_indices(
                call, int(parameters.J), int(context.current_nz)
            )
            if age == int(context.atom.age_index):
                Vd, bg, heq, hcost, dp_arr, bmo, birth_dp, grant = args
                options = kink.reconstruct_tenure_option_values(
                    Vd, bg, heq, hcost, dp_arr, bmo, birth_dp, grant
                )
                reconstructed_choices = np.argmax(options, axis=-1).astype(np.int16)
                mismatch = int(np.count_nonzero(reconstructed_choices != choices))
                maximum = np.max(options, axis=-1)
                errors = np.abs(maximum - values)
                live = np.maximum(maximum, values) > kink.NEGATIVE_INFINITY / 2.0
                live_error = float(np.max(errors[live])) if np.any(live) else 0.0
                full_error = float(np.max(errors))
                if mismatch or live_error > 2.0e-11 or full_error > 2.0e-6:
                    raise RuntimeError(
                        "Conditional option reconstruction differs from the tenure kernel: "
                        f"mismatches={mismatch}, live={live_error:.6g}, full={full_error:.6g}"
                    )
                context.target_slices[income_state] = options
                context.live_value_errors.append(live_error)
                context.full_value_errors.append(full_error)
        return values, choices

    def wrapped_solve_date_policy(**kwargs: Any) -> calendar.PolicyBundle:
        if not context.active:
            return original_solve_date(**kwargs)
        call = context.solve_count
        context.solve_count += 1
        capture_this_call = False
        if call >= context.horizon:
            period = call - context.horizon
            year = pf.CALENDAR_START_YEAR + 4 * period
            capture_this_call = year == int(context.atom.calendar_year)
        if capture_this_call:
            if context.capture_target:
                raise RuntimeError("Nested conditional tenure capture")
            context.capture_target = True
            context.current_period = period
            context.current_parameters = kwargs["P"]
            context.current_nz = len(calendar.model.income_transition_values(kwargs["P"])[0])
            context.kernel_calls = 0
            context.target_slices = {}
            context.live_value_errors = []
            context.full_value_errors = []
        try:
            policy = original_solve_date(**kwargs)
        except Exception:
            context.capture_target = False
            raise
        if not capture_this_call:
            return policy

        parameters = context.current_parameters
        if parameters is None:
            raise RuntimeError("Conditional target solve lost its parameters")
        expected_calls = int(parameters.J) * int(context.current_nz)
        if context.kernel_calls != expected_calls:
            raise RuntimeError(
                f"Captured {context.kernel_calls} tenure calls; expected {expected_calls}"
            )
        options, choices = extract_atom_options_and_choices(
            policy,
            context.target_slices,
            number_of_income_states=context.current_nz,
            atom=context.atom,
        )
        gaps = oracle.validate_mixing_atom_options(
            options,
            choices,
            renter_tenure=context.atom.renter_tenure,
            owner_tenure=context.atom.owner_tenure,
        )
        choices_hash = oracle.array_sha256(policy.tenure_choice)
        values_hash = oracle.array_sha256(policy.V)
        probabilities, probability_record = oracle.probabilities_with_atom_share(
            policy.tenure_choice,
            owner_share=context.owner_share,
            wealth_index=context.atom.wealth_index,
            origin_tenure=context.atom.origin_tenure,
            market=context.atom.market,
            age_index=context.atom.age_index,
            child_state=context.atom.child_state,
            renter_tenure=context.atom.renter_tenure,
            owner_tenure=context.atom.owner_tenure,
        )
        policy.tenure_probs = probabilities
        if choices_hash != oracle.array_sha256(policy.tenure_choice) or values_hash != oracle.array_sha256(policy.V):
            raise RuntimeError("Conditional tenure realization changed Bellman objects")
        context.transformations.append(
            {
                **probability_record,
                "evaluation": context.evaluation_count,
                "period": int(context.current_period),
                "calendar_year": int(context.atom.calendar_year),
                "wealth_index": int(context.atom.wealth_index),
                "origin_tenure": int(context.atom.origin_tenure),
                "market": int(context.atom.market),
                "age_index": int(context.atom.age_index),
                "child_state": int(context.atom.child_state),
                "renter_tenure": int(context.atom.renter_tenure),
                "owner_tenure": int(context.atom.owner_tenure),
                "owner_minus_renter_value_gap_minimum": float(np.min(gaps)),
                "owner_minus_renter_value_gap_maximum": float(np.max(gaps)),
                "owner_minus_renter_value_gap_spread": float(np.ptp(gaps)),
                "complementarity_violation": oracle.complementarity_violation(
                    context.owner_share, gaps
                ),
                "maximum_live_value_reconstruction_error": max(
                    context.live_value_errors, default=0.0
                ),
                "maximum_full_value_reconstruction_error": max(
                    context.full_value_errors, default=0.0
                ),
                "deterministic_choice_sha256": choices_hash,
                "bellman_value_sha256": values_hash,
            }
        )
        context.capture_target = False
        context.current_period = None
        context.current_parameters = None
        context.current_nz = 0
        context.target_slices = {}
        return policy

    def wrapped_evaluate_path(*args: Any, **kwargs: Any) -> person_pf.PersonPathEvaluation:
        prices = np.asarray(kwargs.get("prices", args[0] if args else []), dtype=float)
        context.horizon = int(prices.size)
        context.solve_count = 0
        context.evaluation_count += 1
        context.active = True
        transformation_count = len(context.transformations)
        active_model = calendar.model
        original_kernel = active_model.tenure_choice_kernel
        active_kernel["model"] = active_model
        active_kernel["original"] = original_kernel
        active_model.tenure_choice_kernel = wrapped_tenure_kernel
        try:
            result = original_evaluate_path(*args, **kwargs)
        finally:
            active_model.tenure_choice_kernel = original_kernel
            active_kernel["model"] = None
            active_kernel["original"] = None
            context.active = False
            context.capture_target = False
        if context.solve_count != 2 * context.horizon:
            raise RuntimeError(
                f"Conditional evaluation observed {context.solve_count} policy solves; "
                f"expected {2 * context.horizon}"
            )
        if len(context.transformations) != transformation_count + 1:
            raise RuntimeError("Each path evaluation must transform the atom exactly once")
        transformation = context.transformations[-1]
        mix_rows = [
            row
            for row in result.rows
            if int(row["calendar_year"]) == int(context.atom.calendar_year)
        ]
        if len(mix_rows) != 1:
            raise RuntimeError("Conditional evaluation lacks one mix-year row")
        record = {
            "evaluation": context.evaluation_count,
            "elapsed_seconds": time.perf_counter() - context.started,
            "maximum_market_residual": float(result.maximum_market_residual),
            "maximum_fiscal_residual": max(
                abs(float(row["government_budget_residual"])) for row in result.rows
            ),
            "maximum_policy_reproduction_error": float(
                result.maximum_policy_reproduction_error
            ),
            "maximum_person_identity_error": float(result.maximum_person_identity_error),
            "maximum_head_identity_error": float(result.maximum_head_identity_error),
            "maximum_household_person_head_gap": float(
                result.maximum_household_person_head_gap
            ),
            "maximum_age_head_gap": float(result.maximum_age_head_gap),
            "maximum_feasibility_projection_mass": float(
                result.maximum_feasibility_projection_mass
            ),
            "mix_year_market_residual": float(mix_rows[0]["relative_market_residual"]),
            "mix_year_fiscal_residual": float(mix_rows[0]["government_budget_residual"]),
            "mix_year_owner_rate": float(mix_rows[0]["owner_rate"]),
            **path_array_hashes(result.rows),
            "contract_hashes": context.contract_hashes,
            "transformation": transformation,
        }
        checkpoints = context.output_dir / "conditional_checkpoints"
        checkpoints.mkdir(parents=True, exist_ok=True)
        stem = f"evaluation_{context.evaluation_count:04d}"
        path_checkpoint = checkpoints / f"{stem}.csv"
        json_checkpoint = checkpoints / f"{stem}.json"
        if path_checkpoint.exists() or json_checkpoint.exists():
            raise FileExistsError(f"Refusing to overwrite conditional checkpoint {stem}")
        write_csv_atomic(path_checkpoint, result.rows)
        record["checkpoint_path"] = str(path_checkpoint)
        record["checkpoint_path_sha256"] = file_sha256(path_checkpoint)
        context.evaluation_records.append(record)
        write_json_atomic(json_checkpoint, record)
        write_json_atomic(
            context.output_dir / "latest_conditional_evaluation.json", record
        )
        write_json_atomic(
            context.output_dir / "conditional_evaluation_history.json",
            context.evaluation_records,
        )
        return result

    pf.solve_date_policy = wrapped_solve_date_policy
    person_pf.evaluate_path_at_prices_person_demography = wrapped_evaluate_path

    def restore() -> None:
        if active_kernel["model"] is not None and active_kernel["original"] is not None:
            active_kernel["model"].tenure_choice_kernel = active_kernel["original"]
        pf.solve_date_policy = original_solve_date
        person_pf.evaluate_path_at_prices_person_demography = original_evaluate_path

    return restore


def main() -> None:
    args = parse_args()
    if int(args.horizon) < 2:
        raise ValueError("Horizon must be at least two")
    if int(args.maximum_path_iterations) < 1:
        raise ValueError("Maximum path iterations must be positive")
    allowed_years = {
        pf.CALENDAR_START_YEAR + 4 * period for period in range(int(args.horizon))
    }
    if int(args.mix_calendar_year) not in allowed_years:
        raise ValueError("Mix year must lie on the requested horizon")
    if not math.isfinite(float(args.owner_share)) or not 0.0 <= float(args.owner_share) <= 1.0:
        raise ValueError("Owner share must lie in [0,1]")
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
    atom = AtomSpec(
        calendar_year=int(args.mix_calendar_year),
        wealth_index=int(args.wealth_index),
        origin_tenure=int(args.origin_tenure),
        market=int(args.market),
        age_index=int(args.age_index),
        child_state=int(args.child_state),
        renter_tenure=int(args.renter_tenure),
        owner_tenure=int(args.owner_tenure),
    )
    context = ConditionalCapture(
        output_dir=output_dir,
        atom=atom,
        owner_share=float(args.owner_share),
        contract_hashes={
            "initial_path_sha256": file_sha256(args.initial_path.resolve()),
            "terminal_summary_sha256": file_sha256(args.terminal_summary.resolve()),
            "source_hashes": source_hashes(),
        },
    )
    restore = install_conditional_share(context)
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
        str(int(args.maximum_path_iterations)),
        "--price-damping",
        str(float(args.price_damping)),
        "--transfer-damping",
        str(float(args.transfer_damping)),
        "--maximum-log-price-step",
        str(float(args.maximum_log_price_step)),
        "--maximum-transfer-step",
        str(float(args.maximum_transfer_step)),
        "--psi-persistence",
        str(float(args.psi_persistence)),
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

    production_summary_path = output_dir / "summary.json"
    production_path = output_dir / "transition_path.csv"
    production_summary = json.loads(production_summary_path.read_text(encoding="utf-8"))
    rows = oracle.load_path_rows(production_path)
    if not context.evaluation_records:
        raise RuntimeError("Conditional solver recorded no path evaluation")
    final_record = context.evaluation_records[-1]
    iterations = int(production_summary["iterations"])
    if len(context.evaluation_records) not in {iterations, iterations + 1}:
        raise RuntimeError("Unexpected number of conditional path evaluations")
    accounting_gates = production_summary.get("accounting_gates") or {}
    terminal_gates = production_summary.get("terminal_convergence") or {}
    production_audit = audit_production_packet(
        summary=production_summary,
        rows=rows,
        final_record=final_record,
        expected_horizon=int(args.horizon),
    )
    status = (
        "complete_conditional_path_converged"
        if bool(production_summary.get("path_converged"))
        and all(accounting_gates.values())
        and bool(terminal_gates.get("all_checks_pass"))
        else "complete_conditional_path_unconverged"
    )
    payload = {
        "status": status,
        "interpretation": (
            "A full price/transfer transition solved conditional on one fixed tenure "
            "share. Complementarity also requires the reported final value-gap condition."
        ),
        "case": "rebated-tax2-reform",
        "horizon": int(args.horizon),
        "atom": atom.__dict__,
        "owner_share": float(args.owner_share),
        "initial_path": str(args.initial_path.resolve()),
        "initial_path_sha256": file_sha256(args.initial_path.resolve()),
        "terminal_summary": str(args.terminal_summary.resolve()),
        "terminal_summary_sha256": file_sha256(args.terminal_summary.resolve()),
        "production_summary_sha256": file_sha256(production_summary_path),
        "production_transition_path_sha256": file_sha256(production_path),
        "production_path_status": production_summary.get("path_status"),
        "production_path_converged": production_summary.get("path_converged"),
        "production_accounting_gates": accounting_gates,
        "production_terminal_convergence": terminal_gates,
        "history_reproduction_status": production_summary.get(
            "history_reproduction_status"
        ),
        "conditional_evaluations": len(context.evaluation_records),
        "final_conditional_evaluation": final_record,
        "production_audit": production_audit,
        "source_hashes": source_hashes(),
        "final_state_table_status": (
            "required_separate_exact_oracle_audit_after_conditional_convergence"
        ),
    }
    write_json_atomic(output_dir / "conditional_transition_summary.json", payload)


if __name__ == "__main__":
    main()
