#!/usr/bin/env python3
"""Run one immutable H128 path and expose deterministic tenure value rankings.

This is a default-off diagnostic wrapper around the accepted person-demography
policy driver.  It leaves the production Bellman and forward law unchanged,
captures the conditional tenure values already passed to the exact compiled
tenure kernel, and writes positive-mass decision states at selected calendar
dates.  The wrapper is deliberately restricted to the active one-market model.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Sequence

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


NEGATIVE_INFINITY = -1.0e10
MASS_FLOOR = 1.0e-15


@dataclass
class CaptureContext:
    output_dir: Path
    target_years: set[int]
    horizon: int = 0
    solve_count: int = 0
    active: bool = False
    current_period: int | None = None
    current_year: int | None = None
    current_parameters: SimpleNamespace | None = None
    current_b_grid: np.ndarray | None = None
    current_nz: int = 0
    current_slices: dict[tuple[int, int], dict[str, np.ndarray | float]] = field(
        default_factory=dict
    )
    kernel_calls: int = 0
    validations: list[dict[str, Any]] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terminal-summary", type=Path, required=True)
    parser.add_argument("--initial-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--horizon", type=int, default=128)
    parser.add_argument(
        "--diagnostic-years",
        default="2351,2363,2367",
        help="Comma-separated calendar years on the four-year model grid.",
    )
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
            f"Diagnostic years must lie on the {pf.CALENDAR_START_YEAR}-origin four-year grid "
            f"inside the requested horizon; got {sorted(years)}"
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
    if path.suffix == ".gz":
        handle = gzip.open(path, "wt", newline="", encoding="utf-8")
    else:
        handle = path.open("w", newline="", encoding="utf-8")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interpolate_with_clip(grid: np.ndarray, values: np.ndarray, query: float) -> float:
    bg = np.asarray(grid, dtype=float).reshape(-1)
    vv = np.asarray(values, dtype=float).reshape(-1)
    x = float(query)
    if x <= float(bg[0]):
        return float(vv[0])
    if x >= float(bg[-1]):
        return float(vv[-1])
    lower = int(np.searchsorted(bg, x, side="right") - 1)
    weight = (x - float(bg[lower])) / float(bg[lower + 1] - bg[lower])
    weight = min(max(weight, 0.0), 1.0)
    # Preserve the arithmetic order used by the compiled tenure kernel.
    return float((1.0 - weight) * vv[lower] + weight * vv[lower + 1])


def bellman_slice_indices(call: int, number_of_ages: int, number_of_income_states: int) -> tuple[int, int]:
    """Map the tenure-kernel call order back to Bellman's descending age loop."""

    age = int(number_of_ages) - 1 - int(call) // int(number_of_income_states)
    income_state = int(call) % int(number_of_income_states)
    if age < 0 or age >= int(number_of_ages):
        raise IndexError(f"Tenure-kernel call {call} lies outside the Bellman age loop")
    return age, income_state


def reconstruct_tenure_option_values(
    Vd: np.ndarray,
    b_grid: np.ndarray,
    heq: np.ndarray,
    hcost: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    birth_dp: np.ndarray,
    birth_entry_grant: np.ndarray,
) -> np.ndarray:
    """Mirror ``tenure_choice_kernel`` and return every conditional option value."""

    conditional = np.asarray(Vd, dtype=float)
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    Nb, nt, I, npar, ncs = conditional.shape
    options = np.full((Nb, nt, I, npar, ncs, nt), NEGATIVE_INFINITY)
    for market in range(I):
        for origin in range(nt):
            sale = float(heq[market, origin]) if origin > 0 else 0.0
            for parity in range(npar):
                for child_state in range(ncs):
                    for bb, wealth in enumerate(bg):
                        if origin == 0:
                            options[bb, origin, market, parity, child_state, 0] = conditional[
                                bb, 0, market, parity, child_state
                            ]
                        else:
                            options[bb, origin, market, parity, child_state, 0] = (
                                interpolate_with_clip(
                                    bg,
                                    conditional[:, 0, market, parity, child_state],
                                    wealth + sale,
                                )
                            )
                        for target in range(1, nt):
                            house_cost = float(hcost[market, target])
                            downpayment = float(dp_arr[market, target, parity, child_state])
                            borrowing_floor = float(bmo[market, target, parity, child_state])
                            if origin == target:
                                value = conditional[bb, target, market, parity, child_state]
                            elif origin == 0:
                                post_transaction = float(wealth - house_cost)
                                if bool(birth_dp[parity, child_state, origin, target]):
                                    value = interpolate_with_clip(
                                        bg,
                                        conditional[:, target, market, parity, child_state],
                                        max(post_transaction, borrowing_floor),
                                    )
                                elif float(
                                    birth_entry_grant[market, target, parity, child_state]
                                ) > 0.0:
                                    grant = float(
                                        birth_entry_grant[market, target, parity, child_state]
                                    )
                                    post_grant = post_transaction + grant
                                    value = interpolate_with_clip(
                                        bg,
                                        conditional[:, target, market, parity, child_state],
                                        post_grant,
                                    )
                                    if (
                                        wealth + grant < downpayment
                                        or post_grant < borrowing_floor
                                    ):
                                        value = NEGATIVE_INFINITY
                                else:
                                    value = interpolate_with_clip(
                                        bg,
                                        conditional[:, target, market, parity, child_state],
                                        post_transaction,
                                    )
                                    if (
                                        wealth < downpayment
                                        or post_transaction < borrowing_floor
                                    ):
                                        value = NEGATIVE_INFINITY
                            else:
                                post_transaction = float(wealth + sale - house_cost)
                                value = interpolate_with_clip(
                                    bg,
                                    conditional[:, target, market, parity, child_state],
                                    post_transaction,
                                )
                                if (
                                    wealth < downpayment - sale
                                    or post_transaction < borrowing_floor
                                ):
                                    value = NEGATIVE_INFINITY
                            options[bb, origin, market, parity, child_state, target] = value
    return options


def branch_details(
    *,
    wealth: float,
    origin: int,
    target: int,
    market: int,
    parity: int,
    child_state: int,
    heq: np.ndarray,
    hcost: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    birth_dp: np.ndarray,
    birth_entry_grant: np.ndarray,
) -> tuple[float, bool, str]:
    if target == origin:
        return float(wealth), True, "same_tenure"
    sale = float(heq[market, origin]) if origin > 0 else 0.0
    if target == 0:
        return float(wealth + sale), True, "rent_or_sell"
    house_cost = float(hcost[market, target])
    downpayment = float(dp_arr[market, target, parity, child_state])
    borrowing_floor = float(bmo[market, target, parity, child_state])
    if origin == 0:
        post = float(wealth - house_cost)
        if bool(birth_dp[parity, child_state, origin, target]):
            return max(post, borrowing_floor), True, "birth_downpayment_waiver"
        grant = float(birth_entry_grant[market, target, parity, child_state])
        if grant > 0.0:
            post += grant
            feasible = wealth + grant >= downpayment and post >= borrowing_floor
            return post, bool(feasible), "grant" if feasible else "grant_constraint"
        feasible = wealth >= downpayment and post >= borrowing_floor
        return post, bool(feasible), "purchase" if feasible else "purchase_constraint"
    post = float(wealth + sale - house_cost)
    feasible = wealth >= downpayment - sale and post >= borrowing_floor
    return post, bool(feasible), "rebuy" if feasible else "rebuy_constraint"


def install_capture(context: CaptureContext) -> Callable[[], None]:
    original_evaluate_path = person_pf.evaluate_path_at_prices_person_demography
    original_solve_date = pf.solve_date_policy
    original_evaluate_period = calendar.evaluate_period
    active_kernel: dict[str, Any] = {"model": None, "original": None}

    def wrapped_tenure_kernel(*args: Any) -> tuple[np.ndarray, np.ndarray]:
        original_tenure_kernel = active_kernel["original"]
        if original_tenure_kernel is None:
            raise RuntimeError("Tenure diagnostic kernel wrapper is not active")
        VH, choices = original_tenure_kernel(*args)
        if context.active and context.current_period is not None:
            Vd, bg, heq, hcost, dp_arr, bmo, birth_dp, grant = args
            parameters = context.current_parameters
            if parameters is None:
                raise RuntimeError("Tenure diagnostic lacks current parameters")
            Nz = int(context.current_nz)
            age, income_state = bellman_slice_indices(
                context.kernel_calls, int(parameters.J), Nz
            )
            options = reconstruct_tenure_option_values(
                Vd, bg, heq, hcost, dp_arr, bmo, birth_dp, grant
            )
            reconstructed_choices = np.argmax(options, axis=-1).astype(np.int16)
            mismatch_count = int(np.count_nonzero(reconstructed_choices != choices))
            reconstructed_maximum = np.max(options, axis=-1)
            absolute_error = np.abs(reconstructed_maximum - VH)
            live = np.maximum(reconstructed_maximum, VH) > NEGATIVE_INFINITY / 2.0
            live_value_error = (
                float(np.max(absolute_error[live])) if np.any(live) else 0.0
            )
            full_value_error = float(np.max(absolute_error))
            if (
                mismatch_count != 0
                or live_value_error > 2.0e-11
                or full_value_error > 2.0e-6
            ):
                raise RuntimeError(
                    "Tenure option reconstruction differs from the compiled kernel: "
                    f"mismatches={mismatch_count}, "
                    f"max_live_value_error={live_value_error:.6g}, "
                    f"max_full_value_error={full_value_error:.6g}"
                )
            context.current_slices[(age, income_state)] = {
                "options": options,
                "choices": np.asarray(choices, dtype=np.int16).copy(),
                "heq": np.asarray(heq, dtype=float).copy(),
                "hcost": np.asarray(hcost, dtype=float).copy(),
                "dp_arr": np.asarray(dp_arr, dtype=float).copy(),
                "bmo": np.asarray(bmo, dtype=float).copy(),
                "birth_dp": np.asarray(birth_dp, dtype=bool).copy(),
                "grant": np.asarray(grant, dtype=float).copy(),
                "live_value_error": live_value_error,
                "full_value_error": full_value_error,
            }
            context.kernel_calls += 1
        return VH, choices

    def wrapped_solve_date_policy(**kwargs: Any) -> calendar.PolicyBundle:
        if not context.active:
            return original_solve_date(**kwargs)
        call = context.solve_count
        context.solve_count += 1
        if call >= context.horizon:
            period = call - context.horizon
            year = pf.CALENDAR_START_YEAR + 4 * period
            if year in context.target_years:
                context.current_period = period
                context.current_year = year
                context.current_parameters = kwargs["P"]
                context.current_nz = len(calendar.model.income_transition_values(kwargs["P"])[0])
                context.current_b_grid = np.asarray(kwargs["b_grid"], dtype=float)
                context.current_slices = {}
                context.kernel_calls = 0
        policy = original_solve_date(**kwargs)
        if context.current_period is not None:
            expected = int(context.current_parameters.J) * int(context.current_nz)
            if context.kernel_calls != expected:
                raise RuntimeError(
                    f"Captured {context.kernel_calls} tenure slices; expected {expected}"
                )
        return policy

    def wrapped_evaluate_period(*args: Any, **kwargs: Any) -> calendar.PeriodEvaluation:
        evaluation = original_evaluate_period(*args, **kwargs)
        if context.current_period is not None:
            write_year_diagnostic(context, evaluation, kwargs.get("P", args[2]))
            context.current_period = None
            context.current_year = None
            context.current_parameters = None
            context.current_b_grid = None
            context.current_nz = 0
            context.current_slices = {}
            context.kernel_calls = 0
        return evaluation

    def wrapped_evaluate_path(*args: Any, **kwargs: Any) -> person_pf.PersonPathEvaluation:
        prices = np.asarray(kwargs.get("prices", args[0] if args else []), dtype=float)
        context.horizon = int(prices.size)
        context.solve_count = 0
        context.active = True
        active_model = calendar.model
        original_tenure_kernel = active_model.tenure_choice_kernel
        active_kernel["model"] = active_model
        active_kernel["original"] = original_tenure_kernel
        active_model.tenure_choice_kernel = wrapped_tenure_kernel
        try:
            result = original_evaluate_path(*args, **kwargs)
        finally:
            active_model.tenure_choice_kernel = original_tenure_kernel
            active_kernel["model"] = None
            active_kernel["original"] = None
            context.active = False
        expected_solves = 2 * context.horizon
        if context.solve_count != expected_solves:
            raise RuntimeError(
                f"Diagnostic observed {context.solve_count} date-policy solves; "
                f"expected {expected_solves}"
            )
        return result

    pf.solve_date_policy = wrapped_solve_date_policy
    calendar.evaluate_period = wrapped_evaluate_period
    person_pf.evaluate_path_at_prices_person_demography = wrapped_evaluate_path

    def restore() -> None:
        if active_kernel["model"] is not None and active_kernel["original"] is not None:
            active_kernel["model"].tenure_choice_kernel = active_kernel["original"]
        pf.solve_date_policy = original_solve_date
        calendar.evaluate_period = original_evaluate_period
        person_pf.evaluate_path_at_prices_person_demography = original_evaluate_path

    return restore


def write_year_diagnostic(
    context: CaptureContext,
    evaluation: calendar.PeriodEvaluation,
    parameters: SimpleNamespace,
) -> None:
    if int(parameters.I) != 1:
        raise NotImplementedError("The tenure-kink diagnostic is restricted to I=1")
    if context.current_year is None or context.current_b_grid is None:
        raise RuntimeError("Incomplete tenure diagnostic context")
    bg = context.current_b_grid
    g_post = np.asarray(evaluation.g_post_fertility, dtype=float)
    total_mass = float(np.sum(g_post))
    positive = np.argwhere(g_post > MASS_FLOOR)
    z_grid, _, _ = calendar.model.income_transition_values(parameters)
    rows: list[dict[str, Any]] = []
    max_value_error = 0.0
    mismatch_count = 0
    location_gap = 0.0
    nt = 1 + int(parameters.n_house)
    for bb, origin, market, age, income_state, parity, child_state in positive:
        mass = float(g_post[bb, origin, market, age, income_state, parity, child_state])
        loc_prob = float(
            evaluation.policy.loc_probs[
                bb, origin, market, market, age, income_state, parity, child_state
            ]
        )
        location_gap = max(location_gap, abs(loc_prob - 1.0))
        captured = context.current_slices[(int(age), int(income_state))]
        options = np.asarray(captured["options"], dtype=float)[
            bb, origin, market, parity, child_state, :
        ]
        stored = int(
            evaluation.policy.tenure_choice[
                bb, origin, market, age, income_state, parity, child_state
            ]
        )
        order = np.argsort(-options, kind="stable")
        best = int(order[0])
        second = int(order[1]) if nt > 1 else int(order[0])
        mismatch_count += int(best != stored)
        max_value_error = max(
            max_value_error, float(captured["live_value_error"])
        )
        base: dict[str, Any] = {
            "calendar_year": int(context.current_year),
            "period": int(context.current_period),
            "wealth_index": int(bb),
            "liquid_wealth": float(bg[bb]),
            "origin_tenure": int(origin),
            "market": int(market),
            "age_index": int(age),
            "age_years": float(parameters.age_start) + float(age) * float(parameters.da),
            "income_state": int(income_state),
            "income_multiplier": float(z_grid[income_state]),
            "parity": int(parity),
            "child_state": int(child_state),
            "decision_mass": mass,
            "decision_mass_share": mass / max(total_mass, 1.0e-15),
            "stored_tenure": stored,
            "best_tenure": best,
            "second_tenure": second,
            "best_value": float(options[best]),
            "second_value": float(options[second]),
            "best_minus_second": float(options[best] - options[second]),
            "location_probability": loc_prob,
        }
        for target in range(nt):
            branch_b, feasible, reason = branch_details(
                wealth=float(bg[bb]),
                origin=int(origin),
                target=target,
                market=int(market),
                parity=int(parity),
                child_state=int(child_state),
                heq=np.asarray(captured["heq"]),
                hcost=np.asarray(captured["hcost"]),
                dp_arr=np.asarray(captured["dp_arr"]),
                bmo=np.asarray(captured["bmo"]),
                birth_dp=np.asarray(captured["birth_dp"]),
                birth_entry_grant=np.asarray(captured["grant"]),
            )
            if target == 0:
                housing = interpolate_with_clip(
                    bg,
                    evaluation.policy.hR_pol[
                        :, 0, market, age, income_state, parity, child_state
                    ],
                    branch_b,
                )
            else:
                housing = float(parameters.H_own[target - 1])
            next_wealth = interpolate_with_clip(
                bg,
                evaluation.policy.bp_pol[
                    :, target, market, age, income_state, parity, child_state
                ],
                branch_b,
            )
            tax = (
                max(float(parameters.tau_H), 0.0)
                * float(evaluation.policy.price[market])
                * housing
            )
            prefix = f"tenure_{target}_"
            base[prefix + "value"] = float(options[target])
            base[prefix + "value_minus_best"] = float(options[target] - options[best])
            base[prefix + "feasible"] = bool(feasible)
            base[prefix + "feasibility_reason"] = reason
            base[prefix + "branch_liquid_wealth"] = float(branch_b)
            base[prefix + "housing_services"] = housing
            base[prefix + "property_tax_per_head"] = tax
            base[prefix + "next_liquid_wealth"] = next_wealth
        base["selected_housing_services"] = base[
            f"tenure_{stored}_housing_services"
        ]
        base["selected_property_tax_contribution"] = mass * base[
            f"tenure_{stored}_property_tax_per_head"
        ]
        base["selected_housing_demand_contribution"] = mass * base[
            f"tenure_{stored}_housing_services"
        ]
        base["selected_owner_mass"] = mass * float(stored > 0)
        rows.append(base)
    if mismatch_count != 0 or max_value_error > 2.0e-11 or location_gap > 2.0e-12:
        raise RuntimeError(
            "Tenure state diagnostic validation failed: "
            f"choice_mismatches={mismatch_count}, value_error={max_value_error:.6g}, "
            f"location_probability_gap={location_gap:.6g}"
        )
    path = context.output_dir / f"tenure_value_states_{context.current_year}.csv.gz"
    write_csv(path, rows)
    context.validations.append(
        {
            "calendar_year": int(context.current_year),
            "positive_mass_states": len(rows),
            "decision_mass": total_mass,
            "maximum_compiled_value_error": max_value_error,
            "choice_mismatch_count": mismatch_count,
            "maximum_one_market_location_probability_gap": location_gap,
            "state_csv": str(path),
            "state_csv_sha256": sha256(path),
        }
    )


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
    context = CaptureContext(output_dir=output_dir, target_years=years)
    restore = install_capture(context)
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
        restore()
    if {int(item["calendar_year"]) for item in context.validations} != years:
        raise RuntimeError(
            "Not every requested diagnostic year was captured: "
            f"requested={sorted(years)}, captured={context.validations}"
        )
    production_summary = output_dir / "summary.json"
    if not production_summary.is_file():
        raise RuntimeError("The wrapped production driver did not write summary.json")
    payload = {
        "status": "complete_unpromoted_tenure_value_diagnostic",
        "interpretation": (
            "Default-off readout of the exact deterministic tenure kernel along one "
            "fixed H128 path; it does not change the Bellman problem or forward law."
        ),
        "case": "rebated-tax2-reform",
        "horizon": int(args.horizon),
        "diagnostic_years": sorted(years),
        "initial_path": str(args.initial_path.resolve()),
        "initial_path_sha256": sha256(args.initial_path.resolve()),
        "terminal_summary": str(args.terminal_summary.resolve()),
        "terminal_summary_sha256": sha256(args.terminal_summary.resolve()),
        "source_hashes": {
            "diagnostic_wrapper": sha256(Path(__file__).resolve()),
            "production_driver": sha256(Path(policy_driver.__file__).resolve()),
            "person_demography": sha256(Path(person_pf.__file__).resolve()),
            "perfect_foresight": sha256(Path(pf.__file__).resolve()),
            "solver": sha256(Path(calendar.model.__file__).resolve()),
        },
        "production_summary_sha256": sha256(production_summary),
        "validations": context.validations,
    }
    write_json(output_dir / "tenure_diagnostic_summary.json", payload)


if __name__ == "__main__":
    main()
