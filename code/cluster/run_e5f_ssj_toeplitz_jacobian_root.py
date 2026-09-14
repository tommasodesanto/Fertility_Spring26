"""Bounded sequence-space experiment: measured block-Toeplitz initial Jacobian.

Stage 1 measures, at the original 2007 stationary no-shock economy over a
ten-date horizon, the central finite-difference response of every dated native
residual to one log-coordinate perturbation at the middle date in each of the
three unknown blocks (log price, log pension, log rebate).  Seven exact native
mappings: one baseline plus three central pairs.  Stage 2 reruns the retained
ten-period shocked fixed-terminal root with identical start guess, controls,
endpoint and eight-mapping budget, supplying the assembled Jacobian only as the
Broyden ``initial_jacobian``.  Every economic object, gate, population law and
residual definition is the frozen native one; nothing is redefined here.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from e5f_ssj_toeplitz_jacobian import (assemble_jacobian, central_column, condition_report,
                                       diagonal_default)
from e5f_ssj_scaled_step_root import solve_price_path_scaled

HORIZON = 10
PERTURBED_DATE = 5
STEP = 1e-5
BASELINE_GATE = 2e-4
DRIFT_LIMIT = 1e-5
MAPPING_CAP = 7


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=_json_default) + "\n")
    tmp.replace(path)


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    raise TypeError(type(value).__name__)


def state_gap(actual, reference):
    actual_g = np.asarray(actual.g_pre, dtype=float)
    reference_g = np.asarray(reference.g_pre, dtype=float)
    scale = max(float(reference_g.sum()), 1e-15)
    return dict(
        distribution_relative_l1=float(np.abs(actual_g - reference_g).sum() / scale),
        population_relative_gap=float(abs(actual_g.sum() / scale - 1.0)),
        queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries, dtype=float)
                                              / np.asarray(reference.scheduled_entries, dtype=float) - 1.0))),
        raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries, dtype=float)
                                                  / np.asarray(reference.scheduled_raw_entries, dtype=float) - 1.0))))


def run_derivative_stage(mapping, out, *, horizon, perturbed_date, step, deadline_monotonic,
                         base_coordinates, slope):
    """Seven-mapping measurement loop.  ``mapping(label, u)`` returns a dict with
    ``residual`` (stacked ``(3T,)``), ``gates``, ``seconds`` and, for the baseline,
    ``stationary_drift``.  Pure orchestration so the loop is locally testable."""
    out = Path(out)
    started = time.monotonic()
    evaluations = []

    def call(label, u):
        if len(evaluations) >= MAPPING_CAP or time.monotonic() >= deadline_monotonic:
            raise TimeoutError("Derivative mapping cap or stage budget reached")
        result = mapping(label, np.asarray(u, dtype=float))
        residual = np.asarray(result["residual"], dtype=float)
        if residual.shape != (3 * horizon,) or not np.isfinite(residual).all():
            raise ValueError("Native residual must be finite with shape (3T,)")
        record = dict(label=label, coordinates=np.asarray(u).tolist(), residual=residual.tolist(),
                      max_abs=float(np.max(np.abs(residual))), gates=result["gates"],
                      seconds=float(result["seconds"]))
        if "stationary_drift" in result:
            record["stationary_drift"] = result["stationary_drift"]
        evaluations.append(record)
        save(out / "latest_completed.json", dict(status="mapping_complete", stage="derivative",
             mappings=len(evaluations), latest_label=label, latest_max_abs=record["max_abs"],
             latest_seconds=record["seconds"], elapsed_seconds=time.monotonic() - started))
        return residual

    base = np.asarray(base_coordinates, dtype=float)
    if base.shape != (3, horizon):
        raise ValueError("base coordinates must have shape (3, T)")
    f0 = call("baseline", base)
    baseline_max = float(np.max(np.abs(f0)))
    drift = evaluations[0].get("stationary_drift")
    if drift is None:
        raise RuntimeError("Baseline mapping did not report stationary drift")
    baseline_ok = baseline_max <= BASELINE_GATE and all(v <= DRIFT_LIMIT for v in drift.values())
    if not baseline_ok:
        save(out / "derivative_receipt.json", dict(status="baseline_gate_failed", baseline_max_abs=baseline_max,
             stationary_drift=drift, evaluations=evaluations, elapsed_seconds=time.monotonic() - started))
        raise RuntimeError("Stationary baseline failed the unchanged native gate; no Jacobian built")
    columns = []
    for j, name in enumerate(("log_house_price", "log_pension", "log_rebate")):
        up, down = base.copy(), base.copy()
        up[j, perturbed_date] += step
        down[j, perturbed_date] -= step
        fp = call(f"plus_{name}", up)
        fm = call(f"minus_{name}", down)
        columns.append(central_column(fp, fm, step, horizon))
    jacobian, receipt = assemble_jacobian(columns, horizon, perturbed_date)
    default = diagonal_default(horizon, slope)
    receipt.update(status="complete", baseline_max_abs=baseline_max, baseline_gate=BASELINE_GATE,
                   stationary_drift=drift, drift_limit=DRIFT_LIMIT, step=step,
                   mappings=len(evaluations), evaluations=evaluations,
                   mapping_seconds=[e["seconds"] for e in evaluations],
                   elapsed_seconds=time.monotonic() - started,
                   comparison_to_default=condition_report(jacobian, default),
                   own_date_derivatives=dict(
                       housing_wrt_log_price=float(columns[0][perturbed_date]),
                       paygo_wrt_log_pension=float(columns[1][horizon + perturbed_date]),
                       rebate_wrt_log_rebate=float(columns[2][2 * horizon + perturbed_date]),
                       rebate_wrt_log_price=float(columns[0][2 * horizon + perturbed_date]),
                       housing_wrt_log_rebate=float(columns[2][perturbed_date]),
                       housing_wrt_log_pension=float(columns[1][perturbed_date])),
                   fake_news_derivatives_constructed=False, speedup_claimed=False)
    save(out / "jacobian.json", dict(jacobian=jacobian.tolist(), default_jacobian=default.tolist()))
    save(out / "derivative_receipt.json", receipt)
    return jacobian, receipt


def select_step_rule(step_rule):
    """``clipped`` keeps the retained solver; ``scaled`` swaps in the isolated copy."""
    if step_rule == "clipped":
        return None
    if step_rule == "scaled":
        return solve_price_path_scaled
    raise ValueError("step_rule must be 'clipped' or 'scaled'")


def compare_roots(new_receipt, reference_receipt, horizon):
    def table(receipt):
        rows = []
        for e in receipt.get("history", []):
            if "evaluation" not in e:
                continue
            r = np.asarray(e["residual"], dtype=float)
            rows.append(dict(evaluation=e["evaluation"], phase=e["phase"], score=float(e["score"]),
                             max_housing=float(np.max(np.abs(r[:horizon]))),
                             max_paygo=float(np.max(np.abs(r[horizon:2 * horizon]))),
                             max_rebate=float(np.max(np.abs(r[2 * horizon:]))),
                             evaluation_seconds=float(e.get("evaluation_seconds", float("nan"))),
                             safeguard=e.get("safeguard"), reset_reason=e.get("reset_reason")))
        return rows
    new_rows, ref_rows = table(new_receipt), table(reference_receipt)
    return dict(reference_status=reference_receipt.get("status"), new_status=new_receipt.get("status"),
                reference_converged=bool(reference_receipt.get("converged")),
                new_converged=bool(new_receipt.get("converged")),
                reference_elapsed_seconds=reference_receipt.get("elapsed_seconds"),
                new_elapsed_seconds=new_receipt.get("elapsed_seconds"),
                reference_best_score=(reference_receipt.get("best") or {}).get("score"),
                new_best_score=(new_receipt.get("best") or {}).get("score"),
                reference_history=ref_rows, new_history=new_rows,
                identical_controls_start_endpoint=True, production_eligible=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    m = read(args.manifest)
    for path, digest in m["file_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Changed pinned input: " + path)
    spec = read(m["spec"])
    out = Path(m["output"])
    if (out / "experiment_contract.json").exists():
        raise ValueError("Refusing to overwrite a started experiment")
    sys.path.insert(0, str(Path(spec["batch"]) / "source"))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(m["spec"])
    c.spec_path = Path(m["spec"])
    smoke = read(spec["smoke_summary"])
    if smoke.get("status") != "passed" or smoke["spec_sha256"] != c.driver.sha(m["spec"]):
        raise ValueError("Original native smoke prerequisite is missing")
    endpoint_receipt = read(m["endpoint_receipt"])
    with gzip.open(m["endpoint_pickle"], "rb") as f:
        endpoint = pickle.load(f)
    if (not endpoint.verified or not endpoint_receipt["verified"] or endpoint.receipt != endpoint_receipt
            or endpoint.parameters.psi_child != spec["permanent_psi"]):
        raise ValueError("Verified endpoint does not match its receipt and permanent shock")
    reference = read(m["reference_root_receipt"])
    started_unix = time.time()
    end = started_unix + float(m["seconds"])
    deadline = time.monotonic() + (end - time.time())
    derivative_deadline = time.monotonic() + float(m["derivative_seconds"])
    save(out / "experiment_contract.json", dict(
        manifest_sha256=sha(args.manifest), horizon=HORIZON, perturbed_date=PERTURBED_DATE, step=STEP,
        derivative_mapping_cap=MAPPING_CAP, root_mapping_cap=8, seconds=m["seconds"],
        derivative_seconds=m["derivative_seconds"], started_unix=started_unix, deadline_unix=end,
        afternoon_spec_absolute_deadline_unix=spec["absolute_deadline_unix"],
        afternoon_spec_deadline_applied=False,
        deadline_note="Handoff-authorized two-hour diagnostic cap replaces the expired afternoon batch deadline; all pinned scientific inputs unchanged.",
        population_law=spec["population_law"], no_immigration=True, psi=spec["permanent_psi"],
        comparison_reference=m["reference_root_receipt"], production_eligible=False,
        step_rule=m.get("step_rule", "clipped"), jacobian_source=m.get("jacobian_source"),
        warm_start_receipt=m.get("warm_start_receipt"),
        fake_news_derivatives_constructed=False))
    step_solver = select_step_rule(m.get("step_rule", "clipped"))
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            save(out / "controller_heartbeat.json", dict(remaining_seconds=end - time.time()))
            if time.time() >= end:
                save(out / "controller_failure.json", dict(error="Diagnostic deadline"))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()

    rebated = c.rebated
    q = float(c.packet["evaluation"].policy.price[0])
    pension = float(c.old.parameters.pension)
    rebate = float(c.old.parameters.property_tax_lump_sum_transfer)
    stationary_terminal = NS(parameters=c.old.parameters, policy=c.packet["evaluation"].policy, asset_price=q)
    inherited = rebated.InheritedState(2007, c.old.initial_state)
    base = np.vstack((np.full(HORIZON, np.log(q)), np.full(HORIZON, np.log(pension)),
                      np.full(HORIZON, np.log(rebate))))
    slope = float(c.controls.get("market_slope", 1.63))

    def native_mapping(label, u):
        began = time.monotonic()
        result = c.queue.queue_path(inherited=inherited, old_state=c.old, prices=np.exp(u[0]),
                                    pensions=np.exp(u[1]), transfers=np.exp(u[2]),
                                    psi=float(c.old.parameters.psi_child), terminal=stationary_terminal,
                                    demographics=None)
        if len(result.values) != HORIZON + 1 or len(result.rows) != HORIZON:
            raise ValueError("Native path did not retain T+1 values and T rows")
        blocks = [rebated.dated_residual(demand=row["housing_demand"], supply=row["housing_supply"],
                                         payroll_accounts=row, tax_accounts=row) for row in result.rows]
        residual = rebated.stack_dated_residuals(blocks)
        gates = dict(maximum_mass_accounting_error=float(result.maximum_mass_accounting_error),
                     maximum_policy_reproduction_error=float(result.maximum_policy_reproduction_error),
                     maximum_feasibility_projection_mass=float(result.maximum_feasibility_projection_mass),
                     terminal_g_pre_mass=float(result.person_tail.terminal_state.g_pre.sum()))
        reply = dict(residual=residual, gates=gates, seconds=time.monotonic() - began)
        if label == "baseline":
            reply["stationary_drift"] = state_gap(result.person_tail.terminal_state, c.old.initial_state)
        save(out / "derivative" / "rows" / f"{label}.json", result.rows)
        return reply

    @contextmanager
    def capture(folder):
        native = c.rebated.evaluate_forecast
        timings = []

        def evaluate(**kwargs):
            began = time.monotonic()
            result = native(**kwargs)
            timings.append(time.monotonic() - began)
            case = Path(folder) / "last_sweep"
            c.driver.save(case / "rows.json", result.rows)
            c.driver.save(case / "mapping_timings.json", timings)
            c.driver.save(case / "root_receipt.json", dict(converged=False,
                          status="Last full mapping only; equilibrium acceptance not established"))
            return result
        with patch.object(c.rebated, "evaluate_forecast", evaluate):
            yield

    try:
        with c.queue.original_queue_adapter(), runner.original_receipts(c), \
                c.cache.policy_cache(c.joined.pf, max_bytes=6 * 1024**3):
            if m.get("jacobian_source"):
                jacobian = np.asarray(read(m["jacobian_source"])["jacobian"], dtype=float)
                if jacobian.shape != (3 * HORIZON, 3 * HORIZON) or not np.isfinite(jacobian).all():
                    raise ValueError("Reused Jacobian has the wrong shape")
                receipt = dict(elapsed_seconds=0.0, mapping_seconds=[], reused_from=m["jacobian_source"])
                save(out / "derivative" / "derivative_receipt.json", dict(status="reused", **receipt))
            else:
                jacobian, receipt = run_derivative_stage(
                    native_mapping, out / "derivative", horizon=HORIZON, perturbed_date=PERTURBED_DATE,
                    step=STEP, deadline_monotonic=min(deadline, derivative_deadline),
                    base_coordinates=base, slope=slope)
            save(out / "stage_status.json", dict(derivative="complete", root="starting",
                 elapsed_seconds=time.time() - started_unix))
            native_solve = c.rebated.solve_rebated_forecast
            warm = None
            if m.get("warm_start_receipt"):
                previous = read(m["warm_start_receipt"])
                best = previous["best"]
                if best is None or previous.get("final_reproduction_max_abs") != 0.0:
                    raise ValueError("Warm start requires an exactly reproduced previous best")
                warm = dict(coordinates=np.asarray(best["prices"], dtype=float),
                            jacobian=np.asarray(previous["final_jacobian"], dtype=float),
                            previous_best_score=float(best["score"]))
                if warm["coordinates"].shape != (3 * HORIZON,) or warm["jacobian"].shape != (3 * HORIZON, 3 * HORIZON):
                    raise ValueError("Warm-start receipt has the wrong horizon")
                save(out / "warm_start.json", dict(source=m["warm_start_receipt"], previous_best_score=warm["previous_best_score"],
                     previous_evaluations=previous.get("evaluations"), coordinates=warm["coordinates"].tolist()))

            def solve_with_jacobian(**kwargs):
                controls = dict(kwargs["root_controls"])
                if controls.get("max_evaluations") != 8:
                    raise ValueError("Comparison requires the retained eight-mapping budget")
                if warm is None:
                    controls["initial_jacobian"] = jacobian
                    return native_solve(**dict(kwargs, root_controls=controls))
                controls["initial_jacobian"] = warm["jacobian"]
                x = warm["coordinates"]
                return native_solve(**dict(kwargs, root_controls=controls, initial_prices=x[:HORIZON],
                                           initial_pensions=x[HORIZON:2 * HORIZON], initial_transfers=x[2 * HORIZON:]))
            solver_patch = (patch.object(c.rebated, "_path_root_solver", lambda: step_solver)
                            if step_solver is not None else patch.object(c.rebated, "_path_root_solver",
                                                                         c.rebated._path_root_solver))
            with solver_patch, patch.object(c.rebated, "solve_rebated_forecast", solve_with_jacobian), \
                    capture(out / "transition"):
                result = runner.fixed_terminal_path(c, endpoint, out / "transition", HORIZON,
                                                    spec["permanent_psi"], deadline)
            comparison = compare_roots(result.root_receipt, reference, HORIZON)
            comparison["derivative_stage_seconds"] = receipt["elapsed_seconds"]
            comparison["derivative_mapping_seconds"] = receipt["mapping_seconds"]
            save(out / "comparison.json", comparison)
            save(out / "controller_complete.json", dict(
                finite_root_converged=bool(result.root_receipt["finite_horizon_market_fiscal_converged"]),
                reference_converged=bool(reference.get("converged")),
                new_best_score=comparison["new_best_score"], reference_best_score=comparison["reference_best_score"],
                stationary_endpoint_verified=True, horizon_verified=False, production_eligible=False,
                elapsed_seconds=time.time() - started_unix))
        c.driver.verify_pins(m["file_sha256"])
    except BaseException as exc:
        save(out / "controller_failure.json", dict(error_type=type(exc).__name__, error=str(exc),
             elapsed_seconds=time.time() - started_unix))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
