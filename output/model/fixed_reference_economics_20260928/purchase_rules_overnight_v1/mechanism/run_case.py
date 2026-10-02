"""One bounded, source-authenticated dated financing experiment.

Run separately for each purchase rule, policy duration, and 12/16-date horizon.
Control uses the 80% stationary path; permanent uses a fresh 100% terminal
stationary price/renewal root with the fitted 80% H0 held fixed.  At each dated
path, physical population and both birth-to-entry queues evolve endogenously.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import re
import sys
import time
from pathlib import Path
import numpy as np

from selected_runtime import PACKET, ROOT, construct, read, require, sha

TRANSITION = ROOT / "code/model/experiments/transition_readiness"
sys.path.insert(0, str(ROOT / "code/model/tools"))
sys.path.insert(0, str(TRANSITION / "pinned_tools"))
from integration import map_case
from e5f_ssj_scaled_step_root import solve_price_path_scaled
from e5f_four_shock_acceleration import solve_joint_with_acceleration, extend_measured_jacobian

CONTROLS = ROOT / "output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json"


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    def plain(item):
        if isinstance(item, Path): return str(item)
        if isinstance(item, np.ndarray): return item.tolist()
        if isinstance(item, dict): return {str(k): plain(v) for k,v in item.items()}
        if isinstance(item, (list, tuple)): return [plain(v) for v in item]
        if isinstance(item, np.generic): return plain(item.item())
        if isinstance(item, float) and not math.isfinite(item): return None
        return item
    temp.write_text(json.dumps(plain(value), indent=2, sort_keys=True, allow_nan=False) + "\n")
    temp.replace(path)


def stationary_valid(record, gates):
    return (all(record["gates"].values())
            and abs(record["renewal_residual"]) <= gates["stationary_renewal_tolerance"]
            and math.isfinite(record["population_scale"]) and record["population_scale"] > 0)


def solve_terminal(runtime, folder, plan, deadline):
    """Approved closed-demography terminal contract at phi=1, fixed H0."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    c, gates = plan["endpoint"], plan["gates"]
    q0 = runtime.reference_price
    old_phi = runtime.P.phi.copy()
    latest = {}
    count = 0
    try:
        runtime.P.phi = np.ones_like(old_phi)
        def evaluate(q):
            nonlocal count
            require(time.monotonic() < deadline, "Terminal stage deadline reached")
            count += 1
            packet, record = runtime.stationary(runtime.P.psi_child, float(q[0]), folder / f"point_{count:03d}")
            latest.update(packet=packet, record=record)
            write(folder / "latest_completed.json", record)
            if stationary_valid(record, gates):
                write(folder / "best_so_far.json", record)
            admissible = (record["accounting_valid"] is True
                          and math.isfinite(record["renewal_residual"])
                          and math.isfinite(record["population_scale"])
                          and record["population_scale"] > 0)
            return dict(mapping_valid=admissible,
                        residual=np.array([record["renewal_residual"]]))
        root = solve_price_path_scaled(
            initial_prices=np.array([q0]), evaluate=evaluate,
            project=lambda q: np.clip(q, q0*c["price_bound_ratios"][0], q0*c["price_bound_ratios"][1]),
            slope=c["slope"], market_tolerance=gates["stationary_renewal_tolerance"],
            max_log_step=c["max_log_step"], damping=c["damping"],
            max_evaluations=c["max_evaluations"], deadline_monotonic=deadline,
            max_condition_number=plan["fit"]["max_condition_number"],
            worsening_factor=plan["fit"]["worsening_factor"],
            final_reproduction_tolerance=gates["final_reproduction_tolerance"],
            callback=lambda row: write(folder / "root_progress.json", row))
        write(folder / "root.json", root)
        require(root["converged"],
                "Permanent terminal renewal root or replay failed")
        packet, record = latest["packet"], latest["record"]
        require(stationary_valid(record, gates), "Permanent stationary gate failed")
        require(abs(float(root["final"]["prices"][0])-float(record["price"])) <= 1e-12,
                "Terminal root and saved native packet differ")
        require(float(np.asarray(packet["parameters"].phi).reshape(-1)[0]) == 1.0,
                "Permanent terminal engine did not bind phi=1")
        endpoint = dict(price=record["price"], population_scale=record["population_scale"],
                        phi=1.0, H0=float(runtime.P.H0[0]),
                        stationary_renewal_gap=abs(record["renewal_residual"]))
    finally:
        runtime.P.phi = old_phi
    state = runtime.stationary_state(packet, endpoint["population_scale"])
    native, check = map_case(runtime, kind="permanent", terminal=packet, endpoint=endpoint,
        prices=[endpoint["price"]], pensions=[packet["parameters"].pension],
        folder=folder / "one_step", initial_state=state)
    terminal = runtime.terminal_checks(packet, endpoint, native,
        np.array([runtime.P.psi_child]), tolerance=1e-6, raw_queue_tolerance=1e-6)
    require(terminal["all_checks_pass"] and all(check["gates"].values())
            and max(map(abs, check["market_residual"])) <= gates["market_tolerance"]
            and max(map(abs, check["fiscal_residual"])) <= 1e-6,
            "Permanent terminal one-step native gate failed")
    write(folder / "accepted.json", dict(endpoint=endpoint, terminal=terminal,
                                          root=folder / "root.json", policy_calls=runtime.total_native_calls))
    return packet, endpoint


def measure_control_jacobian(runtime, horizon, folder, plan, deadline):
    from run_e5f_preference_transition import measure_jacobian
    folder = Path(folder)
    endpoint = dict(price=runtime.reference_price, population_scale=runtime.population_scale)
    q = np.full(horizon, runtime.reference_price)
    b = np.full(horizon, runtime.P.pension)
    count = 0
    def evaluate(prices, pensions):
        nonlocal count
        require(time.monotonic() < deadline, "Reference derivative deadline reached")
        count += 1
        native, record = map_case(runtime, kind="control", terminal=runtime.packet,
            endpoint=endpoint, prices=prices, pensions=pensions,
            folder=folder / f"mapping_{count:03d}")
        valid = all(record["gates"].values())
        if count == 1:
            terminal = runtime.terminal_checks(runtime.packet, endpoint, native,
                np.full(horizon, runtime.P.psi_child), tolerance=1e-6,
                raw_queue_tolerance=1e-6)
            valid = (valid and terminal["all_checks_pass"]
                     and max(map(abs, record["market_residual"])) <= plan["gates"]["market_tolerance"]
                     and max(map(abs, record["fiscal_residual"])) <= 1e-6)
            write(folder / "baseline_checks.json", dict(valid=valid, terminal=terminal))
        return dict(mapping_valid=valid, market_residual=record["market_residual"],
                    fiscal_residual=record["fiscal_residual"])
    J = measure_jacobian(evaluate, q, b, 0, 1e-5, folder / "measured",
                         dict(schema="purchase_phi_measured_reference_v1", identity=runtime.identity()))
    require(J.shape == (2*horizon, 2*horizon) and np.isfinite(J).all(),
            "Fresh control Jacobian invalid")
    return J, read(folder / "measured/receipt.json")


def solve_path(runtime, kind, horizon, terminal, endpoint, J, folder, plan, deadline):
    """Original fixed-tax two-block dated price/pension root and all gates."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    q0, b0 = runtime.reference_price, float(runtime.P.pension)
    qT, bT = endpoint["price"], float(terminal["parameters"].pension)
    p, gates = plan["path"], plan["gates"]
    slopes = [float(np.median(np.abs(np.diag(J)[i*horizon:(i+1)*horizon]))) for i in (0,1)]
    require(all(math.isfinite(x) and x > 0 for x in slopes), "Fresh native root slopes absent")
    prices = np.linspace(q0, qT, horizon)
    pensions = np.linspace(b0, bT, horizon)
    count, latest = 0, {}
    def evaluate(q, b):
        nonlocal count
        require(time.monotonic() < deadline, "Dated root stage deadline reached")
        count += 1
        native, record = map_case(runtime, kind=kind, terminal=terminal, endpoint=endpoint,
                                  prices=q, pensions=b, folder=folder / f"mapping_{count:03d}")
        latest.update(native=native, record=record, folder=str(folder / f"mapping_{count:03d}"))
        write(folder / "latest_completed.json", record)
        if all(record["gates"].values()): write(folder / "best_so_far.json", record)
        return dict(mapping_valid=all(record["gates"].values()),
                    market_residual=record["market_residual"], fiscal_residual=record["fiscal_residual"])
    root = solve_joint_with_acceleration(
        closure="fixed_tax", initial_prices=prices, initial_fiscal_values=pensions,
        evaluate=evaluate,
        project_prices=lambda values: np.clip(values, q0*p["price_bound_ratios"][0],
                                               q0*p["price_bound_ratios"][1]),
        fiscal_bounds=[b0*x for x in p["pension_bound_ratios"]],
        market_tolerance=gates["market_tolerance"], fiscal_tolerance=gates["fiscal_tolerance"],
        market_slope=slopes[0], fiscal_slope=slopes[1], max_log_step=p["max_log_step"],
        damping=p["damping"], max_evaluations=p["max_evaluations"],
        deadline_monotonic=deadline, max_condition_number=plan["fit"]["max_condition_number"],
        worsening_factor=plan["fit"]["worsening_factor"],
        final_reproduction_tolerance=gates["final_reproduction_tolerance"],
        initial_jacobian=J, callback=lambda row: write(folder / "root_progress.json", row))
    write(folder / "root.json", root)
    require(root["converged"] and all(root["gates"].values()), "Dated market/fiscal root or replay failed")
    native, record = latest["native"], latest["record"]
    final = root["final"]
    require(np.allclose([row["asset_price"] for row in record["rows"]],
                            final["prices"], rtol=0, atol=1e-12)
            and np.allclose([row["pension_period_units"] for row in record["rows"]],
                            final["fiscal_values"], rtol=0, atol=1e-12),
            "Dated root and saved native path differ")
    terminal_check = runtime.terminal_checks(terminal, endpoint, native,
        np.full(horizon, runtime.P.psi_child), tolerance=gates["terminal_tolerance"],
        raw_queue_tolerance=gates["raw_queue_relative_tolerance"])
    require(terminal_check["all_checks_pass"], "Dated terminal-state gate failed")
    first_births = [dict(period=int(row["period"]), calendar_year=int(row["calendar_year"]),
                         flow=float(np.sum(row["birth_flow_first"])),
                         at_risk_mass=float(np.sum(row["childless_at_risk_mass"])),
                         age_flow=np.asarray(row["birth_flow_first"]).tolist(),
                         age_hazard=np.asarray(row["first_birth_hazard"]).tolist())
                    for row in record["fertility"]]
    dated_state = native.dated_states[4]
    state_path = folder / "accepted_period4_state.pkl.gz"
    runtime.scaffold.dump_checkpoint(state_path, dated_state)
    result = dict(status="passed", kind=kind, horizon=horizon, arm=runtime.arm,
                  phi_path=record["phi_path"], reference_identity=runtime.identity(),
                  endpoint=endpoint, root=root, terminal=terminal_check,
                  first_births=first_births, rows=record["rows"],
                  accepted_mapping=latest["folder"],
                  horizon_overlap_period4_state=dict(path=str(state_path),sha256=sha(state_path),
                                                     period=4,calendar_year=2023,
                                                     role="horizon diagnostic; not historical fit"),
                  policy_calls=runtime.total_native_calls)
    write(folder / "completed.json", result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=("hard", "quarter"), required=True)
    ap.add_argument("--selected-json", type=Path, required=True)
    ap.add_argument("--selected-completed", type=Path, required=True)
    ap.add_argument("--kind", choices=("control", "temporary", "permanent"), required=True)
    ap.add_argument("--horizon", type=int, choices=(1,12,16), required=True)
    ap.add_argument("--smoke-one-date", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--deadline-epoch", type=float, required=True)
    ap.add_argument("--maximum-policy-calls", type=int, required=True)
    args = ap.parse_args()
    require((args.smoke_one_date and args.horizon == 1 and args.kind == "control")
            or (not args.smoke_one_date and args.horizon in (12, 16)),
            "One-date smoke must be a one-period control; production horizons are 12 or 16")
    require(args.maximum_policy_calls > 0 and args.deadline_epoch > time.time()+60,
            "Finite positive policy-call and wall budgets required")
    selected = read(args.selected_json)
    purchase_plan = read(PACKET / "plan.json")
    chain_match = re.fullmatch(r"chain_?(\d+)", args.selected_completed.resolve().parent.parent.name)
    require(selected.get("status") == "postchecked" and selected.get("arm") == args.arm
            and selected.get("target_fingerprint") == purchase_plan["target_fingerprint"]
            and selected.get("weight_fingerprint") == purchase_plan["weight_fingerprint"]
            and chain_match is not None and int(chain_match.group(1)) == int(selected["chain"]),
            "Winner selection or target/weight fingerprint differs")
    selected_receipt = read(args.selected_completed)
    require(selected_receipt.get("status") == "selected_numerically_verified"
            and selected_receipt.get("selected",{}).get("parameters") == selected.get("selected_parameters"),
            "Postchecked selected coordinates differ from native winner")
    plan = read(CONTROLS)
    require(plan["gates"]["stationary_renewal_tolerance"] == 1e-6
            and plan["gates"]["market_tolerance"] == 2e-4
            and plan["gates"]["fiscal_tolerance"] == 2e-5,
            "Approved numerical gates changed")
    start = time.monotonic()
    deadline = start + args.deadline_epoch-time.time()
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / "run_contract.json", dict(arm=args.arm, kind=args.kind, horizon=args.horizon,
        deadline_epoch=args.deadline_epoch, maximum_policy_calls=args.maximum_policy_calls,
        smoke_one_date=args.smoke_one_date,
        selected_json=str(args.selected_json), selected_json_sha256=sha(args.selected_json),
        selected_completed=str(args.selected_completed), selected_sha256=sha(args.selected_completed),
        numerical_controls=str(CONTROLS), numerical_controls_sha256=sha(CONTROLS),
        economic_change="buyer financed share 0.80 to 1.00; all other fitted primitives fixed",
        fixed_H0=True, physical_population="endogenous closed birth-entry queues",
        fiscal_closure="fixed payroll tax; period pension root", transfer_path="zero"))
    try:
        runtime = construct(args.arm, args.selected_completed, args.out / "runtime")
        runtime.arm = args.arm
        with runtime.native_budget(deadline, args.maximum_policy_calls):
            runtime.reconstruct_reference(args.out / "reference")
            terminal = runtime.packet
            endpoint = dict(price=runtime.reference_price, population_scale=runtime.population_scale,
                            phi=0.8, H0=float(runtime.P.H0[0]))
            if args.kind == "permanent":
                terminal, endpoint = solve_terminal(runtime, args.out / "permanent_terminal", plan, deadline)
            if args.kind == "control":
                prices = np.full(args.horizon, runtime.reference_price)
                pensions = np.full(args.horizon, runtime.P.pension)
                native, record = map_case(runtime, kind="control", terminal=terminal,
                    endpoint=endpoint, prices=prices, pensions=pensions, folder=args.out / "control")
                terminal_check = runtime.terminal_checks(terminal, endpoint, native,
                    np.full(args.horizon, runtime.P.psi_child), tolerance=1e-6,
                    raw_queue_tolerance=1e-6)
                require(terminal_check["all_checks_pass"] and all(record["gates"].values())
                        and max(map(abs, record["market_residual"])) <= plan["gates"]["market_tolerance"]
                        and max(map(abs, record["fiscal_residual"])) <= 1e-6,
                        "Control native constant-path gate failed")
                result = dict(status="passed", arm=args.arm,
                    kind="smoke_one_date" if args.smoke_one_date else "control",
                    horizon=args.horizon, phi_path=record["phi_path"],
                    fertility=record["fertility"], rows=record["rows"], terminal=terminal_check,
                    accepted_mapping=str(args.out / "control"),
                    reference_identity=runtime.identity(), policy_calls=runtime.total_native_calls)
                if not args.smoke_one_date:
                    state_path = args.out / "accepted_period4_state.pkl.gz"
                    runtime.scaffold.dump_checkpoint(state_path, native.dated_states[4])
                    result["horizon_overlap_period4_state"] = dict(
                        path=str(state_path),sha256=sha(state_path),period=4,calendar_year=2023,
                        role="horizon diagnostic; not historical fit")
                write(args.out / "completed.json", result)
            else:
                J, derivative_receipt = measure_control_jacobian(runtime, 12, args.out / "derivative", plan, deadline)
                if args.horizon == 16: J = extend_measured_jacobian(derivative_receipt, 16)
                solve_path(runtime, args.kind, args.horizon, terminal, endpoint,
                           J, args.out / "dated_path", plan, deadline)
        write(args.out / "terminal_receipt.json", dict(status="passed", elapsed_seconds=time.monotonic()-start,
            actual_native_policy_calls=runtime.total_native_calls))
    except BaseException as exc:
        write(args.out / "failure.json", dict(status="failed", type=type(exc).__name__,
            message=str(exc), elapsed_seconds=time.monotonic()-start))
        raise


if __name__ == "__main__": main()
