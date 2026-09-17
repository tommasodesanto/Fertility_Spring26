"""Bounded original-household-entry IRFs and historical refits.

This isolated diagnostic keeps the calibrated household parameters and the
four-vintage births/2.1 queue throughout. It never imports the empirical
historical age reweighting or the post-2023 person/headship transition.
"""
from __future__ import annotations

import argparse
import copy
from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"
import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def load_context(spec_path):
    spec = read(spec_path)
    base = read(spec["base_spec"])
    manifest = read(base["runtime_manifest"])
    plan = read(manifest["prior_plan"])
    # The new adapter is importable, but original scientific modules remain
    # pinned to the previously verified cluster snapshot.
    adapter_directory=str(Path(__file__).parent)
    sys.path[:] = [base["helper"],
                    str(Path(plan["source_root"]) / "code/model/tools"),
                    str(Path(plan["source_root"]) / "code/model"), base["cache_directory"],
                    adapter_directory] + [p for p in sys.path if p != adapter_directory]
    import run_e5f_final_rebated_history as driver
    import e5f_rebated_surprises as rebated
    import e5f_exact_policy_cache as cache
    import e5f_original_queue_experiment as queue
    from e5f_balanced_terminal import TerminalAuditControls
    from e5f_rebated_initial_bridge import build_rebated_initial_state
    if Path(driver.__file__).parent.resolve()!=Path(base["helper"]).resolve():
        raise ValueError("Imported driver is not the pinned helper")
    allowed={Path(base["helper"]).resolve(),(Path(plan["source_root"])/"code/model/tools").resolve()}
    if Path(rebated.__file__).parent.resolve() not in allowed:
        raise ValueError("Imported rebated operator is not frozen")
    driver.verify_pins(spec["file_sha256"])
    driver.verify_pins(base["file_sha256"])
    driver.verify_pins(manifest["file_sha256"])
    driver.verify_pins(plan["file_sha256"])
    _, joined, primitive, _, _ = rebated._runtime()
    cache_proof = read(base["cache_proof"])
    if not (cache_proof["exact_mapping_equal"] and cache_proof["household_and_accounting_mapping_valid"]
            and driver.sha(cache.__file__) == base["cache_sha256"] == cache_proof["cache_sha256"]):
        raise ValueError("Exact cache prerequisite failed")
    for pair in read(manifest["kernel_equivalence"])["pairs"]:
        if driver.sha(pair["initial"]) != driver.sha(pair["history"]) or driver.sha(pair["initial"]) != pair["sha256"]:
            raise ValueError("Initial/history kernel mismatch")
    summary = read(manifest["initial_summary"])
    item = summary["checkpoint"]
    cp = item.get("path", item.get("checkpoint"))
    if driver.sha(cp) != item.get("sha256", item.get("checkpoint_sha256")):
        raise ValueError("Initial checkpoint mismatch")
    with gzip.open(cp, "rb") as stream:
        packet = pickle.load(stream)
    raw = read(base["initial_raw_summary"])
    old = build_rebated_initial_state(packet=packet, normalization=raw["normalization"],
        outside_origin_entry_share=plan["outside_origin_entry_share"], preference_change_2023=0.)
    old = queue.initialize_original(old, packet)
    controls = dict(plan["history_root_controls"])
    controls.update(manifest.get("root_controls", {}))
    controls.setdefault("transfer_bounds", [1e-10, 10.])
    audit = TerminalAuditControls(**plan["terminal_template"]["audit_controls"])
    return NS(spec=spec, base=base, manifest=manifest, plan=plan, driver=driver,
        rebated=rebated, cache=cache, queue=queue, joined=joined, primitive=primitive,
        packet=packet, old=old, controls=controls, audit=audit)


def reference(c, evaluation=None, parameters=None):
    e = c.packet["evaluation"] if evaluation is None else evaluation
    P = c.old.parameters if parameters is None else parameters
    return dict(birth_children=float(e.births), adult_population=float(e.g_current.sum()),
        housing_demand=float(e.demand_by_loc[0]), asset_price=float(e.policy.price[0]),
        renter_price=float(P.user_cost_rate) * float(e.policy.price[0]), pension_period=float(P.pension))


def plot_case(c, folder, *, label, status, rows=None):
    folder = Path(folder)
    if rows is not None:
        c.driver.save(folder / "rows.json", rows)
    if not (folder / "rows.json").exists() or not read(folder / "rows.json"):
        return
    c.driver.save(folder / "initial_reference.json", reference(c))
    c.driver.save(folder / "irf_contract.json", dict(label=label, status_label=status,
        shock_description="No preference shock." if "No-shock" in label else "Permanent preference shock.",
        demographic_rule="original_four_vintage_household_queue", birth_to_entry_conversion=1/2.1,
        historical_age_conditioning=False, person_headship_transition=False,
        shock_calendar_year=read(folder/"rows.json")[0].get("calendar_year",2007),
        presentation_results_replaced=False))
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_e5f_stationary_shock_figures.py")),
                    "--case-dir", str(folder)], check=True)


def initial_coordinates(c, count):
    return np.repeat([float(c.packet["evaluation"].policy.price[0]), float(c.old.parameters.pension),
                      float(c.old.parameters.property_tax_lump_sum_transfer)], count+1)


@contextmanager
def capture_paths(c, folder):
    """Keep a readable last sweep even if a later root evaluation times out."""
    native=c.rebated.evaluate_forecast
    folder=Path(folder)
    def evaluate(**kwargs):
        result=native(**kwargs)
        c.driver.save(folder/"last_evaluated_rows.json",result.rows)
        diagnostic=folder/"last_sweep"
        plot_case(c,diagnostic,label=("No-shock: last completed sweep" if kwargs["psi"]==float(c.old.parameters.psi_child)
                                    else "Last completed original-rule sweep"),
            status="Price-path evaluation only; equilibrium acceptance not established.",rows=result.rows)
        return result
    with patch.object(c.rebated,"evaluate_forecast",evaluate):
        yield


@contextmanager
def capture_root_paths(c):
    native=c.driver.solve_forecast
    def solve(**kwargs):
        with capture_paths(c,kwargs["folder"]):
            return native(**kwargs)
    with patch.object(c.driver,"solve_forecast",solve):
        yield


@contextmanager
def original_receipts(c):
    native=c.driver.save
    names={"contract_receipt.json","root_receipt.json","finite_history_complete.json","summary.json"}
    def save(path,value):
        if Path(path).name in names and isinstance(value,dict):
            value.update(c.queue.annotate_original_queue_metadata())
            if value.get("case") in ("A0","A+"):
                value["case_is_internal_compatibility_label"]=True
        return native(path,value)
    with patch.object(c.driver,"save",save):yield


def finite(c, folder, count, psi, deadline, *, max_evaluations=24):
    controls = dict(c.controls, max_evaluations=max_evaluations)
    inherited = c.rebated.InheritedState(2007, c.old.initial_state)
    result, detail = c.driver.solve_forecast(inherited=inherited, old=c.old, demographics=None,
        psi=float(psi), count=count, initial=initial_coordinates(c, count), controls=controls,
        audit=c.audit, deadline=deadline, folder=folder, case="A0")
    result.root_receipt.update(demographic_rule="original_four_vintage_household_queue",
        original_case_label_internal_only=True, historical_age_conditioning=False,
        person_headship_transition=False, permanent_shock=bool(psi!=float(c.old.parameters.psi_child)))
    c.driver.save(Path(folder)/"root_receipt.json", result.root_receipt)
    if detail.get("snapshot"):
        from run_e5f_successive_surprises_overnight import standard_graphs
        standard_graphs(detail["snapshot"], result, Path(folder)/"graphs")
    plot_case(c, folder, label=(f"No-shock replay: {count} periods" if psi==float(c.old.parameters.psi_child)
                               else f"Permanent preference decline: {count} periods"),
        status=("Finite markets and budgets converged; terminal steady state unverified."
                if result.next_state is not None else "Unconverged last evaluated path; diagnostic only."))
    return result, detail


def original_no_shock_smoke(c, folder, deadline):
    P, pf = c.old.parameters, c.joined.pf
    q = float(c.packet["evaluation"].policy.price[0])
    terminal = NS(parameters=P, policy=c.packet["evaluation"].policy, asset_price=q)
    inherited = c.rebated.InheritedState(2007, c.old.initial_state)
    kwargs = dict(inherited=inherited, old_state=c.old, prices=[q]*6,
        pensions=[float(P.pension)]*6, transfers=[float(P.property_tax_lump_sum_transfer)]*6,
        psi=float(P.psi_child), terminal=terminal, demographics=None)
    first = c.queue.queue_path(**kwargs)
    second = c.queue.queue_path(**kwargs)
    if c.driver.clean(first.rows) != c.driver.clean(second.rows):
        raise ValueError("Original-law path replay differs")
    np.testing.assert_array_equal(first.person_tail.terminal_state.g_pre, second.person_tail.terminal_state.g_pre)
    for a, b in zip(first.values, second.values):
        np.testing.assert_array_equal(a, b)
    state = first.person_tail.terminal_state
    initial = c.old.initial_state
    mass = float(initial.g_pre.sum())
    drift = float(np.abs(state.g_pre-initial.g_pre).sum())/mass
    queue_drift = float(np.max(np.abs(np.asarray(state.scheduled_entries)/np.asarray(initial.scheduled_entries)-1)))
    raw_queue_drift = float(np.max(np.abs(np.asarray(state.scheduled_raw_entries)/np.asarray(initial.scheduled_raw_entries)-1)))
    # This additional multi-period drift test accommodates the saved initial
    # fertility normalization's 8.23e-7 error. It does not relax any dated gate.
    if drift > 1e-5 or queue_drift > 1e-5 or raw_queue_drift > 1e-5:
        raise ValueError(f"No-shock stationary drift exceeds diagnostic limit: {drift}, {queue_drift}")
    root, detail = finite(c, Path(folder)/"no_shock_root", 6, float(P.psi_child), deadline,
                          max_evaluations=6)
    if root.next_state is None:
        raise ValueError("Exact original-law no-shock root loop did not pass")
    import e5f_original_queue_terminal as terminal_solver
    terminal_controls=dict(c.controls,max_evaluations=16)
    endpoint=terminal_solver.solve_terminal(old=c.old,psi=float(P.psi_child),audit=c.audit,
        controls=terminal_controls,deadline=deadline,folder=Path(folder)/"original_terminal_replay")
    if not endpoint.verified:
        raise ValueError("Original-preference terminal steady-state replay did not pass")
    coordinates=np.array([q,float(P.pension),float(P.property_tax_lump_sum_transfer)])
    coordinate_gap=float(np.max(np.abs(np.asarray(endpoint.coordinates)/coordinates-1)))
    if coordinate_gap>1e-3:
        raise ValueError("Re-solved original terminal differs materially from calibrated steady state")
    c.driver.save(Path(folder)/"summary.json", dict(status="passed", exact_six_period_replay=True,
        no_shock_distribution_relative_l1=drift, no_shock_queue_relative_drift=queue_drift,
        no_shock_raw_queue_relative_drift=raw_queue_drift,
        original_terminal_verified=True,original_terminal_maximum_coordinate_relative_gap=coordinate_gap,
        no_shock_drift_limit=1e-5, finite_root_reproduced=True,
        demographic_rule="original_four_vintage_household_queue",
        spec_sha256=c.driver.sha(c.spec_path)))


def fixed_terminal_path(c, endpoint, folder, count, psi, deadline):
    from run_e5f_successive_surprises_overnight import standard_graphs
    controls = dict(c.controls)
    for key in ("automatic_fiscal_polish", "fiscal_tolerance", "fiscal_slope", "initial_jacobian"):
        controls.pop(key, None)
    controls["slope"] = controls.pop("market_slope", 1.63)
    controls["max_evaluations"] = 8
    initial = initial_coordinates(c, count).reshape(3, count+1)[:, :-1]
    target = np.asarray(endpoint.coordinates, dtype=float)
    weight = np.linspace(0., 1., count)
    guess = np.exp((1-weight)[None,:]*np.log(initial) + weight[None,:]*np.log(target[:,None]))
    # Do not pass the original queue state to a PersonPFState coercion. Queue
    # terminal distances are checked explicitly below with the same law.
    boundary = NS(parameters=endpoint.parameters, policy=endpoint.policy, asset_price=endpoint.asset_price)
    observations = []; snapshot = {}
    def observe(i, e, P, grid, shared):
        import run_e5f_transition_calibration as fertility
        observations.append(dict(period=i, calendar_year=2007+4*i, **fertility.period_fertility_diagnostics(e,P)))
        if i == 0:
            snapshot.clear(); snapshot.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=c.old.supply_rule)
    def progress(row):
        c.driver.save(Path(folder)/"latest_completed.json",row)
        if row.get("new_best"):c.driver.save(Path(folder)/"best_so_far.json",row)
    result = c.rebated.solve_rebated_forecast(inherited=c.rebated.InheritedState(2007,c.old.initial_state),
        psi=psi, old_state=c.old, terminal=boundary, demographic_primitives=None, count=count,
        initial_prices=guess[0], initial_pensions=guess[1], initial_transfers=guess[2],
        audit_controls=c.audit, root_controls=controls, deadline_monotonic=deadline,
        callback=progress, observer=observe)
    root = result.root_receipt
    c.driver.save(Path(folder)/"root_receipt.json",root)
    if result.path is not None:
        c.driver.save(Path(folder)/"rows.json",result.path.rows)
        c.driver.save(Path(folder)/"fertility.json",observations[-count:])
        actual = result.path.person_tail.terminal_state
        target_state = endpoint.state
        scale = max(float(target_state.g_pre.sum()),1e-15)
        d = dict(distribution_relative_l1=float(np.abs(actual.g_pre-target_state.g_pre).sum())/scale,
            population_relative_gap=abs(float(actual.g_pre.sum())/scale-1),
            queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries)/np.asarray(target_state.scheduled_entries)-1))),
            raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries)/np.asarray(target_state.scheduled_raw_entries)-1))),
            terminal_steady_state_verified=True, finite_path_converged=bool(root.get("finite_horizon_market_fiscal_converged")),
            horizon_comparison_passed=False, production_eligible=False)
        c.driver.save(Path(folder)/"terminal_distance.json",d)
        if snapshot:standard_graphs(snapshot,result,Path(folder)/"graphs")
        plot_case(c,folder,label=f"Permanent shock with solved terminal steady state: {count} periods",
            status="See root and terminal-distance receipts; horizon robustness remains unverified.")
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--spec",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--mode",choices=["smoke","finite","terminal","history"],required=True)
    ap.add_argument("--count",type=int,choices=[6,24,100],default=100)
    ap.add_argument("--seconds",type=float,required=True)
    args=ap.parse_args()
    c=load_context(args.spec);c.spec_path=args.spec
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    if (out/"experiment_contract.json").exists():raise ValueError("Refusing to overwrite experiment")
    remaining=min(args.seconds,float(c.spec["absolute_deadline_unix"])-time.time())
    if remaining < 120:raise TimeoutError("Afternoon budget expired")
    deadline=time.monotonic()+remaining
    c.driver.save(out/"experiment_contract.json",dict(mode=args.mode,count=args.count,
        demographic_rule="original_four_vintage_household_queue",conversion=1/2.1,
        migration=0.,historical_age_conditioning=False,person_headship_transition=False,
        original_parameters_fixed=True,psi_initial=float(c.old.parameters.psi_child),
        psi_permanent=c.spec["permanent_psi"],spec_sha256=c.driver.sha(args.spec),
        presentation_results_replaced=False,production_eligible=False))
    stop=threading.Event()
    def heartbeat():
        while not stop.wait(60):
            c.driver.save(out/"controller_heartbeat.json",dict(remaining_seconds=max(0,deadline-time.monotonic())))
            if time.monotonic()>=deadline:
                c.driver.save(out/"controller_failure.json",dict(error="Afternoon hard deadline"));os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    try:
        if args.mode!="smoke":
            smoke=read(c.spec["smoke_summary"])
            if smoke.get("status")!="passed" or smoke.get("spec_sha256")!=c.driver.sha(args.spec):
                raise ValueError("Exact-loop smoke prerequisite missing")
        with c.queue.original_queue_adapter(), original_receipts(c), capture_root_paths(c), c.cache.policy_cache(c.joined.pf,max_bytes=(12 if args.count==100 else 6)*1024**3):
            if args.mode=="smoke":
                original_no_shock_smoke(c,out,deadline)
            elif args.mode=="finite":
                finite(c,out,args.count,c.spec["permanent_psi"],deadline)
            elif args.mode=="terminal":
                import e5f_original_queue_terminal as terminal
                endpoint=terminal.solve_terminal(old=c.old,psi=c.spec["permanent_psi"],audit=c.audit,
                    controls=c.controls,deadline=min(deadline,time.monotonic()+3600),folder=out/"endpoint")
                if not endpoint.verified:raise RuntimeError("No verified terminal steady state within stage budget")
                with capture_paths(c,out/"transition"):
                    fixed_terminal_path(c,endpoint,out/"transition",args.count,c.spec["permanent_psi"],deadline)
            else:
                # Root mechanics, four empirical targets, scalar preference
                # search and paired policy loop remain the retained driver.
                # Adapter replacements apply to every date, including 2023+.
                if "resume_history" in read(c.spec["history_manifest"]):
                    raise ValueError("Original-law history must be refitted from its own 2007 state")
                c.driver.main(["--manifest",c.spec["history_manifest"],"--case","A0",
                    "--count",str(args.count),"--output",str(out),"--seconds",str(remaining)])
                if [row["year"] for row in read(out/"realized_fit.json")] != [2007,2011,2015,2019]:
                    raise ValueError("Original-law history did not complete all four windows")
                for p in out.glob("policies/*/rows.json"):
                    root=read(p.parent/"root_receipt.json")
                    plot_case(c,p.parent,label=f"Refitted history policy: {p.parent.name}",
                        status="Finite equilibrium; terminal steady state unverified." if root.get("converged") else "Unconverged diagnostic.")
        c.driver.verify_pins(c.spec["file_sha256"])
        c.driver.save(out/"controller_complete.json",dict(mode=args.mode,completed=True,production_eligible=False))
    except BaseException as exc:
        c.driver.save(out/"controller_failure.json",dict(error_type=type(exc).__name__,error=str(exc)))
        raise
    finally:stop.set()


if __name__=="__main__":main()
