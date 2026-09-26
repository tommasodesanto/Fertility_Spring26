#!/usr/bin/env python3
"""Torch-only, bounded three-case estate-receiver stationary diagnostic.

This is an experimental estate diagnostic based on the frozen September 25
low-tax PAYGO runtime.  It is neither a recalibration nor a demographic or
policy transition.  In particular, it holds the selected fixed preference
normalization and reports (rather than closes) the endogenous birth/entry gap.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pickle
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

CASES = ("control", "net_valuation", "net_valuation_transfer")
TAX = 0.08751017424959717
LOCK_SHA = "6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1"
SOURCE_SHA = "237904131d159f775c7ae89d1bbf1e8d1dd70c79ad012c80a36d948658d6f9c6"
OBJECTIVE_SHA = "cce32a5d8208603e9237607da8253995a5f5517ef39fd5a1ee75bee90718ed58"
ROOT = Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1")
SELECTED = ROOT / "results/run_001/oasi_087510/worker_05/point_04"
CHECKPOINT = SELECTED / "initial_state.pkl.gz"
FINAL_RECEIPT = ROOT / "results/run_001/oasi_087510/selected_export/final_pair_receipt.json"
RUNNER = ROOT / "run_pair.py"
MAX_SOLVES = 42
MAX_RECEIVER_ITERATIONS = 40
BUDGET_SECONDS = 3 * 3600


def require_torch() -> None:
    """Refuse model imports and hashed source reads outside the Torch allocation."""
    if not os.environ.get("SLURM_JOB_ID") and os.environ.get("E5F_ESTATE_RECEIVER_TORCH") != "1":
        raise RuntimeError("Torch-only: set by a Torch allocation; local execution is prohibited")
    if "torch" not in socket.gethostname().lower() and not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Torch-only host guard failed")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def table(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def load_runner():
    spec = importlib.util.spec_from_file_location("e5f_sep25_run_pair", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("frozen run_pair loader unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_adapter():
    try:
        import e5f_estate_receiver_adapter as adapter
    except ImportError as exc:
        raise RuntimeError("blocked preflight: lead adapter e5f_estate_receiver_adapter.py is absent") from exc
    required = ("CASES", "install", "configure", "estate_accounts", "solve_case")
    if tuple(adapter.CASES) != CASES or any(not callable(getattr(adapter, key, None)) for key in required[1:]):
        raise RuntimeError("blocked preflight: estate adapter interface differs from reviewed contract")
    return adapter


def verified_runtime(output: Path):
    """Pin frozen source then recreate its reviewed purchase runtime unchanged."""
    runner = load_runner()
    os.environ["EXPECTED_PAIR_LOCK_SHA256"] = LOCK_SHA
    old = runner.load_ancestor()
    lock, objective, _ = runner.read_contract(old)
    if lock["source_manifest_sha256"] != SOURCE_SHA or lock["objective_sha256"] != OBJECTIVE_SHA:
        raise RuntimeError("frozen source or objective pin differs")
    runner.configure(old, "oasi_087510", lock)
    tax, plan, selected, objective = runner.prepare(old, lock)
    if not CHECKPOINT.is_file() or not FINAL_RECEIPT.is_file():
        raise RuntimeError("actual winning checkpoint or final selected receipt is missing")
    final = json.loads(FINAL_RECEIPT.read_text())
    if final.get("selected_case") != str(SELECTED) or final.get("source_manifest_sha256") != SOURCE_SHA:
        raise RuntimeError("final receipt does not identify worker_05/point_04 under the pinned source")
    with gzip.open(CHECKPOINT, "rb") as stream:
        selected_case = pickle.load(stream)
    selected_receipt = json.loads((SELECTED / "receipt.json").read_text())
    if (selected_receipt.get("source_manifest_sha256") != SOURCE_SHA or
            selected_receipt.get("target_system_sha256") != OBJECTIVE_SHA or
            selected_receipt.get("case_checkpoint_sha256") != sha(CHECKPOINT)):
        raise RuntimeError("winning checkpoint fails receipt/source/target verification")
    if (selected_case["parameters"].psi_child != selected_receipt["normalization"]["psi_child"]
            or selected_case["parameters"].tau_pay != TAX
            or selected_case["parameters"].adult_entry_clock != "split_birth_vintage"):
        raise RuntimeError("winning serialized preference, tax or entry clock differs from its receipt")
    runtime = old.setup_runtime(tax, plan, selected, output / "runtime")
    import numpy as np
    np.testing.assert_array_equal(selected_case["b_grid"], selected["b_grid"])
    runtime["original_initial_prices"] = np.asarray(selected["solution"].p_eq).copy()
    return old, tax, plan, selected_case, objective, runtime, selected_receipt


def heartbeat(root: Path, started: float, case: str, status: str, solves: int) -> None:
    write(root / "heartbeat.json", dict(status=status, case=case, elapsed_seconds=time.monotonic()-started,
        native_ge_solves=solves, cap_native_ge_solves=MAX_SOLVES, deadline_seconds=BUDGET_SECONDS,
        updated_epoch=time.time()))


def run_cases(*, adapter, model, native_solver, parameters, b_grid, initial_prices, output: Path,
              deadline_epoch: float, report_case) -> list[dict]:
    """Exact ordered case loop; injectable seams make the Torch smoke deterministic."""
    completed: list[dict] = []
    solves = 0
    def counted_solver(**kwargs):
        nonlocal solves
        if solves >= MAX_SOLVES or time.time() >= deadline_epoch:
            raise TimeoutError("shared native solve or time budget exhausted")
        solves += 1
        write(output / "native_solve_count.json", dict(started=solves, maximum=MAX_SOLVES))
        return native_solver(**kwargs)
    for case in CASES:
        if time.time() >= deadline_epoch:
            raise TimeoutError("three-hour stage deadline exhausted before " + case)
        case_out = output / case
        case_out.mkdir(parents=True, exist_ok=False)
        P = adapter.configure(parameters, case)
        if getattr(P, "estate_receiver", "none") != "none":
            raise RuntimeError("native estate_receiver must be off for isolated diagnostic hooks")
        before = solves
        sol, P, prices, fiscal, estate = adapter.solve_case(
            model=model, native_solver=counted_solver, parameters=P, b_grid=b_grid,
            initial_prices=initial_prices, case=case, deadline_epoch=deadline_epoch,
            output_dir=case_out, max_solves=min(MAX_RECEIVER_ITERATIONS, MAX_SOLVES-solves))
        used = int(estate["native_solves"])
        if used != solves-before:
            raise RuntimeError("adapter and driver native-solve counters disagree")
        accounts = adapter.estate_accounts(sol, P, prices)
        if case == "net_valuation_transfer" and abs(float(accounts["residual"])) > 1e-10:
            raise RuntimeError("estate paid/generated accounting residual exceeds 1e-10")
        report = report_case(case_out, case, sol, P, prices, fiscal, estate, accounts)
        item = dict(case=case, status="completed", native_ge_solves=used, estate_accounts=accounts, report=report)
        completed.append(item)
        write(output / "latest_completed.json", item)
        best = min(completed, key=lambda x: float(x["report"].get("loss", math.inf)))
        write(output / "best_so_far.json", dict(description="descriptive inherited objective only; not selection or adoption", **best))
    return completed


def live_report(old, tax, objective, runtime, selected_receipt, root_started):
    """Reuse the frozen observers/gates; net estate accounts remain supplemental."""
    def report(case_out, case, sol, P, prices, fiscal, estate, accounts):
        import numpy as np
        primitive, model, chain = runtime["primitive"], runtime["model"], runtime["chain"]
        grid = runtime.get("selected", {}).get("b_grid") if "selected" in runtime else None
        if grid is None: raise RuntimeError("runtime lacks frozen wealth grid")
        shared = model.precompute_shared(P, grid); P._fert2_probs = sol.fert2_probs.copy()
        policy = primitive.pf.calendar.policy_from_solution(sol, prices, P, grid, shared)
        pre, reconstruction = primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol, policy, P, grid, shared)
        operator = primitive.pf.transition.operator_gates(sol, policy, pre, P, grid, shared); operator.update(reconstruction)
        for key in ("stationary_post_fertility_nesting_l1", "one_step_constant_path_nesting_l1",
                    "mature_flow_abs_error", "birth_flow_abs_error", "topcode_adjusted_birth_flow_abs_error"):
            if not math.isfinite(float(operator[key])) or abs(operator[key]) > 5e-9:
                raise RuntimeError("native operator gate failed: " + key)
        if (abs(operator["zero_entry_mass_accounting_residual"]) > 2e-8
                or operator["stationary_feasibility_projection_mass"] > 1e-6):
            raise RuntimeError("native mass/projection gate failed")
        supply = primitive.pf.calendar.HousingSupplyRule("static-elastic", float(prices[0]),
            float(P.H0[0]*(P.user_cost_rate*prices[0]/P.r_bar[0])**P.xi_supply[0]), float(P.xi_supply[0]))
        evaluation = primitive.pf.calendar.evaluate_period(prices, pre, P, grid, shared, primitive.pf.calendar.SolveCounter(), supply_rule=supply, supplied_policy=policy)
        if evaluation.relative_market_residual > 2e-4: raise RuntimeError("native market gate failed")
        budget = primitive.dated_budget(evaluation, P, shared, grid, float(P.user_cost_rate*prices[0]))
        purchase = runtime["accounting"].audit_purchase_accounting(evaluation, P, shared, grid, model)
        fiscal = runtime["certify_initial_pension"](evaluation.g_current, P,
            marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
        packet = dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,
                      supply_rule=supply,solution=sol,stationary_g_pre=pre,
                      contract_sha256=OBJECTIVE_SHA,
                      demographic_seed=runtime["selected"].get("demographic_seed"))
        with gzip.open(case_out/"initial_state.pkl.gz", "wb", compresslevel=1) as stream:
            pickle.dump(packet, stream, protocol=5)
        arrays = runtime["audit"].policy_array_audit(dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,solution=sol,stationary_g_pre=pre), case_out)
        if arrays["occupied_negative_steps"] or any(x["nonfinite"] or x["minimum"] < 0 or x["maximum"] > 1 for x in arrays["probabilities"].values()): raise RuntimeError("native value/probability gate failed")
        early = dict(fertility={p: runtime["observe_initial_fertility"](evaluation,P,age_projection=p) for p in ("uniform_birth_time","constant_post_cell")}, housing_wealth=runtime["observe_initial_housing_wealth"](evaluation,P,grid,shared,diagnostic_enabled=True,age_projection="uniform_within_age_cell",diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True))
        recent = runtime["observe_recent_parent_flow"](evaluation,P,diagnostic_enabled=True,snapshot=runtime["SNAPSHOT"],age_projection=runtime["AGE_PROJECTION"],diagnostic_allow_residence_proxy=True,input_provenance={"case_id":case})
        rows = tax.target_rows(objective, early, recent["model_value"], float(chain.extract_moments(sol,P)["tfr"]))
        if len(rows) != 13 or sum(r["weight"] != "" for r in rows) != 12: raise RuntimeError("full frozen 13-row target table failed")
        for row in rows:
            if any(not math.isfinite(float(row[k])) for k in ("target", "model", "gap")):
                raise RuntimeError("Nonfinite target-fit row")
        loss = sum(float(r["loss_contribution"]) for r in rows if r["loss_contribution"] != "")
        if case == "control":
            reference_packet = runtime["selected"]
            for key in ("V", "g", "bp_pol", "c_pol", "hR_pol", "fert_probs", "fert2_probs"):
                np.testing.assert_array_equal(getattr(sol,key),getattr(reference_packet["solution"],key),
                                              err_msg="control replay: "+key)
            np.testing.assert_array_equal(prices, reference_packet["solution"].p_eq)
            np.testing.assert_array_equal(pre, reference_packet["stationary_g_pre"])
            reference = {r["moment"]: r for r in csv.DictReader((SELECTED / "target_fit.csv").open())}
            if any(float(row["model"]) != float(reference[row["moment"]]["model"]) for row in rows):
                raise RuntimeError("control does not reproduce every frozen target model value")
            if loss != float(selected_receipt["loss"]):
                raise RuntimeError("control does not reproduce frozen inherited objective")
        table(case_out / "target_fit.csv", rows)
        baseline = list(csv.DictReader((SELECTED / "parameters.csv").open()))
        actual = tax.actual_parameters(P)
        for row in baseline:
            name = row["parameter"]
            if name in actual and float(row["estimate"]) != float(actual[name]):
                raise RuntimeError("structural parameter changed: " + name)
            if name == "pension_period": row["estimate"] = float(P.pension)
            if name == "psi_child":
                row["estimate"] = float(P.psi_child)
                row["status"] = "held at selected preference; not renormalized in this test"
            elif row["status"] == "experimental free coordinate":
                row["status"] = "held at selected estimate; no recalibration"
        table(case_out / "parameters.csv", baseline + [dict(parameter="estate_transfer_period",estimate=accounts["transfer"],lower="",upper="",near_bound="",status="endogenous balanced transfer in receiver case; zero otherwise")])
        runtime["audit"].standard_diagnostics(dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,solution=sol,stationary_g_pre=pre), case_out, validate_production_young=False)
        if len(list((case_out/"standard_diagnostics").glob("*.png"))) != 17: raise RuntimeError("standard 17-plot packet incomplete")
        receipt = dict(status="completed_experimental_estate_receiver_diagnostic", case=case, loss=loss, payroll_tax=TAX,
            source_manifest_sha256=SOURCE_SHA,target_system_sha256=OBJECTIVE_SHA, fixed_psi_child=float(P.psi_child),
            selected_reference=selected_receipt["case_checkpoint_sha256"],
            fixed_psi_birth_replacement_residual=float(sol.adult_entry_stationary_relative_gap),
            birth_based_adult_entry=float(sol.adult_entry_adjusted_birth_children)/2.1,
            normalized_entrant_flow=float(sol.entry_rate),
            population_closure="fixed normalized entry/age composition; diagnostic, not imposed closed reproduction",
            estate_accounts=accounts, fiscal=fiscal, market_residual=evaluation.relative_market_residual, household_budget=budget, purchase_accounting=purchase, policy_array_gates=arrays,
            experimental_changes="net estate valuation; equal age-45--65 transfer only in transfer case; inherited experimental PAYGO .08751017424959717; not adopted pension rule")
        write(case_out / "receipt.json", receipt); return receipt
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("preflight", "smoke", "run"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(); require_torch()
    if args.output_dir.exists(): raise RuntimeError("output directory must be new")
    args.output_dir.mkdir(parents=True)
    old, tax, plan, selected, objective, runtime, selected_receipt = verified_runtime(args.output_dir)
    adapter = load_adapter(); adapter_receipt = adapter.install(runtime["model"], args.output_dir / "adapter")
    diagnostic_files = [Path(__file__), Path(adapter.__file__),
        Path(__file__).with_name("test_e5f_estate_receiver_adapter.py"),
        Path(__file__).with_name("test_run_e5f_estate_receiver_probe.py")]
    diagnostic_pins = {p.name:sha(p) for p in diagnostic_files}
    preflight = dict(status="source_preflight_passed_no_solve", lock_sha256=LOCK_SHA, source_manifest_sha256=SOURCE_SHA,
        objective_sha256=OBJECTIVE_SHA, selected_checkpoint=str(CHECKPOINT), selected_checkpoint_sha256=sha(CHECKPOINT),
        actual_psi_child=float(selected["parameters"].psi_child), payroll_tax=TAX, adapter=adapter_receipt, native_ge_solves=0,
        selected_chosen_solve_seconds=selected_receipt.get("chosen_solve_seconds"),
        diagnostic_file_sha256=diagnostic_pins,
        indicative_native_seconds_at_reference_speed=MAX_SOLVES*float(selected_receipt["chosen_solve_seconds"]),
        estate_loop_seconds_reference="September18 fixed-preference sandbox:54 minutes; a different model, not a current runtime guarantee",
        experimental_changes="net estate valuation and optional equal age-45--65 transfers; inherited experimental fiscal setting")
    write(args.output_dir / "preflight.json", preflight)
    if args.stage == "preflight": return
    previous_path = os.environ.get("E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT")
    if not previous_path:
        raise RuntimeError("An explicit successful source-preflight receipt is required")
    previous = json.loads(Path(previous_path).read_text())
    for key in ("source_manifest_sha256", "objective_sha256", "selected_checkpoint_sha256", "diagnostic_file_sha256"):
        if previous.get(key) != preflight[key]:
            raise RuntimeError("Preflight fingerprint changed: " + key)
    if previous.get("status") != "source_preflight_passed_no_solve":
        raise RuntimeError("Earlier source preflight failed")
    if args.stage == "smoke":
        # Exercise the exact case-control loop and the actual estate fixed-point
        # adapter on deterministic fixtures before any production native solve.
        tests = [str(p) for p in diagnostic_files if p.name.startswith("test_")]
        with (args.output_dir/"tests.log").open("w") as stream:
            subprocess.run([sys.executable,"-m","pytest","-q",*tests],
                           stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=600)
        write(args.output_dir/"complete.json",dict(status="exact_loop_fixture_smoke_passed",
            native_ge_solves=0, diagnostic_file_sha256=diagnostic_pins,
            source_manifest_sha256=SOURCE_SHA,objective_sha256=OBJECTIVE_SHA,
            selected_checkpoint_sha256=preflight["selected_checkpoint_sha256"],
            scope="driver three-case loop and real adapter fixed point on deterministic fixtures; native baseline replay is first production case"))
        return
    if args.stage == "run":
        smoke_path = os.environ.get("E5F_ESTATE_PROBE_SMOKE_RECEIPT")
        if not smoke_path or not Path(smoke_path).is_file():
            raise RuntimeError("production refuses absent successful exact-loop smoke receipt")
        smoke = json.loads(Path(smoke_path).read_text())
        if (smoke.get("status") != "exact_loop_fixture_smoke_passed" or
                smoke.get("source_manifest_sha256") != SOURCE_SHA or
                smoke.get("objective_sha256") != OBJECTIVE_SHA or
                smoke.get("diagnostic_file_sha256") != diagnostic_pins or
                smoke.get("selected_checkpoint_sha256") != preflight["selected_checkpoint_sha256"]):
            raise RuntimeError("production refuses unsuccessful smoke or mismatched source fingerprint")
    started = time.monotonic(); deadline = time.time() + BUDGET_SECONDS
    runtime["selected"] = selected
    write(args.output_dir/"latest_completed.json",dict(status="none", cases=[]))
    write(args.output_dir/"best_so_far.json",dict(status="none", description="descriptive comparison only; no candidate selection"))
    heartbeat(args.output_dir, started, "starting", "running", 0)
    done = threading.Event()
    def pulse():
        while not done.wait(55):
            count_path = args.output_dir/"native_solve_count.json"
            count = json.loads(count_path.read_text())["started"] if count_path.exists() else 0
            heartbeat(args.output_dir, started, "active_case", "running", count)
    threading.Thread(target=pulse, daemon=True).start()
    try:
        completed = run_cases(adapter=adapter, model=runtime["model"], native_solver=runtime["solve_balanced_initial_equilibrium"], parameters=selected["parameters"], b_grid=selected["b_grid"], initial_prices=runtime["original_initial_prices"], output=args.output_dir, deadline_epoch=deadline, report_case=live_report(old,tax,objective,runtime,selected_receipt,started))
    except Exception as exc:
        write(args.output_dir / "failure.json", dict(status="stopped_preserving_partial_results", error_type=type(exc).__name__, error=str(exc), elapsed_seconds=time.monotonic()-started))
        raise
    finally:
        done.set()
    write(args.output_dir / "complete.json", dict(status="exact_three_case_loop_complete", stage=args.stage, cases=completed, native_ge_solve_cap=MAX_SOLVES, receiver_iteration_cap=MAX_RECEIVER_ITERATIONS, budget_seconds=BUDGET_SECONDS, report_pdf="outstanding: JSON/CSV only until safe frozen PDF assembly is reviewed"))


if __name__ == "__main__": main()
