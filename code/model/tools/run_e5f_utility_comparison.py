#!/usr/bin/env python3
"""Prepare/check a four-arm contract; evaluate only an explicitly approved one.

Use on Torch. Preparation and preflight never solve an equilibrium. This file
does not submit jobs, approve a design, retry failures, or promote a calibration.
"""
from __future__ import annotations

import argparse
import copy
import csv
from dataclasses import asdict
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import time

import numpy as np

import e5f_utility_comparison_runtime as adapter
import e5f_utility_comparison_design as design
import e5f_utility_recovery_policy_v1 as recovery

PARENT_LOCK = "6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1"
INADMISSIBLE_PREFIXES = (
    "Old-steady-state fertility normalization is not bracketed:",
    "Old-steady-state fertility normalization missed tolerance:",
    "Initial housing equilibrium failed its unchanged strict gate",
)
RECOVERY_POLICY = "reviewed_failure_v1"


def recovery_enabled(contract):
    choice = contract.get("candidate_failure_policy")
    if choice not in (None, RECOVERY_POLICY):
        raise RuntimeError("unreviewed candidate failure policy")
    return choice == RECOVERY_POLICY


def recovery_context(contract, contract_path, arm, request, output):
    """Bind structured evidence to this exact externally supervised attempt."""
    if request.get("candidate_id") != output.name:
        raise RuntimeError("candidate identity differs from its output directory")
    runner_hash = contract["files"][Path(__file__).name]["sha256"]
    if request.get("runner_source_sha256") != runner_hash:
        raise RuntimeError("case plan runner source differs from the contract")
    if (request.get("deadline_owner") != "controller"
            or request.get("controller_pid") != os.getppid()):
        raise RuntimeError("reviewed failure policy requires its owning controller")
    return recovery.Context(request["candidate_id"], request["controller_stage"],
        adapter.file_hash(contract_path), runner_hash,
        contract["arms"][arm]["objective"]["sha256"],
        design.canonical_fingerprint(request["point"]))


def authenticated_native_failure_type(contract, runtime):
    """Authenticate the defining solver and its unchanged gates before capture.

    The census is read from the actual exception, never reconstructed from its
    message. This check adds no numerical acceptance rule or model modification.
    """
    model = runtime["model"]
    relative = "code/model/intergen_eqscale_seq_optimized/solver.py"
    expected_path = (Path(contract["reference_root"]) / "source" / relative).resolve()
    manifest_pin = contract["parent_source_inventory"]
    if adapter.file_hash(manifest_pin["path"]) != manifest_pin["sha256"]:
        raise RuntimeError("native failure source manifest changed")
    expected_hash = read(manifest_pin["path"])["files"][relative]
    native_type = model.InfeasibleThetaError
    if (Path(model.__file__).resolve() != expected_path
            or adapter.file_hash(expected_path) != expected_hash
            or sys.modules.get(native_type.__module__) is not model
            or Path(native_type.__init__.__code__.co_filename).resolve() != expected_path
            or native_type.__name__ != "InfeasibleThetaError"
            or model.DEAD_MASS_TOL != 1e-12 or model.DEAD_VALUE_CUTOFF != -1e9):
        raise RuntimeError("native failure class or numerical gate differs from the pinned solver")
    return native_type, dict(native_source_path=str(expected_path),
        native_source_sha256=expected_hash, native_manifest_sha256=manifest_pin["sha256"],
        native_gate_tolerance=model.DEAD_MASS_TOL, native_value_cutoff=model.DEAD_VALUE_CUTOFF)


def failure_status(exc):
    """Only named native gate rejections are inadmissible search proposals.

    This never changes a gate or retries a point. Unclassified failures,
    timeouts and accounting/source/array errors stop further dispatch.
    """
    if isinstance(exc, RuntimeError) and str(exc).startswith(INADMISSIBLE_PREFIXES):
        return "inadmissible_parameter_proposal"
    return "unclassified_failure_stop"


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=adapter.file_hash(path))


def pair_runtime(reference):
    reference = Path(reference).resolve()
    lock = reference / "inputs/launch_lock.json"
    if adapter.file_hash(lock) != PARENT_LOCK:
        raise RuntimeError("only the reviewed immutable September25 nightpair is supported")
    # These expected hashes are from the already reviewed parent's launch lock.
    payload = read(lock)
    for name in ("run_pair.py", "render_pair.py"):
        if adapter.file_hash(reference / name) != payload["runtime_file_sha256"][name]:
            raise RuntimeError("reviewed parent runtime changed: " + name)
    pair = adapter.import_file("four_arm_frozen_pair", reference / "run_pair.py")
    if pair.WORK.resolve() != reference:
        raise RuntimeError("parent source path differs from the frozen execution path")
    os.environ["EXPECTED_PAIR_LOCK_SHA256"] = PARENT_LOCK
    old = pair.load_ancestor()
    lock, objective, bank = pair.read_contract(old)
    return pair, old, lock, objective, bank


def prepare(args):
    if args.output.exists():
        raise RuntimeError("preparation output already exists; no overwrite or automatic retry")
    pair, old, lock, base_objective, _ = pair_runtime(args.reference_root)
    provenance, pension = read(args.provenance), read(args.pension_receipt)
    if pension["results_main"]["rho_2007"] != adapter.PENSION_RATIO:
        raise RuntimeError("pension measurement receipt differs from adopted target")
    by_id = {row["id"]: row for row in provenance["target_rows"]}
    for row in base_objective["target_rows"]:
        source = by_id[row["restriction_id"]]
        if row["target"] != source["target"] or row["actual_weight"] != source["weight"]:
            raise RuntimeError("full target/provenance values or weights differ")
    selected = {}
    for name in ("greaney_179", "oasi_087510"):
        best = read(args.reference_root / "results/run_001" / name / "best_so_far.json")
        receipt_path = Path(best["case_path"]) / "receipt.json"
        receipt = read(receipt_path)
        if receipt["target_system_sha256"] != lock["objective_sha256"] or receipt["point"] != best["point"]:
            raise RuntimeError("parent selected-point receipt differs")
        selected[name] = dict(point=best["point"], receipt=pin(receipt_path))
    ref_receipt = read(selected["oasi_087510"]["receipt"]["path"])
    rent = float(ref_receipt["price"]) * float(ref_receipt["user_cost_rate"])
    adapter.reference_composite_factor(adapter.ALPHA0, rent)
    args.output.mkdir(parents=True)
    target_payload = {k: v for k, v in base_objective.items() if k != "parameter_restrictions"}
    arms = {}
    for arm, (housing, exponent) in adapter.ARMS.items():
        objective = copy.deepcopy(base_objective)
        if housing == "shares":
            objective["parameter_restrictions"] = [r for r in objective["parameter_restrictions"] if r["parameter"] != "h_P"]
            objective["parameter_restrictions"].extend(dict(parameter=k, lower=0., upper=.25, transform="softzero")
                                                       for k in adapter.SHARE_COORDINATES)
        objective_path = args.output / arm / "objective.json"
        write(objective_path, objective)
        seeds = []
        for name, record in selected.items():
            point = copy.deepcopy(record["point"])
            if housing == "shares":
                point.pop("h_P")
                # Starting values, not externally estimated preferences.
                point.update(delta_alpha_jump=.08, delta_alpha=.035)
            seeds.append(dict(label=name+"_structural_seed", point=point))
        arms[arm] = dict(housing=housing, benefit_exponent=exponent, objective=pin(objective_path),
                         free_count=len(objective["parameter_restrictions"]), seeds=seeds)
    tool_dir = Path(__file__).resolve().parent
    files = {name: pin(tool_dir/name) for name in (
        "run_e5f_utility_comparison.py", "e5f_utility_comparison_runtime.py",
        "e5f_utility_comparison_design.py", "run_e5f_utility_comparison_search.py",
        "collect_e5f_utility_comparison.py", "e5f_utility_recovery_policy_v1.py")}
    files["submit_e5f_utility_comparison.sh"] = pin(tool_dir.parent/"submit_e5f_utility_comparison.sh")
    files.update(provenance=pin(args.provenance), pension_receipt=pin(args.pension_receipt),
                 pension_measurement_contract=pin(args.pension_contract), timing_receipt=pin(args.timing_receipt))
    timing = read(args.timing_receipt)
    budget = design.size_search_budget(objective_solve_p90_seconds=timing["p90_seconds"],
        objective_overhead_seconds=args.overhead_seconds, workers_per_arm=args.workers_per_arm,
        total_limit_seconds=args.total_seconds, repeat_reserve_seconds=args.repeat_reserve_seconds,
        export_reserve_seconds=args.export_reserve_seconds, initial_population=40,
        de_generations=3, objective_timeout_seconds=3100)
    catalog = design.define_arms(base_objective["parameter_restrictions"])
    target_contract = design.validate_target_contract(provenance["target_rows"])
    common_seeds = {name:{k:row["point"][k] for k in design.COMMON_NAMES} for name,row in selected.items()}
    floor_seeds = {name:{"h_P":row["point"]["h_P"]} for name,row in selected.items()}
    share_seeds = {name:{"delta_alpha_jump":.08,"delta_alpha":.035} for name in selected}
    bank = design.build_proposal_bank(catalog,shared_seeds=common_seeds,floor_seeds=floor_seeds,
        share_seeds=share_seeds,broad_count=14,medium_count=12,local_count=12,
        medium_unit_scale=.15,local_unit_scale=.04,rng_seed=20260925,budget=budget)
    for name,value in (("budget",budget),("arm_catalog",catalog),("proposal_bank",bank),("target_provenance",target_contract)):
        write(args.output/(name+".json"),value)
        files[name]=pin(args.output/(name+".json"))
    contract = dict(schema="utility_comparison_preparation_v1", status="prepared_awaiting_author_review",
        reference_root=str(args.reference_root.resolve()), parent_lock=pin(pair.LOCK),
        parent_source_inventory=pin(pair.SOURCE_MANIFEST), parent_objective=pin(pair.OBJECTIVE),
        parent_provenance=pin(args.provenance), files=files, arms=arms,
        common_target_fingerprint=adapter.canonical_hash(target_payload),
        pension_ratio=adapter.PENSION_RATIO, selected_seed_receipts=selected,
        reference_rent=rent, reference_rent_source=selected["oasi_087510"]["receipt"],
        normalization=dict(material="A(alpha)=K(alpha0,r*)/K(alpha,r*)",
            K="alpha**alpha*((1-alpha)/r*)**(1-alpha)",
            interpretation="At reference prices, compensated material expenditure is e(m) times the childless level; the floor arm additionally requires rent*h_P.",
            status="experimental cross-child-state utility restriction; mathematical and author review pending",
            fertility_target=2.1, fertility_tolerance=.0005, maximum_stationary_solves_per_objective=23),
        inherited_nonpreference_contract=dict(earnings="B15", entry_wealth="inherited level marginal plus rank coupling",
            adult_entry="half at16 / half at20; births divided by2.1 exactly once",
            dependent_departure_probability=2/9, source_and_numerical_gates="immutable nightpair source and reviewed runtime",
            transition_scope="no policy transitions in this comparison"),
        required_readout=dict(target_rows=13, weighted_rows=12, all_free_and_fixed_parameters=True,
            standard_graph_count=17, repetitions=2, compare_original_to_each_repeat=True, numerical_tolerance=0.),
        unresolved_assumptions=["experimental reference-rent normalization needs author adoption for this comparison",
            "0.86 is a fixed curvature sensitivity, not an externally estimated structural parameter",
            "native first-birth housing proxy does not reproduce the PSID estimator",
            "positive-all-estates model observer differs from child-directed SCF target",
            "entry wealth, estate recipient allocation and mortality mapping remain inherited pending author decisions"],
        proposed_budget=budget, differential_evolution=dict(generations=3,mutation_factor=.7,crossover_probability=1.),
        failure_policy=dict(inadmissible_native_prefixes=list(INADMISSIBLE_PREFIXES),
            inadmissible_proposals_consume_slots=True, maximum_inadmissible_fraction_per_barrier=.5,
            unclassified_or_timeout="stop new dispatch; do not retry", gates="unchanged"),
        launch_permitted=False)
    write(args.output / "contract.json", contract)
    print(json.dumps(dict(status=contract["status"], contract=pin(args.output/"contract.json"),
                          reference_rent=rent, arms=list(arms), native_solve_count=0)))


def verified_contract(path):
    contract = read(path)
    expected = os.environ.get("EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256")
    if expected != adapter.file_hash(path):
        raise RuntimeError("explicit reviewed preparation-contract fingerprint required")
    executing={"run_e5f_utility_comparison.py":__file__,
               "e5f_utility_comparison_runtime.py":adapter.__file__,
               "e5f_utility_comparison_design.py":design.__file__,
               "e5f_utility_recovery_policy_v1.py":recovery.__file__}
    recovery_enabled(contract)
    for name,actual_path in executing.items():
        if Path(actual_path).resolve() != Path(contract["files"][name]["path"]).resolve():
            raise RuntimeError("executing source is not the contract-pinned source: "+name)
    pins = list(contract["files"].values()) + [contract[k] for k in ("parent_lock", "parent_source_inventory", "parent_objective", "reference_rent_source")]
    pins += [record["objective"] for record in contract["arms"].values()]
    for item in pins:
        if adapter.file_hash(item["path"]) != item["sha256"]:
            raise RuntimeError("pinned input changed: " + item["path"])
    for record in contract["arms"].values():
        objective = read(record["objective"]["path"])
        payload = {k:v for k,v in objective.items() if k != "parameter_restrictions"}
        if adapter.canonical_hash(payload) != contract["common_target_fingerprint"]:
            raise RuntimeError("mixed target/weight/measurement contract across arms")
    return contract


def setup(contract, arm, output, *, output_reserved=False):
    pair, old, lock, _, _ = pair_runtime(contract["reference_root"])
    tax, plan, selected, objective, fiscal = adapter.configure_runtime(pair, old, lock, arm,
        contract["arms"][arm]["objective"]["path"], contract["reference_rent"])
    if not output_reserved:
        output.mkdir(parents=True, exist_ok=False)
    os.environ["NUMBA_CACHE_DIR"] = os.environ.get("UTILITY_COMPARISON_NUMBA_CACHE",str(output/"numba_cache"))
    runtime = old.setup_runtime(tax, plan, selected, output)
    return pair, old, lock, tax, selected, objective, runtime, fiscal


def preflight(args):
    contract = verified_contract(args.contract)
    pair,old,lock,tax,selected,objective,runtime,fiscal = setup(contract,args.arm,args.output)
    records = []
    for seed in contract["arms"][args.arm]["seeds"]:
        P = old.apply_point(selected, seed["point"])
        shared = runtime["model"].precompute_shared(P, selected["b_grid"])
        original = runtime["model"].precompute_shared._four_arm_native_precompute(P, selected["b_grid"])
        unchanged = []
        expected_changes = set()
        if contract["arms"][args.arm]["housing"] == "shares":
            expected_changes.add("escale_flat")
        if contract["arms"][args.arm]["benefit_exponent"] != 1.0:
            expected_changes.update(("psi_v","psi_flat","type_map","type_psi","type_cb","type_hb","n_types"))
        for key, value in vars(original).items():
            if key not in expected_changes:
                if not np.array_equal(np.asarray(value), np.asarray(getattr(shared,key))):
                    raise RuntimeError("undeclared native shared-array change: "+key)
                unchanged.append(key)
        for n in range(int(P.n_parity)):
            for m in range(int(P.n_child_states)):
                expected = P.psi_child*float(m)**contract["arms"][args.arm]["benefit_exponent"] if m<=n else 0.
                if shared.psi_v[n,m] != expected:
                    raise RuntimeError("native benefit array differs from current-child reward")
        from e5f_stationary_paygo import bind_initial_balanced_pension
        balanced, accounts = bind_initial_balanced_pension(P, payroll_tax=old.TAX)
        adapter.verify_pension_ratio({"actual_accounts": accounts["predicted_accounts"]})
        if P.adult_entry_clock != "split_birth_vintage" or not np.isfinite(shared.escale_flat).all():
            raise RuntimeError("native preference/entry binding failed")
        records.append(dict(seed=seed["label"], point=seed["point"], payroll_tax=old.TAX,
            pension_period=balanced.pension, benefit=shared.psi_v.tolist(),
            material_multiplier=shared.escale_flat.tolist(), alpha=shared.alpha_flat.tolist(),
            unchanged_native_fields=unchanged))
    receipt = dict(status="native_zero_solve_preflight_passed", arm=args.arm, records=records,
                   fiscal_rule=fiscal, native_solve_count=0, contract_sha256=adapter.file_hash(args.contract))
    write(args.output/"preflight.json",receipt)
    print(json.dumps(dict(status=receipt["status"], arm=args.arm, payroll_tax=old.TAX,
                         native_solve_count=0, receipt=str(args.output/"preflight.json"))))


def evaluate(args):
    contract = verified_contract(args.contract)
    budget=read(contract["files"]["budget"]["path"])
    design.validate_budget(budget)
    if (contract["status"] != "author_approved_frozen_design" or not contract["launch_permitted"]
            or contract.get("approved_budget") != budget or not contract.get("approved_assumptions")):
        raise RuntimeError("no objective launch: frozen design/budget and common assumptions require author review")
    request = read(args.case_plan)
    if request["contract_sha256"] != adapter.file_hash(args.contract) or request["arm"] != args.arm:
        raise RuntimeError("case plan does not belong to the approved frozen contract")
    cap=float(request["objective_cap_seconds"])
    if not np.isfinite(cap) or not 0<cap<=budget["objective_timeout_seconds"]:
        raise ValueError("case wall cap exceeds the finite approved budget")
    if not np.isfinite(float(request["deadline_epoch"])):
        raise ValueError("case requires a finite absolute deadline")
    supervised = recovery_enabled(contract)
    context = recovery_context(contract, args.contract, args.arm, request, args.output) if supervised else None
    # Reserve ownership before the exception path can write a failure receipt.
    # A duplicate request fails here without changing any existing case file.
    args.output.mkdir(parents=True,exist_ok=False)
    if supervised:
        # Attempt metadata stays outside scientific receipts: original/repeat
        # equality must not depend on their distinct IDs or search stages.
        write(args.output / "attempt_provenance.json", dict(context=asdict(context),
            arm=args.arm, candidate_failure_policy=RECOVERY_POLICY,
            case_plan_sha256=adapter.file_hash(args.case_plan)))
    # Absolute deadlines include setup; a shell/Slurm wall limit is also required.
    deadline = min(float(request["deadline_epoch"]), time.time()+cap)
    if deadline <= time.time(): raise TimeoutError("case budget exhausted")
    def alarm(*_): raise TimeoutError("objective wall budget exhausted")
    # Opt-in repair moves wall-time enforcement outside the numerical process.
    # The native per-call absolute deadline and solve-count guards stay intact.
    if not supervised:
        signal.signal(signal.SIGALRM,alarm)
        signal.setitimer(signal.ITIMER_REAL,deadline-time.time())
    start = time.monotonic()
    native_type, native_source = None, None
    try:
        pair,old,lock,tax,selected,objective,runtime,fiscal = setup(contract,args.arm,args.output,output_reserved=True)
        if supervised:
            native_type, native_source = authenticated_native_failure_type(contract, runtime)
        receipt=old.evaluate_point(tax=tax,objective=objective,selected=selected,runtime=runtime,
            point=request["point"],output=args.output/"case",deadline_epoch=deadline,
            graphs=bool(request.get("graphs",False)))
        receipt.update(complete_objective_wall_seconds=time.monotonic()-start,
                       comparison_contract_sha256=adapter.file_hash(args.contract))
        old.write(args.output/"case/receipt.json",receipt)
    except Exception as exc:
        failure = dict(status=failure_status(exc), error_type=type(exc).__name__,
                       error=str(exc), automatic_retry=False)
        if supervised:
            failure.update(arm=args.arm, native_source=native_source,
                           candidate_failure_policy=RECOVERY_POLICY)
            try:
                evidence = recovery.capture_native_failure(exc, expected_native_type=native_type,
                    native_gate_tolerance=(native_source or {}).get("native_gate_tolerance"), context=context)
                payload = asdict(evidence)
                # Invalid diagnostics remain fatal and must not mask the raw
                # exception by failing while serializing NaN or native objects.
                json.dumps(payload, allow_nan=False)
                failure["recovery_evidence"] = payload
            except Exception as capture_error:
                # Preserve the original exception even if its diagnostic is malformed.
                failure["recovery_capture_error"] = type(capture_error).__name__ + ": " + str(capture_error)
        write(args.output/"failure.json", failure)
        raise
    finally:
        if not supervised:
            signal.setitimer(signal.ITIMER_REAL,0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=("prepare","preflight","evaluate"),required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--reference-root",type=Path)
    parser.add_argument("--provenance",type=Path)
    parser.add_argument("--pension-receipt",type=Path)
    parser.add_argument("--pension-contract",type=Path)
    parser.add_argument("--timing-receipt",type=Path)
    parser.add_argument("--workers-per-arm",type=int)
    parser.add_argument("--total-seconds",type=float)
    parser.add_argument("--repeat-reserve-seconds",type=float)
    parser.add_argument("--export-reserve-seconds",type=float)
    parser.add_argument("--overhead-seconds",type=float)
    parser.add_argument("--contract",type=Path)
    parser.add_argument("--arm",choices=tuple(adapter.ARMS))
    parser.add_argument("--case-plan",type=Path)
    args=parser.parse_args()
    {"prepare":prepare,"preflight":preflight,"evaluate":evaluate}[args.stage](args)


if __name__ == "__main__": main()
