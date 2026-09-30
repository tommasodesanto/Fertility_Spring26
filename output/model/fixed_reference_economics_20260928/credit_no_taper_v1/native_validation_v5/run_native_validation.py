#!/usr/bin/env python3
"""Two-case Torch validation for the isolated renter no-taper overlay.

The controller is deliberately a validation harness, not an equilibrium or
recalibration driver.  It starts a fresh interpreter for each case so that the
overlay is registered before any runtime package import.  Case order is fixed:
the exact flag-off control, then the flag-on renter-taper removal.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib
import importlib.abc
import importlib.util
import json
import math
import os
import pickle
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

LABEL = "2007 stationary reference — block0506, September 28 verified export"
ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
PACKET = ROOT / "output/model/fixed_reference_economics_20260928/credit_no_taper_v1"
MANIFEST = ROOT / "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
MANIFEST_SHA = "147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4"
OVERLAY_SHA = "9c6def300f76b2d5ac55c392e8a595fca881b1d78b15c1cba96016a30a3b83b9"
NATIVE = {
    "parameters.py": "66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464",
    "solver.py": "b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1",
    "kernels.py": "639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27",
}
CURRENT_RUNTIME_SHA = "a9625c354909d5f7027d0e2b6c989475ab905c046b4742e6dbe54585bc9367d3"
CASES = ("control", "renter_no_taper")
FLAG = "renter_no_taper_estate_bound"
TOL = 2e-10


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def table(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scientific_plan_fields(plan):
    """Fields that must be identical between Torch preflight and production."""
    keys = ("schema", "reference_label", "reference_manifest_sha256", "cases",
            "maximum_lifecycle_solves", "case_seconds", "total_seconds", "threads",
            "memory_gib", "grid_nodes", "lambda_d", "natural_credit",
            "renormalize_births", "fixed_prices_rents_psi_fiscal", "compiled_mode",
            "required_readout", "overlay_parameters_path", "overlay_parameters_sha256",
            "base_driver_path", "pinned_files")
    return {key: plan[key] for key in keys}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def progress(output, phase, **extra):
    write(output / "progress.json", dict(reference_label=LABEL, phase=phase, time_epoch=time.time(), **extra))


def verify_plan(path):
    require(sys.platform == "linux" and os.environ.get("SLURM_JOB_ID", "").isdigit(),
            "Torch Slurm execution is required")
    plan = read(path)
    require(plan["schema"] == "block0506_renter_no_taper_native_validation_v3", "wrong plan schema")
    require(plan["reference_label"] == LABEL, "wrong reference label")
    require(plan["reference_manifest_sha256"] == MANIFEST_SHA == sha(MANIFEST), "manifest pin differs")
    require(plan["cases"] == list(CASES) and plan["maximum_lifecycle_solves"] == 2, "case/solve contract differs")
    require(plan["case_seconds"] == 360 and plan["total_seconds"] == 1200, "budget differs")
    require(plan["threads"] == 1 and plan["memory_gib"] == 24, "resource contract differs")
    require(plan["grid_nodes"] == 160 and plan["lambda_d"] == 0 and not plan["natural_credit"], "credit/grid contract differs")
    require(plan["renormalize_births"] is False and plan["fixed_prices_rents_psi_fiscal"] is True,
            "economic fixed-input contract differs")
    require(plan["compiled_mode"] is True and plan["required_readout"] == {"fit_rows": 14, "parameter_rows": 31, "standard_plots": 17},
            "compiled/readout contract differs")
    require(plan["preflight_contract"] == {"required": True, "receipt": "controller_mock_tests.json",
            "requires_matching_plan_and_driver": True, "requires_matching_scientific_fields": True},
            "preflight contract differs")
    require(plan["driver_sha256"] == sha(__file__), "driver pin differs")
    require(plan["overlay_parameters_sha256"] == OVERLAY_SHA, "overlay pin differs")
    for name, digest in NATIVE.items():
        require(sha(ROOT / "code/model/intergen_eqscale_seq_optimized" / name) == digest,
                "native source pin differs: " + name)
    for item in plan["pinned_files"]:
        file = Path(item["path"])
        require(file.is_file() and sha(file) == item["sha256"], "pinned file differs: " + str(file))
    return plan


class OverlayParametersFinder(importlib.abc.MetaPathFinder):
    """Serve one authenticated parameter module without changing native source."""
    def __init__(self, overlay):
        self.overlay = Path(overlay).resolve()

    def find_spec(self, fullname, path=None, target=None):
        if fullname == "intergen_eqscale_seq_optimized.parameters":
            return importlib.util.spec_from_file_location(fullname, self.overlay)
        return None


def install_identity_guard_extension(plan):
    """Install the sole overlay exception before the parent authenticator.

    The native guard remains authoritative for the no-overlay case.  With the
    registered overlay, this extension accepts only the exact plan-pinned
    ``parameters`` module; the native solver and every other package module
    must still be inside the current native package.
    """
    overlay = Path(plan["overlay_parameters_path"])
    require(sha(overlay) == OVERLAY_SHA, "effective overlay bytes differ")
    package_name = "intergen_eqscale_seq_optimized"
    parameters_name = package_name + ".parameters"
    package_path = ROOT / "code/model" / package_name
    require(package_name not in sys.modules and parameters_name not in sys.modules,
            "runtime package imported before overlay registration")
    require(sha(package_path / "solver.py") == NATIVE["solver.py"], "native solver hash differs")
    tools = ROOT / "code/model/tools"
    require(sha(tools / "e5f_current_transition_runtime.py") == CURRENT_RUNTIME_SHA,
            "native identity runtime hash differs")
    sys.path.insert(0, str(tools))
    native = importlib.import_module("e5f_current_transition_runtime")
    original_guard = native.require_current_model
    expected_solver = (package_path / "solver.py").resolve()
    expected_package = package_path.resolve()

    def extended_guard(model):
        parameters = sys.modules.get(parameters_name)
        if parameters is None:
            return original_guard(model)
        require(Path(getattr(model, "__file__", "")).resolve() == expected_solver
                and sha(expected_solver) == NATIVE["solver.py"], "Current-source model identity failed")
        require(Path(getattr(parameters, "__file__", "")).resolve() == overlay.resolve()
                and sha(overlay) == OVERLAY_SHA, "overlay parameter identity failed")
        for name, module in list(sys.modules.items()):
            if not name.startswith(package_name) or not getattr(module, "__file__", None):
                continue
            path = Path(module.__file__).resolve()
            if name == parameters_name:
                require(path == overlay.resolve(), "mixed overlay parameter module")
            elif not path.is_relative_to(expected_package):
                raise RuntimeError("Mixed model package: " + name)
        return None

    native.require_current_model = extended_guard
    finder = OverlayParametersFinder(overlay)
    sys.meta_path.insert(0, finder)
    return dict(native=native, finder=finder, overlay=overlay, original_guard=original_guard)


def verify_authenticated_overlay(plan, extension, model):
    """Verify the actual authenticated runtime uses the registered overlay."""
    parameters = sys.modules.get("intergen_eqscale_seq_optimized.parameters")
    require(parameters is not None, "authenticated runtime did not import overlay parameters")
    extension["native"].require_current_model(model)
    require(Path(parameters.__file__).resolve() == extension["overlay"].resolve()
            and sha(parameters.__file__) == OVERLAY_SHA, "authenticated overlay bytes differ")
    require(model.setup_parameters is parameters.setup_parameters,
            "actual native model did not resolve overlay setup_parameters")
    require(model.unsecured_debt_floor is parameters.unsecured_debt_floor,
            "actual native model did not resolve overlay unsecured_debt_floor")


def renter_incidence(ev, P, grid, model, reference_P):
    """Actual occupied renter paths by age, using realised current masses.

    ``policy_mass_branches`` keeps DUE owner-stayer policies separate.  The
    renter slice of every branch is therefore a realised tenure choice, not a
    full-grid location-probability average; in particular no global assertion
    about unused ``loc_probs`` entries is made.
    """
    from e5f_overnight_estate_audit import policy_mass_branches
    branches = policy_mass_branches(ev, P)
    for owner_mass, _, _ in branches[1:]:
        require(float(np.asarray(owner_mass)[:, 0].sum()) == 0.0,
                "owner-stayer branch leaked into realised renter incidence")
    rows = []
    old_weight = np.asarray(reference_P.debt_taper_weights, dtype=float)
    old_cap = np.asarray(reference_P.debt_caps, dtype=float)
    survival = np.asarray(P.survival_probs, dtype=float) if bool(P.use_age_survival) else None
    ages = float(P.age_start) + float(P.da) * np.arange(int(P.J))
    for j, age in enumerate(ages):
        mass = sum(np.asarray(branch[0][:, 0, :, j, :, :, :], dtype=float) for branch in branches)
        # All renter branches use the same policy in the unchanged DUE design;
        # calculate numerator branch-by-branch to avoid a policy average.
        policy_negative = 0.0
        old_bind = old_violation = new_bind = new_violation = 0.0
        current_negative = float(mass[grid < 0.0].sum())
        renter_mass = float(mass.sum())
        old_floor = np.asarray(model.renter_borrowing_floor(reference_P, grid, j), dtype=float)
        death = j == int(P.J) - 1 or (survival is not None and survival[j] < 1.0)
        new_floor = np.zeros_like(grid) if death else np.minimum(grid, 0.0)
        for branch_mass, branch_saving, _ in branches:
            m = np.asarray(branch_mass[:, 0, :, j, :, :, :], dtype=float)
            bnext = np.asarray(branch_saving[:, 0, :, j, :, :, :], dtype=float)
            node = np.asarray(grid).reshape((-1,) + (1,) * 4)
            neg_current = np.broadcast_to(node < 0.0, m.shape)
            neg_save = bnext < 0.0
            policy_negative += float(m[neg_save].sum())
            old_node = old_floor.reshape((-1,) + (1,) * 4)
            new_node = new_floor.reshape((-1,) + (1,) * 4)
            oldb = np.abs(bnext - old_node) <= 1e-8
            newb = np.abs(bnext - new_node) <= 1e-8
            old_bind += float(m[neg_current & oldb].sum())
            new_bind += float(m[neg_current & newb].sum())
            old_violation += float(m[bnext < old_node - 1e-8].sum())
            new_violation += float(m[bnext < new_node - 1e-8].sum())
        rows.append(dict(age_left=float(age), occupied_renter_mass=renter_mass,
            current_negative_renter_assets_mass=current_negative,
            policy_saving_negative_mass=policy_negative,
            old_floor_binding_mass=old_bind, new_floor_binding_mass=new_bind,
            old_floor_violation_mass=old_violation, new_floor_violation_mass=new_violation,
            old_floor=float(old_floor.min()), new_floor_min=float(new_floor.min()),
            death_impossible=not death, old_taper_weight=float(old_weight[j + 1]), old_cap=float(old_cap[j + 1])))
    for row in rows:
        require(row["new_floor_violation_mass"] <= TOL, "operative renter lower-bound violation")
    return rows


def child(args, plan):
    import numpy as np
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    launch = read(out.parent / "launch.json")
    require(launch["plan"] == plan and launch["plan_sha256"] == sha(args.plan), "controller contract absent")
    launch_started, launch_deadline = launch_clock(plan)
    require(launch["started_epoch"] == launch_started and launch["deadline_epoch"] == launch_deadline
            and args.deadline <= launch_deadline, "child launch clock was reset or widened")
    index = CASES.index(args.case)
    if index:
        completed = read(out.parent / "latest_completed.json")["completed"]
        require([x["case"] for x in completed] == list(CASES[:index]), "control not completed first")
        control = out.parent / "control" / "receipt.json"
        require(sha(control) == completed[0]["receipt_sha256"] and read(control)["control"]["status"] == "passed",
                "control receipt changed or failed")
    # This must precede loading the parent authenticator, which imports runtime packages.
    extension = install_identity_guard_extension(plan)
    base = load_module(Path(plan["base_driver_path"]), "_renter_no_taper_reference_helpers")
    progress(out, "authenticate")
    manifest, contract, objective, runtime, prepared, reference = base.authenticate(out)
    verify_authenticated_overlay(plan, extension, prepared.rt["model"])
    P = copy.deepcopy(reference["parameters"])
    grid = np.asarray(reference["b_grid"]).copy()
    require(len(grid) == 160 and P.Nb == 160, "original 160-node grid is required")
    P.native_inherited_distribution_evidence_dir = str(out / "inherited_state_diagnostics")
    original = base.serialized(vars(P))
    if args.case == "renter_no_taper":
        P.renter_no_taper_estate_bound = True
        P.lambda_d = 0.0
        sys.modules["intergen_eqscale_seq_optimized.parameters"].build_debt_caps(P)
    else:
        require(not bool(getattr(P, FLAG, False)), "reference flag unexpectedly active")
    experiment = base.serialized(vars(P))
    changed = {key: dict(reference=original.get(key), experiment=experiment.get(key))
               for key in set(original) | set(experiment) if original.get(key) != experiment.get(key)}
    allowed = {FLAG, "debt_taper_weights", "debt_caps", "mean_labor_income_by_age"}
    require((not index and not changed) or (index and FLAG in changed and "debt_taper_weights" in changed
            and set(changed).issubset(allowed)), "unplanned parameter change")
    for key in ("debt_caps", "mean_labor_income_by_age"):
        if key in changed:
            require(changed[key]["reference"] == changed[key]["experiment"],
                    "rebuilt bookkeeping field changed economically: " + key)
    require(np.array_equal(P.debt_caps, np.zeros_like(P.debt_caps)), "new unsecured credit cap is nonzero")
    require(np.array_equal(P.owner_ltv_multipliers, reference["parameters"].owner_ltv_multipliers),
            "owner LTV rule changed")
    require(float(P.lambda_d) == 0.0 and float(P.psi_child) == float(reference["parameters"].psi_child),
            "credit or preferences changed")
    actual = base.actual_parameters(prepared, P, grid)
    require(len(actual) == 31 and actual == base.actual_parameters(prepared, reference["parameters"], grid),
            "a calibrated/external parameter changed")
    rt, model, cal = prepared.rt, prepared.rt["model"], prepared.rt["primitive"].pf.calendar
    price = np.asarray(reference["solution"].p_eq).copy()
    sd = model.precompute_shared(P, grid)
    progress(out, "one_lifecycle_solve", deadline_epoch=args.deadline)
    require(time.time() < args.deadline, "deadline before solve")
    started = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd, verbose=False, fast_stats=False)
    elapsed = time.monotonic() - started
    require(time.time() < args.deadline, "case deadline exceeded in solve")
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction["stationary_post_fertility_nesting_l1"], 5e-9, "cohort reconstruction")
    runtime.require_abs_gate(reconstruction["stationary_feasibility_projection_mass"], 0.0, "cohort projection")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol, evaluation=ev, stationary_g_pre=pre,
                  supply_rule=supply, demographic_seed=reference.get("demographic_seed"))
    progress(out, "household_and_fiscal_gates", lifecycle_solve_seconds=elapsed)
    # Persist the actual fixed-rule fiscal accounts before the maintained
    # pension certificate can fail.  The subsequent gate is deliberately not
    # weakened or cleared by changing P.pension.
    from e5f_social_security import fiscal_accounts
    write(out / "fiscal_pre_gate.json", cal.jsonable(fiscal_accounts(ev.g_current, P)))
    gates = base.gates(packet, prepared, out, stationary=True)  # preserves any fiscal failure; no clearing adjustment
    rows = renter_incidence(ev, P, grid, model, reference["parameters"])
    table(out / "renter_incidence_by_age.csv", rows)
    fertility = {kind: rt["observe_initial_fertility"](ev, P, age_projection=kind)
                 for kind in ("uniform_birth_time", "constant_post_cell")}
    housing = rt["observe_initial_housing_wealth"](ev, P, grid, sd, diagnostic_enabled=True,
        age_projection="uniform_within_age_cell", diagnostic_allow_family_proxies=True, include_wealth=True,
        include_birth_response=True)
    recent = rt["observe_recent_parent_flow"](ev, P, diagnostic_enabled=True, snapshot=rt["SNAPSHOT"],
        age_projection=rt["AGE_PROJECTION"], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=args.case, reference_checkpoint_sha256=manifest["checkpoint"]["sha256"]))
    completed_fertility = float(rt["chain"].extract_moments(sol, P)["tfr"])
    fits = runtime.score_targets(objective, fertility, housing, recent["model_value"], completed_fertility)
    require(len(fits) == 14, "full fit table missing")
    control = base.exact_control(reference, packet, fits, manifest, prepared, out) if index == 0 else None
    for row in fits:
        if row["moment"] == "initial_normalization":
            row["role"] = "reference benchmark; no fertility renormalization"
    params = copy.deepcopy(manifest["full_parameter_table"])
    for row in params:
        row["status"] = "Fixed inherited reference value; renter lower-bound overlay separately disclosed"
    table(out / "target_fit.csv", fits)
    table(out / "parameters.csv", params)
    write(out / "observers.json", cal.jsonable(dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
    progress(out, "standard_17_plot_rendering")
    rt["audit"].standard_diagnostics(packet, out, validate_production_young=False)
    plots = sorted(x.name for x in (out / "standard_diagnostics").glob("*.png"))
    require(plots == sorted(manifest["standard_diagnostic_names"]), "standard 17 diagnostic set differs")
    if index == 0:
        for name in plots:
            require(sha(out / "standard_diagnostics" / name) == manifest["artifact_hashes"]["standard_diagnostics/" + name],
                    "control diagnostic differs: " + name)
    require(time.time() < args.deadline, "case deadline exceeded in reporting")
    require(sha(args.plan) == launch["plan_sha256"] and sha(__file__) == plan["driver_sha256"], "pin changed during case")
    artifacts = {str(path.relative_to(out)): sha(path) for path in sorted(out.rglob("*")) if path.is_file()}
    receipt = dict(status="passed", reference_label=LABEL, case=args.case, lifecycle_solves=1,
        lifecycle_solve_seconds=elapsed, plan_sha256=sha(args.plan), driver_sha256=sha(__file__),
        overlay_parameters_sha256=sha(Path(plan["overlay_parameters_path"])),
        identity_guard_exception="Only intergen_eqscale_seq_optimized.parameters at the plan-pinned overlay SHA is external; the native solver and every other model module remain current-source pinned.",
        reference_checkpoint=manifest["checkpoint"],
        source_manifest=manifest["source_manifest"], target_weight_fingerprint=manifest["target_weight_fingerprint"],
        parameter_changes=changed, control=control, fixed_prices=price.tolist(), fixed_psi=float(P.psi_child),
        fiscal_residuals=gates["fiscal"], fiscal_certificate=gates["fiscal_certificate"], renter_incidence_rows=len(rows),
        standard_plot_count=len(plots), output_sha256=artifacts, normalization_performed=False,
        interpretation="Fixed-price/rent, fixed-psi and fixed-fiscal conditional validation. It is not a cleared equilibrium, demographic stationary endpoint, or transition.")
    write(out / "receipt.json", cal.jsonable(receipt))
    progress(out, "complete")


def launch_clock(plan):
    started = float(os.environ["NATIVE_VALIDATION_STARTED_EPOCH"])
    deadline = float(os.environ["NATIVE_VALIDATION_DEADLINE_EPOCH"])
    require(math.isfinite(started) and deadline == started + plan["total_seconds"] and time.time() <= deadline,
            "launcher start/deadline missing, reset, or expired")
    return started, deadline


def controller(args, plan, *, child_factory=subprocess.Popen, clock=time.time, sleeper=time.sleep, pin_sha=sha):
    started, deadline = launch_clock(plan)
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / "launch.json", dict(reference_label=LABEL, plan=plan, plan_sha256=sha(args.plan),
        started_epoch=started, deadline_epoch=deadline, slurm_job=os.environ["SLURM_JOB_ID"]))
    completed = []
    write(args.output / "latest_completed.json", dict(reference_label=LABEL, completed=[], lifecycle_solves=0,
        remaining_cases=2))
    for case in CASES:
        require(clock() < deadline and pin_sha(__file__) == plan["driver_sha256"], "global deadline or pin failure")
        case_deadline = min(deadline, clock() + plan["case_seconds"])
        command = [sys.executable, str(Path(__file__).resolve()), "--plan", str(args.plan), "--output",
                   str(args.output / case), "--case", case, "--deadline", str(case_deadline)]
        with (args.output / (case + ".log")).open("w") as log:
            process = child_factory(command, stdout=log, stderr=subprocess.STDOUT,
                env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                         NUMBA_NUM_THREADS="1", MPLBACKEND="Agg"), start_new_session=True)
            try:
                while process.poll() is None:
                    progress(args.output, "case_running", case=case, completed_cases=len(completed),
                             elapsed_seconds=time.time() - started, case_deadline_epoch=case_deadline)
                    if clock() >= case_deadline:
                        raise TimeoutError("case/global deadline reached: " + case)
                    sleeper(2)
                require(process.returncode == 0, "case failed; stop without retry: " + case)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=3)
        receipt = read(args.output / case / "receipt.json")
        require(receipt["status"] == "passed" and receipt["lifecycle_solves"] == 1, "missing case receipt")
        if case == "control":
            require(receipt["control"]["status"] == "passed", "taper removal blocked by failed control")
        completed.append(dict(case=case, receipt_sha256=sha(args.output / case / "receipt.json")))
        write(args.output / "latest_completed.json", dict(reference_label=LABEL, completed=completed,
              lifecycle_solves=len(completed), remaining_cases=2-len(completed)))
    require(clock() < deadline, "global deadline in finalization")
    verify_plan(args.plan)  # authenticate all pins again after both cases
    completion_hashes = {str(p.relative_to(args.output)): sha(p) for p in sorted(args.output.rglob("*")) if p.is_file()}
    write(args.output / "completed.json", dict(status="passed", reference_label=LABEL, completed=completed,
        lifecycle_solves=2, elapsed_seconds=clock()-started, plan_sha256=sha(args.plan),
        completion_hashes=completion_hashes,
        interpretation="Conditional fixed-price validation only; fiscal residuals are reported, not cleared."))


def mock_tests(plan, plan_path):
    """Exercise the production controller with subprocess-only zero-solve children."""
    import shutil
    events = []
    original_started, original_deadline = os.environ.get("NATIVE_VALIDATION_STARTED_EPOCH"), os.environ.get("NATIVE_VALIDATION_DEADLINE_EPOCH")
    root = Path(tempfile.mkdtemp(prefix="native_validation_controller_"))
    try:
        started = time.time(); os.environ["NATIVE_VALIDATION_STARTED_EPOCH"] = str(started)
        os.environ["NATIVE_VALIDATION_DEADLINE_EPOCH"] = str(started + plan["total_seconds"])
        def factory(mode):
            def launch(command, **kwargs):
                case = command[command.index("--case") + 1]; output = command[command.index("--output") + 1]
                script = ("import json,pathlib,sys,time; p=pathlib.Path(sys.argv[1]); p.mkdir(); "
                          + ("time.sleep(0.1); " if mode == "sleep" else "")
                          + ("sys.exit(7)" if mode == "fail" else
                             "json.dump({'status':'passed','lifecycle_solves':1,'control':{'status':'passed'}},open(p/'receipt.json','w'))"))
                return subprocess.Popen([sys.executable, "-c", script, output], **kwargs)
            return launch
        def run(name, mode="pass", **kw):
            ns = argparse.Namespace(plan=plan_path, output=root / name)
            return controller(ns, plan, child_factory=factory(mode), sleeper=lambda _: None, **kw)
        run("success")
        require(read(root / "success" / "completed.json")["lifecycle_solves"] == 2, "actual controller success failed")
        for name, mode, expected in (("control_failure", "fail", RuntimeError), ("timeout", "sleep", TimeoutError)):
            try:
                if mode == "sleep":
                    # Initial controller checks happen before the deadline; on
                    # the first live-child poll the case deadline is reached.
                    # The real child sleeps briefly, so this is no spin poll.
                    ticks = iter((started + 1., started + 1., float(os.environ["NATIVE_VALIDATION_DEADLINE_EPOCH"])))
                    run(name, mode,
                        clock=lambda: next(ticks, float(os.environ["NATIVE_VALIDATION_DEADLINE_EPOCH"])),
                        sleeper=lambda _: time.sleep(.01))
                else: run(name, mode)
            except expected: events.append(name)
            else: raise RuntimeError(name + " did not stop")
        try: run("success")
        except FileExistsError: events.append("duplicate_output_refusal")
        else: raise RuntimeError("duplicate output did not stop")
        try: run("changed_pin", pin_sha=lambda _: "changed")
        except RuntimeError: events.append("changed_pin")
        else: raise RuntimeError("changed pin did not stop")
        # This is the real registration/identity-extension path, not an import toy.
        extension = install_identity_guard_extension(plan)
        solver = importlib.import_module("intergen_eqscale_seq_optimized.solver")
        verify_authenticated_overlay(plan, extension, solver)
        wrong = "intergen_eqscale_seq_optimized.unapproved_external"
        from types import SimpleNamespace
        sys.modules[wrong] = SimpleNamespace(__file__="/tmp/unapproved_external.py")
        try:
            try: extension["native"].require_current_model(solver)
            except RuntimeError: events.append("external_model_module_rejected")
            else: raise RuntimeError("external package module was accepted")
        finally:
            sys.modules.pop(wrong, None)
        # Native-shaped [b,h,I,j,n,m,z] arrays exercise the actual incidence
        # function, including its owner-stayer split and renter slice.
        tools = ROOT / "code/model/tools"
        if str(tools) not in sys.path: sys.path.insert(0, str(tools))
        shape = (3, 2, 1, 2, 1, 1, 1)
        g = np.zeros(shape); g[0, 0, 0, :, 0, 0, 0] = 0.2; g[1, 0, 0, :, 0, 0, 0] = 0.3
        g[1, 1, 0, :, 0, 0, 0] = 0.4  # current owner mass includes stayers
        stay = np.zeros(shape); stay[1, 1, 0, :, 0, 0, 0] = 0.4
        saving = np.zeros(shape); saving[0, 0, 0, 0, 0, 0, 0] = -1.0
        ev = SimpleNamespace(g_current=g, g_stay_distribution=stay,
                             policy=SimpleNamespace(bp_pol=saving, bp_pol_stay=np.zeros(shape), c_pol=np.ones(shape), c_pol_stay=np.ones(shape)))
        fixture_P = SimpleNamespace(native_due_stayer_credit=True, use_age_survival=False, J=2, age_start=20., da=1.,
            debt_taper_weights=np.ones(3), debt_caps=np.zeros(3))
        fixture_model = SimpleNamespace(renter_borrowing_floor=lambda P, b, j: np.minimum(np.asarray(b), 0.0))
        rows = renter_incidence(ev, fixture_P, np.array([-1., 0., 1.]), fixture_model, fixture_P)
        require(len(rows) == 2 and rows[0]["current_negative_renter_assets_mass"] > 0.0,
                "native-shaped renter incidence fixture failed")
        bad_saving = saving.copy(); bad_saving[0, 0, 0, 1, 0, 0, 0] = -1.0
        ev.policy.bp_pol = bad_saving
        try: renter_incidence(ev, fixture_P, np.array([-1., 0., 1.]), fixture_model, fixture_P)
        except RuntimeError: events.append("negative_terminal_renter_saving_rejected")
        else: raise RuntimeError("negative terminal renter saving was accepted")
    finally:
        if original_started is None: os.environ.pop("NATIVE_VALIDATION_STARTED_EPOCH", None)
        else: os.environ["NATIVE_VALIDATION_STARTED_EPOCH"] = original_started
        if original_deadline is None: os.environ.pop("NATIVE_VALIDATION_DEADLINE_EPOCH", None)
        else: os.environ["NATIVE_VALIDATION_DEADLINE_EPOCH"] = original_deadline
        shutil.rmtree(root)
    return dict(status="passed", model_solves=0, actual_controller=True, negative_tests=events,
                plan_sha256=sha(plan_path), driver_sha256=sha(__file__),
                scientific_plan_fields=scientific_plan_fields(plan),
                scientific_plan_sha256=hashlib.sha256(json.dumps(scientific_plan_fields(plan), sort_keys=True,
                    separators=(",", ":")).encode()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=CASES)
    parser.add_argument("--deadline", type=float)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args(); args.plan = args.plan.resolve(); args.output = args.output.resolve()
    plan = verify_plan(args.plan)
    if args.self_test_only:
        require(not args.case and args.deadline is None and not args.output.exists(), "invalid mock-test invocation")
        args.output.mkdir(parents=True)
        write(args.output / "controller_mock_tests.json", mock_tests(plan, args.plan)); return
    require(not args.output.exists(), "Versioned output must not already exist")
    try:
        if args.case:
            require(args.deadline is not None and time.time() < args.deadline <= time.time() + plan["case_seconds"],
                    "child needs a bounded absolute deadline")
            child(args, plan)
        else:
            require(args.deadline is None, "controller owns deadline")
            controller(args, plan)
    except BaseException as error:
        if args.output.exists():
            write(args.output / "failure.json", dict(status="failed", reference_label=LABEL, error_type=type(error).__name__,
                error=str(error), traceback=traceback.format_exc(), retries=0, time_epoch=time.time()))
        raise


if __name__ == "__main__":
    main()
