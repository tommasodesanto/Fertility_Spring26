#!/usr/bin/env python3
"""Minimal authenticated fixed-credit smoke and exact flag-off control."""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import importlib.util
import json
import math
import os
import sys
import time
import types
from pathlib import Path

import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
HERE = Path(__file__).resolve().parent
BASE_PATH = ROOT / "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py"
OVERLAY = HERE / "overlay"
ORIGINAL = ROOT / "code/model/intergen_eqscale_seq_optimized"
MANIFEST = ROOT / "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
LABEL = "2007 stationary reference — block0506, September 28 verified export"
ACTIVE_PLAN_SHA = ""
PINS = {
    "base_authenticator": "96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44",
    "manifest": "147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4",
    "checkpoint": "b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d",
    "parameters.py": "d7b4d23c4bab2153dca0e86a6f189fe6aed6676c72090999d5a1c9105a8d5f8c",
    "solver.py": "cefc1627ff3adb124d223a114122a0a578b5c4dc889635b3deb97c6f2b805a3e",
    "kernels.py": "379d179a8f80477e8d813f60bc526e700a270e37e0f411988aabd9a06fd74352",
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def progress(output, phase, **extra):
    write(output / "progress.json", dict(reference_label=LABEL, phase=phase, time_epoch=time.time(), **extra))


def load_file(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "cannot load " + str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def authenticate(output):
    require(sha(MANIFEST) == PINS["manifest"], "reference manifest hash differs")
    require(sha(BASE_PATH) == PINS["base_authenticator"], "base authenticator hash differs")
    base = load_file(BASE_PATH, "fixed_price_reference_authenticator_v2")
    result = base.authenticate(output)
    manifest = result[0]
    require(manifest["checkpoint"]["sha256"] == PINS["checkpoint"], "checkpoint identity differs")
    return base, result


def overlay_namespace():
    """Load the overlay normally under a fresh package name with frozen fallbacks."""
    for name in ("parameters.py", "solver.py", "kernels.py"):
        require(sha(OVERLAY / name) == PINS[name], "overlay hash differs: " + name)
    package_name = "fixed_credit_runtime_model_v2"
    require(not any(k == package_name or k.startswith(package_name + ".") for k in sys.modules),
            "fresh overlay namespace already imported")
    package = types.ModuleType(package_name)
    package.__package__ = package_name
    package.__path__ = [str(OVERLAY), str(ORIGINAL)]
    package.__spec__ = importlib.util.spec_from_loader(package_name, loader=None, is_package=True)
    sys.modules[package_name] = package
    parameters = importlib.import_module(package_name + ".parameters")
    kernels = importlib.import_module(package_name + ".kernels")
    model = importlib.import_module(package_name + ".solver")
    expected = {"parameters": OVERLAY / "parameters.py", "kernels": OVERLAY / "kernels.py", "solver": OVERLAY / "solver.py"}
    modules = {"parameters": parameters, "kernels": kernels, "solver": model}
    origins = {}
    for name, module in modules.items():
        origin = Path(module.__file__).resolve()
        require(origin == expected[name].resolve(), "wrong module origin: " + name)
        require(sha(origin) == PINS[name + ".py"], "loaded module hash differs: " + name)
        origins[name] = str(origin)
    require(Path(model.__file__).resolve() == (OVERLAY / "solver.py").resolve(), "model.__file__ differs")
    return parameters, kernels, model, origins


def compiled_fixtures(kernels):
    grid = np.array([-2.0, 0.0, 2.0])
    args = (np.array([2.0, 2.0, 2.0]), np.array([2.0, 2.0, 2.0]),
            np.array([[100.0], [0.0], [-100.0]]), np.zeros((3, 1)), 0, grid,
            np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.0]),
            np.array([0.5]), np.array([1.0]), 1.0, 10.0, 1e-6, 0.0, 0.0,
            0.5, 0.5, 0.95, 0.0, 0.0, 0.381966, 0.618034, 1e-5)
    legacy = kernels.full_renter_block_kernel(*args)[1]
    positive = kernels.full_renter_block_kernel(*args, fixed_renter_floor=-1.0)[1]
    zero = kernels.full_renter_block_kernel(*args, fixed_renter_floor=0.0)[1]
    require(legacy[1, 0] >= 0 and positive[1, 0] < -0.5 and zero[1, 0] >= 0, "renter fixture failed")
    b = np.array([-1.0, 0.0, 1.0]); V = np.zeros((3, 2, 1, 1, 1))
    h = np.zeros((1, 2)); dp = np.zeros((1, 2, 1, 1)); bm = np.full((1, 2, 1, 1), -99.0)
    birth = np.zeros((1, 1, 2, 2), dtype=np.bool_); grant = np.zeros((1, 2, 1, 1))
    _, rejected = kernels.tenure_choice_kernel(V, b, np.array([[0.0, 0.5]]), h, dp, bm, birth, grant, V, False, False, True)
    _, boundary = kernels.tenure_choice_kernel(V, b, np.array([[0.0, 1.0]]), h, dp, bm, birth, grant, V, False, False, True)
    require(rejected[0, 1, 0, 0, 0] != 0 and boundary[0, 1, 0, 0, 0] == 0, "sale gate fixture failed")
    require(kernels.NUMBA_AVAILABLE and kernels.full_renter_block_kernel.signatures and kernels.tenure_choice_kernel.signatures,
            "compiled kernels were not executed")
    return {"renter_legacy_at_zero": float(legacy[1, 0]), "renter_positive_cap": float(positive[1, 0]),
            "renter_zero_cap": float(zero[1, 0]), "compiled": True}


def cash_proof(P, grid, model, price):
    sd = model.precompute_shared(P, grid)  # genuine reference parameter/shared-data construction
    b = -0.25581395348837255
    require(np.any(np.isclose(grid, b, rtol=0, atol=1e-14)), "documented entrant wealth is absent")
    rows = []
    for zz in (0, 1):
        z = float(np.asarray(P.z_grid)[zz])
        y = float(model.income_at_state(P, 0, 0, z))
        cash = float(P.R_gross) * b + y
        purchase_cash = b + y / float(P.R_gross)
        rows.append(dict(z_index=zz, z=z, income=y, renter_cash=cash, purchase_cash=purchase_cash))
        require(cash < 0 and purchase_cash < 0, "entrant cell is not the documented infeasible cell")
        expected = (-0.13404978853236, -0.06776975615932)[zz]
        require(math.isclose(cash, expected, rel_tol=0, abs_tol=2e-14), "actual renter cash differs")
        require(math.isclose(purchase_cash, cash / float(P.R_gross), rel_tol=0, abs_tol=2e-14),
                "purchase-income timing is not b + y/Rgross")
    require(float(P.property_tax_lump_sum_transfer) == 0 and float(P.transfer_floor_G0) == 0
            and float(P.transfer_floor_Gn) == 0 and not np.any(sd.gb_flat),
            "cash transfer contract differs")
    require(not np.any(sd.birth_dp) and not np.any(sd.birth_entry_grant)
            and not model.estate_receiver_active(P), "grant/estate-receiver contract differs")
    q = np.asarray(price, float).reshape(-1, 1, 1, 1)
    h = np.asarray(P.H_own, float).reshape(1, -1, 1, 1)
    downpayment = (1 - np.asarray(sd.phi_choice[:, 1:, :, :], float)) * q * h
    require(np.all(q >= 0) and np.all(downpayment >= 0), "actual purchase down payment is negative")
    return {"wealth": b, "gross_return": float(P.R_gross), "cells": rows,
            "minimum_actual_purchase_downpayment": float(np.min(downpayment)),
            "conclusion": "strict D=0 is infeasible for these inherited renters before positive rent, consumption, saving, or any nonnegative-price down payment"}


def array_census(value, prefix="solution", depth=4):
    out = {}
    if isinstance(value, np.ndarray):
        out[prefix] = value
    elif depth and (isinstance(value, dict) or hasattr(value, "__dict__")):
        mapping = value if isinstance(value, dict) else vars(value)
        for key, item in mapping.items():
            if not str(key).startswith("_"):
                out.update(array_census(item, prefix + "." + str(key), depth - 1))
    return out


def compare_solution_arrays(reference, candidate):
    old, new = array_census(reference), array_census(candidate)
    require(set(old) == set(new), "solution array attribute set differs")
    bad, records = {}, {}
    for name in sorted(old):
        a, b = old[name], new[name]
        same_shape = a.shape == b.shape
        exact = same_shape and np.array_equal(a, b)
        row = {"shape": list(a.shape), "candidate_shape": list(b.shape), "exact": bool(exact)}
        if same_shape and a.dtype.kind in "biufc" and b.dtype.kind in "biufc":
            finite = bool(np.isfinite(a).all() and np.isfinite(b).all())
            row["finite"] = finite
            row["max_abs"] = float(np.max(np.abs(a.astype(float) - b.astype(float)), initial=0)) if finite else None
        records[name] = row
        if not exact or row.get("finite") is False:
            bad[name] = row
    require(not bad, "exact control arrays differ: " + json.dumps(dict(list(bad.items())[:12]), sort_keys=True))
    return {"array_count": len(records), "all_exact": True, "arrays": records}


def smoke(output):
    progress(output, "authenticate_reference")
    _, (_, _, _, _, _, reference) = authenticate(output)
    progress(output, "import_overlay_namespace")
    _, kernels, model, origins = overlay_namespace()
    P = copy.deepcopy(reference["parameters"]); grid = np.asarray(reference["b_grid"]).copy()
    require(getattr(P, "unsecured_credit_limit", None) is None, "reference unexpectedly enables scalar credit")
    floors = {}
    for credit in (0.0, 2.5):
        P.unsecured_credit_limit = credit
        floors[str(credit)] = model.renter_borrowing_floor(P, np.array([-4.0, -1.0, 2.0]), 0).tolist()
    require(floors == {"0.0": [0.0, 0.0, 0.0], "2.5": [-2.5, -2.5, -2.5]}, "scalar helper floors differ")
    P.unsecured_credit_limit = 0.0
    result = dict(status="passed", lifecycle_solves=0, reference_label=LABEL, module_origins=origins,
                  source_hashes=PINS, driver_sha256=sha(__file__), plan_sha256=ACTIVE_PLAN_SHA,
                  helper_floors=floors, compiled_fixtures=compiled_fixtures(kernels),
                  strict_zero_cash_proof=cash_proof(P, grid, model, reference["solution"].p_eq),
                  strict_lifecycle_run_performed=False, checkpoint_sha256=PINS["checkpoint"])
    write(output / "receipt.json", result); progress(output, "complete", status="passed")


def control(output):
    deadline = float(os.environ["CASE_DEADLINE_EPOCH"])
    progress(output, "authenticate_reference")
    base, (manifest, _, objective, runtime, prepared, reference) = authenticate(output)
    require(time.time() < deadline, "control deadline reached during authentication")
    _, _, model, origins = overlay_namespace()
    P = copy.deepcopy(reference["parameters"]); grid = np.asarray(reference["b_grid"]).copy()
    require(getattr(P, "unsecured_credit_limit", None) is None, "control scalar must be absent/None")
    P.native_inherited_distribution_evidence_dir = str(output / "inherited_distribution_evidence")
    sd = model.precompute_shared(P, grid); price = np.asarray(reference["solution"].p_eq).copy()
    progress(output, "one_control_lifecycle", deadline_epoch=deadline)
    started = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd, verbose=False, fast_stats=False)
    require(time.time() < deadline, "control deadline exceeded")
    arrays = compare_solution_arrays(reference["solution"], sol)
    P._fert2_probs = sol.fert2_probs.copy()
    cal = prepared.rt["primitive"].pf.calendar
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction["stationary_post_fertility_nesting_l1"], 5e-9, "reconstruction")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol, evaluation=ev,
                  stationary_g_pre=pre, supply_rule=supply, demographic_seed=reference.get("demographic_seed"))
    gate = base.gates(packet, prepared, output, stationary=True)
    rt = prepared.rt
    fertility = {p: rt["observe_initial_fertility"](ev, P, age_projection=p)
                 for p in ("uniform_birth_time", "constant_post_cell")}
    housing = rt["observe_initial_housing_wealth"](ev, P, grid, sd, diagnostic_enabled=True,
        age_projection="uniform_within_age_cell", diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = rt["observe_recent_parent_flow"](ev, P, diagnostic_enabled=True,
        snapshot=rt["SNAPSHOT"], age_projection=rt["AGE_PROJECTION"], diagnostic_allow_residence_proxy=True,
        input_provenance={"case_id": "control", "reference_checkpoint_sha256": PINS["checkpoint"]})
    completed = float(rt["chain"].extract_moments(sol, P)["tfr"])
    fits = runtime.score_targets(objective, fertility, housing, recent["model_value"], completed)
    params = base.actual_parameters(prepared, P, grid)
    require(len(fits) == 14 and len(params) == 31, "complete fit/parameter readout unavailable")
    require(all(float(params[r["parameter"]]) == float(r["estimate"]) for r in manifest["full_parameter_table"]),
            "31-parameter control differs")
    base.exact_control(reference, packet, fits, manifest, prepared, output)
    base.table(output / "target_fit.csv", fits)
    base.table(output / "parameters.csv", manifest["full_parameter_table"])
    require(time.time() < deadline, "control reporting deadline exceeded")
    write(output / "receipt.json", dict(status="passed", lifecycle_solves=1,
          lifecycle_solve_seconds=time.monotonic() - started, module_origins=origins,
          driver_sha256=sha(__file__), plan_sha256=ACTIVE_PLAN_SHA,
          array_comparison=arrays, fit_rows=14, parameter_rows=31, gates=gate,
          reference_label=LABEL, checkpoint_sha256=PINS["checkpoint"]))
    progress(output, "complete", status="passed")


def main():
    global ACTIVE_PLAN_SHA
    parser = argparse.ArgumentParser(); parser.add_argument("mode", choices=("smoke", "control"))
    parser.add_argument("--output", type=Path, required=True); parser.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    require(sys.platform == "linux" and os.environ.get("SLURM_JOB_ID", "").isdigit(), "Torch Slurm required")
    require(all(os.environ.get(k) == "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")), "one-thread contract differs")
    plan = json.loads(args.plan.read_text())
    ACTIVE_PLAN_SHA = sha(args.plan)
    require(plan["schema"] == "block0506_fixed_credit_runtime_validation_v2" and plan["reference_manifest_sha256"] == PINS["manifest"], "plan contract differs")
    require(plan["driver_sha256"] == sha(__file__) and plan["base_authenticator_sha256"] == PINS["base_authenticator"], "driver/authenticator plan pin differs")
    require(plan["overlay_sha256"] == {k: PINS[k] for k in ("parameters.py", "solver.py", "kernels.py")}, "overlay plan pins differ")
    require(plan["threads"] == 1 and plan["modes"]["smoke"] == {"lifecycle_solves": 0, "minutes": 5, "memory_gib": 16}, "smoke budget differs")
    require(plan["modes"]["control"] == {"lifecycle_solves": 1, "case_seconds": 300, "total_seconds": 900, "memory_gib": 24}, "control budget differs")
    args.output.mkdir(parents=True, exist_ok=False)
    started = float(os.environ["LAUNCH_STARTED_EPOCH"])
    case_deadline = float(os.environ["CASE_DEADLINE_EPOCH"])
    total_deadline = float(os.environ["TOTAL_DEADLINE_EPOCH"])
    require(case_deadline == started + 300 and total_deadline == started + 900, "launcher clocks differ")
    write(args.output / "launch.json", dict(status="started", mode=args.mode,
          started_epoch=started, case_deadline_epoch=case_deadline, total_deadline_epoch=total_deadline,
          plan_sha256=ACTIVE_PLAN_SHA, slurm_job=os.environ["SLURM_JOB_ID"],
          maximum_lifecycle_solves=0 if args.mode == "smoke" else 1))
    try:
        (smoke if args.mode == "smoke" else control)(args.output)
    except Exception as error:
        write(args.output / "failure.json", {"status": "failed", "mode": args.mode, "error": repr(error), "time_epoch": time.time()})
        progress(args.output, "failed", error=repr(error)); raise


if __name__ == "__main__":
    main()
