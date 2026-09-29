#!/usr/bin/env python3
"""Fail-closed preparation controller for the e5f preference estimator.

This module deliberately uses only the standard library.  In particular, plan
inspection and plan creation never import the numerical model.  Numerical work
is delegated to an explicitly supplied Torch-only smoke harness or to the
already-gated estimator driver.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import time


IDENTIFICATION_SHA256 = "68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf"
REFERENCE_MANIFEST_SHA256 = "147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4"
SOURCE_NAMES = (
    "run_e5f_preference_transition.py", "e5f_four_shock_acceleration.py",
    "e5f_exact_policy_cache.py", "e5f_social_security_root.py",
    "e5f_matched_pf_path_root.py", "e5f_ssj_scaled_step_root.py",
    "e5f_ssj_toeplitz_jacobian.py", "e5f_preference_shock_fit.py",
    "run_e5f_preference_estimation.py",
)


def fail(message):
    raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def source_pins(source_dir):
    source_dir = Path(source_dir)
    pins = {name: sha(source_dir / name) for name in SOURCE_NAMES}
    if set(pins) != set(SOURCE_NAMES):
        fail("Incomplete estimator source pins")
    return pins


def target_contract(blocks, annual):
    with Path(blocks).open(newline="") as stream:
        block_rows = list(csv.DictReader(stream))
    with Path(annual).open(newline="") as stream:
        annual_rows = list(csv.DictReader(stream))
    result = []
    for period, year in enumerate((2007, 2011, 2015, 2019)):
        block = [row for row in block_rows if int(row["decision_year"]) == year]
        sample = [row for row in annual_rows if year < int(row["year"]) <= year + 4]
        if len(block) != 1 or len(sample) != 4 or len({row["year"] for row in sample}) != 4:
            fail("Exactly four unique annual observations per fertility window required")
        row = block[0]
        value = sum(float(item["period_tfr_births_per_woman"]) for item in sample) / 4
        if (int(row["birth_year_start"]) != year + 1 or int(row["birth_year_end"]) != year + 4
                or abs(value - float(row["period_tfr_arithmetic_mean"])) > 1e-12):
            fail("Historical target construction changed")
        if any(item["status"] != "verified_published_final" for item in sample):
            fail("Unverified annual data")
        result.append({
            "moment": f"period_tfr_{year + 1}_{year + 4}", "decision_year": year,
            "period": period, "birth_year_start": year + 1, "birth_year_end": year + 4,
            "target": value, "estimator": "equal-weight arithmetic mean of four published NCHS annual TFRs",
            "sample": "US annual birth-registration rates; published female exposure",
            "fixed_effects": "not applicable", "clustering": "not applicable", "uncertainty": None,
            "sources": sorted({item["source_url"] for item in sample}),
        })
    return {
        "schema": "retained_nchs_four_windows_v1", "rows": result,
        "blocks": {"path": str(Path(blocks)), "sha256": sha(blocks)},
        "annual": {"path": str(Path(annual)), "sha256": sha(annual)},
        "model_measurement": "period_tfr_topcode_adjusted; sum of four-year birth flows divided by own-age household mass",
        "caveat": "Retained household-rate analogue, not literal female-exposure TFR; no new age-support target adopted",
    }


def maximum_policy_calls(kind, seed_horizon, fit_evaluations, endpoint_evaluations, path_evaluations, horizons):
    """Mirror the driver’s deliberately conservative call ceiling, including the measured seed."""
    stages = 4 if kind == "four_successive" else 1
    return 10 * seed_horizon + stages * fit_evaluations * (
        endpoint_evaluations + 2 + sum(2 * horizon * path_evaluations for horizon in horizons)
    ) + 8


def require_positive_budgets(budget):
    for name, value in budget.items():
        if (not isinstance(value, (int, float)) or isinstance(value, bool)
                or not math.isfinite(value) or value <= 0):
            fail(f"Finite positive budget required: {name}")


def require_positive_cache_budget(value):
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        fail("Cache memory budget must be a positive integer byte count")


def check_readiness(receipt_path, pins):
    receipt_path = Path(receipt_path)
    receipt = read_json(receipt_path)
    if receipt.get("status") != "PASS" or not receipt.get("tests_passed"):
        fail("Readiness receipt is not PASS")
    if not receipt.get("native_endpoint_and_fertility_smoke_passed"):
        fail("Native endpoint/path/handoff readiness did not pass")
    if receipt.get("preference_changes") or receipt.get("historical_fit"):
        fail("Readiness must be unchanged-preference and test-only")
    if receipt.get("estimator_sources") != pins:
        fail("Readiness source pins do not match the staged estimator")
    return {"path": str(receipt_path), "sha256": sha(receipt_path)}


def load_staged_driver(source_dir):
    """Load the staged stdlib-only driver without retaining import path changes."""
    source_dir = Path(source_dir).resolve()
    names = ("run_e5f_preference_transition", "run_e5f_preference_estimation")
    prior = {name: sys.modules.get(name) for name in names}
    try:
        for name in names:
            spec = importlib.util.spec_from_file_location(name, source_dir / (name + ".py"))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
        return sys.modules[names[-1]]
    finally:
        for name in names:
            if prior[name] is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = prior[name]


def preflight_plan(plan, source_dir):
    """Check a completed launch-shaped plan with the staged driver's gates."""
    driver = load_staged_driver(source_dir)
    launch_plan = copy.deepcopy(plan)
    launch_plan["execution_enabled"] = True
    driver.validate_plan(launch_plan, launching=True)
    driver.validate_launch_inputs(launch_plan)


def plan_from_args(args):
    require_positive_cache_budget(args.cache_max_bytes)
    if args.kind not in ("one_permanent", "four_successive"):
        fail("kind must be one_permanent or four_successive")
    horizons = list(args.horizon)
    if len(horizons) < 2 or any(h < 6 for h in horizons) or any(a >= b for a, b in zip(horizons, horizons[1:])):
        fail("At least two increasing horizons of six or more periods required")
    if args.seed_horizon < 6 or not 0 < args.perturbed_date < args.seed_horizon - 1 or not 0 < args.log_step < .01:
        fail("Invalid measured-Jacobian settings")
    if args.reference_manifest_sha256 != REFERENCE_MANIFEST_SHA256:
        fail("Wrong fixed block0506 reference manifest")
    if args.fit_evaluations < 5 or args.endpoint_evaluations < 2 or args.path_evaluations < 2:
        fail("Plan must reserve the driver's minimum root and replay evaluations")
    pins = source_pins(args.source_dir)
    readiness = check_readiness(args.readiness, pins)
    seed_receipt = None
    if args.jacobian_receipt is not None:
        seed_receipt = {"path": str(args.jacobian_receipt), "sha256": sha(args.jacobian_receipt)}
    budget = {
        "total_seconds": args.total_seconds, "candidate_seconds": args.candidate_seconds,
        "endpoint_seconds": args.endpoint_seconds, "mapping_seconds": args.mapping_seconds,
        "path_seconds": args.path_seconds, "jacobian_seconds": args.jacobian_seconds,
    }
    require_positive_budgets(budget)
    calls = maximum_policy_calls(args.kind, args.seed_horizon, args.fit_evaluations,
                                 args.endpoint_evaluations, args.path_evaluations, horizons)
    budget["maximum_policy_calls"] = calls
    plan = {
        "schema": "block0506_surprise_estimation_v1", "execution_enabled": bool(args.enable_execution),
        "kind": args.kind, "reference_manifest_sha256": args.reference_manifest_sha256,
        "housing": args.housing, "credit": "saved_reference_unchanged",
        "expectations": "current_shock_permanent_until_next_surprise", "fiscal": "fixed_payroll_tax_endogenous_pension",
        "outside_entry": 0, "retention": 1, "property_rebate": 0, "initial_level": "saved_reference_level",
        "search_bound_ratios": [.01, 2.], "target_contract": target_contract(args.empirical_blocks, args.annual),
        "source_pins": pins, "readiness_receipt": readiness,
        "acceleration": {"seed_horizon": args.seed_horizon, "perturbed_date": args.perturbed_date,
                           "log_step": args.log_step, "seed_receipt": seed_receipt},
        "horizons": horizons, "horizon_comparison_periods": 1 if args.kind == "four_successive" else 4,
        "horizon_relative_tolerance": 1e-3, "budget": budget,
        "fit": {"max_evaluations": args.fit_evaluations, "log_difference_step": .01, "fertility_tolerance": .005,
                "max_log_step": .15, "damping": .7, "max_condition_number": 1e8, "worsening_factor": 1.5,
                "reproduction_tolerance": 1e-8},
        "endpoint": {"max_evaluations": args.endpoint_evaluations, "price_bound_ratios": [.05, 20.], "max_log_step": .15,
                     "damping": .7, "slope": 1., "renewal_tolerance": 1e-6, "reproduction_tolerance": 1e-10},
        "path": {"max_evaluations": args.path_evaluations, "price_bound_ratios": [.05, 20.], "pension_bound_ratios": [.05, 20.],
                 "cache_max_bytes": args.cache_max_bytes, "market_tolerance": 2e-4, "fiscal_tolerance": 1e-6,
                 "market_slope": 1., "fiscal_slope": 1., "max_log_step": .15, "damping": .7,
                 "final_reproduction_tolerance": 1e-10,
                 "terminal_tolerances": {key: 1e-3 for key in ("population_relative_gap", "normalized_distribution_l1", "birth_queue_maximum_relative_gap", "asset_price_relative_gap", "renter_price_relative_gap", "psi_absolute_gap")},
                 "raw_queue_relative_tolerance": 1e-3},
        "labels": {"estate": "provisional inherited estate settlement remains outstanding"},
    }
    preflight_plan(plan, args.source_dir)
    return plan


def create_plan(args):
    write_json(args.output, plan_from_args(args))


def run_harness(args):
    """Run a lead-staged test harness; it owns native imports and writes the receipt."""
    command = [args.python, str(args.harness), *args.harness_arg]
    started = time.monotonic()
    result = subprocess.run(command, timeout=args.wall_seconds)
    if result.returncode:
        fail(f"Smoke harness failed with exit status {result.returncode}")
    pins = source_pins(args.source_dir)
    check_readiness(args.readiness, pins)
    write_json(args.summary, {"status": "PASS", "mode": "native-readiness", "seconds": time.monotonic() - started,
                              "readiness": {"path": str(args.readiness), "sha256": sha(args.readiness)}, "source_pins": pins,
                              "historical_fit": False, "preference_changes": False})


def run_pure_tests(args):
    result = subprocess.run([args.python, "-m", "unittest", "-v", *args.test_module], cwd=args.test_cwd,
                            timeout=args.wall_seconds)
    if result.returncode:
        fail(f"Pure controller tests failed with exit status {result.returncode}")


def execute(args):
    plan = read_json(args.plan)
    if not plan.get("execution_enabled"):
        fail("Execution plan is disabled")
    if sha(args.plan) != args.plan_sha256:
        fail("Pinned execution plan changed")
    pins = source_pins(args.source_dir)
    readiness = check_readiness(plan["readiness_receipt"]["path"], pins)
    if readiness != plan["readiness_receipt"] or plan.get("source_pins") != pins:
        fail("Plan/readiness/source identity mismatch")
    command = [args.python, str(Path(args.source_dir) / "run_e5f_preference_estimation.py"),
               "--plan", str(args.plan), "--plan-sha256", args.plan_sha256, "--output", str(args.output), "--execute"]
    raise SystemExit(subprocess.run(command).returncode)


def parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--python", default=sys.executable)
    common.add_argument("--wall-seconds", type=float, default=1200.)
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="mode", required=True)
    pure = sub.add_parser("pure-tests", parents=[common]); pure.add_argument("--test-cwd", required=True)
    pure.add_argument("--test-module", action="append", required=True); pure.set_defaults(func=run_pure_tests)
    ready = sub.add_parser("native-readiness", parents=[common]); ready.add_argument("--harness", type=Path, required=True)
    ready.add_argument("--harness-arg", action="append", default=[]); ready.add_argument("--source-dir", type=Path, required=True)
    ready.add_argument("--readiness", type=Path, required=True); ready.add_argument("--summary", type=Path, required=True); ready.set_defaults(func=run_harness)
    create = sub.add_parser("create-plan")
    create.add_argument("--kind", required=True); create.add_argument("--source-dir", type=Path, required=True)
    create.add_argument("--readiness", type=Path, required=True); create.add_argument("--empirical-blocks", type=Path, required=True)
    create.add_argument("--annual", type=Path, required=True); create.add_argument("--output", type=Path, required=True)
    create.add_argument("--reference-manifest-sha256", required=True); create.add_argument("--housing", choices=("fixed_stock", "elastic_reference"), default="fixed_stock")
    create.add_argument("--horizon", type=int, action="append", required=True); create.add_argument("--seed-horizon", type=int, required=True)
    create.add_argument("--perturbed-date", type=int, required=True); create.add_argument("--log-step", type=float, default=1e-5)
    create.add_argument("--cache-max-bytes", type=int, default=2 * 1024 ** 3)
    create.add_argument("--jacobian-receipt", type=Path)
    create.add_argument("--fit-evaluations", type=int, required=True); create.add_argument("--endpoint-evaluations", type=int, required=True); create.add_argument("--path-evaluations", type=int, required=True)
    for name in ("total", "candidate", "endpoint", "mapping", "path", "jacobian"):
        create.add_argument("--" + name + "-seconds", type=float, required=True)
    create.add_argument("--enable-execution", action="store_true")
    create.set_defaults(func=create_plan)
    launch = sub.add_parser("execute", parents=[common]); launch.add_argument("--plan", type=Path, required=True)
    launch.add_argument("--plan-sha256", required=True); launch.add_argument("--source-dir", type=Path, required=True)
    launch.add_argument("--output", type=Path, required=True); launch.set_defaults(func=execute)
    return root


if __name__ == "__main__":
    arguments = parser().parse_args()
    arguments.func(arguments)
