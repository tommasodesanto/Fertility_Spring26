#!/usr/bin/env python3
"""Prepare a pinned four-cell earnings-process by entrant-wealth battery.

This packages the frozen V5 source and emits smoke and worker plans. It does
not solve the model, copy files to a cluster, use SSH, or submit jobs.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import random
import shutil
import tempfile
from pathlib import Path
from typing import Any

SCHEMA = "e5f_earnings_entry_battery_manifest_v1"
V5_OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
ADAPTER_ARM = "literature_income_purchase"
EXPECTED_SOURCE_COUNT = 641
WORKERS_PER_CELL = 10
PROPOSALS_PER_WORKER = 6
CASE_SECONDS = 3200
NATIVE_SECONDS = 3100
WRAPPER_SECONDS = 3150
TOTAL_SECONDS = 3600

CELLS: dict[str, dict[str, Any]] = {
    "A": {
        "label": "single_zero_assets", "earnings_specification": "single_persistent_ar1",
        "entry_rule": "zero_assets",
        "income": {"rho_period": 0.7345934905942886,
                   "persistent_innovation_sd_period": 0.4838308245314463,
                   "transitory_sd_period": 0.0, "n_persistent": 7, "n_iid": 1},
    },
    "B": {
        "label": "single_fixed_reference_entry", "earnings_specification": "single_persistent_ar1",
        "entry_rule": "fixed_reference_marginal",
        "income": {"rho_period": 0.7345934905942886,
                   "persistent_innovation_sd_period": 0.4838308245314463,
                   "transitory_sd_period": 0.0, "n_persistent": 7, "n_iid": 1},
    },
    "C": {
        "label": "two_shock_zero_assets", "earnings_specification": "persistent_plus_transitory",
        "entry_rule": "zero_assets",
        "income": {"rho_period": 0.7761446698586146,
                   "persistent_innovation_sd_period": 0.43743545312066157,
                   "transitory_sd_period": 0.1649906105257306,
                   "n_persistent": 7, "n_iid": 3},
    },
    "D": {
        "label": "two_shock_fixed_reference_entry", "earnings_specification": "persistent_plus_transitory",
        "entry_rule": "fixed_reference_marginal",
        "income": {"rho_period": 0.7761446698586146,
                   "persistent_innovation_sd_period": 0.43743545312066157,
                   "transitory_sd_period": 0.1649906105257306,
                   "n_persistent": 7, "n_iid": 3},
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fingerprint(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def read(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def write(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def validate_cells() -> None:
    combos = {(v["earnings_specification"], v["entry_rule"]) for v in CELLS.values()}
    expected = {
        ("single_persistent_ar1", "zero_assets"),
        ("single_persistent_ar1", "fixed_reference_marginal"),
        ("persistent_plus_transitory", "zero_assets"),
        ("persistent_plus_transitory", "fixed_reference_marginal"),
    }
    if len(CELLS) != 4 or combos != expected:
        raise ValueError("cells must form the declared 2-by-2 earnings/entry design")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pinned search controller: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def candidate_pool(plan: dict[str, Any], size: int = 60) -> list[dict[str, float]]:
    """One anchor plus native, deterministic local proposals shared by cells."""
    if size != WORKERS_PER_CELL * PROPOSALS_PER_WORKER:
        raise ValueError("the production pool must contain exactly 60 matched proposals")
    anchor = {k: float(v) for k, v in plan["starting_structural_parameters"].items()}
    if set(anchor) != set(plan["parameter_bounds"]):
        raise ValueError("anchor and frozen nine-parameter bound contract differ")
    controller = load_module(Path(plan["files"]["search_controller"]["path"]),
                             "pinned_e5f_battery_search_controller")
    seen = {fingerprint(anchor)}
    proposals = controller.proposal_batch(
        anchor, plan, random.Random(20260922), size - 1, seen, 0.35)
    points = [anchor, *proposals]
    if len(points) != size or len({fingerprint(x) for x in points}) != size:
        raise ValueError("native proposal_batch did not produce 60 unique shared points")
    for point in points:
        for name, value in point.items():
            low, high = map(float, plan["parameter_bounds"][name])
            if not low <= float(value) <= high:
                raise ValueError(f"proposal outside frozen bound: {name}")
    return points


def smoke_parameters(plan: dict[str, Any]) -> dict[str, float]:
    point = {k: float(v) for k, v in plan["starting_structural_parameters"].items()}
    point["beta_annual"] = 0.985
    low, high = map(float, plan["parameter_bounds"]["beta_annual"])
    if not low <= point["beta_annual"] <= high:
        raise ValueError("beta=.985 smoke point is outside the frozen bound")
    return point


def _cell_plan(base: dict[str, Any], cell_id: str, structural: dict[str, float],
               case_id: str) -> dict[str, Any]:
    spec = CELLS[cell_id]
    plan = copy.deepcopy(base)
    plan["battery_cell_id"] = cell_id
    plan["battery_cell_label"] = spec["label"]
    plan["structural_parameters"] = dict(structural)
    plan["entry_specification"] = {
        **copy.deepcopy(base.get("entry_specification", {})),
        "rule": spec["entry_rule"],
        "classification": "externally_fixed_diagnostic",
        "status": "authorized_four_cell_diagnostic_only",
        "basis": ("Diagnostic comparison: zero wealth at age 18, independent of current income."
                  if spec["entry_rule"] == "zero_assets" else
                  "Diagnostic comparison: retain the inherited checkpoint-P model wealth marginal and rank coupling. This is not a denominator-aligned empirical estimate."),
        "source_contract": ("zero assets at age 18, independent of current income" if spec["entry_rule"] == "zero_assets"
                            else "inherited checkpoint-P model wealth marginal and rank coupling; not denominator-aligned empirical estimate"),
        "required_runtime_receipts": [
            "source age/sample and ratio field/denominator", "model-income conversion",
            "frontier-censored mass", "marginal mean/L1/max gap", "rank-coupling declaration",
            "persistent/iid independence",
        ],
        "reference_marginal_receipt": {
            "source_checkpoint": "inherited checkpoint P, ages 18-24, 5 income-ratio nodes",
            "empirical_ratio": "NETWORTH2R / INCFAMR (nonhousing net worth / family income); denominator differs from model-income conversion",
            "baseline_mean": 0.18651967924681825,
            "prezero_marginal_gap_max_abs": 1.11e-16,
            "historical_replay_frontier_censored_mass": 0.0,
            "historical_grid_clip_mass": 0.0,
            "historical_raw_wealth_range": [-8.21846, 11.47594],
            "historical_grid_range": [-12.0, 30.0],
            "frontier_censored_mass": "must be measured and reported by each new run; historical zero is not assumed",
            "model_income_conversion": "wealth = empirical ratio * P.income[0,0] * old_z / (4*(1-0.179)); P.income is four-year after-payroll-tax model earnings",
        },
    }
    plan["entry_specification"].pop("entry_wealth", None)
    prior = copy.deepcopy(base["income_specification"])
    prior.update({
        "author_decision": "approved_diagnostic",
        "authorization_basis": "Author authorized a four-cell earnings-by-entry-wealth battery with a common age profile; both entry rules remain diagnostic and no paper specification is adopted.",
        "mapping": "direct_period",
        "constructor_arguments": copy.deepcopy(spec["income"]),
        "max_relative_discrete_level_covariance_error": 0.15,
        "resolution_classification": "coarse_resolution_diagnostic_not_certified",
        "resolution_note": "Seven persistent states are used in every cell. Resolution is not certified; repeat any selected point with the 15-state persistent grid before interpretation or adoption.",
    })
    prior["entry_wealth"] = ("Zero assets at age 18 independent of current income; diagnostic cell"
                              if spec["entry_rule"] == "zero_assets" else
                              "Inherited checkpoint-P model wealth marginal and rank coupling; diagnostic cell, not denominator-aligned empirical estimate")
    plan["income_specification"] = prior
    plan["income_entry_battery"] = {
        "earnings_specification": spec["earnings_specification"],
        "income_mapping": "direct_period", "persistent_states": spec["income"]["n_persistent"],
        "iid_states": spec["income"]["n_iid"], "entry_rule": spec["entry_rule"],
        "income_parameters": copy.deepcopy(spec["income"]),
        "coarse_covariance_gate": 0.15, "old_v5_gate_changed": False,
        "grid_resolution_certified": False,
    }
    plan["cases"] = [{"id": case_id, "arm": ADAPTER_ARM,
                      "native_seconds": NATIVE_SECONDS,
                      "wrapper_seconds": WRAPPER_SECONDS,
                      "seconds": CASE_SECONDS, "repetitions": 1}]
    plan["search_arm"] = ADAPTER_ARM
    plan["case_id"] = case_id
    plan["adapter_path"] = plan["files"]["adapter"]["path"]
    plan["adapter_sha256"] = plan["files"]["adapter"]["sha256"]
    if float(plan["initial_psi"]) != float(base["initial_psi"]):
        raise ValueError("initial psi differs across battery cells")
    return plan


def _candidate_record(bundle: Path, execution_root: Path, base: dict[str, Any],
                      cell_id: str, structural: dict[str, float], case_id: str,
                      proposal_id: str) -> dict[str, Any]:
    local = bundle / "inputs" / "earnings_entry_battery" / case_id / "plan.json"
    plan = _cell_plan(base, cell_id, structural, case_id)
    write(local, plan)
    remote = execution_root.resolve() / local.relative_to(bundle.resolve())
    return {"proposal_id": proposal_id, "case_id": case_id,
            "plan_path": str(remote), "plan_sha256": sha256(local),
            "objective_canonical_sha256": plan["objective_canonical_sha256"],
            "shared_parameters": dict(structural)}


def build(plan_path: Path, bundle: Path, execution_root: Path, python_path: str,
          period_income_source: Path, minimum_next_proposal_seconds: float,
          runner_source: Path | None = None, launcher_source: Path | None = None) -> dict[str, Any]:
    validate_cells()
    plan_path, bundle, execution_root = plan_path.resolve(), bundle.resolve(), execution_root.resolve()
    period_income_source = period_income_source.resolve()
    if bundle.exists():
        raise ValueError(f"bundle must be a new immutable directory: {bundle}")
    if not 1 <= float(minimum_next_proposal_seconds) <= CASE_SECONDS:
        raise ValueError("minimum_next_proposal_seconds must be in [1, 3200]")
    source_plan = read(plan_path)
    if source_plan.get("status") != "validated_by_lead" or source_plan.get("config_fixed") is not True:
        raise ValueError("input must be the validated, fixed-configuration V5 plan")
    if source_plan.get("objective_canonical_sha256") != V5_OBJECTIVE:
        raise ValueError("frozen V5 objective fingerprint mismatch")
    wealth = source_plan.get("wealth_grid_specification", {})
    if (int(wealth.get("original_nodes", -1)) + int(wealth.get("extra_nodes", 0)) != 160
            or float(wealth.get("upper", float("nan"))) != 3000.0):
        raise ValueError("battery must retain the frozen V5 160-node, 3000 upper grid")
    if not period_income_source.is_file():
        raise ValueError(f"missing pure-AR-capable constructor: {period_income_source}")
    # The added runtime constructor is separately pinned. Frozen source files,
    # including the original 641-file source inventory, remain untouched.
    staged_plan = copy.deepcopy(source_plan)
    staged_plan["files"]["period_income"] = {
        "path": str(period_income_source), "sha256": sha256(period_income_source)}
    staged_plan["bundle_role"] = "shared_runtime_payload_only"
    with tempfile.TemporaryDirectory(prefix="e5f_earnings_entry_plan_") as temp:
        temp_plan = Path(temp) / "plan.json"
        write(temp_plan, staged_plan)
        from prepare_e5f_earnings_wealth_run import prepare
        prepare(temp_plan, bundle, execution_root, python_path)
    base = read(bundle / "plan.json")
    if base.get("source_manifest", {}).get("file_count") != EXPECTED_SOURCE_COUNT:
        raise ValueError("frozen V5 source inventory must contain exactly 641 files")
    if base["files"]["period_income"]["sha256"] != sha256(period_income_source):
        raise ValueError("pure-AR constructor hash changed during packaging")
    if base["objective_canonical_sha256"] != V5_OBJECTIVE:
        raise ValueError("V5 objective changed during packaging")
    copied_support = []
    for source in (runner_source, launcher_source):
        if source is None:
            continue
        source = source.resolve()
        destination = bundle / "tools" / source.name
        shutil.copy2(source, destination)
        copied_support.append({"path": str(destination), "sha256": sha256(destination)})
    objective = read(source_plan["files"]["run_contract"]["path"])
    objective_file = Path(objective["working_objective"]["path"])
    objective_payload = read(objective_file)
    target_payload = {k: v for k, v in objective_payload.items() if k != "parameter_restrictions"}
    target_sha = fingerprint(target_payload)
    base["target_system_sha256"] = target_sha
    source_sha = fingerprint(base["source_manifest"]["files"])
    points = candidate_pool(source_plan)
    beta_smoke = smoke_parameters(source_plan)
    production: list[dict[str, Any]] = []
    smoke: list[dict[str, Any]] = []
    for cell_id in CELLS:
        smoke_case = f"{cell_id}_smoke_beta0985"
        smoke_proposal = _candidate_record(bundle, execution_root, base, cell_id,
            beta_smoke, smoke_case, "smoke_beta0985")
        smoke.append({"cell_id": cell_id, "task_id": 1, "proposals": [smoke_proposal]})
        for worker in range(WORKERS_PER_CELL):
            proposals = []
            for within in range(PROPOSALS_PER_WORKER):
                draw_index = worker * PROPOSALS_PER_WORKER + within
                case_id = f"{cell_id}_worker{worker+1:02d}_proposal{within+1:02d}"
                proposals.append(_candidate_record(bundle, execution_root, base, cell_id,
                    points[draw_index], case_id, f"proposal_{within+1:02d}"))
            production.append({"cell_id": cell_id, "task_id": worker + 1,
                               "proposals": proposals})
    cell_receipts = {}
    for cell_id, spec in CELLS.items():
        cell_receipts[cell_id] = {
            "label": spec["label"], "earnings_specification": spec["earnings_specification"],
            "persistent_states": spec["income"]["n_persistent"],
            "iid_states": spec["income"]["n_iid"], "income_mapping": "direct_period",
            "entry_rule": spec["entry_rule"], "income_parameters": spec["income"],
            "source_inventory_sha256": source_sha, "source_file_count": EXPECTED_SOURCE_COUNT,
            "period_income_runtime_sha256": sha256(period_income_source),
            "income_resolution": "7 persistent states; 1 iid state for single AR(1), 3 for two-shock; exploratory only",
            "coarse_discrete_level_covariance_error_gate": 0.15,
            "old_v5_gate_changed": False,
        }
    manifest = {
        "schema": SCHEMA, "status": "ready", "objective_canonical_sha256": V5_OBJECTIVE,
        "target_system_sha256": target_sha, "source_inventory_sha256": source_sha,
        "source_file_count": EXPECTED_SOURCE_COUNT, "scored_moment_count": 13,
        "parameter_row_count": 17, "workers_per_cell": WORKERS_PER_CELL,
        "worker_count": 40, "proposals_per_worker": PROPOSALS_PER_WORKER,
        "case_seconds": CASE_SECONDS, "native_seconds": NATIVE_SECONDS,
        "wrapper_seconds": WRAPPER_SECONDS, "total_seconds": TOTAL_SECONDS,
        "minimum_next_proposal_seconds": float(minimum_next_proposal_seconds),
        "maximum_production_objectives": 240,
        "maximum_production_stationary_solves": 1920,
        "maximum_smoke_objectives": 4, "maximum_smoke_stationary_solves": 32,
        "maximum_total_stationary_solves": 1952,
        "estimated_seconds_per_objective": {
            "reference_45_income_states_seconds": [2000, 2800],
            "seven_state_single_process_seconds": "unmeasured; smoke establishes runtime",
            "21_state_two_shock_seconds": "unmeasured; smoke establishes runtime"},
        "estimated_single_case_rss_gib": "unmeasured; verify within worker allocation during preflight",
        "maximum_income_level_covariance_error": 0.15,
        "resolution_certified": False,
        "selected_candidate_repeat_requirement": "Repeat the selected point on 15 persistent states before interpretation or adoption; the seven-state battery does not certify grid resolution.",
        "wealth_grid": {"nodes": 160, "upper": 3000,
                        "source_plan_sha256": sha256(plan_path),
                        "unchanged_from_v5": True},
        "initial_psi": float(base["initial_psi"]),
        "proposal_generation": {"seed": 20260922, "native_method": "proposal_batch",
            "scale": 0.35, "anchor_in_pool": True,
            "pairing": "Each worker/proposal parameter vector is shared across A/B/C/D."},
        "cells": cell_receipts,
        "smoke": {"cases": smoke}, "production": {"cases": production},
        "economic_scope": "Four-cell earnings-process by entry-wealth diagnostic. Targets, weights, source inventory, nine structural coordinates, initial psi, housing preferences, and solver gates are unchanged. No adoption claim.",
    }
    if runner_source is not None:
        manifest["runner"] = {"path": str(execution_root / "tools" / runner_source.name),
                              "sha256": sha256(runner_source)}
    if launcher_source is not None:
        manifest["launcher"] = {"path": str(execution_root / "tools" / launcher_source.name),
                                "sha256": sha256(launcher_source)}
    write(bundle / "manifest.json", manifest)
    hash_manifest = read(bundle / "hash_manifest.json")
    hash_manifest["files"].append({"label": "pure_ar_period_income_constructor",
        "path": str(bundle / "tools" / period_income_source.name),
        "source": str(period_income_source), "sha256": sha256(bundle / "tools" / period_income_source.name)})
    for item in copied_support:
        hash_manifest["files"].append({"label": Path(item["path"]).name,
            "path": item["path"], "source": item["path"], "sha256": item["sha256"]})
    for row in [*smoke, *production]:
        for proposal in row["proposals"]:
            local_plan = Path(proposal["plan_path"].replace(str(execution_root), str(bundle)))
            hash_manifest["files"].append({"label": "case_plan:" + proposal["case_id"],
                "path": str(local_plan), "source": str(local_plan), "sha256": sha256(local_plan)})
    hash_manifest["battery_manifest_sha256"] = sha256(bundle / "manifest.json")
    write(bundle / "hash_manifest.json", hash_manifest)
    write(bundle / "battery_receipt.json", {
        "schema": "e5f_earnings_entry_battery_receipt_v1", "status": "prepared_hash_verified",
        "manifest": str(bundle / "manifest.json"), "manifest_sha256": sha256(bundle / "manifest.json"),
        "production_worker_rows": len(production), "smoke_worker_rows": len(smoke),
        "production_plan_count": sum(len(x["proposals"]) for x in production),
        "smoke_plan_count": sum(len(x["proposals"]) for x in smoke),
        "source_file_count": EXPECTED_SOURCE_COUNT, "source_inventory_sha256": source_sha,
        "objective_canonical_sha256": V5_OBJECTIVE, "target_system_sha256": target_sha,
        "runtime_constructor_sha256": sha256(period_income_source),
        "grid_resolution_certified": False,
    })
    return {"status": "prepared_hash_verified", "bundle": str(bundle),
            "manifest": str(bundle / "manifest.json"), "production_plan_count": 240,
            "smoke_plan_count": 4, "source_file_count": EXPECTED_SOURCE_COUNT}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument("--period-income-source", type=Path,
                        default=Path(__file__).with_name("build_period_earnings_process.py"))
    parser.add_argument("--minimum-next-proposal-seconds", type=float, required=True,
                        help="Fixed minimum reserve; production also reserves 1.1 times its verified smoke and latest proposal runtimes")
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--launcher", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.plan, args.bundle, args.execution_root, args.python,
                           args.period_income_source, args.minimum_next_proposal_seconds,
                           args.runner, args.launcher), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
