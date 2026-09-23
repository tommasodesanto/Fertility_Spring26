#!/usr/bin/env python3
"""Prepare four hash-pinned E5F utility templates; no model solve or submission.

The only generated adapter source change is the entry-rule guard.  The existing
earnings adapter still implements the inherited marginal and rank coupling.
"""
from __future__ import annotations

import argparse
import copy
import json
import shutil
import sys
import tempfile
from pathlib import Path

import prepare_e5f_earnings_entry_battery as battery
from prepare_e5f_earnings_wealth_run import prepare as prepare_earnings

SCHEMA = "e5f_utility_overnight_manifest_v1"
CELLS = ("B_floor", "B_shares", "D_floor", "D_shares")
COUNTS = {"B_floor": 18, "B_shares": 18, "D_floor": 2, "D_shares": 2}
MODES = {"floor": "floor_control", "shares": "child_dependent_shares"}
# Retained reference tilts; the overnight driver does not depend on the older
# unsubmitted zero-entry comparison launcher.
SHARED_NAMES = ("H0", "beta_annual", "chi", "first_birth_fixed_cost",
                "kappa_fert", "kappa_fert_continuation", "theta0", "theta1")
ANCHOR_TILTS = {"delta_alpha_jump": 0.020344, "delta_alpha": 0.016919}
ADAPTER_GUARD_OLD = "if plan.get('entry_specification', {}).get('rule') != 'zero_assets':"
ADAPTER_GUARD_NEW = "if plan.get('entry_specification', {}).get('rule') != 'fixed_reference_marginal':"


def rebase(obj, source: Path, target: Path):
    source, target = str(source.resolve()), str(target.resolve())
    if isinstance(obj, dict):
        return {k: rebase(v, Path(source), Path(target)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [rebase(v, Path(source), Path(target)) for v in obj]
    if isinstance(obj, str) and (obj == source or obj.startswith(source + "/")):
        return target + obj[len(source):]
    return obj


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise ValueError(f"adapter guard occurrence count {text.count(old)}; source review required")
    return text.replace(old, new, 1)


def make_local_plan(base: dict, bundle: Path, remote: Path, cell: str,
                    selected: dict, work: Path) -> dict:
    earnings_cell = cell[0]
    plan = battery._cell_plan(base, earnings_cell, selected, f"{cell}_template")
    spec = battery.CELLS[earnings_cell]
    n_persistent = 15 if earnings_cell == "B" else 7
    n_iid = 1 if earnings_cell == "B" else 3
    args = plan["income_specification"]["constructor_arguments"]
    args.update(n_persistent=n_persistent, n_iid=n_iid)
    plan["income_entry_battery"].update(persistent_states=n_persistent, iid_states=n_iid)
    plan["income_specification"]["resolution_note"] = (
        f"{n_persistent} persistent by {n_iid} iid states; experimental comparison, "
        "not household-grid convergence certification")
    plan["entry_specification"]["status"] = "authorized_utility_comparison_only"
    plan["preference_specification"] = {
        "mapping": MODES[cell.split("_")[1]], "author_decision": "approved_diagnostic"}
    if cell.endswith("shares"):
        structural = {k: float(v) for k, v in selected.items() if k != "h_P"}
        structural.update(ANCHOR_TILTS)
        bounds = dict(plan["parameter_bounds"])
        bounds.pop("h_P")
        bounds.update(delta_alpha_jump=[0.0, 0.25], delta_alpha=[0.0, 0.25])
        plan["parameter_bounds"] = bounds
    else:
        structural = dict(selected)
    plan["structural_parameters"] = dict(structural)
    plan["starting_structural_parameters"] = dict(structural)
    plan["comparison_draw"] = {
        "shared_parameters": {k: structural[k] for k in SHARED_NAMES},
        "preference_parameters": {k: structural[k] for k in structural if k not in SHARED_NAMES}}
    plan["cases"] = [{"id": f"{cell}_template", "arm": battery.ADAPTER_ARM,
                      "native_seconds": 6900, "wrapper_seconds": 7100,
                      "seconds": 7200, "repetitions": 1}]
    plan["case_id"] = f"{cell}_template"
    plan = rebase(plan, remote, bundle)
    run = rebase(battery.read(Path(plan["files"]["run_contract"]["path"])), remote, bundle)
    run_path = work / "run_contract.local.json"
    battery.write(run_path, run)
    plan["files"]["run_contract"] = {"path": str(run_path), "sha256": battery.sha256(run_path)}
    plan["source_root"] = str(bundle / "source")
    return plan


def build(plan_path: Path, selected_path: Path, bundle: Path, execution_root: Path,
          python_path: str, period_income_source: Path, adapter_source: Path,
          runner_source: Path, launcher_source: Path) -> dict:
    plan_path, selected_path, bundle = map(lambda p: p.resolve(), (plan_path, selected_path, bundle))
    execution_root = execution_root.resolve()
    if bundle.exists():
        raise ValueError("bundle must be new")
    source = battery.read(plan_path)
    selected = battery.read(selected_path)
    if selected.get("entry_specification", {}).get("rule") != "fixed_reference_marginal":
        raise ValueError("selected source must be inherited-entry case B")
    if source.get("objective_canonical_sha256") != battery.V5_OBJECTIVE:
        raise ValueError("source objective fingerprint changed")
    if set(selected["structural_parameters"]) != set(source["parameter_bounds"]):
        raise ValueError("selected structural dimension changed")
    if selected["structural_parameters"]["h_P"] > 2.3 or source["parameter_bounds"]["beta_annual"][1] > .99:
        raise ValueError("frozen upper bound changed")
    staged = copy.deepcopy(source)
    staged["files"]["period_income"] = {"path": str(period_income_source.resolve()),
                                         "sha256": battery.sha256(period_income_source)}
    with tempfile.TemporaryDirectory(prefix="e5f_utility_source_") as scratch:
        path = Path(scratch) / "plan.json"
        battery.write(path, staged)
        prepare_earnings(path, bundle, execution_root, python_path)
    base = battery.read(bundle / "plan.json")
    if base["source_manifest"]["file_count"] != 641:
        raise ValueError("source inventory changed")
    adapter = bundle / "tools" / adapter_source.name
    patched = replace_once(adapter_source.read_text(), ADAPTER_GUARD_OLD, ADAPTER_GUARD_NEW)
    patched = replace_once(patched,
        "raise ValueError('This matched experiment authorizes only explicit zero-asset entry')",
        "raise ValueError('This overnight comparison requires inherited heterogeneous entrant wealth')")
    adapter.write_text(patched)
    for supporting in (runner_source, launcher_source):
        shutil.copy2(supporting, bundle / "tools" / supporting.name)
    if str(bundle / "tools") not in sys.path:
        sys.path.insert(0, str(bundle / "tools"))
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("overnight_preference_adapter", adapter)
    assert spec and spec.loader
    pref_adapter = module_from_spec(spec)
    spec.loader.exec_module(pref_adapter)
    templates = {}
    selected_vector = {k: float(v) for k, v in selected["structural_parameters"].items()}
    for cell in CELLS:
        work = bundle / "inputs" / "utility_templates" / cell
        work.mkdir(parents=True)
        local = make_local_plan(base, bundle, execution_root, cell, selected_vector, work)
        local["files"]["adapter"] = {"path": str(adapter), "sha256": battery.sha256(adapter)}
        local["adapter_path"] = str(adapter)
        local["adapter_sha256"] = battery.sha256(adapter)
        generated = work / "generated"
        derived = pref_adapter.prepare_plan(local, generated)
        contract_path = Path(derived["files"]["preference_contract"]["path"])
        contract = battery.read(contract_path)
        contract["unchanged"] = ["power equivalence scale", "earnings process for declared cell",
                                 "inherited heterogeneous entrant wealth", "current-income purchase eligibility",
                                 "target rows and weights", "numerical gates"]
        battery.write(contract_path, contract)
        # Rebase generated contract links while the local files still exist.
        for key in ("run_contract", "preference_contract"):
            p = Path(derived["files"][key]["path"])
            battery.write(p, rebase(battery.read(p), bundle, execution_root))
        for key, item in derived["files"].items():
            p = Path(item["path"])
            if p.is_relative_to(bundle):
                item["path"] = str(execution_root / p.relative_to(bundle))
                item["sha256"] = battery.sha256(p)
        derived["source_root"] = str(execution_root / "source")
        derived["adapter_path"] = str(execution_root / "tools" / adapter.name)
        derived["adapter_sha256"] = battery.sha256(adapter)
        template = work / "template_plan.json"
        battery.write(template, derived)
        templates[cell] = {"path": str(execution_root / template.relative_to(bundle)),
            "sha256": battery.sha256(template),
            "objective_canonical_sha256": derived["objective_canonical_sha256"],
            "target_system_sha256": derived["target_system_sha256"],
            "bounds": derived["parameter_bounds"], "seed_parameters": derived["structural_parameters"],
            "entry_rule": derived["entry_specification"]["rule"],
            "income_grid": {"persistent": 15 if cell[0] == "B" else 7,
                            "iid": 1 if cell[0] == "B" else 3}}
    hashes = {v["target_system_sha256"] for v in templates.values()}
    if len(hashes) != 1:
        raise ValueError("target system differs across utility cells")
    manifest = {"schema": SCHEMA, "status": "ready", "templates": templates,
        "target_system_sha256": hashes.pop(), "source_file_count": 641,
        "source_manifest_sha256": battery.fingerprint(base["source_manifest"]["files"]),
        "frozen_v5_plan_sha256": battery.sha256(plan_path),
        "selected_B_plan_sha256": battery.sha256(selected_path),
        "entry_rule": "fixed_reference_marginal", "wealth_grid": {"nodes": 160, "upper": 3000},
        "workers": COUNTS, "total_workers": 40, "paired_seed_count": 18,
        "smoke_case_seconds": 7200, "case_seconds": 7200,
        "production_stop_utc": "2026-09-23T12:00:00Z",
        "verification_stop_utc": "2026-09-23T13:00:00Z",
        "report_deadline_utc": "2026-09-23T14:00:00Z",
        "maximum_nested_solves_per_objective": 8,
        "native_seconds_per_objective": 6900,
        "observed_seconds_per_objective": {"B_15x1": 1187, "D_7x3": 1699},
        "maximum_proposals_per_worker": 18,
        "maximum_objectives": {"smoke": 8, "production": 720, "verification": 8, "total": 736},
        "maximum_stationary_solves": 5888,
        "proposal_scale_design": {"paired_primary_chains": "seed IDs 1-6 local, 7-12 medium, 13-18 broad",
            "supplementary_chains": "seed ID 1 local, 2 medium",
            "initial_scales": {"local": 0.8, "medium": 1.6, "broad": 3.2},
            "selection": "each persistent worker proposes around its own best valid point; shared random innovations pair seed IDs but parameter vectors can diverge"},
        "adapter_source_change": {"original_sha256": battery.sha256(adapter_source),
            "bundled_sha256": battery.sha256(adapter), "old": ADAPTER_GUARD_OLD,
            "new": ADAPTER_GUARD_NEW, "scientific_meaning": "accept fixed reference wealth marginal and rank coupling; no zero fallback"},
        "runner": {"path": str(execution_root / "tools" / runner_source.name),
                   "sha256": battery.sha256(runner_source)},
        "launcher": {"path": str(execution_root / "tools" / launcher_source.name),
                     "sha256": battery.sha256(launcher_source)}}
    battery.write(bundle / "manifest.json", manifest)
    return {"bundle": str(bundle), "manifest": str(bundle / "manifest.json"),
            "manifest_sha256": battery.sha256(bundle / "manifest.json"),
            "target_system_sha256": manifest["target_system_sha256"], "templates": len(templates)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--selected-b-plan", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument("--period-income-source", type=Path, required=True)
    parser.add_argument("--adapter", type=Path, default=Path(__file__).with_name("run_e5f_preference_share_candidate.py"))
    parser.add_argument("--runner", type=Path, default=Path(__file__).with_name("run_e5f_utility_overnight.py"))
    parser.add_argument("--launcher", type=Path, default=Path(__file__).parents[2] / "cluster" / "submit_e5f_utility_overnight.sh")
    args = parser.parse_args()
    print(json.dumps(build(args.plan, args.selected_b_plan, args.bundle, args.execution_root,
                           args.python, args.period_income_source, args.adapter,
                           args.runner, args.launcher), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
