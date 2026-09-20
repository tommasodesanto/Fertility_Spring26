"""Bounded adapter for the PSID persistent-plus-transitory income candidate.

The frozen initial probe and scorer remain untouched.  This adapter adds an
explicit, hash-pinned income contract and applies its override only after the
probe's parenthood candidate has been bound.  ``preflight`` performs no model
solve; ``pilot`` runs one cold native candidate and the frozen 13-row scorer.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import os
import runpy
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
PARAMETERS = ("beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
              "theta0", "theta1", "first_birth_fixed_cost", "h_P")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                  allow_nan=False).encode()).hexdigest()


def read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _artifact(plan: dict[str, Any], key: str) -> Path:
    path = Path(plan[key])
    if path.exists():
        return path
    staged = Path(__file__).resolve().parent / path.name
    if staged.exists():
        return staged
    return path


def _inventory(root: Path) -> set[str]:
    return {str(p.relative_to(root)) for p in (root / "code/model").rglob("*")
            if p.is_file() and "__pycache__" not in p.parts
            and p.suffix not in (".pyc", ".nbc", ".nbi")}


def validate_plan(plan: dict[str, Any], *, require_source: bool = True) -> dict[str, Any]:
    if plan.get("schema") != "e5f_income_candidate_calibration_v1":
        raise ValueError("wrong candidate calibration schema")
    if plan.get("objective_canonical_sha256") != OBJECTIVE:
        raise ValueError("frozen target/weight objective hash changed")
    if plan.get("free_parameters") != list(PARAMETERS) or plan.get("free_parameter_count") != 9:
        raise ValueError("candidate must retain all nine structural coordinates")
    if plan.get("fertility_normalization") != 2.1 or plan.get("payroll_tax") != 0.179:
        raise ValueError("normalization or fiscal contract changed")
    candidate = read(_artifact(plan, "candidate_json"))
    if candidate.get("status") != "diagnostic_candidate_only" or candidate.get("state_count") != 15:
        raise ValueError("candidate payload is not the declared 15-state diagnostic")
    if candidate.get("permanent_types") is not False:
        raise ValueError("permanent-income types must be explicitly disabled")
    if sha(_artifact(plan, "candidate_json")) != plan["candidate_json_sha256"]:
        raise ValueError("candidate JSON fingerprint mismatch")
    if fingerprint(candidate) != plan["candidate_payload_fingerprint"]:
        raise ValueError("candidate payload fingerprint mismatch")
    if sha(_artifact(plan, "constructor_path")) != plan["constructor_sha256"]:
        raise ValueError("candidate constructor fingerprint mismatch")
    if sha(_artifact(plan, "adapter_path")) != plan["adapter_sha256"]:
        raise ValueError("adapter fingerprint mismatch")
    source_sha = plan.get("source_sha256")
    if source_sha is None:
        manifest = read(Path(plan["source_manifest_path"])); provenance = manifest.get("source_provenance", {})
        source_sha = provenance.get("observation_snapshot", {}).get("source_sha256")
        if not isinstance(source_sha, dict):
            raise ValueError("source manifest lacks the frozen observation snapshot")
    if require_source:
        source = Path(plan["source_root"])
        if not source.is_dir():
            raise RuntimeError("BLOCKED: pinned 641-file native source root is unavailable")
        expected = set(source_sha)
        actual = _inventory(source)
        if actual != expected:
            raise RuntimeError(f"BLOCKED: source inventory mismatch; expected {len(expected)}, found {len(actual)}")
        for name, pin in source_sha.items():
            if sha(source / name) != pin:
                raise RuntimeError(f"BLOCKED: frozen source hash mismatch: {name}")
    return candidate


def candidate_overrides(candidate: dict[str, Any]) -> dict[str, Any]:
    """Build the solver payload from the audited constructor, never defaults."""
    # Caller installs the frozen source root before this function; retain the
    # active checkout only as a local construction fallback.
    if str(ROOT / "code/model") not in sys.path:
        sys.path.append(str(ROOT / "code/model"))
    annual = candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
    mapping = candidate["four_year_diagnostic_mapping"]
    rho = float(annual["rho_annual"])
    variance = float(annual["persistent_variance"])
    innovation_sd = (max((1.0 - rho * rho) * variance, 0.0)) ** 0.5
    constructor_path = Path(__file__).resolve().parent / "build_persistent_transitory_income_candidate.py"
    if not constructor_path.exists():
        constructor_path = ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py"
    constructor = _load(constructor_path,
                        "candidate_constructor")
    overrides, metadata = constructor.build_persistent_transitory_income_candidate(
        rho_annual=rho,
        persistent_innovation_sd_annual=innovation_sd,
        transitory_log_sd_period=float(mapping["transitory_log_sd_period"]),
        period_years=float(candidate["period_years"]), persistent_states=5,
    )
    if len(overrides["z_grid"]) != 15:
        raise ValueError("constructor did not produce exactly 15 income states")
    if candidate.get("state_count") != len(overrides["z_grid"]):
        raise ValueError("candidate metadata/state payload mismatch")
    return overrides


def _load(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def run_wrapped_probe(plan: dict[str, Any], contract: Path, output: Path) -> None:
    """Invoke the exact native probe with an in-memory income binding."""
    source = Path(plan["source_root"])
    sys.path[:0] = [str(source / "code/model/tools"), str(source / "code/model")]
    parent = __import__("e5f_parenthood_utility")
    original = parent.bind_parenthood_utility
    candidate = read(Path(plan["candidate_json"]))
    overrides = candidate_overrides(candidate)

    def bind(base: Any, structural: dict[str, Any]) -> Any:
        bound = original(base, structural)
        for key, value in overrides.items():
            setattr(bound, key, copy.deepcopy(value))
        from intergen_eqscale_seq_optimized.parameters import build_debt_caps
        bound = build_debt_caps(bound)
        bound.income_candidate_id = plan["candidate_id"]
        bound.income_candidate_fingerprint = plan["candidate_payload_fingerprint"]
        return bound

    parent.bind_parenthood_utility = bind
    probe = source / "code/model/tools/run_e5f_initial_revision_probe.py"
    old_argv = sys.argv
    try:
        sys.argv = [str(probe), "--contract", str(contract), "--contract-sha256", sha(contract),
                    "--case", "new_balanced", "--output", str(output)]
        runpy.run_path(str(probe), run_name="__main__")
    finally:
        sys.argv = old_argv


def run_scored_pilot(plan: dict[str, Any], output: Path, *, preflight_only: bool = False,
                     parameters: dict[str, float] | None = None,
                     initial_psi: float | None = None, repetitions: int = 1,
                     case_id: str = "income_candidate_pilot") -> dict[str, Any]:
    """Run one exact scored candidate through the frozen wrapper.

    The wrapper's child process is redirected only to this adapter's probe
    entry point; its validator, observer, scorer, gates, and 17-graph packet
    remain the frozen implementations.
    """
    template = Path(plan["template_dir"])
    inputs = template / "inputs" if (template / "inputs").is_dir() else template
    wrapper_path = Path(plan.get("wrapper_path", inputs / "run_scored_candidate.py"))
    for path in (Path(plan["initial_contract_path"]), Path(plan["run_contract_path"]), wrapper_path):
        if not path.exists():
            raise RuntimeError(f"BLOCKED: frozen native template missing {path}")
    sys.path[:0] = [str(Path(plan["source_root"]) / "code/model/tools"),
                    str(Path(plan["source_root"]) / "code/model")]
    initial = read(Path(plan["initial_contract_path"]))
    candidate = read(_artifact(plan, "candidate_json"))
    overrides = candidate_overrides(candidate)
    if repetitions not in (1, 2):
        raise ValueError("repetitions must be 1 or 2")
    structural = copy.deepcopy(parameters if parameters is not None else plan["pilot_parameters"])
    psi = float(plan["pilot_initial_psi"] if initial_psi is None else initial_psi)
    initial.update(case_id=case_id, repetitions=repetitions,
                   structural_candidate=structural,
                   initial_psi=psi,
                   income_candidate_id=plan["candidate_id"],
                   income_candidate_payload_fingerprint=plan["candidate_payload_fingerprint"],
                   income_candidate_json_sha256=plan["candidate_json_sha256"],
                   income_constructor_sha256=plan["constructor_sha256"],
                   income_adapter_sha256=plan["adapter_sha256"],
                   income_candidate_overrides={k: (v.tolist() if hasattr(v, "tolist") else v)
                                               for k, v in overrides.items()},
                   cold_solve_required=True)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    contract = copy.deepcopy(read(Path(plan["run_contract_path"])))
    contract["case_id"] = case_id
    contract["source_root"] = plan["source_root"]
    contract["initial_solve_contract"] = {"path": str(out / "initial_contract.json"),
                                          "sha256": ""}
    (out / "initial_contract.json").write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")
    contract["initial_solve_contract"]["sha256"] = sha(out / "initial_contract.json")
    (out / "run_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    wrapper = _load(wrapper_path, "frozen_candidate_wrapper")
    if preflight_only:
        wrapper.preflight(out / "run_contract.json", sha(out / "run_contract.json"))
        receipt = {"status": "candidate_wrapper_preflight_passed", "solves": 0,
                   "run_contract": str(out / "run_contract.json"),
                   "initial_contract": str(out / "initial_contract.json")}
        (out / "candidate_preflight.json").write_text(json.dumps(receipt, indent=2) + "\n")
        return receipt
    original_child = wrapper.run_child
    adapter = Path(__file__).resolve()

    def child(command, **kwargs):
        if any(str(x).endswith("run_e5f_initial_revision_probe.py") for x in command):
            command = [sys.executable, str(adapter), "--mode", "probe-child", "--plan", str(Path(plan["plan_path"])),
                       "--contract", command[command.index("--contract") + 1],
                       "--output", command[command.index("--output") + 1]]
        return original_child(command, **kwargs)

    wrapper.run_child = child
    wrapper.run(str(out / "run_contract.json"), sha(out / "run_contract.json"), out / "evaluation")
    summary = read(out / "evaluation/summary.json")
    return {"status": "pilot_completed", "summary": summary,
            "output": str(out / "evaluation")}


def assemble_pilot_contract(plan: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    """Return the reviewable native pilot contract without starting a solve."""
    manifest = read(Path(plan["source_manifest_path"]))
    source_sha = manifest["source_provenance"]["observation_snapshot"]["source_sha256"]
    return {
        "schema": "e5f_income_candidate_native_pilot_v1",
        "status": "assembled_review_required",
        "candidate_id": plan["candidate_id"],
        "candidate_payload_fingerprint": plan["candidate_payload_fingerprint"],
        "source_root": plan["source_root"],
        "source_manifest_fingerprint": fingerprint(source_sha),
        "objective_canonical_sha256": OBJECTIVE,
        "structural_parameter_count": 9,
        "structural_parameter_bounds": {"beta_annual": [0.94, 0.99]},
        "target_rows": 12,
        "normalization_target": 2.1,
        "payroll_tax": 0.179,
        "housing_supply_elasticity": 0.63,
        "income_state_count": candidate["state_count"],
        "cold_solve_required": True,
        "cached_checkpoint_policy": "reference-only; no policy/distribution reuse",
        "worker_count": 1,
        "thread_count": 1,
        "native_seconds": 420,
        "wrapper_seconds": 600,
        "exact_repetitions": 1,
        "pilot_output_contract": "full native checkpoint, 13-row target fit, 9-parameter table, gates",
        "launch_command": "python3 code/model/tools/run_e5f_income_candidate_calibration.py --mode pilot --plan calibration_plan.json --output pilot",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--mode", choices=("preflight", "pilot-contract", "prepare-only", "pilot", "probe-child"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--plan-path", dest="plan_path", type=Path)
    parser.add_argument("--parameters-json", type=Path)
    parser.add_argument("--initial-psi", type=float)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--case-id", default="income_candidate_pilot")
    args = parser.parse_args()
    plan = read(args.plan)
    plan["plan_path"] = str(args.plan.resolve())
    candidate = validate_plan(plan, require_source=args.mode == "pilot")
    if args.mode == "preflight":
        print(json.dumps({"status": "preflight_passed", "solves": 0,
                          "candidate_id": plan["candidate_id"], "state_count": candidate["state_count"]}))
        return
    if args.mode == "pilot-contract":
        print(json.dumps(assemble_pilot_contract(plan, candidate), indent=2, sort_keys=True))
        return
    if args.mode == "prepare-only":
        if args.output is None:
            raise ValueError("prepare-only requires --output")
        print(json.dumps(run_scored_pilot(plan, args.output, preflight_only=True), sort_keys=True))
        return
    if args.mode == "probe-child":
        if args.contract is None or args.output is None:
            raise ValueError("probe-child requires --contract and --output")
        run_wrapped_probe(plan, args.contract, args.output)
        return
    if args.output is None:
        raise ValueError("pilot requires --output")
    if args.mode == "pilot":
        if not plan.get("pilot_parameters") or "pilot_initial_psi" not in plan:
            raise ValueError("pilot parameters and normalized initial psi are required")
        parameters = read(args.parameters_json) if args.parameters_json else None
        print(json.dumps(run_scored_pilot(
            plan, args.output, parameters=parameters, initial_psi=args.initial_psi,
            repetitions=args.repetitions, case_id=args.case_id), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
