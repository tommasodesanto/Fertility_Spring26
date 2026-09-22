"""Three explicitly contracted diagnostic arms on the frozen native scorer.

No model runs by import. Source files are immutable; economic runtime changes
are separately hashed and recorded, never represented as unchanged behavior.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import runpy
import sys
from pathlib import Path

import numpy as np

import build_literature_period_income as income
import build_period_earnings_process as period_income
import e5f_earnings_wealth_contract as accounting

OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
ARMS = ("reference", "literature_income", "literature_income_purchase")


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_plan(plan):
    if plan.get("objective_canonical_sha256") != OBJECTIVE:
        raise ValueError("unexpected complete target system")
    required = {"adapter", "accounting", "income", "wrapper", "initial_contract", "run_contract"}
    if not required.issubset(plan["files"]):
        raise ValueError("missing adapter/input pins")
    for label, item in plan["files"].items():
        if digest(item["path"]) != item["sha256"]:
            raise ValueError(f"changed contracted file: {label}")
    for label, actual in (("adapter", __file__), ("accounting", accounting.__file__),
                          ("income", income.__file__)):
        if Path(plan["files"][label]["path"]).resolve() != Path(actual).resolve():
            raise ValueError(f"loaded a different {label} module")
    if plan["income_specification"]["mapping"] == "direct_period":
        item = plan["files"].get("period_income")
        if item is None or Path(item["path"]).resolve() != Path(period_income.__file__).resolve():
            raise ValueError("direct-period constructor path/hash must be pinned")


def candidate(plan):
    spec = plan["income_specification"]
    if spec.get("author_decision") != "approved_diagnostic":
        raise ValueError("income-period approximation decision is pending")
    if spec["mapping"] == "conventional_endpoint":
        overrides, metadata = income.build_conventional_endpoint_income_adapter(**spec["constructor_arguments"])
        continuous = np.asarray(metadata["continuous_endpoint_level_covariances"])
    elif spec["mapping"] == "direct_period":
        overrides, metadata = period_income.build_period_earnings_process(**spec["constructor_arguments"])
        continuous = np.asarray(metadata["continuous_level_covariances"])
    else:
        raise ValueError("unsupported income-period mapping; no automatic fallback")
    discrete = np.asarray(metadata["discrete_level_covariances"])
    relative_error = float(np.max(np.abs(discrete / continuous - 1.)))
    if (not metadata["stationary"] or not metadata["iid_transition_independent"]
            or abs(float(overrides["z_weights"] @ overrides["z_grid"]) - 1.) > 1e-12
            or relative_error > spec["max_relative_discrete_level_covariance_error"]):
        raise ValueError("income grid failed the explicit distribution-approximation gate")
    metadata["maximum_relative_discrete_level_covariance_error"] = relative_error
    return overrides, metadata


def case_for(plan, arm):
    matches = [c for c in plan["cases"] if c["arm"] == arm]
    if arm not in ARMS or len(matches) != 1:
        raise ValueError("arm must have exactly one explicit case")
    return matches[0]


def run_probe(plan, initial_path, output):
    initial = read(initial_path)
    arm = initial["earnings_wealth_arm"]
    source = Path(plan["source_root"])
    sys.path[:0] = [str(source / "code/model/tools"), str(source / "code/model")]
    if arm != "reference":
        overrides, metadata = candidate(plan)
        from intergen_eqscale_seq_optimized import solver as model
        from intergen_eqscale_seq_optimized.parameters import build_debt_caps
        import e5f_parenthood_utility as parent
        import run_e5f_matched_pf_smoke as primitive
        bind_original = parent.bind_parenthood_utility
        accounting.install_fixed_entry(model)

        def bind(base, structural):
            bound = bind_original(base, structural)
            grid = model.make_grid(bound)
            components = metadata["components"]
            entry_rule = plan.get("entry_specification", {}).get("rule", "fixed_reference_marginal")
            if entry_rule == "fixed_reference_marginal":
                conditional, entry_receipt = accounting.rank_coupled_entry(
                    model, bound, grid, np.asarray(components["persistent_weights"]),
                    np.asarray(components["iid_weights"]))
            elif entry_rule == "zero_assets":
                zero = np.flatnonzero(np.asarray(grid) == 0.)
                if len(zero) != 1:
                    raise ValueError("zero-asset entry requires an exact zero grid node")
                conditional = np.zeros((len(grid), len(overrides["z_grid"])))
                conditional[zero[0], :] = 1.
                entry_receipt = dict(rule="zero_assets", entry_age=float(bound.age_start),
                    externally_fixed=True, empirical_joint_distribution_estimated=False,
                    candidate_wealth_mean=0., candidate_wealth_marginal=conditional[:, 0].tolist())
            else:
                raise ValueError("unsupported entry-wealth rule; no fallback")
            for key, value in overrides.items():
                setattr(bound, key, copy.deepcopy(value))
            bound.fixed_reference_entry_grid = grid.copy()
            bound.fixed_reference_entry_conditional = conditional
            bound.earnings_wealth_arm = arm
            bound = build_debt_caps(bound)
            write(output.parent / "income_process.json", metadata)
            write(output.parent / "entry_wealth.json", entry_receipt)
            return bound

        parent.bind_parenthood_utility = bind
        if arm == "literature_income_purchase":
            accounting.install_purchase_income(model, output.parent / "purchase_source.diff")
            budget_original = primitive.dated_budget
            audits = []

            def budget(evaluation, P, shared, grid, rent):
                native = budget_original(evaluation, P, shared, grid, rent)
                receipt = accounting.audit_purchase_accounting(evaluation, P, shared, grid, model)
                audits.append(receipt)
                write(output.parent / "purchase_accounting.json", audits)
                return {**native, "purchase_accounting": receipt}

            primitive.dated_budget = budget
    probe = source / "code/model/tools/run_e5f_initial_revision_probe.py"
    sys.argv = [str(probe), "--contract", str(initial_path), "--contract-sha256", digest(initial_path),
                "--case", "new_balanced", "--output", str(output)]
    runpy.run_path(str(probe), run_name="__main__")


def run_case(plan_path, plan, arm, output, repetitions, preflight=False):
    case = case_for(plan, arm)
    if repetitions != case["repetitions"]:
        raise ValueError("requested repetitions differ from contracted case")
    if arm != "reference":
        candidate(plan)  # Fail before any household solve if the decision is pending.
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    initial = copy.deepcopy(read(plan["files"]["initial_contract"]["path"]))
    initial.update(case_id=case["id"], repetitions=repetitions,
                   structural_candidate=plan["structural_parameters"], initial_psi=plan["initial_psi"],
                   earnings_wealth_arm=arm, earnings_wealth_plan_sha256=digest(plan_path),
                   earnings_wealth_adapter_files=plan["files"],
                   scope="Fixed-parameter earnings/entry/purchase diagnostic; no calibration adoption",
                   seconds=case["native_seconds"])
    write(output / "initial_contract.json", initial)
    contract = copy.deepcopy(read(plan["files"]["run_contract"]["path"]))
    contract.update(case_id=case["id"], source_root=plan["source_root"],
                    initial_solve_contract={"path": str(output / "initial_contract.json"),
                                            "sha256": digest(output / "initial_contract.json")},
                    seconds=case["wrapper_seconds"])
    write(output / "run_contract.json", contract)
    runtime_contract = dict(arm=arm, plan_sha256=digest(plan_path),
        frozen_source_manifest_sha256=hashlib.sha256(json.dumps(initial["source_sha256"],
            sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        additional_runtime_files={k: plan["files"][k] for k in ("adapter", "accounting", "income", "period_income") if k in plan["files"]},
        changed_economic_objects=([] if arm == "reference" else
            ["income process", "entry wealth: " + plan.get("entry_specification", {}).get("rule", "fixed_reference_marginal_rank_coupling")] +
            (["current-income purchase eligibility", "ordinary transaction wealth map",
              "end-of-period owner mortgage floor"] if arm.endswith("purchase") else [])),
        interpretation="Native score source fingerprints identify the frozen base; this additional contract identifies economic runtime changes.")
    write(output / "runtime_contract.json", runtime_contract)
    wrapper = load(plan["files"]["wrapper"]["path"], "frozen_earnings_wealth_wrapper")
    if preflight:
        wrapper.preflight(output / "run_contract.json", digest(output / "run_contract.json"))
        receipt = {"status": "preflight_passed", "arm": arm, "household_solves": 0}
        write(output / "preflight.json", receipt)
        return receipt
    original_child = wrapper.run_child

    def child(command, **kwargs):
        if any(str(x).endswith("run_e5f_initial_revision_probe.py") for x in command):
            command = [sys.executable, str(Path(__file__).resolve()), "--mode", "probe-child",
                       "--plan", str(plan_path), "--contract", str(output / "initial_contract.json"),
                       "--output", command[command.index("--output") + 1]]
        return original_child(command, **kwargs)

    wrapper.run_child = child
    result = wrapper.run(str(output / "run_contract.json"), digest(output / "run_contract.json"),
                         output / "evaluation")
    if arm == "reference" and abs(result["loss"] - plan["reference_loss"]) > 1e-8:
        raise RuntimeError("reference objective failed to reproduce; subsequent arms blocked")
    if arm.endswith("purchase"):
        audits = read(output / "evaluation/purchase_accounting.json")
        if len(audits) != repetitions:
            raise RuntimeError("missing purchase-accounting repetition audit")
        if repetitions == 2 and audits[0] != audits[1]:
            raise RuntimeError("purchase-accounting repetitions differ")
        runtime_contract["generated_solver_sha256"] = digest(output / "evaluation/purchase_source.generated.py")
        runtime_contract["generated_tenure_support_sha256"] = digest(output / "evaluation/purchase_source.tenure.py")
        runtime_contract["source_diff_sha256"] = digest(output / "evaluation/purchase_source.diff")
    runtime_contract["status"] = "verified_diagnostic"
    write(output / "runtime_contract.json", runtime_contract)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("pilot", "probe-child"), default="pilot")
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument("--repetitions", type=int)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    plan = read(plan_path)
    verify_plan(plan)
    if args.preflight:
        receipts = [run_case(plan_path, plan, c["arm"], args.output / c["id"],
                             c["repetitions"], True) for c in plan["cases"]]
        write(args.output / "receipt.json", {"status": "preflight_passed", "cases": receipts, "solves": 0})
    elif args.mode == "probe-child":
        run_probe(plan, args.contract, args.output)
    else:
        print(json.dumps(run_case(plan_path, plan, args.arm, args.output, args.repetitions)))


if __name__ == "__main__":
    main()
