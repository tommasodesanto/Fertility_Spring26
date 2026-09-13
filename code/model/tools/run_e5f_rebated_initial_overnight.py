"""Run the pinned initial objective with an equal property-tax rebate.

The saved capped-beta packet supplies the source snapshot, candidate, objective,
weights, and original numerical gates.  This wrapper changes one closure only:
inside every initial stationary solve it roots the scalar equal household rebate,
while the saved routine continues to balance PAYGO pensions separately.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys
import time

SCHEMA = "e5f_rebated_initial_overnight_v1"
ANNUAL_PROPERTY_TAX = 0.01
PERIOD_PROPERTY_TAX = 0.04
REBATE_RELATIVE_TOLERANCE = 1e-6
PENSION_RELATIVE_TOLERANCE = 1e-6
ROOMS_TARGET = 0.7202462623815278
NORMALIZATION_TARGET = 2.1
DEFAULT_WORKERS = 6
MAXIMUM_POOL = 18
MAXIMUM_WALL_SECONDS = 3 * 60 * 60
MAXIMUM_REBATE_SOLVES_PER_NORMALIZATION_CALL = 20
_ACTUAL_SOLVES = 0
_ACTUAL_SOLVE_LIMIT = None


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(base, path):
    path = Path(path)
    return path.resolve() if path.is_absolute() else (Path(base) / path).resolve()


def load_module(name, path):
    path = Path(path).resolve()
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def saved_packet(template):
    """Resolve templates from the saved resume receipt, never a default old target."""
    template = Path(template).resolve()
    plan_path = template / "plan_capped_beta_099.json"
    plan = read(plan_path)
    if sha(template / "run_capped_beta.py") != plan["controller_sha256"]:
        raise ValueError("Saved capped-beta controller fingerprint changed")
    seed_score = Path(plan["resume_score_path"]).resolve()
    if sha(seed_score) != plan["resume_score_sha256"]:
        raise ValueError("Saved seed score fingerprint changed")
    evaluation = seed_score.parent.parent
    case = evaluation.parent
    initial_path = case / "initial_contract.json"
    run_path = case / "run_contract.json"
    initial = read(initial_path)
    run = read(run_path)
    base = run_path.parent
    scorer_path = resolve(base, run["scorer"]["path"])
    scored_wrapper = scorer_path.with_name("run_scored_candidate.py")
    if sha(scored_wrapper) != run["wrapper_sha256"]:
        raise ValueError("Saved scored-candidate wrapper fingerprint changed")
    objective_path = resolve(base, run["working_objective"]["path"])
    objective = read(objective_path)
    source_root = Path(run["source_root"]).resolve()
    return dict(plan=plan, plan_path=plan_path, seed_score_path=seed_score,
                seed_score=read(seed_score), initial_path=initial_path,
                initial=initial, run_path=run_path, run=run,
                scored_wrapper=scored_wrapper, objective=objective,
                objective_path=objective_path, source_root=source_root)


def validate_scientific_contract(packet):
    initial, objective, run = packet["initial"], packet["objective"], packet["run"]
    rows = objective["target_rows"]
    scored = [row for row in rows if row["role"] != "normalization_separate_from_scored_objective"]
    normalizations = [row for row in rows if row["restriction_id"] == "initial_normalization"]
    rooms = [row for row in rows if row["restriction_id"] == "first_birth_rooms"]
    restrictions = {row["parameter"]: dict(row) for row in objective["parameter_restrictions"]}
    if len(scored) != 12 or len(normalizations) != 1 or normalizations[0]["target"] != NORMALIZATION_TARGET:
        raise ValueError("Saved full twelve-row objective or separate 2.1 normalization changed")
    if len(rooms) != 1 or rooms[0]["target"] != ROOMS_TARGET:
        raise ValueError("Saved first-birth rooms target changed")
    if len(restrictions) != 9 or set(initial["structural_candidate"]) != set(restrictions):
        raise ValueError("Saved candidate is not the complete nine-coordinate problem")
    restrictions["beta_annual"]["upper"] = 0.99
    beta = float(initial["structural_candidate"]["beta_annual"])
    if not restrictions["beta_annual"]["lower"] <= beta <= 0.99:
        raise ValueError("Saved annual beta violates the approved [.94,.99] cap")
    if (initial["fertility_normalization"] != NORMALIZATION_TARGET
            or initial["payroll_tax"] != 0.179 or not initial["normalize"]
            or not initial["observe_early"] or run["source_root"] != str(packet["source_root"])):
        raise ValueError("Saved initial normalization, PAYGO, observer, or source contract changed")
    return restrictions


def _relative(residual, *flows):
    return abs(float(residual)) / max(*(abs(float(x)) for x in flows), 1e-12)


def solve_equal_rebate(original, *, model, parameters, b_grid, initial_prices,
                       payroll_tax, marginal_tolerance, fiscal_tolerance):
    """Root equal rebates outside the unchanged housing/PAYGO equilibrium."""
    global _ACTUAL_SOLVES

    def evaluate(transfer):
        global _ACTUAL_SOLVES
        _ACTUAL_SOLVES += 1
        if _ACTUAL_SOLVE_LIMIT is not None and _ACTUAL_SOLVES > _ACTUAL_SOLVE_LIMIT:
            raise TimeoutError("Explicit actual stationary-solve ceiling exhausted")
        candidate = copy.deepcopy(parameters)
        if not math.isclose(float(candidate.tau_H), PERIOD_PROPERTY_TAX, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Initial property-tax primitive is not the saved 1% annual/4% period rate")
        candidate.property_tax_lump_sum_transfer = float(round(float(transfer), 12))
        return original(model=model, parameters=candidate, b_grid=b_grid,
                        initial_prices=initial_prices, payroll_tax=payroll_tax,
                        marginal_tolerance=marginal_tolerance,
                        fiscal_tolerance=fiscal_tolerance)

    def ledger(result):
        sol, P, _, pension = result
        revenue = float(sol.property_tax_revenue)
        outlays = float(sol.property_tax_transfer_outlays)
        residual = float(sol.property_tax_budget_residual)
        mass = float(outlays / P.property_tax_lump_sum_transfer) if P.property_tax_lump_sum_transfer else float(sol.g.sum())
        if not all(math.isfinite(x) for x in (revenue, outlays, residual, mass)) or revenue <= 0 or mass <= 0:
            raise RuntimeError("Nonfinite or nonpositive property-tax accounting")
        return residual, revenue, outlays, mass, pension

    lower = 0.0
    low_result = evaluate(lower)
    f_low, revenue, _, mass, _ = ledger(low_result)
    if f_low < 0:
        raise RuntimeError("Equal-rebate root has negative revenue residual at zero transfer")
    upper = max(1e-8, 1.25 * revenue / mass)
    best = low_result
    high_result = evaluate(upper)
    f_high, _, _, _, _ = ledger(high_result)
    if abs(f_high) < abs(f_low):
        best = high_result
    while f_high > 0 and upper < 8.0:
        upper *= 2.0
        high_result = evaluate(upper)
        f_high, _, _, _, _ = ledger(high_result)
        if abs(f_high) < abs(ledger(best)[0]):
            best = high_result
    if f_high > 0:
        raise RuntimeError("Could not bracket the equal property-tax rebate")
    for _ in range(MAXIMUM_REBATE_SOLVES_PER_NORMALIZATION_CALL - 2):
        residual, revenue, outlays, mass, _ = ledger(best)
        if _relative(residual, revenue, outlays) <= REBATE_RELATIVE_TOLERANCE:
            break
        proposal = lower - f_low * (upper - lower) / (f_high - f_low)
        width = upper - lower
        proposal = min(upper - 0.1 * width, max(lower + 0.1 * width, proposal))
        trial = evaluate(proposal)
        trial_f = ledger(trial)[0]
        if trial_f > 0:
            lower, f_low = proposal, trial_f
        else:
            upper, f_high = proposal, trial_f
        if abs(trial_f) < abs(ledger(best)[0]):
            best = trial
    residual, revenue, outlays, mass, pension = ledger(best)
    sol, P, price, pension = best
    relative = _relative(residual, revenue, outlays)
    if relative > REBATE_RELATIVE_TOLERANCE:
        raise RuntimeError(f"Equal-rebate accounting gate failed: {relative:.9g}")
    if not math.isclose(outlays, float(P.property_tax_lump_sum_transfer) * mass,
                        rel_tol=1e-10, abs_tol=1e-10):
        raise RuntimeError("Equal per-household rebate outlays do not match household mass")
    receipt = dict(pension)
    receipt["property_tax_rebate"] = dict(
        fiscal_convention="balanced_budget_equal_rebate",
        annual_property_tax_rate=ANNUAL_PROPERTY_TAX,
        period_property_tax_rate=PERIOD_PROPERTY_TAX,
        transfer_period_units=float(P.property_tax_lump_sum_transfer),
        property_tax_revenue=revenue, transfer_outlays=outlays,
        budget_residual=residual, relative_budget_residual=relative,
        actual_stationary_solves=_ACTUAL_SOLVES)
    return sol, P, price, receipt


def raw_mode(source_root, driver_path, driver_args):
    global _ACTUAL_SOLVES, _ACTUAL_SOLVE_LIMIT
    source_root = Path(source_root).resolve()
    driver_path = Path(driver_path).resolve()
    if not driver_path.is_relative_to(source_root):
        raise ValueError("Raw driver escapes the saved source root")
    sys.path[:0] = [str(source_root / "code/model/tools"), str(source_root / "code/model")]
    driver = load_module("rebated_initial_raw_driver", driver_path)
    try:
        contract_path = driver_args[driver_args.index("--contract") + 1]
    except (ValueError, IndexError) as exc:
        raise ValueError("Raw rebated run requires the saved initial contract") from exc
    contract = read(contract_path)
    _ACTUAL_SOLVES = 0
    _ACTUAL_SOLVE_LIMIT = (int(contract["maximum_stationary_solves_per_repetition"])
                           * int(contract["repetitions"])
                           * MAXIMUM_REBATE_SOLVES_PER_NORMALIZATION_CALL)
    original = driver.solve_balanced_initial_equilibrium
    driver.solve_balanced_initial_equilibrium = lambda **kwargs: solve_equal_rebate(original, **kwargs)
    sys.argv = [str(driver_path), *driver_args]
    driver.main()


def accounting_receipts(evaluation_dir, source_root):
    source_root = Path(source_root).resolve()
    sys.path[:0] = [str(source_root / "code/model/tools"), str(source_root / "code/model")]
    from intergen_eqscale_seq_optimized import solver as model
    receipts = []
    raw = Path(evaluation_dir) / "raw"
    for repetition in sorted(raw.glob("repetition_*")):
        summary = read(repetition / "summary.json")
        with gzip.open(repetition / "initial_state.pkl.gz", "rb") as stream:
            packet = pickle.load(stream)
        P, evaluation = packet["parameters"], packet["evaluation"]
        revenue = float(model.property_tax_revenue_from_distribution(
            evaluation.g_current, evaluation.policy.hR_pol, evaluation.policy.price, P))
        outlays = float(P.property_tax_lump_sum_transfer) * float(evaluation.g_current.sum())
        residual = revenue - outlays
        pension = summary["fiscal"]["actual_accounts"]
        pension_relative = abs(float(pension["pension_budget_residual"])) / max(
            abs(float(pension["payroll_tax_revenue"])), abs(float(pension["pension_outlays"])), 1e-12)
        rebate_relative = _relative(residual, revenue, outlays)
        if (not math.isclose(float(P.tau_H), PERIOD_PROPERTY_TAX, rel_tol=0, abs_tol=1e-12)
                or float(P.property_tax_lump_sum_transfer) <= 0
                or rebate_relative > REBATE_RELATIVE_TOLERANCE
                or pension_relative > PENSION_RELATIVE_TOLERANCE):
            raise RuntimeError("Final independent pension/rebate accounting gate failed")
        receipts.append(dict(repetition=repetition.name,
            checkpoint=str(repetition / "initial_state.pkl.gz"),
            checkpoint_sha256=summary["checkpoint_sha256"],
            transfer_period_units=float(P.property_tax_lump_sum_transfer),
            property_tax_revenue=revenue, rebate_outlays=outlays,
            rebate_residual=residual, rebate_relative_gap=rebate_relative,
            pension_relative_gap=pension_relative))
    if not receipts:
        raise RuntimeError("No rebated checkpoint was produced")
    return receipts


def candidate(packet, restrictions, item, output):
    controller = load_module("saved_capped_beta_controller", Path(packet["plan_path"]).parent / "run_capped_beta.py")
    controller.validate_proposal(item, restrictions)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    initial = copy.deepcopy(packet["initial"])
    initial.update(structural_candidate=item["parameters"], initial_psi=item["initial_psi"],
                   repetitions=item.get("repetitions", 1),
                   equal_property_tax_rebate=True,
                   rebate_relative_tolerance=REBATE_RELATIVE_TOLERANCE,
                   fiscal_closure="equal property-tax rebate; PAYGO pension separately balanced")
    write(output / "initial_contract.json", initial)
    run = copy.deepcopy(packet["run"])
    old_base = packet["run_path"].parent
    for key in ("working_objective", "scorer", "validator"):
        run[key]["path"] = str(resolve(old_base, packet["run"][key]["path"]))
    for key, entry in run.get("objective_source_files", {}).items():
        entry["path"] = str(resolve(old_base, packet["run"]["objective_source_files"][key]["path"]))
    run["case_id"] = item["case_id"]
    run["initial_solve_contract"] = dict(path=str(output / "initial_contract.json"),
                                          sha256=sha(output / "initial_contract.json"))
    write(output / "run_contract.json", run)
    scored = load_module("saved_rebated_scored_wrapper", packet["scored_wrapper"])
    original_run_child = scored.run_child

    def intercept(command, **kwargs):
        if len(command) > 1 and Path(command[1]).name == "run_e5f_initial_revision_probe.py":
            command = [sys.executable, str(Path(__file__).resolve()), "raw",
                       "--source-root", str(packet["source_root"]),
                       "--driver", str(command[1]), "--", *command[2:]]
        return original_run_child(command, **kwargs)

    scored.run_child = intercept
    evaluation = output / "evaluation"
    try:
        summary = scored.run(output / "run_contract.json", sha(output / "run_contract.json"), evaluation)
        receipts = accounting_receipts(evaluation, packet["source_root"])
        write(evaluation / "rebate_accounting.json", dict(status="verified", receipts=receipts))
        score = read(evaluation / "scored_repetition_01/score.json")
        controller.validate_score(score, item, restrictions)
        result = dict(case_id=item["case_id"], status="verified", loss=score["loss"],
                      output=str(evaluation), proposal=item, score=score,
                      receipt_status=summary["status"], accounting=receipts)
        if item.get("repetitions") == 2:
            second = read(evaluation / "scored_repetition_02/score.json")
            controller.validate_score(second, item, restrictions)
            result["second_signature_equal"] = (summary.get("exact_loss_equality") is True
                and controller.numeric_signature(score) == controller.numeric_signature(second)
                and receipts[0] == {**receipts[1], "repetition": receipts[0]["repetition"],
                                    "checkpoint": receipts[0]["checkpoint"],
                                    "checkpoint_sha256": receipts[0]["checkpoint_sha256"]})
        return result
    except Exception as exc:
        write(output / "branch_failure.json", dict(status="hard_error", error_type=type(exc).__name__, error=str(exc)))
        return dict(case_id=item["case_id"], status="failed", output=str(evaluation),
                    proposal=item, failure_detail=str(exc))


def write_selected_tables(out, result):
    out = Path(out)
    score = result["score"]
    for name, rows in (("selected_target_fit.csv", score["target_fit"]),
                       ("selected_parameters.csv", score["parameters"])):
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (out / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
    write(out / "selected_checkpoint.json", result["accounting"][-1])


def run_smoke(template, output):
    packet = saved_packet(template)
    restrictions = validate_scientific_contract(packet)
    proposal = copy.deepcopy(packet["plan"]["resume_proposal"])
    proposal["case_id"] = "rebated_initial_smoke"
    proposal["repetitions"] = 1
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    write(out / "launch_contract.json", dict(schema=SCHEMA, mode="smoke",
        template=str(Path(template).resolve()), source_root=str(packet["source_root"]),
        saved_plan_sha256=sha(packet["plan_path"]), workers=1,
        wall_seconds=2100, annual_property_tax_rate=ANNUAL_PROPERTY_TAX,
        complete_scored_moments=12, free_parameters=9, beta_upper=0.99,
        normalization=NORMALIZATION_TARGET, rooms_target=ROOMS_TARGET))
    write(out / "latest_completed.json", dict(status="running", case_id=proposal["case_id"]))
    result = candidate(packet, restrictions, proposal, out / "case")
    if result["status"] != "verified":
        write(out / "summary.json", dict(status="failed", result=result)); raise SystemExit(2)
    write_selected_tables(out, result)
    write(out / "best_so_far.json", dict(status="verified", case_id=result["case_id"],
                                           loss=result["loss"], output=result["output"]))
    write(out / "latest_completed.json", dict(status="verified", case_id=result["case_id"],
                                                loss=result["loss"], output=result["output"]))
    summary = dict(status="verified_rebated_initial_smoke", loss=result["loss"],
                   checkpoint=result["accounting"][-1], source_root=str(packet["source_root"]))
    write(out / "summary.json", summary)
    print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    raw = sub.add_parser("raw")
    raw.add_argument("--source-root", type=Path, required=True)
    raw.add_argument("--driver", type=Path, required=True)
    raw.add_argument("driver_args", nargs=argparse.REMAINDER)
    run = sub.add_parser("run")
    run.add_argument("--template", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--smoke", action="store_true", required=True)
    args = parser.parse_args()
    if args.mode == "raw":
        values = args.driver_args[1:] if args.driver_args[:1] == ["--"] else args.driver_args
        raw_mode(args.source_root, args.driver, values)
    else:
        run_smoke(args.template, args.output)


if __name__ == "__main__":
    main()
