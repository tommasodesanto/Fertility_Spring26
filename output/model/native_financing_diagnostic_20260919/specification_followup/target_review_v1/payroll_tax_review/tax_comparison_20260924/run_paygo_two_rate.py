#!/usr/bin/env python3
"""Two fixed-parameter stationary PAYGO tax diagnostics on the selected B-floor checkpoint.

The dry run exercises both case slots and output bookkeeping without a model solve.
Production requires the original selected checkpoint, frozen source, and external
1800-second process timeout. This is an experiment, not recalibration.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import os
import pickle
import sys
import threading
import time
from pathlib import Path

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"
os.environ["NUMBA_DISABLE_JIT"] = "0"

RATES = (("current_179", 0.179), ("proposal_087510", 0.08751017424959717))
PSI = 0.14245465246024056
CHECKPOINT_SHA = "83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee"
PLAN_SHA = "9cad55d0b8ec186b95d3cb4e9f866add69ba78d5e894b0cdf780ce8191f1889c"
SOURCE_SHA = "76406fcc10206d9e30bcc29d4e18accdf7fffdd9219e50e6e11c1c3360f01336"
TARGET_SHA = "10d80bcfca64aad511058b4a1d427283c3c2f548803213f21688e091f7457107"
PARAMETERS = {
    "H0": 8.967215659504557, "beta_annual": 0.9762497404726131,
    "chi": 0.8879849680171324, "first_birth_fixed_cost": 0.35326910501721936,
    "h_P": 2.044201723718249, "kappa_fert": 0.2092758576051789,
    "kappa_fert_continuation": 0.48206876668272913,
    "theta0": 0.08816920777433607, "theta1": 0.0875688104118386,
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def finite_json(value):
    """Preserve unavailable legacy diagnostics as null, without changing arrays."""
    if isinstance(value, dict):
        return {k: finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def source_preflight(source_root, selected_plan):
    require(sha(selected_plan) == PLAN_SHA, "selected plan hash differs")
    plan = json.loads(Path(selected_plan).read_text())
    require(plan["case_id"] == "worker09_proposal16", "selected case differs")
    require(plan["target_system_sha256"] == TARGET_SHA, "target fingerprint differs")
    require(plan["structural_parameters"] == PARAMETERS, "selected parameters differ")
    require(plan["preference_specification"]["mapping"] == "floor_control", "utility differs")
    require(plan["income_specification"]["constructor_arguments"]["n_persistent"] == 15,
            "persistent-income grid differs")
    require(plan["entry_specification"]["rule"] == "fixed_reference_marginal", "entry rule differs")
    require(plan["wealth_grid_specification"]["upper"] == 3000, "wealth grid differs")
    files = plan["source_manifest"]["files"]
    require(len(files) == plan["source_manifest"]["file_count"] == 641, "source inventory differs")
    import hashlib as _hashlib
    canonical = json.dumps(files, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    require(_hashlib.sha256(canonical).hexdigest() == SOURCE_SHA, "source manifest fingerprint differs")
    root = Path(source_root).resolve()
    for relative, expected in files.items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root) and path.is_file() and sha(path) == expected,
                "frozen source differs: " + relative)
    return plan


def actual_parameters(P):
    return {
        "H0": float(P.H0[0]), "beta_annual": float(P.beta) ** .25,
        "chi": float(P.chi), "first_birth_fixed_cost": float(P.first_birth_fixed_cost),
        "h_P": float(P.hbar_first_child_jump + P.hbar_child_rooms),
        "kappa_fert": float(P.kappa_fert),
        "kappa_fert_continuation": float(P.kappa_fert_continuation),
        "theta0": float(P.theta0), "theta1": float(P.theta1),
    }


def check_checkpoint(packet):
    P = packet["parameters"]
    require(float(P.psi_child) == PSI, "selected psi differs")
    require(float(P.tau_pay) == RATES[0][1], "selected payroll rate differs")
    require(float(P.hbar_child_rooms) == 0., "selected floor utility decomposition differs")
    for key, value in PARAMETERS.items():
        require(abs(actual_parameters(P)[key] - value) <= 2e-12 * max(1, abs(value)),
                "selected checkpoint structural parameter differs: " + key)
    require(int(P.J) == 17 and int(P.I) == 1 and int(P.Nb) == 160, "selected grid dimensions differ")
    require(float(packet["b_grid"][-1]) == 3000., "selected wealth upper bound differs")
    require(float(P.xi_supply[0]) == .63, "selected housing supply elasticity differs")


def target_rows(objective, early, recent_value, completed_fertility):
    projection = objective["cps_projection"]
    mapping = {
        "cps_childlessness": ("fertility", projection, "moments", "childless_rate_40_44"),
        "cps_exactly_one": ("fertility", projection, "moments", "exactly_one_among_mothers_40_44"),
        "nchs_mean_age": ("fertility", projection, "moments", "period_mean_age_first_birth"),
        "nchs_share30": ("fertility", projection, "moments", "period_share_first_births_age30plus"),
        "wealth_earnings": ("housing_wealth", "moments", "aggregate_wealth_to_annual_gross_labor_earnings"),
        "bequest_wealth": ("housing_wealth", "moments", "annual_bequest_flow_to_aggregate_wealth"),
        "old_dispersion": ("housing_wealth", "moments", "old_total_wealth_to_annual_income_p90_p50_7684"),
        "mean_rooms": ("housing_wealth", "moments", "aggregate_mean_occupied_rooms_capped9_18_85"),
        "ownership_30_55": ("housing_wealth", "moments", "own_rate_30_55"),
        "first_birth_rooms": ("housing_wealth", "moments", "housing_increment_0to1"),
        "family_rooms": ("housing_wealth", "moments", "prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9"),
    }
    rows = []
    for original in objective["target_rows"]:
        name = original["restriction_id"]
        if name == "initial_normalization":
            value = completed_fertility
        elif name == "recent_parent_ownership":
            value = recent_value
        else:
            value = early
            for part in mapping[name]:
                value = value[part]
        value = float(value)
        target = float(original["target"])
        weight = original["actual_weight"]
        rows.append(dict(moment=name, target=target, model=value, gap=value-target,
                         weight="" if weight is None else weight,
                         loss_contribution="" if weight is None else float(weight)*(value-target)**2,
                         descriptive_fixed_psi=True))
    require(len(rows) == 13, "full frozen target table unavailable")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--selected-plan", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--objective", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget-seconds", type=int, default=1800)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    require(args.budget_seconds == 1800, "shared computation budget must be 1800 seconds")
    plan = source_preflight(args.source_root, args.selected_plan)
    require(sha(args.objective) == "be21e426f6b2f67b220ab42a65d928b77d8a86354702c6510ad136f02f09ce1b",
            "frozen target/objective file differs")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    cases = [dict(case=case, payroll_tax=rate, status="pending") for case, rate in RATES]
    write(out / "latest_completed.json", dict(status="none", cases=cases))
    if args.dry_run:
        if args.checkpoint is not None:
            require(str(args.source_root.resolve()) == plan["source_root"] and
                    sha(args.checkpoint) == CHECKPOINT_SHA,
                    "remote dry run requires actual selected checkpoint and frozen source")
        for item in cases:
            item["status"] = "dry_run_case_slot_verified_no_solve"
            write(out / "latest_completed.json", dict(status="dry_run", cases=cases))
        write(out / "dry_run.json", dict(status="passed", model_solves=0, cases=cases,
                                          selected_plan_sha256=PLAN_SHA, source_manifest_sha256=SOURCE_SHA))
        return
    require(args.checkpoint is not None and sha(args.checkpoint) == CHECKPOINT_SHA,
            "actual selected checkpoint missing or hash differs")
    sys.dont_write_bytecode = True
    os.environ["NUMBA_CACHE_DIR"] = str(out / "numba_cache")
    (out / "numba_cache").mkdir()
    require(str(args.source_root.resolve()) == plan["source_root"], "remote frozen source path differs")
    bundle_tools = args.source_root.parent / "tools"
    require(sha(bundle_tools / "e5f_earnings_wealth_contract.py") ==
            plan["files"]["accounting"]["sha256"], "reviewed runtime accounting helper differs")
    require(sha(bundle_tools / "run_e5f_preference_share_candidate.py") ==
            plan["files"]["adapter"]["sha256"], "reviewed floor adapter differs")
    sys.path[:0] = [str(args.source_root / "code/model/tools"), str(args.source_root / "code/model")]
    with gzip.open(args.checkpoint, "rb") as stream:
        selected = pickle.load(stream)
    check_checkpoint(selected)
    import run_e5f_matched_pf_smoke as primitive
    sys.path.insert(0, str(bundle_tools))
    import e5f_earnings_wealth_contract as accounting
    import run_e5f_independent_numerical_audit as audit
    from e5f_stationary_paygo import solve_balanced_initial_equilibrium, certify_initial_pension
    from e5f_social_security import fiscal_accounts
    from e5f_initial_fertility_observer import observe_initial_fertility
    from e5f_initial_housing_observer import observe_initial_housing_wealth
    from e5f_recent_parent_flow_observer import observe_recent_parent_flow, SNAPSHOT, AGE_PROJECTION
    chain, model = primitive.pf.transition.configure_sequential_model()
    accounting.install_fixed_entry(model)
    accounting.install_explicit_transaction_grid(model)
    runtime_diff = out / "reviewed_purchase_runtime.diff"
    accounting.install_purchase_income(model, runtime_diff)
    expected_runtime = {
        ".diff": "592a3c90bbc187ba6aecd5fe0b6491ef2f90aa965f4043eae0101baceac23e4f",
        ".generated.py": "900b8bac964a2f91e8f9f47c566cd89c25e316e5db0afbe0975b58cf1a52bd28",
        ".tenure.py": "01ff6c5402cf8317d4b927916f0ac1233026fece2716138163a44990b13403a0",
        ".allocation.py": "86c56d2666f69f757dd69a5b87647686bbd5c4205a8a07cf20651317065e4d8a",
    }
    for suffix, expected in expected_runtime.items():
        require(sha(runtime_diff.with_suffix(suffix)) == expected,
                "reviewed B-floor generated runtime differs: " + suffix)
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    primitive.pf.transition.calendar.model = model
    objective = json.loads(args.objective.read_text())
    require(objective["source_fingerprints"] and plan["objective_canonical_sha256"] ==
            "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1",
            "frozen objective differs")
    start = time.monotonic()
    done = threading.Event()
    active = {"case": None, "status": "started", "elapsed_seconds": 0.}
    def heartbeat():
        while not done.wait(60):
            active["elapsed_seconds"] = time.monotonic()-start
            write(out / "heartbeat.json", active)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        for item in cases:
            require(time.monotonic()-start < args.budget_seconds, "shared 1800-second budget exhausted")
            case, rate = item["case"], item["payroll_tax"]
            case_out = out / case
            case_out.mkdir()
            active.update(case=case, status="solving", elapsed_seconds=time.monotonic()-start)
            write(out / "heartbeat.json", active)
            try:
                P = copy.deepcopy(selected["parameters"])
                require(float(P.psi_child) == PSI, "psi changed before solve")
                sol, P, price, fiscal = solve_balanced_initial_equilibrium(
                    model=model, parameters=P, b_grid=selected["b_grid"],
                    initial_prices=selected["solution"].p_eq, payroll_tax=rate,
                    marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
                require(float(P.psi_child) == PSI and float(P.tau_pay) == rate, "psi/tax changed in solve")
                write(case_out / "solve_complete.json", finite_json(primitive.pf.calendar.jsonable(dict(
                    status="equilibrium_solve_completed", case=case, rate=rate,
                    psi_child=float(P.psi_child), price=float(price[0]), fiscal=fiscal,
                    best_eq_error=float(sol.timings["best_eq_error"]),
                    strict_converged=bool(sol.timings["strict_converged"]),
                    elapsed_seconds=time.monotonic()-start))))
                grid = selected["b_grid"]
                shared = model.precompute_shared(P, grid)
                P._fert2_probs = sol.fert2_probs.copy()
                policy = primitive.pf.calendar.policy_from_solution(sol, price, P, grid, shared)
                pre, reconstruction = primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol, policy, P, grid, shared)
                operator = primitive.pf.transition.operator_gates(sol, policy, pre, P, grid, shared)
                operator.update(reconstruction)
                for name in ("stationary_post_fertility_nesting_l1", "one_step_constant_path_nesting_l1",
                             "mature_flow_abs_error", "birth_flow_abs_error", "topcode_adjusted_birth_flow_abs_error"):
                    require(abs(operator[name]) <= 5e-9, "operator gate failed: " + name)
                require(abs(operator["zero_entry_mass_accounting_residual"]) <= 2e-8 and
                        operator["stationary_feasibility_projection_mass"] <= 1e-6, "mass/feasibility gate failed")
                supply = primitive.pf.calendar.HousingSupplyRule("static-elastic", float(price[0]),
                    float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),
                    float(P.xi_supply[0]))
                evaluation = primitive.pf.calendar.evaluate_period(price, pre, P, grid, shared,
                    primitive.pf.calendar.SolveCounter(), supply_rule=supply, supplied_policy=policy)
                require(evaluation.relative_market_residual <= 2e-4, "dated market gate failed")
                budget = primitive.dated_budget(evaluation, P, shared, grid,
                                                float(P.user_cost_rate*price[0]))
                purchase = accounting.audit_purchase_accounting(evaluation, P, shared, grid, model)
                fiscal = certify_initial_pension(evaluation.g_current, P,
                    marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
                packet = dict(parameters=P, b_grid=grid, evaluation=evaluation, shared=shared,
                    supply_rule=supply, solution=sol, stationary_g_pre=pre,
                    demographic_seed=selected.get("demographic_seed"),
                    contract_sha256=selected.get("contract_sha256"))
                arrays = audit.policy_array_audit(packet, case_out)
                require(arrays["occupied_negative_steps"] == 0, "occupied value monotonicity gate failed")
                require(all(not x["nonfinite"] and x["minimum"] >= 0 and x["maximum"] <= 1
                            for x in arrays["probabilities"].values()), "probability gate failed")
                fertility = {projection: observe_initial_fertility(evaluation, P, age_projection=projection)
                             for projection in ("uniform_birth_time", "constant_post_cell")}
                housing = observe_initial_housing_wealth(evaluation, P, grid, shared,
                    diagnostic_enabled=True, age_projection="uniform_within_age_cell",
                    diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True)
                early = dict(fertility=fertility, housing_wealth=housing)
                case_checkpoint = case_out / "initial_state.pkl.gz"
                with gzip.open(case_checkpoint, "wb", compresslevel=1) as stream:
                    pickle.dump(packet, stream, protocol=5)
                case_checkpoint_sha = sha(case_checkpoint)
                recent = observe_recent_parent_flow(evaluation, P, diagnostic_enabled=True,
                    snapshot=SNAPSHOT, age_projection=AGE_PROJECTION,
                    diagnostic_allow_residence_proxy=True,
                    input_provenance=dict(case_id=case, checkpoint_sha256=case_checkpoint_sha))
                moments = chain.extract_moments(sol, P)
                rows = target_rows(objective, early, recent["model_value"], float(moments["tfr"]))
                accounts = fiscal_accounts(evaluation.g_current, P)
                gross_annual = accounts["payroll_tax_base_period"] / 4
                worker_mass = accounts["worker_household_mass"]
                main = dict(case=case, rate=rate, psi_child=float(P.psi_child), price=float(price[0]),
                    fiscal=accounts, gross_annual_per_worker=gross_annual/worker_mass,
                    disposable_annual_per_worker=(1-rate)*gross_annual/worker_mass,
                    aggregate_wealth=housing["moments"]["aggregate_wealth_to_annual_gross_labor_earnings"]*gross_annual,
                    legacy_moments=moments, early_moments=early, recent_parent=recent,
                    operator_gates=operator, household_budget=budget,
                    purchase_accounting=purchase, policy_array_gates=arrays,
                    target_system_sha256=TARGET_SHA, source_manifest_sha256=SOURCE_SHA,
                    selected_checkpoint_sha256=CHECKPOINT_SHA,
                    case_checkpoint_sha256=case_checkpoint_sha, calibrated_smm=False)
                write(case_out / "summary.json", finite_json(primitive.pf.calendar.jsonable(main)))
                table(case_out / "all_target_moments.csv", rows)
                parameter_rows = [dict(parameter=name, estimate=actual_parameters(P)[name],
                    lower=plan["parameter_bounds"][name][0], upper=plan["parameter_bounds"][name][1],
                    near_bound=min(actual_parameters(P)[name]-plan["parameter_bounds"][name][0],
                        plan["parameter_bounds"][name][1]-actual_parameters(P)[name]) <=
                        .01*(plan["parameter_bounds"][name][1]-plan["parameter_bounds"][name][0]),
                    status="selected experimental structural coordinate, held fixed") for name in PARAMETERS]
                parameter_rows += [dict(parameter=name, estimate=value, lower="", upper="",
                    near_bound="", status=status) for name, value, status in (
                        ("psi_child", P.psi_child, "selected normalized value, held fixed"),
                        ("payroll_tax", P.tau_pay, "case rate"),
                        ("pension_period", P.pension, "endogenous PAYGO balance"))]
                table(case_out / "parameters.csv", parameter_rows)
                audit.standard_diagnostics(packet, case_out, validate_production_young=False)
                require(len(list((case_out / "standard_diagnostics").glob("*.png"))) == 17,
                        "standard diagnostic packet incomplete")
                item.update(status="completed", price=float(price[0]),
                            pension_period=float(P.pension), elapsed_seconds=time.monotonic()-start)
                write(out / "latest_completed.json", dict(status="running", cases=cases))
            except Exception as exc:
                item.update(status="failed", error_type=type(exc).__name__, error=str(exc),
                            elapsed_seconds=time.monotonic()-start)
                write(case_out / "failure.json", item)
                write(out / "latest_completed.json", dict(status="failed", cases=cases))
                raise
        write(out / "complete.json", dict(status="completed", cases=cases,
            elapsed_seconds=time.monotonic()-start, solves=2, calibrated_smm=False))
    finally:
        done.set()


if __name__ == "__main__":
    main()
