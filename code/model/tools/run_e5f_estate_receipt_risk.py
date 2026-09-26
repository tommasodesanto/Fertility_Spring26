#!/usr/bin/env python3
"""Bounded Torch-only, fixed-price conditional-mean versus receipt-risk test.

The reference is the retained net-valuation estate case, not the newly adopted
pension contract. Three household solves; no recalibration or price search.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import math
import pickle
import shutil
from pathlib import Path
import subprocess
import sys
import threading
import time

import run_e5f_estate_receiver_probe as probe

TASK = Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/estate_receiver_probe_20260925_v1")
REFERENCE = TASK / "results/v2/run/net_valuation"
PROFILE = TASK / "results/recipient_evidence/v1/inheritance_receipt_profile_mapped_4y.csv"
BUDGET_SECONDS = 25 * 60


def run_cases(cases, solve, report, root, deadline):
    """One ordered solve per case; the control must pass before either change."""
    completed = []
    for name in cases:
        if time.time() >= deadline:
            raise TimeoutError("Receipt-risk stage exhausted its 25-minute budget")
        probe.write(root / "progress.json", dict(status="solving", case=name,
                                                 completed=[r["case"] for r in completed]))
        tick = time.monotonic()
        sol, P = solve(name)
        receipt = report(name, sol, P)
        receipt["wall_seconds"] = time.monotonic() - tick
        completed.append(receipt)
        probe.write(root / name / "receipt.json", receipt)
        probe.write(root / "latest_completed.json", receipt)
        probe.write(root / "best_so_far.json", dict(
            description="descriptive inherited fit only; not specification selection",
            **min(completed, key=lambda row: row["loss"])))
    return completed


def loop_fixture():
    """Exercise the exact ordered loop and stop-on-failure without a solve."""
    import tempfile
    from types import SimpleNamespace
    seen = []
    def solve(name):
        seen.append(name)
        return SimpleNamespace(), SimpleNamespace()
    def report(name, sol, P):
        if name == "conditional_mean":
            raise RuntimeError("intentional fixture gate")
        return dict(case=name, loss=1.)
    with tempfile.TemporaryDirectory() as folder:
        try:
            run_cases(("no_receipt", "conditional_mean", "receipt_lottery"),
                      solve, report, Path(folder), time.time() + 10.)
        except RuntimeError as exc:
            if str(exc) != "intentional fixture gate":
                raise
        else:
            raise AssertionError("The exact case loop did not stop at a failed gate")
    if seen != ["no_receipt", "conditional_mean"]:
        raise AssertionError("Receipt-risk case ordering differs")
    return dict(status="passed", ordered_loop=True, fail_closed=True)


def check_control(sol, pre, reference):
    import numpy as np
    for key in ("V", "g", "bp_pol", "c_pol", "hR_pol", "fert_probs", "fert2_probs"):
        np.testing.assert_array_equal(getattr(sol, key), getattr(reference["solution"], key),
                                      err_msg="zero-receipt native reproduction: " + key)
    np.testing.assert_array_equal(sol.p_eq, reference["solution"].p_eq)
    np.testing.assert_array_equal(pre, reference["stationary_g_pre"])


def lifecycle_rows(P, pre, prices, grid, ledger):
    import numpy as np
    rows = []
    for j in range(int(P.J)):
        mass = pre[:, :, :, j]
        by_wealth = mass.sum(axis=tuple(range(1, mass.ndim)))
        total = float(by_wealth.sum())
        cdf = np.cumsum(by_wealth) / total
        quantiles = {f"financial_wealth_p{int(p*100)}": float(grid[min(np.searchsorted(cdf, p), len(grid)-1)])
                     for p in (.1, .5, .9)}
        mean = float(np.dot(grid, by_wealth) / total)
        row = dict(age=float(P.age_start) + j * float(P.da), mass=total,
                   ownership_before_current_choices=float(mass[:, 1:].sum() / total),
                   financial_wealth_before_current_choices=mean,
                   financial_wealth_sd=float(np.sqrt(np.dot((grid - mean)**2, by_wealth) / total)),
                   receipt_probability=float(P.estate_receipt_risk_profile[j]["probability"])
                   if P.estate_receipt_risk_case == "receipt_lottery" else
                   float(P.estate_receipt_risk_case == "conditional_mean" and
                         P.estate_receipt_risk_profile[j]["relative_mean"] > 0),
                   mean_receipt=ledger["by_age"][j]["expected_receipt_flow"] / total,
                   **quantiles)
        rows.append(row)
    return rows


def report_case(name, sol, P, *, output, runtime, tax, objective, reference,
                retained_receipt, scale, pins, adapter, risk):
    import numpy as np
    case_out = output / name
    case_out.mkdir(parents=True, exist_ok=True)
    primitive, model, chain = runtime["primitive"], runtime["model"], runtime["chain"]
    grid = reference["b_grid"]
    prices = np.asarray(reference["solution"].p_eq).copy()
    shared = model.precompute_shared(P, grid)
    P._fert2_probs = sol.fert2_probs.copy()
    ledger = risk.normalized_ledger(P, sol)
    policy = primitive.pf.calendar.policy_from_solution(sol, prices, P, grid, shared)
    pre, reconstruction = primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol, policy, P, grid, shared)
    operator = primitive.pf.transition.operator_gates(sol, policy, pre, P, grid, shared)
    operator.update(reconstruction)
    for key in ("stationary_post_fertility_nesting_l1", "one_step_constant_path_nesting_l1",
                "mature_flow_abs_error", "birth_flow_abs_error", "topcode_adjusted_birth_flow_abs_error"):
        if not math.isfinite(float(operator[key])) or abs(operator[key]) > 5e-9:
            raise RuntimeError("Native operator gate failed: " + key)
    if (abs(operator["zero_entry_mass_accounting_residual"]) > 2e-8
            or operator["stationary_feasibility_projection_mass"] > 1e-6):
        raise RuntimeError("Native mass/projection gate failed")
    if name == "no_receipt":
        check_control(sol, pre, reference)
    supply = primitive.pf.calendar.HousingSupplyRule("static-elastic", float(prices[0]),
        float(P.H0[0] * (P.user_cost_rate * prices[0] / P.r_bar[0])**P.xi_supply[0]), float(P.xi_supply[0]))
    evaluation = primitive.pf.calendar.evaluate_period(prices, pre, P, grid, shared,
        primitive.pf.calendar.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    # Fixed-price attribution deliberately reports excess demand; it is not GE.
    budget = primitive.dated_budget(evaluation, P, shared, grid, float(P.user_cost_rate * prices[0]))
    purchase = runtime["accounting"].audit_purchase_accounting(evaluation, P, shared, grid, model)
    fiscal = runtime["certify_initial_pension"](evaluation.g_current, P,
        marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    packet = dict(parameters=P, b_grid=grid, evaluation=evaluation, shared=shared,
                  supply_rule=supply, solution=sol, stationary_g_pre=pre,
                  contract_sha256=probe.OBJECTIVE_SHA,
                  demographic_seed=reference.get("demographic_seed"))
    with gzip.open(case_out / "initial_state.pkl.gz", "wb", compresslevel=1) as stream:
        pickle.dump(packet, stream, protocol=5)
    arrays = runtime["audit"].policy_array_audit(packet, case_out)
    if arrays["occupied_negative_steps"] or any(row["nonfinite"] or row["minimum"] < 0 or row["maximum"] > 1
                                               for row in arrays["probabilities"].values()):
        raise RuntimeError("Native value/probability gate failed")
    early = dict(fertility={p: runtime["observe_initial_fertility"](evaluation, P, age_projection=p)
                           for p in ("uniform_birth_time", "constant_post_cell")},
                 housing_wealth=runtime["observe_initial_housing_wealth"](evaluation, P, grid, shared,
                     diagnostic_enabled=True, age_projection="uniform_within_age_cell",
                     diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True))
    recent = runtime["observe_recent_parent_flow"](evaluation, P, diagnostic_enabled=True,
        snapshot=runtime["SNAPSHOT"], age_projection=runtime["AGE_PROJECTION"],
        diagnostic_allow_residence_proxy=True, input_provenance={"case_id": name})
    rows = tax.target_rows(objective, early, recent["model_value"], float(chain.extract_moments(sol, P)["tfr"]))
    if len(rows) != 13 or sum(row["weight"] != "" for row in rows) != 12:
        raise RuntimeError("Full frozen target table failed")
    if any(not math.isfinite(float(row[key])) for row in rows for key in ("model", "target", "gap")):
        raise RuntimeError("Nonfinite target-fit row")
    loss = sum(float(row["loss_contribution"]) for row in rows if row["loss_contribution"] != "")
    if name == "no_receipt" and reference.get("wealth_grid_subdivision", 1) == 1:
        retained_rows = {r["moment"]: r for r in csv.DictReader((REFERENCE / "target_fit.csv").open())}
        if loss != retained_receipt["loss"] or any(float(r["model"]) != float(retained_rows[r["moment"]]["model"]) for r in rows):
            raise RuntimeError("Zero-receipt control does not reproduce every retained target")
    probe.table(case_out / "target_fit.csv", rows)
    parameters = list(csv.DictReader((REFERENCE / "parameters.csv").open()))
    actual = tax.actual_parameters(P)
    for row in parameters:
        if row["parameter"] in actual and float(row["estimate"]) != float(actual[row["parameter"]]):
            raise RuntimeError("Fixed structural parameter changed: " + row["parameter"])
    parameters.append(dict(parameter="receipt_profile_scale", estimate=scale, lower="", upper="", near_bound="",
                           status="fixed to retained net estate pool; empirical relative age profile; experimental"))
    probe.table(case_out / "parameters.csv", parameters)
    estate = adapter.estate_accounts(sol, P, prices)
    account = dict(generated_net_period=estate["generated_net_period"],
                   generated_gross_period=estate["generated_gross_period"],
                   paid_period=ledger["paid_period"],
                   paid_minus_generated=ledger["paid_period"] - estate["generated_net_period"],
                   interpretation="fixed receipt scale; funding residual reported, not closed")
    probe.write(case_out / "receipt_transport.json", ledger)
    probe.write(case_out / "operator_gates.json", operator)
    probe.table(case_out / "lifecycle_and_receipts.csv", lifecycle_rows(P, pre, prices, grid, ledger))
    runtime["audit"].standard_diagnostics(packet, case_out, validate_production_young=False)
    if len(list((case_out / "standard_diagnostics").glob("*.png"))) != 17:
        raise RuntimeError("Standard 17-figure diagnostic packet incomplete")
    return dict(status="completed_experimental_fixed_price_case", case=name, loss=loss,
                payroll_tax=float(P.tau_pay), price=prices.tolist(), receipt_scale=scale,
                source_pins=pins, target_system_sha256=probe.OBJECTIVE_SHA,
                relative_market_residual=float(evaluation.relative_market_residual),
                market_clearing_required=False, fiscal=fiscal, household_budget=budget,
                purchase_accounting=purchase, policy_array_gates=arrays, estate_accounts=account,
                operator_gates=operator, clipped_wealth_loss=ledger["clipped_wealth_loss"],
                fixed_preference_birth_replacement_gap=float(sol.adult_entry_stationary_relative_gap),
                zero_receipt_exact_reproduction=name == "no_receipt",
                control_reference="retained net-valuation arrays and targets" if reference.get("wealth_grid_subdivision", 1) == 1
                else "fresh unpatched net-valuation arrays and calendar distribution on the finer grid",
                wealth_grid_subdivision=reference.get("wealth_grid_subdivision", 1),
                experimental_changes="age-pooled receipt profile, restricted age support, IID zero/positive risk versus certain conditional mean, start-period wealth timing; fixed net-valuation price and fixed scale; no recalibration; inherited experimental tax")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--subdivide-wealth-grid", type=int, choices=(1, 2), default=1,
                        help="2 inserts each interval midpoint; four solves including an unpatched numerical control")
    args = parser.parse_args()
    probe.require_torch()
    if args.output_dir.exists():
        raise FileExistsError("Refusing to replace an earlier receipt-risk run")
    output = args.output_dir
    output.mkdir(parents=True)
    started = time.monotonic()
    deadline = time.time() + BUDGET_SECONDS
    done = threading.Event()
    def pulse():
        while not done.wait(55):
            progress = json.loads((output / "progress.json").read_text()) if (output / "progress.json").exists() else {}
            probe.write(output / "heartbeat.json", dict(elapsed_seconds=time.monotonic() - started,
                budget_seconds=BUDGET_SECONDS, updated_epoch=time.time(), **progress))
    thread = threading.Thread(target=pulse, daemon=True)
    thread.start()
    try:
        old, tax, plan, selected, objective, runtime, selected_receipt = probe.verified_runtime(output)
        import numpy as np
        import e5f_estate_receipt_risk_adapter as risk
        adapter = probe.load_adapter()
        adapter.install(runtime["model"], output / "net_adapter")
        hook_names = ("solve_bellman_full_markov_income", "forward_distribution_markov_income",
                      "advance_cohort_one_period_markov_income")
        unpatched_functions = {name: getattr(runtime["model"], name) for name in hook_names}
        installed = risk.install(runtime["model"], output / "risk_adapter")
        with gzip.open(REFERENCE / "initial_state.pkl.gz", "rb") as stream:
            reference = pickle.load(stream)
        retained = json.loads((REFERENCE / "receipt.json").read_text())
        if (retained["source_manifest_sha256"] != probe.SOURCE_SHA
                or retained["target_system_sha256"] != probe.OBJECTIVE_SHA
                or retained["payroll_tax"] != probe.TAX or retained["case"] != "net_valuation"):
            raise RuntimeError("Retained net-valuation reference differs from contract")
        P0, grid = reference["parameters"], reference["b_grid"]
        np.testing.assert_array_equal(grid, selected["b_grid"])
        runtime["selected"] = selected
        ages = float(P0.age_start) + np.arange(int(P0.J)) * float(P0.da)
        profile = risk.pooled_profile(PROFILE, ages)
        base_mass = np.asarray(reference["solution"].g).sum(axis=(0, 1, 2, 4, 5, 6))
        mean_profile = np.asarray([row["relative_mean"] for row in profile])
        D0 = adapter.estate_accounts(reference["solution"], P0, reference["solution"].p_eq)["generated_net_period"]
        scale = D0 / float(np.dot(base_mass, mean_profile))
        if not np.isfinite(scale) or scale <= 0:
            raise RuntimeError("Invalid baseline estate/profile scale")
        if args.subdivide_wealth_grid == 2:
            refined = np.sort(np.r_[grid, .5 * (grid[:-1] + grid[1:])])
            P0 = risk.embed_reference_entry_on_refined_grid(P0, grid, refined)
            grid = refined
        solve_count = 3 + int(args.subdivide_wealth_grid == 2)
        files = [Path(__file__), Path(risk.__file__), Path(__file__).with_name("e5f_estate_receipt_jump.py"),
                 Path(__file__).with_name("test_e5f_estate_receipt_jump.py"),
                 Path(__file__).with_name("test_e5f_estate_receipt_risk_adapter.py")]
        pins = dict(frozen_source_manifest=probe.SOURCE_SHA, objective=probe.OBJECTIVE_SHA,
                    selected_original_checkpoint=probe.sha(probe.CHECKPOINT),
                    retained_net_checkpoint=probe.sha(REFERENCE / "initial_state.pkl.gz"),
                    receipt_profile=probe.sha(PROFILE), diagnostic_files={p.name: probe.sha(p) for p in files})
        (output / "source").mkdir()
        for path in files:
            shutil.copyfile(path, output / "source" / path.name)
        probe.write(output / "contract.json", dict(source_pins=pins, installed_hooks=installed,
            cases=list(risk.CASES), scale=scale, profile=profile, fixed_price=reference["solution"].p_eq.tolist(),
            fixed_reference_estate_pool=D0, maximum_household_solves=solve_count,
            wealth_grid_subdivision=args.subdivide_wealth_grid, wealth_grid_nodes=int(len(grid)),
            wealth_grid_bounds=[float(grid[0]), float(grid[-1])],
            entry_grid_rule="original point masses and income conditionals retained exactly; zero mass at inserted knots",
            budget_seconds=BUDGET_SECONDS, estimated_minutes=10,
            support="exact age nodes26..78; experimental zero receipts at18,22,82",
            rank_mapping="age-only pooled usual-income groups50/40/10; no labor-income rank mapping",
            scientific_scope="fixed-price/fixed-scale risk attribution; unbalanced estates and housing reported; no GE or recalibration",
            adopted=False))
        tests = [str(p) for p in files if p.name.startswith("test_")]
        with (output / "tests.log").open("w") as stream:
            subprocess.run([sys.executable, "-m", "pytest", "-q", *tests], stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=180)
        probe.write(output / "smoke.json", loop_fixture())
        if args.subdivide_wealth_grid == 2:
            # Independent numerical control: solve with the original net
            # Bellman and original native/cohort distribution functions.
            # Source code and the original saved reference are untouched.
            probe.write(output / "progress.json", dict(status="solving", case="unpatched_finer_grid_reference", completed=[]))
            model = runtime["model"]
            patched_functions = {name: getattr(model, name) for name in hook_names}
            try:
                for name, function in unpatched_functions.items():
                    setattr(model, name, function)
                reference_P = copy.deepcopy(P0)
                reference_sol = model.solve_markov_income_at_prices(
                    np.asarray(reference["solution"].p_eq).copy(), reference_P, grid, fast_stats=False)
                shared = model.precompute_shared(reference_P, grid)
                reference_P._fert2_probs = reference_sol.fert2_probs.copy()
                calendar = runtime["primitive"].pf.calendar
                policy = calendar.policy_from_solution(reference_sol, reference_sol.p_eq, reference_P, grid, shared)
                reference_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
                    reference_sol, policy, reference_P, grid, shared)
            finally:
                for name, function in patched_functions.items():
                    setattr(model, name, function)
            reference = dict(reference, b_grid=grid, parameters=reference_P, solution=reference_sol,
                             stationary_g_pre=reference_pre, wealth_grid_subdivision=2)
            probe.write(output / "finer_grid_reference.json", dict(status="unpatched_reference_solved",
                reconstruction=reconstruction, wealth_grid_nodes=int(len(grid)),
                original_price_and_receipt_scale_fixed=True, economic_changes=[]))
        def solve(name):
            P = risk.configure(P0, name, profile, scale, grid)
            sol = runtime["model"].solve_markov_income_at_prices(
                np.asarray(reference["solution"].p_eq).copy(), P, grid, fast_stats=False)
            return sol, P
        def report(name, sol, P):
            return report_case(name, sol, P, output=output, runtime=runtime, tax=tax,
                objective=objective, reference=reference, retained_receipt=retained,
                scale=scale, pins=pins, adapter=adapter, risk=risk)
        completed = run_cases(risk.CASES, solve, report, output, deadline)
        paid = [row["estate_accounts"]["paid_period"] for row in completed[1:]]
        if max(abs(value - D0) for value in paid) > 1e-10:
            raise RuntimeError("Paired age-conditional means no longer preserve common total expected receipts")
        probe.write(output / "complete.json", dict(status="completed", cases=completed,
            elapsed_seconds=time.monotonic() - started, household_solves=solve_count, source_pins=pins,
            decision_scope="receipt-risk attribution only; economic specification remains experimental"))
    except Exception as exc:
        failure = dict(status="failed", error_type=type(exc).__name__, error=str(exc),
                       elapsed_seconds=time.monotonic() - started)
        if hasattr(exc, "account"):
            failure["account"] = exc.account
        probe.write(output / "failed.json", failure)
        raise
    finally:
        done.set()
        thread.join(timeout=2)


if __name__ == "__main__":
    main()
