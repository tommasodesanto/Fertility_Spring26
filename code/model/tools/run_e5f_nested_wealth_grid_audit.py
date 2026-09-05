#!/usr/bin/env python3
"""Bounded nested-grid impact check, conditional on the certified 2023 measure.

Run inside the frozen production snapshot. All economic parameters, inherited
state probabilities, endpoints, and housing choices are held fixed. Only wealth
interpolation resolution changes; global saving is used at every resolution.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import pickle
import threading
import time
from pathlib import Path

import numpy as np

import run_e5f_global_saving_quantification as saving

audit = saving.audit
model, baseline, calendar, transition = saving.model, saving.baseline, saving.calendar, saving.transition
AUDIT_SHA = "ba0f94f52705f43fd0bf7a18ea8c2053f9897674abeb20b030eaf37c57de8623"
SAVING_SHA = "0e89d55a6cd3d0432dcff39b24683a7f6233ac0820a03b47a94eda8e2a924586"
CHECKPOINT_SHA = "bbe10a21a843facaf2bceed56e89281e00992d632cddab640ff0c572d3eb494f"
MARKET_TOL = 1e-5
POLICIES = ("baseline", "supply-plus-20", "dependent-child-ltv95")


def nested_grid(bg, subdivisions):
    bg = np.asarray(bg, dtype=float)
    if subdivisions not in (1, 2, 4) or bg.ndim != 1 or np.any(np.diff(bg) <= 0):
        raise ValueError("Only sorted original grids subdivided 1, 2, or 4 times are allowed")
    fine = np.empty((bg.size - 1) * subdivisions + 1)
    for k in range(subdivisions):
        fine[k:-1:subdivisions] = bg[:-1] + (k / subdivisions) * np.diff(bg)
    fine[-1] = bg[-1]
    old = np.arange(bg.size) * subdivisions
    if not np.array_equal(fine[old], bg) or np.any(np.diff(fine) <= 0):
        raise RuntimeError("Nested-grid exact-node gate failed")
    return fine, old


def embed_measure(g, bg, fine, old):
    result = np.zeros((fine.size,) + g.shape[1:], dtype=g.dtype)
    result[old] = g
    added = np.ones(fine.size, dtype=bool)
    added[old] = False
    if not np.array_equal(result[old], g) or np.any(result[added] != 0):
        raise RuntimeError("Exact inherited-measure embedding failed")
    original_marginals, fine_marginals = g.sum(axis=0), result.sum(axis=0)
    shape = (-1,) + (1,) * (g.ndim - 1)
    wealth_before = float(np.sum(g * bg.reshape(shape)))
    wealth_after = float(np.sum(result * fine.reshape(shape)))
    receipt = dict(original_nodes=int(bg.size), new_nodes=int(fine.size),
        old_nodes_identical=bool(np.array_equal(fine[old], bg)),
        old_node_mass_identical=bool(np.array_equal(result[old], g)),
        mass_on_added_nodes=float(result[added].sum()),
        original_mass=float(g.sum()), new_mass=float(result.sum()),
        max_other_state_marginal_gap=float(np.max(np.abs(original_marginals-fine_marginals))),
        original_liquid_wealth=wealth_before, new_liquid_wealth=wealth_after,
        wealth_gap=wealth_after-wealth_before, lower_endpoint=float(fine[0]), upper_endpoint=float(fine[-1]))
    if abs(receipt["original_mass"]-receipt["new_mass"]) > 2e-12 or abs(receipt["wealth_gap"]) > 2e-12 or receipt["max_other_state_marginal_gap"] > 2e-12:
        raise RuntimeError("Inherited marginal accounting failed")
    return result, receipt


def economic_fingerprint(P):
    values = {k:v for k,v in vars(P).items() if k != "Nb" and not k.startswith("_")}
    return hashlib.sha256(json.dumps(calendar.jsonable(values), sort_keys=True).encode()).hexdigest()


def run_case(packet, subdivisions, name, out, references, *, clear=True):
    saving.operator_contract()
    saving.set_method("global")
    P = copy.deepcopy(packet["parameters"])
    audit.policy.apply_policy(P, audit.policy.POLICIES[name])
    economic_before = economic_fingerprint(P)
    fine, old = nested_grid(packet["b_grid"], subdivisions)
    P.Nb = fine.size
    if economic_fingerprint(P) != economic_before:
        raise RuntimeError("An economic parameter changed with grid resolution")
    rule = audit.policy.policy_supply_rule(packet["supply_rule"], audit.policy.POLICIES[name])
    state = copy.copy(packet["state"])
    state.g_pre, embedding = embed_measure(state.g_pre, packet["b_grid"], fine, old)
    state.initial_policy = None
    state.price_guess = float(references[name]["asset_price"])
    label = f"n{fine.size}_{name}_{'cleared' if clear else 'reproduction'}"
    case_out = out / "cases" / label
    case_out.mkdir(parents=True)
    audit.save_json(case_out / "embedding_receipt.json", embedding)
    started = time.perf_counter()
    counter = calendar.SolveCounter()
    audit.progress(out, "case_started", case=label, subdivisions=subdivisions, nodes=fine.size)
    if clear:
        e, shared, fallback = baseline.evaluate_state(state, P, fine, rule, counter, MARKET_TOL, 30)
    else:
        shared = model.precompute_shared(P, fine)
        e = calendar.evaluate_period(np.asarray([state.price_guess]), state.g_pre, P, fine, shared, counter, rule)
        fallback = False
    health = calendar.distribution_health(dict(pre=e.g_pre, post=e.g_post_fertility, current=e.g_current))
    case_packet = dict(packet, parameters=P, b_grid=fine, evaluation=e, state=state,
        shared=shared, supply_rule=rule, diagnostic_case=dict(method="global", policy=name,
        subdivisions=subdivisions, interpretation="common original pre-choice measure; no history re-estimation"))
    # Save failure evidence too: production floors and any monotonicity flags
    # are measured before checking acceptance. No numerical gate is weakened.
    budget = audit.budget_audit(case_packet, case_out)
    arrays = audit.policy_array_audit(case_packet, case_out)
    audit.standard_diagnostics(case_packet, case_out, validate_production_young=False)
    checkpoint = case_out / "case_state.pkl.gz"
    with gzip.open(checkpoint, "wb", compresslevel=1) as stream:
        pickle.dump(case_packet, stream, protocol=5)
    quantities = saving.quantities(e, P)
    common_mass_gap = abs(e.g_post_fertility.sum() - packet["evaluation"].g_post_fertility.sum())
    pre_gate_change = float(np.max(np.abs(e.g_pre-state.g_pre)))
    current_mass_gap = abs(e.g_current.sum() - e.g_post_fertility.sum())
    probability_ok = all(p["nonfinite"] == 0 and p["minimum"] >= -1e-7 and p["maximum"] <= 1+1e-7 for p in arrays["probabilities"].values())
    gates = dict(finite_nonnegative_mass=health["nonfinite_distribution_count"] == 0 and health["min_distribution_mass"] >= -1e-14,
        common_initial_mass=bool(common_mass_gap <= 2e-10 and pre_gate_change <= 2e-10),
        current_mass_accounting=bool(current_mass_gap <= 2e-10),
        production_feasibility=bool(e.feasibility_projection_mass <= 1e-6),
        production_market=bool(e.relative_market_residual <= 2e-4),
        requested_market=bool(e.relative_market_residual <= MARKET_TOL),
        probability_ranges=probability_ok,
        occupied_wealth_value_monotonicity=arrays["occupied_negative_steps"] == 0,
        standard_graphs=len(list((case_out / "standard_diagnostics").glob("*.png"))) == 17)
    reproduction_gap = None
    if subdivisions == 1:
        reproduction_gap = max(abs(quantities[k]-float(references[name][k])) for k in
            ("births_per_household", "ownership", "rooms_per_household", "asset_price", "explicit_births", "adjusted_births"))
        gates["coarse_exact_reproduction"] = reproduction_gap <= 2e-10
    row = dict(case=label, nodes=int(fine.size), subdivisions=subdivisions, policy=name,
        method="global", market_cleared=clear, **quantities,
        requested_market_tolerance=MARKET_TOL, all_gates_passed=all(gates.values()),
        current_mass_gap=float(current_mass_gap), common_mass_gap=float(common_mass_gap),
        pre_gate_max_change=pre_gate_change, coarse_reproduction_gap=reproduction_gap,
        budget_excess_share=float(budget["budget_excess_share"]),
        weighted_positive_budget_gap=float(budget["weighted_positive_budget_gap"]),
        occupied_negative_value_steps=arrays["occupied_negative_steps"],
        bellman_solves=counter.bellman, elapsed_seconds=time.perf_counter()-started,
        grid_fallback=fallback, case_checkpoint_sha256=audit.digest(checkpoint),
        diagnostic_directory=str(case_out))
    audit.save_json(case_out / "case_receipt.json", dict(row=row, gates=gates, embedding=embedding,
        policy_arrays=arrays, budget=budget, economic_parameters_excluding_Nb_sha256=economic_before,
        inherited_checkpoint_sha256=CHECKPOINT_SHA,
        interpretation="Budget-floor excess is reported separately; its nonzero value is not declared exactly feasible. Grid effects are conditional on the unchanged original inherited history."))
    audit.save_json(out / "latest_completed_case.json", row)
    audit.save_json(out / "best_so_far.json", dict(status="diagnostic_only_no_calibration_selection", latest_admissible_case=row if all(gates.values()) else None))
    audit.progress(out, "case_complete", **row)
    if not all(gates.values()):
        raise RuntimeError(f"Case acceptance failed; preserved evidence: {gates}")
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--reference-csv", type=Path, required=True)
    parser.add_argument("--reference-csv-sha256", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--stage", choices=("smoke", "fine-policies", "finer"), required=True)
    parser.add_argument("--previous-summary", type=Path)
    parser.add_argument("--previous-summary-sha256")
    args = parser.parse_args()
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    for path, expected in ((args.checkpoint, CHECKPOINT_SHA), (Path(audit.__file__), AUDIT_SHA),
                           (Path(saving.__file__), SAVING_SHA), (args.reference_csv, args.reference_csv_sha256)):
        if audit.digest(path) != expected:
            raise RuntimeError(f"Pinned input/helper differs: {path}")
    previous = None
    if args.stage != "smoke":
        if not args.previous_summary or audit.digest(args.previous_summary) != args.previous_summary_sha256:
            raise RuntimeError("Previous-stage receipt missing or mismatched")
        previous = json.loads(args.previous_summary.read_text())
        required_stage = "smoke" if args.stage == "fine-policies" else "fine-policies"
        if previous.get("status") != "complete" or previous.get("stage") != required_stage or previous.get("driver_sha256") != audit.digest(__file__):
            raise RuntimeError("Required exact-loop smoke/stage did not pass this driver")
    transition.configure_sequential_model()
    saving.install_sequential_operators()
    packet = audit.load_checkpoint(args.checkpoint)
    packet["checkpoint_sha256"] = CHECKPOINT_SHA
    _, configured = transition.configure_sequential_model()
    live = baseline.calibration.code_fingerprint_contract(configured)["bundle_sha256"]
    if live != audit.BUNDLE or packet["input_hashes"]["code_bundle_sha256"] != live:
        raise RuntimeError("Frozen production source differs")
    with args.reference_csv.open() as stream:
        references = {r["policy"]:r for r in csv.DictReader(stream) if r["method"] == "global" and r["pricing"] == "separately_cleared"}
    if set(references) != set(POLICIES):
        raise RuntimeError("Three completed global-saving reference cases required")
    started = time.perf_counter()
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            audit.save_json(out / "runtime_heartbeat.json", dict(elapsed_seconds=time.perf_counter()-started,
                utc=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())))
    threading.Thread(target=heartbeat, daemon=True).start()
    original_solve = calendar.solve_policy
    def traced_solve(price, P, bg, shared, counter):
        t = time.perf_counter()
        result = original_solve(price, P, bg, shared, counter)
        audit.progress(out, "bellman_complete", nodes=len(bg), price=np.asarray(price).tolist(), seconds=time.perf_counter()-t, count=counter.bellman)
        return result
    calendar.solve_policy = traced_solve
    contract = dict(stage=args.stage, source_bundle_sha256=live, driver_sha256=audit.digest(__file__),
        audit_helper_sha256=AUDIT_SHA, global_saving_helper_sha256=SAVING_SHA,
        inherited_checkpoint_sha256=CHECKPOINT_SHA, reference_csv_sha256=args.reference_csv_sha256,
        previous_summary_sha256=args.previous_summary_sha256, operator_contract=saving.operator_contract(),
        estimated_parameters="all eleven fixed at task_010", external_restrictions="eta=0.63; tenure kappa=0.005",
        population="exact original 2023 pre-choice household measure, embedded without redistribution",
        entry="inherited queue unchanged; no forward population transition evaluated",
        fiscal="existing supply and family-credit policies; tax and rebate objects unchanged",
        geography="unchanged single-market aggregate closure", resolution="explicit diagnostic override of wealth grid only",
        outstanding="history on refined grid, re-estimation, upper-bound expansion, and long-run policy effects are not evaluated")
    audit.save_json(out / "run_contract.json", contract)
    audit.save_json(out / "latest_completed_case.json", dict(status="no_case_complete"))
    audit.save_json(out / "best_so_far.json", dict(status="diagnostic_only_no_calibration_selection"))
    cases = [(1, "baseline", False), (2, "baseline", True)] if args.stage == "smoke" else (
        [(2, name, True) for name in POLICIES[1:]] if args.stage == "fine-policies" else [(4, name, True) for name in POLICIES])
    rows = []
    try:
        for subdivisions, name, clear in cases:
            rows.append(run_case(packet, subdivisions, name, out, references, clear=clear))
            baseline.write_csv(out / "case_summary.csv", rows)
        audit.save_json(out / "summary.json", dict(contract, status="complete", rows=rows,
            elapsed_seconds=time.perf_counter()-started, all_case_gates_passed=True))
        audit.progress(out, "stage_complete", cases=len(rows), elapsed_seconds=time.perf_counter()-started)
    except Exception as exc:
        audit.save_json(out / "failure.json", dict(contract, status="failed", error=repr(exc), completed_rows=rows,
            elapsed_seconds=time.perf_counter()-started))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
