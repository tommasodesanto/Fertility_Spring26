#!/usr/bin/env python3
"""Bounded fixed-state price trace for the rejected tax-policy date.

This diagnoses a saved failure; it does not replace a policy path or calibration.
Run against the immutable scientific snapshot on a compute node.
"""
import argparse
import concurrent.futures as futures
import gzip
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import sys
import time

ROOT = Path(os.environ.get("E5F_DIAGNOSTIC_SOURCE_ROOT", Path(__file__).resolve().parents[3]))
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]
import numpy as np
import run_e5f_joint_nested_finalize as finalizer
import run_e5f_joint_overnight_case as adapter

STATE = None


def initialize(checkpoint, progress, checkpoint_sha):
    global STATE
    adapter.verify(checkpoint, checkpoint_sha)
    finalizer.configure_policy_model()
    assert finalizer.calibration.code_fingerprint_contract(finalizer.solver)["bundle_sha256"] == adapter.BUNDLE
    with gzip.open(checkpoint, "rb") as stream:
        packet = pickle.load(stream)
    e, P, grid, shared = [packet[k] for k in ("evaluation", "parameters", "b_grid", "shared")]
    row = adapter.read_csv(progress)[-1]
    assert int(row["calendar_year"]) == 2047 and math.isclose(P.tau_H / P.period_years, .02)
    assert e.policy.price[0] == float(row["asset_price"])
    # Exactly the transition and entrant replacement in advance_from_evaluation.
    g, _, deaths, _ = finalizer.policy.transition.advance_sequential_calendar_distribution(
        e, np.zeros(int(P.I)), P, grid, shared)
    entrants = np.asarray(P.entry_shares, dtype=float).reshape(-1).copy()
    entrants /= entrants.sum()
    entrants *= float(row["entrant_flow_next"])
    g[:, :, :, 0, :, :, :] = finalizer.policy.calendar.entrant_cohort(entrants, P, grid)
    expected = float(e.g_post_fertility.sum()) - deaths + float(row["entrant_flow_next"])
    assert abs(float(g.sum()) - expected) < 2e-10 and np.isfinite(g).all() and g.min() >= 0
    STATE = dict(packet=packet, g=g, previous=e.policy, expected_mass=expected)


def evaluate(price):
    start = time.time()
    p = STATE["packet"]
    e = finalizer.policy.calendar.evaluate_period(
        np.array([price]), STATE["g"], p["parameters"], p["b_grid"], p["shared"],
        finalizer.policy.calendar.SolveCounter(), supply_rule=p["supply_rule"])
    demand, supply = float(e.demand_by_loc[0]), float(e.supply_by_loc[0])
    return dict(price=price, demand=demand, supply=supply, excess=demand-supply,
                relative_residual=e.relative_market_residual, births=e.births,
                mass=float(e.g_current.sum()), projection=e.feasibility_projection_mass,
                elapsed_seconds=time.time()-start)


def endpoint_replay(args):
    trace = adapter.read_json(args.endpoint_trace)
    assert trace["status"] == "bounded_trace_complete"
    assert trace["checkpoint_sha256"] == args.checkpoint_sha256
    initialize(args.checkpoint, args.progress, args.checkpoint_sha256)
    p = STATE["packet"]; P = p["parameters"]; results = []; evaluations = []
    for label, row in zip(("lower", "upper"), trace["bracket"]):
        folder = args.outdir / label; folder.mkdir(parents=True, exist_ok=True)
        if args.reuse_endpoint_checkpoints:
            with gzip.open(folder/"dated_state.pkl.gz", "rb") as stream:
                saved = pickle.load(stream)
            e = saved["evaluation"]
            assert np.array_equal(e.inherited_g_pre, STATE["g"])
            assert np.array_equal(saved["b_grid"], p["b_grid"])
            assert float(e.policy.price[0]) == row["price"]
        else:
            e = finalizer.policy.calendar.evaluate_period(
                np.array([row["price"]]), STATE["g"], P, p["b_grid"], p["shared"],
                finalizer.policy.calendar.SolveCounter(), supply_rule=p["supply_rule"])
        assert e.relative_market_residual == row["relative_residual"]
        packet = dict(p, evaluation=e)
        if args.reuse_endpoint_checkpoints:
            budget = adapter.read_json(folder/"budget_summary.json")
            arrays = adapter.read_json(folder/"policy_array_summary.json")
            assert len(list((folder/"standard_diagnostics").glob("*.png"))) == 17
        else:
            with gzip.open(folder/"dated_state.pkl.gz", "wb", compresslevel=1) as stream:
                pickle.dump(packet, stream, protocol=5)
            finalizer.audit.standard_diagnostics(packet, folder, validate_production_young=False)
            budget = finalizer.audit.budget_audit(packet, folder)
            arrays = finalizer.audit.policy_array_audit(packet, folder)
        assert budget["budget_excess_mass"] <= 2e-10 and arrays["occupied_negative_steps"] == 0
        component = [float(np.sum(e.g_current[:,0] * e.policy.hR_pol[:,0]))]
        component += [float(np.sum(e.g_current[:,k+1]))*float(h) for k,h in enumerate(P.H_own)]
        assert abs(sum(component)-float(e.demand_by_loc.sum())) < 2e-10
        results.append(dict(label=label, price=row["price"], housing_by_product=component,
            budget_excess_mass=budget["budget_excess_mass"], occupied_value_drops=arrays["occupied_negative_steps"],
            checkpoint_sha256=adapter.digest(folder/"dated_state.pkl.gz")))
        evaluations.append(e)
    low, high = evaluations
    diff = np.max(np.abs(high.policy.joint_choice.probabilities-low.policy.joint_choice.probabilities),axis=(-2,-1))*STATE["g"]
    largest = []
    for flat in np.argsort(diff.ravel())[-12:][::-1]:
        index = np.unravel_index(int(flat), diff.shape)
        largest.append(dict(index=[int(x) for x in index], wealth=float(p["b_grid"][index[0]]),
            age=float(P.age_start+index[3]*P.da), inherited_mass=float(STATE["g"][index]),
            weighted_maximum_plan_probability_change=float(diff[index]),
            lower_probabilities=low.policy.joint_choice.probabilities[index].tolist(),
            upper_probabilities=high.policy.joint_choice.probabilities[index].tolist(),
            lower_products=low.policy.joint_choice.products[index].tolist(),
            upper_products=high.policy.joint_choice.products[index].tolist()))
    result=dict(status="exact_endpoint_replay", checkpoint_sha256=args.checkpoint_sha256,
        progress_sha256=adapter.digest(args.progress), trace_sha256=adapter.digest(args.endpoint_trace),
        fixed_population_mass=float(STATE["g"].sum()), owner_products=np.asarray(P.H_own).tolist(),
        phi=np.asarray(P.phi).tolist(), endpoints=results, largest_occupied_plan_changes=largest,
        interpretation="These are non-clearing diagnostic endpoints, not equilibrium solutions.", production_promoted=False)
    adapter.write_json(args.outdir/"endpoint_replay.json", result)
    print(json.dumps(result), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--checkpoint-sha256", required=True)
    ap.add_argument("--progress", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-rounds", type=int, default=10)
    ap.add_argument("--seconds", type=int, default=1500)
    ap.add_argument("--endpoint-trace", type=Path, help="Replay a completed trace's two boundary prices without searching.")
    ap.add_argument("--reuse-endpoint-checkpoints", action="store_true", help="Read and verify existing endpoint packets; perform no new solve or overwrite.")
    args = ap.parse_args()
    out = args.outdir; out.mkdir(parents=True, exist_ok=True)
    if (out / "trace.csv").exists():
        raise RuntimeError("Refusing to overwrite a diagnostic trace")
    if args.endpoint_trace:
        if (out/"endpoint_replay.json").exists():
            raise RuntimeError("Refusing to overwrite an endpoint replay")
        endpoint_replay(args); return
    assert 1 <= args.workers <= 8 and 1 <= args.max_rounds <= 10 and args.seconds <= 1500
    start = time.time(); rows = []
    def save(status, **extra):
        finalizer.policy.baseline.write_csv(out / "trace.csv", rows)
        adapter.write_json(out / "latest.json", dict(status=status, elapsed_seconds=time.time()-start,
            completed_price_evaluations=len(rows), minimum_residual=min((r["relative_residual"] for r in rows), default=None),
            checkpoint_sha256=args.checkpoint_sha256, fixed_failed_year=2051, production_promoted=False, **extra))
    # Two endpoint solves smoke-test the exact worker, collection and artifact loop.
    previous = float(adapter.read_csv(args.progress)[-1]["asset_price"])
    with futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"),
            initializer=initialize, initargs=(args.checkpoint, args.progress, args.checkpoint_sha256)) as pool:
        def batch(prices, round_number):
            submitted = [pool.submit(evaluate, float(p)) for p in prices]
            found = []
            for f in futures.as_completed(submitted):
                r = f.result(); r["round"] = round_number
                assert r["projection"] <= 1e-6 and all(math.isfinite(float(v)) for v in r.values())
                rows.append(r); found.append(r); save("running", round=round_number)
                print(json.dumps(r), flush=True)
            return sorted(found, key=lambda r: r["price"])
        endpoints = batch([previous / 1.18, previous], 0)
        if endpoints[0]["excess"] * endpoints[1]["excess"] > 0:
            save("unbracketed_diagnostic_interval"); return
        lo, hi = endpoints
        smoke_elapsed = time.time()-start
        remaining_forecast = args.max_rounds * max(r["elapsed_seconds"] for r in endpoints) * 1.5
        adapter.write_json(out / "loop_smoke.json", dict(status="pass", endpoints=2,
            worker_count=args.workers, elapsed_seconds=smoke_elapsed,
            maximum_remaining_evaluations=7*args.max_rounds,
            remaining_walltime_forecast_seconds=remaining_forecast,
            interpretation="One parallel seven-price round per refinement; no parameter search."))
        if smoke_elapsed + remaining_forecast > args.seconds:
            save("smoke_forecast_exceeds_budget"); return
        for round_number in range(1, args.max_rounds+1):
            if time.time()-start + max(r["elapsed_seconds"] for r in endpoints)*1.5 > args.seconds:
                save("time_budget", bracket=[lo,hi]); break
            grid = np.linspace(lo["price"], hi["price"], 9)[1:-1]
            middle = batch(grid, round_number)
            points = [lo,*middle,hi]
            brackets = [(a,b) for a,b in zip(points,points[1:]) if a["excess"]*b["excess"] <= 0]
            if not brackets:
                save("lost_sign_bracket", points=points); break
            lo, hi = min(brackets,key=lambda ab: ab[1]["price"]-ab[0]["price"])
            passed = [r for r in rows if r["relative_residual"] <= 2e-4]
            save("running", bracket=[lo,hi], points_meeting_market_gate=len(passed))
        else:
            save("bounded_trace_complete", bracket=[lo,hi],
                 points_meeting_market_gate=sum(r["relative_residual"]<=2e-4 for r in rows))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax, subset, label in [(axes[0],rows,"All evaluated prices"),
            (axes[1],[r for r in rows if r["round"]>=max(x["round"] for x in rows)-1],"Final two refinement rounds")]:
        ax.scatter([r["price"] for r in subset],[r["excess"]/r["supply"] for r in subset],s=24)
        ax.axhline(0,color="black",linewidth=.6)
        ax.axhline(2e-4,color="gray",linestyle="--"); ax.axhline(-2e-4,color="gray",linestyle="--")
        ax.set_title(label);ax.set_xlabel("House price");ax.set_ylabel("Signed excess demand / supply")
        ax.ticklabel_format(axis="x",useOffset=True)
    fig.suptitle("Supplemental diagnostic: fixed inherited tax-policy state, 2051")
    fig.savefig(out/"price_trace.png",dpi=180);plt.close(fig)


if __name__ == "__main__":
    main()
