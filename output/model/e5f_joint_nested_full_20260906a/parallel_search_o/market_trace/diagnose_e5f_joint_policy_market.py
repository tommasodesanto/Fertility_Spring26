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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--checkpoint-sha256", required=True)
    ap.add_argument("--progress", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-rounds", type=int, default=10)
    ap.add_argument("--seconds", type=int, default=1500)
    args = ap.parse_args()
    out = args.outdir; out.mkdir(parents=True, exist_ok=True)
    if (out / "trace.csv").exists():
        raise RuntimeError("Refusing to overwrite a diagnostic trace")
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
