"""Validate surrogate/BO proposals against the REAL model solver.

Two checks, both using the actual ``run_model_cp_dt`` solve path at the
diagnostic grid (J=16, Nb=60, Nz=5 income states, n_house=5, max_iter_eq=3):

1. Reproduction check: re-solve a few existing data records at their stored theta
   and confirm we reproduce their stored moments. This proves the override
   construction here matches the pipeline that generated the training data.

2. Proposal validation: solve the BO/exploit proposals and compare the surrogate
   *predicted* loss/moments to the *actual* solved loss/moments. This confirms
   the emulator is calibrated and that its picks beat the data-set incumbent.

Each solve takes ~40-50s, so this is meant to be run on a small number of points
(``--limit``) and is happy in the background.
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from intergen_housing_fertility.solver import run_model_cp_dt
from intergen_housing_fertility.calibration import (
    base_overrides, extract_moments, diagnostic_loss,
)
from intergen_housing_fertility.local_panel import income_process_overrides

from . import data as D

OUT_ROOT = os.path.join(D.REPO_ROOT, "output/model/intergen_surrogate_calibration")


def solve_theta(theta: dict, J=16, Nb=60, n_house=5, max_iter_eq=3,
                income_states=5) -> dict:
    overrides = {
        **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
        **income_process_overrides(income_states),
        **{p: float(theta[p]) for p in D.THETA_PARAMS},
    }
    t0 = time.time()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    return {
        "moments": {k: (float(v) if np.isfinite(_f(v)) else None)
                    for k, v in moments.items()},
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "p_eq": [float(x) for x in np.atleast_1d(p_eq)],
        "elapsed_sec": time.time() - t0,
    }


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def reproduction_check(ds: D.Dataset, recs: list[dict], n: int = 2, seed: int = 0):
    """Re-solve a few records and compare to their stored moments."""
    rng = np.random.default_rng(seed)
    # Choose a couple of low-loss records (most relevant region).
    order = np.argsort(ds.loss)
    pick = [int(order[0]), int(order[len(order) // 2])][:n]
    # Map dataset rows back to records: rebuild theta from ds.X.
    results = []
    for row in pick:
        theta = {p: float(ds.X[row, k]) for k, p in enumerate(ds.theta_params)}
        solved = solve_theta(theta)
        stored = {m: float(ds.M[m][row]) for m in ds.moment_names}
        recomputed = {m: solved["moments"].get(m) for m in ds.moment_names}
        diffs = {m: (None if recomputed[m] is None
                     else abs(recomputed[m] - stored[m])) for m in ds.moment_names}
        max_abs = max((d for d in diffs.values() if d is not None), default=None)
        results.append({
            "row": row, "theta": theta, "stored": stored,
            "recomputed": recomputed, "max_abs_moment_diff": max_abs,
            "elapsed_sec": solved["elapsed_sec"],
            "market_residual": solved["market_residual"],
        })
        print(f"  repro row {row}: max abs moment diff = "
              f"{max_abs:.2e} ({solved['elapsed_sec']:.0f}s)")
    return results


def validate(ds: D.Dataset, proposals: list[dict], tag: str):
    rows = []
    for i, prop in enumerate(proposals):
        theta = prop["theta"]
        solved = solve_theta(theta)
        m = solved["moments"]
        m_full = dict(m); m_full["market_residual"] = solved["market_residual"]
        actual_loss = float(diagnostic_loss(m_full, targets=ds.targets,
                                            weights=ds.weights))
        pred_loss = prop.get("pred_loss")
        rows.append({
            "tag": tag, "i": i,
            "pred_loss": pred_loss,
            "actual_loss": actual_loss,
            "market_residual": solved["market_residual"],
            "elapsed_sec": solved["elapsed_sec"],
            "actual_moments": {mm: m.get(mm) for mm in ds.moment_names},
            "pred_moments": prop.get("pred_moments"),
            "theta": theta,
        })
        print(f"  {tag}[{i}] pred_loss={pred_loss:.3f}  "
              f"actual_loss={actual_loss:.3f}  "
              f"resid={solved['market_residual']:.1e}  "
              f"({solved['elapsed_sec']:.0f}s)")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-set", default="candidate_no_timing_v0")
    ap.add_argument("--limit", type=int, default=4,
                    help="max EI-batch proposals to solve")
    ap.add_argument("--repro", type=int, default=2,
                    help="number of reproduction-check solves")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    outdir = os.path.join(OUT_ROOT, args.target_set)
    recs = D.load_records()
    ds = D.build_dataset(args.target_set, recs)

    print(f"[validate] target_set={args.target_set}  "
          f"data-best loss={ds.loss.min():.3f}")

    out = {"target_set": args.target_set, "data_best_loss": float(ds.loss.min())}

    if args.repro > 0:
        print(f"[validate] reproduction check ({args.repro}) ...")
        out["reproduction"] = reproduction_check(ds, recs, n=args.repro,
                                                 seed=args.seed)

    prop_path = os.path.join(outdir, "bo_proposals.json")
    bo = json.load(open(prop_path))
    ei_batch = bo.get("global_ei_batch", bo.get("ei_batch", []))[: args.limit]
    local_batch = bo.get("local_batch", [])[: args.limit]
    exploit = bo["exploit_point"]

    print(f"[validate] solving exploit point ...")
    out["exploit_validation"] = validate(ds, [exploit], tag="exploit")
    print(f"[validate] solving {len(ei_batch)} GLOBAL EI proposals ...")
    out["ei_validation"] = validate(ds, ei_batch, tag="global_ei")
    if local_batch:
        print(f"[validate] solving {len(local_batch)} LOCAL trust-region proposals ...")
        out["local_validation"] = validate(ds, local_batch, tag="local")

    # Summary: did any proposal beat the incumbent? emulator calibration.
    all_val = (out["exploit_validation"] + out["ei_validation"]
               + out.get("local_validation", []))
    actual = np.array([r["actual_loss"] for r in all_val], float)
    pred = np.array([r["pred_loss"] for r in all_val], float)

    def _mode_best(rows):
        a = np.array([r["actual_loss"] for r in rows], float)
        return float(np.nanmin(a)) if a.size else None

    out["summary"] = {
        "incumbent_loss": float(ds.loss.min()),
        "best_proposal_actual_loss": float(np.nanmin(actual)),
        "n_beat_incumbent": int(np.sum(actual < ds.loss.min())),
        "pred_vs_actual_mae": float(np.nanmean(np.abs(pred - actual))),
        "pred_vs_actual_corr": (float(np.corrcoef(pred, actual)[0, 1])
                                if np.sum(np.isfinite(pred + actual)) > 2 else None),
        "best_global_ei_actual": _mode_best(out["ei_validation"]),
        "best_local_actual": _mode_best(out.get("local_validation", [])),
        "exploit_actual": _mode_best(out["exploit_validation"]),
    }
    print(f"[validate] best proposal actual loss = "
          f"{out['summary']['best_proposal_actual_loss']:.3f} "
          f"(incumbent {ds.loss.min():.3f}); "
          f"{out['summary']['n_beat_incumbent']} beat incumbent")

    json.dump(out, open(os.path.join(outdir, "validation.json"), "w"), indent=1,
              default=lambda o: o.item() if hasattr(o, "item") else str(o))

    # Markdown.
    md = [f"# Real-solver validation: `{args.target_set}`\n"]
    md.append(f"- Incumbent (data-set best) loss: **{ds.loss.min():.3f}**\n")
    md.append(f"- Best proposal actual loss: "
              f"**{out['summary']['best_proposal_actual_loss']:.3f}**\n")
    md.append(f"- Proposals beating incumbent: "
              f"**{out['summary']['n_beat_incumbent']}/{len(all_val)}**\n")
    md.append(f"- Surrogate predicted-vs-actual loss MAE: "
              f"**{out['summary']['pred_vs_actual_mae']:.3f}**, "
              f"corr: **{out['summary']['pred_vs_actual_corr']}**\n")
    if "reproduction" in out:
        md.append("\n## Reproduction check\n")
        md.append("| row | max abs moment diff vs stored | residual |\n|---|---:|---:|\n")
        for r in out["reproduction"]:
            md.append(f"| {r['row']} | {r['max_abs_moment_diff']:.2e} | "
                      f"{r['market_residual']:.1e} |\n")
    md.append("\n## Proposal validation\n")
    md.append("| tag | pred loss | actual loss | residual |\n|---|---:|---:|---:|\n")
    for r in all_val:
        md.append(f"| {r['tag']}[{r['i']}] | {r['pred_loss']:.3f} | "
                  f"{r['actual_loss']:.3f} | {r['market_residual']:.1e} |\n")
    with open(os.path.join(outdir, "VALIDATION.md"), "w") as fh:
        fh.write("".join(md))
    print(f"[validate] wrote {os.path.join(outdir, 'validation.json')}")


if __name__ == "__main__":
    main()
