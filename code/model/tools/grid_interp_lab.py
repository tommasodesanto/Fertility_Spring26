#!/usr/bin/env python3
"""Grid + interpolation lab: diagnose how a solution uses the liquid-wealth grid,
and print ready-to-run commands to re-solve on a resized grid / different
interpolation method and re-plot a policy slice.

Usage (from repo root):

  # 1. Diagnose grid usage of the current solution
  code/model/.venv/bin/python code/model/tools/grid_interp_lab.py --diagnose

  # 2. Re-solve on a dense-core + sparse-buffer grid (and/or cubic interp), then
  #    re-plot the poor-renter slice. The script prints the exact two commands;
  #    run them, or pass --run to execute them for you.
  code/model/.venv/bin/python code/model/tools/grid_interp_lab.py \
      --b-min -8 --b-max 25 --b-core-lo -8 --b-core-hi 12 --b-frac-core 0.82 \
      --interp-method monotone_cubic --name dense_cubic --run
"""
from __future__ import annotations
import argparse
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
MODEL = REPO / "code/model"
VENV = MODEL / ".venv/bin/python"
PY = str(VENV if VENV.exists() else sys.executable)
CURRENT = REPO / "output/model/intergen_model_run_current"
SOURCE = REPO / "output/model/intergen_room_distribution_current_best_20260623/summary.json"


def diagnose(cache_path: Path) -> None:
    d = pickle.load(open(cache_path, "rb"))
    sol = d["baseline"]["sol"]
    bg = np.asarray(sol.b_grid, float).reshape(-1)
    g = np.asarray(sol.g, float)
    mb = g.sum(axis=tuple(range(1, g.ndim)))
    occ = mb > 1e-9
    cm = np.cumsum(mb) / mb.sum()
    print(f"grid: [{bg[0]:.2f}, {bg[-1]:.2f}], Nb={bg.size}")
    print(f"occupied (mass>1e-9): {int(occ.sum())}/{bg.size}  range [{bg[occ][0]:.2f}, {bg[occ][-1]:.2f}]")
    qs = {q: float(bg[min(np.searchsorted(cm, q), bg.size - 1)]) for q in (0.001, 0.01, 0.5, 0.99, 0.999)}
    print("wealth quantiles:  " + "  ".join(f"{int(q*1000)/10}%={v:.2f}" for q, v in qs.items()))
    # per-region spacing
    diffs = np.diff(bg)
    regions = [("[-8,-3)", -8, -3), ("[-3,0)", -3, 0), ("[0,3)", 0, 3),
               ("[3,6)", 3, 6), ("[6,10)", 6, 10), ("[10,bmax]", 10, bg[-1] + 1)]
    print("median node spacing by region:")
    for name, lo, hi in regions:
        sel = (bg[:-1] >= lo) & (bg[:-1] < hi)
        if sel.any():
            print(f"   {name:>11}: dB={np.median(diffs[sel]):.3f}  ({int(sel.sum())} gaps)")
    # how much mass sits in the COARSE upper-middle vs the dense core
    frac_above6 = float(mb[bg >= 6].sum() / mb.sum())
    print(f"mass at b>=6 (currently coarse 1.56-spaced region): {frac_above6:.3f}")
    # suggested bounds: a hair below deepest reachable mortgage; a buffer above 99.9%ile
    p = float(np.asarray(d['baseline']['p_eq']).reshape(-1)[0])
    phi = float(np.asarray(d['baseline']['P'].phi).reshape(-1)[0])
    Hmax = float(np.max(np.asarray(d['baseline']['P'].H_own)))
    deepest = -(1 - phi) * 0 - phi * p * Hmax  # = -phi*p*Hmax (max-LTV owner)
    print(f"deepest reachable b' (max mortgage -phi*p*Hmax): {deepest:.2f}")
    print(f"SUGGESTED dense core ~[{deepest-1.5:.0f}, {qs[0.999]+3:.0f}], sparse buffer to ~{max(25, qs[0.999]*2):.0f}")


def resolve_cmds(args) -> tuple[list[str], list[str]]:
    outdir = REPO / "output/model" / f"intergen_grid_lab_{args.name}"
    packet = [PY, str(MODEL / "tools/build_intergen_mechanics_packet.py"),
              "--source", str(SOURCE), "--outdir", str(outdir),
              "--target-set", "candidate_replacement_young_old_roomgap_v1",
              "--J", "16", "--Nb", str(args.nb), "--income-states", "5",
              "--H-own", "2,4,6,8,10", "--hR-max", "6.0", "--max-iter-eq", str(args.max_iter_eq),
              "--skip-contact-sheet"]
    for flag, val in [("--b-min", args.b_min), ("--b-max", args.b_max),
                      ("--b-core-lo", args.b_core_lo), ("--b-core-hi", args.b_core_hi),
                      ("--b-mid-hi", args.b_mid_hi), ("--b-frac-low", args.b_frac_low),
                      ("--b-frac-core", args.b_frac_core), ("--b-frac-mid", args.b_frac_mid),
                      ("--interp-method", args.interp_method)]:
        if val is not None:
            packet += [flag, str(val)]
    plot = [PY, str(MODEL / "tools/plot_intergen_policy_slice.py"),
            "--cache", str(outdir / "solution_cache.pkl"), "--outdir", str(outdir),
            "--age", str(args.slice_age), "--z", str(args.slice_z),
            "--tenure", "renter", "--parity", "0", "--child-state", "0",
            "--b-min", "-1.5", "--b-max", "4", "--name", f"{args.name}_slice"]
    return packet, plot, outdir


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diagnose", action="store_true", help="Just report grid usage of --cache and exit.")
    ap.add_argument("--cache", type=Path, default=CURRENT / "solution_cache.pkl")
    ap.add_argument("--name", default="resized")
    ap.add_argument("--nb", type=int, default=60)
    ap.add_argument("--max-iter-eq", type=int, default=10)
    ap.add_argument("--b-min", type=float, default=None)
    ap.add_argument("--b-max", type=float, default=None)
    ap.add_argument("--b-core-lo", type=float, default=None)
    ap.add_argument("--b-core-hi", type=float, default=None)
    ap.add_argument("--b-mid-hi", type=float, default=None)
    ap.add_argument("--b-frac-low", type=float, default=None)
    ap.add_argument("--b-frac-core", type=float, default=None)
    ap.add_argument("--b-frac-mid", type=float, default=None)
    ap.add_argument("--interp-method", default=None, choices=["linear", "monotone_cubic"])
    ap.add_argument("--slice-age", type=float, default=30.0)
    ap.add_argument("--slice-z", type=float, default=1.0)
    ap.add_argument("--run", action="store_true", help="Execute the re-solve + re-plot (slow) instead of just printing.")
    args = ap.parse_args()

    print("=== grid usage of:", args.cache, "===")
    diagnose(args.cache)
    if args.diagnose:
        return
    packet, plot, outdir = resolve_cmds(args)
    print("\n=== re-solve command ===\n" + " ".join(packet))
    print("\n=== re-plot slice command ===\n" + " ".join(plot))
    if args.run:
        print("\n[running re-solve...]")
        subprocess.run(packet, cwd=str(MODEL), check=True)
        print("\n[re-plotting slice...]")
        subprocess.run(plot, cwd=str(MODEL), check=True)
        print("\n=== grid usage of the RESIZED solution ===")
        diagnose(outdir / "solution_cache.pkl")


if __name__ == "__main__":
    main()
