#!/usr/bin/env python3
"""Report-only initial-state group decomposition for native financing cases.

The script reuses the saved native policy arrays and calendar evaluator.  It
never solves a household problem.  Groups are initial tenure (renter versus
owner) crossed with four common, pooled wealth quartiles among fertile-age,
childless households.  Housing outcomes are therefore mapped from each
initial subgroup, rather than inferred from a realized-tenure policy average.
"""
from __future__ import annotations

import argparse, csv, gzip, json, pickle, sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from run_e5f_native_financing_diagnostic import (
    CASES, CHECKPOINT_SHA256, POLICY_NAMES, arm_parameters, install_paths,
    policy_arrays, validate_contract,
)

CASE_ORDER = ("baseline", "mortgage_only", "unsecured_only", "both")

AGE_LO, AGE_HI = 26.0, 38.0
ROOM_THRESHOLD = 6.0
ATOL = 3e-10


def _age(meta: dict[str, Any], j: int) -> float:
    return float(meta["age_start"] + j * meta["da"])


def _case_dirs(root: Path) -> dict[str, Path]:
    found = {}
    for p in sorted((root / "cases").glob("*/arrays.npz")):
        for case in CASE_ORDER:
            if p.parent.name == case or p.parent.name.startswith(case + "_"):
                found.setdefault(case, p.parent)
    return found


def _load_case(path: Path) -> dict[str, np.ndarray]:
    rec = json.loads((path / "receipt.json").read_text())
    if rec.get("status") != "completed":
        raise ValueError(f"case receipt is not completed: {rec.get('status')}")
    if rec.get("contract", {}).get("checkpoint_sha256") != CHECKPOINT_SHA256:
        raise ValueError("case checkpoint hash mismatch")
    with np.load(path / "arrays.npz", allow_pickle=False) as z:
        required = {"g_pre", "g_post_fertility", "g_current", "births", *POLICY_NAMES}
        missing = sorted(required - set(z.files))
        if missing: raise ValueError(f"missing mandatory arrays: {missing}")
        return {k: np.asarray(z[k]) for k in z.files}


def _weighted_quantiles(values: np.ndarray, weights: np.ndarray, probs=(.25, .5, .75)) -> np.ndarray:
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    order = np.argsort(values, kind="mergesort"); values, weights = values[order], weights[order]
    total = float(weights.sum())
    if total <= 0: raise ValueError("eligible childless pool has zero mass")
    cdf = np.cumsum(weights) / total
    return np.asarray([values[np.searchsorted(cdf, p, side="left")] for p in probs])


def _age_sets(P: Any, meta: dict[str, Any], J: int) -> tuple[np.ndarray, np.ndarray]:
    # Native convention is one-based A_f_start/A_f_end, represented by Python
    # indices [A_f_start-1, A_f_end).
    lo = int(getattr(P, "A_f_start")); hi = int(getattr(P, "A_f_end"))
    fertile = np.arange(max(0, lo - 1), min(J, hi), dtype=int)
    age2638 = np.asarray([j for j in fertile if AGE_LO <= _age(meta, j) <= AGE_HI], dtype=int)
    return fertile, age2638


def _eligible_mask(shape: tuple[int, ...], ages: np.ndarray) -> np.ndarray:
    m = np.zeros(shape, dtype=bool)
    for j in np.asarray(ages, dtype=int): m[:, :, :, j, :, 0, :] = True
    return m


def group_definition(g: np.ndarray, b_grid: np.ndarray, fertile: np.ndarray, age2638: np.ndarray, *, age_start: float = 18.0, da: float = 4.0) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    eligible = _eligible_mask(g.shape, fertile)
    pooled = np.take(g, fertile, axis=3)[:, :, :, :, :, 0, :].sum(axis=(1, 2, 3, 4, 5))
    # b_grid values are repeated across all non-wealth eligible states.
    cuts = _weighted_quantiles(np.asarray(b_grid), pooled)
    groups = []
    for tenure, label in ((0, "renter"), (1, "owner")):
        for q in range(4):
            lo = -np.inf if q == 0 else cuts[q - 1]
            hi = np.inf if q == 3 else cuts[q]
            w = np.asarray(b_grid)
            wmask = (w >= lo) & ((w <= hi) if q == 3 else (w < hi))
            mask = np.zeros_like(g, dtype=bool)
            ts = tenure if tenure == 0 else slice(1, g.shape[1])
            for j in np.asarray(fertile, dtype=int): mask[wmask, ts, :, j, :, 0, :] = True
            groups.append({"group": f"{label}_q{q+1}", "initial_tenure": label, "wealth_quartile": q+1, "mask": mask, "age2638_mask": _eligible_mask(g.shape, age2638) & mask})
    # The owner slice is a single group across all owner products; assert this
    # is exhaustive/disjoint over the eligible pool.
    union = np.zeros_like(g, dtype=bool)
    for x in groups:
        if np.any(union & x["mask"]): raise ValueError("group masks overlap")
        union |= x["mask"]
    if not np.array_equal(union, eligible): raise ValueError("group masks do not exhaust eligible pool")
    return groups, {"fertile_indices": fertile.tolist(), "age26_38_indices": age2638.tolist(), "age26_38": [float(age_start + da*j) for j in age2638], "wealth_quartile_cuts": cuts.tolist(), "duplicate_cuts": bool(np.any(np.diff(cuts) == 0)), "eligible_mass": float(g[eligible].sum())}


def _first(gpre: np.ndarray, gpost: np.ndarray, ages: np.ndarray) -> float:
    pre = np.take(gpre, ages, axis=3)[:, :, :, :, :, 0, :]
    post = np.take(gpost, ages, axis=3)[:, :, :, :, :, 0, :]
    loss = float(pre.sum() - post.sum())
    if loss < -3e-10:
        raise ValueError(f"n=0 mass increases at selected ages: {loss:.17g}")
    return max(0.0, loss)


def _housing(gc: np.ndarray, h: np.ndarray, P: Any) -> dict[str, float]:
    renter = float(gc[:, 0, ...].sum()); owner = float(gc[:, 1:, ...].sum())
    rooms_r = float(np.sum(gc[:, 0, ...] * h[:, 0, ...]))
    H = np.asarray(getattr(P, "H_own"))
    rooms_o = float(sum(gc[:, t, ...].sum() * H[t-1] for t in range(1, gc.shape[1])))
    ge6_r = float(np.sum(gc[:, 0, ...] * (h[:, 0, ...] >= ROOM_THRESHOLD)))
    ge6_o = float(sum(gc[:, t, ...].sum() * (H[t-1] >= ROOM_THRESHOLD) for t in range(1, gc.shape[1])))
    return {"renter_mass": renter, "owner_mass": owner, "rooms_total": rooms_r + rooms_o,
            "mass_renting_ge6_rooms": ge6_r, "mass_owning_ge6_rooms": ge6_o}


def _native_map(subset: np.ndarray, case: str, a: dict[str, np.ndarray], packet: dict[str, Any], source_root: Path) -> Any:
    install_paths(source_root)
    primitive = __import__("run_e5f_matched_pf_smoke")
    model = __import__("intergen_eqscale_seq_optimized.solver", fromlist=["x"])
    P = arm_parameters(packet["parameters"], case)
    grid = np.asarray(packet["b_grid"]); price = np.asarray(a["price"])
    primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    shared = model.precompute_shared(P, grid)
    sol = SimpleNamespace(**{k: a[k] for k in POLICY_NAMES if k != "price"})
    P._fert2_probs = np.asarray(a["fert2_probs"]).copy()
    policy = primitive.pf.calendar.policy_from_solution(sol, price, P, grid, shared)
    return primitive.pf.calendar.evaluate_period(price, subset, P, grid, shared, primitive.pf.calendar.SolveCounter(), supply_rule=packet["supply_rule"], supplied_policy=policy), P


def analyze_case(case: str, a: dict[str, np.ndarray], packet: dict[str, Any], source_root: Path, meta: dict[str, Any], groups: list[dict[str, Any]], fertile: np.ndarray, age2638: np.ndarray) -> list[dict[str, Any]]:
    rows = []
    for x in groups:
        mask = x["mask"]; sub = np.where(mask, a["g_pre"], 0.0)
        ev, P = _native_map(sub, case, a, packet, source_root)
        if not np.array_equal(np.asarray(ev.g_pre), sub): raise ValueError(f"{case}/{x['group']} mapper changed subgroup g_pre")
        gc = np.asarray(ev.g_current); mass = float(sub.sum()); current = float(gc.sum())
        h = np.asarray(a["hR_pol"]); renter = float(gc[:, 0, ...].sum()); owner = float(gc[:, 1:, ...].sum())
        rooms_r = float(np.sum(gc[:, 0, ...] * h[:, 0, ...]))
        H = np.asarray(getattr(P, "H_own")); rooms_o = float(sum(gc[:, t, ...].sum() * H[t-1] for t in range(1, gc.shape[1])))
        ge6_r = float(np.sum(gc[:, 0, ...] * (h[:, 0, ...] >= ROOM_THRESHOLD)))
        ge6_o = float(sum(gc[:, t, ...].sum() * (H[t-1] >= ROOM_THRESHOLD) for t in range(1, gc.shape[1])))
        first = _first(sub, np.asarray(ev.g_post_fertility), fertile)
        first2638 = _first(sub, np.asarray(ev.g_post_fertility), age2638)
        rows.append({"case":case,"group":x["group"],"initial_tenure":x["initial_tenure"],"wealth_quartile":x["wealth_quartile"],"initial_mass":mass,"current_mass":current,"birth_flow":float(np.asarray(ev.births).sum()),"firstbirth_flow":first,"firstbirth_rate":first/mass if mass else np.nan,"firstbirth_26_38":first2638,"mean_rooms":(rooms_r+rooms_o)/current if current else np.nan,"ownership_rate":owner/current if current else np.nan,"owner_mass":owner,"rooms_total":rooms_r+rooms_o,"mass_renting_ge6_rooms":ge6_r,"mass_owning_ge6_rooms":ge6_o,"response_difference_baseline":np.nan,"firstbirth_rate_difference_baseline":np.nan,"aggregate_response_contribution":np.nan})
    return rows


def build(args: argparse.Namespace) -> None:
    args.output.mkdir(parents=True, exist_ok=True)
    contract = validate_contract(args.checkpoint, args.replay, args.source_root)
    install_paths(args.source_root)
    with gzip.open(args.checkpoint, "rb") as f: packet = pickle.load(f)
    P = packet["parameters"]; meta = {"age_start":float(P.age_start),"da":float(P.da)}
    dirs = _case_dirs(args.input); all_rows=[]; definition=None; aggregate={}
    for case in CASE_ORDER:
        if case not in dirs: raise ValueError(f"missing completed case: {case}")
        a = _load_case(dirs[case])
        if not np.array_equal(a["g_pre"], packet["stationary_g_pre"]): raise ValueError(f"{case} g_pre differs from checkpoint exactly")
        fertile, age2638 = _age_sets(P, meta, a["g_pre"].shape[3])
        groups, definition = group_definition(a["g_pre"], packet["b_grid"], fertile, age2638, age_start=meta["age_start"], da=meta["da"])
        rows = analyze_case(case,a,packet,args.source_root,meta,groups,fertile,age2638); all_rows.extend(rows)
        elig = _eligible_mask(a["g_pre"].shape, fertile); pooled=np.where(elig,a["g_pre"],0.)
        ev,_ = _native_map(pooled,case,a,packet,args.source_root)
        hs = _housing(np.asarray(ev.g_current), np.asarray(a["hR_pol"]), arm_parameters(packet["parameters"], case))
        aggregate[case]={"eligible_mass":float(pooled.sum()),"current_mass":float(np.asarray(ev.g_current).sum()),"birth_flow":float(np.asarray(ev.births).sum()),"firstbirth_flow":_first(pooled,np.asarray(ev.g_post_fertility),fertile),"firstbirth_26_38":_first(pooled,np.asarray(ev.g_post_fertility),age2638), **hs}
        rr=[r for r in rows]; checks={k:float(sum(r[k] for r in rr)) for k in ("current_mass","birth_flow","firstbirth_flow","firstbirth_26_38")}
        checks.update({k:float(sum(r[k] for r in rr)) for k in ("owner_mass","rooms_total","mass_renting_ge6_rooms","mass_owning_ge6_rooms")})
        for k,v in checks.items():
            if not np.isclose(v, aggregate[case][k], rtol=0, atol=2e-8): raise ValueError(f"grouped {k} does not reconcile for {case}: {v} vs {aggregate[case][k]}")
        aggregate[case]["max_group_reconciliation_gap"] = max(abs(v-aggregate[case][k]) for k,v in checks.items())
        # The saved full native array provides the independent first-birth check.
        full_native = _first(np.asarray(a["g_pre"]), np.asarray(a["g_post_fertility"]), fertile)
        if not np.isclose(full_native, aggregate[case]["firstbirth_flow"], rtol=0, atol=2e-8):
            raise ValueError(f"{case} eligible first-birth pool does not match saved full-array native flow")
    by={(r["case"],r["group"]):r for r in all_rows}
    base={r["group"]:r for r in all_rows if r["case"]=="baseline"}
    for r in all_rows:
        d=r["firstbirth_flow"]-base[r["group"]]["firstbirth_flow"]; r["response_difference_baseline"]=d
        r["firstbirth_rate_difference_baseline"] = r["firstbirth_rate"]-base[r["group"]]["firstbirth_rate"]
        r["aggregate_response_contribution"]=d / max(aggregate[r["case"]]["eligible_mass"], 1e-15)
    # Age-specific first-birth table is calculated directly from native g_pre/g_post.
    age_rows=[]
    for case in CASE_ORDER:
        a = _load_case(dirs[case]); fertile, _ = _age_sets(P, meta, a["g_pre"].shape[3])
        for j in fertile:
            den=float(np.asarray(a["g_pre"][:, :, :, j, :, 0, :]).sum()); flow=_first(a["g_pre"], a["g_post_fertility"], np.asarray([j]))
            age_rows.append({"case":case,"age_index":int(j),"age":float(_age(meta,j)),"initial_mass":den,"firstbirth_flow":flow,"firstbirth_rate":flow/den if den else np.nan})
    with (args.output/"age_firstbirth.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(age_rows[0]),lineterminator="\n"); w.writeheader(); w.writerows(age_rows)
    with (args.output/"groups.csv").open("w",newline="") as f:
        fields=list(all_rows[0]); w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n"); w.writeheader(); w.writerows(all_rows)
    (args.output/"metadata.json").write_text(json.dumps({"contract":contract,"groups":definition,"aggregate_reconciliation":aggregate,"room_threshold":ROOM_THRESHOLD,"scope":"fertile-age initial childless pool; age 26-38 subtotal reported separately"},indent=2)+"\n")
    import matplotlib.pyplot as plt
    fig, ax=plt.subplots(1,2,figsize=(9,4),constrained_layout=True)
    for case in CASE_ORDER:
        vals=[next(r["firstbirth_rate_difference_baseline"] for r in all_rows if r["case"]==case and r["group"]==f"renter_q{q}") for q in range(1,5)]
        ax[0].plot(range(1,5),vals,marker="o",label=f"{case} renter")
        vals=[next(r["firstbirth_rate_difference_baseline"] for r in all_rows if r["case"]==case and r["group"]==f"owner_q{q}") for q in range(1,5)]
        ax[0].plot(range(1,5),vals,marker="s",linestyle="--",label=f"{case} owner")
    for ten in ("renter","owner"):
        vals=[next(r["mean_rooms"] for r in all_rows if r["case"]=="baseline" and r["group"]==f"{ten}_q{q}") for q in range(1,5)]
        ax[1].plot(range(1,5),vals,marker="o",label=ten)
    ax[0].set_title("First-birth-rate response by initial group"); ax[1].set_title("Baseline mean physical rooms")
    for axy in ax: axy.set_xlabel("pooled wealth quartile"); axy.legend()
    fig.savefig(args.output/"eligible_groups_summary.png",dpi=160); plt.close(fig)


def main(argv=None):
    p=argparse.ArgumentParser(); p.add_argument("--input",type=Path,required=True); p.add_argument("--output",type=Path,required=True); p.add_argument("--checkpoint",type=Path,required=True); p.add_argument("--replay",type=Path,required=True); p.add_argument("--source-root",type=Path,required=True); build(p.parse_args(argv))

if __name__ == "__main__": main()
