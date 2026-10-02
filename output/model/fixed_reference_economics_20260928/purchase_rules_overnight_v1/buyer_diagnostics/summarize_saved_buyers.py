"""Read-only, zero-solve buyer diagnostics from a native solution_arrays.npz.

The asset state is net liquid wealth. Thus (Q-A)/Q is an *implied net funding
ratio*, not a separately observed mortgage loan-to-value ratio.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    return float(v[np.searchsorted(np.cumsum(w), probability * np.sum(w), side="left")])


def summarize(values: np.ndarray, weights: np.ndarray) -> dict:
    positive = weights > 0
    values = values[positive]
    weights = weights[positive]
    if not len(values):
        return {"mass": 0.0}
    total = float(np.sum(weights))
    return {
        "mass": total,
        "mean": float(np.sum(weights * values) / total),
        "p10": weighted_quantile(values, weights, 0.10),
        "p25": weighted_quantile(values, weights, 0.25),
        "p50": weighted_quantile(values, weights, 0.50),
        "p75": weighted_quantile(values, weights, 0.75),
        "p90": weighted_quantile(values, weights, 0.90),
        "share_above_80_percent": float(np.sum(weights[values > 0.8]) / total),
        "share_above_90_percent": float(np.sum(weights[values > 0.9]) / total),
        "share_above_100_percent": float(np.sum(weights[values > 1.0]) / total),
        "mass_above_100_percent_equivalently_negative_A": float(np.sum(weights[values > 1.0])),
        "cdf": [
            {"net_ltv": q, "share_at_or_below": float(np.sum(weights[values <= q]) / total)}
            for q in (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2)
        ],
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path, required=True, help="Native selected stage solution_arrays.npz")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--owner-rungs", required=True, help="Comma-separated physical owner unit sizes")
    ap.add_argument("--sale-cost", type=float, required=True, help="Fractional sale cost psi")
    ap.add_argument("--age-start", type=float, default=18)
    ap.add_argument("--age-step", type=float, default=4)
    ap.add_argument("--young-min", type=float, default=25)
    ap.add_argument("--young-max", type=float, default=34)
    ap.add_argument("--first-birth-flow", type=Path,
                    help="Optional exact pre-tenure first-birth flow, key first_birth_flow, shape of g")
    ap.add_argument("--pre-birth", type=Path,
                    help="Optional saved pre-fertility distribution NPZ, key g_pre, shape of g")
    ap.add_argument("--fecundity", type=float,
                    help="Required with --pre-birth: conception success probability (constant across ages)")
    ap.add_argument("--fecundity-by-age", help="Comma-separated native get_fecundity_by_age(P) vector")
    ap.add_argument("--hard-phi", type=float,
                    help="With exact first-birth flow, report strict down-payment ineligibility at this phi")
    args = ap.parse_args()
    rungs = np.array([float(x) for x in args.owner_rungs.split(",")], dtype=float)
    if not len(rungs) or np.any(rungs <= 0) or np.any(np.diff(rungs) <= 0):
        raise ValueError("Owner rungs must be positive and strictly increasing")
    if not 0 <= args.sale_cost < 1:
        raise ValueError("Sale cost must be in [0,1)")
    with np.load(args.stage, allow_pickle=False) as saved:
        g = saved["g_beginning_distribution"]
        probs = saved["tenure_probs"]
        fert_probs = saved["fert_probs"]
        b = saved["b_grid"]
        price = saved["p_eq"]
        if "shared.phi_choice" in saved:
            phi_choice = saved["shared.phi_choice"]
        else:
            phi_choice = None
    if g.ndim != 7 or probs.shape != g.shape + (g.shape[1],):
        raise ValueError("Unexpected distribution/tenure probability shapes")
    if len(rungs) != g.shape[1] - 1 or len(b) != g.shape[0] or len(price) != g.shape[2]:
        raise ValueError("Rungs, wealth grid, price or solution dimensions disagree")
    if np.any(g < -1e-13) or not np.isfinite(g).all():
        raise ValueError("Invalid pre-tenure distribution")
    if np.any(probs < -1e-6) or not np.isfinite(probs).all():
        raise ValueError("Invalid tenure probabilities")
    if np.max(np.abs(np.sum(probs, axis=-1)[g > 1e-12] - 1)) > 2e-5:
        raise ValueError("Occupied tenure probabilities do not sum to one")
    ages = args.age_start + args.age_step * np.arange(g.shape[3])
    young = (ages >= args.young_min) & (ages <= args.young_max)
    if not np.any(young):
        raise ValueError("Young age window contains no model cell")
    groups = {"renter_to_owner_buyers": ([], []), "owner_switchers": ([], []),
              "all_purchasers": ([], [])}
    for origin in range(g.shape[1]):
        old_q = 0.0 if origin == 0 else price[None, :] * rungs[origin - 1]
        for dest in range(1, g.shape[1]):
            if origin == dest:
                continue  # Stayers do not close a new purchase.
            weight = (g[:, origin, :, :, :, :, :][:, :, young]
                      * probs[:, origin, :, :, :, :, :, dest][:, :, young])
            cash = b[:, None] if origin == 0 else b[:, None] + (1 - args.sale_cost) * old_q
            q = price[None, :] * rungs[dest - 1]
            ratio = np.maximum(0.0, (q - cash) / q)
            ratio = np.broadcast_to(ratio[:, :, None, None, None, None], weight.shape)
            values = ratio.ravel()
            masses = weight.ravel()
            if np.sum(masses) <= 0:
                continue
            label = "renter_to_owner_buyers" if origin == 0 else "owner_switchers"
            groups[label][0].append(values)
            groups[label][1].append(masses)
            groups["all_purchasers"][0].append(values)
            groups["all_purchasers"][1].append(masses)
    summaries = {
        key: summarize(np.concatenate(v[0]), np.concatenate(v[1])) if v[0] else {"mass": 0.0}
        for key, v in groups.items()
    }
    all_age_negative = {"renter_to_owner_buyers": {"buyer_mass": 0.0, "negative_A_mass": 0.0},
                        "owner_switchers": {"buyer_mass": 0.0, "negative_A_mass": 0.0}}
    for origin in range(g.shape[1]):
        label = "renter_to_owner_buyers" if origin == 0 else "owner_switchers"
        old_q = 0.0 if origin == 0 else price[None, :] * rungs[origin - 1]
        cash = b[:, None] if origin == 0 else b[:, None] + (1 - args.sale_cost) * old_q
        negative = (cash < 0)[:, :, None, None, None, None]
        for dest in range(1, g.shape[1]):
            if origin == dest:
                continue
            flow = g[:, origin] * probs[:, origin, :, :, :, :, :, dest]
            all_age_negative[label]["buyer_mass"] += float(np.sum(flow))
            all_age_negative[label]["negative_A_mass"] += float(np.sum(flow * negative))
    for values in all_age_negative.values():
        values["negative_A_share"] = (
            values["negative_A_mass"] / values["buyer_mass"] if values["buyer_mass"] else None)
    out = {
        "source": str(args.stage.resolve()), "source_sha256": sha256(args.stage),
        "definition": "implied net closing funding ratio max(0,Q-A)/Q; A=b+(1-psi)*old owner asset value; Q=price*new owner rung",
        "gross_mortgage_ltv_identified": False,
        "young_age_cells": [float(x) for x in ages[young]],
        "young_age_label": f"model state ages {args.young_min:g}–{args.young_max:g}",
        "groups": summaries,
        "all_age_buyers_negative_closing_cash": all_age_negative,
        "renter_to_owner_definition": "origin tenure renter, destination owner; lifetime first ownership is not observed",
        "owner_switcher_definition": "origin owner, different destination owner rung",
        "stayers_excluded": True,
        "new_parent_unable_to_buy": {
            "status": "requires_exact_first_birth_flow_and_rule_specific_eligibility",
            "reason": "g_beginning_distribution is post-fertility; the n=1,m=1 stock includes prior first births. Renting is not proof of credit ineligibility.",
        },
    }
    if phi_choice is not None:
        out["phi_choice_range"] = [float(np.min(phi_choice)), float(np.max(phi_choice))]
    if args.first_birth_flow is not None and args.pre_birth is not None:
        raise ValueError("Supply either exact first-birth flow or pre-birth distribution")
    birth = None
    if args.first_birth_flow is not None:
        with np.load(args.first_birth_flow, allow_pickle=False) as saved:
            birth = saved["first_birth_flow"]
    if args.pre_birth is not None:
        if (args.fecundity is None) == (args.fecundity_by_age is None):
            raise ValueError("--pre-birth requires exactly one authenticated fecundity input")
        if args.fecundity_by_age is not None:
            fec = np.array([float(x) for x in args.fecundity_by_age.split(",")])
        else:
            fec = np.full(g.shape[3], args.fecundity)
        if fec.shape != (g.shape[3],) or np.any((fec < 0) | (fec > 1)) or not np.isfinite(fec).all():
            raise ValueError("Fecundity vector must be one valid probability per age cell")
        with np.load(args.pre_birth, allow_pickle=False) as saved:
            key = "distribution.g_pre" if "distribution.g_pre" in saved else "g_pre"
            pre = saved[key]
        if pre.shape != g.shape:
            raise ValueError("Pre-birth distribution shape disagrees with solution")
        if float(np.sum(pre[..., 0, 1:])) > 1e-10:
            raise ValueError("Childless readiness states are occupied; this extractor needs a separate birth-flow array")
        birth = np.zeros_like(g)
        if fert_probs.shape != g.shape[:-2] + (g.shape[-2],):
            raise ValueError("Fertility probabilities have unexpected shape")
        birth[..., 1, 1] = pre[..., 0, 0] * fert_probs[..., 1] * fec[None, None, None, :, None]
    if birth is not None:
        if birth.shape != g.shape or np.any(birth < -1e-12) or np.any(birth - g > 1e-8):
            raise ValueError("First-birth flow must be post-birth pre-tenure mass contained in g")
        if args.hard_phi is not None:
            # Strict closing cash test for an origin renter. This measures
            # down-payment exclusion only; utility and ending-estate constraints
            # can make more owner choices unavailable.
            first = birth[:, 0]
            min_q = price * rungs[0]
            blocked = b[:, None] < (1 - args.hard_phi) * min_q[None, :]
            blocked = blocked[:, :, None, None, None, None]
            denom = float(np.sum(first))
            out["new_parent_unable_to_buy"] = {
                "status": "strict_downpayment_ineligibility_only",
                "population": "first-birth flow from origin renters, all ages",
                "flow_mass": denom,
                "share_ineligible_for_every_owner_rung": float(np.sum(first * blocked) / denom) if denom else None,
                "caveat": "Only the closing cash screen is evaluated. Other feasibility and preference channels are excluded.",
            }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
