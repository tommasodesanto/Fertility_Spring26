#!/usr/bin/env python3
"""Audit active tenure value gaps at occupied wealth-grid nodes.

This is a read-only diagnostic for saved one-market intergenerational solution
caches. It reconstructs the conditional tenure-branch values implied by the
stored policies and next-period value function, then applies the same transaction
wealth and feasibility logic used by the tenure-choice kernel. The goal is to
make the aggregate consumption kinks inspectable without changing the solver.
"""

from __future__ import annotations

import argparse
import csv
import math
import pickle
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.solver import (  # noqa: E402
    apply_child_aging,
    bequest_utility_vec,
    get_completed_fertility,
    income_at_state,
    income_transition_values,
    precompute_shared,
    pti_adjusted_downpayment,
)


NEG_INF = -1.0e10
DEFAULT_WEALTHS = "0,0.146514,4.2487"
DEFAULT_CACHE = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "best_de_g044_i022_packet/solution_cache.pkl"
)


def main() -> None:
    args = parse_args()
    cache_path = args.cache.resolve()
    outdir = args.outdir.resolve() if args.outdir else cache_path.parent / "active_tenure_value_audit"
    outdir.mkdir(parents=True, exist_ok=True)

    payload = load_cache(cache_path)
    base = payload["baseline"]
    sol = base["sol"]
    P = base["P"]
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    price = owner_asset_price_vector(sol, base)
    rent = rent_vector(sol, P, price)
    shared = precompute_shared(P, b_grid)
    hcost, heq, dp_arr, bmo = transaction_arrays(P, price, b_grid, shared)

    wealth_targets = parse_float_list(args.wealths)
    selected_states = select_states(
        sol,
        P,
        wealth_targets=wealth_targets,
        top_per_wealth=int(args.top_per_wealth),
        extra_top_mass=int(args.extra_top_mass),
    )
    if not selected_states:
        raise RuntimeError("No positive-mass states selected")

    cache_vd: dict[tuple[int, int, int], np.ndarray] = {}
    state_rows: list[dict[str, Any]] = []
    option_rows: list[dict[str, Any]] = []
    for state_id, state in enumerate(selected_states, start=1):
        bb, to, i, j, zz, nn, cs = state["index"]
        key = (int(i), int(j), int(zz))
        if key not in cache_vd:
            cache_vd[key] = conditional_values_for_market_age_income(sol, P, int(i), int(j), int(zz))
        vd = cache_vd[key]
        option_data = evaluate_options(
            sol=sol,
            P=P,
            b_grid=b_grid,
            price=price,
            shared=shared,
            hcost=hcost,
            heq=heq,
            dp_arr=dp_arr,
            bmo=bmo,
            vd=vd,
            state=state,
        )
        state_summary = summarize_state(state_id, state, option_data, sol, P)
        state_rows.append(state_summary)
        for row in option_data:
            row["state_id"] = state_id
            for key_state, val_state in state_summary.items():
                if key_state in {
                    "state_id",
                    "best_option",
                    "second_option",
                    "gap_best_minus_second",
                    "stored_owner_probability",
                    "computed_owner_probability",
                    "max_abs_probability_diff",
                }:
                    continue
                row[key_state] = val_state
            option_rows.append(row)

    write_csv(outdir / "state_summary.csv", state_rows)
    write_csv(outdir / "active_tenure_value_gaps.csv", option_rows)
    plot_value_gaps(option_rows, state_rows, outdir / "selected_state_value_gaps.png")
    write_readme(
        outdir,
        cache_path=cache_path,
        payload=payload,
        P=P,
        price=price,
        rent=rent,
        state_rows=state_rows,
        option_rows=option_rows,
    )
    print(outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="Path to solution_cache.pkl.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output folder.")
    parser.add_argument(
        "--wealths",
        default=DEFAULT_WEALTHS,
        help="Comma-separated current liquid wealth nodes to inspect; nearest grid points are used.",
    )
    parser.add_argument(
        "--top-per-wealth",
        type=int,
        default=3,
        help="For each selected wealth node, inspect this many highest-mass state cells.",
    )
    parser.add_argument(
        "--extra-top-mass",
        type=int,
        default=3,
        help="Also inspect this many highest-mass state cells overall.",
    )
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} is not an intergen solution cache")
    for key in ("sol", "P"):
        if key not in payload["baseline"]:
            raise ValueError(f"{path} baseline is missing {key!r}")
    return payload


def parse_float_list(text: str) -> list[float]:
    out: list[float] = []
    for item in str(text).split(","):
        item = item.strip()
        if item:
            out.append(float(item))
    return out


def maybe_vector_value(values: np.ndarray, idx: int) -> float:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        return math.nan
    return float(arr[min(max(int(idx), 0), arr.size - 1)])


def interp_scalar(grid: np.ndarray, values: np.ndarray, query: float) -> float:
    x = np.asarray(grid, dtype=float).reshape(-1)
    y = np.asarray(values, dtype=float).reshape(-1)
    if x.size < 2 or y.size != x.size or not math.isfinite(query):
        return math.nan
    q = float(np.clip(query, x[0], x[-1]))
    return float(np.interp(q, x, y))


def crra_utility(consumption: np.ndarray | float, services: np.ndarray | float, alpha: float, oms: float) -> np.ndarray:
    c = np.maximum(np.asarray(consumption, dtype=float), 1.0e-10)
    h = np.maximum(np.asarray(services, dtype=float), 1.0e-10)
    composite = (c**alpha) * (h ** (1.0 - alpha))
    if abs(oms) < 1.0e-8:
        return np.log(composite)
    return composite**oms / oms


def owner_asset_price_vector(sol: Any, base: dict[str, Any]) -> np.ndarray:
    fallback = np.asarray(base.get("p_eq", getattr(sol, "p_eq", [math.nan])), dtype=float).reshape(-1)
    values = getattr(sol, "owner_asset_price", getattr(sol, "p_eq", fallback))
    arr = np.asarray(values, dtype=float).reshape(-1)
    return arr if arr.size else fallback


def rent_vector(sol: Any, P: Any, price: np.ndarray) -> np.ndarray:
    fallback = float(getattr(P, "user_cost_rate", math.nan)) * np.asarray(price, dtype=float).reshape(-1)
    values = getattr(sol, "owner_user_cost", fallback)
    arr = np.asarray(values, dtype=float).reshape(-1)
    return arr if arr.size else np.asarray(fallback, dtype=float).reshape(-1)


def tenure_label(P: Any, tenure_index: int) -> str:
    if int(tenure_index) <= 0:
        return "rent"
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    if int(tenure_index) - 1 < h_own.size:
        return f"own_H{h_own[int(tenure_index) - 1]:g}"
    return f"own_{int(tenure_index)}"


def liquidated_housing_value(P: Any, price: np.ndarray, tenure_index: int, market_index: int = 0) -> float:
    if int(tenure_index) <= 0:
        return 0.0
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float).reshape(-1)
    if int(tenure_index) - 1 >= h_own.size:
        return math.nan
    p = maybe_vector_value(price, market_index)
    if not math.isfinite(p):
        return math.nan
    return float((1.0 - float(P.psi)) * p * h_own[int(tenure_index) - 1])


def transaction_arrays(
    P: Any,
    price: np.ndarray,
    b_grid: np.ndarray,
    shared: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    I = int(P.I)
    nt = 1 + int(P.n_house)
    npar = int(P.n_parity)
    ncs = int(P.n_child_states)
    hcost = np.zeros((I, nt))
    heq = np.zeros((I, nt))
    dp_arr = np.zeros((I, nt, npar, ncs))
    bmo = np.zeros((I, nt, npar, ncs))
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    for i in range(I):
        p = maybe_vector_value(price, i)
        for ten in range(1, nt):
            house_cost = p * float(h_own[ten - 1])
            hcost[i, ten] = house_cost
            heq[i, ten] = (1.0 - float(P.psi)) * house_cost
            for nn in range(npar):
                for cs in range(ncs):
                    phi_ncs = float(shared.phi_choice[i, ten, nn, cs])
                    dp_arr[i, ten, nn, cs] = (1.0 - phi_ncs) * house_cost
                    bmo[i, ten, nn, cs] = -phi_ncs * house_cost
    if bool(getattr(P, "use_pti_constraint", False)):
        z_grid, _, _ = income_transition_values(P)
        income = np.array([income_at_state(P, i, 0, float(z_grid[0])) for i in range(I)], dtype=float)
        dp_arr = pti_adjusted_downpayment(dp_arr, hcost, income, P, b_grid)
    return hcost, heq, dp_arr, bmo


def conditional_values_for_market_age_income(sol: Any, P: Any, i: int, j: int, zz: int) -> np.ndarray:
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    V = np.asarray(sol.V, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    price = owner_asset_price_vector(sol, {"p_eq": getattr(sol, "p_eq", np.array([math.nan]))})
    rent = rent_vector(sol, P, price)
    shared = precompute_shared(P, b_grid)
    z_grid, _, Pi_z = income_transition_values(P)

    Nb = b_grid.size
    nt = 1 + int(P.n_house)
    npar = int(P.n_parity)
    ncs = int(P.n_child_states)
    beta = float(P.beta)
    Rg = float(P.R_gross)
    alpha = float(P.alpha_cons)
    oms = 1.0 - float(P.sigma)
    yj = income_at_state(P, i, j, float(z_grid[zz]))
    Rv = Rg * b_grid + float(yj)
    ri = maybe_vector_value(rent, i)

    if j == int(P.J) - 1:
        Vnr_i = np.zeros((Nb, nt, 1, npar, ncs))
        for ten in range(nt):
            hv = maybe_vector_value(price, i) * float(P.H_own[ten - 1]) if ten > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vnr_i[:, ten, 0, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)
    else:
        Vnr_i = np.zeros((Nb, nt, 1, npar, ncs))
        for znext in range(len(z_grid)):
            Vnr_i[:, :, 0, :, :] += float(Pi_z[zz, znext]) * V[:, :, i, j + 1, znext, :, :]
    Vc_i = apply_child_aging(Vnr_i, P, Nb, nt, 1, npar, ncs)[:, :, 0, :, :]

    out = np.full((Nb, nt, npar, ncs), NEG_INF, dtype=float)
    cb = np.asarray(shared.c_bar, dtype=float)
    hb = np.asarray(shared.h_bar, dtype=float)
    psi_child = np.asarray(shared.psi_v, dtype=float)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    owner_service_premium = max(float(getattr(P, "chi", 1.0)), 1.0e-8)
    owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
    owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
    owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))
    p = maybe_vector_value(price, i)

    for nn in range(npar):
        for cs in range(ncs):
            cb_c = float(cb[nn, cs])
            hb_c = float(hb[nn, cs])
            psi_c = float(psi_child[nn, cs])
            bp = bp_pol[:, 0, i, j, zz, nn, cs]
            surplus = Rv - cb_c - ri * hb_c - bp
            cap = ri * (float(P.hR_max) - hb_c) / max(1.0 - alpha, 1.0e-12)
            unc = surplus <= cap
            c_t = np.where(unc, alpha * np.maximum(surplus, 1.0e-10), Rv - cb_c - ri * float(P.hR_max) - bp)
            h_t = np.where(unc, (1.0 - alpha) / ri * np.maximum(surplus, 1.0e-10), float(P.hR_max) - hb_c)
            continuation = np.array([interp_scalar(b_grid, Vc_i[:, 0, nn, cs], x) for x in bp], dtype=float)
            val = crra_utility(c_t, h_t, alpha, oms) + psi_c + beta * continuation
            val[surplus <= 1.0e-10] = NEG_INF
            out[:, 0, nn, cs] = val

            for ten in range(1, nt):
                hs = float(h_own[ten - 1])
                extra_size_cost = owner_size_cost * p * max(hs - owner_size_cost_ref, 0.0) ** owner_size_cost_power
                owner_cost = (float(P.delta) + float(P.tau_H)) * p * hs + extra_size_cost
                bp_o = bp_pol[:, ten, i, j, zz, nn, cs]
                c_t_o = Rv - owner_cost - cb_c - bp_o
                h_t_o = owner_service_premium * max(hs - owner_h_bar_scale * hb_c, 1.0e-10)
                continuation_o = np.array([interp_scalar(b_grid, Vc_i[:, ten, nn, cs], x) for x in bp_o], dtype=float)
                val_o = crra_utility(c_t_o, h_t_o, alpha, oms) + psi_c + beta * continuation_o
                val_o[c_t_o <= 1.0e-10] = NEG_INF
                out[:, ten, nn, cs] = val_o
    return out


def branch_info(
    *,
    b: float,
    origin_tenure: int,
    target_tenure: int,
    market: int,
    parity: int,
    child_state: int,
    P: Any,
    hcost: np.ndarray,
    heq: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    shared: Any,
) -> dict[str, Any]:
    to = int(origin_tenure)
    tn = int(target_tenure)
    if tn == to:
        return {
            "branch_liquid_wealth": float(b),
            "required_current_liquid": 0.0,
            "downpayment": 0.0,
            "borrowing_floor": float(bmo[market, tn, parity, child_state]) if tn > 0 else math.nan,
            "feasible": True,
            "feasibility_reason": "same_tenure",
        }
    sale = float(heq[market, to]) if to > 0 else 0.0
    if tn <= 0:
        return {
            "branch_liquid_wealth": float(max(b + sale, 0.0)),
            "required_current_liquid": 0.0,
            "downpayment": 0.0,
            "borrowing_floor": math.nan,
            "feasible": True,
            "feasibility_reason": "sell_or_rent",
        }

    house_cost = float(hcost[market, tn])
    dpn = float(dp_arr[market, tn, parity, child_state])
    bmn = float(bmo[market, tn, parity, child_state])
    if to <= 0:
        branch = float(b - house_cost)
        grant = float(shared.birth_entry_grant[market, tn, parity, child_state])
        if bool(shared.birth_dp[parity, child_state, to, tn]):
            return {
                "branch_liquid_wealth": max(branch, bmn),
                "required_current_liquid": 0.0,
                "downpayment": dpn,
                "borrowing_floor": bmn,
                "feasible": True,
                "feasibility_reason": "birth_downpayment_waiver",
            }
        if grant > 0.0:
            branch_g = branch + grant
            feasible = (b + grant) >= dpn and branch_g >= bmn
            return {
                "branch_liquid_wealth": branch_g,
                "required_current_liquid": max(dpn - grant, 0.0),
                "downpayment": dpn,
                "borrowing_floor": bmn,
                "feasible": bool(feasible),
                "feasibility_reason": "grant" if feasible else "grant_constraint",
            }
        feasible = b >= dpn and branch >= bmn
        return {
            "branch_liquid_wealth": branch,
            "required_current_liquid": dpn,
            "downpayment": dpn,
            "borrowing_floor": bmn,
            "feasible": bool(feasible),
            "feasibility_reason": "purchase" if feasible else "downpayment_or_floor",
        }

    branch = float(b + sale - house_cost)
    required = dpn - sale
    feasible = b >= required and branch >= bmn
    return {
        "branch_liquid_wealth": branch,
        "required_current_liquid": required,
        "downpayment": dpn,
        "borrowing_floor": bmn,
        "feasible": bool(feasible),
        "feasibility_reason": "rebuy" if feasible else "rebuy_constraint",
    }


def evaluate_options(
    *,
    sol: Any,
    P: Any,
    b_grid: np.ndarray,
    price: np.ndarray,
    shared: Any,
    hcost: np.ndarray,
    heq: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    vd: np.ndarray,
    state: dict[str, Any],
) -> list[dict[str, Any]]:
    bb, to, i, j, zz, nn, cs = state["index"]
    b = float(b_grid[bb])
    nt = vd.shape[1]
    raw_vals: list[float] = []
    rows: list[dict[str, Any]] = []
    stored_probs = target_probabilities(sol, bb, to, i, j, zz, nn, cs)
    for tn in range(nt):
        info = branch_info(
            b=b,
            origin_tenure=to,
            target_tenure=tn,
            market=i,
            parity=nn,
            child_state=cs,
            P=P,
            hcost=hcost,
            heq=heq,
            dp_arr=dp_arr,
            bmo=bmo,
            shared=shared,
        )
        branch_b = float(info["branch_liquid_wealth"])
        if info["feasible"]:
            if tn == to:
                value = float(vd[bb, tn, nn, cs])
            else:
                value = interp_scalar(b_grid, vd[:, tn, nn, cs], branch_b)
        else:
            value = NEG_INF
        raw_vals.append(value)
        c_val = interp_scalar(b_grid, np.asarray(sol.c_pol)[:, tn, i, j, zz, nn, cs], branch_b)
        bp_val = interp_scalar(b_grid, np.asarray(sol.bp_pol)[:, tn, i, j, zz, nn, cs], branch_b)
        if tn <= 0:
            h_val = interp_scalar(b_grid, np.asarray(sol.hR_pol)[:, 0, i, j, zz, nn, cs], branch_b)
        else:
            h_val = float(np.asarray(P.H_own, dtype=float).reshape(-1)[tn - 1])
        rows.append(
            {
                "target_tenure": tn,
                "target_label": tenure_label(P, tn),
                "feasible": bool(info["feasible"]),
                "feasibility_reason": info["feasibility_reason"],
                "value": value,
                "branch_liquid_wealth": branch_b,
                "post_tenure_total_wealth": branch_b + liquidated_housing_value(P, price, tn, i),
                "required_current_liquid": float(info["required_current_liquid"]),
                "downpayment": float(info["downpayment"]),
                "borrowing_floor": float(info["borrowing_floor"]) if math.isfinite(float(info["borrowing_floor"])) else math.nan,
                "stored_probability": float(stored_probs[tn]) if tn < stored_probs.size else math.nan,
                "computed_probability": math.nan,
                "branch_consumption": c_val,
                "branch_housing_services": h_val,
                "branch_next_liquid_wealth": bp_val,
            }
        )
    probs = logit_probs(np.asarray(raw_vals, dtype=float), float(getattr(P, "tenure_choice_kappa", 0.0)))
    best = float(np.max(raw_vals))
    for row, pr in zip(rows, probs):
        row["computed_probability"] = float(pr)
        row["value_minus_best"] = float(row["value"] - best) if row["value"] > NEG_INF / 2 else math.nan
    return rows


def logit_probs(values: np.ndarray, kappa: float) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    out = np.zeros(vals.size, dtype=float)
    finite = vals > NEG_INF / 2
    if not np.any(finite):
        return out
    best = float(np.max(vals[finite]))
    if kappa <= 1.0e-12:
        winners = np.isclose(vals, best)
        out[winners] = 1.0 / max(float(np.sum(winners)), 1.0)
        return out
    ev = np.zeros(vals.size, dtype=float)
    ev[finite] = np.exp((vals[finite] - best) / kappa)
    denom = float(np.sum(ev))
    if denom > 0.0:
        out = ev / denom
    return out


def target_probabilities(sol: Any, bb: int, to: int, i: int, j: int, zz: int, nn: int, cs: int) -> np.ndarray:
    nt = np.asarray(sol.g).shape[1]
    tp = getattr(sol, "tenure_probs", None)
    if tp is None:
        out = np.zeros(nt, dtype=float)
        out[int(sol.tenure_choice[bb, to, i, j, zz, nn, cs])] = 1.0
        return out
    return np.asarray(tp[bb, to, i, j, zz, nn, cs, :], dtype=float).reshape(-1)


def select_states(
    sol: Any,
    P: Any,
    *,
    wealth_targets: Iterable[float],
    top_per_wealth: int,
    extra_top_mass: int,
) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    selected: dict[tuple[int, int, int, int, int, int, int], dict[str, Any]] = {}

    def add_indices(indices: list[tuple[int, int, int, int, int, int, int]], reason: str) -> None:
        for idx in indices:
            mass = float(g[idx])
            if mass <= 0.0:
                continue
            bb, to, i, j, zz, nn, cs = idx
            key = tuple(int(x) for x in idx)
            selected[key] = {
                "index": key,
                "selection_reason": reason,
                "state_mass": mass,
                "current_liquid_wealth": float(b_grid[bb]),
                "current_total_wealth": float(b_grid[bb] + liquidated_housing_value(P, owner_asset_price_vector(sol, {"p_eq": getattr(sol, "p_eq", [math.nan])}), to, i)),
                "current_tenure": int(to),
                "current_label": tenure_label(P, to),
                "market": int(i),
                "age_index": int(j),
                "age": float(age_grid[j]),
                "income_state": int(zz),
                "parity": int(nn),
                "child_state": int(cs),
            }

    for target in wealth_targets:
        bb = int(np.argmin(np.abs(b_grid - float(target))))
        slice_mass = g[bb, :, :, :, :, :, :]
        flat_order = np.argsort(slice_mass.reshape(-1))[::-1]
        indices: list[tuple[int, int, int, int, int, int, int]] = []
        for flat in flat_order[: max(top_per_wealth * 4, top_per_wealth)]:
            if float(slice_mass.reshape(-1)[flat]) <= 0.0:
                continue
            rest = np.unravel_index(int(flat), slice_mass.shape)
            indices.append((bb, *[int(x) for x in rest]))
            if len(indices) >= top_per_wealth:
                break
        add_indices(indices, f"nearest_b_to_{target:g}")

    if extra_top_mass > 0:
        flat_order = np.argsort(g.reshape(-1))[::-1]
        indices = []
        for flat in flat_order:
            if float(g.reshape(-1)[flat]) <= 0.0:
                continue
            indices.append(tuple(int(x) for x in np.unravel_index(int(flat), g.shape)))
            if len(indices) >= extra_top_mass:
                break
        add_indices(indices, "top_mass_overall")

    out = list(selected.values())
    out.sort(key=lambda r: (-float(r["state_mass"]), float(r["current_liquid_wealth"]), int(r["current_tenure"])))
    return out


def summarize_state(state_id: int, state: dict[str, Any], option_rows: list[dict[str, Any]], sol: Any, P: Any) -> dict[str, Any]:
    values = np.asarray([r["value"] for r in option_rows], dtype=float)
    probs_stored = np.asarray([r["stored_probability"] for r in option_rows], dtype=float)
    probs_computed = np.asarray([r["computed_probability"] for r in option_rows], dtype=float)
    finite = values > NEG_INF / 2
    if np.any(finite):
        order = np.argsort(values)[::-1]
        best_idx = int(order[0])
        second_idx = int(order[1]) if order.size > 1 else best_idx
        gap = float(values[best_idx] - values[second_idx]) if second_idx != best_idx and values[second_idx] > NEG_INF / 2 else math.inf
    else:
        best_idx = 0
        second_idx = 0
        gap = math.nan
    bb, to, i, j, zz, nn, cs = state["index"]
    stored_choice = int(sol.tenure_choice[bb, to, i, j, zz, nn, cs])
    row = {
        "state_id": state_id,
        **state,
        "stored_choice": stored_choice,
        "stored_choice_label": tenure_label(P, stored_choice),
        "best_option": int(best_idx),
        "best_option_label": tenure_label(P, best_idx),
        "second_option": int(second_idx),
        "second_option_label": tenure_label(P, second_idx),
        "gap_best_minus_second": gap,
        "stored_owner_probability": float(np.sum(probs_stored[1:])),
        "computed_owner_probability": float(np.sum(probs_computed[1:])),
        "max_abs_probability_diff": float(np.nanmax(np.abs(probs_stored - probs_computed))),
    }
    row.pop("index", None)
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_value_gaps(option_rows: list[dict[str, Any]], state_rows: list[dict[str, Any]], path: Path) -> None:
    state_ids = [int(r["state_id"]) for r in state_rows[:8]]
    if not state_ids:
        return
    fig, axes = plt.subplots(len(state_ids), 1, figsize=(9, max(2.2, 1.7 * len(state_ids))), squeeze=False)
    for ax, state_id in zip(axes[:, 0], state_ids):
        rows = [r for r in option_rows if int(r["state_id"]) == state_id]
        labels = [r["target_label"] for r in rows]
        vals = [max(float(r.get("value_minus_best", math.nan)), -0.25) if math.isfinite(float(r.get("value_minus_best", math.nan))) else -0.25 for r in rows]
        colors = ["0.1" if abs(v) < 1.0e-10 else "0.55" for v in vals]
        ax.bar(labels, vals, color=colors)
        st = next(r for r in state_rows if int(r["state_id"]) == state_id)
        ax.set_title(
            f"state {state_id}: b={float(st['current_liquid_wealth']):.3g}, "
            f"{st['current_label']}, age {float(st['age']):.0f}, mass={float(st['state_mass']):.3g}",
            fontsize=9,
        )
        ax.set_ylabel("V - max")
        ax.set_ylim(min(-0.25, min(vals) * 1.1), 0.02)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_readme(
    outdir: Path,
    *,
    cache_path: Path,
    payload: dict[str, Any],
    P: Any,
    price: np.ndarray,
    rent: np.ndarray,
    state_rows: list[dict[str, Any]],
    option_rows: list[dict[str, Any]],
) -> None:
    max_prob_diff = max(float(r["max_abs_probability_diff"]) for r in state_rows)
    min_gap = min(float(r["gap_best_minus_second"]) for r in state_rows if math.isfinite(float(r["gap_best_minus_second"])))
    with (outdir / "README.md").open("w", encoding="utf-8") as fh:
        fh.write("# Active Tenure Value-Gap Audit\n\n")
        fh.write(f"- Source cache: `{cache_path}`\n")
        fh.write(f"- Target set in cache: `{payload.get('target_set', '')}`\n")
        fh.write(f"- Owner asset price: `{np.asarray(price, dtype=float).reshape(-1).tolist()}`\n")
        fh.write(f"- Flow rent/user cost: `{np.asarray(rent, dtype=float).reshape(-1).tolist()}`\n")
        fh.write(f"- Tenure logit kappa: `{float(getattr(P, 'tenure_choice_kappa', 0.0)):.8g}`\n")
        fh.write("- Method: reconstruct conditional branch values from saved policies and next-period value function, then apply the tenure-kernel transaction/feasibility accounting.\n")
        fh.write(f"- Selected states: `{len(state_rows)}`; option rows: `{len(option_rows)}`\n")
        fh.write(f"- Largest stored-vs-computed tenure-probability difference: `{max_prob_diff:.6g}`\n")
        fh.write(f"- Smallest finite best-minus-second value gap among selected states: `{min_gap:.6g}`\n\n")
        fh.write("## State Summary\n\n")
        fh.write("| state | b | current | age | mass | stored choice | best | second | gap | Pr(owner) |\n")
        fh.write("|---:|---:|---|---:|---:|---|---|---|---:|---:|\n")
        for row in state_rows:
            fh.write(
                f"| {int(row['state_id'])} | {float(row['current_liquid_wealth']):.6g} | "
                f"{row['current_label']} | {float(row['age']):.0f} | {float(row['state_mass']):.6g} | "
                f"{row['stored_choice_label']} | {row['best_option_label']} | {row['second_option_label']} | "
                f"{float(row['gap_best_minus_second']):.6g} | {float(row['stored_owner_probability']):.4f} |\n"
            )
        fh.write("\n## Files\n\n")
        fh.write("- `state_summary.csv`\n")
        fh.write("- `active_tenure_value_gaps.csv`\n")
        fh.write("- `selected_state_value_gaps.png`\n")


if __name__ == "__main__":
    main()
