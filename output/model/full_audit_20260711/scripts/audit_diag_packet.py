#!/usr/bin/env python3
"""Audit-only diagnostic packet for a July-10 combined-spec candidate.

Rebuilds the exact combined-spec environment (same merge order as
tmp/overnight_combined_20260710/run_wave2_lateral.py and
output/model/full_audit_20260711/scripts/reproduce_candidate.py):
{base_overrides, **extra(production profile + fixed spec), **income(rouwenhorst),
**theta} -> intergen_housing_fertility.solver.run_model_cp_dt, then:

1. writes the package STANDARD diagnostic packet via
   intergen_housing_fertility.diagnostics.write_diagnostics (no new graph set);
2. runs numerical pathology checks (value monotonicity, saving/consumption,
   fertility probabilities, tenure by age, boundary states, market clearing,
   stationary mass accounting) and writes them to numeric_checks.json.

Writes ONLY into output/model/full_audit_20260711/. Does not touch production.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.diagnostics import write_diagnostics
from intergen_housing_fertility.solver import run_model_cp_dt
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0
VALID_V = -1e9  # infeasible nodes carry -1e10 (kernels.py:606 NEG_INF = -1e10)


def build_overrides(nb: int, max_iter_eq: int, theta: dict) -> tuple[dict, dict, dict]:
    targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    income = lp.income_process_overrides(
        5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
    )
    extra = production_profile_overrides()
    extra.update(
        {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": np.array([1.75]),
            "normalize_bequest_utility": True,
        }
    )
    extra["max_iter_eq"] = int(max_iter_eq)
    overrides = {
        **lp.base_overrides(J=PRODUCTION_J, Nb=nb, n_house=5, max_iter_eq=int(max_iter_eq)),
        **extra,
        **income,
        **theta,
    }
    return overrides, targets, weights


def jnum(x):
    if isinstance(x, dict):
        return {str(k): jnum(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jnum(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nb", type=int, default=120)
    ap.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    args = ap.parse_args()

    theta = dict(json.loads(args.record.read_text())["theta"])
    overrides, targets, weights = build_overrides(args.nb, args.max_iter_eq, theta)

    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    wall = time.perf_counter() - t0

    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- standard package diagnostic packet (existing graph set) ----
    write_diagnostics(sol, P, outdir)

    # ---- numeric pathology checks ----
    g = np.asarray(sol.g)  # (Nb, nt, I, J, Nz, npar, ncs)
    V = np.asarray(sol.V)
    c = np.asarray(sol.c_pol)
    bp = np.asarray(sol.bp_pol)
    hR = np.asarray(sol.hR_pol)
    fp = np.asarray(sol.fert_probs)
    tc = np.asarray(sol.tenure_choice)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    Nb, nt, I, J, Nz, npar, ncs = g.shape
    total_mass = float(np.sum(g))
    checks: dict = {
        "wall_seconds": wall,
        "p_eq": p_eq,
        "shape": {"Nb": Nb, "nt": nt, "I": I, "J": J, "Nz": Nz, "n_parity": npar, "n_child_states": ncs},
        "b_grid": {"b_min": float(b_grid[0]), "b_max": float(b_grid[-1])},
    }

    # 1. value monotonicity in wealth on valid nodes
    valid = V > VALID_V
    viol_count = 0
    viol_populated_mass = 0.0
    worst = []
    by_age = np.zeros(J, dtype=int)
    for ten in range(nt):
        for i in range(I):
            for j in range(J):
                for zz in range(Nz):
                    for nn in range(npar):
                        for cs in range(ncs):
                            v = V[:, ten, i, j, zz, nn, cs]
                            m = valid[:, ten, i, j, zz, nn, cs]
                            idx = np.flatnonzero(m)
                            if idx.size < 2:
                                continue
                            dv = np.diff(v[idx])
                            bad = np.flatnonzero(dv < -1e-9)
                            if bad.size:
                                viol_count += int(bad.size)
                                by_age[j] += int(bad.size)
                                mass_here = float(np.sum(g[idx[bad], ten, i, j, zz, nn, cs]))
                                viol_populated_mass += mass_here
                                worst.append(
                                    {
                                        "state": {"ten": ten, "i": i, "j": j, "z": zz, "n": nn, "cs": cs},
                                        "n_drops": int(bad.size),
                                        "max_drop": float(np.min(dv[bad])),
                                        "mass_at_drop_nodes": mass_here,
                                    }
                                )
    worst.sort(key=lambda d: d["max_drop"])
    checks["value_monotonicity"] = {
        "pairwise_decreases_valid_nodes": viol_count,
        "violations_by_age_index": by_age.tolist(),
        "mass_at_violation_nodes": viol_populated_mass,
        "mass_share": viol_populated_mass / max(total_mass, 1e-14),
        "worst_10": worst[:10],
        "n_invalid_nodes": int(np.sum(~valid)),
        "invalid_share_of_nodes": float(np.mean(~valid)),
        "populated_invalid_mass": float(np.sum(g[~valid])),
    }

    # 2. saving / consumption
    pop = g > 1e-14
    c_pop = c[pop]
    checks["consumption"] = {
        "min_c_on_populated_states": float(np.min(c_pop)),
        "n_populated_states_c_le_0": int(np.sum(c_pop <= 0.0)),
        "mass_c_le_0": float(np.sum(g[pop & (c <= 0.0)])) if np.any(pop & (c <= 0.0)) else 0.0,
        "min_c_anywhere_valid": float(np.min(c[valid])) if np.any(valid) else np.nan,
    }
    bp_pop = bp[pop]
    checks["savings_policy"] = {
        "min_bp_populated": float(np.min(bp_pop)),
        "max_bp_populated": float(np.max(bp_pop)),
        "share_populated_bp_at_bmin": float(np.mean(np.isclose(bp_pop, b_grid[0]))),
        "share_populated_bp_at_bmax": float(np.mean(np.isclose(bp_pop, b_grid[-1]))),
        "mass_bp_at_bmax": float(np.sum(g[pop & np.isclose(bp, b_grid[-1])])),
        "mass_bp_at_bmin": float(np.sum(g[pop & np.isclose(bp, b_grid[0])])),
    }
    mass_bnode = np.sum(g, axis=(1, 2, 3, 4, 5, 6))
    checks["wealth_distribution"] = {
        "mass_at_lowest_b_node": float(mass_bnode[0]),
        "mass_at_highest_b_node": float(mass_bnode[-1]),
        "share_at_lowest_b_node": float(mass_bnode[0] / max(total_mass, 1e-14)),
        "share_at_highest_b_node": float(mass_bnode[-1] / max(total_mass, 1e-14)),
        "share_top_3_b_nodes": float(np.sum(mass_bnode[-3:]) / max(total_mass, 1e-14)),
        "share_negative_wealth": float(np.sum(mass_bnode[b_grid < 0]) / max(total_mass, 1e-14)),
    }

    # 3. fertility probabilities
    fp_min = float(np.min(fp))
    fp_max = float(np.max(fp))
    rowsum = np.sum(fp, axis=-1)
    fw = range(P.A_f_start - 1, P.A_f_end + 1)
    p_birth = 1.0 - fp[..., 0]  # prob of >=1 birth for childless
    pzero_cells = 0
    pzero_mass = 0.0
    fert_cells = 0
    for j in fw:
        pj = p_birth[:, :, :, j, :]
        gj = g[:, :, :, j, :, 0, 0]
        fert_cells += pj.size
        z = pj <= 1e-12
        pzero_cells += int(np.sum(z))
        pzero_mass += float(np.sum(gj[z]))
    childless_fertile_mass = float(sum(np.sum(g[:, :, :, j, :, 0, 0]) for j in fw))
    checks["fertility_probs"] = {
        "min": fp_min,
        "max": fp_max,
        "rowsum_min": float(np.min(rowsum)),
        "rowsum_max": float(np.max(rowsum)),
        "fertile_window_age_indices": [int(j) for j in fw],
        "p_birth_zero_cells": pzero_cells,
        "p_birth_cells_total": fert_cells,
        "p_birth_zero_cell_share": pzero_cells / max(fert_cells, 1),
        "childless_mass_at_p_birth_zero": pzero_mass,
        "childless_fertile_mass": childless_fertile_mass,
        "childless_mass_share_at_p_birth_zero": pzero_mass / max(childless_fertile_mass, 1e-14),
    }

    # 4. tenure / ownership by age
    own_by_age = np.asarray(sol.own_by_age, dtype=float).reshape(-1)
    ages = P.age_start + np.arange(J) * P.da
    jump_j = int(np.argmax(np.diff(own_by_age))) if J > 1 else 0
    checks["ownership_by_age"] = {
        "ages": ages.tolist(),
        "own_by_age": own_by_age.tolist(),
        "largest_jump_between_ages": [float(ages[jump_j]), float(ages[jump_j + 1])],
        "largest_jump_size": float(np.max(np.diff(own_by_age))),
        "own_rate": float(sol.own_rate),
    }
    # tenure-state mass by age
    mass_ten_age = np.sum(g, axis=(0, 2, 4, 5, 6))  # (nt, J)
    checks["tenure_mass_by_age"] = {
        "renter_share_by_age": (mass_ten_age[0] / np.maximum(np.sum(mass_ten_age, axis=0), 1e-14)).tolist(),
        "top_rung_share_by_age": (mass_ten_age[-1] / np.maximum(np.sum(mass_ten_age, axis=0), 1e-14)).tolist(),
        "owner_mass_by_rung": np.sum(mass_ten_age[1:], axis=1).tolist(),
        "H_own": np.asarray(P.H_own, dtype=float).tolist(),
    }

    # 5. boundary: renter cap bunching, top rung usage
    hR_cap = float(P.hR_max)
    renter_mass = g[:, 0]  # (Nb, I, J, Nz, npar, ncs)
    hr_choice = hR[:, 0]
    at_cap = np.isclose(hr_choice, hR_cap, atol=1e-9)
    renter_total = float(np.sum(renter_mass))
    checks["renter_cap"] = {
        "hR_max": hR_cap,
        "renter_mass_total": renter_total,
        "renter_mass_at_cap": float(np.sum(renter_mass[at_cap])),
        "renter_share_at_cap": float(np.sum(renter_mass[at_cap]) / max(renter_total, 1e-14)),
    }
    owner_mass = float(np.sum(mass_ten_age[1:]))
    checks["top_owner_rung"] = {
        "top_rung_mass": float(np.sum(mass_ten_age[-1])),
        "top_rung_share_of_owners": float(np.sum(mass_ten_age[-1]) / max(owner_mass, 1e-14)),
        "bottom_rung_share_of_owners": float(np.sum(mass_ten_age[1]) / max(owner_mass, 1e-14)),
        "owner_rung_shares": (np.sum(mass_ten_age[1:], axis=1) / max(owner_mass, 1e-14)).tolist(),
    }

    # 6. market clearing
    demand = np.asarray(sol.rental_demand_by_market) + np.asarray(sol.owner_demand_by_market)
    supply = np.asarray(sol.housing_supply)
    checks["market"] = {
        "p_eq": np.asarray(p_eq).tolist(),
        "demand_by_market": demand.tolist(),
        "supply_by_market": supply.tolist(),
        "rel_residual_by_market": ((demand - supply) / np.maximum(supply, 1e-12)).tolist(),
        "best_max_abs_rel_excess": float(sol.best_max_abs_rel_excess),
        "tol_eq": float(getattr(P, "tol_eq", np.nan)),
        "aggregate_rental_demand": float(sol.aggregate_rental_demand),
        "aggregate_owner_demand": float(sol.aggregate_owner_demand),
        "aggregate_supply": float(sol.aggregate_housing_supply),
    }

    # 7. mass accounting by age
    age_mass = np.sum(g, axis=(0, 1, 2, 4, 5, 6))
    entry_total = float(np.sum(np.asarray(P.entry_by_loc, dtype=float)))
    checks["mass_accounting"] = {
        "total_mass": total_mass,
        "sol_total_mass_attr": float(getattr(sol, "total_mass", np.nan)),
        "age_mass": age_mass.tolist(),
        "entry_mass_total_per_cohort": entry_total,
        "age_mass_over_entry": (age_mass / max(entry_total, 1e-14)).tolist(),
        "max_abs_dev_from_entry": float(np.max(np.abs(age_mass - entry_total))),
        "normalize_population_mass": bool(getattr(P, "normalize_population_mass", True)),
        "mass_negative_nodes": int(np.sum(g < -1e-14)),
        "min_g": float(np.min(g)),
    }

    # cross-check moments and loss against the record
    moments = lp.extract_moments(sol, P)
    loss = lp.diagnostic_loss(moments, targets=targets, weights=weights)
    checks["loss_crosscheck"] = {"rank_loss": float(loss)}

    (outdir / "numeric_checks.json").write_text(json.dumps(jnum(checks), indent=2, sort_keys=True))
    print(json.dumps(jnum({k: checks[k] for k in ("loss_crosscheck", "market")}), indent=1))
    print(f"packet written to {outdir}, wall={wall:.1f}s")


if __name__ == "__main__":
    main()
