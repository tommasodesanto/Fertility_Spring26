#!/usr/bin/env python3
"""Extension probes: (1) current-bound candidate — entry-node feasibility,
zero-rowsum fertility rows by age, and the policies followed on populated
infeasible nodes; (2) light numeric pass + standard packet for the relaxed
(housing) candidate, solved in PE at its stored p_eq. Audit-only."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from audit_diag_packet import build_overrides, jnum, VALID_V  # noqa: E402

from intergen_housing_fertility.diagnostics import write_diagnostics  # noqa: E402
from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402

BASE = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711")


def solve_pe(record: Path, p_fixed: float):
    theta = dict(json.loads(record.read_text())["theta"])
    overrides, targets, weights = build_overrides(120, 10, theta)
    overrides["solve_mode"] = "pe"
    overrides["p_fixed"] = np.array([float(p_fixed)])
    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    return sol, P, time.perf_counter() - t0


def light_checks(sol, P) -> dict:
    g = np.asarray(sol.g)
    V = np.asarray(sol.V)
    fp = np.asarray(sol.fert_probs)
    hR = np.asarray(sol.hR_pol)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    Nb, nt, I, J, Nz, npar, ncs = g.shape
    total = float(np.sum(g))
    invalid = V <= VALID_V
    # monotonicity count (valid nodes)
    viol = 0
    for ten in range(nt):
        for i in range(I):
            for j in range(J):
                for zz in range(Nz):
                    for nn in range(npar):
                        for cs in range(ncs):
                            v = V[:, ten, i, j, zz, nn, cs]
                            m = ~invalid[:, ten, i, j, zz, nn, cs]
                            idx = np.flatnonzero(m)
                            if idx.size >= 2:
                                viol += int(np.sum(np.diff(v[idx]) < -1e-9))
    p_birth = 1.0 - fp[..., 0]
    pz_mass = 0.0
    c_mass = 0.0
    for j in range(P.A_f_start - 1, P.A_f_end):  # operative KFE window j..A_f_end-1
        gj = g[:, :, :, j, :, 0, 0]
        pz_mass += float(np.sum(gj[p_birth[:, :, :, j, :] <= 1e-12]))
        c_mass += float(np.sum(gj))
    renter_mass = g[:, 0]
    at_cap = np.isclose(hR[:, 0], float(P.hR_max), atol=1e-9)
    age_mass = np.sum(g, axis=(0, 1, 2, 4, 5, 6))
    demand = np.asarray(sol.rental_demand_by_market) + np.asarray(sol.owner_demand_by_market)
    supply = np.asarray(sol.housing_supply)
    return {
        "own_rate": float(sol.own_rate),
        "own_by_age": np.asarray(sol.own_by_age, dtype=float).tolist(),
        "value_monotonicity_pairwise_decreases": viol,
        "populated_invalid_mass": float(np.sum(g[invalid])),
        "populated_invalid_share": float(np.sum(g[invalid]) / total),
        "min_c_populated": float(np.min(np.asarray(sol.c_pol)[g > 1e-14])),
        "p_birth_zero_childless_mass_share_operative_window": pz_mass / max(c_mass, 1e-14),
        "renter_share_at_cap": float(np.sum(renter_mass[at_cap]) / max(float(np.sum(renter_mass)), 1e-14)),
        "share_at_lowest_b_node": float(np.sum(g[0]) / total),
        "share_at_highest_b_node": float(np.sum(g[-1]) / total),
        "age_mass_max_abs_dev": float(np.max(np.abs(age_mass - age_mass[0]))),
        "market_rel_residual": ((demand - supply) / np.maximum(supply, 1e-12)).tolist(),
        "total_mass": total,
    }


def main() -> None:
    # ---- Part 1: current-bound extras ----
    sol, P, wall = solve_pe(BASE / "repro/current_repro_exact.json", 0.6897624769313723)
    g = np.asarray(sol.g)
    V = np.asarray(sol.V)
    fp = np.asarray(sol.fert_probs)
    c = np.asarray(sol.c_pol)
    bp = np.asarray(sol.bp_pol)
    tc = np.asarray(sol.tenure_choice)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    Nb, nt, I, J, Nz, npar, ncs = g.shape
    invalid = V <= VALID_V
    out: dict = {"wall_seconds": wall}

    entry_idx = int(np.argmax(b_grid >= float(getattr(P, "b_entry_fixed", 0.0))))
    out["entry_node"] = {
        "b_entry_fixed": float(getattr(P, "b_entry_fixed", 0.0)),
        "entry_b_value": float(b_grid[entry_idx]),
        "V_at_entry_by_z_childless_renter_j0": V[entry_idx, 0, 0, 0, :, 0, 0].tolist(),
        "invalid_at_entry_by_z": invalid[entry_idx, 0, 0, 0, :, 0, 0].tolist(),
        "entry_mass_by_z": g[entry_idx, 0, 0, 0, :, 0, 0].tolist(),
    }

    rowsum = np.sum(fp, axis=-1)
    zr_by_age = []
    for j in range(J):
        gj = g[:, :, :, j, :, 0, 0]
        zr_by_age.append(float(np.sum(gj[rowsum[:, :, :, j, :] < 1e-12])))
    out["zero_rowsum_childless_mass_by_age"] = zr_by_age

    pi_mask = invalid & (g > 1e-14)
    out["invalid_populated_policies"] = {
        "n_nodes": int(np.sum(pi_mask)),
        "mass": float(np.sum(g[pi_mask])),
        "c_pol_min": float(np.min(c[pi_mask])) if np.any(pi_mask) else None,
        "c_pol_max": float(np.max(c[pi_mask])) if np.any(pi_mask) else None,
        "bp_min": float(np.min(bp[pi_mask])) if np.any(pi_mask) else None,
        "bp_max": float(np.max(bp[pi_mask])) if np.any(pi_mask) else None,
        "tenure_choice_values": np.unique(tc[pi_mask]).tolist() if np.any(pi_mask) else [],
    }
    # j0 breakdown of invalid populated mass by (z, parity)
    j0 = (g[:, :, :, 0, :, :, :] * invalid[:, :, :, 0, :, :, :])
    out["invalid_populated_j0"] = {
        "by_z": np.sum(j0, axis=(0, 1, 2, 4, 5)).tolist(),
        "by_parity": np.sum(j0, axis=(0, 1, 2, 3, 5)).tolist(),
        "total": float(np.sum(j0)),
    }
    (BASE / "diagnostics/current_bound/extra_checks.json").write_text(json.dumps(jnum(out), indent=2, sort_keys=True))
    print("current extras:", json.dumps(jnum(out["entry_node"]), indent=1))
    print("zero_rowsum_by_age:", [round(x, 5) for x in zr_by_age])
    print("invalid_populated_policies:", json.dumps(jnum(out["invalid_populated_policies"]), indent=1))
    print("invalid_populated_j0:", json.dumps(jnum(out["invalid_populated_j0"]), indent=1))

    # ---- Part 2: relaxed candidate light pass (PE at stored p_eq) ----
    rec = BASE / "repro/housing_repro_exact.json"
    p_eq = float(json.loads(rec.read_text())["p_eq"][0])
    sol2, P2, wall2 = solve_pe(rec, p_eq)
    outdir = BASE / "diagnostics/relaxed_light"
    outdir.mkdir(parents=True, exist_ok=True)
    write_diagnostics(sol2, P2, outdir)
    lc = light_checks(sol2, P2)
    lc["wall_seconds"] = wall2
    lc["note"] = "PE solve at stored p_eq=%r (light pass, not a GE re-solve)" % p_eq
    (outdir / "numeric_checks_light.json").write_text(json.dumps(jnum(lc), indent=2, sort_keys=True))
    print("relaxed light:", json.dumps(jnum({k: lc[k] for k in lc if k not in ("own_by_age",)}), indent=1))


if __name__ == "__main__":
    main()
