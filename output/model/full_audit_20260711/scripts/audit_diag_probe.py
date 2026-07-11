#!/usr/bin/env python3
"""Follow-up probes on flagged pathologies, solved once at fixed p_eq (PE mode).

Probes: (a) where the populated-infeasible mass (g>0 on V=-1e10 nodes) lives;
(b) magnitude/location of huge V drops among 'valid' nodes; (c) fert_probs rows
summing to zero and their populated mass; (d) decomposition of p_birth=0 mass.
Writes probe_checks.json next to the main packet. Audit-only.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from audit_diag_packet import build_overrides, jnum, VALID_V  # noqa: E402

from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402

RECORD = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711/repro/current_repro_exact.json")
OUT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711/diagnostics/current_bound")
P_EQ = 0.6897624769313723


def main() -> None:
    theta = dict(json.loads(RECORD.read_text())["theta"])
    overrides, targets, weights = build_overrides(120, 10, theta)
    overrides["solve_mode"] = "pe"
    overrides["p_fixed"] = np.array([P_EQ])
    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    wall = time.perf_counter() - t0

    g = np.asarray(sol.g)
    V = np.asarray(sol.V)
    fp = np.asarray(sol.fert_probs)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    Nb, nt, I, J, Nz, npar, ncs = g.shape
    total = float(np.sum(g))
    out: dict = {"wall_seconds": wall, "p_fixed": P_EQ,
                 "residual": float(sol.best_max_abs_rel_excess)}

    # (a) populated infeasible mass: where
    invalid = V <= VALID_V
    pi = g * invalid
    out["populated_invalid"] = {
        "total_mass": float(np.sum(pi)),
        "share": float(np.sum(pi) / total),
        "by_tenure": np.sum(pi, axis=(0, 2, 3, 4, 5, 6)).tolist(),
        "by_age": np.sum(pi, axis=(0, 1, 2, 4, 5, 6)).tolist(),
        "by_parity": np.sum(pi, axis=(0, 1, 2, 3, 4, 6)).tolist(),
        "by_child_state": np.sum(pi, axis=(0, 1, 2, 3, 4, 5)).tolist(),
        "by_z": np.sum(pi, axis=(0, 1, 2, 3, 5, 6)).tolist(),
        "wealth_range_of_mass": [
            float(b_grid[np.min(np.flatnonzero(np.sum(pi, axis=(1, 2, 3, 4, 5, 6)) > 1e-12))]),
            float(b_grid[np.max(np.flatnonzero(np.sum(pi, axis=(1, 2, 3, 4, 5, 6)) > 1e-12))]),
        ] if np.any(np.sum(pi, axis=(1, 2, 3, 4, 5, 6)) > 1e-12) else None,
    }

    # (b) V value distribution among valid nodes
    valid = ~invalid
    vv = V[valid]
    out["V_valid_distribution"] = {
        "min": float(np.min(vv)),
        "p0001": float(np.quantile(vv, 1e-4)),
        "p01": float(np.quantile(vv, 0.01)),
        "median": float(np.median(vv)),
        "max": float(np.max(vv)),
        "n_below_minus_1e6": int(np.sum(vv < -1e6)),
        "populated_mass_below_minus_1e6": float(np.sum(g[valid & (V < -1e6)])),
    }

    # (c) fert prob rows summing to zero
    rowsum = np.sum(fp, axis=-1)  # (Nb, nt, I, J, Nz)
    zero_rows = rowsum < 1e-12
    # childless mass at those rows within fertile window
    zmass = 0.0
    cmass = 0.0
    for j in range(P.A_f_start - 1, P.A_f_end + 1):
        gj = g[:, :, :, j, :, 0, 0]
        zmass += float(np.sum(gj[zero_rows[:, :, :, j, :]]))
        cmass += float(np.sum(gj))
    out["fert_rowsum_zero"] = {
        "cells": int(np.sum(zero_rows)),
        "cell_share": float(np.mean(zero_rows)),
        "childless_fertile_mass_on_zero_rows": zmass,
        "childless_fertile_mass": cmass,
    }

    # (d) p_birth = 0 decomposition (fertile window, childless)
    p_birth = 1.0 - fp[..., 0]
    dec = {"by_z": np.zeros(Nz), "by_tenure": np.zeros(nt), "by_age": np.zeros(J),
           "valid_node_mass": 0.0, "invalid_node_mass": 0.0}
    for j in range(P.A_f_start - 1, P.A_f_end + 1):
        z0 = p_birth[:, :, :, j, :] <= 1e-12  # (Nb, nt, I, Nz)
        gj = g[:, :, :, j, :, 0, 0]
        m = gj * z0
        dec["by_z"] += np.sum(m, axis=(0, 1, 2))
        dec["by_tenure"] += np.sum(m, axis=(0, 2, 3))
        dec["by_age"][j] = float(np.sum(m))
        vj = valid[:, :, :, j, :, 0, 0]
        dec["valid_node_mass"] += float(np.sum(m[vj]))
        dec["invalid_node_mass"] += float(np.sum(m[~vj]))
    # wealth profile of the p_birth=0 childless mass vs all childless
    wz = np.zeros(Nb)
    wall_c = np.zeros(Nb)
    for j in range(P.A_f_start - 1, P.A_f_end + 1):
        z0 = p_birth[:, :, :, j, :] <= 1e-12
        gj = g[:, :, :, j, :, 0, 0]
        wz += np.sum(gj * z0, axis=(1, 2, 3))
        wall_c += np.sum(gj, axis=(1, 2, 3))
    def wmean(w):
        s = float(np.sum(w))
        return float(np.sum(w * b_grid) / s) if s > 1e-14 else np.nan
    dec["mean_wealth_pbirth0"] = wmean(wz)
    dec["mean_wealth_all_childless_fertile"] = wmean(wall_c)
    out["p_birth_zero_decomposition"] = dec

    (OUT / "probe_checks.json").write_text(json.dumps(jnum(out), indent=2, sort_keys=True))
    print(json.dumps(jnum(out), indent=1))


if __name__ == "__main__":
    main()
