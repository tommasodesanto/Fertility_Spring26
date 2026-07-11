#!/usr/bin/env python3
"""Locate exactly which j=0 states carry populated-infeasible mass, and print
fertility probabilities at those nodes. Audit-only."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from audit_diag_packet import build_overrides, jnum, VALID_V  # noqa: E402

from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402

BASE = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711")


def main() -> None:
    theta = dict(json.loads((BASE / "repro/current_repro_exact.json").read_text())["theta"])
    overrides, _, _ = build_overrides(120, 10, theta)
    overrides["solve_mode"] = "pe"
    overrides["p_fixed"] = np.array([0.6897624769313723])
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)

    g = np.asarray(sol.g)
    V = np.asarray(sol.V)
    fp = np.asarray(sol.fert_probs)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    invalid = V <= VALID_V

    out = {"rows": []}
    for j in (0, 1):
        gi = g[:, :, :, j] * invalid[:, :, :, j]  # (Nb, nt, I, Nz, npar, ncs)
        nzr = np.argwhere(gi > 1e-10)
        for (b, ten, i, zz, nn, cs) in nzr[:60]:
            out["rows"].append(
                {
                    "j": int(j),
                    "b": float(b_grid[b]),
                    "b_idx": int(b),
                    "ten": int(ten),
                    "z": int(zz),
                    "n": int(nn),
                    "cs": int(cs),
                    "mass": float(gi[b, ten, i, zz, nn, cs]),
                    "V": float(V[b, ten, i, j, zz, nn, cs]),
                    "fert_probs": fp[b, ten, i, j, zz, :].tolist() if nn == 0 else None,
                }
            )
    # mass at j=0 by b node (all states)
    mj0 = np.sum(g[:, :, :, 0], axis=(1, 2, 3, 4, 5))
    out["j0_mass_by_bnode_nonzero"] = {str(float(b_grid[k])): float(mj0[k]) for k in np.flatnonzero(mj0 > 1e-12)}
    # V at entry node for parents cs=1
    k0 = int(np.argmax(b_grid >= 0.0))
    out["V_entry_node_parents"] = {
        f"z{zz}_n{nn}": float(V[k0, 0, 0, 0, zz, nn, 1]) for zz in range(g.shape[4]) for nn in (1, 2)
    }
    (BASE / "diagnostics/current_bound/locate_invalid_j0.json").write_text(json.dumps(jnum(out), indent=2))
    print(json.dumps(jnum(out), indent=1)[:6000])


if __name__ == "__main__":
    main()
