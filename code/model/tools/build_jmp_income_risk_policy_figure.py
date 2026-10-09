"""First-birth policy function, baseline vs lower income risk, for the JMP slides (saved solutions, no solve).

Object: same as first_birth_policy.pdf (left panel), i.e. birth_count_realized_probs[:, renter, loc0, j=2 (age 26-29), z_idx, n=0, m=0, :] @ arange(4)
= P(attempt) x P(conception) for a childless renter, against liquid wealth b >= 0, at the SAME income-state INDEX (middle state, 4) in both arms.
Baseline  : output/model/production/2007/solution (loaded with production.storage.load_case, as in build_jmp_policy_function_figure.py).
Lower risk: "risk A" (sigma_eps 0.484 -> 0.35, grid renormalized to mean 1), saved cell
            output/model/reconciliation_14p40_20261007/precautionary/cells/riskA_fixed/solution_arrays.npz (fixed price, rebate T at base).
The z value of the middle state differs across arms because the risk-A grid is renormalized; both are printed and written to the CSV.

Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
     output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_income_risk_policy_figure.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]
from build_jmp_policy_function_figure import AGE_FOR_LEFT, OUT, line, load  # noqa: E402

RISKA = ROOT / "output/model/reconciliation_14p40_20261007/precautionary/cells/riskA_fixed/solution_arrays.npz"
Z_IDX = 4
X_MAX = 22.0


def main():
    sol, P, b, z, realized, attempt, valid, kvec = load()
    j = int(round((AGE_FOR_LEFT - P.age_start) / P.da))
    assert j == 2
    d = np.load(RISKA)
    bA, zA = d["b_grid"].reshape(-1), d["type_values"].reshape(-1)
    assert np.array_equal(b, bA)
    realA, attA, validA = d["birth_count_realized_probs"], d["fert_probs"], d["V"] > -1e9
    y0, a0 = line(realized, attempt, valid, kvec, j, Z_IDX)
    y1, a1 = line(realA, attA, validA, kvec, j, Z_IDX)
    print(f"middle state index {Z_IDX}: z baseline = {z[Z_IDX]:.4f}, z lower risk = {zA[Z_IDX]:.4f}")
    # verification against the published first_birth_policy values (middle earnings, left panel)
    ref = {}
    with (OUT / "first_birth_policy_values.csv").open() as f:
        for r in csv.DictReader(f):
            if r["panel"] == "left_by_income" and r["income_state"] == "middle earnings":
                ref[round(float(r["liquid_wealth"]), 6)] = float(r["first_birth_prob_realized"])
    m = (b >= 0) & (b <= X_MAX) & np.isfinite(y0)
    diff = max(abs(ref[round(float(x), 6)] - float(v)) for x, v in zip(b[m], y0[m]))
    print(f"max abs diff, baseline curve vs first_birth_policy middle-earnings curve: {diff:.3e} over {int(m.sum())} nodes")
    m1 = (b >= 0) & (b <= X_MAX) & np.isfinite(y1)
    print("nodes with b>=0 valid: baseline", int(m.sum()), "lower risk", int(m1.sum()))
    print(f"p at b=0: baseline {y0[b == 0][0]:.4f}, lower risk {y1[b == 0][0]:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "font.family": "sans-serif"})
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(b[m], y0[m], color="#2171b5", lw=2.4, label="baseline income risk")
    ax.plot(b[m1], y1[m1], color="#d95f0e", lw=2.4, label="lower income risk")
    ax.set_xlim(0, X_MAX)
    ax.set_ylim(0, 1)
    ax.set_xlabel("liquid wealth (model units)")
    ax.set_ylabel("probability of a first birth\nin the 4-year period")
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "income_risk_policy.pdf")
    fig.savefig(OUT / "income_risk_policy.png", dpi=200)
    with (OUT / "income_risk_policy_values.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arm", "z_state_index", "z", "age", "liquid_wealth", "first_birth_prob_realized", "attempt_prob"])
        for arm, zz, bb, yy, aa, mm in [("baseline", z[Z_IDX], b, y0, a0, m), ("lower income risk (risk A)", zA[Z_IDX], b, y1, a1, m1)]:
            for x, yv, av in zip(bb[mm], yy[mm], aa[mm]):
                w.writerow([arm, Z_IDX, f"{zz:.4f}", 26, f"{x:.6f}", f"{yv:.6f}", f"{av:.6f}"])
    print("wrote", OUT / "income_risk_policy.pdf")


if __name__ == "__main__":
    main()
