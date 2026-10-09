"""First-birth policy-function figure for the JMP slides, from the saved 2007 solution (no solve).

Quantity: realized-first-birth probability of a childless renter in one 4-year model period, i.e. exactly what the
standard diagnostics label "expected children" / "Fertility choice":
    birth_count_realized_probs[:, tenure=0, loc=0, age_idx, z_idx, n=0, m=0, :] @ arange(n_parity)
(engine diagnostics.plot_policy_childless_renter; readiness gate is off, so child-state index 0).
With the one-birth cap only outcomes 0 and 1 have mass, so this is P(attempt) x P(conception).
The CSV also stores the pre-fecundity attempt probability (fert_probs) for reference.

Run:  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
      output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
      code/model/tools/build_jmp_policy_function_figure.py [--verify]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]
from production.storage import load_case  # noqa: E402

SOLUTION = ROOT / "output/model/production/2007/solution"
OUT = ROOT / "output/model/jmp_slides_policy_functions_20261007"
TARGET_Z = {"low earnings": 0.47, "middle earnings": 0.78, "high earnings": 1.29}
X_LIMITS = (0.0, 22.0)       # renters cannot borrow: b >= 0; liquid wealth, model units; 21.97 = 99.5th pctile of beginning wealth, curves nearly flat above ~15
AGE_FOR_LEFT = 26            # decision age index j = (26-18)/4 = 2, i.e. the 26-29 band
TENURE_RENTER, LOC, N0, M0 = 0, 0, 0, 0


def load():
    result, _ = load_case(SOLUTION)
    sol, P = result.solution, result.P
    assert not bool(getattr(P, "readiness_gate_enabled", False)) and bool(P.birth_count_choice_enabled)
    b = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    z = np.asarray(sol.type_values, dtype=float).reshape(-1)
    realized = np.asarray(sol.birth_count_realized_probs)   # (b, tenure, loc, age, z, n, m, k)
    attempt = np.asarray(sol.fert_probs)                    # (b, tenure, loc, age, z, k)
    valid = np.asarray(sol.V) > -1e9                        # (b, tenure, loc, age, z, n, m)
    kvec = np.arange(int(P.n_parity))
    return sol, P, b, z, realized, attempt, valid, kvec


def line(realized, attempt, valid, kvec, j, zi):
    ix = (slice(None), TENURE_RENTER, LOC, j, zi)
    y = realized[ix + (N0, M0, slice(None))] @ kvec
    a = attempt[ix + (slice(None),)] @ kvec
    ok = valid[ix + (N0, M0)]
    return np.where(ok, y, np.nan), np.where(ok, a, np.nan)


def verify(sol, P, z, realized, kvec):
    """Replay the engine's by-age weighted line and compare with the numbers read off the saved PNG."""
    pre = np.asarray(sol.birth_count_pre_distribution)
    for zi, png in [(3, 0.094), (4, 0.449)]:
        for j, label in [(3, 30)]:
            w = pre[:, :, :, j, zi, 0, 0]
            probs = realized[:, :, :, j, zi, 0, 0, :]
            val = float(np.sum(w * (probs @ kvec)) / np.sum(w))
            print(f"weighted by-age line z={z[zi]:.3f} age {label}: {val:.4f} (PNG ~{png})")
    j, zi = 3, 4
    y, _ = line(realized, np.asarray(sol.fert_probs), np.asarray(sol.V) > -1e9, kvec, j, zi)
    print(f"policy_childless_renter_age30 z=0.778 at b=3000 (right edge): {y[-1]:.4f} (PNG plateau ~0.84)")


def main(do_verify):
    sol, P, b, z, realized, attempt, valid, kvec = load()
    ages = P.age_start + P.da * np.arange(P.J)
    j26 = int(round((AGE_FOR_LEFT - P.age_start) / P.da))
    zidx = {k: int(np.argmin(np.abs(z - v))) for k, v in TARGET_Z.items()}
    if do_verify:
        verify(sol, P, z, realized, kvec)
    ib0 = int(np.searchsorted(b, 0.0))
    assert b[ib0] == 0.0
    ages_f = [j for j in range(int(P.A_f_end)) if np.any(realized[:, 0, 0, j, :, 0, 0, :].sum(-1) > 0)]
    rows = []
    colors = {"low earnings": "#9ecae1", "middle earnings": "#2171b5", "high earnings": "#08306b"}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "font.family": "sans-serif"})
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 3.6))
    for name, zi in zidx.items():
        lw = 2.4 if name == "middle earnings" else 1.8
        lab = f"{name} (z = {z[zi]:.2f})"
        y, a = line(realized, attempt, valid, kvec, j26, zi)
        m = (b >= X_LIMITS[0]) & (b <= X_LIMITS[1]) & np.isfinite(y)
        axl.plot(b[m], y[m], color=colors[name], lw=lw, label=lab)
        for bb, yy, aa in zip(b[m], y[m], a[m]):
            rows.append(["left_by_income", name, f"{z[zi]:.4f}", 26, f"{bb:.6f}", f"{yy:.6f}", f"{aa:.6f}"])
        ys, as_ = [], []
        for j in ages_f:
            yj, aj = line(realized, attempt, valid, kvec, j, zi)
            ys.append(yj[ib0]); as_.append(aj[ib0])
            rows.append(["right_by_age", name, f"{z[zi]:.4f}", int(ages[j]), "0.000000", f"{yj[ib0]:.6f}", f"{aj[ib0]:.6f}"])
        axr.plot(ages[ages_f], ys, color=colors[name], lw=lw, marker="o", ms=4.5, label=lab)
    axl.set_xlim(*X_LIMITS)
    axl.set_xlabel("liquid wealth (model units)")
    axl.set_ylabel("probability of a first birth\nin the 4-year period")
    axl.set_title("By income (renters, age 26-29)", loc="left", fontsize=11)
    axr.set_xlabel("age (start of 4-year period)")
    axr.set_title("By age (renters, liquid wealth = 0)", loc="left", fontsize=11)
    axr.set_xticks(ages[ages_f])
    for ax in (axl, axr):
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.25, lw=0.6)
    axl.legend(frameon=False, loc="lower right", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "first_birth_policy.pdf")
    fig.savefig(OUT / "first_birth_policy.png", dpi=200)
    with (OUT / "first_birth_policy_values.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "income_state", "z", "age", "liquid_wealth", "first_birth_prob_realized", "attempt_prob"])
        w.writerows(rows)
    print("wrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    main(ap.parse_args().verify)
