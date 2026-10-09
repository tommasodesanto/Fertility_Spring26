"""Owning and first-birth policy functions, financed share 80% vs 95%, for the JMP slides (saved solutions, no solve).

Household: childless renter (n=0, m=0, tenure 0, loc 0), age band 26-29 (j=2), middle earnings (z index 4), liquid wealth b in [-3, 15], valid nodes only.
Left : end-of-period owning probability, no-child branch: tenure_probs[b, 0, 0, 2, 4, 0, 0, 1:].sum()   (as in build_jmp_child_policy_figure.py).
Right: first-birth probability in the period: birth_count_realized_probs[b, 0, 0, 2, 4, 0, 0, :] @ arange(4) (as in build_jmp_policy_function_figure.py).
Baseline (financed share 0.80): output/model/rental_menu_precaution_20261007/cells/A6/solution_arrays.npz
Credit   (financed share 0.95): .../cells/A6_LTV95/solution_arrays.npz (fixed price, everything else equal).
Check: the 80% curves are compared with the production solution (output/model/production/2007/solution) loaded as in the other two scripts.

Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
     output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_credit_policy_figure.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]
from build_jmp_policy_function_figure import OUT, load  # noqa: E402

CELLS = ROOT / "output/model/rental_menu_precaution_20261007/cells"
J_IDX, Z_IDX = 2, 4
X_LIMITS = (-3.0, 15.0)
IX = (slice(None), 0, 0, J_IDX, Z_IDX)


def curves(V, tp, rp):
    valid = V[IX + (0, 0)] > -1e9
    own = tp[IX + (0, 0, slice(1, None))].sum(axis=-1)
    birth = rp[IX + (0, 0, slice(None))] @ np.arange(rp.shape[-1])
    return np.where(valid, own, np.nan), np.where(valid, birth, np.nan), valid


def cell(name):
    d = np.load(CELLS / name / "solution_arrays.npz")
    return d["b_grid"].reshape(-1), curves(d["V"], d["tenure_probs"], d["birth_count_realized_probs"])


def main():
    b0, (own0, bir0, v0) = cell("A6")
    b1, (own1, bir1, v1) = cell("A6_LTV95")
    assert np.array_equal(b0, b1)
    b = b0
    sol, P, bp, z, realized, attempt, valid, kvec = load()
    assert np.array_equal(b, bp)
    ownp, birp, vp = curves(np.asarray(sol.V), np.asarray(sol.tenure_probs), realized)
    for nm, a, c in [("owning", own0, ownp), ("first birth", bir0, birp)]:
        both = np.isfinite(a) & np.isfinite(c)
        print(f"A6 vs production, {nm}: max abs diff {np.max(np.abs(a[both] - c[both])):.3e} over {int(both.sum())} nodes; "
              f"valid-mask equal: {bool(np.array_equal(np.isfinite(a), np.isfinite(c)))}")
    sel = (b >= X_LIMITS[0]) & (b <= X_LIMITS[1])
    m0, m1 = sel & v0, sel & v1
    print("valid nodes in range: 80%", int(m0.sum()), "95%", int(m1.sum()))
    both = m0 & m1
    for nm, a, c in [("owning", own0, own1), ("first birth", bir0, bir1)]:
        gap = np.where(both, c - a, np.nan)
        i = int(np.nanargmax(np.abs(gap)))
        print(f"largest |gap| {nm}: {gap[i]:+.4f} at b={b[i]:.3f} (80%={a[i]:.4f}, 95%={c[i]:.4f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "font.family": "sans-serif"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, (t, y0, y1) in zip(axes, [("Probability of owning", own0, own1), ("Probability of a first birth", bir0, bir1)]):
        ax.plot(b[m0], y0[m0], color="#9ecae1", lw=2.4, label="financed share 80%")
        ax.plot(b[m1], y1[m1], color="#08306b", lw=2.4, label="financed share 95%")
        ax.set_xlim(*X_LIMITS)
        ax.set_ylim(0, 1)
        ax.set_xlabel("liquid wealth (model units)")
        ax.set_title(t, loc="left", fontsize=11)
        ax.grid(alpha=0.25, lw=0.6)
    axes[0].legend(frameon=False, loc="upper left", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "credit_policy.pdf")
    fig.savefig(OUT / "credit_policy.png", dpi=200)
    with (OUT / "credit_policy_values.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["financed_share", "liquid_wealth", "owner_prob_no_child", "first_birth_prob_realized"])
        for sh, mm, o, bb in [(0.80, m0, own0, bir0), (0.95, m1, own1, bir1)]:
            for i in np.nonzero(mm)[0]:
                w.writerow([sh, f"{b[i]:.6f}", f"{o[i]:.6f}", f"{bb[i]:.6f}"])
    print("wrote", OUT / "credit_policy.pdf")


if __name__ == "__main__":
    main()
