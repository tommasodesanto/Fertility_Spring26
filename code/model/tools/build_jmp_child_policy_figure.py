"""Policy functions of a childless renter with and without a first child, for the JMP slides (saved 2007 solution, no solve).

Household: inherited renter, age band 26-29 (j=2), middle earnings (z index 4, z=0.78), liquid wealth b.
Two branches at the same b: no child = family state (n=0, m=0); first child born this period = (n=1, m=1).
Objects (same as engine diagnostics.plot_policy_childless_renter, which the standard plots use), all at [b, tenure=0, loc=0, j, z, n, m]:
  rooms       : H_own[tenure_choice-1] if the argmax tenure choice is an owner size, else hR_pol   ("housing after tenure choice")
  consumption : c_pol                                                                                 (consumption if renting)
  owner prob  : tenure_probs[..., 1:].sum()                                                           ("owner-entry policy")
Nodes with V <= -1e9 are masked.  --verify also prints the tenure-probability-weighted (expected) rooms/consumption used in
output/model/rental_menu_precaution_20261007/analyze_branches.py; those expected columns are written to the CSV for reference.

Fourth panel: value of trying for a first child, pi_a * G(b), G = V^H(b,1,1) - xi - V^H(b,0,0), where V^H is the housing-stage
value after the birth outcome (engine household.py: VI, count_active branch). Only the menu value V (birth_count.birth_count_menu,
cap=1) is saved, so V^H is recovered exactly: V(3,3)=VI(3,3); V(n,n)=kappa_c*log(e^{VI(n,n)/kappa_c}+e^{(pi VI(n+1,n+1)+(1-pi)VI(n,n))/kappa_c})
for n=2,1 (kappa_c=kappa_fert_continuation, no fixed cost); V(0,0)=kappa_1*log(e^{VI00/kappa_1}+e^{(pi(VI11-xi)+(1-pi)VI00)/kappa_1}).
The attempt probability is then sigmoid(pi*G/kappa_1); it is checked against fert_probs[...,1] (tolerance 1e-6) before anything is plotted.
Default plots tenure-weighted (expected) rooms/consumption; --standard uses the argmax-tenure / renting-branch objects.

Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
     output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_child_policy_figure.py [--verify]
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
J_IDX, Z_IDX = 2, 4          # age 26-29, z = 0.7776
X_LIMITS = (-3.0, 15.0)     # author (Oct 9): show the negative-wealth region; renters there are buying this period
RENTAL_CAP = 6.0
BRANCHES = {"no child": (0, 0), "first child born this period": (1, 1)}


def interp(b, arr, x):
    i = int(np.clip(np.searchsorted(b, x) - 1, 0, len(b) - 2))
    t = float(np.clip((x - b[i]) / (b[i + 1] - b[i]), 0, 1))
    return (1 - t) * arr[i] + t * arr[i + 1]


def expected(sol, P, price, ib, ten, j, zi, n, m):
    """Tenure-probability-weighted (owner share, rooms, consumption), as in analyze_branches.branch (renter origin: ten=0)."""
    b = np.asarray(sol.b_grid, float).reshape(-1)
    H = np.asarray(P.H_own, float)
    Hval = np.concatenate([[0.0], price * H])
    tp = np.asarray(sol.tenure_probs)[ib, ten, 0, j, zi, n, m].astype(float)
    tp = tp / tp.sum() if tp.sum() > 0 else tp
    out = np.zeros(3)
    for d in range(len(tp)):
        if tp[d] < 1e-10:
            continue
        stay = d == ten
        ba = b[ib] if stay else b[ib] + ((1 - P.psi) * Hval[ten] - Hval[d]) / P.R_gross
        c = interp(b, np.asarray(sol.c_pol_stay if (stay and d > 0) else sol.c_pol)[:, d, 0, j, zi, n, m], ba)
        rooms = interp(b, np.asarray(sol.hR_pol)[:, 0, 0, j, zi, n, m], ba) if d == 0 else H[d - 1]
        out += tp[d] * np.array([d > 0, rooms, c])
    return out


def _lse(a, c, k):
    m = np.maximum(a, c)
    return m + k * np.log(np.exp((a - m) / k) + np.exp((c - m) / k))


def _bisect(target, f, lo_shift=60.0):
    """Solve f(x)=target for increasing f on [target-lo_shift, target]."""
    lo, hi = target - lo_shift, target.copy()
    for _ in range(200):
        x = 0.5 * (lo + hi)
        big = f(x) > target
        hi, lo = np.where(big, x, hi), np.where(big, lo, x)
    return 0.5 * (lo + hi)


def trying_value(sol, P, b):
    """Recover V^H and G(b) for the state; verify the engine's attempt logit against saved fert_probs. Returns dict."""
    V = np.asarray(sol.V)
    ap, rp, fp = (np.asarray(getattr(sol, k)) for k in ("birth_count_action_probs", "birth_count_realized_probs", "fert_probs"))
    k1, kc, xi = float(P.kappa_fert), float(P.kappa_fert_continuation), float(P.first_birth_fixed_cost)
    assert int(P.birth_count_choice_cap) == 1
    v = lambda n, m: V[:, 0, 0, J_IDX, Z_IDX, n, m].astype(float)
    ok = np.all([v(n, n) > -1e8 for n in range(4)], axis=0)
    a1, r1 = ap[:, 0, 0, J_IDX, Z_IDX, 0, 0, 1], rp[:, 0, 0, J_IDX, Z_IDX, 0, 0, 1]
    pis = r1[a1 > 1e-8] / a1[a1 > 1e-8]
    assert np.ptp(pis) < 1e-9, "conception probability not constant across nodes"
    pi = float(pis[0])
    v33 = np.where(ok, v(3, 3), 0.0)
    vi22 = _bisect(np.where(ok, v(2, 2), 0.0), lambda x: _lse(x, pi * v33 + (1 - pi) * x, kc))
    vi11 = _bisect(np.where(ok, v(1, 1), 0.0), lambda x: _lse(x, pi * vi22 + (1 - pi) * x, kc))
    vi00 = _bisect(np.where(ok, v(0, 0), 0.0), lambda x: _lse(x, pi * (vi11 - xi) + (1 - pi) * x, k1))
    G = np.where(ok, vi11 - xi - vi00, np.nan)
    with np.errstate(over="ignore"):
        att = 1.0 / (1.0 + np.exp(-pi * G / k1))
        att1 = 1.0 / (1.0 + np.exp(-pi * (vi22 - vi11) / kc))
    err = float(np.nanmax(np.abs(att - fp[:, 0, 0, J_IDX, Z_IDX, 1])[ok]))
    err_n1 = float(np.nanmax(np.abs(att1 - ap[:, 0, 0, J_IDX, Z_IDX, 1, 1, 1])[ok]))
    print(f"attempt-probability reproduction: max abs err {err:.2e} at {int(ok.sum())} valid nodes "
          f"(n=1 second-attempt check {err_n1:.2e}); pi_a={pi:.5f}, kappa_1={k1:.5f}, kappa_c={kc:.5f}, xi={xi:.5f}")
    if not (err < 1e-6 and err_n1 < 1e-6):
        raise RuntimeError("attempt probability not reproduced; refusing to plot the value panel")
    cross = [b[i] - G[i] * (b[i + 1] - b[i]) / (G[i + 1] - G[i]) for i in range(len(b) - 1)
             if ok[i] and ok[i + 1] and G[i] < 0 <= G[i + 1]]
    print("G crosses zero (upward) at b =", [f"{c:.4f}" for c in cross])
    return dict(G=G, piG=pi * G, pi=pi, ok=ok, cross=cross)


def main(do_verify, use_expected=True, out=None):
    global OUT
    OUT = Path(out) if out else OUT
    result, _ = load_case(SOLUTION)
    sol, P = result.solution, result.P
    price = float(result.price)
    b = np.asarray(sol.b_grid, float).reshape(-1)
    z = np.asarray(sol.type_values, float).reshape(-1)
    H = np.asarray(P.H_own, float)
    V, tc, tp = np.asarray(sol.V), np.asarray(sol.tenure_choice), np.asarray(sol.tenure_probs)
    c_pol, hR = np.asarray(sol.c_pol), np.asarray(sol.hR_pol)
    ix = lambda n, m: (slice(None), 0, 0, J_IDX, Z_IDX, n, m)

    def standard(n, m):
        valid = V[ix(n, m)] > -1e9
        t = tc[ix(n, m)]
        rooms = np.where(t <= 0, hR[ix(n, m)], H[np.maximum(t - 1, 0)])
        own = tp[ix(n, m) + (slice(1, None),)].sum(axis=-1)
        return [np.where(valid, a, np.nan) for a in (rooms, c_pol[ix(n, m)], own)] + [valid]

    if do_verify:
        verify(sol, P, price, b, z)

    if use_expected:   # tenure-probability-weighted rooms and consumption (what tables_branches.md reports)
        def standard(n, m, _std=standard):
            rooms, c, own, valid = _std(n, m)
            e = np.array([expected(sol, P, price, i, 0, J_IDX, Z_IDX, n, m) if valid[i] else [np.nan] * 3
                          for i in range(len(b))])
            return [e[:, 1], e[:, 2], own, valid]

    val = trying_value(sol, P, b)
    sel = (b >= X_LIMITS[0]) & (b <= X_LIMITS[1])
    data = {k: standard(*nm) for k, nm in BRANCHES.items()}
    for k, d in data.items():
        print(f"{k}: valid at all {sel.sum()} nodes in [{X_LIMITS[0]},{X_LIMITS[1]}]? {bool(d[3][sel].all())}; "
              f"first valid b = {b[sel & d[3]].min():.3f}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "font.family": "sans-serif"})
    colors = {"no child": "#9ecae1", "first child born this period": "#08306b"}
    titles = ["Rooms", "Consumption", "Probability of owning", "Value of trying for a child"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.4))
    axes = axes.ravel()
    rows = []
    for (k, d) in data.items():
        n, m = BRANCHES[k]
        ok = sel & d[3]
        for p, ax in enumerate(axes):
            ax.plot(b[ok], d[p][ok], color=colors[k], lw=2.2, label=k)
        for i in np.nonzero(ok)[0]:
            e = expected(sol, P, price, i, 0, J_IDX, Z_IDX, n, m)
            rows.append([k, n, m, f"{b[i]:.6f}", f"{d[0][i]:.6f}", f"{d[1][i]:.6f}", f"{d[2][i]:.6f}",
                         f"{e[1]:.6f}", f"{e[2]:.6f}", f"{e[0]:.6f}",
                         f"{val['G'][i]:.6f}" if k.startswith("first") else "", f"{val['piG'][i]:.6f}" if k.startswith("first") else ""])
    okv = sel & val["ok"]
    axes[3].plot(b[okv], val["piG"][okv], color=colors["first child born this period"], lw=2.2)
    axes[3].axhline(0.0, color="0.35", ls="--", lw=1)
    near = okv & (b >= -1.5)   # y-range from the economically relevant region; the tail near the floor is clipped
    lo, hi = np.nanmin(val["piG"][near]), np.nanmax(val["piG"][near])
    axes[3].set_ylim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
    axes[3].set_ylabel(r"$\pi_a\,[V^H(1,1)-\xi-V^H(0,0)]$", fontsize=10)
    axes[0].axhline(RENTAL_CAP, color="0.35", ls="--", lw=1)
    axes[0].text(X_LIMITS[1], RENTAL_CAP + 0.2, "rental cap", ha="right", va="bottom", fontsize=10, color="0.25")
    axes[2].set_ylim(0, 1)
    for ax, t in zip(axes, titles):
        ax.set_xlim(*X_LIMITS)
        ax.set_xlabel("liquid wealth (model units)")
        ax.set_title(t, loc="left", fontsize=11)
        ax.grid(alpha=0.25, lw=0.6)
    axes[1].legend(frameon=False, loc="lower right", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "child_policy.pdf")
    fig.savefig(OUT / "child_policy.png", dpi=200)
    with (OUT / "child_policy_values.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["branch", "n", "m", "liquid_wealth", "rooms_after_argmax_tenure", "consumption_if_renting",
                    "owner_prob", "rooms_expected_over_tenure", "consumption_expected_over_tenure", "owner_prob_check", "G_value_of_birth_less_cost", "pi_times_G"])
        w.writerows(rows)
    print("wrote", OUT)


def verify(sol, P, price, b, z):
    """Population-weighted b=0, z=0.78 renter row at ages 22-33, as in tables_branches.md (wait / birth)."""
    g = np.asarray(sol.g)[:, :, 0, :, :, 0, 0]
    ib = int(np.nonzero(b == 0.0)[0][0])
    acc = {k: np.zeros(3) for k in BRANCHES}
    mass = 0.0
    for j in (1, 2, 3):
        w = g[ib, 0, j, Z_IDX]
        mass += w
        for k, (n, m) in BRANCHES.items():
            acc[k] += w * expected(sol, P, price, ib, 0, j, Z_IDX, n, m)
    for k in acc:
        o, r, c = acc[k] / mass
        print(f"weighted b=0, z=0.78, ages 22-33, {k}: owner {o:.2f}, rooms {r:.2f}, consumption {c:.3f}")
    print("table row: owner 0.34 / 0.07, rooms 4.38 / 6.00, consumption 1.773 / 1.609")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--expected", action="store_true", help="tenure-weighted rooms and consumption (the default)")
    ap.add_argument("--standard", action="store_true", help="plot argmax-tenure rooms and renting-branch consumption instead of tenure-weighted")
    ap.add_argument("--out", help="output directory (default: the slides folder)")
    a = ap.parse_args()
    main(a.verify, not a.standard, a.out)
