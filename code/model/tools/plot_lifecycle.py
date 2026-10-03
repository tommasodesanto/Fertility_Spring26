"""Reproduce the slide-deck six-panel lifecycle figure from Python sol.

Equivalent of MATLAB's diag11_lifecycle.png. Computes everything from sol.g
(the steady-state distribution) so it works for any solved theta.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dt_cp_model.evaluate import solve_theta
from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.theta import apply_theta
from dt_cp_model.solver import run_model_cp_dt


def lifecycle_arrays(sol: SimpleNamespace, P: SimpleNamespace, b_grid: np.ndarray) -> dict:
    g = sol.g  # (Nb, nt, I, J, npar, ncs)
    Nb, nt, I, J, npar, ncs = g.shape

    pop_by_age_loc = g.sum(axis=(0, 1, 4, 5))  # (I, J)
    pop_by_age = pop_by_age_loc.sum(axis=0)  # (J,)

    own_mass_by_age_loc = g[:, 1:, :, :, :, :].sum(axis=(0, 1, 4, 5))  # (I, J)
    pop_by_loc_age = g.sum(axis=(0, 1, 4, 5))  # (I, J)
    own_by_age_loc = np.where(pop_by_loc_age > 0, own_mass_by_age_loc / pop_by_loc_age, 0.0)

    parity_mass_by_age = g.sum(axis=(0, 1, 2, 5))  # (J, npar)
    parity_share_by_age = parity_mass_by_age / np.where(pop_by_age[:, None] > 0, pop_by_age[:, None], 1.0)

    mean_parity_by_age = (parity_share_by_age * np.arange(npar)[None, :]).sum(axis=1)
    tfr_lifecycle = 2.0 * mean_parity_by_age  # one-shot model: 2x mean parity

    b_axis = b_grid
    b_mass_by_age = g.sum(axis=(1, 2, 4, 5))  # (Nb, J)
    mean_b_by_age = (b_axis[:, None] * b_mass_by_age).sum(axis=0) / np.where(
        b_mass_by_age.sum(axis=0) > 0, b_mass_by_age.sum(axis=0), 1.0
    )

    loc_share_by_age = pop_by_age_loc / np.where(pop_by_age[None, :] > 0, pop_by_age[None, :], 1.0)

    own_by_parity = sol.own_by_parity

    return {
        "own_by_age_loc": own_by_age_loc,
        "tfr_lifecycle": tfr_lifecycle,
        "parity_share_by_age": parity_share_by_age,
        "mean_b_by_age": mean_b_by_age,
        "loc_share_by_age": loc_share_by_age,
        "own_by_parity": own_by_parity,
    }


def plot_lifecycle(arrays: dict, P: SimpleNamespace, sol: SimpleNamespace, out_path: Path, title: str) -> None:
    J = P.J
    age_axis = np.arange(J) + P.age_start

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title, fontsize=14)

    own_by_age_loc = arrays["own_by_age_loc"]
    own_overall = float(sol.own_rate * 100.0)
    ax = axes[0, 0]
    ax.plot(age_axis, own_by_age_loc[1, :], color="tab:red", label="Center", linewidth=2)
    ax.plot(age_axis, own_by_age_loc[0, :], color="tab:blue", label="Periphery", linewidth=2)
    ax.set_xlabel("Age")
    ax.set_ylabel("Own rate")
    ax.set_title(f"Ownership (agg={own_overall:.0f}%)")
    ax.set_ylim(0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    tfr_lc = arrays["tfr_lifecycle"]
    tfr_final = float(2.0 * sol.mean_parity)
    ax = axes[0, 1]
    ax.plot(age_axis, tfr_lc, color="tab:blue", linewidth=2)
    ax.set_xlabel("Age")
    ax.set_ylabel("TFR (2 x mean parity)")
    ax.set_title(f"Fertility lifecycle (TFR={tfr_final:.2f})")
    ax.grid(True, alpha=0.3)

    parity_share = arrays["parity_share_by_age"]
    npar = parity_share.shape[1]
    ax = axes[0, 2]
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:purple"][:npar]
    cum = np.zeros(J)
    for n in range(npar):
        ax.fill_between(age_axis, cum, cum + parity_share[:, n], color=colors[n], label=f"n={n}", alpha=0.85)
        cum += parity_share[:, n]
    ax.set_xlabel("Age")
    ax.set_ylabel("Share")
    ax.set_title("Parity composition")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="center right")

    mean_b = arrays["mean_b_by_age"]
    ax = axes[1, 0]
    ax.plot(age_axis, mean_b, color="tab:blue", linewidth=2)
    ax.set_xlabel("Age")
    ax.set_ylabel("Mean b (financial)")
    ax.set_title("Financial wealth lifecycle")
    ax.grid(True, alpha=0.3)

    loc_share = arrays["loc_share_by_age"]
    c_share_total = float(sol.pop_share[1] * 100.0)
    ax = axes[1, 1]
    ax.fill_between(age_axis, 0.0, loc_share[0, :], color="tab:blue", label="Periphery", alpha=0.85)
    ax.fill_between(age_axis, loc_share[0, :], 1.0, color="tab:orange", label="Center", alpha=0.85)
    ax.set_xlabel("Age")
    ax.set_ylabel("Share")
    ax.set_title(f"Location (C share={c_share_total:.0f}%)")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="center right")

    own_by_par = arrays["own_by_parity"]
    p0 = float(own_by_par[0])
    p_max = float(own_by_par.max())
    gap = p_max - p0
    ax = axes[1, 2]
    ax.bar(np.arange(len(own_by_par)), own_by_par, color="tab:blue")
    ax.set_xlabel("Parity")
    ax.set_ylabel("Own rate")
    ax.set_title(f"Ownership by parity (gap={gap:.2f})")
    ax.set_ylim(0, 1.0)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup", default="fast")
    parser.add_argument("--theta", default="x0", help="comma-separated theta or 'x0'")
    parser.add_argument("--max-iter-eq", type=int, default=120)
    parser.add_argument("--out", type=Path, default=Path("benchmarks/diag11_lifecycle_python.png"))
    args = parser.parse_args()

    setup = build_calibration_setup(args.setup)
    if args.theta.strip().lower() == "x0":
        theta = setup.x0
    else:
        theta = np.array([float(x) for x in args.theta.split(",")])

    P = apply_theta(setup.P_base, theta, setup.names)
    P.max_iter_eq = args.max_iter_eq
    sol, P, p_eq = run_model_cp_dt(P, verbose=False)

    from dt_cp_model.utils import make_grid

    b_grid = make_grid(P)

    arrays = lifecycle_arrays(sol, P, b_grid)
    title = f"Fig 11: Lifecycle Profiles (Python, setup={args.setup}, theta={args.theta})"
    plot_lifecycle(arrays, P, sol, args.out, title)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
