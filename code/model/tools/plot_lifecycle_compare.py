"""Side-by-side comparison of Python vs MATLAB lifecycle arrays."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.theta import apply_theta
from dt_cp_model.utils import make_grid


def python_arrays(theta: np.ndarray, overrides: dict) -> tuple[dict, SimpleNamespace, SimpleNamespace]:
    setup = build_calibration_setup("benchmark")
    P = apply_theta(setup.P_base, theta, setup.names)
    for k, v in overrides.items():
        setattr(P, k, v if not isinstance(v, list) else np.array(v))
    P.max_iter_eq = 200
    sol, P, p_eq = run_model_cp_dt(P, verbose=False)
    b_grid = make_grid(P)
    g = sol.g
    Nb, nt, I, J, npar, ncs = g.shape
    pop_by_age_loc = g.sum(axis=(0, 1, 4, 5))
    pop_by_age = pop_by_age_loc.sum(axis=0)
    own_mass_by_age_loc = g[:, 1:, :, :, :, :].sum(axis=(0, 1, 4, 5))
    own_by_age_loc = np.where(pop_by_age_loc > 0, own_mass_by_age_loc / pop_by_age_loc, 0.0)
    parity_mass_by_age = g.sum(axis=(0, 1, 2, 5))
    parity_share_by_age = parity_mass_by_age / np.where(pop_by_age[:, None] > 0, pop_by_age[:, None], 1.0)
    mean_parity_by_age = (parity_share_by_age * np.arange(npar)[None, :]).sum(axis=1)
    tfr_lc = 2.0 * mean_parity_by_age
    b_mass_by_age = g.sum(axis=(1, 2, 4, 5))
    mean_b = (b_grid[:, None] * b_mass_by_age).sum(axis=0) / np.where(b_mass_by_age.sum(axis=0) > 0, b_mass_by_age.sum(axis=0), 1.0)
    loc_share = pop_by_age_loc / np.where(pop_by_age[None, :] > 0, pop_by_age[None, :], 1.0)
    return {
        "own_by_age_loc": own_by_age_loc,
        "tfr_lc": tfr_lc,
        "parity_share_by_age": parity_share_by_age,
        "mean_b_by_age": mean_b,
        "loc_share_by_age": loc_share,
        "own_by_parity": np.array(sol.own_by_parity),
        "own_rate": float(sol.own_rate),
        "mean_parity": float(sol.mean_parity),
        "pop_share": np.array(sol.pop_share),
    }, sol, P


def matlab_arrays(json_path: Path) -> dict:
    d = json.loads(json_path.read_text())
    own_by_age_loc = np.array(d["own_by_age_loc"])  # (J, I) or (I, J) depending on MATLAB column-major
    if own_by_age_loc.shape[0] == 2:
        # (I, J)
        pass
    else:
        own_by_age_loc = own_by_age_loc.T
    parity_share_by_age = np.array(d["parity_share_by_age"])  # (J, npar)
    if parity_share_by_age.shape[1] != 4 and parity_share_by_age.shape[0] == 4:
        parity_share_by_age = parity_share_by_age.T
    loc_share_by_age = np.array(d["loc_share_by_age"])
    if loc_share_by_age.shape[0] != 2:
        loc_share_by_age = loc_share_by_age.T
    fert_by_age = np.array(d["fert_by_age"])
    # cumulative fertility lifecycle (sum of births)
    tfr_lc = np.cumsum(fert_by_age) * 2.0
    return {
        "own_by_age_loc": own_by_age_loc,
        "tfr_lc": tfr_lc,
        "parity_share_by_age": parity_share_by_age,
        "loc_share_by_age": loc_share_by_age,
        "own_by_parity": np.array(d["own_by_parity"]),
        "own_rate": float(d["own_rate"]),
        "mean_parity": float(d["mean_parity"]),
        "pop_share": np.array(d["pop_share"]),
    }


def make_panel(ax, title, x, py, ml=None, py_label="Python", ml_label="MATLAB", ylabel=""):
    ax.plot(x, py, color="tab:blue", linewidth=2, label=py_label)
    if ml is not None:
        ax.plot(x, ml, color="tab:red", linewidth=1.5, linestyle="--", label=ml_label)
    ax.set_title(title)
    ax.set_xlabel("Age")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta", required=True)
    parser.add_argument("--matlab-json", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--E-loc", default="0,0.1318")
    parser.add_argument("--r-bar", default="0.04,0.0762")
    parser.add_argument("--H-own", default="5.3725,5.84,6.3075,6.875,7.625,8.68")
    parser.add_argument("--hR-max", type=float, default=5.1)
    args = parser.parse_args()

    theta = np.array([float(x) for x in args.theta.split(",")])
    overrides = {
        "E_loc": [float(x) for x in args.E_loc.split(",")],
        "r_bar": [float(x) for x in args.r_bar.split(",")],
        "H_own": [float(x) for x in args.H_own.split(",")],
        "hR_max": args.hR_max,
    }
    py, sol, P = python_arrays(theta, overrides)
    ml = matlab_arrays(args.matlab_json)

    J = P.J
    age_axis = np.arange(J) + P.age_start

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    title = (
        f"Lifecycle Profiles: Python vs MATLAB at bridge-bench theta\n"
        f"Python own={py['own_rate']*100:.1f}% TFR={2*py['mean_parity']:.2f}  |  "
        f"MATLAB own={ml['own_rate']*100:.1f}% TFR={2*ml['mean_parity']:.2f}"
    )
    fig.suptitle(title, fontsize=12)

    ax = axes[0, 0]
    ax.plot(age_axis, py["own_by_age_loc"][1, :], color="tab:red", label="Python Center", linewidth=2)
    ax.plot(age_axis, py["own_by_age_loc"][0, :], color="tab:blue", label="Python Periphery", linewidth=2)
    ax.plot(age_axis, ml["own_by_age_loc"][1, :], color="tab:red", linestyle="--", label="MATLAB Center", alpha=0.7)
    ax.plot(age_axis, ml["own_by_age_loc"][0, :], color="tab:blue", linestyle="--", label="MATLAB Periphery", alpha=0.7)
    ax.set_title("Ownership by age x location")
    ax.set_xlabel("Age"); ax.set_ylabel("Own rate"); ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)

    ax = axes[0, 1]
    ax.plot(age_axis, py["tfr_lc"], color="tab:blue", label="Python", linewidth=2)
    ax.plot(age_axis, ml["tfr_lc"], color="tab:red", linestyle="--", label="MATLAB", alpha=0.8)
    ax.set_title("Fertility lifecycle (TFR cumulative)")
    ax.set_xlabel("Age"); ax.set_ylabel("TFR")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    ax = axes[0, 2]
    width = 0.35
    n_par = py["own_by_parity"].size
    x = np.arange(n_par)
    ax.bar(x - width/2, py["own_by_parity"], width, label="Python", color="tab:blue")
    ax.bar(x + width/2, ml["own_by_parity"], width, label="MATLAB", color="tab:red")
    ax.set_title("Ownership by parity")
    ax.set_xlabel("Parity"); ax.set_ylabel("Own rate"); ax.set_ylim(0, 1)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis="y")

    ax = axes[1, 0]
    n_par = py["parity_share_by_age"].shape[1]
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:purple"][:n_par]
    cum = np.zeros(J)
    for n in range(n_par):
        ax.fill_between(age_axis, cum, cum + py["parity_share_by_age"][:, n], color=colors[n], alpha=0.6, label=f"Py n={n}")
        cum += py["parity_share_by_age"][:, n]
    cum = np.zeros(J)
    for n in range(n_par):
        ax.plot(age_axis, cum + ml["parity_share_by_age"][:, n], color=colors[n], linestyle="--", linewidth=1.2)
        cum += ml["parity_share_by_age"][:, n]
    ax.set_title("Parity composition (fill=Python, dashed=MATLAB)")
    ax.set_xlabel("Age"); ax.set_ylabel("Share")
    ax.set_ylim(0, 1); ax.legend(fontsize=7, loc="center right")

    ax = axes[1, 1]
    ax.plot(age_axis, py["loc_share_by_age"][1, :], color="tab:orange", label="Python Center", linewidth=2)
    ax.plot(age_axis, ml["loc_share_by_age"][1, :], color="tab:orange", linestyle="--", label="MATLAB Center", alpha=0.7)
    ax.set_title("Center share by age")
    ax.set_xlabel("Age"); ax.set_ylabel("Share"); ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    ax = axes[1, 2]
    diff_own_age = py["own_by_age_loc"] - ml["own_by_age_loc"]
    ax.plot(age_axis, diff_own_age[1, :], color="tab:red", label="Center: Py - ML")
    ax.plot(age_axis, diff_own_age[0, :], color="tab:blue", label="Periphery: Py - ML")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("Ownership diff (Python - MATLAB)")
    ax.set_xlabel("Age"); ax.set_ylabel("diff")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.out}")
    print(f"Python:  own={py['own_rate']:.4f} TFR={2*py['mean_parity']:.4f} pop_share={py['pop_share']}")
    print(f"MATLAB:  own={ml['own_rate']:.4f} TFR={2*ml['mean_parity']:.4f} pop_share={ml['pop_share']}")
    py_own_age = py["own_by_age_loc"]
    ml_own_age = ml["own_by_age_loc"]
    diff = py_own_age - ml_own_age
    print(f"own_by_age_loc max abs diff: {np.max(np.abs(diff)):.4f}")
    print(f"parity_share max abs diff:   {np.max(np.abs(py['parity_share_by_age'] - ml['parity_share_by_age'])):.4f}")
    print(f"loc_share max abs diff:      {np.max(np.abs(py['loc_share_by_age'] - ml['loc_share_by_age'])):.4f}")
    print(f"own_by_parity max abs diff:  {np.max(np.abs(py['own_by_parity'] - ml['own_by_parity'])):.4f}")


if __name__ == "__main__":
    main()
