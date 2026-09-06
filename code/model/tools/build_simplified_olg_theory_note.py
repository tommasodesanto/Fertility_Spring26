#!/usr/bin/env python3
"""Build the integrated illustrative-theory note and analytical transition figure.

Run from any directory. This uses derivatives of the unchanged household and
equilibrium equations, not finite simulated reforms or a calibrated economy.
The main slides and existing proof/example files are not rebuilt.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import schur

from verify_simplified_olg_mixed_transition import anchor, linearization, values
from verify_simplified_olg_local_transition import young_choices, complex_jacobian

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"
TMP = ROOT / "tmp/pdfs/simplified_olg_polish"
SOURCE = ROOT / "latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex"
PDF = ROOT / "output/pdf/simplified_olg_amendment_proposal.pdf"
BLUE, RED, GRAY = "#26527a", "#b0493d", "#7a7a7a"


def transition_data():
    p, Z, hh = anchor(4)
    lin = linearization(p, Z, hh)
    v = Z[[1, 2]]
    Fv = complex_jacobian(lambda z: values(Z, z, p)[0], v)
    Gv = complex_jacobian(lambda z: values(Z, z, p)[1], v)
    Ft = complex_jacobian(
        lambda z: values(Z, v, dict(p, theta=z[0]))[0], [p["theta"]]).ravel()
    Gt = complex_jacobian(
        lambda z: values(Z, v, dict(p, theta=z[0]))[1], [p["theta"]]).ravel()
    Dtheta = Gt - Gv @ np.linalg.solve(Fv, Ft)
    stationary_theta = np.linalg.solve(np.eye(6) - lin["J"], Dtheta)
    S, E, count = schur(lin["J"], output="real",
                        sort=lambda re, im: np.hypot(re, im) < 1)
    assert count == 4
    Es, stable = E[:, :count], S[:count, :count]

    def response(endpoint):
        c = np.linalg.solve(lin["B"] @ Es, -lin["B"] @ endpoint)
        return np.array([endpoint + Es @ np.linalg.matrix_power(stable, t) @ c
                         for t in range(42)])

    # These are directional derivatives as epsilon -> 0 with delta=.5*epsilon.
    # The policy-date derivative is the baseline plus the phi response. Their
    # interaction is second order. No size of a finite admissible shock is claimed.
    baseline = -response(stationary_theta)
    credit = response(lin["derivative"])
    policy_date, relative_reform = 1, .5
    policy = baseline.copy()
    policy[policy_date:] += relative_reform * credit[:-policy_date]
    end0 = -stationary_theta
    end1 = end0 + relative_reform * lin["derivative"]
    fertility0 = np.diff(baseline[:, 2]) / (p["nu"] * Z[2])
    fertility1 = np.diff(policy[:, 2]) / (p["nu"] * Z[2])
    base_error = np.max(abs(baseline[1:] - baseline[:-1] @ lin["J"].T + Dtheta))
    post_error = np.max(abs(policy[policy_date+1:] -
                           policy[policy_date:-1] @ lin["J"].T + Dtheta -
                           relative_reform * lin["Dphi"]))
    boundary_error = np.max(abs(lin["B"] @ (policy[policy_date] - baseline[policy_date])))
    assert max(base_error, post_error, boundary_error) < 1e-11
    assert np.max(abs(policy[:policy_date] - baseline[:policy_date])) == 0
    assert abs(fertility0[0] - fertility1[0]) < 1e-13
    assert fertility0[0] < 0 and fertility1[policy_date] > fertility0[policy_date]
    assert sum(end0[2:4]) < 0 and sum(end1[2:4]) > sum(end0[2:4])
    # At the chosen intervention date, the baseline cohort/population is falling.
    assert sum(baseline[2, 2:4]) < sum(baseline[1, 2:4]) < 0

    def stationary_schedule(z):
        P, N, theta, phi = z
        T = p["q"] * p["tau"] * P * p["Hbar"] / N
        choices = young_choices([P] * 3, [T] * 2, dict(p, theta=theta, phi=phi))
        pi = choices["pi"]
        old_h = pi * choices["owner"]["z"][5] + (1-pi) * choices["renter"]["z"][5]
        return np.array([N / 2 * (choices["housing"] + old_h) - p["Hbar"],
                         choices["fertility"]])

    point = np.array([Z[0], Z[2]+Z[3], p["theta"], p["phi"]])
    J = complex_jacobian(stationary_schedule, point)
    finite_errors = []
    for k in range(4):
        step = np.zeros(4); step[k] = 1e-5
        finite_errors.append(np.max(abs((stationary_schedule(point+step) -
                                        stationary_schedule(point-step))/2e-5 - J[:, k])))
    assert max(finite_errors) < 1e-7
    curves = []
    for name, theta, phi, endpoint in (("Initial", 0., 0., np.zeros(6)),
                                      ("Baseline", -1., 0., end0),
                                      ("Policy", -1., relative_reform, end1)):
        shift = theta * J[:, 2] + phi * J[:, 3]
        NP, Ns = -J[0, 0]/J[0, 1], -shift[0]/J[0, 1]
        nP, ns = J[1, 0] + J[1, 1]*NP, shift[1]+J[1, 1]*Ns
        Pstar = -ns/nP
        Nstar = NP*Pstar+Ns
        assert max(abs(Pstar-endpoint[0]), abs(Nstar-sum(endpoint[2:4]))) < 1e-10
        curves.append(dict(name=name, N_P=NP, N_shift=Ns, n_P=nP, n_shift=ns))
    return p, Z, baseline, policy, end0, end1, fertility0, fertility1, curves, dict(
        interpretation="Analytical directional derivatives at the mixed-tenure economy, delta=.5 epsilon and policy_date=1; not finite reforms or calibrated paths.",
        surprise_boundary="Pre-policy choices use baseline price expectations. The actual-state recurrence is imposed from the announcement onward, not across the surprise at policy_date-1.",
        policy_date=policy_date, relative_reform=relative_reform,
        baseline_map_error=float(base_error), post_policy_map_error=float(post_error),
        common_inherited_boundary_error=float(boundary_error),
        stationary_schedule_central_difference_error=float(max(finite_errors)),
        baseline_states=baseline.tolist(), policy_states=policy.tolist(),
        baseline_fertility=fertility0.tolist(), policy_fertility=fertility1.tolist(),
        baseline_endpoint=end0.tolist(), policy_endpoint=end1.tolist(),
        curves=curves, stationary_schedule_jacobian=J.tolist(),
        parameters={k: float(v) for k, v in p.items()})


def draw_transition():
    p, Z, base, policy, end0, end1, n0, n1, curves, receipt = transition_data()
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": 14, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42,
                         "legend.frameon": False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    fig.subplots_adjust(left=.07, right=.98, bottom=.19, top=.84, wspace=.30)
    Ngrid, Pgrid = np.linspace(-5.6, .5, 250), np.linspace(-3.6, .5, 250)
    for curve, color in zip(curves, (GRAY, BLUE, RED)):
        axes[0].plot(Ngrid, (Ngrid-curve["N_shift"])/curve["N_P"],
                     color=color, lw=1.5, ls="--", alpha=.85, label=curve["name"])
        axes[1].plot(Pgrid, curve["n_P"]*Pgrid+curve["n_shift"],
                     color=color, lw=1.5, ls="--", alpha=.85)

    # Arrows compare the initial, policy-date, and limiting equilibria. They
    # do not trace every date or claim monotone adjustment between these points.
    def arrows(ax, points, color):
        for a, b in zip(points[:-1], points[1:]):
            ax.annotate("", xy=b, xytext=a,
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.45,
                                        shrinkA=2, shrinkB=2, mutation_scale=10))

    for panel, ax in enumerate(axes):
        def position(states, fertility):
            return np.column_stack((states[:-1, 2]+states[:-1, 3], states[:-1, 0])) if panel == 0 else np.column_stack((states[:-1, 0], fertility))
        b = position(base, n0); c = position(policy, n1)
        old = np.array([0., 0.])
        terminal0 = np.array([sum(end0[2:4]), end0[0]]) if panel == 0 else np.array([end0[0], 0.])
        terminal1 = np.array([sum(end1[2:4]), end1[0]]) if panel == 0 else np.array([end1[0], 0.])
        arrows(ax, np.vstack((old, b[:2], terminal0)), BLUE)
        arrows(ax, np.vstack((b[1], c[1], terminal1)), RED)
        ax.plot(*b[0], ".", ms=5, color=BLUE, zorder=5)
        ax.annotate(r"$t=0$", b[0], xytext=(-31, -5 if panel == 0 else -16),
                    textcoords="offset points", color=BLUE, fontsize=12)
        for point, label, color, offset in (
            (old, r"$S_-$", GRAY, (5, 8)),
            (terminal0, r"$S_0$", BLUE, (-15, -20)),
            (terminal1, r"$S_1$", RED, (6, 7)),
            (b[1], r"$I^0$", BLUE, (5, -17)),
            (c[1], r"$I^1$", RED, (6, 4)),
        ):
            ax.plot(*point, "o", ms=4.5, color=color, zorder=5)
            ax.annotate(label, point, xytext=offset, textcoords="offset points",
                        color=color, fontsize=14, zorder=6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(length=0)
    axes[0].set(xlim=(-5.6, .5), ylim=(-3.6, .6),
                xlabel=r"Adult households, $N_{\rm hh}$", ylabel=r"House price, $P$",
                title="(a) Housing and population")
    axes[1].set(xlim=(-3.6, .5), ylim=(-.98, .39), xlabel=r"House price, $P$",
                ylabel=r"Mean fertility, $\bar n$", title="(b) Prices and fertility")
    axes[1].axhline(0, color=GRAY, lw=.8, ls=":")
    axes[1].set_yticks([0], [r"$1/\nu$"])
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=3,
               bbox_to_anchor=(.53, 1.02), title="Stationary schedules", fontsize=12,
               title_fontsize=12)
    for ext in ("pdf", "png"):
        fig.savefig(OUT/f"combined_transition_figure.{ext}", dpi=190, bbox_inches="tight")
    plt.close(fig)
    return receipt


def compile_note():
    engine = shutil.which("pdflatex") or "/Library/TeX/texbin/pdflatex"
    env = dict(os.environ)
    env["PATH"] = str(Path(engine).parent)+os.pathsep+env.get("PATH", "")
    for run in (1, 2):
        with (TMP/f"build-pass{run}.txt").open("w") as log:
            subprocess.run([engine, "-interaction=nonstopmode", "-halt-on-error",
                            "-file-line-error", f"-output-directory={TMP}", str(SOURCE)],
                           cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    log = (TMP/f"{SOURCE.stem}.log").read_text()
    bad = [line for line in log.splitlines() if any(word in line for word in
           ("Overfull", "undefined", "multiply defined", "LaTeX Error"))]
    if bad:
        raise RuntimeError("LaTeX checks failed: "+"; ".join(bad))
    shutil.copy2(TMP/f"{SOURCE.stem}.pdf", PDF)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); TMP.mkdir(parents=True, exist_ok=True)
    receipt = draw_transition()
    if not args.figure_only:
        compile_note()
    sources = (Path(__file__), SOURCE,
               Path(__file__).with_name("verify_simplified_olg_mixed_transition.py"),
               Path(__file__).with_name("verify_simplified_olg_local_transition.py"))
    receipt["source_hashes"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in sources}
    (OUT/"integrated_note_figure_checks.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("Analytical transition checks passed; figure written.")
    if not args.figure_only:
        print(f"Compiled twice without overflow or undefined references: {PDF}")


if __name__ == "__main__":
    main()
