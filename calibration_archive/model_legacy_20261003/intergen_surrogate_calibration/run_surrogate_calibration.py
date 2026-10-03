"""Driver: train the moment emulator, run identification + BO, write a packet.

Usage (from code/model, using the model venv):

    .venv/bin/python -m intergen_surrogate_calibration.run_surrogate_calibration \
        --target-set candidate_no_timing_v0 --n-keep 1500 --restarts 2 --iter 80

Outputs go to output/model/intergen_surrogate_calibration/<target_set>/:

    cv_summary.json, cv_loss_scatter.png, cv_moment_r2.png
    identification.json, surrogate_jacobian_<anchor>.png, fd_vs_surrogate_sv.png
    ard_relevance.png
    bo_proposals.json, pareto_frontier.json, pareto_frontier.png
    REPORT.md

This is a diagnostic methods prototype, not a production SMM calibration.
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import data as D
from . import identification as ID
from . import bo as BO
from .emulator import MomentEmulator

OUT_ROOT = os.path.join(D.REPO_ROOT, "output/model/intergen_surrogate_calibration")

# Conflicting moments for the frontier probe. The first three present in the
# active target set are used (the June-18 ledger flags ownership vs old-age vs
# housing-cost as jointly hard; the replacement set substitutes old wealth and
# owner room separation).
CONFLICT_MOMENT_PREFERENCE = [
    "own_rate", "old_age_own_rate", "housing_user_cost_share",
    "old_nonhousing_wealth_to_income_6575",
    "prime30_55_childless_owner_share_rooms_ge6",
    "old_age_parent_childless_gap", "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "tfr",
]


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-set", default="candidate_no_timing_v0")
    ap.add_argument("--n-keep", type=int, default=1500)
    ap.add_argument("--frac-best", type=float, default=0.5)
    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--iter", type=int, default=80)
    ap.add_argument("--cv-folds", type=int, default=5)
    ap.add_argument("--cv-n", type=int, default=900,
                    help="subset size used for cross-validation (speed).")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip-cv", action="store_true")
    args = ap.parse_args()

    outdir = os.path.join(OUT_ROOT, args.target_set)
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()

    print(f"[1/6] loading records ...")
    recs = D.load_records()
    ds = D.build_dataset(args.target_set, recs)
    print(f"      {ds.X.shape[0]} usable rows; data-set best loss {ds.loss.min():.3f}")

    keep = D.curate_indices(ds.Xs, ds.loss, n_keep=args.n_keep,
                            frac_best=args.frac_best, seed=args.seed)
    print(f"      curated training subset: {keep.size} points")

    report = []
    report.append(f"# Surrogate calibration packet: `{args.target_set}`\n")
    report.append(f"Generated: surrogate-BO prototype. Diagnostic, not a "
                  f"production SMM calibration.\n")
    report.append(f"- Tier-1 records used: **{ds.X.shape[0]}** "
                  f"(target moments all finite)\n")
    report.append(f"- Data-set best loss under this target system: "
                  f"**{ds.loss.min():.3f}**\n")
    report.append(f"- Curated training subset: **{keep.size}** "
                  f"(best {args.frac_best:.0%} by loss + farthest-point fill)\n")
    report.append(f"- Loss = model `diagnostic_loss` (weighted SSE) over "
                  f"{len(ds.moment_names)} moments\n")

    # ---- cross-validation ----------------------------------------------
    if not args.skip_cv:
        print(f"[2/6] cross-validation ({args.cv_folds}-fold) ...")
        cv_idx = keep if keep.size <= args.cv_n else \
            D.curate_indices(ds.Xs, ds.loss, n_keep=args.cv_n,
                             frac_best=args.frac_best, seed=args.seed + 7)
        emu_cv = MomentEmulator(ds, n_restarts=args.restarts,
                                max_iter=args.iter, seed=args.seed, verbose=False)
        cv = emu_cv.cross_validate(cv_idx, k=args.cv_folds, seed=args.seed)
        json.dump(_jsonable({k: v for k, v in cv.items()
                             if k not in ("oof_loss", "true_loss", "idx")}),
                  open(os.path.join(outdir, "cv_summary.json"), "w"), indent=1)
        print(f"      loss R^2={cv['loss_r2']:.3f}  "
              f"rank-corr={cv['loss_rank_corr']:.3f}")

        # Plots.
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(cv["true_loss"], cv["oof_loss"], s=8, alpha=0.4)
        lim = [0, np.nanpercentile(cv["true_loss"], 99)]
        ax.plot(lim, lim, "k--", lw=1)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("true loss"); ax.set_ylabel("out-of-fold predicted loss")
        ax.set_title(f"{args.target_set}\nloss R²={cv['loss_r2']:.3f}, "
                     f"rank-corr={cv['loss_rank_corr']:.3f}")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "cv_loss_scatter.png"),
                                        dpi=120); plt.close(fig)

        names = ds.moment_names
        r2v = [cv["moment_r2"][n] for n in names]
        fig, ax = plt.subplots(figsize=(7, 5))
        order = np.argsort(r2v)
        ax.barh([names[i] for i in order], [r2v[i] for i in order])
        ax.set_xlabel("out-of-fold R²"); ax.set_xlim(min(0, min(r2v)), 1)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_title("Per-moment emulator accuracy (CV)")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "cv_moment_r2.png"),
                                        dpi=120); plt.close(fig)

        report.append("\n## Emulator accuracy (cross-validated)\n")
        report.append(f"- Held-out loss R²: **{cv['loss_r2']:.3f}**, "
                      f"rank correlation: **{cv['loss_rank_corr']:.3f}**\n")
        report.append("- Per-moment out-of-fold R² (worst→best):\n\n")
        report.append("| moment | R² |\n|---|---:|\n")
        for i in np.argsort(r2v):
            report.append(f"| `{names[i]}` | {r2v[i]:.3f} |\n")

    # ---- fit full emulator ---------------------------------------------
    print(f"[3/6] fitting full emulator on {keep.size} points ...")
    emu = MomentEmulator(ds, n_restarts=args.restarts, max_iter=args.iter,
                         seed=args.seed, verbose=True)
    emu.fit(keep)

    # ---- identification -------------------------------------------------
    print(f"[4/6] identification vs finite-difference audit ...")
    fd = ID.load_fd_audit()
    id_out = {"anchors": {}}
    sv_compare = {}
    # The FD audit was computed under one specific target system; only do the
    # head-to-head when the active target set matches its moment list.
    fd_target_set = fd.get("config", {}).get("target_set")
    fd_comparable = (fd_target_set == args.target_set)
    if not fd_comparable:
        print(f"      FD audit target set is '{fd_target_set}', not "
              f"'{args.target_set}'; skipping FD head-to-head, "
              f"reporting surrogate-only identification.")
    for label in (fd["points"] if fd_comparable else []):
        cmp = ID.compare_to_fd(emu, fd, label)
        id_out["anchors"][label] = {
            "surrogate_rank": cmp["surrogate_rank"],
            "surrogate_cond": cmp["surrogate_cond"],
            "surrogate_smallest_sv": cmp["surrogate_smallest_sv"],
            "fd_rank": cmp["fd_rank"],
            "fd_cond": cmp["fd_cond"],
            "fd_vs_surrogate_pearson": cmp.get("fd_vs_surrogate_pearson"),
            "fd_vs_surrogate_sign_agreement": cmp.get("fd_vs_surrogate_sign_agreement"),
        }
        # Heatmap of the surrogate scaled Jacobian.
        S = cmp["surrogate_scaled_matrix"]
        fig, ax = plt.subplots(figsize=(8, 6))
        vmax = np.nanpercentile(np.abs(S), 98) or 1.0
        im = ax.imshow(S, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(ds.theta_params)))
        ax.set_xticklabels(ds.theta_params, rotation=90, fontsize=7)
        ax.set_yticks(range(len(ds.moment_names)))
        ax.set_yticklabels(ds.moment_names, fontsize=7)
        ax.set_title(f"Surrogate scaled Jacobian @ {label}")
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"surrogate_jacobian_{label}.png"),
                    dpi=120); plt.close(fig)
        # Singular-value spectra: surrogate.
        sv_compare[label] = {"surrogate_sv":
                             np.linalg.svd(S, compute_uv=False).tolist()}

    # ARD relevance.
    rel = ID.ard_relevance(emu)
    fig, ax = plt.subplots(figsize=(7, 4))
    order = np.argsort(rel["param_relevance"])
    ax.barh([rel["param_names"][i] for i in order],
            [rel["param_relevance"][i] for i in order])
    ax.set_xlabel("mean ARD relevance (1/lengthscale, unit cube)")
    ax.set_title("Parameter relevance across moments (emulator)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "ard_relevance.png"),
                                    dpi=120); plt.close(fig)
    id_out["ard_relevance"] = {
        "param_names": rel["param_names"],
        "param_relevance": rel["param_relevance"].tolist(),
        "lengthscales": rel["lengthscales"].tolist(),
    }
    json.dump(_jsonable(id_out),
              open(os.path.join(outdir, "identification.json"), "w"), indent=1)

    report.append("\n## Identification: surrogate vs finite-difference audit\n")
    if not fd_comparable:
        report.append(f"_FD audit was computed under `{fd_target_set}`; this run "
                      f"uses `{args.target_set}`, so no head-to-head. "
                      f"ARD relevance below is surrogate-only._\n")
    else:
        report.append("| anchor | FD rank | FD cond | surrogate rank | surrogate cond | "
                      "Jacobian sign-agree | Pearson |\n|---|---:|---:|---:|---:|---:|---:|\n")
        for label, a in id_out["anchors"].items():
            def f(x):
                return "—" if x is None else (f"{x:.3g}" if isinstance(x, float) else str(x))
            report.append(f"| {label} | {f(a['fd_rank'])} | {f(a['fd_cond'])} | "
                          f"{a['surrogate_rank']} | {f(a['surrogate_cond'])} | "
                          f"{f(a['fd_vs_surrogate_sign_agreement'])} | "
                          f"{f(a['fd_vs_surrogate_pearson'])} |\n")
    pr = id_out["ard_relevance"]
    rorder = np.argsort(pr["param_relevance"])
    report.append("\n- ARD relevance (least→most identified parameter): "
                  + ", ".join(f"`{pr['param_names'][i]}`"
                              for i in rorder) + "\n")

    # ---- BO + frontier --------------------------------------------------
    print(f"[5/6] Bayesian-optimization proposals + Pareto frontier ...")
    f_best = float(ds.loss[keep].min())
    batch = BO.propose_batch(emu, f_best=f_best, q=args.batch, seed=args.seed)
    exploit = BO.propose_exploit(emu, seed=args.seed + 3)
    local = BO.propose_local(emu, keep, q=args.batch, seed=args.seed + 9)
    bo_out = {
        "f_best_in_training": f_best,
        "global_ei_batch": batch,
        "exploit_point": exploit,
        "local_batch": local,
        # Back-compat key consumed by validate_proposals (global EI batch).
        "ei_batch": batch,
    }
    json.dump(_jsonable(bo_out),
              open(os.path.join(outdir, "bo_proposals.json"), "w"), indent=1)

    # Save the fitted emulator so refinement/validation never re-fits.
    emu.save(os.path.join(outdir, "emulator.pkl"))

    conflict_moments = [m for m in CONFLICT_MOMENT_PREFERENCE
                        if m in ds.moment_names][:3]
    pareto = BO.pareto_front(emu, conflict_moments, seed=args.seed + 5)
    json.dump(_jsonable(pareto),
              open(os.path.join(outdir, "pareto_frontier.json"), "w"), indent=1)

    # Pareto plot (pairwise projections of the conflict objectives).
    om = pareto["objective_moments"]
    nd = np.array(pareto["nondominated_objs"])
    if nd.ndim == 2 and nd.shape[0] > 0:
        pairs = [(0, 1), (0, 2), (1, 2)]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, (i, j) in zip(axes, pairs):
            ax.scatter(nd[:, i], nd[:, j], s=10, alpha=0.5)
            ax.set_xlabel(f"|{om[i]} − target|")
            ax.set_ylabel(f"|{om[j]} − target|")
        fig.suptitle("Emulated Pareto front of ledger-conflict moments")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "pareto_frontier.png"),
                                        dpi=120); plt.close(fig)

    report.append("\n## Bayesian-optimization proposals\n")
    report.append(f"- Best loss in training pool (incumbent): **{f_best:.3f}**\n")
    report.append("- Two proposal modes. **Global EI** searches the whole box; in "
                  "13-D the data cloud is sparse and the smooth emulator "
                  "extrapolates (and discrete owner-room moments snap to a "
                  "different rung than predicted), so global predicted losses are "
                  "optimistic and must be confirmed. **Local** proposals perturb "
                  "the best incumbents inside the data manifold (trust region) -- "
                  "the reliable, realistic use with real-solver confirmation.\n")
    report.append(f"- Surrogate exploit point predicted loss: "
                  f"**{exploit['pred_loss']:.3f}** (global argmin; expect "
                  f"over-optimism).\n\n")
    report.append("Global EI batch:\n\n| # | EI | predicted loss | nn-dist to data |\n"
                  "|---:|---:|---:|---:|\n")
    for b in batch:
        report.append(f"| {b['rank']} | {b['ei']:.3g} | {b['pred_loss']:.3f} | — |\n")
    report.append("\nLocal trust-region batch:\n\n"
                  "| # | EI | predicted loss | nn-dist to data |\n|---:|---:|---:|---:|\n")
    for b in local:
        report.append(f"| {b['rank']} | {b['ei']:.3g} | {b['pred_loss']:.3f} | "
                      f"{b['nn_dist_to_data']:.3f} |\n")

    report.append("\n## Frontier probe: are the ledger-conflict moments jointly reachable?\n")
    report.append(f"Objectives: {', '.join('`'+m+'`' for m in om)} "
                  f"(absolute deviation from target).\n\n")
    report.append("| moment | marginal best achievable | at best-joint point | target |\n")
    report.append("|---|---:|---:|---:|\n")
    for m in om:
        report.append(f"| `{m}` | {pareto['marginal_min'][m]:.3f} | "
                      f"{pareto['best_joint_objs'][m]:.3f} | "
                      f"{ds.targets[m]:.3f} |\n")
    report.append(f"\nNondominated emulated points: {pareto['n_nondominated']} "
                  f"of {pareto['n_pool']} sampled.\n")

    # ---- write report ---------------------------------------------------
    report.append(f"\n---\n_Total wall time: {time.time()-t0:.0f}s. "
                  f"Emulator: {len(emu.gps)} ARD-GPs on {keep.size} points._\n")
    with open(os.path.join(outdir, "REPORT.md"), "w") as fh:
        fh.write("".join(report))

    print(f"[6/6] done in {time.time()-t0:.0f}s. Packet: {outdir}")
    print(f"      report: {os.path.join(outdir, 'REPORT.md')}")


if __name__ == "__main__":
    main()
