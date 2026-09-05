#!/usr/bin/env python3
"""Build the bounded-refinement readout from verified saved tables, without solves."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[3]
NAMES = {
    "tfr": "Completed fertility", "childless_rate": "Childlessness",
    "mean_age_first_birth": "Mean age at first birth",
    "share_first_births_age30plus": "First births at age 30+",
    "housing_increment_0to1": "Rooms response to first birth",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Rooms, 3+ minus 1-2 children",
    "own_family_gap": "Family ownership gap", "own_rate": "Ownership, ages 30-55",
    "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms",
    "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / annual earnings",
    "annual_bequest_flow_to_aggregate_wealth": "Bequest flow / wealth",
    "old_total_wealth_to_annual_income_p90_p50_7684": "Older wealth/income p90 / p50",
}
PARAMETERS = {
    "beta_annual": r"Annual patience $\beta$", "kappa_fert": r"First-birth dispersion $\kappa_1$",
    "kappa_fert_continuation": r"Later-birth dispersion $\kappa_C$",
    "chi": r"Owner-service premium $\chi$", "H0": r"Supply intercept $H_0$",
    "theta0": r"Bequest strength $\theta_0$", "theta1": r"Bequest shifter $\theta_1$",
    "hbar_child_rooms": "Per-child room floor", "hbar_first_child_jump": "First-child room intercept",
    "first_birth_fixed_cost": "First-birth fixed cost", "psi_child_change_2023": "2007-2023 child-value change",
    "tenure_choice_kappa": r"Tenure dispersion $\kappa$", "housing_supply_elasticity": r"Dated supply elasticity $\eta$",
}


def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def localize(path, run):
    marker = run.name + "/"
    if marker not in str(path):
        raise RuntimeError("Result path is outside the named run")
    return run / str(path).split(marker, 1)[1]


def fit_table(rows):
    lines = ["| Moment | Target | Model | Gap | Weight | Loss |", "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append("| " + " | ".join([NAMES[r["moment"]], f"{float(r['target']):.6f}",
            f"{float(r['model']):.6f}", f"{float(r['gap']):+.6f}", f"{float(r['weight']):.6g}",
            f"{float(r['loss_contribution']):.4f}"]) + " |")
    return "\n".join(lines)


def coordinate_identification(run):
    """Read-only finite differences; these matrices never generate candidates."""
    import numpy as np
    stage = run / "round1"
    summary = read_json(stage / "report/summary.json")
    if summary["status"] != "complete" or summary["completed_cases"] != 23:
        raise RuntimeError("Identification diagnostic requires all 23 verified cases")
    values = {}
    for i in range(1, 24):
        task = stage / f"task_{i:03d}"
        s = read_json(task / "summary.json")
        receipt = read_json(task / "case_receipt.json")
        for name in ("summary.json", "target_fit_long.csv"):
            if sha(task / name) != receipt["artifact_sha256"][name]:
                raise RuntimeError("Coordinate evidence hash mismatch")
        fit = read_csv(task / "target_fit_long.csv")
        contract = [(r["moment"], float(r["target"]), float(r["weight"])) for r in fit]
        values[i] = (np.asarray(s["panel_design"]["unit_vector"]),
                     np.asarray([float(r["standardized_gap"]) for r in fit]), contract, s)
    u0, r0, contract, anchor = values[1]
    J = np.empty((12, 11)); sides = []
    for j, spec in enumerate(anchor["panel_design"]["domain"]):
        um, rm, cm, _ = values[2+2*j]; up, rp, cp, _ = values[3+2*j]
        if cm != contract or cp != contract:
            raise RuntimeError("Mixed coordinate target contracts")
        if not np.array_equal(np.delete(um,j),np.delete(u0,j)) or not np.array_equal(np.delete(up,j),np.delete(u0,j)):
            raise RuntimeError("Not a single-coordinate comparison")
        forward, backward = (rp-r0)/(up[j]-u0[j]), (r0-rm)/(u0[j]-um[j])
        J[:,j] = (rp-rm)/(up[j]-um[j])
        norm = np.linalg.norm(forward)*np.linalg.norm(backward)
        sides.append(dict(parameter=spec["name"], forward_backward_cosine=float(np.dot(forward,backward)/norm) if norm else None,
                          central_norm=float(np.linalg.norm(J[:,j]))))
    _, singular, vt = np.linalg.svd(J, full_matrices=False)
    diagnostic = dict(scope="Local weighted finite differences, not global identification or standard errors",
        algebraic_rank=int(np.linalg.matrix_rank(J)), relative_rank_1e_minus3=int(np.sum(singular >= singular[0]*1e-3)),
        singular_values=singular.tolist(), weakest_squared_loadings={spec["name"]:float(v*v) for spec,v in zip(anchor["panel_design"]["domain"],vt[-1])},
        step_agreement=sides, parameter_count=11, moment_count=12)
    out = run / "report"; out.mkdir(exist_ok=True)
    (out / "coordinate_identification.json").write_text(json.dumps(diagnostic,indent=2)+"\n")
    np.savetxt(out / "weighted_coordinate_jacobian.csv", J, delimiter=",")
    return diagnostic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--coordinate-diagnostics-only", action="store_true")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    identification = coordinate_identification(run)
    if args.coordinate_diagnostics_only:
        print(json.dumps(identification))
        return
    if args.source is None:
        raise RuntimeError("A final Markdown source path is required")
    stages = [(name, read_json(run / name / "report/summary.json")) for name in ("smoke", "round1", "round2", "repeats")]
    if any(s["status"] != "complete" for _, s in stages):
        raise RuntimeError("Final readout requires all planned stages collected completely")
    best = read_json(run / "best_so_far.json")["best"]
    selected_path = localize(best["summary"], run)
    selected = read_json(selected_path)
    if sha(selected_path) != best["summary_sha256"]:
        raise RuntimeError("Selected summary hash mismatch")
    task = selected_path.parent
    for i in (1, 2):
        receipt = read_json(run / f"repeats/task_{i:03d}/case_receipt.json")
        if not receipt.get("reference", {}).get("exact_twelve_row_fit"):
            raise RuntimeError("Final exact repeat was not certified")
        if receipt["reference"]["reference_sha256"] != best["summary_sha256"]:
            raise RuntimeError("Repeats do not belong to the cross-round best candidate")
    production = ROOT / "output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010"
    production_summary = read_json(production / "summary.json")
    production_fit, selected_fit = read_csv(production / "target_fit_long.csv"), read_csv(task / "target_fit_long.csv")
    if len(production_fit) != 12 or len(selected_fit) != 12:
        raise RuntimeError("The complete twelve-row target fit is required")
    for old, new in zip(production_fit, selected_fit):
        if (old["moment"], old["target"], old["weight"]) != (new["moment"], new["target"], new["weight"]):
            raise RuntimeError("Target or weight changed")
    lookup = lambda rows: {r["moment"]: float(r["model"]) for r in rows}
    old, new = lookup(production_fit), lookup(selected_fit)
    old_loss = float(production_summary["best_candidate"]["transition_loss"])
    loss = float(selected["best_candidate"]["transition_loss"])
    param = read_csv(task / "parameter_table.csv")
    old_param = {r["parameter"]: r for r in read_csv(production / "parameter_table.csv")}
    if sum(r["is_free_parameter"] == "True" for r in param) != 11:
        raise RuntimeError("Free parameter count changed")
    improvement = 100*(1-loss/old_loss)
    first_target = float(next(r["target"] for r in selected_fit if r["moment"] == "housing_increment_0to1"))
    graphdir = run / "report"
    graphdir.mkdir(exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(8, 4.7), constrained_layout=True)
    pos = np.arange(12)
    ax.barh(pos+.17, [float(r["loss_contribution"]) for r in production_fit], .32, label="Retained production", color="#9CAAB3")
    ax.barh(pos-.17, [float(r["loss_contribution"]) for r in selected_fit], .32, label="Verified new candidate", color="#197C80")
    ax.set_yticks(pos, [NAMES[r["moment"]] for r in selected_fit], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("Contribution to the unchanged twelve-moment loss")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
    fig.savefig(graphdir / "loss_contributions.pdf"); fig.savefig(graphdir / "loss_contributions.png", dpi=180); plt.close(fig)
    plot = (graphdir / "loss_contributions.pdf").relative_to(ROOT)
    counts = sum(s["completed_cases"] for _,s in stages)
    ptable = ["| Parameter | Production | New candidate | Bounds / restriction | Near bound? |", "|---|---:|---:|---|---|"]
    for r in param:
        name = r["parameter"]
        if name not in PARAMETERS:
            continue
        free = r["is_free_parameter"] == "True"
        bounds = f"[{float(r['lower_bound']):g}, {float(r['upper_bound']):g}]" if free else "Externally fixed"
        near = ("Raw scale only" if name == "theta1" and r["near_bound"] == "True" else "Yes" if r["near_bound"] == "True" else "No") if free else "Not estimated"
        ptable.append(f"| {PARAMETERS[name]} | {float(old_param[name]['value']):.6f} | {float(r['value']):.6f} | {bounds} | {near} |")
    progress = ["| Evaluation stage | Completed cases | Best loss in stage |", "|---|---:|---:|"]
    labels = {"smoke": "Inherited-candidate exact repeats", "round1": "Small coordinate moves", "round2": "Joint parameter moves", "repeats": "Final exact repeats"}
    for name,s in stages:
        progress.append(f"| {labels[name]} | {s['completed_cases']} | {float(s['best']['loss']):.6f} |")
    flagged = []
    for name, _ in stages:
        for case_dir in sorted((run / name).glob("task_*")):
            receipt = read_json(case_dir / "case_receipt.json")
            for filename, expected in receipt["artifact_sha256"].items():
                if sha(case_dir / filename) != expected:
                    raise RuntimeError(f"Report input changed: {case_dir / filename}")
            screen = receipt["policy_array_diagnostic"]
            if screen["occupied_negative_steps"]:
                flagged.append({"case": str(case_dir.relative_to(run)), **screen})
    selected_screen = read_json(task / "case_receipt.json")["policy_array_diagnostic"]
    max_flag_mass = max((r["share_pre_choice_mass_at_negative_steps"] for r in flagged), default=0.)
    selected_budget = read_json(task / "case_receipt.json")["budget_diagnostic"]
    spike_path = task / "supplemental_spike_exposure.json"
    spike_note = ""
    if spike_path.exists():
        spike = read_json(spike_path)
        if spike["checkpoint_sha256"] != sha(task / "dated_state.pkl.gz"):
            raise RuntimeError("Supplemental spike diagnostic uses another checkpoint")
        spike_note = (f"The standard age-30 first-birth panel also shows a probability spike "
            f"at {spike['state_count']} negative-wealth renter state. Its settled-childless "
            f"pre-choice mass is {spike['pre_choice_mass']:.6g}. This is a check of current "
            "exposure, not a proof of optimal policies outside occupied states. The original "
            "graphs and their dense income-state labels are preserved.")
    standard_names = [
        "fertility_by_age", "ownership_by_age", "housing_market", "housing_prices",
        "market_clearing_by_market", "market_clearing_residuals", "owner_rungs", "tenure_services",
        "fertility_policy_by_age_income_state", "housing_by_age_income_state",
        "ownership_by_age_income_state", "liquid_wealth_by_age_income_state", "income_state_outcomes",
        "policy_childless_renter_age30", "wealth_dist_childless_renter_age30",
        "policy_childless_renter_age42", "wealth_dist_childless_renter_age42",
    ]
    if {p.stem for p in (task / "standard_diagnostics").glob("*.png")} != set(standard_names):
        raise RuntimeError("The standard seventeen-graph set changed")
    graph_appendix = []
    for i, name in enumerate(standard_names):
        if i % 2 == 0:
            graph_appendix.append("\\clearpage\n\n# Standard solution diagnostics\n")
        path = (task / "standard_diagnostics" / f"{name}.png").relative_to(ROOT)
        graph_appendix.append(f"![{name.replace('_', ' ').capitalize()}.]({path}){{width=88%}}\n")
    doc = f'''---
title: "What the bounded calibration refinement achieved"
subtitle: "Same model, targets, weights, and external restrictions"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The result for discussion

The best verified candidate has loss **{loss:.6f}**, compared with **{old_loss:.6f}**
for the retained production calibration: a **{improvement:.2f}% reduction** in the
same twelve-row objective. Two fresh executions reproduce every fit row,
parameter, preference normalization, and saved numerical history entry exactly.
This is a bounded local improvement, with weak parameter identification still
unresolved. It is recorded for review; production has not been silently replaced.

The first-birth housing response moves from **{old['housing_increment_0to1']:.4f}**
to **{new['housing_increment_0to1']:.4f} rooms**, against **{first_target:.4f}** targeted.
Mean occupied rooms move from **{old['aggregate_mean_occupied_rooms_18_85']:.4f}**
to **{new['aggregate_mean_occupied_rooms_18_85']:.4f}**, against **{float(next(r['target'] for r in selected_fit if r['moment']=='aggregate_mean_occupied_rooms_18_85')):.4f}** targeted.
The first-birth response explains **{100*new['housing_increment_0to1']/first_target:.1f}%**
of its target. A lower total objective must be assessed alongside this remaining gap.

The largest gain comes from the rooms contrast between families with three or more
children and families with one or two: its model value changes from
**{old['prime30_55_parent_3plus_minus_1to2_mean_rooms']:.4f}** to
**{new['prime30_55_parent_3plus_minus_1to2_mean_rooms']:.4f} rooms**, against
**{float(next(r['target'] for r in selected_fit if r['moment']=='prime30_55_parent_3plus_minus_1to2_mean_rooms')):.4f}** targeted.
Mean occupied rooms remain too high and become slightly higher at this candidate.

{chr(10).join(progress)}

All **{counts}** completed evaluations retain the original source and target
fingerprints, twelve moments, eleven free parameters, the 120-point wealth grid,
dated housing-supply elasticity 0.63, and tenure-choice dispersion 0.005. They reconstruct
the old state and full 2007-2023 path. Every case has a terminal checkpoint,
complete target and parameter tables, and the unchanged seventeen-figure packet.
No policy path or long-run response was recalculated in this exercise.

\\clearpage

# Where the fit changes

![Supplemental comparison of every objective contribution. Both series use identical target values and weights.]( {plot} ){{width=100%}}

The first round tests small separate parameter moves. The second combines
successful moves and tests a child-space rebalancing: increase the first-child
intercept while reducing the per-child floor. With
$\\bar h(m)=\\bar h_J+m\\bar h_C$ for children at home $m>0$, the proposed changes
$d\\bar h_J=d$ and $d\\bar h_C=-d/3$ preserve the three-child floor and raise the
first-child floor by $2d/3$. Actual occupied rooms still depend on prices, saving,
tenure choice, and fertility selection; the floor identity is not an identity
for observed housing use.

All ranking uses the complete unchanged objective. No target was removed,
reweighted, or added. The earlier ridge planner's rejection is retained; these
direct candidate evaluations do not invert its unstable Jacobian or establish
precise estimates of the weakly informed bequest parameter.

\\clearpage

# Complete target fit: verified candidate

Gaps equal model minus target. Contributions are weight times squared gap.
The weights retain their existing empirical or synthetic interpretations;
this loss is not a specification-test statistic.

{fit_table(selected_fit)}

# Complete target fit: retained production

{fit_table(production_fit)}

\\clearpage

# Complete parameter comparison

{chr(10).join(ptable)}

The old fertility-preference intercept is re-normalized to replacement for
every parameter vector. New candidate values are
$\\psi_{{2007}}={float(selected['old_psi_child']):.8f}$ and
$\\psi_{{2023}}={float(selected['best_candidate']['new_psi_child']):.8f}$.
The stored near-bound flag uses raw parameter units. For the bequest shifter,
interpret it alongside its position in the logarithmic search interval; a
raw-scale flag alone does not establish a binding search constraint.

The supply restriction 0.63 applies to dated markets. The inherited initialization
retains elasticity 1.75 before the dated supply intercept is normalized. Both
values and that order are unchanged across this comparison; whether this is the
intended economic contract remains a discussion item in the separate code audit.

The complete smaller-step panel has algebraic Jacobian rank
**{identification['algebraic_rank']} of 11** and relative rank
**{identification['relative_rank_1e_minus3']} of 11** at the declared $10^{{-3}}$
threshold. These are local, step-dependent diagnostics, not standard errors or
proof of global identification. The saved receipt includes every singular value,
weak-direction loading, and forward/backward derivative comparison.
The weakest direction has **{100*identification['weakest_squared_loadings']['theta1']:.1f}%**
of its squared loading on the bequest shifter. The supply-intercept derivative
has forward/backward cosine **{next(r['forward_backward_cosine'] for r in identification['step_agreement'] if r['parameter']=='H0'):.3f}**;
opposite directional responses make a smooth local inverse unreliable here.

\\clearpage

# Numerical checks and their limits

All completed cases pass the unchanged production acceptance checks. Those gates
do not certify every policy function. The separate adjacent-wealth screen flags
occupied value decreases in **{len(flagged)}** evaluations; the largest exposed
pre-choice mass share is **{max_flag_mass:.3g}**. These flags are retained and are
consistent with the already documented need to check the local saving optimizer;
this search does not establish their cause through a new global-optimizer solve.
The selected candidate has **{selected_screen['occupied_negative_steps']}** flagged
occupied value decreases. Its reporting-floor budget-excess share is
**{selected_budget['budget_excess_share']:.3g}**. Probability arrays are finite
and remain inside $[0,1]$. The seventeen standard graphs follow this note unchanged.

{spike_note}

The concurrent source-code audit also reproduces an incomplete policy-storage
interface: later-birth probabilities reside in a mutable parameter field instead
of the saved calendar policy. Reuse after an intervening solve can mix policies.
No affected result has been demonstrated in this historical calibration: the
nonzero preference drift clears warm starts at each later date. The interface
requires a separate repair and targeted replay before relying on vulnerable
policy-reuse paths. Exact reproduction alone would not detect this defect.

# What remains to discuss

The first-birth rooms contrast remains a provisional empirical-to-model mapping.
The young-ownership diagnostic is still excluded from estimation; model and
empirical age windows need a common contract, including existing prime-age and
old-wealth statistics. The overnight evidence on near-universal late-life
ownership and sensitivity to rental opportunities remains relevant. A better
in-sample objective does not independently resolve those issues.

Supply and credit results in the earlier morning audit belong to the retained
production calibration. This refinement supplies no revised policy magnitudes,
welfare ranking, or resolution of the household-entry transition.

The complete per-case fit and parameter tables are in each stage's
`report/all_target_fits.csv` and `report/all_parameters.csv` under
`output/model/e5f_morning_refinement_20260905a/`. Its README provides stage plans,
execution receipts, checkpoints, standard plots, and reproduction commands.
The selected case is **{selected_path.relative_to(run)}**.

{chr(10).join(graph_appendix)}
'''
    # Standard Markdown image syntax; avoid spaces inside a link target.
    doc = doc.replace(f"]( {plot} )", f"]({plot})")
    args.source.parent.mkdir(parents=True, exist_ok=True)
    args.source.write_text(doc.rstrip() + "\n")
    receipt = {"selected_summary_sha256": sha(selected_path), "source_sha256": sha(args.source),
        "loss": loss, "production_loss": old_loss, "improvement_percent": improvement,
        "completed_evaluations": counts, "target_count": 12, "free_parameter_count": 11,
        "flagged_value_screens": flagged,
        "fit_input_sha256": {str(p): sha(p) for p in (task / "target_fit_long.csv", task / "parameter_table.csv", production / "target_fit_long.csv", production / "parameter_table.csv")}}
    (graphdir / "report_data_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(args.source)


if __name__ == "__main__":
    main()
