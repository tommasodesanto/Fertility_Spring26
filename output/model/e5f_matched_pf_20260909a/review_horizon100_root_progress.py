"""Validate collected root receipts and rebuild the provisional advisor tables.

Run from any directory with Python and NumPy. No model solves or source edits.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys

import numpy as np

BASE = Path(__file__).resolve().parent
REPO = BASE.parents[2]
sys.path.insert(0, str(REPO / "tmp/e5f_matched_pf/code/model/tools"))
from collect_e5f_matched_pf_price_jacobian import gates_pass

ROOT = BASE / "meeting_receipts/historical_root_h100_01/sequential"
LABELS = {
    "tfr": "Completed fertility",
    "childless_rate": "Childless share",
    "mean_age_first_birth": "Mean age at first birth",
    "share_first_births_age30plus": "First births at age 30+ (share)",
    "housing_increment_0to1": "First-birth housing response (rooms)",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Rooms gap: 3+ versus 1–2 children, ages 30–55",
    "own_family_gap": "Parent ownership gap (share units)",
    "own_rate": "Ownership share",
    "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms, ages 18–85",
    "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / annual gross labor earnings",
    "annual_bequest_flow_to_aggregate_wealth": "Annual bequests / wealth",
    "old_total_wealth_to_annual_income_p90_p50_7684": "Wealth / income dispersion, ages 76–84 (p90/p50)",
}


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def fmt(value):
    return "—" if value == "" else f"{float(value):.8g}"


def main():
    trials = []
    reference = None
    for folder in sorted(ROOT.glob("evaluation_*")):
        if not (folder / "summary.json").exists():
            continue
        summary = json.loads((folder / "summary.json").read_text())
        contract = json.loads((folder / "contract.json").read_text())
        assert summary["mapping_valid"] is True
        gates_pass(summary["gates"])
        for name, expected in summary["artifact_sha256"].items():
            artifact = folder / name
            if artifact.exists():
                assert hashlib.sha256(artifact.read_bytes()).hexdigest() == expected
        fits = read_csv(folder / "target_fit.csv")
        assert len(fits) == 12 and {r["moment"] for r in fits} == set(LABELS)
        identity = ([tuple(r[k] for k in ("moment", "target", "weight")) for r in fits],
                    (folder / "parameters.csv").read_bytes(), contract["source_sha256"],
                    contract["contract_sha256"])
        if reference is None:
            reference = identity
        assert identity == reference
        contributions = []
        for row in fits:
            gap = float(row["model"]) - float(row["target"])
            contribution = float(row["weight"]) * gap ** 2
            assert np.isclose(gap, float(row["gap"]), rtol=1e-12, atol=1e-12)
            assert np.isclose(contribution, float(row["loss_contribution"]), rtol=1e-12, atol=1e-12)
            contributions.append(contribution)
        assert np.isclose(sum(contributions), summary["loss"], rtol=1e-12, atol=1e-12)
        path = read_csv(folder / "transition_path.csv")
        assert len(path) == 100
        assert [int(r["calendar_year"]) for r in path] == summary["years"]
        residual = np.array([(float(r["housing_demand"]) - float(r["housing_supply"]))
                             / float(r["housing_supply"]) for r in path])
        assert np.array_equal(residual, summary["residual"])
        assert np.array_equal([float(r["asset_price"]) for r in path], summary["prices"])
        score = float(np.max(np.abs(residual)))
        assert np.isclose(score, summary["maximum_market_residual"], rtol=0, atol=1e-12)
        trials.append(dict(evaluation=int(folder.name.split("_")[-1]),
                           elapsed_seconds=summary["elapsed_seconds"],
                           maximum_market_residual=score, mapping_gates_pass=True,
                           market_gate_pass=score <= 2e-4,
                           full_fit_and_parameter_directory=str(folder),
                           available_artifacts_verified=True,
                           terminal_distance=summary["terminal_distance"]))
    assert trials
    best = min(trials, key=lambda r: r["maximum_market_residual"])
    folder = Path(best["full_fit_and_parameter_directory"])
    fits = read_csv(folder / "target_fit.csv")
    parameters = read_csv(folder / "parameters.csv")
    summary = json.loads((folder / "summary.json").read_text())
    review = dict(status="collected_root_trials_validated_not_final_certification",
                  verified_trials=trials, unchanged_parameters_targets_weights_and_sources=True,
                  best_collected_evaluation=best["evaluation"],
                  qualification="This collector verifies individual trials. Final root replay, horizon stability, calibration and policy certification require separate evidence.")
    (BASE / "horizon100_root_progress_review.json").write_text(json.dumps(review, indent=2) + "\n")
    lines = ["# 100-date perfect-foresight solution: provisional readout", "",
             "This is a price-solver diagnostic at inherited parameters, not a new calibration. "
             "The tables below use the collected trial with the smallest market residual. "
             "Final reproduction and horizon stability are not certified by this document.", "",
             "| Trial | Maximum market gap | Mapping checks | Seconds |",
             "|---|---:|---|---:|"]
    lines += [f"| {r['evaluation']} | {100*r['maximum_market_residual']:.6f}% | Pass | {r['elapsed_seconds']:.1f} |" for r in trials]
    lines += ["", "The market tolerance is **0.02%**. Parameters, empirical targets and weights are identical across these trials.", "",
              f"## Full target fit — trial {best['evaluation']}", "",
              f"Objective at these provisional prices: **{summary['loss']:.8f}**. Shares remain in fraction units. This objective is not a calibrated-equilibrium loss.", "",
              "| Moment | Target | Model | Gap | Weight | Loss contribution |",
              "|---|---:|---:|---:|---:|---:|"]
    lines += ["| " + " | ".join([LABELS[r["moment"]]] + [fmt(r[k]) for k in ("target", "model", "gap", "weight", "loss_contribution")]) + " |" for r in fits]
    lines += ["", "Four ACS targets pool 2012–2023 while the current model observer uses 2023; two parent/child group definitions remain unresolved. The approved childbirth event-study target is unchanged. See overnight_target_mapping_review.md for authoritative sources and dates.", "",
              "## Complete parameter and restriction table", "",
              "The eleven free coordinates below are inherited inputs to this price solve. No parameter has been re-estimated in this run. Near-bound flags reproduce the existing parameter receipt.", "",
              "| Parameter | Value | Lower | Upper | Free coordinate | Near bound | Restriction/status |",
              "|---|---:|---:|---:|---|---|---|"]
    lines += ["| " + " | ".join([r["parameter"], fmt(r["value"]), fmt(r["lower_bound"]), fmt(r["upper_bound"]), r["is_free_parameter"], r["near_bound"], r["status"].replace("_", " ")]) + " |" for r in parameters]
    tail = best["terminal_distance"]
    lines += ["", "## Terminal-distance checks for this trial", "",
              "These endpoint distances do not replace a comparison of historical prices and moments across longer horizons.", "",
              "| Distance | Value | Tolerance | Pass |", "|---|---:|---:|---|"]
    lines += [f"| {k.replace('_', ' ')} | {v:.8g} | {tail['tolerances'][k]:.8g} | {tail['checks'][k]} |" for k, v in tail["metrics"].items()]
    lines += ["", "No new matched policy path is available. The main intended property-tax comparison keeps equal household rebates in both tax regimes. Market convergence, replay, horizon stability, empirical alignment and matched re-estimation remain separate requirements.", "",
              f"Source receipts: `{folder}`. Regenerate with `python3 {Path(__file__).name}` after collecting complete trials. No model solve is performed.", ""]
    (BASE / "HORIZON100_PROGRESS.md").write_text("\n".join(lines))
    print(f"Verified {len(trials)} trials; best market gap {best['maximum_market_residual']:.9g}. Rebuilt full tables.")


if __name__ == "__main__":
    main()
