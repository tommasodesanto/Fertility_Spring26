"""Read-only reduction of the two saved fixed-price credit solutions."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
OUT = ROOT / "output/model/experiments/birth_count_choice/credit_at_binary_winner_v1"
sys.path[:0] = [str(HERE), str(ROOT / "code/model/tools")]

from run_cap2_at_binary_winner import inputs_and_receipt
from model.inputs import load_inputs
from model.estate_contract import apply_experiment_flags, experiment_flags
from model.reporting import build_context
from model.credit import bind_engine_credit
from model.engine import solver
from audit_e5f_estate_resource_account import signed_accounts
from e5f_overnight_estate_audit import policy_mass_branches


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")


def reconstruct(label, phi, source, point, h0):
    arm = OUT / label
    source_arrays = arm / "solution_arrays.npz"
    solve_receipt = json.loads((arm / "solve_completed.json").read_text())
    if sha(source_arrays) != solve_receipt["arrays_sha256"]:
        raise RuntimeError(f"{label}: saved array hash drift")
    P, grid = load_inputs(parameters=point, external_inputs={"H0": [h0], "phi": [phi] * 4})
    apply_experiment_flags(P, experiment_flags(1))
    report = Path(tempfile.mkdtemp(prefix="read_only_reconstruction_", dir=arm))
    context = build_context(P, grid, report, price_start=source["selected"]["price"],
                            deadline=time.time() + 300, max_lifecycle=32, closure="fixed_h0")
    Q = bind_engine_credit(P, "corrected", float(P.unsecured_credit_limit))
    shared = solver.precompute_shared(Q, grid)
    with np.load(source_arrays, allow_pickle=False) as archive:
        solution = SimpleNamespace(**{key: archive[key].copy() for key in archive.files
                                      if not key.startswith("shared.")})
    Q._fert2_probs = solution.fert2_probs.copy()
    Q.birth_count_realized_probs = solution.birth_count_realized_probs.copy()
    Q.birth_count_action_probs = solution.birth_count_action_probs.copy()
    cal = context["prepared"].rt["primitive"].pf.calendar
    price = np.asarray([source["selected"]["price"]])
    policy = cal.policy_from_solution(solution, price, Q, grid, shared)
    pre, checks = cal.reconstruct_stationary_pre_fertility(solution, policy, Q, grid, shared)
    if checks["stationary_post_fertility_nesting_l1"] > 5e-9 or checks["stationary_feasibility_projection_mass"] != 0:
        raise RuntimeError(f"{label}: saved distribution reconstruction fails inherited numerical gate")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(Q.H0[0] * (Q.user_cost_rate * price[0] / Q.r_bar[0]) ** Q.xi_supply[0]),
        float(Q.xi_supply[0]))
    evaluation = cal.evaluate_period(price, pre, Q, grid, shared, cal.SolveCounter(),
                                     supply_rule=supply, supplied_policy=policy)
    ledger = context["prepared"].estate.audit(evaluation, Q, grid)
    native_gate_ledger_equal = None
    if label == "phi_080" or OUT.name.endswith("_v2"):
        saved = json.loads((arm / "reporting/phase_b_ge" / label / "gates.json").read_text())["estate"]
        native_gate_ledger_equal = ledger == saved
        if not native_gate_ledger_equal:
            raise RuntimeError(f"{label}: read-only ledger differs from native PE gate ledger")
    death = np.r_[1 - np.asarray(Q.survival_probs), 1.]
    branches = policy_mass_branches(evaluation, Q)
    tenure = []
    for ten in range(1 + len(Q.H_own)):
        house_value = 0. if ten == 0 else float(price[0] * Q.H_own[ten - 1])
        accounts = [signed_accounts(mass[:, ten:ten + 1], saving[:, ten:ten + 1],
                                    death, np.asarray([house_value]), float(Q.psi))
                    for mass, saving, _ in branches]
        tenure.append(dict(tenure_index=ten, rooms=0. if ten == 0 else float(Q.H_own[ten - 1]),
            net_negative=sum(a["totals"]["net_negative"] for a in accounts),
            negative_estate_death_mass=sum(a["totals"]["negative_estate_death_mass"] for a in accounts),
            net_positive=sum(a["totals"]["net_positive"] for a in accounts)))
    totals = ledger["estate"]["totals"]
    for key in ("net_negative", "negative_estate_death_mass", "net_positive"):
        if abs(sum(t[key] for t in tenure) - totals[key]) > 1e-12:
            raise RuntimeError(f"{label}: tenure estate split differs from audited ledger: {key}")
    write(arm / "estate_from_saved_arrays.json", dict(
        status="read_only_reconstruction_matches_native" if native_gate_ledger_equal else "read_only_reconstruction_unaccepted_arm",
        native_gate_ledger_exact_match=native_gate_ledger_equal, source_arrays_sha256=sha(source_arrays),
        reconstruction=checks, ledger=ledger, by_tenure=tenure,
        gate_net_negative_limit=1e-10,
        production_gate_passes=totals["net_negative"] <= 1e-10))
    pre_counts = np.asarray(solution.birth_count_pre_distribution).sum(axis=(0, 1, 2, 4, 6))
    post_counts = np.asarray(solution.birth_count_post_distribution).sum(axis=(0, 1, 2, 4, 6))
    if pre_counts.shape != (int(Q.J), 4) or post_counts.shape != pre_counts.shape:
        raise RuntimeError("Saved children-ever-born distribution shape drift")
    shares_pre = pre_counts / pre_counts.sum(axis=1, keepdims=True)
    shares_post = post_counts / post_counts.sum(axis=1, keepdims=True)
    age25 = .125 * shares_pre[1] + .875 * shares_post[1]
    if abs(age25.sum() - 1.) > 2e-14:
        raise RuntimeError("Age-25 projected count shares do not sum to one")
    if label == "phi_080":
        observed = json.loads((arm / "reporting/phase_b_ge/phi_080/observers.json").read_text())
        target = observed["fertility"]["uniform_birth_time"]["moments"]["mean_children_ever_born_capped3_age25"]
        if abs(float(age25 @ np.arange(4)) - float(target)) > 1e-12:
            raise RuntimeError("Age-25 projection differs from accepted fertility observer")
    aggregate = dict(label=label, phi=phi, source_arrays_sha256=sha(source_arrays),
        age25_shares=age25.tolist(), age25_cdf=np.cumsum(age25).tolist(),
        age25_motherhood=float(1 - age25[0]), age25_mean_children_ever_born=float(age25 @ np.arange(4)),
        age25_children_conditional_on_motherhood=float(age25 @ np.arange(4) / (1 - age25[0])),
        ownership_all_households=float(solution.g[:, 1:].sum() / solution.g.sum()),
        ownership_age22_25=float(solution.own_by_age[1]),
        ownership_age26_29=float(solution.own_by_age[2]),
        expected_births=float(np.asarray(solution.birth_count_expected_children_by_age).sum()),
        births_by_order=np.asarray(solution.birth_count_births_by_order_by_age).sum(axis=1).tolist(),
        estate_negative_period=totals["net_negative"],
        estate_negative_death_mass=totals["negative_estate_death_mass"],
        estate_total_death_mass=totals["death_mass"],
        accepted_native_pe_gates=label == "phi_080" or OUT.name.endswith("_v2"),
        full_target_report=label == "phi_080")
    age_rows = []
    for j in range(int(Q.J)):
        for timing, shares in (("pre_birth", shares_pre[j]), ("post_birth", shares_post[j])):
            age_rows.append(dict(arm=label, age_start=float(Q.age_start + j * Q.da), timing=timing,
                **{f"share_n_{n if n < 3 else '3plus'}": float(shares[n]) for n in range(4)}))
    return aggregate, age_rows


def main():
    global OUT
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-version", choices=("v1", "v2"), default="v1")
    args = parser.parse_args()
    OUT = ROOT / f"output/model/experiments/birth_count_choice/credit_at_binary_winner_{args.output_version}"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from model import calibration
    source, point, h0, _, _, target_pin, weight_pin = inputs_and_receipt()
    prepared = json.loads((OUT / "prepared.json").read_text())
    if prepared["target_fingerprint"] != target_pin or prepared["weight_fingerprint"] != weight_pin:
        raise RuntimeError("Prepared target/weight identity drift")
    arms = {}
    rows = []
    for label, phi in (("phi_080", .8), ("phi_100", 1.)):
        arms[label], age_rows = reconstruct(label, phi, source, point, h0)
        rows.extend(age_rows)
    with (OUT / "children_ever_born_by_age.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    base, relaxed = arms["phi_080"], arms["phi_100"]
    values = np.arange(4)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for arm, color in ((base, "#2864a4"), (relaxed, "#bc6442")):
        label = f"financed share {arm['phi']:.1f}"
        if arm["phi"] == 1. and args.output_version == "v1":
            label += " (provisional; estate check failed)"
        elif arm["phi"] == 1.:
            label += " (PE gates passed; fixed price)"
        ax.step(values, arm["age25_cdf"], where="post", marker="o", color=color,
                label=label)
    ax.set(xlabel="children ever born by completed interview age 25", ylabel="cumulative household share",
           xticks=values, ylim=(0, 1.03), title="Age-25 children ever born, fixed price")
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(OUT / "age25_children_cdf.png", dpi=180); plt.close(fig)
    result = dict(status="read_only_postprocess_complete", source_search_sha256=prepared["source_search_sha256"],
        postprocessor_sha256=sha(__file__), source_code_pins=prepared["source_sha256"],
        price=source["selected"]["price"], fixed_H0=h0,
        age25_projection="0.125 pre-birth + 0.875 post-birth at completed interview age 25",
        arms=arms,
        change=dict(expected_births_percent=100 * (relaxed["expected_births"] / base["expected_births"] - 1),
            ownership_percentage_points=100 * (relaxed["ownership_all_households"] - base["ownership_all_households"]),
            age25_motherhood_percentage_points=100 * (relaxed["age25_motherhood"] - base["age25_motherhood"]),
            age25_children_conditional_on_motherhood=relaxed["age25_children_conditional_on_motherhood"] - base["age25_children_conditional_on_motherhood"]),
        relaxed_acceptance=("native_PE_gates_passed; fixed_price_not_GE" if args.output_version == "v2"
            else "failed_negative_estate_gate; diagnostics provisional, not accepted counterfactual"),
        no_new_lifecycle_solve=True)
    write(OUT / "diagnostic_summary.json", result)
    print(json.dumps(dict(status=result["status"], change=result["change"],
                          relaxed_net_negative=relaxed["estate_negative_period"]), indent=2))


if __name__ == "__main__":
    main()
