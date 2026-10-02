"""Assemble results.csv for the selling-cost dose probe (read-only; no solves)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARMS = ["psi06", "psi04", "psi02", "psi00"]
LABELS = []
for a in ARMS:
    LABELS += [a + "_phi08", a + "_phi10"]


def model_of(target_fit, moment):
    for row in target_fit:
        if row["moment"] == moment:
            return float(row["model"])
    raise KeyError(moment)


def main():
    results = json.loads((HERE / "case_results.json").read_text())
    rows = {}

    def put(name, vals):
        rows[name] = vals

    for lab in LABELS:
        r = results[lab]
        assert r.get("status") == "passed", lab
        rep = Path(r["report"])
        tf = r["target_fit"]
        ex = r["extra"]
        young = list(csv.DictReader(
            (rep / "young_ownership_age_measurement.csv").open()))[0]
        bf = ex["birth_flows_uniform_birth_time"]
        cap = ex["renter_cap_by_children_home"]
        ownm = ex["ownership_by_children_home"]
        ent = ex["entrant_cell_j0"]
        vals = {
            "first_birth_flow": bf["first_birth_flow"],
            "second_birth_flow": bf["second"],
            "third_birth_flow": bf["third_bin"],
            "completed_fertility_TFR": r["closure"]["completed_fertility_tfr"],
            "childlessness_40_44": model_of(tf, "cps_childlessness"),
            "ownership_30_55": model_of(tf, "ownership_30_55"),
            "own_rate_25_34": float(young["own_rate_2534"]),
            "first_birth_rooms": model_of(tf, "first_birth_rooms"),
            "renter_cap_share_m1": cap["1"]["cap_share_of_renters"],
            "renter_cap_share_m2": cap["2"]["cap_share_of_renters"],
            "renter_cap_share_m3": cap["3"]["cap_share_of_renters"],
            "entrant_first_birth_prob": ent["first_birth_probability"],
            "entrant_ownership": ent["ownership_rate"],
            "gates_status": r["closure"]["gates_status"],
            "renewal_residual": r["closure"]["renewal_residual"],
            "housing_residual": r["closure"]["absolute_housing_residual"],
            "paygo_residual": r["closure"]["actual_paygo_residual"],
            "lifecycle_solve_seconds": r["lifecycle_solve_seconds"],
        }
        for k, v in vals.items():
            rows.setdefault(k, {})[lab] = v
    order = list(rows)
    fin_cols = []
    for a in ARMS:
        fin_cols.append("fin_effect_%s(phi10-phi08)" % a)
    with (HERE / "results.csv").open("w", newline="") as s:
        w = csv.writer(s)
        w.writerow(["row", *LABELS, *fin_cols])
        for k in order:
            v = rows[k]
            if k == "gates_status":
                w.writerow([k, *[v[lab] for lab in LABELS],
                            *["" for _ in fin_cols]])
                continue
            eff = {a: v[a + "_phi10"] - v[a + "_phi08"] for a in ARMS}
            w.writerow([k, *[v[lab] for lab in LABELS],
                        *[eff[a] for a in ARMS]])
    # params_used.json: exact case parameter contracts.
    ref = json.loads((HERE / "reference_point.json").read_text()) \
        if (HERE / "reference_point.json").exists() else {"point": {}}
    params_used = {}
    for lab in LABELS:
        r = results[lab]
        params_used[lab] = dict(
            free_coordinates=ref.get("point", {}), price=r["price"], H0=r["H0"],
            financed_share_phi=r["phi"], phi_vector=[r["phi"]] * 4,
            hbar_first_child_jump=r["hbar_first_child_jump"],
            hbar_child_rooms=r["hbar_child_rooms"],
            selling_cost_psi=r.get("selling_cost_psi"),
            tenure_choice_kappa=r.get("tenure_choice_kappa"),
            H_own=r.get("H_own"), n_house=r.get("n_house"),
            purchase_saving_fraction=r["purchase_saving_fraction"],
            unsecured_credit_d_bar=0.0, entry_arm="nonnegative_mean",
            wealth_grid_nodes=120, income_states=9,
            hbar_reporting_accommodation=r.get("hbar_reporting_accommodation", False),
            report_dir=r["report"])
    (HERE / "params_used.json").write_text(
        json.dumps(params_used, indent=2, sort_keys=True) + "\n")
    print("wrote results.csv with", len(order), "rows")


if __name__ == "__main__":
    main()
