"""Assemble results.csv for the per-child need probe (read-only; no solves)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
LABELS = ["need0_phi08", "need0_phi10", "need1_phi08", "need1_phi10"]


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
        rep = Path(r["report"])
        tf = r["target_fit"]
        ex = r["extra"]
        young = list(csv.DictReader(
            (rep / "young_ownership_age_measurement.csv").open()))[0]
        lc = {x["moment"]: x for x in tf}
        bf = ex["birth_flows_uniform_birth_time"]
        cap = ex["renter_cap_by_children_home"]
        ent = ex["entrant_cell_j0"]
        vals = {
            "first_birth_flow": bf["first_birth_flow"],
            "second_birth_flow": bf["second"],
            "third_birth_flow": bf["third_bin"],
            "second_birth_flow_gcheck": ex["birth_flows_g_arrays"]["second"],
            "third_birth_flow_gcheck": ex["birth_flows_g_arrays"]["third_bin"],
            "completed_fertility_TFR": r["closure"]["completed_fertility_tfr"],
            "completed_fertility_targetfit": model_of(tf, "initial_normalization"),
            "childlessness_40_44": model_of(tf, "cps_childlessness"),
            "ownership_30_55": model_of(tf, "ownership_30_55"),
            "own_rate_25_34": float(young["own_rate_2534"]),
            "first_birth_rooms": model_of(tf, "first_birth_rooms"),
            "mean_rooms": model_of(tf, "mean_rooms"),
            "renter_cap_share_m1": cap["1"]["cap_share_of_renters"],
            "renter_cap_share_m2": cap["2"]["cap_share_of_renters"],
            "renter_cap_share_m3": cap["3"]["cap_share_of_renters"],
            "renter_cap_share_m0_memo": cap["0"]["cap_share_of_renters"],
            "entrant_first_birth_prob": ent["first_birth_probability"],
            "entrant_ownership": ent["ownership_rate"],
            "renewal_residual": r["closure"]["renewal_residual"],
            "housing_residual": r["closure"]["absolute_housing_residual"],
            "paygo_residual": r["closure"]["actual_paygo_residual"],
            "lifecycle_solve_seconds": r["lifecycle_solve_seconds"],
        }
        for k, v in vals.items():
            rows.setdefault(k, {})[lab] = v
    order = list(rows)
    with (HERE / "results.csv").open("w", newline="") as s:
        w = csv.writer(s)
        w.writerow(["row", *LABELS, "fin_effect_need0(phi10-phi08)",
                    "fin_effect_need1(phi10-phi08)", "DiD(need1-Need0)"])
        for k in order:
            v = rows[k]
            e0 = v["need0_phi10"] - v["need0_phi08"]
            e1 = v["need1_phi10"] - v["need1_phi08"]
            w.writerow([k, v["need0_phi08"], v["need0_phi10"],
                        v["need1_phi08"], v["need1_phi10"], e0, e1, e1 - e0])
    # params_used.json: exact case parameter contracts.
    ref = json.loads((HERE / "reference_point.json").read_text())
    params_used = {}
    for lab in LABELS:
        r = results[lab]
        params_used[lab] = dict(
            free_coordinates=ref["point"], price=r["price"], H0=r["H0"],
            financed_share_phi=r["phi"], phi_vector=[r["phi"]] * 4,
            hbar_first_child_jump=r["hbar_first_child_jump"],
            hbar_child_rooms=r["hbar_child_rooms"],
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
