#!/usr/bin/env python3
"""Reproduce saved-B late-age fertility accounting and cached NCHS comparison.

This is a read-only receipt builder: it does not run CPS extraction, solve the
model, or access raw NCHS natality files. Run from any working directory.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
OUT = Path(__file__).resolve().parent
B_PATH = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/native_summary.json"
NCHS_PATH = ROOT / "code/data/nchs_natality_timing/first_birth_counts_year_age.csv"
OPUS_PATH = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/phase3_empirical_adjudication/final.md"
OBSERVER_PATH = ROOT / "code/model/tools/e5f_initial_fertility_observer.py"
PERIOD_PATH = ROOT / "code/model/tools/run_e5f_transition_calibration.py"
TOPCODE_PATH = ROOT / "code/model/tools/run_e5f_open_population_transition.py"
NCHS_README_PATH = ROOT / "code/data/nchs_natality_timing/README.md"
TOP = 3.602359422009
YEARS = range(2003, 2007)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def mean(mass: list[float], top: float) -> float:
    total = sum(mass)
    return (mass[1] + 2 * mass[2] + top * mass[3]) / total


def main() -> None:
    summary = json.loads(B_PATH.read_text())
    nchs = list(csv.DictReader(NCHS_PATH.open(newline="")))
    nchs = [r for r in nchs if int(r["year"]) in YEARS]
    nchs_total = sum(int(r["n_first_births"]) for r in nchs)
    nchs_42plus = sum(int(r["n_first_births"]) for r in nchs if int(r["age"]) >= 42)
    nchs_by_age: dict[int, int] = {}
    for r in nchs:
        age = int(r["age"])
        nchs_by_age[age] = nchs_by_age.get(age, 0) + int(r["n_first_births"])

    fertility = summary["early_measurement"]["fertility"]
    rows: list[dict[str, object]] = []
    checks: dict[str, dict[str, float]] = {}
    for variant, obj in fertility.items():
        a = obj["accounting"]
        starts = a["age_cell_start"]
        mids = a["age_cell_midpoint"]
        pre = a["pre_parity_mass_by_age"]
        post = a["post_parity_mass_by_age"]
        flows = a["parity_birth_flows_by_age"]
        max_identity_error = 0.0
        max_topcode_flow_error = 0.0
        for i, start in enumerate(starts):
            p, q, f = pre[i], post[i], flows[i]
            mass = sum(p)
            pre_mean, post_mean = mean(p, TOP), mean(q, TOP)
            adjusted_flow = f[0] + f[1] + (TOP - 2.0) * f[2]
            identity_error = abs((post_mean - pre_mean) - adjusted_flow / mass)
            max_identity_error = max(max_identity_error, identity_error)
            # Pre/post identity also checks parity transition accounting by category.
            for k in range(3):
                inferred = sum(p[: k + 1]) - sum(q[: k + 1])
                max_topcode_flow_error = max(max_topcode_flow_error, abs(inferred - f[k]))
            rows.append({
                "record_type": "model_age_cell",
                "variant": variant,
                "age_cell_start": start,
                "age_cell_midpoint": mids[i],
                "age_cell_interval": f"[{int(start)},{int(start + 4)})",
                "age_cell_population_mass_pre_fertility": mass,
                "pre_mass_n0": p[0], "pre_mass_n1": p[1], "pre_mass_n2": p[2], "pre_mass_n3plus": p[3],
                "post_mass_n0": q[0], "post_mass_n1": q[1], "post_mass_n2": q[2], "post_mass_n3plus": q[3],
                "pre_mean_children_topbin_3_602359422009": pre_mean,
                "post_mean_children_topbin_3_602359422009": post_mean,
                "same_cell_stock_mean_change_topbin_adjusted": post_mean - pre_mean,
                "first_birth_flow": f[0], "second_birth_flow": f[1],
                "third_birth_entry_flow": f[2], "explicit_birth_flow_sum_first_second_third": sum(f),
                "topcode_adjusted_birth_child_flow": adjusted_flow,
                "first_birth_share_of_all_model_first_birth_flow": None,
                "first_birth_flow_rate_per_age_cell_mass": f[0] / mass if mass else 0.0,
                "third_birth_entry_share_of_age_cell_mass": f[2] / mass if mass else 0.0,
                "third_birth_entry_share_of_pre_n0_n1_n2_risk_mass": f[2] / sum(p[:3]) if sum(p[:3]) else 0.0,
                "same_cell_stock_change_identity_error": identity_error,
                "denominator_note": "stationary model household reproductive-member mass; flows and masses are model-weighted, not CPS/NCHS counts",
            })
        first_total = float(a["first_birth_flow"])
        first_42plus = sum(flows[i][0] for i, start in enumerate(starts) if start >= 42)
        explicit_total = float(a["explicit_birth_flow"])
        flow_total = sum(sum(f) for f in flows)
        checks[variant] = {
            "sum_first_flow_vector_minus_saved_total": sum(f[0] for f in flows) - first_total,
            "sum_explicit_flow_vector_minus_saved_total": flow_total - explicit_total,
            "max_same_cell_topbin_stock_flow_identity_abs_error": max_identity_error,
            "max_parity_transition_stock_flow_identity_abs_error": max_topcode_flow_error,
            "first_birth_share_age42plus": first_42plus / first_total,
            "first_birth_flow_age42plus": first_42plus,
            "third_birth_entry_share_cell42_mass": flows[6][2] / sum(pre[6]),
            "age42_cell_topbin_mean_pre": mean(pre[6], TOP),
            "age42_cell_topbin_mean_post": mean(post[6], TOP),
            "age42_cell_topbin_mean_change": mean(post[6], TOP) - mean(pre[6], TOP),
        }

    # Share denominator is all model first-birth flow. Use explicit labels so
    # this cannot be mistaken for a completed-fertility cross-age stock.
    for row in rows:
        if row["record_type"] == "model_age_cell":
            variant = str(row["variant"])
            total = float(fertility[variant]["accounting"]["first_birth_flow"])
            row["first_birth_share_of_all_model_first_birth_flow"] = float(row["first_birth_flow"]) / total if total else 0.0

    nchs_rows = []
    for age in sorted(nchs_by_age):
        nchs_rows.append({"record_type": "nchs_exact_age", "year_start": 2003, "year_end": 2006,
                          "exact_age": age, "model_cell_start_after_boundary_collapse": 18 if age < 22 else (42 if age >= 42 else 18 + 4 * ((age - 18) // 4)),
                          "first_birth_count": nchs_by_age[age], "share_of_all_nchs_first_births_ages12_49": nchs_by_age[age] / nchs_total})
    nchs_rows.append({"record_type": "nchs_collapsed_42plus", "year_start": 2003, "year_end": 2006,
                      "exact_age": "42+", "model_cell_start_after_boundary_collapse": 42,
                      "first_birth_count": nchs_42plus, "share_of_all_nchs_first_births_ages12_49": nchs_42plus / nchs_total})

    fieldnames = list(rows[0].keys())
    with (OUT / "late_fertility_by_age.csv").open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    with (OUT / "late_fertility_nchs_first_births.csv").open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=list(nchs_rows[0].keys()))
        w.writeheader(); w.writerows(nchs_rows)

    receipt = {
        "status": "complete_saved_artifact_diagnostic",
        "method": "Read existing B native_summary and existing pooled NCHS first-birth age cache; no checkpoint load, model solve, CPS rerun, raw NCHS access, or download.",
        "model_population_and_phase": fertility["uniform_birth_time"]["metadata"]["population_approximation"] + "; " + fertility["uniform_birth_time"]["metadata"]["stock_phase"],
        "age_cell_definition": "four-year interval [start,start+4); last fertile cell [42,46), represented by midpoint 44",
        "topbin_value_children": TOP,
        "nchs_period": "2003-2006 pooled, exact ages as stored, first live birth only",
        "nchs_first_birth_counts_total_ages12_49": nchs_total,
        "nchs_first_birth_counts_age42plus": nchs_42plus,
        "nchs_share_age42plus": nchs_42plus / nchs_total,
        "model_variants": checks,
        "nchs_all_birth_orders": {"available": False, "reason": "The checked-in cache is first-birth-only; its builder/README defines order==1. It cannot establish all-birth flow at ages 42-45."},
        "all_order_future_source": "If needed, reconstruct from archived CDC/NCHS Natality public-use detail files for 2003-2006, retaining every live-birth-order record and age of mother (NCHS README documents mager and lbo_rec for these years); aggregate all known-order births by maternal age using an explicitly stated unknown-order policy. The local project cache has no such all-order-by-age table.",
        "interpretation_limits": [
            "Pre/post are the same evaluation's within-cell fertility-choice distributions; their difference is validated against contemporaneous modeled transition flows. They are not adjacent-age stocks or a cohort trajectory.",
            "The model transition vector explicitly counts first, second, and third births entering 3+. The top-code-adjusted child flow weights third-bin entry by T-2=1.602359422009; because the model has no separate 4+ state/flow, this is a top-coded measurement convention, not observed fourth-or-later births.",
            "The 42+ NCHS denominator is all cached first births at ages 12-49 in 2003-2006; the model denominator is total stationary model first-birth flow. Similar age labels do not make the populations, period, or selection mechanism identical.",
            "The 40-44 CPS stock projection is a separate observer statistic. The terminal-cell flow is not a same-population CPS completed-fertility comparison, and no fundamental inconsistency or loss lower bound follows from this diagnostic.",
        ],
        "primary_source_references": [
            "NCHS Natality public-use files (2003-2006) as locally documented in code/data/nchs_natality_timing/README.md; cache has only first-birth counts.",
            "CDC/NCHS 2003 Natality public-use documentation: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat2003doc.pdf",
        ],
        "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in [B_PATH, NCHS_PATH, OPUS_PATH, OBSERVER_PATH, PERIOD_PATH, TOPCODE_PATH, NCHS_README_PATH]},
    }
    (OUT / "late_fertility_diagnostic.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"out": str(OUT), "nchs_42plus_share": receipt["nchs_share_age42plus"], "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
