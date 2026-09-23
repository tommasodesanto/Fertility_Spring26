#!/usr/bin/env python3
"""Reproduce the saved CPS fertility summaries and NCHS timing convention.

Read-only with respect to all source data, model builders, and calibration
targets. CPS fixed-width records are streamed from the local IPUMS gzip file.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[6]
CPS_GZ = Path("/Users/tommasodesanto/Desktop/Projects/Datasets/CPS/extract3/cps_00003.dat.gz")
LOADER = CPS_GZ.with_name("loader.do")
NCHS = ROOT / "code/data/nchs_natality_timing/first_birth_counts_year_age.csv"
MODEL_B = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/native_summary.json"
MODEL_B_INITIAL_CONTRACT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/initial_contract.json"
MODEL_B_EVALUATION_RECEIPT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/verified_evaluation_receipt.json"
FERTILITY_RECEIPT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/fertility_receipt.json"
REFERENCE_CPS_BUILDER = ROOT / "output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/extract_fertility_availability.py"
OBSERVER_SOURCE = ROOT / "code/model/tools/e5f_initial_fertility_observer.py"
NORMALIZATION_SOURCE = ROOT / "code/model/tools/run_e5f_transition_calibration.py"
CALIBRATION_SOURCE = ROOT / "code/model/intergen_eqscale_seq_optimized/calibration.py"
SOLVER_SOURCE = ROOT / "code/model/intergen_eqscale_seq_optimized/solver.py"
CHAIN_SOURCE = ROOT / "code/model/intergen_eqscale_seq_optimized/run_e1_chain.py"
TRANSITION_SOURCE = ROOT / "code/model/tools/run_e5f_open_population_transition.py"
AUDIT_CHAIN_SOURCE = ROOT / "code/model/tools/audit_closed_reproductive_closure.py"
TOP_BIN = 3.602359422009
WIDTH = 261


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cps_records() -> tuple[dict[int, list[tuple[int, float]]], dict, str]:
    rows: dict[int, list[tuple[int, float]]] = {2004: [], 2006: []}
    partition_hash = {year: hashlib.sha256() for year in rows}
    count_codes = {year: Counter() for year in rows}
    june_rows = Counter()
    last_key = (0, 0)
    with gzip.open(CPS_GZ, "rb") as src:
        while True:
            rec = src.read(WIDTH)
            if not rec:
                break
            if len(rec) != WIDTH or not rec.endswith(b"\n"):
                raise ValueError("Malformed fixed-width record or record missing terminal newline")
            year, month = int(rec[0:4]), int(rec[9:11])
            key = (year, month)
            if key < last_key:
                raise ValueError("CPS records are not sorted by year/month; partition identity is unsafe")
            last_key = key
            if year not in rows or month != 6:
                continue
            partition_hash[year].update(rec)
            june_rows[year] += 1
            sex = int(rec[148:149])
            age = int(rec[146:148])
            if sex != 2 or not 40 <= age <= 44:
                continue
            children = int(rec[238:241])
            weight = int(rec[250:260]) / 10000.0
            count_codes[year][children] += 1
            if 0 <= children <= 20 and weight > 0:
                rows[year].append((children, weight))
    audit = {
        str(year): {
            "june_partition_rows": june_rows[year],
            "june_partition_sha256_uncompressed": partition_hash[year].hexdigest(),
            "female_age_40_44_frever_code_counts": dict(count_codes[year]),
            "valid_weighted_sample_n": len(rows[year]),
        }
        for year in rows
    }
    return rows, audit, sha256(CPS_GZ)


def cps_moments(label: str, data: list[tuple[int, float]]) -> dict:
    sw = sum(w for _, w in data)
    weighted_counts = {k: sum(w for n, w in data if n == k) for k in (0, 1, 2)}
    weighted_counts["3+"] = sum(w for n, w in data if n >= 3)
    shares = {str(k): v / sw for k, v in weighted_counts.items()}
    p0, p1, p2, p3 = (shares["0"], shares["1"], shares["2"], shares["3+"])
    mothers_weight = sw - weighted_counts[0]
    exact_one_weight = weighted_counts[1]
    top_weight = weighted_counts["3+"]
    m_uncap = sum(n * w for n, w in data) / sw
    m_cap3 = sum(min(n, 3) * w for n, w in data) / sw
    m_cap5 = sum(min(n, 5) * w for n, w in data) / sw
    mu3_uncap = sum(n * w for n, w in data if n >= 3) / top_weight
    mu3_cap5 = sum(min(n, 5) * w for n, w in data if n >= 3) / top_weight
    # The target-review contract's 3.602 is a top-bin value applied to this
    # same sample's literal weighted 0/1/2/3+ shares; it is a diagnostic.
    m_topbin_3602 = p1 + 2 * p2 + TOP_BIN * p3
    s1_mothers = exact_one_weight / mothers_weight
    mu2plus_observed = (sum(n * w for n, w in data if n >= 2) / sw) / (p2 + p3)
    mu2plus_required_for_2p1 = (2.1 / (1 - p0) - s1_mothers) / (1 - s1_mothers)
    p3_fraction_among_2plus = p3 / (p2 + p3)
    required_p3_fraction_under_T = (mu2plus_required_for_2p1 - 2) / (TOP_BIN - 2)
    return {
        "window": label,
        "n_unweighted": len(data),
        "sum_supplement_weights": sw,
        "weighted_shares_0_1_2_3plus": shares,
        "weighted_childlessness_share": p0,
        "weighted_exactly_one_among_mothers": s1_mothers,
        "weighted_mean_children_ever_born_uncapped": m_uncap,
        "weighted_mean_children_ever_born_capped_at_3": m_cap3,
        "weighted_mean_children_ever_born_capped_at_5": m_cap5,
        "weighted_3plus_conditional_mean_uncapped": mu3_uncap,
        "weighted_3plus_conditional_mean_capped_at_5": mu3_cap5,
        "weighted_mean_with_3plus_replaced_by_3_602": m_topbin_3602,
        "conditional_2plus_mean_observed_same_population": mu2plus_observed,
        "3plus_share_among_2plus_observed_same_population": p3_fraction_among_2plus,
        "same_population_identity_diagnostic_only": {
            "hypothetical_completed_mean": 2.1,
            "required_mean_among_2plus": mu2plus_required_for_2p1,
            "assumed_3plus_group_value": TOP_BIN,
            "required_3plus_share_among_2plus_if_other_groups_are_exactly_2_and_3_602": required_p3_fraction_under_T,
            "interpretation": "Arithmetic applied only to this single CPS population. It does not compare the age-40-44 stock observer to the model's separate completed-fertility normalization."
        },
        "weighted_denominators": {
            "all_valid_sample": sw,
            "mothers": mothers_weight,
            "exactly_one": exact_one_weight,
            "3plus": top_weight,
        },
    }


def nchs_moments() -> tuple[list[dict], dict, str]:
    with NCHS.open(newline="") as f:
        records = list(csv.DictReader(f))
    counts: dict[int, Counter] = {year: Counter() for year in range(2003, 2007)}
    for r in records:
        y = int(r["year"])
        if y in counts:
            counts[y][int(r["age"])] += int(r["n_first_births"])

    def cell_midpoint(age: int) -> int:
        if age < 22:
            return 20
        if age >= 42:
            return 44
        return 20 + 4 * ((age - 18) // 4)

    def mapping_by_bracket(age_counts: Counter, total: int) -> dict:
        brackets = [
            ("under_18", lambda age: age < 18),
            ("18_21", lambda age: 18 <= age <= 21),
            ("22_41", lambda age: 22 <= age <= 41),
            ("42_45", lambda age: 42 <= age <= 45),
            ("over_45", lambda age: age > 45),
        ]
        parts = {}
        for name, contains in brackets:
            n_births = sum(n for age, n in age_counts.items() if contains(age))
            delta = sum((cell_midpoint(age) - (age + 0.5)) * n
                        for age, n in age_counts.items() if contains(age)) / total
            parts[name] = {
                "first_birth_count": n_births,
                "share": n_births / total,
                "mapping_contribution_years": delta,
            }
        return parts

    byyear = []
    pooled_counts: Counter = Counter()
    for year, age_counts in counts.items():
        pooled_counts.update(age_counts)
        total = sum(age_counts.values())
        raw_age = sum(age * n for age, n in age_counts.items()) / total
        exact = sum((age + 0.5) * n for age, n in age_counts.items()) / total
        mapped = sum(cell_midpoint(age) * n for age, n in age_counts.items()) / total
        parts = mapping_by_bracket(age_counts, total)
        byyear.append({
            "year": year,
            "first_birth_count": total,
            "mean_recorded_single_age": raw_age,
            "mean_single_age_interval_midpoint_age_plus_0_5": exact,
            "mean_exact_single_age_midpoint_age_plus_0_5": exact,
            "age_plus_half_is_single_year_interval_midpoint_assumption": True,
            "mean_current_cell_midpoint_mapping": mapped,
            "single_age_bin_midpoint_effect_years": exact - raw_age,
            "mapping_effect_beyond_exact_age_midpoint_years": mapped - exact,
            "mapping_effect_beyond_single_age_interval_midpoint_years": mapped - exact,
            "difference_vs_recorded_single_age_total_years": mapped - raw_age,
            "first_bin_total_mapping_contribution_years_ages_below_22": parts["under_18"]["mapping_contribution_years"] + parts["18_21"]["mapping_contribution_years"],
            "mapping_effect_by_age_bracket": parts,
        })
    total = sum(pooled_counts.values())
    raw_age = sum(age * n for age, n in pooled_counts.items()) / total
    exact = sum((age + 0.5) * n for age, n in pooled_counts.items()) / total
    mapped = sum(cell_midpoint(age) * n for age, n in pooled_counts.items()) / total
    parts = mapping_by_bracket(pooled_counts, total)
    out = {
        "period": "2003-2006 pooled annual first-birth counts, ages 12-49",
        "first_birth_count": total,
        "mean_recorded_single_age": raw_age,
        "mean_single_age_interval_midpoint_age_plus_0_5": exact,
        "mean_exact_single_age_midpoint_age_plus_0_5": exact,
        "age_plus_half_is_single_year_interval_midpoint_assumption_not_exact_continuous_age": True,
        "mean_current_cell_midpoint_mapping": mapped,
        "single_age_bin_midpoint_effect_years": exact - raw_age,
        "mapping_effect_beyond_exact_age_midpoint_years": mapped - exact,
        "difference_vs_recorded_single_age_total_years": mapped - raw_age,
        "first_bin_total_mapping_contribution_years_ages_below_22": parts["under_18"]["mapping_contribution_years"] + parts["18_21"]["mapping_contribution_years"],
        "mapping_effect_by_age_bracket": parts,
        "share_first_births_age_30plus_exact_age": sum(n for age, n in pooled_counts.items() if age >= 30) / total,
        "mapping": "age<22 -> 20; 22-25 -> 24; 26-29 -> 28; 30-33 -> 32; 34-37 -> 36; 38-41 -> 40; age>=42 -> 44",
        "operator": "Pooled period first-birth counts; no female exposure denominator; not cohort completed fertility.",
    }
    return byyear, out, sha256(NCHS)


def model_b_readout() -> dict:
    d = json.loads(MODEL_B.read_text())
    initial_contract = json.loads(MODEL_B_INITIAL_CONTRACT.read_text())
    evaluation_receipt = json.loads(MODEL_B_EVALUATION_RECEIPT.read_text())
    fertility = d["early_measurement"]["fertility"]
    norm = d["normalization"]
    uniform_accounting = fertility["uniform_birth_time"]["accounting"]
    age_starts = uniform_accounting["age_cell_start"]
    norm_index = age_starts.index(42.0)
    norm_post_mass = uniform_accounting["post_parity_mass_by_age"][norm_index]
    norm_total = sum(norm_post_mass)
    observer_g_post_shares = [mass / norm_total for mass in norm_post_mass]
    post_mass_by_age = uniform_accounting["post_parity_mass_by_age"]
    age46plus_mass = [
        sum(post_mass_by_age[j][k] for j, start_age in enumerate(age_starts) if start_age >= 46.0)
        for k in range(4)
    ]
    age46plus_total = sum(age46plus_mass)
    observer_g_post_age46plus_shares = [mass / age46plus_total for mass in age46plus_mass]
    legacy = d["legacy_stationary_moments"]
    legacy_shares = [
        legacy["parity_share_0"], legacy["parity_share_1"],
        legacy["parity_share_2plus"] - legacy["parity_share_3plus"],
        legacy["parity_share_3plus"],
    ]
    return {
        "source": str(MODEL_B),
        "source_sha256": sha256(MODEL_B),
        "frozen_run_provenance": {
            "source_commit": initial_contract.get("source_commit"),
            "initial_contract_sha256": sha256(MODEL_B_INITIAL_CONTRACT),
            "evaluation_receipt_sha256": sha256(MODEL_B_EVALUATION_RECEIPT),
            "input_sha256": evaluation_receipt.get("input_sha256"),
            "source_fingerprints": evaluation_receipt.get("source_fingerprints"),
            "normalization_path_frozen_source_sha256": {
                path: value for path, value in initial_contract.get("source_sha256", {}).items()
                if any(term in path for term in [
                    "intergen_eqscale_seq_optimized/calibration.py",
                    "intergen_eqscale_seq_optimized/run_e1_chain.py",
                    "intergen_eqscale_seq_optimized/solver.py",
                    "tools/audit_closed_reproductive_closure.py",
                    "tools/run_e5f_open_population_transition.py",
                    "tools/run_e5f_transition_calibration.py",
                ])
            },
        },
        "normalization": {
            **norm,
            "implementation": "solve_old_steady_state varies psi_child until chain.extract_moments(...)[tfr] is within tolerance of 2.1. The chain is loaded by run_e5f_open_population_transition.configure_sequential_model from audit_closed_reproductive_closure.load_chain, which reloads run_e1_chain; that chain imports intergen_eqscale_seq_optimized.calibration.extract_moments. Under literal_topcode, it values parity_dist with P.tfr_top_bin_weight. The frozen solver computes parity_dist from g starting at P.A_f_end; model ages start at 18 in four-year cells, and A_f_end=7 puts this normalization on post-fertility age cells starting at 46.",
            "normalization_source": str(NORMALIZATION_SOURCE),
            "normalization_source_sha256": sha256(NORMALIZATION_SOURCE),
            "calibration_source": str(CALIBRATION_SOURCE),
            "calibration_source_sha256": sha256(CALIBRATION_SOURCE),
            "solver_source": str(SOLVER_SOURCE),
            "solver_source_sha256": sha256(SOLVER_SOURCE),
            "chain_source": str(CHAIN_SOURCE),
            "chain_source_sha256": sha256(CHAIN_SOURCE),
            "transition_source": str(TRANSITION_SOURCE),
            "transition_source_sha256": sha256(TRANSITION_SOURCE),
            "chain_loader_source": str(AUDIT_CHAIN_SOURCE),
            "chain_loader_source_sha256": sha256(AUDIT_CHAIN_SOURCE),
            "observer_source": str(OBSERVER_SOURCE),
            "observer_source_sha256": sha256(OBSERVER_SOURCE),
        },
        "post_fertility_ages_46plus_normalization_population": {
            "age_cell_starts_included": [age for age in age_starts if age >= 46.0],
            "shares_from_saved_stationary_summary_0_1_2_3plus": legacy_shares,
            "top_bin_value": TOP_BIN,
            "mean_using_top_bin_value_3_602359422009": legacy_shares[1] + 2 * legacy_shares[2] + TOP_BIN * legacy_shares[3],
            "mean_if_top_bin_is_literal_3": legacy_shares[1] + 2 * legacy_shares[2] + 3 * legacy_shares[3],
            "saved_model_normalization_statistic": norm["completed_fertility"],
            "independent_observer_g_post_fertility_shares_ages_46plus_0_1_2_3plus": observer_g_post_age46plus_shares,
            "independent_observer_g_post_fertility_shares_age42_cell_0_1_2_3plus": observer_g_post_shares,
            "interpretation": "The stationary solver forms parity_dist from g over age indices starting at P.A_f_end. With 18-start, four-year periods and A_f_end=7 (fertility through the cell starting at 42), this is the post-fertility 46+ population. Applying the calibrated top-bin value reproduces the 2.1 normalization; literal 3 gives the saved mean_completed_fertility. The independent initial-fertility observer reads g_post_fertility and is retained separately. Neither is the CPS [40,45) age-projected observer.",
        },
        "uniform_birth_time_age_40_44": fertility["uniform_birth_time"]["parity_shares_40_44"],
        "constant_post_cell_age_40_44": fertility["constant_post_cell"]["parity_shares_40_44"],
        "uniform_birth_time_window_parity_mass_age_40_44": uniform_accounting["window_parity_mass"],
        "same_evaluation_full_age_pre_fertility_parity_mass_by_age": uniform_accounting["pre_parity_mass_by_age"],
        "same_evaluation_full_age_post_fertility_parity_mass_by_age": uniform_accounting["post_parity_mass_by_age"],
        "same_evaluation_full_age_post_fertility_parity_shares_by_age": [
            [value / sum(row) for value in row] if sum(row) > 0 else None
            for row in uniform_accounting["post_parity_mass_by_age"]
        ],
        "legacy_stationary_moments": {
            k: d["legacy_stationary_moments"].get(k)
            for k in ["mean_completed_fertility", "tfr", "parity_share_0", "parity_share_1", "parity_share_2plus", "parity_share_3plus"]
        },
        "observer_warning": "The 40-44 shares are projected stock moments; the separately recorded normalization uses the post-fertility age-46+ parity distribution from the stationary solver. The saved initial-fertility observer uses g_pre/g_post_fertility and is also distinct. Do not infer incompatibility or a loss bound across these observers.",
        "fertility_receipt_sha256": sha256(FERTILITY_RECEIPT),
    }


def main() -> None:
    rows, partition_info, cps_source_sha = cps_records()
    cps = [cps_moments(str(y), rows[y]) for y in (2004, 2006)]
    cps.append(cps_moments("2004+2006 pooled weighted records", rows[2004] + rows[2006]))
    nchs_byyear, nchs_pooled, nchs_sha = nchs_moments()
    receipt = {
        "schema": "independent_empirical_fertility_reproduction_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "read_only": True,
        "cps": {
            "source": str(CPS_GZ),
            "compressed_size_bytes": CPS_GZ.stat().st_size,
            "compressed_sha256": cps_source_sha,
            "loader": str(LOADER),
            "loader_sha256": sha256(LOADER),
            "reference_builder": str(REFERENCE_CPS_BUILDER),
            "reference_builder_sha256": sha256(REFERENCE_CPS_BUILDER),
            "sample": "June 2004/2006; female (SEX=2); age 40-44; valid FREVER 0-20; FRSUPPWT>0; exclude 999/NIU; supplied supplement weight divided by 10,000.",
            "record_width_bytes": WIDTH,
            "partitions": partition_info,
            "moments": cps,
        },
        "nchs": {
            "source": str(NCHS),
            "source_sha256": nchs_sha,
            "year_rows": nchs_byyear,
            "pooled": nchs_pooled,
        },
        "model_b_local_receipt": model_b_readout(),
        "scope_limits": [
            "No raw-data standard errors were recomputed.",
            "No target, weight, builder, or model code was edited.",
            "2.1 is shown only as a hypothetical identity applied within one weighted CPS sample; it is kept separate from the model's full-population normalization and age-40-44 projected stock observer.",
            "The current NCHS age mapping is compared with exact single-age midpoints. The midpoint differences are measurement arithmetic, not evidence of an economic target error.",
        ],
    }
    (OUT / "empirical_reproduction.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    with (OUT / "cps_weighted_moments.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["window", "n_unweighted", "sum_supplement_weights", "share_0", "share_1", "share_2", "share_3plus", "childless_share", "exactly_one_among_mothers", "mean_uncapped", "mean_capped3", "mean_capped5", "mean_3plus_uncapped", "mean_3plus_capped5", "mean_topbin_3_602359422009", "observed_3plus_share_among_2plus", "required_3plus_share_among_2plus_if_same_population_mean_2_1"])
        writer.writeheader()
        for r in cps:
            s = r["weighted_shares_0_1_2_3plus"]
            writer.writerow({
                "window": r["window"], "n_unweighted": r["n_unweighted"], "sum_supplement_weights": r["sum_supplement_weights"],
                "share_0": s["0"], "share_1": s["1"], "share_2": s["2"], "share_3plus": s["3+"],
                "childless_share": r["weighted_childlessness_share"], "exactly_one_among_mothers": r["weighted_exactly_one_among_mothers"],
                "mean_uncapped": r["weighted_mean_children_ever_born_uncapped"], "mean_capped3": r["weighted_mean_children_ever_born_capped_at_3"],
                "mean_capped5": r["weighted_mean_children_ever_born_capped_at_5"],
                "mean_3plus_uncapped": r["weighted_3plus_conditional_mean_uncapped"], "mean_3plus_capped5": r["weighted_3plus_conditional_mean_capped_at_5"],
                "mean_topbin_3_602359422009": r["weighted_mean_with_3plus_replaced_by_3_602"],
                "observed_3plus_share_among_2plus": r["3plus_share_among_2plus_observed_same_population"],
                "required_3plus_share_among_2plus_if_same_population_mean_2_1": r["same_population_identity_diagnostic_only"]["required_3plus_share_among_2plus_if_other_groups_are_exactly_2_and_3_602"],
            })
    with (OUT / "nchs_midpoint_comparison.csv").open("w", newline="") as f:
        fields = ["year", "first_birth_count", "mean_recorded_single_age", "mean_single_age_interval_midpoint_age_plus_0_5", "age_plus_half_is_single_year_interval_midpoint_assumption", "mean_current_cell_midpoint_mapping", "single_age_bin_midpoint_effect_years", "mapping_effect_beyond_single_age_interval_midpoint_years", "difference_vs_recorded_single_age_total_years", "first_bin_total_mapping_contribution_years_ages_below_22"]
        for bracket in ("under_18", "18_21", "22_41", "42_45", "over_45"):
            fields.extend([f"{bracket}_birth_count", f"{bracket}_share", f"{bracket}_mapping_contribution_years"])
        fields.append("share_first_births_age_30plus_exact_age")
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in nchs_byyear:
            flat = {k: v for k, v in row.items() if k in fields}
            for bracket, result in row["mapping_effect_by_age_bracket"].items():
                flat[f"{bracket}_birth_count"] = result["first_birth_count"]
                flat[f"{bracket}_share"] = result["share"]
                flat[f"{bracket}_mapping_contribution_years"] = result["mapping_contribution_years"]
            writer.writerow(flat)
    print(json.dumps({"output": str(OUT), "cps": cps, "nchs_pooled": nchs_pooled}, indent=2))


if __name__ == "__main__":
    main()
