from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

TOOLS = Path(__file__).resolve().parent
REPOSITORY = TOOLS.parents[2]
NCHS_SOURCE = REPOSITORY / "code" / "data" / "nchs_natality_timing"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import build_e5f_transition_historical_validation as builder


YEARS = builder.CALENDAR_YEARS
AGES = tuple(range(18, 43, 4))
MODEL_AGES = tuple(range(18, 83, 4))
TARGETS = tuple(
    zip(
        builder.ACTIVE_TARGET_SYSTEM.moment_names,
        builder.ACTIVE_TARGET_SYSTEM.target_values,
        builder.ACTIVE_TARGET_SYSTEM.weights,
    )
)
NCHS_FILES = (
    "first_birth_counts_year_age.csv",
    "first_birth_counts_manifest.csv",
    "order_exclusion_shares_by_year.csv",
    "cohort_truncation.csv",
    "timing_targets.csv",
    "timing_targets_exact_age.csv",
    "timing_targets_model_comparable.csv",
    "timing_target_contract.csv",
    "timing_age_bin_counts.csv",
    "build_first_birth_timing.R",
    "test_first_birth_timing_targets.R",
    "timing_target_metadata.json",
)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def target_fingerprint() -> str:
    payload = [[name, target, weight] for name, target, weight in TARGETS]
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def make_fixture(tmp_path: Path) -> argparse.Namespace:
    if target_fingerprint() != builder.EXPECTED_TARGET_FINGERPRINT:
        raise AssertionError("Test target ledger no longer matches the pinned system")
    case = tmp_path / "task" / "cases" / "task_001"
    case.mkdir(parents=True)
    population = {year: 1.0 + 0.03 * period for period, year in enumerate(YEARS)}

    terminal_models: dict[str, float] = {}
    terminal_losses: dict[str, float] = {}
    transition_rows: list[dict[str, object]] = []
    moment_rows: list[dict[str, object]] = []
    fertility_rows: list[dict[str, object]] = []
    fertility_age_rows: list[dict[str, object]] = []
    for period, year in enumerate(YEARS):
        transition_rows.append(
            {
                "scenario": "task_001",
                "period": period,
                "years_from_start": year - 2007,
                "asset_price_index": 1.0 + 0.02 * period,
                "housing_user_cost": 0.1 * (1.0 + 0.02 * period),
                "population_index": population[year],
                "population_target_index": population[year],
                "population_target_gap": 0.0,
                "owner_rate": 0.60 + 0.005 * period,
                "relative_market_residual": 1e-7 * (period + 1),
                "mass_accounting_residual": -1e-14 * period,
                "census_age_bridge_maximum_target_gap": 1e-15 * period,
                "policy_case": "none",
                "policy_active": False,
                "closure": "open_birth_vintage",
            }
        )
        for target_index, (name, target, weight) in enumerate(TARGETS):
            model = (
                5.5 + 0.1 * period
                if name == builder.ROOM_MOMENT
                else target + (4 - period) * 0.001 * (target_index + 1)
            )
            terminal = year == 2023
            if terminal:
                terminal_models[name] = model
                terminal_losses[name] = weight * (model - target) ** 2
            moment_rows.append(
                {
                    "candidate": "task_001",
                    "period": period,
                    "calendar_year": float(year),
                    "moment": name,
                    "model": model,
                    "role": builder.TERMINAL_ROLE if terminal else builder.DIAGNOSTIC_ROLE,
                    "target": target if terminal else "",
                    "weight": weight if terminal else "",
                    "gap": model - target if terminal else "",
                    "loss_contribution": weight * (model - target) ** 2 if terminal else "",
                }
            )
        explicit_rates = [(1.9 - 0.1 * period) / 7.0] * 7
        adjusted_rates = [rate + 0.01 for rate in explicit_rates]
        first_rates = [0.01 * (index + 1) * (1 + 0.02 * period) for index in range(7)]
        fertility_rows.append(
            {
                "period": period,
                "calendar_year": year,
                "period_tfr_explicit_states": sum(explicit_rates),
                "period_tfr_topcode_adjusted": sum(adjusted_rates),
                "period_first_birth_mean_age": 999.0,
                "period_first_birth_share_age30plus": 0.999,
                "role": builder.PERIOD_FERTILITY_ROLE,
                "denominator": "model pre-fertility adult-household mass by age",
            }
        )
        for age_index, age in enumerate(MODEL_AGES):
            age_mass = 0.1 + 0.01 * age_index
            fertile = age_index < len(AGES)
            first_rate = first_rates[age_index] if fertile else 0.0
            explicit_rate = explicit_rates[age_index] if fertile else 0.0
            adjusted_rate = adjusted_rates[age_index] if fertile else 0.0
            fertility_age_rows.append(
                {
                    "period": period,
                    "calendar_year": year,
                    "age_index": age_index,
                    "age_cell_start": age,
                    "age_cell_midpoint": age + 2,
                    "age": age + 2,
                    "age_mass": age_mass,
                    "birth_flow_first": first_rate * age_mass,
                    "birth_flow_explicit": explicit_rate * age_mass,
                    "birth_flow_topcode_adjusted": adjusted_rate * age_mass,
                    "age_specific_birth_rate_explicit": explicit_rate,
                    "age_specific_birth_rate_topcode_adjusted": adjusted_rate,
                }
            )
    write_csv(case / "transition_path.csv", transition_rows)
    write_csv(case / "dated_model_moments.csv", moment_rows)
    write_csv(case / "dated_period_fertility.csv", fertility_rows)
    write_csv(case / "dated_period_fertility_by_age.csv", fertility_age_rows)

    entry = 0.05
    outside_share = 0.169
    mature = entry
    outside = outside_share * entry
    retention = (1.0 - outside_share) * entry / mature
    summary = {
        "status": "complete_transition_calibration_panel_task",
        "source_sha256": "a" * 64,
        "code_fingerprints": {
            "bundle_sha256": "b" * 64,
            "files": {"driver": "c" * 64},
        },
        "target_fingerprint": builder.EXPECTED_TARGET_FINGERPRINT,
        "target_set": builder.EXPECTED_TARGET_SET,
        "target_count": 12,
        "transition_free_parameter_count": 10,
        "target_measurements": {
            "tfr": "fixture cohort TFR",
            "childless_rate": "fixture cohort childlessness",
            "mean_age_first_birth": "fixture cohort timing",
            "share_first_births_age30plus": "fixture cohort timing share",
            "housing_increment_0to1": "fixture event response",
            "remaining_targets": "fixture 2023 cross-section",
        },
        "population_bridge": {
            "status": "externally_normalized_not_estimated",
            "target_indices": {
                str(period): population[year] for period, year in enumerate(YEARS)
            },
        },
        "population_validation_status": "matched by construction",
        "renewal_accounting_contract": {
            "replacement_fertility": 2.1,
            "effective_birth_to_household_conversion": 1.0 / 2.1,
            "birth_to_entry_effect_lag_years": 20.0,
            "birth_to_entry_effect_lag_dates": 5,
            "birth_vintage_queue_waiting_slots": 4,
            "birth_measure": "topcode_adjusted_birth_children",
        },
        "renewal_accounting_old_state": {
            "old_entry_flow_E": entry,
            "old_queue_mature_flow_B": mature,
            "old_queue_B_over_E": mature / entry,
            "outside_flow_M": outside,
            "old_renewal_residual": 0.0,
        },
        "renewal_retention": retention,
        "outside_origin_entry_share": outside_share,
        "old_model_completed_fertility": 2.1,
        "old_completed_fertility_reference": 2.1,
        "old_fertility_normalization": {
            "status": "derived_intercept",
            "target": 2.1,
            "completed_fertility": 2.1,
        },
        "model_profile": {
            "name": builder.EXPECTED_MODEL_PROFILE,
            "profile_id": builder.EXPECTED_MODEL_PROFILE_ID,
            "income_state_count": 15,
            "first_birth_fixed_cost_semantics": (
                "one-time utility cost paid only when the first child arrives"
            ),
        },
        "policy_case": "none",
        "post_2023_periods": 0,
        "best_candidate": {
            "candidate": "task_001",
            "policy_case": "none",
            "post_2023_periods": 0,
            "transition_loss": sum(terminal_losses.values()),
            "terminal_tfr": terminal_models["tfr"],
            "terminal_childless_rate": terminal_models["childless_rate"],
            "terminal_mean_age_first_birth": terminal_models["mean_age_first_birth"],
            "terminal_share_first_births_age30plus": terminal_models[
                "share_first_births_age30plus"
            ],
        },
    }
    write_json(case.parent.parent / "summary.json", summary)

    acs_dir = tmp_path / "acs"
    acs_housing = acs_dir / "national_householder_housing_path.csv"
    housing_rows: list[dict[str, object]] = []
    exposure_rows: list[dict[str, object]] = []
    for period, year in enumerate(YEARS):
        sample = year * 100 + 1
        base_weight = 1000.0 + year
        ownership = 0.62 + 0.005 * period
        rooms = 5.6 + 0.05 * period
        housing_rows.append(
            {
                "calendar_year": year,
                "sample_code": sample,
                "base_head_records": 100,
                "base_head_weight": base_weight,
                "tenure_valid_records": 100,
                "tenure_valid_weight": base_weight,
                "owner_records": 60,
                "owner_weight": ownership * base_weight,
                "ownership_rate": ownership,
                "rooms_valid_records": 100,
                "rooms_valid_weight": base_weight,
                "rooms_weighted_sum": rooms * base_weight,
                "mean_rooms_literal": rooms,
                "ownership_rate_index_2007": ownership / 0.62,
                "mean_rooms_literal_index_2007": rooms / 5.6,
            }
        )
        for age in AGES:
            exposure_rows.append(
                {
                    "calendar_year": year,
                    "sample_code": sample,
                    "age_lower": age,
                    "age_upper": age + 3,
                    "female_records": 100,
                    "female_perwt_mass": 1000,
                    "female_exposure_share_18_45": 1 / 7,
                }
            )
    write_csv(acs_housing, housing_rows)
    acs_exposure = acs_dir / "national_female_exposure_path.csv"
    write_csv(acs_exposure, exposure_rows)
    write_json(
        acs_dir / "metadata.json",
        {
            "contract_sha256": "e" * 64,
            "source": {"sha256": "f" * 64},
            "builder": {"sha256": "1" * 64},
            "output_files": {
                acs_housing.name: {"sha256": digest(acs_housing)},
                acs_exposure.name: {"sha256": digest(acs_exposure)},
            },
        },
    )
    write_json(
        acs_dir / "verification_report.json",
        {
            "status": "pass",
            "canonical_csv_byte_exact": True,
            "canonical_female_exposure_csv_byte_exact": True,
            "canonical_metadata_byte_exact": True,
            "self_checks": [{"check": "fixture", "passed": True}],
        },
    )

    historical_dir = tmp_path / "historical"
    raw_dir = historical_dir / "observed_fred_raw"
    raw_dir.mkdir(parents=True)
    raw_levels = {
        "consumer_price": ("CPIAUCSL", [100 + 10 * p for p in range(5)]),
        "fertility": ("SPDYNTFRTINUSA", [2.1 - 0.1 * p for p in range(5)]),
        "house_price": ("CSUSHPINSA", [100 + 15 * p for p in range(5)]),
        "rent": ("CUUR0000SEHA", [100 + 10 * p for p in range(5)]),
    }
    source_receipts: dict[str, dict[str, object]] = {}
    for source_name, (series_id, levels) in raw_levels.items():
        raw_path = raw_dir / f"{series_id}.csv"
        write_csv(
            raw_path,
            [
                {"observation_date": f"{year}-01-01", series_id: levels[period]}
                for period, year in enumerate(YEARS)
            ],
        )
        source_receipts[source_name] = {
            "id": series_id,
            "raw_file": str(raw_path.resolve()),
            "raw_sha256": digest(raw_path),
            "download_url": f"https://example.test/{series_id}.csv",
            "fred_page": f"https://example.test/{series_id}",
        }
    historical_metadata = historical_dir / "observed_data_metadata.json"
    write_json(
        historical_metadata,
        {
            "aggregation": "arithmetic mean of nonmissing observations within calendar year",
            "series": source_receipts,
        },
    )
    households = {year: 100.0 * population[year] for year in YEARS}
    historical_rows: list[dict[str, object]] = []
    for period, year in enumerate(YEARS):
        cpi = raw_levels["consumer_price"][1][period] / 100.0
        price = raw_levels["house_price"][1][period] / 100.0
        rent = raw_levels["rent"][1][period] / 100.0
        historical_rows.append(
            {
                "calendar_year": year,
                "observed_population_index": population[year],
                "observed_nominal_house_price_index": price,
                "observed_nominal_rent_index": rent,
                "observed_real_house_price_index": price / cpi,
                "observed_real_rent_index": rent / cpi,
                "observed_price_to_rent_index": price / rent,
                "observed_tfr": raw_levels["fertility"][1][period],
            }
        )
    historical_comparison = historical_dir / "historical_comparison.csv"
    write_csv(historical_comparison, historical_rows)
    household_receipt = {
        "title": "Table HH-3 total households, restricted to fixture heads",
        "url": "https://www2.census.gov/fixture/hh3.xls",
        "revision_note": "The revised 2011 row is used.",
        "model_age_totals_thousands": {
            str(year): households[year] for year in YEARS
        },
        "published_hh3_totals_thousands": {
            str(year): households[year] + 5.0 for year in YEARS
        },
    }
    write_json(
        historical_dir / "historical_comparison_summary.json",
        {
            "normalization": "annual averages indexed to 2007=1",
            "number_of_model_dates": 5,
            "start_year": 2007,
            "last_overlapping_year": 2023,
            "population_units": (
                "Census HH-3 households with heads ages 18--85 versus model "
                "adult-household mass"
            ),
            "sources": {**source_receipts, "households": household_receipt},
            "first_row": historical_rows[0],
            "last_row": historical_rows[-1],
        },
    )

    return argparse.Namespace(
        case_dir=case,
        acs_housing_csv=acs_housing,
        acs_female_exposure_csv=acs_exposure,
        historical_comparison_csv=historical_comparison,
        historical_metadata_json=historical_metadata,
        nchs_counts_csv=NCHS_SOURCE / "first_birth_counts_year_age.csv",
        nchs_manifest_csv=NCHS_SOURCE / "first_birth_counts_manifest.csv",
        nchs_metadata_json=NCHS_SOURCE / "timing_target_metadata.json",
        output_dir=tmp_path / "packet",
    )


def copy_nchs_contract(destination: Path) -> None:
    destination.mkdir(parents=True)
    for name in NCHS_FILES:
        shutil.copy2(NCHS_SOURCE / name, destination / name)


def scalar_fields(row: dict[str, str]) -> list[str]:
    return [
        row["gap_2007_model_minus_data"],
        row["gap_2011_model_minus_data"],
        row["gap_2023_model_minus_data"],
        row["change_gap_model_minus_data"],
        row["rmse_five_dates"],
        row["rmse_evaluation_window"],
        row["interval_sign_agreement_fraction"],
        row["interval_sign_agreement_count"],
        row["interval_count"],
    ]


class HistoricalValidationBuilderTests(unittest.TestCase):
    def test_builder_writes_audited_packet(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            result = builder.build_packet(args)
            self.assertTrue(
                result["validation_gates"]["target_count_fingerprint_and_loss_reconciled"]
            )
            self.assertTrue(
                result["validation_gates"]["nchs_metadata_hash_contract_verified"]
            )
            self.assertIn("mixed-scope", result["scope_warning"])
            for stem in (
                "imposed_inputs_matched_by_construction",
                "untargeted_fertility_housing_validation",
                "untargeted_price_rent_validation",
                "numerical_residual_audit",
            ):
                self.assertGreater((args.output_dir / f"{stem}.png").stat().st_size, 0)
                self.assertGreater((args.output_dir / f"{stem}.pdf").stat().st_size, 0)

            comparisons = read_csv(args.output_dir / "five_date_model_data_comparison.csv")
            timing = next(
                row
                for row in comparisons
                if row["series_id"] == "period_first_birth_mean_age"
                and row["calendar_year"] == "2007"
            )
            expected_model_mean = sum(
                midpoint * (index + 1)
                for index, midpoint in enumerate(range(20, 45, 4))
            ) / sum(range(1, 8))
            self.assertAlmostEqual(float(timing["model_value"]), expected_model_mean)
            self.assertAlmostEqual(float(timing["data_value"]), 25.83978425331493)
            timing_share = next(
                row
                for row in comparisons
                if row["series_id"] == "period_first_birth_share_age30plus"
                and row["calendar_year"] == "2007"
            )
            self.assertAlmostEqual(
                float(timing_share["data_value"]), 0.23592839449896857
            )
            rooms = [
                row
                for row in comparisons
                if row["series_id"] == "mean_rooms_literal_18_85"
            ]
            self.assertEqual(
                rooms[0]["classification"],
                "untargeted_holdout_cross_section_only_coding_break",
            )
            self.assertEqual(
                {row["classification"] for row in rooms[1:]},
                {"untargeted_holdout_longitudinal_2011_2023"},
            )

            fit = {
                row["series_id"]: row
                for row in read_csv(args.output_dir / "five_date_fit_statistics.csv")
            }
            for series_id in (
                "period_tfr_explicit_states",
                "period_tfr_topcode_adjusted_sensitivity",
                "price_to_rent_index",
            ):
                self.assertEqual(
                    fit[series_id]["scalar_statistics_status"],
                    "withheld_noncomparable_measurement",
                )
                self.assertEqual(scalar_fields(fit[series_id]), [""] * 9)
            room_fit = fit["mean_rooms_literal_18_85"]
            self.assertEqual(room_fit["evaluation_years"], "2011;2015;2019;2023")
            self.assertEqual(room_fit["n_common_dates"], "4")
            self.assertEqual(room_fit["gap_2007_model_minus_data"], "")
            self.assertEqual(room_fit["change_start_year"], "2011")
            self.assertEqual(room_fit["rmse_five_dates"], "")
            self.assertNotEqual(room_fit["rmse_evaluation_window"], "")

    def test_fails_closed_on_renewal_target_and_loss_corruption(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            summary_path = args.case_dir.parent.parent / "summary.json"
            summary = json.loads(summary_path.read_text())
            summary["renewal_accounting_contract"][
                "effective_birth_to_household_conversion"
            ] = 0.5
            write_json(summary_path, summary)
            with self.assertRaisesRegex(ValueError, "conversion"):
                builder.build_packet(args)

        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            moment_path = args.case_dir / "dated_model_moments.csv"
            rows = read_csv(moment_path)
            for row in rows:
                if row["calendar_year"] == "2023.0" and row["moment"] == "tfr":
                    model = float(row["model"])
                    target = float(row["target"]) + 0.001
                    weight = float(row["weight"])
                    row["target"] = target
                    row["gap"] = model - target
                    row["loss_contribution"] = weight * (model - target) ** 2
            write_csv(moment_path, rows)
            with self.assertRaisesRegex(ValueError, "fingerprint"):
                builder.build_packet(args)

        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            summary_path = args.case_dir.parent.parent / "summary.json"
            summary = json.loads(summary_path.read_text())
            summary["best_candidate"]["transition_loss"] += 1.0
            write_json(summary_path, summary)
            with self.assertRaisesRegex(ValueError, "selected loss"):
                builder.build_packet(args)

    def test_fails_closed_on_historical_receipt_and_population_corruption(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            rows = read_csv(args.historical_comparison_csv)
            rows[-1]["observed_tfr"] = float(rows[-1]["observed_tfr"]) + 0.01
            write_csv(args.historical_comparison_csv, rows)
            with self.assertRaisesRegex(ValueError, "does not reproduce raw receipts"):
                builder.build_packet(args)

        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            transition_path = args.case_dir / "transition_path.csv"
            rows = read_csv(transition_path)
            rows[1]["population_index"] = float(rows[1]["population_index"]) + 0.001
            rows[1]["population_target_index"] = float(
                rows[1]["population_target_index"]
            ) + 0.001
            write_csv(transition_path, rows)
            summary_path = args.case_dir.parent.parent / "summary.json"
            summary = json.loads(summary_path.read_text())
            summary["population_bridge"]["target_indices"]["1"] += 0.001
            write_json(summary_path, summary)
            with self.assertRaisesRegex(ValueError, "observed household population"):
                builder.build_packet(args)

    def test_fails_closed_on_nchs_hash_and_operator_corruption(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = make_fixture(root)
            nchs = root / "nchs"
            copy_nchs_contract(nchs)
            args.nchs_counts_csv = nchs / "first_birth_counts_year_age.csv"
            args.nchs_manifest_csv = nchs / "first_birth_counts_manifest.csv"
            args.nchs_metadata_json = nchs / "timing_target_metadata.json"
            with args.nchs_counts_csv.open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "fingerprint changed"):
                builder.build_packet(args)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = make_fixture(root)
            nchs = root / "nchs"
            copy_nchs_contract(nchs)
            args.nchs_counts_csv = nchs / "first_birth_counts_year_age.csv"
            args.nchs_manifest_csv = nchs / "first_birth_counts_manifest.csv"
            args.nchs_metadata_json = nchs / "timing_target_metadata.json"
            metadata = json.loads(args.nchs_metadata_json.read_text())
            metadata["recommended_operator"]["boundary_treatment"] = "corrupted"
            write_json(args.nchs_metadata_json, metadata)
            with mock.patch.object(
                builder,
                "EXPECTED_NCHS_METADATA_SHA256",
                digest(args.nchs_metadata_json),
            ):
                with self.assertRaisesRegex(ValueError, "boundary treatment"):
                    builder.build_packet(args)

    def test_fails_closed_on_legacy_start_labeled_model_ages(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = make_fixture(Path(directory))
            path = args.case_dir / "dated_period_fertility_by_age.csv"
            rows = read_csv(path)
            for row in rows:
                row.pop("age_cell_start")
                row.pop("age_cell_midpoint")
                row["age"] = float(row["age"]) - 2.0
            write_csv(path, rows)
            with self.assertRaisesRegex(ValueError, "missing required columns"):
                builder.build_packet(args)


if __name__ == "__main__":
    unittest.main()
