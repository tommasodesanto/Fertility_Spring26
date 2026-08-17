from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path

import finalize_first_birth_rooms_target_receipt as receipt


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture(root: Path) -> tuple[Path, Path, Path]:
    output = root / "output"
    output.mkdir(parents=True)
    source = root / "source.dta"
    source.write_bytes(b"source")
    do_file = root / "builder.do"
    do_file.write_text("clear all\n", encoding="utf-8")
    (output / "sa_rooms_first_birth_household_aligned_v1.log").write_text(
        "completed\n", encoding="utf-8"
    )
    event_rows = []
    for event_time in range(-7, 12):
        coefficient_name = (
            f"F{abs(event_time)}event" if event_time < 0 else f"L{event_time}event"
        )
        estimate = -0.05 if event_time < 0 else 0.05
        standard_error = 0.08
        if event_time == -2:
            estimate = 0.0
            standard_error = 0.0
        if event_time == 3:
            estimate = 0.4
            standard_error = 0.1
        event_rows.append(
            {
                "coefficient_name": coefficient_name,
                "relative_time": event_time,
                "b": estimate,
                "se": standard_error,
                "ci_lo": estimate - 1.96 * standard_error,
                "ci_hi": estimate + 1.96 * standard_error,
            }
        )
    write_csv(output / "event_study_estimates.csv", event_rows)
    write_csv(
        output / "target_receipt.csv",
        [
            {
                "moment": "housing_increment_0to1",
                "estimate": 0.45,
                "standard_error": 0.12,
                "contrast_start_time": -1,
                "contrast_end_time": 3,
                "regression_omitted_time": -2,
                "component_l3": 0.4,
                "component_f1": -0.05,
                "covariance_l3_f1": 0.001,
                "input_observations": 120,
                "input_individuals": 30,
                "treated_individuals": 20,
                "never_treated_observations": 40,
                "estimation_observations": 100,
                "estimation_individuals": 25,
                "estimation_household_years": 100,
                "est_never_obs": 35,
                "est_confirmed_never_ids": 10,
                "excl_multi_fu_hhyears": 5,
                "excl_multi_fu_women": 8,
                "excl_untimed_parent_ids": 2,
                "excl_unknown_history_ids": 3,
                "excl_bad_room_gap_rows": 4,
                "pre_event_mean_rooms": 5.0,
                "max_women_before_dedup": 2,
                "runtime_seconds": 10.0,
                "estimator": "eventstudyinteract",
                "sample": "current women age 18+, ref/spouse, positive IW, single-FID dwelling, one woman per HH-year",
                "fixed_effects": "person and survey-year fixed effects; age and education controls",
                "clustering": "individual ID",
                "weighting": "PSID longitudinal pweight IW",
                "control_group": "confirmed zero-child women from full relationship history",
                "rooms_alignment": "ACTUALROOMS_ shifted forward one observed interview within individual",
                "source_file": "PSID/PSIDSHELF_MOBILITY.dta",
                "fertility_timing": "first biological child across RELCHI1-20 TYPE/BYEAR records",
                "single_fu_only": 1,
                "status": "corrected_primary_target",
            }
        ],
    )
    return output, source, do_file


def test_valid_fixture_and_corruption_rejection() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        output, source, do_file = fixture(Path(temporary))
        contract = receipt.build_contract(output, source, do_file)
        assert contract["target"]["estimate"] == 0.45
        assert math.isclose(
            contract["target"]["inverse_variance_weight"], 1.0 / 0.12**2, abs_tol=1e-12
        )
        rows = receipt.read_rows(output / "target_receipt.csv")
        rows[0]["estimation_household_years"] = "99"
        write_csv(output / "target_receipt.csv", rows)
        try:
            receipt.build_contract(output, source, do_file)
        except RuntimeError as error:
            assert "one row per household-year" in str(error)
        else:
            raise AssertionError("Corrupt household-year receipt was accepted")

        output, source, do_file = fixture(Path(temporary) / "second")
        rows = receipt.read_rows(output / "event_study_estimates.csv")
        rows = [row for row in rows if row["coefficient_name"] != "F6event"]
        write_csv(output / "event_study_estimates.csv", rows)
        try:
            receipt.build_contract(output, source, do_file)
        except RuntimeError as error:
            assert "event-time grid" in str(error)
        else:
            raise AssertionError("Receipt without F6event was accepted")


if __name__ == "__main__":
    test_valid_fixture_and_corruption_rejection()
    print("FIRST_BIRTH_ROOMS_RECEIPT_TESTS_PASS tests=3")
