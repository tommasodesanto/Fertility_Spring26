from __future__ import annotations

import csv
from pathlib import Path

import compile_empirical_roundup_v1 as v1


SCRIPT_DIR = Path(__file__).resolve().parent
WORKDIR = SCRIPT_DIR.parents[1]
OUTDIR = WORKDIR / "output" / "empirical_roundup_v2"

SPRING26 = WORKDIR.parent
PSID_OUTPUT = SPRING26 / "code" / "data" / "psid_followup_mar2026" / "output"
MMS_OUTPUT = SPRING26 / "code" / "data" / "mms_center_periphery" / "output_middle_center"

NEW_OUT_INCOME = WORKDIR / "output" / "empirical_roundup_income_elasticity_v1"
NEW_OUT_WEALTH = WORKDIR / "output" / "empirical_roundup_first_birth_by_wealth_v1"


def _row(rows: list[dict[str, str]], **kwargs: str) -> dict[str, str]:
    for row in rows:
        if all(row.get(k) == v for k, v in kwargs.items()):
            return row
    raise KeyError(f"Missing row with filters {kwargs}")


def _f(x: float, digits: int = 3) -> str:
    return f"{x:.{digits}f}"


def _monotone(values: list[float]) -> bool:
    return all(a <= b for a, b in zip(values, values[1:])) or all(
        a >= b for a, b in zip(values, values[1:])
    )


def build_roundup() -> tuple[list[dict[str, str]], str]:
    location_rows = v1.read_csv_rows(MMS_OUTPUT / "mms_location_by_parent_summary.csv")
    move_reg_rows = v1.read_csv_rows(MMS_OUTPUT / "mms_move_regressions_allparent.csv")
    first_birth_rooms = v1.read_csv_rows(
        PSID_OUTPUT / "sa_rooms_first_birth_one_variant_v1" / "rooms_f_c_y_all_summary.csv"
    )[0]
    first_birth_one_two = v1.read_csv_rows(
        PSID_OUTPUT / "rooms_first_birth_one_vs_two_horizon_v1" / "rooms_first_birth_one_vs_two_horizon_v1.csv"
    )
    second_birth_rooms = v1.read_csv_rows(
        PSID_OUTPUT
        / "sa_rooms_second_birth_with_onechild_controls_v1"
        / "rooms_s_c_y_all_summary.csv"
    )[0]
    med = v1.parse_mediation_log(PSID_OUTPUT / "fertility_wealth_v1" / "fertility_wealth_v9.log")
    within_wealth_gap = v1.weighted_gap_within_wealth(
        PSID_OUTPUT / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"
    )
    move_for_size_design = v1.read_csv_rows(
        PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_design_summary_v1.csv"
    )[0]
    moved_to_own_state = v1.read_csv_rows(
        PSID_OUTPUT / "moved_to_own_state_change_v1" / "moved_to_own_state_change_summary_v1.csv"
    )
    income_rows = v1.read_csv_rows(NEW_OUT_INCOME / "income_elasticity_summary_v1.csv")
    wealth_rows = v1.read_csv_rows(NEW_OUT_WEALTH / "first_birth_by_wealth_summary_v1.csv")
    wealth_bins = v1.read_csv_rows(NEW_OUT_WEALTH / "wealth_quintile_summary_v1.csv")

    nonparents = v1.row_by_value(location_rows, "parent_status", "Non-Parents")
    newparents = v1.row_by_value(location_rows, "parent_status", "New Parents")
    move_center = v1.row_by_value(
        move_reg_rows, "outcome", "Center destination among within-CBSA movers"
    )
    interstate_post3 = v1.row_by_two_values(
        moved_to_own_state, "sample", "post3", "measure", "share_interstate_geo_weighted"
    )
    one_kid_by2 = v1.row_by_value(first_birth_one_two, "horizon_group", "one_kid_by2")
    two_plus_by2 = v1.row_by_value(first_birth_one_two, "horizon_group", "two_plus_by2")
    ln_rooms = _row(income_rows, outcome="ln_rooms")
    own = _row(income_rows, outcome="own")

    rooms_p3 = [
        v1.to_float(_row(wealth_rows, outcome="rooms", wealth_q=str(q))["coef_p3"])
        for q in range(1, 6)
    ]
    own_p3 = [
        v1.to_float(_row(wealth_rows, outcome="own", wealth_q=str(q))["coef_p3"])
        for q in range(1, 6)
    ]
    move_p3 = [
        v1.to_float(_row(wealth_rows, outcome="moved_for_size", wealth_q=str(q))["coef_p3"])
        for q in range(1, 6)
    ]
    base_rooms = [
        v1.to_float(_row(wealth_rows, outcome="rooms", wealth_q=str(q))["pre_event_mean"])
        for q in range(1, 6)
    ]
    base_own = [
        v1.to_float(_row(wealth_rows, outcome="own", wealth_q=str(q))["pre_event_mean"])
        for q in range(1, 6)
    ]

    q1_med = v1.to_float(_row(wealth_bins, wealth_q="1")["median_nw_pre"])
    q5_med = v1.to_float(_row(wealth_bins, wealth_q="5")["median_nw_pre"])

    hazard_tercile = v1.read_csv_rows(
        PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_tercile.csv"
    )
    hazard_quintile = v1.read_csv_rows(
        PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_quintile.csv"
    )
    hazard_broken = all(v1.to_float(r["first_birth_rate"]) == 0.0 for r in hazard_tercile) and all(
        v1.to_float(r["first_birth_rate"]) == 0.0 for r in hazard_quintile
    )

    rows = [
        {
            "prediction_id": "P1a",
            "prediction": "Parents sort away from expensive central locations",
            "status": "tested",
            "exportability": "slide_ready",
            "headline": (
                f"Center share falls from {_f(v1.to_float(nonparents['center_share']))} "
                f"for non-parents to {_f(v1.to_float(newparents['center_share']))} for new parents."
            ),
            "key_files": "; ".join(
                [
                    str(MMS_OUTPUT / "mms_location_by_parent_summary.csv"),
                    str(MMS_OUTPUT / "mms_location_by_parent.png"),
                ]
            ),
            "notes": "This is the clean ACS sorting fact. The exact conditional-on-income rent-slope regression is still not isolated separately.",
        },
        {
            "prediction_id": "P1b",
            "prediction": "Birth shifts families toward cheaper/peripheral housing markets",
            "status": "partially_tested",
            "exportability": "exportable_indirect",
            "headline": (
                f"Among within-CBSA movers, parents are {_f(v1.to_float(move_center['estimate_pp']), 1)} pp "
                f"less likely to choose a center destination; PSID moved-for-size post3 mean = "
                f"{_f(v1.to_float(move_for_size_design['moved_for_size_post3']))}."
            ),
            "key_files": "; ".join(
                [
                    str(MMS_OUTPUT / "mms_move_regressions_allparent.csv"),
                    str(PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_design_summary_v1.csv"),
                ]
            ),
            "notes": "The spatial piece is established in ACS and the family-space move margin is established in PSID, but there is still no direct PSID cheaper-destination event study.",
        },
        {
            "prediction_id": "P2a",
            "prediction": "Parents have lower housing income elasticity",
            "status": "tested_new",
            "exportability": "slide_ready",
            "headline": (
                f"FE log-rooms elasticity: childless {_f(v1.to_float(ln_rooms['childless_effect']))}, "
                f"parents {_f(v1.to_float(ln_rooms['parent_effect']))}; "
                f"FE ownership semi-elasticity: childless {_f(v1.to_float(own['childless_effect']))}, "
                f"parents {_f(v1.to_float(own['parent_effect']))}."
            ),
            "key_files": "; ".join(
                [
                    str(NEW_OUT_INCOME / "income_elasticity_summary_v1.csv"),
                    str(NEW_OUT_INCOME / "income_elasticity_models_v1.tex"),
                ]
            ),
            "notes": "This now supports the Stone-Geary attenuation prediction directly in PSID. The parent interaction is negative and precisely estimated in both outcomes.",
        },
        {
            "prediction_id": "P2b",
            "prediction": "Parents have lower housing price elasticity",
            "status": "not_run_current_workflow",
            "exportability": "missing",
            "headline": "Still not identified cleanly in the current PSID workflow because there is no usable metro-rent panel in the implemented shelf data.",
            "key_files": "",
            "notes": "This remains a real gap, not a skipped line in the existing code.",
        },
        {
            "prediction_id": "P3a",
            "prediction": "Wealth should mute the first-birth housing adjustment / family flight margin",
            "status": "proxy_run_no_clean_gradient",
            "exportability": "diagnostic_not_headline",
            "headline": (
                f"Pre-birth housing is strongly wealth-graded: baseline rooms rise from {_f(base_rooms[0])} to {_f(base_rooms[-1])}, "
                f"and baseline ownership from {_f(base_own[0])} to {_f(base_own[-1])}; "
                f"but k=+3 event-study effects are not monotone across wealth quintiles."
            ),
            "key_files": "; ".join(
                [
                    str(NEW_OUT_WEALTH / "first_birth_by_wealth_summary_v1.csv"),
                    str(NEW_OUT_WEALTH / "rooms_overlay_v1.png"),
                    str(NEW_OUT_WEALTH / "own_overlay_v1.png"),
                    str(NEW_OUT_WEALTH / "moved_for_size_overlay_v1.png"),
                ]
            ),
            "notes": (
                f"Executed proxy test says wealth clearly shifts the baseline housing stock, but not in a stable monotone post-birth jump. "
                f"Monotone at k=+3: rooms={str(_monotone(rooms_p3)).lower()}, own={str(_monotone(own_p3)).lower()}, move_for_size={str(_monotone(move_p3)).lower()}."
            ),
        },
        {
            "prediction_id": "L3a",
            "prediction": "First birth causes a discrete housing-demand jump",
            "status": "tested",
            "exportability": "slide_ready",
            "headline": (
                f"First-birth rooms are +{_f(v1.to_float(first_birth_rooms['coef_p3']))} at k=+3; "
                f"one-kid by +2 = +{_f(v1.to_float(one_kid_by2['mean_change']))}, "
                f"two-plus by +2 = +{_f(v1.to_float(two_plus_by2['mean_change']))}."
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "sa_rooms_first_birth_one_variant_v1" / "rooms_f_c_y_all_summary.csv"),
                    str(PSID_OUTPUT / "rooms_first_birth_one_vs_two_horizon_v1" / "rooms_first_birth_one_vs_two_horizon_v1.csv"),
                ]
            ),
            "notes": "Still one of the strongest PSID facts in the repo.",
        },
        {
            "prediction_id": "L3a_ext",
            "prediction": "Second birth also raises housing demand",
            "status": "tested",
            "exportability": "slide_ready",
            "headline": f"Second-birth rooms are +{_f(v1.to_float(second_birth_rooms['coef_p3']))} at k=+3.",
            "key_files": str(
                PSID_OUTPUT
                / "sa_rooms_second_birth_with_onechild_controls_v1"
                / "rooms_s_c_y_all_summary.csv"
            ),
            "notes": "Useful extension of the family-space mechanism.",
        },
        {
            "prediction_id": "L3b",
            "prediction": "Low slack should sharply depress first-birth hazard",
            "status": "still_missing",
            "exportability": "missing",
            "headline": "No executed slack-hazard design yet.",
            "key_files": "",
            "notes": "This is now the main distinctive theory prediction that still lacks a direct empirical counterpart.",
        },
        {
            "prediction_id": "L3c",
            "prediction": "Ownership / family-space relaxation raises fertility",
            "status": "strong_support",
            "exportability": "slide_ready",
            "headline": (
                f"Within wealth bins, the weighted childlessness gap between stay-renters and become-owners is "
                f"{_f(within_wealth_gap)}; mediation log gives ownership coef = {_f(med['owner_controlled_coef'] or float('nan'))}."
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"),
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "fertility_wealth_v9.log"),
                ]
            ),
            "notes": "The current repo evidence keeps pointing to ownership as the operative margin rather than liquid wealth by itself.",
        },
        {
            "prediction_id": "QC1",
            "prediction": "Existing first-birth hazard-by-wealth exports are usable",
            "status": "broken",
            "exportability": "do_not_cite",
            "headline": f"All-zero hazard artifact = {str(hazard_broken).lower()}.",
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_tercile.csv"),
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_quintile.csv"),
                ]
            ),
            "notes": "These should stay out of slides until debugged.",
        },
        {
            "prediction_id": "SUPP1",
            "prediction": "Housing adjustment around birth is mostly local rather than interstate",
            "status": "tested",
            "exportability": "supporting",
            "headline": (
                f"Among move-to-own transitions in post3, the interstate share is "
                f"{_f(v1.to_float(interstate_post3['value']))}."
            ),
            "key_files": str(
                PSID_OUTPUT / "moved_to_own_state_change_v1" / "moved_to_own_state_change_summary_v1.csv"
            ),
            "notes": "Helpful as a support fact, but not a substitute for a true price or center-periphery design inside PSID.",
        },
        {
            "prediction_id": "NEW1",
            "prediction": "Executed wealth-stratified first-birth event studies",
            "status": "completed",
            "exportability": "diagnostic",
            "headline": (
                f"Pre-birth median net worth ranges from {_f(q1_med, 0)} in Q1 to {_f(q5_med, 0)} in Q5; "
                f"all three outcomes now have by-quintile event-study figures on disk."
            ),
            "key_files": "; ".join(
                [
                    str(NEW_OUT_WEALTH / "wealth_quintile_summary_v1.csv"),
                    str(NEW_OUT_WEALTH / "rooms_overlay_v1.png"),
                    str(NEW_OUT_WEALTH / "own_overlay_v1.png"),
                    str(NEW_OUT_WEALTH / "moved_for_size_overlay_v1.png"),
                ]
            ),
            "notes": "Useful for the paper appendix or a robustness slide, but the resulting pattern is not yet the clean wealth-threshold result from the theory.",
        },
    ]

    md_lines = [
        "# Empirical Roundup v2",
        "",
        "This roundup reflects the Stata runs executed in this turn, not just the legacy outputs already on disk.",
        "",
        "## Main Assessment",
        "",
        "- ACS still handles the spatial sorting side cleanly.",
        "- PSID still handles the discrete family-space activation and ownership margin cleanly.",
        "- The new FE regression directly supports lower housing income elasticity for parents.",
        "- The new wealth-stratified first-birth event studies show large baseline housing differences by wealth, but not a clean monotone differential birth response.",
        "- The biggest remaining empirical hole is the low-slack hazard / nonlinear squeeze result; price elasticity is also still missing under the current data layout.",
        "",
        "## Audit Table",
        "",
        "| ID | Prediction | Status | Exportability | Headline | Notes |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        md_lines.append(
            "| {prediction_id} | {prediction} | {status} | {exportability} | {headline} | {notes} |".format(
                **{k: str(v).replace("|", "\\|") for k, v in row.items()}
            )
        )

    md_lines.extend(
        [
            "",
            "## Executed This Turn",
            "",
            f"- [empirical_roundup_income_elasticity_v1.do]({(SCRIPT_DIR / 'empirical_roundup_income_elasticity_v1.do').as_posix()})",
            f"- [empirical_roundup_first_birth_by_wealth_v1.do]({(SCRIPT_DIR / 'empirical_roundup_first_birth_by_wealth_v1.do').as_posix()})",
            f"- [income_elasticity_summary_v1.csv]({(NEW_OUT_INCOME / 'income_elasticity_summary_v1.csv').as_posix()})",
            f"- [first_birth_by_wealth_summary_v1.csv]({(NEW_OUT_WEALTH / 'first_birth_by_wealth_summary_v1.csv').as_posix()})",
        ]
    )

    return rows, "\n".join(md_lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, markdown = build_roundup()

    csv_path = OUTDIR / "empirical_roundup_v2.csv"
    md_path = OUTDIR / "empirical_roundup_v2.md"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "prediction_id",
                "prediction",
                "status",
                "exportability",
                "headline",
                "key_files",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    md_path.write_text(markdown, encoding="utf-8")
    print(csv_path)
    print(md_path)


if __name__ == "__main__":
    main()
