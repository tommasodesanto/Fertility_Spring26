from __future__ import annotations

import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
WORKDIR = SCRIPT_DIR.parents[1]
SPRING26 = WORKDIR.parent
FERTILITY_ROOT = SPRING26.parent

PSID_OUTPUT = SPRING26 / "code" / "data" / "psid_followup_mar2026" / "output"
MMS_OUTPUT = SPRING26 / "code" / "data" / "mms_center_periphery" / "output_middle_center"
OUTDIR = WORKDIR / "output" / "empirical_roundup_v1"


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def row_by_value(rows: list[dict[str, str]], key: str, value: str) -> dict[str, str]:
    for row in rows:
        if row.get(key) == value:
            return row
    raise KeyError(f"Missing row where {key} == {value!r}")


def row_by_two_values(
    rows: list[dict[str, str]], key1: str, value1: str, key2: str, value2: str
) -> dict[str, str]:
    for row in rows:
        if row.get(key1) == value1 and row.get(key2) == value2:
            return row
    raise KeyError(f"Missing row where {key1} == {value1!r} and {key2} == {value2!r}")


def to_float(value: str | None) -> float:
    if value is None or value == "":
        return float("nan")
    return float(value)


def fmt(x: float, digits: int = 3) -> str:
    return f"{x:.{digits}f}"


def parse_mediation_log(path: Path) -> dict[str, float | None]:
    text = path.read_text(encoding="utf-8")

    wealth_matches = re.findall(
        r"(?m)^total_n~2530 \|\s+([-.0-9]+)\s+[.0-9]+\s+[-.0-9]+\s+([.0-9]+)",
        text,
    )
    owner_matches = re.findall(
        r"(?m)^ owner_by_35 \|\s+([-.0-9]+)\s+[.0-9]+\s+[-.0-9]+\s+([.0-9]+)",
        text,
    )

    out: dict[str, float | None] = {
        "wealth_total_coef": None,
        "wealth_total_p": None,
        "wealth_first_stage_coef": None,
        "wealth_first_stage_p": None,
        "wealth_controlled_coef": None,
        "wealth_controlled_p": None,
        "owner_controlled_coef": None,
        "owner_controlled_p": None,
    }

    if len(wealth_matches) >= 1:
        out["wealth_total_coef"] = float(wealth_matches[0][0])
        out["wealth_total_p"] = float(wealth_matches[0][1])
    if len(wealth_matches) >= 2:
        out["wealth_first_stage_coef"] = float(wealth_matches[1][0])
        out["wealth_first_stage_p"] = float(wealth_matches[1][1])
    if len(wealth_matches) >= 3:
        out["wealth_controlled_coef"] = float(wealth_matches[2][0])
        out["wealth_controlled_p"] = float(wealth_matches[2][1])
    if owner_matches:
        out["owner_controlled_coef"] = float(owner_matches[-1][0])
        out["owner_controlled_p"] = float(owner_matches[-1][1])

    return out


def weighted_gap_within_wealth(path: Path) -> float:
    rows = read_csv_rows(path)
    stayed = {
        int(row["wealth_q"]): (to_float(row["childless_by_45"]), to_float(row["n_obs"]))
        for row in rows
        if row["owner_by_35"] == "0"
    }
    became = {
        int(row["wealth_q"]): (to_float(row["childless_by_45"]), to_float(row["n_obs"]))
        for row in rows
        if row["owner_by_35"] == "1"
    }

    numer = 0.0
    denom = 0.0
    for q in sorted(set(stayed) & set(became)):
        stayed_rate, stayed_n = stayed[q]
        became_rate, became_n = became[q]
        weight = stayed_n + became_n
        numer += (stayed_rate - became_rate) * weight
        denom += weight
    return numer / denom


def build_roundup() -> tuple[list[dict[str, str]], str]:
    location_rows = read_csv_rows(MMS_OUTPUT / "mms_location_by_parent_summary.csv")
    move_reg_rows = read_csv_rows(MMS_OUTPUT / "mms_move_regressions_allparent.csv")
    origin_rows = read_csv_rows(MMS_OUTPUT / "mms_origin_transition_summary_allparent.csv")
    four_way_rows = read_csv_rows(MMS_OUTPUT / "mms_four_way_move_shares_allparent.csv")

    first_birth_rooms = read_csv_rows(
        PSID_OUTPUT / "sa_rooms_first_birth_one_variant_v1" / "rooms_f_c_y_all_summary.csv"
    )[0]
    first_birth_one_two = read_csv_rows(
        PSID_OUTPUT / "rooms_first_birth_one_vs_two_horizon_v1" / "rooms_first_birth_one_vs_two_horizon_v1.csv"
    )
    first_birth_group_one = read_csv_rows(
        PSID_OUTPUT / "sa_rooms_first_birth_grouped_v1" / "rooms_f_c_y_one_kid_by3_summary.csv"
    )[0]
    first_birth_group_two = read_csv_rows(
        PSID_OUTPUT / "sa_rooms_first_birth_grouped_v1" / "rooms_f_c_y_two_plus_by3_summary.csv"
    )[0]
    second_birth_rooms = read_csv_rows(
        PSID_OUTPUT
        / "sa_rooms_second_birth_with_onechild_controls_v1"
        / "rooms_s_c_y_all_summary.csv"
    )[0]
    within_wealth_gap = weighted_gap_within_wealth(
        PSID_OUTPUT / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"
    )
    dual_rates = read_csv_rows(PSID_OUTPUT / "fertility_wealth_v1" / "dual_rates_quintile_v9.csv")
    move_for_size_design = read_csv_rows(
        PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_design_summary_v1.csv"
    )[0]
    move_for_size_f = read_csv_rows(
        PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_fstats_v1.csv"
    )
    moved_to_own_state = read_csv_rows(
        PSID_OUTPUT / "moved_to_own_state_change_v1" / "moved_to_own_state_change_summary_v1.csv"
    )
    wealth_terciles = read_csv_rows(PSID_OUTPUT / "wealth_v1" / "wealth_tercile_outcomes_v1.csv")
    med = parse_mediation_log(PSID_OUTPUT / "fertility_wealth_v1" / "fertility_wealth_v9.log")

    first_birth_hazard_tercile = read_csv_rows(
        PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_tercile.csv"
    )
    first_birth_hazard_quintile = read_csv_rows(
        PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_quintile.csv"
    )
    hazard_broken = all(to_float(row["first_birth_rate"]) == 0.0 for row in first_birth_hazard_tercile) and all(
        to_float(row["first_birth_rate"]) == 0.0 for row in first_birth_hazard_quintile
    )

    nonparents = row_by_value(location_rows, "parent_status", "Non-Parents")
    newparents = row_by_value(location_rows, "parent_status", "New Parents")

    center_center_np = row_by_value(
        origin_rows, "origin_label", "Center"
    )
    # The file has two center-origin rows and two periphery-origin rows. Split manually.
    center_origin_rows = [row for row in origin_rows if row["origin_label"] == "Center"]
    periphery_origin_rows = [row for row in origin_rows if row["origin_label"] == "Periphery"]

    def origin_share(rows: list[dict[str, str]], parent_group: str, dest: str) -> float:
        row = next(
            row
            for row in rows
            if row["parent_compare_all"] == parent_group and row["dest_label"] == dest
        )
        return to_float(row["share"])

    move_across = row_by_value(move_reg_rows, "outcome", "Across-CBSA mover")
    move_center = row_by_value(move_reg_rows, "outcome", "Center destination among within-CBSA movers")

    origin_reg_rows = read_csv_rows(MMS_OUTPUT / "mms_origin_regressions_allparent.csv")
    origin_center = row_by_value(
        origin_reg_rows,
        "outcome",
        "Center destination among within-CBSA movers, center-origin mass",
    )
    origin_periphery = row_by_value(
        origin_reg_rows,
        "outcome",
        "Center destination among within-CBSA movers, periphery-origin mass",
    )

    one_kid_by2 = row_by_value(first_birth_one_two, "horizon_group", "one_kid_by2")
    two_plus_by2 = row_by_value(first_birth_one_two, "horizon_group", "two_plus_by2")
    f_twins = row_by_value(move_for_size_f, "stat", "F_twins_first_stage")
    f_samesex = row_by_value(move_for_size_f, "stat", "F_samesex_first_stage")
    interstate_post3 = row_by_two_values(
        moved_to_own_state,
        "sample",
        "post3",
        "measure",
        "share_interstate_geo_weighted",
    )

    rows = [
        {
            "prediction_id": "P1a",
            "prediction": "Parent share falls with local rent, conditional on income",
            "status": "partially_tested",
            "exportability": "broader_sorting_fact_slide_ready",
            "headline": (
                f"Center share: non-parents {fmt(to_float(nonparents['center_share']))} vs "
                f"new parents {fmt(to_float(newparents['center_share']))}"
            ),
            "key_files": "; ".join(
                [
                    str(MMS_OUTPUT / "mms_location_by_parent_summary.csv"),
                    str(MMS_OUTPUT / "mms_location_by_parent.png"),
                ]
            ),
            "notes": "Broader sorting fact is established; exact conditional rent-slope regression not found on disk.",
        },
        {
            "prediction_id": "P1b",
            "prediction": "First birth triggers relocation toward lower-price areas",
            "status": "partially_tested",
            "exportability": "indirect_support_exportable",
            "headline": (
                f"Parents are {fmt(to_float(move_center['estimate_pp']), 1)} pp less likely to choose a "
                f"center destination among within-CBSA movers; legacy PSID move-for-size mean = "
                f"{fmt(to_float(move_for_size_design['moved_for_size_post3']))}"
            ),
            "key_files": "; ".join(
                [
                    str(MMS_OUTPUT / "mms_move_regressions_allparent.csv"),
                    str(MMS_OUTPUT / "mms_origin_transition_summary_allparent.csv"),
                    str(FERTILITY_ROOT / "Outputs" / "Graphs" / "mv_s_f_c_y_all.png"),
                ]
            ),
            "notes": "Strong indirect support via ACS within-metro destination sorting and PSID moves for size; no direct cheaper-area birth event study.",
        },
        {
            "prediction_id": "P2a",
            "prediction": "Parents have lower housing income elasticity",
            "status": "untested_but_feasible",
            "exportability": "missing",
            "headline": "Current PSID workflow contains income, rooms, and ownership; no elasticity script existed before this turn.",
            "key_files": "",
            "notes": "Implemented in this turn as a new Stata script using household FE and housing proxies.",
        },
        {
            "prediction_id": "P2b",
            "prediction": "Parents have lower housing price elasticity",
            "status": "untested_not_feasible_current_workflow",
            "exportability": "missing",
            "headline": "No live PSID metro panel or cleaned rent merge exists in the current workflow.",
            "key_files": "",
            "notes": "Do not treat this as an easy extension under the current data layout.",
        },
        {
            "prediction_id": "P3a",
            "prediction": "Family flight is decreasing in pre-birth wealth",
            "status": "untested_spatially",
            "exportability": "indirect_wealth_margin_exportable",
            "headline": (
                f"Weighted childlessness gap between stay-renters and become-owners within wealth bins = "
                f"{fmt(within_wealth_gap, 3)}"
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"),
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "dual_rates_quintile_v9.csv"),
                ]
            ),
            "notes": "Wealth clearly matters for ownership/fertility margins, but not yet for spatial relocation at birth.",
        },
        {
            "prediction_id": "P3b",
            "prediction": "Wealth threshold above which parents stay central",
            "status": "untested",
            "exportability": "missing",
            "headline": "No direct spatial wealth-threshold exercise is on disk.",
            "key_files": "",
            "notes": "Requires a spatial measure that the current PSID workflow does not have.",
        },
        {
            "prediction_id": "L3a",
            "prediction": "First child causes a discrete housing-demand jump",
            "status": "tested",
            "exportability": "slide_ready",
            "headline": (
                f"First-birth rooms: +{fmt(to_float(first_birth_rooms['coef_p3']))} at k=+3; "
                f"one kid by +2 = +{fmt(to_float(one_kid_by2['mean_change']))}, "
                f"two-plus by +2 = +{fmt(to_float(two_plus_by2['mean_change']))}"
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "sa_rooms_first_birth_one_variant_v1" / "rooms_f_c_y_all_summary.csv"),
                    str(PSID_OUTPUT / "rooms_first_birth_one_vs_two_horizon_v1" / "rooms_first_birth_one_vs_two_horizon_v1.csv"),
                ]
            ),
            "notes": "One of the cleanest PSID facts in the repo.",
        },
        {
            "prediction_id": "L3a_ext",
            "prediction": "Second birth raises housing demand again",
            "status": "tested",
            "exportability": "slide_ready",
            "headline": f"Second-birth rooms: +{fmt(to_float(second_birth_rooms['coef_p3']))} at k=+3",
            "key_files": str(
                PSID_OUTPUT
                / "sa_rooms_second_birth_with_onechild_controls_v1"
                / "rooms_s_c_y_all_summary.csv"
            ),
            "notes": "Useful extension of the discrete family-space activation fact.",
        },
        {
            "prediction_id": "L3b",
            "prediction": "Low slack makes the first-child housing cost sharply nonlinear",
            "status": "untested",
            "exportability": "missing",
            "headline": "No finished slack-hazard or low-slack spline design is on disk.",
            "key_files": "",
            "notes": "Implemented in this turn only as a recommended next Stata path, not as a finished executed result.",
        },
        {
            "prediction_id": "L3c",
            "prediction": "Relaxing the ownership/family-space constraint raises fertility",
            "status": "strongly_supported_descriptive",
            "exportability": "slide_ready",
            "headline": (
                f"Mediation log: wealth p(total) = {fmt(med['wealth_total_p'] or float('nan'), 3)}, "
                f"wealth p(controlled) = {fmt(med['wealth_controlled_p'] or float('nan'), 3)}, "
                f"ownership coef(controlled) = {fmt(med['owner_controlled_coef'] or float('nan'), 3)}"
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"),
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "fertility_wealth_v9.log"),
                ]
            ),
            "notes": "Strong reduced-form support for the ownership constraint margin, though not a final causal claim.",
        },
        {
            "prediction_id": "QC1",
            "prediction": "First-birth hazard by wealth artifact is usable",
            "status": "broken",
            "exportability": "not_exportable",
            "headline": f"Hazard CSVs all-zero = {str(hazard_broken).lower()}",
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_tercile.csv"),
                    str(PSID_OUTPUT / "fertility_wealth_v1" / "first_birth_hazard_by_quintile.csv"),
                ]
            ),
            "notes": "Exists on disk but should not be cited until debugged.",
        },
        {
            "prediction_id": "QC2",
            "prediction": "Post-birth wealth-transition summaries are core evidence",
            "status": "exploratory_only",
            "exportability": "exploratory",
            "headline": (
                f"Top pre-birth liquid-wealth tercile own_post3 = {fmt(to_float(wealth_terciles[2]['own_post3']))}; "
                f"middle tercile own_post3 = {fmt(to_float(wealth_terciles[1]['own_post3']))}"
            ),
            "key_files": "; ".join(
                [
                    str(PSID_OUTPUT / "wealth_v1" / "wealth_tercile_outcomes_v1.csv"),
                    str(PSID_OUTPUT / "wealth_v1" / "wealth_transition_regressions_v1.tex"),
                ]
            ),
            "notes": "Interesting, but not yet clean enough for the headline narrative.",
        },
        {
            "prediction_id": "SUPP1",
            "prediction": "Birth raises moves for space",
            "status": "tested",
            "exportability": "supported_but_split_across_folders",
            "headline": (
                f"Moved-for-size post3 mean = {fmt(to_float(move_for_size_design['moved_for_size_post3']))}; "
                f"twins first-stage F = {fmt(to_float(f_twins['value']), 1)}"
            ),
            "key_files": "; ".join(
                [
                    str(FERTILITY_ROOT / "Outputs" / "Graphs" / "mv_s_f_c_y_all.png"),
                    str(PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_design_summary_v1.csv"),
                    str(PSID_OUTPUT / "moved_for_size_iv_alltenure_v1" / "moved_for_size_fstats_v1.csv"),
                ]
            ),
            "notes": "Core event-study fact is legacy; IV follow-up is robustness/mechanism support.",
        },
        {
            "prediction_id": "SUPP2",
            "prediction": "Move-to-own transitions are mostly not interstate around first birth",
            "status": "tested",
            "exportability": "exportable",
            "headline": (
                f"Among moved-to-own events in post3, interstate_geo share = "
                f"{fmt(to_float(interstate_post3['value']))}"
            ),
            "key_files": str(
                PSID_OUTPUT / "moved_to_own_state_change_v1" / "moved_to_own_state_change_summary_v1.csv"
            ),
            "notes": "Useful support for the idea that the housing adjustment margin is relatively local, but not a substitute for a spatial price design.",
        },
    ]

    md_lines = [
        "# Empirical Roundup v1",
        "",
        "This file consolidates the current ACS and PSID evidence already on disk.",
        "",
        "## Main Read",
        "",
        "- ACS already carries the spatial sorting claims.",
        "- PSID already carries the family-space activation, ownership, and wealth-margin claims.",
        "- The still-missing direct tests are price elasticity attenuation, slack nonlinearity, and wealth-moderated spatial flight.",
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
            "## Highest-Value Next Steps",
            "",
            "- Run first-birth housing event studies by pre-birth wealth quintile for rooms, ownership, and move-for-size.",
            "- Run parent-vs-childless income-elasticity regressions using PSID rooms and ownership as housing proxies.",
            "- Do not treat MSA-rent merges as part of the current PSID workflow unless a real metro identifier is first verified.",
        ]
    )

    return rows, "\n".join(md_lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, markdown = build_roundup()

    csv_path = OUTDIR / "empirical_roundup_v1.csv"
    md_path = OUTDIR / "empirical_roundup_v1.md"

    fieldnames = [
        "prediction_id",
        "prediction",
        "status",
        "exportability",
        "headline",
        "key_files",
        "notes",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    md_path.write_text(markdown, encoding="utf-8")

    print(csv_path)
    print(md_path)


if __name__ == "__main__":
    main()
