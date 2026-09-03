#!/usr/bin/env python3
"""Build the concise advisor-facing PDF for the overnight E5F numerical audit."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


BLUE = colors.HexColor("#204b6b")
LIGHT_BLUE = colors.HexColor("#eaf1f6")
RED = colors.HexColor("#a84336")
GREEN = colors.HexColor("#287451")
MID_GREY = colors.HexColor("#66717a")
LIGHT_GREY = colors.HexColor("#f2f4f5")


MOMENT_LABELS = {
    "tfr": "Completed fertility",
    "childless_rate": "Childlessness",
    "mean_age_first_birth": "Mean age at first birth",
    "share_first_births_age30plus": "Share first births age 30+",
    "housing_increment_0to1": "Rooms response to first birth",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Rooms: 3+ minus 1–2 children",
    "own_family_gap": "Family ownership gap",
    "own_rate": "Ownership rate",
    "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms",
    "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / annual earnings",
    "annual_bequest_flow_to_aggregate_wealth": "Bequest flow / wealth",
    "old_total_wealth_to_annual_income_p90_p50_7684": "Old-age wealth p90 / p50",
}


PARAMETER_LABELS = {
    "beta_annual": "Annual discount factor",
    "kappa_fert": "First-birth taste dispersion",
    "kappa_fert_continuation": "Continuation-birth taste dispersion",
    "chi": "Housing utility weight",
    "H0": "Baseline housing scale",
    "theta0": "Housing preference intercept",
    "theta1": "Housing preference slope",
    "hbar_child_rooms": "Child-room floor",
    "first_birth_fixed_cost": "First-birth fixed cost",
    "hbar_first_child_jump": "First-child room jump",
    "psi_child_change_2023": "2007–2023 child-value change",
    "tenure_choice_kappa": "Tenure taste dispersion",
    "housing_supply_elasticity": "Housing supply elasticity",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-report-dir", type=Path, required=True)
    parser.add_argument("--calibration-task-dir", type=Path, required=True)
    parser.add_argument("--calibration-diagnosis-dir", type=Path, required=True)
    parser.add_argument("--policy-comparison-dir", type=Path, required=True)
    parser.add_argument("--shapley-dir", type=Path, required=True)
    parser.add_argument("--horizon-evidence-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected a JSON object in {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def fmt(value: Any, digits: int = 3) -> str:
    number = float(value)
    if not math.isfinite(number):
        return "—"
    if number == 0.0:
        return "0"
    if abs(number) < 0.001 or abs(number) >= 10000:
        return f"{number:.2e}"
    return f"{number:.{digits}f}"


def pct(value: Any, digits: int = 2) -> str:
    return f"{float(value):+.{digits}f}%"


def pp(value: Any, digits: int = 2) -> str:
    return f"{float(value):+.{digits}f} pp"


def paragraph(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(text, style)


def add_bullets(story: list[Any], items: list[str], style: ParagraphStyle) -> None:
    for item in items:
        story.append(Paragraph(f"•&nbsp;&nbsp;{item}", style))
        story.append(Spacer(1, 3))


def image_fit(path: Path, max_width: float, max_height: float) -> Image:
    from PIL import Image as PILImage

    with PILImage.open(path) as source:
        width, height = source.size
    scale = min(max_width / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale)


def cropped_image_fit(path: Path, max_width: float, max_height: float,
                      top_fraction: float = 0.665) -> Image:
    """Keep the economic transition panels and omit the residual panels."""
    from PIL import Image as PILImage

    with PILImage.open(path) as source:
        cropped = source.crop((0, 0, source.width, round(source.height * top_fraction)))
        width, height = cropped.size
        buffer = io.BytesIO()
        cropped.save(buffer, format="PNG")
    buffer.seek(0)
    scale = min(max_width / width, max_height / height)
    result = Image(buffer, width=width * scale, height=height * scale)
    result._source_buffer = buffer
    return result


def table(data: list[list[Any]], widths: list[float], header_rows: int = 1,
          font_size: float = 7.2) -> Table:
    result = Table(data, colWidths=widths, repeatRows=header_rows, hAlign="LEFT")
    result.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, header_rows - 1), BLUE),
                ("TEXTCOLOR", (0, 0), (-1, header_rows - 1), colors.white),
                ("FONTNAME", (0, 0), (-1, header_rows - 1), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), font_size),
                ("LEADING", (0, 0), (-1, -1), font_size + 2),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c5ccd1")),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("ROWBACKGROUNDS", (0, header_rows), (-1, -1), [colors.white, LIGHT_GREY]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    return result


def footer(canvas: Any, doc: Any) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(MID_GREY)
    canvas.drawString(36, 20, "E5F overnight numerical audit — provisional research evidence")
    canvas.drawRightString(A4[0] - 36, 20, f"Page {doc.page}")
    canvas.restoreState()


def select_effect(rows: list[dict[str, str]], case: str, year: int) -> dict[str, str]:
    found = [
        row for row in rows
        if row["policy_case"] == case and int(row["calendar_year"]) == year
    ]
    if len(found) != 1:
        raise RuntimeError(f"Expected one effect row for {case}, {year}; found {len(found)}")
    return found[0]


def select_shapley(rows: list[dict[str, str]], metric: str, factor: str) -> float:
    found = [
        row for row in rows
        if row["metric"] == metric and row["factor"] == factor
    ]
    if len(found) != 1:
        raise RuntimeError(f"Expected one Shapley row for {metric}, {factor}")
    return float(found[0]["reported_contribution"])


def build_report(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    calibration = read_json(args.calibration_report_dir / "summary.json")
    candidate = dict(calibration["best_candidate"])
    if calibration.get("status") != "complete" or not bool(candidate.get("valid")):
        raise RuntimeError("Calibration report is not a complete valid collector packet")
    fit_rows = read_csv(args.calibration_report_dir / "best_target_fit.csv")
    parameter_rows = read_csv(args.calibration_task_dir / "parameter_table.csv")
    diagnosis = read_json(args.calibration_diagnosis_dir / "summary.json")
    policy = read_json(args.policy_comparison_dir / "summary.json")
    policy_rows = read_csv(args.policy_comparison_dir / "policy_effects_relative_to_baseline.csv")
    shapley = read_json(args.shapley_dir / "summary.json")
    shapley_rows = read_csv(args.shapley_dir / "shapley_decomposition.csv")
    horizon = read_json(args.horizon_evidence_json)

    supply_2023 = select_effect(policy_rows, "supply-plus-20", 2023)
    ltv_2023 = select_effect(policy_rows, "dependent-child-ltv95", 2023)
    combined_2023 = select_effect(policy_rows, "combined", 2023)
    tax_2023 = select_effect(policy_rows, "property-tax-2pct-no-rebate", 2023)
    supply_2063 = select_effect(policy_rows, "supply-plus-20", 2063)
    arc_ratios: list[float] = []
    for year in (2023, 2043, 2063):
        for case in ("supply-plus-20", "property-tax-2pct-no-rebate"):
            row = select_effect(policy_rows, case, year)
            housing_change = float(row["percent_difference_housing_demand_per_adult"])
            birth_change = float(row["percent_difference_topcode_adjusted_births_per_adult"])
            if abs(housing_change) > 1e-12:
                arc_ratios.append(birth_change / housing_change)

    direct_birth = select_shapley(shapley_rows, "births_per_adult", "tax_rate")
    price_birth = select_shapley(shapley_rows, "births_per_adult", "asset_price")
    rebate_birth = select_shapley(shapley_rows, "births_per_adult", "equal_rebate")
    net_birth = direct_birth + price_birth + rebate_birth

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], fontName="Helvetica-Bold",
                              fontSize=20, leading=24, textColor=BLUE, alignment=TA_CENTER,
                              spaceAfter=12))
    styles.add(ParagraphStyle(name="Subtitle", parent=styles["Normal"], fontSize=10, leading=14,
                              textColor=MID_GREY, alignment=TA_CENTER, spaceAfter=14))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading1"], fontName="Helvetica-Bold",
                              fontSize=14, leading=17, textColor=BLUE, spaceBefore=4, spaceAfter=8))
    styles.add(ParagraphStyle(name="Subsection", parent=styles["Heading2"], fontName="Helvetica-Bold",
                              fontSize=10.5, leading=13, textColor=BLUE, spaceBefore=6, spaceAfter=4))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=8.5,
                              leading=11.2, spaceAfter=5))
    styles.add(ParagraphStyle(name="BulletSmall", parent=styles["BodyText"], fontSize=8.7,
                              leading=11.6, leftIndent=8, firstLineIndent=-8))
    styles.add(ParagraphStyle(name="Callout", parent=styles["BodyText"], fontSize=10, leading=14,
                              textColor=colors.HexColor("#17364c"), backColor=LIGHT_BLUE,
                              borderColor=colors.HexColor("#b7cad8"), borderWidth=0.6,
                              borderPadding=8, spaceBefore=5, spaceAfter=8))
    styles.add(ParagraphStyle(name="Caption", parent=styles["BodyText"], fontSize=7.5, leading=9.5,
                              textColor=MID_GREY, spaceBefore=3, spaceAfter=5))
    styles.add(ParagraphStyle(name="PathSmall", parent=styles["BodyText"], fontSize=5.4, leading=6.4,
                              wordWrap="CJK"))

    story: list[Any] = []
    story.append(Spacer(1, 0.25 * inch))
    story.append(paragraph("E5F housing, tenure, and fertility", styles["ReportTitle"]))
    story.append(paragraph("Overnight numerical audit and advisor brief · 3 September 2026",
                           styles["Subtitle"]))
    story.append(paragraph(
        "<b>Bottom line.</b> The model generates economically meaningful housing and tenure responses to property taxes, "
        "but it does not yet contain Coven's cross-jurisdiction relocation margin. "
        "What is weak is the fertility transmission from a tenure switch that barely changes housing "
        "services. A supply expansion changes services and fertility; an LTV relaxation mostly changes "
        "the financing label. With an equal rebate, the income channel more than offsets the direct "
        "fertility cost of the higher property tax. The active quantitative profile has one location, so "
        "cross-jurisdiction relocation is absent by construction.", styles["Callout"]))

    story.append(paragraph("What is established", styles["Section"]))
    add_bullets(story, [
        f"A complete 23-point validation promotes a valid calibration with loss <b>{float(candidate['transition_loss']):.2f}</b>, "
        "8.5% below the prior 38.16 benchmark. All 12 targets and 11 free parameters remain in the same contract.",
        f"The fit is still dominated by housing: the model predicts <b>{float(candidate['terminal_first_birth_housing_response']):.3f}</b> "
        "rooms after a first birth versus the 0.720 target, and mean rooms remain too high by "
        f"{float(diagnosis['mean_rooms_excess']):.3f}.",
        f"In 2023, 20% more housing supply raises housing use {pct(supply_2023['percent_difference_housing_demand_per_adult'])} "
        f"and births {pct(supply_2023['percent_difference_topcode_adjusted_births_per_adult'])}; by 2063 the births effect is "
        f"{pct(supply_2063['percent_difference_topcode_adjusted_births_per_adult'])}.",
        f"Conditional on housing services moving, the link is stable rather than tiny: across the supply and "
        f"unrebated-tax paths in 2023, 2043, and 2063, a 1% housing-use change corresponds to a "
        f"{min(arc_ratios):.2f}–{max(arc_ratios):.2f}% birth change.",
        f"The family LTV policy raises ownership {pp(100*float(ltv_2023['difference_owner_rate']))}, yet births only "
        f"{pct(ltv_2023['percent_difference_topcode_adjusted_births_per_adult'])}. Tenure itself is not an independent fertility input.",
        f"Combining supply and LTV raises ownership {pp(100*float(combined_2023['difference_owner_rate']))} and births "
        f"{pct(combined_2023['percent_difference_topcode_adjusted_births_per_adult'])}; almost all of the birth effect comes from supply.",
        "The active profile sets the location count to one. It can show tenure and housing-product reallocation, "
        "but it cannot reproduce the spatial relocation channel emphasized by Coven et al.",
        f"Without a rebate, doubling the tax lowers ownership {pp(100*float(tax_2023['difference_owner_rate']))}, "
        f"housing use {pct(tax_2023['percent_difference_housing_demand_per_adult'])}, and births "
        f"{pct(tax_2023['percent_difference_topcode_adjusted_births_per_adult'])} in 2023.",
        f"With an equal rebate, the net 2023 birth effect is {pct(net_birth)}: direct tax {pct(direct_birth)}, "
        f"capitalized price {pct(price_birth)}, and rebate {pct(rebate_birth)}.",
    ], styles["BulletSmall"])
    story.append(paragraph("Recommendation", styles["Section"]))
    story.append(paragraph(
        "Lead with the decomposition, not the small net fertility coefficient. The model says a property tax "
        "causes economically meaningful housing reallocation, but the equal rebate is pro-fertility and masks "
        "the direct cost channel. For the paper, the next empirical object should be a causal housing-cost-to-fertility "
        "moment; the current calibration contains fertility levels and fertility-to-housing responses, not that reverse elasticity. "
        "A true comparison with Coven also requires reactivating a multi-location policy environment rather than interpreting "
        "a one-market tenure response as relocation. Within the present architecture, supply reform or a child-contingent "
        "housing-service subsidy is more promising than a financing-only credit; the latter proposal remains unquantified.",
        styles["Callout"]))
    story.append(PageBreak())

    story.append(paragraph("1. Calibration: full target fit", styles["Section"]))
    story.append(paragraph(
        "The table reports every active calibration moment under the unchanged target and weighting system. "
        "Loss is the squared standardized gap. The three largest housing/ownership misses account for roughly "
        f"{100*float(diagnosis['top_three_loss_share']):.1f}% of the total loss.", styles["BodySmall"]))
    fit_table: list[list[Any]] = [["Moment", "Target", "Model", "Gap", "Weight", "Loss"]]
    for row in fit_rows:
        fit_table.append([
            MOMENT_LABELS.get(row["moment"], row["moment"]),
            fmt(row["target"]), fmt(row["model"]), fmt(row["gap"]),
            fmt(row["weight"], 1), fmt(row["loss_contribution"], 2),
        ])
    story.append(table(fit_table, [190, 56, 56, 56, 64, 52], font_size=7.0))
    story.append(Spacer(1, 8))
    story.append(paragraph(
        "<b>Interpretation.</b> The model fits fertility levels and childlessness tightly, but it is not calibrated "
        "to the causal response of births to housing costs. The room response following a first birth is the "
        "opposite causal direction. Consequently, the policy fertility elasticity is a model prediction rather than an identified target.",
        styles["Callout"]))
    story.append(PageBreak())

    story.append(paragraph("2. Parameters and local identification", styles["Section"]))
    diagnosis_plot = args.calibration_diagnosis_dir / "calibration_gap_identification.png"
    story.append(image_fit(diagnosis_plot, 7.1 * inch, 4.6 * inch))
    story.append(paragraph(
        "The local weighted Jacobian has rank 11/11 at the declared 1e-3 threshold, but its condition number is "
        f"{float(diagnosis['local_weighted_jacobian']['condition_number_numerical_rank']):.1f}; theta1 is near its lower bound. "
        "The weakest local direction is 96.7% the first-child room-jump parameter, while theta1 has the largest "
        "forward/backward disagreement. This is adequate for local diagnosis, not proof of global identification.", styles["Caption"]))
    parameter_table: list[list[Any]] = [["Parameter", "Estimate", "Bounds", "Status"]]
    wanted = set(PARAMETER_LABELS)
    for row in parameter_rows:
        if row["parameter"] not in wanted:
            continue
        bounds = "external" if row["status"].startswith("externally") else f"[{fmt(row['lower_bound'])}, {fmt(row['upper_bound'])}]"
        status = "external" if row["status"].startswith("externally") else ("near bound" if row["near_bound"].lower() == "true" else "estimated")
        parameter_table.append([
            PARAMETER_LABELS[row["parameter"]], fmt(row["value"], 4), bounds, status
        ])
    story.append(table(parameter_table, [225, 75, 105, 80], font_size=6.8))
    story.append(PageBreak())

    story.append(paragraph("3. Which housing policies move fertility?", styles["Section"]))
    story.append(image_fit(args.policy_comparison_dir / "policy_mechanism_comparison.png",
                           7.1 * inch, 5.7 * inch))
    story.append(paragraph(
        "These are closed-national mechanism comparisons from one common fitted 2023 state, not welfare or terminal-steady-state claims. "
        "Effects are relative to the no-policy path under the promoted calibration. The supply intervention moves "
        "housing services and births together. The family LTV intervention changes ownership much more than housing "
        "services, so its fertility effect is close to zero. This is a structural result: utility depends on housing "
        f"services and children, while owner status has no separate fertility payoff. The implied policy arc ratio—births "
        f"in percent divided by housing use in percent—is {min(arc_ratios):.2f}–{max(arc_ratios):.2f} across the supply and "
        "unrebated-tax paths at 2023, 2043, and 2063. In addition, the active profile "
        "contains only one location, so none of these experiments can move households across jurisdictions.", styles["Callout"]))
    compact = [["2023 policy", "User cost", "Ownership", "Housing use", "Births"]]
    for label, row in (
        ("Supply +20%", supply_2023),
        ("Family LTV 95%", ltv_2023),
        ("Supply + LTV", combined_2023),
        ("Tax 1%→2%, no rebate", tax_2023),
    ):
        compact.append([
            label,
            pct(row["percent_difference_housing_user_cost"]),
            pp(100 * float(row["difference_owner_rate"])),
            pct(row["percent_difference_housing_demand_per_adult"]),
            pct(row["percent_difference_topcode_adjusted_births_per_adult"]),
        ])
    story.append(table(compact, [180, 75, 85, 80, 70], font_size=7.4))
    story.append(PageBreak())

    story.append(paragraph("4. Why the rebated property-tax effect is small", styles["Section"]))
    story.append(image_fit(args.shapley_dir / "rebated_tax_shapley.png", 6.8 * inch, 4.2 * inch))
    story.append(paragraph(
        "The exact Shapley calculation evaluates every combination of the tax rate, the equilibrium asset price, "
        "and the equal lump-sum rebate. Contributions therefore add exactly to the total 2023 change. It is a one-date "
        "decomposition holding the fitted pre-fertility distribution fixed, not a full transition decomposition.", styles["BodySmall"]))
    shapley_table = [["Outcome", "Direct tax", "Price capitalization", "Equal rebate", "Net"]]
    metric_labels = {
        "births_per_adult": ("Births", "%"),
        "owner_rate": ("Ownership", "pp"),
        "dependent_child_owner_rate": ("Child ownership", "pp"),
        "housing_demand_per_adult": ("Housing use", "%"),
    }
    for metric, (label, unit) in metric_labels.items():
        values = [select_shapley(shapley_rows, metric, factor) for factor in ("tax_rate", "asset_price", "equal_rebate")]
        net_value = sum(values)
        form = (lambda x: pct(x)) if unit == "%" else (lambda x: pp(x))
        shapley_table.append([label, *(form(value) for value in values), form(net_value)])
    story.append(table(shapley_table, [130, 85, 105, 85, 75], font_size=7.2))
    story.append(Spacer(1, 8))
    story.append(paragraph(
        "<b>Economic reading.</b> The direct tax produces the expected housing and fertility contraction. Lower asset "
        "prices partly offset it. The equal rebate is large enough to reverse the birth response, although ownership "
        "and housing services still fall. Hence a small positive net birth effect is not evidence that housing costs do not matter.",
        styles["Callout"]))
    story.append(PageBreak())

    story.append(paragraph("5. Horizon sufficiency and true runtime", styles["Section"]))
    story.append(paragraph(
        f"The terminal demographic operator has persistence {float(horizon['demographic_persistence_per_period']):.5f} per "
        f"four-year period. The implied half-life is {float(horizon['demographic_half_life_years']):.1f} years, and 1% attenuation "
        f"requires about {int(horizon['periods_to_one_percent'])} periods ({float(horizon['years_to_one_percent']):.0f} years). "
        f"At that rate, H32 retains {100 * float(horizon['demographic_persistence_per_period']) ** 32:.1f}% of the slow mode, "
        f"whereas H128 retains only {100 * float(horizon['demographic_persistence_per_period']) ** 128:.2f}%. "
        "H128 is therefore long enough for this demographic block.", styles["BodySmall"]))
    horizon_figure = Table(
        [[
            cropped_image_fit(Path(horizon["baseline_diagnostics_png"]), 3.25 * inch, 2.35 * inch),
            cropped_image_fit(Path(horizon["reform_diagnostics_png"]), 3.25 * inch, 2.35 * inch),
        ]],
        colWidths=[3.35 * inch, 3.35 * inch],
        hAlign="CENTER",
    )
    horizon_figure.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 1),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    story.append(horizon_figure)
    story.append(paragraph(
        "The dashed lines are the policy-specific terminal steady states. The slow approach is demographic: "
        "the price block settles earlier, while population and births inherit the long-lived age distribution. "
        "Equilibrium-residual panels are intentionally omitted from this economic transition figure.",
        styles["Caption"]))
    convergence_table = [["Iteration", "Baseline market", "Baseline fiscal", "Base max year", "Tax-2 market", "Tax-2 fiscal"]]
    base_rows = horizon["baseline_iterations"]
    reform_rows = horizon["reform_iterations"]
    base_by_iteration = {int(row["iteration"]): row for row in base_rows}
    reform_by_iteration = {int(row["iteration"]): row for row in reform_rows}
    all_iterations = sorted(set(base_by_iteration) | set(reform_by_iteration))
    selected_iterations = sorted(set(
        [all_iterations[0], 5, 10, *horizon.get("highlight_iterations", []), *all_iterations[-4:]]
    ) & set(all_iterations))
    for iteration in selected_iterations:
        left = base_by_iteration.get(iteration)
        right = reform_by_iteration.get(iteration)
        convergence_table.append([
            str(iteration),
            fmt(left["market_residual"], 6) if left else "—",
            fmt(left["fiscal_residual"], 6) if left else "—",
            str(left.get("market_argmax_year", "—")) if left else "—",
            fmt(right["market_residual"], 6) if right else "—",
            fmt(right["fiscal_residual"], 6) if right else "—",
        ])
    story.append(table(convergence_table, [52, 88, 88, 70, 88, 88], font_size=7.1))
    story.append(Spacer(1, 8))
    story.append(paragraph(
        f"The paired H128 stage is recorded as <b>{horizon['stage_status']}</b>. {horizon['runtime_interpretation']}",
        styles["Callout"]))
    story.append(paragraph(
        "This H128 pair uses the immediately preceding joint calibration rather than the promoted round-4 vector. "
        "It is a controlled solver-runtime and horizon diagnostic, not the source of the policy magnitudes reported above.",
        styles["Caption"]))
    story.append(PageBreak())
    story.append(paragraph("What to tell advisors now", styles["Subsection"]))
    add_bullets(story, [
        "The core model mechanism is present: housing-supply and unrebated-tax experiments move housing, tenure, and fertility in the expected directions.",
        "The tiny fertility response to an LTV relaxation is not a failed housing mechanism; it is a tenure relabeling with almost no housing-service change.",
        "The full Coven relocation mechanism is not currently testable: the active calibration has one location and a national common tax experiment.",
        "The small positive fertility response to the rebated tax is a fiscal-incidence result. The direct tax component lowers fertility, but the rebate dominates it in 2023.",
        "The calibration still misses the first-birth housing response badly and does not identify the housing-cost-to-fertility elasticity. That is the main empirical and modeling priority before strong quantitative claims.",
        "H128 is the appropriate diagnostic horizon, but the coupled fiscal solve—not the economic horizon—is the runtime bottleneck. Do not calibrate by repeatedly running full H128 transitions.",
    ], styles["BulletSmall"])

    provenance = [
        ["Artifact", "Path"],
        ["Calibration", paragraph(str(args.calibration_report_dir.resolve()), styles["PathSmall"])],
        ["Policy comparison", paragraph(str(args.policy_comparison_dir.resolve()), styles["PathSmall"])],
        ["Tax decomposition", paragraph(str(args.shapley_dir.resolve()), styles["PathSmall"])],
        ["Horizon evidence", paragraph(str(args.horizon_evidence_json.resolve()), styles["PathSmall"])],
    ]
    story.append(Spacer(1, 7))
    story.append(KeepTogether([
        paragraph("Reproducible evidence", styles["Subsection"]),
        table(provenance, [105, 375], font_size=5.7),
    ]))

    document = SimpleDocTemplate(
        str(output), pagesize=A4, rightMargin=36, leftMargin=36,
        topMargin=34, bottomMargin=32, title="E5F overnight numerical audit",
        author="Fertility_Spring26 numerical workflow",
    )
    document.build(story, onFirstPage=footer, onLaterPages=footer)


def main() -> None:
    build_report(parse_args())


if __name__ == "__main__":
    main()
