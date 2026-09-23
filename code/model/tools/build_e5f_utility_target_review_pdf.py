#!/usr/bin/env python3
"""Build a review PDF from the collector's E5F utility overnight readout.

The builder formats collected files. It does not score cases or infer scientific
verification from a favorable objective value.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (BaseDocTemplate, Frame, HRFlowable, Image,
                                PageBreak, PageTemplate, Paragraph, Spacer,
                                Table, TableStyle)

CELLS = ("B_floor", "B_shares", "D_floor", "D_shares")
SHARED = ("H0", "beta_annual", "chi", "first_birth_fixed_cost",
          "kappa_fert", "kappa_fert_continuation", "theta0", "theta1")
EXPECTED_TARGETS = 13
EXPECTED_GRAPHS = 17
PARAM_COUNTS = {"B_floor": 17, "D_floor": 17, "B_shares": 19, "D_shares": 19}

styles = getSampleStyleSheet()
STY = {
    "Title": ParagraphStyle("UTitle", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=18, leading=22, textColor=colors.HexColor("#17324d"), spaceAfter=10),
    "H1": ParagraphStyle("UH1", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=13, leading=16, textColor=colors.HexColor("#17324d"), spaceBefore=5, spaceAfter=7),
    "H2": ParagraphStyle("UH2", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=10, leading=12, textColor=colors.HexColor("#2c526e"), spaceBefore=5, spaceAfter=4),
    "Body": ParagraphStyle("UBody", parent=styles["BodyText"], fontName="Helvetica", fontSize=9.5, leading=13, textColor=colors.HexColor("#222222"), spaceAfter=6),
    "Small": ParagraphStyle("USmall", parent=styles["BodyText"], fontName="Helvetica", fontSize=8, leading=9.5, textColor=colors.HexColor("#333333"), spaceAfter=4),
    "Head": ParagraphStyle("UHead", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=7.5, leading=8.5, textColor=colors.white),
    "Cell": ParagraphStyle("UCell", parent=styles["BodyText"], fontName="Helvetica", fontSize=7.5, leading=8.8),
    "Cap": ParagraphStyle("UCap", parent=styles["BodyText"], fontName="Helvetica-Oblique", fontSize=8, leading=9.5, textColor=colors.HexColor("#555555"), spaceBefore=2, spaceAfter=5),
}


def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def ascii_text(value):
    replacements = {"\u2013": "-", "\u2014": "-", "\u2212": "-", "\u2011": "-",
                    "\u2264": "<=", "\u2265": ">=", "\u03b2": "beta", "\u03c1": "rho",
                    "\u03c3": "sigma", "\u03b8": "theta", "\u03ba": "kappa",
                    "\u03c8": "psi", "\u2018": "'", "\u2019": "'",
                    "\u201c": '"', "\u201d": '"'}
    result = str(value if value is not None else "")
    for old, new in replacements.items():
        result = result.replace(old, new)
    return result.encode("ascii", "ignore").decode("ascii")


def fmt(value, digits=5):
    if value is None or value == "":
        return "-"
    try:
        number = float(value)
        return f"{number:,.2f}" if abs(number) >= 1000 else f"{number:.{digits}g}"
    except (ValueError, TypeError):
        return ascii_text(value)


def para(value, style="Body", trusted_markup=False):
    text = ascii_text(value)
    return Paragraph(text if trusted_markup else escape(text), STY[style])


def make_table(rows, widths, repeat=1):
    data = []
    for row_index, row in enumerate(rows):
        data.append([value if isinstance(value, Paragraph) else para(value, "Head" if row_index == 0 else "Cell") for value in row])
    table = Table(data, colWidths=widths, repeatRows=repeat, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2f607c")),
        ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#b8c5ce")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f1f5f7")]),
    ]))
    return table


class Doc(BaseDocTemplate):
    def __init__(self, filename, title):
        super().__init__(filename, pagesize=letter, leftMargin=.55*inch, rightMargin=.55*inch,
                         topMargin=.58*inch, bottomMargin=.5*inch, title=ascii_text(title))
        self.addPageTemplates([PageTemplate(id="main", frames=[Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="normal")], onPage=self.footer)])

    def footer(self, canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#c8d2d8"))
        canvas.line(.55*inch, .4*inch, 7.95*inch, .4*inch)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#52616b"))
        canvas.drawString(.55*inch, .24*inch, "E5F utility comparison - target review")
        canvas.drawRightString(7.95*inch, .24*inch, str(doc.page))
        canvas.restoreState()


def narrative(path):
    data = read_json(path)
    for key in ("title", "date", "summary", "recommendations", "limitations"):
        if key not in data:
            raise ValueError(f"narrative missing required key: {key}")
    for key in ("summary", "recommendations", "limitations"):
        if not isinstance(data[key], list):
            raise ValueError(f"narrative {key} must be a list of paragraphs")
    for section in data.get("sections", []):
        if not isinstance(section, dict) or not section.get("title"):
            raise ValueError("each narrative section needs a title")
        if not isinstance(section.get("paragraphs", []), list):
            raise ValueError("section paragraphs must be a list")
        table = section.get("table")
        if table and (not isinstance(table.get("headers"), list) or not isinstance(table.get("rows"), list)):
            raise ValueError("section table needs headers and rows")
    for source in data.get("sources", []):
        if not source.get("label") or not str(source.get("url", "")).startswith(("https://", "http://")):
            raise ValueError("sources need a label and http(s) URL")
    return data


def moment_name(row):
    return "Initial fertility normalization" if row.get("restriction_id") == "initial_normalization" else row.get("label", row.get("restriction_id", ""))


def fingerprint_and_plan(base, summary_hash):
    score = read_json(base / "score.json")
    plan = read_json(base / "plan.json")
    receipt_path = base / "verified_evaluation_receipt.json"
    if not receipt_path.is_file():
        raise ValueError(f"{base}: collector-verified evaluation receipt missing")
    receipt = read_json(receipt_path)
    target_hash = plan.get("target_system_sha256")
    if not target_hash or target_hash != summary_hash:
        raise ValueError(f"{base}: target/weight fingerprint differs from collection summary")
    if receipt.get("checkpoint_sha256") != score.get("checkpoint_sha256"):
        raise ValueError(f"{base}: score and verification receipt checkpoint hashes differ")
    if len(score.get("target_fit", [])) != EXPECTED_TARGETS:
        raise ValueError(f"{base}: expected 13 scored target rows")
    if len(score.get("parameters", [])) != PARAM_COUNTS[base.parent.name]:
        raise ValueError(f"{base}: unexpected estimated-parameter row count")
    if not math.isfinite(float(score.get("loss", float("nan")))):
        raise ValueError(f"{base}: score has no finite loss")
    return score, plan, receipt


def validate_plan_dimensions(plan, cell):
    structural = plan.get("structural_parameters", {})
    expected = set(SHARED) | ({"h_P"} if cell.endswith("floor") else {"delta_alpha_jump", "delta_alpha"})
    if set(structural) != expected:
        raise ValueError(f"{cell}: structural coordinates differ from the declared floor/share design")
    if "psi" in structural:
        raise ValueError(f"{cell}: fertility scale psi must remain separately normalized")


def validate_target_csv(base, score):
    rows = read_csv(base / "target_fit.csv")
    by_id = {row.get("restriction_id"): row for row in rows}
    scored = {row.get("restriction_id"): row for row in score["target_fit"]}
    if len(rows) != EXPECTED_TARGETS or len(by_id) != EXPECTED_TARGETS or set(by_id) != set(scored):
        raise ValueError(f"{base}: target_fit.csv identities/count differ from score.json")
    for key, item in scored.items():
        csv_row = by_id[key]
        for field in ("target", "model", "gap", "actual_weight", "loss_contribution"):
            expected = item.get(field)
            actual = csv_row.get(field, "")
            if expected is None:
                if actual not in (None, ""):
                    raise ValueError(f"{base}: CSV {field} should be empty for {key}")
            elif float(actual) != float(expected):
                raise ValueError(f"{base}: CSV {field} differs from score.json for {key}")


def validate_targets(score, cell, reference):
    rows = score["target_fit"]
    by_id = {row.get("restriction_id"): row for row in rows}
    if len(by_id) != EXPECTED_TARGETS:
        raise ValueError(f"{cell}: duplicated/missing target identities")
    empirical = [row for row in rows if row.get("restriction_id") != "initial_normalization"]
    norm = by_id.get("initial_normalization")
    if len(empirical) != 12 or not norm or norm.get("scored") is not False:
        raise ValueError(f"{cell}: expected 12 scored targets and one separate normalization")
    if norm.get("actual_weight") is not None or norm.get("loss_contribution") is not None:
        raise ValueError(f"{cell}: separate normalization must not enter the weighted loss")
    contributions = []
    for row in empirical:
        if not row.get("scored", True):
            raise ValueError(f"{cell}: an empirical target is not scored")
        if not math.isclose(float(row["model"]) - float(row["target"]), float(row["gap"]), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"{cell}: gap does not reconcile for {row['restriction_id']}")
        weight, gap = float(row["actual_weight"]), float(row["gap"])
        contribution = float(row["loss_contribution"])
        if not all(math.isfinite(value) for value in (weight, gap, contribution)) or weight <= 0:
            raise ValueError(f"{cell}: nonfinite or nonpositive target weight for {row['restriction_id']}")
        if not math.isclose(weight * gap * gap, contribution, rel_tol=1e-12, abs_tol=1e-9):
            raise ValueError(f"{cell}: weighted squared-gap contribution does not reconcile for {row['restriction_id']}")
        contributions.append(contribution)
    if not math.isclose(math.fsum(contributions), float(score["loss"]), rel_tol=1e-12, abs_tol=1e-9):
        raise ValueError(f"{cell}: 12 target contributions do not sum to the scored loss")
    ids = set(by_id)
    if reference is not None and ids != set(reference):
        raise ValueError(f"{cell}: target identities differ across utility cells")
    if reference is not None:
        for key in ids:
            left, right = by_id[key], reference[key]
            if float(left["target"]) != float(right["target"]):
                raise ValueError(f"{cell}: target differs across utility cells for {key}")
            if key != "initial_normalization" and float(left["actual_weight"]) != float(right["actual_weight"]):
                raise ValueError(f"{cell}: actual weight differs across utility cells for {key}")
    return by_id


def selected_comparison(scores):
    ref = next(iter(scores.values()))
    ids = list(ref)
    rows = [["Target moment", "Target", *CELLS]]
    for key in ids:
        rows.append([moment_name(ref[key]), fmt(ref[key]["target"]), *[
            fmt(scores[cell][key]["model"]) if cell in scores else "Unavailable" for cell in CELLS]])
    return rows


def common_smoke_parameters(plans):
    if set(plans) != set(CELLS):
        return None
    vectors = {cell: plan.get("structural_parameters", {}) for cell, plan in plans.items()}
    shared = [name for name in SHARED if all(name in vectors[cell] for cell in CELLS)]
    if tuple(shared) != SHARED:
        raise ValueError("smoke plans do not contain the eight declared shared coordinates")
    for name in shared:
        values = [float(vectors[cell][name]) for cell in CELLS]
        if not all(value == values[0] for value in values[1:]):
            raise ValueError(f"smoke anchor differs in shared structural coordinate {name}")
    if any("psi" in vectors[cell] for cell in CELLS):
        raise ValueError("psi must be separately normalized, not a structural search coordinate")
    return shared


def counts_story(readout):
    inventory = read_json(readout / "inventory.json")
    counts = {}
    for cell in CELLS:
        counts[cell] = {state: sum(row.get("cell") == cell and row.get("status") == state for row in inventory)
                        for state in ("unrun", "verified_scored", "rejected", "failed", "running", "incomplete", "collection_rejected")}
    return counts


def build_story(narr, readout):
    story = [para(narr["title"], "Title"), para(narr["date"], "Small"),
             HRFlowable(width="100%", color=colors.HexColor("#2f607c"), thickness=1), Spacer(1, 8),
             para("Summary", "H1")]
    for item in narr["summary"]:
        story.append(para(item))
    if readout is None:
        story += [para("Collection status", "H1"), para("Overnight calibration results remain uncollected. The empirical review and frozen target definitions are available below; new model moments, estimated parameters, numerical verification and diagnostic plots are unavailable.")]
    else:
        summary_path = readout / "collection_summary.json"
        inventory_path = readout / "inventory.json"
        if not summary_path.is_file() or not inventory_path.is_file():
            raise ValueError("readout must contain collection_summary.json and inventory.json")
        summary = read_json(summary_path)
        if not summary.get("target_system_sha256"):
            raise ValueError("collection summary does not pin the target/weight fingerprint")
        status_counts = counts_story(readout)
        selected = {}
        smoke = {}
        selected_plans = {}
        smoke_plans = {}
        common_parameter_block = None
        selection = {item["cell"]: item for item in read_json(readout / "selection.json")}
        for cell in CELLS:
            selected_base = readout / cell / "selected"
            if selected_base.is_dir():
                score, plan, _ = fingerprint_and_plan(selected_base, summary["target_system_sha256"])
                validate_plan_dimensions(plan, cell)
                validate_target_csv(selected_base, score)
                chosen = selection.get(cell)
                if not chosen:
                    raise ValueError(f"{cell}: selected files are present without a collector selection record")
                if not math.isclose(float(chosen.get("loss", float("nan"))), float(score["loss"]), rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f"{cell}: selection record loss differs from copied score")
                if chosen.get("checkpoint_sha256") != score.get("checkpoint_sha256"):
                    raise ValueError(f"{cell}: selection record checkpoint differs from copied score")
                selected[cell] = validate_targets(score, cell, next(iter(selected.values()), None))
                selected_plans[cell] = plan
            smoke_base = readout / cell / "smoke_anchor"
            if smoke_base.is_dir():
                score, plan, _ = fingerprint_and_plan(smoke_base, summary["target_system_sha256"])
                validate_plan_dimensions(plan, cell)
                validate_target_csv(smoke_base, score)
                smoke[cell] = validate_targets(score, cell, next(iter(smoke.values()), None))
                smoke_plans[cell] = plan
        if selected:
            story += [para("Selected points: all 13 target moments", "H1"),
                      para("B uses one persistent earnings process; D adds iid risk. All cells retain heterogeneous entry wealth. Each value is from the lowest-loss verified smoke or production case. Repeats are excluded from selection. These finite-search results are not a convergence certificate.")]
            story.append(make_table(selected_comparison(selected), [2.45*inch, .72*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch]))
        if smoke:
            story += [PageBreak(), para("Common smoke anchors: all 13 target moments", "H1"),
                      para("These are the cell-specific smoke_01 anchor cases copied by the collector. Missing cells remain unavailable. A smoke anchor is a checked starting point, not a search result.")]
            story.append(make_table(selected_comparison(smoke), [2.45*inch, .72*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch]))
        if len(smoke_plans) == len(CELLS):
            names = common_smoke_parameters(smoke_plans)
            common_rows = [["Shared coordinate", "Smoke anchor value"]]
            for name in names:
                common_rows.append([name, fmt(smoke_plans[CELLS[0]]["structural_parameters"][name])])
            common_parameter_block = [para("Shared smoke parameters verified from plans", "H2"),
                                      para("The eight listed structural coordinates are identical in all four smoke plans. Utility specifications still differ in their remaining structural dimensions: floor cells search h_P, while child-dependent-share cells search two share tilts. Fertility utility scale psi is normalized separately and is not included among these fixed coordinates."),
                                      make_table(common_rows, [2.2*inch, 1.2*inch])]
        else:
            common_parameter_block = [para("Shared smoke parameters", "H2"), para("Unavailable: all four verified smoke_anchor plans were not present, so shared-parameter equality was not asserted.")]
        total_counts = {state: sum(item[state] for item in status_counts.values()) for state in next(iter(status_counts.values()), {})}
        story += [PageBreak(), para("Collection and verification status", "H1"),
                  para(f"Collector inventory: {len(read_json(readout / 'inventory.json'))} planned cases; {summary.get('attempted_cases', 'unavailable')} attempted; {summary.get('selected_cells', 'unavailable')} selected cells. This status describes collection, not economic adequacy.")]
        count_rows = [["Cell", "Verified scores", "Rejected", "Failed", "Running/started unresolved", "Incomplete", "Unrun", "Collection rejected"]]
        for cell in CELLS:
            c = status_counts[cell]
            count_rows.append([cell, c["verified_scored"], c["rejected"], c["failed"], c["running"], c["incomplete"], c["unrun"], c["collection_rejected"]])
        count_rows.append(["Total"] + [total_counts.get(k, 0) for k in ("verified_scored", "rejected", "failed", "running", "incomplete", "unrun", "collection_rejected")])
        story.append(make_table(count_rows, [.82*inch, .78*inch, .58*inch, .48*inch, 1.05*inch, .82*inch, .45*inch, .9*inch]))
        repeat_rows = [["Cell", "Runner-reported completed", "Collector-verified repeats", "Finalizer receipt status"]]
        for cell in CELLS:
            receipt_record = next((r for r in read_json(readout / "verification_receipts.json") if r.get("cell") == cell), {})
            repeats = receipt_record.get("repeats", [])
            completed = sum(item.get("status") == "completed" for item in repeats)
            checked = int(receipt_record.get("verified_repeat_count", 0))
            final = receipt_record.get("scientific_receipt_status", "not supplied by collector")
            if final == "not_present":
                final = "not present"
            elif final.startswith("raw_receipt_reports_exact_twice"):
                final = "raw receipt reports exact twice; selected binding matches; not independently revalidated"
            elif final.startswith("raw_receipt_status_requires_review"):
                final = "raw receipt status needs review; selected binding matches"
            elif final.startswith("not_applicable_"):
                final = "requires review: " + final.removeprefix("not_applicable_").replace("_", " ")
            repeat_rows.append([cell, completed, checked, final])
        verification_text = narr.get("scientific_verification_status", "Unavailable; repeat completion alone is not a scientific comparison.")
        story += [Spacer(1, 8), para("Selected-point repeat receipts", "H2"),
                  para("Runner-reported completions, collector-scored repeats, and raw finalizer receipts are separate evidence. A raw finalizer statement is not independently revalidated here. File-based running or started states are unresolved without separate scheduler or process confirmation."),
                  make_table(repeat_rows, [1.05*inch, 1.2*inch, 1.25*inch, 2.5*inch])]
        story += [Spacer(1, 5), para("Scientific verification status", "H2"), para(verification_text)]
        story += [Spacer(1, 8), *common_parameter_block]
    story += [PageBreak(), para("Recommendations", "H1")]
    for item in narr["recommendations"]:
        story.append(para(item))
    story.append(para("Limitations", "H1"))
    for item in narr["limitations"]:
        story.append(para(item))
    for section in narr.get("sections", []):
        story += [PageBreak(), para(section["title"], "H1")]
        for item in section.get("paragraphs", []):
            story.append(para(item))
        if section.get("table"):
            table = section["table"]
            widths = [7.35*inch/max(1, len(table["headers"]))] * len(table["headers"])
            story.append(make_table([table["headers"], *table["rows"]], widths))
    if readout is not None:
        story.append(PageBreak())
        for cell in CELLS:
            base = readout / cell / "selected"
            if not base.is_dir():
                story += [para(f"Appendix: {cell} selected readout unavailable", "H1"), para("No verified selected case was collected for this cell. No target-fit, parameter, or figure values are supplied."), PageBreak()]
                continue
            score, plan, receipt = fingerprint_and_plan(base, summary["target_system_sha256"])
            validate_plan_dimensions(plan, cell)
            targets = score["target_fit"]
            params = read_csv(base / "parameters_actual_bounds.csv")
            if len(params) != PARAM_COUNTS[cell]:
                raise ValueError(f"{cell}: actual-bounds table row count differs from {PARAM_COUNTS[cell]}")
            target_rows = [["Target moment", "Target", "Model", "Gap", "Actual weight", "Loss contribution"]]
            for row in targets:
                weight = row.get("actual_weight") if row.get("actual_weight") is not None else None
                contribution = row.get("loss_contribution") if row.get("loss_contribution") is not None else None
                target_rows.append([moment_name(row), fmt(row.get("target")), fmt(row.get("model")), fmt(row.get("gap")), fmt(weight), fmt(contribution)])
            story += [para(f"Appendix: {cell} full selected target fit", "H1"),
                      para(f"Collector-selected stage: {selection[cell]['stage']}; case: {selection[cell]['case_id']}; loss: {fmt(score.get('loss'), 9)}. Collection verification receipt checkpoint SHA-256: {receipt.get('checkpoint_sha256', 'unavailable')}.") ,
                      make_table(target_rows, [2.15*inch, .75*inch, .75*inch, .75*inch, 1.0*inch, 1.0*inch]),
                      Spacer(1, 7), para(f"All {PARAM_COUNTS[cell]} parameter estimates, bounds, and status", "H2"),
                      para("Bounds are taken from the collector's parameters_actual_bounds.csv. 'externally_fixed_or_derived' indicates the parameter is not a searched coordinate in this plan; a missing bound is shown as unavailable.", "Small")]
            parameter_rows = [["Parameter", "Estimate", "Actual lower", "Actual upper", "Near bound", "Status"]]
            for row in params:
                status = row.get("restriction_type", "")
                if status == "searched":
                    status = "searched coordinate"
                parameter_rows.append([row.get("parameter", ""), fmt(row.get("estimate")), fmt(row.get("actual_lower")), fmt(row.get("actual_upper")), fmt(row.get("near_actual_bound")), status])
            story.append(make_table(parameter_rows, [1.55*inch, .8*inch, .85*inch, .85*inch, .8*inch, 1.55*inch]))
            story.append(PageBreak())
            plots = sorted((base / "standard_diagnostics").glob("*.png"))
            if len(plots) != EXPECTED_GRAPHS:
                raise ValueError(f"{cell}: expected 17 original standard plots, found {len(plots)}")
            for index in range(0, len(plots), 2):
                story.append(para(f"Appendix: {cell} original standard diagnostics ({index+1}-{min(index+2, len(plots))} of 17)", "H1"))
                for image_path in plots[index:index+2]:
                    story.append(para(image_path.stem.replace("_", " "), "H2"))
                    image = Image(str(image_path))
                    image._restrictSize(7.1*inch, 3.55*inch)
                    story += [image, para("Original diagnostic image retained unchanged.", "Cap")]
                if index + 2 < len(plots):
                    story.append(PageBreak())
            story.append(PageBreak())
    if narr.get("sources"):
        if not story or not isinstance(story[-1], PageBreak):
            story.append(PageBreak())
        story.append(para("Sources", "H1"))
        for source in narr["sources"]:
            label = escape(ascii_text(source["label"]))
            url = escape(str(source["url"]), {'"': "&quot;"})
            story.append(Paragraph(f'<link href="{url}" color="#2f607c">{label}</link> - {escape(ascii_text(source["url"]))}', STY["Body"]))
    return story


def build(readout_path, narrative_path, output_path):
    narr = narrative(narrative_path)
    readout = Path(readout_path) if readout_path else None
    if readout is not None and not readout.is_dir():
        raise FileNotFoundError(f"readout directory unavailable: {readout}")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    story = build_story(narr, readout)
    Doc(str(output), narr["title"]).build(story)
    from pypdf import PdfReader
    return {"output": str(output), "pages": len(PdfReader(str(output)).pages),
            "readout_available": readout is not None,
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readout", type=Path, help="collector output directory; omit when unavailable")
    parser.add_argument("--narrative", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.readout, args.narrative, args.output), indent=2))


if __name__ == "__main__":
    main()
