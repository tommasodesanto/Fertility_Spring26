#!/usr/bin/env python3
"""Fail-closed Torch-only PDF builder for the three-case estate diagnostic.

This reader never imports the model or unpickles saved states.  It consumes only
the completed report artifacts written by ``run_e5f_estate_receiver_probe.py``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import socket
from pathlib import Path

CASES = ("control", "net_valuation", "net_valuation_transfer")
COMPLETE_STATUS = "exact_three_case_loop_complete"
CASE_STATUS = "completed_experimental_estate_receiver_diagnostic"
PAGE = (792, 612)  # landscape letter


def require_torch() -> None:
    """Keep hashing and PDF work inside a Torch allocation."""
    if not os.environ.get("SLURM_JOB_ID") and os.environ.get("E5F_ESTATE_REPORT_TORCH") != "1":
        raise RuntimeError("Torch-only: run inside an allocation or set E5F_ESTATE_REPORT_TORCH=1 on Torch")
    if "torch" not in socket.gethostname().lower() and not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Torch-only host guard failed")


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"cannot read JSON: {path}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON object required: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
    except OSError as exc:
        raise RuntimeError(f"cannot read CSV: {path}") from exc
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def fmt(value: object) -> str:
    value = number(value)
    if value is None:
        return ""
    if value and abs(value) < 0.0005:
        return f"{value:.2e}"
    return f"{value:.3f}"


def pick(mapping: dict, *names: str) -> object:
    for name in names:
        if name in mapping:
            return mapping[name]
    return ""


def wrap(value: object, styles, style="BodyText"):
    from reportlab.platypus import Paragraph
    from xml.sax.saxutils import escape
    return Paragraph(escape(str(value)), styles[style])


def validate(run_root: Path, figure_root: Path | None = None) -> tuple[dict, dict[str, dict], list[str], dict[str, str]]:
    complete = read_json(run_root / "complete.json")
    if complete.get("status") != COMPLETE_STATUS:
        raise RuntimeError("complete.json does not certify the exact three-case loop")
    declared = complete.get("cases")
    if not isinstance(declared, list) or [item.get("case") for item in declared] != list(CASES):
        raise RuntimeError("complete.json must list control, net_valuation, net_valuation_transfer in order")
    packets: dict[str, dict] = {}
    fingerprints: dict[str, str] = {}
    if figure_root is not None:
        formatting = read_json(figure_root / "figure_formatting.json")
        if formatting.get("status") != "complete_same_figures_legend_formatting" or formatting.get("native_ge_solves") != 0:
            raise RuntimeError("completed report-only legend formatting required")
        fingerprints["figure_formatting.json"] = sha256(figure_root / "figure_formatting.json")
    diagnostic_names: list[str] | None = None
    required_target = {"moment", "target", "model", "gap", "weight", "loss_contribution"}
    required_parameter = {"parameter", "estimate", "lower", "upper", "near_bound", "status"}
    for case in CASES:
        root = run_root / case
        receipt_path, target_path, parameter_path = root / "receipt.json", root / "target_fit.csv", root / "parameters.csv"
        iterations_path = root / "estate_iterations.json"
        receipt = read_json(receipt_path)
        if receipt.get("status") != CASE_STATUS or receipt.get("case") != case:
            raise RuntimeError(f"{case}: receipt is not a completed matching case")
        targets, parameters = read_csv(target_path), read_csv(parameter_path)
        try:
            iterations = json.loads(iterations_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"{case}: missing or invalid estate_iterations.json") from exc
        if not iterations or iterations[-1].get("status") != "completed":
            raise RuntimeError(f"{case}: final estate iteration is incomplete")
        prices = iterations[-1].get("price", [])
        if len(prices) != 1 or number(prices[0]) is None or float(prices[0]) <= 0.:
            raise RuntimeError(f"{case}: valid one-market price is required")
        if len(targets) != 13 or not required_target.issubset(targets[0]):
            raise RuntimeError(f"{case}: require exactly 13 complete target rows")
        if len(parameters) != 26 or not required_parameter.issubset(parameters[0]):
            raise RuntimeError(f"{case}: parameter rows/bounds/restrictions are incomplete")
        if len({r["moment"] for r in targets}) != 13 or any(
                number(r[key]) is None for r in targets for key in ("target", "model", "gap")):
            raise RuntimeError(f"{case}: target names or numeric values are invalid")
        accounts = receipt.get("estate_accounts", {})
        if any(number(accounts.get(key)) is None for key in (
                "generated_net_period", "generated_gross_period", "paid_period",
                "transfer", "recipient_mass", "period_years", "residual")):
            raise RuntimeError(f"{case}: required estate accounting fields are missing")
        if float(accounts["period_years"]) != 4.:
            raise RuntimeError(f"{case}: period length differs from the four-year reference")
        entrant_flow = number(receipt.get("normalized_entrant_flow"))
        if entrant_flow is None or entrant_flow <= 0.:
            raise RuntimeError(f"{case}: positive household entrant flow is required")
        figures = (figure_root / case if figure_root is not None else root) / "standard_diagnostics"
        if figure_root is not None and formatting["cases"][case]["checkpoint_sha256"] != sha256(root / "initial_state.pkl.gz"):
            raise RuntimeError("formatted figures came from a different checkpoint")
        names = sorted(item.name for item in figures.glob("*.png"))
        if len(names) != 17 or len(set(names)) != 17:
            raise RuntimeError(f"{case}: require exactly 17 diagnostic PNGs")
        if diagnostic_names is None:
            diagnostic_names = names
        elif names != diagnostic_names:
            raise RuntimeError("standard diagnostic names differ across cases")
        pins = {key: receipt.get(key) for key in ("source_manifest_sha256", "target_system_sha256", "selected_reference")}
        if any(not value for value in pins.values()):
            raise RuntimeError(f"{case}: source, target, or selected-reference pin is absent")
        packets[case] = dict(root=root, figures=figures, figure_legends_formatted=figure_root is not None,
                            receipt=receipt, targets=targets, parameters=parameters, price=prices[0])
        fingerprints.update({f"{case}/{receipt_path.name}": sha256(receipt_path),
                             f"{case}/{target_path.name}": sha256(target_path),
                             f"{case}/{parameter_path.name}": sha256(parameter_path),
                             f"{case}/{iterations_path.name}": sha256(iterations_path)})
        for name in names:
            fingerprints[f"{case}/standard_diagnostics/{name}"] = sha256(figures / name)
    for pin in ("source_manifest_sha256", "target_system_sha256", "selected_reference"):
        if len({packets[case]["receipt"][pin] for case in CASES}) != 1:
            raise RuntimeError(f"case packets disagree on {pin}")
    target_contract = [(r["moment"], r["target"], r["weight"]) for r in packets[CASES[0]]["targets"]]
    parameter_names = [r["parameter"] for r in packets[CASES[0]]["parameters"]]
    if len(set(parameter_names)) != 26:
        raise RuntimeError("duplicate parameter rows")
    for case in CASES[1:]:
        if [(r["moment"], r["target"], r["weight"]) for r in packets[case]["targets"]] != target_contract:
            raise RuntimeError("target definitions or weights differ across cases")
        if [r["parameter"] for r in packets[case]["parameters"]] != parameter_names:
            raise RuntimeError("parameter names differ across cases")
    fingerprints["complete.json"] = sha256(run_root / "complete.json")
    return complete, packets, diagnostic_names or [], fingerprints


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColorRGB(.28, .33, .38)
    canvas.drawString(36, 20, "Estate-recipient diagnostic | experimental report-only comparison")
    canvas.drawRightString(PAGE[0] - 36, 20, f"Page {doc.page}")
    canvas.restoreState()


def render_and_check(output: Path, packets: dict) -> dict:
    """Render on Torch and check the PDF's numeric table columns against CSVs."""
    import pymupdf
    from PIL import Image, ImageDraw
    qa = output.with_name(output.stem + "_qa")
    qa.mkdir(exist_ok=False)
    verified = {"target_numeric_cells": 0, "parameter_numeric_cells": 0}
    with pymupdf.open(output) as doc:
        if len(doc) != 64 or sum(bool(p.get_images()) for p in doc) != 51:
            raise RuntimeError("PDF must contain 13 text/table pages and 51 figure pages")
        def column(page, left, right):
            values = []
            for word in page.get_text("words", sort=True):
                if left < (word[0] + word[2]) / 2 < right and number(word[4]) is not None:
                    values.append(word[4])
            return values
        for case in CASES:
            fit_pages = [p for p in doc if f"{case}: complete target fit" in p.get_text()]
            if len(fit_pages) != 1:
                raise RuntimeError("target table must occupy one page: " + case)
            for key, left, right in (("target",306,391), ("model",391,476),
                    ("gap",476,561), ("weight",561,641), ("loss_contribution",641,736)):
                expected = [fmt(r[key]) for r in packets[case]["targets"] if r[key] != ""]
                if column(fit_pages[0], left, right) != expected:
                    raise RuntimeError(f"PDF target numeric column differs: {case}/{key}")
                verified["target_numeric_cells"] += len(expected)
            param_pages = [p for p in doc if f"{case}: all parameters, bounds, and restrictions" in p.get_text()]
            if len(param_pages) != 2:
                raise RuntimeError("parameter tables must occupy two pages: " + case)
            for key, left, right in (("estimate",221,296), ("lower",296,361), ("upper",361,426)):
                expected = [fmt(r[key]) for r in packets[case]["parameters"] if r[key] != ""]
                actual = [value for page in param_pages for value in column(page, left, right)]
                if actual != expected:
                    raise RuntimeError(f"PDF parameter numeric column differs: {case}/{key}")
                verified["parameter_numeric_cells"] += len(expected)
        sheets = []
        for start in range(0, len(doc), 8):
            sheet = Image.new("RGB", (1584, 656), "white")
            draw = ImageDraw.Draw(sheet)
            for offset, index in enumerate(range(start, min(start+8, len(doc)))):
                page = doc[index]
                for word in page.get_text("words"):
                    if word[0] < 30 or word[2] > 762 or word[1] < 15 or word[3] > 602:
                        raise RuntimeError(f"PDF text exceeds page margins: page {index+1}")
                pix = page.get_pixmap(matrix=pymupdf.Matrix(1.5,1.5), alpha=False)
                frame = Image.frombytes("RGB", (pix.width,pix.height), pix.samples)
                if not page.get_images():
                    frame.save(qa / f"page_{index+1:02d}.png")
                frame.thumbnail((396,306))
                x,y = (offset%4)*396, (offset//4)*328
                sheet.paste(frame, (x,y+20)); draw.text((x+4,y+3),f"Page {index+1}",fill="black")
            path = qa / f"contact_{start//8+1:02d}.png"
            sheet.save(path); sheets.append(path.name)
    return dict(**verified, text_margins_passed=True, contact_sheets=sheets,
                visual_review="pending human/agent inspection", qa_directory=str(qa))


def story_for(complete: dict, packets: dict[str, dict], names: list[str]):
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import Spacer, PageBreak, Table, TableStyle, Image

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=7.2, leading=8.5))
    styles.add(ParagraphStyle(name="TableHeader", parent=styles["Small"], textColor=colors.white))
    styles["Title"].fontSize, styles["Title"].leading = 20, 24
    styles["Heading1"].fontSize, styles["Heading1"].leading = 14, 17
    styles["BodyText"].fontSize, styles["BodyText"].leading = 9, 12
    out = [wrap("Estate-recipient diagnostic", styles, "Title"), Spacer(1, 7)]
    out += [wrap("Purpose and experimental design", styles, "Heading1"),
            wrap("The reference is the September 25 low-tax selected calibration. Its control, net-valuation, and net-valuation-plus-transfer cases hold preferences, B15 earnings, entrant wealth, other transfers/floors, housing supply and targets fixed. These experiments retain normalized population composition; they are not recalibrated or adopted specifications.", styles),
            wrap("The control values estates as max(b' + qh, 0); both treatments use net = max(b' + (1 - psi) qh, 0). Financial wealth b' already includes debt. Selling-cost psi is distinct from fertility preference psi_child. All cases retain the same entrant wealth distribution and fixed preferences. There are no age-18 transfers. In the receiving case, anticipated equal payments go to the four-year grid ages 46, 50, 54, 58, and 62.", styles),
            wrap("The inherited experimental payroll-tax rate is 8.751%, not a newly adopted fiscal rule. Population composition is normalized; the birth-entry gap is reported below and is not a certification of closed reproduction. The gross estate-flow observer and its target definition are unchanged across cases; net accounting is supplemental.", styles), Spacer(1, 8)]

    summary = [[wrap(h, styles, "TableHeader") for h in ("Case", "Loss", "Net estates / period", "Paid / period", "Transfer / recipient / period", "Transfer / year", "Recipient mass", "Gross estates / year", "Market residual", "Birth-entry gap", "House price")]]
    for case in CASES:
        receipt, targets = packets[case]["receipt"], packets[case]["targets"]
        accounts = receipt.get("estate_accounts", {})
        if not isinstance(accounts, dict):
            raise RuntimeError(f"{case}: estate_accounts must be an object")
        period_years = float(accounts["period_years"])
        transfer = float(accounts["transfer"])
        gross = float(accounts["generated_gross_period"]) / period_years
        summary.append([wrap(x, styles, "Small") for x in (case, fmt(receipt.get("loss")),
            fmt(accounts["generated_net_period"]), fmt(accounts["paid_period"]),
            fmt(transfer), fmt(transfer / period_years), fmt(accounts["recipient_mass"]),
            fmt(gross), fmt(receipt.get("market_residual")), fmt(receipt.get("fixed_psi_birth_replacement_residual")), fmt(packets[case]["price"]))])
    out += [wrap("Compact case summary", styles, "Heading1"), Table(summary, colWidths=[78, 52, 66, 54, 62, 58, 57, 58, 61, 62, 56], repeatRows=1, hAlign="LEFT"), Spacer(1, 4),
            wrap("Estate flows are per normalized household population; the transfer is per eligible recipient. All money amounts use model income units. Annual equivalents divide four-year amounts by four. Gross-flow measurement is unchanged; its level can respond to household choices. Equal payments exhaust net estates only in the receiving case. The other two cases retain an estate outflow.", styles, "Small"), PageBreak()]
    out[-4].setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), colors.HexColor("#1d4e70")), ("TEXTCOLOR", (0,0), (-1,0), colors.white), ("GRID", (0,0), (-1,-1), .25, colors.HexColor("#9aa6b2")), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3), ("TOPPADDING", (0,0), (-1,-1), 3), ("BOTTOMPADDING", (0,0), (-1,-1), 3)]))

    def table_pages(case: str, title: str, headers: tuple[str, ...], rows: list[list[str]], widths: list[float], chunk: int):
        for start in range(0, len(rows), chunk):
            heading = f"{case}: {title}" if case else title
            out.append(wrap(heading + ("" if start == 0 else " (continued)"), styles, "Heading1"))
            table = Table([[wrap(h, styles, "TableHeader") for h in headers]] + [[wrap(v, styles, "Small") for v in row] for row in rows[start:start+chunk]], colWidths=widths, repeatRows=1)
            table.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), colors.HexColor("#1d4e70")), ("TEXTCOLOR", (0,0), (-1,0), colors.white), ("GRID", (0,0), (-1,-1), .25, colors.HexColor("#aeb8c2")), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3), ("TOPPADDING", (0,0), (-1,-1), 3), ("BOTTOMPADDING", (0,0), (-1,-1), 3)]))
            out.extend([table, PageBreak()])
    comparison = []
    for index, control_row in enumerate(packets[CASES[0]]["targets"]):
        valuation_row = packets[CASES[1]]["targets"][index]
        receiving_row = packets[CASES[2]]["targets"][index]
        comparison.append([control_row["moment"], fmt(control_row["target"]),
            fmt(control_row["model"]), fmt(valuation_row["model"]), fmt(receiving_row["model"]),
            fmt(float(receiving_row["model"]) - float(valuation_row["model"]))])
    table_pages("", "All target moments: common units across cases",
        ("Moment", "Target", "Control", "Net valuation", "Net + receipts", "Receipt effect vs net valuation"),
        comparison, [220, 80, 80, 95, 95, 110], 13)
    for case in CASES:
        table_pages(case, "complete target fit", ("Moment", "Target", "Model", "Gap", "Weight", "Loss"),
                    [[r["moment"], fmt(r["target"]), fmt(r["model"]), fmt(r["gap"]), fmt(r["weight"]), fmt(r["loss_contribution"])] for r in packets[case]["targets"]], [250, 85, 85, 85, 80, 95], 13)
        table_pages(case, "all parameters, bounds, and restrictions", ("Parameter", "Estimate", "Lower", "Upper", "Near bound", "Restriction / status"),
                    [[r["parameter"], fmt(r["estimate"]), fmt(r["lower"]), fmt(r["upper"]), r["near_bound"], r["status"]] for r in packets[case]["parameters"]], [165, 75, 65, 65, 70, 240], 18)
    legend_note = ("Income-state legends are compactly formatted from the saved solutions, with plotted data preserved."
                   if packets[CASES[0]]["figure_legends_formatted"] else "Original figure legends are retained.")
    out += [wrap("Provenance and caveats", styles, "Heading1"),
            wrap("All values are read from completed receipts, target-fit CSVs, parameter CSVs, and saved standard diagnostics. Retained measurement approximations include the native first-birth housing observer versus the empirical event-study estimator, the recent-parent residence proxy, and all-positive estates versus the child-directed empirical target. The inherited loss does not decide which estate closure is appropriate. Transfers are deterministic and pooled; there is no parent-child matching or inheritance risk. These are fixed-preference effects under normalized entry composition, not recalibrated policy effects.", styles),
            wrap("Observed entry wealth is a stock containing prior saving and any prior assistance; later inheritance is a separate flow. This experiment funds the later receipts from net deaths, while entry wealth remains externally supplied. It does not close the full intergenerational wealth account or adopt the separate common-scale entry candidate. Modeling a pre-entry grant would require separating and funding that component of observed entry wealth, rather than adding it twice. A future transition must update the estate pool and recipient mass each date; changing cohort entry also changes aggregate entrant resources even with a fixed per-entrant wealth law.", styles),
            wrap(f"Run completion status: {complete.get('status')}. The stable figures are reproduced one image per landscape page. {legend_note} Input fingerprints and report counts are recorded in the adjacent JSON sidecar; visual review is recorded separately.", styles)]
    for key in ("source_manifest_sha256", "target_system_sha256", "selected_reference"):
        out.append(wrap(f"{key}: {packets[CASES[0]]['receipt'][key]}", styles, "Small"))
    out.append(PageBreak())
    scale_rows = []
    for case in CASES:
        receipt = packets[case]["receipt"]
        accounts = receipt["estate_accounts"]
        entrants = float(receipt["normalized_entrant_flow"])
        scale_rows.append([case, fmt(entrants),
            fmt(float(accounts["generated_gross_period"]) / entrants),
            fmt(float(accounts["generated_net_period"]) / entrants)])
    out.extend([wrap("Mechanical scale of an entry-only allocation", styles, "Heading1"),
        wrap("Dividing the saved period estate pool by the saved flow of entering households gives the one-time grant each entrant would receive if that pool were allocated entirely at entry. The grant is a wealth stock in units of average annual gross working earnings. It must not be divided by four. This calculation holds estate generation and entrant flow fixed; changing the recipient rule would also change saving, fertility and prices. No entry-only equilibrium is solved here.", styles), Spacer(1, 8)])
    allocation_table = Table(
        [[wrap(h, styles, "TableHeader") for h in ("Case supplying the pool", "Entrants / period", "Gross pool / entrant", "Net pool / entrant")]]
        + [[wrap(v, styles, "Small") for v in row] for row in scale_rows],
        colWidths=[230, 140, 150, 150], repeatRows=1)
    allocation_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#1d4e70")),
        ("GRID", (0,0), (-1,-1), .25, colors.HexColor("#aeb8c2")),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("TOPPADDING", (0,0), (-1,-1), 6), ("BOTTOMPADDING", (0,0), (-1,-1), 6)]))
    out.extend([allocation_table, PageBreak()])
    for case in CASES:
        for name in names:
            image_path = packets[case]["figures"] / name
            out.extend([wrap(f"{case}: {name}", styles, "Heading1"), Spacer(1, 4), Image(str(image_path), width=9.65*inch, height=6.5*inch, kind="proportional"), PageBreak()])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new PDF path; overwrite is refused")
    parser.add_argument("--figure-root", type=Path, help="completed display-only legend re-export")
    args = parser.parse_args()
    require_torch()
    run_root, output = args.run_root.resolve(), args.output.resolve()
    sidecar = output.with_suffix(".json")
    if output.exists() or sidecar.exists():
        raise RuntimeError("refusing to overwrite PDF or JSON sidecar")
    complete, packets, names, fingerprints = validate(run_root, args.figure_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    from reportlab.lib.pagesizes import landscape, letter
    from reportlab.platypus import SimpleDocTemplate
    document = SimpleDocTemplate(str(output), pagesize=landscape(letter), leftMargin=36, rightMargin=36, topMargin=34, bottomMargin=34)
    document.build(story_for(complete, packets, names), onFirstPage=footer, onLaterPages=footer)
    qa = render_and_check(output, packets)
    from pypdf import PdfReader
    parameter_table_pages = sum(math.ceil(len(packets[case]["parameters"]) / 18) for case in CASES)
    sidecar.write_text(json.dumps({"status": "report_built_unverified_visual_qa", "run_root": str(run_root), "pdf": str(output), "input_fingerprints_sha256": fingerprints, "case_count": 3, "target_table_pages": 3, "comparative_target_table_count": 1, "parameter_table_pages": parameter_table_pages, "target_parameter_table_page_count": 3 + parameter_table_pages, "supplemental_entry_allocation_table_count": 1, "figure_count": 51, "diagnostics_per_case": 17, "page_count": len(PdfReader(str(output)).pages), "complete_status": complete["status"], "pdf_qa":qa}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
