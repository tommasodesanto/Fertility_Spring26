"""Build the saved-array credit distribution report; no model solve is run.

Usage (from repository root):
  python output/model/experiments/credit_distribution_report_v1/build_pdf.py
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A3, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph


ROOT = Path(__file__).resolve().parents[4]
WORK = Path(__file__).resolve().parent
MANIFEST = WORK / "manifest.json"
OUTPUT = ROOT / "output/pdf/credit_and_fertility_distributions.pdf"
PAGE_W, PAGE_H = landscape(A3)
LEFT = 55
RIGHT = PAGE_W - 55
NAVY = colors.HexColor("#173047")
TEAL = colors.HexColor("#007F82")
AMBER = colors.HexColor("#BB6D27")
PALE = colors.HexColor("#ECF3F3")
FAINT = colors.HexColor("#F5F7F8")
GREY = colors.HexColor("#506171")
LIGHT = colors.HexColor("#CBD5D9")
RED = colors.HexColor("#AB4D43")


def path_from(raw: str | Path) -> Path:
    p = Path(raw)
    return p if p.is_absolute() else ROOT / p


def section_of(item: dict) -> str:
    val = str(item.get("section", item.get("population", item.get("group", "reference")))).lower()
    return "low" if "low" in val else "reference"


class Report:
    def __init__(self, manifest: dict):
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        self.c = canvas.Canvas(str(OUTPUT), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
        self.c.setTitle("Credit, housing and family size")
        self.c.setAuthor("Fertility Spring 2026 research project")
        self.manifest = manifest
        self.page = 0
        self.toc: list[tuple[str, int]] = []
        self.chapter = ""

    def paragraph(self, text: str, x: float, top: float, width: float, size=11, leading=None,
                  color=GREY, font="Helvetica", max_height=500):
        style = ParagraphStyle(
            "body", fontName=font, fontSize=size, leading=leading or size * 1.42,
            textColor=color, alignment=TA_LEFT, spaceAfter=0,
        )
        para = Paragraph(text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"), style)
        _, height = para.wrap(width, max_height)
        if height > max_height + .01:
            raise ValueError(f"Paragraph does not fit: {text[:90]}")
        para.drawOn(self.c, x, top - height)
        return top - height

    def line(self, y, x=LEFT, right=RIGHT, color=LIGHT, width=.6):
        self.c.setStrokeColor(color)
        self.c.setLineWidth(width)
        self.c.line(x, y, right, y)

    def start(self, title: str, chapter: str, kicker="RESULTS", book=False):
        if self.page:
            self.c.showPage()
        self.page += 1
        self.chapter = chapter
        self.c.setFillColor(NAVY)
        self.c.rect(0, PAGE_H - 9, PAGE_W, 9, fill=1, stroke=0)
        self.c.setFont("Helvetica-Bold", 8.5)
        self.c.setFillColor(TEAL)
        self.c.drawString(LEFT, PAGE_H - 47, kicker.upper())
        self.c.setFont("Helvetica-Bold", 22)
        self.c.setFillColor(NAVY)
        self.c.drawString(LEFT, PAGE_H - 78, title)
        self.line(PAGE_H - 95)
        self.line(52)
        self.c.setFont("Helvetica", 8)
        self.c.setFillColor(GREY)
        self.c.drawString(LEFT, 34, "CREDIT, HOUSING AND FAMILY SIZE  /  SAVED FIXED-PRICE COMPARISON")
        self.c.drawRightString(RIGHT, 34, f"{self.page:02d}")
        if book:
            self.toc.append((title, self.page))
        return PAGE_H - 116

    def badge(self, x, y, text, fill=PALE, ink=TEAL, pad=9):
        w = stringWidth(text, "Helvetica-Bold", 9) + 2 * pad
        self.c.setFillColor(fill)
        self.c.roundRect(x, y - 4, w, 24, 7, fill=1, stroke=0)
        self.c.setFont("Helvetica-Bold", 9)
        self.c.setFillColor(ink)
        self.c.drawString(x + pad, y + 5, text)
        return x + w

    def image(self, path: Path, x: float, y_top: float, max_w: float, max_h: float):
        if not path.is_file():
            raise FileNotFoundError(path)
        im = PILImage.open(path)
        iw, ih = im.size
        scale = min(max_w / iw, max_h / ih)
        w, h = iw * scale, ih * scale
        self.c.drawImage(str(path), x + (max_w - w) / 2, y_top - h, w, h,
                         preserveAspectRatio=True, mask="auto")
        return y_top - h

    def cover(self):
        self.start("Credit, housing and family size", "Cover", "research diagnostic")
        self.c.setFillColor(PALE)
        self.c.roundRect(LEFT, 372, RIGHT - LEFT, 312, 16, fill=1, stroke=0)
        self.c.setFont("Helvetica-Bold", 12)
        self.c.setFillColor(TEAL)
        self.c.drawString(LEFT + 34, 642, "FOUR SAVED HOUSEHOLD CASES")
        self.c.setFont("Helvetica-Bold", 31)
        self.c.setFillColor(NAVY)
        self.c.drawString(LEFT + 34, 590, "Financing changes at held prices")
        self.paragraph(
            "The main comparison changes the financed share of a home purchase from 80% to 100% "
            "in the current working chain-13 population. A separate extension repeats the comparison "
            "when every household remains at the lowest productivity state.",
            LEFT + 34, 545, RIGHT - LEFT - 68, size=16, leading=23, color=NAVY,
        )
        x = LEFT + 34
        x = self.badge(x, 418, "80% vs 100% financed") + 12
        x = self.badge(x, 418, "Two populations") + 12
        self.badge(x, 418, "17 standard views per case")
        self.paragraph("October 3, 2026  |  Prepared from saved native solution arrays and plotted diagnostics",
                       LEFT, 336, RIGHT - LEFT, 11, color=GREY)
        self.line(305)
        self.paragraph(
            "Interpretation boundary: these are fixed-price household solutions. A 100% financed share "
            "removes the modeled purchase down payment; it does not make all borrowing unrestricted. "
            "The report does not claim market clearing, a recalibration, or a policy general equilibrium.",
            LEFT, 276, RIGHT - LEFT, 12, leading=18, color=NAVY,
        )
        self.c.setFont("Helvetica-Bold", 10)
        self.c.setFillColor(TEAL)
        self.c.drawString(LEFT, 181, "CONTENTS")
        self.line(169)
        self.paragraph(
            "Scope and definitions  /  Current productivity distributions  /  Permanent lowest-productivity "
            "extension  /  Numeric comparison tables  /  Sources and limitations  /  Seventeen four-case diagnostic views",
            LEFT, 151, RIGHT - LEFT, 10.7, leading=16, color=NAVY,
        )

    def guide(self):
        y = self.start("How to read this report", "Guide", "scope and definitions", True)
        self.c.setFont("Helvetica-Bold", 13)
        self.c.setFillColor(NAVY)
        self.c.drawString(LEFT, y - 15, "Experiment matrix")
        rows = [
            ("Current productivity process", "80%", "100%", "Main comparison"),
            ("Permanent lowest productivity", "80%", "100%", "Separate extension"),
        ]
        cols = [LEFT, LEFT + 425, LEFT + 585, LEFT + 735]
        top = y - 47
        for i, row in enumerate(rows):
            by = top - i * 52
            self.c.setFillColor(PALE if i % 2 == 0 else FAINT)
            self.c.roundRect(LEFT, by - 34, RIGHT - LEFT, 43, 6, fill=1, stroke=0)
            for j, cell in enumerate(row):
                self.c.setFont("Helvetica-Bold" if j == 0 else "Helvetica", 11)
                self.c.setFillColor(NAVY if j == 0 else GREY)
                self.c.drawString(cols[j] + 12, by - 15, cell)
        y = top - 133
        definitions = [
            ("Financed share", "The fraction of a home purchase financed by mortgage debt. The down payment uses 1 minus this share."),
            ("Children ever born", "A stock of births up to the model age. It differs from the number of children currently at home."),
            ("Ownership share", "The population mass choosing owner tenure divided by total mass at the age or group shown."),
            ("Age", "The start of a four-year model cell for the reproductive household member, not a single-year survey interview age."),
            ("Fixed price", "The saved house price is held while household choices and distributions change; housing markets need not clear."),
        ]
        for term, definition in definitions:
            self.c.setFont("Helvetica-Bold", 11)
            self.c.setFillColor(TEAL)
            self.c.drawString(LEFT + 4, y, term.upper())
            self.paragraph(definition, LEFT + 180, y + 6, RIGHT - LEFT - 188, 11, color=NAVY)
            y -= 62
        self.paragraph("The appendix preserves all 17 standard diagnostic plot types for each of the four saved cases. "
                       "Raw policy lines can include unoccupied numerical income and wealth states. "
                       "Unvalidated 100% financing consumption and next-asset aggregates are omitted from the main comparison.",
                       LEFT, y - 5, RIGHT - LEFT, 10.5)
        self.line(166)
        self.paragraph(
            "Executed inputs: house price 0.776057 in all four cases. The lower-productivity extension "
            "sets permanent productivity to 0.103468, retains the entrant financial-wealth marginal, "
            "keeps the payroll-tax rate at 0.080281 and rebalances the four-year pension to 0.094961.",
            LEFT, 151, RIGHT - LEFT, 9.5, leading=13, color=NAVY, max_height=77,
        )

    def main_figure_page(self, items: list[dict], heading: str, chapter: str, idx: int):
        y = self.start(heading, chapter, f"{chapter}  /  comparison {idx}", idx == 1)
        w = (RIGHT - LEFT - 32) / 2
        if len(items) == 1:
            fig = items[0]
            title = str(fig.get("title", Path(fig["path"]).stem.replace("_", " ")))
            self.c.setFillColor(TEAL if chapter == "Reference" else AMBER)
            self.c.setFont("Helvetica-Bold", 13)
            self.c.drawString(LEFT, y - 11, title[:76])
            image_bottom = self.image(path_from(fig["path"]), LEFT, y - 29,
                                      RIGHT - LEFT, 490)
            cap = str(fig.get("caption", ""))
            if cap:
                self.paragraph(cap, LEFT, image_bottom - 14, RIGHT - LEFT,
                               11, leading=15.2, max_height=94)
            return
        for j, fig in enumerate(items):
            x = LEFT + j * (w + 32)
            title = str(fig.get("title", Path(fig["path"]).stem.replace("_", " ")))
            self.c.setFillColor(TEAL if chapter == "Reference" else AMBER)
            self.c.setFont("Helvetica-Bold", 13)
            self.c.drawString(x, y - 11, title[:76])
            image_bottom = self.image(path_from(fig["path"]), x, y - 28, w, 455)
            cap = str(fig.get("caption", ""))
            if cap:
                self.paragraph(cap, x, min(355, image_bottom - 15), w, 10.5,
                               leading=14.6, max_height=115)

    def note_page(self, title: str, chapter: str, notes: list[tuple[str, str]]):
        y = self.start(title, chapter, "interpretation and checks", True)
        for label, body in notes:
            self.c.setFillColor(PALE)
            self.c.roundRect(LEFT, y - 98, RIGHT - LEFT, 94, 10, fill=1, stroke=0)
            self.c.setFont("Helvetica-Bold", 12)
            self.c.setFillColor(TEAL)
            self.c.drawString(LEFT + 18, y - 25, label.upper())
            self.paragraph(body, LEFT + 18, y - 37, RIGHT - LEFT - 36, 11,
                           leading=15, color=NAVY, max_height=59)
            y -= 115
            if y < 135:
                y = self.start(title + " (continued)", chapter, "interpretation and checks")

    def table_pages(self, entry: dict):
        src = path_from(entry["path"])
        with src.open(newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            rows = list(reader)
        if not rows:
            return
        title = entry.get("title", src.stem.replace("_", " ").title())
        headers, data = rows[0], rows[1:]
        n = len(headers)
        if n == 0:
            return
        # Wide tables remain complete; small type is reserved for dense numerical source tables.
        col_w = (RIGHT - LEFT) / n
        fsize = min(9, max(5.8, 105 / max(12, n)))
        row_h = 20 if n <= 12 else 18
        head_h = 50
        per_page = max(1, int((PAGE_H - 185 - head_h) // row_h))
        if not data:
            data = [[]]
        for start in range(0, len(data), per_page):
            end = min(start + per_page, len(data))
            y = self.start(title + (f"  ({start+1}-{end})" if len(data) > per_page else ""),
                           "Source tables", "selected source values", start == 0)
            self.paragraph(f"Source: {src.relative_to(ROOT)}. More granular age, wealth and room tables "
                           "are retained as companion CSVs in the same analysis folder.", LEFT, y, RIGHT - LEFT, 9.5)
            y -= 57
            self.c.setFillColor(NAVY)
            self.c.roundRect(LEFT, y - head_h + 2, RIGHT - LEFT, head_h, 4, fill=1, stroke=0)
            for j, h in enumerate(headers):
                box_x = LEFT + j * col_w + 4
                words = str(h).replace("_", " ").split()
                line1 = ""
                lines = []
                for word in words:
                    candidate = f"{line1} {word}".strip()
                    if stringWidth(candidate, "Helvetica-Bold", fsize) > col_w - 9 and line1:
                        lines.append(line1)
                        line1 = word
                    else:
                        line1 = candidate
                if line1:
                    lines.append(line1)
                self.c.setFillColor(colors.white)
                self.c.setFont("Helvetica-Bold", fsize)
                for k, ln in enumerate(lines[:4]):
                    self.c.drawString(box_x, y - 13 - k * (fsize + 2), ln)
            yy = y - head_h
            for i, row in enumerate(data[start:end]):
                self.c.setFillColor(FAINT if i % 2 == 0 else colors.white)
                self.c.rect(LEFT, yy - row_h + 2, RIGHT - LEFT, row_h, stroke=0, fill=1)
                for j in range(n):
                    val = str(row[j]) if j < len(row) else ""
                    if not val:
                        val = "n/a"
                    if val.startswith("unavailable_"):
                        val = "Unavailable"
                    if val:
                        # Retain all precision in source CSV; show a readable table value.
                        try:
                            v = float(val)
                            val = f"{v:.6g}" if math.isfinite(v) else val
                        except ValueError:
                            pass
                        while stringWidth(val, "Helvetica", fsize) > col_w - 8 and fsize >= 5.8:
                            fsize -= .2
                        if stringWidth(val, "Helvetica", fsize) > col_w - 8:
                            val = val[:max(1, int((col_w - 8) / (fsize * .52)) - 1)] + "..."
                    self.c.setFillColor(NAVY)
                    self.c.setFont("Helvetica", fsize)
                    self.c.drawString(LEFT + j * col_w + 4, yy - row_h + 7, val)
                yy -= row_h
            self.line(yy)

    def appendix(self, case_dirs: list[tuple[str, Path]]):
        standard_names = [
            "fertility_by_age", "ownership_by_age", "tenure_services", "owner_rungs",
            "income_state_outcomes", "housing_by_age_income_state", "ownership_by_age_income_state",
            "liquid_wealth_by_age_income_state", "fertility_policy_by_age_income_state",
            "wealth_dist_childless_renter_age30", "wealth_dist_childless_renter_age42",
            "policy_childless_renter_age30", "policy_childless_renter_age42",
            "housing_prices", "housing_market", "market_clearing_by_market",
            "market_clearing_residuals",
        ]
        labels = {
            "fertility_by_age": "Children and births by age",
            "ownership_by_age": "Homeownership by age",
            "tenure_services": "Tenure and housing services",
            "owner_rungs": "Owner housing products",
            "income_state_outcomes": "Outcomes by income state",
            "housing_by_age_income_state": "Housing by age and income state",
            "ownership_by_age_income_state": "Ownership by age and income state",
            "liquid_wealth_by_age_income_state": "Financial wealth by age and income state",
            "fertility_policy_by_age_income_state": "Fertility choices by age and income state",
            "wealth_dist_childless_renter_age30": "Wealth: childless renters at age 30",
            "wealth_dist_childless_renter_age42": "Wealth: childless renters at age 42",
            "policy_childless_renter_age30": "Policies: childless renters at age 30",
            "policy_childless_renter_age42": "Policies: childless renters at age 42",
            "housing_prices": "Held housing prices",
            "housing_market": "Housing quantities at the held price",
            "market_clearing_by_market": "Market residuals by market",
            "market_clearing_residuals": "Market residuals",
        }
        box_w = (RIGHT - LEFT - 24) / 2
        box_h = 278
        for i, name in enumerate(standard_names):
            y = self.start(labels[name], "Standard diagnostics", f"appendix  /  view {i+1:02d} of 17", i == 0)
            for j, (label, directory) in enumerate(case_dirs):
                col = j % 2
                row = j // 2
                x = LEFT + col * (box_w + 24)
                top = y - row * 304
                self.c.setFillColor(PALE if row == 0 else colors.HexColor("#F9F2EA"))
                self.c.roundRect(x, top - box_h, box_w, box_h, 7, fill=1, stroke=0)
                self.c.setFont("Helvetica-Bold", 10)
                self.c.setFillColor(TEAL if row == 0 else AMBER)
                self.c.drawString(x + 10, top - 17, label)
                self.image(directory / f"{name}.png", x + 6, top - 25, box_w - 12, box_h - 36)
            note = "Fixed-price diagnostic; plotted market residuals are not equilibrium clearing tests."
            if name in {"fertility_policy_by_age_income_state", "policy_childless_renter_age30",
                        "policy_childless_renter_age42", "income_state_outcomes",
                        "housing_by_age_income_state", "ownership_by_age_income_state",
                        "liquid_wealth_by_age_income_state"}:
                note += " Raw policies may show unoccupied numerical states."
            note += " Both 100% aggregate policy reporters failed their unchanged value gate; raw diagnostics are shown, not certified aggregate consumption or saving."
            self.paragraph(note, LEFT, 101, RIGHT - LEFT, 8.5, leading=11, max_height=41)

    def finish(self):
        self.c.save()
        return OUTPUT, self.page


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Data manifest not ready: {MANIFEST}")
    manifest = json.loads(MANIFEST.read_text())
    report = Report(manifest)
    report.cover()
    report.guide()
    reference = [x for x in manifest.get("figures", []) if section_of(x) == "reference"]
    low = [x for x in manifest.get("figures", []) if section_of(x) == "low"]
    for section, figures, heading in [
        ("Reference", reference, "Current productivity process"),
        ("Lower productivity", low, "Permanent lowest productivity"),
    ]:
        for i, figure in enumerate(figures, 1):
            report.main_figure_page([figure], heading, section, i)
    report.note_page("What the saved cases establish", "Interpretation", [
        ("Current productivity process",
         "The all-age ownership share rises from 67.927% under 80% financing to 77.404% under 100% financing. At ages 42-45, the childless share is 18.888% and 19.611%, respectively. Credit relaxation changes ownership substantially while the children-ever-born comparison is smaller and has a different sign. These are saved household distributions at the same housing price."),
        ("Permanent lowest productivity",
         "The all-age ownership share rises from 0.077% to 3.792%, while age-22, age-30 and age-42 ownership remains near 0.08% in both financing arms. Virtually all households have no children ever born at every saved age; the largest share with any children is 1.73e-12. Productivity and the balanced pension also differ from the current population."),
        ("Numerical qualification",
         "Both 100% financed native distributions pass their distribution checks, but the unchanged aggregate policy reporter stops at infeasible buyer-value cells. The flagged population masses are 5.26e-38 in the current-productivity case and 2.68e-15 in the lower-productivity case. Their consumption and next-asset aggregates are therefore unavailable. The raw distribution still supports the children, ownership, room and financial-wealth figures."),
        ("Economic scope",
         "The four cases retain the working post-interest chain-13 parameters, current birth menu, no Estate A and fixed house prices. There is no price adjustment, new calibration, general-equilibrium claim or population forecast."),
    ])
    for table in manifest.get("tables", []):
        report.table_pages(table)
    report.note_page("Sources, reproduction and limits", "Provenance", [
        ("Saved case identities",
         "Main reference: output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/. Lower productivity: output/model/experiments/low_productivity_credit_v1/. Each case retains its executed parameters, solution arrays and standard diagnostics."),
        ("No new solve",
         "The accompanying output/model/experiments/credit_distribution_report_v1/manifest.json records each included figure and table. Granular CSVs remain beside the figures. Regenerate this PDF with build_pdf.py in the same folder, using the bundled Python environment. Figure and table preparation reads saved arrays only."),
        ("Comparability",
         "Within each population, the 80% and 100% financing arms hold the other executed household inputs and the house price fixed. Across populations, productivity, entrant wealth-income allocation and balanced pension differ. The latter comparison therefore cannot be attributed to financing alone."),
        ("Raw diagnostic scope",
         "The appendix reproduces the standard diagnostic images for auditability. Several plot titles come from a general-equilibrium diagnostic routine, but their market residual panels here are evaluations at held prices. Income-state policy lines in the lower-productivity population include eight unoccupied numerical nodes."),
    ])
    refbase = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation"
    lowbase = ROOT / "output/model/experiments/low_productivity_credit_v1"
    dirs = [
        ("Current productivity  /  80% financed", refbase / "phi_080/standard_diagnostics"),
        ("Current productivity  /  100% financed", refbase / "phi_100/standard_diagnostics"),
        ("Lowest productivity  /  80% financed", lowbase / "phi_08/standard_diagnostics"),
        ("Lowest productivity  /  100% financed", lowbase / "phi_10/standard_diagnostics"),
    ]
    report.appendix(dirs)
    output, pages = report.finish()
    print(json.dumps({"pdf": str(output), "pages": pages, "figures": len(reference) + len(low),
                      "tables": len(manifest.get("tables", []))}, indent=2))


if __name__ == "__main__":
    main()
