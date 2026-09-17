#!/usr/bin/env python3
"""Build the bounded E5F recovery review figures and PDF from a pinned JSON input.

This is an artist/reporting driver only.  It never imports the model and never
solves an equilibrium.  ``--mode plots`` writes PNGs and CSVs; ``--mode report``
assembles those artifacts and any available standard diagnostic PNGs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

DEFAULT_DIR = Path("output/model/e5f_original_queue_20260913a/long_successive_refit/recovery_review")
BLUE, ORANGE, GREY, RED = "#1f5fa6", "#d97815", "#777777", "#c73e3a"


def _f(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def _seq(x):
    if isinstance(x, (list, tuple)):
        return [_f(v) for v in x]
    return []


def _get(d, *names, default=float("nan")):
    for n in names:
        if isinstance(d, dict) and n in d and d[n] not in (None, ""):
            return _f(d[n], default)
    return default


def _input(path):
    with path.open() as fh:
        return json.load(fh)


def _rows(data):
    # The pinned recovery input deliberately keeps the native macro and
    # fertility panels as separate top-level arrays.
    return list(data.get("rows", []) or []), list(data.get("fertility", []) or [])


def _year(row, i):
    return _get(row, "window_start_year", "calendar_year", "year", "observation_year", default=2007 + 4 * i)


def _series(rows, *names):
    return [_get(r, *names) for r in rows]


def _finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def _initial_age_rates(data):
    f = data.get("initial", {}).get("fertility", {})
    rates = _seq(f.get("age_specific_birth_rate_topcode_adjusted"))
    if rates:
        return rates
    flow, mass = _seq(f.get("birth_flow_topcode_adjusted")), _seq(f.get("age_mass"))
    return [a / b if b else 0.0 for a, b in zip(flow, mass)]


def _completed(rows, data):
    """Reconstruct children ever born using the age-cell stock identity."""
    previous = []
    total = 0.0
    for rate in _initial_age_rates(data):
        total += rate
        previous.append(total)
    out = []
    for row in rows:
        rates = _seq(row.get("age_specific_birth_rate_topcode_adjusted"))
        if not rates:
            flow, mass = _seq(row.get("birth_flow_topcode_adjusted")), _seq(row.get("age_mass"))
            rates = [a / b if b else 0.0 for a, b in zip(flow, mass)]
        if not rates:
            out.append(float("nan")); continue
        current = [0.0] + previous[:-1]
        current = [a + b for a, b in zip(current, rates)]
        ages = _seq(row.get("age_cell_start"))
        try:
            idx = ages.index(42.0)
        except ValueError:
            idx = 6
        out.append(current[idx] if idx < len(current) else float("nan"))
        previous = current
    return out


def _save_csv(path, rows, fields=None):
    if not rows: return
    fields = fields or list(rows[0])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _source_hashes(data, input_path):
    result = {"review_input.json": hashlib.sha256(input_path.read_bytes()).hexdigest()}
    supplied = data.get("source_hashes", data.get("source_sha256", data.get("source", [])))
    if isinstance(supplied, dict):
        result.update({str(k): str(v) for k, v in supplied.items()})
        return result
    for item in supplied if isinstance(supplied, list) else []:
        if isinstance(item, dict): result.update({str(k): str(v) for k, v in item.items()})
    return result


def build_plots(data, out, input_path):
    import numpy as np
    import matplotlib.pyplot as plt

    out.mkdir(parents=True, exist_ok=True)
    macro, fert = _rows(data)
    years = np.array([_year(r, i) for i, r in enumerate(macro)])
    if not len(years): raise ValueError("review_input.json contains no native macro rows")
    fert_years = np.array([_year(r, i) for i, r in enumerate(fert or macro)])
    pf = _series(fert, "period_tfr_topcode_adjusted")
    completed = _completed(fert, data)
    initial_tfr=float(data['initial']['fertility']['period_tfr_topcode_adjusted'])
    initial_rates=np.asarray(_initial_age_rates(data))
    rates=np.asarray([f['age_specific_birth_rate_topcode_adjusted'] for f in fert])
    ages=np.asarray(fert[0]['age_cell_start']); k=int(np.flatnonzero(ages==42)[0])
    diagonal=np.asarray([sum(initial_rates[j] if t-k+j<0 else rates[t-k+j,j] for j in range(k+1)) for t in range(len(fert))])
    np.testing.assert_allclose(completed,diagonal,atol=2e-12,rtol=0)
    np.testing.assert_allclose(pf,rates.sum(axis=1),atol=2e-12,rtol=0)

    slide = data.get("slide_rows", []) or []
    slide_year = [_get(r, "window_start_year", "year", default=2007 + 4*i) for i, r in enumerate(slide)]

    def checked_plot(ax,x,y,*args,**kwargs):
        line,=ax.plot(x,y,*args,**kwargs)
        np.testing.assert_array_equal(line.get_xdata(),np.asarray(x))
        np.testing.assert_array_equal(line.get_ydata(),np.asarray(y))
        return line

    def style(ax):
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.grid(alpha=.15)

    def index(rows, names, base):
        return 100 * np.array([_get(r, *names) / base for r in rows])
    initial_q = data.get("initial", {}).get("quantities", {})
    terminal = data.get("terminal", {})
    adult0 = _get(initial_q, "adult_population", "population_households")
    housing0 = _get(initial_q, "housing_demand")
    price0 = _get(initial_q, "asset_price")
    term_tfr = _get(terminal, "birth_children_topcode_adjusted") / _get(terminal, "entry_flow")

    # 1. Four target points are observations; every later orange point is unfitted.
    targets = data.get("targets", []) or []
    ty = [_get(t, "year", "observation_year") for t in targets]
    tv = [_get(t, "target", "value") for t in targets]
    fig, ax = plt.subplots(figsize=(10, 5.8))
    if any(_finite(v) for v in tv): checked_plot(ax, ty, tv, "s--", color=BLUE, lw=2, ms=6, label="Observed target")
    checked_plot(ax, np.r_[2003,2007,fert_years], np.r_[initial_tfr,initial_tfr,pf], "o-", color=ORANGE, lw=2.2, ms=3.5, label="Unfitted forecast after first shock")
    if slide:
        checked_plot(ax, slide_year, _series(slide, "period_fertility"), "-", color=GREY, lw=1.8, label="Slide transition curve")
    ax.axhline(initial_tfr, color=BLUE, lw=1, ls=":", label="Pre-shock steady state (2.1)")
    ax.axvline(2007, color="black", lw=.8)
    ax.set(xlim=(2003,2065), xlabel="Start of four-year period", ylabel="Period fertility", title="Fertility fit and first-shock recovery")
    ax.legend(frameon=False, fontsize=9); style(ax); fig.tight_layout(); fig.savefig(out/"fertility_fit_2007_2063.png", dpi=180); plt.close(fig)

    # 2. The full panel is exclusively the 104-period native path.
    specs = [(np.array(pf), "Period fertility", term_tfr),
             (index(macro,("adult_population",),adult0), "Household heads (initial = 100)", 100*_get(terminal,"population_households")/adult0),
             (index(macro,("housing_demand",),housing0), "Total housing (initial = 100)", 100*_get(terminal,"housing_demand")/housing0),
             (index(macro,("asset_price",),price0), "House price (initial = 100)", 100*_get(terminal,"asset_price")/price0)]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2)); axes = axes.ravel()
    for ax,(vals,label,tvv) in zip(axes,specs):
        checked_plot(ax, fert_years if label == "Period fertility" else years, vals, "-", color=ORANGE, lw=1.3, label="Native transition")
        if _finite(tvv):
            checked_plot(ax, [2423], [tvv], "D", color=RED, ms=5, label="Terminal reference (2423)")
            ax.axhline(tvv, color=RED, lw=.7, ls=":")
        ax.set_ylabel(label); ax.set_xlabel("Start of four-year period"); style(ax)
    axes[0].legend(frameon=False, fontsize=8); fig.suptitle("Full transition after the first shock", y=.99); fig.tight_layout(); fig.savefig(out/"full_transition_four_panel.png", dpi=180); plt.close(fig)

    # 3. Short 2x2 inspection panel, with native levels normalized at 2007.
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True); axes=axes.ravel()
    short = years <= 2063; fyshort=fert_years <= 2063
    panels=[[(index(macro,("adult_population",),adult0),BLUE,"Households"),(index(macro,("housing_demand",),housing0),ORANGE,"Housing demand"),(100*index(macro,("housing_demand",),housing0)/index(macro,("adult_population",),adult0),GREY,"Housing per household")],
            [(index(macro,("asset_price",),price0),BLUE,"House price"),(index(macro,("renter_price",),_get(initial_q,"renter_price")),ORANGE,"Renter price")],
            [(np.array(pf),BLUE,"Period fertility"),(np.array(completed),ORANGE,"Completed fertility, age 42–45")],
            [(index(macro,("housing_demand",),housing0),BLUE,"Housing demand"),(index(macro,("housing_supply",),housing0),ORANGE,"Housing supply")]]
    for ax, lines in zip(axes,panels):
        for values,color,label in lines:
            mask=fyshort if len(values)==len(fert_years) else short; xx=fert_years if len(values)==len(fert_years) else years
            checked_plot(ax, xx[mask],values[mask],"o-",color=color,lw=1.5,ms=2.5,label=label)
        ax.set(xlim=(2007,2063)); ax.legend(frameon=False,fontsize=8); style(ax)
    axes[0].set_ylabel("Initial steady state = 100"); axes[1].set_ylabel("Initial steady state = 100"); axes[2].set_ylabel("Children"); axes[3].set_ylabel("Initial steady state = 100")
    axes[2].set_xlabel("Start of four-year period"); axes[3].set_xlabel("Start of four-year period"); fig.suptitle("Quantities, prices, and fertility, 2007–2063"); fig.tight_layout(); fig.savefig(out/"quantities_prices_2007_2063.png", dpi=180); plt.close(fig)

    # 4. Fiscal and market residuals; units are explicit and no common threshold asserted.
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for names,color,label in [(('relative_market_residual',),BLUE,"Market residual (%)"), (('scaled_pension_budget_residual',),ORANGE,"Pension budget residual (%)"), (('scaled_government_budget_residual',),GREY,"Government budget residual (%)")]:
        vals = [_get(r,*names) * 100 for r in macro]
        if any(_finite(v) for v in vals): checked_plot(ax, years, vals, "o-", ms=2.5, lw=1.5, color=color, label=label)
    ax.axhline(0,color="black",lw=.7); ax.set(xlabel="Start of four-year period", ylabel="Percent", title="Fiscal and market residuals"); ax.legend(frameon=False, fontsize=8); style(ax); fig.tight_layout(); fig.savefig(out/"fiscal_market_residuals.png", dpi=180); plt.close(fig)

    # 5. Lifecycle fertility: native rates only; no fabricated 2023 housing/wealth profiles.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    initial_rates = _initial_age_rates(data)
    post_rates = _seq(fert[0].get("age_specific_birth_rate_topcode_adjusted")) if fert else []
    rates_2023 = _seq(next((r.get("age_specific_birth_rate_topcode_adjusted") for r in fert if int(_year(r, 0)) == 2023), []))
    ages = _seq(data.get("initial",{}).get("fertility",{}).get("age_cell_start"))
    for label, vals, color in [("2007 pre-shock", initial_rates, BLUE), ("2007 post-shock", post_rates, ORANGE), ("2023", rates_2023, RED)]:
        if vals: checked_plot(ax, ages[:len(vals)], np.array(vals)/4*1000, "o-", ms=3, lw=1.5, label=label, color=color)
    ax.set(xlim=(18,46), xlabel="Age", ylabel="Births per 1,000 women per year (household analogue)", title="Lifecycle fertility profiles"); ax.legend(frameon=False, fontsize=9); style(ax); fig.tight_layout(); fig.savefig(out/"lifecycle_fertility_profiles.png", dpi=180); plt.close(fig)

    # Machine-readable artist data, preserving full active rows.
    fit_rows = []
    target_map={int(y):v for y,v in zip(ty,tv)}
    for i, y in enumerate(fert_years):
        target=target_map.get(int(y),""); cand=pf[i]
        fit_rows.append({"year":int(y),"target":target,"candidate_forecast":cand,"gap":cand-target if target != "" else "","fitted":"false","status":"future forecast; not fitted" if int(y)!=2007 else "candidate verified; fit missed"})
    _save_csv(out/"forecast.csv", fit_rows, ["year","target","candidate_forecast","gap","fitted","status"])
    _save_csv(out/"shock_fit.csv", [r for r in fit_rows if r["target"] != ""])
    _save_csv(out/"target_fit.csv", data.get("calibration_fit", []) or [])
    _save_csv(out/"parameters.csv", data.get("parameters", []) or [])
    # Independently check the cohort-diagonal recurrence and plotted artists.
    assert len(completed) == len(fert) and all(_finite(v) for v in completed)
    assert all(np.allclose(a, b, atol=2e-12, rtol=0) for a,b in [(np.asarray(pf), np.asarray(_series(fert,"period_tfr_topcode_adjusted"))), (np.asarray(index(macro,("adult_population",),adult0)), 100*np.asarray(_series(macro,"adult_population"))/adult0)])
    verification = {"status":"passed", "no_model_solves":True, "data_artist_checks": {"macro_rows":len(macro),"fertility_rows":len(fert),"completed_fertility_reconstructed":len(completed),"cohort_diagonals_validated":True,"artist_arrays_equal_native":True}, "source_hashes":_source_hashes(data,input_path), "input_hash":hashlib.sha256(input_path.read_bytes()).hexdigest(), "figure_list":["fertility_fit_2007_2063.png","full_transition_four_panel.png","quantities_prices_2007_2063.png","fiscal_market_residuals.png","lifecycle_fertility_profiles.png"]}
    (out/"verification.json").write_text(json.dumps(verification, indent=2)+"\n")
    return verification


def build_report(data, out):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak, KeepTogether
    from pypdf import PdfReader
    from PIL import Image as PILImage
    from xml.sax.saxutils import escape
    pdf = out / "recovery_review.pdf"
    width, height = landscape(A4)
    usable = width - 64
    doc = SimpleDocTemplate(str(pdf), pagesize=(width,height), rightMargin=26, leftMargin=26, topMargin=28, bottomMargin=28)
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small",fontName="Helvetica",fontSize=9,leading=12))
    styles.add(ParagraphStyle(name="HeadCell",fontName="Helvetica-Bold",fontSize=9,leading=12,textColor=colors.white))
    story=[]
    def text(value,style="BodyText"):
        return Paragraph(escape(str(value).replace("−","-")),styles[style])
    def title(value):
        story.extend([text(value,"Heading1"),Spacer(1,8)])
    def table(headers, rows, widths):
        body=[[text(v,"HeadCell") for v in headers]]
        body.extend([[text(v,"Small") for v in row] for row in rows])
        t=Table(body,colWidths=widths,repeatRows=1,hAlign="LEFT")
        t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor(BLUE)),("VALIGN",(0,0),(-1,-1),"TOP"),("BOTTOMPADDING",(0,0),(-1,-1),6),("TOPPADDING",(0,0),(-1,-1),6),("LINEBELOW",(0,0),(-1,0),.5,colors.grey),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f1f4f7")])]))
        story.append(t)
    def num(v):
        if v in (None,""):return "--"
        try:return format(float(v),".5g")
        except (ValueError,TypeError):return str(v)
    def picture(path,max_w,max_h):
        with PILImage.open(path) as im:w,h=im.size
        scale=min(max_w/w,max_h/h)
        return Image(str(path),width=w*scale,height=h*scale)
    title("Verified first-shock candidate: review packet")
    story.append(text("Saved results only. The 104-period forecast passes finite-horizon market and fiscal checks and reproduces exactly. It is not a completed historical fit: no shock meets both equilibrium and fertility-fit acceptance."))
    story.append(Spacer(1,9))
    story.append(text("Preferences fall from 0.1489153 to 0.1289153 in 2007 and are expected to remain there forever. The remaining historical shocks have not been estimated. Later orange observations are forecasts conditional on that single shock."))
    story.append(Spacer(1,12))
    fr={int(r["calendar_year"]):r for r in data["fertility"]}
    rows=[]
    for t in data["targets"]:
        v=fr[int(t["year"])]["period_tfr_topcode_adjusted"]
        rows.append([t["observation_years"],num(t["target"]),num(v),num(v-t["target"]),t["status"]])
    table(["Observed years","Data","Candidate forecast","Gap","Fit status"],rows,[92,68,110,75,usable-345])
    story.extend([Spacer(1,10),text("The grey comparator is the curve in the saved presentation PDF: four announced preference levels, still unconverged. Its old permanent-shock caption was inconsistent with the plotted data. The two experiments therefore differ in both shocks and expectations.","Small"),PageBreak()])
    figures=[("Fertility: data, slide curve and recovered candidate","fertility_fit_2007_2063.png"),("Full transition and verified stationary endpoint","full_transition_four_panel.png"),("Quantities, prices and reconstructed completed fertility","quantities_prices_2007_2063.png"),("Market and fiscal residuals","fiscal_market_residuals.png"),("Fertility over the lifecycle","lifecycle_fertility_profiles.png")]
    for heading,name in figures:
        title(heading);story.append(picture(out/name,usable,430))
        story.append(text("Candidate forecast after one shock; later shocks remain unfitted. Four-year model periods. Source: verified native rows and fertility observations.","Small"));story.append(PageBreak())
    title("Initial parameter vector: complete empirical review")
    story.append(text("All 13 rows are retained. These describe the unchanged initial state, not a new calibration. The current provenance review labels the scored restrictions proposed rather than certified SMM estimates. The 2.1 normalization is separate and unscored; -- denotes an unavailable/not-applicable weight or contribution.","Small"));story.append(Spacer(1,8))
    fits=data["calibration_fit"]
    table(["Moment","Target","Model","Gap","Weight","Contribution"],[[r["label"],num(r["target"]),num(r["model"]),num(r["gap"]),num(r["actual_weight"]),num(r["loss_contribution"])] for r in fits],[272,75,75,75,105,usable-602])
    story.append(PageBreak())
    params=data["parameters"]
    for offset in range(0,len(params),9):
        title("Parameters and restrictions"+(" (continued)" if offset else ""))
        table(["Parameter","Value","Lower","Upper","Near bound","Role / restriction"],[[r["parameter"],num(r["estimate"]),num(r["lower"]),num(r["upper"]),r["near_bound"],r["interpretation"]] for r in params[offset:offset+9]],[165,75,65,65,72,usable-442])
        story.append(text("Annual beta is estimated subject to a 0.99 upper bound. Full source status, raw bounds and provenance are preserved in parameters.csv; this packet does not promote these values as new estimates.","Small"));story.append(PageBreak())
    title("Validation status and saved-data coverage")
    td=data["terminal_distance"]
    table(["Check / object","Result"],[["Finite market/fiscal equilibrium","Passed for this first candidate; exact final reproduction."],["Historical shock fit","Zero of four accepted. First gap exceeds 0.005 tolerance."],["Verified terminal stationary equilibrium","Passed. This does not certify the path has reached it."],["Terminal household mass gap",f"{100*td['population_relative_gap']:.3f}%"],["Terminal distribution relative L1 distance",f"{100*td['distribution_relative_l1']:.3f}%"],["Horizon robustness","Not established."],["2023 fertility by age","Available from the candidate forecast."],["2023 housing, ownership and wealth by age; intergenerational allocation","Not saved for this candidate; would require additional household-state extraction/replay."],["New policy comparison","Not run: historical fit did not reach 2023."]],[275,usable-275])
    story.append(PageBreak())
    diag=sorted((out/"source/recovered_candidate/round_01/graphs/standard_diagnostics").glob("*.png"))
    assert len(diag)==17,"Expected unchanged standard 17-plot diagnostic set"
    for i in range(0,len(diag),2):
        title("Household diagnostics: 2007 after the first shock")
        story.append(text("Saved standard plots. The legacy file lifecycle_2023.csv belongs to this 2007 state. Expected-choice statistics in the legacy fertility panels are distinct from the measured all-birth age rates above; summary.json's tfr field measures completed fertility.","Small"));story.append(Spacer(1,14))
        cells=[[text(p.stem.replace("_"," "),"Small"),picture(p,(usable-30)/2,310)] for p in diag[i:i+2]]
        if len(cells)==1:cells.append("")
        story.append(Table([cells],colWidths=[usable/2,usable/2],style=[("VALIGN",(0,0),(-1,-1),"TOP")]))
        if i+2<len(diag):story.append(PageBreak())
    def footer(canvas,document):
        canvas.setFont("Helvetica",8);canvas.drawRightString(width-32,15,f"{document.page}")
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    pages=len(PdfReader(str(pdf)).pages)
    if pages>25:raise ValueError(f"Unexpected report length: {pages}")
    return pdf,pages


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--mode", choices=("plots","report","all"), required=True); ap.add_argument("--input", type=Path); ap.add_argument("--output-dir", type=Path, default=DEFAULT_DIR); args=ap.parse_args()
    out=args.output_dir; inp=args.input or out/"review_input.json"
    data=_input(inp)
    if args.mode in ("plots","all"):
        build_plots(data,out,inp)
        if args.mode == "all":
            import subprocess
            report_python=Path.home()/".cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"
            subprocess.run([str(report_python),str(Path(__file__).resolve()),"--mode","report","--input",str(inp.resolve()),"--output-dir",str(out.resolve())],check=True)
    else:
        build_report(data,out)
        # report mode remains non-plotting and records its own verification marker.
        meta=out/"verification.json"; old=json.loads(meta.read_text()) if meta.exists() else {}
        old.update({"report_built":True,"report_path":"recovery_review.pdf","no_model_solves":True}); meta.write_text(json.dumps(old,indent=2)+"\n")
    print(out)


if __name__ == "__main__": main()
