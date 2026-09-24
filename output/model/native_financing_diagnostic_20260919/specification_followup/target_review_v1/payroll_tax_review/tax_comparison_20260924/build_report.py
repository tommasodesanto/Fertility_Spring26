#!/usr/bin/env python3
"""Build the reviewed two-rate diagnostic PDF from collected compact artifacts."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_first"
PDF = Path(__file__).resolve().parents[7] / "output/pdf/paygo_tax_comparison_20260924.pdf"
W, H = 792, 612
M = 42
INK = colors.HexColor("#152238")
BLUE = colors.HexColor("#175c9e")
GRAY = colors.HexColor("#526275")
LINE = colors.HexColor("#d9e1e8")


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def num(value, *, percent=False):
    if value in (None, ""):
        return "-"
    x = float(value)
    if percent:
        return f"{100*x:.3f}%"
    if x and abs(x) < .001:
        return f"{x:.2e}"
    return f"{x:,.3f}"


def text(c, x, y, value, *, size=9, color=INK, font="Helvetica"):
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawString(x, y, value)


def right(c, x, y, value, *, size=9, color=INK, font="Helvetica"):
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawRightString(x, y, value)


def footer(c, page):
    c.setStrokeColor(LINE)
    c.line(M, 30, W-M, 30)
    text(c, M, 17, "Experimental fixed-parameter PAYGO comparison | September 24, 2026", size=7.5, color=GRAY)
    right(c, W-M, 17, f"{page}", size=7.5, color=GRAY)


def heading(c, title, subtitle=None):
    text(c, M, H-44, title, size=20, font="Helvetica-Bold")
    if subtitle:
        text(c, M, H-62, subtitle, size=9, color=GRAY)
    c.setStrokeColor(LINE)
    c.line(M, H-74, W-M, H-74)


def table(c, x, y_top, widths, header, body, row_h=20, sizes=None):
    total = sum(widths)
    c.setFillColor(BLUE)
    c.roundRect(x, y_top-row_h, total, row_h, 3, stroke=0, fill=1)
    xx=x
    for i,h in enumerate(header):
        text(c,xx+6,y_top-row_h+6,h,size=8,color=colors.white,font="Helvetica-Bold")
        xx+=widths[i]
    y=y_top-row_h
    for j,row in enumerate(body):
        y-=row_h
        if j%2==0:
            c.setFillColor(colors.HexColor("#f5f8fb"))
            c.rect(x,y,total,row_h,stroke=0,fill=1)
        xx=x
        for i,value in enumerate(row):
            value=str(value)
            if i==0:
                text(c,xx+6,y+6,value,size=(sizes or {}).get(i,8))
            else:
                right(c,xx+widths[i]-6,y+6,value,size=(sizes or {}).get(i,8))
            xx+=widths[i]
        c.setStrokeColor(LINE)
        c.line(x,y,x+total,y)
    return y


def draw_fitted_image(c, path, x, y, box_w, box_h):
    image=ImageReader(str(path))
    iw,ih=image.getSize()
    scale=min(box_w/iw,box_h/ih)
    dw,dh=iw*scale,ih*scale
    c.drawImage(image,x+(box_w-dw)/2,y+(box_h-dh)/2,width=dw,height=dh,preserveAspectRatio=True,mask='auto')


def main():
    complete=json.loads((RESULTS/"complete.json").read_text())
    assert complete["status"]=="completed" and complete["solves"]==2
    comparison={r["moment"]:r for r in rows(ROOT/"comparison.csv")}
    targets=rows(ROOT/"full_target_comparison.csv")
    params_a=rows(RESULTS/"current_179/full_parameters_actual_bounds.csv")
    params_b=rows(RESULTS/"proposal_087510/full_parameters_actual_bounds.csv")
    assert len(targets)==13 and len(params_a)==len(params_b)==17
    graphs_a={p.name:p for p in (RESULTS/"current_179/standard_diagnostics").glob("*.png")}
    graphs_b={p.name:p for p in (RESULTS/"proposal_087510/standard_diagnostics").glob("*.png")}
    assert len(graphs_a)==len(graphs_b)==17 and set(graphs_a)==set(graphs_b)
    PDF.parent.mkdir(parents=True,exist_ok=True)
    c=canvas.Canvas(str(PDF),pagesize=(W,H),pageCompression=1)
    c.setTitle("PAYGO payroll tax: fixed-parameter comparison")
    c.setAuthor("Fertility_Spring26 diagnostic")
    page=1
    heading(c,"PAYGO payroll tax: two-rate diagnostic","Same selected parameters and utility scale; endogenous balanced pension and housing equilibrium")
    keys=[
        ("Payroll rate","payroll_tax_rate","percent"),
        ("Pension / mean annual gross earnings","pension_annual_over_mean_working_gross","percent"),
        ("Payroll revenue (four-year units)","payroll_revenue_period",None),
        ("Pension outlays (four-year units)","pension_outlays_period",None),
        ("Mean working disposable income","mean_working_annual_disposable",None),
        ("Aggregate wealth","aggregate_wealth",None),
        ("Ownership, ages 30-55","ownership_30_55","percent"),
        ("Mean occupied rooms","mean_rooms_capped9_18_85",None),
        ("House price","house_price",None),
        ("Childlessness, ages 40-44","childlessness_40_44","percent"),
        ("Mean age at first birth","first_birth_mean_age",None),
        ("Top-coded children-ever-born stock","model_tfr",None),
        ("Raw stationary children-state mean","mean_completed_fertility_legacy",None),
    ]
    body=[]
    for label,key,kind in keys:
        row=comparison[key]
        a=float(row["current_179"]); b=float(row["proposal_087510"])
        body.append([label,num(a,percent=kind=="percent"),num(b,percent=kind=="percent"),
                     f"{100*(b-a):+.3f} pp" if kind=="percent" else f"{b-a:+,.3f}"])
    y=table(c,M,H-91,[310,128,145,125],["Outcome","17.9% rate","8.751% rate","Change"],body,row_h=26,sizes={0:8.5})
    text(c,M,y-24,"Finding: the proposal halves the pension and raises wealth, ownership and fertility at fixed preferences.",size=8.6,font="Helvetica-Bold")
    text(c,M,y-40,"It cannot select a tax rate: later adopted targets are not applied, and demographic entry is fixed.",size=8.5,color=GRAY)
    text(c,M,y-55,"The 2.1 top-coded stock moves to 2.252. It equals adjusted birth flow / fixed entry flow here.",size=8.5,color=GRAY)
    text(c,M,y-70,"This is not a new demographic steady state, period total fertility rate, or cohort-completion measure.",size=8.2,color=GRAY)
    footer(c,page);c.showPage();page+=1

    heading(c,"Complete frozen target comparison","September 23 frozen target/weight contract; later adopted targets not applied; all 13 rows")
    body=[]
    for r in targets:
        body.append([r["moment"].replace("_"," "),num(r["target"]),num(r["current_model"]),
                     num(r["proposal_model"]),num(r["proposal_gap"]),num(r["weight"]),num(r["proposal_contribution"])])
    y=table(c,M,H-93,[205,77,77,77,80,90,97],
        ["Moment","Target","Current","Proposal","Proposal gap","Weight","Contribution"],body,row_h=29,sizes={0:8,1:7.5,2:7.5,3:7.5,4:7.5,5:7.5,6:7.5})
    text(c,M,y-22,"The separate normalization row has no weight. Weighted contributions do not choose the policy rate.",size=8.2,color=GRAY)
    text(c,M,y-37,"Full-precision values, current gaps and current contributions are in full_target_comparison.csv.",size=8.2,color=GRAY)
    footer(c,page);c.showPage();page+=1

    heading(c,"Complete parameter and restriction table","Nine selected structural coordinates are held fixed; payroll rate changes and pension balances endogenously")
    body=[]
    for a,b in zip(params_a,params_b):
        assert a["parameter"]==b["parameter"]
        low=a["actual_lower"]
        high=a["actual_upper"]
        body.append([a["parameter"],num(a["estimate"]),num(b["estimate"]),num(low),num(high),
                     "yes" if a["near_actual_bound"]=="True" else "no" if a["near_actual_bound"]=="False" else "-"])
    y=table(c,M,H-91,[234,106,106,88,88,86],
        ["Parameter","Current","Proposal","Lower","Upper","Near bound"],body,row_h=23,sizes={0:8})
    text(c,M,y-22,"Bounds are the selected run's actual search restrictions; blank means externally fixed or derived.",size=8.2,color=GRAY)
    text(c,M,y-37,"The housing supply elasticity is 0.630. The September 24 maturation/entry decision is not implemented.",size=8.2,color=GRAY)
    footer(c,page);c.showPage();page+=1

    names=sorted(graphs_a)
    for start in range(0,len(names),2):
        pair=names[start:start+2]
        heading(c,"Standard diagnostic figures",f"Matched current and proposed cases | figure pairs {start+1}-{start+len(pair)} of 17")
        for index,name in enumerate(pair):
            y_top=H-103-index*225
            text(c,M,y_top,name.replace("_"," ").removesuffix(".png"),size=9,font="Helvetica-Bold")
            left_box=(M,y_top-190,344,177)
            right_box=(M+364,y_top-190,344,177)
            text(c,left_box[0]+4,y_top-16,"Current: 17.9%",size=8,color=BLUE)
            text(c,right_box[0]+4,y_top-16,"Proposal: 8.751%",size=8,color=BLUE)
            draw_fitted_image(c,graphs_a[name],left_box[0],left_box[1],left_box[2],left_box[3]-20)
            draw_fitted_image(c,graphs_b[name],right_box[0],right_box[1],right_box[2],right_box[3]-20)
        footer(c,page);c.showPage();page+=1
    c.save()
    print(PDF)


if __name__=="__main__":
    main()
