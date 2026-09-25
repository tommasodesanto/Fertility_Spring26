from pathlib import Path
import csv, hashlib, json
import pymupdf as fitz
from PIL import Image, ImageDraw

bundle = Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/commute_calibration_20260924_v1")
run = bundle / "results/run_001"
case = run / "worker_03/point_02"
original = run / "selected_export/experimental_commute_calibration.pdf"
outdir = Path("/scratch/td2248/commute_report_review_20260924")
out = outdir / "commute_calibration_review_20260924.pdf"
fits = list(csv.DictReader((case / "target_fit.csv").open()))
params = list(csv.DictReader((case / "parameters.csv").open()))
receipt = json.loads((case / "receipt.json").read_text())
assert len(fits) == 13 and len(params) == 24
assert sum(bool(r["weight"]) for r in fits) == 12
assert receipt["loss"] == 381.0485741441105
assert next(r for r in params if r["parameter"] == "beta_annual")["upper"] == "0.9995"

W, H = 792, 612
LEFT, RIGHT = 38, 754
doc = fitz.open()
ink = (0.10, 0.15, 0.21)
muted = (0.32, 0.37, 0.42)
red = (0.62, 0.12, 0.12)
blue = (0.10, 0.29, 0.48)
line = (0.77, 0.80, 0.83)
pale = (0.94, 0.96, 0.98)

def page(title, subtitle):
    p = doc.new_page(width=W, height=H)
    p.draw_rect(fitz.Rect(0,0,W,9), color=blue, fill=blue)
    p.insert_text((LEFT,39), title, fontname="helv", fontsize=17, color=ink)
    p.insert_text((LEFT,59), subtitle, fontname="helv", fontsize=8.5, color=muted)
    p.draw_line((LEFT,69),(RIGHT,69),color=line,width=0.7)
    return p

def put(p, x, y, s, size=10, color=ink):
    p.insert_text((x,y), str(s), fontname="helv", fontsize=size, color=color)

def paragraph(p, y, heading, lines, color=ink):
    put(p, LEFT, y, heading, 11, blue)
    y += 20
    for s in lines:
        put(p, LEFT+10, y, s, 9.3, color)
        y += 16
    return y+10

def footer(p, n):
    p.draw_line((LEFT,579),(RIGHT,579),color=line,width=0.7)
    put(p,LEFT,596,"Experimental diagnostic | frozen numerical run; corrected report pages",8,muted)
    put(p,RIGHT-38,596,f"{n} / 21",8,muted)

p = page("Experimental commute calibration - reviewed report",
         "Four corrected report pages; 17 selected-case diagnostic pages copied unchanged from frozen export")
y = 95
y = paragraph(p,y,"Numerical result",[
    "Selected weighted loss 381.0485741441105; 12 scored moments, 8 free coordinates,",
    "and a separate completed-fertility normalization to 2.1.",
    "23 objective cases scored; 1 additional case timed out inside a native solve; 17 planned cases unrun.",
    "143 native stationary solves completed; 1 began but did not complete. Export met the one-hour deadline."
])
y = paragraph(p,y,"Target choice pending clarification",[
    "The frozen experiment scored a first-birth room response target of 1.465.",
    "This was our interpretation of the author's phrase 'the new one'; the earlier household target is 1.025.",
    "The author has questioned that interpretation. No target edit, rescore, or numerical rerun is implied."
],red)
y = paragraph(p,y,"Search-bound review failure",[
    "Actual frozen annual discount-factor (beta) search bound: [0.94, 0.9995]. Requested upper cap: 0.99.",
    "This is a preflight/review error. All 23 scored betas were <= 0.9762497404726131,",
    "including the selected beta 0.9762497404726131; the contract itself did not enforce 0.99."
],red)
y = paragraph(p,y,"Measurement and specification limits",[
    "First-birth housing uses an unmatched stationary model proxy for the PSID panel response.",
    "The model bequest row measures gross positive estates, including childless households;",
    "the accepted SCF target is child-directed. The adopted 16/20 birth-entry queue is unimplemented.",
    "PAYGO tax 0.08751017424959717 is experimental, not adopted."
])
put(p,LEFT,y+1,"Plot note: owner demand concentrates at the 10-room maximum in the source figure; cause unproven.",8.5,muted)
put(p,LEFT,y+17,"Income-state labels in the unchanged source figures remain crowded.",8.5,muted)
footer(p,1)

def fmt(v):
    if v is None or v == "":
        return ""
    try:
        return format(float(v), ".12g")
    except (ValueError, TypeError):
        return str(v)

def draw_table(p, headers, widths, rows, y0, row_h, font=8.3, status_col=None):
    x0=LEFT
    total=sum(widths)
    p.draw_rect(fitz.Rect(x0,y0,x0+total,y0+25),color=blue,fill=blue)
    x=x0
    for h,w in zip(headers,widths):
        rc=p.insert_textbox(fitz.Rect(x+5,y0+4,x+w-3,y0+23),h,
                            fontname="helv",fontsize=9,color=(1,1,1))
        assert rc>=0,(h,rc)
        x+=w
    for j,row in enumerate(rows):
        top=y0+25+j*row_h
        if j%2==0: p.draw_rect(fitz.Rect(x0,top,x0+total,top+row_h),color=None,fill=pale)
        p.draw_line((x0,top+row_h),(x0+total,top+row_h),color=line,width=0.4)
        x=x0
        for col,(val,w) in enumerate(zip(row,widths)):
            f=7.7 if col==status_col else font
            rect=fitz.Rect(x+5,top+4,x+w-4,top+row_h-2)
            rc=p.insert_textbox(rect,str(val),fontname="helv",fontsize=f,color=ink)
            assert rc>=0,(j,col,val,rc)
            x+=w
    p.draw_rect(fitz.Rect(x0,y0,x0+total,y0+25+len(rows)*row_h),color=line,width=0.7)

p=page("Complete target fit","Exact frozen 13-row objective; positive weights retained on all 12 scored rows")
frows=[[r["moment"],fmt(r["target"]),fmt(r["model"]),fmt(r["gap"]),
        fmt(r["weight"]),fmt(r["loss_contribution"])] for r in fits]
draw_table(p,["Moment","Target","Model","Gap","Weight","Loss contribution"],
           [164,105,105,105,100,137],frows,88,31,8.1)
put(p,LEFT,550,"Gap = model - target. The fertility-normalization row has no objective weight.",8.5,muted)
footer(p,2)

p=page("Free parameters and normalized inputs",
       "Actual frozen search bounds shown without rounding; near bound uses 1% of the stated interval")
put(p,LEFT,88,"Annual beta upper bound 0.9995 was erroneous; requested cap 0.99. Selected beta is below both.",9,red)
prows=[]
for r in params[:12]:
    prows.append([r["parameter"],fmt(r["estimate"]),fmt(r["lower"]),fmt(r["upper"]),
                  r["near_bound"],r["status"]])
draw_table(p,["Parameter","Estimate","Lower","Upper","Near bound","Status"],
           [178,120,80,80,75,183],prows,103,34,8.4,status_col=5)
footer(p,3)

p=page("Retained settings and annual-to-period inputs",
       "Four-year depreciation is compounded; property tax follows the source's four-year linear convention")
prows=[]
for r in params[12:]:
    prows.append([r["parameter"],fmt(r["estimate"]),fmt(r["lower"]),fmt(r["upper"]),
                  r["near_bound"],r["status"]])
draw_table(p,["Parameter","Estimate","Lower","Upper","Near bound","Status"],
           [178,120,80,80,75,183],prows,91,34,8.4,status_col=5)
footer(p,4)

original_doc=fitz.open(original)
assert len(original_doc)==21
doc.insert_pdf(original_doc,from_page=4,to_page=20)
assert len(doc)==21
doc.save(out,garbage=4,deflate=True)
doc.close()

new=fitz.open(out)
assert len(new)==21
target_text=new[1].get_text()
parameter_text=new[2].get_text()+new[3].get_text()
assert all(r["moment"] in target_text for r in fits)
assert all(r["parameter"] in parameter_text for r in params)
assert all(fmt(r[k]) in target_text for r in fits for k in ("target","model","gap","weight","loss_contribution") if r[k])
assert all(fmt(r[k]) in parameter_text for r in params for k in ("estimate","lower","upper") if r[k])
assert "0.9995" in new[0].get_text() and "1.025" in new[0].get_text()
figures_identical=True
for i in range(17):
    a=original_doc[i+4].get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False).samples
    b=new[i+4].get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False).samples
    if a!=b: figures_identical=False;break
assert figures_identical
qa=outdir/"qa"
qa.mkdir(exist_ok=True)
thumbs=[]
for i in range(4):
    pix=new[i].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(qa/f"page_{i+1:02d}.png")
    im.thumbnail((500,385))
    tile=Image.new("RGB",(520,420),"white")
    tile.paste(im,((520-im.width)//2,5))
    ImageDraw.Draw(tile).text((8,395),f"Corrected page {i+1}",fill="black")
    thumbs.append(tile)
canvas=Image.new("RGB",(1040,840),"#dddddd")
for i,im in enumerate(thumbs): canvas.paste(im,((i%2)*520,(i//2)*420))
canvas.save(qa/"report_contact_sheet.png")
sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
review=dict(status="reviewed_report_only",pdf_path=str(out),pdf_sha256=sha(out),
    original_pdf_sha256=sha(original),figure_pages_pixel_identical=figures_identical,
    pages=21,corrected_report_pages=4,unchanged_selected_figure_pages=17,
    target_rows_checked=len(fits),parameter_rows_checked=len(params),
    selected_loss=receipt["loss"],objective_sha256=receipt["target_system_sha256"],
    selected_checkpoint_sha256=receipt["selected_checkpoint_sha256"],
    reference_checkpoint_sha256=receipt["selected_checkpoint_sha256"],
    final_case_checkpoint_sha256=receipt["case_checkpoint_sha256"],
    source_manifest_sha256=receipt["source_manifest_sha256"],
    beta_frozen_upper=0.9995,beta_requested_upper=0.99,
    beta_selected=receipt["point"]["beta_annual"],
    first_birth_target_frozen=1.465,first_birth_choice_status="pending author clarification",
    numeric_run_unchanged=True)
(outdir/"report_review_receipt.json").write_text(json.dumps(review,indent=2,sort_keys=True)+"\n")
print(json.dumps(review,sort_keys=True))
