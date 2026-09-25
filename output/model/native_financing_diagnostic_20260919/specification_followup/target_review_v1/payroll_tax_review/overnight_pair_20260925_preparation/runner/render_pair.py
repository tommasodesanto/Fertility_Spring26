#!/usr/bin/env python3
"""Report-only selected export for one arm; no household or GE solve."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,os,pickle,sys,tempfile
from pathlib import Path
import pymupdf as fitz
import numpy as np
import run_pair

def scientific_checkpoint(case,chain,lock,old):
    receipt=json.loads((case/"receipt.json").read_text())
    checkpoint=case/"initial_state.pkl.gz"
    if (receipt["case_checkpoint_sha256"]!=old.sha(checkpoint)
        or receipt["target_system_sha256"]!=lock["objective_sha256"]
        or receipt["source_manifest_sha256"]!=lock["source_manifest_sha256"]):
        raise RuntimeError("saved selected/repeat source, target, or checkpoint pin differs")
    with gzip.open(checkpoint,"rb") as f: packet=pickle.load(f)
    if packet["parameters"].adult_entry_clock!="split_birth_vintage":
        raise RuntimeError("saved selected/repeat checkpoint has old adult entry clock")
    solution=packet["solution"]
    return dict(price=receipt["price"],native_price=np.asarray(solution.p_eq),
        native_V=np.asarray(solution.V),native_g=np.asarray(solution.g),
        V=np.asarray(packet["evaluation"].policy.V),
        g=np.asarray(packet["evaluation"].g_current),
        stationary_g_pre=np.asarray(packet["stationary_g_pre"]),
        moments=chain.extract_moments(packet["solution"],packet["parameters"]),
        psi=float(packet["parameters"].psi_child),
        normalization={k:v for k,v in receipt["normalization"].items()
                       if k!="stationary_solve_seconds"},
        loss=receipt["loss"])

def same_numerics(a,b,tolerance):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(
            same_numerics(a[k],b[k],tolerance) for k in a)
    if isinstance(a,(list,tuple,np.ndarray)):
        try:
            left=np.asarray(a,dtype=float);right=np.asarray(b,dtype=float)
            return left.shape==right.shape and bool(np.allclose(left,right,
                               rtol=0,atol=tolerance,equal_nan=True))
        except (TypeError,ValueError):return list(a)==list(b)
    if isinstance(a,(int,float,np.number)) and isinstance(b,(int,float,np.number)):
        return bool(np.isclose(a,b,rtol=0,atol=tolerance,equal_nan=True))
    return a==b

def fmt(value):
    if value is None or value=="": return ""
    try: x=float(value)
    except (TypeError,ValueError): return str(value)
    if x and abs(x)<.0005: return f"{x:.3e}"
    if abs(x)>=1000: return f"{x:,.3f}"
    return f"{x:.3f}"

def table_page(doc,title,subtitle,headers,widths,rows,row_height):
    p=doc.new_page(width=792,height=612)
    p.draw_rect(fitz.Rect(0,0,792,9),color=(.10,.29,.48),fill=(.10,.29,.48))
    p.insert_text((38,39),title,fontsize=17,color=(.1,.15,.21))
    p.insert_text((38,59),subtitle,fontsize=8.5,color=(.3,.35,.4))
    x0=38;y0=86
    p.draw_rect(fitz.Rect(x0,y0,x0+sum(widths),y0+25),color=(.10,.29,.48),fill=(.10,.29,.48))
    x=x0
    for label,w in zip(headers,widths):
        assert p.insert_textbox(fitz.Rect(x+4,y0+4,x+w-2,y0+23),label,
            fontsize=8.5,color=(1,1,1))>=0
        x+=w
    for i,row in enumerate(rows):
        y=y0+25+i*row_height
        if i%2==0:p.draw_rect(fitz.Rect(x0,y,x0+sum(widths),y+row_height),
            color=None,fill=(.94,.96,.98))
        x=x0
        for value,w in zip(row,widths):
            assert p.insert_textbox(fitz.Rect(x+4,y+3,x+w-3,y+row_height-2),
                str(value),fontsize=8,color=(.1,.15,.21))>=0,(title,i,value)
            x+=w
    return p

def render(arm,run_root):
    old=run_pair.load_ancestor()
    lock,objective,bank=run_pair.read_contract(old)
    run_pair.configure(old,arm,lock)
    tax,plan,selected,objective=run_pair.prepare(old,lock)
    arm_root=run_root/arm
    best=json.loads((arm_root/"best_so_far.json").read_text())
    repeat_receipts=[]; repeat_missing=[]
    for i in (1,2):
        complete_path=arm_root/f"repeat_{i:02d}/complete.json"
        if not complete_path.exists():
            repeat_missing.append(i);continue
        complete=json.loads(complete_path.read_text())
        assert complete["status"]=="exact_selected_repeat_complete"
        assert complete["selected_case"]==best["case_path"]
        r=json.loads((arm_root/f"repeat_{i:02d}/selected_repeat/receipt.json").read_text())
        assert r["point"]==best["point"] and r["target_system_sha256"]==lock["objective_sha256"]
        repeat_receipts.append((i,r))
    case=Path(best["case_path"])
    selected_receipt=json.loads((case/"receipt.json").read_text())
    assert selected_receipt["target_system_sha256"]==lock["objective_sha256"]
    assert selected_receipt["source_manifest_sha256"]==lock["source_manifest_sha256"]
    assert selected_receipt["point"]==best["point"]
    with tempfile.TemporaryDirectory(prefix="repeat_signature_",dir=arm_root) as scratch:
        runtime=old.setup_runtime(tax,plan,selected,Path(scratch))
        signature=scientific_checkpoint(case,runtime["chain"],lock,old)
        scientific_verified=[]
        for i,_ in repeat_receipts:
            comparison=scientific_checkpoint(arm_root/f"repeat_{i:02d}/selected_repeat",
                                             runtime["chain"],lock,old)
            if not same_numerics(signature,comparison,float(lock["repeat_model_absolute_tolerance"])):
                raise RuntimeError(f"selected scientific checkpoint differs in repeat {i}")
            scientific_verified.append(i)
        graph_count=len(list((case/"standard_diagnostics").glob("*.png")))
        if graph_count not in (0,17):
            raise RuntimeError("selected case has a partial standard diagnostics gallery")
        if graph_count==0:
            with gzip.open(case/"initial_state.pkl.gz","rb") as f:packet=pickle.load(f)
            runtime["audit"].standard_diagnostics(packet,case,validate_production_young=False)
        if len(list((case/"standard_diagnostics").glob("*.png")))!=17:
            raise RuntimeError("selected saved case lacks 17 standard diagnostic graphs")
    fits=list(csv.DictReader((case/"target_fit.csv").open()))
    params=list(csv.DictReader((case/"parameters.csv").open()))
    assert len(fits)==13 and sum(bool(r["weight"]) for r in fits)==12 and len(params)>=25
    base_rows={r["moment"]:r for r in fits}
    max_model_gap=0.;max_loss_gap=0.;max_parameter_gap=0.
    for i,rep in repeat_receipts:
        max_loss_gap=max(max_loss_gap,abs(float(rep["loss"])-float(best["loss"])))
    # Every economic table cell must reproduce; timing and case-specific paths may differ.
    for i,_ in repeat_receipts:
        rep_rows=list(csv.DictReader((arm_root/f"repeat_{i:02d}/selected_repeat/target_fit.csv").open()))
        assert len(rep_rows)==13
        for row in rep_rows:
            assert row["moment"] in base_rows
            base=base_rows[row["moment"]]
            for field in ("target","model","gap","weight","loss_contribution"):
                if row[field]==base[field]=="":continue
                max_model_gap=max(max_model_gap,abs(float(row[field])-float(base[field])))
        repeated_params=list(csv.DictReader((arm_root/f"repeat_{i:02d}/selected_repeat/parameters.csv").open()))
        assert len(repeated_params)==len(params)
        for base,row in zip(params,repeated_params):
            assert row["parameter"]==base["parameter"] and row["status"]==base["status"]
            for field in ("estimate","lower","upper"):
                if row[field]==base[field]=="":continue
                max_parameter_gap=max(max_parameter_gap,abs(float(row[field])-float(base[field])))
    if max_model_gap>float(lock["repeat_model_absolute_tolerance"]) or max_loss_gap>float(lock["repeat_loss_absolute_tolerance"]):
        raise RuntimeError("selected exact repeats exceed reviewed tolerances")
    if max_parameter_gap>float(lock["repeat_model_absolute_tolerance"]):
        raise RuntimeError("selected repeat parameter table differs")
    # Reuse the saved-case plot exporter; it never solves the model.
    old.write(arm_root/"deadline.json",json.loads((run_root/"deadline.json").read_text()))
    old.export_selected(arm_root)
    original=arm_root/"selected_export/experimental_commute_calibration.pdf"
    source_pdf=fitz.open(original)
    assert len(source_pdf)==21
    doc=fitz.open()
    p=doc.new_page(width=792,height=612)
    p.draw_rect(fitz.Rect(0,0,792,9),color=(.10,.29,.48),fill=(.10,.29,.48))
    p.insert_text((38,42),"Paired overnight calibration - "+arm,fontsize=18)
    lines=[
        f"PAYGO payroll tax: {run_pair.RATES[arm]:.12g}; pension balances endogenously.",
        f"Selected weighted loss: {best['loss']:.6g}; 12 scored moments, 8 free coordinates, psi normalized to 2.1.",
        "Both tax arms use identical targets, bounds, source, common seeds, and proposal bank.",
        "This is a recalibrated-arm comparison, not a fixed-parameter causal tax effect.",
        "First-birth rooms target: exact 1.465; model observer is an unmatched stationary proxy.",
        "SCF bequest target is child-directed; model observer sums gross positive estates.",
        "Adult entry: split 16/20 birth-vintage queue, once divided by 2.1; closed renewal gate passed.",
        f"Exact repeats verified: {len(repeat_receipts)}/2; missing or failed: {repeat_missing}.",
        f"Available-repeat max table-cell gap {max_model_gap:.3e}, loss gap {max_loss_gap:.3e}.",
        f"Objective SHA256: {lock['objective_sha256']}",
        f"Source inventory SHA256: {lock['source_manifest_sha256']}",
        "17 standard selected-case figure pages follow the complete tables."
    ]
    y=82
    for line in lines:
        assert p.insert_textbox(fitz.Rect(44,y,748,y+30),line,fontsize=10)>=0,line
        y+=35
    table_page(doc,"Complete target fit","Exact source values remain in target_fit.csv",
        ["Moment","Target","Model","Gap","Weight","Loss"],
        [164,105,105,105,100,137],
        [[r["moment"],fmt(r["target"]),fmt(r["model"]),fmt(r["gap"]),
          fmt(r["weight"]),fmt(r["loss_contribution"])] for r in fits],31)
    for index,start in enumerate((0,13),1):
        table_page(doc,f"Parameters and restrictions ({index}/2)","Exact bounds remain in parameters.csv",
            ["Parameter","Estimate","Lower","Upper","Near","Status"],
            [178,120,80,80,75,183],
            [[r["parameter"],fmt(r["estimate"]),fmt(r["lower"]),fmt(r["upper"]),
              r["near_bound"],r["status"]] for r in params[start:start+13]],32)
    doc.insert_pdf(source_pdf,from_page=4,to_page=20)
    assert len(doc)==21
    output=arm_root/"selected_export"/f"paired_{arm}_review.pdf"
    doc.save(output,garbage=4,deflate=True)
    internal=arm_root/"selected_export/INTERNAL_NOT_FOR_DELIVERY"
    internal.mkdir(exist_ok=True)
    original.rename(internal/original.name)
    old_receipt=arm_root/"selected_export/final_receipt.json"
    if old_receipt.exists():old_receipt.rename(internal/old_receipt.name)
    reviewed=fitz.open(output)
    assert len(reviewed)==21
    assert all(r["moment"] in reviewed[1].get_text() for r in fits)
    assert all(r["parameter"] in reviewed[2].get_text()+reviewed[3].get_text() for r in params)
    for i in range(17):
        assert source_pdf[i+4].get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False).samples == reviewed[i+4].get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False).samples
    ledgers=list(arm_root.rglob("stationary_solves.json"))
    attempted_cases={p.parent for p in ledgers}
    attempted_cases.update(p.parent for p in arm_root.rglob("receipt.json")
                           if "selected_export" not in p.parts)
    attempted_cases.update(p.parent for p in arm_root.rglob("failure.json"))
    solve_events=[event for path in ledgers for event in json.loads(path.read_text())]
    native_started=sum(event.get("started_epoch") is not None for event in solve_events)
    native_completed=sum(event.get("status")=="completed" for event in solve_events)
    objective_scored=sum((case/"receipt.json").exists() for case in attempted_cases)
    objective_attempted=len(attempted_cases)
    final=dict(status="paired_selected_export_complete",arm=arm,loss=best["loss"],
        selected_case=best["case_path"],objective_sha256=lock["objective_sha256"],
        source_manifest_sha256=lock["source_manifest_sha256"],pdf_path=str(output),
        pdf_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),target_rows=13,
        weighted_rows=12,parameter_rows=len(params),figure_pages=17,
        repeats_verified=len(repeat_receipts),repeats_missing_or_failed=repeat_missing,
        exact_repeat_claim=len(scientific_verified)==2,
        scientific_checkpoint_repeats_verified=scientific_verified,
        repeat_model_max_absolute_gap=max_model_gap,repeat_parameter_max_absolute_gap=max_parameter_gap,
        repeat_loss_max_absolute_gap=max_loss_gap,
        repeat_absolute_tolerance=lock["repeat_model_absolute_tolerance"],
        objective_planned=364,objective_attempted=objective_attempted,
        objective_scored=objective_scored,objective_incomplete=objective_attempted-objective_scored,
        objective_unrun=364-objective_attempted,native_started=native_started,
        native_completed=native_completed,native_incomplete=native_started-native_completed)
    old.write(arm_root/"selected_export/final_pair_receipt.json",final)
    print(json.dumps(final,sort_keys=True))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--arm",choices=tuple(run_pair.RATES),required=True)
    p.add_argument("--run-root",type=Path,required=True)
    a=p.parse_args()
    render(a.arm,a.run_root.resolve())
if __name__=="__main__":main()
