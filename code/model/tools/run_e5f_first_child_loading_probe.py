#!/usr/bin/env python3
"""Torch-only four-point diagnostic: first-child housing share, no housing floor.

Reuse the authenticated comparison runtime and all scientific gates. Curvature
is fixed, one-child flow benefit preserved, no calibration or normalization.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import threading
import time

import run_e5f_soft_housing_probe as common

CASES = (('floor_concave', None), ('shares_zero', 0.), ('shares_010', .1), ('shares_020', .2))
LABELS = {'control':'Floor, linear benefit','floor_concave':'Floor, concave benefit',
          'shares_zero':'No floor, no loading','shares_010':'No floor, +10 points','shares_020':'No floor, +20 points'}


def configure(P, name, loading):
    out = copy.deepcopy(P)
    if name == 'control': return out
    out.utility_child_benefit_exponent = .86
    out.utility_comparison_arm = 'floor_concave' if loading is None else 'shares_concave'
    if loading is not None:
        out.child_room_floor = False
        out.hbar_first_child_jump = out.hbar_child_rooms = 0.
        out.delta_alpha_jump, out.delta_alpha = float(loading), 0.
    return out


def parameters(reference, P, name, loading):
    rows = list(csv.DictReader((reference/'parameters.csv').open()))
    for row in rows:
        if row['parameter'] == 'h_P' and loading is not None:
            row.update(estimate=0., lower='', upper='', near_bound='',status='housing floor removed in this experiment')
        if row['parameter'] == 'child_benefit_exponent':
            row.update(estimate=P.utility_child_benefit_exponent,status='fixed diagnostic restriction; not estimated')
    for key, value, status in (
        ('child_benefit_curvature', 1.-P.utility_child_benefit_exponent, 'fixed trial value; distinct from fertility taste-shock scales'),
        ('child_benefit_CRRA_coefficient', P.psi_child*P.utility_child_benefit_exponent, 'psi in psi*m^(1-kappa)/(1-kappa); preserves one-child benefit'),
        ('delta_alpha_jump', P.delta_alpha_jump, 'fixed first-child housing share increment; no later-child increment'),
        ('delta_alpha', P.delta_alpha, 'fixed at zero; no additional housing share loading for subsequent children')):
        rows.append(dict(parameter=key,estimate=value,lower='',upper='',near_bound='',status=status))
    return rows


def smoke(model, P, grid):
    import numpy as np
    import e5f_utility_comparison_runtime as adapter
    checks = []
    for name, loading in (('control',None),)+CASES:
        Q = configure(P,name,loading)
        sd = model.precompute_shared(Q,grid)
        native = model.precompute_shared._four_arm_native_precompute(Q,grid)
        exponent = 1. if name == 'control' else .86
        for n in range(Q.n_parity):
            for m in range(n+1):
                ix = n+Q.n_parity*m
                assert sd.psi_flat.reshape(-1)[ix] == Q.psi_child*m**exponent
                alpha = Q.alpha_cons-(float(loading or 0.) if m>0 else 0.)
                assert sd.alpha_flat.reshape(-1)[ix] == alpha
                floor = P.hbar_first_child_jump if loading is None and m>0 else 0.
                assert sd.hb_flat.reshape(-1)[ix] == floor
                if loading is not None:
                    factor = float(adapter.reference_composite_factor(alpha,Q.utility_reference_rent))
                    k = alpha**alpha*((1-alpha)/Q.utility_reference_rent)**(1-alpha)
                    a = Q.alpha_cons
                    k0 = a**a*((1-a)/Q.utility_reference_rent)**(1-a)
                    np.testing.assert_allclose(factor*k,k0,rtol=1e-14,atol=0)
                    np.testing.assert_allclose(sd.escale_flat.reshape(-1)[ix],native.escale_flat.reshape(-1)[ix]*factor**(1-Q.sigma),rtol=1e-14)
        if name == 'control':
            for key in ('alpha_flat','hb_flat','cb_flat','escale_flat','psi_flat'):
                np.testing.assert_array_equal(getattr(sd,key),getattr(native,key))
        checks.append(dict(case=name,loading=loading,exponent=exponent,one_child_benefit=Q.psi_child))
    return dict(status='passed',scope='Exact case configuration loop, native arrays, common-price normalization and unchanged control',cases=checks)


def render(output, records, failures):
    """Readable summary and full tables followed by the unchanged graph sets."""
    from reportlab.pdfgen import canvas
    from reportlab.lib.utils import ImageReader
    from reportlab.lib import colors
    from reportlab.platypus import Table, TableStyle, Paragraph
    from reportlab.lib.styles import getSampleStyleSheet
    import matplotlib.pyplot as plt
    import numpy as np
    out = output/'first_child_loading.pdf'
    cv = canvas.Canvas(str(out),pagesize=(842,595))
    page = 0
    styles = getSampleStyleSheet()
    small = styles['BodyText']; small.fontSize=8; small.leading=10
    def heading(title,subtitle=''):
        nonlocal page
        page += 1
        cv.setFont('Helvetica-Bold',17); cv.drawString(32,561,title)
        cv.setFont('Helvetica',9); cv.drawString(32,543,subtitle)
        cv.setFont('Helvetica',8); cv.drawRightString(810,18,str(page))
    def draw_table(data, widths, top=515, size=8):
        tbl=Table(data,colWidths=widths,repeatRows=1)
        tbl.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),size),('LEADING',(0,0),(-1,-1),size+2),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E9EEF3')),('BOTTOMPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.5,colors.grey)]))
        _,height=tbl.wrap(778,500); tbl.drawOn(cv,32,top-height)
    def fmt(x):
        if x=='':return ''
        if isinstance(x,str):
            try:x=float(x)
            except ValueError:return x
        if isinstance(x,bool):return str(x)
        return f'{x:.3g}' if abs(x)>1000 or (0<abs(x)<.001) else f'{x:.3f}'
    heading('First-child housing loading without Stone-Geary','Bounded fixed-preference diagnostic; no recalibration or production adoption')
    data=[['Case','Parent housing weight','Fertility','Birth rooms','Childless (%)','Loss']]
    for r in records:
        moments={x['moment']:x['model'] for x in r['target_fit']}
        data.append([LABELS[r['case']],fmt(.267+(r['loading'] or 0.)),fmt(moments['initial_normalization']),fmt(moments['first_birth_rooms']),fmt(100*moments['cps_childlessness']),fmt(r['loss'])])
    draw_table(data,[165,115,90,110,100,100])
    notes=[
        'One-child benefit held fixed. New cases fix child-benefit curvature at 0.140; original control is linear.',
        'Positive loading raises housing weight only when children are at home. Later children add no further share loading.',
        'Earnings, initial wealth/income, mortality, bequests, credit, shocks, supply and fiscal rule are unchanged.',
        'Changing housing weights uses the inherited fixed-reference-rent utility normalization; no estimated scale added.',
        'Entry composition is held fixed. Reproduction gaps are reported rather than closed by refitting child benefits.',
        'Birth rooms is the inherited stationary proxy, not a matched event study. Frozen old ACS targets retained for comparison.',
        'The original floor response also misses the 1.465-room target. These fixed points cannot rank optimized specifications.']
    y=330
    for line in notes:
        cv.setFont('Helvetica',9); cv.drawString(32,y,line); y-=19
    for failure in failures:
        cv.drawString(32,y,('FAILED '+failure['case']+': '+failure['error'])[:150]); y-=18
    cv.showPage()
    heading('Specification and interpretation','All symbols below describe economic objects; numerical choices are separate')
    formulas=[r'$u(C,m)=\frac{C^{1-\sigma}}{1-\sigma}+\psi\frac{m^{1-\kappa}}{1-\kappa}$',
              r'$C=A(m)\,\frac{c^{\alpha(m)}s^{1-\alpha(m)}}{e(m)}$',
              r'$\alpha(m)=\alpha_0-\lambda\,\mathbf{1}\{m>0\},\qquad e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}$']
    fig,ax=plt.subplots(figsize=(10,2.2));ax.axis('off')
    for i,line in enumerate(formulas):ax.text(.02,.85-i*.34,line,fontsize=17)
    fig.savefig(output/'equations.png',dpi=170,bbox_inches='tight');plt.close(fig)
    cv.drawImage(ImageReader(str(output/'equations.png')),32,325,width=730,height=190,preserveAspectRatio=True,anchor='c')
    lines=[
        'c: nonhousing consumption; s: housing services, including the unchanged owner premium; m: children at home.',
        'sigma: material-consumption curvature (2 here); kappa: child-benefit curvature (0.140 in new cases).',
        'psi weights enjoyment of children relative to material consumption; C is the household material-consumption index.',
        'alpha_0 = 0.733 is the childless consumption weight; lambda is the positive first-child housing loading.',
        'e(m) is the household needs scale: 1.000 without children, about 1.234 with one and 1.450 with two.',
        'A(m) is a substantive fixed normalization, not an estimated parameter. For uncapped renters at the reference rent,',
        'equal expenditure attains the same material index before household needs. This does not extend to all owners.',
        'The benefit code uses b*m^(1-kappa). Displayed psi=(1-kappa)*b, so changing curvature preserves b at m=1.',
        'Housing shares increase from 0.267 without children to 0.367 or 0.467 with children in the positive-loading cases.']
    for i,line in enumerate(lines):cv.setFont('Helvetica',9);cv.drawString(32,298-i*21,line)
    cv.showPage()
    heading('Fertility incentives and demographic closure','Positive gap means trying now is preferred to waiting before the current taste draw')
    data=[['Case','First births at\npositive gaps (%)','Later births at\npositive gaps (%)','Zero first-birth\nattempts (% risk)','Replacement\ngap (%)']]
    for r in records:
        first,later=r['incentives'][:2]
        data.append([LABELS[r['case']],fmt(100*first['birth_share_positive_interior']),fmt(100*later['birth_share_positive_interior']),fmt(100*first['risk_share_zero_try']),fmt(100*r['birth_replacement_relative_gap'])])
    draw_table(data,[185,150,150,150,130])
    cv.setFont('Helvetica',10)
    cv.drawString(32,320,'Continuation values retain future taste uncertainty. This is not a zero-shock equilibrium or a lifetime childlessness comparison.')
    cv.drawString(32,300,'Zero probabilities combine infeasibility and underflow. Changes in weighted shares also reflect population composition.')
    cv.drawString(32,280,'Replacement gaps are reported at fixed child benefit and fixed normalized entry; no closed-renewal condition is imposed.')
    cv.showPage()
    for r in records:
        folder=output/r['case']
        heading(LABELS[r['case']]+': complete target fit','Same target definitions and weights across cases; normalization row is descriptive at fixed child benefit')
        data=[['Moment','Target','Model','Gap','Weight','Loss contribution']]
        data += [[row['moment']]+[fmt(row[k]) for k in ('target','model','gap','weight','loss_contribution')] for row in r['target_fit']]
        draw_table(data,[208,96,96,96,120,145]);cv.showPage()
        params=list(csv.DictReader((folder/'parameters.csv').open()))
        for start in range(0,len(params),16):
            heading(LABELS[r['case']]+': parameters and restrictions','All values fixed in this diagnostic; inherited estimation bounds are shown where applicable')
            data=[['Parameter','Value','Lower','Upper','Near bound','Status']]
            for row in params[start:start+16]:
                data.append([Paragraph(row['parameter'].replace('_',' '),small),fmt(row['estimate']),fmt(row['lower']),fmt(row['upper']),row['near_bound'],Paragraph(row['status'],small)])
            draw_table(data,[175,72,55,55,65,356],size=8);cv.showPage()
    # Preserve every standard diagnostic with identical names/layouts, paired by case.
    names=sorted(p.name for p in (output/records[0]['case']/'standard_diagnostics').glob('*.png'))
    for name in names:
        for start in range(0,len(records),2):
            heading(name.replace('_',' ').replace('.png',''),'Standard diagnostic; same graph set for every completed case')
            for j,r in enumerate(records[start:start+2]):
                cv.setFont('Helvetica-Bold',10);cv.drawString(32+395*j,517,LABELS[r['case']])
                cv.drawImage(ImageReader(str(output/r['case']/'standard_diagnostics'/name)),32+395*j,75,width=385,height=420,preserveAspectRatio=True,anchor='c')
            cv.showPage()
    cv.save()
    common.write(output/'pdf_receipt.json',dict(pages=page,case_count=len(records),standard_graphs_per_case=len(names),path=str(out)))
    import fitz
    from PIL import Image, ImageDraw
    qa=output/'pdf_qa';qa.mkdir(exist_ok=True)
    doc=fitz.open(out);outside=[]
    for start in range(0,len(doc),12):
        sheet=Image.new('RGB',(1140,1160),'#dddddd');draw=ImageDraw.Draw(sheet)
        for k in range(start,min(start+12,len(doc))):
            page=doc[k]
            for block in page.get_text('blocks'):
                if block[0]<-1 or block[1]<-1 or block[2]>843 or block[3]>596:outside.append((k+1,block[:4]))
            pix=page.get_pixmap(matrix=fitz.Matrix(.44,.44),alpha=False)
            image=Image.frombytes('RGB',(pix.width,pix.height),pix.samples)
            x=((k-start)%3)*380;y=((k-start)//3)*290
            sheet.paste(image,(x,y+20));draw.text((x+6,y+4),str(k+1),fill='black')
            if k<5:page.get_pixmap(matrix=fitz.Matrix(1.2,1.2),alpha=False).save(str(qa/f'page_{k+1}.png'))
        sheet.save(qa/f'contact_{start//12+1}.png')
    assert not outside,outside
    common.write(qa/'receipt.json',dict(pages=len(doc),text_outside_page=outside,contact_sheets=(len(doc)+11)//12))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--render-only',action='store_true');args=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Torch allocation required')
    if args.render_only:
        completed=json.loads((args.output/'complete.json').read_text());render(args.output,completed['records'],completed['failures']);return
    import numpy as np
    import collect_e5f_utility_comparison as collector
    args.output.mkdir(parents=True,exist_ok=False)
    os.environ[collector.ENV_PIN]=common.PIN
    contract,runner=collector._load_contract(common.BASE/'launch_v1/contract.json')
    _,_,_,tax,_,objective,runtime,_=runner.setup(contract,'floor_linear',args.output/'runtime')
    choice=collector.read_json(common.BASE/'results/run_001/floor_linear/selected.json')
    selected=collector.scientific_checkpoint(choice['original_case_output'],contract,'floor_linear',runtime,tax)
    packet=selected['packet'];P=packet['parameters'];grid=packet['b_grid'];model=runtime['model']
    started=time.monotonic();state={'phase':'smoke','elapsed_seconds':0.};done=threading.Event()
    def heartbeat():
        while not done.wait(30):
            state['elapsed_seconds']=time.monotonic()-started;common.write(args.output/'heartbeat.json',state)
    threading.Thread(target=heartbeat,daemon=True).start()
    common.write(args.output/'smoke.json',smoke(model,P,grid))
    common.write(args.output/'design.json',dict(cases=CASES,source_pins=selected['pins'],target_contract=common.PIN,
        restrictions='Four diagnostic fixed points: floor/concave, no floor at first-child share loadings 0,.1,.2. No later-child share loading. Curvature .14 fixed, first-child benefit fixed. No refit.',
        normalization='Inherited compensated-expenditure normalization at reference rent; substantive trial restriction, no extra estimated parameter.',
        unchanged='Earnings, entry wealth/income, timing, mortality, bequests, credit, shocks, housing products, supply, fiscal rule, targets, weights and numerical gates.',
        budget='Four new GE solves; 25 minute allocation including charts. Reference stationary solve about135 seconds; roughly9 minutes plus reporting. No retries or search.',
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(common.__file__))}))
    records=[];failures=[]
    for name,loading in (('control',None),)+CASES:
        state['phase']=name;common.write(args.output/'heartbeat.json',state)
        case=args.output/name;case.mkdir();Q=configure(P,name,loading)
        try:
            if name=='control':current=packet
            else:
                sol,Q,prices,_=runtime['solve_balanced_initial_equilibrium'](model=model,parameters=Q,b_grid=grid,initial_prices=np.asarray(packet['solution'].p_eq).copy(),payroll_tax=float(P.tau_pay),marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
                sd=model.precompute_shared(Q,grid);Q._fert2_probs=sol.fert2_probs.copy();cal=runtime['primitive'].pf.calendar
                policy=cal.policy_from_solution(sol,prices,Q,grid,sd)
                pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,Q,grid,sd)
                assert reconstruction['stationary_post_fertility_nesting_l1']<5e-9
                supply=cal.HousingSupplyRule('static-elastic',float(prices[0]),float(Q.H0[0]*(Q.user_cost_rate*prices[0]/Q.r_bar[0])**Q.xi_supply[0]),float(Q.xi_supply[0]))
                evaluation=cal.evaluate_period(prices,pre,Q,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=policy)
                current=dict(parameters=Q,b_grid=grid,solution=sol,shared=sd,stationary_g_pre=pre,evaluation=evaluation,supply_rule=supply,ancestry_contract_sha256=common.PIN,demographic_seed=packet.get('demographic_seed'))
                with gzip.open(case/'initial_state.pkl.gz','wb',compresslevel=1) as stream:pickle.dump(current,stream,protocol=5)
            record=common.report(current,runtime,tax,objective,selected['case'],case,name,parameter_rows=parameters(selected['case'],Q,name,loading),old_housing_floor=P.hbar_first_child_jump)
            record.update(loading=loading,child_curvature=1-Q.utility_child_benefit_exponent)
            common.write(case/'receipt.json',record);records.append(record)
            common.write(args.output/'latest_completed.json',record)
            common.write(args.output/'best_so_far.json',dict(interpretation='Descriptive fixed-point loss, not calibrated winner',case=min(records,key=lambda r:r['loss'])['case']))
            print(json.dumps(dict(case=name,status='completed',elapsed_seconds=time.monotonic()-started,loss=record['loss'])),flush=True)
        except Exception as exc:
            failure=dict(case=name,status='case_failed_requires_classification',error=str(exc),type=type(exc).__name__)
            common.write(case/'failure.json',failure);failures.append(failure)
            print(json.dumps(failure),flush=True)
            # Continue only to the other predeclared points; no altered gates or fallback.
    common.write(args.output/'complete.json',dict(status='completed_fixed_design' if not failures else 'completed_with_failures',records=records,failures=failures,elapsed_seconds=time.monotonic()-started,job=os.environ['SLURM_JOB_ID']))
    state['phase']='rendering';render(args.output,records,failures);done.set()
    print(json.dumps(dict(status='finished',elapsed_seconds=time.monotonic()-started)),flush=True)


if __name__=='__main__':main()
