#!/usr/bin/env python3
"""Bounded Torch-only PE test of housing affordability and fertility preferences.

No production equilibrium, preference adoption, or broad calibration. Frozen
native runtime, fixed reference normalization, common prices and fiscal inputs.
"""
import argparse, copy, csv, gzip, hashlib, json, os, pickle, threading, time
from pathlib import Path
import numpy as np
import run_e5f_soft_housing_probe as common
import run_e5f_first_child_loading_probe as loading

ROOT=common.BASE.parent
PROBE=ROOT/'first_child_loading_probe_20260926/run_001'
ARMS=(('floor_concave',None),('shares_010',.1),('shares_020',.2))
LABEL={'floor_concave':'Floor','shares_010':'Share +10 points','shares_020':'Share +20 points'}
B_BRACKET=(.001,.400)
TARGET=2.1
TOL=.002
MAX_B_EVAL=12


def save(path,obj):
    with gzip.open(path,'wb',compresslevel=1) as f:pickle.dump(obj,f,protocol=5)


def setup(out):
    import collect_e5f_utility_comparison as collector
    os.environ[collector.ENV_PIN]=common.PIN
    contract,runner=collector._load_contract(common.BASE/'launch_v1/contract.json')
    _,_,_,tax,_,objective,rt,_=runner.setup(contract,'floor_linear',out/'runtime')
    choice=collector.read_json(common.BASE/'results/run_001/floor_linear/selected.json')
    selected=collector.scientific_checkpoint(choice['original_case_output'],contract,'floor_linear',rt,tax)
    return rt,tax,objective,selected


def loop_smoke(P,grid,model):
    result=loading.smoke(model,P,grid)
    cases=[]
    for arm,lam in ARMS:
        for factor in (1.,1.1):
            for b in (P.psi_child,)+B_BRACKET:
                Q=loading.configure(P,arm,lam);Q.psi_child=b
                sd=model.precompute_shared(Q,grid)
                assert Q.utility_reference_rent==P.utility_reference_rent
                assert Q.pension==P.pension and Q.tau_pay==P.tau_pay
                assert np.isfinite(sd.psi_flat).all()
                cases.append([arm,factor,b])
    result.update(exact_loop_configurations=cases,bracket=B_BRACKET,tolerance=TOL,max_b_evaluations=MAX_B_EVAL)
    return result


def solve(rt,P,grid,price):
    model=rt['model'];cal=rt['primitive'].pf.calendar
    sd=model.precompute_shared(P,grid)
    sol=model.solve_markov_income_at_prices(price,P,grid,SD=sd)
    P._fert2_probs=sol.fert2_probs.copy()
    policy=cal.policy_from_solution(sol,price,P,grid,sd)
    pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,sd)
    ev=cal.evaluate_period(price,pre,P,grid,sd,cal.SolveCounter(),supplied_policy=policy)
    assert reconstruction['stationary_post_fertility_nesting_l1']<5e-9
    return dict(parameters=P,b_grid=grid,solution=sol,shared=sd,stationary_g_pre=pre,evaluation=ev)


def validate(packet,rt,out):
    P=packet['parameters'];sol=packet['solution'];ev=packet['evaluation'];grid=packet['b_grid'];sd=packet['shared']
    gates=rt['primitive'].pf.transition.operator_gates(sol,ev.policy,packet['stationary_g_pre'],P,grid,sd)
    for key in ('one_step_constant_path_nesting_l1','mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
        assert abs(gates[key])<5e-9,(key,gates[key])
    assert abs(gates['zero_entry_mass_accounting_residual'])<2e-8
    assert ev.feasibility_projection_mass<1e-6
    budget=rt['primitive'].dated_budget(ev,P,sd,grid,float(P.user_cost_rate*sol.p_eq[0]))
    assert budget['budget_excess_mass']<=2e-10
    arrays=rt['audit'].policy_array_audit(packet,out)
    assert not arrays['occupied_negative_steps']
    assert not any(v['nonfinite'] or v['minimum']<0 or v['maximum']>1 for v in arrays['probabilities'].values())
    return dict(operator=gates,budget=budget,feasibility_projection_mass=float(ev.feasibility_projection_mass),market_residual=float(ev.relative_market_residual),market_clearing_required=False)


def common_states(packet,control,rt,out):
    """Same initial state weights, optimizing fertility then location/tenure."""
    P=packet['parameters'];grid=packet['b_grid'];policy=packet['evaluation'].policy
    cal=rt['primitive'].pf.calendar;sd=packet['shared'];g=control['stationary_g_pre']
    rows=[]
    for tenure,label in ((0,'current_renters'),(1,'current_owners')):
        for low in (False,True):
            weights=np.zeros_like(g)
            for j in range(P.J):
                age=P.age_start+j*P.da
                if age>34 or not P.A_f_start<=j+1<=P.A_f_end:continue
                if tenure==0:weights[:,0,:,j,:,0,0]=g[:,0,:,j,:,0,0]
                else:weights[:,1:,:,j,:,0,0]=g[:,1:,:,j,:,0,0]
            if low:weights[grid>1.]=0
            ev=cal.evaluate_period(policy.price,weights,P,grid,sd,cal.SolveCounter(),supplied_policy=policy)
            mass=float(weights.sum());current=ev.g_current
            assert ev.feasibility_projection_mass<1e-6
            own=float(current[:,1:].sum());h=float(cal.housing_demand_by_location(current,policy.hR_pol,P).sum())
            c=float((current*policy.c_pol).sum())
            rows.append(dict(group=label+('_wealth_le_1' if low else ''),mass=mass,birth_probability=float(ev.births)/mass,realized_ownership=own/mass,realized_housing=h/mass,realized_nonhousing=c/mass,projected_mass=float(ev.feasibility_projection_mass)))
    common.table(out/'common_states.csv',rows)
    # Common-grid conditional policies and gaps; no occupancy conditioning.
    rows=[];n_z=int(g.shape[4]);zs=sorted(set([0,n_z//2,n_z-1]))
    for j in range(P.J):
        age=P.age_start+j*P.da
        if age not in (22,26,30):continue
        for z in zs:
            for n,m in ((0,0),(1,1)):
                pair=policy.fert_probs[:,0,0,j,z,:2] if n==0 else policy.fert2_probs[:,0,0,j,z,:,n-1,m]
                for ib,b in enumerate(grid):
                    p0,p1=map(float,pair[ib]);interior=p0>0 and p1>0
                    rows.append(dict(age=age,income_index=z,birth_number=n+1,wealth=float(b),common_mass=float(g[ib,0,0,j,z,n,m]),attempt_probability=p1,gap_over_shock_scale=float(np.log(p1)-np.log(p0)) if interior else '',choice_status='finite' if interior else ('unavailable' if p0==p1==0 else 'endpoint'),conditional_renter_housing=float(policy.hR_pol[ib,0,0,j,z,n,m]),conditional_renter_consumption=float(policy.c_pol[ib,0,0,j,z,n,m])))
    common.table(out/'common_grid.csv',rows)
    return rows


def report(packet,rt,tax,objective,reference,out,arm,lam,phase,factor,control):
    P=packet['parameters'];ev=packet['evaluation'];sol=packet['solution'];grid=packet['b_grid'];sd=packet['shared']
    gates=validate(packet,rt,out)
    early=dict(fertility={p:rt['observe_initial_fertility'](ev,P,age_projection=p) for p in ('uniform_birth_time','constant_post_cell')},housing_wealth=rt['observe_initial_housing_wealth'](ev,P,grid,sd,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True))
    recent=rt['observe_recent_parent_flow'](ev,P,diagnostic_enabled=True,snapshot=rt['SNAPSHOT'],age_projection=rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance={'case_id':out.name})
    tfr=float(rt['chain'].extract_moments(sol,P)['tfr'])
    fits=tax.target_rows(objective,early,recent['model_value'],tfr);common.table(out/'target_fit.csv',fits)
    params=loading.parameters(reference,P,arm,lam)
    actual=tax.actual_parameters(P)
    for row in params:
        if row['parameter'] in actual:row['estimate']=actual[row['parameter']]
        if row['parameter']=='pension_period':row['estimate']=P.pension
        if row['parameter']=='psi_child':
            row['estimate']=P.psi_child
            row['status']='Only adjusted coordinate in matched baseline; frozen under price shock' if phase=='matched' else 'Fixed inherited child benefit'
        elif row['status']=='experimental free coordinate':row['status']='Inherited search coordinate; fixed in this diagnostic'
        if row['parameter']=='psi_child' and phase=='matched':
            row.update(lower=B_BRACKET[0],upper=B_BRACKET[1],near_bound=min(P.psi_child-B_BRACKET[0],B_BRACKET[1]-P.psi_child)<=.02*(B_BRACKET[1]-B_BRACKET[0]))
        if row['parameter']=='child_benefit_CRRA_coefficient':row['estimate']=P.psi_child*.86
    common.table(out/'parameters.csv',params)
    inc=common.incentives(packet,out);common_states(packet,control,rt,out)
    # All 17 standard panels retained for stationary PE solutions, labeled PE.
    rt['audit'].standard_diagnostics(packet,out,validate_production_young=False)
    assert len(list((out/'standard_diagnostics').glob('*.png')))==17
    result=dict(case=out.name,arm=arm,phase=phase,price_factor=factor,price=float(sol.p_eq[0]),rent=float(P.user_cost_rate*sol.p_eq[0]),b=float(P.psi_child),fertility=tfr,target_fit=fits,incentives=inc,gates=gates,replacement_gap=float(sol.adult_entry_stationary_relative_gap),closure='Permanent fixed-price partial equilibrium; fixed fiscal inputs and normalized entry; market, fiscal and demographic imbalances not cleared')
    common.write(out/'receipt.json',result);save(out/'state.pkl.gz',packet)
    return result


def render(out):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from reportlab.pdfgen import canvas
    from reportlab.platypus import Table,TableStyle,Paragraph
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib import colors
    from reportlab.lib.utils import ImageReader
    complete=json.loads((out/'complete.json').read_text());records=complete['records']
    # Report-only wording: fiscal primitives are frozen in these PE cases.
    for record in records:
        path=out/record['case']/'parameters.csv';params=list(csv.DictReader(path.open()))
        for row in params:
            if row['parameter']=='pension_period':row['status']='Held at authenticated baseline pension; no rebalancing'
            if row['parameter']=='payroll_tax':row['status']='Held at authenticated baseline tax rate'
            if row['parameter']=='child_benefit_CRRA_coefficient' and record['phase']=='matched':row['status']='psi=(1-kappa)*b; matched b held fixed under price shock'
        common.table(path,params)
    figures=[]
    data=list(csv.DictReader((out/records[0]['case']/'common_grid.csv').open()))
    zs=(json.loads((out/'plot_income_states.json').read_text())['indices'] if (out/'plot_income_states.json').exists() else sorted({int(v['income_index']) for v in data}))
    income_labels=dict(zip(zs,('Low income (p10)','Middle income (p50)','High income (p90)')))
    fig,axs=plt.subplots(1,3,figsize=(12,3.7),constrained_layout=True)
    for ax,z in zip(axs,zs):
        for age in (22,26,30):
            vals=[v for v in data if int(v['birth_number'])==1 and float(v['age'])==age and int(v['income_index'])==z and float(v['wealth'])<=5]
            ax.plot([float(v['wealth']) for v in vals],[float(v['common_mass']) for v in vals],label=f'Age {age}')
        ax.set_title(income_labels[z]);ax.set_xlabel('Beginning liquid wealth');ax.set_ylabel('Fixed control pre-choice mass')
    axs[0].legend();fig.suptitle('Common occupied weights: current renters without prior children')
    path=out/'common_occupied_weights.png';fig.savefig(path,dpi=140);plt.close(fig);figures.append(path)
    for phase in ('fixed','matched'):
        selected=[r for r in records if r['phase']==phase]
        if not selected:continue
        fig,axs=plt.subplots(3,3,figsize=(12,9),sharex=True,constrained_layout=True)
        for row_i,(arm,lam) in enumerate(ARMS):
            for r in selected:
                if r['arm']!=arm:continue
                data=list(csv.DictReader((out/r['case']/'common_grid.csv').open()));z=zs[1]
                for col,age in enumerate((22,26,30)):
                    vals=[v for v in data if int(v['birth_number'])==1 and float(v['age'])==age and int(v['income_index'])==z and float(v['wealth'])<=5 and v['gap_over_shock_scale']!='']
                    axs[row_i,col].plot([float(v['wealth']) for v in vals],[float(v['gap_over_shock_scale']) for v in vals],label='Baseline' if r['price_factor']==1 else 'Prices +10%',ls='-' if r['price_factor']==1 else '--')
                    axs[row_i,col].axhline(0,color='gray',lw=.5);axs[row_i,col].set_ylim(-20,20);axs[row_i,col].set_title(f'{LABEL[arm]}, age {age}')
                    axs[row_i,col].set_xlabel('Beginning liquid wealth');axs[row_i,col].set_ylabel('Attempt-minus-wait / shock scale')
        for row_i,(arm,_) in enumerate(ARMS):
            if not any(r['arm']==arm for r in selected):
                for ax in axs[row_i]:ax.text(.5,.5,'Not completed',ha='center',transform=ax.transAxes);ax.set_title(LABEL[arm])
        axs[0,0].legend(fontsize=8);fig.suptitle(f'{phase.title()} benefit: first-birth incentives by age, middle income state')
        path=out/f'{phase}_gaps_by_age.png';fig.savefig(path,dpi=140);plt.close(fig);figures.append(path)
        for birth in (1,2):
            fig,axs=plt.subplots(3,3,figsize=(12,9),sharex=True,constrained_layout=True)
            for row_i,(arm,lam) in enumerate(ARMS):
                for r in selected:
                    if r['arm']!=arm:continue
                    data=list(csv.DictReader((out/r['case']/'common_grid.csv').open()))
                    for col,z in enumerate(zs):
                        vals=[v for v in data if int(v['birth_number'])==birth and float(v['age'])==26 and int(v['income_index'])==z and float(v['wealth'])<=5 and v['gap_over_shock_scale']!='']
                        axs[row_i,col].plot([float(v['wealth']) for v in vals],[float(v['gap_over_shock_scale']) for v in vals],label='Baseline' if r['price_factor']==1 else 'Housing prices +10%',ls='-' if r['price_factor']==1 else '--')
                        axs[row_i,col].axhline(0,color='gray',lw=.5);axs[row_i,col].set_ylim(-20,20);axs[row_i,col].set_title(f'{LABEL[arm]}, {income_labels[z]}')
                        axs[row_i,col].set_xlabel('Beginning liquid wealth');axs[row_i,col].set_ylabel('Attempt-minus-wait / shock scale')
            for row_i,(arm,_) in enumerate(ARMS):
                if not any(r['arm']==arm for r in selected):
                    for ax in axs[row_i]:ax.text(.5,.5,'Not completed',ha='center',transform=ax.transAxes);ax.set_title(LABEL[arm])
            axs[0,0].legend(fontsize=8);fig.suptitle(f'{phase.title()} child benefit: child {birth}, current renters, age 26\nCommon state grid; endpoint gaps omitted, full values/classification in CSV')
            path=out/f'{phase}_birth_{birth}_gaps.png';fig.savefig(path,dpi=140);plt.close(fig);figures.append(path)
        for measure,label in (('attempt_probability','Attempt probability'),('conditional_renter_housing','Conditional renter housing'),('conditional_renter_consumption','Conditional renter consumption')):
            fig,axs=plt.subplots(1,3,figsize=(12,3.7),constrained_layout=True)
            for ax,(arm,lam) in zip(axs,ARMS):
                for r in selected:
                    if r['arm']!=arm:continue
                    data=list(csv.DictReader((out/r['case']/'common_grid.csv').open()));z=zs[1]
                    for birth in ((1,) if measure=='attempt_probability' else (1,2)):
                        vals=[v for v in data if int(v['birth_number'])==birth and float(v['age'])==26 and int(v['income_index'])==z and 0<=float(v['wealth'])<=5 and v['choice_status']!='unavailable' and (measure=='attempt_probability' or (float(v['conditional_renter_consumption'])>0 and float(v['conditional_renter_housing'])>0))]
                        ax.plot([float(v['wealth']) for v in vals],[float(v[measure]) for v in vals],label=f'{"Base" if r["price_factor"]==1 else "+10%"}, {birth-1} children',ls='-' if r['price_factor']==1 else '--')
                ax.set_title(LABEL[arm]);ax.set_xlabel('Beginning liquid wealth');ax.set_ylabel(label)
                if not any(r['arm']==arm for r in selected):ax.text(.5,.5,'Not completed',ha='center',transform=ax.transAxes)
            axs[0].legend(fontsize=7);fig.suptitle(f'{phase.title()} benefit, age 26, middle income state\nChild 1: current (n,m)=(0,0); child 2: (1,1). Housing is family-state policy, not event response.')
            path=out/f'{phase}_{measure}.png';fig.savefig(path,dpi=140);plt.close(fig);figures.append(path)
    cv=canvas.Canvas(str(out/'housing_fertility_cost_diagnostic.pdf'),pagesize=(842,595));pages=0
    styles=getSampleStyleSheet();style=styles['BodyText'];style.fontSize=7;style.leading=9
    def page(title,sub=''):
        nonlocal pages
        pages+=1;cv.setFont('Helvetica-Bold',16);cv.drawString(28,563,title);cv.setFont('Helvetica',9);cv.drawString(28,545,sub);cv.drawRightString(814,18,str(pages))
    def table(rows,widths,top=520,size=8):
        t=Table(rows,colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('FONTSIZE',(0,0),(-1,-1),size),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eef3')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));_,h=t.wrap(780,500);t.drawOn(cv,28,top-h)
    def fmt(v):
        if v=='':return ''
        try:v=float(v)
        except (TypeError,ValueError):return str(v)
        return f'{v:.3g}' if abs(v)>=10000 or 0<abs(v)<.001 else f'{v:.3f}'
    all_matched=all(v['status']=='matched' for v in complete['matching_status'].values())
    page('Housing costs, fertility and the first-child commitment',('All three matched comparisons complete' if all_matched else 'Partial diagnostic: larger-loading match unfinished')+'; no production adoption or policy forecast')
    for i,line in enumerate([
        'Three preference arms share control prices, earnings, pensions, taxes, estate rules, shocks, entry and grids.',
        'Both house asset prices and rents rise 10%; rent = fixed user-cost rate x house price. Credit primitives stay fixed.',
        'Current renters can optimize tenure: results separately show conditional renter branches and realized choices.',
        'Existing owners receive asset-price effects; their response is a separate diagnostic, not pure affordability.',
        'Child utility is b*m^0.86. Fixed-benefit trials hold b; matched trials change only b to baseline fertility 2.1.',
        'The share utility reference-rent normalizer remains fixed at its inherited value throughout.',
        'The floor is ongoing housing need. A separate inherited first-birth utility cost of 0.507 stays fixed.',
        'Gap = attempt now minus wait, excluding current taste draw, retaining future shock-inclusive values.',
        'Negative gaps do not mean children have no direct benefit or that lifetime childlessness is preferred.',
        'Stationary PE distributions adjust to prices. Market/fiscal/replacement imbalances are reported, not cleared.',
        'The birth-room statistic remains an unmatched stationary proxy; historical target system retained for comparison.',
        'Full target/parameter tables and all 17 standard stationary diagnostics per reported solution follow.',
        ('All three positive child-benefit matches meet the declared fertility tolerance.' if all_matched else 'The time budget stopped +20-point matching before acceptance. No missing result is imputed.')]):
        cv.setFont('Helvetica',10);cv.drawString(28,513-i*27,line)
    cv.showPage()
    page('Direct comparison of housing-cost sensitivity','Permanent house-price and rent increase of 10%; matched trials hold their newly matched benefit under the shock')
    rows=[['Benefit treatment','Arm','Baseline fertility','Fertility change (%)','Birth-age change']]
    for base in records:
        if base['price_factor']!=1:continue
        shocked=next((r for r in records if r['phase']==base['phase'] and r['arm']==base['arm'] and r['price_factor']==1.1),None)
        if shocked is None:continue
        bm={r['moment']:r['model'] for r in base['target_fit']};sm={r['moment']:r['model'] for r in shocked['target_fit']}
        rows.append([base['phase'].title(),LABEL[base['arm']],fmt(base['fertility']),fmt(100*(shocked['fertility']/base['fertility']-1)),fmt(sm['nchs_mean_age']-bm['nchs_mean_age'])])
    table(rows,[150,160,150,160,160])
    for i,line in enumerate([
        'Fixed-benefit comparisons combine the housing mechanism with different starting fertility levels.',
        'Matching mean fertility removes that difference, but does not match every hazard or household distribution.',
        'Larger price sensitivity alone is not empirical validation; no target for this elasticity is imposed here.',
        'Higher mean birth age is an endpoint timing shift, not an identified transition path of postponed births.',
        f'Common baseline house price {records[0]["price"]:.3f}; rent {records[0]["rent"]:.3f}; both multiplied by 1.100.']):
        cv.setFont('Helvetica',10);cv.drawString(28,280-i*26,line)
    cv.showPage()
    if (out/'common_state_changes.csv').exists():
        page('Matched benefits: the same households face higher housing costs','Childless ages 18-34; common control weights. Owners also receive asset revaluation; tenure is optimized.')
        rows=[['Arm','Initial tenure','Birth probability\nchange (pp)','Owner probability\nchange (pp)','Housing\nchange','Nonhousing\nchange']]
        for s in csv.DictReader((out/'common_state_changes.csv').open()):
            if s['phase']!='matched' or s['group'] not in ('current_renters','current_owners'):continue
            rows.append([LABEL[s['arm']],s['group'].replace('current_',''),fmt(100*float(s['birth_probability_change'])),fmt(100*float(s['realized_ownership_change'])),fmt(s['realized_housing_change']),fmt(s['realized_nonhousing_change'])])
        table(rows,[150,105,140,140,120,135])
        for i,line in enumerate(['These differences hold beginning states fixed, separating choice responses from stationary composition.',
            'Normalized gap averages can be dominated by tiny probabilities and exclude endpoint states.',
            'Use realized birth responses to compare exposure; gap plots describe conditional choice incentives.',
            'Income/wealth plots use ages 22, 26 and 30; this common-weight table uses the full eligible 18-34 group.']):
            cv.setFont('Helvetica',10);cv.drawString(28,265-i*26,line)
        cv.showPage()
    for phase in ('fixed','matched'):
        selected=[r for r in records if r['phase']==phase]
        if not selected:continue
        page(phase.title()+' child benefit: price sensitivity','Rows are stationary partial equilibria; mean age is the inherited birth-timing observer')
        rows=[['Arm','Price factor','b','Fertility','Childless %','Birth age','Rooms jump','Market gap %']]
        for r in selected:
            m={v['moment']:v['model'] for v in r['target_fit']}
            rows.append([LABEL[r['arm']],fmt(r['price_factor']),fmt(r['b']),fmt(r['fertility']),fmt(100*m['cps_childlessness']),fmt(m['nchs_mean_age']),fmt(m['first_birth_rooms']),fmt(100*r['gates']['market_residual'])])
        table(rows,[135,90,80,85,100,95,90,105]);cv.showPage()
        page(phase.title()+': identical young renter states','Ages 18-34, no prior children; fixed control weights. Tenure is fully optimized.')
        rows=[['Arm','Price factor','Birth probability','Realized owner %','Housing','Nonhousing']]
        for r in selected:
            s=next(v for v in csv.DictReader((out/r['case']/'common_states.csv').open()) if v['group']=='current_renters')
            rows.append([LABEL[r['arm']],fmt(r['price_factor']),fmt(s['birth_probability']),fmt(100*float(s['realized_ownership'])),fmt(s['realized_housing']),fmt(s['realized_nonhousing'])])
        table(rows,[150,110,130,130,110,150]);cv.showPage()
        page(phase.title()+': direct child benefit and current incentives','Case-specific stationary weights; positive-gap birth shares and risk-weighted median gaps are different objects')
        rows=[['Arm','Price factor','b','First births at\npositive gaps (%)','Later births at\npositive gaps (%)','First risk median\ngap / shock scale']]
        for r in selected:
            first,later=r['incentives'][:2]
            rows.append([LABEL[r['arm']],fmt(r['price_factor']),fmt(r['b']),fmt(100*first['birth_share_positive_interior']),fmt(100*later['birth_share_positive_interior']),fmt(first['gap_over_kappa_median'])])
        table(rows,[145,100,85,150,150,150]);cv.showPage()
    page('Budget, authentication and limitations','Source pins, solve history, loop smoke and receipts are retained beside this report')
    lines=[f'Native evaluations: {complete["evaluations"]}; elapsed before reporting: {complete["elapsed_seconds"]:.1f} seconds.',f'Matching: positive b bracket {B_BRACKET}; tolerance {TOL}; original cap 12; continuation cap {complete.get("max_b_evaluations",MAX_B_EVAL)}.']
    lines += [f'{LABEL[arm]}: {value["status"]}; baseline fertility {value.get("fertility",value.get("last_fertility",float("nan"))):.3f}.' for arm,value in complete['matching_status'].items()]
    lines += ['All reported cases passed unchanged household probability, budget, value and transition checks.', 'The control is numerically replayed before new experiments. No source model edits or relaxed scientific gates.', 'Fiscal inputs remain frozen; actual fiscal residuals are measured in the final verification receipt.', 'Liquid wealth <=1 is an illustrative model-unit subset, not an empirical bottom quantile.', 'Plots clip normalized gaps visually to [-20,20]; endpoint probabilities are classified, never inverted.', 'Stationary endpoints do not identify a transition delay or separate all cohort timing effects.', 'Matched mean fertility does not match the first-birth hazard, childlessness or full wealth distribution.']
    for i,line in enumerate(lines):cv.setFont('Helvetica',9);cv.drawString(28,512-i*27,line)
    cv.showPage()
    for fig in figures:page(fig.stem.replace('_',' '));cv.drawImage(ImageReader(str(fig)),28,45,width=780,height=480,preserveAspectRatio=True,anchor='c');cv.showPage()
    for r in records:
        page(r['case']+': full target fit','Descriptive inherited objective, not a recalibrated or market-clearing fit')
        rows=[['Moment','Target','Model','Gap','Weight','Contribution']]+[[v['moment']]+[fmt(v[k]) for k in ('target','model','gap','weight','loss_contribution')] for v in r['target_fit']]
        table(rows,[205,105,105,105,120,140]);cv.showPage()
        params=list(csv.DictReader((out/r['case']/'parameters.csv').open()))
        for start in range(0,len(params),16):
            page(r['case']+': parameter restrictions')
            rows=[['Parameter','Value','Lower','Upper','Near bound','Status']]
            for p in params[start:start+16]:rows.append([Paragraph(p['parameter'].replace('_',' '),style),fmt(p['estimate']),fmt(p['lower']),fmt(p['upper']),p['near_bound'],Paragraph(p['status'],style)])
            table(rows,[180,75,55,55,70,345],size=7);cv.showPage()
    for r in records:
        figs=sorted((out/r['case']/'standard_diagnostics').glob('*.png'))
        for start in range(0,len(figs),2):
            page(r['case']+': stationary PE diagnostics','Standard graph set; prices are externally fixed, market clearing is not imposed')
            for j,f in enumerate(figs[start:start+2]):cv.drawImage(ImageReader(str(f)),28+395*j,45,width=385,height=475,preserveAspectRatio=True,anchor='c')
            cv.showPage()
    cv.save()
    import fitz
    from PIL import Image,ImageDraw
    qa=out/'pdf_qa';qa.mkdir(exist_ok=True);doc=fitz.open(out/'housing_fertility_cost_diagnostic.pdf');outside=[]
    for start in range(0,len(doc),16):
        sheet=Image.new('RGB',(1360,1040),'#ddd');draw=ImageDraw.Draw(sheet)
        for k in range(start,min(start+16,len(doc))):
            p=doc[k]
            for b in p.get_text('blocks'):
                if b[0]<0 or b[1]<0 or b[2]>843 or b[3]>596:outside.append((k+1,b[:4]))
            pix=p.get_pixmap(matrix=fitz.Matrix(.39,.39),alpha=False);im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);x=(k-start)%4*340;y=(k-start)//4*260;sheet.paste(im,(x,y+20));draw.text((x+4,y+4),str(k+1),fill='black')
            if k<6:p.get_pixmap(matrix=fitz.Matrix(1.2,1.2),alpha=False).save(str(qa/f'page_{k+1}.png'))
        sheet.save(qa/f'contact_{start//16+1}.png')
    common.write(out/'pdf_receipt.json',dict(pages=len(doc),outside_text=outside,sha256=hashlib.sha256((out/'housing_fertility_cost_diagnostic.pdf').read_bytes()).hexdigest()))
    assert not outside,outside


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--stage',choices=['smoke','run','render'],required=True);a=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Torch allocation required')
    if a.stage=='render':render(a.output);return
    a.output.mkdir(parents=True,exist_ok=False);start=time.monotonic();rt,tax,obj,selected=setup(a.output)
    control=selected['packet'];P=control['parameters'];grid=control['b_grid'];model=rt['model']
    common.write(a.output/'smoke.json',loop_smoke(P,grid,model))
    common.write(a.output/'design.json',dict(arms=ARMS,bracket=B_BRACKET,target=TARGET,tolerance=TOL,max_b_evaluations=MAX_B_EVAL,contract=common.PIN,source_pins=selected['pins'],driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),frozen_reference_rent=P.utility_reference_rent,baseline_price=np.asarray(control['solution'].p_eq).tolist(),fixed_pension=P.pension,fixed_payroll_tax=P.tau_pay,job=os.environ['SLURM_JOB_ID']))
    if a.stage=='smoke':print('EXACT LOOP SMOKE PASSED',flush=True);return
    state={'phase':'control replay'};done=threading.Event()
    def heartbeat():
        while not done.wait(30):common.write(a.output/'heartbeat.json',dict(state,elapsed_seconds=time.monotonic()-start))
    threading.Thread(target=heartbeat,daemon=True).start()
    replay=solve(rt,copy.deepcopy(P),grid,np.asarray(control['solution'].p_eq));deltas={}
    for key in ('V','c_pol','hR_pol','bp_pol','fert_probs','fert2_probs','g'):
        x=np.asarray(getattr(replay['solution'],key));y=np.asarray(getattr(control['solution'],key));deltas[key]=float(np.max(np.abs(x-y)));np.testing.assert_allclose(x,y,atol=1e-10,rtol=1e-12)
    common.write(a.output/'control_replay.json',deltas);del replay
    records=[];history=[];matches={};count=1;baseprice=np.asarray(control['solution'].p_eq)
    def evaluate(arm,lam,b,factor,label,full=False):
        nonlocal count
        if time.monotonic()-start>1950:raise TimeoutError('Stop before report reserve')
        assert count<48
        state.update(phase=label,evaluations=count);out=a.output/label;out.mkdir()
        Q=loading.configure(P,arm,lam);Q.psi_child=float(b)
        packet=solve(rt,Q,grid,baseprice*factor);count+=1;tfr=float(rt['chain'].extract_moments(packet['solution'],Q)['tfr'])
        entry=dict(case=label,arm=arm,b=float(b),factor=factor,fertility=tfr,elapsed_seconds=time.monotonic()-start);history.append(entry);common.table(a.output/'solve_history.csv',history);common.write(a.output/'latest_completed.json',entry)
        print(json.dumps(entry),flush=True)
        if full:
            phase='matched' if label.startswith('matched') else 'fixed'
            rec=report(packet,rt,tax,obj,selected['case'],out,arm,lam,phase,factor,control);records.append(rec);common.write(a.output/'comparison_summary.json',records)
        return packet,tfr,out
    for arm,lam in ARMS:
        for factor in (1.,1.1):evaluate(arm,lam,P.psi_child,factor,f'fixed_{arm}_{factor:.1f}',True)
    for arm,lam in ARMS:
        lo,hi=B_BRACKET;p_lo,tlo,_=evaluate(arm,lam,lo,1.,f'matchsearch_{arm}_00');del p_lo
        p_hi,thi,_=evaluate(arm,lam,hi,1.,f'matchsearch_{arm}_01');del p_hi
        if not tlo<=TARGET<=thi:matches[arm]={'status':'unbracketed','low':tlo,'high':thi};continue
        matched=None
        for it in range(2,MAX_B_EVAL):
            b=(lo+hi)/2;packet,tfr,folder=evaluate(arm,lam,b,1.,f'matchsearch_{arm}_{it:02d}')
            if abs(tfr-TARGET)<=TOL:matched=(b,packet,folder);break
            if tfr<TARGET:lo=b
            else:hi=b
        if matched is None:matches[arm]={'status':'budget_exhausted','last_fertility':tfr};continue
        b,packet,folder=matched;out=a.output/f'matched_{arm}_1.0';out.mkdir()
        records.append(report(packet,rt,tax,obj,selected['case'],out,arm,lam,'matched',1.,control))
        matches[arm]={'status':'matched','b':b,'fertility':tfr,'search_evaluations':it+1}
        evaluate(arm,lam,b,1.1,f'matched_{arm}_1.1',True)
    common.write(a.output/'complete.json',dict(status='completed',records=records,matching_status=matches,evaluations=count,elapsed_seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID']))
    state['phase']='render';render(a.output);done.set();print('COMPLETE',flush=True)

if __name__=='__main__':main()
