#!/usr/bin/env python3
"""Verified joint-choice policy paths from the exact fitted 2023 population.

The historical calibration uses its pinned Census bridge. Post-2023 paths use
the maintained finite-horizon closed-national experiment (M=0, rho=1), with
current prices treated as permanent at each date. No perfect-foresight claim.
"""
from __future__ import annotations
import argparse, copy, gzip, json, math, pickle, sys, time, traceback, textwrap
from pathlib import Path
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import numpy as np
import run_e5f_joint_overnight_case as adapter
import run_e5f_independent_numerical_audit as audit
import run_e5f_post2023_policy_mechanisms as policy
import run_e5f_transition_calibration as calibration
from intergen_eqscale_seq_optimized import solver

CASES=('baseline','supply-plus-20','dependent-child-ltv95','property-tax-2pct-no-rebate')


def selected_path(path):
    value=adapter.read_json(path)
    if 'best_candidate' in value:return Path(path)
    return Path(value['best']['summary'])


def verify_selected(path,contract):
    s=adapter.read_json(path);out=path.parent;receipt=adapter.read_json(out/'case_receipt.json')
    planpath=out.parent/'plan.json';plan=adapter.load_plan(planpath,receipt['plan_sha256'])
    case=next(x for x in plan['cases'] if x['id']==receipt['case_id'])
    adapter.validate_result(out,plan,case)
    for rel,sha in receipt['artifact_sha256'].items():adapter.verify(out/rel,sha)
    if s['code_fingerprints']['bundle_sha256'] != contract['code_bundle_sha256']:
        raise RuntimeError('Selected policy-path calibration has different scientific code')
    if calibration.code_fingerprint_contract(solver)['bundle_sha256'] != contract['code_bundle_sha256']:
        raise RuntimeError('Policy path scientific source has drifted')
    for name,sha in contract['base_plan']['helper_sha256'].items():adapter.verify(ROOT/'code/model/tools'/name,sha)
    if s['target_fingerprint'] != adapter.TARGET:raise RuntimeError('Target contract drift')
    return s,receipt


def prepare(path,s):
    packet=audit.load_checkpoint(path.parent/'dated_state.pkl.gz')
    P=packet['parameters'];e=packet['evaluation']
    if not P.joint_nested_choice or e.policy.joint_choice is None:
        raise RuntimeError('Selected checkpoint has no joint household policies')
    inherited = getattr(e, 'inherited_g_pre', None)
    if inherited is None:
        raise RuntimeError('Checkpoint lacks the original inherited population; rebuild it')
    if inherited.shape != e.g_pre.shape or not np.isfinite(inherited).all() or inherited.min() < 0:
        raise RuntimeError('Invalid original inherited population')
    replay, projected = policy.calendar.gate_pre_fertility_distribution(
        inherited, e.policy, P, packet['b_grid'], packet['shared'])
    if not np.array_equal(replay, e.g_pre) or projected != e.feasibility_projection_mass:
        raise RuntimeError('Original inherited population does not reproduce the fitted feasibility gate')
    rows=adapter.read_csv(path.parent/'cases'/s['best_candidate']['candidate']/'transition_path.csv')
    if len(rows)!=5 or [int(float(x['calendar_year'])) for x in rows] != [2007,2011,2015,2019,2023]:
        raise RuntimeError('Expected complete five-date fitted history')
    # A path row saves the queue AFTER that date advances. The 2019 row
    # therefore contains exactly the queue inherited at the start of 2023.
    queue=json.loads(rows[3]['birth_queue_scheduled_flows'])
    raw=json.loads(rows[3]['birth_queue_raw_state_scheduled_flows'])
    if len(queue)!=4 or len(raw)!=4 or min(queue+raw)<0:raise RuntimeError('Invalid inherited birth queue')
    state=policy.baseline.DynamicState(inherited.copy(),queue,raw,float(e.policy.price[0]),None)
    initial_mass=float(rows[0]['adult_population'])
    prepared=SimpleNamespace(b_grid=packet['b_grid'],supply_rule=packet['supply_rule'],initial_mass_2007=initial_mass)
    return packet,prepared,state,rows


def report(out,selected,contract,smoke=False):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    s=adapter.read_json(selected);fit=adapter.read_csv(selected.parent/'target_fit_long.csv')
    params=adapter.read_csv(selected.parent/'parameter_table.csv')
    receipt=adapter.read_json(out/'equilibrium_receipt.json')
    paths={name:adapter.read_csv(out/name/'policy_path.csv') for name in CASES if (out/name/'policy_path.csv').exists()}
    effects=[]
    if 'baseline' in paths:
        base=paths['baseline']
        for name,rows in paths.items():
            if name=='baseline':continue
            for k in (0,len(base)-1):
                if len(rows)!=len(base):raise RuntimeError('Mismatched policy horizons')
                b,r=base[k],rows[k]
                effects.append(dict(policy=name,year=int(float(r['calendar_year'])),
                    births_percent=100*(float(r['birth_children_topcode_adjusted'])/float(b['birth_children_topcode_adjusted'])-1),
                    ownership_pp=100*(float(r['owner_rate'])-float(b['owner_rate'])),
                    rooms_percent=100*(float(r['housing_demand_per_adult'])/float(b['housing_demand_per_adult'])-1)))
    policy.baseline.write_csv(out/'policy_effects.csv',effects)
    worst=max(fit,key=lambda x:float(x['loss_contribution']))
    first=next(x for x in fit if x['moment']=='housing_increment_0to1')
    searchroot=Path(contract['output_root'])/'search'
    exact=adapter.read_json(searchroot/'final_verification.json') if (searchroot/'final_verification.json').exists() else {}
    jac=adapter.read_json(searchroot/'jacobian_diagnostics.json') if (searchroot/'jacobian_diagnostics.json').exists() else {}
    message=[f"Experimental joint-choice model: objective {s['best_candidate']['transition_loss']:.3f} on the original twelve targets.",
        f"First-birth housing response: {float(first['model']):.3f} rooms; target {float(first['target']):.3f}.",
        f"Largest loss contribution: {worst['moment']} ({float(worst['loss_contribution']):.3f}).",
        f"Two exact final repetitions: {'passed' if exact.get('status')=='pass' else 'not yet certified'}. Local identification: {jac.get('status','not yet assessed')}.",
        f"Equilibrium paths: {receipt['status']}. Markets clear separately at each date; households treat current prices as permanent.",
        "Joint shocks are simultaneous. Tenure is committed before conception; housing size, consumption and saving may respond within that tenure.",
        "One common nesting parameter applies to all birth orders. This is an experimental restriction, and the production calibration remains unchanged."]
    if smoke:message.insert(0,'SMOKE TEST: these are verification cases, not the result of the overnight search.')
    lines=['# Simultaneous-choice calibration and equilibrium paths','',*message,'',
        'Full target-fit table: target_fit.csv. Full parameter table: parameters.csv.',
        'Policy paths retain the closed-national finite-horizon benchmark: no outside entry after 2023, retention one, inherited four-slot birth queue and 2.1 replacement conversion. Supply elasticity remains 0.63; property-tax revenue is discarded; no rebates or grants.',
        'Every completed date has the standard seventeen diagnostic graphs and a saved checkpoint.','',
        'Regenerate this readout with the same command and --report-only.']
    (out/'MORNING_SUMMARY.md').write_text('\n\n'.join(lines)+'\n')
    (out/'target_fit.csv').write_bytes((selected.parent/'target_fit_long.csv').read_bytes())
    (out/'parameters.csv').write_bytes((selected.parent/'parameter_table.csv').read_bytes())
    def figpage(pdf,fig):
        pdf.savefig(fig);plt.close(fig)
    with PdfPages(out/'joint_nested_readout.pdf') as pdf:
        fig=plt.figure(figsize=(11.7,8.3));fig.text(.06,.92,'Simultaneous-choice calibration',fontsize=22,weight='bold')
        y=.84
        for text in message:
            wrapped=textwrap.fill(text,110);fig.text(.07,y,wrapped,fontsize=11,va='top');y-=.040*(wrapped.count('\n')+1)+.025
        fig.text(.07,.04,'Experimental research readout | All displayed results are subject to the saved verification receipts.',fontsize=9)
        figpage(pdf,fig)
        for title,rows,columns,labels in [
            ('Complete target fit',fit,['moment','target','model','gap','weight','loss_contribution'],['Moment','Target','Model','Gap','Weight','Loss']),
            ('Parameters and restrictions',params,['parameter','value','lower_bound','upper_bound','is_free_parameter','near_bound'],['Parameter','Estimate','Lower','Upper','Estimated','Near bound'])]:
            fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title(title,fontsize=18,pad=20)
            body=[]
            for row in rows:
                values=[]
                for k in columns:
                    v=row[k]
                    if k==columns[0]:v=textwrap.fill(v.replace('_',' '),39)
                    else:
                        try:
                            x=float(v);v=f'{x:.6g}' if math.isfinite(x) else '—'
                        except (TypeError,ValueError):pass
                    values.append(v)
                body.append(values)
            table=ax.table(cellText=body,colLabels=labels,cellLoc='right',colWidths=[.43,.115,.115,.115,.115,.105],loc='center')
            table.auto_set_font_size(False);table.set_fontsize(8.2)
            for (r,c),cell in table.get_celld().items():
                cell.set_height(.052 if title.startswith('Parameters') else .061)
                if c==0:cell.get_text().set_ha('left')
                if r==0:cell.set_facecolor('#dce8ec');cell.get_text().set_weight('bold')
            fig.subplots_adjust(left=.06,right=.94,top=.88,bottom=.06);figpage(pdf,fig)
        if paths:
            fig,axes=plt.subplots(2,3,figsize=(11.7,8.3),constrained_layout=True)
            fields=[('asset_price','House price'),('housing_demand_per_adult','Rooms per household'),('owner_rate','Owner share'),
                    ('birth_children_topcode_adjusted','Births (top-code adjusted)'),('topcode_adjusted_births_per_adult','Births per household'),('population_index_2023','Population, 2023 = 1')]
            for ax,(field,label) in zip(axes.flat,fields):
                for name,rows in paths.items():ax.plot([float(r['calendar_year']) for r in rows],[float(r[field]) for r in rows],label=name.replace('-',' '))
                ax.set_title(label);ax.grid(alpha=.2)
            axes.flat[0].legend(fontsize=7);fig.suptitle('Market-clearing paths: closed-national experiment',fontsize=17);figpage(pdf,fig)
        if (searchroot/'jacobian_supplemental.png').exists():
            fig=plt.figure(figsize=(11.7,8.3));ax=fig.add_axes([.03,.06,.94,.88]);ax.imshow(plt.imread(searchroot/'jacobian_supplemental.png'));ax.axis('off');figpage(pdf,fig)
    return effects


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--selected-summary',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--contract',type=Path,required=True);ap.add_argument('--smoke',action='store_true');ap.add_argument('--report-only',action='store_true')
    a=ap.parse_args();out=a.outdir.resolve();out.mkdir(parents=True,exist_ok=True)
    contract=adapter.read_json(a.contract);selected=selected_path(a.selected_summary)
    summary,selected_receipt=verify_selected(selected,contract)
    if a.report_only:report(out,selected,contract,a.smoke);return
    if (out/'equilibrium_receipt.json').exists():raise RuntimeError('Refusing to overwrite completed equilibrium paths')
    policy.transition.configure_sequential_model()
    policy.calendar.apply_fertility=policy.transition.apply_sequential_fertility
    policy.calendar.advance_calendar_distribution=policy.transition.advance_sequential_calendar_distribution
    packet,prepared,state,history=prepare(selected,summary)
    adapter.write_json(out/'inherited_state_verification.json', dict(
        status='exact_feasibility_replay', source_summary_sha256=adapter.digest(selected),
        inherited_population_sha256=policy.baseline.array_sha256(state.g_pre),
        fitted_gated_population_sha256=policy.baseline.array_sha256(packet['evaluation'].g_pre),
        inherited_to_gated_l1=float(np.abs(state.g_pre-packet['evaluation'].g_pre).sum()),
        fitted_projection_mass=float(packet['evaluation'].feasibility_projection_mass),
        birth_queue_source_year=2019, scheduled_entries=state.scheduled_entries,
        scheduled_raw_entries=state.scheduled_raw_entries))
    post=1 if a.smoke else 10;start=time.time();receipts={};failures={}
    for name in CASES:
        folder=out/name;folder.mkdir(parents=True,exist_ok=True);count=0
        original=policy.baseline.evaluate_state
        def observed(*args,**kwargs):
            nonlocal count
            evaluation,shared,fallback=original(*args,**kwargs)
            P=args[1];year=2023+4*count;target=folder/f'date_{year}';target.mkdir(parents=True,exist_ok=True)
            current=dict(parameters=P,b_grid=prepared.b_grid,evaluation=evaluation,shared=shared,supply_rule=policy.policy_supply_rule(prepared.supply_rule,policy.POLICIES[name]))
            if name=='baseline' and count==0:
                gap=float(np.abs(evaluation.g_current-packet['evaluation'].g_current).sum())
                if gap>2e-10:raise RuntimeError(f'Baseline does not reproduce selected 2023 state: {gap}')
            audit.standard_diagnostics(current,target,validate_production_young=False)
            arrays=audit.policy_array_audit(current,target);budget=audit.budget_audit(current,target)
            if arrays['occupied_negative_steps'] != 0:raise RuntimeError('Occupied value monotonicity failed on policy path')
            for bounds in arrays['probabilities'].values():
                if bounds['nonfinite'] or bounds['minimum']<0 or bounds['maximum']>1:raise RuntimeError('Policy probability gate failed')
            if budget['budget_excess_mass']>2e-10:raise RuntimeError('Material occupied budget violation on policy path')
            with gzip.open(target/'dated_state.pkl.gz','wb',compresslevel=1) as stream:pickle.dump(current,stream,protocol=5)
            count+=1
            return evaluation,shared,fallback
        policy.baseline.evaluate_state=observed
        try:
            rows,gates=policy.run_policy_path(prepared,state,packet['parameters'],float(state.g_pre.sum()),policy.POLICIES[name],
                post_2023_periods=post,market_tol=2e-4,market_max_iter=60,progress_dir=folder)
            if gates['maximum_market_residual']>2e-4 or gates['maximum_mass_residual']>2e-10:
                raise RuntimeError('Equilibrium-path market or mass gate failed')
            policy.baseline.write_csv(folder/'policy_path.csv',rows);policy.make_case_figure(rows,folder)
            receipts[name]=dict(status='complete',dates=len(rows),gates=gates,source_summary_sha256=adapter.digest(selected))
            adapter.write_json(folder/'receipt.json',receipts[name])
        except Exception as error:
            failures[name]=dict(error=str(error),traceback=traceback.format_exc())
            adapter.write_json(folder/'failure.json',failures[name])
            if 'Housing market did not clear' not in str(error) and 'market gate failed' not in str(error):raise
        finally:policy.baseline.evaluate_state=original
    receipt=dict(status='complete' if not failures else 'partial_policy_failures',smoke=a.smoke,elapsed_seconds=time.time()-start,
        inherited_state_verification_sha256=adapter.digest(out/'inherited_state_verification.json'),
        cases=receipts,failures=failures,selected_summary=str(selected),selected_summary_sha256=adapter.digest(selected),
        scientific_bundle=contract['code_bundle_sha256'],target_fingerprint=adapter.TARGET,
        expectations='temporary equilibrium: each current price treated as permanent',population_closure='closed national: M=0, rho=1',
        replacement_conversion=1/2.1,supply_elasticity=.63,fiscal_closure='tax revenue discarded; no rebate or grant',production_promoted=False)
    adapter.write_json(out/'equilibrium_receipt.json',receipt);report(out,selected,contract,a.smoke)
    if a.smoke and failures:raise RuntimeError('Policy-loop smoke has rejected cases')
    print(json.dumps(dict(status=receipt['status'],pdf=str(out/'joint_nested_readout.pdf'))),flush=True)
if __name__=='__main__':main()
