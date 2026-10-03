"""No-solve paired Estate-A readout, with complete fit/parameter and estate diagnostics."""
from __future__ import annotations
import argparse,copy,csv,json,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
PROJECT=Path(__file__).resolve().parents[4]
OUTPUT=PROJECT/'output/model/experiments/birth_count_choice/estate_a_v1'
BASELINE=PROJECT/'output/model/local_solution/cases/20261003T175652812716Z_b1c72f13'
CLAUDE_A=PROJECT/'output/model/experiments/estate_valuation_20261003/cases/net_sale'
BASELINE_MULTIPLE=PROJECT/'output/model/experiments/birth_count_choice/current_params_v1/cases/20261003T203001912131Z_812b878e'


def estate_diagnostics(result):
    from model.engine.distribution import add_aggregate_wealth_bequest_flow_moments
    from model.engine.utils import weighted_quantile
    P=copy.deepcopy(result.P); sol=result.solution
    asset=np.asarray(sol.g_beginning_distribution); death=np.asarray(sol.g)
    prices=np.full(int(P.I),result.price)
    rows=[]
    for name,net in [('gross_housing',False),('net_selling_cost',True)]:
        Q=copy.deepcopy(P); Q.estate_flow_net_of_selling_cost=net; Q.estate_receiver='none'
        if bool(getattr(Q,'native_due_stayer_credit',False)):
            Q._g_stay_distribution=sol.g_stay_distribution; Q._bp_pol_stay=sol.bp_pol_stay
        stats=SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(stats,asset,death,sol.bp_pol,Q,result.b_grid,prices)
        rows += [dict(measure='annual_positive_death_estate_flow',definition=name,value=stats.annual_bequest_flow),
                 dict(measure='annual_positive_death_estate_flow_to_living_gross_wealth',definition=name,value=stats.annual_bequest_flow_to_aggregate_wealth)]
        vals,wts=[],[]; signed=positive=negative=mass=0.
        for j in range(int(P.J)):
            age=float(P.age_start)+j*float(P.da)
            if age<65: continue
            hazard=(1-float(P.survival_probs[j])) if j<int(P.J)-1 and P.use_age_survival else (1. if j==int(P.J)-1 else 0.)
            for ten in range(death.shape[1]):
                h=prices[0]*float(P.H_own[ten-1])*(1-float(P.psi) if net else 1.) if ten else 0.
                v=sol.bp_pol[:,ten,:,j]+h; w=death[:,ten,:,j]
                if bool(getattr(P,'native_due_stayer_credit',False)):
                    stay=np.asarray(sol.g_stay_distribution[:,ten,:,j]); sv=sol.bp_pol_stay[:,ten,:,j]+h
                    w=w-stay; vv=np.r_[v.ravel(),sv.ravel()]; ww=np.r_[w.ravel(),stay.ravel()]
                else: vv=v.ravel(); ww=w.ravel()
                if ww.min() < -1e-10: raise RuntimeError('Invalid death mixture')
                ww=np.maximum(ww,0)*hazard
                use=ww>0; vals.append(vv[use]); wts.append(ww[use])
                signed+=float(np.sum(vv*ww)); positive+=float(np.sum(np.maximum(vv,0)*ww)); negative+=float(np.sum(np.minimum(vv,0)*ww)); mass+=float(ww.sum())
        values=np.concatenate(vals); weights=np.concatenate(wts)
        for measure,value in [('old65plus_signed_death_estate_flow',signed/float(P.period_years)),('old65plus_positive_death_estate_flow',positive/float(P.period_years)),('old65plus_negative_death_estate_flow',negative/float(P.period_years)),('old65plus_mean_estate_per_death',signed/mass),('old65plus_median_estate_per_death',float(weighted_quantile(values,weights,.5)))]:
            rows.append(dict(measure=measure,definition=name,value=value))
    age_index=int(round((82-float(P.age_start))/float(P.da)))
    if float(P.age_start)+age_index*float(P.da)!=82: raise RuntimeError('Age 82 is not a saved model age')
    weight=np.asarray(death[:,:,:,age_index]); mass=float(weight.sum())
    ownership=float(weight[:,1:].sum())/mass
    saving=np.asarray(sol.bp_pol[:,:,:,age_index]); total=float(np.sum(weight*saving))
    if bool(getattr(P,'native_due_stayer_credit',False)):
        stay=np.asarray(sol.g_stay_distribution[:,:,:,age_index])
        stay_saving=np.asarray(sol.bp_pol_stay[:,:,:,age_index])
        total+=float(np.sum(stay*(stay_saving-saving)))
    rows += [dict(measure='ownership_age82',definition='actual_post_tenure_living_mass',value=ownership),
             dict(measure='mean_post_saving_bp_age82',definition='actual_post_tenure_stayer_corrected',value=total/mass)]
    return rows


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,default=BASELINE)
    parser.add_argument('--baseline-multiple',type=Path,default=BASELINE_MULTIPLE)
    parser.add_argument('--single',type=Path,default=OUTPUT/'single/latest')
    parser.add_argument('--multiple',type=Path,default=OUTPUT/'multiple/latest')
    parser.add_argument('--out',type=Path,default=OUTPUT)
    args=parser.parse_args(argv);sys.path.insert(0,str(Path(__file__).resolve().parent))
    from summarize import read_case,write_csv,markdown_table
    from model.storage import load_case
    from model.estate_contract import rescore_rows,contract,experiment_flags
    cases={};allfits=[];allparams=[];estates=[]
    newcontract,tp,wp,_=contract()
    coordinates=None
    for arm,path in [('baseline',args.baseline),('baseline_multiple',args.baseline_multiple),('single',args.single),('multiple',args.multiple)]:
        case,meta,oldfits,params,_,_=read_case(path)
        if coordinates is None: coordinates=meta['parameters']
        if meta['parameters']!=coordinates: raise RuntimeError('Fixed parameter coordinates differ')
        if arm in ('single','multiple'):
            flags=json.loads((case/'input_contract.json').read_text())['experiment_flags']
            if flags!=experiment_flags(1 if arm=='single' else 3): raise RuntimeError('Saved Estate-A flags differ')
            for directory,count in [('policy_plots',8),('aggregate_plots',7)]:
                if len(list((case/directory).glob('*.png')))!=count: raise RuntimeError('Incomplete cached plots')
        fits,residual=rescore_rows(oldfits)
        result,_=load_case(case)
        diagnostics=estate_diagnostics(result)
        cases[arm]=dict(case=str(case),birth_cap=(3 if arm in ('baseline_multiple','multiple') else 1),estate_a=(arm in ('single','multiple')),effective_estate_flags={key:getattr(result.P,key,False) for key in ('bequest_net_of_selling_cost','estate_flow_net_of_selling_cost')},price=meta['price'],new_contract_loss=float(residual@residual),old_contract_loss=sum(float(r['loss_contribution']) for r in oldfits if r['role']=='scored'),closure=meta['closure'])
        allfits += [dict(arm=arm,**r) for r in fits]
        allparams += [dict(arm=arm,**r) for r in params]
        estates += [dict(arm=arm,**r) for r in diagnostics]
    args.out.mkdir(parents=True,exist_ok=True)
    write_csv(args.out/'paired_target_fit_new_contract.csv',allfits)
    write_csv(args.out/'paired_parameters.csv',allparams)
    write_csv(args.out/'paired_estate_diagnostics.csv',estates)
    (args.out/'comparison.json').write_text(json.dumps(dict(status='completed_fixed_parameters_no_recalibration',cases=cases,target_fingerprint=tp,weight_fingerprint=wp,ten_parameters_exactly_matched=True),indent=2)+'\n')
    claude_note=[]
    if CLAUDE_A.is_dir():
        claude_case,claude_meta,claude_fit,_,_,_=read_case(CLAUDE_A)
        claude_result,_=load_case(claude_case)
        claude_estates=estate_diagnostics(claude_result)
        # Claude changed utility only; replace its gross-flow model observer with
        # independently recomputed net flow before comparing that one moment.
        recomputed=next(r['value'] for r in claude_estates if r['measure']=='annual_positive_death_estate_flow_to_living_gross_wealth' and r['definition']=='net_selling_cost')
        for row in claude_fit:
            if row['moment']=='bequest_wealth':
                row['model']=str(recomputed);row['gap']=str(recomputed-float(row['target']))
                row['loss_contribution']=str(float(row['weight'])*float(row['gap'])**2)
        revised,_=rescore_rows(claude_fit)
        own=[r for r in allfits if r['arm']=='single']
        comparisons=[dict(moment=a['moment'],target=a['target'],weight=a['weight'],role=a['role'],claude_model=a['model'],estate_a_single_model=b['model'],gap_model_difference=float(b['model'])-float(a['model'])) for a,b in zip(revised,own)]
        write_csv(args.out/'claude_a_comparison.csv',comparisons)
        write_csv(args.out/'claude_a_estate_diagnostics.csv',claude_estates)
        claude_note=['','## Independent Claude-A comparison','',
            'Claude-A changed utility only. Its gross bequest-flow report is replaced here by the net flow recomputed from its saved distribution and saving policies; other saved model moments are retained. Both are rescored with the same new wealth target.',
            f"Claude price: {claude_meta['price']}; Estate-A single price: {cases['single']['price']}.",
            '[All moment comparisons](claude_a_comparison.csv), [Claude estate and age-82 diagnostics](claude_a_estate_diagnostics.csv).','']
        claude_note+=markdown_table(['Moment','Target','Claude A','Estate A single','Model difference'],[[r[k] for k in ('moment','target','claude_model','estate_a_single_model','gap_model_difference')] for r in comparisons])
    lines=['# Estate-A fixed-parameter comparison' ,'',
        'The old single and old multiple baselines isolate Estate-A within each menu; source paths and actual saved estate flags are recorded in comparison.json. Both experimental arms use post-saving net estates W=bp+(1-psi)*P*h, with no extra interest multiplier and no recipient mapping. The single arm caps intended births at one; the multiple arm caps at three. The new wealth target is 4.45838713455674; all other empirical values and weights, including the bequest target, are retained. This is a fixed-parameter GE comparison, not recalibration.','',
        'Old native fit tables remain diagnostic and unchanged. The authoritative common new-target comparison is below. Living old-age wealth continues to use beginning b+P*h; death-estate diagnostics use saving bp and chosen housing. The retained SCF bequest target has an outstanding wealth-scope/recipient mapping mismatch and remains provisional. Descriptive count hazards use pre-birth exposure, unlike baseline post-birth descriptions.','']
    lines+=markdown_table(['Arm','Birth cap','Estate A','Price','Old diagnostic loss','Common new-target loss'],[(a,v['birth_cap'],v['estate_a'],v['price'],v['old_contract_loss'],v['new_contract_loss']) for a,v in cases.items()])
    for arm in cases:
        lines+=['',f'## {arm}: complete target fit','']
        lines+=markdown_table(['Moment','Target','Model','Gap','Weight','Loss contribution','Role'],[[r[k] for k in ('moment','target','model','gap','weight','loss_contribution','role')] for r in allfits if r['arm']==arm])
        lines+=['',f'## {arm}: all parameter restrictions','']
        lines+=markdown_table(['Parameter','Estimate','Lower','Upper','Near bound','Status'],[[r[k] for k in ('parameter','estimate','lower','upper','near_bound','status')] for r in allparams if r['arm']==arm])
    lines+=['','## Estate diagnostics','']+markdown_table(['Arm','Measure','Estate definition','Value'],[[r[k] for k in ('arm','measure','definition','value')] for r in estates])
    lines+=['','[Complete fit CSV](paired_target_fit_new_contract.csv), [complete parameters CSV](paired_parameters.csv), [estate diagnostics CSV](paired_estate_diagnostics.csv).','']
    lines+=claude_note
    (args.out/'RESULTS.md').write_text('\n'.join(lines))
    print(args.out/'RESULTS.md')

if __name__=='__main__':main()
