"""Cached local collateral envelope: age22 inherited renters only. No solve."""
import csv,json,argparse
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from extract_joint_budget import HERE,ROOT,REF,ALT,load,owner4,continuation,digest,get_fecundity_by_age

def main():
    pins=json.loads((ALT/'runtime_after.json').read_text())['imported_modules_after'];source=[]
    for name in ('production.engine.household','production.engine.parameters','production.engine.shared','production.engine.kernels'):
        pin=pins[name]
        if digest(ROOT/pin['path'])!=pin['sha256']:raise RuntimeError('Production source mismatch: '+name)
        source.append(dict(module=name,**pin))
    public0=json.loads((REF/'executed_P.json').read_text());public1=json.loads((ALT/'executed_P.json').read_text())
    changes=[k for k in public0 if not k.startswith('_') and public0[k]!=public1[k]]
    if changes!=['phi']:raise RuntimeError('Only phi may differ: '+str(changes))
    P,a=load(REF);T,t=load(ALT)
    if not np.array_equal(a['b_grid'],t['b_grid']) or not np.array_equal(a['p_eq'],t['p_eq']):raise RuntimeError('Paired grid/price mismatch')
    if list(P.H_own)!=[2,4,6,8,10]:raise RuntimeError('Unexpected housing menu')
    fec=get_fecundity_by_age(P);pi=fec[None,None,None,:,None]
    attempt=a['fert_probs'][...,1];wait=a['fert_probs'][...,0]
    fertile=(np.arange(P.J)>=P.A_f_start-1)&(np.arange(P.J)<P.A_f_end)
    G=a['g_beginning_distribution'][...,0,0].copy()
    den=1-pi*attempt
    if np.any(den[:,:,:,fertile,:]<=0):raise RuntimeError('Noninvertible first birth map')
    G[:,:,:,fertile,:]/=den[:,:,:,fertile,:]
    flow=(G*pi*attempt).sum(axis=(0,1,2,4))
    native=np.asarray(public0['_first_births_by_age'])
    flow_error=float(np.max(np.abs(flow[fertile]-native[fertile])))
    if flow_error>1e-10 or np.min(G)<-1e-12:raise RuntimeError('Prebirth first-flow reconstruction fails')
    W=G*pi*pi*attempt*wait/P.kappa_fert;W[:,:,:,~fertile,:]=0
    oldW=G*pi*pi*attempt*(1-attempt)/P.kappa_fert;oldW[:,:,:,~fertile,:]=0
    j=1;selectedW=W[:,0,0,j,:];selectedG=G[:,0,0,j,:]
    totalW=float(W.sum());sliceW=float(selectedW.sum())
    cuts=json.loads((HERE.parent/'phi095_decomposition/receipt.json').read_text())['income_cuts']
    earnings=np.asarray(P.income)[0,j]*np.asarray(P.z_grid)/(P.period_years*(1-P.tau_pay))
    rows=[];products={};maxbudget=0.;minsurplus=np.inf;maxnorm=0.
    for z in range(P.Nz):
        cv=continuation(P,a,z)
        for bi,b in enumerate(a['b_grid']):
            weight=float(selectedW[bi,z])
            if weight<=0:continue
            row=dict(b=float(b),income_state=z+1,annual_gross_earnings=float(earnings[z]),income_group=['low','mid','high'][int(np.searchsorted(cuts,earnings[z],side='left'))],exposure=float(selectedG[bi,z]),W=weight)
            for label,n,m in [('wait',0,0),('success',1,1)]:
                probs=np.asarray(a['tenure_probs'][bi,0,0,j,z,n,m],dtype=float)
                if not np.all(np.isfinite(probs)) or np.min(probs)<0 or np.max(probs)>1:raise RuntimeError('Invalid tenure probability')
                maxnorm=max(maxnorm,abs(float(probs.sum())-1))
                gain=0.;binding=0.;unsupported_q=0.;missing=[];qzero=[]
                for ten,house in enumerate(P.H_own,1):
                    q=float(probs[ten]);key=(label,int(house))
                    if key not in products:products[key]=dict(branch=label,rooms=int(house),nodes=0,supported_nodes=0,q_zero_nodes=0,q_zero_multiplier=None,unsupported_positive_q_nodes=0,W=0.,supported_W=0.,q_zero_W=0.,unsupported_positive_q_W=0.,qW=0.,supported_qW=0.,uncovered_qW=0.,selected_binding_qW=0.,floor_derivative_qW=0.,failure_reasons={})
                    r=products[key];r['nodes']+=1;r['W']+=weight;r['qW']+=weight*q
                    if q==0:
                        r['q_zero_nodes']+=1;r['q_zero_W']+=weight;qzero.append(str(int(house)))
                        continue # mu is unmeasured, not assigned zero; fixed-menu q is exactly zero.
                    try:
                        o=owner4(P,a,z,float(b),n,m,ten=ten,cv_all=cv)
                    except RuntimeError as exc:
                        reason=str(exc);missing.append(str(int(house))+':'+reason);unsupported_q+=q
                        r['unsupported_positive_q_nodes']+=1;r['unsupported_positive_q_W']+=weight;r['uncovered_qW']+=weight*q
                        r['failure_reasons'][reason]=r['failure_reasons'].get(reason,0)+1
                        continue
                    r['supported_nodes']+=1;r['supported_W']+=weight;r['supported_qW']+=weight*q
                    local=q*o['conditional_current_phi_derivative'];floor=q*o['floor_binding_mixture_weight']
                    gain+=local;binding+=floor;r['selected_binding_qW']+=weight*floor;r['floor_derivative_qW']+=weight*local
                    maxbudget=max(maxbudget,abs(o['mixture_budget_error']))
                    for endpoint in o['endpoint_controls']:
                        maxbudget=max(maxbudget,abs(endpoint['budget_error']));minsurplus=min(minsurplus,endpoint['c']-o['committed_consumption'])
                row.update({label+'_current_derivative_supported':gain,label+'_owner_q':float(probs[1:].sum()),label+'_selected_floor_q_supported':binding,label+'_uncovered_positive_q':unsupported_q,label+'_zero_q_products':','.join(qzero),label+'_unsupported_products':';'.join(missing)})
            row['positive_q_complete']=row['wait_uncovered_positive_q']==0 and row['success_uncovered_positive_q']==0
            row['all_products_measured']=row['positive_q_complete'] and not row['wait_zero_q_products'] and not row['success_zero_q_products']
            row['success_minus_wait_current_derivative']=row['success_current_derivative_supported']-row['wait_current_derivative_supported'] if row['positive_q_complete'] else None
            p0=np.asarray(a['fert_probs'][bi,0,0,j,z,:2]);p1=np.asarray(t['fert_probs'][bi,0,0,j,z,:2])
            valid=bool(np.all(p0>0)&np.all(p1>0)&np.all(np.isfinite(p0))&np.all(np.isfinite(p1)))
            row['permanent_gap_supported']=valid
            row['permanent_success_minus_wait_delta']=float(P.kappa_fert/fec[j]*(np.log(p1[1])-np.log(p1[0])-np.log(p0[1])+np.log(p0[0]))) if valid else None
            rows.append(row)
    def summarize(label,rr):
        full=sum(r['W'] for r in rr);strict=[r for r in rr if r['positive_q_complete']];strictW=sum(r['W'] for r in strict);common=[r for r in strict if r['permanent_gap_supported']];commonW=sum(r['W'] for r in common)
        result=dict(group=label,nodes=len(rr),baseline_W=full,positive_q_complete_W=strictW,positive_q_complete_share=strictW/full if full else None,all_products_measured_W=sum(r['W'] for r in rr if r['all_products_measured']),permanent_common_support_W=commonW,uncovered_permanent_common_W=strictW-commonW,current_local_first_birth_flow_derivative=sum(r['W']*r['success_minus_wait_current_derivative'] for r in strict),current_weighted_success_minus_wait_derivative=sum(r['W']*r['success_minus_wait_current_derivative'] for r in strict)/strictW if strictW else None,permanent_weighted_success_minus_wait_delta=sum(r['W']*r['permanent_success_minus_wait_delta'] for r in common)/commonW if commonW else None)
        for family in ('wait','success'):
            qW=sum(r['W']*r[family+'_owner_q'] for r in strict);floorW=sum(r['W']*r[family+'_selected_floor_q_supported'] for r in strict)
            result[family+'_selected_owner_probability']=qW/strictW if strictW else None
            result[family+'_selected_floor_probability']=floorW/strictW if strictW else None
            result[family+'_floor_share_given_owner']=floorW/qW if qW else None
            result[family+'_uncovered_qW']=sum(r['W']*r[family+'_uncovered_positive_q'] for r in rr)
            result[family+'_any_uncovered_positive_q_W']=sum(r['W'] for r in rr if r[family+'_uncovered_positive_q']>0)
            result[family+'_any_zero_q_product_W']=sum(r['W'] for r in rr if r[family+'_zero_q_products'])
        return result
    summaries=[summarize('all age22 inherited renters',rows)]+[summarize(g,[r for r in rows if r['income_group']==g]) for g in ('low','mid','high')]
    with (HERE/'age22_renter_envelope.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=dict(scope='Age22 inherited renters; all baseline W>0 b,z nodes; five owner products; baseline current collateral envelope holding future continuation and menu/screens fixed',formula='sum_h q_h Q_h mu_h SUCCESS minus WAIT; W=Gpre*pi^2*saved_attempt*saved_wait/kappa; rent direct derivative zero because its floor is unchanged',all_ages_all_beginning_tenure_first_birth_W=totalW,age22_renter_W=sliceW,age22_renter_share_of_all_first_birth_W=sliceW/totalW,strict_support_share_of_all_first_birth_W=summaries[0]['positive_q_complete_W']/totalW,income_cuts=cuts,income_definition='Fixed baseline pooled fertile childless-exposure annual gross labor earnings bins; current Markov state',summaries=summaries,product_coverage=list(products.values()),checks=dict(first_birth_native_flow_max_error=flow_error,max_budget_error=maxbudget,minimum_consumption_surplus=float(minsurplus),maximum_saved_tenure_probability_sum_error=maxnorm,old_a_one_minus_a_W= float(oldW.sum()),saved_both_probability_W_minus_old=float(W.sum()-oldW.sum())),source_pins=source,sources=[dict(path=str(path),arrays_sha256=digest(path/'solution_arrays.npz'),P_sha256=digest(path/'executed_P.json')) for path in (REF,ALT)],uncovered_convention='q=0 has unmeasured mu; fixed numerical menu contribution is zero, not evidence against menu-opening effects. Positive-q unsupported derivatives remain missing; complete-node aggregates do not impute them.',caveat='Local borrowing-capacity envelope is not a permanent-policy causal decomposition, welfare measure, or forward simulated birth-flow change; no new solve or optimization call.')
    result['extraction_driver_sha256']=digest(Path(__file__))
    result['conditional_control_driver_sha256']=digest(HERE/'extract_joint_budget.py')
    (HERE/'age22_renter_envelope.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(nodes=len(rows),summary=summaries[0],checks=result['checks']),indent=2))

def regenerate_standard_graphs(case=ALT,outdir=HERE/'standard_graph_replay',check_only=False):
    """Hydrate saved arrays/scalars, then reuse the unchanged canonical plotter."""
    from production.engine.diagnostics import write_diagnostics,_summary,_dumps_json
    case=Path(case);outdir=Path(outdir)
    pin=json.loads((ALT/'runtime_after.json').read_text())['imported_modules_after']['production.engine.diagnostics']
    if digest(ROOT/pin['path'])!=pin['sha256']:raise RuntimeError('Canonical diagnostics source mismatch')
    saved=(case/'standard_diagnostics/summary.json').read_text()
    values=json.loads(saved)
    sol=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in values.items()})
    with np.load(case/'solution_arrays.npz',allow_pickle=False) as z:
        for k in z.files:
            if not k.startswith('shared.'):setattr(sol,k,z[k].copy())
    p=json.loads((case/'executed_P.json').read_text());P=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in p.items()})
    if _dumps_json(_summary(sol,P))!=saved:raise RuntimeError('Cached diagnostic summary does not reconstruct exactly')
    if check_only:return dict(status='cached_diagnostic_summary_exact',case=str(case),expected_standard_graphs=17,plotter_source_sha256=pin['sha256'])
    if not outdir.resolve().is_relative_to(HERE.resolve()):raise RuntimeError('Graph replay output must remain in owned joint_budget folder')
    write_diagnostics(sol,P,outdir)
    count=len(list(outdir.glob('*.png')))
    if count!=17:raise RuntimeError('Canonical graph packet does not contain17 PNGs')
    return dict(status='canonical_standard_graphs_regenerated',output=str(outdir),standard_graphs=count)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--standard-graphs',choices=['check','write']);parser.add_argument('--case',type=Path,default=ALT);parser.add_argument('--output',type=Path,default=HERE/'standard_graph_replay');args=parser.parse_args()
    if args.standard_graphs:print(json.dumps(regenerate_standard_graphs(args.case,args.output,check_only=args.standard_graphs=='check'),indent=2))
    else:main()
