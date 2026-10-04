"""Age-22 cached owner-four-room controls and current-floor shadow values.

No solve and no optimization call. Reconstructs continuation using the same
income expectation and child-aging function used by the saved Bellman solve.
Numerical transaction interpolation is shown as its two-node mixture.
"""
from pathlib import Path
from types import SimpleNamespace
import json,sys,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/'code/model'))
from production.engine.household import apply_child_aging, validate_native_solvency_mode
from production.engine.parameters import parent_age_maturation_active, get_fecundity_by_age
from production.engine.shared import income_transition_values, get_phi_choice_tensor, DEAD_VALUE_CUTOFF
REF=ROOT/'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/phi_080'
ALT=HERE.parent/'phi_095_run1'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def load(path):
    d=json.loads((path/'executed_P.json').read_text());P=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in d.items()})
    names=['V','c_pol','bp_pol','tenure_probs','fert_probs','fert_value','g_beginning_distribution','b_grid','p_eq','income_transition','shared.cb_flat','shared.hb_flat','shared.alpha_flat','shared.escale_flat','shared.psi_flat']
    with np.load(path/'solution_arrays.npz',allow_pickle=False) as z:a={k:z[k].copy() for k in names}
    if parent_age_maturation_active(P):raise RuntimeError('Successful-birth controls require newborn-exempt capture, absent from saved arrays')
    if P.survival_probs[1]!=1 or P.transfer_floor_G0!=0 or P.transfer_floor_Gn!=0 or P.owner_size_cost!=0:
        raise RuntimeError('This illustration requires the exact retained age22/no-transfer/no-size-cost contract')
    if validate_native_solvency_mode(P) or not P.native_purchase_income or not P.native_exact_allocation_output or not P.exhaustive_saving_control:
        raise RuntimeError('Unsupported continuation or savings runtime contract')
    if float(P.H_own[1])!=4 or P.I!=1 or P.chi<=0 or P.readiness_gate_enabled:
        raise RuntimeError('Unsupported illustration dimensions or preferences')
    raw_transition=np.array(P.Pi_z,copy=True)
    _,_,executed_transition=income_transition_values(P)
    if not np.array_equal(executed_transition,a['income_transition']):raise RuntimeError('Canonical executed income transition identity mismatch')
    a['transition_normalization_max_change']=float(np.max(np.abs(raw_transition-executed_transition)))
    return P,a

def continuation(P,a,z):
    v=np.zeros((P.Nb,1+P.n_house,P.I,P.n_parity,P.n_child_states))
    for zz in range(P.Nz):
        if a['income_transition'][z,zz]>0:v+=a['income_transition'][z,zz]*a['V'][:,:,:,2,zz,:,:]
    return apply_child_aging(v,P,P.Nb,1+P.n_house,P.I,P.n_parity,P.n_child_states,age_index=1)

def owner4(P,a,z,b,n,m,ten=2,cv_all=None):
    house=float(P.H_own[ten-1]);grid=a['b_grid'];price=float(a['p_eq'][0]);Q=house*price;K=(P.delta+P.tau_H)*Q;y=float(P.income[0,1]*P.z_grid[z]);x=b-Q/P.R_gross
    lo=int(np.searchsorted(grid,x,side='right')-1);lo=max(0,min(lo,len(grid)-2));weight=(x-grid[lo])/(grid[lo+1]-grid[lo])
    if not 0<=weight<=1:raise RuntimeError('Transaction leaves saved support')
    cv=(continuation(P,a,z) if cv_all is None else cv_all)[:,ten,0,n,m]
    state=n+P.n_parity*m;cb=float(a['shared.cb_flat'][0,state]);hb=float(a['shared.hb_flat'][0,state]);alpha=float(a['shared.alpha_flat'][0,state]);es=float(a['shared.escale_flat'][0,state])
    hsurplus=house-P.owner_h_bar_scale*hb
    if hsurplus<=0:raise RuntimeError('Infeasible owner housing surplus')
    Ko=(P.chi*hsurplus)**((1-alpha)*(1-P.sigma));floor=-float(get_phi_choice_tensor(P)[0,ten,n,m])*Q
    endpoint=[]
    for index,w in [(lo,1-weight),(lo+1,weight)]:
        if w==0:continue
        c=float(a['c_pol'][index,ten,0,1,z,n,m]);bp=float(a['bp_pol'][index,ten,0,1,z,n,m]);slack=bp-floor;budget=P.R_gross*grid[index]+y-K-bp-c
        if abs(budget)>1e-9 or c-cb<=1e-10 or slack < -1e-9:raise RuntimeError('Unsupported conditional control/budget')
        ix=int(np.searchsorted(grid,bp,side='right')-1);ix=max(0,min(ix,len(grid)-2))
        slope=(cv[ix+1]-cv[ix])/(grid[ix+1]-grid[ix]);uc=es*Ko*alpha*(c-cb)**(alpha*(1-P.sigma)-1)
        value=es*Ko*(c-cb)**(alpha*(1-P.sigma))/(1-P.sigma)+float(a['shared.psi_flat'][0,state])+P.beta*np.interp(bp,grid,cv)
        if not np.isfinite(value) or value<=DEAD_VALUE_CUTOFF:raise RuntimeError('Dead conditional owner value')
        binds=abs(slack)<=1e-9;mu=uc-P.beta*slope if binds else 0.
        if binds and mu < -1e-8:raise RuntimeError('Negative floor multiplier beyond roundoff')
        at_knot=bool(abs(bp-grid[ix])<=1e-10 and ix>0)
        left_slope=(cv[ix]-cv[ix-1])/(grid[ix]-grid[ix-1]) if at_knot else slope
        if not binds and not (P.beta*slope-1e-8<=uc<=P.beta*left_slope+1e-8):
            raise RuntimeError('Interior saving first-order/kink inequalities fail')
        endpoint.append(dict(index=index,weight=float(w),c=c,bp=bp,floor_slack=slack,budget_error=float(budget),floor_binding=binds,uc=float(uc),beta_continuation_right_slope=float(P.beta*slope),beta_continuation_left_slope=float(P.beta*left_slope),saving_grid_knot=at_knot,mu=float(mu)))
    mix=lambda key:sum(e['weight']*e[key] for e in endpoint)
    bi=int(np.where(grid==b)[0][0]);prob=float(a['tenure_probs'][bi,0,0,1,z,n,m,ten])
    return dict(phi=float(P.phi[n]),income_state=z+1,beginning_b=b,family='WAIT' if n==0 else 'SUCCESS',four_year_aftertax_income=y,purchase_cost=Q,owner_carrying_cost=K,nominal_downpayment=(1-float(P.phi[n]))*Q,purchase_screen_slack=P.R_gross*b+y-(1-float(P.phi[n]))*Q,owner4_choice_probability=prob,ending_asset_floor=floor,committed_consumption=cb,mixed_consumption_surplus=mix('c')-cb,consumption_ceiling_at_floor=P.R_gross*b+y-Q-K-floor,mixed_consumption=mix('c'),mixed_ending_assets=mix('bp'),mixed_floor_slack=mix('floor_slack'),floor_binding_mixture_weight=sum(e['weight'] for e in endpoint if e['floor_binding']),conditional_current_phi_derivative=Q*mix('mu'),owner4_inclusive_phi_contribution=prob*Q*mix('mu'),mixture_budget_error=P.R_gross*b+y-Q-K-mix('c')-mix('bp'),endpoint_controls=endpoint)

def main():
    pins=json.loads((ALT/'runtime_after.json').read_text())['imported_modules_after']
    source_pins=[]
    for name in ['production.engine.household','production.engine.parameters','production.engine.shared','production.engine.kernels']:
        pin=pins[name];actual=digest(ROOT/pin['path'])
        if actual!=pin['sha256']:raise RuntimeError('Executed production source changed: '+name)
        source_pins.append(dict(module=name,**pin))
    P,a=load(REF);Q,c=load(ALT);rows=[];gains=[]
    for z in (3,4,5):
        for b in (0.,1.):
            for label,R,d in [('baseline',P,a),('relaxed',Q,c)]:
                for n,m in [(0,0),(1,1)]:rows.append(dict(case=label,**owner4(R,d,z,b,n,m)))
            bi=int(np.where(a['b_grid']==b)[0][0]);pi=float(get_fecundity_by_age(P)[1]);k=P.kappa_fert
            def branch(d):
                probs=d['fert_probs'][bi,0,0,1,z,:2];I=d['fert_value'][bi,0,0,1,z]
                wait=I+k*np.log(probs[0]);D=k*(np.log(probs[1])-np.log(probs[0]))/pi
                return wait,D
            w0,D0=branch(a);w1,D1=branch(c)
            gains.append(dict(income_state=z+1,beginning_b=b,wait_credit_gain=float(w1-w0),success_credit_gain=float(w1-w0+D1-D0),success_minus_wait_credit_gain=float(D1-D0)))
    checks=dict(illustrations=len(rows),endpoints=sum(len(r['endpoint_controls']) for r in rows),max_endpoint_budget_error=max(abs(e['budget_error']) for r in rows for e in r['endpoint_controls']),max_mixture_budget_error=max(abs(r['mixture_budget_error']) for r in rows),minimum_consumption_surplus=min(e['c']-r['committed_consumption'] for r in rows for e in r['endpoint_controls']),saving_optimality='Floor and interior/kink first-order inequalities pass at all endpoints; no reoptimization performed')
    result=dict(contract='age22 inherited renter, common b/income; owner4 fixed-house current-floor shadow illustration, not aggregate causal decomposition',checks=checks,source_pins=source_pins,transition_normalization_max_change=[a['transition_normalization_max_change'],c['transition_normalization_max_change']],newborn_exemption_active=False,newborn_exemption_reason='Executed P lacks child_maturation_mode; exact runtime helper defaults constant; VI_ex is not called',rows=rows,full_permanent_credit_branch_gains=gains,sources=[dict(path=str(p),arrays_sha256=digest(p/'solution_arrays.npz'),P_sha256=digest(p/'executed_P.json')) for p in (REF,ALT)],caveat='Q*mu holds continuation and current menu fixed locally; total permanent credit changes future continuation, other housing choices and screens. Interpolated controls are the implemented two-grid-node mixture; not exact continuous-state reoptimization.')
    (HERE/'joint_budget.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'example':[r for r in rows if r['income_state']==5 and r['beginning_b']==0],'gains':[g for g in gains if g['income_state']==5 and g['beginning_b']==0]},indent=2))

if __name__=='__main__':main()
