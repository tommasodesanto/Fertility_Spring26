"""Reviewed frozen native forward operators on saved policies; no solves."""
import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace
import numpy as np

PRE_HASH = 'caf760cb04c6a941cb6890510a832d5c0b61b2744f4fa4c16b89e1c8f3317903'
CASES = ('00_reference_p1.00','04_lifetime_repayment_only_p1.00')
ORDERS = ('first_births','second_births','third_bin_entries')
FIELDS = ('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs',
          'fert_probs','fert_value','fert2_probs','g_beginning_distribution','entry_by_loc','p_eq','b_grid')

def require(ok, message):
    if not ok:
        raise ValueError(message)

def forbidden(*args, **kwargs):
    raise RuntimeError('Household/GE solving is forbidden in saved-policy reduction')

def install_solve_traps(*modules):
    blocked=[]
    for module in modules:
        for name in dir(module):
            if name.startswith('solve') and callable(getattr(module,name)):
                setattr(module,name,forbidden)
                blocked.append(module.__name__+'.'+name)
    return blocked

def validate_flows(ev, saved):
    values=[float(ev.g_post_fertility[...,n:,:].sum()-ev.g_pre[...,n:,:].sum()) for n in range(1,4)]
    require(abs(sum(values)-float(ev.births)) < 2e-10, 'Birth-order sum')
    for name,value in zip(ORDERS,values):
        require(abs(value-float(saved[name])) < 2e-10, 'Retained native '+name+' differs')
    require(abs(float(ev.births)-float(saved['births'])) < 2e-10, 'Retained native total differs')
    return dict(zip(ORDERS,values),births=float(ev.births))

def grouped_flows(native_apply, pre, policy, P, grid):
    """Mask origin family only; native birth operator retains wealth and age."""
    result={}
    for n in range(3):
        for m in range(n+1):
            masked=np.zeros_like(pre)
            masked[...,n,m]=pre[...,n,m]
            post,births,_=native_apply(masked,policy.fert_probs,P,policy.fert2_probs)
            # Only origin n can have a birth: crossing n+1 measures that flow.
            state_flow=post[...,n+1:,:].sum(axis=(-2,-1))-masked[...,n+1:,:].sum(axis=(-2,-1))
            require(abs(float(state_flow.sum())-float(births)) < 2e-10,'Masked birth-order total')
            require(float(state_flow.min()) >= -2e-12,'Negative grouped birth flow')
            # Remaining dimensions: wealth, inherited tenure, location, age, income.
            for j in range(P.J):
                for negative in (False,True):
                    selected=(grid < 0) if negative else (grid >= 0)
                    flow=float(state_flow[selected,:,:,j,:].sum())
                    exposure=float(masked[selected,:,:,j,:,:,:].sum())
                    if flow != 0 or exposure != 0:
                        key=(j,n,m,negative)
                        result[key]=dict(flow=flow,exposure=exposure)
    return result

def three_corner(rows):
    output=[]
    for key in sorted(set().union(*(set(x) for x in rows))):
        values=[x.get(key,dict(flow=0.,exposure=0.)) for x in rows]
        j,n,m,negative=key
        a,b,c=[x['flow'] for x in values]
        output.append(dict(age=18+4*j,children_ever_born_before=n,children_at_home_before=m,
                           negative_wealth=negative,baseline_flow=a,credit_at_baseline_flow=b,credit_cohort_flow=c,
                           policy_change=b-a,composition_change=c-b,total_change=c-a,
                           baseline_exposure=values[0]['exposure'],credit_cohort_exposure=values[2]['exposure']))
    return output

def main(packet, root, work):
    started=time.monotonic()
    require(not work.exists(),'Refusing existing analysis output directory')
    work.mkdir(parents=True)
    sys.path.insert(0,str(packet))
    import fixed_price_responses as driver
    auth=driver.authenticate_candidate(work/'runtime_auth')
    cal=auth['context']['prepared'].rt['primitive'].pf.calendar
    transition=auth['context']['prepared'].rt['primitive'].pf.transition
    model=auth['context']['prepared'].rt['model']
    require(cal.apply_fertility is transition.apply_sequential_fertility,'Native sequential operator not installed')
    blocked=install_solve_traps(auth['solver'],model,cal)
    with np.load(root/'q0_reference_inherited_states.npz',allow_pickle=False) as arrays:
        retained_pre=arrays['g_pre']
    require(hashlib.sha256(retained_pre.tobytes()).hexdigest()==PRE_HASH,'Retained PRE identity')
    grid=auth['grid']; price=np.asarray([driver.BINDING['candidate_price']])
    require(auth['P'].age_start==18 and auth['P'].da==4 and auth['P'].J==17,'Frozen age clock')
    require(not getattr(auth['P'],'readiness_gate_enabled',False),'Origin childless-state convention')
    corners=[];checks=[];regimes=[];baseline_pre=None
    for index,case in enumerate(CASES):
        folder=root/case
        closure=json.loads((folder/'closure.json').read_text())
        receipt=json.loads((folder/'receipt.json').read_text())
        require(driver.sha(folder/'parameters.csv')==receipt['parameters_sha256'],'Saved parameter receipt')
        P=copy.deepcopy(auth['P'] if index==0 else auth['natural'])
        public_before=auth['base'].serialized({k:v for k,v in vars(P).items() if not k.startswith('_')})
        with np.load(folder/'solution_arrays.npz',allow_pickle=False) as arrays:
            data={name:arrays[name] for name in FIELDS}
            for name in ('bp_pol_stay','c_pol_stay'):
                if name in arrays.files:data[name]=arrays[name]
        require(np.array_equal(data.pop('b_grid'),grid) and np.array_equal(data['p_eq'],price),'Saved grid/price')
        sol=SimpleNamespace(**data)
        P._fert2_probs=sol.fert2_probs.copy()
        shared=auth['solver'].precompute_shared(P,grid)
        policy=cal.policy_from_solution(sol,price,P,grid,shared)
        pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
        require(reconstruction['stationary_post_fertility_nesting_l1'] <= 5e-9,'Retained POST nesting')
        require(reconstruction['stationary_feasibility_projection_mass']==0,'Reconstruction projection')
        if index==0:
            require(np.array_equal(pre,retained_pre) and hashlib.sha256(pre.tobytes()).hexdigest()==PRE_HASH,'Exact baseline PRE replay')
            baseline_pre=pre
        ev=cal.evaluate_period(price,pre,P,grid,shared,cal.SolveCounter(),supplied_policy=policy)
        require(ev.feasibility_projection_mass==0,'Evaluation projection')
        checks.append(dict(case=case,reconstruction=reconstruction,flows=validate_flows(ev,closure['cohort_summary'])))
        fertility=auth['context']['prepared'].rt['observe_initial_fertility'](ev,P,age_projection='uniform_birth_time')
        old=json.loads((folder/'observers.json').read_text())['fertility']['uniform_birth_time']['moments']['childless_rate_40_44']
        require(abs(fertility['moments']['childless_rate_40_44']-old)<2e-10,'Retained age40–44 childlessness')
        adjusted=transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),P)['topcode_adjusted_birth_children']
        require(abs(adjusted-closure['adjusted_births'])<2e-10,'Adjusted births')
        require(abs(adjusted/float(closure['actual_entry_rate'])-closure['cohort_summary']['completed_fertility'])<2e-10,'Saved completed fertility reconciliation')
        grouped=grouped_flows(cal.apply_fertility,pre,policy,P,grid)
        require(abs(sum(x['flow'] for x in grouped.values())-float(ev.births))<2e-10,'Grouped cohort flow sum')
        if index==0:
            corners.append(grouped)
        else:
            impact=cal.evaluate_period(price,baseline_pre,P,grid,shared,cal.SolveCounter(),supplied_policy=policy)
            require(impact.feasibility_projection_mass==0,'Impact projection')
            checks.append(dict(case=case+'_baseline_states',flows=validate_flows(impact,closure['baseline_state_impact'])))
            corners.append(grouped_flows(cal.apply_fertility,baseline_pre,policy,P,grid))
            corners.append(grouped)
        require(auth['base'].serialized({k:v for k,v in vars(P).items() if not k.startswith('_')})==public_before,
                'Forward operator changed public economic parameters')
        regimes.append(dict(case=case,pre_sha256=hashlib.sha256(pre.tobytes()).hexdigest(),
                            public_parameters_unchanged=True,
                            occupied_negative_wealth_mass=float(pre[grid<0].sum()),pre_mass=float(pre.sum()),
                            natural_support_certified=False,existing_support_status=receipt['support_status']))
        del ev,sol,data
    output=dict(status='native_saved_policy_forward_completed',household_solves=0,ge_solves=0,
                script_sha256=globals().get('EXECUTED_SOURCE_SHA256') or driver.sha(__file__),
                rows=three_corner(corners),checks=checks,regimes=regimes,blocked_solver_entries=blocked,
                baseline_pre_sha256=PRE_HASH,flow_units='raw birth events per normalized cross section per four-year period',
                order='policy first at baseline PRE, composition at credit policy; reverse corner not evaluated',
                value_recovery_exclusions='None: native birth operator uses all saved probabilities including endpoints; inherited support/projection gates remain strict',
                natural_credit_limitation='Original support-limited diagnostic; no full unoccupied support or grid convergence certificate',
                elapsed_seconds=time.monotonic()-started,maximum_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(output,indent=2,allow_nan=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    a=p.parse_args();main(a.packet,a.root,a.work)
