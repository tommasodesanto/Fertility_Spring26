#!/usr/bin/env python3
"""Three native fixed-price preference arms; no price root or recalibration."""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, os, signal, sys, time, traceback
from pathlib import Path
from types import SimpleNamespace
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
ECON=HERE.parent
MECHANISM=ECON/'utility_floor_psi_v1/mechanism_responses_v1'
OLD=ECON/'utility_floor_round2_v1'
BASE=ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner'
ARMS=('with_A','no_A','constant_alpha')
Q=.8198539089139171
OLD_TABLE=ECON/'utility_calibration_round1_v1/deployment/attempt2/comparison_requested/with_A_parameters.csv'

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=str)+'\n');tmp.replace(path)
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module
def rows(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def table(path,data):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)

class LifecycleBudget:
    def __init__(self,solver):self.solver=solver;self.used=0
    def __getattr__(self,name):return getattr(self.solver,name)
    def solve_markov_income_at_prices(self,*args,**kwargs):
        require(self.used<3,'Three native lifecycle-attempt cap exceeded');self.used+=1
        return self.solver.solve_markov_income_at_prices(*args,**kwargs)


def initialize(output):
    """Authenticate frozen reporting exactly once, validate every arm before solving."""
    import numpy as np
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    # Apply the same reviewed overlay before loading the authored private runtime.
    sys.path.insert(0,str(MECHANISM));load('runtime_overlay',HERE/'runtime_overlay.py')
    mechanism=load('share_A_fixed_price_mechanics',MECHANISM/'fixed_price_responses.py')
    sys.path.insert(0,str(OLD));sys.modules.pop('inputs',None)
    native=load('share_A_native_reporting',OLD/'runner.py')
    inputs=load('share_A_entry_inputs',ECON/'entry_calibration_pilot_v1/inputs.py')
    sys.path.insert(0,str(BASE))
    base=load('share_A_base_workflow',BASE/'run_comparison.py')
    ge=load('share_A_ge_validation',ECON/'entry_calibration_pilot_v1/phase_b_pilot.py')
    sys.path.insert(0,str(ROOT/'code/model'))
    sys.path.insert(0,str(ECON/'credit_no_taper_v1/small_credit_v1/source'))
    from small_credit_lab import credit
    from refactor_lab.engine import solver
    native.verify_sources()
    (output/'runtime').mkdir(parents=True,exist_ok=True)
    ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=output/'runtime'))
    native.install_reporter_on_authored(base.authored);base.authored.authenticate_frozen(ctx)
    helper=load('share_A_aggregate_helpers',ECON/'elasticity_v1/source_v2/run_credit.py')
    old_rows=rows(OLD_TABLE);expected0={r['parameter']:float(r['estimate']) for r in old_rows}
    require(len(old_rows)==31 and len(expected0)==31,'Old reference must have exactly 31 parameters')
    reference=json.loads((HERE/'reference/with_A_closure.json').read_text())
    require(float(reference['price'])==Q,'Old with-A price changed')
    point={k:expected0[k] for k in inputs.PARAMETERS};_,bounds,_=inputs.seed_and_bounds()
    require(expected0['delta_alpha_jump']==.1303383207216736,'Old birth-share loading changed')
    P0,grid=inputs.proposal();P0,entry=inputs.entry(P0,grid,'nonnegative_mean');P0=inputs.bind(P0,point,bounds)
    from small_credit_lab.engine import shared as es, child_preferences as ec, kernels as ek
    from refactor_lab.engine import shared as cs, child_preferences as cc, kernels as ck
    for current,executed in ((cs,es),(cc,ec),(ck,ek)):
        require(sha(current.__file__)==sha(executed.__file__),'Executed frozen primitive differs from native refactor: '+current.__name__)
    budget=LifecycleBudget(solver);auths={};arrays={};fingerprints={}
    for arm in ARMS:
        P=copy.deepcopy(P0);P.child_room_floor=False;P.hbar_first_child_jump=0.;P.hbar_child_rooms=0.;P.delta_alpha=0.
        P.delta_alpha_jump=0. if arm=='constant_alpha' else expected0['delta_alpha_jump']
        P.compensated_child_housing_shares=arm=='with_A'
        credit.bind_engine_credit(P,'corrected',0.)
        expected=dict(expected0);expected['delta_alpha_jump']=P.delta_alpha_jump
        context=dict(ctx);context.update(P=P,b_grid=grid,selected_d_bar=0.,reference_psi=float(P.psi_child),expected_dimensions={'wealth_grid_nodes':120,'income_states':9},expected_parameters=expected,free_coordinates=[],price_start=Q)
        actual=context['fp'].actual_parameters(context['prepared'],P,grid)
        ge.validate_parameter_estimates(context,old_rows,actual)
        require(len(actual)==31 and all(abs(float(actual[k])-v)<=2e-12*max(1.,abs(v)) for k,v in expected.items()),'Effective old full31 mismatch: '+arm)
        sd=solver.precompute_shared(P,grid)
        from small_credit_lab.engine import solver as executed
        live=executed.precompute_shared(P,grid)
        arm_arrays={k:np.asarray(v) for k,v in vars(sd).items() if isinstance(v,np.ndarray) and v.dtype!=object}
        for k in ('h_bar','c_bar','g_bar','alpha_flat','psi_v','escale_flat'):np.testing.assert_array_equal(getattr(sd,k),getattr(live,k))
        for k,value in arm_arrays.items():arrays[arm+'__'+k]=value
        fingerprints[arm]={k:dict(shape=list(v.shape),dtype=str(v.dtype),sha256=hashlib.sha256(v.tobytes()).hexdigest()) for k,v in arm_arrays.items()}
        param_rows=copy.deepcopy(old_rows)
        for row in param_rows:
            key=row['parameter'];row['reference_estimate']=row['estimate'];row['estimate']=actual[key]
            row['status']='Fixed old reference value; experimental preference arm'
            if key=='delta_alpha_jump' and arm=='constant_alpha':row['status']='Experimental restriction: constant consumption share';row['lower']=row['upper']=0.;row['near_bound']=True
        auths[arm]=dict(native=native,solver=budget,credit=credit,ge=ge,report_helpers=helper,context=context,P=P,natural=P,grid=grid,entry=entry,params_rows=param_rows,actual_parameters=actual)
    from refactor_lab.engine.child_preferences import apply_child_preferences
    np.testing.assert_array_equal(arrays['with_A__alpha_flat'],arrays['no_A__alpha_flat'])
    require(np.all(arrays['constant_alpha__alpha_flat']==float(auths['constant_alpha']['P'].alpha_cons)), 'Constant alpha arm is not constant')
    # Verify the compensator by invoking the existing native preference switch on
    # a copy of the uncompensated material multiplier, with no model solve.
    alpha=arrays['no_A__alpha_flat'].copy()
    material=arrays['no_A__escale_flat'].copy()
    benefit=arrays['no_A__psi_v'].copy()
    probe=copy.deepcopy(auths['with_A']['P']);probe.child_benefit_curvature=0.
    apply_child_preferences(probe,alpha,benefit,material)
    np.testing.assert_array_equal(material,arrays['with_A__escale_flat'])
    for arm in ARMS:
        require(np.all(arrays[arm+'__h_bar']==0.),'Nonzero housing floor')
        np.testing.assert_array_equal(arrays[arm+'__psi_v'],arrays['with_A__psi_v'])
    np.testing.assert_array_equal(arrays['no_A__escale_flat'],arrays['constant_alpha__escale_flat'])
    np.savez_compressed(output/'native_precomputed_arrays.npz',**arrays)
    # Native arrays, rather than a second implementation of utility, are the evidence.
    comparison=[]
    for name in ('alpha_flat','escale_flat','psi_v','h_bar'):
        values={arm:arrays[arm+'__'+name].reshape(-1).tolist() for arm in ARMS}
        comparison.append(dict(array=name,values=values))
    write(output/'utility_precompute_comparison.json',dict(arrays=comparison,fingerprints=fingerprints,entry=entry,old_table_sha256=sha(OLD_TABLE),price=Q,lifecycle_solves=budget.used))
    mechanism.BINDING=dict(mechanism.BINDING,candidate_price=Q,candidate_loss=18.128926)
    mechanism.CONTRACT=dict(mechanism.CONTRACT,candidate_case='old_with_A_fixed_q_decomposition',candidate_label='old with-A parameter reference; fixed-q utility decomposition')
    require(budget.used==0,'Initializer unexpectedly solved a lifecycle')
    write(output/'runtime_parameter_fingerprints.json',{arm:dict(compensated_child_housing_shares=bool(auth['P'].compensated_child_housing_shares),actual_parameters=auth['actual_parameters'],fixed_entry_sha256=entry['conditional_sha256'],rental_menu_cap=float(auth['P'].hR_max),owner_room_menu=np.asarray(auth['P'].H_own).tolist(),parameters_csv_sha256=sha(OLD_TABLE)) for arm,auth in auths.items()})
    write(output/'initializer_receipt.json',dict(status='passed_zero_lifecycle',arms=list(ARMS),full_parameter_rows=31,lifecycle_solves=0,price=Q,entry=entry))
    return dict(auths=auths,mechanism=mechanism,budget=budget)


def evaluate(state,arm,out,deadline):
    """Reuse the unchanged fixed-price household, accounting, PAYGO and observer gates."""
    auth=state['auths'][arm];mechanism=state['mechanism']
    audit=auth['context']['prepared'].rt['audit'];original=audit.standard_diagnostics
    def diagnostics(packet,*args,**kwargs):
        result=original(packet,*args,**kwargs)
        import numpy as np
        g=np.asarray(packet['evaluation'].g_current);h=np.asarray(packet['solution'].hR_pol)
        P=packet['parameters'];require(g.shape==h.shape and g.ndim==7,'Supplemental native housing dimensions changed')
        profile=[]
        for j in range(g.shape[3]):
            for m in range(int(P.n_child_states)):
                mask=np.asarray([[min(n,cs,int(P.n_child_states)-1)==m for cs in range(int(P.n_child_states))] for n in range(int(P.n_parity))])
                renter=g[:,0,:,j]*mask;rm=float(renter.sum());rrooms=float((renter*np.where(renter>0,h[:,0,:,j],0.)).sum())
                owner_mass=0.;owner_rooms=0.
                for t,size in enumerate(np.asarray(P.H_own).reshape(-1),1):
                    mass=float((g[:,t,:,j]*mask).sum());owner_mass+=mass;owner_rooms+=mass*float(size)
                total=rm+owner_mass
                profile.append(dict(age_start=float(P.age_start)+j*float(P.da),children_at_home=m,household_mass=total,ownership=owner_mass/total if total else None,mean_rooms=(rrooms+owner_rooms)/total if total else None,renter_mean_rooms=rrooms/rm if rm else None,owner_mean_rooms=owner_rooms/owner_mass if owner_mass else None,mapped_unit_rent=Q*float(P.user_cost_rate),mean_renter_spending=Q*float(P.user_cost_rate)*rrooms/rm if rm else None))
        table(out/'supplemental_parent_housing_rent_profiles.csv',profile)
        return result
    audit.standard_diagnostics=diagnostics
    try:result=mechanism.make_price_cell(auth,'reference',1.,out,deadline)
    finally:audit.standard_diagnostics=original
    result.update(arm=arm,economic_changes=[] if arm=='with_A' else ['Compensated housing-share utility factor disabled']+(['Birth consumption-share jump set to zero'] if arm=='constant_alpha' else []),candidate_base_loss=18.128926,baseline_state_impact_scope='Own-arm reconstructed stationary pre-fertility state; not a common with-A inherited distribution')
    write(out/'closure.json',result)
    auth['native'].table(out/'parameters.csv',auth['params_rows'])
    restrictions=[dict(parameter=r['parameter'],value=r['estimate'],restriction='held fixed throughout fixed-q decomposition',compensation_enabled=arm=='with_A') for r in auth['params_rows']]
    table(out/'parameter_restrictions.csv',restrictions)
    # Saved native housing policies and distributions support supplementary views.
    import numpy as np
    with np.load(out/'solution_arrays.npz',allow_pickle=False) as saved:
        profiles={key:saved[key] for key in saved.files if any(token in key.lower() for token in ('housing','h_r','h_o','hr','ho','rent','g_beginning','p_owner','fert2'))}
    profiles['wealth_grid']=auth['grid'];profiles['income_grid']=np.asarray(auth['P'].z_grid);profiles['mapped_rent']=np.asarray([Q*float(auth['P'].user_cost_rate)])
    np.savez_compressed(out/'supplemental_parent_housing_rent_profiles.npz',**profiles)
    write(out/'supplemental_profile_index.json',dict(status='native_solved_arrays_no_additional_lifecycle',arrays={k:list(v.shape) for k,v in profiles.items()},mapped_rent=Q*float(auth['P'].user_cost_rate),interpretation='Conditional native policies; use realized tenure and beginning distribution for occupied parent aggregates.'))
    return result


def write_incremental_comparison(output,records):
    fit=[];parameters=[];compact=[]
    for record in records:
        if record['status']!='completed':continue
        arm=record['arm'];case_rows=rows(output/arm/'target_fit.csv')
        fit.extend(dict(arm=arm,**r) for r in case_rows)
        parameters.extend(dict(arm=arm,**r) for r in rows(output/arm/'parameters.csv'))
        compact.extend(dict(arm=arm,Moment=r['moment'],Target=r['target'],Model=r['model']) for r in case_rows)
    if fit:table(output/'comparison_target_fit.csv',fit)
    if parameters:table(output/'comparison_parameters.csv',parameters)
    if compact:table(output/'comparison_moment_target_model.csv',compact)
    write(output/'comparison_summary.json',dict(status='experimental_fixed_price_not_adopted',cases=[dict(arm=r['arm'],cohort_summary=r['closure'].get('cohort_summary')) for r in records if r['status']=='completed']))

def run(output,deadline_epoch,*,state=None,evaluator=None,mock=False):
    started=time.time();deadline=min(float(deadline_epoch),started+1200.);output=Path(output);output.mkdir(parents=True,exist_ok=False)
    require(deadline>started,'Deadline has expired')
    state=initialize(output/'initialization') if state is None else state;evaluator=evaluate if evaluator is None else evaluator
    records=[];attempts=0
    write(output/'launch.json',dict(started_epoch=started,deadline_epoch=deadline,global_budget_seconds=1200,case_budget_seconds=600,lifecycle_cap=3,cpu_threads=1,arms=list(ARMS),price=Q,mocked=mock))
    for arm in ARMS:
        if time.time()>=deadline:
            records.append(dict(arm=arm,status='not_started_global_deadline'));break
        out=output/arm;out.mkdir();case_deadline=min(deadline,time.time()+600.)
        attempts+=1;require(attempts<=3,'Three-case controller cap exceeded')
        write(out/'latest.json',dict(status='running',arm=arm,deadline_epoch=case_deadline))
        old=signal.getsignal(signal.SIGALRM)
        def alarm(*_):raise TimeoutError('Declared fixed-q case deadline')
        signal.signal(signal.SIGALRM,alarm);timer=signal.setitimer(signal.ITIMER_REAL,max(.001,case_deadline-time.time()))
        try:
            result=evaluator(state,arm,out,case_deadline)
            record=dict(arm=arm,status='completed',closure=result,closure_path=str(out/'closure.json'))
        except BaseException as exc:
            record=dict(arm=arm,status='failed',error=repr(exc),traceback=traceback.format_exc());records.append(record)
            write(out/'failure.json',record);write(output/'completed.json',dict(status='failed_no_retry',cases=records,case_attempts=attempts));raise
        finally:
            signal.setitimer(signal.ITIMER_REAL,*timer);signal.signal(signal.SIGALRM,old)
        records.append(record)
        if not mock:write_incremental_comparison(output,records)
        write(output/'latest_completed.json',dict(cases=records,case_attempts=attempts,lifecycle_solves=0 if mock else state['budget'].used))
    summary=[];fit=[];parameters=[]
    for record in records:
        if record['status']!='completed':continue
        arm=record['arm'];closure=record['closure'];summary.append(dict(arm=arm,price=Q,cohort_summary=closure.get('cohort_summary'),renewal_residual=closure.get('renewal_residual_reported_not_imposed'),market_residual=closure.get('relative_market_residual_reported_not_imposed')))
        if not mock:
            fit.extend(dict(arm=arm,**r) for r in rows(output/arm/'target_fit.csv'))
            parameters.extend(dict(arm=arm,**r) for r in rows(output/arm/'parameters.csv'))
    if fit:table(output/'comparison_target_fit.csv',fit)
    if parameters:table(output/'comparison_parameters.csv',parameters)
    write(output/'comparison_summary.json',dict(status='experimental_fixed_price_not_adopted',cases=summary))
    write(output/'completed.json',dict(status='completed' if len(records)==3 and all(r['status']=='completed' for r in records) else 'incomplete',cases=records,case_attempts=attempts,lifecycle_solves=0 if mock else state['budget'].used,production_adoption=False))
    return records


def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=('initialize-only','smoke-loop','run'));p.add_argument('--initialize-only',action='store_true');p.add_argument('--smoke-loop',action='store_true');p.add_argument('--output','--out',dest='output',type=Path,required=True);p.add_argument('--deadline-epoch',type=float);p.add_argument('--go-reviewed',action='store_true');args=p.parse_args()
    args.mode=args.mode or ('initialize-only' if args.initialize_only else 'smoke-loop' if args.smoke_loop else 'run')
    if args.mode=='initialize-only':initialize(args.output);return
    if args.mode=='smoke-loop':
        seen=[]
        def mocked(state,arm,out,deadline):
            seen.append(arm);result=dict(lifecycle_seconds=0.,cohort_summary={'completed_fertility':2.1},renewal_residual_reported_not_imposed=0.,relative_market_residual_reported_not_imposed=.1);write(out/'closure.json',result);return result
        run(args.output,time.time()+60,state={},evaluator=mocked,mock=True);require(seen==list(ARMS),'Mock actual three-case controller drift');return
    require(sys.platform.startswith('linux') and bool(os.environ.get('SLURM_JOB_ID')),'Numerical execution is restricted to Torch Slurm');require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Numerical execution requires one allocated CPU');require(args.go_reviewed,'Numerical execution requires --go-reviewed');require(args.deadline_epoch is not None,'Numerical execution requires explicit deadline');run(args.output,args.deadline_epoch)
if __name__=='__main__':main()
