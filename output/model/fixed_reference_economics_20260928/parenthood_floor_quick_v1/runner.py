"""Isolated fixed-parameter parenthood-floor diagnostic. No search or adoption."""
from __future__ import annotations
import argparse, copy, importlib.util, json, os, signal, sys, time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PILOT=ROOT/'output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1'
sys.path.insert(0,str(PILOT))
import inputs
spec=importlib.util.spec_from_file_location('parenthood_pilot_helpers',PILOT/'runner.py')
pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)
PLAN=json.loads((HERE/'plan.json').read_text())
require=inputs.require
write=pilot.write
BASE_BIND=inputs.bind
ARM="floor"

def floor_bind(P,point,bounds):
    Q=BASE_BIND(P,point,bounds)
    Q.child_room_floor=ARM=="floor"
    Q.hbar_first_child_jump=float(PLAN['h_P']) if ARM=='floor' else 0.
    Q.hbar_child_rooms=0.
    Q.delta_alpha_jump=float(point["delta_alpha_jump"]) if ARM=="no_A" else 0.
    Q.delta_alpha=0.
    Q.compensated_child_housing_shares=False
    require(Q.owner_h_bar_scale==1.,'Owner physical-room floor convention drift')
    return Q

def expected_parameters(point):
    x={r['parameter']:float(r['estimate']) for r in PLAN['control_parameters']}
    x.update(point,delta_alpha_jump=float(point['delta_alpha_jump']) if ARM=='no_A' else 0.,delta_alpha=0.,h_P=float(PLAN['h_P']) if ARM=='floor' else 0.)
    x['child_benefit_CRRA_coefficient']=(1-point['child_benefit_curvature'])*x['psi_child']
    return x

def verify():
    for rel,digest in json.loads((HERE/'source_pins.json').read_text()).items():
        require(inputs.sha(ROOT/rel)==digest,'Diagnostic source drift: '+rel)
    receipt=json.loads((HERE/'authenticated_control.json').read_text())
    for name,digest in receipt['compact_control_hashes'].items():
        require(inputs.sha(HERE/'reference_control'/name)==digest,'Control table drift '+name)
    require(receipt['exact_two_independent_full_ge_repeats'] and len(receipt['standard_plot_hashes'])==17,'Control not authenticated')
    return receipt

def utility_checks(P,grid,out):
    # Shared arrays and kernel configuration are the exact arrays used by the native solve.
    sys.path.insert(0,str(ROOT/'code/model'))
    from refactor_lab.engine import solver
    Q=floor_bind(P,PLAN['seed'],PLAN['bounds']);sd=solver.precompute_shared(Q,grid)
    require(Q.preference_spec=='eqscale' and Q.eqscale_form=='power' and Q.child_state_mode=='independent_count','Wrong utility architecture')
    require(not Q.compensated_child_housing_shares and Q.delta_alpha==Q.hbar_child_rooms==0.,'Compensation/later-child loading retained')
    for n in range(Q.n_parity):
        for m in range(min(n+1,Q.n_child_states)):
            expected_floor=PLAN['h_P'] if m>0 and ARM=='floor' else 0.
            require(sd.h_bar[n,m]==expected_floor,'Floor is not first-child-only')
            require(sd.c_bar[n,m]==sd.g_bar[n,m]==0.,'Nonhousing floor introduced')
            alpha=float(sd.alpha_flat.reshape((Q.n_parity,Q.n_child_states),order='F')[n,m])
            expected_alpha=Q.alpha_cons-Q.delta_alpha_jump if m>0 and ARM=='no_A' else Q.alpha_cons
            require(alpha==expected_alpha,'Child share differs')
            benefit=Q.psi_child*m**(1-Q.child_benefit_curvature) if m else 0.
            require(sd.psi_v[n,m]==benefit,'Nonlinear benefit differs')
            e=((2+.7*m)/2)**.7
            require(abs(sd.escale_flat.reshape((Q.n_parity,Q.n_child_states),order="F")[n,m]-e**(Q.sigma-1))<2e-15,'Equivalence scale changed')
            if m:
                c,h=1.7,4.2
                native=sd.escale_flat.reshape((Q.n_parity,Q.n_child_states),order="F")[n,m]*(c**alpha*(h-sd.h_bar[n,m])**(1-alpha))**(1-Q.sigma)/(1-Q.sigma)+sd.psi_v[n,m]
                requested=((c**alpha*(h-expected_floor)**(1-alpha))/e)**(1-Q.sigma)/(1-Q.sigma)+benefit
                require(abs(native-requested)<2e-15,'CRRA floor formula differs')
    require(Q.hR_max>PLAN['h_P'],'Rental cap makes parent housing infeasible')
    # Owner kernel rejects physical h<=floor before multiplying residual rooms by chi.
    from refactor_lab.engine import household,kernels
    import inspect
    kernel=inspect.getsource(kernels.full_owner_block_kernel.py_func)
    require('if strict_hbar_feasibility and ht_c <= 0.0:' in kernel and 'ht_c = owner_service_premium * ht_c' in kernel,'Native strict owner-floor gate missing')
    household_src=inspect.getsource(household)
    require('strict_owner_hbar_feasibility = int(' in household_src,'Household floor gate binding missing')
    write(out/'utility_verification.json',dict(status='passed_zero_solve',arm=ARM,floor=Q.hbar_first_child_jump,compensation=False,child_share_loading=Q.delta_alpha_jump,nonlinear_child_benefit=True,first_child_only_floor=True,owner_services='chi*(physical_rooms-h_P)',renter_services='physical_rooms-h_P',strict_owner_feasibility=True,lifecycle_solves=0))
    return Q

def install_reporter_adapter():
    # Patch only this process's authenticated observer, leaving its source untouched.
    sys.path.insert(0,str(pilot.BASE))
    spec=importlib.util.spec_from_file_location('parenthood_native_workflow',pilot.BASE/'run_comparison.py')
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    authored=base.authored;original=authored.authenticate_frozen
    def authenticate(ctx):
        result=original(ctx);fp=ctx['fp']
        if not getattr(fp,'_parenthood_floor_reporting_adapter',False):
            native=fp.actual_parameters
            def actual(prepared,P,grid):
                values=native(prepared,P,grid)
                require(float(P.hbar_child_rooms)==0.,'Unexpected later-child floor')
                values['h_P']=float(P.hbar_first_child_jump) if P.child_room_floor else 0.
                return values
            fp.actual_parameters=actual;fp._parenthood_floor_reporting_adapter=True
        return result
    authored.authenticate_frozen=authenticate
    return base

def validate_native_reporter(Q,grid,out):
    from types import SimpleNamespace
    base=install_reporter_adapter()
    auth_out=out/'native_reporter_validation';auth_out.mkdir()
    ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=auth_out))
    base.authored.authenticate_frozen(ctx)
    from small_credit_lab import credit
    from small_credit_lab.engine import solver
    credit.bind_engine_credit(Q,'corrected',0.)
    live_sd=solver.precompute_shared(Q,grid)
    # Native execution and component check must use byte-identical material preference sources.
    import inspect
    from small_credit_lab.engine import shared as executed_shared,child_preferences as executed_child,kernels as executed_kernels
    from refactor_lab.engine import shared as checked_shared,child_preferences as checked_child,kernels as checked_kernels
    for checked,executed in [(checked_shared,executed_shared),(checked_child,executed_child),(checked_kernels,executed_kernels)]:
        require(inputs.sha(Path(checked.__file__))==inputs.sha(Path(executed.__file__)),'Executed preference/kernel source differs')
    expected=expected_parameters(PLAN['seed'])
    actual=ctx['fp'].actual_parameters(ctx['prepared'],Q,grid)
    import phase_b_pilot as ge
    ge.validate_parameter_estimates(dict(expected_parameters=expected),PLAN['control_parameters'],actual)
    floor=PLAN['h_P'] if ARM=='floor' else 0.
    require(float(live_sd.h_bar[1,1])==floor and float(live_sd.h_bar[1,0])==0.,'Executed shared floor differs')
    write(out/'native_reporter_verification.json',dict(status='passed_zero_solve',actual_parameters=actual,expected_h_P=floor,source_identity_authenticated=True,lifecycle_solves=0))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['preflight','run'],required=True)
    parser.add_argument('--arm',choices=PLAN['arms'],required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--deadline-seconds',type=float,default=1200)
    global ARM
    args=parser.parse_args();ARM=args.arm;args.out.mkdir(parents=True,exist_ok=False)
    deadline=time.time()+args.deadline_seconds
    signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('Diagnostic global deadline')))
    signal.setitimer(signal.ITIMER_REAL,args.deadline_seconds)
    try:
        control=verify();P,grid=inputs.proposal();P,entry=inputs.entry(P,grid,'nonnegative_mean')
        Q=utility_checks(P,grid,args.out)
        if args.mode=="preflight":
            validate_native_reporter(Q,grid,args.out)
        else:
            # Preflight authenticates in a separate interpreter. Native run must
            # authenticate once through native_evaluator; its candidate validator
            # checks all 31 effective parameters before the first lifecycle call.
            install_reporter_adapter()
        write(args.out/'input_contract.json',dict(plan=PLAN,entry=entry,control=control,closure='price clears actual birth renewal; population clears absolute housing supply',refit=False,experimental_not_adopted=True))
        inputs.bind=floor_bind;pilot.expected_parameters=expected_parameters
        if args.mode=='preflight':
            # Run the exact inherited price-loop smoke with all numerical gates unchanged.
            import phase_b_pilot as ge
            ge_dir=args.out/'ge_loop_mock';ge_dir.mkdir()
            fake=dict(P=Q,q_ref=.8,out=ge_dir,selected_d_bar=0.,fp=type('FP',(),{'write':staticmethod(write)})(),prepared=object(),manifest={},objective={},runtime=object(),reference={})
            budget=type('Budget',(),{'remaining_lifecycle':20,'stage_deadline_seconds':300,'deadline_epoch':time.time()+3600})()
            result=ge.smoke_phase_b(fake,budget)
            write(args.out/'completed.json',dict(status='preflight_passed',lifecycle_solves=0,price_loop=result));return
        require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_CPUS_PER_TASK','1')=='1','Native run requires one-core Torch')
        evaluate=pilot.native_evaluator(args.out,'nonnegative_mean',Q,grid,deadline)
        import phase_b_pilot as ge
        native_observe=ge.observe_price
        def observe_fixed(ctx,live,label,*,final=False):
            ctx['free_coordinates']=[]
            result=native_observe(ctx,live,label,final=final)
            if final:
                folder=Path(ctx['out'])/'phase_b_ge'/label
                rows=pilot.readtable(folder/'parameters.csv')
                for row in rows:
                    row['near_bound']=''
                    if row['parameter']=='h_P':row['status']='experimental fixed historical parenthood physical-room floor' if ARM=='floor' else 'fixed zero housing floor'
                    elif row['parameter']=='delta_alpha_jump':row['status']='fixed common pilot-selected share loading' if ARM=='no_A' else 'experimental fixed zero child share loading'
                    elif row['parameter'] in inputs.PARAMETERS:row['status']='fixed common pilot-selected parameter; no refit'
                    elif row['parameter']=='utility_reference_rent':row['status']='retained inactive reference rent; compensation off'
                pilot.table(folder/'parameters.csv',rows)
            return result
        ge.observe_price=observe_fixed
        results=[]
        for label in ('000_baseline','001_'+ARM+'_full_ge_repeat'):
            write(args.out/'latest.json',dict(status='running',label=label,deadline_epoch=deadline,refit=False))
            r=evaluate(label,PLAN['seed']);require(r['status']=='passed','Floor GE did not pass');results.append(r)
            write(args.out/'latest_completed.json',r)
        repeated=pilot.compare_repeated(Path(results[0]['report']),Path(results[1]['report']))
        write(args.out/'completed.json',dict(status='fixed_parameter_variant_verified',arm=ARM,control_reused=True,control_loss=PLAN['control_loss'],variant_loss=float(np.asarray(results[0]['residual'])@np.asarray(results[0]['residual'])),results=results,exact_full_ge_repeat=repeated,refit=False,experimental_not_adopted=True))
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),no_auto_retry=True));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
