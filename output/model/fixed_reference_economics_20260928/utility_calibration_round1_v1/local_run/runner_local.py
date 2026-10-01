"""Bounded iterative full-GE calibration: three four-hour lanes."""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, math, os, signal, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import inputs
HERE=Path(__file__).resolve().parent.parent
ROOT=HERE.parents[3]
BASE=ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner'
PLAN=json.loads((HERE/'plan.json').read_text())
require=inputs.require

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=str)+'\n');temp.replace(path)

def table(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def readtable(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def target_identity(rows):return [{k:r[k] for k in ('moment','target','weight','role')} for r in rows]

def residual(rows):
    require(len(rows)==14 and target_identity(rows)==PLAN['target_contract'],'Target-and-weight contract drift')
    scored=[r for r in rows if r['role']=='scored']
    require(len(scored)==10,'Ten scored targets required for nine free coordinates')
    r=np.asarray([math.sqrt(float(x['weight']))*float(x['gap']) for x in scored])
    require(np.isfinite(r).all(),'Nonfinite scored residual')
    require(abs(float(r@r)-sum(float(x['loss_contribution']) for x in scored))<1e-8,'Native target loss inconsistent')
    return r

def expected_parameters(point,dimensions=(120,9),arm='floor'):
    x={r['parameter']:float(r['estimate']) for r in PLAN['reference_parameter_table']}
    x.update(point,wealth_grid_nodes=dimensions[0],income_states=dimensions[1],delta_alpha=0.,h_P=float(point['h_P']) if arm=='floor' else 0.,delta_alpha_jump=float(point['delta_alpha_jump']) if arm=='no_A' else 0.)
    x['child_benefit_CRRA_coefficient']=(1-point['child_benefit_curvature'])*x['psi_child']
    return x

def verify_sources():
    pins=json.loads((HERE/'source_pins.json').read_text())
    for rel,digest in pins.items():require(inputs.sha(ROOT/rel)==digest,'Pilot source drift: '+rel)


def compare_repeated(a,b):
    for name,count in (('target_fit.csv',14),('parameters.csv',31)):
        x,y=readtable(a/name),readtable(b/name)
        require(len(x)==len(y)==count and x==y,'Repeated full GE differs: '+name)
    require(json.loads((a/'closure.json').read_text())==json.loads((b/'closure.json').read_text()),'Repeated full GE closure differs')
    aa={p.name:inputs.sha(p) for p in (a/'standard_diagnostics').glob('*.png')}
    bb={p.name:inputs.sha(p) for p in (b/'standard_diagnostics').glob('*.png')}
    require(len(aa)==17 and aa==bb,'Repeated full GE standard plots differ')
    return dict(status='exact_full_ge_repeat_passed',target_rows=14,parameter_rows=31,standard_plot_hashes=aa)


def step_sizes(seed,bounds):
    return {k:max(1e-6,min(.02*max(abs(seed[k]),.01),.005*(bounds[k][1]-bounds[k][0]))) for k in seed}


def reserve_seconds(lane,observed):
    size=inputs.LANES[lane]['size']
    return max(float(PLAN['final_reserve_seconds'][size]),3*observed,710+2*1.3*observed)


def search(out,seed,bounds,evaluate,deadline,*,mock=False,lane='floor_s0'):
    """Exact same iterative controller for zero-solve mock and native run."""
    out=Path(out);cases=[];best=None;started=time.time();rounds=[];coordinates=inputs.parameters(lane);dimension=len(coordinates)
    size=inputs.LANES[lane]['size'];cap=PLAN['maximum_full_ge'][size]
    minimum=PLAN['minimum_ge_seconds'][size]
    # Freeze the reserve before exploration; increase after measured slower solves.
    observed=minimum;reserve=reserve_seconds(lane,observed)
    def run(label,point,kind):
        nonlocal best,observed,reserve
        inputs.check_point(point,bounds)
        require(len(cases)<cap,'Full-GE cap reached')
        require(time.time()<deadline,'Global four-hour deadline reached')
        t=time.monotonic()
        evaluation_deadline=deadline if kind=='selected_repeat' else deadline-reserve
        write(out/'latest.json',dict(status='running_full_GE',label=label,kind=kind,parameters=point,completed_full_ge=len(cases),deadline_epoch=deadline,evaluation_deadline_epoch=evaluation_deadline,final_reserve_seconds=reserve))
        result=evaluate(label,point,evaluation_deadline)
        row=dict(label=label,kind=kind,parameters=dict(point),seconds=time.monotonic()-t,**result)
        cases.append(row)
        if row['status']=='passed':
            rr=np.asarray(row['residual']);require(rr.size==10 and np.isfinite(rr).all(),'Wrong residual dimension')
            row['loss']=float(rr@rr)
            observed=max(observed,row['seconds']);reserve=reserve_seconds(lane,observed)
            if kind not in ('baseline_repeat','selected_repeat') and (best is None or row['loss']<best['loss']):best=row
        write(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(cases),elapsed_seconds=time.time()-started))
        write(out/'best_so_far.json',dict(status='provisional_until_two_final_full_GE_repeats',best=best))
        write(out/'cases.json',cases)
        return row
    def enough(required=1):
        estimate=max(minimum,1.3*observed)
        return len(cases)+required+2<=cap and (mock or deadline-time.time()>reserve+max(710.,required*estimate))
    require(enough(2),'Insufficient total budget for baseline/repeat and verification reserve')
    baseline=run('000_baseline',seed,'baseline');require(baseline['status']=='passed','Baseline full GE failed; search blocked')
    twin=run('001_baseline_repeat',seed,'baseline_repeat');require(twin['status']=='passed','Baseline repeated full GE failed; search blocked')
    require(baseline['residual']==twin['residual'],'Baseline repeated residual differs')
    if not mock:write(out/'baseline_repeat.json',compare_repeated(Path(baseline['report']),Path(twin['report'])))
    stop='maximum_rounds'
    for round_index in range(PLAN['maximum_rounds'][size]):
        if not enough():stop='budget_before_fresh_jacobian';break
        center=dict(best);steps=step_sizes(center['parameters'],bounds)
        probe_rows=[];probe_steps=[];budget_hit=False
        identification=dict(round_index=round_index,center_label=center['label'],center_parameters=center['parameters'],free_coordinates=list(coordinates),scored_targets=10,valid_probes=0,fresh_at_selected=False,status='incomplete_jacobian')
        for k in coordinates:
            if not enough():budget_hit=True;break
            point=dict(center['parameters']);h=steps[k]
            if point[k]+h>bounds[k][1]:h=-h
            require(bounds[k][0]<=point[k]+h<=bounds[k][1],'Probe outside bounds')
            point[k]+=h;probe=run(f'{len(cases):03d}_r{round_index}_probe_{k}',point,'finite_difference')
            if probe['status']=='passed':probe_rows.append(probe);probe_steps.append(h)
            else:
                budget_hit=probe['status']=='budget_exhausted';break
        identification['valid_probes']=len(probe_rows)
        if len(probe_rows)==dimension:
            spans=np.asarray([bounds[k][1]-bounds[k][0] for k in coordinates])
            J=np.column_stack([(np.asarray(row['residual'])-center['residual'])/(h/spans[i]) for i,(row,h) in enumerate(zip(probe_rows,probe_steps))])
            sv=np.linalg.svd(J,compute_uv=False);rank=int(np.linalg.matrix_rank(J))
            identification.update(status='round_center_local_jacobian' if rank==dimension else 'underidentified_local_jacobian',rank=rank,singular_values=sv.tolist(),normalized_coordinate_jacobian=J.tolist())
            if rank==dimension:
                rr=np.asarray(center['residual']);ridge=max(float(sv[0]**2)*1e-4,1e-10)
                delta=np.linalg.solve(J.T@J+ridge*np.eye(dimension),-J.T@rr)
                trust=np.asarray([3*steps[k]/spans[i] for i,k in enumerate(coordinates)])
                delta=np.clip(delta,-trust,trust)
                identification.update(ridge=ridge,normalized_trust=trust.tolist(),normalized_step=delta.tolist())
                for damping in (.5,.2,1.):
                    if not enough():budget_hit=True;break
                    point={k:float(np.clip(center['parameters'][k]+damping*delta[i]*spans[i],*bounds[k])) for i,k in enumerate(coordinates)}
                    proposal=run(f'{len(cases):03d}_r{round_index}_gn_{damping}',point,'damped_Gauss_Newton')
                    if proposal['status']=='budget_exhausted':budget_hit=True;break
        identification['budget_stop']=budget_hit
        rounds.append(identification);write(out/'rounds.json',rounds);write(out/'identification.json',identification)
        if len(probe_rows)!=dimension:stop='incomplete_jacobian_budget' if budget_hit else 'incomplete_jacobian';break
        if budget_hit:stop='search_budget_exhausted';break
        if identification['rank']!=dimension:stop='rank_deficient';break
        if best['label']==center['label']:stop='no_local_improvement';break
    identification=rounds[-1] if rounds else dict(status='no_jacobian_budget',fresh_at_selected=False,valid_probes=0)
    write(out/'identification.json',identification)
    require(best is not None,'No valid calibration point')
    selected=dict(best);write(out/'selected_provisional.json',selected)
    repeats=[]
    for i in range(2):
        if not mock and deadline-time.time()<=710:break
        row=run(f'{len(cases):03d}_selected_repeat_{i}',selected['parameters'],'selected_repeat')
        require(row['status']=='passed' and row['residual']==selected['residual'],'Selected full GE repeat failed')
        if not mock:write(out/f'selected_repeat_{i}.json',compare_repeated(Path(selected['report']),Path(row['report'])))
        repeats.append(row)
    status='selected_verified' if len(repeats)==2 else 'provisional_budget_exhausted'
    result=dict(status=status,lane=lane,selected=selected,baseline=baseline,selected_repeats=len(repeats),completed_full_ge=len(cases),maximum_full_ge=cap,rounds_completed=len(rounds),search_stop_reason=stop,identification=identification,lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in cases),elapsed_seconds=time.time()-started,deadline_epoch=deadline,experimental_not_adopted=True)
    write(out/'completed.json',result)
    return result


def is_search_budget_exit(exc,budget,global_deadline):
    if not budget.deadline_epoch<global_deadline:return False
    remaining=budget.deadline_epoch-time.time()
    guards={'No time reserve for selected reporting and exact repeat':budget.stage_deadline_seconds+400,
            'No time for exact repeat':budget.stage_deadline_seconds,
            'No GE solve and repeat reserve':0.,'Price or total deadline exceeded':0.}
    return str(exc) in guards and remaining<=guards[str(exc)]


def native_evaluator(out,lane,P,grid,deadline,price_start=None):
    coordinates=inputs.parameters(lane)
    arm=inputs.LANES[lane]['arm'];dims=inputs.LANES[lane]['dimensions']
    # Import the same indexed full-GE integration, never the default CLI.
    sys.path.insert(0,str(BASE))
    spec=importlib.util.spec_from_file_location('pilot_matched_workflow',BASE/'run_comparison.py')
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    sys.path.insert(0,str(HERE))
    import phase_b_pilot as ge
    ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=out))
    install_reporter_on_authored(base.authored)
    base.authored.authenticate_frozen(ctx)
    ctx.update(P=P,b_grid=grid,selected_d_bar=inputs.ARMS[arm],reference_psi=float(P.psi_child),expected_dimensions={'wealth_grid_nodes':dims[0],'income_states':dims[1]},free_coordinates=list(coordinates),price_start=float(price_start) if price_start is not None else float(ctx['q_ref']),phase_b_max_new_lifecycle=32)
    from small_credit_lab.engine.shared import annual_gross_income_at_state
    from small_credit_lab.engine import solver
    sd=solver.precompute_shared(P,grid)
    require(float(sd.cb_flat[0,0])==float(sd.hb_flat[0,0])==float(sd.gb_flat[0,0])==0.,'Necessary cash preflight missing childless floors')
    actual_income=np.asarray([annual_gross_income_at_state(P,0,0,float(z)) for z in P.z_grid])
    np.testing.assert_array_equal(actual_income,P.income[0,0]*P.z_grid/P.period_years/(1-P.tau_pay))
    cal=ctx['prepared'].rt['primitive'].pf.calendar
    cohort=cal.entrant_cohort(np.asarray([1.]),P,grid)
    np.testing.assert_allclose(cohort.sum(axis=(1,2,4,5)),P.fixed_reference_entry_conditional*P.z_weights[None,:],rtol=0,atol=2e-16)
    require(abs(cohort.sum()-1)<2e-12,'Entrant mass changed')
    write(out/'native_entry_verification.json',dict(calendar_joint_maximum_error=float(abs(cohort.sum(axis=(1,2,4,5))-P.fixed_reference_entry_conditional*P.z_weights[None,:]).max()),income_formula_exact=True,childless_cash_floors_zero=True,lifecycle_solves=0))
    seed,bounds,_=inputs.seed_and_bounds(lane)
    for row in ctx['manifest']['full_parameter_table']:
        if row['parameter'] in bounds:
            row['lower'],row['upper']=map(str,bounds[row['parameter']])
    install_observer_metadata(ge,arm,bounds)
    from small_credit_lab import credit
    credit.bind_engine_credit(P,'corrected',0.)
    actual=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)
    expected=expected_parameters(seed,dims,arm)
    ge.validate_parameter_estimates(dict(expected_parameters=expected),PLAN['reference_parameter_table'],actual)
    live_sd=solver.precompute_shared(P,grid)
    from refactor_lab.engine import solver as checked_solver
    checked_sd=checked_solver.precompute_shared(P,grid)
    for key in ('h_bar','c_bar','g_bar','alpha_flat','psi_v','escale_flat'):
        np.testing.assert_array_equal(getattr(live_sd,key),getattr(checked_sd,key))
    require(float(live_sd.h_bar[1,1])==expected['h_P'] and float(live_sd.h_bar[1,0])==0.,'Executed physical floor differs')
    import inspect
    from small_credit_lab.engine import shared as executed_shared,child_preferences as executed_child,kernels as executed_kernels
    from refactor_lab.engine import shared as checked_shared,child_preferences as checked_child,kernels as checked_kernels
    for checked,executed in [(checked_shared,executed_shared),(checked_child,executed_child),(checked_kernels,executed_kernels)]:
        require(inputs.sha(Path(checked.__file__))==inputs.sha(Path(executed.__file__)),'Executed preference source differs')
    write(out/'native_initializer_verification.json',dict(status='passed_zero_solve',actual_parameters=actual,expected_parameters=expected,authenticated_once=True,lifecycle_solves=0))
    def evaluate(label,point,evaluation_deadline):
        directory=out/label;directory.mkdir()
        candidate=dict(ctx);candidate.update(P=inputs.bind(P,point,bounds,arm),out=directory,expected_parameters=expected_parameters(point,dims,arm),deadline_epoch=evaluation_deadline)
        # No direct field may silently fail to propagate through the imported engine.
        from small_credit_lab import credit
        credit.bind_engine_credit(candidate['P'],'corrected',inputs.ARMS[arm])
        actual=candidate['fp'].actual_parameters(candidate['prepared'],candidate['P'],grid)
        ge.validate_parameter_estimates(candidate,PLAN['reference_parameter_table'],actual)
        budget=base.ArmBudget(directory,evaluation_deadline);budget.max_lifecycle=32
        write(directory/'proposed_parameters.json',dict(free=point,actual=actual,credit=inputs.ARMS[arm],starting_price=candidate['price_start'],target_contract_sha256=inputs.canonical(PLAN['target_contract'])))
        try:
            solved=ge.run_phase_b(candidate,dict(selected_d_bar=inputs.ARMS[arm]),budget)
            if solved['status']!='passed':
                exhausted=solved['status']=='uncomputed_bounded_budget' or solved.get('price_search',{}).get('termination_reason')=='budget_or_repeat_reserve'
                return dict(status='budget_exhausted' if exhausted else 'inadmissible_numerical',reason=solved['status'],lifecycle_solves=budget.used_lifecycle,price_search=solved.get('price_search'))
        except TimeoutError as exc:
            if evaluation_deadline<deadline and time.time()>=evaluation_deadline:
                return dict(status='budget_exhausted',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            raise
        except RuntimeError as exc:
            if is_search_budget_exit(exc,budget,deadline):
                write(directory/'search_budget_exhausted.json',dict(reason=str(exc),lifecycle_solves=budget.used_lifecycle))
                return dict(status='budget_exhausted',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            # Only explicit failure to locate a root can reject an exploratory proposal.
            # Accounting, feasibility, source and parameter failures halt the arm.
            if label!='000_baseline' and str(exc).startswith('Renewal root unbracketed'):
                write(directory/'inadmissible.json',dict(reason=str(exc),lifecycle_solves=budget.used_lifecycle))
                return dict(status='inadmissible_numerical',reason=str(exc),lifecycle_solves=budget.used_lifecycle)
            raise
        report=directory/'phase_b_ge/selected_root'
        rows=readtable(report/'target_fit.csv');rr=residual(rows)
        # The native selected-price repeat verifies arrays/tables; verify actual PNGs too.
        repeat=directory/'phase_b_ge/selected_repeat_final'
        write(directory/'native_selected_repeat.json',compare_repeated(report,repeat))
        return dict(status='passed',residual=rr.tolist(),report=str(report),lifecycle_solves=budget.used_lifecycle,price=solved['selected_price'],starting_price=candidate['price_start'],population=solved['selected']['population_scale'])
    return evaluate


def install_reporter_on_authored(authored):
    original=authored.authenticate_frozen
    def authenticate(ctx):
        result=original(ctx);fp=ctx['fp']
        if not getattr(fp,'_utility_floor_reporting_adapter',False):
            native=fp.actual_parameters
            def actual(prepared,P,grid):
                values=native(prepared,P,grid)
                require(float(P.hbar_child_rooms)==0.,'Unexpected later-child floor')
                values['h_P']=float(P.hbar_first_child_jump) if P.child_room_floor else 0.
                return values
            fp.actual_parameters=actual;fp._utility_floor_reporting_adapter=True
        return result
    authored.authenticate_frozen=authenticate


def install_observer_metadata(ge,arm,bounds):
    original=ge.observe_price
    def observe(ctx,live,label,*,final=False):
        result=original(ctx,live,label,final=final)
        if final:
            path=Path(ctx['out'])/'phase_b_ge'/label/'parameters.csv';rows=readtable(path)
            for row in rows:
                key=row['parameter']
                if key in bounds:
                    row['lower'],row['upper']=map(str,bounds[key]);row['status']='free in experimental utility calibration'
                    lo,hi=bounds[key];x=float(row['estimate']);row['near_bound']=str(min(x-lo,hi-x)<=.01*(hi-lo))
                elif key in ('h_P','delta_alpha_jump','delta_alpha'):
                    row['status']='fixed zero under experimental utility contract';row['lower']=row['upper']=row['near_bound']=''
                elif key=='utility_reference_rent':row['status']='retained inactive normalization; compensation off'
            table(path,rows)
        return result
    ge.observe_price=observe


def utility_checks(P,grid,lane,out):
    arm=inputs.LANES[lane]['arm'];seed,bounds,_=inputs.seed_and_bounds(lane)
    Q=inputs.bind(P,seed,bounds,arm)
    sys.path.insert(0,str(ROOT/'code/model'))
    from refactor_lab.engine import solver
    sd=solver.precompute_shared(Q,grid)
    require(Q.preference_spec=='eqscale' and Q.eqscale_form=='power' and Q.child_state_mode=='independent_count','Wrong utility architecture')
    require(not Q.compensated_child_housing_shares and Q.delta_alpha==Q.hbar_child_rooms==0.,'Utility compensation retained')
    for n in range(Q.n_parity):
        for m in range(min(n+1,Q.n_child_states)):
            floor=seed['h_P'] if m>0 and arm=='floor' else 0.
            require(sd.h_bar[n,m]==floor and sd.c_bar[n,m]==sd.g_bar[n,m]==0.,'Physical-room floor differs')
            alpha=float(sd.alpha_flat.reshape((Q.n_parity,Q.n_child_states),order='F')[n,m])
            expected_alpha=Q.alpha_cons-Q.delta_alpha_jump if m>0 and arm=='no_A' else Q.alpha_cons
            require(alpha==expected_alpha,'Housing share differs')
            require(sd.psi_v[n,m]==(Q.psi_child*m**(1-Q.child_benefit_curvature) if m else 0.),'Nonlinear benefit differs')
            e=((2+.7*m)/2)**.7
            require(abs(sd.escale_flat.reshape((Q.n_parity,Q.n_child_states),order='F')[n,m]-e**(Q.sigma-1))<2e-15,'Equivalence scale changed')
    require(Q.hR_max>seed.get('h_P',0.),'Rental cap below floor')
    write(out/'utility_verification.json',dict(status='passed_zero_solve',arm=arm,h_P=Q.hbar_first_child_jump,delta_alpha_jump=Q.delta_alpha_jump,compensation=False,nonlinear_child_benefit=True,owner_services='chi*(physical_rooms-h_P)',renter_services='physical_rooms-h_P',lifecycle_solves=0))
    return Q


def mocked_evaluator(out,seed,bounds,lane):
    coordinates=inputs.parameters(lane);dim=len(coordinates)
    matrix=np.vstack((np.eye(dim),np.ones((10-dim,dim))*.2))
    spans=np.asarray([bounds[k][1]-bounds[k][0] for k in coordinates]);P,_=inputs.proposal(lane);arm=inputs.LANES[lane]['arm']
    def evaluate(label,point,end):
        Q=inputs.bind(P,point,bounds,arm)
        require(Q.beta==point['beta_annual']**Q.period_years and Q.eps_fert==point['kappa_fert'],'Mock transformed binding failed')
        for k in coordinates:
            if k=='beta_annual':continue
            actual=Q.hbar_first_child_jump if k=='h_P' else getattr(Q,k)
            require(actual==point[k],'Mock parameter binding failed '+k)
        r=matrix@np.asarray([(point[k]-seed[k])/spans[i] for i,k in enumerate(coordinates)])-np.linspace(.1,.2,10)
        return dict(status='passed',residual=r.tolist(),report=str(out/label),lifecycle_solves=0,price=.8,population=1.)
    return evaluate


def source_fingerprint():
    verify_sources()
    return inputs.canonical(json.loads((HERE/'source_pins.json').read_text()))


def report_hashes(folder):
    folder=Path(folder)
    files=[folder/'target_fit.csv',folder/'parameters.csv',folder/'closure.json',*sorted((folder/'standard_diagnostics').glob('*.png'))]
    require(len(files)==20,'Smoke full reporting incomplete')
    return {str(p.relative_to(folder)):inputs.sha(p) for p in files}


def validate_smoke_receipt(path,arm):
    path=Path(path);r=json.loads(path.read_text());require(r['status']=='utility_smoke_verified' and r['arm']==arm,'Wrong smoke gate')
    require(r['source_fingerprint']==source_fingerprint(),'Smoke source fingerprint drift')
    require(r['target_contract_sha256']==inputs.canonical(PLAN['target_contract']),'Smoke target fingerprint drift')
    require(r['main_seed']==inputs.LANES[arm+'_s0']['seed'],'Smoke seed drift')
    require(r['independent_full_ge_repeats']==1 and r['full_ge_count']==2,'Smoke repeats missing')
    for report in r['reports']:
        require(report_hashes(path.parent/report['relative_path'])==report['file_sha256'],'Smoke receipt artifact drift')
    price=float(r['selected_price']);require(math.isfinite(price) and price>0,'Invalid smoke root')
    return r


def smoke(out,lane,evaluate,deadline,*,mock=False):
    seed,bounds,_=inputs.seed_and_bounds(lane);results=[]
    for label in ('000_baseline','001_baseline_repeat'):
        write(out/'latest.json',dict(status='running_full_GE',label=label,parameters=seed,deadline_epoch=deadline))
        row=evaluate(label,seed,deadline);require(row['status']=='passed','Smoke full GE failed; this arm search blocked')
        results.append(row);write(out/'latest_completed.json',row);write(out/'best_so_far.json',row);write(out/'cases.json',results)
    require(results[0]['residual']==results[1]['residual'] and results[0]['price']==results[1]['price'],'Smoke repeat differs')
    if not mock:compare_repeated(Path(results[0]['report']),Path(results[1]['report']))
    receipt=dict(status='utility_smoke_verified' if not mock else 'mock_smoke_verified',arm=inputs.LANES[lane]['arm'],source_fingerprint=source_fingerprint(),target_contract_sha256=inputs.canonical(PLAN['target_contract']),main_seed=seed,starting_price=results[0].get('starting_price'),selected_price=results[0]['price'],independent_full_ge_repeats=1,full_ge_count=2,lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in results),reports=[])
    if not mock:
        receipt['reports']=[dict(relative_path=str(Path(r['report']).relative_to(out)),file_sha256=report_hashes(r['report'])) for r in results]
    write(out/'smoke_receipt.json',receipt);write(out/'completed.json',receipt);return receipt


def alarm(*_):raise TimeoutError('Utility calibration hard deadline reached')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=('preflight','smoke','run'),required=True);ap.add_argument('--arm',choices=('floor','no_A','constant_alpha'),required=True);ap.add_argument('--start',type=int,choices=range(6),default=0);ap.add_argument('--stage',choices=('smoke','search'),default='search');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--smoke-receipt',type=Path);ap.add_argument('--deadline-seconds',type=float);ap.add_argument('--deadline-epoch',type=float);args=ap.parse_args()
    require(not args.out.exists(),'Refusing existing output');args.out.mkdir(parents=True)
    lane=f'{args.arm}_s{args.start}';is_smoke=args.mode=='smoke' or args.mode=='preflight' and args.stage=='smoke'
    require(not is_smoke or args.start==0,'Smoke must use exact main seed')
    started=time.time();seconds=1800 if is_smoke else 5400
    deadline=min(started+min(args.deadline_seconds or seconds,seconds),PLAN['common_deadline_epoch'],args.deadline_epoch or float('inf'))
    require(deadline>started,'Absolute deadline already reached')
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,deadline-started)
    try:
        verify_sources();P,grid=inputs.proposal(lane);P,entry=inputs.entry(P,grid,'nonnegative_mean');seed,bounds,_=inputs.seed_and_bounds(lane)
        Q=utility_checks(P,grid,lane,args.out)
        price_start=None;receipt=None
        if args.mode=='run':
            require(args.smoke_receipt is not None,'Search requires verified smoke receipt');receipt=validate_smoke_receipt(args.smoke_receipt,args.arm);price_start=receipt['selected_price']
        write(args.out/'input_contract.json',dict(lane=lane,entry=entry,seed=seed,bounds=bounds,free_coordinates=inputs.parameters(lane),target_contract=PLAN['target_contract'],target_contract_sha256=inputs.canonical(PLAN['target_contract']),fixed_H0=float(P.H0[0]),fixed_psi=float(P.psi_child),economic_changes=PLAN['economic_changes'][args.arm],starting_guesses_not_fixed_changes=True,starting_price=price_start,smoke_receipt_sha256=inputs.sha(args.smoke_receipt) if receipt else None,closure='Price clears birth renewal; population clears absolute housing supply',experimental_not_adopted=True))
        np.savez_compressed(args.out/'initial_entry_distribution.npz',b_grid=grid,z_grid=P.z_grid,z_weights=P.z_weights,conditional=P.fixed_reference_entry_conditional)
        if args.mode=='preflight':
            # Exact native initializer and reporter authentication, zero lifecycle calls.
            initialize=args.out/'native_initializer';initialize.mkdir()
            native_evaluator(initialize,lane,Q,grid,deadline,price_start)
            import phase_b_pilot as ge
            fake_out=args.out/'ge_loop_mock';fake_out.mkdir()
            fake=dict(P=Q,q_ref=.8,price_start=.8,out=fake_out,selected_d_bar=0.,fp=SimpleNamespace(write=write),prepared=object(),manifest={},objective={},runtime=object(),reference={},phase_b_max_new_lifecycle=32)
            budget=SimpleNamespace(remaining_lifecycle=32,stage_deadline_seconds=300,deadline_epoch=time.time()+3600)
            mock_result=ge.smoke_phase_b(fake,budget)
            expected=expected_parameters(seed,arm=args.arm);rows=PLAN['reference_parameter_table']
            for row in rows:
                actual=dict(expected);actual[row['parameter']]+=1.
                try:ge.validate_parameter_estimates(dict(expected_parameters=expected),rows,actual)
                except RuntimeError:pass
                else:raise RuntimeError('Actual parameter drift accepted '+row['parameter'])
            write(args.out/'phase_b_mock_checks.json',dict(price_loop=mock_result,rejected_parameter_drifts=31,lifecycle_solves=0))
            evaluate=mocked_evaluator(args.out,seed,bounds,lane)
        else:
            require(sys.platform=='darwin' and all(os.environ.get(k)=='2' for k in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')),'Explicit local authorization requires two threads per arm')
            evaluate=native_evaluator(args.out,lane,Q,grid,deadline,price_start)
        if is_smoke:result=smoke(args.out,lane,evaluate,deadline,mock=args.mode=='preflight')
        else:result=search(args.out,seed,bounds,evaluate,deadline,mock=args.mode=='preflight',lane=lane)
        if args.mode=='preflight':require(result['lifecycle_solves']==0,'Preflight called lifecycle solver')
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-started,no_auto_retry=True));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
