"""Bounded iterative full-GE calibration: three four-hour lanes."""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, math, os, signal, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import inputs
HERE=Path(__file__).resolve().parent
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

def expected_parameters(point,dimensions=(120,9)):
    x={r['parameter']:float(r['estimate']) for r in PLAN['reference_parameter_table']}
    x.update(point,wealth_grid_nodes=dimensions[0],income_states=dimensions[1])
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
    return {k:max(1e-6,min(.02*max(abs(seed[k]),.01),.005*(bounds[k][1]-bounds[k][0]))) for k in inputs.PARAMETERS}


def reserve_seconds(lane,observed):
    size=inputs.LANES[lane]['size']
    return max(float(PLAN['final_reserve_seconds'][size]),3*observed,710+2*1.3*observed)


def search(out,seed,bounds,evaluate,deadline,*,mock=False,lane='nonnegative_mean_120x9'):
    """Exact same iterative controller for zero-solve mock and native run."""
    out=Path(out);cases=[];best=None;started=time.time();rounds=[]
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
        identification=dict(round_index=round_index,center_label=center['label'],center_parameters=center['parameters'],free_coordinates=list(inputs.PARAMETERS),scored_targets=10,valid_probes=0,fresh_at_selected=False,status='incomplete_jacobian')
        for k in inputs.PARAMETERS:
            if not enough():budget_hit=True;break
            point=dict(center['parameters']);h=steps[k]
            if point[k]+h>bounds[k][1]:h=-h
            require(bounds[k][0]<=point[k]+h<=bounds[k][1],'Probe outside bounds')
            point[k]+=h;probe=run(f'{len(cases):03d}_r{round_index}_probe_{k}',point,'finite_difference')
            if probe['status']=='passed':probe_rows.append(probe);probe_steps.append(h)
            else:
                budget_hit=probe['status']=='budget_exhausted';break
        identification['valid_probes']=len(probe_rows)
        if len(probe_rows)==9:
            spans=np.asarray([bounds[k][1]-bounds[k][0] for k in inputs.PARAMETERS])
            J=np.column_stack([(np.asarray(row['residual'])-center['residual'])/(h/spans[i]) for i,(row,h) in enumerate(zip(probe_rows,probe_steps))])
            sv=np.linalg.svd(J,compute_uv=False);rank=int(np.linalg.matrix_rank(J))
            identification.update(status='round_center_local_jacobian' if rank==9 else 'underidentified_local_jacobian',rank=rank,singular_values=sv.tolist(),normalized_coordinate_jacobian=J.tolist())
            if rank==9:
                rr=np.asarray(center['residual']);ridge=max(float(sv[0]**2)*1e-4,1e-10)
                delta=np.linalg.solve(J.T@J+ridge*np.eye(9),-J.T@rr)
                trust=np.asarray([3*steps[k]/spans[i] for i,k in enumerate(inputs.PARAMETERS)])
                delta=np.clip(delta,-trust,trust)
                identification.update(ridge=ridge,normalized_trust=trust.tolist(),normalized_step=delta.tolist())
                for damping in (.5,.2,1.):
                    if not enough():budget_hit=True;break
                    point={k:float(np.clip(center['parameters'][k]+damping*delta[i]*spans[i],*bounds[k])) for i,k in enumerate(inputs.PARAMETERS)}
                    proposal=run(f'{len(cases):03d}_r{round_index}_gn_{damping}',point,'damped_Gauss_Newton')
                    if proposal['status']=='budget_exhausted':budget_hit=True;break
        identification['budget_stop']=budget_hit
        rounds.append(identification);write(out/'rounds.json',rounds);write(out/'identification.json',identification)
        if len(probe_rows)!=9:stop='incomplete_jacobian_budget' if budget_hit else 'incomplete_jacobian';break
        if budget_hit:stop='search_budget_exhausted';break
        if identification['rank']!=9:stop='rank_deficient';break
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


def native_evaluator(out,lane,P,grid,deadline):
    arm=inputs.LANES[lane]['arm'];dims=inputs.LANES[lane]['dimensions']
    # Import the same indexed full-GE integration, never the default CLI.
    sys.path.insert(0,str(BASE))
    spec=importlib.util.spec_from_file_location('pilot_matched_workflow',BASE/'run_comparison.py')
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    sys.path.insert(0,str(HERE))
    import phase_b_pilot as ge
    ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=out))
    base.authored.authenticate_frozen(ctx)
    ctx.update(P=P,b_grid=grid,selected_d_bar=inputs.ARMS[arm],reference_psi=float(P.psi_child),expected_dimensions={'wealth_grid_nodes':dims[0],'income_states':dims[1]},free_coordinates=list(inputs.PARAMETERS))
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
    def evaluate(label,point,evaluation_deadline):
        directory=out/label;directory.mkdir()
        candidate=dict(ctx);candidate.update(P=inputs.bind(P,point,bounds),out=directory,expected_parameters=expected_parameters(point,dims),deadline_epoch=evaluation_deadline)
        # No direct field may silently fail to propagate through the imported engine.
        from small_credit_lab import credit
        credit.bind_engine_credit(candidate['P'],'corrected',inputs.ARMS[arm])
        actual=candidate['fp'].actual_parameters(candidate['prepared'],candidate['P'],grid)
        ge.validate_parameter_estimates(candidate,PLAN['reference_parameter_table'],actual)
        budget=base.ArmBudget(directory,evaluation_deadline)
        write(directory/'proposed_parameters.json',dict(free=point,actual=actual,credit=inputs.ARMS[arm],target_contract_sha256=inputs.canonical(PLAN['target_contract'])))
        try:
            solved=ge.run_phase_b(candidate,dict(selected_d_bar=inputs.ARMS[arm]),budget)
            if solved['status']!='passed':return dict(status='inadmissible_numerical',reason='GE did not converge in bounded budget',lifecycle_solves=budget.used_lifecycle)
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
        return dict(status='passed',residual=rr.tolist(),report=str(report),lifecycle_solves=budget.used_lifecycle,price=solved['selected_price'],population=solved['selected']['population_scale'])
    return evaluate


def mocked_evaluator(out,seed,bounds,lane):
    # A deterministic dense ten-by-nine residual system exercises the actual loop.
    matrix=np.vstack((np.eye(9),np.ones((1,9))*.2))
    spans=np.asarray([bounds[k][1]-bounds[k][0] for k in inputs.PARAMETERS])
    P,_=inputs.proposal(lane)
    def evaluate(label,point,evaluation_deadline):
        require(set(point)==set(seed),'Mock coordinate binding failed')
        Q=inputs.bind(P,point,bounds)
        require(Q.beta==point['beta_annual']**Q.period_years and Q.eps_fert==point['kappa_fert'] and Q.rho==Q.rho_hat==1/Q.beta-1,'Mock transformed binding failed')
        require(all(float(getattr(Q,k))==point[k] for k in inputs.PARAMETERS if k!='beta_annual'),'Mock direct binding failed')
        rr=matrix@((np.asarray([point[k]-seed[k] for k in inputs.PARAMETERS]))/spans)-np.linspace(.1,.2,10)
        return dict(status='passed',residual=rr.tolist(),report=str(out/label),lifecycle_solves=0,price=.8,population=1.)
    return evaluate


def phase_b_mock_checks(out, P, lane, deadline):
    arm=inputs.LANES[lane]['arm'];dims=inputs.LANES[lane]['dimensions']
    # The copied PhaseB imports only the at-price function; engine imports are lazy.
    matched=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
    sys.path.insert(0,str(matched))
    import phase_b_pilot as ge
    fake=dict(P=P,q_ref=.8,out=out,selected_d_bar=inputs.ARMS[arm],fp=SimpleNamespace(write=write),prepared=object(),manifest={},objective={},runtime=object(),reference={})
    # Mock GE has a logical3600s budget; the real preflight alarm remains300s.
    budget=SimpleNamespace(remaining_lifecycle=20,stage_deadline_seconds=300,deadline_epoch=time.time()+3600)
    result=ge.smoke_phase_b(fake,budget)
    seed,_,rows=inputs.seed_and_bounds(lane);expected=expected_parameters(seed,dims)
    ctx=dict(expected_parameters=expected)
    ge.validate_parameter_estimates(ctx,rows,expected)
    rejected=0
    for row in rows:
        actual=dict(expected);actual[row['parameter']]+=1
        try:ge.validate_parameter_estimates(ctx,rows,actual)
        except RuntimeError:rejected+=1
        else:raise RuntimeError('Actual parameter drift accepted '+row['parameter'])
    require(rejected==31,'Wrong stale parameter negative-test count')
    result['rejected_actual_parameter_drifts']=rejected
    write(Path(out)/'phase_b_mock_checks.json',result)
    return result


def alarm(_signum,_frame):raise TimeoutError('Hard four-hour round deadline reached')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=('preflight','run'),required=True);ap.add_argument('--lane',choices=tuple(inputs.LANES),required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-seconds',type=float,default=14400);ap.add_argument('--deadline-epoch',type=float);args=ap.parse_args()
    require(not args.out.exists(),'Refusing existing pilot output');args.out.mkdir(parents=True)
    started=time.time();deadline=min(started+min(args.deadline_seconds,14400),args.deadline_epoch if args.deadline_epoch is not None else float('inf'))
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.01,deadline-started))
    try:
        verify_sources();P,grid=inputs.proposal(args.lane);arm=inputs.LANES[args.lane]['arm'];P,entry=inputs.entry(P,grid,arm)
        seed,bounds,rows=inputs.seed_and_bounds(args.lane)
        write(args.out/'input_contract.json',dict(entry=entry,seed=seed,bounds=bounds,target_contract=PLAN['target_contract'],target_contract_sha256=inputs.canonical(PLAN['target_contract']),fixed_H0=float(P.H0[0]),fixed_psi=float(P.psi_child),annual_interest=(1+P.q)**(1/P.period_years)-1,closure='birth-renewal price; population clears absolute housing supply',economic_changes='Two author-approved entry/credit specifications, three numerical lanes; same 2% saving/borrowing rate, homeowner financing and fixed benefit; no silent feasibility relocation; five-bin rather than raw-survey censoring'))
        np.savez_compressed(args.out/'initial_entry_distribution.npz',b_grid=grid,z_grid=P.z_grid,z_weights=P.z_weights,conditional=P.fixed_reference_entry_conditional)
        if args.mode=='run':require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Native calibration requires one-core Torch Slurm')
        if args.mode=='preflight':phase_b_mock_checks(args.out/'ge_loop_mock',P,args.lane,deadline)
        evaluator=mocked_evaluator(args.out,seed,bounds,args.lane) if args.mode=='preflight' else native_evaluator(args.out,args.lane,P,grid,deadline)
        result=search(args.out,seed,bounds,evaluator,deadline,mock=args.mode=='preflight',lane=args.lane)
        if args.mode=='preflight':require(result['completed_full_ge']==PLAN['maximum_full_ge'][inputs.LANES[args.lane]['size']] and result['lifecycle_solves']==0 and result['selected_repeats']==2,'Exact mocked iterative search loop incomplete')
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-started,no_auto_retry=True));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
