"""Bounded cluster experiment: terminal -> native smoke -> surprise fitting -> policies."""
from pathlib import Path
import argparse,copy,csv,gzip,hashlib,json,os,pickle,subprocess,sys,threading,time
import numpy as np

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(d,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x),indent=2)+'\n');q.replace(p)

def standard_graphs(snapshot,result,folder):
    import run_e5f_independent_numerical_audit as diagnostic_writer
    from unittest.mock import patch
    original=diagnostic_writer.write_diagnostics
    actual_rent=float(result.path.rows[0]['renter_price'])
    def writer(stats,Q,destination):stats.owner_user_cost=np.array([actual_rent]);return original(stats,Q,destination)
    with patch.object(diagnostic_writer,'write_diagnostics',writer):diagnostic_writer.standard_diagnostics(snapshot,folder,validate_production_young=False)
    if len(list((folder/'standard_diagnostics').glob('*.png')))!=17:raise RuntimeError('Missing standard graphs')

def next_psi(trials, initial, seed_step, bound):
    """Bracketed scalar proposals; do not extrapolate outside explicit bounds."""
    if not trials:return initial+seed_step
    valid=sorted([(p,e) for p,e in trials if np.isfinite(e)])
    brackets=[(a,b) for a,b in zip(valid,valid[1:]) if a[1]*b[1]<0]
    if brackets:
        a,b=min(brackets,key=lambda pair:pair[1][0]-pair[0][0]);x=a[0]-a[1]*(b[0]-a[0])/(b[1]-a[1]);return float(np.clip(x,a[0]+.15*(b[0]-a[0]),b[0]-.15*(b[0]-a[0])))
    if len(valid)>=2:
        a,b=sorted(valid,key=lambda x:abs(x[1]))[:2]
        if abs(a[1]-b[1])>1e-10:return float(np.clip(a[0]-a[1]*(a[0]-b[0])/(a[1]-b[1]),*bound))
    p,e=min(valid,key=lambda x:abs(x[1]))
    return float(np.clip(p-(abs(seed_step) if e>0 else -abs(seed_step)),*bound))

def large_owner_observation(e,P):
    """Small saved sufficient statistics; no policy or equilibrium modifications."""
    if P.child_state_mode!='independent_count':raise ValueError('Dependent-count observer requires independent_count')
    rows=[]
    slots=[k for k,h in enumerate(P.H_own,start=1) if h>=6]
    for j in range(P.J):
        g=e.g_current[:,:,:,j,:,:,:]
        td=g.sum(axis=(0,2,3,4))
        rows.append(dict(age=float(P.age_start+j*P.da),age_width=float(P.da),
            without_children=float(td[slots,0].sum()),with_children=float(td[slots,1:].sum())))
    return rows

def carry_finite_diagnostic(result, *, enabled, inherited, old, demographics, psi, module):
    """Allow experimental state propagation while retaining failed horizon flags.

    Production admission is unchanged. The exact first-period replay is mandatory.
    """
    if not enabled or result.next_state is not None:return result
    receipt=result.root_receipt
    if not receipt.get('finite_horizon_market_fiscal_converged'):return result
    f=receipt['final']
    next_state=module.first_period_state(inherited=inherited,old_state=old,
        demographics=demographics,path=result.path,prices=f['prices'],
        pensions=f['fiscal_values'],psi=psi)
    receipt['diagnostic_finite_horizon_state_carry']=True
    receipt['production_eligible']=False
    receipt['horizon_verified']=False
    realized=dict(result.path.rows[0],forecast_vintage_year=inherited.year,
        expected_next_asset_price=float(f['prices'][1]),expected_constant_psi=float(psi))
    return module.SurpriseResult(result.path,receipt,next_state,realized)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--arm',type=int,required=True);args=ap.parse_args();plan=read(args.plan)
    finite_sequence=bool(plan.get('finite_sequence_diagnostic',False))
    if finite_sequence and not plan.get('skip_policies',False):raise ValueError('Finite-horizon sequence diagnostic cannot launch policies')
    root=Path(plan['source_root']);out=Path(plan['output_root'])/f'arm_{args.arm}';out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();deadline=start+plan['total_seconds'];stop=threading.Event();state={'phase':'preflight'}
    def heartbeat():
        while not stop.wait(60):
            save(out/'heartbeat.json',dict(state,elapsed=time.monotonic()-start))
            if time.monotonic()>deadline:save(out/'failure.json',dict(error='Total budget exhausted',**state));os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    for path,digest in plan['file_sha256'].items():
        if sha(path)!=digest:raise ValueError('Pinned input/source changed: '+path)
    sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import run_e5f_candidate_terminal as td
    import e5f_successive_surprises as surprise
    import e5f_balanced_terminal as terminal_module
    import run_e5f_transition_calibration as fertility
    from e5f_approved_initial_state import build_approved_initial_state
    score=read(plan['initial_score_path'])
    if score['contract_sha256']!=plan['target_fingerprint'] or len(score['target_fit'])!=13:raise ValueError('Initial target/weight fingerprint mismatch')
    if not any(r.get('target')==.7202462623815278 for r in score['target_fit']):raise ValueError('Retained rooms target missing')
    c=plan['terminal_template'];td.validate_contract(c);td.verify_sources(c)
    for name in ('initial_checkpoint','initial_summary','initial_contract'):td.verify(c[name]['path'],c[name]['sha256'])
    with gzip.open(c['initial_checkpoint']['path'],'rb') as f:packet=pickle.load(f)
    initial_summary=read(c['initial_summary']['path']);initial_contract=read(c['initial_contract']['path'])
    td.validate_initial_receipts(c,packet,initial_summary,initial_contract)
    old=build_approved_initial_state(packet=packet,normalization=initial_summary['normalization'],outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.,fertility_tolerance=5e-4)
    demographics=packet['demographic_seed'];initial_psi=float(old.parameters.psi_child)
    inherited=surprise.InheritedState(2007,old.initial_state);del packet
    if plan.get('stationary_restart_2019'):
        spec=plan['stationary_restart_2019'];static_plan=read(spec['fit_plan']);selected=read(spec['fit_summary'])
        if (static_plan['terminal_template']['initial_checkpoint']!=c['initial_checkpoint'] or
            selected['status']!='stationary_fertility_target_matched' or selected['decision_year']!=2015):
            raise ValueError('2019 stationary restart must use the same pinned structural seed and the window ending2019')
        point=selected['selected'];td.verify(point['checkpoint'],point['checkpoint_sha256'])
        with gzip.open(point['checkpoint'],'rb') as f:stationary=pickle.load(f)
        Q=stationary['parameters'];sol=stationary['solution']
        if not np.array_equal(stationary['b_grid'],old.b_grid):raise ValueError('Stationary restart grid changed')
        from run_e5f_matched_pf_history import pf
        births=fertility.closure.topcode_consistent_renewal_accounting(sol,Q)
        adjusted=float(births['topcode_adjusted_birth_children'])/2.1;raw=float(sol.total_births_kfe)/2.1
        start_state=pf.PFInitialState(stationary['stationary_g_pre'].copy(),[adjusted]*4,[raw]*4)
        inherited=surprise.InheritedState(2019,start_state)
        save(out/'restart_2019.json',dict(source=spec,checkpoint_sha256=point['checkpoint_sha256'],psi_2019=float(Q.psi_child),
            household_mass=float(start_state.g_pre.sum()),raw_birth_queue=start_state.scheduled_raw_entries,
            adjusted_birth_queue=start_state.scheduled_entries,distribution_reset_2023=False,
            interpretation='Conditional stationary household distribution dated2019; birth queues from its constant birth flows. Original supply curve and external2023person anchor retained; no claim of full demographic stationarity.'))
    targets=list(csv.DictReader(Path(plan['empirical_blocks']).open()));realized=[];smoked=False;trial_index=0;warm=None;warm_year=None
    if plan.get('stationary_restart_2019'):targets=[t for t in targets if int(t['decision_year'])>=2019]
    if plan.get('verified_native_smoke'):
        native=read(plan['verified_native_smoke'])
        if (not native['finite_horizon_market_fiscal_converged'] or native['start_year']!=2007
                or abs(native['psi']-(initial_psi+plan['seed_steps'][args.arm]))>1e-14):
            raise ValueError('Pinned native smoke does not match this preference and initial date')
        smoked=True
    audit=terminal_module.TerminalAuditControls(**c['audit_controls'])
    fit_deadline=deadline-plan['policy_reserve_seconds']
    try:
        for target in targets:
            year=int(target['decision_year']);desired=float(target['period_tfr_arithmetic_mean']);trials=[];winner=None;seen=[]
            if inherited.year!=year:raise RuntimeError('Inherited shock clock differs from empirical window')
            for attempt in range(plan['maximum_trials_per_window']):
                if time.monotonic()+1800>fit_deadline:raise TimeoutError('Historical fitting budget reached; policy reserve preserved')
                center=initial_psi if not realized else realized[-1]['psi']
                seed_step=plan.get('seed_steps_by_year',{}).get(str(year),plan['seed_steps'][args.arm])
                proposal_step=plan.get('proposal_step',seed_step) if trials else seed_step
                psi=next_psi(trials,center,proposal_step,(initial_psi-.20,initial_psi+.02))
                if not trials and not seen and str(year) in plan.get('initial_psi_by_year',{}):psi=float(plan['initial_psi_by_year'][str(year)])
                if any(abs(psi-v)<1e-6 for v in seen):
                    psi=float(np.clip(center+(attempt+1)*plan['seed_steps'][args.arm],initial_psi-.20,initial_psi+.02))
                    if any(abs(psi-v)<1e-6 for v in seen):break
                seen.append(psi);folder=out/f'trial_{trial_index:02d}_{year}';trial_index+=1;folder.mkdir()
                state.update(phase='terminal',year=year,psi=psi,trial=trial_index);save(out/'latest_stage.json',state)
                tc=copy.deepcopy(c);tc['psi_change_from_initial']=psi-initial_psi
                tc.update(plan.get('terminal_contract_overrides',{}));save(folder/'terminal_contract.json',tc)
                terminal_driver=Path(plan.get('terminal_driver',root/'code/model/tools/run_e5f_candidate_terminal.py'))
                with (folder/'terminal.log').open('w') as log:
                    status=subprocess.run([sys.executable,'-B',str(terminal_driver),'--contract',str(folder/'terminal_contract.json'),'--contract-sha256',sha(folder/'terminal_contract.json'),'--output',str(folder/'terminal')],stdout=log,stderr=subprocess.STDOUT,timeout=1860)
                if status.returncode:
                    save(folder/'rejected.json',dict(stage='terminal',exit=status.returncode));continue
                with gzip.open(folder/'terminal/terminal_state.pkl.gz','rb') as f:t=pickle.load(f)
                tr=read(folder/'terminal/root_receipt.json');payload=tr['final']['payload']
                terminal=terminal_module.BalancedTerminalEndpoint(t['parameters'],t['b_grid'],t['policy'],t['endpoint'],t['social_security'],payload['diagnostics'],payload['household_gates'])
                pstart=float(plan['warm_price_2007']);bstart=float(plan['warm_pension_2007']);previous=None;result=None;observed={}
                counts=plan.get('forecast_counts',([6] if not smoked else [])+[28,56])
                if counts!=[6] and not smoked and counts[0]!=6:raise ValueError('Long diagnostic requires its pinned native smoke first')
                for count in counts:
                    if time.monotonic()+900>fit_deadline:break
                    prices=np.linspace(pstart,float(terminal.policy.price[0]),count);pensions=np.linspace(bstart,float(terminal.parameters.pension),count)
                    if plan.get('initialization_receipt'):
                        init=read(plan['initialization_receipt']);initial_coordinates=init.get('final') or init['best'];f=initial_coordinates;n=min(count,len(f['prices']))
                        prices[:n]=np.asarray(f['prices'][:n])*plan.get('initial_price_multiplier',1.)
                        pensions[:n]=f['fiscal_values'][:n]
                        if count>n:
                            prices[n:]=np.linspace(prices[n-1],float(terminal.policy.price[0]),count-n+1)[1:]
                            pensions[n:]=np.linspace(pensions[n-1],float(terminal.parameters.pension),count-n+1)[1:]
                        if plan.get('verified_native_smoke') and plan.get('initialization_native_prefix',True):
                            q=read(plan['verified_native_smoke'])['final'];k=min(count,len(q['prices']))
                            prices[:k]=q['prices'][:k];pensions[:k]=q['fiscal_values'][:k]
                    if finite_sequence and warm is not None:
                        f=warm.get('final') or warm['best'];offset=(year-warm_year)//4
                        for arr,key,tail in ((prices,'prices',float(terminal.policy.price[0])),(pensions,'fiscal_values',float(terminal.parameters.pension))):
                            source=np.asarray(f[key])[offset:];take=min(count,len(source));arr[:take]=source[:take]
                            if take<count:arr[take:]=np.linspace(arr[take-1],tail,count-take+1)[1:]
                    if previous is not None:
                        f=previous.root_receipt.get('final') or previous.root_receipt.get('best')
                        if f is not None:
                            take=min(count,len(f['prices']));prices[:take]=f['prices'][:take];pensions[:take]=f['fiscal_values'][:take]
                    for continuation in range(3):
                        rc=copy.deepcopy(plan['history_root_controls']);rc['initial_jacobian']=None
                        if continuation==0 and plan.get('initialization_receipt') and len(initial_coordinates['prices'])==count:
                            rc['initial_jacobian']=init.get('final_jacobian')
                        if finite_sequence and warm is not None:
                            rc['initial_jacobian']=warm.get('final_jacobian') if warm_year==year else None
                        if continuation and result is not None:
                            f=result.root_receipt.get('final') or result.root_receipt.get('best')
                            if f is None:break
                            prices=np.asarray(f['prices']);pensions=np.asarray(f['fiscal_values']);rc['initial_jacobian']=result.root_receipt.get('final_jacobian')
                        stage=folder/f'forecast_{count}_{continuation}';stage.mkdir();state.update(phase='forecast',dates=count,continuation=continuation);save(out/'latest_stage.json',state)
                        snapshot={};measurements=[];dated_allocation={}
                        def observe(i,e,P,grid,shared):
                            if i==0:measurements.clear();dated_allocation.clear();snapshot.clear();snapshot.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=old.supply_rule)
                            if year+4*i==2023:
                                dated_allocation.update(calendar_year=2023,large_owner_age_cells=large_owner_observation(e,P),interpretation='Expected2023allocation within this forecast vintage, not automatically realized history')
                            measurements.append(dict(calendar_year=year+4*i,**fertility.period_fertility_diagnostics(e,P)))
                            save(stage/'latest_date.json',dict(year=year+4*i));state['last_completed_year']=year+4*i
                        def progress(record):
                            save(stage/'latest_completed.json',record)
                            if record.get('new_best'):save(stage/'best_so_far.json',record)
                        if continuation==0 and plan.get('initialization_native_prefix') is False and not finite_sequence:
                            if (len(initial_coordinates['prices'])!=count or
                                not np.array_equal(prices,np.asarray(initial_coordinates['prices'])) or
                                not np.array_equal(pensions,np.asarray(initial_coordinates['fiscal_values']))):
                                raise ValueError('Continuation must preserve the complete pinned path coordinates')
                            save(stage/'initialization_verified.json',dict(receipt=plan['initialization_receipt'],sha256=sha(plan['initialization_receipt']),full_path_preserved=True,jacobian_reused=rc['initial_jacobian'] is not None))
                        result=surprise.solve_surprise(inherited=inherited,psi=psi,old_state=old,terminal=terminal,terminal_root_receipt=tr,
                            demographic_primitives=demographics,terminal_demographic_primitives=t['demographic_seed'],count=count,initial_prices=prices,initial_pensions=pensions,
                            audit_controls=audit,root_controls=rc,deadline_monotonic=min(fit_deadline,time.monotonic()+(1800 if count==6 else 7200)),pension_tail_tolerance=.01,callback=progress,observer=observe)
                        result=carry_finite_diagnostic(result,enabled=finite_sequence,inherited=inherited,old=old,demographics=demographics,psi=psi,module=surprise)
                        if dated_allocation:save(stage/'allocation_2023.json',dict(dated_allocation,finite_converged=result.root_receipt['finite_horizon_market_fiscal_converged'],horizon_verified=False,forecast_vintage_year=year))
                        surprise.persist_episode(stage/'vintage',year,result,provenance={'plan_sha256':sha(args.plan),'terminal_contract_sha256':sha(folder/'terminal_contract.json'),'target_fingerprint':plan['target_fingerprint']})
                        save(stage/'fertility.json',measurements);previous=result
                        if finite_sequence and (result.root_receipt.get('final') or result.root_receipt.get('best')):
                            warm=result.root_receipt;warm_year=year
                        if plan.get('forecast_diagnostic_only',False) and result.root_receipt['finite_horizon_market_fiscal_converged']:
                            standard_graphs(snapshot,result,out/'native_graphs')
                            save(out/'summary.json',dict(status='finite_forecast_diagnostic_only',dates=count,year=year,psi=psi,data=desired,
                                model=measurements[0]['period_tfr_topcode_adjusted'],gap=measurements[0]['period_tfr_topcode_adjusted']-desired,
                                finite_converged=True,terminal_distance_passed=result.root_receipt['terminal_distance_passed'],
                                historical_fit_complete=False,horizon_verified=False,production_eligible=False,forecast_folder=str(stage)))
                            return
                        if count==6:
                            if not result.root_receipt['finite_horizon_market_fiscal_converged']:
                                if continuation==2:raise RuntimeError('Native exact-loop smoke failed after bounded continuations; no long new-timing solve')
                                continue
                            smoked=True;save(out/'native_smoke.json',dict(passed=True,year=year,psi=psi,folder=str(stage)))
                            if plan.get('native_smoke_only',False):
                                standard_graphs(snapshot,result,out/'native_graphs')
                                save(out/'summary.json',dict(status='native_smoke_only',year=year,psi=psi,data=desired,
                                    model=measurements[0]['period_tfr_topcode_adjusted'],
                                    gap=measurements[0]['period_tfr_topcode_adjusted']-desired,
                                    finite_converged=True,historical_fit_complete=False,horizon_verified=False,production_eligible=False,
                                    forecast_folder=str(stage),measurement_approximation_retained=True))
                                return
                            break
                        if result.next_state is not None:break
                    if result is not None and result.next_state is not None:break
                if result is None or result.next_state is None:
                    save(folder/'rejected.json',dict(stage='forecast_or_terminal_distance'));continue
                standard_graphs(snapshot,result,folder/'accepted_graphs')
                value=measurements[0]['period_tfr_topcode_adjusted'];error=value-desired;trials.append((psi,error))
                save(folder/'fit.json',dict(finite_horizon_diagnostic=finite_sequence,terminal_distance_passed=result.root_receipt['terminal_distance_passed'],year=year,psi=psi,data=desired,model=value,gap=error,measurement='retained household-rate analogue of published femaleTFR'))
                if winner is None or abs(error)<winner['error_abs']:
                    winner=dict(error_abs=abs(error),psi=psi,model=value,data=desired,year=year,folder=str(stage),finite_horizon_diagnostic=finite_sequence,terminal_distance_passed=result.root_receipt['terminal_distance_passed']);winner_state=result.next_state;winner_receipt=result.root_receipt
                    save(out/'best_so_far.json',dict(realized=realized,current=winner))
                    with gzip.open(folder/'first_period_diagnostics.pkl.gz','wb',compresslevel=1) as f:pickle.dump(snapshot,f,protocol=5)
                if abs(error)<=plan['fertility_fit_tolerance']:break
            if winner is None:raise RuntimeError('No accepted full forecast for this inherited state')
            if winner['error_abs']>plan['fertility_fit_tolerance']:raise RuntimeError('Shock fit did not meet declared tolerance; do not carry inaccurate fit forward')
            realized.append(winner);inherited=winner_state;warm=winner_receipt;warm_year=year;save(out/'realized_fit.json',realized)
            with gzip.open(out/f'realized_state_{inherited.year}.pkl.gz','wb',compresslevel=1) as f:pickle.dump(inherited,f,protocol=5)
        if finite_sequence:
            save(out/'finite_horizon_sequence_fit.json',dict(realized=realized,historical_fit_complete=False,finite_horizon_sequence_fit_complete=True,horizon_verified=False,production_eligible=False,forecast_counts=plan['forecast_counts']))
            save(out/'summary.json',dict(status='finite_horizon_sequence_fit_complete',realized=realized,horizon_verified=False,production_eligible=False,policies_launched=False))
            return
        save(out/'historical_fit_complete.json',dict(realized=realized,measurement_approximation_retained=True,horizon_verified=False))
        state['phase']='policy'
        from run_e5f_successive_surprise_policy import run_policies
        run_policies(inherited=inherited,psi=realized[-1]['psi'],old=old,demographics=demographics,plan=plan,out=out/'policies',deadline=deadline)
        save(out/'summary.json',dict(status='historical_fit_and_policy_attempts_complete',realized=realized,production_eligible=False))
    except Exception as exc:
        save(out/'failure.json',dict(state,error_type=type(exc).__name__,error=str(exc),realized=realized));raise
    finally:stop.set()
if __name__=='__main__':main()
