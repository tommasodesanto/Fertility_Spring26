"""Reproduce a rejected first root step at its boundary only, without a history solve."""
import argparse,copy,gzip,hashlib,json,os,pickle,sys,time
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[key]='1'
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--case-dir',type=Path,required=True);ap.add_argument('--helper',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    read=lambda p:json.loads(Path(p).read_text())
    manifest=read(args.manifest);plan=read(manifest['prior_plan']);root=Path(plan['source_root'])
    sys.path[:0]=[str(args.helper),str(root/'code/model/tools'),str(root/'code/model')]
    from run_e5f_final_rebated_history import verify_pins,save,initial_seed_step
    import e5f_rebated_surprises as rebated
    import e5f_closed_finite_boundary as closed
    from e5f_rebated_initial_bridge import build_rebated_initial_state
    from e5f_balanced_terminal import TerminalAuditControls
    from e5f_matched_pf_path_root import solve_price_path
    from run_e5f_successive_surprises_overnight import next_psi
    verify_pins(manifest['file_sha256']);verify_pins(plan['file_sha256'])
    _,joined,primitive,_,rent_domain=rebated._runtime()
    summary=read(manifest['initial_summary']);item=summary['checkpoint']
    cp=Path(item.get('path',item.get('checkpoint')))
    digest=item.get('sha256',item.get('checkpoint_sha256'))
    contract=read(args.case_dir/'contract_receipt.json')
    assert sha(cp)==digest==contract['initial_checkpoint_sha256']
    assert sha(args.manifest)==contract['manifest_sha256']
    with gzip.open(cp,'rb') as f:packet=pickle.load(f)
    raw=read(cp.parent/'summary.json')
    old=build_rebated_initial_state(packet=packet,normalization=raw['normalization'],
        outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.)
    count=contract['count'];width=count+1
    saved_path=args.case_dir/'window_2007/trial_00/latest_completed.json';saved=read(saved_path)
    assert saved['evaluation']==1 and saved['mapping_valid'] is True
    p0=np.asarray(saved['prices']);residual=np.asarray(saved['residual'])
    controls=dict(plan['history_root_controls']);controls.update(manifest['root_controls'])
    slope=float(controls['market_slope']);bounds=[controls[k] for k in ('price_bounds','pension_bounds','transfer_bounds')]
    psi=next_psi([],float(old.parameters.psi_child),initial_seed_step(plan,manifest),
        (float(old.parameters.psi_child)-.20,float(old.parameters.psi_child)+.02))
    def project(values):
        x=np.asarray(values).reshape(3,width).copy()
        for j,(lo,hi) in enumerate(bounds):x[j]=np.clip(x[j],lo,hi)
        from types import SimpleNamespace
        endpoint=SimpleNamespace(parameters=old.parameters,asset_price=float(x[0,-1]))
        x[0,:-1]=rent_domain.project_price_path_to_positive_rents(x[0,:-1],terminal=endpoint,minimum_rent_share=1e-6)[0]
        if np.any(x[0]>bounds[0][1]):raise ValueError('Projected price out of bounds')
        return x.ravel()
    captured=[]
    class Captured(BaseException):pass
    def replay_first_then_capture(prices):
        if not captured:
            np.testing.assert_array_equal(prices,p0);captured.append(prices.copy())
            return dict(residual=residual,mapping_valid=True)
        captured.append(prices.copy());raise Captured()
    rc={k:controls[k] for k in ('market_tolerance','max_log_step','damping','max_evaluations','max_condition_number','worsening_factor','final_reproduction_tolerance')}
    try:
        solve_price_path(initial_prices=p0,evaluate=replay_first_then_capture,project=project,slope=slope,
            default_jacobian=np.diag(np.r_[np.full(width,-slope),np.full(2*width,-200.)]),
            deadline_monotonic=time.monotonic()+60,**rc)
    except Captured:pass
    assert len(captured)==2
    step=np.clip(controls['damping']*residual/np.r_[np.full(width,slope),np.full(2*width,200.)],-controls['max_log_step'],controls['max_log_step'])
    np.testing.assert_array_equal(captured[1],project(np.exp(np.log(p0)+step)))
    args.out.mkdir(parents=True,exist_ok=False)
    inputs=[args.manifest,Path(manifest['prior_plan']),cp,args.case_dir/'contract_receipt.json',saved_path,Path(__file__)]
    save(args.out/'inputs.json',dict(source_sha256={str(p):sha(p) for p in inputs},psi=psi,count=count,
        initial_coordinates=p0,failed_proposal=captured[1],boundary_only=True,no_history_solved=True))
    model=primitive.model;native=primitive.dated_budget;results=[];current={}
    def budget(evaluation,P,shared,grid,rent):
        try:return native(evaluation,P,shared,grid,rent)
        except RuntimeError as exc:
            if not str(exc).startswith('Dated budget gate failed:'):raise
            rows=[];p=evaluation.policy;g=evaluation.g_current
            for age in range(P.J):
                for tenure in range(g.shape[1]):
                    for zz,z in enumerate(P.z_grid):
                        for n in range(P.n_parity):
                            for child in range(P.n_child_states):
                                idx=(slice(None),tenure,0,age,zz,n,child);mass=g[idx]
                                if mass.sum()<=0:continue
                                income=model.income_at_state(P,0,age,float(z));resources=P.R_gross*grid+income
                                grant=float(shared.gb_flat.reshape(-1)[n+P.n_parity*child])
                                resources+=np.clip(grant-(P.R_gross*np.maximum(grid,0)+income),0,grant)
                                if tenure==0:cost=rent*p.hR_pol[idx]
                                else:
                                    h=P.H_own[tenure-1];cost=(P.delta+P.tau_H)*p.price[0]*h
                                    cost+=getattr(P,'owner_size_cost',0)*p.price[0]*max(h-getattr(P,'owner_size_cost_ref',6),0)**getattr(P,'owner_size_cost_power',2)
                                gap=p.c_pol[idx]+cost+p.bp_pol[idx]-resources
                                for b in np.flatnonzero((gap>1e-9)&(mass>0)):
                                    rows.append(dict(wealth_index=int(b),wealth=grid[b],age_index=age,tenure=tenure,
                                        income_state=zz,children_ever_born=n,child_state=child,mass=mass[b],gap=gap[b],
                                        consumption=p.c_pol[idx][b],next_wealth=p.bp_pol[idx][b],resources=resources[b],
                                        housing_cost=float(cost[b]) if np.ndim(cost) else float(cost)))
            current['violations']=dict(count=len(rows),total_mass=sum(r['mass'] for r in rows),
                largest=sorted(rows,key=lambda r:r['mass'],reverse=True)[:30],error=str(exc))
            raise
    primitive.dated_budget=budget
    try:
        half_step=np.clip(.5*controls['damping']*residual/np.r_[np.full(width,slope),np.full(2*width,200.)],
            -controls['max_log_step'],controls['max_log_step'])
        half_coordinates=project(np.exp(np.log(p0)+half_step))
        for label,coordinates in [('initial',p0),('full_step',captured[1]),
                ('half_step',half_coordinates),
                ('full_step_repeat',captured[1])]:
            current=dict(label=label,started_unix=time.time());blocks=coordinates.reshape(3,width)
            Q=copy.deepcopy(old.parameters);Q.psi_child=float(psi)
            current['boundary_coordinates']=[float(x[-1]) for x in blocks]
            try:
                value=closed.boundary_evaluation(parameters=Q,g_pre=old.initial_state.g_pre,grid=old.b_grid,
                    supply_rule=old.supply_rule,price=blocks[0,-1],pension=blocks[1,-1],transfer=blocks[2,-1],
                    audit_controls=TerminalAuditControls(**plan['terminal_template']['audit_controls']),
                    deadline_monotonic=time.monotonic()+600)
                current.update(status='boundary_household_gates_pass',gates=value.gates,residuals=value.residuals)
            except Exception as exc:current.update(status='rejected',error_type=type(exc).__name__,error=str(exc))
            current['seconds']=time.time()-current['started_unix'];results.append(current)
            save(args.out/'latest_completed.json',current);save(args.out/'summary.json',dict(results=results,boundary_only=True,no_history_solved=True))
            print(json.dumps(dict(label=label,status=current['status'],seconds=current['seconds'])),flush=True)
    finally:primitive.dated_budget=native
    verify_pins(manifest['file_sha256'])

if __name__=='__main__':main()
