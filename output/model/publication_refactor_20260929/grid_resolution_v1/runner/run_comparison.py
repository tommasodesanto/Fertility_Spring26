"""Pinned two-arm full renewal-price/population GE; no recalibration."""
from __future__ import annotations
import argparse, copy, csv, hashlib, json, os, subprocess, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
ROOT=PACKET.parents[3]
MATCHED=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
sys.path[:0]=[str(MATCHED),str(MATCHED/'source'),str(ROOT/'code/model')]
import driver as authored
import single_price
import phase_b_grid as ge
from refactor_lab.inputs import decode, serialized, sha256_file
from refactor_lab.engine.utils import make_grid
write=authored.write_json
D=.14
CHANGED={'Nb','Nz','z_grid','z_weights','Pi_z','earnings_transaction_grid','fixed_reference_entry_grid','fixed_reference_entry_conditional'}
REMOVED={'permanent_income_group_index','permanent_income_base_state_index'}

def require(ok,message):
    if not ok: raise RuntimeError(message)

def verify_sources():
    pins=json.loads((HERE/'source_hashes.json').read_text())
    for rel,digest in pins.items():
        require(sha256_file(ROOT/rel)==digest,'Runner source drift: '+rel)
    prep=json.loads((PACKET/'preflight.json').read_text())
    for rel,digest in prep['source_pins'].items():
        require(sha256_file(ROOT/rel)==digest,'Preparation source drift: '+rel)
    require(sha256_file(PACKET/'proposed_120x9/bundle.json')==prep['proposed_bundle_sha256'],'Proposal bundle drift')
    return prep

def apply_proposal(context,prep):
    original=context['P']; folder=PACKET/'proposed_120x9'
    meta=json.loads((folder/'bundle.json').read_text())
    require(meta['parent_bundle_sha256']==authored.BUNDLE_SHA,'Proposal parent drift')
    require(sha256_file(folder/'arrays.npz')==meta['arrays_sha256'],'Proposal arrays drift')
    with np.load(folder/'arrays.npz',allow_pickle=False) as arrays:
        fields={k:decode(v,arrays) for k,v in meta['parameters'].items()}
        grid=arrays['b_grid'].copy()
    require(set(vars(original))-set(fields)==REMOVED,'Proposal removed fields differ')
    require(not set(fields)-set(vars(original)),'Unexpected proposal fields')
    for key in set(fields)-CHANGED:
        require(serialized(fields[key])==serialized(getattr(original,key)),'Economic input drift: '+key)
    P=SimpleNamespace(**fields)
    require(P.Nb==120 and P.Nz==9 and P.fixed_reference_entry_conditional.shape==(120,9),'Proposal dimensions differ')
    np.testing.assert_array_equal(make_grid(P),grid)
    oldmass=np.einsum("ij,j->i",original.fixed_reference_entry_conditional,original.z_weights,optimize=False)
    newmass=np.einsum("ij,j->i",P.fixed_reference_entry_conditional,P.z_weights,optimize=False)
    restored=np.zeros(len(context['b_grid']))
    lookup={float(v):i for i,v in enumerate(context['b_grid'])}
    for v,m in zip(grid,newmass): restored[lookup[float(v)]]=m
    np.testing.assert_allclose(restored,oldmass,rtol=0,atol=3e-16)
    require(not P.permanent_income_levels_enabled,'Dormant income groups activated')
    context.update(P=P,b_grid=grid,proposal_bundle_sha256=prep['proposed_bundle_sha256'])

class ArmBudget(authored.Budget):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.max_lifecycle=6


def policy_snapshot(live):
    values={'b_grid':live['b_grid'],'z_grid':live['P'].z_grid}
    for key in ('V','bp_pol','c_pol','hR_pol'):
        value=getattr(live['sol'],key,None)
        require(isinstance(value,np.ndarray),'Missing required policy snapshot: '+key)
        values[key]=value
    tp=live['sol'].tenure_probs
    require(isinstance(tp,np.ndarray) and tp.ndim==8 and tp.shape[:-1]==values['V'].shape,'Tenure choice probability layout differs')
    require(np.isfinite(tp).all() and np.min(tp)>=0 and np.max(tp)<=1,'Tenure choice probability bounds fail')
    values['owner_choice_probability']=tp[...,1:].sum(axis=-1)
    return values


def arm(args,name,prep,deadline,mock=False):
    out=args.out/name;out.mkdir(parents=True)
    context=authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=out))
    context['deadline_epoch']=deadline
    if not mock: authored.authenticate_frozen(context) # authenticate original BEFORE grid replacement
    if name=='proposal_120x9':apply_proposal(context,prep)
    dims={'wealth_grid_nodes':len(context['b_grid']),'income_states':len(context['P'].z_grid)}
    context['expected_dimensions']=dims
    context['selected_d_bar']=D;context['reference_psi']=float(context['P'].psi_child)
    write(out/'effective_input_contract.json',dict(reference=context['loaded'].identity,dimensions=dims,unsecured_credit_limit=D,credit_status='previously_authorized_experiment_not_adopted',psi_fixed=context['reference_psi'],proposal_bundle_sha256=context.get('proposal_bundle_sha256'),projection_disclosure=prep['disclosure'],entry_moments=prep['entry_moments'],closure='price clears actual birth renewal; population clears fixed physical housing supply',recalibration=False))
    budget=ArmBudget(out,deadline,smoke=mock)
    native_observe=ge.observe_price
    def observe_snapshot(ctx,live,label,*,final=False):
        receipt=native_observe(ctx,live,label,final=final)
        if final:
            values=policy_snapshot(live)
            np.savez_compressed(out/'phase_b_ge'/label/'common_support_policies.npz',**values)
        return receipt
    if not mock:ge.observe_price=observe_snapshot
    if mock:
        context.update(fp=SimpleNamespace(write=write),prepared=object(),manifest={},objective={},runtime=object(),reference={})
        (out/'phase_b_ge').mkdir()
        result=ge.smoke_phase_b(context,budget)
    else:
        seed=single_price.solve_fixed_price(context,D,context['q_ref'],budget,'phase_a_fixed_d14',out/'seed')
        result=ge.run_phase_b(context,dict(selected_d_bar=D,selected_live=seed),budget)
        require(result['status']=='passed','Arm GE not certified within budget')
    ge.observe_price=native_observe
    write(out/'completed.json',dict(result=result,lifecycle_solves=budget.used_lifecycle,dimensions=dims))
    return result


def compare(out):
    folders=[out/n/'phase_b_ge/selected_root' for n in ('control_160x15','proposal_120x9')]
    for filename,key in [('target_fit.csv','moment'),('parameters.csv','parameter')]:
        rows=[]
        for p in folders:
            with (p/filename).open(newline='') as stream: rows.append(list(csv.DictReader(stream)))
        require(len(rows[0])==len(rows[1])==(14 if filename.startswith('target') else 31),'Table row mismatch')
        require([r[key] for r in rows[0]]==[r[key] for r in rows[1]],'Table keys differ')
        write(out/('comparison_'+filename+'.json'),[dict(control=a,proposal=b) for a,b in zip(*rows)])
    policy_started=time.monotonic();policy_comparison={}
    with np.load(folders[0]/'common_support_policies.npz') as a, np.load(folders[1]/'common_support_policies.npz') as b:
        wealth=b['b_grid'];income=b['z_grid'];oldincome=a['z_grid']
        index=np.searchsorted(a['b_grid'],wealth)
        np.testing.assert_array_equal(a['b_grid'][index],wealth)
        require(income[0]>=oldincome[0] and income[-1]<=oldincome[-1],'Income supports differ')
        high=np.clip(np.searchsorted(oldincome,income),1,len(oldincome)-1);low=high-1
        weight=(income-oldincome[low])/(oldincome[high]-oldincome[low])
        w=weight.reshape((1,1,1,1,len(income),1,1))
        vlo=a['V'][index][:,:,:,:,low];vhi=a['V'][index][:,:,:,:,high];newV=b['V']
        require(newV.ndim==7 and vlo.shape==newV.shape,'Value state layout differs')
        child_ok=np.arange(newV.shape[6])[None,:]<=np.arange(newV.shape[5])[:,None]
        child_ok=child_ok.reshape((1,1,1,1,1)+child_ok.shape)
        old_feasible=(vlo>-1e9)&(vhi>-1e9)&np.isfinite(vlo)&np.isfinite(vhi)&child_ok
        new_feasible=(newV>-1e9)&np.isfinite(newV)&child_ok
        common=old_feasible&new_feasible
        policy_comparison['feasibility']=dict(common_state_count=int(common.sum()),old_only_state_count=int((old_feasible&~new_feasible).sum()),new_only_state_count=int((new_feasible&~old_feasible).sum()),invalid_child_states_excluded=True,dead_value_cutoff=-1e9,interpolation_stencil_requires_both_old_income_states_feasible=True)
        for key in ('V','bp_pol','c_pol','hR_pol','owner_choice_probability'):
            av,bv=a[key],b[key]
            require(av.ndim==7 and av.shape[0]==160 and av.shape[4]==15 and bv.shape==newV.shape,'Policy state layout differs: '+key)
            mapped=(1-w)*av[index][:,:,:,:,low]+w*av[index][:,:,:,:,high]
            gap=bv-mapped
            require(np.isfinite(gap[common]).all(),'Nonfinite mapped feasible policy')
            def maximum(values,mask):
                return float(abs(values[mask]).max()) if np.any(mask) else None
            policy_comparison[key]=dict(status='subset_wealth_linear_income_common_feasible_support',max_abs=maximum(gap,common),mean_abs=float(abs(gap[common]).mean()) if np.any(common) else None,poorest_max_abs=maximum(gap[0],common[0]),richest_max_abs=maximum(gap[-1],common[-1]),youngest_max_abs=maximum(gap[:,:,:,0],common[:,:,:,0]),oldest_max_abs=maximum(gap[:,:,:,-1],common[:,:,:,-1]),renter_max_abs=maximum(gap[:,0],common[:,0]),owner_max_abs=maximum(gap[:,1:],common[:,1:]),no_children_max_abs=maximum(gap[:,:,:,:,:,0],common[:,:,:,:,:,0]),high_children_max_abs=maximum(gap[:,:,:,:,:,-1],common[:,:,:,:,:,-1]),interpretation='Interpolation diagnostic on mutually feasible states, not exact equality or occupied-mass accuracy')
    policy_comparison['workflow_seconds']=time.monotonic()-policy_started
    write(out/'comparison_common_support_policies.json',policy_comparison)
    closures=[json.loads((p/'closure.json').read_text()) for p in folders]
    write(out/'comparison_closure.json',dict(control=closures[0],proposal=closures[1],difference={k:closures[1][k]-closures[0][k] for k in closures[0] if isinstance(closures[0][k],(float,int)) and isinstance(closures[1].get(k),(float,int))}))


def negative_parameter_preflight():
    with (PACKET/'reference_parameters.csv').open(newline='') as stream: rows=list(csv.DictReader(stream))
    actual={r['parameter']:float(r['estimate']) for r in rows}
    attempts=0
    for dims in ({'wealth_grid_nodes':160,'income_states':15},{'wealth_grid_nodes':120,'income_states':9}):
        ctx={'expected_dimensions':dims};correct=dict(actual,**dims)
        ge.validate_parameter_estimates(ctx,rows,correct)
        for row in rows:
            wrong=dict(correct);key=row['parameter'];wrong[key]+=1
            try:ge.validate_parameter_estimates(ctx,rows,wrong)
            except RuntimeError:attempts+=1
            else:raise RuntimeError('Parameter drift accepted: '+key)
        stale=dict(correct,wealth_grid_nodes=120 if dims['wealth_grid_nodes']==160 else 160)
        try:ge.validate_parameter_estimates(ctx,rows,stale)
        except RuntimeError:attempts+=1
        else:raise RuntimeError('Stale reference/proposal dimensions accepted')
    return dict(valid_dimension_arms=2,rejected_parameter_drifts=attempts,lifecycle_solves=0)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['preflight','full','arm']);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float);ap.add_argument('--mock-child',action='store_true');ap.add_argument('--arm-name',choices=['control_160x15','proposal_120x9'])
    args=ap.parse_args()
    if args.mode!='arm':
        require(not args.out.exists(),'Refusing existing output');args.out.mkdir(parents=True)
    else:require(args.arm_name and args.out.is_dir(),'Internal arm requires existing parent output')
    if args.mode=='full' or (args.mode=='arm' and not args.mock_child):
        require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Paired production batch requires Torch Slurm')
        require(int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'One CPU per task required')
    started=time.time();deadline=min(args.deadline_epoch or started+2400,started+2400)
    try:
        prep=verify_sources()
        parameter_negative_tests=negative_parameter_preflight()
        write(args.out/'parameter_gate_preflight.json',parameter_negative_tests)
        if args.mode=='arm':
            arm(args,args.arm_name,prep,deadline,mock=args.mock_child)
            return
        results={}
        for name in ('control_160x15','proposal_120x9'):
            remaining=deadline-time.time();require(remaining>0,'Global deadline reached')
            command=[sys.executable,str(__file__),'arm','--arm-name',name,'--out',str(args.out),'--deadline-epoch',str(deadline)]
            if args.mode=='preflight':command.append('--mock-child')
            child_env=dict(os.environ,NUMBA_CACHE_DIR=str(args.out/name/'numba_cache'))
            arm_started=time.monotonic()
            process=subprocess.run(command,timeout=remaining,env=child_env)
            require(process.returncode==0,'Fresh-process arm failed: '+name)
            results[name]=json.loads((args.out/name/'completed.json').read_text())['result']
            write(args.out/name/'external_workflow_timing.json',dict(seconds=time.monotonic()-arm_started,includes='subprocess imports, authentication, compilation, solve, native observers, reporting, serialization',cold_cache_directory=child_env['NUMBA_CACHE_DIR']))
        if args.mode=='full':compare(args.out)
        write(args.out/'completed.json',dict(status='mock_exact_loop_zero_solves' if args.mode=='preflight' else 'full_passed',results=results,total_workflow_seconds=time.time()-started,maximum_lifecycle_solves=12,total_seconds=2400,lifecycle_solves=sum(json.loads((args.out/name/'completed.json').read_text())['lifecycle_solves'] for name in results),not_recalibrated=True,limitations='Mock preflight authenticates bundles and pins, not checkpoint observer runtime or numerical GE. Full mode authenticates observer before replacing inputs; reporting remains subject to unchanged gates.'))
    except BaseException as exc:
        write((args.out/args.arm_name/'failure.json') if args.mode=='arm' else (args.out/'failure.json'),dict(type=type(exc).__name__,message=str(exc),total_workflow_seconds=time.time()-started,no_auto_retry=True));raise
if __name__=='__main__':main()
