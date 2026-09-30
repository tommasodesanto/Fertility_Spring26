"""Isolated five-ratio/current-income entry comparison; no estimation or adoption."""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, os, subprocess, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASE=ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner'
sys.path.insert(0,str(BASE))
spec=importlib.util.spec_from_file_location('matched_grid_workflow',BASE/'run_comparison.py')
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
authored,ge,write=base.authored,base.ge,base.write
from small_credit_lab.engine.distribution import entry_wealth_grid_weights, entry_wealth_ratio_distribution
from small_credit_lab.engine.shared import annual_gross_income_at_state, income_at_state
from small_credit_lab.engine import solver
from refactor_lab.inputs import serialized, sha256_file
ARMS=('control_fixed_reference','candidate_five_ratios')
D=.53
require=base.require

def context(out):
    return authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=out))

def verify_sources():
    pins=json.loads((HERE/'source_pins.json').read_text())
    for rel,digest in pins.items(): require(sha256_file(ROOT/rel)==digest,'Source drift: '+rel)
    # This also authenticates the complete reused engine and reporting source inventory.
    for rel,digest in json.loads((BASE/'source_hashes.json').read_text()).items():
        require(sha256_file(ROOT/rel)==digest,'Matched source drift: '+rel)

def apply_entry(ctx,name):
    P,grid=ctx['P'],ctx['b_grid']; original=copy.deepcopy(P)
    require(P.Nb==160 and P.Nz==15 and P.native_fixed_reference_entry is True,'Retained dimensions/mode differ')
    ratios,wr=entry_wealth_ratio_distribution(P)
    require(len(ratios)==5 and np.isfinite(ratios).all() and (wr>0).all(),'Exactly five empirical ratio nodes required')
    income=np.asarray([annual_gross_income_at_state(P,0,0,float(z)) for z in P.z_grid])
    points=ratios[:,None]*income[None,:]
    draw_mass=wr[:,None]*np.asarray(P.z_weights)[None,:]
    proposal=copy.deepcopy(P); proposal.native_fixed_reference_entry=False
    C=np.zeros((160,15))
    for zi,z in enumerate(P.z_grid):
        ix,wt=entry_wealth_grid_weights(grid,proposal,i=0,j=0,z_value=float(z)); C[ix,zi]=wt
    np.testing.assert_allclose(C.sum(axis=0),1.,rtol=0,atol=2e-15)
    if name==ARMS[1]: P.fixed_reference_entry_conditional=C
    require(set(vars(P))==set(vars(original)),'Unexpected parameter field')
    for key in vars(original):
        if name==ARMS[1] and key=='fixed_reference_entry_conditional': continue
        require(serialized(getattr(P,key))==serialized(getattr(original,key)),'Unapproved economic/input change: '+key)
    # Verify the operative fixed-matrix path equals the native ratio path for candidate.
    effective=np.zeros_like(C)
    for zi,z in enumerate(P.z_grid):
        ix,wt=entry_wealth_grid_weights(grid,P,i=0,j=0,z_value=float(z)); effective[ix,zi]=wt
    np.testing.assert_array_equal(effective,P.fixed_reference_entry_conditional)
    if name==ARMS[1]: np.testing.assert_array_equal(effective,C)
    marginal=effective@P.z_weights
    clip=points<grid[0]; high=points>grid[-1]
    return dict(arm=name,ratio_nodes=ratios.tolist(),ratio_weights=wr.tolist(),ratio_source=P.entry_wealth_ratio_source,
        ratio_draw_independent_of_current_income=True if name==ARMS[1] else False,
        law='native ratio*current annual gross entrant income, projected to retained grid' if name==ARMS[1] else 'retained fixed_reference_entry_conditional',
        runtime_mode='native_fixed_reference_entry=True with verified conditional matrix',
        conditional_sha256=hashlib.sha256(effective.tobytes()).hexdigest(),
        mean_wealth=float(grid@marginal),mean_annual_entry_income=float(income@P.z_weights),
        candidate_projection=dict(lower_clipped_draw_mass=float(draw_mass[clip].sum()),upper_clipped_draw_mass=float(draw_mass[high].sum()),
            total_clipped_draw_mass=float(draw_mass[clip|high].sum()),raw_mean_wealth=float((points*draw_mass).sum()),
            projected_mean_wealth=float(grid@(C@P.z_weights)),minimum_raw_wealth=float(points.min()),grid_minimum=float(grid[0]),
            disclosure='Native finite-support approximation clips out-of-grid raw wealth before linear projection. This is distinct from prohibited forward feasibility relocation; raw debt support is not exactly preserved.'),
        changed_input_fields=['fixed_reference_entry_conditional'] if name==ARMS[1] else [])

def entry_preflight(out):
    rows={};arrays={}
    for name in ARMS:
        ctx=context(out);rows[name]=apply_entry(ctx,name);P=ctx['P'];grid=ctx['b_grid']
        sd=solver.precompute_shared(P,grid)
        require(float(sd.gb_flat[0,0])==0,'Childless entrant transfer floor changed')
        slack=[float(P.R_gross*grid[bi]+income_at_state(P,0,0,float(P.z_grid[zi]))+D-float(sd.cb_flat[0,0])-P.user_cost_rate*ctx['q_ref']*float(sd.hb_flat[0,0])) for bi,zi in zip(*np.nonzero(P.fixed_reference_entry_conditional>0))]
        require(min(slack)>1e-6,'Necessary entrant current budget fails at fixed D53')
        rows[name]['necessary_current_budget_minimum_slack']=min(slack)
        arrays[name]=P.fixed_reference_entry_conditional.copy()
        arrays['b_grid']=grid;arrays['z_grid']=P.z_grid;arrays['z_weights']=P.z_weights
    np.savez_compressed(out/'initial_entry_distributions.npz',**arrays)
    write(out/'entry_preflight.json',dict(arms=rows,lifecycle_solves=0,credit=D,dimensions=[160,15]))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(13,3.8))
    for name,label in zip(ARMS,('Retained conditional mapping','Five ratio nodes × current income')):
        C=arrays[name];m=C@arrays['z_weights'];grid=arrays['b_grid'];y=np.asarray([annual_gross_income_at_state(ctx['P'],0,0,float(z)) for z in arrays['z_grid']])
        ax[0].step(grid,np.cumsum(m),where='post',label=label);ax[1].plot(y,grid@C,label=label)
    ax[0].set(xlabel='Entrant net financial wealth',ylabel='Cumulative entrant probability',xlim=(-1,4))
    ax[1].set(xlabel='Current annual gross entrant income',ylabel='Conditional mean entrant wealth')
    delta=arrays[ARMS[1]]-arrays[ARMS[0]];im=ax[2].imshow(delta,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-abs(delta).max(),vmax=abs(delta).max());fig.colorbar(im,ax=ax[2],label='Conditional probability difference')
    ax[2].set(xlabel='Income-state index',ylabel='Wealth-grid index',title='Candidate minus control')
    ax[0].legend(fontsize=8);fig.suptitle('Initial entrant distributions');fig.tight_layout();fig.savefig(out/'initial_entry_distribution.png',dpi=150);plt.close(fig)
    return rows

def seed_fit(ctx,live):
    """Native observers at prescribed reference prices; zero extra lifecycle calls."""
    fp,prepared,manifest,objective,runtime,_=ge._observer_context(ctx)
    cal=prepared.rt['primitive'].pf.calendar;P,grid,sd,sol=(live[k] for k in ('P','b_grid','sd','sol'))
    policy=cal.policy_from_solution(sol,live['price'],P,grid,sd);pre,_=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,sd)
    supply=cal.HousingSupplyRule('static-elastic',float(live['price'][0]),float(P.H0[0]*(P.user_cost_rate*live['price'][0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
    ev=cal.evaluate_period(live['price'],pre,P,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=policy)
    fertility={k:prepared.rt['observe_initial_fertility'](ev,P,age_projection=k) for k in ('uniform_birth_time','constant_post_cell')}
    housing=prepared.rt['observe_initial_housing_wealth'](ev,P,grid,sd,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
    recent=prepared.rt['observe_recent_parent_flow'](ev,P,diagnostic_enabled=True,snapshot=prepared.rt['SNAPSHOT'],age_projection=prepared.rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance=dict(case_id='prescribed_reference_price_seed',reference_checkpoint_sha256=manifest['checkpoint']['sha256']))
    (Path(ctx['out'])/'seed_reference_price').mkdir(parents=True,exist_ok=True)
    fits=runtime.score_targets(objective,fertility,housing,recent['model_value'],float(prepared.rt['chain'].extract_moments(sol,P)['tfr']))
    require(len(fits)==14,'Reference-price seed needs all14 target rows');fp.table(Path(ctx['out'])/'seed_reference_price/target_fit.csv',fits)
    write(Path(ctx['out'])/'seed_reference_price/interpretation.json',dict(price=float(live['price'][0]),lifecycle_solves_added=0,interpretation='Prescribed reference prices, fixed benefit and native cohort distribution; not a renewed demographic equilibrium.'))


def arm(args,deadline):
    name=args.arm_name;out=args.out/name;require(not (out/'completed.json').exists(),'Refusing completed arm');out.mkdir(parents=True,exist_ok=True)
    ctx=context(out);ctx['deadline_epoch']=deadline
    if not args.mock_child: authored.authenticate_frozen(ctx)
    entry=apply_entry(ctx,name);ctx.update(selected_d_bar=D,reference_psi=float(ctx['P'].psi_child),expected_dimensions={'wealth_grid_nodes':160,'income_states':15})
    write(out/'effective_input_contract.json',dict(entry=entry,credit=D,reference=ctx['loaded'].identity,psi_fixed=ctx['reference_psi'],closure='price clears actual birth renewal; population clears fixed physical supply',experimental=True,recalibration=False))
    budget=base.ArmBudget(out,deadline,smoke=args.mock_child)
    if args.mock_child:
        ctx.update(fp=SimpleNamespace(write=write),prepared=object(),manifest={},objective={},runtime=object(),reference={});(out/'phase_b_ge').mkdir();result=ge.smoke_phase_b(ctx,budget)
    else:
        # Check the actual frozen calendar entrant cohort before any lifecycle call.
        cal=ctx['prepared'].rt['primitive'].pf.calendar;cohort=cal.entrant_cohort(np.asarray([1.]),ctx['P'],ctx['b_grid'])
        expected=ctx['P'].fixed_reference_entry_conditional*ctx['P'].z_weights[None,:]
        require(cohort.ndim==6 and cohort.shape[0]==160 and cohort.shape[3]==15,'Calendar entrant axes differ')
        np.testing.assert_allclose(cohort.sum(axis=(1,2,4,5)),expected,rtol=0,atol=2e-16)
        write(out/'calendar_entry_verification.json',dict(maximum_joint_error=float(abs(cohort.sum(axis=(1,2,4,5))-expected).max()),lifecycle_solves=0))
        require(abs(float(cohort.sum())-1)<2e-12,'Observer entrant mass differs')
        original=ge.observe_price
        def observe(c,live,label,*,final=False):
            require(float(getattr(live['P'],'_entry_censored_mass',0.))==0.,'Forward entry relocation is prohibited')
            receipt=original(c,live,label,final=final)
            if label=='phase_a_selected_qref': seed_fit(c,live)
            if final: np.savez_compressed(out/'phase_b_ge'/label/'common_support_policies.npz',**base.policy_snapshot(live))
            return receipt
        ge.observe_price=observe
        try:
            seed=base.single_price.solve_fixed_price(ctx,D,ctx['q_ref'],budget,'phase_a_fixed_d53',out/'seed')
            result=ge.run_phase_b(ctx,dict(selected_d_bar=D,selected_live=seed),budget)
            require(result['status']=='passed','GE uncomputed within unchanged budget')
        finally:ge.observe_price=original
    write(out/'completed.json',dict(result=result,lifecycle_solves=budget.used_lifecycle,dimensions=[160,15]))


def compare(out):
    for location in ('phase_b_ge/selected_root','seed_reference_price'):
        pair=[]
        for name in ARMS:
            with (out/name/location/'target_fit.csv').open(newline='') as s:pair.append(list(csv.DictReader(s)))
        require(len(pair[0])==len(pair[1])==14,'Full target rows absent')
        rows=[]
        for a,b in zip(*pair):
            require(all(a[k]==b[k] for k in ('moment','target','weight','role')),'Target/weight contract differs')
            rows.append(dict(moment=a['moment'],target=a['target'],weight=a['weight'],role=a['role'],control_model=a['model'],candidate_model=b['model'],candidate_minus_control=float(b['model'])-float(a['model']) if a['model'] and b['model'] else '',control_gap=a['gap'],candidate_gap=b['gap'],control_loss_contribution=a['loss_contribution'],candidate_loss_contribution=b['loss_contribution']))
        path=out/('comparison_target_fit.csv' if location.startswith('phase') else 'comparison_reference_price_target_fit.csv')
        with path.open('w',newline='') as s:w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    folders=[out/n/'phase_b_ge/selected_root' for n in ARMS]
    params=[]
    for f in folders:
        with (f/'parameters.csv').open(newline='') as s:params.append(list(csv.DictReader(s)))
    require(len(params[0])==31 and params[0]==params[1],'31parameter rows must be identical')
    write(out/'comparison_parameters.json',params[0])
    closures=[json.loads((f/'closure.json').read_text()) for f in folders];write(out/'comparison_closure.json',dict(control=closures[0],candidate=closures[1]))
    changes={}
    with np.load(folders[0]/'common_support_policies.npz') as a,np.load(folders[1]/'common_support_policies.npz') as b:
        for k in ('b_grid','z_grid'):np.testing.assert_array_equal(a[k],b[k])
        common=(a['V']>-1e9)&(b['V']>-1e9)
        child_ok=np.arange(a['V'].shape[6])[None,:]<=np.arange(a['V'].shape[5])[:,None];common &= child_ok.reshape((1,1,1,1,1)+child_ok.shape)
        for k in ('V','bp_pol','c_pol','hR_pol','owner_choice_probability'):
            diff=b[k]-a[k];require(np.isfinite(diff[common]).all(),'Nonfinite common policy');changes[k]=dict(common_feasible_max_abs=float(abs(diff[common]).max()),common_feasible_mean_abs=float(abs(diff[common]).mean()),interpretation='Unweighted common feasible states on identical grid; not occupied-mass accuracy')
    write(out/'comparison_policies.json',changes)
    # Keep every standard plot and verify both exact repeats, in addition to native table/array gates.
    plot_receipts={}
    for name in ARMS:
        a=out/name/'phase_b_ge/selected_root/standard_diagnostics';b=out/name/'phase_b_ge/selected_repeat_final/standard_diagnostics'
        hashes={p.name:sha256_file(p) for p in a.glob('*.png')};twins={p.name:sha256_file(p) for p in b.glob('*.png')}
        require(len(hashes)==17 and hashes==twins,'17standard PNG exact repeat differs');plot_receipts[name]=hashes
    write(out/'standard_plot_repeat_hashes.json',plot_receipts)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=('preflight','full','arm'));ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float);ap.add_argument('--arm-name',choices=ARMS);ap.add_argument('--mock-child',action='store_true');args=ap.parse_args()
    if args.mode!='arm':require(not args.out.exists(),'Refusing existing output');args.out.mkdir(parents=True)
    started=time.time();deadline=min(args.deadline_epoch or started+2400,started+2400)
    try:
        verify_sources()
        if args.mode=='arm':
            if not args.mock_child:require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Native arm requires single-core Torch Slurm')
            arm(args,deadline);return
        entry_preflight(args.out)
        if args.mode=='full':require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Full batch requires single-core Torch Slurm')
        for name in ARMS:
            command=[sys.executable,str(__file__),'arm','--arm-name',name,'--out',str(args.out),'--deadline-epoch',str(deadline)]
            if args.mode=='preflight':command+=['--mock-child']
            env=dict(os.environ,NUMBA_CACHE_DIR=str(args.out/name/'numba_cache'));t=time.monotonic();p=subprocess.run(command,env=env,timeout=max(.01,deadline-time.time()));require(p.returncode==0,'Arm failed: '+name)
            write(args.out/name/'workflow_timing.json',dict(seconds=time.monotonic()-t,includes='imports/authentication/compilation/GE/observers/standard plots/exact repeat'))
        if args.mode=='full':compare(args.out)
        solves=sum(json.loads((args.out/n/'completed.json').read_text())['lifecycle_solves'] for n in ARMS)
        require(solves<=(0 if args.mode=='preflight' else 40),'Solve cap violated')
        write(args.out/'completed.json',dict(status='mock_exact_loop_zero_solves' if args.mode=='preflight' else 'full_passed',lifecycle_solves=solves,total_workflow_seconds=time.time()-started,maximum_lifecycle_solves=40,budget_seconds=2400,recalibration=False,no_auto_retry=True))
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),no_auto_retry=True));raise
if __name__=='__main__':main()
