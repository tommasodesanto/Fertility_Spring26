"""Replay existing dated matched observer with saved policies; forbid solves."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[key]='1'
import pathlib,sys,json,csv,hashlib,signal,time,resource
from types import SimpleNamespace
import numpy as np
OUT=pathlib.Path(__file__).resolve().parent
ROOT=OUT.parents[4]
BASE=OUT.parent
SRC=BASE/'chain_11/postcheck_results/selected_postcheck'
SNAP=BASE/'chain_11/postcheck_results/runtime_auth/runtime_preparation/native_preparation/source_snapshot'
START=time.time()
def deadline(*args):raise TimeoutError('Focused matched cohort 10-minute budget exhausted')
signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,600)
sys.path[:0]=[str(SNAP/'code/model/tools'),str(SNAP/'code/model'),str(ROOT/'code/model'),str(BASE.parent/'utility_floor_round2_v1')]
import inputs
P,grid=inputs.proposal('floor_s0');P,entry=inputs.entry(P,grid,'nonnegative_mean')
point=json.loads((SRC/'proposed_parameters.json').read_text())['free']
bounds=dict(inputs.LANES['floor_s0']['bounds'],psi_child=(.01,.5))
P=inputs.bind(P,point,bounds,'floor');P.unsecured_credit_limit=0.
import run_dynamic_population_transition as cal
import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as measure
from intergen_eqscale_seq_optimized import solver as model
cal.model=model;cal.apply_fertility=transition.apply_sequential_fertility
# These guards make any accidental solve an immediate error.
def no_solve(*args,**kwargs):raise RuntimeError('A model/lifecycle solve is forbidden in saved-policy observer replay')
cal.solve_policy=no_solve
for name in ('solve_markov_income_at_prices','solve_markov_income_equilibrium','solve_bellman_full_markov_income'):
 if hasattr(model,name):setattr(model,name,no_solve)
arr=SRC/'phase_b_ge/selected_repeat/stage/solution_arrays.npz'
need=('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs','bp_pol_stay','c_pol_stay','g_beginning_distribution','entry_by_loc','p_eq','b_grid')
with np.load(arr,allow_pickle=False)as a:
 sol=SimpleNamespace(**{k:a[k]for k in need if k not in ('b_grid',)})
 sd=SimpleNamespace(**{k[7:]:a[k]for k in a.files if k.startswith('shared.')})
 assert np.array_equal(grid,a['b_grid'])
sd.nc=P.n_parity*P.n_child_states
assert P.child_state_mode=='independent_count'
assert P.beta==point['beta_annual']**P.period_years
assert sd.h_bar[1,1]==point['h_P'] and sd.h_bar[1,0]==0
assert P.psi_child==point['psi_child']
for key,value in point.items():
 if key not in ('beta_annual','h_P'):assert getattr(P,key)==value
price=sol.p_eq;P._fert2_probs=sol.fert2_probs.copy()
policy=cal.policy_from_solution(sol,price,P,grid,sd)
pre,recon=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,sd)
assert recon['stationary_post_fertility_nesting_l1']<5e-9,recon
assert recon['stationary_feasibility_projection_mass']==0,recon
counter=cal.SolveCounter();ev=cal.evaluate_period(price,pre,P,grid,sd,counter,supplied_policy=policy)
branch=measure.begin_dated_first_birth_housing_branch(ev,P,grid,sd,origin_period=0)
captured={}
def profiler(frame,event,arg):
 if event=='return'and frame.f_code is measure.finish_dated_first_birth_housing_branch.__code__:
  captured['treated']=frame.f_locals['treated_current'];captured['control']=frame.f_locals['control_current']
sys.setprofile(profiler)
try: result=measure.finish_dated_first_birth_housing_branch(branch,ev,P,grid,sd,destination_period=1)
finally:sys.setprofile(None)
obs=json.loads((SRC/'phase_b_ge/selected_root/observers.json').read_text())['housing_wealth']
saved=next(r['branch']for r in obs['rows']if 'branch'in r)
for key in ('control_mean_housing','treated_mean_housing','origin_mass','destination_mass','housing_response','treated_continuation_births'):
 assert result[key]==saved[key],(key,result[key],saved[key])
assert counter.total==0
rows=[]
for label,g in captured.items():
 mass=float(g.sum());rent=g[:,0];h=np.where(rent>0,policy.hR_pol[:,0],0);rm=float(rent.sum());rr=float((rent*h).sum());om=float(g[:,1:].sum());orr=sum(float(g[:,t].sum())*float(size)for t,size in enumerate(P.H_own,1));cap=float((rent*(np.abs(h-P.hR_max)<=1e-8)).sum())
 row=dict(branch=label,mass=mass,renter_fraction=rm/mass,owner_fraction=om/mass,renter_mean_rooms=rr/rm,owner_mean_rooms=orr/om,all_mean_rooms=(rr+orr)/mass,renter_cap_fraction=cap/rm,cap_mass_fraction_of_branch=cap/mass,rental_contribution_to_branch_mean=rr/mass,owner_contribution_to_branch_mean=orr/mass)
 assert abs(row['all_mean_rooms']-result[label+'_mean_housing'])<=2e-15,(row,result)
 rows.append(row)
with(OUT/'matched_birth_housing.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
pins={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [pathlib.Path(measure.__file__),pathlib.Path(cal.__file__),pathlib.Path(transition.__file__),pathlib.Path(model.__file__),arr,SRC/'proposed_parameters.json',SRC/'native_selected_repeat.json',ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/proposed_120x9/bundle.json',ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1/proposed_120x9/arrays.npz']}
peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
assert peak<8*1024**3,peak
receipt=dict(status='exact_saved_observer_reproduction',point='verified local chain11; original-weight loss114.13604148137114; not newer Torch109 point',elapsed_seconds=time.time()-START,budget_seconds=600,threads=1,memory_budget_gib=8,observed_peak_rss_bytes=peak,free_parameters_verified=True,saved_room_floor_exact=True,reconstruction=recon,exact_saved_observer_fields=['control_mean_housing','treated_mean_housing','origin_mass','destination_mass','housing_response','treated_continuation_births'],result=result,source_pins=pins,model_solves=counter.total,capture='Python return-event profiler, no source edits',switches='Unavailable: origin-destination tenure tags not preserved in current branch distribution',purchase_affordability='Unavailable in existing observer; requires separate exact state-cash rule measurement')
(OUT/'matched_birth_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(rows,indent=2));print('EXACT',receipt['status'],'SECONDS',receipt['elapsed_seconds'])
