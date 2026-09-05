from pathlib import Path
import sys,json,pickle,platform,csv
from types import SimpleNamespace
import numpy as np
r=Path.cwd();f=r/'tmp/e5f_overnight_20260905a/frozen_production/code/model';sys.path[:0]=[str(f/'tools'),str(f)]
from intergen_eqscale_seq_optimized import solver as m,kernels as k,parameters as pm
import run_dynamic_population_transition as cal
import run_e5f_open_population_transition as tr
cal.model=m
checks={}
# Exact financed-share threshold: wealth below .4 cannot buy a two-unit house,
# with p=1 and phi=.8; wealth .4+epsilon can.
bg=np.array([-2.,0.,.399999,.400001,1.,3.]); V=np.zeros((6,2,1,4,4));V[:,1]=10.
hc=np.array([[0.,2.]]); he=.94*hc;phi=.8
bm=np.zeros((1,2,4,4));bm[:,1]=-phi*2.; dp=np.zeros_like(bm);dp[:,1]=(1-phi)*2.
_,tc,pr=k.tenure_logit_kernel(V,bg,he,hc,dp,bm,np.zeros((4,4,2,2),bool),np.zeros_like(dp),.005)
assert tc[2,0,0,0,0]==0 and tc[3,0,0,0,0]==1
assert pr[2,0,0,0,0,1]==0 and pr[3,0,0,0,0,1]==1
checks['downpayment_threshold']={'pH':2.,'phi':.8,'threshold':.4,'below_buys':False,'above_buys':True}
# Full one-period saving/tenure/income/child-maturation operator and its Bellman adjoint.
P=pm.apply_overrides(pm.setup_parameters(),{'J':3,'J_R':3,'A_f_end':2,'n_house':1,'H_own':np.array([2.]),'n_parity':4,'child_state_mode':'independent_count','sequential_births':True,'z_grid':np.array([.8,1.2]),'z_weights':np.array([.5,.5]),'Pi_z':np.array([[.8,.2],[.2,.8]]),'income_type_transition':'markov','use_income_types':True,'use_numba_scatter':False})
bg=np.array([-1.,0.,1.,2.]);SD=m.precompute_shared(P,bg);shape=(4,2,1,2,4,4); rng=np.random.default_rng(20260905);g=rng.random(shape)
for n in range(4):g[...,n,n+1:]=0
g/=g.sum();v=rng.normal(size=shape)
polshape=(4,2,1,3,2,4,4);bp=np.broadcast_to(bg[:,None,None,None,None,None,None],polshape).copy(); tc=np.broadcast_to(np.arange(2)[None,:,None,None,None,None,None],polshape).copy();lp=np.ones((4,2,1,1,3,2,4,4))
mp=cal.build_transition_maps(np.array([1.]),P,bg,SD)
gn=m.advance_cohort_one_period_markov_income(g,0,lp,tc,None,bp,P,bg,SD,mp.lmm_idx,mp.lmm_wt,mp.tmx_idx,mp.tmx_wt,True,P.Pi_child,P.Pi_z)
ev=np.empty_like(v)
for z in range(2):
 mix=sum(P.Pi_z[z,zn]*v[:,:,:,zn,:,:] for zn in range(2));ev[:,:,:,z,:,:]=m.apply_child_aging(mix,P,4,2,1,4,4,age_index=0)
gap=float(np.sum(gn*v)-np.sum(g*ev)); assert abs(gap)<1e-14 and abs(gn.sum()-g.sum())<1e-14
checks['bellman_forward_adjoint']={'gap':gap,'mass_error':float(gn.sum()-g.sum()),'includes':['identity tenure and saving','nontrivial income transition','binomial child maturation']}
# At most one birth per period: new first births cannot join same-period second pool.
gpre=np.zeros(polshape);gpre[1,0,0,0,0,0,0]=1.;gpre[1,0,0,0,0,1,1]=1.;fp=np.zeros((4,2,1,3,2,4));fp[...,1]=1.;P._fert2_probs=np.zeros((4,2,1,3,2,2,2,4));P._fert2_probs[...,1,:,:]=1.
gpost,births,_=tr.apply_sequential_fertility(gpre,fp,P)
assert births==2 and gpost[...,1,1].sum()==1 and gpost[...,2,2].sum()==1 and gpost[...,3,:].sum()==0
checks['one_birth_per_period']={'births':births,'parity1_mass':float(gpost[...,1,1].sum()),'parity2_mass':float(gpost[...,2,2].sum()),'parity3_mass':float(gpost[...,3,:].sum())}
# Birth-vintage impulse reaches entrants in next-period states 4,8,...,20 years.
q=[0.]*4;due=[]
for t in range(6):
 x,q=tr.advance_birth_vintage_queue(q,2.1 if t==0 else 0.,1/2.1);due.append(x)
assert due==[0.,0.,0.,0.,1.,0.]
checks['queue_impulse']={'next_state_years':[4,8,12,16,20,24],'due':due}

out={"python":sys.version,"numpy":np.__version__,"numba_available":k.NUMBA_AVAILABLE,"numba_jit_disabled":__import__("os").environ.get("NUMBA_DISABLE_JIT"),"source":"frozen production","checks":checks}
(r/"output/model/e5f_code_architecture_audit_20260905a/invariant_checks.json").write_text(json.dumps(out,indent=2)+"\n");print(json.dumps(out,indent=2))
