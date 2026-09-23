#!/usr/bin/env python3
"""No-solve reconstruction of inherited level-marginal versus July ratio entry."""
import copy,gzip,json,sys,types,pickle
from pathlib import Path
import numpy as np
repo=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
bundle=repo/'tmp/utility_overnight_20260923_v1'
out=repo/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry'
# V5 checkpoint was serialized with Python's private pathlib module path.
import pathlib
pathlib.__path__=[]
local=types.ModuleType('pathlib._local'); local.Path=pathlib.Path; local.PosixPath=pathlib.PosixPath; local.PurePath=pathlib.PurePath
sys.modules['pathlib._local']=local
# Checkpoint constructors import plot helpers, but the object read below contains
# no plots; stub only the module initialization so this remains an inspection.
mat=types.ModuleType('matplotlib'); mat.__path__=[]; mat.use=lambda *a,**k:None
plt=types.ModuleType('matplotlib.pyplot'); mat.pyplot=plt
sys.modules['matplotlib']=mat; sys.modules['matplotlib.pyplot']=plt
sys.path[:0]=[str((bundle/'source/code/model/tools').resolve()),str((bundle/'source/code/model').resolve()),str((bundle/'source/code/model/intergen_eqscale_seq_optimized').resolve()),str((bundle/'tools').resolve())]
from intergen_eqscale_seq_optimized import solver as model
import build_period_earnings_process as period_income
import e5f_earnings_wealth_contract as contract
plan=json.loads((bundle/'inputs/utility_templates/B_floor/template_plan.json').read_text())
with gzip.open(bundle/'inputs/normalized_checkpoint/normalized_old.pkl.gz','rb') as f:
    checkpoint=pickle.load(f)
old=checkpoint['old']; P=copy.deepcopy(old.parameters); grid=np.asarray(old.b_grid,dtype=float)
spec=plan['income_specification']['constructor_arguments']
overrides,metadata=period_income.build_period_earnings_process(**spec)
old_cond=np.zeros((len(grid),len(P.z_grid)))
old_z,old_w,_=model.income_transition_values(P)
for k,z in enumerate(old_z):
    ii,ww=model.entry_wealth_grid_weights(grid,P,i=0,j=0,z_value=float(z))
    old_cond[ii,k]=ww
old_marg=old_cond@old_w
# The B15 candidate uses the same P.income[0,0], four-year period and payroll
# gross-up; only the persistent z grid is replaced by the frozen B15 input.
P15=copy.deepcopy(P)
for k,v in overrides.items(): setattr(P15,k,copy.deepcopy(v))
P15.Nz=len(P15.z_grid)
new_z,new_w,_=model.income_transition_values(P15)
new_cond,receipt=contract.rank_coupled_entry(model,P,grid,metadata['components']['persistent_weights'],metadata['components']['iid_weights'])
new_annual=np.array([model.annual_gross_income_at_state(P15,0,0,float(z)) for z in new_z])
ratio_joint=new_cond*new_w[None,:]
ratio_values=grid[:,None]/new_annual[None,:]
keep=ratio_joint>0
x=ratio_values[keep]; w=ratio_joint[keep]
def wmean(x,w): return float(np.sum(x*w)/np.sum(w))
def wq(x,w,p):
    o=np.argsort(x,kind='stable'); xx=x[o]; ww=w[o]
    return float(xx[np.searchsorted(np.cumsum(ww),p*np.sum(ww),side='left')])
def bins(x,w,n=5):
    o=np.argsort(x,kind='stable'); xx=x[o]; ww=w[o]
    b=np.minimum(n,np.maximum(1,np.ceil(np.cumsum(ww)/np.sum(ww)*n).astype(int)))
    return [dict(bin=int(k),weight=float(ww[b==k].sum()/ww.sum()),mean=float(np.sum(xx[b==k]*ww[b==k])/ww[b==k].sum())) for k in range(1,n+1)]
def weighted_spearman_tied_midranks(x,y,w):
    """Weighted Pearson correlation of weighted midpoint ranks; ties share rank."""
    def midrank(a,ww):
        _,inverse=np.unique(a,return_inverse=True)
        mass=np.bincount(inverse,weights=ww)
        rank=(np.cumsum(mass)-.5*mass)/mass.sum()
        return rank[inverse]
    rx=midrank(x,w); ry=midrank(y,w); sw=w.sum()
    mx=float(np.sum(w*rx)/sw); my=float(np.sum(w*ry)/sw)
    vx=float(np.sum(w*(rx-mx)**2)/sw); vy=float(np.sum(w*(ry-my)**2)/sw)
    cov=float(np.sum(w*(rx-mx)*(ry-my))/sw)
    return cov/np.sqrt(vx*vy)
wealth_joint=np.broadcast_to(grid[:,None],ratio_joint.shape)[keep]
income_joint=np.broadcast_to(new_annual[None,:],ratio_joint.shape)[keep]
imposed_spearman=weighted_spearman_tied_midranks(wealth_joint,income_joint,w)
# Exact saved wealth marginal and new income-state conditional wealth means/ratios.
conditional=[]
for j,z in enumerate(new_z):
    m=float(grid@new_cond[:,j]); y=float(new_annual[j])
    conditional.append(dict(state=j+1,z=float(z),weight=float(new_w[j]),annual_gross_income=y,
      conditional_mean_wealth=m,implied_mean_wealth_income_ratio=m/y,
      july_ratio_node_conditional_mean=float(P.entry_wealth_ratio_nodes@P.entry_wealth_ratio_weights),
      conditional_wealth_mean_july=float(y*(P.entry_wealth_ratio_nodes@P.entry_wealth_ratio_weights))))
# Extend grid only as in the frozen numerical-domain contract. All source entry
# wealth points lie within old support, so this does not alter their values.
summary={
 'schema':'entry_units_model_mapping_v1',
 'inputs':{
   'utility_bundle_manifest_sha256':__import__('hashlib').sha256((bundle/'manifest.json').read_bytes()).hexdigest(),
   'B_floor_template_sha256':plan['files']['income']['sha256'],
   'template_plan_sha256':__import__('hashlib').sha256((bundle/'inputs/utility_templates/B_floor/template_plan.json').read_bytes()).hexdigest(),
   'checkpoint_sha256':__import__('hashlib').sha256((bundle/'inputs/normalized_checkpoint/normalized_old.pkl.gz').read_bytes()).hexdigest(),
   'entry_reference_audit_sha256':__import__('hashlib').sha256((repo/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/entry_reference_audit.json').read_bytes()).hexdigest(),
   'selected_B_checkpoint_audit_sha256':__import__('hashlib').sha256((repo/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/selected_B_checkpoint_audit.json').read_bytes()).hexdigest(),
   'mapping':'B15 period-process input and saved pre-solve inherited checkpoint; no solver/equilibrium invocation'},
 'definition':{'wealth':'beginning/pre-choice liquid wealth, model units','income':'age-18 annual gross income = P.income[0,0]*z/(period_years*(1-tau_pay))','period_years':float(P15.period_years),'tau_pay':float(P15.tau_pay),'P_income_0_0':float(P15.income[0,0]),'age18_income_profile':float(P15.income_age_profile[0]),'rank_rule':'old total-income rank intervals overlapped with new B persistent-rank intervals; no perfect-correlation claim'},
 'B15':{'z_grid':new_z.tolist(),'z_weights':new_w.tolist(),'weighted_mean_z':float(new_w@new_z),'annual_gross_income_by_state':new_annual.tolist(),'weighted_mean_annual_gross_income':float(new_w@new_annual),'annual_gross_income_min':float(new_annual.min()),'annual_gross_income_max':float(new_annual.max())},
 'inherited_reference_level_marginal':{'mean_wealth':float(grid@old_marg),'min_supported_grid_wealth':float(grid[old_marg>0].min()),'max_supported_grid_wealth':float(grid[old_marg>0].max()),'raw_reference_conditional_shape':list(old_cond.shape),'new_rank_coupled_conditional_shape':list(new_cond.shape),'marginal_max_abs_gap':float(np.max(np.abs(new_cond@new_w-old_marg))),'marginal_L1_gap':float(np.sum(np.abs(new_cond@new_w-old_marg)))},
 'imposed_entry_wealth_income_rank_dependence':{'weighted_spearman':float(imposed_spearman),'definition':'Weighted Pearson correlation of weighted midpoint ranks; exact ties share their pooled weight rank. Weights are B15 entry joint probabilities over beginning/pre-choice wealth-grid nodes and age-18 annual gross income states.','sample':'Model age-18 entrants, not the empirical age-18-24 PSID family-year sample.','interpretation':'Generated by the inherited wealth-level marginal plus interval-overlap rank coupling to B15 persistent income ranks; it is imposed and not estimated.'},
 'implied_ratio_if_keep_level_marginal':{'weighted_mean':wmean(x,w),'weighted_median':wq(x,w,.5),'weighted_quintile_means':bins(x,w),'pooled_distribution_states':int(len(x))},
 'july_ratio_node_mapping':{'nodes':np.asarray(P.entry_wealth_ratio_nodes).tolist(),'weights':np.asarray(P.entry_wealth_ratio_weights).tolist(),'weighted_mean':float(P.entry_wealth_ratio_nodes@P.entry_wealth_ratio_weights),'weighted_median':wq(np.asarray(P.entry_wealth_ratio_nodes),np.asarray(P.entry_wealth_ratio_weights),.5),'weighted_quintile_means':bins(np.asarray(P.entry_wealth_ratio_nodes),np.asarray(P.entry_wealth_ratio_weights)),'conditional_mean_ratio_each_state':float(P.entry_wealth_ratio_nodes@P.entry_wealth_ratio_weights),'pooled_entry_wealth_mean':float((P.entry_wealth_ratio_nodes@P.entry_wealth_ratio_weights)*(new_w@new_annual))},
 'state_conditional':conditional,
 'rank_mapping_receipt':{k:v for k,v in receipt.items() if k not in ('reference_wealth_marginal','candidate_wealth_marginal')}
}
(out/'model_entry_mapping.json').write_text(json.dumps(summary,indent=2)+'\n')
with (out/'model_entry_mapping_by_state.csv').open('w') as f:
    f.write('state,z,weight,annual_gross_income,conditional_mean_wealth,implied_mean_wealth_income_ratio,july_ratio_node_conditional_mean,conditional_wealth_mean_july\n')
    for r in conditional: f.write(','.join(str(r[k]) for k in ['state','z','weight','annual_gross_income','conditional_mean_wealth','implied_mean_wealth_income_ratio','july_ratio_node_conditional_mean','conditional_wealth_mean_july'])+'\n')
print(json.dumps({k:summary[k] for k in ['definition','B15','inherited_reference_level_marginal','implied_ratio_if_keep_level_marginal','july_ratio_node_mapping']},indent=2))
