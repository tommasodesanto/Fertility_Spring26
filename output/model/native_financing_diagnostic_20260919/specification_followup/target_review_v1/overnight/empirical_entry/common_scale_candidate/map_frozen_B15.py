#!/usr/bin/env python3
"""Map aggregate 3x5 empirical entry nodes to frozen B15 z-rank intervals; no solve."""
from pathlib import Path
import csv, json, hashlib, math, re
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/common_scale_candidate'
PREV=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry'
BUNDLE=ROOT/'tmp/utility_overnight_20260923_v1'
B15=json.loads((PREV/'model_entry_mapping.json').read_text())
GRID=json.loads((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/wealth_grid.json').read_text())
LIFECYCLE=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/native_lifecycle_2023.csv'
PARAMS=BUNDLE/'source/code/model/intergen_eqscale_seq/parameters.py'
actual=json.loads((HERE/'actual_frozen_parameters.json').read_text())['arms']['B_floor']
period=actual['period_years']; age_start=actual['working_age_nodes'][0]
J_R=len(actual['working_age_nodes']); da=period
breaks=actual['income_age_breaks']; values=actual['income_age_values']
profile_norm=actual['income_age_profile'][:J_R]
profile_mean=sum(profile_norm)/J_R
# Actual saved parameters, with frozen B15 and utility bindings: never source defaults.
with LIFECYCLE.open() as f: life=list(csv.DictReader(f))
worker_life=[r for r in life if float(r['age_node'])<age_start+J_R*da]
if len(worker_life)!=J_R: raise ValueError('Saved B lifecycle does not cover all expected working-age cells')
age_mass=actual['working_age_mass']
if max(abs(a-float(r['mass'])) for a,r in zip(age_mass,worker_life))>1e-12:
 raise ValueError('Actual B15 analytic working-age masses disagree with saved B lifecycle')
if min(age_mass)<=0: raise ValueError('Nonpositive saved working-age mass')
# Model has one location; use B15 invariant z distribution (mean one) at every working age.
zgrid=[float(v) for v in B15['B15']['z_grid']]; zp=[float(v) for v in B15['B15']['z_weights']]
mean_z=sum(a*b for a,b in zip(zgrid,zp))
if abs(sum(zp)-1)>1e-12: raise ValueError('B15 weights do not sum to one')
p0=actual['working_aftertax_period_income'][0]; tau=actual['tau_pay']
p0_norm=profile_norm[0]
# Reconstruct P.income[j] from the saved age-18 after-tax period income and frozen normalized age profile.
P_income=actual['working_aftertax_period_income']
model_work_mean=sum(m*(inc*mean_z/(period*(1-tau))) for m,inc in zip(age_mass,P_income))/sum(age_mass)
if abs(model_work_mean-actual['mean_annual_gross_working_income'])>1e-12:
 raise ValueError('Independent mean-income formula disagrees with native income_at_state')
# Actual cell probabilities and conditional means from the corrected empirical 3x5 table.
with (HERE/'entry_nodes_3x5.csv').open() as f: cells=list(csv.DictReader(f))
pi=[float(x['probability']) for x in cells]; omega=[float(x['conditional_mean_omega']) for x in cells]
if abs(sum(pi)-1)>1e-12: raise ValueError(f'Cell probabilities sum to {sum(pi)}')
pt=[sum(pi[t*5:(t+1)*5]) for t in range(3)]
qcond=[[pi[t*5+q]/pt[t] for q in range(5)] for t in range(3)]
te=[0.]
for v in pt: te.append(te[-1]+v)
ze=[0.]
for v in zp: ze.append(ze[-1]+v)
overlap=[[max(0.,min(te[t+1],ze[j+1])-max(te[t],ze[j])) for j in range(len(zp))] for t in range(3)]
# P(t,q,z)=P(q|t)*overlap(t,z): preserve both exact marginals.
J=[]
for t in range(3):
 for q in range(5): J.append([qcond[t][q]*overlap[t][j] for j in range(len(zp))])
if any(abs(sum(J[k])-pi[k])>2e-12 for k in range(15)): raise ValueError('Wealth node marginal not preserved')
if any(abs(sum(J[k][j] for k in range(15))-zp[j])>2e-12 for j in range(len(zp))): raise ValueError('B15 income marginal not preserved')
b=[v*model_work_mean for v in omega]
# Native age-18 resources are R_gross*b + P.income[0,0]*z (after-tax four-year flow).
R=actual['R_gross']
y_net_period=[p0*z for z in zgrid]
wealth_grid=[float(v) for v in GRID['grid']]
with (HERE/'discretization_support.csv').open() as f: support={ (int(r['earnings_tercile']),int(r['wealth_quintile']),int(r['z_state'])):r for r in csv.DictReader(f)}
with (HERE/'candidate_metrics.csv').open() as f: empirical_metrics={r['metric']:float(r['value']) for r in csv.DictReader(f)}
resources=[[R*b[k]+y_net_period[j] for j in range(len(zgrid))] for k in range(15)]
negative=[[resources[k][j]<=0 for j in range(len(zgrid))] for k in range(15)]
outside=[v<min(wealth_grid) or v>max(wealth_grid) for v in b]
def weighted_rank(vals,weights):
 keys=sorted(set(vals)); masses={x:sum(w for v,w in zip(vals,weights) if v==x) for x in keys}; lower=0.; rank={}
 for x in keys: rank[x]=(lower+.5*masses[x])/sum(weights); lower+=masses[x]
 return [rank[x] for x in vals]
def wcov(a,c,w):
 ma=sum(x*z for x,z in zip(a,w))/sum(w); mc=sum(x*z for x,z in zip(c,w))/sum(w)
 return sum(z*(x-ma)*(y-mc) for x,y,z in zip(a,c,w))/sum(w)
wealth=[]; income=[]; weight=[]
for k in range(15):
 for j in range(len(zgrid)): wealth.append(b[k]);income.append(zgrid[j]);weight.append(J[k][j])
rw=weighted_rank(wealth,weight); rz=weighted_rank(income,weight)
sp=wcov(rw,rz,weight)/math.sqrt(wcov(rw,rw,weight)*wcov(rz,rz,weight))
outside_mass=sum(sum(J[k]) for k in range(15) if outside[k])
negative_mass=sum(J[k][j] for k in range(15) for j in range(len(zgrid)) if negative[k][j])
profile_source=BUNDLE/'source/code/model/intergen_eqscale_seq/parameters.py'
summary={
 'schema':'common_scale_candidate_B15_mapping_v1',
 'inputs':{'empirical_nodes_sha256':hashlib.sha256((HERE/'entry_nodes_3x5.csv').read_bytes()).hexdigest(),'prior_frozen_B15_mapping_sha256':hashlib.sha256((PREV/'model_entry_mapping.json').read_bytes()).hexdigest(),'selected_B_grid_sha256':hashlib.sha256((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/wealth_grid.json').read_bytes()).hexdigest(),'selected_B_checkpoint_audit_sha256':hashlib.sha256((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/selected_B_checkpoint_audit.json').read_bytes()).hexdigest(),'actual_parameters_sha256':hashlib.sha256((HERE/'actual_frozen_parameters.json').read_bytes()).hexdigest(),'income_age_profile_source_sha256':hashlib.sha256(profile_source.read_bytes()).hexdigest(),'saved_lifecycle_sha256':hashlib.sha256(LIFECYCLE.read_bytes()).hexdigest(),'utility_manifest_sha256':hashlib.sha256((BUNDLE/'manifest.json').read_bytes()).hexdigest(),'utility_hash_manifest_sha256':hashlib.sha256((BUNDLE/'hash_manifest.json').read_bytes()).hexdigest()},
 'method':{'rank_overlap':'Actual earnings-tercile probability intervals overlapped with B15 discrete z-state rank intervals; within each empirical tercile, actual conditional omega-quintile probabilities distribute the overlap. Both empirical node and B15 z marginals are preserved exactly.','wealth_scale':'b_model = omega * actual frozen working-age mean annual gross earnings. Actual serialized seed and both utility bindings inspected without solving; native income and analytic stationary age weights agree with saved B lifecycle.','cash_only_resource_sign':'R_gross*b + Y_age18 <= 0 retained only as a descriptive cash-only statistic; withdrawn as a feasibility or necessary-resource test because age-18 unsecured debt rollover is allowed. See renter_slack_correction.json for source-traced renter-floor gate.','return':'Actual inherited R_gross=1.02^4; source defaults differ and are not used; no equilibrium solve.','grid_support':'Compared exact scaled node values with selected B wealth-grid min/max; no clipping, censoring, or interpolation.'},
 'model_working_age_scale':{'annual_gross_mean_per_working_household':model_work_mean,'period_years':period,'J_R':J_R,'age_nodes':[float(r['age_node']) for r in worker_life],'saved_age_mass':[float(v) for v in age_mass],'actual_profile_breaks':breaks,'actual_profile_values':values,'profile_normalized_working_mean':profile_mean,'profile_normalized':profile_norm,'age18_normalized_profile':p0_norm,'P_income_age18_aftertax_period':p0,'tau_pay':tau,'z_weighted_mean':mean_z,'z_weights_source':'Frozen B15 distribution; invariant mean z preserved by rank overlap coupling.'},
 'empirical_comparison':{'family_years':empirical_metrics['July family-years'],'persons':empirical_metrics['July persons'],'E_omega':empirical_metrics['mean omega (weighted)'],'pooled_ratio_of_weighted_sums_NW2_over_EARNINDRRC':empirical_metrics['pooled ratio of weighted sums NETWORTH2R/EARNINDRRC'],'raw_weighted_spearman_omega_gross_earnings':empirical_metrics['raw weighted Spearman omega/gross earnings'],'same_sample_node_imputed_weighted_spearman':empirical_metrics['node-imputed weighted Spearman omega/gross earnings'],'inherited_reference_level_law_imposed_spearman':B15['imposed_entry_wealth_income_rank_dependence']['weighted_spearman'],'age_sample_note':'PSID childless-renter family-years ages18-24; mapped model is age-18 entrant distribution.'},
 'working_age_mass_source':{'source':'Saved selected B seven-state native_lifecycle_2023.csv, 12 working age cells, used only for age weights; all masses equal to numerical precision.','why_usable_for_illustrative_mapping':'Actual inherited survival is enabled, equals one through working ages, and declines in retirement. B15 analytic age masses computed from actual P equal the saved B07 working masses; no source-default survival assumption is used.','resolution_limit':'No B15 stationary equilibrium/checkpoint was solved or validated; this uses the common demographic timing/profile and is a unit mapping only.'},
 'B15':{'states':len(zgrid),'probability_sum':sum(zp),'annual_gross_income_age18_by_state':B15['B15']['annual_gross_income_by_state'],'aftertax_period_income_age18_by_state':y_net_period,'R_gross':R},
 'discretization':{'candidate_product_cells_absent_from_empirical_rank_overlap':sum(1 for r in support.values() if r['candidate_combination_absent_in_ranked_microdata']=='TRUE'),'candidate_mass_on_absent_combinations':sum(float(r['candidate_product_probability']) for r in support.values() if r['candidate_combination_absent_in_ranked_microdata']=='TRUE'),'interpretation':'Product-within-tercile mapping can assign probability to omega-quintile/z-rank combinations with zero support in the finer microdata rank-overlap table.'},
 'candidate':{'wealth_nodes':len(b),'min_b':min(b),'max_b':max(b),'mean_b':sum(x*z for x,z in zip(b,pi)),'mean_omega':sum(x*z for x,z in zip(omega,pi)),'model_imposed_weighted_spearman_wealth_z':sp,'grid_min':min(wealth_grid),'grid_max':max(wealth_grid),'outside_grid_mass':outside_mass,'outside_grid_support_pairs':sum(1 for k in range(15) for j in range(len(zgrid)) if J[k][j]>0 and outside[k]),'nonpositive_native_resources_mass':negative_mass,'nonpositive_support_pairs':sum(1 for k in range(15) for j in range(len(zgrid)) if J[k][j]>0 and negative[k][j]),'minimum_R_b_plus_aftertax_Y':min(min(row) for row in resources),'maximum_R_b_plus_aftertax_Y':max(max(row) for row in resources)}}
(HERE/'B15_mapping_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'mapped_joint_wealth_income.csv').open('w',newline='') as f:
 w=csv.writer(f); w.writerow(['earnings_tercile','omega_quintile','b_model_units','z_state_1based','z_level','aftertax_period_income_age18','joint_probability','empirical_rank_overlap_probability','candidate_combination_absent_in_ranked_microdata','R_b_plus_Y_period_units','cashflow_positive','within_selected_B_grid_support'])
 for t in range(3):
  for q in range(5):
   k=t*5+q
   for j,z in enumerate(zgrid):
    prob=J[k][j];res=resources[k][j];srow=support[(t+1,q+1,j+1)]
    w.writerow([t+1,q+1,f'{b[k]:.16g}',j+1,f'{z:.16g}',f'{y_net_period[j]:.16g}',f'{prob:.16g}',srow['empirical_rank_overlap_probability'],srow['candidate_combination_absent_in_ranked_microdata'],f'{res:.16g}',res>0,min(wealth_grid)<=b[k]<=max(wealth_grid)])
# Compact decomposition of nonpositive native resources on unchanged mapped support.
rows=[]
for j,z in enumerate(zgrid):
 total=sum(J[k][j] for k in range(15)); fail=sum(J[k][j] for k in range(15) if negative[k][j]); rows.append(['income_state',str(j+1),z,total,fail,fail/total if total else None])
for t in range(3):
 total=sum(J[k][j] for k in range(t*5,(t+1)*5) for j in range(len(zgrid))); fail=sum(J[k][j] for k in range(t*5,(t+1)*5) for j in range(len(zgrid)) if negative[k][j]); rows.append(['earnings_tercile',str(t+1),None,total,fail,fail/total if total else None])
for k in range(15):
 total=sum(J[k]);fail=sum(J[k][j] for j in range(len(zgrid)) if negative[k][j]);rows.append(['omega_node',f'{k//5+1}:{k%5+1}',omega[k],total,fail,fail/total if total else None])
with (HERE/'cashflow_failure_decomposition.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['dimension','group','z_or_omega_value','candidate_probability_mass','nonpositive_resource_mass','nonpositive_share_within_group']);w.writerows(rows)
print(json.dumps(summary['model_working_age_scale'],indent=2));print(json.dumps(summary['discretization'],indent=2));print(json.dumps(summary['candidate'],indent=2))
