#!/usr/bin/env python3
"""Re-evaluate candidate entry nodes using frozen age-18 renter floor; no solve/raw read."""
from pathlib import Path
import csv, json, math, hashlib
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/common_scale_candidate'
GRID=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/wealth_grid.json'
BUNDLE=ROOT/'tmp/utility_overnight_20260923_v1'
SOLVER=BUNDLE/'source/code/model/intergen_eqscale_seq/solver.py'
PARAM=BUNDLE/'source/code/model/intergen_eqscale_seq/parameters.py'

def brackets(grid, x):
    if x < grid[0] or x > grid[-1]: raise ValueError('wealth node outside frozen grid')
    lo=0
    hi=len(grid)-1
    while hi-lo>1:
        mid=(lo+hi)//2
        if grid[mid] <= x: lo=mid
        else: hi=mid
    if x == grid[-1]: return [(len(grid)-1,1.0)]
    w=(x-grid[lo])/(grid[hi]-grid[lo])
    return [(lo,1-w),(hi,w)]

summary=json.loads((HERE/'B15_mapping_summary.json').read_text())
# Keep the earlier R*b+Y table, but label it as a cash-only statistic.
old_cash=HERE/'cashflow_failure_decomposition.csv'
new_cash=HERE/'cash_only_resource_sign_decomposition.csv'
if old_cash.exists():
 old_cash.rename(new_cash)
if new_cash.exists():
 text=new_cash.read_text().replace('nonpositive_resource_mass','nonpositive_Rb_plus_Y_cash_only_mass_descriptive').replace('nonpositive_share_within_group','cash_only_nonpositive_share_within_group_descriptive')
 new_cash.write_text(text)
grid=[float(x) for x in json.loads(GRID.read_text())['grid']]
R=float(summary['B15']['R_gross'])
y=[float(x) for x in summary['B15']['aftertax_period_income_age18_by_state']]
source=PARAM.read_text(); solver=SOLVER.read_text()
# Frozen source asserts exactly the defaults/overrides applicable to candidate B templates.
checks={
 'lambda_d_default_zero':'P.lambda_d = 0.0' in source,
 'taper_starts_42_ends_62':'P.debt_taper_start_age = 42.0' in source and 'P.debt_taper_end_age = 62.0' in source,
 'transfer_floors_default_zero':'P.transfer_floor_G0 = 0.0' in source and 'P.transfer_floor_Gn = 0.0' in source,
 'age18_next_taper_one':'ages[middle]' in source and 'ages >= end' in source,
 'renter_floor_clipped_to_grid':'np.maximum(renter_borrowing_floor(P, b_grid, j), b_grid[0])' in solver,
 'native_unsecured_floor_formula':'np.minimum(float(s_next) * np.minimum(u, 0.0), -float(D_next))' in source,
 'childless_utility_cbar_zero':'"c_bar_0": 0.0' in (BUNDLE/'source/code/model/tools/e5f_parenthood_utility.py').read_text(),
 'childless_utility_eqscale_zeros_hbar':'c_bar[nn, cs] = 0.0' in solver and 'h_bar[nn, cs] = 0.0' in solver,
 'gate_required_flow_definition':'required_flow = float(SD.c_bar[nn, cs]) + float(r_hat[i]) * float(SD.h_bar[nn, cs])' in solver,
 'gate_slack_definition':'"slack": float(resources - required_flow - floor)' in solver,
}
if not all(checks.values()): raise RuntimeError(f'Frozen-source check failed: {checks}')
# Both B_floor and B_shares templates use the same childless zero cbar/hbar. Their
# parenthood-specific variation applies only when children are present.
rows=[]; gridmass={}; total=0.; continuous_min=math.inf; interp_min=math.inf
continuous_nonpositive=0.; interp_nonpositive=0.; exact_nodes=0
with (HERE/'mapped_joint_wealth_income.csv').open() as f:
 for r in csv.DictReader(f):
  prob=float(r['joint_probability'])
  if prob<=0: continue
  b=float(r['b_model_units']); zix=int(r['z_state_1based'])-1
  floor=max(min(b,0.0),grid[0]) # s_next=1, D_next=0 at age18
  res=R*b+y[zix]
  slack=res-floor # transfer=0; childless cbar=hbar=0
  continuous_min=min(continuous_min,slack); continuous_nonpositive += prob*(slack<=0)
  support=brackets(grid,b)
  for idx,w in support:
   mass=prob*w; bg=grid[idx]
   floorg=max(min(bg,0.0),grid[0])
   slackg=R*bg+y[zix]-floorg
   interp_min=min(interp_min,slackg); interp_nonpositive += mass*(slackg<=0)
   gridmass[(idx,zix)]=gridmass.get((idx,zix),0.0)+mass
   rows.append([int(r['earnings_tercile']),int(r['omega_quintile']),zix+1,prob,b,slack,idx,bg,w,mass,slackg])
  total+=prob; exact_nodes+=1
if abs(total-1)>1e-10: raise ValueError(f'joint mass {total}')
summary2={
 'schema':'age18_childless_renter_slack_correction_v1',
 'scope':'Frozen-code algebra and saved candidate aggregates only; no solve and no raw-person reread.',
 'source_files':{str(PARAM.relative_to(ROOT)):hashlib.sha256(PARAM.read_bytes()).hexdigest(),str(SOLVER.relative_to(ROOT)):hashlib.sha256(SOLVER.read_bytes()).hexdigest()},
 'frozen_source_checks':checks,
 'configuration':{
  'lambda_d':0.0,'debt_taper_start_age':42.0,'debt_taper_end_age':62.0,
  'entry_age':18,'next_age_node':22,'next_period_taper_s':1.0,'next_period_debt_cap_D':0.0,
  'age18_renter_floor':'max(min(b,0), b_grid[0])',
  'transfer_floor_G0':0.0,'transfer_floor_Gn':0.0,
  'childless_cbar_both_overnight_utilities':0.0,'childless_hbar_both_overnight_utilities':0.0,
  'current_resources':'R_gross*b + P.income[0,0]*z',
  'required_flow':'cbar + r_hat*hbar = 0 at childless entry under both frozen utility variants',
  'slack':'R_gross*b + P.income[0,0]*z + transfer - required_flow - renter_floor',
  'R_gross':R,'P_income_age18_aftertax_period':float(summary['model_working_age_scale']['P_income_age18_aftertax_period']),
  'utility_templates':['B_floor (eqscale childless cbar/hbar zero)','B_shares (childless cbar/hbar zero)']},
 'candidate_support':{
  'joint_candidate_mass':total,'positive_node_state_pairs':exact_nodes,
  'continuous_candidate_node_min_slack':continuous_min,
  'continuous_candidate_nonpositive_slack_mass':continuous_nonpositive,
  'interpolation_grid_state_support_points':len(gridmass),
  'interpolation_grid_state_min_slack':interp_min,
  'interpolation_grid_support_nonpositive_slack_mass':interp_nonpositive,
  'grid_min_max':[min(grid),max(grid)],
  'interpretation':'Within the current-period resource/floor gate used by frozen solver.py, all candidate mass passes strictly if recorded slack is positive. This is not Bellman continuation viability, utility optimality, or an equilibrium result. The previously reported Rb+Y sign shares are cash-only descriptive statistics and are withdrawn as feasibility/necessity measures.'},
 'limitations':[
  'No Bellman or equilibrium solve was run; future value and optimal housing/consumption choices are not checked.',
  'All-age/parent states can have positive child-dependent housing floors; this diagnostic concerns childless entrant state only.',
  'Interpolation support applies the age-18 floor at neighboring frozen wealth-grid nodes; it does not certify all future policy transitions.',
  'Candidate common-scale entry law is still not adopted or equilibrium-validated.']}
with (HERE/'renter_slack_grid_support.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['earnings_tercile','omega_quintile','z_state_1based','candidate_joint_probability','continuous_b','continuous_slack_at_age18_renter_floor','grid_index_0based','grid_b','interpolation_weight','interpolated_probability_mass','grid_support_slack']);w.writerows(rows)
(HERE/'renter_slack_correction.json').write_text(json.dumps(summary2,indent=2)+'\n')
# Replace erroneous feasibility labels in the prior candidate summary, keeping raw metrics.
summary['candidate']['cash_only_nonpositive_Rb_plus_Y_mass_descriptive_withdrawn_as_feasibility']=summary['candidate'].pop('nonpositive_native_resources_mass', summary['candidate'].get('cash_only_nonpositive_Rb_plus_Y_mass_descriptive_withdrawn_as_feasibility', 0.017955065869277066))
summary['candidate']['cash_only_nonpositive_support_pairs_descriptive_withdrawn_as_feasibility']=summary['candidate'].pop('nonpositive_support_pairs', summary['candidate'].get('cash_only_nonpositive_support_pairs_descriptive_withdrawn_as_feasibility', 5))
summary['candidate']['cash_only_min_Rb_plus_Y_descriptive']=summary['candidate'].pop('minimum_R_b_plus_aftertax_Y', summary['candidate'].get('cash_only_min_Rb_plus_Y_descriptive', -0.5850845130625311))
summary['candidate']['cash_only_max_Rb_plus_Y_descriptive']=summary['candidate'].pop('maximum_R_b_plus_aftertax_Y', summary['candidate'].get('cash_only_max_Rb_plus_Y_descriptive', 28.527974763456807))
summary['candidate']['actual_childless_renter_slack_correction']='renter_slack_correction.json'
(HERE/'B15_mapping_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary2['candidate_support'],indent=2))
