#!/usr/bin/env python3
"""Map aggregate 3x5 empirical entry nodes to frozen B15 z-rank intervals, no solve."""
from pathlib import Path
import csv, json, hashlib
import math
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/common_scale_candidate'
PREV=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry'
B15=json.loads((PREV/'model_entry_mapping.json').read_text())
GRID=json.loads((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/wealth_grid.json').read_text())
with (HERE/'entry_nodes_3x5.csv').open() as f: cells=list(csv.DictReader(f))
# Actual cell probabilities from the 1,835-row July sample (ties retained in bins).
with (HERE/'entry_nodes_3x5.csv').open() as f: cells=list(csv.DictReader(f))
pi=[float(x['probability']) for x in cells]
omega=[float(x['conditional_mean_omega']) for x in cells]
if abs(sum(pi)-1)>1e-12: raise ValueError(f'Cell probabilities sum to {sum(pi)}')
pt=[sum(pi[t*5:(t+1)*5]) for t in range(3)]
qcond=[[pi[t*5+q]/pt[t] for q in range(5)] for t in range(3)]
zp=[float(v) for v in B15['B15']['z_weights']]; y=[float(v) for v in B15['B15']['annual_gross_income_by_state']]
if abs(sum(zp)-1)>1e-12: raise ValueError('B15 weights do not sum to 1')
# Overlap old tercile rank intervals with the B15 state rank intervals.
def edges(p):
 e=[0.]
 for v in p: e.append(e[-1]+v)
 return e
te=edges(pt); ze=edges(zp); overlap=[[max(0.,min(te[t+1],ze[j+1])-max(te[t],ze[j])) for j in range(len(zp))] for t in range(3)]
# P(t,q,z)=P(q|t)*overlap(t,z); verifies entrant-node and B15 state marginals.
J=[]
for t in range(3):
 for q in range(5): J.append([qcond[t][q]*overlap[t][j] for j in range(len(zp))])
if any(abs(sum(J[k])-pi[k])>2e-12 for k in range(15)): raise ValueError('Wealth node marginal not preserved')
if any(abs(sum(J[k][j] for k in range(15))-zp[j])>2e-12 for j in range(len(zp))): raise ValueError('B15 income marginal not preserved')
mean_y=sum(a*b for a,b in zip(zp,y)); b=[v*mean_y for v in omega]
# Model return uses frozen source parameter q=(1.04^period_years)-1, R_gross=1+q.
period_years=float(B15['definition']['period_years']); R=1.04**period_years
wealth_grid=[float(v) for v in GRID['grid']]
resources=[[R*b[k]+y[j] for j in range(len(y))] for k in range(15)]
negative=[[resources[k][j]<=0 for j in range(len(y))] for k in range(15)]
outside=[v<min(wealth_grid) or v>max(wealth_grid) for v in b]
# Weighted tie-midrank Spearman, same rule as earlier audit.
def weighted_rank(vals,weights):
 keys=sorted(set(vals)); masses={x:sum(w for v,w in zip(vals,weights) if v==x) for x in keys}; lower=0.; rank={}
 for x in keys: rank[x]=(lower+.5*masses[x])/sum(weights); lower+=masses[x]
 return [rank[x] for x in vals]
def wcov(a,c,w):
 ma=sum(x*z for x,z in zip(a,w))/sum(w); mc=sum(x*z for x,z in zip(c,w))/sum(w)
 return sum(z*(x-ma)*(y-mc) for x,y,z in zip(a,c,w))/sum(w)
wealth=[]; income=[]; weight=[]
for k in range(15):
 for j in range(len(y)):
  wealth.append(b[k]);income.append(y[j]);weight.append(J[k][j])
rw=weighted_rank(wealth,weight); ry=weighted_rank(income,weight)
sp=wcov(rw,ry,weight)/math.sqrt(wcov(rw,rw,weight)*wcov(ry,ry,weight))
outside_mass=sum(sum(J[k]) for k in range(15) if outside[k])
negative_mass=sum(J[k][j] for k in range(15) for j in range(len(y)) if negative[k][j])
summary={
 'schema':'common_scale_candidate_B15_mapping_v1',
 'inputs':{'empirical_nodes_sha256':hashlib.sha256((HERE/'entry_nodes_3x5.csv').read_bytes()).hexdigest(),'prior_frozen_B15_mapping_sha256':hashlib.sha256((PREV/'model_entry_mapping.json').read_bytes()).hexdigest(),'selected_B_grid_sha256':hashlib.sha256((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/wealth_grid.json').read_bytes()).hexdigest(),'selected_B_checkpoint_audit_sha256':hashlib.sha256((ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/selected_B_checkpoint_audit.json').read_bytes()).hexdigest()},
 'method':{'rank_overlap':'Actual earnings-tercile probability intervals overlapped with B15 discrete z-state rank intervals; within each empirical tercile, actual conditional wealth-quintile probabilities distribute the overlap. No correlation or exact quantile mass imposed.','wealth_scale':'b_model = omega * B15 weighted mean age-18 annual gross earnings per working household; candidate omega is nonhousing NETWORTH2R divided by same-wave weighted EARNINDRRC mean among working-age RP.','annual_gross_return':'R=1.04^period_years; period_years read from frozen B15 mapping; no equilibrium solve.','cashflow_test':'R*b + annual gross Y_age18 > 0 evaluated on every positive-mass node pair.','grid_support':'Compared exact scaled node values with selected B wealth-grid min/max; no clipping, censoring, or interpolation.'},
 'B15':{'states':len(zp),'probability_sum':sum(zp),'weighted_mean_annual_gross_earnings':mean_y,'min_annual_gross_earnings':min(y),'max_annual_gross_earnings':max(y),'R_gross':R,'source':B15['definition']},
 'candidate':{'wealth_nodes':len(b),'min_b':min(b),'max_b':max(b),'mean_b':sum(x*z for x,z in zip(b,pi)),'mean_omega':sum(x*z for x,z in zip(omega,pi)),'model_imposed_weighted_spearman':sp,'grid_min':min(wealth_grid),'grid_max':max(wealth_grid),'outside_grid_mass':outside_mass,'outside_grid_support_pairs':sum(1 for k in range(15) for j in range(len(y)) if J[k][j]>0 and outside[k]),'R_b_plus_Y_nonpositive_mass':negative_mass,'nonpositive_support_pairs':sum(1 for k in range(15) for j in range(len(y)) if J[k][j]>0 and negative[k][j]),'minimum_R_b_plus_Y':min(min(row) for row in resources),'maximum_R_b_plus_Y':max(max(row) for row in resources)}}
(HERE/'B15_mapping_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'mapped_joint_wealth_income.csv').open('w',newline='') as f:
 w=csv.writer(f); w.writerow(['earnings_tercile','wealth_quintile','b_model_units','z_state_1based','annual_gross_income_model_units','joint_probability','R_b_plus_Y','cashflow_positive','within_selected_B_grid_support'])
 for t in range(3):
  for q in range(5):
   k=t*5+q
   for j in range(len(y)):
    p=J[k][j];res=resources[k][j];w.writerow([t+1,q+1,f'{b[k]:.16g}',j+1,f'{y[j]:.16g}',f'{p:.16g}',f'{res:.16g}',res>0,min(wealth_grid)<=b[k]<=max(wealth_grid)])
print(json.dumps(summary['candidate'],indent=2))
