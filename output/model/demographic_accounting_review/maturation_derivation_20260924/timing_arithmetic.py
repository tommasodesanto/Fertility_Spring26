"""Exact timing-kernel arithmetic and compact saved-receipt reuse; no model solve."""
from fractions import Fraction as F
from math import sqrt
from pathlib import Path
import json, hashlib

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
p=F(2,9); q=1-p
pulse=[]
for k in [1,2,3,4,5,7,10]:
    pulse.append({'years':4*k,'geometric_cdf':float(1-q**k),
                  'geometric_point_mass':float(p*q**(k-1)),
                  'queue_cdf':int(k>=5),
                  'step_new_entry_fraction_of_long_run_increment':float(1-q**k),
                  'step_queue_entry_fraction_of_long_run_increment':int(k>=5),
                  'step_cumulative_extra_entrants_normalized':float(F(k)-q*(1-q**k)/p),
                  'step_queue_cumulative_extra_entrants_normalized':max(k-4,0)})
# Exercise exactly the source queue's pop-append-next-date pattern with4slots.
queue=[F(0)]*4
queue_entries=[]
for t in range(11):
    due=queue.pop(0)
    queue.append(F(1) if t==0 else F(0))
    queue_entries.append({'entry_year':4*(t+1),'entry':float(due)})
assert [x['entry_year'] for x in queue_entries if x['entry']]==[20]
assert sum(x['entry'] for x in queue_entries)==1
# These are probabilities/normalized response kernels, not observed shares.
assert sum(p*q**(k-1) for k in range(1,101))+q**100==1
assert 1-q**4==F(4160,6561)
assert F(4)/p==18
assert F(16)*q/p**2==252
median=next(k for k in range(1,30) if 1-q**k>=F(1,2))*4
assert median==12
hist_rel=Path('output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/child_accounting/initial_parent_death.json')
hist_path=REPO/hist_rel
hist=json.loads(hist_path.read_text())
survivals={int(r['age']):1-float(r['exit_probability']) for r in hist['age_rows']}
assert min(a for a,s in survivals.items() if s<1)==66
assert hist['last_fertile_age']==42
assert 66-42>20
age_examples=[]
for birth_age in [18,30,42]:
    pending=1.; normal=0.; death=0.; mean=0.; departures=[]
    for a in range(birth_age,83,4):
        s=survivals[a]
        normal_k=pending*s*float(p)
        death_k=pending*(1-s)
        years=a-birth_age+4
        mean+=years*(normal_k+death_k)
        departures.append({'years':years,'ordinary':normal_k,'death':death_k})
        normal+=normal_k; death+=death_k
        pending*=s*float(q)
    assert abs(normal+death-1)<1e-14
    assert pending==0
    # Every death release is later than20yr; before20 no mortality adjustment.
    assert all(r['years']>20 or r['death']==0 for r in departures)
    assert abs(sum(r['ordinary']+r['death'] for r in departures if r['years']<=20)-float(1-q**5))<1e-14
    age_examples.append({'birth_parent_age':birth_age,'ordinary_route_probability':normal,
                         'death_route_probability':death,'mean_entry_delay_with_release_years':mean,
                         'earliest_possible_death_release_delay_years':66-birth_age+4})
# Exact compressed forward-cohort moment: ages do not affect household policy
# when policy conditions only on existing state and thinning is symmetric.
for m in [1,2,3]:
    for tagged_count in range(m+1):
        from itertools import product
        from math import comb
        for mn in range(m+1):
            enumerated=F(0)
            for stays in product((0,1),repeat=m):
                if sum(stays)==mn:
                    enumerated+=q**mn*p**(m-mn)*sum(stays[:tagged_count])
            K=F(comb(m,mn))*q**mn*p**(m-mn)
            assert enumerated==K*F(mn,m)*tagged_count
sources={
 str(hist_rel):[[1,29]],
 'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/child_accounting/README.md':[[1,39]],
 'code/model/tools/run_e5f_open_population_transition.py':[[231,238],[351,370],[867,921],[1228,1247],[1479,1482]],
 'code/model/tools/run_e5f_perfect_foresight_transition.py':[[650,685]],
}
record={
 'schema':'maturation_timing_assessment_v1','status':'passed_pure_arithmetic_not_model_experiment',
 'scope':{'model_solves':0,'equilibrium_solves':0,'checkpoint_loads':0,'source_changes':0,'new_empirical_estimates':0},
 'mu':float(p),'period_years':4,'queue_waiting_slots':4,'queue_effect_lag_years':20,
 'geometric_delay_no_parent_death':{'mean_years':18.,'variance_years_squared':252.,'sd_years':sqrt(252),'median_years':median,
 'fraction_enter_before18':float(1-q**4),'fraction_enter_after20':float(q**5)},
 'normalized_pulse_and_step':pulse,'queue_impulse_check':queue_entries,
 'historical_literal_departure_flows':{'M':hist['maturation'],'R':hist['dependents_losing_parent'],'B':hist['births'],
 'R_over_M_plus_R':hist['dependents_losing_parent']/(hist['maturation']+hist['dependents_losing_parent']),
 'R_over_dependent_stock':hist['dependents_losing_parent']/hist['dependents'],
 'R_div_2p1':hist['dependents_losing_parent']/2.1,
 'label':'Historical stationary flow share, not current cohort probability or difference versus queued entry.'},
 'historical_age_specific_illustrations':age_examples,
 'age_specific_limits':'Actual historical survival schedule; no birth-age weights used. These are conditional cohort illustrations, not aggregate estimates or population bounds.',
 'hybrid_death_entry':{'result_under_reviewed_support':'zero additional entries if already-entered offspring are excluded','latest_birth_parent_age':42,'first_mortality_parent_age':66,'minimum_child_age_at_start_of_death_period':24,'minimum_delay_to_next_date_death_release':28,'assumes_queue_lag_years':20,'latest_selected_parameter_reverification':False},
 'cohort_moment_derivation':{'formula':'a_prime(k+1,xprime) += s * choice_transport * K(mprime|m) * (mprime/m) * a(k,x); inject births at k=0; remove death-trigger early entries from matching future queue slot','holds_if':'Identity-symmetric thinning and existing state-only household policies. It is a forward accounting statistic, not an age-restricted dependence law.'},
 'sources':{p:{'sha256':hashlib.sha256((REPO/p).read_bytes()).hexdigest(),'lines':lines} for p,lines in sources.items()},
 'verification':['Exact geometric pulse sums and moments','Exact four-slot pop-append queue impulse produces20yr next-date lag','Step response equals cumulative pulse kernel','Age-specific death/ordinary routes add to1 with terminal release','Mortality adjustment is zero through20yr under reviewed fertility/survival support','Enumerated child subsets reproduce attached-cohort first-moment formula'],
}
record['artifacts']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['timing_assessment.md','timing_arithmetic.py'] if (HERE/n).exists()}
(HERE/'timing_assessment.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':record['status'],'death_flow_share':record['historical_literal_departure_flows']['R_over_M_plus_R'],'age_examples':age_examples},indent=2))
