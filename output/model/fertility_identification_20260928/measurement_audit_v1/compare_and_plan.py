"""Torch-only arithmetic on small saved tables; no model imports or solves."""
import os
import sys
import csv
import json
import hashlib
from pathlib import Path
import numpy as np

assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
B = ROOT/'output/model/fertility_identification_20260928'
OUT = B/'measurement_audit_v1'
INPUT = OUT/'inputs'
SOURCE = B/'resume_v1/selected_export/primary'
read = lambda p: json.loads(Path(p).read_text())
rows = lambda p: list(csv.DictReader(Path(p).open()))


def write(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def table(name, records):
    with (OUT/name).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)


observer = read(SOURCE/'observers.json')['fertility']['uniform_birth_time']
frozen = read(INPUT/'fixed_reference_manifest.json')
assert frozen['checkpoint']['sha256'] == read(SOURCE/'receipt.json')['case_checkpoint_sha256']
assert frozen['source_manifest']['sha256'] == read(SOURCE/'receipt.json')['source_manifest_sha256']
top_weight = float(frozen['actual_serialized_parameters']['tfr_top_bin_weight'])
count_weights = np.array([0, 1, 2, top_weight])
a = observer['accounting']
pre = np.array(a['pre_parity_mass_by_age']); post = np.array(a['post_parity_mass_by_age'])
ages = np.array(a['age_cell_start']); counts = np.arange(4)
data25 = read(INPUT/'early_fertility_target.json')['estimates']['25']
model_shares = np.array([observer['ever_born_shares_age25'][key] for key in ('0','1','2','3plus')])
sd = data25['share_with_any_birth']; sm = 1-model_shares[0]
ed = data25['mean_children_ever_born_capped3']; em = model_shares@counts
kd, km = ed/sd, em/sm
extensive = (sm-sd)*(km+kd)/2
intensive = (km-kd)*(sm+sd)/2
assert abs(extensive+intensive-(em-ed)) < 1e-12
decomposition = dict(data_mother_share=sd, model_mother_share=sm,
    data_capped_children_given_mother=kd, model_capped_children_given_mother=km,
    data_early_fertility=ed, model_early_fertility=em, gap=em-ed,
    symmetric_mother_share_component=extensive,
    symmetric_children_given_mother_component=intensive,
    fraction_gap_from_children_given_mother=intensive/(em-ed),
    model_age25_count_shares=model_shares.tolist(),
    data_above3_share=data25['weighted_share_frever_above3'],
    data_raw_mean=data25['mean_children_ever_born'],
    interpretation='Arithmetic decomposition, not a causal attribution or target-reachability claim')
write('early_fertility_decomposition.json', decomposition)

bycell = np.zeros(7); raw_total=0.; raw_age_sum=0.; under18=0.
for row in rows(INPUT/'first_birth_counts_year_age.csv'):
    if 2003 <= int(row['year']) <= 2006:
        age=int(row['age']); value=float(row['n_first_births'])
        cell=0 if age<=21 else 6 if age>=42 else (age-18)//4
        bycell[cell]+=value; raw_total+=value; raw_age_sum+=age*value
        if age<18: under18+=value
assert raw_total == 6611269
emp_first = bycell/raw_total
model_first = np.array(a['parity_birth_flows_by_age'])[:7,0]/a['first_birth_flow']
assert abs(model_first.sum()-1)<1e-12
table('first_birth_age_cells.csv', [dict(age_lower=int(ages[i]),age_upper=int(ages[i]+3),
    empirical_share=float(emp_first[i]),model_share=float(model_first[i]),
    gap=float(model_first[i]-emp_first[i]),
    empirical_tail_mapping='ages12-21' if i==0 else 'ages42-49' if i==6 else 'same cell',
    role='untargeted individual shares; mean targeted, age30plus validation') for i in range(7)])
fit={r['moment']:r for r in rows(SOURCE/'target_fit.csv')}
assert abs(emp_first@(ages[:7]+2)-float(fit['nchs_mean_age']['target']))<1e-10
assert abs(model_first@(ages[:7]+2)-float(fit['nchs_mean_age']['model']))<1e-10
write('first_birth_timing.json', dict(empirical_counts=raw_total,
    empirical_midpoint_mean=float(emp_first@(ages[:7]+2)),
    empirical_completed_integer_age_mean=raw_age_sum/raw_total,
    empirical_share_under18=under18/raw_total,
    model_mean=float(model_first@(ages[:7]+2)),
    conditional_claude_ceiling=float((1-float(fit['cps_childlessness']['target']))*(1.875*emp_first[0]+.875*emp_first[1])),
    ceiling_assumptions=['Fix empirical first-birth shares','Treat age40-44 motherhood as eventual motherhood',
        'Convert stationary period shares to lifetime probabilities','At most one birth per model cell'],
    active_target_infeasibility_proved=False))

profile=[]
for row in rows(INPUT/'age_profile_candidates.csv'):
    if row['window']!='pooled': continue
    lo=int(row['age_lower']); hi=int(row['age_upper'])+1
    left=np.maximum(ages,lo); right=np.minimum(ages+4,hi)
    overlap=np.maximum(right-left,0)/4
    post_weight=((left+right)/2-ages)/4
    mass=(overlap[:,None]*((1-post_weight[:,None])*pre+post_weight[:,None]*post)).sum(axis=0)
    shares=mass/mass.sum()
    empirical=np.array([float(row[k]) for k in ('share_0','share_1','share_2','share_3plus')])
    model_cap=float(shares@counts); data_cap=float(empirical@counts)
    profile.append(dict(age_lower=lo,age_upper=hi-1,data_capped3=data_cap,model_capped3=model_cap,
        gap_capped3=model_cap-data_cap,data_mother_share=float(1-empirical[0]),
        model_mother_share=float(1-shares[0]),data_given_mother=data_cap/(1-empirical[0]),
        model_given_mother=model_cap/(1-shares[0]),data_3plus=float(empirical[3]),model_3plus=float(shares[3]),
        data_uncapped=float(row['mean_CEB_uncapped']),data_coded3602=float(row['mean_model_coded_CEB']),
        model_coded3602=float(shares@count_weights),role='untargeted lifecycle diagnostic'))
table('fertility_lifecycle_matched_windows.csv',profile)
last=profile[-1]
closed=read(SOURCE/'receipt.json')['normalization']['completed_fertility']
write('stationary_approximation_comparison.json',dict(
    imposed_replacement_normalization=closed,
    actual_saved_top_bin_weight=top_weight,
    empirical_age40_44_uncapped=last['data_uncapped'],
    empirical_age40_44_coded3602=last['data_coded3602'],
    normalization_minus_empirical_uncapped=closed-last['data_uncapped'],
    normalization_minus_empirical_coded3602=closed-last['data_coded3602'],
    model_age40_44_capped3=last['model_capped3'],data_age40_44_capped3=last['data_capped3'],
    model_age40_44_coded3602=last['model_coded3602'],
    model_age40_44_3plus=last['model_3plus'],data_age40_44_3plus=last['data_3plus'],
    actual_model_terminal_post_count_shares=(post[6]/post[6].sum()).tolist(),
    actual_model_terminal_post_capped3=float(post[6]@counts/post[6].sum()),
    actual_model_terminal_post_coded3602=float(post[6]@count_weights/post[6].sum()),
    caution='Normalization and observed age40-44 stock have different timing; capped-three and model-coded means additionally use different top-bin weights. Gaps quantify approximation, not a model error or causal early-gap attribution.'))

# Reconstruct full-precision trial coordinates; predictions remain unsolved.
jac = {(r['parameter'],r['moment']):r for r in rows(INPUT/'jacobian.csv')}
svd=read(INPUT/'jacobian_scaled_svd.json'); params=svd['parameters']; moments=svd['moments']
anchor=np.array([svd['column_scales'][p] for p in params])
weights=np.array([float(fit[m]['weight']) for m in moments])
gaps=np.array([float(fit[m]['gap']) for m in moments])
J=np.array([[float(jac[p,m]['full_derivative'])*anchor[i] for i,p in enumerate(params)] for m in moments])
Jh=np.array([[float(jac[p,m]['half_derivative'])*anchor[i] for i,p in enumerate(params)] for m in moments])
d=np.linalg.solve(J.T@(weights[:,None]*J)+10*np.eye(len(params)),-J.T@(weights*gaps))
alpha=min(.5,.15/max(abs(d)))
step=alpha*d
trial=anchor*np.exp(step)
restrictions={r['parameter']:r for r in rows(SOURCE/'parameters.csv')}
assert all(float(restrictions[p]['lower']) <= trial[i] <= float(restrictions[p]['upper']) for i,p in enumerate(params))
psi=float(restrictions['psi_child']['estimate'])
psi_prediction=psi+sum(float(jac[p,'psi_child']['full_derivative'])*anchor[i]*step[i] for i,p in enumerate(params))
psi_half=psi+sum(float(jac[p,'psi_child']['half_derivative'])*anchor[i]*step[i] for i,p in enumerate(params))
predicted=[]
for moment,r in fit.items():
    jf=np.array([float(jac[p,moment]['full_derivative'])*anchor[i] for i,p in enumerate(params)])
    jh=np.array([float(jac[p,moment]['half_derivative'])*anchor[i] for i,p in enumerate(params)])
    m=float(r['model'])+jf@step; mh=float(r['model'])+jh@step
    wt=float(r['weight']) if r['weight'] else 0
    predicted.append(dict(moment=moment,role=r['role'],target=float(r['target']),
        reference=float(r['model']),predicted_full=float(m),predicted_half=float(mh),
        full_minus_half=float(m-mh),weight=wt,predicted_loss=float(wt*(m-float(r['target']))**2),
        status='linear prediction only; no candidate solve'))
table('proposed_step_predictions.csv',predicted)
table('proposed_step_parameters.csv',[dict(parameter=p,anchor=float(anchor[i]),trial=float(trial[i]),
    log_step=float(step[i]),lower=float(restrictions[p]['lower']),upper=float(restrictions[p]['upper'])) for i,p in enumerate(params)])
write('bounded_experiment_plan.json',dict(status='prepared_not_launched',specification_changes=[],
    economic_changes=[dict(object='ten estimated preference and housing-supply coordinates',
        status='proposed diagnostic recalibration',values='trial_point below; no adoption'),
        dict(object='child-benefit level psi',status='derived within proposed diagnostic',
        rule='renormalize to replacement fertility2.1 and enforce renewal in both arms')],
    reference_label='2007 stationary reference — block0506, September 28 verified export',
    proposal_method='full-step log-coordinate ridge10 Gauss-Newton; uniform damping alpha<=.5 and max logstep<=.15',
    original_undamped_linear_loss=float((gaps+J@d)@(weights*(gaps+J@d))),
    damping=float(alpha),max_absolute_log_step=float(max(abs(step))),
    predicted_full_loss=float((gaps+J@step)@(weights*(gaps+J@step))),
    predicted_half_loss=float((gaps+Jh@step)@(weights*(gaps+Jh@step))),
    trial_point=dict(zip(params,map(float,trial))),
    predicted_initial_psi_full=float(psi_prediction),predicted_initial_psi_half=float(psi_half),
    normalization_rule='Re-solve psi to2.1 and enforce demographic renewal for each evaluation',
    paired_arms=['retained initial psi .14281100340255604','full-step derivative-predicted initial psi'],
    initial_bracket_step=.005,price_initialization='identical retained initialization in both arms',
    objective_count=2,maximum_stationary_solves_per_objective=23,maximum_stationary_solves_total=46,
    objective_cap_seconds=1800,whole_experiment_cap_seconds=4200,workers=2,
    unchanged=['income','initial wealth/income','timing','transfers/floors','preferences functional form',
        'target/weight fingerprint','parameter bounds','all scientific gates','17 standard plots'],
    stop=['after pair','unknown or integrity failure','deadline; no automatic extension/retry'],
    launch_prerequisites=['Resolve any saved-income mismatch before candidate inference',
        'Pin isolated opt-in normalization-start adapter and all inputs',
        'Torch synthetic exact dispatch/receipt/failure/timeout smoke',
        'Implement the explicit paired-consistency screens below before seeing results',
        'Retain complete tables and17 standard plots; no automatic promotion'],
    expected_runtime_basis='Reference6 stationary solves970.676seconds; no assumed warm-start saving',
    paired_consistency=dict(scored_max_abs_difference_in_working_scale=.01,
        validation_absolute_difference_rule='.01*max(1,abs(reference),abs(target))',
        psi_absolute_difference=1e-4,price_relative_difference=1e-4,loss_absolute_difference=.05,
        status='proposed numerical equivalence screens, additional to unchanged scientific gates',
        failure_action='report sensitivity; do not claim equivalent faster normalization or promote'),
    retained_scientific_gates=dict(completed_fertility_absolute_gap=5e-4,
        positive_psi=True,demographic_renewal='original gate unchanged',
        market='original solver and audit gates unchanged',paygo_scaled_residual=1e-6,
        all_other_gates='unchanged full runtime budget, estate, transaction, probability and occupied-value checks'),
    promotion='requires separate final exact-repeat verification; none budgeted or launched by this plan',
    identification='10scored moments/10free coordinates plus separate normalization; numerical full rank is not statistical identification',
    policy_transition_reestimation=False))
write('small_table_input_hashes.json',{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUT.iterdir()) if p.is_file()})
print(json.dumps(dict(decomposition=decomposition,lifecycle=profile,
    plan=read(OUT/'bounded_experiment_plan.json'))))
