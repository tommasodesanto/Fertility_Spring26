"""Independent aggregate arithmetic after actual-seed correction; no raw read/solve."""
from pathlib import Path
import csv, hashlib, json, math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[7]
def read(name):
    return json.loads((HERE / name).read_text())
def rows(name):
    with (HERE / name).open() as f:
        return list(csv.DictReader(f))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

actual = read('actual_frozen_parameters.json')
P = actual['arms']['B_floor']
nodes = rows('entry_nodes_3x5.csv')
joint = rows('mapped_joint_wealth_income.csv')
support = rows('renter_slack_grid_support.csv')
summary = read('B15_mapping_summary.json')
slack = read('renter_slack_correction.json')
pi = np.array([float(r['probability']) for r in nodes])
omega = np.array([float(r['conditional_mean_omega']) for r in nodes])
J = np.zeros((15, 15)); b = np.zeros(15)
for r in joint:
    k = (int(r['earnings_tercile'])-1)*5 + int(r['omega_quintile'])-1
    z = int(r['z_state_1based'])-1
    J[k, z] = float(r['joint_probability']); b[k] = float(r['b_model_units'])
z = np.array(P['z_grid']); zp = np.array(P['z_weights'])
mass = np.array(P['working_age_mass'])
gross = np.array(P['annual_gross_income_by_age_state'])
mean_earn = float(mass @ gross @ zp / mass.sum())
y = P['working_aftertax_period_income'][0] * z
slacks = P['R_gross'] * b[:, None] + y[None, :] - np.minimum(b[:, None], 0.)
occupied = J > 0
grid_slack_gap = 0.; supported = set(); support_mass = 0.
for r in support:
    bg = float(r['grid_b']); zi = int(r['z_state_1based'])-1
    expected = P['R_gross']*bg + y[zi] - min(bg, 0.)
    grid_slack_gap = max(grid_slack_gap, abs(expected-float(r['grid_support_slack'])))
    if float(r['interpolated_probability_mass']) > 0:
        supported.add((int(r['grid_index_0based']), zi))
    support_mass += float(r['interpolated_probability_mass'])
checks = dict(
    both_utility_bindings_same_income_and_credit=all(
        v == actual['arms']['B_shares'][k] for k,v in P.items() if k!='template_sha256'),
    actual_checkpoint_hash_matches=all(sha(ROOT / f)==h for f,h in actual['source_sha256'].items()),
    model_mean_annual_working_gross_is_one=abs(mean_earn-1.)<1e-12,
    node_mapping_uses_actual_model_gross_mean=np.max(abs(b-omega*mean_earn))<1e-12,
    empirical_node_marginal_preserved=np.max(abs(J.sum(axis=1)-pi))<1e-12,
    income_marginal_preserved=np.max(abs(J.sum(axis=0)-zp))<1e-12,
    mass_one=abs(J.sum()-1.)<1e-12 and abs(support_mass-1.)<1e-12,
    native_annual_return_two_percent=abs(P['R_gross']-1.02**4)<1e-15,
    no_credit_or_transfer_override=P['lambda_d']==P['next_period_debt_cap_D']==P['transfer_floor_G0']==P['transfer_floor_Gn']==0.,
    inherited_debt_rollover_at_entry=P['next_period_taper_s']==1.,
    working_survival_one=np.array_equal(np.array(P['survival_probs'])[:11],np.ones(11)),
    all_occupied_continuous_current_slacks_positive=bool(np.all(slacks[occupied]>0)),
    reported_continuous_min_correct=abs(slacks[occupied].min()-slack['candidate_support']['continuous_candidate_node_min_slack'])<1e-12,
    neighboring_grid_slacks_correct=grid_slack_gap<1e-12,
    all_neighboring_grid_current_slacks_positive=all(float(r['grid_support_slack'])>0 for r in support),
    positive_grid_support_count_correct=len(supported)==slack['candidate_support']['interpolation_grid_state_support_points'],
)
if not all(checks.values()):
    raise ValueError(checks)
receipt = dict(schema='entry_common_scale_lead_review_v2', status='verified_saved_aggregate_candidate',
    checks={k:bool(v) for k,v in checks.items()},
    findings=dict(model_mean_annual_working_gross=mean_earn, mean_entry_wealth=float(pi@b),
        empirical_node_marginal_max_gap=float(np.max(abs(J.sum(axis=1)-pi))),
        income_marginal_max_gap=float(np.max(abs(J.sum(axis=0)-zp))),
        min_continuous_current_slack=float(slacks[occupied].min()),
        min_interpolation_current_slack=min(float(r['grid_support_slack']) for r in support),
        positive_node_state_pairs=int(occupied.sum()), positive_grid_state_pairs=len(supported)),
    limitations=['Current childless renter resource check only; no household, continuation, later-child-state or equilibrium validation.',
        'Joint 3x5 node/rank coupling is a proposed discretization; ages18-24 proxy and cross-wave normalization remain assumptions.',
        'Raw microdata construction was reviewed earlier; this pass checks saved aggregate arithmetic and actual serialized parameters.'],
    source_trace=['Actual normalized_old checkpoint, not source defaults, supplies return/profile/survival/credit/transfers.',
        'Native run rebases only supply, initializes utility, binds listed structural coordinates, replaces z, rebuilds debt caps, then analytically balances pension.',
        'Extraction reproduces the relevant binding operations in both arms without model solves.'],
    inputs_sha256={f:sha(HERE/f) for f in ('actual_frozen_parameters.json','entry_nodes_3x5.csv',
        'mapped_joint_wealth_income.csv','renter_slack_grid_support.csv','B15_mapping_summary.json',
        'renter_slack_correction.json')})
(HERE/'lead_review.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(status=receipt['status'],checks=len(checks),findings=receipt['findings']),indent=2))
