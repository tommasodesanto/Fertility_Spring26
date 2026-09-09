"""One dated Bellman replay plus saved-policy KFE; no changed numerical gates."""
import argparse, copy, gzip, json, os, pickle, sys, time
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
p = argparse.ArgumentParser()
p.add_argument('--contract', type=Path, required=True)
p.add_argument('--contract-sha256', required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
root = Path(json.loads(a.contract.read_text())['source_root'])
sys.path[:0] = [str(root/'code/model'), str(root/'code/model/tools')]
import numpy as np
import run_e5f_matched_pf_smoke as smoke
c = smoke.load_contract(a.contract, a.contract_sha256)
assert str(root) == '/scratch/td2248/projects/Fertility_Spring26_matched_pf_20260909e'
a.output.mkdir(parents=True, exist_ok=False)
start = time.monotonic()
smoke.transition.configure_sequential_model()
smoke.calendar.model = smoke.model
smoke.calendar.apply_fertility = smoke.transition.apply_sequential_fertility
with gzip.open(c['checkpoint'], 'rb') as f:
    packet = pickle.load(f)
P = copy.deepcopy(packet['parameters'])
P.joint_nested_choice = P.fertility_nest_choice = True
P.two_shock_choice = False
P.exhaustive_saving_control = True
P.property_tax_lump_sum_transfer = 0.
saved_path = root/'output/model/e5f_matched_pf_20260909a/pilot_04/nested/stationary_arrays.npz'
saved_sha = smoke.digest(saved_path)
with np.load(saved_path, allow_pickle=False) as z:
    saved = {k: z[k] for k in z.files}
grid, pre = saved['wealth_grid'], saved['g_pre']
price = float(packet['evaluation'].policy.price[0])
rent = float(smoke.pf.rents_from_asset_prices([price], price, P)[0])
shared = smoke.model.precompute_shared(P, grid)
result = dict(source_root=str(root), contract_sha256=a.contract_sha256,
              saved_arrays_sha256=saved_sha, scope='one Bellman replay and saved-policy KFE; no gates changed')
def report(stage):
    result.update(stage=stage, elapsed_seconds=time.monotonic()-start)
    smoke.pf.write_json(a.output/'diagnosis.json', result)
    print(stage, flush=True)
report('inputs_verified')
constant = smoke.pf.solve_date_policy(price=price, rent=rent, P=P, b_grid=grid,
                                     shared=shared, continuation_V=saved['V'])
all_actual = smoke.policy_arrays(constant)
result['native_array_gaps'] = {k: float(np.max(np.abs(v-saved[k]))) for k,v in all_actual.items() if k != 'tenure_probs'}
reference = replace(constant, **{k: saved[k] for k in smoke.FIELDS if k in saved},
                    joint_choice=SimpleNamespace(**{k: saved['joint_'+k] for k in smoke.JOINT_FIELDS}))
entry = np.array([float(pre[:,:,:,0].sum())])
supply, _ = smoke.calendar.normalize_date0_housing_supply(pre, reference, P, grid, shared, 'static-elastic')
ev = smoke.calendar.evaluate_period(np.array([price]), pre, P, grid, shared,
    smoke.calendar.SolveCounter(), supply_rule=supply, supplied_policy=constant)
nxt, _, deaths, mass_residual = smoke.transition.advance_sequential_calendar_distribution(ev, entry, P, grid, shared)
diff = np.abs(ev.policy.tenure_probs - saved['tenure_probs'])
state_diff = diff.max(axis=-1)
bad = state_diff > 2e-10
post = ev.g_post_fertility
weighted = post[...,None] * diff
worst = tuple(int(i) for i in np.unravel_index(diff.argmax(), diff.shape))
result['conditioning'] = dict(max_absolute=float(diff.max()), cells_above_2e10=int((diff>2e-10).sum()),
    post_states_above_2e10=int(bad.sum()), post_mass_on_mismatches=float(post[bad].sum()),
    maximum_post_mass_on_mismatches=float(post[bad].max()) if bad.any() else 0.,
    weighted_product_mass_l1=float(weighted.sum()), maximum_weighted_product_mass=float(weighted.max()),
    worst_index=worst, worst_post_mass=float(post[worst[:-1]]),
    worst_reference=float(saved['tenure_probs'][worst]), worst_dated=float(ev.policy.tenure_probs[worst]),
    support=[dict(min_post_mass=threshold, count=int((post>threshold).sum()),
                  max_probability_gap=float(state_diff[post>threshold].max()) if np.any(post>threshold) else 0.,
                  bad_count=int((bad & (post>threshold)).sum())) for threshold in [0.,1e-16,1e-14,1e-12,1e-10,1e-8]])
result['dated_invariants'] = dict(current_l1=float(np.abs(ev.g_current-saved['g_current']).sum()),
    invariant_l1=float(np.abs(nxt-pre).sum()), births=float(ev.births), mass_residual=float(mass_residual),
    feasibility_projection_mass=float(ev.feasibility_projection_mass), deaths=float(deaths))
result['dated_budget'] = smoke.dated_budget(ev,P,shared,grid,rent)
# Same-pool reference factorization is the meaningful native-policy comparison.
ref_post, ref_effective, ref_births, _, _ = smoke.calendar.factor_joint_distribution(ev.g_pre, reference, P)
result['same_pool_reference'] = dict(post_max=float(np.abs(ref_post-post).max()),
    effective_max=float(np.abs(ref_effective-ev.policy.tenure_probs).max()),
    births_gap=abs(float(ref_births.sum())-float(ev.births)))
report('dated_replay_complete')
# Reconstruct the original stationary KFE, without another Bellman optimization.
P._joint_choice = reference.joint_choice
P._fert2_probs = reference.fert2_probs
kfe_tp = saved['tenure_probs'].copy()
g_kfe, stats = smoke.model.forward_distribution_markov_income(
    reference.bp_pol, reference.hR_pol, reference.tenure_choice,
    reference.loc_probs, reference.fert_probs, reference.V, np.array([rent]), np.array([price]),
    P, grid, shared, fast_stats=False, tenure_probs=kfe_tp)
result['saved_policy_kfe'] = dict(current_l1_to_saved=float(np.abs(g_kfe-saved['g_current']).sum()),
    effective_max_to_saved=float(np.abs(kfe_tp-saved['tenure_probs']).max()),
    births=float(stats.total_births_kfe), births_gap_to_dated=abs(float(stats.total_births_kfe)-float(ev.births)),
    post_l1_to_dated=float(np.abs(stats.g_beginning_distribution-post).sum()),
    post_max_to_dated=float(np.abs(stats.g_beginning_distribution-post).max()))
assert smoke.digest(saved_path) == saved_sha
smoke.load_contract(a.contract,a.contract_sha256)
result['inputs_unchanged'] = True
report('complete')
