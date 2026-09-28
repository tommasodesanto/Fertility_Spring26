"""Torch-only fixed-policy spacing diagnostic; NOT a model solution/calibration.

Reuse the authenticated reference choices and native cohort propagation. After
one successful birth, optionally apply the saved subsequent-birth probability
at the destination state once more. Never chain three new births. No Bellman,
price, benefit-normalization, target or shared-source changes occur.
"""
import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
import numpy as np

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE = ROOT / 'output/model/fertility_identification_20260928'
OUT = BASE / 'two_births_v1'
REFERENCE = BASE / 'resume_v1/selected_export/primary'
sys.path.insert(0, str(ROOT / 'code/model/tools'))
import run_e5f_fertility_identification as driver
import e5f_evening_calibration_runtime as runtime


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def table(name, rows):
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fertility_step(pre, first_probability, continuation, pi, extra):
    """State axes b, tenure, location, income, children ever born, at home.

    All initial flows are computed from pre, not the changing output. A success
    from (n,m) moves to (n+1,m+1); an optional extra success moves to (n+2,m+2).
    continuation axes after fixing age: b, tenure, location, income, action,
    initial number of children minus one, children at home.
    """
    post = pre.copy()
    flows = np.zeros(3)
    extra_flows = np.zeros(3)
    for n in range(3):
        for m in range(n + 1):
            prob = first_probability if n == 0 else continuation[..., 1, n - 1, m]
            success = pre[..., n, m] * prob * pi
            second = np.zeros_like(success)
            if extra and n <= 1:
                second = success * continuation[..., 1, n, m + 1] * pi
            assert np.all(second >= 0) and np.all(second <= success + 1e-15)
            post[..., n, m] -= success
            post[..., n + 1, m + 1] += success - second
            flows[n] += float(success.sum())
            if n <= 1:
                post[..., n + 2, m + 2] += second
                flows[n + 1] += float(second.sum())
                extra_flows[n + 1] += float(second.sum())
    assert post.min() >= -1e-14 and np.isfinite(post).all()
    assert abs(post.sum() - pre.sum()) < 2e-12
    counts = np.arange(4)
    before = pre.sum(axis=(0, 1, 2, 3, 5))
    after = post.sum(axis=(0, 1, 2, 3, 5))
    assert abs((after - before) @ counts - flows.sum()) < 2e-12
    assert np.max(np.abs(np.cumsum(before - after)[:3] - flows)) < 2e-12
    return post, flows, extra_flows


def toy_checks():
    shape = (1, 1, 1, 1, 4, 4)
    first = np.ones(shape[:4])
    continuation = np.ones(shape[:4] + (2, 2, 4))
    for n in range(4):
        pre = np.zeros(shape)
        pre[..., n, n] = 1.
        for extra in (False, True):
            post, flows, more = fertility_step(pre, first, continuation, 1., extra)
            destination = min(n + (2 if extra else 1), 3)
            assert post[..., destination, destination].item() == 1.
            assert flows.sum() == destination - n
        post, flows, more = fertility_step(pre, first, continuation, 0., True)
        assert np.array_equal(post, pre) and flows.sum() == 0.
    pre = np.zeros(shape); pre[..., 0, 0] = 1.
    post, flows, more = fertility_step(pre, first, continuation, .5, True)
    assert np.array_equal(post.sum(axis=(0, 1, 2, 3, 5)), [.5, .25, .25, 0.])
    first[...] = 0.
    post, flows, more = fertility_step(pre, first, continuation, 1., True)
    assert np.array_equal(post, pre) and flows.sum() == 0.


def count_moments(counts):
    shares = counts / counts.sum()
    mean = float(shares @ np.arange(4))
    mothers = float(1 - shares[0])
    return dict(children_capped3=mean, mother_share=mothers,
                children_given_mother=mean / mothers, shares=shares.tolist())


def window_counts(pre, post, ages, lower, upper):
    left = np.maximum(ages, lower)
    right = np.minimum(ages + 4., upper)
    weights = np.maximum(right - left, 0.) / 4.
    fraction = np.clip(((left + right) / 2. - ages) / 4., 0., 1.)
    return (weights[:, None] * ((1-fraction[:, None])*pre + fraction[:, None]*post)).sum(axis=0)


started = time.monotonic()
toy_checks()
print('Toy mass, birth-count, no-third-newbirth, zero-success and no-attempt checks PASS', flush=True)
contract, objectives = driver.verify(BASE / 'contract_v1/contract.json')
hashes = json.loads((REFERENCE / 'artifact_hashes.json').read_text())
for name, digest in hashes.items():
    assert sha(REFERENCE / name) == digest, name
receipt = json.loads((REFERENCE / 'receipt.json').read_text())
assert sha(REFERENCE / 'initial_state.pkl.gz') == receipt['case_checkpoint_sha256']
evaluator = runtime.setup(dict(contract, objective=contract['lanes']['primary']['objective']),
                          objectives['primary'], OUT / ('runtime_' + os.environ['SLURM_JOB_ID']))
with gzip.open(REFERENCE / 'initial_state.pkl.gz', 'rb') as stream:
    packet = pickle.load(stream)
P = copy.deepcopy(packet['parameters'])
sol = packet['solution']; grid = packet['b_grid']; shared = packet['shared']
model = evaluator.rt['model']; calendar = evaluator.rt['primitive'].pf.calendar
assert P.sequential_births and not P.joint_nested_choice
assert model.independent_child_maturation_active(P)
assert not model.parent_age_maturation_active(P) and not model.readiness_gate_active(P)
assert P.n_parity == P.n_child_states == 4 and P.period_years == 4.
policy = calendar.policy_from_solution(sol, sol.p_eq, P, grid, shared)
assert policy.bp_pol_stay is not None and P.native_due_stayer_credit
for probabilities in (policy.fert_probs, policy.fert2_probs):
    assert np.isfinite(probabilities).all()
    assert probabilities.min() >= 0. and probabilities.max() <= 1.
_, _, Pi_z = model.income_transition_values(P)
fec = model.get_fecundity_by_age(P)
ages = P.age_start + np.arange(P.J) * P.da
saved_pre = packet['evaluation'].g_pre
saved_post = sol.g_beginning_distribution
results = {}

for name, extra in [('control', False), ('two_births_fixed_policy', True)]:
    current = saved_pre[:, :, :, 0].copy()
    pre_counts = []; post_counts = []; birth_flows = []; extra_flows = []
    max_pre_gap = max_post_gap = max_mass_gap = max_dead_mass = 0.
    for j in range(P.J):
        assert time.monotonic() - started < 840, '14-minute internal budget'
        pre_counts.append(current.sum(axis=(0, 1, 2, 3, 5)))
        if name == 'control':
            max_pre_gap = max(max_pre_gap, float(np.max(np.abs(current - saved_pre[:, :, :, j]))))
        if P.A_f_start <= j + 1 <= P.A_f_end:
            post, flows, more = fertility_step(current, policy.fert_probs[:, :, :, j, :, 1],
                policy.fert2_probs[:, :, :, j], float(fec[j]), extra)
        else:
            post = current.copy(); flows = np.zeros(3); more = np.zeros(3)
        post_counts.append(post.sum(axis=(0, 1, 2, 3, 5)))
        birth_flows.append(flows); extra_flows.append(more)
        dead_mass = float(post[policy.V[:, :, :, j] <= model.DEAD_VALUE_CUTOFF].sum())
        max_dead_mass = max(max_dead_mass, dead_mass)
        assert dead_mass <= model.DEAD_MASS_TOL, (name, j, 'dead-state exposure', dead_mass)
        if name == 'control':
            max_post_gap = max(max_post_gap, float(np.max(np.abs(post - saved_post[:, :, :, j]))))
        print(f'{name} age {ages[j]:g}: count={count_moments(post_counts[-1])["children_capped3"]:.6f}', flush=True)
        if j + 1 < P.J:
            survival = float(P.survival_probs[j]) if P.use_age_survival else 1.
            current = model.advance_cohort_one_period_markov_income(
                survival * post, j, policy.loc_probs, policy.tenure_choice, policy.tenure_probs,
                policy.bp_pol, P, grid, shared, policy.maps.lmm_idx, policy.maps.lmm_wt,
                policy.maps.tmx_idx, policy.maps.tmx_wt, bool(P.use_stochastic_aging),
                P.Pi_child if P.use_stochastic_aging else None, Pi_z,
                bp_pol_stay=policy.bp_pol_stay)
            mass_gap = abs(float(current.sum()) - survival * float(post.sum()))
            max_mass_gap = max(max_mass_gap, mass_gap)
            assert mass_gap < 2e-10, (name, j, 'cohort mass', mass_gap)
    pre = np.array(pre_counts); post = np.array(post_counts)
    flows = np.array(birth_flows); more = np.array(extra_flows)
    if name == 'control':
        assert max_pre_gap < 5e-9 and max_post_gap < 5e-9
        assert abs(flows.sum() - packet['evaluation'].births) < 5e-9
        for column, field in enumerate(('_first_births_by_age', '_second_births_by_age', '_third_births_by_age')):
            assert np.max(np.abs(flows[:, column] - np.asarray(getattr(P, field)))) < 5e-9
    else:
        assert np.max(np.abs(flows[:, 0] - np.array(results['control']['birth_flows'])[:, 0])) < 5e-9
        assert np.max(np.abs(pre[:, 0] - np.array(results['control']['pre_counts'])[:, 0])) < 5e-9
    terminal = count_moments(post[6])
    terminal['model_coded_children'] = float(np.array(terminal['shares']) @ [0, 1, 2, P.tfr_top_bin_weight])
    adjusted_births = model.adjusted_births(float(flows.sum()), float(flows[:, 2].sum()), P.tfr_top_bin_weight)
    entry = float(pre[0].sum())
    potential_entry = model.potential_entry_households(adjusted_births)
    results[name] = dict(age25=count_moments(window_counts(pre, post, ages, 25., 26.)),
        ages40_44=count_moments(window_counts(pre, post, ages, 40., 45.)),
        terminal_post_fertility=terminal,
        mapped_mean_first_birth_age=float(flows[:, 0] @ (ages + 2.) / flows[:, 0].sum()),
        first_birth_cell_shares=(flows[:, 0]/flows[:, 0].sum()).tolist(),
        pre_counts=pre.tolist(), post_counts=post.tolist(), birth_flows=flows.tolist(),
        extra_birth_flows=more.tolist(), max_control_pre_gap=max_pre_gap,
        max_control_post_gap=max_post_gap, max_cohort_mass_gap=max_mass_gap,
        max_dead_state_exposure=max_dead_mass,
        renewal_diagnostic=dict(entry_households=entry, potential_entry_households=potential_entry,
            adjusted_births=adjusted_births, births_per_entry=adjusted_births/entry,
            fertility_gap=adjusted_births/entry-2.1, entry_residual=entry-potential_entry,
            potential_entry_percent_above_fixed_entry=100*(potential_entry/entry-1)))
    write('latest_completed.json', dict(completed=name, results=results))

original = results['control']; diagnostic = results['two_births_fixed_policy']
assert abs(original['age25']['children_capped3'] - .5354261243026106) < 5e-9
assert abs(original['mapped_mean_first_birth_age'] - 25.93278056306455) < 5e-9
rows = []
for lower in [20, 25, 30, 35, 40]:
    row = dict(age_window=f'{lower}-{lower+4}')
    for name, result in results.items():
        moments = count_moments(window_counts(np.array(result['pre_counts']), np.array(result['post_counts']), ages, lower, lower+5))
        row.update({f'{name}_{key}':value for key,value in moments.items() if key != 'shares'})
    rows.append(row)
table('lifecycle.csv', rows)
write('result.json', dict(status='passed_fixed_policy_diagnostic_not_equilibrium',
    label='2007 stationary reference — block0506, September 28 verified export',
    reference_checkpoint_sha256=receipt['case_checkpoint_sha256'],
    contract_sha256=sha(BASE/'contract_v1/contract.json'), script_sha256=sha(__file__),
    standard_plots_preserved_and_hash_verified=17, model_solves=0, cohort_replays=2,
    slurm_job=os.environ['SLURM_JOB_ID'], elapsed_seconds=time.monotonic()-started,
    economic_changes=['Experimental extra subsequent-birth opportunity after one successful birth; at most two new births per cell.'],
    retained=['Saved birth/tenure/location/saving choices, income, entry, prices, preferences and psi.',
              'Same age-window projection; no new within-cell birth dates.', 'All calibration targets/weights/bounds and original 17 plots.'],
    limitations=['No household reoptimization or market/fiscal/renewal clearing.',
        'Completed fertility 2.1 is not imposed in this fixed-policy accounting diagnostic; it remains mandatory for a calibration.',
        'Age25 stock uses the inherited linear pre/post interpolation; distinct within-cell birth dates are unspecified.',
        'No full-fit loss is scored because this is not a model solution.'], results=results))
print('PASS: authenticated baseline replay and experimental fixed-policy replay; zero model solves', flush=True)
