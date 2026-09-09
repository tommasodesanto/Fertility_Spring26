"""Bounded calendar-time forecast at a supplied stationary policy, not a PF root."""
import argparse
import gzip
import hashlib
import json
import pickle
from pathlib import Path
import sys
import time
from types import SimpleNamespace

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--contract', type=Path, required=True)
parser.add_argument('--contract-sha256', required=True)
parser.add_argument('--anchor-summary', type=Path, required=True)
parser.add_argument('--anchor-summary-sha256', required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--periods', type=int, required=True)
parser.add_argument('--seconds', type=int, required=True)
parser.add_argument('--initial-state', choices=('anchor', 'stationary'), default='anchor')
a = parser.parse_args()
assert 1 <= a.periods <= 100 and 1 <= a.seconds <= 240
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert digest(a.contract) == a.contract_sha256
c = json.loads(a.contract.read_text())
sys.path[:0] = [str(Path(c['source_root'])/'code/model/tools'),
               str(Path(c['source_root'])/'code/model')]
import numpy as np
import run_e5f_matched_pf_baseline as baseline
import run_e5f_perfect_foresight_person_demography as person
import run_e5f_perfect_foresight_person_demography_policy as checks
person.transition.configure_sequential_model()
person.calendar.apply_fertility = person.transition.apply_sequential_fertility
person.calendar.advance_calendar_distribution = person.transition.advance_sequential_calendar_distribution
person.calendar.distribution_rows = person.transition.independent_child_distribution_rows
baseline.joined.load_smoke_contract(a.contract, a.contract_sha256, 'sequential', maximum_seconds=7200)
assert digest(a.anchor_summary) == a.anchor_summary_sha256
anchor = json.loads(a.anchor_summary.read_text())
state_file = a.anchor_summary.parent/'final_state.npz'
assert anchor['arm'] == 'sequential' and anchor['mapping_valid']
assert digest(state_file) == anchor['artifact_sha256']['final_state.npz']
assert digest(c['terminal_checkpoint']) == c['terminal_checkpoint_sha256']
with gzip.open(c['terminal_checkpoint'], 'rb') as stream:
    terminal = pickle.load(stream)
P, grid, policy = terminal['parameters'], terminal['b_grid'], terminal['policy']
assert P.tau_H == .04 and P.property_tax_lump_sum_transfer == 0.
assert int(P.period_years) == 4 and not P.joint_nested_choice
z = np.load(state_file)
assert np.array_equal(z['wealth_grid'], grid)
state = person.PersonPFState(z['g_pre'], person.CohortState(int(z['year']), z['persons'], z['heads']))
state.persons.validated()
assert state.persons.year == 2119
demography = terminal['demographics']
price = float(policy.price[0])
reference = SimpleNamespace(asset_price=price, renter_price=P.user_cost_rate*price,
    equal_transfer=0., psi_child=P.psi_child,
    state=person.PersonPFState(terminal['fixed_point'].g_pre, terminal['fixed_point'].persons))
if a.initial_state == 'stationary':
    state = person.PersonPFState(reference.state.g_pre.copy(),
        person.CohortState(2119, reference.state.persons.persons.copy(), reference.state.persons.heads.copy()))
calendar, transition = person.calendar, person.transition
shared = calendar.model.precompute_shared(P, grid)
counter = calendar.SolveCounter()
shares = np.asarray(P.entry_shares, dtype=float).copy()
template = calendar.entrant_cohort(shares/shares.sum(), P, grid)
if a.output.exists() and any(a.output.iterdir()):
    raise FileExistsError(a.output)
a.output.mkdir(parents=True, exist_ok=True)
started = time.monotonic()
rows = []
def save(name, value):
    person.pf.write_json(a.output/name, value)
save('contract.json', dict(parent_contract_sha256=a.contract_sha256,
    anchor_summary_sha256=a.anchor_summary_sha256, final_state_sha256=digest(state_file),
    terminal_checkpoint_sha256=c['terminal_checkpoint_sha256'], script_sha256=digest(__file__),
    periods=a.periods, seconds=a.seconds, initial_state=a.initial_state,
    fixed_policy=True, market_clearing=False))
for step in range(a.periods + 1):
    view = SimpleNamespace(terminal_state=state, prices=[price], rents=[reference.renter_price],
                           rows=[{'equal_transfer_period_units': 0.}])
    distance = checks.terminal_convergence_diagnostics(view, terminal=reference, psi_path=[P.psi_child])
    row = dict(step=step, year=state.persons.year, elapsed_seconds=time.monotonic()-started,
               distance=distance, bellman_solves=counter.total)
    rows.append(row)
    save('latest_completed.json', row)
    save('history.json', rows)
    if a.initial_state == 'stationary':
        assert max(distance['metrics'][k] for k in (
            'resident_persons_relative_gap', 'household_heads_relative_gap',
            'normalized_household_distribution_l1', 'normalized_person_age_sex_l1',
            'normalized_head_age_sex_l1')) <= 1e-8
    if ((distance['all_checks_pass'] and a.initial_state != 'stationary')
            or step == a.periods or time.monotonic()-started >= a.seconds):
        break
    e = calendar.evaluate_period(np.array([price]), state.g_pre, P, grid, shared, counter,
        supply_rule=terminal['supply_rule'], supplied_policy=policy)
    births = transition.calendar_topcode_birth_accounting(e.g_pre, e.g_post_fertility, float(e.births), P)
    raw, _, _, raw_residual = transition.advance_sequential_calendar_distribution(e, np.zeros(int(P.I)), P, grid, shared)
    sex, survival, migration = demography.block_inputs(state.persons.year, 4)
    g, people, ledger = person.advance_household_person_block(raw, state.persons,
        model_period_births=float(births['topcode_adjusted_birth_children']), period_years=4,
        birth_sex_shares=sex, survival=survival, net_migration=migration,
        headship_rates=demography.headship_rates, model_age_start=int(P.age_start),
        model_age_cell_width=int(P.da), number_of_model_age_cells=int(P.J), age_axis=3,
        empty_age_templates={0: template})
    health = calendar.distribution_health(dict(pre=e.g_pre, post=e.g_post_fertility,
        current=e.g_current, raw=raw, next=g))
    assert counter.total == 0 and e.feasibility_projection_mass == 0.
    assert abs(float(raw_residual)) <= 2e-10
    assert abs(float(e.g_pre.sum()-e.g_post_fertility.sum())) <= 2e-10
    assert abs(float(e.g_post_fertility.sum()-e.g_current.sum())) <= 2e-10
    assert health['nonfinite_distribution_count'] == 0 and health['min_distribution_mass'] >= -1e-13
    assert ledger.person.person_identity_max_abs <= 2e-10
    assert ledger.person.head_identity_max_abs <= 2e-10
    assert abs(ledger.household_person_head_gap) <= 2e-10
    state = person.PersonPFState(g, people)
save('summary.json', dict(status='completed_bounded_fixed_policy_forecast',
    first_threshold_year=state.persons.year if distance['all_checks_pass'] else None,
    completed_periods=rows[-1]['step'], elapsed_seconds=time.monotonic()-started,
    bellman_solves=counter.total, final_distance=distance, equilibrium_certified=False))
