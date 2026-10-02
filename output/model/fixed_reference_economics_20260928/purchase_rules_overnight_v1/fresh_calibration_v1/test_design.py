"""Zero-native exact-loop and 24-slot source/seed contract check."""
import hashlib
import json
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
D = json.loads((HERE / 'design.json').read_text())
P = json.loads((PACKET / 'plan.json').read_text())
assert sha(PACKET / 'plan.json') == D['original_plan_sha256']
assert sha(PACKET / 'restart_controller_v2/controller.py') == D['reviewed_censor_source_sha256']
assert P['target_fingerprint'] == D['target_fingerprint']
assert P['weight_fingerprint'] == D['weight_fingerprint']
assert P['bounds'] == D['full_approved_bounds']
assert P['free_coordinates'] == D['free_coordinates']
assert len(D['starts']) == 24 and len(D['centers']) == 6
assert D['maximum_objective_calls_per_chain'] == 250
assert D['wall_seconds_per_chain'] == 14400 and D['final_reserve_seconds'] == 900
assert [s['slot'] for s in D['starts']] == list(range(24))
seen = set()
for s in D['starts']:
    slot = s['slot']
    arm = 'hard' if slot < 12 else 'quarter'
    chain = slot if slot < 12 else slot + 12
    assert (s['arm'], s['center_arm'], s['original_chain_index']) == (arm, arm, chain)
    assert s['weight_profile'] == 'base_control'
    assert any(c['arm'] == arm and c['chain'] == s['center_chain'] and c['origin'] == s['center_origin'] for c in D['centers'])
    assert set(s['parameters']) == set(P['free_coordinates'])
    assert all(P['bounds'][k][0] <= v <= P['bounds'][k][1] for k,v in s['parameters'].items())
    key = tuple(s['parameters'][k] for k in P['free_coordinates'])
    assert key not in seen
    seen.add(key)

sys.path.insert(0, str(PACKET))
sys.argv = ['run_psi.py', '--chain', '0']
import run_psi
from restart_controller_v2.controller import install_censored_optimizer

install_censored_optimizer(run_psi, 8)
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / 'search'
    calls = [0]
    def evaluate(label, point, deadline):
        calls[0] += 1
        if calls[0] == 2:
            return dict(status='budget_exhausted', reason='uncomputed_bounded_budget', lifecycle_solves=31)
        return dict(status='passed', residual=[float(calls[0])] + [0.] * 9,
                    lifecycle_solves=1, report=str(out / label))
    result = run_psi.optimize(out, {'x':.5}, {'x':(0.,1.)}, ('x',), evaluate,
                              time.time()+3600, toy=True)
    cases = json.loads((out / 'cases.json').read_text())
    assert len(cases) >= 3 and cases[0]['status'] == 'passed'
    assert cases[1]['optimizer_only_censor'] and cases[1]['loss'] is None
    assert cases[1]['computed_valid_loss'] is False
    assert cases[2]['status'] == 'passed' and result['selected']['status'] == 'passed'
    late = run_psi.optimize(Path(tmp) / 'deadline', {'x':.5}, {'x':(0.,1.)},
                            ('x',), evaluate, time.time()+1, toy=True)
    assert late['objective_calls'] == 0
    assert late['search_stop_reason'] == 'four_hour_actual_start_final_reserve'
print(json.dumps(dict(status='passed_zero_native_solves',starts=24,centers=6,
                      arms={'hard':12,'quarter':12},censored_case_continued=True,
                      deadline_stopped=True)))
