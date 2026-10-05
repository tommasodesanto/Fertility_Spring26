"""Execution-only checks for the isolated fixed-psi candidate driver."""
from __future__ import annotations
import contextlib
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import pytest

SOURCE = Path(__file__).with_name('parallel_transition_candidate.py')
spec = importlib.util.spec_from_file_location('parallel_transition_candidate', SOURCE)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_no_arbitrage_initial_path_matches_recursion():
    P = SimpleNamespace(R_gross=1.08243216, delta=.05545379079326218,
                        tau_H=.042393443095490375)
    for h in (24, 32):
        q = runner.no_arbitrage_prices(P, .65, h)
        A = P.R_gross + P.delta + P.tau_H
        assert len(q) == h
        assert np.allclose(A*q[:-1]-runner.RENT, q[1:], rtol=0, atol=2e-16)
        assert q[-1] == pytest.approx((.65+runner.RENT)/A)
    with pytest.raises(ValueError, match='24 or 32'):
        runner.no_arbitrage_prices(P, .65, 6)


def fake_contract(tmp_path, monkeypatch, *, native_gate=True, bad_queue=False, bad_source=False,
                  fail_second_mapping=False, fail_first_mapping=False):
    source = tmp_path/'frozen'; source.mkdir(parents=True)
    for name in ('two_shock.py','two_shock_runtime.py'):
        (source/name).write_text(name)
    evidence = tmp_path/'evidence.json'
    evidence.write_text(json.dumps(dict(accounting_valid=True,gates=dict(mass=True,projection=True,
                                         policy_reproduction=True,dated_audits=True))))
    state = SimpleNamespace(g_pre=np.array([.4,.6]),scheduled_entries=np.array([.1,.2]),
                            scheduled_raw_entries=np.array([.1,.2]))
    def state_hash(s, queue_values):
        arrays=[np.asarray(s.g_pre),np.asarray(queue_values(s.scheduled_entries)),
                np.asarray(queue_values(s.scheduled_raw_entries))]
        if bad_queue:
            raise ValueError('Invalid raw queue')
        return runner.hashlib.sha256(b''.join(a.tobytes() for a in arrays)).hexdigest()
    driver = SimpleNamespace(state_hash=state_hash,pinned=runner.pinned)
    interrupted_mapping_number=1 if fail_first_mapping else 2 if fail_second_mapping else None
    class RootSolver:
        @staticmethod
        def solve_joint_with_acceleration(**kwargs):
            q=np.asarray(kwargs['initial_prices']); b=np.asarray(kwargs['initial_fiscal_values'])
            assert kwargs['max_evaluations']==12
            assert kwargs['market_tolerance']==2e-4 and kwargs['fiscal_tolerance']==2e-5
            assert kwargs['final_reproduction_tolerance']==1e-10
            assert kwargs['initial_jacobian'].shape==(2*len(q),2*len(q))
            assert np.array_equal(q,runner.no_arbitrage_prices(P,.65,len(q)))
            budget_result=dict(status='time_or_evaluation_budget',converged=False,
                gates=dict(market_replay=False,fiscal_replay=False),
                final_reproduction_max_abs=None,final=None)
            if fail_first_mapping:
                try:
                    kwargs['evaluate'](q,b)
                except TimeoutError:
                    return budget_result
            else:
                kwargs['evaluate'](q,b)
            if fail_second_mapping:
                try:
                    kwargs['evaluate'](q,b) # interrupted reserved fresh replay
                except TimeoutError:
                    return budget_result
            else:
                kwargs['evaluate'](q,b) # original reserved fresh replay
            return dict(converged=True,gates=dict(market_replay=True,fiscal_replay=True),
                        final_reproduction_max_abs=0.,final=dict(prices=q,fiscal_values=b,
                        market_residual=np.zeros(len(q)),fiscal_residual=np.zeros(len(q))))
        @staticmethod
        def extend_measured_jacobian(seed,horizon):
            assert seed['mapping_count']==5 and seed['horizon']==12
            return np.eye(2*horizon)
    monkeypatch.setitem(sys.modules,'e5f_four_shock_acceleration',RootSolver)
    original=SimpleNamespace(inner=SimpleNamespace(write=runner.write_full))
    retained=SimpleNamespace(original_modules=lambda:(original,None))
    P=SimpleNamespace(R_gross=1.08243216,delta=.05545379079326218,
                      tau_H=.042393443095490375)
    class BaseAdapter:
        def __init__(self,rt,plan,inherited_state,start_year,initialization=None):
            self.rt,self.plan,self.inherited_state,self.start_year=rt,plan,inherited_state,start_year
            self.calls=0;self.q=.7;self.pension=.2
        def identity(self):
            return dict(self.rt.identity(),stage_start_year=2007,
                        inherited_state_sha256=state_hash(self.inherited_state,lambda x:x))
        def _endpoint(self,psi,folder,deadline):
            self.calls+=1;self.rt.total_native_calls+=1
            folder.mkdir(parents=True,exist_ok=True)
            runner.write_full(folder/'endpoint.json',dict(psi=psi,price=.65,pension=.2))
            return dict(parameters=SimpleNamespace(pension=.2)),dict(price=.65,stationary_pass=True,
                                                                      stationary_renewal_gap=0.)
        def select_warm(self,horizon,psi): return None,'fresh_measured_seed_initialization'
        def retain_warm(self,*args): pass
        def _mapping(self,terminal,endpoint,q,b,psi_path,folder,deadline):
            if interrupted_mapping_number == self.calls:
                self.rt.total_native_calls+=1 # native entry began, but _account never observed completion
                raise TimeoutError('mapping interrupted by budget')
            self.calls+=1;self.rt.total_native_calls+=1
            folder.mkdir(parents=True,exist_ok=True)
            h=len(q);rows=[dict(calendar_year=2007+4*i,asset_price=float(q[i]),
                              pension_period_units=float(b[i])) for i in range(h)]
            fertility=[dict(period_tfr_topcode_adjusted=1.8) for _ in range(h)]
            record=dict(accounting_valid=True,gates=dict(mass=native_gate,projection=True,
                policy_reproduction=True,dated_audits=True),market_residual=np.zeros(h),
                fiscal_residual=np.zeros(h),rows=rows,fertility=fertility)
            runner.write_full(folder/'native_record.json',record)
            dated={i:dict(state=state) for i in range(h)}
            native_reply=SimpleNamespace(dated_states=dated,values=[np.ones(2) for _ in range(h+1)],
                floor_runtime_paths=dict(prices=np.asarray(q),pensions=np.asarray(b),
                                         psi_path=np.asarray(psi_path),start_year=2007))
            return native_reply,record
    class FakeRT:
        def __init__(self):self.total_native_calls=0;self.P=P;self.initial_state=state;self.housing='static-elastic'
        def identity(self):return dict(source_pins={'fake':'sha'},reference_sha256='fresh')
        def terminal_checks(self,*args,**kwargs):return dict(all_checks_pass=False)
    class FakeRuntime:
        def __init__(self,plan,output):
            self.plan=plan;self.rt=FakeRT();self.adapters={};self.seeds={};self.queue_values=lambda x:x
            self.fitter=SimpleNamespace(fit_one=lambda **kw:pytest.fail('scalar fitter called'))
        def bind_budget(self,deadline,heartbeat):self.heartbeat=heartbeat
        @contextlib.contextmanager
        def budget_context(self):yield
        def prepare_reference(self,deadline,folder):
            self.rt.total_native_calls+=2;self.heartbeat('native_call',actual_native_calls=2)
            return dict(policy_calls=2,reference_checkpoint=runner.pin(evidence),
                        reconstruction_receipt=runner.pin(evidence))
        def initial_state(self):return state
        def measure_seed(self,**kwargs):
            self.rt.total_native_calls+=5;self.heartbeat('native_call',actual_native_calls=7)
            self.adapters[0]=native.StageAdapter(self.rt,self.plan,state,2007)
            seed=dict(mapping_count=5,horizon=12,identity=self.adapters[0].identity(),
                      source_evidence=[runner.pin(evidence)]*5,accounting_valid=True,
                      policy_calls=5)
            self.seeds[0]=seed;return seed
        def evaluate_stage(self,**kwargs):
            adapter=self.adapters[0]
            reply=adapter.evaluate(psi=kwargs['psi'],start_year=2007,horizon=kwargs['horizon'],
                seed=self.seeds[0],gates=self.plan['gates'],budget=self.plan['budget'],
                endpoint_controls={},path_controls=self.plan['path'],deadline=kwargs['deadline'],
                folder=kwargs['folder'])
            raw=reply['native_reply'];paths=raw.floor_runtime_paths
            reply.update(prices=paths['prices'],pensions=paths['pensions'],values=raw.values,
                dated_states=raw.dated_states,source_evidence=[runner.pin(kwargs['folder']/'root.json')]+[
                    runner.pin(p) for p in sorted(kwargs['folder'].rglob('native_record.json'))])
            return reply
        def export_state(self,reply,*,index,year,folder):
            folder.mkdir(parents=True,exist_ok=True)
            return dict(calendar_year=2023,period_index=4,exact_native_state=True,
                        reconstructed_or_rescaled=False,path=str(folder/'actual_2023.pkl.gz'))
        def render_standard(self,reply,folder,deadline):
            plots={}
            directory=folder/'date_000'/'standard_diagnostics';directory.mkdir(parents=True)
            for name in self.plan['standard_plot_names']:
                path=directory/name;path.write_bytes(b'fake PNG')
                plots[name]=runner.sha(path)
            return dict(sampled_dates=[dict(period=0,plots=plots)])
    native=SimpleNamespace(StageAdapter=BaseAdapter,NativeRuntime=FakeRuntime,retained=retained,
                           write=runner.write_full,pin=runner.pin)
    plan=dict(identity=dict(source_pins={'fake':'sha'},reference_sha256='fresh'),
              seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),
              gates=dict(market_tolerance=2e-4,fiscal_tolerance=2e-5,
                         final_reproduction_tolerance=1e-10,stationary_renewal_tolerance=1e-6),
              budget=dict(seed_seconds=1800,candidate_seconds=7200,path_seconds=6000,
                          endpoint_seconds=1800,mapping_seconds=1800,render_seconds=900,
                          maximum_policy_calls=20000),
              path=dict(max_evaluations=12,price_bound_ratios=[.05,20.],
                        pension_bound_ratios=[.05,20.],max_log_step=.15,damping=.7),
              fit=dict(max_condition_number=1e8,worsening_factor=2.),mode='diagnostic',
              standard_plot_names=[f'plot_{i}.png' for i in range(17)])
    if bad_source: plan['identity']['source_pins']={'wrong':'sha'}
    manifest=tmp_path/'fit_manifest.json';manifest.write_text('{}')
    monkeypatch.setattr(runner,'load_frozen',lambda config:(plan,driver,native,manifest,source))
    return dict(driver=runner.pin(SOURCE),fit_manifest=runner.pin(manifest),
                frozen_package=str(source),task_id='c00_h24',psi=.14736308634876963,
                horizon=24,budget=dict(external_seconds=10800,internal_seconds=10700,
                                       maximum_policy_calls=20000))


def test_full_fake_loop_preserves_order_and_nonacceptance(tmp_path,monkeypatch):
    config=fake_contract(tmp_path,monkeypatch,native_gate=False)
    result=runner.run(config,tmp_path/'run')
    assert result['status']=='candidate_unaccepted' # terminal diagnostic is nongating;
    assert result['root_pass'] and result['replay_pass'] and not result['accounting_valid']
    assert result['native_calls']==10 and result['seed_native_calls']==5
    assert result['scalar_fitter_called'] is False
    assert result['own_vintage_2015']['local_index']==2
    assert result['own_vintage_2015']['selected_stage_handoff'] is False
    assert len(result['first_two_historical_fertility_rows'])==2
    assert len(result['full_dated_rows'])==24
    assert (tmp_path/'run/candidate/warm_start.json').is_file()
    warm=json.loads((tmp_path/'run/candidate/warm_start.json').read_text())
    assert warm['kind']=='no_arbitrage_backward_from_fresh_endpoint'
    assert (tmp_path/'run/candidate/root_full.json').is_file()
    assert (tmp_path/'run/latest_completed.json').is_file()
    assert not (tmp_path/'run/best_so_far.json').exists()


def test_full_fake_loop_accepts_only_original_root_gates(tmp_path,monkeypatch):
    config=fake_contract(tmp_path,monkeypatch)
    result=runner.run(config,tmp_path/'accepted')
    assert result['status']=='root_gate_passed_pending_horizon'
    assert result['accepted_candidate'] is False and result['pending_horizon_comparison'] is True
    assert result['terminal_diagnostic_pass'] is False
    assert result['terminal_diagnostic_gating'] is False
    assert len(result['standard_diagnostics']['sampled_dates'][0]['plots'])==17
    assert (tmp_path/'accepted/best_so_far.json').is_file()
    assert result['own_vintage_2015']['selected_stage_handoff'] is False
    progress=json.loads((tmp_path/'accepted/progress.json').read_text())
    assert progress['accepted_candidate'] is False
    assert progress['root_gate_passed'] is True
    assert progress['pending_horizon_comparison'] is True


def test_interrupted_replay_pins_latest_successful_mapping_and_preserves_root_outcome(tmp_path,monkeypatch):
    config=fake_contract(tmp_path,monkeypatch,fail_second_mapping=True)
    result=runner.run(config,tmp_path/'interrupted')
    assert result['status']=='candidate_unaccepted'
    assert result['root_pass'] is False
    assert result['replay_pass'] is False
    assert result['root']['status']=='time_or_evaluation_budget'
    assert result['root']['final'] is None
    assert result['path_evaluations']==2  # attempted path maps
    assert result['path_mappings_completed']==1
    assert result['latest_completed_map_number']==1
    pinned=runner.pinned(result['final_mapping_pin'])
    assert pinned == (tmp_path/'interrupted/candidate/map_001/native_record.json').resolve()
    assert result['replay_maximum_gap'] is None
    assert result['candidate_native_calls']==3 # endpoint, successful map, interrupted native entry
    assert result['candidate_completed_operation_calls']==2
    assert result['candidate_unfinished_native_calls']==1
    saved=json.loads((tmp_path/'interrupted/latest_completed.json').read_text())
    assert saved['root']['status']=='time_or_evaluation_budget'
    assert saved['root']['final'] is None
    assert saved['final_mapping_pin']==result['final_mapping_pin']


def test_interrupted_first_mapping_still_fails_without_latest_completed_map(tmp_path,monkeypatch):
    config=fake_contract(tmp_path,monkeypatch,fail_first_mapping=True)
    with pytest.raises(ValueError,match='Original root produced no completed native mapping'):
        runner.run(config,tmp_path/'no_completed_map')
    failure=json.loads((tmp_path/'no_completed_map/failure.json').read_text())
    assert failure['status']=='failed'


def test_fake_queue_and_source_failures_leave_receipt(tmp_path,monkeypatch):
    config=fake_contract(tmp_path,monkeypatch,bad_source=True)
    with pytest.raises(ValueError,match='seed identity'):
        runner.run(config,tmp_path/'bad_source')
    assert (tmp_path/'bad_source/failure.json').is_file()
    assert not (tmp_path/'bad_source/result.json').exists()
    config=fake_contract(tmp_path/'queue',monkeypatch,bad_queue=True)
    with pytest.raises(ValueError,match='Invalid raw queue'):
        runner.run(config,tmp_path/'bad_queue')
    assert (tmp_path/'bad_queue/failure.json').is_file()


def test_pin_mutation_rejected(tmp_path):
    path=tmp_path/'source.py';path.write_text('original')
    proof=runner.pin(path)
    path.write_text('changed')
    with pytest.raises(ValueError,match='Pinned file differs'):
        runner.pinned(proof)
