"""Deterministic exact-loop controller tests; zero native model solves."""
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
import difflib
from types import SimpleNamespace
from contextlib import contextmanager
from unittest.mock import patch

import one_shock_floor as runner


def stationary_gates():
    return dict(household_budget={},purchase={},estate={},policy_arrays={},stationary_operator={},
        housing_market_clearing_required=False,feasibility_projection_mass=0.,
        fiscal_certificate=dict(marginal_gate=True,fiscal_gate=True,marginal_tolerance=1e-9,fiscal_tolerance=1e-6))


def plan():
    original,_ = runner.original_modules()
    root = runner.HERE.parents[3]
    blocks = root/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/inputs/empirical_blocks.csv'
    annual = root/'output/model/e5f_matched_pf_20260909a/path_pilot_20260910/fertility_data/annual_fertility_2007_2023.csv'
    def pin(path):
        return dict(path=str(path),sha256=runner.sha(path))
    identity = {key:'a'*64 for key in runner.IDENTITY_KEYS-{'source_pins'}}
    identity['reference_sha256']=runner.sha(__file__)
    identity['source_pins'] = {'test_runtime':'a'*64}
    return dict(schema='current_floor_one_permanent_v1',kind='one_permanent',start_year=2007,
        fiscal_relaxation_authorized=True,gates=runner.GATES.copy(),mode='diagnostic',horizons=[6,8],identity=identity,handoff=pin(Path(__file__)),
        source_files=dict(controller=pin(runner.HERE/'one_shock_floor.py'),
            original_estimator=pin(runner.PINNED/'run_e5f_preference_estimation.py'),
            original_fitter=pin(runner.PINNED/'e5f_preference_shock_fit.py'),runtime=pin(Path(__file__))),
        target_contract=original.target_contract(blocks,annual),
        budget=dict(total_seconds=30,seed_seconds=2,candidate_seconds=5,endpoint_seconds=2,
            mapping_seconds=1,path_seconds=2,render_seconds=1,maximum_policy_calls=200),
        seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),initial_psi=.2,psi_bound_ratios=[.01,2.],
        standard_plot_names=[f'plot_{i:02d}.png' for i in range(17)],
        fit=dict(max_evaluations=6,log_difference_step=.01,fertility_tolerance=.005,max_log_step=.15,
            damping=.7,max_condition_number=1e8,worsening_factor=1.5,reproduction_tolerance=1e-8),
        path=dict(max_evaluations=2,price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],max_log_step=.15,damping=.7),
        endpoint=dict(max_evaluations=2,price_bound_ratios=[.05,20.],max_log_step=.15,damping=.7,slope=1.))


class FakeRuntime:
    def __init__(self,p):
        self.p = p
        self.calls = []
        self.accounting_valid = True
        self.seed_identity = p['identity']
        self.bad_export = False

    def identity(self):
        return self.p['identity']

    def measure_seed(self,**kwargs):
        self.calls.append(('seed',kwargs['horizon']))
        return dict(identity=self.seed_identity,horizon=12,mapping_count=5,
            unknown_blocks=['log_house_price','log_period_pension'],native_measured=True,
            matrix=[[1]],accounting_valid=True,policy_calls=60)

    def evaluate(self,**kwargs):
        psi,H = kwargs['psi'],kwargs['horizon']
        self.calls.append(('evaluate',psi,H))
        rows = [dict(calendar_year=2007+4*i,asset_price=1.,renter_price=.1,adult_population=1.,
            birth_children=1.,housing_demand=1.,pension_period_units=.1) for i in range(H)]
        target = self.p['target_contract']['rows'][3]['target']
        fertility = [dict(period_tfr_topcode_adjusted=target+.1*math.log(psi/.2)) for _ in range(H)]
        return dict(identity=self.identity(),psi=psi,horizon=H,accounting_valid=self.accounting_valid,policy_calls=1,
            root_and_terminal_pass=True,root_pass=True,replay_pass=True,terminal_pass=True,stationary_pass=True,stationary_renewal_gap=0.,
            market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,rows=rows,fertility=fertility,
            reference_manifest_sha256=self.identity()['reference_sha256'],source_pins=self.identity()['source_pins'],
            housing='static-elastic',shock_contract=dict(psi=psi,start_year=2007))

    def export_2023(self,reply,folder):
        self.calls.append(('export',reply['psi'],reply['horizon']))
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
        path=folder/'actual_2023.pkl.gz';path.write_bytes(b'fake native state for output-only controller test')
        return dict(calendar_year=2023,period_index=4,exact_native_state=True,
            reconstructed_or_rescaled=self.bad_export,queue_lags=[16,20],forecast_and_continuation_saved=True,
            path=str(path),sha256=runner.sha(path))

    def render_standard(self,reply,folder,deadline):
        folder.mkdir(parents=True)
        paths = [folder/name for name in self.p['standard_plot_names']]
        for path in paths:
            path.write_bytes(b'fake controller smoke plot, not native visual evidence')
        return paths


class Tests(unittest.TestCase):
    def test_stationary_audit_zeros_and_descriptive_false_are_healthy(self):
        record=dict(accounting_valid=True,gates=stationary_gates())
        self.assertFalse(all(record['gates'].values()))
        self.assertTrue(runner.stationary_mapping_valid(record))
        for field in ('marginal_gate','fiscal_gate'):
            bad=copy.deepcopy(record);bad['gates']['fiscal_certificate'][field]=False
            self.assertFalse(runner.stationary_mapping_valid(bad))
        bad=copy.deepcopy(record);bad['gates']['feasibility_projection_mass']=1e-12
        self.assertFalse(runner.stationary_mapping_valid(bad))
        bad=copy.deepcopy(record);bad['accounting_valid']=False
        self.assertFalse(runner.stationary_mapping_valid(bad))
        bad=copy.deepcopy(record);bad['gates']['fiscal_certificate']['fiscal_tolerance']=2e-5
        with self.assertRaisesRegex(ValueError,'Original native stationary'):
            runner.stationary_mapping_valid(bad)
    def compatibility_fixture(self,directory):
        directory=Path(directory);p=plan()
        old=runner.HERE.parents[3]/'output/model/transition_readiness_v1/current_floor_preparation/controller/controller_before_checkpoint.py'
        snapshot=directory/'historical_generator.py';snapshot.write_bytes(old.read_bytes())
        pin=lambda path:dict(path=str(path),sha256=runner.sha(path))
        prep={k:copy.deepcopy(p[k]) for k in ('identity','source_files','target_contract','seed','gates')}
        prep['source_files']['controller']['sha256']=runner.sha(snapshot)
        original=directory/'original_preparation.json';runner.write(original,prep)
        p['prepared_native_inputs']=pin(original)
        old_pin=prep['source_files']['controller'];consumer=p['source_files']['controller']
        diff=directory/'reviewed.diff'
        diff.write_text(''.join(difflib.unified_diff(snapshot.read_text().splitlines(True),Path(consumer['path']).read_text().splitlines(True),fromfile=old_pin['path'],tofile=consumer['path'])))
        scopes=['output_checkpoint','explicit_generation_consumer_validation','stationary_gate_boolean_normalization']
        audit=directory/'reviewed_audit.json'
        runner.write(audit,dict(schema='controller_checkpoint_consumer_source_delta_audit_v1',accepted_numerics_unchanged=True,
            lead_review_accepted=True,original_preparation=pin(original),old_snapshot=pin(snapshot),old_source=old_pin,
            new_source=consumer,diff=pin(diff),change_scopes=scopes,
            checks=runner.compatibility_ast_checks(snapshot.read_text(),Path(consumer['path']).read_text())))
        receipt=dict(schema='current_floor_prepared_consumer_compatibility_v1',compatibility_approved=True,
            lead_review_accepted=True,original_preparation=pin(original),original_generator_controller=old_pin,
            original_generator_snapshot=pin(snapshot),consumer_controller=consumer,allowed_source_roles=['controller'],
            reviewed_audit=pin(audit),reviewed_diff=pin(diff),change_scopes=scopes)
        receipt_path=directory/'compatibility.json';runner.write(receipt_path,receipt)
        p['prepared_consumer_compatibility']=pin(receipt_path)
        return p,prep,receipt,receipt_path

    def test_explicit_approved_consumer_compatibility_keeps_original_generator_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            p,prep,receipt,path=self.compatibility_fixture(directory)
            original_digest=runner.sha(p['prepared_native_inputs']['path'])
            runner.validate_consumer_compatibility(p,prep)
            self.assertEqual(runner.sha(p['prepared_native_inputs']['path']),original_digest)
            self.assertNotEqual(prep['source_files']['controller']['sha256'],runner.sha(prep['source_files']['controller']['path']))

    def test_consumer_compatibility_fails_without_approval_or_with_stale_diff(self):
        with tempfile.TemporaryDirectory() as directory:
            p,prep,receipt,path=self.compatibility_fixture(directory)
            receipt['compatibility_approved']=False;runner.write(path,receipt)
            p['prepared_consumer_compatibility']['sha256']=runner.sha(path)
            with self.assertRaisesRegex(ValueError,'explicit lead-reviewed'):
                runner.validate_consumer_compatibility(p,prep)
            receipt['compatibility_approved']=True;runner.write(path,receipt)
            p['prepared_consumer_compatibility']['sha256']=runner.sha(path)
            Path(receipt['reviewed_diff']['path']).write_text('stale diff')
            with self.assertRaisesRegex(ValueError,'Missing or changed pin'):
                runner.validate_consumer_compatibility(p,prep)

    def test_compatibility_rejects_multiple_source_changes_and_missing_pin(self):
        with tempfile.TemporaryDirectory() as directory:
            p,prep,receipt,path=self.compatibility_fixture(directory)
            original=copy.deepcopy(prep)
            prep.update(schema='current_floor_native_reference_seed_v1',reference_repeat_verified=True,fresh_seed_baseline_verified=True)
            p.pop('prepared_consumer_compatibility')
            with self.assertRaisesRegex(ValueError,'compatibility pin required'):
                runner.validate_preparation(p,prep)
            prep['source_files']['runtime']['sha256']='b'*64
            with self.assertRaisesRegex(ValueError,'Only the reviewed controller'):
                runner.validate_preparation(p,prep)

    def test_candidate_checkpoint_exports_once_without_extra_evaluation(self):
        p=plan();rt=FakeRuntime(p)
        with tempfile.TemporaryDirectory() as directory:
            c=runner.Controller(p,rt,directory);c.prepare();reply=c.evaluate(.2)
            self.assertEqual(len([r for r in rt.calls if r[0]=='evaluate']),2)
            self.assertEqual(len([r for r in rt.calls if r[0]=='export']),1)
            checkpoint=json.loads(runner.pinned(reply['checkpoint_2023']).read_text())
            self.assertFalse(checkpoint['shock_fit_complete'])
            self.assertTrue(checkpoint['one_shock_measurement_certified'])
            self.assertFalse(checkpoint['state_experiment_ready'])
            self.assertEqual(json.loads((Path(directory)/'best_so_far.json').read_text())['checkpoint_2023'],reply['checkpoint_2023'])

    def test_checkpoint_timeout_never_marks_partial_candidate_accepted(self):
        p=plan()
        class PartialExport(FakeRuntime):
            def export_2023(self,reply,folder):
                super().export_2023(reply,folder)
                raise TimeoutError('checkpoint interrupted after partial output')
        rt=PartialExport(p)
        with tempfile.TemporaryDirectory() as directory:
            c=runner.Controller(p,rt,directory);c.prepare()
            with self.assertRaisesRegex(TimeoutError,'checkpoint interrupted'):
                c.evaluate(.2)
            self.assertFalse((Path(directory)/'latest_completed.json').exists())
            self.assertFalse((Path(directory)/'candidate_0001/state_2023_checkpoint/checkpoint_receipt.json').exists())
            self.assertEqual(len([r for r in rt.calls if r[0]=='evaluate']),2)

    def test_horizon_unstable_checkpoint_has_no_measurement_certificate(self):
        p=plan()
        class Unstable(FakeRuntime):
            def evaluate(self,**kwargs):
                reply=super().evaluate(**kwargs)
                if kwargs['horizon']==8:reply['rows'][4]['asset_price']=1.01
                return reply
        rt=Unstable(p)
        with tempfile.TemporaryDirectory() as directory:
            c=runner.Controller(p,rt,directory);c.prepare();reply=c.evaluate(.2)
            self.assertFalse(reply['certified'])
            checkpoint=json.loads(runner.pinned(reply['payload']['checkpoint_2023']).read_text())
            self.assertTrue(checkpoint['horizon_unstable'])
            self.assertFalse(checkpoint['one_shock_measurement_certified'])
            self.assertFalse((Path(directory)/'latest_completed.json').exists())

    def test_different_psi_cannot_replace_same_psi_fresh_replay_initializer(self):
        import numpy as np
        p=plan();H=6
        class Runtime:
            packet={'parameters':SimpleNamespace(pension=.1)};reference_price=1.;P=SimpleNamespace(psi_child=.2,pension=.1)
            housing='static-elastic';calls=[]
            def identity(self):return p['identity']
            @contextmanager
            def native_budget(self,*args):yield
            def mapping(self,terminal,endpoint,q,b,psi,folder,**kwargs):
                self.calls.append((float(psi[0]),q.copy(),b.copy()))
                record=FakeRuntime(p).evaluate(psi=float(psi[0]),horizon=len(q))
                record.update(accounting_valid=True,policy_calls=1,gates={'native':True},
                    market_residual=[0.]*len(q),fiscal_residual=[0.]*len(q))
                return SimpleNamespace(),record
            def terminal_checks(self,*args,**kwargs):return {'all_checks_pass':True}
        rt=Runtime();adapter=runner.NativeAdapter(rt,p)
        def root(q,J):
            return dict(converged=True,gates={'mapping':True,'housing':True,'social_security':True,'market_replay':True,'fiscal_replay':True},
                final=dict(prices=np.full(H,q),fiscal_values=np.full(H,.1)),final_jacobian=np.eye(2*H)*J,
                final_reproduction_max_abs=0.)
        a=root(1.1,-1.);b=root(1.2,-2.)
        adapter.retain_warm(H,.2,a,'verified_psi_a')
        adapter.retain_warm(H,.21,b,'verified_psi_b')
        a['final']['prices'][:]=9. # Stored coordinates must be independent copies.
        selected,kind=adapter.select_warm(H,.2)
        np.testing.assert_array_equal(selected['prices'],np.full(H,1.1))
        self.assertEqual(kind,'same_psi_verified_root_initialization')
        self.assertEqual(adapter.select_warm(H,.22)[1],'latest_other_psi_initialization')
        adapter._endpoint=lambda *args:(rt.packet,dict(price=1.,population_scale=.922,stationary_pass=True,stationary_renewal_gap=0.))
        def fresh_root(**kwargs):
            np.testing.assert_array_equal(kwargs['initial_prices'],np.full(H,1.1))
            np.testing.assert_array_equal(kwargs['initial_jacobian'],-np.eye(2*H))
            kwargs['evaluate'](kwargs['initial_prices'],kwargs['initial_fiscal_values'])
            kwargs['evaluate'](kwargs['initial_prices'],kwargs['initial_fiscal_values'])
            return root(1.1,-1.)
        with tempfile.TemporaryDirectory() as directory,patch('e5f_four_shock_acceleration.extend_measured_jacobian',return_value=-np.eye(2*H)),patch('e5f_four_shock_acceleration.solve_joint_with_acceleration',side_effect=fresh_root):
            adapter.evaluate(psi=.2,start_year=2007,horizon=H,seed={},gates=p['gates'],budget=p['budget'],
                endpoint_controls=p['endpoint'],path_controls=p['path'],deadline=__import__('time').monotonic()+5,folder=Path(directory))
        self.assertEqual(len(rt.calls),2)
        self.assertEqual([r[0] for r in rt.calls],[.2,.2])
        p['identity']['engine_sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'snapshot identity'):
            adapter.select_warm(H,.2)

    def test_diagnostic_fit_can_report_terminal_nonpass_without_full_certificate(self):
        p=plan()
        class DiagnosticFake(FakeRuntime):
            def evaluate(self,**kwargs):
                reply=super().evaluate(**kwargs)
                reply.update(terminal_pass=False,root_and_terminal_pass=False)
                return reply
            def compare_2023(self,*args):
                return dict(passed=True,available=True,state_value_integrity=True,state_gaps={'normalized_distribution_l1':0.})
        with tempfile.TemporaryDirectory() as directory:
            result=runner.Controller(p,DiagnosticFake(p),directory).run()
            self.assertTrue(result['one_shock_fit_certified'])
            self.assertFalse(result['scientific_validation'])
            self.assertFalse(result['numerical_path_certified'])
            self.assertFalse(result['state_experiment_ready'])
            self.assertTrue(result['state_physical_horizon_stable'])
            self.assertTrue(result['state_value_integrity'])
            self.assertFalse(result['full_path_certified'])
            self.assertTrue(result['exploratory_state_available'])
            self.assertEqual(result['terminal_passes'],[False,False])

    def test_diagnostic_comparison_includes_actual_2023_macro_row(self):
        p=plan();rt=FakeRuntime(p)
        a=rt.evaluate(psi=.2,horizon=6);b=rt.evaluate(psi=.2,horizon=8)
        b['rows'][4]['asset_price']=1.01
        result=runner.diagnostic_horizon_comparison(a,b)
        self.assertFalse(result['passed'])
        self.assertEqual(result['early_macro_years'][-1],2023)

    def test_production_candidate_still_rejects_terminal_nonpass(self):
        p=plan();p.update(mode='production',horizons=[104,128])
        class TerminalFail(FakeRuntime):
            def evaluate(self,**kwargs):
                reply=super().evaluate(**kwargs)
                reply.update(root_and_terminal_pass=False,terminal_pass=False)
                return reply
        with tempfile.TemporaryDirectory() as directory:
            c=runner.Controller(p,TerminalFail(p),directory);c.prepare()
            self.assertFalse(c.evaluate(.2)['certified'])

    def test_actual_2023_state_queues_gate_separately_from_reported_values(self):
        import numpy as np
        def reply(H):
            state=SimpleNamespace(g_pre=np.array([.4,.6]),scheduled_entries=np.array([.1,.2]),scheduled_raw_entries=np.array([.2,.3]))
            native=SimpleNamespace(dated_states={4:dict(state=state)},values=[np.array([1.,2.]) for _ in range(H+1)],
                floor_runtime_paths=dict(prices=np.ones(H),pensions=np.full(H,.1),psi_path=np.full(H,.2)))
            return dict(native_reply=native)
        a,b=reply(6),reply(8)
        b['native_reply'].values[4]=np.array([5.,10.])
        result=runner.compare_2023_states(a,b,lambda q:q)
        self.assertTrue(result['passed'])
        self.assertTrue(result['state_value_integrity'])
        self.assertGreater(result['values']['current_2023_V']['occupied_weighted_relative_gap'],0)
        self.assertFalse(result['values']['current_2023_V']['gating'])
        b['native_reply'].dated_states[4]['state'].scheduled_raw_entries[0]=.3
        self.assertFalse(runner.compare_2023_states(a,b,lambda q:q)['passed'])
        b['native_reply'].values[5]=np.array([float('nan'),2.])
        self.assertFalse(runner.compare_2023_states(a,b,lambda q:q)['state_value_integrity'])
        b['native_reply'].values[5]=np.ones((1,2))
        self.assertFalse(runner.compare_2023_states(a,b,lambda q:q)['state_value_integrity'])

    def test_failed_callback_after_native_increment_reports_actual_call(self):
        p=plan()
        class FailingRuntime:
            packet={};reference_price=1.;P=SimpleNamespace(psi_child=.2,pension=.1)
            total_native_calls=0
            @contextmanager
            def native_budget(self,*args):yield
            def mapping(self,*args,**kwargs):
                self.total_native_calls+=1
                raise RuntimeError('callback failed after counted native call')
        rt=FailingRuntime();adapter=runner.NativeAdapter(rt,p)
        with tempfile.TemporaryDirectory() as directory:
            controller=runner.Controller(p,adapter,directory)
            try:
                adapter._mapping({}, {}, [1.],[.1],[.2],Path(directory),__import__('time').monotonic()+5)
            except RuntimeError as exc:
                receipt=runner.failure_receipt(exc,controller)
                runner.write(Path(directory)/'failure.json',receipt)
            else:self.fail('Expected failed callback')
            saved=json.loads((Path(directory)/'failure.json').read_text())
            self.assertEqual(saved['actual_policy_calls'],1)
            self.assertEqual(saved['control_accounted_policy_calls'],0)
            self.assertTrue(saved['policy_call_count_mismatch'])
            self.assertEqual(saved['policy_call_count_source'],'runtime_total_native_calls')

    def test_exhausted_native_budget_blocks_before_mapping_call(self):
        p=plan()
        class CappedRuntime:
            packet={};reference_price=1.;P=SimpleNamespace(psi_child=.2,pension=.1)
            mapping_called=False
            @contextmanager
            def native_budget(self,deadline,remaining_calls):
                self.remaining_calls=remaining_calls
                if remaining_calls<=0:raise TimeoutError('Native solve budget exhausted before LC')
                yield
            def mapping(self,*args,**kwargs):
                self.mapping_called=True
                raise AssertionError('Must not reach a native call after cap')
        rt=CappedRuntime();adapter=runner.NativeAdapter(rt,p)
        adapter.calls=p['budget']['maximum_policy_calls']
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(TimeoutError,'before LC'):
                adapter._mapping({}, {}, [1.], [.1], [.2],Path(directory),__import__('time').monotonic()+5)
        self.assertEqual(rt.remaining_calls,0)
        self.assertFalse(rt.mapping_called)

    def test_full_numerical_path_does_not_certify_outstanding_policy_contract(self):
        p=plan();p.update(mode='production',horizons=[104,128])
        flags=runner.production_flags(p)
        self.assertTrue(flags['numerical_path_certified'])
        self.assertFalse(flags['production_ready'])
        p['policy_contract_closed']=True
        self.assertTrue(runner.production_flags(p)['production_ready'])
        p['mode']='diagnostic'
        self.assertFalse(runner.production_flags(p)['production_ready'])

    def test_native_adapter_exact_seed_and_joint_loop_with_fake_mapping(self):
        import numpy as np
        p=plan()
        p['endpoint'].update(price_bound_ratios=[.05,20.],slope=1.,max_log_step=.15,damping=.7,max_evaluations=4)
        p['path'].update(price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],max_log_step=.15,damping=.7,max_evaluations=4)
        class NativeFake:
            reference_price=1.;population_scale=.922;housing='static-elastic'
            P=SimpleNamespace(psi_child=.2,pension=.1)
            packet={'parameters':P}
            initial_states=[]
            guard_calls=[]
            @contextmanager
            def native_budget(self,deadline,remaining_calls):
                self.guard_calls.append((deadline,remaining_calls))
                self.guarded=True
                try:yield
                finally:self.guarded=False
            def identity(self):return p['identity']
            def stationary(self,psi,price,folder):
                assert self.guarded
                return self.packet,dict(price=price,psi_child=psi,pension=.1,population_scale=.922,
                    renewal_residual=price-1.,gates=stationary_gates(),accounting_valid=True,policy_calls=1)
            def stationary_state(self,packet,population_scale):return ('endpoint_state',population_scale)
            def terminal_checks(self,*args,**kwargs):return {'all_checks_pass':True}
            def mapping(self,terminal,endpoint,q,b,psi,folder,initial_state,start_year):
                assert self.guarded
                self.initial_states.append(initial_state)
                H=len(q)
                record=FakeRuntime(p).evaluate(psi=float(psi[0]),horizon=H)
                record.update(market_residual=(-np.log(q)).tolist(),fiscal_residual=(-np.log(b/.1)).tolist(),
                    gates={'native':True},accounting_valid=True,policy_calls=H)
                return {'fake_native':True},record
        rt=NativeFake();adapter=runner.NativeAdapter(rt,p)
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            seed=adapter.measure_seed(**p['seed'],gates=p['gates'],budget=p['budget'],deadline=__import__('time').monotonic()+5,folder=folder/'seed')
            self.assertEqual(seed['matrix'].shape,(24,24))
            np.testing.assert_allclose(np.diag(seed['matrix']),-1.,atol=1e-10)
            self.assertEqual(seed['mapping_count'],5)
            result=adapter.evaluate(psi=.2,start_year=2007,horizon=6,seed=seed,gates=p['gates'],budget=p['budget'],
                endpoint_controls=p['endpoint'],path_controls=p['path'],deadline=__import__('time').monotonic()+5,folder=folder/'candidate')
            self.assertTrue(result['root_and_terminal_pass'])
            self.assertEqual(result['replay_maximum_gap'],0.)
            self.assertIn(('endpoint_state',.922),rt.initial_states)
            self.assertNotIn(('endpoint_state',1.),rt.initial_states)
            self.assertEqual(result['policy_calls'],15)
            self.assertEqual(rt.guard_calls[0][1],200)
            self.assertEqual(rt.guard_calls[1][1],188)
            self.assertEqual(rt.guard_calls[-1][1],131)

    def test_short_smoke_terminal_nonpass_does_not_hide_native_root_pass(self):
        reply=dict(accounting_valid=True,stationary_pass=True,root_pass=True,replay_pass=True,
            market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,
            terminal_pass=False,root_and_terminal_pass=False)
        receipt=runner.smoke_readiness(reply)
        self.assertTrue(receipt['native_setup_verified'])
        self.assertFalse(receipt['diagnostic_terminal_pass'])
        self.assertFalse(receipt['production_ready'])

    def test_prepared_native_inputs_changed_identity_rejected_no_fallback(self):
        p=plan()
        prepared=dict(schema='current_floor_native_reference_seed_v1',reference_repeat_verified=True,
            fresh_seed_baseline_verified=True,identity={})
        with self.assertRaisesRegex(ValueError,'Prepared native inputs differ: identity'):
            runner.validate_preparation(p,prepared)

    def test_prepared_native_reference_seed_roundtrip_has_zero_reconstruction_calls(self):
        import numpy as np
        p=plan()
        class PreparedFake(FakeRuntime):
            def prepare_native_reference(self,deadline,folder):
                folder.mkdir(parents=True)
                checkpoint=folder/'selected_native_packet.pkl.gz'
                checkpoint.write_bytes(b'mocked already verified native packet')
                receipt=folder/'reference_reconstruction.json'
                runner.write(receipt,dict(status='passed',identity=self.identity(),checkpoint_sha256=runner.sha(checkpoint),policy_calls=2))
                self.calls.append(('reconstruct',2))
                return dict(accounting_valid=True,policy_calls=2,
                    reference_checkpoint=dict(path=str(checkpoint),sha256=runner.sha(checkpoint)),
                    reconstruction_receipt=dict(path=str(receipt),sha256=runner.sha(receipt)))
            def measure_seed(self,**kwargs):
                folder=kwargs['folder'];original,_=runner.original_modules()
                H=kwargs['horizon']
                def evaluate(q,b):
                    return dict(mapping_valid=True,market_residual=-np.log(q),fiscal_residual=-np.log(b/.1))
                matrix=original.inner.measure_jacobian(evaluate,np.ones(H),np.full(H,.1),kwargs['perturbed_date'],kwargs['log_step'],folder/'measured',dict(identity=self.identity(),native_measured=True))
                runner.write(folder/'baseline_checks.json',dict(valid=True,terminal={'all_checks_pass':True}))
                receipt=json.loads((folder/'measured/receipt.json').read_text())
                self.calls.append(('seed',5))
                return dict(receipt,matrix=matrix,accounting_valid=True,policy_calls=60)
            def restore_native_reference(self,preparation,folder):
                runner.pinned(preparation['reference_checkpoint'])
                self.calls.append(('restore',0))
                return dict(accounting_valid=True,policy_calls=0)
        with tempfile.TemporaryDirectory() as directory:
            cold=PreparedFake(p)
            cold_controller=runner.Controller(p,cold,Path(directory)/'cold');cold_controller.prepare()
            self.assertEqual(cold_controller.policy_calls,62)
            prep=Path(directory)/'cold/native_preparation.json'
            resumed=copy.deepcopy(p);resumed['prepared_native_inputs']=dict(path=str(prep),sha256=runner.sha(prep))
            resumed['mode']='production';resumed['horizons']=[104,128]
            warm=PreparedFake(resumed);warm_controller=runner.Controller(resumed,warm,Path(directory)/'reuse')
            warm_controller.prepare()
            self.assertEqual(warm.calls,[('restore',0)])
            self.assertEqual(warm_controller.policy_calls,0)
            np.testing.assert_array_equal(cold_controller.seed['matrix'],warm_controller.seed['matrix'])

    def test_complete_exact_loop_diagnostic_has_no_production_claim(self):
        p = plan();rt = FakeRuntime(p)
        with tempfile.TemporaryDirectory() as directory:
            result = runner.Controller(p,rt,directory).run()
            self.assertFalse(result['production_ready'])
            self.assertEqual(result['fitted_parameter']['estimate'],.2)
            self.assertGreaterEqual(len([r for r in rt.calls if r[0]=='evaluate']),8)
            self.assertTrue((Path(directory)/'best_so_far.json').is_file())
            self.assertTrue((Path(directory)/'heartbeat.json').is_file())
            self.assertEqual(len([r for r in rt.calls if r[0]=='export']),len([r for r in rt.calls if r[0]=='evaluate'])//2)
            import csv
            with (Path(directory)/'fertility_fit.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r['role'] for r in rows],['validation']*3+['fitted'])
            self.assertEqual([float(r['weight']) for r in rows],[0.,0.,0.,1.])

    def test_accounting_failure_stops_before_fertility_score(self):
        p = plan();rt = FakeRuntime(p);rt.accounting_valid=False
        with tempfile.TemporaryDirectory() as directory:
            c = runner.Controller(p,rt,directory);c.prepare()
            with self.assertRaisesRegex(ValueError,'accounting failed'):
                c.evaluate(.2)
            self.assertFalse((Path(directory)/'latest_completed.json').exists())

    def test_old_seed_identity_rejected(self):
        p = plan();rt = FakeRuntime(p);rt.seed_identity={}
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,'Fresh measured'):
                runner.Controller(p,rt,directory).prepare()

    def test_short_horizons_cannot_claim_production(self):
        p = plan();p['mode']='production'
        with self.assertRaisesRegex(ValueError,'104/128'):
            runner.preflight(p)

    def test_unapproved_extra_fiscal_relaxation_rejected(self):
        p = plan();p['gates']['fiscal_tolerance']=1e-4
        with self.assertRaisesRegex(ValueError,'2e-5'):
            runner.preflight(p)

    def test_unbounded_time_rejected(self):
        p = plan();p['budget']['total_seconds']=float('inf')
        with self.assertRaisesRegex(ValueError,'finite budget'):
            runner.preflight(p)

    def test_target_drift_rejected_before_native(self):
        p = plan();p['target_contract']['rows'][3]['target']+=.1;rt=FakeRuntime(p)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,'fingerprint changed'):
                runner.Controller(p,rt,directory).prepare()
            self.assertEqual(rt.calls,[])

    def test_rescaled2023_export_rejected(self):
        p=plan();rt=FakeRuntime(p);rt.bad_export=True
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,'Exact index4'):
                runner.Controller(p,rt,directory).run()


if __name__ == '__main__':
    unittest.main()
