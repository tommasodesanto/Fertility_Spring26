"""Deterministic exact-loop controller tests; zero native model solves."""
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from contextlib import contextmanager

import one_shock_floor as runner


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
            root_and_terminal_pass=True,stationary_pass=True,stationary_renewal_gap=0.,
            market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,rows=rows,fertility=fertility,
            reference_manifest_sha256=self.identity()['reference_sha256'],source_pins=self.identity()['source_pins'],
            housing='static-elastic',shock_contract=dict(psi=psi,start_year=2007))

    def export_2023(self,reply,folder):
        return dict(calendar_year=2023,period_index=4,exact_native_state=True,
            reconstructed_or_rescaled=self.bad_export,queue_lags=[16,20],forecast_and_continuation_saved=True)

    def render_standard(self,reply,folder,deadline):
        folder.mkdir(parents=True)
        paths = [folder/name for name in self.p['standard_plot_names']]
        for path in paths:
            path.write_bytes(b'fake controller smoke plot, not native visual evidence')
        return paths


class Tests(unittest.TestCase):
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
                    renewal_residual=price-1.,gates={'native':True},accounting_valid=True,policy_calls=1)
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
