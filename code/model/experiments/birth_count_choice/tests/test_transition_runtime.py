"""No-solve checks for the current saved-case transition adapter."""
import copy
import csv
import gzip
import importlib
import json
import pickle
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
r = importlib.import_module('experiments.birth_count_choice.transition_runtime')


def handoff(folder):
    names = ('metadata.json', 'native_result.npz', 'input_contract.json', 'target_fit.csv',
             'parameters.csv', 'native/phase_b_ge/selected_root/closure.json')
    sources = list((r.ROOT/'code/model/experiments/birth_count_choice/model').rglob('*.py'))
    sources += [Path(r.__file__), Path(r.retained.__file__)]
    record = dict(schema=r.SCHEMA, saved_case=str(r.CASE),
        saved_files={name:r.sha(r.CASE/name) for name in names},
        source_pins={str(path.relative_to(r.ROOT)):r.sha(path) for path in sources})
    path = Path(folder)/'handoff.json'; path.write_text(json.dumps(record))
    return dict(path=str(path), sha256=r.sha(path))


def report_fixture(folder):
    folder = Path(folder); (folder/'standard_diagnostics').mkdir(parents=True)
    for i in range(17): (folder/'standard_diagnostics'/f'{i}.png').write_bytes(b'original-render')
    tables = {'target_fit.csv':[dict(moment=str(i),target='1',model='1',gap='0',weight='2',loss_contribution='0',role='scored') for i in range(14)],
              'parameters.csv':[dict(parameter=str(i),estimate='.18',lower='0',upper='1',status='fixed') for i in range(31)]}
    for name,rows in tables.items():
        with (folder/name).open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    r.write(folder/'closure.json',dict(population_scale=1.,numeric=[.2,3],status='passed'))
    return folder


class CurrentTransitionTests(unittest.TestCase):
    def test_constructor_restores_explicit_thread_request_after_frozen_import(self):
        import numba
        original=r.reporting.build_context
        def one_thread_import(*args, **kwargs):
            context=original(*args, **kwargs)
            r.os.environ['NUMBA_NUM_THREADS']='1'
            return context
        with tempfile.TemporaryDirectory() as temp:
            pin=handoff(temp)
            with patch.dict(r.os.environ, {'NUMBA_NUM_THREADS':'8'}), \
                 patch.object(numba.config,'NUMBA_NUM_THREADS',8), \
                 patch.object(numba,'set_num_threads') as mask, \
                 patch.object(numba,'get_num_threads',return_value=8), \
                 patch.object(r.reporting,'build_context',one_thread_import):
                r.CurrentEstateARuntime.from_handoff(pin,Path(temp)/'runtime')
                mask.assert_called_once_with(8)
                self.assertEqual(r.os.environ['NUMBA_NUM_THREADS'],'8')
                receipt=r.read(Path(temp)/'runtime/constructor.json')
                self.assertEqual(receipt['numba_threads']['requested'],8)
                self.assertEqual(receipt['numba_threads']['actual'],8)
            with patch.dict(r.os.environ, {'NUMBA_NUM_THREADS':str(numba.config.NUMBA_NUM_THREADS+1)}):
                with self.assertRaisesRegex(RuntimeError,'initialized pool limit'):
                    r.CurrentEstateARuntime.from_handoff(pin,Path(temp)/'invalid')

    def test_disclosed_platform_criteria_value_controls_support_and_occupied_states(self):
        saved=SimpleNamespace(V=np.array([100.,-1e10,1.]),birth_count_pre_distribution=np.array([0.,0.,.2]),
            c_pol=np.array([100.,3.,1.]),bp_pol=np.array([100.,3.,1.]),
            c_pol_stay=np.array([100.,3.,1.]),bp_pol_stay=np.array([100.,3.,1.]),
            tenure_probs=np.array([.2,0.,.8],dtype=np.float32),other=np.array([1.,2.,3.]))
        fresh=copy.deepcopy(saved);fresh.V[1]+=1e-6
        for field in r.CONTROL_FIELDS:
            getattr(fresh,field)[0]+=1e-7
            getattr(fresh,field)[1]+=1.
            getattr(fresh,field)[2]+=5e-11
        fresh.tenure_probs[0]=np.nextafter(fresh.tenure_probs[0],np.float32(1.))
        result=r.verify_solution_arrays(saved,fresh)
        self.assertTrue(result['all_applicable_criteria_passed'])
        self.assertFalse(result['criteria']['optimizer_identity_certified'])
        self.assertEqual(result['arrays']['c_pol']['details']['dead_audit']['absolute_maximum_gap'],1.)
        mutations=(lambda x:x.V.__setitem__(0,100.+2e-10),lambda x:x.V.__setitem__(1,0.),
            lambda x:x.c_pol.__setitem__(0,100.+1e-4),lambda x:x.c_pol.__setitem__(2,1.+2e-10),
            lambda x:x.birth_count_pre_distribution.__setitem__(0,.2),
            lambda x:x.tenure_probs.__setitem__(0,saved.tenure_probs[0]+np.float32(5*np.spacing(saved.tenure_probs[0]))),
            lambda x:x.tenure_probs.__setitem__(1,np.nextafter(np.float32(0.),np.float32(1.))),
            lambda x:x.other.__setitem__(0,1.+2e-10))
        for mutate in mutations:
            wrong=copy.deepcopy(fresh);mutate(wrong)
            with self.assertRaises(RuntimeError):r.verify_solution_arrays(saved,wrong)

    def test_solution_bridge_checks_complete_inventory_masks_and_absolute_gaps(self):
        saved = SimpleNamespace(V=np.array([1.,np.nan,np.inf,-np.inf]), discrete=np.array([1],dtype=np.int16))
        fresh = copy.deepcopy(saved);fresh.V[0]+=5e-11
        self.assertLessEqual(r.verify_solution_arrays(saved,fresh)['maximum_gap'],1e-10)
        for mutate in (lambda x:delattr(x,'discrete'),lambda x:setattr(x,'discrete',x.discrete.astype(np.int32)),
                       lambda x:x.V.__setitem__(0,1.+2e-10),lambda x:x.V.__setitem__(1,0.),
                       lambda x:x.discrete.__setitem__(0,2)):
            wrong = copy.deepcopy(saved); mutate(wrong)
            with self.assertRaises(RuntimeError):r.verify_solution_arrays(saved,wrong)

    def test_platform_report_bridge_allows_only_derived_roundoff_and_rendering(self):
        with tempfile.TemporaryDirectory() as temp:
            left=report_fixture(Path(temp)/'saved');right=report_fixture(Path(temp)/'fresh')
            original=(right/'target_fit.csv').read_text()
            (right/'target_fit.csv').write_text(original.replace('0,1,1,0,2,0,scored','0,1,1.00000000005,0.00000000005,2,0,scored',1))
            (right/'standard_diagnostics/0.png').write_bytes(b'fresh-render')
            result=r.compare_saved_platform_reports(left,right)
            self.assertEqual(len(result['render_only_hash_differences']),1)
            with self.assertRaises(RuntimeError):r.compare_reports(left,right)
            for replacement in ('0,1.00000000001,1,0,2,0,scored','0,1,1,0,2.00000000001,0,scored',
                                '0,1,1.000000001,0,2,0,scored','0,1,1,0,2,0,validation'):
                (right/'target_fit.csv').write_text(original.replace('0,1,1,0,2,0,scored',replacement,1))
                with self.assertRaises(RuntimeError):r.compare_saved_platform_reports(left,right)
            (right/'target_fit.csv').write_text(original)
            r.write(right/'closure.json',dict(population_scale=1.+2e-10,numeric=[.2,3],status='passed'))
            with self.assertRaises(RuntimeError):r.compare_saved_platform_reports(left,right)

    def test_saved_bridge_rejects_missing_stale_or_incomplete_array_proof(self):
        with tempfile.TemporaryDirectory() as temp:
            runtime=r.CurrentEstateARuntime();runtime.report=report_fixture(Path(temp)/'saved')
            current=report_fixture(Path(temp)/'repeat_0/phase_b_ge/selected_root')
            runtime.case=Path(temp)/'case';runtime.case.mkdir();(runtime.case/'native_result.npz').write_bytes(b'archive')
            runtime.saved=SimpleNamespace(solution=SimpleNamespace(**{f'x{i}':np.array([float(i)]) for i in range(78)}))
            runtime.identity=lambda:dict(source='verified')
            with self.assertRaisesRegex(RuntimeError,'requires verified'):runtime.compare_repeated(runtime.report,current)
            checkpoint=Path(temp)/'native.pkl.gz'
            with gzip.open(checkpoint,'wb') as stream:
                pickle.dump(dict(identity=runtime.identity(),solution=copy.deepcopy(runtime.saved.solution)),stream)
            arrays=r.verify_solution_arrays(runtime.saved.solution,copy.deepcopy(runtime.saved.solution))
            proof=dict(schema='current_estate_a_saved_array_bridge_v1',status='passed',identity=runtime.identity(),
                saved_archive_sha256=r.sha(runtime.case/'native_result.npz'),current_report=str(current.resolve()),
                native_checkpoint=dict(path=str(checkpoint),sha256=r.sha(checkpoint)),arrays=arrays)
            proof_path=current.parent.parent/'saved_array_verification.json';r.write(proof_path,proof)
            runtime.compare_repeated(runtime.report,current)
            for mutate in (lambda x:x.update(identity=dict(source='old')),lambda x:x['arrays']['arrays'].pop('x0'),
                           lambda x:x['arrays'].update(maximum_gap=2e-10),
                           lambda x:x['arrays']['criteria'].update(value_eps_multiplier=32.),
                           lambda x:x['arrays']['arrays']['x0'].update(criterion='looser'),
                           lambda x:x['arrays']['arrays']['x0'].update(status='failed'),
                           lambda x:x['native_checkpoint'].update(sha256='0'*64)):
                wrong=copy.deepcopy(proof);mutate(wrong);r.write(proof_path,wrong)
                with self.assertRaises(RuntimeError):runtime.compare_repeated(runtime.report,current)

    def test_complete_native_api_requires_actual_objects(self):
        runtime = r.CurrentEstateARuntime()
        runtime.model = r.reporting.production_model_facade()
        with runtime.native_api_bindings() as callbacks:
            self.assertEqual(len(callbacks), 34)
            for name in ('_gate_dead_mass_at_age', 'estate_housing_value', 'realize_stayer_cross_section'):
                if name != 'estate_housing_value':
                    self.assertIn(name, callbacks)
        original = runtime.model._gate_dead_mass_at_age
        def substitute(*args, **kwargs):
            return original(*args, **kwargs)
        substitute.__module__ = original.__module__
        runtime.model._gate_dead_mass_at_age = substitute
        with self.assertRaisesRegex(RuntimeError, 'Native compatibility object differs'):
            with runtime.native_api_bindings():
                pass

    def test_owned_count_snapshots_survive_parameter_changes_and_cache_pickle(self):
        import pickle
        P = SimpleNamespace(birth_count_action_probs=np.array([[.2,.8]]),
                            birth_count_realized_probs=np.array([[.6,.4]]))
        p = r.attach_owned_count_policy(SimpleNamespace(V=np.zeros(1)), P)
        P.birth_count_action_probs[:] = 0
        P.birth_count_realized_probs[:] = 0
        q = pickle.loads(pickle.dumps(p))
        np.testing.assert_array_equal(q.birth_count_action_probs, [[.2,.8]])
        np.testing.assert_array_equal(q.birth_count_realized_probs, [[.6,.4]])
        with self.assertRaisesRegex(RuntimeError, 'unavailable'):
            r.attach_owned_count_policy(SimpleNamespace(V=np.zeros(2)), P)

    def test_handoff_source_or_case_drift_fails_before_context(self):
        with tempfile.TemporaryDirectory() as folder:
            pin = handoff(folder)
            record = json.loads(Path(pin['path']).read_text())
            record['saved_files']['input_contract.json'] = '0'*64
            Path(pin['path']).write_text(json.dumps(record)); pin['sha256'] = r.sha(pin['path'])
            with self.assertRaisesRegex(RuntimeError, 'Saved case file differs'):
                r.authenticate_handoff(pin)

    def test_actual_context_bootstrap_has_no_solve_and_preserves_native_state(self):
        def forbidden(*args, **kwargs):
            raise AssertionError('Bootstrap attempted a numerical solve')
        forbidden.__module__ = r.household.__name__
        with tempfile.TemporaryDirectory() as folder:
            pin = handoff(folder)
            with patch.object(r.household, 'solve_bellman_full_markov_income', forbidden), \
                 patch.object(r.solver, 'solve_bellman_full_markov_income', forbidden):
                runtime = r.CurrentEstateARuntime.from_handoff(pin, Path(folder)/'runtime')
                receipt = runtime.bootstrap_saved_reference(Path(folder)/'bootstrap')
            # Restore the intercepted callee on the shared facade, retaining
            # the object owned by reporting and its guarded recent observer.
            runtime.model.solve_bellman_full_markov_income = r.household.solve_bellman_full_markov_income
            recent = next(cell.cell_contents for cell in runtime.rt['observe_recent_parent_flow'].__closure__
                          if hasattr(cell.cell_contents, 'observe_recent_parent_flow'))
            self.assertIs(runtime.model, runtime.rt['model'])
            self.assertIs(recent.observe_recent_parent_flow.__globals__['_production_model_facade'], runtime.model)
            self.assertEqual(runtime.total_native_calls, 0)
            status = runtime.native_status()
            self.assertEqual(status['housing'], 'static-elastic')
            self.assertEqual(status['grid'], [120, 9])
            self.assertEqual(status['population_scale'], 1.0009489339264241)
            self.assertTrue(status['reconstruction_pending'])
            self.assertFalse(runtime.reference_verified)
            self.assertEqual(receipt['policy_calls'], 0)
            self.assertAlmostEqual(receipt['actual_initial_population'], 1.0009489339264241, places=12)
            self.assertFalse(np.array_equal(receipt['adjusted_queue'], receipt['raw_queue']))
            self.assertEqual(runtime.P.birth_count_choice_cap, 1)
            self.assertFalse(runtime.P.joint_nested_choice)
            self.assertTrue(runtime.P.bequest_net_of_selling_cost)
            self.assertTrue(runtime.P.estate_flow_net_of_selling_cost)
            before = runtime.pf.calendar.model
            with runtime.native_bindings():
                self.assertIs(runtime.pf.calendar.model, runtime.model)
                self.assertIs(recent.observe_recent_parent_flow.__globals__['_production_model_facade'],
                              runtime.pf.calendar.model)
                self.assertIs(runtime.pf.transition.apply_sequential_fertility, r.apply_count_fertility)
                with self.assertRaisesRegex(RuntimeError, 'reconfiguration forbidden'):
                    runtime.pf.transition.configure_sequential_model()
            self.assertIs(runtime.pf.calendar.model, before)
            P = copy.deepcopy(runtime.P); P.joint_nested_choice = True
            runtime.P = P
            with self.assertRaisesRegex(RuntimeError, 'Joint-choice'):
                r.validate_contract(P, runtime.saved.closure, r.read(runtime.case/'input_contract.json'))


if __name__ == '__main__':
    unittest.main()
