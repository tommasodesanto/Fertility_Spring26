"""Historical-root provenance tests; no model or cluster evaluations."""
from __future__ import annotations

import copy
from contextlib import ExitStack
import gzip
import json
from pathlib import Path
import pickle
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

import collect_e5f_matched_pf_price_jacobian as collector
import run_e5f_matched_pf_historical_root as driver


class RootRuntimeContractTest(unittest.TestCase):
    def test_long_root_budget_reaches_input_verification_but_excess_budget_fails(self):
        class InputVerificationReached(Exception):
            pass
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'contract.json'
            for seconds in (18000, driver.MAXIMUM_ROOT_SECONDS, driver.MAXIMUM_ROOT_SECONDS + 1):
                path.write_text(json.dumps(dict(schema=driver.joined.HISTORY_SMOKE_SCHEMA,
                    seconds=seconds, checkpoint='/frozen/checkpoint', checkpoint_sha256='a' * 64)))
                with mock.patch.object(driver.primitive, 'verify',
                        side_effect=[None, InputVerificationReached()]):
                    expected = InputVerificationReached if seconds <= driver.MAXIMUM_ROOT_SECONDS else ValueError
                    with self.subTest(seconds=seconds), self.assertRaises(expected):
                        driver.joined.load_smoke_contract(path, 'b' * 64, 'sequential',
                            maximum_seconds=driver.MAXIMUM_ROOT_SECONDS)


class JacobianPacketValidationTest(unittest.TestCase):
    def setUp(self):
        self.evaluator = 'code/model/tools/run_e5f_matched_pf_baseline.py'
        self.model_source = 'code/model/intergen_eqscale_seq_optimized/solver.py'
        self.bridge_source = 'code/model/tools/run_e5f_matched_pf_history.py'
        shared = dict(
            checkpoint_sha256='1' * 64, selected_summary_sha256='2' * 64,
            normalized_checkpoint_sha256='3' * 64,
            normalized_summary_sha256='4' * 64,
            normalized_contract_sha256='5' * 64,
            terminal_checkpoint_sha256='6' * 64,
            terminal_summary_sha256='7' * 64,
            terminal_contract_sha256='8' * 64,
            demographic_sources={name: {'path': '/data/' + filename, 'sha256': '9' * 64}
                for name, filename in dict(population_mid='population_mid.csv',
                    births_mid='births_mid.csv', survival='survival.csv',
                    vintage_2025='vintage_2025_age_sex.csv',
                    acs_headship='acs_headship_profiles.csv').items()},
            target_fingerprint='a' * 64, path_date_count=6,
            initial_price_rule='log_old_to_selected_2023_then_terminal',
            terminal_preference_rule='hold_normalized_2023_intercept',
            probe_log_step=.01,
            source_sha256={self.evaluator: 'b' * 64, self.model_source: 'c' * 64,
                           self.bridge_source: 'd' * 64},
        )
        self.packet = dict(
            schema=collector.SCHEMA,
            status='complete_validated_finite_difference_jacobian',
            arm='sequential', target_fingerprint='a' * 64,
            prices=[1.] * 6, residual=[.1] * 6,
            years=list(range(2007, 2031, 4)), jacobian=(-np.eye(6)).tolist(),
            shared_contract=shared,
            provenance=dict(anchor={'directory': '/probes/anchor'},
                columns=[{'directory': f'/probes/column_{j}', 'coordinate': j} for j in range(6)]),
        )
        self.contract = copy.deepcopy(shared)
        self.contract['arm'] = 'sequential'
        self.contract['source_sha256'][self.evaluator] = 'e' * 64
        self.contract['source_sha256']['code/model/tools/run_e5f_matched_pf_historical_root.py'] = 'f' * 64
        self.contract['reviewed_evaluator_driver_change'] = dict(
            path=self.evaluator, from_sha256='b' * 64, to_sha256='e' * 64,
            scope='explicit supplied-price hook and optional dated-state checkpoint; unchanged economic evaluator')

    def validate(self, *, packet=None, contract=None, arm='sequential', repeated=None):
        p = self.packet if packet is None else packet
        c = self.contract if contract is None else contract
        with mock.patch.object(collector, 'collect', return_value=p if repeated is None else repeated) as recollect:
            result = driver.validate_jacobian_packet(p, c, arm)
        return result, recollect

    def test_full_source_pin_chain_and_explicit_hook_change_pass(self):
        p_before, c_before = copy.deepcopy(self.packet), copy.deepcopy(self.contract)
        result, recollect = self.validate()
        self.assertIsNone(result)
        recollect.assert_called_once_with('/probes/anchor', [f'/probes/column_{j}' for j in range(6)])
        self.assertEqual(self.packet, p_before)
        self.assertEqual(self.contract, c_before)

    def test_nested_arm_can_use_its_own_matching_packet(self):
        self.packet['arm'] = self.contract['arm'] = 'nested'
        self.validate(arm='nested')

    def test_changed_model_or_bridge_source_is_rejected(self):
        for source in (self.model_source, self.bridge_source):
            c = copy.deepcopy(self.contract)
            c['source_sha256'][source] = '0' * 64
            with self.subTest(source=source), self.assertRaisesRegex(ValueError, 'economic source changed'):
                self.validate(contract=c)

    def test_missing_existing_source_pin_is_rejected(self):
        del self.contract['source_sha256'][self.model_source]
        with self.assertRaisesRegex(ValueError, 'economic source changed'):
            self.validate()

    def test_reviewed_hook_requires_exact_from_to_and_scope(self):
        for name in ('path', 'from_sha256', 'to_sha256', 'scope'):
            c = copy.deepcopy(self.contract)
            c['reviewed_evaluator_driver_change'][name] = 'not the reviewed change'
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'reviewed'):
                self.validate(contract=c)

    def test_old_hook_scope_cannot_authorize_new_checkpoint_source(self):
        self.contract['reviewed_evaluator_driver_change']['scope'] = (
            'explicit supplied-price hook only; unchanged economic evaluator')
        with self.assertRaisesRegex(ValueError, 'reviewed'):
            self.validate()

    def runtime_change(self):
        name = 'code/model/tools/run_e5f_matched_pf_historical_root.py'
        self.packet['shared_contract']['source_sha256'][name] = '1' * 64
        self.contract['reviewed_root_runtime_changes'] = dict(
            scope='root runtime ceiling and explicit provenance checks only; unchanged economic evaluator and numerical root',
            files={name: dict(from_sha256='1' * 64, to_sha256='f' * 64)})
        return name

    def test_explicit_runtime_change_preserves_scientific_source_guards(self):
        self.runtime_change()
        self.validate()
        self.contract['source_sha256'][self.model_source] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'economic source changed'):
            self.validate()

    def test_runtime_change_requires_exact_pins_scope_and_existing_file(self):
        name = self.runtime_change()
        for mutation in ('missing', 'from', 'to', 'scope', 'extra', 'deleted_source'):
            c = copy.deepcopy(self.contract)
            review = c['reviewed_root_runtime_changes']
            if mutation == 'missing':
                del c['reviewed_root_runtime_changes']
            elif mutation in ('from', 'to'):
                review['files'][name][mutation + '_sha256'] = '0' * 64
            elif mutation == 'scope':
                review['scope'] = 'arbitrary numerical change'
            elif mutation == 'extra':
                review['files'][self.model_source] = dict(from_sha256='c'*64, to_sha256='0'*64)
            else:
                del c['source_sha256'][name]
            with self.subTest(mutation=mutation), self.assertRaisesRegex(ValueError, 'reviewed root runtime'):
                self.validate(contract=c)

    def test_changed_checkpoint_summary_or_demography_is_rejected(self):
        for name in ('checkpoint_sha256', 'selected_summary_sha256',
                     'normalized_checkpoint_sha256', 'normalized_summary_sha256',
                     'normalized_contract_sha256', 'terminal_checkpoint_sha256',
                     'terminal_summary_sha256', 'terminal_contract_sha256',
                     'demographic_sources'):
            c = copy.deepcopy(self.contract)
            c[name] = {} if name == 'demographic_sources' else '0' * 64
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'scientific input contracts differ'):
                self.validate(contract=c)

    def test_wrong_packet_arm_and_target_are_rejected(self):
        for name, value in (('arm', 'nested'), ('target_fingerprint', '0' * 64)):
            p = copy.deepcopy(self.packet)
            p[name] = value
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.validate(packet=p)

    def test_current_root_contract_requires_explicit_matching_arm(self):
        for arm in (None, 'nested'):
            c = copy.deepcopy(self.contract)
            if arm is None:
                del c['arm']
            else:
                c['arm'] = arm
            with self.subTest(arm=arm), self.assertRaises(ValueError):
                self.validate(contract=c)

    def test_recollection_failure_is_propagated_before_root(self):
        with mock.patch.object(collector, 'collect', side_effect=ValueError('failed numerical reproduction')) as recollect:
            with self.assertRaisesRegex(ValueError, 'failed numerical reproduction'):
                driver.validate_jacobian_packet(self.packet, self.contract, 'sequential')
        self.assertEqual(recollect.call_count, 1)

    def test_changed_recollected_jacobian_is_rejected(self):
        repeated = copy.deepcopy(self.packet)
        repeated['jacobian'][0][0] += .001
        with self.assertRaisesRegex(ValueError, 'does not reproduce'):
            self.validate(repeated=repeated)

    def test_changed_recollected_provenance_is_rejected(self):
        repeated = copy.deepcopy(self.packet)
        repeated['provenance']['columns'][0]['summary_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'does not reproduce'):
            self.validate(repeated=repeated)

    def test_dates_shape_and_nonfinite_matrix_rejected(self):
        variants = []
        p = copy.deepcopy(self.packet); p['years'][-1] += 4; variants.append(p)
        p = copy.deepcopy(self.packet); p['jacobian'] = [[1.]]; variants.append(p)
        p = copy.deepcopy(self.packet); p['jacobian'][0][0] = float('nan'); variants.append(p)
        for p in variants:
            with self.subTest(years=p['years'], matrix_shape=np.asarray(p['jacobian']).shape), self.assertRaisesRegex(ValueError, 'horizon or finite matrix'):
                self.validate(packet=p)

    def test_changed_horizon_rule_and_probe_step_rejected(self):
        for name, value in (('path_date_count', 7), ('probe_log_step', .02),
                            ('initial_price_rule', 'different'),
                            ('terminal_preference_rule', 'different')):
            c = copy.deepcopy(self.contract); c[name] = value
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'scientific input contracts differ'):
                self.validate(contract=c)


class RestartReceiptTest(unittest.TestCase):
    def setUp(self):
        fixture = JacobianPacketValidationTest()
        fixture.setUp()
        self.packet, self.contract = fixture.packet, fixture.contract
        self.contract['jacobian_packet_sha256'] = '7' * 64
        self.previous = copy.deepcopy(self.contract)
        self.previous['source_sha256'][fixture.evaluator] = '0' * 64
        best = dict(prices=[1.] * 6, residual=[.01] * 6, score=.01,
                    mapping_valid=True, payload={'directory': '/old/evaluation_005'})
        self.history = dict(best=best, final=copy.deepcopy(best), evaluations=6,
            final_jacobian=(-np.eye(6)).tolist(), converged=False,
            status='evaluation_budget_exhausted', final_reproduction_max_abs=0.)
        self.summary = dict(arm='sequential', best=copy.deepcopy(best), evaluations=6,
            finite_horizon_market_converged=False, final_reproduction_max_abs=0.)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)

    def save_receipts(self):
        for name, filename, data in (
            ('restart_contract', 'contract.json', self.previous),
            ('restart_history', 'root_history.json', self.history),
            ('restart_summary', 'summary.json', self.summary)):
            path = Path(self.directory.name)/filename
            path.write_text(json.dumps(data))
            self.contract[name] = str(path)
            self.contract[name+'_sha256'] = driver.primitive.digest(path)

    def load(self):
        return driver.load_restart(self.contract, self.packet, 'sequential')

    def test_no_restart_uses_original_panel(self):
        prices, matrix = self.load()
        self.assertEqual(prices, self.packet['prices'])
        self.assertEqual(matrix, self.packet['jacobian'])

    def test_valid_reproduced_unfinished_root_supplies_prices_and_preconditioner(self):
        self.history['best']['prices'] = self.history['final']['prices'] = [1.2] * 6
        self.summary['best'] = copy.deepcopy(self.history['best'])
        self.history['final_jacobian'] = (-2*np.eye(6)).tolist()
        self.save_receipts()
        prices, matrix = self.load()
        np.testing.assert_array_equal(prices, [1.2] * 6)
        np.testing.assert_array_equal(matrix, -2*np.eye(6))

    def test_partial_and_hash_only_receipts_are_rejected(self):
        for key in ('restart_contract', 'restart_contract_sha256'):
            with self.subTest(key=key):
                self.contract[key] = 'incomplete'
                with self.assertRaises(ValueError):
                    self.load()
                del self.contract[key]

    def test_changed_scientific_inputs_or_kernel_are_rejected(self):
        for key in ('checkpoint_sha256', 'target_fingerprint', 'arm',
                    'jacobian_packet_sha256', 'path_date_count', 'demographic_sources'):
            original = self.previous[key]
            self.previous[key] = 'changed'
            self.save_receipts()
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'scientific inputs'):
                self.load()
            self.previous[key] = original
        source = 'code/model/intergen_eqscale_seq_optimized/solver.py'
        self.previous['source_sha256'][source] = 'changed'
        self.save_receipts()
        with self.assertRaisesRegex(ValueError, 'economic source'):
            self.load()

    def test_changed_receipt_content_fails_hash_verification(self):
        self.save_receipts()
        Path(self.contract['restart_history']).write_text('{}')
        with self.assertRaises((ValueError, RuntimeError)):
            self.load()

    def test_receipts_require_one_directory_and_canonical_filenames(self):
        self.save_receipts()
        for value in ('relative/summary.json', str(Path(self.directory.name)/'renamed.json'),
                      str(Path(self.directory.name)/'other'/'summary.json')):
            self.contract['restart_summary'] = value
            with self.subTest(path=value), self.assertRaises(ValueError):
                self.load()

    def test_failed_unreproduced_or_already_converged_roots_rejected(self):
        history, summary = copy.deepcopy(self.history), copy.deepcopy(self.summary)
        for variant in ('missing_final', 'invalid_best', 'invalid_final', 'already_converged',
                        'failed_replay', 'nonfinite_replay', 'different_summary_best', 'different_counts'):
            self.history, self.summary = copy.deepcopy(history), copy.deepcopy(summary)
            if variant == 'missing_final': self.history['final'] = None
            elif variant == 'invalid_best': self.history['best']['mapping_valid'] = False
            elif variant == 'invalid_final': self.history['final']['mapping_valid'] = False
            elif variant == 'already_converged': self.summary['finite_horizon_market_converged'] = True
            elif variant == 'failed_replay': self.summary['final_reproduction_max_abs'] = .01
            elif variant == 'nonfinite_replay': self.summary['final_reproduction_max_abs'] = float('nan')
            elif variant == 'different_summary_best': self.summary['best']['score'] = .02
            elif variant == 'different_counts': self.summary['evaluations'] = 5
            self.save_receipts()
            with self.subTest(variant=variant), self.assertRaises(ValueError):
                self.load()

    def test_matrix_price_residual_shapes_and_replay_values_rejected(self):
        history = copy.deepcopy(self.history)
        for variant in ('matrix_shape', 'matrix_nan', 'negative_price', 'changed_final_price',
                        'final_residual_scalar', 'changed_final_residual', 'wrong_score'):
            self.history = copy.deepcopy(history)
            if variant == 'matrix_shape': self.history['final_jacobian'] = [[1.]]
            elif variant == 'matrix_nan': self.history['final_jacobian'][0][0] = float('nan')
            elif variant == 'negative_price': self.history['best']['prices'][0] = -1.
            elif variant == 'changed_final_price': self.history['final']['prices'][0] = 1.1
            elif variant == 'final_residual_scalar': self.history['final']['residual'] = .01
            elif variant == 'changed_final_residual': self.history['final']['residual'][0] += 1e-7
            elif variant == 'wrong_score': self.history['best']['score'] = .02
            self.summary['best'] = copy.deepcopy(self.history['best'])
            self.save_receipts()
            with self.subTest(variant=variant), self.assertRaises(ValueError):
                self.load()


class Initial2023CheckpointTest(unittest.TestCase):
    """Run the actual probe orchestration with model evaluations replaced by mocks."""

    def run_probe(self, output, *, save_state=True, gate_failure=False):
        baseline = driver.baseline
        import e5f_matched_pf_moments as moments
        import run_e5f_perfect_foresight_person_demography_policy as terminal_checks
        old_parameters = SimpleNamespace(psi_child=.3, user_cost_rate=.1,
            tau_H=.04, property_tax_lump_sum_transfer=0., joint_nested_choice=False)
        people = SimpleNamespace(year=2023, persons=np.array([2.]), heads=np.array([1.]))
        demographics = SimpleNamespace(initial_person_state=people)
        old = SimpleNamespace(parameters=old_parameters, psi_path=np.linspace(.3, -.02, 5),
            supply_rule=SimpleNamespace(initial_price=1., initial_housing_stock=2., eta=.63),
            stationary_g_pre=np.array([1.]), b_grid=np.array([0., 1.]),
            shared=None, policy=SimpleNamespace(price=np.array([1.])),
            diagnostics={'normalization': {}}, initial_state=object(), historical_conditioning=object())
        terminal = dict(parameters=old_parameters, policy=SimpleNamespace(price=[1.], V=np.array([3.])),
            fixed_point=SimpleNamespace(g_pre=np.array([5.]), persons=people))
        history_pre = np.array([11., 12.])
        final = SimpleNamespace(g_pre=np.array([91., 92.]),
            persons=SimpleNamespace(year=2031, persons=np.array([4.]), heads=np.array([3.])))
        result = SimpleNamespace(history=SimpleNamespace(terminal_state=SimpleNamespace(g_pre=history_pre)),
            person_tail=SimpleNamespace(terminal_state=final), bellman_solves=12,
            rows=[dict(calendar_year=y, housing_demand=2., housing_supply=2.) for y in range(2007, 2031, 4)])
        report = dict(target_fit_rows=[dict(moment='test', target=1., model=1.)],
            parameter_rows=[dict(parameter='test', value=1.)], loss=0.)
        seed = SimpleNamespace(price=1., selected={'panel_design': {'domain': [
            dict(name='test', lower=0., upper=1., transform='linear')]},
            'best_candidate': {'theta': [1.]}})
        contract = dict(path_date_count=6, probe_coordinate=-1, probe_log_step=.01,
            terminal_preference_rule='hold_normalized_2023_intercept',
            initial_price_rule='log_old_to_selected_2023_then_terminal', target_fingerprint='target',
            save_initial_2023_state=save_state)
        args = SimpleNamespace(arm='sequential', contract_sha256='contract-hash')
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(baseline, 'load_normalized', return_value=(old, demographics)))
            stack.enter_context(mock.patch.object(baseline, 'load_terminal', return_value=terminal))
            stack.enter_context(mock.patch.object(driver.rent_domain, 'project_price_path_to_positive_rents',
                side_effect=lambda prices, **kwargs: (np.asarray(prices), {'adjusted_period_count': 0})))
            stack.enter_context(mock.patch.object(baseline.pf, 'rents_from_asset_prices', return_value=np.ones(6)))
            stack.enter_context(mock.patch.object(moments.measurement, 'TRANSITION_SEARCH_DOMAIN', ()))
            stack.enter_context(mock.patch.object(moments.measurement, 'e5_target_system_for_profile',
                return_value=SimpleNamespace(fingerprint='target')))
            stack.enter_context(mock.patch.object(moments.measurement, 'first_birth_accounting_by_age', return_value={}))
            stack.enter_context(mock.patch.object(baseline.pf.calendar, 'evaluate_period', return_value=object()))
            observer = stack.enter_context(mock.patch.object(moments, 'HistoricalMomentObserver'))
            observer.return_value.report.return_value = report
            evaluator = stack.enter_context(mock.patch.object(baseline.joined, 'evaluate_history_and_person_tail',
                return_value=result))
            stack.enter_context(mock.patch.object(baseline.joined, 'check_smoke_gates', return_value={'passed': True},
                side_effect=RuntimeError('failed full path gate') if gate_failure else None))
            stack.enter_context(mock.patch.object(terminal_checks, 'terminal_convergence_diagnostics',
                return_value={'all_checks_pass': False}))
            summary = baseline.run_history_probe(seed, contract, args, output, {},
                lambda name, data: baseline.pf.write_json(output/name, data), time.monotonic(),
                prices_override=np.ones(6))
        return summary, old, people, history_pre, evaluator.call_args.kwargs

    def test_checkpoint_is_2023_prechoice_not_final_state_and_not_equilibrium(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary, old, people, history_pre, inputs = self.run_probe(output)
            with gzip.open(output/'initial_2023.pkl.gz', 'rb') as stream:
                saved = pickle.load(stream)
            np.testing.assert_array_equal(saved['initial_state'].g_pre, history_pre)
            np.testing.assert_array_equal(saved['initial_state'].persons.persons, people.persons)
            np.testing.assert_array_equal(saved['initial_state'].persons.heads, people.heads)
            self.assertEqual(saved['initial_state'].persons.year, 2023)
            self.assertEqual(saved['parameters'].psi_child, old.psi_path[-1])
            self.assertEqual(old.parameters.psi_child, .3)
            self.assertEqual(vars(saved['supply_rule']), vars(old.supply_rule))
            self.assertEqual(saved['contract_sha256'], 'contract-hash')
            self.assertFalse(saved['equilibrium_certified'])
            self.assertFalse(saved['policy_announcement_included'])
            self.assertIn('before 2023 fertility/tenure choices', saved['state_timing'])
            self.assertEqual(summary['artifact_sha256']['initial_2023.pkl.gz'],
                driver.primitive.digest(output/'initial_2023.pkl.gz'))
            np.testing.assert_array_equal(inputs['transfer_path'], np.zeros(6))

    def test_disabled_optional_checkpoint_is_not_written(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary, *_ = self.run_probe(output, save_state=False)
            self.assertFalse((output/'initial_2023.pkl.gz').exists())
            self.assertNotIn('initial_2023.pkl.gz', summary['artifact_sha256'])

    def test_failed_full_path_gate_cannot_write_inherited_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, 'failed full path gate'):
                self.run_probe(output, gate_failure=True)
            self.assertFalse((output/'initial_2023.pkl.gz').exists())
            self.assertFalse((output/'summary.json').exists())


if __name__ == '__main__':
    unittest.main()
