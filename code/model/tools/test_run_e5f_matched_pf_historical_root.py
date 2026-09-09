"""Historical-root provenance tests; no model or cluster evaluations."""
from __future__ import annotations

import copy
import unittest
from unittest import mock

import numpy as np

import collect_e5f_matched_pf_price_jacobian as collector
import run_e5f_matched_pf_historical_root as driver


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
            scope='explicit supplied-price hook only; unchanged economic evaluator')

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


if __name__ == '__main__':
    unittest.main()
