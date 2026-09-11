"""Pure scorer checks. Unit-test weights below are synthetic, never adopted."""
import copy
import csv
from pathlib import Path
import unittest

import score_initial as score

ROOT = Path(__file__).resolve().parent


def fixture():
    # Reuse the complete saved metadata; only the required sample declaration
    # and weights are synthetic unit-test inputs, not a proposed calibration.
    with (ROOT.parent/'initial_fit_readout/target_fit.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row['target'] = float(row['target'])
        row['sample'] = 'Synthetic unit fixture; authoritative sample supplied by caller in real use'
        row['actual_weight'] = None if row['restriction_id'] == 'initial_normalization' else 1.
        if row['restriction_id'] == 'recent_parent_ownership':
            row['model_observation'] = score.RECENT_MOMENT
    with (ROOT.parent/'initial_fit_readout/parameters.csv').open() as stream:
        parameter_rows = list(csv.DictReader(stream))
    restrictions = [{key: (float(row[key]) if key in ('lower', 'upper') else row[key])
                     for key in ('parameter', 'lower', 'upper', 'transform')}
                    for row in parameter_rows if row['parameter'] in score.PARAMETERS]
    approximation = {'approximation_id': 'unit_fixture_maintained_approximation',
                     'maintained': True, 'exact_acs': False,
                     'declaration': 'Synthetic certification fixture; no empirical equivalence claim',
                     'observer_id': score.RECENT_OBSERVER, 'moment': score.RECENT_MOMENT,
                     'snapshot': 'synchronized_post_fertility_snapshot',
                     'age_projection': 'uniform_within_age_cell',
                     'diagnostic_allow_residence_proxy': True}
    contract = {'schema': score.SCHEMA, 'contract_id': 'synthetic_unit_test_only',
                'objective_name': 'synthetic_known_score', 'cps_projection': 'uniform_birth_time',
                'source_fingerprints': {'initial_science': 'a'*64, 'recent_observer': 'b'*64},
                'target_rows': rows, 'parameter_restrictions': restrictions,
                'near_bound_fraction': .01, 'normalization_tolerance': 5e-4,
                'recent_parent_approximation': approximation}
    early = {'calibrated_smm': False, 'weights_assigned': False,
             'empirical_target_contract_activated': False,
             'fertility': {p: {'metadata': {'age_projection': p}, 'moments': {}}
                           for p in ('uniform_birth_time', 'constant_post_cell')},
             'housing_wealth': {'moments': {'recent_parent_minus_no_resident_child_ownership_30_55': None}}}
    for row in rows:
        name = row['restriction_id']
        if name in ('initial_normalization', 'recent_parent_ownership'):
            continue
        for projection in ('uniform_birth_time', 'constant_post_cell'):
            path = score.MOMENT_PATHS[name].format(projection=projection).split('.')
            dest = early
            for key in path[:-1]:
                dest = dest[key]
            dest[path[-1]] = row['target']
    checkpoint = 'c'*64
    recent = {'observer_id': score.RECENT_OBSERVER, 'moment': score.RECENT_MOMENT,
              'available': True, 'status': 'diagnostic_only_measurement_approximations_unresolved',
              'production_eligible': False, 'target_contract_activated': False,
              'model_value': score.RECENT_TARGET,
              'metadata': {'snapshot': approximation['snapshot'], 'age_projection': approximation['age_projection'],
                           'diagnostic_allow_residence_proxy': True, 'age_interval': [30., 56.],
                           'policy_input_provenance': {'checkpoint_sha256': checkpoint}},
              'groups': {'selected_birth': {'denominator': 1., 'owner_numerator': score.RECENT_TARGET,
                                             'ownership_rate': score.RECENT_TARGET},
                         'current_empty': {'denominator': 1., 'owner_numerator': 0., 'ownership_rate': 0.}}}
    normalization = {'target': 2.1, 'completed_fertility': 2.1001,
                     'absolute_gap': abs(2.1001-2.1),
                     'psi_child': float(next(r['estimate'] for r in parameter_rows if r['parameter'] == 'psi_child'))}
    receipt = {'schema': 'e5f_initial_score_receipt_v1', 'status': 'verified',
               'source_fingerprints': copy.deepcopy(contract['source_fingerprints']),
               'checkpoint_sha256': checkpoint, 'numerical_gates_verified': True,
               'recent_parent_certified': True, 'recent_parent_approximation_id': approximation['approximation_id']}
    inputs = {'early_measurement': early, 'recent_parent_observation': recent,
              'normalization': normalization, 'parameters': parameter_rows}
    receipt['input_sha256'] = {key: score.fingerprint(value) for key, value in inputs.items()}
    return contract, inputs, receipt


class ScoreInitialTests(unittest.TestCase):
    def setUp(self):
        self.contract, self.inputs, self.receipt = fixture()

    def call(self, *, repin_inputs=False, pin=None):
        if repin_inputs:
            self.receipt['input_sha256'] = {k: score.fingerprint(v) for k,v in self.inputs.items()}
        return score.score_initial(self.contract, expected_contract_sha256=pin or score.fingerprint(self.contract),
                                   evaluation_receipt=self.receipt, **self.inputs)

    def target(self, name):
        return next(row for row in self.contract['target_rows'] if row['restriction_id'] == name)

    def test_exact_known_score_and_all_rows_retained(self):
        housing = self.inputs['early_measurement']['housing_wealth']['moments']
        housing['annual_bequest_flow_to_aggregate_wealth'] += .5
        housing['aggregate_mean_occupied_rooms_capped9_18_85'] += 2.
        self.target('bequest_wealth')['actual_weight'] = 4.
        self.target('mean_rooms')['actual_weight'] = 3.
        before = copy.deepcopy((self.contract, self.inputs))
        result = self.call(repin_inputs=True)
        self.assertEqual(result['loss'], 13.)  # 4*(.5)^2 + 3*(2)^2
        self.assertEqual(len(result['target_fit']), 13)
        self.assertEqual(len(result['parameters']), 17)
        self.assertEqual(result['scored_moment_count'], 12)
        self.assertEqual((self.contract, self.inputs), before)
        self.assertFalse(result['identification_established'])
        self.assertFalse(result['benchmark_certified'])

    def test_normalization_never_scored(self):
        first = self.call()
        self.inputs['normalization'].update(completed_fertility=2.1004, absolute_gap=abs(2.1004-2.1))
        second = self.call(repin_inputs=True)
        self.assertEqual(first['loss'], second['loss'])
        self.assertEqual(first['loss'], 0.)
        row = next(x for x in second['target_fit'] if x['restriction_id'] == 'initial_normalization')
        self.assertIsNone(row['loss_contribution']); self.assertIsNone(row['actual_weight'])

    def test_normalization_failure_is_not_a_score_penalty(self):
        self.inputs['normalization'].update(completed_fertility=2.102, absolute_gap=.002)
        with self.assertRaisesRegex(score.ScoreContractError, 'normalization failed'):
            self.call(repin_inputs=True)

    def test_contract_mutation_rejects_original_fingerprint(self):
        original_pin = score.fingerprint(self.contract)
        self.target('mean_rooms')['actual_weight'] *= 2
        with self.assertRaisesRegex(score.ScoreContractError, 'Contract fingerprint mismatch'):
            self.call(pin=original_pin)

    def test_observation_mutation_rejects_receipt_fingerprint(self):
        self.inputs['early_measurement']['housing_wealth']['moments']['housing_increment_0to1'] += .1
        with self.assertRaisesRegex(score.ScoreContractError, 'input fingerprint mismatch'):
            self.call()

    def test_missing_duplicate_and_unknown_rows_rejected(self):
        for kind in ('missing', 'duplicate', 'unknown', 'fewer_than_nine'):
            with self.subTest(kind=kind):
                self.setUp()
                if kind == 'missing': self.contract['target_rows'].pop()
                elif kind == 'duplicate': self.contract['target_rows'].append(copy.deepcopy(self.contract['target_rows'][-1]))
                elif kind == 'unknown': self.contract['target_rows'][-1]['restriction_id'] = 'substituted_target'
                else: self.contract['target_rows'] = self.contract['target_rows'][:8]
                with self.assertRaises(score.ScoreContractError): self.call()

    def test_missing_provenance_rejected(self):
        for field in score.PROVENANCE_FIELDS:
            with self.subTest(field=field):
                self.setUp(); self.target('family_rooms')[field] = ''
                with self.assertRaisesRegex(score.ScoreContractError, field): self.call()

    def test_nonpositive_missing_and_nonfinite_weights_rejected(self):
        for value in (0., -1., None, True, float('nan'), float('inf')):
            with self.subTest(value=value):
                self.setUp(); self.target('bequest_wealth')['actual_weight'] = value
                with self.assertRaises(score.ScoreContractError): self.call()

    def test_missing_or_null_model_is_not_silently_omitted(self):
        for remove in (True, False):
            self.setUp(); moments = self.inputs['early_measurement']['housing_wealth']['moments']
            if remove: del moments['housing_increment_0to1']
            else: moments['housing_increment_0to1'] = None
            with self.assertRaises(score.ScoreContractError): self.call(repin_inputs=True)

    def test_source_and_recent_certificate_required(self):
        for field, value in [('source_fingerprints', {'initial_science': 'd'*64}),
                             ('recent_parent_certified', False), ('numerical_gates_verified', False),
                             ('recent_parent_approximation_id', 'another_approximation')]:
            with self.subTest(field=field):
                self.setUp(); self.receipt[field] = value
                with self.assertRaises(score.ScoreContractError): self.call()

    def test_recent_target_and_named_approximation_cannot_be_substituted(self):
        for change in ('target', 'exact_acs', 'not_maintained', 'legacy_observer', 'unavailable'):
            with self.subTest(change=change):
                self.setUp()
                if change == 'target': self.target('recent_parent_ownership')['target'] = .17
                elif change == 'exact_acs': self.contract['recent_parent_approximation']['exact_acs'] = True
                elif change == 'not_maintained': self.contract['recent_parent_approximation']['maintained'] = False
                elif change == 'legacy_observer': self.inputs['recent_parent_observation']['moment'] = 'old_any_dependent_gap'
                else: self.inputs['recent_parent_observation']['available'] = False
                with self.assertRaises(score.ScoreContractError): self.call(repin_inputs=True)

    def test_recent_checkpoint_and_group_arithmetic_checked(self):
        for change in ('checkpoint', 'ratio', 'difference', 'empty_group'):
            with self.subTest(change=change):
                self.setUp(); recent = self.inputs['recent_parent_observation']
                if change == 'checkpoint': recent['metadata']['policy_input_provenance']['checkpoint_sha256'] = 'e'*64
                elif change == 'ratio': recent['groups']['selected_birth']['owner_numerator'] += .01
                elif change == 'difference': recent['model_value'] += .01
                else: recent['groups']['selected_birth']['denominator'] = 0.
                with self.assertRaises(score.ScoreContractError): self.call(repin_inputs=True)

    def test_explicit_projection_selects_only_its_named_observer(self):
        self.inputs['early_measurement']['fertility']['constant_post_cell']['moments']['childless_rate_40_44'] += .25
        self.assertEqual(self.call(repin_inputs=True)['loss'], 0.)
        self.contract['cps_projection'] = 'constant_post_cell'
        for row in self.contract['target_rows']:
            row['model_observation'] = row['model_observation'].replace('uniform_birth_time', 'constant_post_cell')
        self.assertEqual(self.call()['loss'], .0625)
        self.contract['cps_projection'] = ['uniform_birth_time', 'constant_post_cell']
        with self.assertRaisesRegex(score.ScoreContractError, 'one explicit CPS projection'): self.call()

    def test_wrong_mapping_and_projection_metadata_rejected(self):
        self.target('cps_exactly_one')['model_observation'] = score.MOMENT_PATHS['cps_childlessness'].format(projection='uniform_birth_time')
        with self.assertRaisesRegex(score.ScoreContractError, 'mapping'): self.call()
        self.setUp(); self.inputs['early_measurement']['fertility']['uniform_birth_time']['metadata']['age_projection'] = 'constant_post_cell'
        with self.assertRaisesRegex(score.ScoreContractError, 'metadata mismatch'): self.call(repin_inputs=True)

    def test_parameters_bounds_duplicates_and_nonfinite_rejected(self):
        for change in ('missing', 'duplicate', 'outside', 'bound', 'nonfinite'):
            with self.subTest(change=change):
                self.setUp(); parameters = self.inputs['parameters']
                if change == 'missing': parameters.pop(0)
                elif change == 'duplicate': parameters.append(copy.deepcopy(parameters[0]))
                elif change == 'outside': parameters[0]['estimate'] = '1.0'
                elif change == 'bound': parameters[0]['lower'] = '.90'
                else: parameters[0]['estimate'] = 'nan'
                with self.assertRaises(score.ScoreContractError): self.call(repin_inputs=True)

    def test_current_model_receipt_replaces_stale_baseline_metadata(self):
        result = self.call()
        recent = next(row for row in result['target_fit'] if row['restriction_id'] == 'recent_parent_ownership')
        self.assertTrue(recent['model_available'])
        self.assertIsNone(recent['model_source_path'])
        self.assertEqual(recent['model_checkpoint_sha256'], self.receipt['checkpoint_sha256'])
        self.assertEqual(recent['evaluation_input'], 'recent_parent_observation')
        self.assertEqual(recent['model_source_json_location'], 'model_value')
        self.assertEqual(recent['definition'], self.target('recent_parent_ownership')['definition'])

    def test_normalized_preference_matches_full_parameter_receipt(self):
        self.inputs['normalization']['psi_child'] += .01
        with self.assertRaisesRegex(score.ScoreContractError, 'Normalized preference differs'):
            self.call(repin_inputs=True)

    def test_overflow_is_rejected_instead_of_infinite_objective(self):
        self.target('wealth_earnings')['target'] = -1e308
        self.inputs['early_measurement']['housing_wealth']['moments']['aggregate_wealth_to_annual_gross_labor_earnings'] = 1e308
        with self.assertRaises(score.ScoreContractError): self.call(repin_inputs=True)


if __name__ == '__main__':
    unittest.main()
