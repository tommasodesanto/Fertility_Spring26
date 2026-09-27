"""Pure preparation checks: no model imports or solves."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import prepare_e5f_supervised_night as prep


def tables():
    objective = {'target_rows': [
        {'restriction_id': 'initial_normalization', 'target': 2.1, 'actual_weight': None},
        {'restriction_id': 'early_fertility', 'target': .81, 'actual_weight': 100.},
        {'restriction_id': 'mean_rooms', 'target': 5.73, 'actual_weight': 128.}],
        'parameter_restrictions': [{'parameter': 'x', 'lower': 0., 'upper': 1.}]}
    provenance = {'target_rows': [dict(id=r['restriction_id'], target=r['target'],
                                     weight=r['actual_weight']) for r in objective['target_rows']]}
    return objective, provenance


class PreparationTests(unittest.TestCase):
    def test_primary_preserves_every_input(self):
        o, p = tables(); before = copy.deepcopy((o, p))
        prep.weighting_update(o, p, 'primary')
        self.assertEqual((o, p), before)

    def test_early_experiment_changes_only_one_weight_and_label(self):
        o, p = tables(); before = copy.deepcopy((o, p))
        prep.weighting_update(o, p, 'early_fertility_3000')
        for rows, old, key, name in [(o['target_rows'], before[0]['target_rows'], 'actual_weight', 'restriction_id'),
                                    (p['target_rows'], before[1]['target_rows'], 'weight', 'id')]:
            for row, previous in zip(rows, old):
                if row[name] == 'early_fertility':
                    self.assertEqual(row[key], 3000.)
                    row = {k: v for k, v in row.items() if k not in (key, 'weight_status')}
                    previous = {k: v for k, v in previous.items() if k != key}
                self.assertEqual(row, previous)
        self.assertEqual(o['parameter_restrictions'], before[0]['parameter_restrictions'])

    def test_identity_keeps_normalization_unscored(self):
        o, p = tables(); prep.weighting_update(o, p, 'identity')
        self.assertEqual([r['actual_weight'] for r in o['target_rows']], [None, 1., 1.])
        self.assertEqual([r['weight'] for r in p['target_rows']], [None, 1., 1.])

    def test_missing_early_provenance_fails(self):
        o, p = tables(); p['target_rows'].pop(1)
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            prep.weighting_update(o, p, 'early_fertility_3000')

    def case(self, root):
        checkpoint = root / 'initial_state.pkl.gz'; checkpoint.write_bytes(b'fixture only')
        c = {'objective': {'sha256': 'objective'}, 'source_manifest': {'sha256': 'source'}}
        receipt = dict(status='verified_provisional_calibration_point', target_system_sha256='objective',
                       source_manifest_sha256='source', case_checkpoint_sha256=prep.sha(checkpoint), point={'x': .5})
        (root / 'receipt.json').write_text(json.dumps(receipt))
        return c, receipt

    def test_verified_start_retains_point_and_pins_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); c, receipt = self.case(root)
            point, pins = prep.initial_case(c, tables()[0], root)
            self.assertEqual(point, receipt['point'])
            self.assertEqual(pins['initial_case_receipt'], prep.pin(root / 'receipt.json'))

    def test_changed_checkpoint_and_invalid_bound_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); c, receipt = self.case(root)
            (root / 'initial_state.pkl.gz').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'provenance/checkpoint'):
                prep.initial_case(c, tables()[0], root)
            c, receipt = self.case(root); receipt['point']['x'] = 2.
            (root / 'receipt.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'outside bounds'):
                prep.initial_case(c, tables()[0], root)


if __name__ == '__main__':
    unittest.main()
