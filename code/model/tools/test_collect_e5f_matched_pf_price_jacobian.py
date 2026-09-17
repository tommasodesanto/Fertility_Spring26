"""Synthetic file receipts test the collector; no model computation."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import collect_e5f_matched_pf_price_jacobian as collector


class JacobianCollectorTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.A = np.array([[-2., .3], [.2, -1.]])
        self.F = np.array([.1, -.2])
        self.anchor = self.case(-1)
        self.probes = [self.case(j) for j in range(2)]

    def case(self, coordinate):
        folder = self.root / str(coordinate)
        folder.mkdir()
        prices = np.array([.6, .7])
        if coordinate >= 0:
            prices[coordinate] *= np.exp(.01)
        F = self.F.copy() if coordinate < 0 else self.F + .01 * self.A[:, coordinate]
        c = dict(arm='sequential', years=[2007, 2011], anchor_prices=[.6, .7],
            target_fingerprint='target1', path_date_count=2,
            normalized_checkpoint_sha256='a'*64, terminal_checkpoint_sha256='b'*64,
            checkpoint_sha256='c'*64, selected_summary_sha256='d'*64,
            demographic_sources={name: {'path': '/' + name + '.csv', 'sha256': 'e'*64}
                for name in ('population_mid','births_mid','survival','vintage_2025','acs_headship')},
            source_sha256={'solver.py': 'f'*64}, initial_price_rule='explicit_initial',
            terminal_preference_rule='explicit_terminal', psi_path=[.3, .2], transfer_path=[0., 0.],
            terminal_price=.8, supply_rule={'initial_price': .6, 'elasticity': .63},
            probe_log_step=.01, probe_coordinate=coordinate, prices=prices.tolist(),
            contract_sha256=str(coordinate+2)*64)
        with (folder/'transition_path.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['calendar_year','asset_price','housing_demand','housing_supply'])
            writer.writeheader()
            for year, price, residual in zip(c['years'], prices, F):
                writer.writerow(dict(calendar_year=year, asset_price=price, housing_demand=1+residual, housing_supply=1.))
        s = dict(status=collector.STATUS, arm=c['arm'], years=c['years'],
            prices=c['prices'], anchor_prices=c['anchor_prices'], residual=F.tolist(),
            maximum_market_residual=float(np.max(np.abs(F))), mapping_valid=True,
            gates={'mass': {'value': 0., 'tolerance': 1e-9, 'passed': True}},
            target_fingerprint=c['target_fingerprint'], contract_sha256=c['contract_sha256'],
            artifact_sha256={'transition_path.csv': collector.digest(folder/'transition_path.csv')})
        (folder/'contract.json').write_text(json.dumps(c))
        (folder/'summary.json').write_text(json.dumps(s))
        return folder

    def mutate(self, folder, filename, change):
        path = folder / filename
        obj = json.loads(path.read_text())
        change(obj)
        path.write_text(json.dumps(obj))

    def test_valid_coupled_signed_jacobian(self):
        result = collector.collect(self.anchor, list(reversed(self.probes)))
        np.testing.assert_allclose(result['jacobian'], self.A, atol=2e-14)
        np.testing.assert_allclose(result['residual'], self.F, atol=1e-15)
        self.assertAlmostEqual(result['condition_number'], np.linalg.cond(self.A))
        self.assertEqual([x['coordinate'] for x in result['provenance']['columns']], [0, 1])

    def test_mixed_fingerprint_rejected(self):
        for filename in ('contract.json', 'summary.json'):
            self.mutate(self.probes[0], filename, lambda x: x.update(target_fingerprint='other'))
        with self.assertRaisesRegex(ValueError, 'target_fingerprint'):
            collector.collect(self.anchor, self.probes)

    def test_missing_and_duplicate_coordinate_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Missing probe'):
            collector.collect(self.anchor, self.probes[:1])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            collector.collect(self.anchor, [self.probes[0], self.probes[0]])

    def test_projection_that_changes_other_price_rejected(self):
        for filename in ('contract.json', 'summary.json'):
            self.mutate(self.probes[0], filename, lambda x: x['prices'].__setitem__(1, .701))
        with self.assertRaisesRegex(ValueError, 'single-coordinate'):
            collector.collect(self.anchor, self.probes)

    def test_unsigned_residual_is_rejected_against_hashed_csv(self):
        self.mutate(self.anchor, 'summary.json', lambda x: x.update(residual=np.abs(x['residual']).tolist()))
        with self.assertRaisesRegex(ValueError, 'signed .*residual'):
            collector.collect(self.anchor, self.probes)

    def test_changed_transition_csv_hash_rejected(self):
        with (self.anchor/'transition_path.csv').open('a') as stream:
            stream.write('\n')
        with self.assertRaisesRegex(ValueError, 'hash'):
            collector.collect(self.anchor, self.probes)

    def test_claimed_pass_with_bad_gate_value_rejected(self):
        self.mutate(self.anchor, 'summary.json', lambda x: x['gates']['mass'].update(value=.1))
        with self.assertRaisesRegex(ValueError, 'Failed numerical gate'):
            collector.collect(self.anchor, self.probes)

    def test_mixed_source_and_missing_demographic_pin_rejected(self):
        self.mutate(self.probes[0], 'contract.json', lambda x: x['source_sha256'].update({'solver.py': 'a'*64}))
        with self.assertRaisesRegex(ValueError, 'source_sha256'):
            collector.collect(self.anchor, self.probes)
        self.mutate(self.anchor, 'contract.json', lambda x: x['demographic_sources'].pop('survival'))
        with self.assertRaisesRegex(ValueError, 'five demographic'):
            collector.collect(self.anchor, self.probes)


if __name__ == '__main__':
    unittest.main()
