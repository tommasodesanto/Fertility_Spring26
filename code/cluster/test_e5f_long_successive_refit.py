"""Small routing/fit tests; no household or equilibrium model solves."""
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from types import SimpleNamespace as NS

import numpy as np

import run_e5f_long_successive_refit as refit


class RefitContractTests(unittest.TestCase):
    def test_wall_budget_keeps_space_for_final_mapping_and_artifacts(self):
        manifest = dict(mapping_seconds_budget_104=3900, artifact_reserve_seconds=900)
        self.assertEqual(refit.budgeted_root_evaluations(9000, 104, 8, manifest), 2)
        self.assertEqual(refit.budgeted_root_evaluations(8000, 104, 8, manifest), 1)
        self.assertEqual(refit.budgeted_root_evaluations(40000, 104, 8, manifest), 8)
        self.assertEqual(refit.budgeted_root_evaluations(500, 104, 8, manifest), 0)

    def test_recovery_rejects_changed_native_coordinates_and_calendar(self):
        with TemporaryDirectory() as directory:
            p = Path(directory)
            rows = [dict(calendar_year=2007, psi_child=.12, asset_price=.7,
                         pension_period_units=2., equal_transfer_period_units=.18)]
            root = dict(calendar_year=2007, count=1, psi=.12,
                        best=dict(mapping_valid=True, prices=[.7, 2., .18]))
            (p/'rows.json').write_text(json.dumps(rows)); (p/'root.json').write_text(json.dumps(root))
            source = dict(psi=.12, native_rows=str(p/'rows.json'), root_receipt=str(p/'root.json'))
            self.assertEqual(refit.recovery_root(source, 2007, 1), root)
            with self.assertRaises(ValueError): refit.recovery_root(source, 2011, 1)
            rows[0]['asset_price'] = .71
            (p/'rows.json').write_text(json.dumps(rows))
            with self.assertRaises(ValueError): refit.recovery_root(source, 2007, 1)

    def test_insufficient_second_round_preserves_first_round_receipt(self):
        with TemporaryDirectory() as directory:
            manifest = dict(max_rounds=4, max_root_evaluations=8, reserve_verification_time=True,
                            mapping_seconds_budget_104=3900, artifact_reserve_seconds=900)
            c = NS(manifest=manifest, driver=NS())
            terminal = NS(coordinates=np.array([.7,2.,.18]))
            result = NS(path=None, next_state=None, root_receipt={'finite_horizon_market_fiscal_converged':False})
            root = dict(finite_horizon_market_fiscal_converged=False,
                        best=dict(prices=[.7]*104+[2.]*104+[.18]*104))
            saved = {}
            with patch.object(refit, 'endpoint', return_value=terminal), \
                 patch.object(refit, 'run_mapping', return_value=(result, [], root)) as mapping, \
                 patch.object(refit, 'budgeted_root_evaluations', side_effect=[2,1]), \
                 patch.object(refit, 'save', side_effect=lambda driver,path,value:saved.update({Path(path).name:value})):
                _,_,record,_ = refit.solve_candidate(c,None,None,inherited=None,psi=.12,year=2007,
                    count=104,target=1.975,folder=directory,deadline=refit.time.monotonic()+10000,endpoint_start=None)
            self.assertEqual(mapping.call_count,1)
            self.assertEqual(record['root_receipt'],root)
            self.assertEqual(record['rounds'],1)
            self.assertIn('candidate.json',saved)

    def test_same_terminal_boundary_leaves_matched_policy_tail(self):
        self.assertEqual([refit.horizon(y) for y in (2007, 2011, 2015, 2019, 2023)],
                         [104, 103, 102, 101, 100])
        self.assertEqual(refit.horizon(2019) - 1, refit.horizon(2023))
        with self.assertRaises(ValueError):
            refit.horizon(2008)
        with self.assertRaises(ValueError):
            refit.horizon(2423)

    def test_initial_seed_and_bounds(self):
        self.assertAlmostEqual(refit.next_psi([], .12, (-.05, .17), []), .12)
        self.assertAlmostEqual(refit.next_psi([], .25, (-.05, .17), []), .17)

    def test_bracketed_fit_stays_inside_bracket(self):
        trials = [dict(psi=.08, gap=-.4), dict(psi=.14, gap=.001)]
        x = refit.next_psi(trials, .12, (-.05, .17), [.08, .14])
        self.assertGreaterEqual(x, .08 + .15 * .06 - 1e-12)
        self.assertLessEqual(x, .14 - .15 * .06 + 1e-12)

    def test_failed_candidate_is_not_a_fertility_observation(self):
        trials = [dict(psi=.10, gap=None), dict(psi=.11, gap=float('nan'))]
        x = refit.next_psi(trials, .12, (-.05, .17), [.10, .11, .12])
        self.assertTrue(np.isfinite(x))
        self.assertNotIn(x, [.10, .11, .12])

    def test_scaled_solver_patch_keeps_factory_interface(self):
        original = lambda: 'original'
        actual = lambda **kwargs: kwargs
        c = NS(rebated=NS(_path_root_solver=original))
        with refit.scaled_root_context(c, NS(solve_price_path_scaled=actual)):
            self.assertIs(c.rebated._path_root_solver(), actual)
        self.assertIs(c.rebated._path_root_solver, original)

    def test_unconverged_forecast_cannot_be_accepted(self):
        result = NS(path=NS(), next_state=NS(), root_receipt={
            'finite_horizon_market_fiscal_converged': False})
        ok, gap = refit.accepted(result, [{'period_tfr_topcode_adjusted': 1.974875}], 1.974875)
        self.assertFalse(ok)
        self.assertIsNone(gap)

    def test_valid_wrong_fit_is_available_to_scalar_search(self):
        result = NS(path=NS(), next_state=NS(), root_receipt={
            'finite_horizon_market_fiscal_converged': True})
        ok, gap = refit.accepted(result, [{'period_tfr_topcode_adjusted': 1.9}], 1.8)
        self.assertFalse(ok)
        self.assertAlmostEqual(gap, .1)


if __name__ == '__main__':
    unittest.main()
