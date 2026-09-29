#!/usr/bin/env python3
"""Pure tests for the fixed-price pension-accounting diagnostic."""
import importlib.util
import math
import numpy as np
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

PATH = Path(__file__).resolve().parents[2] / 'cluster' / 'run_e5f_preference_budget_diagnostic.py'
SPEC = importlib.util.spec_from_file_location('budget_diagnostic', PATH)
driver = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(driver)
TRANSITION_PATH = Path(__file__).with_name('run_e5f_preference_transition.py')
TRANSITION_SPEC = importlib.util.spec_from_file_location('budget_diagnostic_transition', TRANSITION_PATH)
transition = importlib.util.module_from_spec(TRANSITION_SPEC); TRANSITION_SPEC.loader.exec_module(transition)


class WriterTests(unittest.TestCase):
    def test_numpy_values_roundtrip_and_invalid_values_reject(self):
        path = Path(tempfile.mkdtemp()) / 'receipt.json'
        value = dict(nested=[np.array([1., 2.]), dict(float=np.float64(.5), integer=np.int64(4), boolean=np.bool_(True))])
        driver.write(path, value)
        self.assertEqual(driver.read(path), {'nested': [[1., 2.], {'boolean': True, 'float': .5, 'integer': 4}]})
        with self.assertRaises(ValueError): driver.write(path, {'nonfinite': np.float64(np.nan)})
        with self.assertRaises(TypeError): driver.write(path, {'unsupported': object()})


class BudgetUpdateTests(unittest.TestCase):
    def test_exact_balance_and_nonmutation(self):
        rows = [dict(payroll_tax_revenue=12., pension_outlays=3., pension_period_units=2., asset_price=1.)]
        self.assertEqual(driver.budget_update(rows), [8.]); self.assertEqual(rows[0]['pension_period_units'], 2.)
        self.assertAlmostEqual(driver.budget_update(rows, .5)[0], 4.)

    def test_invalid_units_rejected(self):
        good = dict(payroll_tax_revenue=1., pension_outlays=1., pension_period_units=1., asset_price=1.)
        for key, value in [('payroll_tax_revenue', 0.), ('pension_outlays', -1.), ('asset_price', math.nan)]:
            row = good.copy(); row[key] = value
            with self.assertRaises(ValueError): driver.budget_update([row])


class IdentityTests(unittest.TestCase):
    def test_full_baseline_identity_guards(self):
        ref = dict(q=1., psi=.1, pension=.2)
        fiscal = 1.01787e-6
        revenue = 1.
        outlays = 1. - fiscal
        calculated_fiscal = (revenue - outlays) / max(abs(revenue), abs(outlays))
        rows = [dict(calendar_year=2007+4*t, asset_price=1., psi_child=.1, pension_period_units=.2,
                     housing_demand=1., housing_supply=1., scaled_pension_budget_residual=calculated_fiscal,
                     payroll_tax_revenue=revenue, pension_outlays=outlays) for t in range(104)]
        record = dict(backward_forward_policy_error=0., cache={}, diagnostic_packets=[], fertility=[{}]*104,
                      final_mass=1., fiscal_residual=[calculated_fiscal]*104,
                      gates=dict(mass=True, policy_reproduction=True, projection=True, dated_audits=True), initial_mass=1.,
                      market_residual=[0.]*104, mass_error=0., population_l1=0., projection_mass=0., rows=rows, seconds=1.)
        driver.validate_baseline(record, 104, ref)
        record['rows'][0]['calendar_year'] = 2008
        with self.assertRaises(ValueError): driver.validate_baseline(record, 104, ref)

    def test_baseline_rejects_missing_gates_changed_price_or_residual_mismatch(self):
        ref = dict(q=1., psi=.1, pension=.2)
        row = dict(calendar_year=2007, asset_price=1., psi_child=.1, pension_period_units=.2, housing_demand=1., housing_supply=1.,
                   payroll_tax_revenue=1., pension_outlays=1., scaled_pension_budget_residual=0.)
        record = dict(backward_forward_policy_error=0., cache={}, diagnostic_packets=[], fertility=[{}], final_mass=1., fiscal_residual=[0.],
                      gates=dict(mass=True, policy_reproduction=True, projection=True, dated_audits=True), initial_mass=1., market_residual=[0.],
                      mass_error=0., population_l1=0., projection_mass=0., rows=[row], seconds=1.)
        for mutate in (lambda r: r['gates'].pop('mass'), lambda r: r['rows'][0].update(asset_price=2.),
                       lambda r: r.update(fiscal_residual=[1e-8])):
            bad = {k: (v.copy() if isinstance(v, dict) else list(v) if isinstance(v, list) else v) for k, v in record.items()}
            bad['rows'] = [row.copy()]; bad['gates'] = record['gates'].copy(); mutate(bad)
            with self.assertRaises(ValueError): driver.validate_baseline(bad, 1, ref)

    def test_mode_horizon_and_required_fields(self):
        config = dict(mode='smoke', reference_manifest_sha256='147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4',
            source_pins={name:'a'*64 for name in driver.SOURCE_NAMES}, baseline_mapping=None,
            endpoint_receipt={'path':'nope','sha256':'a'*64}, horizon=7, cache_max_bytes=0, mapping_seconds=1, total_seconds=1,
            numerical_threads=None, market_tolerance=2e-4, fiscal_tolerance=1e-6, replay_tolerance=1e-10, smoke_pension_factor=1.00001)
        config.pop('baseline_mapping')
        with patch.object(driver, 'pinned', return_value=Path('x')):
            with self.assertRaises(ValueError): driver.validate_config(config)


class ExactLoopTests(unittest.TestCase):
    def config(self):
        return dict(mode='smoke', horizon=6, total_seconds=10, mapping_seconds=5, cache_max_bytes=0,
                    smoke_pension_factor=1.00001)

    def row(self, pension=1.):
        return dict(calendar_year=2007, asset_price=1., psi_child=.1, pension_period_units=pension,
                    payroll_tax_revenue=2., pension_outlays=1., housing_demand=1., housing_supply=1.,
                    scaled_pension_budget_residual=0.)

    def record(self, fiscal=0.):
        rows=[self.row() for _ in range(6)]
        for row in rows: row['scaled_pension_budget_residual']=fiscal
        return dict(rows=rows, fertility=[dict(period_tfr_topcode_adjusted=.5)]*6, market_residual=[0.]*6, fiscal_residual=[fiscal]*6,
                    cache={'actual_solves':1,'hits':2}, gates={'all':True})

    def make(self):
        d=driver.Diagnostic(self.config(), Path(tempfile.mkdtemp()), SimpleNamespace(plain=transition.plain))
        d.setup=lambda: setattr(d, 'reference', dict(q=1., psi=.1, pension=1.))
        d.evaluator=SimpleNamespace(rt={'primitive':SimpleNamespace(pf=SimpleNamespace(birth_queue_values=lambda x: np.asarray(x)))})
        return d

    def test_fresh_replay_is_required_at_same_coordinates(self):
        d=self.make(); calls=[]
        def mapping(name, pensions):
            calls.append((name, list(pensions))); record=self.record(); terminal={'all_checks_pass':True}
            summary=dict(accepted=True, max_market_error=0., max_fiscal_error=0., residual_ratios={}, cache={}, gates={}, terminal=terminal)
            state=type('State', (), dict(g_pre=np.array([[0.]]), scheduled_entries=[0.], scheduled_raw_entries=[0.]))()
            return type('Result', (), dict(terminal_state=state))(),record,terminal,summary
        d.mapping=mapping
        result=d.run()
        self.assertTrue(result['numerical_certified']); self.assertEqual([x[0] for x in calls], ['baseline','trial','fresh_replay'])
        self.assertEqual(calls[1][1], calls[2][1])

    def test_trial_fiscal_failure_stops_without_damping_search(self):
        d=self.make(); calls=[]
        def mapping(name, pensions):
            calls.append(name); record=self.record(2e-6 if name == 'trial' else 0.); terminal={'all_checks_pass':True}
            summary=dict(accepted=name != 'trial', max_market_error=0., max_fiscal_error=abs(record['fiscal_residual'][0]), residual_ratios={}, cache={}, gates={}, terminal=terminal)
            state=type('State', (), dict(g_pre=np.array([[0.]]), scheduled_entries=[0.], scheduled_raw_entries=[0.]))()
            return type('Result', (), dict(terminal_state=state))(),record,terminal,summary
        d.mapping=mapping
        result=d.run()
        self.assertFalse(result['numerical_certified']); self.assertEqual(calls, ['baseline','trial'])

    def test_real_mapping_writes_numpy_terminal_and_replay_receipts(self):
        d=self.make(); calls=[]
        def native_mapping(*_args, **_kwargs):
            state=SimpleNamespace(g_pre=np.zeros((7,)), scheduled_entries=np.arange(7.), scheduled_raw_entries=np.arange(7.))
            record=self.record()
            for row in record['rows']:
                row['payroll_tax_revenue']=1.; row['pension_outlays']=1.
            return SimpleNamespace(terminal_state=state), record
        def terminal_checks(*_args, **_kwargs):
            calls.append(True)
            return dict(all_checks_pass=True, terminal_birth_queue=np.arange(7.), stationary_birth_queue=np.arange(7.),
                        native_flag=np.bool_(True), count=np.int64(7), scalar=np.float64(.25))
        d.inner.mapping=native_mapping; d.inner.terminal_checks=terminal_checks
        result=d.run()
        self.assertTrue(result['numerical_certified']); self.assertEqual(len(calls), 3)
        latest=driver.read(d.output / 'latest_completed.json')
        self.assertEqual(latest['terminal']['terminal_birth_queue'], list(np.arange(7.)))
        self.assertEqual(latest['terminal']['stationary_birth_queue'], list(np.arange(7.)))
        self.assertEqual(latest['terminal']['native_flag'], True)
        self.assertEqual(latest['terminal']['count'], 7)
        self.assertEqual(len(driver.read(d.output / 'best_so_far.json')['terminal']['terminal_birth_queue']), 7)
        driver.write(d.output / 'completed_result.json', result)
        self.assertEqual(driver.read(d.output / 'completed_result.json')['outcome'], 'certified')


if __name__ == '__main__': unittest.main()
