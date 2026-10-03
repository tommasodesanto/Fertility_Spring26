"""Synthetic, zero-solve estate-A audit adapter checks."""
from pathlib import Path
import copy
import importlib
import importlib.util
import sys
import unittest
from types import SimpleNamespace
import numpy as np

EXPERIMENT = Path(__file__).resolve().parents[1]
ROOT = EXPERIMENT.parents[3]
sys.path.insert(0, str(EXPERIMENT.parent))
sys.path.insert(0, str(ROOT/'code/model/tools'))
adapter = importlib.import_module('birth_count_choice.model.estate_audit_adapter')


def load_audit(relative):
    spec = importlib.util.spec_from_file_location('estate_a_test_frozen', ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(**overrides):
    P = SimpleNamespace(I=1, J=2, H_own=np.array([100.]), period_years=4.,
        psi=.06, use_age_survival=True, survival_probs=np.ones(1),
        estate_receiver='none', bequest_net_of_selling_cost=True,
        estate_flow_net_of_selling_cost=True)
    vars(P).update(overrides)
    g = np.zeros((1, 2, 1, 2, 1, 1, 1))
    g[0, 0, 0, 0, 0, 0, 0] = 1.
    g[0, 1, 0, 1, 0, 0, 0] = 1.
    bp = np.zeros_like(g)
    bp[0, 1, 0, 1, 0, 0, 0] = -40.
    ev = SimpleNamespace(g_current=g, g_pre=g.copy(),
        policy=SimpleNamespace(bp_pol=bp, price=np.ones(1)))
    return ev, P, np.array([5.])


class EstateAAuditTests(unittest.TestCase):
    def setUp(self):
        self.module = load_audit('code/model/tools/e5f_overnight_estate_audit.py')
        self.audit, self.receipt = adapter.adapt_audit(self.module.audit)

    def test_net_estate_accounting_unchanged_metadata_and_no_mutation(self):
        ev, P, bg = fixture()
        params = copy.deepcopy(vars(P)); mass = ev.g_current.copy(); bp = ev.policy.bp_pol.copy()
        result = self.audit(ev, P, bg)
        self.assertEqual(result['available_estates_period'], 54.)
        self.assertTrue(result['donor_bequest_valuation_changes'])
        self.assertEqual(result['estate_a_identity'], adapter.IDENTITY)
        self.assertIn('no extra interest', result['donor_estate_valuation'])
        self.assertFalse(result['adult_transfers'])
        self.assertFalse(result['certifies_counterparty_or_physical_housing_settlement'])
        self.assertFalse(any('donor utility remains unchanged' in c for c in result['caveats']))
        self.assertTrue(self.receipt['all_accounting_arithmetic_and_numerical_gates_unchanged'])
        np.testing.assert_array_equal(mass, ev.g_current)
        np.testing.assert_array_equal(bp, ev.policy.bp_pol)
        for key in params:
            np.testing.assert_array_equal(params[key], getattr(P, key))
        P.bequest_net_of_selling_cost = P.estate_flow_net_of_selling_cost = False
        baseline = self.module.audit(ev, P, bg)
        for key in ('estate','entry','available_estates_period','funding_gate_tolerance',
                    'net_residual_after_entry_period','residual_sink_period'):
            self.assertEqual(result[key], baseline[key])
        self.assertEqual(self.audit(ev, P, bg), baseline)

    def test_shortfall_retains_dedicated_exception_and_labels_ledger(self):
        ev, P, _ = fixture()
        with self.assertRaises(self.module.EstateFundingShortfall) as error:
            self.audit(ev, P, np.array([60.]))
        ledger = error.exception.audit
        self.assertEqual(ledger['funding_shortfall_period'], 6.)
        self.assertIsNone(ledger['residual_sink_period'])
        self.assertTrue(ledger['donor_bequest_valuation_changes'])
        self.assertEqual(ledger['estate_a_identity'], adapter.IDENTITY)

    def test_partial_flags_transfers_and_other_gates_rejected(self):
        for name, value in (('estate_receiver','ages_45_65'), ('estate_lump_sum_transfer',.1),
                            ('estate_probe_transfer',.1), ('estate_tax_rate',.1),
                            ('bequest_net_of_selling_cost',False),
                            ('estate_flow_net_of_selling_cost',False),
                            ('estate_flow_net_of_selling_cost',1),
                            ('use_postdecision_current_distribution',False)):
            ev, P, bg = fixture(**{name:value})
            with self.assertRaises(ValueError): self.audit(ev, P, bg)
        ev, P, bg = fixture()
        ev.g_pre[0, 0, 0, 0, 0, 0, 0] = 0.
        ev.g_pre[0, 1, 0, 0, 0, 0, 0] = 1.
        with self.assertRaises(ValueError): self.audit(ev, P, bg)

    def test_actual_pinned_audit_supported_and_default_exact(self):
        frozen = load_audit('tmp/e5f_overnight_local_20260927/portable/tools_v4/e5f_overnight_estate_audit.py')
        audit, receipt = adapter.adapt_audit(frozen.audit)
        ev, P, bg = fixture()
        self.assertEqual(audit(ev, P, bg)['available_estates_period'], 54.)
        self.assertEqual(receipt['original_file_sha256'], '31bcdbeca73036da9306dcb2f3ef28a2c3ee43d878f60b75aa044e08aba5f962')
        P.bequest_net_of_selling_cost = P.estate_flow_net_of_selling_cost = False
        self.assertEqual(audit(ev,P,bg), frozen.audit(ev,P,bg))


if __name__ == '__main__':
    unittest.main()
