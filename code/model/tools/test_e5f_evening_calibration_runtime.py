"""Run on Torch only: zero-solve contract/binding fixtures for evening runtime."""
import copy
import unittest
from types import SimpleNamespace as NS
from pathlib import Path
from unittest.mock import patch
import numpy as np
import e5f_evening_calibration_runtime as r

IDS=['initial_normalization','cps_childlessness','cps_exactly_one','nchs_mean_age','nchs_share30',
     'early_fertility','wealth_earnings','bequest_wealth','old_dispersion','mean_rooms',
     'ownership_30_55','first_birth_rooms','family_rooms','recent_parent_ownership']
def objective():
    return dict(target_rows=[dict(restriction_id=k,target=2.1 if k=='initial_normalization' else 1.,
           actual_weight=None if k=='initial_normalization' else 0. if k in r.VALIDATION else 1.) for k in IDS],
           parameter_restrictions=[dict(parameter=k,lower=.001 if k=='tenure_choice_kappa' else 0.,upper=.1 if k=='tenure_choice_kappa' else 10.) for k in r.FREE])

class EveningRuntime(unittest.TestCase):
    def test_stale_namespace_rejected(self):
        wrong=NS(__file__='/frozen/tools/e5f_overnight_estate_audit.py')
        with self.assertRaisesRegex(RuntimeError,'not current'):
            r.require_current_tool(wrong,'e5f_overnight_estate_audit')
        current=NS(__file__=str(r.ROOT/'code/model/tools/e5f_overnight_estate_audit.py'))
        self.assertIs(r.require_current_tool(current,'e5f_overnight_estate_audit'),current)

    def test_finite_scalar_gates(self):
        r.require_abs_gate(1e-11,2e-10,'budget')
        for value in (float('nan'),float('inf'),-float('inf'),1e-8):
            with self.assertRaises(RuntimeError):r.require_abs_gate(value,2e-10,'budget')

    def test_exact_counts_and_explicit_validation(self):
        o=objective();r.validate_objective(o)
        for key in r.VALIDATION:
            bad=copy.deepcopy(o)
            next(x for x in bad['target_rows'] if x['restriction_id']==key)['actual_weight']=1.
            with self.assertRaises(ValueError):r.validate_objective(bad)
        bad=copy.deepcopy(o);bad['parameter_restrictions'].pop()
        with self.assertRaises(ValueError):r.validate_objective(bad)

    def test_zero_weights_do_not_drop_rows(self):
        rows=[dict(moment=k,target=1.,model=3.,gap=2.,weight=1.,loss_contribution=4.) for k in IDS]
        with patch.object(r.base,'score_targets',return_value=rows):
            actual=r.score_targets(objective(),{}, {},0.,2.1)
        self.assertEqual(len(actual),14)
        self.assertEqual(sum(x['loss_contribution'] for x in actual if x['loss_contribution']!=''),40.)
        for row in actual:
            if row['moment'] in r.VALIDATION:self.assertEqual(row['loss_contribution'],0.)

    def test_binding_preserves_primitives_and_aliases(self):
        P=NS(H0=np.array([6.]),period_years=4,sigma=2.,alpha_cons=.733,child_room_floor=False,
             hbar_first_child_jump=0.,hbar_child_rooms=0.,delta_alpha=0.,phi=np.array([.8]),xi_supply=np.array([.63]),
             utility_reference_rent=.15,q=.08,tenure_choice_kappa=.005)
        e=object.__new__(r.EveningObjective);e.obj=objective();e.prepared={'parameters':P};e.c={'fixed':{'reference_rent':.15}}
        e.case_output=Path('/tmp/only-a-string-no-write');e.due=True
        point={k:.5 for k in r.FREE};point.update(beta_annual=.96,H0=7.,tenure_choice_kappa=.01)
        new=e.bind(point)
        self.assertEqual(new.beta,.96**4);self.assertEqual(new.rho,new.rho_hat)
        self.assertEqual(new.eps_fert,new.kappa_fert);self.assertTrue(new.native_due_stayer_credit)
        self.assertEqual(new.q,P.q);self.assertEqual(P.tenure_choice_kappa,.005)
        np.testing.assert_array_equal(P.H0,[6.]);np.testing.assert_array_equal(new.H0,[7.])
        point['tenure_choice_kappa']=.2
        with self.assertRaises(ValueError):e.bind(point)

if __name__=='__main__':unittest.main()
