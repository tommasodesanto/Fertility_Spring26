"""Optional initial step preserves the normalizer cache and default path."""
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import run_e5f_transition_calibration as calibration


class NormalizationStepTests(unittest.TestCase):
    def solve(self, step=None):
        calls=[]
        def solve_ge(chain, overrides):
            psi=overrides['psi_child']; calls.append(psi)
            return SimpleNamespace(tfr=2.1+psi-.07), SimpleNamespace(), np.array([1.]), 3.
        chain=SimpleNamespace(extract_moments=lambda s,p:{'tfr':s.tfr})
        options={} if step is None else {'initial_step':step}
        with patch.object(calibration,'closure',SimpleNamespace(solve_ge=solve_ge)):
            result=calibration.solve_old_steady_state(chain,{},initial_psi=.14,
                completed_fertility_target=2.1,completed_fertility_tolerance=5e-4,normalize=True,**options)
        return result,calls

    def test_default_and_explicit_default_match_candidate_sequence_and_counts(self):
        implicit,a=self.solve(); explicit,b=self.solve(.25)
        self.assertEqual(a,b)
        self.assertEqual(a[:2],[.14,.14-.25])
        self.assertEqual(implicit[4],explicit[4])
        self.assertEqual(implicit[4]['stationary_solves'],len(a))
        self.assertEqual(implicit[4]['stationary_solve_seconds'],3.*len(a))

    def test_warm_step_is_explicit_and_meets_same_target(self):
        result,calls=self.solve(.05)
        self.assertEqual(calls[:3],[.14,.14-.05,.14-.1])
        self.assertLessEqual(result[4]['absolute_gap'],5e-4)
        self.assertEqual(result[4]['initial_psi'],.14)
        self.assertEqual(result[4]['initial_step'],.05)
        self.assertEqual(len(calls),len(set(calls)))

    def test_invalid_step_fails_before_ge(self):
        for value in [0.,-.1,float('nan'),float('inf'),25.]:
            with patch.object(calibration,'closure',SimpleNamespace(solve_ge=lambda *a:self.fail('no solve allowed'))):
                with self.assertRaises(ValueError):
                    calibration.solve_old_steady_state(None,{},initial_psi=.1,
                        completed_fertility_target=2.1,completed_fertility_tolerance=5e-4,
                        normalize=True,initial_step=value)

    def test_missed_tolerance_is_not_silently_accepted(self):
        calls=[]
        def solve_ge(chain, overrides):
            psi=overrides['psi_child']; calls.append(psi)
            return SimpleNamespace(tfr=2. if psi<.065 else 2.2),None,np.array([1.]),1.
        chain=SimpleNamespace(extract_moments=lambda s,p:{'tfr':s.tfr})
        with patch.object(calibration,'closure',SimpleNamespace(solve_ge=solve_ge)):
            with self.assertRaisesRegex(RuntimeError,'missed tolerance'):
                calibration.solve_old_steady_state(chain,{},initial_psi=.14,
                    completed_fertility_target=2.1,completed_fertility_tolerance=5e-4,
                    normalize=True,initial_step=.05)
        self.assertEqual(len(calls),len(set(calls)))


if __name__=='__main__': unittest.main()
