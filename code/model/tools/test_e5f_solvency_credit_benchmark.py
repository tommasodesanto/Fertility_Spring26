import unittest
from types import SimpleNamespace
import numpy as np
import e5f_solvency_credit_benchmark as s

class SolvencyTests(unittest.TestCase):
    def test_default_noop(self):
        obj = SimpleNamespace(marker=object())
        before = vars(obj).copy()
        self.assertTrue(s.install(obj, '/nonexistent/no/write')['baseline_noop'])
        self.assertEqual(vars(obj), before)

    def test_terminal_net_liquidation(self):
        self.assertEqual(s.net_estate_value(-9.5, 10, .05), 0)
        self.assertLess(s.net_estate_value(-10, 10, .05), 0)

    def test_tiny_positive_risk_still_disqualifies(self):
        result = s.strict_expectation(np.array([[1., 2.], [s.DEAD, 3.]]), [1-1e-14, 1e-14])
        self.assertEqual(result[0], s.DEAD)
        self.assertGreater(result[1], 0)

    def test_zero_probability_does_not_disqualify(self):
        result = s.strict_expectation(np.array([[1.], [s.DEAD]]), [1., 0.])
        self.assertEqual(result[0], 1)

    def test_upper_interval_and_dead_column(self):
        bg = np.array([-3., -1., 0., 2.])
        v = np.array([[s.DEAD,s.DEAD], [s.DEAD,s.DEAD], [-2.,s.DEAD], [-1.,s.DEAD]])
        f, dead = s.natural_support_floor(v,bg)
        np.testing.assert_equal(f,[0.,2.]); np.testing.assert_equal(dead,[False,True])

    def test_hole_rejected(self):
        with self.assertRaises(ValueError):
            s.natural_support_floor(np.array([[1.],[s.DEAD],[2.]]),[-1.,0.,1.])

    def test_classifier_limitation_is_explicit(self):
        # Finite but extremely negative valid utility cannot be certified by
        # this draft: it deliberately inherits the native cutoff classification.
        f, dead = s.natural_support_floor(np.array([[-2e9],[-2.]]),[-1.,0.])
        self.assertEqual(f[0],0.)

    def test_independent_death_audit(self):
        r=s.audit_solvency_arrays(np.array([.4,.6]),[-10.,-9.],10.,[1e-14,1.],.05,[-20.,20.])
        self.assertEqual(r['negative_estate_exposure_mass'],.4)
        self.assertGreater(r['negative_estate_liability'],0)

    def test_survival_transition_length_and_terminal_death(self):
        P = SimpleNamespace(J=3, use_age_survival=True, survival_probs=np.array([1., .8]))
        np.testing.assert_allclose(s.death_probabilities(P), [0., .2, 1.])
        P.use_age_survival = False
        np.testing.assert_equal(s.death_probabilities(P), [0., 0., 1.])

    def test_bad_survival_support_rejected(self):
        for survival in ([1.], [1., 1., 0.], [1., np.nan], [1., 1.1]):
            with self.subTest(survival=survival), self.assertRaises(ValueError):
                s.death_probabilities(SimpleNamespace(J=3, use_age_survival=True, survival_probs=survival))

    def test_owner_wrapper_argument_positions(self):
        captured=[]
        def native(*args):
            captured.append(args)
            return (np.ones((3,1)),np.ones((3,1)),np.ones((3,1)))
        args=[0]*28
        args[2]=np.array([[s.DEAD],[1.],[2.]])
        args[5]=np.array([-2.,-1.,0.])
        s._kernel_wrapper(native,False)(*args)
        self.assertEqual(captured[0][12][0],-1.)
        self.assertEqual(captured[0][21],0)
        self.assertEqual(captured[0][22],0)

if __name__ == '__main__': unittest.main()
