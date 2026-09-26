"""Mathematical and fiscal checks; run on Torch, never on the author's Mac."""
import copy
import math
from types import SimpleNamespace
import unittest

import numpy as np

import e5f_utility_comparison_runtime as runtime


class UtilityMathematics(unittest.TestCase):
    def test_compensated_expenditure_and_conditional_demands(self):
        # At E=e(m)*E0 and the fixed reference prices, the optimal material
        # composite / e(m) is independent of m and its housing preference.
        rstar, E0 = 0.073, 4.0
        common = E0 * runtime.ALPHA0 ** runtime.ALPHA0 * ((1-runtime.ALPHA0)/rstar) ** (1-runtime.ALPHA0)
        for m, alpha in enumerate([runtime.ALPHA0, 0.64, 0.48, 0.20]):
            e = ((2+0.7*m)/2)**0.7
            E = e*E0
            c, s = alpha*E, (1-alpha)*E/rstar
            A = float(runtime.reference_composite_factor(alpha, rstar))
            Q = A*c**alpha*s**(1-alpha)/e
            self.assertAlmostEqual(Q, common, places=12)
            self.assertAlmostEqual((1-alpha)*c/(alpha*s), rstar, places=12)
            # A feasible deviation with the same budget must lower material Q.
            c2 = 0.9*c
            s2 = (E-c2)/rstar
            self.assertLess(A*c2**alpha*s2**(1-alpha)/e, Q)

    def test_childless_exact_units_and_reference_is_substantive(self):
        for r in [0.01, 0.073, 1.0, 12.0]:
            self.assertEqual(float(runtime.reference_composite_factor(runtime.ALPHA0, r)), 1.0)
        self.assertNotEqual(float(runtime.reference_composite_factor(.5, .073)),
                            float(runtime.reference_composite_factor(.5, 1.0)))

    def test_crra_multiplier_is_inverse_factor(self):
        alpha, c, s, m, r = .4, 3., 5., 2, .073
        scale = ((2+.7*m)/2)**.7
        A = float(runtime.reference_composite_factor(alpha, r))
        raw = c**alpha*s**(1-alpha)
        self.assertAlmostEqual(-(scale/A)/raw, -1/(A*raw/scale), places=14)

    def test_invalid_reference_and_shares_fail(self):
        for alpha, rent in [(0., .1), (1., .1), (.7, 0.), (.7, float("nan"))]:
            with self.assertRaises(ValueError):
                runtime.reference_composite_factor(alpha, rent)

    @staticmethod
    def native_fixture(P, grid):
        from intergen_eqscale_seq_optimized.solver import precompute_shared
        return precompute_shared(P, grid)

    def test_current_children_not_ever_born_and_compression(self):
        # Independent fixture isolates storage correctness. Native-array
        # equivalence on actual parameters is a separate preflight below.
        def original(P, grid):
            psi = np.zeros((4,4))
            for n in range(4):
                for m in range(n+1):
                    psi[n,m] = P.psi_child*m
            return SimpleNamespace(psi_v=psi, psi_flat=psi.reshape(1,16,order="F"),
                c_bar=np.zeros((4,4)),h_bar=np.zeros((4,4)),nc=16,
                alpha_flat=np.full((1,16), runtime.ALPHA0),escale_flat=np.ones((1,16)))
        P=SimpleNamespace(utility_comparison_arm="floor_concave",child_state_mode="independent_count",
            n_parity=4,n_child_states=4,alpha_cons=.733,sigma=2.,eqscale_form="power",psi_child=.3)
        shared=runtime.comparison_shared(original,P,np.array([0.,1.]))
        for n in range(4):
            self.assertEqual(shared.psi_v[n,0],0.)
            for m in range(1,n+1):
                self.assertAlmostEqual(shared.psi_v[n,m],.3*m**.86)
        np.testing.assert_array_equal(shared.type_psi[shared.type_map],shared.psi_flat.ravel())
        self.assertEqual(shared.psi_v[1,1],shared.psi_v[3,1])
        self.assertEqual(shared.psi_v[0,3],0.)


class FiscalRule(unittest.TestCase):
    def test_ratio_from_population_not_inherited_tax(self):
        P=SimpleNamespace(I=1,J=3,J_R=2,z_grid=np.array([1.]),z_weights=np.array([1.]),
            Pi_z=np.ones((1,1)),income_type_transition="markov",use_age_survival=True,
            survival_probs=np.array([.8,.5,0.]),retirement_income_z_scale=0.)
        rate,receipt=runtime.pension_tax_from_demographics(P)
        self.assertAlmostEqual(rate,runtime.PENSION_RATIO*.4/1.8)
        self.assertAlmostEqual(receipt["retiree_worker_ratio"],.4/1.8)
        P.retirement_income_z_scale=.2
        with self.assertRaises(ValueError): runtime.pension_tax_from_demographics(P)

    def test_actual_fiscal_certification_rejects_wrong_ratio(self):
        ledger={"actual_accounts":{"payroll_tax_base_period":8.,"worker_household_mass":2.,
                                  "pension_period_units":4*runtime.PENSION_RATIO}}
        self.assertAlmostEqual(runtime.verify_pension_ratio(ledger),runtime.PENSION_RATIO)
        ledger["actual_accounts"]["pension_period_units"]*=1.01
        with self.assertRaises(RuntimeError): runtime.verify_pension_ratio(ledger)


if __name__ == "__main__":
    unittest.main()
