"""Small diagnostic-prefix checks; no Bellman or equilibrium evaluations.

Run with NUMBA_DISABLE_JIT=1 and the tools/model directories on PYTHONPATH.
"""
import copy
import unittest
from types import SimpleNamespace

import numpy as np

from e5f_historical_paygo_prefix import (
    YEARS, predict_historical_paygo_prefix, compare_historical_paygo_prefix,
)
from e5f_social_security import bind_social_security_income, fiscal_accounts


def fixture():
    P = SimpleNamespace(I=1, J=4, J_R=2, period_years=4.,
        scale_flows_to_period=True, tau_pay=.179, sequential_births=True,
        income_type_transition='markov', z_grid=np.array([.5, 1.5]),
        z_weights=np.array([2., 1.]), Pi_z=np.array([[.8, .2], [.1, .9]]),
        use_age_survival=True, survival_probs=np.array([.7, .8, .9]),
        w_hat=np.array([2.]), income_age_profile=np.array([1., 2., 1., 1.]),
        pension=1., retirement_income_z_scale=.5,
        n_house=1, n_parity=2, n_child_states=2, n_child_stages=1,
        child_state_mode='independent_count', use_numba_scatter=False)
    bind_social_security_income(P)
    initial = np.array([[1., 0.], [0., 1.], [.5, 0.], [0., .25]])
    ages = np.array([[1., 1., .5, .25], [2., 3., 1., 2.],
        [1., 4., 2., 3.], [3., 2., 4., 1.], [4., 3., 2., 1.]])
    return P, initial, ages


def predict(P, initial, ages, years=YEARS):
    return predict_historical_paygo_prefix(initial_age_income_mass=initial,
        parameters=P, observed_age_masses=ages, years=years)


def actual_ledger(P, marginal, benefit):
    dated = bind_social_security_income(copy.copy(P), pension_period=benefit)
    # Full seven-dimensional mass spread over wealth and tenure.
    shares = np.array([[.1, .2], [.3, .4]])
    g = shares[:, :, None, None, None, None, None] * marginal[None, None, None, :, :, None, None]
    return fiscal_accounts(g, dated)


class HistoricalPaygoPrefixTests(unittest.TestCase):
    def test_markov_orientation_entry_and_actual_small_transport_projection(self):
        from intergen_eqscale_seq_optimized import solver
        from intergen_eqscale_seq_optimized.parameters import make_independent_child_count_transition_matrix

        P, initial, ages = fixture()
        result = predict(P, initial, ages)
        M = result['age_income_marginals']
        np.testing.assert_allclose(M.sum(axis=2), ages, rtol=1e-14)
        np.testing.assert_allclose(M[1, 1], ages[1, 1] * np.array([.8, .2]))
        np.testing.assert_allclose(M[1, 2], ages[1, 2] * np.array([.1, .9]))
        np.testing.assert_allclose(M[1:, 0] / ages[1:, 0, None],
                                   np.tile([2./3., 1./3.], (4, 1)))

        # Compare one age of the actual transport, with saving, tenure and
        # child transitions active, then the same within-age reweighting.
        grid = np.array([0., 1.])
        shape = (2, 2, 1, 2, 2, 2)
        cohort = np.zeros(shape)
        weights = np.arange(1., 17.).reshape(2, 2, 1, 2, 2)
        weights /= weights.sum()
        for zz in range(2):
            cohort[:, :, :, zz, :, :] = initial[1, zz] * weights
        loc = np.ones((2, 2, 1, 1, 4, 2, 2, 2))
        choices = np.zeros((2, 2, 1, 4, 2, 2, 2), dtype=np.int16)
        tenure = np.empty(choices.shape + (2,))
        tenure[..., 0], tenure[..., 1] = .25, .75
        bp = np.broadcast_to(grid[:, None, None, None, None, None, None], choices.shape)
        idx, wt = solver.interp_indices(grid, grid)
        lidx, lwt = np.broadcast_to(idx, (1, 2, 2)), np.broadcast_to(wt, (1, 2, 2))
        tidx = np.broadcast_to(idx, (1, 2, 2, 2, 2, 2))
        twt = np.broadcast_to(wt, tidx.shape)
        child = make_independent_child_count_transition_matrix(2., 2)
        advanced = solver.advance_cohort_one_period_markov_income(
            P.survival_probs[1] * cohort, 1, loc, choices, tenure, bp, P,
            grid, SimpleNamespace(nc=4), lidx, lwt, tidx, twt, True, child, P.Pi_z)
        advanced *= ages[1, 2] / advanced.sum()
        np.testing.assert_allclose(advanced.sum(axis=(0, 1, 2, 4, 5)), M[1, 2], atol=1e-14)
        before = cohort.sum(axis=(0, 1, 2, 4, 5))
        values = np.ones_like(cohort)
        values[0] = -1e100
        self.assertGreater(solver._censor_entry_dead_mass(cohort, values), 0.)
        np.testing.assert_allclose(cohort.sum(axis=(0, 1, 2, 4, 5)), before)
        self.assertEqual(float(cohort[0].sum()), 0.)

    def test_period_units_retirement_boundary_and_date_scale_invariance(self):
        P, initial, ages = fixture()
        original_income = P.income.copy()
        result = predict(P, initial, ages)
        # Age index J_R is retired. Worker base is 4*(2*.5 + 2*2*1.5).
        # Retiree exposure is .5*.75 + .25*1.25; benefit is already four-year.
        expected = .179 * 28. / .6875
        self.assertAlmostEqual(result['pensions_period'][0], expected)
        self.assertAlmostEqual(result['predicted_ledgers'][0]['payroll_tax_base_period'], 28.)
        self.assertAlmostEqual(result['predicted_ledgers'][0]['retiree_benefit_exposure'], .6875)
        scales = np.array([7., .25, 3., 2., 10.])
        scaled = predict(P, initial * scales[0], ages * scales[:, None])
        np.testing.assert_allclose(scaled['pensions_period'], result['pensions_period'], rtol=1e-14)
        np.testing.assert_allclose(scaled['age_income_marginals'],
                                   result['age_income_marginals'] * scales[:, None, None])
        for i in range(5):
            self.assertAlmostEqual(scaled['predicted_ledgers'][i]['payroll_tax_revenue'],
                                   result['predicted_ledgers'][i]['payroll_tax_revenue'] * scales[i])
        np.testing.assert_array_equal(P.income, original_income)
        self.assertEqual(P.pension, 1.)
        np.testing.assert_array_equal(initial, fixture()[1])

    def test_actual_ledger_comparison_separates_prediction_from_trial_balance(self):
        P, initial, ages = fixture()
        result = predict(P, initial, ages)
        ledgers = [actual_ledger(P, M, 1.2 * benefit) for M, benefit in
                   zip(result['age_income_marginals'], result['pensions_period'])]
        receipt = compare_historical_paygo_prefix(prediction=result,
            actual_ledgers=ledgers, actual_years=YEARS)
        self.assertTrue(receipt['all_comparisons_pass'])
        self.assertFalse(receipt['actual_trial_budgets_pass'])
        self.assertFalse(receipt['historical_equilibrium_certified'])
        self.assertFalse(receipt['actual_marginals_verified'])
        changed = result['age_income_marginals'][0].copy()
        changed[0] = [0., 1.]
        ledgers[0] = actual_ledger(P, changed, result['pensions_period'][0])
        failed = compare_historical_paygo_prefix(prediction=result,
            actual_ledgers=ledgers, actual_years=YEARS)
        self.assertFalse(failed['all_comparisons_pass'])
        self.assertFalse(failed['rows'][0]['pension_prediction_pass'])
        self.assertFalse(failed['rows'][0]['ledger_factors_pass'])

    def test_invalid_inputs_and_looser_or_malformed_comparisons_fail(self):
        P, initial, ages = fixture()
        for attr, value in [('I', 2), ('tau_pay', .18), ('period_years', 1.),
            ('sequential_births', False), ('income_type_transition', 'fixed'),
            ('Pi_z', P.Pi_z.T), ('z_weights', np.array([0., 0.])),
            ('survival_probs', np.array([0., 1., 1.])), ('J_R', 4)]:
            with self.subTest(attr=attr), self.assertRaises(ValueError):
                changed = copy.copy(P)
                setattr(changed, attr, value)
                predict(changed, initial, ages)
        for invalid in [initial * np.nan, -initial, np.zeros_like(initial), initial * 2.]:
            with self.assertRaises(ValueError):
                predict(P, invalid, ages)
        for invalid in [ages * np.nan, -ages, np.zeros_like(ages), ages[:4]]:
            with self.assertRaises(ValueError):
                predict(P, initial, invalid)
        with self.assertRaises(ValueError):
            predict(P, initial, ages, years=(2007, 2011, 2015, 2019, 2027))
        result = predict(P, initial, ages)
        for tolerance in [1.01e-6, 0., np.nan]:
            with self.assertRaises(ValueError):
                compare_historical_paygo_prefix(prediction=result,
                    actual_ledgers=result['predicted_ledgers'], actual_years=YEARS,
                    fiscal_tolerance=tolerance)
        inconsistent = copy.deepcopy(result['predicted_ledgers'])
        inconsistent[0]['pension_outlays'] *= 2.
        with self.assertRaises(ValueError):
            compare_historical_paygo_prefix(prediction=result,
                actual_ledgers=inconsistent, actual_years=YEARS)


if __name__ == '__main__':
    unittest.main()
