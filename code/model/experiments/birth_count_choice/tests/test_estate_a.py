"""No-solve estate-A checks; production sources and economics stay untouched."""
from pathlib import Path
import importlib
import sys
import unittest
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parents[1]))
parameters = importlib.import_module('birth_count_choice.model.engine.parameters')
distribution = importlib.import_module('birth_count_choice.model.engine.distribution')
household = importlib.import_module('birth_count_choice.model.engine.household')
production_parameters = importlib.import_module('production.engine.parameters')
production_distribution = importlib.import_module('production.engine.distribution')


def toy_case(**overrides):
    # One terminal-age renter cell and owner cell; owner bp already nets debt.
    P = SimpleNamespace(J=1, J_R=1, I=1, z_grid=np.ones(1),
        period_years=4., income=np.array([[100.]]), tau_pay=.2,
        H_own=np.array([10.]), psi=.06, R=1.08243216, phi=.8,
        estate_receiver='none', estate_lump_sum_transfer=100.,
        age_start=50., da=4., theta0=1., theta1=1., theta_n=.2, sigma=2.)
    for name, value in overrides.items():
        setattr(P, name, value)
    mass = np.zeros((1, 2, 1, 1, 1, 1, 1))
    mass[0, 0, 0, 0, 0, 0, 0] = 2.
    mass[0, 1, 0, 0, 0, 0, 0] = 3.
    bp = np.zeros_like(mass)
    bp[0, 0, 0, 0, 0, 0, 0] = 7.
    bp[0, 1, 0, 0, 0, 0, 0] = -40.
    return P, mass, bp, np.array([40.]), np.array([10.])


def flow_stats(case, module=distribution):
    P, mass, bp, bg, ph = case
    stats = SimpleNamespace()
    module.add_aggregate_wealth_bequest_flow_moments(stats, mass, mass, bp, P, bg, ph)
    return stats


class EstateATests(unittest.TestCase):
    def test_default_off_matches_production_exactly(self):
        for receiver in ('none', 'ages_45_65'):
            for explicit_off in (False, True):
                case = toy_case(estate_receiver=receiver)
                if explicit_off:
                    case[0].estate_flow_net_of_selling_cost = False
                new, old = flow_stats(case), flow_stats(case, production_distribution)
                self.assertEqual(vars(new).keys(), vars(old).keys())
                for name in vars(old):
                    np.testing.assert_array_equal(getattr(new, name), getattr(old, name))
                for accounting in (False, True):
                    for price, rooms in ((10., 10.), (.73, 2.3), (0., 10.)):
                        self.assertEqual(
                            parameters.estate_housing_value(case[0], price, rooms, for_accounting=accounting),
                            production_parameters.estate_housing_value(case[0], price, rooms, for_accounting=accounting))

    def test_owner_100_debt_40_net_estate_54_without_interest(self):
        for interest in (1., 1.08243216, 2.):
            case = toy_case(R=interest, bequest_net_of_selling_cost=True,
                            estate_flow_net_of_selling_cost=True)
            P = case[0]
            utility_housing = parameters.estate_housing_value(P, 10., 10., for_accounting=False)
            flow_housing = parameters.estate_housing_value(P, 10., 10., for_accounting=True)
            self.assertEqual(utility_housing, 94.)
            self.assertEqual(-40. + utility_housing, 54.)
            self.assertEqual(utility_housing, flow_housing)
            self.assertEqual(flow_stats(case).annual_bequest_flow, (2. * 7. + 3. * 54.) / 4.)

    def test_actual_bequest_utility_uses_same_estate(self):
        P = toy_case(bequest_net_of_selling_cost=True,
                     estate_flow_net_of_selling_cost=True)[0]
        estate = np.array([7., -40. + parameters.estate_housing_value(
            P, 10., 10., for_accounting=False)])
        for children in (0, 1, 3):
            value = household.bequest_utility_vec(estate, children, P)
            np.testing.assert_array_equal(value, -(1. + .2 * children) / (1. + np.array([7., 54.])))

    def test_flow_and_utility_switches_are_independent(self):
        P = toy_case(bequest_net_of_selling_cost=True)[0]
        self.assertEqual(parameters.estate_housing_value(P, 10., 10., for_accounting=False), 94.)
        self.assertEqual(parameters.estate_housing_value(P, 10., 10., for_accounting=True), 100.)
        case = toy_case(estate_flow_net_of_selling_cost=True)
        self.assertEqual(parameters.estate_housing_value(case[0], 10., 10., for_accounting=False), 100.)
        self.assertEqual(flow_stats(case).annual_bequest_flow, 44.)

    def test_renter_estate_unchanged_and_receiver_stays_none(self):
        P = toy_case(bequest_net_of_selling_cost=True,
                     estate_flow_net_of_selling_cost=True)[0]
        self.assertFalse(parameters.estate_receiver_active(P))
        self.assertEqual(P.estate_receiver, 'none')
        for j in range(8):
            self.assertEqual(parameters.estate_transfer_at_age(P, j), 0.)
        for accounting in (False, True):
            self.assertEqual(parameters.estate_housing_value(P, 10., 0., for_accounting=accounting), 0.)

    def test_owner_stayer_uses_its_post_saving_wealth(self):
        case = toy_case(bequest_net_of_selling_cost=True,
                        estate_flow_net_of_selling_cost=True, native_due_stayer_credit=True)
        P, mass, bp, _, _ = case
        P._g_stay_distribution = np.zeros_like(mass)
        P._g_stay_distribution[0, 1, 0, 0, 0, 0, 0] = 1.
        P._bp_pol_stay = bp.copy()
        P._bp_pol_stay[0, 1, 0, 0, 0, 0, 0] = -20.
        # One of three owners stays: 2*54 + 1*74, plus two renter estates of 7.
        self.assertEqual(flow_stats(case).annual_bequest_flow, (2. * 7. + 2. * 54. + 74.) / 4.)
        P.estate_flow_net_of_selling_cost = False
        old = flow_stats(case, production_distribution)
        self.assertEqual(flow_stats(case).annual_bequest_flow, old.annual_bequest_flow)

    def test_current_financed_share_leaves_nonnegative_net_estate(self):
        P = toy_case(estate_flow_net_of_selling_cost=True)[0]
        for gross_house_value in (20., 40., 100.):
            owner_floor = -P.phi * gross_house_value
            net_housing = parameters.estate_housing_value(P, gross_house_value, 1., for_accounting=True)
            self.assertGreaterEqual(owner_floor + net_housing, 0.)
            self.assertAlmostEqual(owner_floor + net_housing, .14 * gross_house_value)


if __name__ == '__main__':
    unittest.main()
