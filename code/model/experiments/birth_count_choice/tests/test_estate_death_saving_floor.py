"""No-solve checks of the experimental net-death-estate saving constraint."""
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.engine.household import (death_possible_at_age, net_estate_saving_floor,
                                    renter_borrowing_floor, owner_borrowing_floor)
from model.engine.kernels import full_owner_block_kernel, full_renter_block_kernel

PRICE = .7811670615311468
HOUSE = 2.
PSI = .06
GRID = np.array([-2., -1., 0., 1., 2.])
VALUES = np.zeros((len(GRID), 1))
ZERO = np.zeros(1)
ONE = np.ones(1)
GOLDEN = (.3819660112501051, .6180339887498949, 1e-5)


def parameters(net=True, survival=True, fixed_credit=0.):
    return SimpleNamespace(J=3, use_age_survival=survival,
        survival_probs=np.array([1., .9, 1.]), psi=PSI,
        estate_flow_net_of_selling_cost=net, estate_receiver='none',
        unsecured_credit_limit=fixed_credit, debt_taper_weights=np.ones(4),
        debt_caps=np.ones(4), native_purchase_income=True,
        owner_ltv_multipliers=np.ones(4))


def owner_policy(phi, estate_floor, *, stay=False, exhaustive=False):
    return full_owner_block_kernel(
        np.full(5, 3.), np.full(5, 3.), VALUES, np.zeros_like(VALUES), 0, GRID,
        ZERO, ZERO, ZERO, ZERO, np.full(1, .5), ONE,
        np.array([-phi * PRICE * HOUSE]), .1, HOUSE, 0., 1., 1e-10,
        .5, -1., .95, 0., 0., *GOLDEN,
        exhaustive_saving=int(exhaustive), due_stayer=stay,
        # The preexisting DUE stayer already had its own death floor.
        due_death_floor=-(1. - PSI) * PRICE * HOUSE if stay else -np.inf,
        net_estate_death_floor=estate_floor)[1]


def renter_policy(estate_floor, *, exhaustive=False):
    return full_renter_block_kernel(
        np.full(5, 3.), np.full(5, 3.), VALUES, np.zeros_like(VALUES), 0, GRID,
        ZERO, ZERO, ZERO, ZERO, np.full(1, .5), ONE,
        1., 5., 1e-10, 0., 0., .5, -1., .95, 1., 1., *GOLDEN,
        exhaustive_saving=int(exhaustive), net_estate_death_floor=estate_floor)[1]


class DeathEstateSavingFloorTests(unittest.TestCase):
    def test_mortality_and_terminal_only_when_net_estate_active(self):
        p = parameters()
        self.assertFalse(death_possible_at_age(p, 0))
        self.assertTrue(death_possible_at_age(p, 1))
        self.assertTrue(death_possible_at_age(p, 2))
        expected = -(1. - PSI) * PRICE * HOUSE
        self.assertEqual(net_estate_saving_floor(p, 1, PRICE, HOUSE), expected)
        self.assertEqual(net_estate_saving_floor(p, 2, PRICE, HOUSE), expected)
        self.assertEqual(net_estate_saving_floor(p, 1, PRICE, 0.), 0.)
        self.assertEqual(net_estate_saving_floor(p, 0, PRICE, HOUSE), -np.inf)
        p.estate_flow_net_of_selling_cost = False
        self.assertEqual(net_estate_saving_floor(p, 1, PRICE, HOUSE), -np.inf)
        p.use_age_survival = False
        self.assertTrue(death_possible_at_age(p, 2))
        self.assertFalse(death_possible_at_age(p, 1))

    def test_owner_buyer_stayer_both_search_paths_and_phi(self):
        floor = -(1. - PSI) * PRICE * HOUSE
        for exhaustive in (False, True):
            for stay in (False, True):
                for phi in (.8, 1.):
                    active = owner_policy(phi, floor, stay=stay, exhaustive=exhaustive)
                    legacy = owner_policy(phi, -np.inf, stay=stay, exhaustive=exhaustive)
                    self.assertGreaterEqual(float(active.min()), floor - 1e-12)
                    if phi == .8:
                        np.testing.assert_array_equal(active, legacy)
                    if phi == 1. and not stay:
                        self.assertLess(float(legacy.min()), floor - .01)
                        self.assertAlmostEqual(float(active.min()), floor, places=10)
        # The non-full path also maximizes its preexisting owner floor with the
        # same continuous death floor before optimizing, including changed tenure.
        p = parameters()
        for phi in (.8, 1.):
            collateral = -phi * PRICE * HOUSE
            old = owner_borrowing_floor(p, GRID, collateral, 1)
            effective = np.maximum(old, floor)
            self.assertTrue(np.all(effective >= floor))
            self.assertEqual(float(effective.min()), max(collateral, floor))

    def test_renter_full_and_fast_paths(self):
        p = parameters(fixed_credit=None)
        self.assertEqual(float(renter_borrowing_floor(p, np.array([-1.]), 1)[0]), 0.)
        self.assertEqual(float(renter_borrowing_floor(p, np.array([-1.]), 0)[0]), -1.)
        for exhaustive in (False, True):
            active = renter_policy(0., exhaustive=exhaustive)
            legacy = renter_policy(-np.inf, exhaustive=exhaustive)
            self.assertGreaterEqual(float(active.min()), -1e-12)
            self.assertLess(float(legacy.min()), -.1)
        p.estate_flow_net_of_selling_cost = False
        self.assertEqual(float(renter_borrowing_floor(p, np.array([-1.]), 1)[0]), -1.)
        p.unsecured_credit_limit = 1.
        self.assertEqual(float(renter_borrowing_floor(p, np.array([-1.]), 1)[0]), 0.)


if __name__ == '__main__':
    unittest.main()
