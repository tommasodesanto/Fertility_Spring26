"""Torch-only synthetic tests of actual AHS rooms versus retained capped rows.

Run in a Torch allocation from code/model/tools:
    python -m unittest test_e5f_ahs_rooms_observer
No household model import, solve, target or weight modification is needed.
"""
import os
import sys

if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Run AHS observer tests only in a Torch allocation")

import json
from types import SimpleNamespace
import unittest

import numpy as np

from e5f_initial_housing_observer import (
    AGE_PROJECTION, AHS_MEAN_ROOMS, MOMENT_NAMES, observe_initial_housing_wealth,
)


def fixture():
    P = SimpleNamespace(
        J=17, da=4., period_years=4., age_start=18.,
        n_house=1, I=1, n_parity=4, n_child_states=4,
        H_own=np.array([11.]), hR_max=6., z_grid=np.array([.5, 1.5]),
        child_state_mode="independent_count", child_bin_high_cutoff=3,
        chi=8.,
    )
    g = np.zeros((2, 2, 1, 17, 2, 4, 4))
    h_r = np.full_like(g, np.nan)
    ev = SimpleNamespace(g_current=g, policy=SimpleNamespace(hR_pol=h_r))
    return P, ev


def observe(P, ev, **kwargs):
    return observe_initial_housing_wealth(
        ev, P, None, None, diagnostic_enabled=True, age_projection=AGE_PROJECTION,
        include_wealth=False, **kwargs)


def row(result, name):
    return next(item for item in result["rows"] if item["moment"] == name)


class AhsRoomsObserverTests(unittest.TestCase):
    def test_disabled_stays_disabled_and_old_moment_positions_survive(self):
        result = observe_initial_housing_wealth(None, None, None, None)
        self.assertIsNone(result["moments"][AHS_MEAN_ROOMS])
        self.assertEqual(MOMENT_NAMES[0], "aggregate_mean_occupied_rooms_capped9_18_85")
        self.assertEqual(MOMENT_NAMES[3], "prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9")
        self.assertEqual(MOMENT_NAMES[9], "housing_increment_0to1")

    def test_realized_tenure_counts_owners_as_rooms_not_renter_policy_or_premium(self):
        P, ev = fixture()
        # Same underlying state: 25% rent four rooms, 75% own eleven rooms.
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = .25
        ev.g_current[0, 1, 0, 3, 0, 0, 0] = .75
        ev.policy.hR_pol[0, 0, 0, 3, 0, 0, 0] = 4.
        ev.policy.hR_pol[0, 1, 0, 3, 0, 0, 0] = 999.
        before_g, before_h = ev.g_current.copy(), ev.policy.hR_pol.copy()
        result = observe(P, ev)
        self.assertEqual(result["moments"][AHS_MEAN_ROOMS], 9.25)
        self.assertEqual(result["moments"][MOMENT_NAMES[0]], 7.75)
        self.assertEqual(result["moments"]["own_rate_30_55"], .75)
        detail = row(result, AHS_MEAN_ROOMS)
        self.assertTrue(detail["configured_room_support_within_ahs_topcode"])
        self.assertEqual(detail["realized_mass_above_ahs_topcode"], 0.)
        self.assertIsNone(detail["rooms_capped_at"])
        self.assertFalse(detail["model_topcode_applied"])
        self.assertFalse(detail["production_eligible"])
        np.testing.assert_array_equal(ev.g_current, before_g)
        np.testing.assert_array_equal(ev.policy.hR_pol, before_h)
        json.dumps(result, allow_nan=False)

    def test_uncapped_renters_preserve_income_state_weights_and_legacy_cap(self):
        P, ev = fixture()
        P.hR_max = 12.
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = .25
        ev.g_current[0, 0, 0, 3, 1, 0, 0] = .75
        ev.policy.hR_pol[0, 0, 0, 3, 0, 0, 0] = 4.
        ev.policy.hR_pol[0, 0, 0, 3, 1, 0, 0] = 12.
        result = observe(P, ev)
        self.assertEqual(result["moments"][AHS_MEAN_ROOMS], 10.)
        self.assertEqual(result["moments"][MOMENT_NAMES[0]], 7.75)

    def test_family_gap_keeps_cap9_and_current_dependent_groups(self):
        P, ev = fixture()
        P.hR_max = 12.
        ev.g_current[0, 0, 0, 3, 0, 1, 1] = 1.
        ev.policy.hR_pol[0, 0, 0, 3, 0, 1, 1] = 12.
        ev.g_current[0, 1, 0, 3, 0, 3, 3] = 1.
        result = observe(P, ev, diagnostic_allow_family_proxies=True)
        self.assertEqual(result["moments"][AHS_MEAN_ROOMS], 11.5)
        self.assertEqual(result["moments"][MOMENT_NAMES[0]], 9.)
        self.assertEqual(result["moments"][MOMENT_NAMES[3]], 0.)
        self.assertEqual(row(result, MOMENT_NAMES[3])["rooms_capped_at"], 9.)

    def test_all_ages_keep_household_mass_weights_including_endpoint_cells(self):
        P, ev = fixture()
        ev.g_current[0, 0, 0, 0, 0, 0, 0] = 2.
        ev.policy.hR_pol[0, 0, 0, 0, 0, 0, 0] = 3.
        ev.g_current[0, 1, 0, 16, 0, 0, 0] = 1.
        result = observe(P, ev)
        self.assertAlmostEqual(result["moments"][AHS_MEAN_ROOMS], 17. / 3.)
        self.assertEqual(row(result, AHS_MEAN_ROOMS)["denominator"], 3.)
        self.assertEqual(result["age_overlap_weights"]["18_85"], [1.] * 17)

    def test_future_support_above21_is_reported_and_never_silently_topcoded(self):
        P, ev = fixture()
        P.H_own = np.array([24.])
        ev.g_current[0, 1, 0, 3, 0, 0, 0] = 1.
        result = observe(P, ev)
        self.assertEqual(result["moments"][AHS_MEAN_ROOMS], 24.)
        self.assertEqual(result["moments"][MOMENT_NAMES[0]], 9.)
        detail = row(result, AHS_MEAN_ROOMS)
        self.assertFalse(detail["configured_room_support_within_ahs_topcode"])
        self.assertEqual(detail["realized_mass_above_ahs_topcode"], 1.)

    def test_zero_mass_is_unavailable_and_occupied_invalid_renter_is_rejected(self):
        P, ev = fixture()
        result = observe(P, ev)
        self.assertIsNone(result["moments"][AHS_MEAN_ROOMS])
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = 1.
        with self.assertRaisesRegex(ValueError, "Occupied renters"):
            observe(P, ev)


if __name__ == "__main__":
    unittest.main()
