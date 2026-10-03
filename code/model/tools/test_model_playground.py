"""No-solve checks for the interactive model parameter binding."""
from __future__ import annotations

import copy
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from model_playground import ModelPlayground, PARAMETER_ORDER, bind_parameters


class ParameterBindingTests(unittest.TestCase):
    def setUp(self):
        self.values = dict(zip(PARAMETER_ORDER, (0.97, 1.1, 0.3, 0.12, 0.36, 0.1, 2.5, 0.1, 0.02, 0.18)))
        self.P = SimpleNamespace(
            period_years=4.0, beta=0.96 ** 4, rho=0.96 ** -4 - 1,
            rho_hat=0.96 ** -4 - 1, kappa_fert=0.2, eps_fert=0.2,
            child_room_floor=False, hbar_first_child_jump=0.0,
            hbar_child_rooms=0.0, chi=1.0,
        )

    def test_native_transform_and_primitive_mapping(self):
        original = copy.deepcopy(vars(self.P))
        Q = bind_parameters(self.P, self.values)
        self.assertEqual(Q.beta, self.values["beta_annual"] ** 4)
        self.assertEqual(Q.rho, 1 / Q.beta - 1)
        self.assertEqual(Q.rho_hat, Q.rho)
        self.assertEqual(Q.eps_fert, self.values["kappa_fert"])
        self.assertTrue(Q.child_room_floor)
        self.assertEqual(Q.hbar_first_child_jump, self.values["h_P"])
        self.assertEqual(Q.hbar_child_rooms, 0.0)
        for name in PARAMETER_ORDER:
            if name not in {"beta_annual", "h_P"}:
                self.assertEqual(getattr(Q, name), self.values[name])
        self.assertEqual(vars(self.P), original)

    def test_experiments_are_not_clamped_to_calibration_bounds(self):
        values = dict(self.values, chi=100.0)
        self.assertEqual(bind_parameters(self.P, values).chi, 100.0)

    def test_native_overrides_are_preserved_and_ten_parameters_take_precedence(self):
        edited = copy.deepcopy(self.P)
        edited.utility_reference_rent = 0.42
        edited.chi = 7.0
        Q = bind_parameters(edited, self.values)
        self.assertEqual(Q.utility_reference_rent, 0.42)
        self.assertEqual(Q.chi, self.values["chi"])

    def test_rejects_wrong_parameter_set_and_nonfinite_values(self):
        with self.assertRaises(ValueError):
            bind_parameters(self.P, {key: value for key, value in self.values.items() if key != "h_P"})
        with self.assertRaises(ValueError):
            bind_parameters(self.P, dict(self.values, chi=np.nan))


class ResultIsolationTests(unittest.TestCase):
    def test_mock_solve_keeps_result_and_live_parameters_independent(self):
        point = dict(zip(PARAMETER_ORDER, (0.97, 1.1, 0.3, 0.12, 0.36, 0.1, 2.5, 0.1, 0.02, 0.18)))
        P = SimpleNamespace(
            period_years=4.0, beta=0.97 ** 4, rho=1 / (0.97 ** 4) - 1,
            rho_hat=1 / (0.97 ** 4) - 1, kappa_fert=0.12, eps_fert=0.12,
            child_room_floor=True, hbar_first_child_jump=2.5,
            hbar_child_rooms=0.0, utility_reference_rent=0.33,
        )
        model = ModelPlayground.__new__(ModelPlayground)
        model.P = P
        model._authenticated_base_P = copy.deepcopy(P)
        model.params = dict(point)
        model._reference_params = dict(point)
        model.b_grid = np.array([-1.0, 0.0, 1.0])
        model.solver = object()
        model.reference_price = 0.7
        model.last_result = None
        candidates = []

        def fake_solve(candidate, *_args):
            candidates.append(candidate)
            return SimpleNamespace()

        with patch("model_playground.solve_at_price", side_effect=fake_solve):
            result = model.solve()
        self.assertIs(model.P, P)
        self.assertIsNot(result.P, P)
        self.assertEqual(candidates[0].chi, point["chi"])
        model.P.utility_reference_rent = 0.42
        model.params["chi"] = 1.3
        self.assertEqual(result.P.utility_reference_rent, 0.33)
        self.assertEqual(result.P.chi, point["chi"])
        model.reset()
        self.assertIs(model.P, P)
        self.assertEqual(model.P.utility_reference_rent, 0.33)
        self.assertEqual(model.params, point)


if __name__ == "__main__":
    unittest.main()
