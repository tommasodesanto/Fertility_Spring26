from __future__ import annotations

import unittest

import numpy as np

from intergen_seq_fertility.parameters import apply_overrides, setup_parameters


class CombinedSpecificationTests(unittest.TestCase):
    def test_scalar_h0_override_is_coerced_to_one_market_vector(self) -> None:
        params = apply_overrides(setup_parameters(), {"H0": 4.25})
        np.testing.assert_array_equal(params.H0, np.array([4.25]))

    def test_finance_rates_are_four_year_objects(self) -> None:
        q = (1.0 + 0.02) ** 4 - 1.0
        delta = 1.0 - (1.0 - 0.011) ** 4
        params = apply_overrides(setup_parameters(), {"q": q, "delta": delta})
        self.assertAlmostEqual(params.q, 0.08243215999999998, places=14)
        self.assertAlmostEqual(params.delta, 0.04327930935900004, places=14)


if __name__ == "__main__":
    unittest.main()
