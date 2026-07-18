from __future__ import annotations

import numpy as np
import unittest

from intergen_seq_fertility.parameters import apply_overrides, setup_parameters


class AgeSurvivalTests(unittest.TestCase):
    def test_age_survival_default_is_inert(self) -> None:
        P = setup_parameters()
        self.assertFalse(P.use_age_survival)
        np.testing.assert_array_equal(P.survival_probs, np.ones(P.J - 1))

    def test_age_survival_override_is_external_and_validated(self) -> None:
        P = setup_parameters()
        schedule = np.ones(P.J - 1)
        schedule[-4:] = [0.9391263063710125, 0.9184976343249724, 0.8849521927812863, 0.8300468061015381]
        P = apply_overrides(P, {"use_age_survival": True, "survival_probs": schedule})
        np.testing.assert_array_equal(P.survival_probs, schedule)

        with self.assertRaisesRegex(ValueError, "exactly J-1"):
            apply_overrides(setup_parameters(), {"survival_probs": np.ones(P.J)})
        with self.assertRaisesRegex(ValueError, r"\[0, 1\]"):
            apply_overrides(setup_parameters(), {"survival_probs": np.full(P.J - 1, 1.01)})


if __name__ == "__main__":
    unittest.main()
