from __future__ import annotations

import math
import unittest
from types import SimpleNamespace

from e5f_post2023_tenure_sensitivity import apply_post2023_tenure_choice_kappa


class Post2023TenureSensitivityTest(unittest.TestCase):
    @staticmethod
    def prepared(value: float = 0.0) -> SimpleNamespace:
        return SimpleNamespace(
            base_overrides={"tenure_choice_kappa": value},
            parameters=SimpleNamespace(tenure_choice_kappa=value),
        )

    def test_default_is_noop(self) -> None:
        prepared = self.prepared()
        contract = apply_post2023_tenure_choice_kappa(prepared, None)
        self.assertFalse(contract["active"])
        self.assertEqual(prepared.base_overrides["tenure_choice_kappa"], 0.0)
        self.assertEqual(prepared.parameters.tenure_choice_kappa, 0.0)

    def test_positive_value_updates_both_parameter_views(self) -> None:
        prepared = self.prepared()
        contract = apply_post2023_tenure_choice_kappa(prepared, 0.010017488787185433)
        self.assertTrue(contract["active"])
        self.assertEqual(
            prepared.base_overrides["tenure_choice_kappa"], 0.010017488787185433
        )
        self.assertEqual(
            prepared.parameters.tenure_choice_kappa, 0.010017488787185433
        )

    def test_nonpositive_or_nonfinite_values_fail(self) -> None:
        for value in (0.0, -0.01, math.inf, math.nan):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    apply_post2023_tenure_choice_kappa(self.prepared(), value)


if __name__ == "__main__":
    unittest.main()
