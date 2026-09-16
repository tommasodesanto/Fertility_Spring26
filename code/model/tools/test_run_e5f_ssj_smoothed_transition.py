"""Pure tests for the smoothed-transition driver's fit table and scale override; no native model."""
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "cluster"))
import run_e5f_ssj_smoothed_transition as drv  # noqa: E402


class TestSmoothedTransitionPieces(unittest.TestCase):
    def test_fit_table_rows_and_markdown(self):
        a = dict(asset_price=0.70, owner_rate=0.566, adjusted_births=0.13, ignored="x")
        b = dict(asset_price=0.63, owner_rate=0.52, adjusted_births=0.13)
        rows, md = drv.fit_table(a, b, 0.005, 0.05)
        self.assertEqual([r["moment"] for r in rows], ["asset_price", "owner_rate", "adjusted_births"])
        self.assertAlmostEqual(rows[0]["relative_change"], -0.1)
        self.assertIn("| owner_rate | 0.566 | 0.52 | -8.13% |", md)

    def test_with_scale_copies_parameters_only(self):
        old = NS(parameters=NS(tenure_choice_kappa=0.005, psi_child=0.1), policy="P", initial_state="S")
        probe = drv.with_scale(old, 0.05)
        self.assertEqual(probe.parameters.tenure_choice_kappa, 0.05)
        self.assertEqual(old.parameters.tenure_choice_kappa, 0.005)
        self.assertIs(probe.policy, old.policy)


if __name__ == "__main__":
    unittest.main()
