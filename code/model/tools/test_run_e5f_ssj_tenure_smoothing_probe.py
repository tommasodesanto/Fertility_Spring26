"""Pure test of the tenure-smoothing probe summary; no native model."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "cluster"))
import run_e5f_ssj_tenure_smoothing_probe as probe  # noqa: E402


class TestProbeSummary(unittest.TestCase):
    def test_jumps_and_residual_blocks(self):
        T = 104
        own = np.linspace(0.57, 0.81, T); own[71:] += 0.017  # one 1.7pp jump between dates 70 and 71
        rows = [dict(owner_rate=float(own[t]), birth_children=0.05, adult_population=0.5, property_tax_revenue=0.04) for t in range(T)]
        residual = np.zeros(3 * T); residual[2 * T + 70] = -0.98; residual[5] = 1e-3
        s = probe.summarize_rows(rows, residual, T)
        self.assertEqual(s["owner_jump_dates_above_threshold"][0][0], 71)
        self.assertAlmostEqual(s["max_abs_owner_jump_pp"], 1.7 + 100 * (own[1] - own[0]), places=6)
        self.assertEqual(s["rebate_argmax_date"], 70)
        self.assertEqual(s["max_abs_housing"], 1e-3)
        self.assertEqual(s["rebate_spikes_above_0p15"], [(70, -0.98)])


if __name__ == "__main__":
    unittest.main()
