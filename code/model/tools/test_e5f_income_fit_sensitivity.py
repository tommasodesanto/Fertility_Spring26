import gzip
import json
import pickle
import tempfile
import unittest
import importlib.util
from pathlib import Path

import run_e5f_income_fit_sensitivity as sensitivity


class FitSensitivityTests(unittest.TestCase):
    def setUp(self):
        self.anchor = {
            "beta_annual": .99,
            "kappa_fert": .19777115505824105,
            "kappa_fert_continuation": .34649369327384283,
            "chi": .9530594713531557,
            "H0": 10.15538448099704,
            "theta0": .17078349399905368,
            "theta1": .07348722988916863,
            "first_birth_fixed_cost": .1673816654593897,
            "h_P": 2.3,
        }
        self.bounds = {
            "beta_annual": (.94, .99), "kappa_fert": (.02, 50.),
            "kappa_fert_continuation": (.02, 50.), "chi": (.1, 5.),
            "H0": (.2, 80.), "theta0": (0., 8.), "theta1": (.02, 16.),
            "first_birth_fixed_cost": (0., 8.), "h_P": (.1, 2.3),
        }

    def test_probe_count_and_one_sided_bounds(self):
        specs = sensitivity.build_probe_specs(self.anchor, self.bounds)
        self.assertEqual(len(specs), 24)
        self.assertEqual(sum(row["group"] == "full" for row in specs), 16)
        self.assertEqual(sum(row["group"] == "half" for row in specs), 8)
        self.assertNotIn("full_beta_annual_plus", {row["probe_id"] for row in specs})
        self.assertNotIn("full_h_P_plus", {row["probe_id"] for row in specs})
        for row in specs:
            value = row["parameters"][row["coordinate"]]
            lo, hi = self.bounds[row["coordinate"]]
            self.assertGreaterEqual(value, lo)
            self.assertLessEqual(value, hi)

    def test_rejection_whitelist_does_not_swallow_contract_failures(self):
        self.assertTrue(sensitivity.expected_rejection(RuntimeError("stationary root did not converge")))
        self.assertFalse(sensitivity.expected_rejection(RuntimeError("target fingerprint mismatch")))
        self.assertFalse(sensitivity.expected_rejection(RuntimeError("budget gate failed")))

    @unittest.skipUnless(importlib.util.find_spec("numpy"), "numpy is supplied by the Torch runtime")
    def test_native_array_signature_changes_with_array_bytes(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as directory:
            evaluation = Path(directory) / "evaluation"
            raw = evaluation / "raw/repetition_01"
            raw.mkdir(parents=True)
            (raw / "policy_array_summary.json").write_text(json.dumps({"maximum_occupied_value_drop": 0.0}))
            packet = {"b_grid": np.array([0., 1.]), "stationary_g_pre": np.array([.4, .6]),
                      "evaluation": {"policy": {"births": np.array([.1, .2])}}}
            with gzip.open(raw / "initial_state.pkl.gz", "wb") as stream:
                pickle.dump(packet, stream)
            first = sensitivity.artifact_signature(evaluation, require_files=True)
            packet["stationary_g_pre"][0] = .5
            with gzip.open(raw / "initial_state.pkl.gz", "wb") as stream:
                pickle.dump(packet, stream)
            second = sensitivity.artifact_signature(evaluation, require_files=True)
            self.assertNotEqual(first["raw/repetition_01/initial_state.pkl.gz"], second["raw/repetition_01/initial_state.pkl.gz"])


if __name__ == "__main__":
    unittest.main()
