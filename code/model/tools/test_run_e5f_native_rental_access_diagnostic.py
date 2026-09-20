import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_native_rental_access_diagnostic as d


class TestNativeRentalAccess(unittest.TestCase):
    def test_two_by_two_arm_contract(self):
        self.assertEqual(d.ARMS["capped_phi08"], ("capped", 0.8))
        self.assertEqual(d.ARMS["capped_phi1"], ("capped", 1.0))
        self.assertEqual(d.ARMS["open_phi08"], ("open", 0.8))
        self.assertEqual(d.ARMS["open_phi1"], ("open", 1.0))

    def test_open_arm_only_changes_cap_and_finance_controls(self):
        base = SimpleNamespace(phi=np.array([0.8]), H_own=np.array([4.0, 8.0]), hR_max=6.0)
        # No finance rebuild is needed for the phi=.8 arm in this structural check.
        arm = d.arm_parameters(base, "open", 0.8)
        self.assertEqual(arm.hR_max, 8.0)
        self.assertEqual(set(d.changed(base, arm)), {"hR_max"})

    def test_policy_schema_is_required(self):
        bad = SimpleNamespace(V=np.ones(1))
        with self.assertRaisesRegex(ValueError, "missing mandatory"):
            d.policy_arrays(bad)

    def test_exact_control_tolerance(self):
        p = SimpleNamespace(**{n: np.array([1.0]) for n in d.POLICY_NAMES})
        q = SimpleNamespace(**{n: np.array([1.0]) for n in d.POLICY_NAMES})
        d.compare_exact(p, q)
        q.price = np.array([1.0 + 2e-10])
        with self.assertRaises(AssertionError): d.compare_exact(p, q)

    def test_contract_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "checkpoint contract"):
            d.validate_contract(Path(__file__), Path(__file__).parent, Path(__file__).parent)

    def test_probability_and_finite_mass_gates(self):
        ev = SimpleNamespace(
            g_pre=np.array([1.0]), g_post_fertility=np.array([1.0]),
            g_current=np.array([1.0]), births=np.array([0.0]),
            policy=SimpleNamespace(**{n: np.array([1.0]) for n in d.POLICY_NAMES}),
        )
        ev.policy.fert_probs = np.array([1.01])
        with self.assertRaisesRegex(ValueError, "probability gate"):
            d.gates(ev)
        ev.policy.fert_probs = np.array([1.0])
        ev.g_current = np.array([np.nan])
        with self.assertRaisesRegex(ValueError, "nonfinite mass"):
            d.gates(ev)

    def test_synthetic_sequence_order_and_failure_stop(self):
        seen, progress, compared = [], [], []

        def runner(case):
            seen.append(case)
            if case == "rental_access": raise RuntimeError("synthetic failure")

        with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
            d.drive_case_order(runner, lambda i, case: progress.append((i, case)), lambda: compared.append(True))
        self.assertEqual(seen, ["baseline_control", "rental_access"])
        self.assertEqual(progress, [(1, "baseline_control")])
        self.assertEqual(compared, [])


if __name__ == "__main__": unittest.main()
