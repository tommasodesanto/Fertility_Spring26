#!/usr/bin/env python3
"""Torch-only targeted tests for the applied renter no-taper overlay; no solves."""
from __future__ import annotations

import ast
import hashlib
import importlib
import importlib.util
import inspect
import json
import os
import sys
import unittest
from pathlib import Path

import numpy as np

PACKET = Path(__file__).resolve().parent
ROOT = Path(os.environ["FROZEN_PROJECT_ROOT"]).resolve()
OVERLAY = Path(os.environ["RENTER_NO_TAPER_OVERLAY"]).resolve()
REL = Path("code/model/intergen_eqscale_seq_optimized")
FLAG = "renter_no_taper_estate_bound"
REFERENCE_LABEL = "2007 stationary reference — block0506, September 28 verified export"
MANIFEST = json.loads((ROOT / "output/model/fertility_identification_20260928/fixed_reference_manifest.json").read_text())


def load_parameters(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def load_overlay_solver():
    """Register the changed module, then import the real native solver path."""
    sys.path.insert(0, str(ROOT / "code/model"))
    importlib.import_module("intergen_eqscale_seq_optimized")
    name = "intergen_eqscale_seq_optimized.parameters"
    sys.modules.pop(name, None)
    overlay = load_parameters(OVERLAY / REL / "parameters.py", name)
    sys.modules[name] = overlay
    sys.modules.pop("intergen_eqscale_seq_optimized.solver", None)
    solver = importlib.import_module("intergen_eqscale_seq_optimized.solver")
    return overlay, solver


class RenterNoTaper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = load_parameters(ROOT / REL / "parameters.py", "frozen_parameters")
        cls.params, cls.solver = load_overlay_solver()

    def make(self, *, active: bool, survival: np.ndarray | None = None, lam: float = 0.0):
        P = self.params.setup_parameters()
        P.lambda_d = lam
        P.renter_no_taper_estate_bound = active
        if survival is not None:
            P.use_age_survival = True
            P.survival_probs = np.asarray(survival, dtype=float)
        self.params.build_debt_caps(P)
        return P

    def test_absent_flag_matches_frozen_builder(self):
        frozen = self.base.setup_parameters()
        changed = self.params.setup_parameters()
        self.base.build_debt_caps(frozen)
        self.params.build_debt_caps(changed)
        np.testing.assert_array_equal(changed.debt_taper_weights, frozen.debt_taper_weights)
        np.testing.assert_array_equal(changed.debt_caps, frozen.debt_caps)
        np.testing.assert_array_equal(changed.owner_ltv_multipliers, frozen.owner_ltv_multipliers)

    def test_default_off_identity_at_nonzero_lambda(self):
        frozen = self.base.setup_parameters()
        changed = self.params.setup_parameters()
        frozen.lambda_d = changed.lambda_d = 0.25
        self.base.build_debt_caps(frozen)
        self.params.build_debt_caps(changed)
        np.testing.assert_array_equal(changed.debt_taper_weights, frozen.debt_taper_weights)
        np.testing.assert_array_equal(changed.debt_caps, frozen.debt_caps)

    def test_all_ages_zero_and_negative_debt(self):
        P = self.make(active=True)
        for j in range(P.J):
            floor = self.solver.renter_borrowing_floor(P, np.array([-3.0, 0.0, 2.0]), j)
            expected = np.array([0.0, 0.0, 0.0]) if j == P.J - 1 else np.array([-3.0, 0.0, 0.0])
            np.testing.assert_array_equal(floor, expected)
            self.assertEqual(float(P.debt_caps[j + 1]), 0.0)

    def test_first_mortality_and_terminal_rule(self):
        template = self.params.setup_parameters()
        survival = np.ones(template.J - 1)
        survival[0] = 0.999
        P = self.make(active=True, survival=survival)
        self.assertEqual(float(P.debt_taper_weights[1]), 0.0)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, 0)), 0.0)
        self.assertEqual(float(P.debt_taper_weights[-1]), 0.0)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, P.J - 1)), 0.0)

    def test_later_first_mortality_preserves_prior_rollover(self):
        template = self.params.setup_parameters()
        survival = np.ones(template.J - 1)
        survival[3] = 0.999
        P = self.make(active=True, survival=survival)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, 2)), -2.0)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, 3)), 0.0)

    def test_serialized_reference_survival_and_native_mortality_timing(self):
        saved = MANIFEST["actual_serialized_parameters"]
        self.assertEqual(float(saved["lambda_d"]), 0.0)
        self.assertTrue(saved["use_age_survival"])
        survival = np.asarray(saved["survival_probs"], dtype=float)
        template = self.params.setup_parameters()
        self.assertEqual(template.J, int(saved["J"]))
        self.assertEqual(survival.shape, (template.J - 1,))
        P = self.make(active=True, survival=survival)
        first_death = int(np.flatnonzero(survival < 1.0)[0])
        self.assertEqual(first_death, 12)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, first_death - 1)), -2.0)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, first_death)), 0.0)
        self.assertEqual(float(self.solver.renter_borrowing_floor(P, -2.0, P.J - 1)), 0.0)

    def test_positive_new_credit_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "lambda_d == 0"):
            self.make(active=True, lam=1e-6)

    def test_flag_survives_native_override_rebuild(self):
        P = self.params.setup_parameters()
        P = self.params.apply_overrides(P, {FLAG: True})
        self.assertTrue(P.renter_no_taper_estate_bound)
        np.testing.assert_array_equal(P.debt_caps, np.zeros(P.J + 1))

    def test_native_buyer_and_incumbent_owner_floors_are_unchanged(self):
        for active in (False, True):
            P = self.make(active=active)
            P.native_purchase_income = True
            P.native_due_stayer_credit = True
            buyer_floor = self.solver.owner_borrowing_floor(P, 1.0, -1.0, 0)
            incumbent_floor = self.solver.owner_borrowing_floor(P, 1.0, -1.0, 0, stay_on=True)
            self.assertEqual(float(buyer_floor), -1.0)
            self.assertEqual(float(incumbent_floor), -1.0)

    def test_overlay_module_is_the_module_used_by_native_floor_path(self):
        imported = importlib.import_module("intergen_eqscale_seq_optimized.parameters")
        self.assertIs(imported, self.params)
        self.assertEqual(Path(inspect.getsourcefile(self.params)).resolve(), OVERLAY / REL / "parameters.py")
        self.assertEqual(self.solver.renter_borrowing_floor.__module__, "intergen_eqscale_seq_optimized.solver")

    def test_native_paths_continue_to_pass_the_same_scalars(self):
        source = (ROOT / REL / "solver.py").read_text()
        tree = ast.parse(source)
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and
                 isinstance(n.func, ast.Name) and n.func.id == "_savings_stage"]
        self.assertGreaterEqual(len(calls), 2)
        for call in calls:
            args = [ast.unparse(x) for x in call.args]
            self.assertIn("s_next", args)
            self.assertIn("D_next", args)
            self.assertIn("renter_floor", args)
        self.assertIn("rollover_floor = s_next *", (ROOT / REL / "kernels.py").read_text())


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RenterNoTaper)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    receipt = {
        "status": "PASS" if result.wasSuccessful() else "FAIL",
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "model_solves": 0,
        "standard_plots_changed": 0,
        "reference": REFERENCE_LABEL,
        "test_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (OVERLAY / "test_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    raise SystemExit(not result.wasSuccessful())
