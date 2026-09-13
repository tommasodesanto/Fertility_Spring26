from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

import numpy as np

import run_e5f_rebated_initial_overnight as r


class RebateRootTests(unittest.TestCase):
    def parameters(self):
        return SimpleNamespace(tau_H=0.04, property_tax_lump_sum_transfer=0.0)

    def original(self, **kwargs):
        P = kwargs["parameters"]
        transfer = float(P.property_tax_lump_sum_transfer)
        mass = 4.0
        revenue = 2.0 - 0.2 * transfer
        outlays = mass * transfer
        sol = SimpleNamespace(
            property_tax_revenue=revenue,
            property_tax_transfer_outlays=outlays,
            property_tax_budget_residual=revenue - outlays,
            g=np.ones(4),
        )
        pension = {"actual_accounts": {"budget_residual": 0.0}}
        return sol, P, np.array([1.0]), pension

    def solve(self, original=None):
        return r.solve_equal_rebate(
            original or self.original, model=object(), parameters=self.parameters(),
            b_grid=np.array([0.0]), initial_prices=np.array([1.0]),
            payroll_tax=0.179, marginal_tolerance=1e-9, fiscal_tolerance=1e-6,
        )

    def test_roots_equal_rebate_and_preserves_pension_receipt(self):
        sol, P, _, receipt = self.solve()
        self.assertLessEqual(
            abs(sol.property_tax_budget_residual)
            / max(sol.property_tax_revenue, sol.property_tax_transfer_outlays),
            r.REBATE_RELATIVE_TOLERANCE,
        )
        self.assertAlmostEqual(sol.property_tax_transfer_outlays,
                               4.0 * P.property_tax_lump_sum_transfer)
        self.assertIn("actual_accounts", receipt)
        self.assertEqual(receipt["property_tax_rebate"]["fiscal_convention"],
                         "balanced_budget_equal_rebate")

    def test_rejects_changed_property_tax_primitive(self):
        P = self.parameters(); P.tau_H = 0.08
        with self.assertRaisesRegex(ValueError, "property-tax primitive"):
            r.solve_equal_rebate(
                self.original, model=object(), parameters=P,
                b_grid=np.array([0.0]), initial_prices=np.array([1.0]),
                payroll_tax=0.179, marginal_tolerance=1e-9,
                fiscal_tolerance=1e-6,
            )

    def test_unbracketed_rebate_is_hard_failure(self):
        def bad(**kwargs):
            P = kwargs["parameters"]
            sol = SimpleNamespace(property_tax_revenue=2.0,
                property_tax_transfer_outlays=float(P.property_tax_lump_sum_transfer),
                property_tax_budget_residual=1.0,
                g=np.ones(1))
            return sol, P, np.array([1.0]), {"actual_accounts": {}}
        with self.assertRaisesRegex(RuntimeError, "bracket"):
            self.solve(bad)


class ContractTests(unittest.TestCase):
    def objective(self):
        rows = [dict(restriction_id="initial_normalization", target=2.1,
                     role="normalization_separate_from_scored_objective")]
        rows.extend(dict(restriction_id=f"m{i}", target=0.0,
                         role="proposed_scored_restriction") for i in range(11))
        rows.append(dict(restriction_id="first_birth_rooms",
                         target=r.ROOMS_TARGET, role="proposed_scored_restriction"))
        restrictions = [dict(parameter=name, lower=0.0, upper=2.0)
                        for name in ("beta_annual", "kappa_fert",
                                     "kappa_fert_continuation", "chi", "H0",
                                     "theta0", "theta1", "first_birth_fixed_cost", "h_P")]
        restrictions[0].update(lower=0.94, upper=0.9995)
        return dict(target_rows=rows, parameter_restrictions=restrictions)

    def packet(self):
        objective = self.objective()
        candidate = {row["parameter"]: (0.99 if row["parameter"] == "beta_annual" else 1.0)
                     for row in objective["parameter_restrictions"]}
        return dict(objective=objective,
            initial=dict(structural_candidate=candidate,
                         fertility_normalization=2.1, payroll_tax=0.179,
                         normalize=True, observe_early=True),
            run=dict(source_root="/saved/source"),
            source_root=Path("/saved/source"))

    def test_full_objective_nine_coordinates_and_beta_cap(self):
        restrictions = r.validate_scientific_contract(self.packet())
        self.assertEqual(len(restrictions), 9)
        self.assertEqual(restrictions["beta_annual"]["upper"], 0.99)

    def test_target_change_is_rejected(self):
        packet = self.packet()
        next(row for row in packet["objective"]["target_rows"]
             if row["restriction_id"] == "first_birth_rooms")["target"] = 0.7
        with self.assertRaisesRegex(ValueError, "rooms target"):
            r.validate_scientific_contract(packet)

    def test_saved_packet_follows_resume_receipt_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            template = root / "template"; template.mkdir()
            controller = template / "run_capped_beta.py"; controller.write_text("x=1\n")
            case = root / "arbitrary_saved_case"; evaluation = case / "evaluation"
            scored = evaluation / "scored_repetition_01"; scored.mkdir(parents=True)
            seed = scored / "score.json"; seed.write_text("{}\n")
            initial = case / "initial_contract.json"; initial.write_text("{}\n")
            inputs = root / "explicit_inputs"; inputs.mkdir()
            scorer = inputs / "score_initial.py"; scorer.write_text("x=1\n")
            wrapper = inputs / "run_scored_candidate.py"; wrapper.write_text("x=1\n")
            objective = inputs / "working_contract.json"; objective.write_text('{"target_rows": []}\n')
            validator = inputs / "validator.py"; validator.write_text("x=1\n")
            run = case / "run_contract.json"
            run.write_text(json.dumps(dict(source_root=str(root / "saved_source"),
                scorer=dict(path=str(scorer)), validator=dict(path=str(validator)),
                working_objective=dict(path=str(objective)), wrapper_sha256=r.sha(wrapper))) + "\n")
            plan = dict(controller_sha256=r.sha(controller), resume_score_path=str(seed),
                        resume_score_sha256=r.sha(seed))
            (template / "plan_capped_beta_099.json").write_text(json.dumps(plan) + "\n")
            packet = r.saved_packet(template)
            self.assertEqual(packet["initial_path"], initial.resolve())
            self.assertEqual(packet["run_path"], run.resolve())
            self.assertNotIn("md_exact_loop", str(packet["run_path"]))


if __name__ == "__main__":
    unittest.main()
