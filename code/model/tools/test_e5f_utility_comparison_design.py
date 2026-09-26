"""Focused stdlib design checks; run on Torch, never importing model runtime.

    python3 code/model/tools/test_e5f_utility_comparison_design.py

All timing fixtures here are explicit hypothetical inputs, not adopted budgets.
"""
from __future__ import annotations

import copy
import json
import math
import unittest

import e5f_utility_comparison_design as design


def restrictions():
    return [
        {"parameter": "beta_annual", "lower": .94, "upper": .99, "transform": "discount"},
        {"parameter": "kappa_fert", "lower": .02, "upper": 50., "transform": "log"},
        {"parameter": "kappa_fert_continuation", "lower": .02, "upper": 50., "transform": "log"},
        {"parameter": "chi", "lower": .1, "upper": 5., "transform": "log"},
        {"parameter": "H0", "lower": .2, "upper": 80., "transform": "log"},
        {"parameter": "theta0", "lower": 0., "upper": 8., "transform": "softzero"},
        {"parameter": "first_birth_fixed_cost", "lower": 0., "upper": 8., "transform": "softzero"},
        {"parameter": "h_P", "lower": .1, "upper": 2.3, "transform": "log"},
    ]


def budget(**changes):
    arguments = dict(objective_solve_p90_seconds=100., objective_overhead_seconds=20.,
                     workers_per_arm=2, total_limit_seconds=10000.,
                     repeat_reserve_seconds=300., export_reserve_seconds=60.,
                     initial_population=8, de_generations=2, objective_timeout_seconds=300.)
    arguments.update(changes)
    return design.size_search_budget(**arguments)


def seeds():
    shared = {
        "a": dict(H0=7.5, beta_annual=.986, chi=1., first_birth_fixed_cost=.3,
                  kappa_fert=.29, kappa_fert_continuation=.45, theta0=.04),
        "b": dict(H0=6.3, beta_annual=.965, chi=2., first_birth_fixed_cost=0.,
                  kappa_fert=.08, kappa_fert_continuation=2., theta0=.6),
    }
    floor = {"a": {"h_P": 2.3}, "b": {"h_P": .1}}
    share = {"a": {"delta_alpha_jump": .020344, "delta_alpha": .016919},
             "b": {"delta_alpha_jump": 0., "delta_alpha": .25}}
    return shared, floor, share


def bank(**changes):
    shared, floor, share = seeds()
    arguments = dict(shared_seeds=shared, floor_seeds=floor, share_seeds=share,
                     broad_count=2, medium_count=2, local_count=2,
                     medium_unit_scale=.2, local_unit_scale=.04,
                     rng_seed=321, budget=budget())
    arguments.update(changes)
    return design.build_proposal_bank(design.define_arms(restrictions()), **arguments)


def target_rows():
    # Synthetic provenance has the same required structure as supplied frozen
    # rows.  This test makes no empirical claim about the placeholder targets.
    return [dict(id=name, target=2.1 if name == "initial_normalization" else 1.,
                 weight=None if name == "initial_normalization" else 2.,
                 sample="fixture sample", definition="fixture definition",
                 model_observation="fixture.observer", uncertainty_status="fixture only",
                 mapping_warning={"unmatched": True},
                 source=dict(path="fixture.json", builder="fixture.py",
                             record_id=name, contract_id="fixture"))
            for name in design.TARGET_IDS]


class ArmAndTransformTests(unittest.TestCase):
    def test_frozen_full_domain_and_sensitivity_contract(self):
        arms = design.define_arms(restrictions())
        self.assertEqual(tuple(arms), design.ARM_NAMES)
        for name, arm in arms.items():
            self.assertEqual(arm["child_reward_exponent"], .86 if name.endswith("concave") else 1.)
            self.assertEqual(arm["free_parameter_count"], 8 if name.startswith("floor") else 9)
            self.assertNotIn("theta1", arm["free_parameter_names"])
            self.assertNotIn("psi_child", arm["free_parameter_names"])
        changed = restrictions()
        changed[0]["upper"] = .9995
        with self.assertRaises(ValueError):
            design.define_arms(changed)
        with self.assertRaises(ValueError):
            design.define_arms(restrictions() + [restrictions()[0]])

    def test_transform_endpoints_interior_and_explicit_share_softzero(self):
        arms = design.define_arms(restrictions())
        rows = arms["shares_linear"]["parameter_restrictions"] + restrictions()
        for row in rows:
            for unit in (0., .01, .25, .5, .9, 1.):
                value = design.unit_to_physical(unit, row)
                self.assertGreaterEqual(value, row["lower"])
                self.assertLessEqual(value, row["upper"])
                self.assertAlmostEqual(design.physical_to_unit(value, row), unit, places=12)
            self.assertEqual(design.unit_to_physical(0., row), row["lower"])
            self.assertEqual(design.unit_to_physical(1., row), row["upper"])
        share_row = next(row for row in rows if row["parameter"] == "delta_alpha")
        self.assertEqual(design.unit_to_physical(.5, share_row), .0625)
        self.assertEqual(design.physical_to_unit(.0625, share_row), .5)
        beta = restrictions()[0]
        self.assertAlmostEqual(design.unit_to_physical(.5, beta), .9775)
        for invalid in (-.1, 1.1, math.nan, True):
            with self.assertRaises(ValueError):
                design.unit_to_physical(invalid, beta)
        with self.assertRaises(ValueError):
            design.physical_to_unit(.995, beta)


class BudgetTests(unittest.TestCase):
    def test_equal_counts_and_stationary_upper_bound(self):
        plan = budget(objective_solve_p90_seconds=1000., objective_overhead_seconds=100.,
                      workers_per_arm=10, total_limit_seconds=30000.,
                      repeat_reserve_seconds=2000., export_reserve_seconds=1000.,
                      initial_population=40, de_generations=4, objective_timeout_seconds=3100.)
        self.assertEqual(plan["search_waves"], 20)
        self.assertEqual(plan["total_objective_attempts"], 812)
        self.assertEqual(plan["total_attempts_per_arm"], 203)
        self.assertEqual(plan["search_population_evaluations_per_arm"], 200)
        self.assertEqual(plan["reused_smoke_seed_cases_per_arm"], 1)
        self.assertEqual(plan["stationary_calls_upper_bound"], 812 * 23)
        self.assertEqual(plan["planned_seconds"], 27200.)
        self.assertEqual(plan["smoke_estimate_seconds"], 2200.)
        self.assertIn("sequential within each arm", plan["smoke_schedule"])
        self.assertFalse(plan["timing_is_completion_guarantee"])
        design.validate_budget(json.loads(json.dumps(plan)))

    def test_generation_barrier_rounding_and_derived_maximum(self):
        plan = budget(initial_population=9, workers_per_arm=4, de_generations=2)
        self.assertEqual(plan["search_waves"], 9)
        self.assertNotEqual(plan["search_waves"], math.ceil(27 / 4))
        automatic = budget(total_limit_seconds=1800., repeat_reserve_seconds=240.,
                           export_reserve_seconds=120., de_generations=None)
        self.assertEqual(automatic["de_generations"], 1)
        self.assertLessEqual(automatic["planned_seconds"], 1800.)
        with self.assertRaises(ValueError):
            budget(total_limit_seconds=1800., repeat_reserve_seconds=240.,
                   export_reserve_seconds=120., de_generations=2)

    def test_invalid_or_forged_finite_budget_fails(self):
        for changes in ({"workers_per_arm": 0}, {"workers_per_arm": True},
                        {"de_generations": 1.5}, {"de_generations": -1},
                        {"repeat_reserve_seconds": 119.}, {"total_limit_seconds": 100.},
                        {"objective_overhead_seconds": -1.}, {"objective_solve_p90_seconds": math.nan},
                        {"objective_timeout_seconds": 110.}, {"initial_population": 3}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                budget(**changes)
        forged = budget()
        forged["total_objective_attempts"] += 1
        with self.assertRaises(ValueError):
            bank(budget=forged)


class ProposalTests(unittest.TestCase):
    def test_deterministic_common_draws_and_complete_coordinate_moves(self):
        first = bank()
        self.assertEqual(first, bank())
        shared, floor, share = seeds()
        for arm, points in first["arms"].items():
            names = design.FLOOR_NAMES if arm.startswith("floor") else design.SHARE_NAMES
            extra = floor if arm.startswith("floor") else share
            self.assertEqual(len(points), 8)
            self.assertEqual(points[0]["parameters"], {**shared["a"], **extra["a"]})
            self.assertEqual(points[1]["parameters"], {**shared["b"], **extra["b"]})
            for point in points:
                self.assertEqual(set(point["unit"]), set(names))
                self.assertTrue(all(0. <= value <= 1. for value in point["unit"].values()))
                reference = first["arms"]["floor_linear"][point["slot"]]
                for name in design.COMMON_NAMES:
                    self.assertEqual(point["unit"][name], reference["unit"][name])
                    self.assertEqual(point["parameters"][name], reference["parameters"][name])
                if point["stage"] in {"medium", "local"}:
                    center = {**shared[point["center_seed"]], **extra[point["center_seed"]]}
                    self.assertTrue(all(point["parameters"][name] != center[name] for name in names))
        self.assertNotEqual(first["arms"], bank(rng_seed=322)["arms"])

    def test_seed_and_count_contract_errors(self):
        shared, floor, share = seeds()
        for changes in ({"broad_count": 0}, {"local_count": 3},
                        {"local_unit_scale": .2}, {"medium_unit_scale": 1.1},
                        {"share_seeds": {"a": share["a"]}}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                bank(**changes)
        shared["a"]["theta1"] = .01
        with self.assertRaises(ValueError):
            bank(shared_seeds=shared)
        _, floor, _ = seeds()
        floor["a"]["h_P"] = 2.31
        with self.assertRaises(ValueError):
            bank(floor_seeds=floor)


class DifferentialEvolutionTests(unittest.TestCase):
    def setUp(self):
        self.arms = design.define_arms(restrictions())
        self.initial = bank()["arms"]
        self.plan = budget()

    def generation(self, arm="floor_linear", population=None, scores=None, **changes):
        population = self.initial[arm] if population is None else population
        scores = {row["id"]: 100. for row in population} if scores is None else scores
        arguments = dict(generation=1, budget=self.plan, rng_seed=42, mutation_factor=.7)
        arguments.update(changes)
        return design.make_de_generation(self.arms, arm, population, scores, **arguments)

    def test_complete_finite_trials_and_common_donor_streams(self):
        reference = self.generation()
        self.assertEqual(reference, self.generation())
        for arm in design.ARM_NAMES:
            trials = self.generation(arm)
            self.assertEqual(len(trials), 8)
            for index, trial in enumerate(trials):
                self.assertEqual(trial["donor_slots"], reference[index]["donor_slots"])
                self.assertNotIn(index, trial["donor_slots"])
                self.assertEqual(len(set(trial["donor_slots"])), 3)
                self.assertEqual(trial["crossover_probability"], 1.)
                for name in design.COMMON_NAMES:
                    self.assertEqual(trial["unit"][name], reference[index]["unit"][name])

    def test_completed_failures_consume_slots_and_selection_advances_barrier(self):
        parents = self.initial["floor_linear"]
        parent_scores = {row["id"]: 100. for row in parents}
        parent_scores[parents[2]["id"]] = None
        trials = self.generation(scores=parent_scores)
        trial_scores = {row["id"]: 100. for row in trials}
        trial_scores.update({trials[0]["id"]: 50., trials[1]["id"]: None, trials[2]["id"]: 30.})
        survivors, scores = design.select_de_generation(parents, parent_scores, trials, trial_scores)
        self.assertEqual([row["id"] for row in survivors[:4]],
                         [trials[0]["id"], parents[1]["id"], trials[2]["id"], parents[3]["id"]])
        self.assertTrue(all(row["selection_generation"] == 1 for row in survivors))
        self.assertTrue(all(row["selection_generation"] == 0 for row in parents))
        self.assertEqual(len(self.generation(population=survivors, scores=scores, generation=2)), 8)
        with self.assertRaises(ValueError):
            self.generation(population=survivors, scores=scores, generation=3)

    def test_missing_receipts_bad_pairs_and_skipped_barriers_stop(self):
        parents = self.initial["floor_linear"]
        scores = {row["id"]: 100. for row in parents}
        missing = dict(scores)
        missing.pop(parents[0]["id"])
        with self.assertRaises(ValueError):
            self.generation(scores=missing)
        with self.assertRaises(ValueError):
            self.generation(scores={key: None for key in scores})
        with self.assertRaises(ValueError):
            self.generation(generation=2)
        altered = copy.deepcopy(parents)
        altered[0]["parameters"]["H0"] += 1.
        with self.assertRaises(ValueError):
            self.generation(population=altered)
        trials = self.generation()
        trial_scores = {row["id"]: 99. for row in trials}
        trials[0]["parent_id"] = parents[1]["id"]
        with self.assertRaises(ValueError):
            design.select_de_generation(parents, scores, trials, trial_scores)


class TargetContractTests(unittest.TestCase):
    def test_complete_rows_and_source_metadata_preserved(self):
        source = target_rows()
        receipt = design.validate_target_contract(source)
        self.assertEqual(receipt["target_rows"], source)
        self.assertEqual(receipt["positive_weight_rows"], 12)
        self.assertEqual(receipt["normalization"]["target"], 2.1)
        self.assertEqual(receipt["free_parameter_counts"], {"floor": 8, "shares": 9})
        source[1]["weight"] = 3.
        self.assertNotEqual(design.validate_target_contract(source)["target_rows_sha256"],
                            receipt["target_rows_sha256"])
        self.assertEqual(receipt["target_rows"][1]["weight"], 2.)

    def test_missing_or_demoted_target_and_source_fail(self):
        for mode in ("missing", "duplicate", "zero_weight", "normalization", "source", "sample"):
            rows = target_rows()
            if mode == "missing":
                rows.pop()
            elif mode == "duplicate":
                rows[-1] = copy.deepcopy(rows[1])
            elif mode == "zero_weight":
                rows[1]["weight"] = 0.
            elif mode == "normalization":
                rows[0]["weight"] = 1.
            elif mode == "source":
                rows[1]["source"].pop("builder")
            else:
                rows[1]["sample"] = ""
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                design.validate_target_contract(rows)


if __name__ == "__main__":
    unittest.main()
