#!/usr/bin/env python3
"""Pure controller tests; intentionally do not run the model."""
import math
import random
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_joint_nested_long_search as search

class ControllerTests(unittest.TestCase):
    def test_failure_classification(self):
        self.assertEqual(search.classify_failure("Housing market did not clear"), "market_nonconvergence")
        self.assertIsNone(search.classify_failure("mass accounting gate failed"))
        self.assertIsNone(search.classify_failure("fertility accounting tolerance failed"))
        self.assertEqual(search.classify_failure("Old-steady-state fertility normalization missed tolerance: x"), "fertility_normalization")
        self.assertEqual(search.classify_failure("x", "InfeasibleThetaError"), "infeasible_theta")

    def test_population_changes_all_coordinates_and_is_bounded(self):
        domain=[{"name":f"x{i}","lower":.1,"upper":10.,"transform":"log"} for i in range(11)]
        domain[1]["name"]="tenure_choice_kappa"; domain[2]["name"]="joint_nest_lambda"; domain[9]["name"]="hbar_first_child_jump"
        center=[.5]*11; pop=search.initial_population(center,domain,random.Random(20260906))
        self.assertEqual(len(pop),32)
        self.assertTrue(all(len(u)==11 and all(0<=x<=1 for x in u) for u in pop))
        self.assertTrue(all(any(x != .5 for x in u) for u in pop[1:]))

    def test_budget_and_valid_incumbent_logic(self):
        obj=object.__new__(search.Search); obj.completed=358; obj.c={"max_histories":360,"max_workers":12,"case_timeout_seconds":3600}; obj.finish=1e20; obj.search_finish=1e20
        self.assertTrue(obj.can_fit(2, final=True)); self.assertFalse(obj.can_fit(3, final=True)); self.assertFalse(obj.can_fit(2))
        valid=[{"loss":3.}, {"loss":2.}]
        rejected={"rejection_type":"market_nonconvergence", "error":"no calibrated loss"}
        self.assertEqual(search.best_completed(valid)["loss"],2.)
        self.assertNotIn("loss", rejected)  # A rejected proposal cannot become an incumbent.

if __name__ == "__main__": unittest.main()
