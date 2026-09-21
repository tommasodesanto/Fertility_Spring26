import tempfile
import unittest
from pathlib import Path

from run_e5f_income_candidate_search import OBJECTIVE, fingerprint
from run_e5f_income_local_refinement import PARAMETERS, proposals, run_refinement


class LocalRefinementTest(unittest.TestCase):
    def setUp(self):
        self.params = dict(zip(PARAMETERS, [.99, .3377, .3978, 1.0496, 8.11,
                                            .081, .085, .265, 2.3]))
        rows = [{"restriction_id": f"moment_{i:02d}", "target": float(i + 1),
                 "model": float(i + 1), "gap": 0., "actual_weight": 1.,
                 "loss_contribution": 1., "scored": i < 12} for i in range(13)]
        self.plan = {"candidate_payload_fingerprint": "candidate-pin",
                     "parameter_bounds": {n: [0., 10.] for n in PARAMETERS},
                     "target_signature": fingerprint([(r["restriction_id"], r["target"],
                                                        r["actual_weight"], r["scored"]) for r in rows])}
        self.plan["parameter_bounds"]["beta_annual"] = [.94, .99]
        self.plan["parameter_bounds"]["h_P"] = [0., 2.3]
        self.rows = rows

    def receipt(self, loss, params=None, *, exact=True):
        params = params or self.params
        structural = [{"parameter": n, "estimate": params[n], "structural_coordinate": True}
                      for n in PARAMETERS]
        structural += [{"parameter": f"nuisance_{i}", "estimate": float(i),
                        "structural_coordinate": False} for i in range(8)]
        return {"schema": "e5f_initial_minimum_distance_result_v1", "contract_sha256": OBJECTIVE,
                "loss": loss, "objective": loss, "candidate_payload_fingerprint": "candidate-pin",
                "parameters": structural, "target_fit": self.rows,
                "normalization": {"target": 2.1, "absolute_gap": 0.},
                "_summary": {"status": "verified_scored_candidate", "objective_canonical_sha256": OBJECTIVE,
                             "loss": loss, "repetitions": 2 if exact else 1,
                             "exact_loss_equality": exact}}

    def test_design_is_twenty_unique_and_bounded(self):
        points = proposals(self.params, self.plan)
        self.assertEqual(len(points), 20)
        self.assertEqual(len({tuple(p[n] for n in PARAMETERS) for p in points}), 20)
        for point in points:
            self.assertTrue(.94 <= point["beta_annual"] <= .99)
            self.assertTrue(0. <= point["h_P"] <= 2.3)
            self.assertNotEqual(point, self.params)

    def test_smoke_gate_and_provisional_selection(self):
        calls = []
        anchor = self.receipt(326.9831988727637)

        def fake(point, case, **kwargs):
            calls.append(kwargs)
            case.mkdir(parents=True, exist_ok=True)
            if kwargs["repetitions"] == 2:
                return self.receipt(326.9831988727637, point)
            return self.receipt(400., point, exact=False)

        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "out"
            result = run_refinement(self.plan, Path(folder) / "plan.json", out,
                                    anchor, fake, anchor_psi=.23949950404168222,
                                    seconds=2., workers=10)
            self.assertIn(result["status"], ("provisional", "verified_selection"))
            self.assertEqual(calls[0]["repetitions"], 2)
            self.assertTrue((out / "design.json").exists())

    def test_bad_anchor_smoke_fails_closed(self):
        anchor = self.receipt(326.9831988727637)

        def fake(point, case, **kwargs):
            return self.receipt(326.99, point)

        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ValueError):
                run_refinement(self.plan, Path(folder) / "plan.json", Path(folder) / "out",
                               anchor, fake, anchor_psi=.2, seconds=2., workers=2)


if __name__ == "__main__":
    unittest.main()
