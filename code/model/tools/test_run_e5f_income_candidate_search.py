import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from run_e5f_income_candidate_search import (
    ContractError, OBJECTIVE, ORIGINAL_PSI, _kill_process_group,
    run_search, proposals, validate_parameters,
)

PARAMETERS = ("beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
              "theta0", "theta1", "first_birth_fixed_cost", "h_P")
EXTRA_PARAMETERS = tuple(f"nuisance_{i:02d}" for i in range(8))
ALL_PARAMETER_ROWS = PARAMETERS + EXTRA_PARAMETERS


class BoundedSearchTest(unittest.TestCase):
    def setUp(self):
        from run_e5f_income_candidate_search import STOP
        STOP.clear()

    def plan(self):
        pilot = dict(zip(PARAMETERS, [.99, .3377, .3978, 1.0496, 8.11,
                                      .081, .085, .265, 2.3]))
        bounds = {name: [0., 10.] for name in PARAMETERS}
        bounds.update(beta_annual=[.94, .99], h_P=[0., 2.3])
        rows = self.receipt(1.)["target_fit"]
        return {"objective_canonical_sha256": OBJECTIVE,
                "candidate_payload_fingerprint": "candidate-pin",
                "candidate_json_sha256": "candidate-json-pin",
                "constructor_sha256": "constructor-pin",
                "parameter_bounds": bounds, "pilot_parameters": pilot,
                "target_signature": self.target_signature(rows)}

    @staticmethod
    def target_signature(rows):
        import run_e5f_income_candidate_search as search
        return search.fingerprint([(x["restriction_id"], x["target"],
                                    x.get("actual_weight"), x.get("scored"))
                                   for x in rows])

    def receipt(self, loss, *, parameters=None, timestamp="2026-09-19T00:00:00Z",
                provenance="fixture"):
        parameters = parameters or dict(zip(PARAMETERS,
                                             [.99, .3377, .3978, 1.0496, 8.11,
                                              .081, .085, .265, 2.3]))
        targets = []
        for i in range(13):
            scored = i < 12
            targets.append({"restriction_id": f"moment_{i:02d}", "target": float(i + 1),
                            "model": float(i + 1) + (loss / 1000 if scored else 0.),
                            "gap": loss / 1000 if scored else 0., "actual_weight": 1.0,
                            "loss_contribution": loss / 12 if scored else 0., "scored": scored})
        params = [{"parameter": name, "estimate": float(parameters.get(name, i + 1)),
                   "structural_coordinate": name in PARAMETERS}
                  for i, name in enumerate(ALL_PARAMETER_ROWS)]
        return {"schema": "e5f_initial_minimum_distance_result_v1",
                "contract_sha256": OBJECTIVE, "loss": float(loss), "objective": float(loss),
                "normalization": {"target": 2.1, "absolute_gap": 0.0,
                                   "completed_fertility": 2.1, "psi_child": .1},
                "candidate_payload_fingerprint": "candidate-pin", "parameters": params,
                "target_fit": targets,
                "_summary": {"status": "verified_scored_candidate",
                              "objective_canonical_sha256": OBJECTIVE,
                              "objective": float(loss), "loss": float(loss),
                              "repetitions": 2, "exact_loss_equality": True},
                "timestamp": timestamp, "provenance": provenance}

    @staticmethod
    def add_mock_diagnostics(evaluation, repetitions=(1, 2)):
        for repetition in repetitions:
            path = Path(evaluation) / "raw" / f"repetition_{repetition:02d}" / "standard_diagnostics"
            path.mkdir(parents=True, exist_ok=True)
            for i in range(17):
                (path / f"diagnostic_{i:02d}.png").write_bytes(b"PNG fixture")

    def test_order_has_all_first_directions_and_at_most_sixteen(self):
        plan = self.plan()
        values = proposals(plan["pilot_parameters"], plan)
        self.assertEqual(len(values), 16)
        self.assertEqual(values[0]["beta_annual"], .985)
        self.assertEqual(values[1]["kappa_fert"], .3377 * 1.2)
        self.assertNotEqual(values[8], plan["pilot_parameters"])

    def test_fake_evaluator_exercises_exact_sixteen_batch_cases_failure_best_and_repeat(self):
        plan, incumbent, calls = self.plan(), self.receipt(100.), []
        best_score = [100.]

        def fake(params, case, **kwargs):
            calls.append(kwargs)
            case.mkdir(parents=True, exist_ok=True)
            if kwargs["repetitions"] == 2:
                evaluation = case / "evaluation"
                self.add_mock_diagnostics(evaluation, (2,))
                return dict(self.receipt(best_score[0], parameters=params), _evaluation=str(evaluation))
            if abs(params["beta_annual"] - .985) < 1e-12:
                raise RuntimeError("numerical candidate failure")
            score = 90. - len(calls)
            best_score[0] = min(best_score[0], score)
            evaluation = case / "evaluation"
            self.add_mock_diagnostics(evaluation, (1,))
            return dict(self.receipt(score, parameters=params), _evaluation=str(evaluation))

        with tempfile.TemporaryDirectory() as folder:
            result = run_search(plan, Path(folder) / "search", fake, incumbent=incumbent,
                                search_seconds=2100, max_proposals=16, workers=4)
            self.assertEqual(result["proposal_count"], 16)
            self.assertTrue(any(row["status"] == "failed" for row in result["cases"]))
            self.assertEqual(result["verification"]["status"], "verified")
            self.assertEqual(sum(x["repetitions"] == 2 for x in calls), 1)
            self.assertEqual(len(list((Path(folder) / "search" /
                                       "selected_standard_diagnostics").glob("*.png"))), 17)

    def test_incumbent_wins_and_verifies_at_original_psi(self):
        plan, incumbent, calls = self.plan(), self.receipt(10.), []

        def fake(params, case, **kwargs):
            calls.append(kwargs)
            score = 10. if kwargs["repetitions"] == 2 else 20.
            return self.receipt(score, parameters=params)

        with tempfile.TemporaryDirectory() as folder:
            result = run_search(plan, Path(folder) / "search", fake, incumbent=incumbent,
                                search_seconds=2100, max_proposals=16, workers=4)
        self.assertEqual(result["selected"]["status"], "incumbent")
        self.assertEqual(result["verification"]["status"], "verified")
        self.assertEqual([x["psi"] for x in calls[-1:]], [ORIGINAL_PSI])

    def test_budget_stop_before_dispatch(self):
        plan, incumbent, calls = self.plan(), self.receipt(100.), []
        with tempfile.TemporaryDirectory() as folder:
            result = run_search(plan, Path(folder) / "search",
                                lambda *a, **k: calls.append(k), incumbent=incumbent,
                                search_seconds=499)
            self.assertEqual(len(calls), 1)
            self.assertEqual(result["proposal_count"], 0)
            self.assertEqual(json.loads((Path(folder) / "search/stop.json").read_text())["status"],
                             "budget_stop")

    def test_contract_mismatch_fails_closed(self):
        plan, incumbent, bad = self.plan(), self.receipt(100.), self.receipt(90.)
        bad["candidate_payload_fingerprint"] = "wrong-candidate"
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ContractError):
                run_search(plan, Path(folder) / "search", lambda *a, **k: bad,
                           incumbent=incumbent, search_seconds=2100, max_proposals=16, workers=4)

    def test_timestamp_and_provenance_are_ignored_by_numeric_comparison(self):
        plan, incumbent, calls = self.plan(), self.receipt(100.), []

        def fake(params, case, **kwargs):
            calls.append(kwargs)
            return self.receipt(50., parameters=params,
                                timestamp="later" if kwargs["repetitions"] == 2 else "earlier",
                                provenance="rerun" if kwargs["repetitions"] == 2 else "candidate")

        with tempfile.TemporaryDirectory() as folder:
            result = run_search(plan, Path(folder) / "search", fake, incumbent=incumbent,
                                search_seconds=2100, max_proposals=16, workers=4)
        self.assertEqual(result["verification"]["status"], "verified")
        self.assertTrue(result["verification"]["numeric_fit_equal"])

    def test_validate_parameters_requires_exact_structural_estimates(self):
        expected = dict(zip(PARAMETERS, [.99, .3377, .3978, 1.0496, 8.11,
                                          .081, .085, .265, 2.3]))
        receipt = self.receipt(1., parameters=expected)
        validate_parameters(receipt, expected)
        receipt["parameters"][0]["estimate"] += 1e-11
        with self.assertRaises(ContractError):
            validate_parameters(receipt, expected)

    def test_nested_process_cleanup_starts_new_session_and_kills_descendant(self):
        script = ("import subprocess,time,sys; "
                  "subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True); "
                  "time.sleep(30)")
        proc = subprocess.Popen([sys.executable, "-c", script], start_new_session=True)
        try:
            time.sleep(.15)
            self.assertEqual(os.getpgid(proc.pid), proc.pid)
            import psutil
            children = psutil.Process(proc.pid).children(recursive=True)
            self.assertTrue(children)
            _kill_process_group(proc)
            self.assertIsNotNone(proc.poll())
            for child in children:
                if psutil.pid_exists(child.pid):
                    self.assertIn(psutil.Process(child.pid).status(), ("zombie", "dead"))
        finally:
            if proc.poll() is None:
                _kill_process_group(proc)


if __name__ == "__main__":
    unittest.main()
