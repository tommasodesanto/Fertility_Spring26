"""Zero-model-solve contract checks for the isolated financing path."""
from __future__ import annotations

import ast
import difflib
import hashlib
import json
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "code/model/tools"))
import dated_phi
from integration import bind_phi_path, financing_path
import run_case
from selected_runtime import authenticate_selected


class PhiPathTests(unittest.TestCase):
    def test_fresh_startup_defers_legacy_import_until_after_authentication_point(self):
        script = (
            f"import sys; sys.path[:0]={[str(Path(__file__).parent), str(ROOT / 'code/model/tools'), str(ROOT / 'code/model/experiments/transition_readiness/pinned_tools')]!r}; "
            "import run_case; "
            "assert run_case.map_case is None; "
            "assert not any(n.startswith('intergen_eqscale_seq_optimized') for n in sys.modules); "
            "import integration; "
            "assert any(n.startswith('intergen_eqscale_seq_optimized') for n in sys.modules)"
        )
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_authenticated_compact_saved_repeat(self):
        with tempfile.TemporaryDirectory() as tmp:
            postcheck = Path(tmp) / "chain_9/postcheck"
            root = postcheck / "selected_postcheck/phase_b_ge/selected_root"
            repeat = postcheck / "selected_postcheck/phase_b_ge/selected_repeat"
            root.mkdir(parents=True)
            (repeat / "stage").mkdir(parents=True)
            (root / "closure.json").write_text(json.dumps({"price": .7, "population_scale": 1.,
                                                             "report_only": 17}))
            compact = repeat / "closure.json"
            compact.write_text(json.dumps({"price": .7, "population_scale": 1.}))
            arrays = repeat / "stage/solution_arrays.npz"
            arrays.write_bytes(b"native-checkpoint-fixture")
            (postcheck / "input_contract.json").write_text(json.dumps(dict(
                purchase_rule="hard", normalized_population=1., owner_financed_share=.8,
                entry=dict(arm="nonnegative_mean"), weight_contract_sha256="weight",
                parameter_contract_sha256="parameters")))
            completed = postcheck / "completed.json"
            completed.write_text(json.dumps(dict(status="selected_numerically_verified",
                selected=dict(parameters={"p": .1}, weight_fingerprint="weight"),
                selected_postcheck=dict(status="passed", parameters={"p": .1},
                                        weight_fingerprint="weight"))))
            self.assertEqual(authenticate_selected("hard", completed)[2:],
                             (root.resolve(), repeat.resolve()))
            arrays.unlink()
            with self.assertRaisesRegex(RuntimeError, "state arrays missing"):
                authenticate_selected("hard", completed)
            arrays.write_bytes(b"native-checkpoint-fixture")
            compact.write_text(json.dumps({"price": .8, "population_scale": 1.}))
            with self.assertRaisesRegex(RuntimeError, "repeat closure differs"):
                authenticate_selected("hard", completed)

    def test_paths_and_terminal_values(self):
        self.assertEqual(financing_path("control", 3).tolist(), [.8,.8,.8])
        self.assertEqual(financing_path("temporary", 3).tolist(), [1.,.8,.8])
        self.assertEqual(financing_path("permanent", 3).tolist(), [1.,1.,1.])
        with self.assertRaises(ValueError): financing_path("temporary", 0)

    def test_backward_uses_date_specific_phi_and_exact_policy_function(self):
        calls = []
        old_model = dated_phi.calendar.model
        fake_model = SimpleNamespace(precompute_shared=lambda P, grid: None)
        def fake_policy(**kw):
            calls.append((float(kw["P"].phi[0]), float(kw["P"].psi_child)))
            return SimpleNamespace(V=np.array([len(calls)],float))
        try:
            dated_phi.calendar.model = fake_model
            with patch.object(dated_phi.base, "solve_date_policy", fake_policy):
                values, count = dated_phi.backward_value_path(
                    prices=np.ones(3), rents=np.ones(3), psi_path=np.full(3,.2),
                    phi_path=np.array([1.,.8,.8]), terminal_V=np.array([0.]),
                    base_parameters=SimpleNamespace(phi=np.full(4,.8),psi_child=.2),
                    b_grid=np.array([0.]))
        finally:
            dated_phi.calendar.model = old_model
        self.assertEqual(count, 3)
        self.assertEqual(calls, [(.8,.2),(.8,.2),(1.,.2)])
        self.assertEqual(values[0].tolist(),[3.])

    def test_binding_restores_unmodified_engine(self):
        sentinel = lambda **kwargs: kwargs
        pf = SimpleNamespace(evaluate_path_at_prices=sentinel)
        with patch.object(dated_phi,"base",pf):
            with bind_phi_path(pf, [1.,.8]):
                self.assertIsNot(pf.evaluate_path_at_prices, sentinel)
                with self.assertRaises(RuntimeError): pf.evaluate_path_at_prices(prices=[1.])
        self.assertIs(pf.evaluate_path_at_prices, sentinel)

    def test_dated_copy_has_only_declared_model_changes(self):
        original = (ROOT / "code/model/tools/run_e5f_perfect_foresight_transition.py").read_text()
        copied = (Path(__file__).parent / "dated_phi.py").read_text()
        source_ast, copy_ast = ast.parse(original), ast.parse(copied)
        expected = {
            "backward_value_path": "021c5d8cbdee0c1b51b071c77799846a848030fb67defb6789f51edd0b0b6f97",
            "evaluate_path_at_prices": "96c779d7645477070ff8f08cf99db8d9a072e99a217ad7ed0fba05e733723259",
        }
        for name in expected:
            a = next(x for x in source_ast.body if isinstance(x,ast.FunctionDef) and x.name==name)
            b = next(x for x in copy_ast.body if isinstance(x,ast.FunctionDef) and x.name==name)
            changed = "\n".join(difflib.unified_diff(
                ast.get_source_segment(original,a).splitlines(),
                ast.get_source_segment(copied,b).splitlines(),n=0))
            self.assertEqual(hashlib.sha256(changed.encode()).hexdigest(), expected[name])

    def test_terminal_and_dated_root_loops_without_model_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            P = SimpleNamespace(phi=np.full(4,.8),psi_child=.2,H0=np.array([6.]),pension=.3)
            runtime = SimpleNamespace(P=P,reference_price=.7,arm="hard",total_native_calls=0,
                stationary_state=lambda packet,scale: SimpleNamespace(g_pre=np.ones((1,))),
                terminal_checks=lambda *args,**kwargs: dict(all_checks_pass=True),
                identity=lambda: dict(test=True))
            runtime.scaffold=SimpleNamespace(dump_checkpoint=lambda path,obj:Path(path).write_bytes(b"mock-state"))
            observed_phi=[]
            def stationary(psi,price,out):
                observed_phi.append(float(runtime.P.phi[0]))
                Path(out).mkdir(parents=True)
                return (dict(parameters=SimpleNamespace(phi=runtime.P.phi.copy(),pension=.3)),
                        dict(price=price,population_scale=1.2,renewal_residual=0.,
                             accounting_valid=True,gates=dict(native=True)))
            runtime.stationary=stationary
            def fake_map(runtime, **kw):
                Path(kw["folder"]).mkdir(parents=True)
                h=len(kw["prices"])
                fertility=[dict(period=i,calendar_year=2007+4*i,birth_flow_first=np.array([.1]),
                                childless_at_risk_mass=np.array([.2]),first_birth_hazard=np.array([.5]))
                           for i in range(h)]
                return (SimpleNamespace(rows=[dict(period=i) for i in range(h)],dated_states={4:dict(mock=True)}),
                        dict(gates=dict(native=True),market_residual=[0.]*h,
                             fiscal_residual=[0.]*h,fertility=fertility,
                             phi_path=financing_path(kw["kind"],h).tolist(),
                             rows=[dict(period=i,asset_price=float(kw["prices"][i]),
                                        pension_period_units=float(kw["pensions"][i])) for i in range(h)]))
            plan=dict(endpoint=dict(price_bound_ratios=[.05,20.],slope=1.,max_log_step=.15,
                                    damping=.7,max_evaluations=18),
                      path=dict(price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],
                                max_log_step=.15,damping=1.,max_evaluations=12),
                      fit=dict(max_condition_number=1e8,worsening_factor=1.5),
                      gates=dict(stationary_renewal_tolerance=1e-6,final_reproduction_tolerance=1e-10,
                                 market_tolerance=2e-4,fiscal_tolerance=2e-5,terminal_tolerance=1e-3,
                                 raw_queue_relative_tolerance=1e-3))
            def endpoint_root(**kw):
                kw["evaluate"](np.array([.71]));kw["evaluate"](np.array([.72]))
                return dict(converged=True,final=dict(prices=[.72]))
            def path_root(**kw):
                kw["evaluate"](kw["initial_prices"],kw["initial_fiscal_values"])
                return dict(converged=True,gates=dict(replay=True),final=dict(
                    prices=kw["initial_prices"].tolist(),
                    fiscal_values=kw["initial_fiscal_values"].tolist()))
            with (patch.object(run_case,"solve_price_path_scaled",endpoint_root),
                  patch.object(run_case,"solve_joint_with_acceleration",path_root),
                  patch.object(run_case,"map_case",fake_map)):
                terminal, endpoint=run_case.solve_terminal(runtime,folder/"terminal",plan,time.monotonic()+60)
                self.assertEqual(observed_phi,[1.,1.])
                self.assertEqual(runtime.P.phi.tolist(),[.8]*4)
                result=run_case.solve_path(runtime,"permanent",12,terminal,endpoint,
                                           np.eye(24),folder/"path",plan,time.monotonic()+60)
            self.assertEqual(result["first_births"][0]["flow"],.1)
            self.assertEqual(result["phi_path"],[1.]*12)
            self.assertTrue((folder/"terminal/accepted.json").is_file())
            self.assertTrue((folder/"path/completed.json").is_file())

    def test_real_scalar_root_accepts_offroot_admissible_mapping(self):
        observations=[]
        def evaluate(prices):
            residual=float(prices[0]-.72)
            observations.append(residual)
            return dict(mapping_valid=True,residual=np.array([residual]))
        root=run_case.solve_price_path_scaled(
            initial_prices=np.array([.7]),evaluate=evaluate,
            project=lambda q: np.clip(q,.1,2.),slope=1.,market_tolerance=1e-6,
            max_log_step=.15,damping=.7,max_evaluations=18,
            deadline_monotonic=time.monotonic()+60,max_condition_number=1e8,
            worsening_factor=1.5,final_reproduction_tolerance=1e-10)
        self.assertTrue(root["converged"])
        self.assertNotIn("gates",root)
        self.assertTrue(any(abs(x)>1e-6 for x in observations))
        self.assertLessEqual(abs(observations[-1]),1e-6)


if __name__ == "__main__": unittest.main()
