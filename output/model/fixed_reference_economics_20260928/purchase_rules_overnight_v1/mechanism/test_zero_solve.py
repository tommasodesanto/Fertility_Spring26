"""Zero-model-solve contract checks for the isolated financing path."""
from __future__ import annotations

import ast
import copy
import difflib
import hashlib
import json
import symtable
import builtins
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
sys.path.insert(0, str(ROOT / "code/model/experiments/transition_readiness"))
import dated_phi
from integration import bind_phi_path, financing_path
import run_case
from selected_runtime import authenticate_selected
import selected_runtime


class PhiPathTests(unittest.TestCase):
    def test_full_saved_array_gate_rejects_all_material_differences(self):
        values={f"field_{i:03d}":np.array([float(i)]) for i in range(93)}
        values["choice"]=np.array([1],dtype=np.int64)
        values["nonfinite"]=np.array([float("nan"),1.])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"arrays.npz"
            np.savez(path,**values)
            with patch.object(selected_runtime,"full_native_arrays",return_value=values):
                self.assertEqual(selected_runtime.compare_full_native_arrays(path,{})["fields"],95)
            for name,replacement,error in (
                ("choice",np.array([2],dtype=np.int64),"discrete"),
                ("field_000",np.array([1e-5]),"numeric"),
                ("nonfinite",np.array([0.,1.]),"nonfinite")):
                altered=dict(values);altered[name]=replacement
                with patch.object(selected_runtime,"full_native_arrays",return_value=altered):
                    with self.assertRaisesRegex(RuntimeError,error):
                        selected_runtime.compare_full_native_arrays(path,{})
            omitted=dict(values);omitted.pop("field_000")
            with patch.object(selected_runtime,"full_native_arrays",return_value=omitted):
                with self.assertRaisesRegex(RuntimeError,"key set"):
                    selected_runtime.compare_full_native_arrays(path,{})

    def test_renderer_exception_only_for_different_artifact_renderers(self):
        from PIL import Image, PngImagePlugin
        import floor_runtime
        with tempfile.TemporaryDirectory() as tmp:
            original,current=[Path(tmp)/name for name in ("original","current")]
            for root,version in ((original,"3.7.1"),(current,"3.10.0")):
                (root/"standard_diagnostics").mkdir(parents=True)
                for i in range(17):
                    info=PngImagePlugin.PngInfo()
                    info.add_text("Software","Matplotlib version"+version)
                    Image.new("RGB",(1,1),(i,0,0)).save(root/"standard_diagnostics"/f"plot_{i:02d}.png",pnginfo=info)
            with (patch.object(selected_runtime,"compare_full_native_arrays",return_value={"fields":95}),
                  patch.object(floor_runtime,"compare_normalized_reference_reports",
                               side_effect=RuntimeError("Normalized reference standard plots differ"))):
                result=selected_runtime.compare_renderer_aware_reference(None,original,current,{},Path(tmp)/"arrays.npz")
                self.assertEqual(result["reason"],"cross_renderer_full_arrays_and_reports")
                for plot in (current/"standard_diagnostics").glob("*.png"):
                    info=PngImagePlugin.PngInfo();info.add_text("Software","Matplotlib version3.7.1")
                    Image.new("RGB",(1,1),(0,0,0)).save(plot,pnginfo=info)
                with self.assertRaisesRegex(RuntimeError,"standard plots differ"):
                    selected_runtime.compare_renderer_aware_reference(None,original,current,{},Path(tmp)/"arrays.npz")
                (current/"standard_diagnostics/plot_00.png").unlink()
                with self.assertRaisesRegex(RuntimeError,"plot count differs"):
                    selected_runtime.compare_renderer_aware_reference(None,original,current,{},Path(tmp)/"arrays.npz")
            # A changed calibration table cannot be classified as renderer drift.
            for plot in (current/"standard_diagnostics").glob("*.png"):
                info=PngImagePlugin.PngInfo();info.add_text("Software","Matplotlib version3.10.0")
                Image.new("RGB",(1,1)).save(plot,pnginfo=info)
            info=PngImagePlugin.PngInfo();info.add_text("Software","Matplotlib version3.10.0")
            Image.new("RGB",(1,1)).save(current/"standard_diagnostics/plot_00.png",pnginfo=info)
            with (patch.object(selected_runtime,"compare_full_native_arrays",return_value={"fields":95}),
                  patch.object(floor_runtime,"compare_normalized_reference_reports",
                               side_effect=RuntimeError("Normalized reference numeric value differs: target_fit.csv"))):
                with self.assertRaisesRegex(RuntimeError,"target_fit.csv"):
                    selected_runtime.compare_renderer_aware_reference(None,original,current,{},Path(tmp)/"arrays.npz")

    @staticmethod
    def valid_stationary_record():
        return dict(price=.7, population_scale=1.2, renewal_residual=0.,
                    accounting_valid=True, absolute_housing_demand=6.,
                    absolute_housing_supply=6., gates=dict(
                        household_budget={}, purchase={}, estate={}, policy_arrays={},
                        stationary_operator={}, housing_market_clearing_required=False,
                        feasibility_projection_mass=0., fiscal_certificate=dict(
                            marginal_gate=True, fiscal_gate=True,
                            marginal_tolerance=1e-9, fiscal_tolerance=1e-6)))

    def test_stationary_audit_interprets_descriptive_false_and_rejects_failures(self):
        controls=dict(stationary_renewal_tolerance=1e-6,market_tolerance=2e-4)
        base=self.valid_stationary_record()
        self.assertTrue(run_case.stationary_valid(base,controls))
        for change in (
            lambda r:r.update(accounting_valid=False),
            lambda r:r['gates']['fiscal_certificate'].update(fiscal_gate=False),
            lambda r:r['gates'].update(feasibility_projection_mass=1e-8),
            lambda r:r.update(renewal_residual=2e-6),
            lambda r:r.update(absolute_housing_supply=5.),
            lambda r:r.update(absolute_housing_demand=.01001,absolute_housing_supply=.01),
            lambda r:r.update(absolute_housing_demand=float('nan')),
        ):
            bad=copy.deepcopy(base);change(bad)
            self.assertFalse(run_case.stationary_valid(bad,controls))

    def test_copied_path_globals_are_bound(self):
        source = (Path(__file__).parent / "dated_phi.py").read_text()
        table = symtable.symtable(source, "dated_phi.py", "exec")
        module_names = set(dated_phi.__dict__) | set(dir(builtins))
        for name in ("backward_value_path", "evaluate_path_at_prices"):
            function = next(child for child in table.get_children() if child.get_name() == name)
            missing = sorted(symbol.get_name() for symbol in function.get_symbols()
                             if symbol.is_global() and symbol.is_referenced()
                             and symbol.get_name() not in module_names)
            self.assertEqual(missing, [], name)
        self.assertIs(dated_phi.validate_entry_queues, dated_phi.base.validate_entry_queues)

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
                        dict(self.valid_stationary_record(),price=price))
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
