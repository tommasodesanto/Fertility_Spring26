"""Small cached-workflow tests; fixtures never import or solve the native model."""
from pathlib import Path
import runpy
import sys
import json
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pytest

from . import equilibrium, inputs, workflow
from .inputs import load_inputs
from .storage import StoredResult, load_case, load_latest, publish_latest, reserve_case, save_case


def _fixture_result():
    shape = (2, 6, 1, 2, 2, 2, 2)
    solution = SimpleNamespace(
        b_grid=np.array([-1.0, 2.0]),
        V=np.arange(np.prod(shape), dtype=float).reshape(shape),
        c_pol=np.full(shape, 0.4), hR_pol=np.full(shape, 2.0),
        bp_pol=np.full(shape, 0.1), c_pol_stay=np.full(shape, 0.3),
        bp_pol_stay=np.full(shape, 0.2),
        tenure_probs=np.full(shape + (6,), 1 / 6),
        fert_probs=np.full(shape[:-1] + (2,), 0.5),
        fert2_probs=np.full(shape[:-1] + (2,), 0.5),
        g_beginning_distribution=np.full(shape, 1 / np.prod(shape)),
        g=np.full(shape, 1 / np.prod(shape)),
        g_stay_distribution=np.full(shape, 1 / np.prod(shape)),
        type_values=np.array([[1.0, np.nan], [2.0, 3.0]]),
        nested={"values": (np.array([1, 2], dtype=np.int16), float("nan"))},
    )
    P = SimpleNamespace(H_own=np.array([2, 4, 6, 8, 10]), age_start=18,
                        period_years=4, R_gross=1.08, psi=0.06,
                        purchase_timing="inherited_only")
    return StoredResult(solution, P, solution.b_grid.copy(), 0.776,
                        parameters={"beta": 0.96}, label="fixture")


def _fake_core(monkeypatch):
    fixture = _fixture_result()

    def load_inputs(parameters=None, external_inputs=None, native_overrides=None):
        return fixture.P, fixture.b_grid.copy()

    def solve_stationary_ge(P, grid, *, out, price_start, budget_seconds, closure):
        out = Path(out)
        out.mkdir(parents=True, exist_ok=False)
        report = out / "selected_root"
        report.mkdir()
        (report / "target_fit.csv").write_text("moment,target,model,gap,weight,loss\n" +
                                              "m,1,1,0,1,0\n" * 14)
        (report / "parameters.csv").write_text("parameter,value\n" + "x,1\n" * 31)
        plots = report / "standard_diagnostics"
        plots.mkdir()
        for index in range(17):
            (plots / f"diagnostic_{index:02d}.png").write_bytes(b"fixture")
        return {"solution": fixture.solution, "P": fixture.P,
                "b_grid": fixture.b_grid, "price": 0.776,
                "closure": {"name": "fixed_h0", "status": "converged",
                            "renewal_residual": 0.0, "population_scale": 1.1,
                            "implied_N1_H0": 6.4},
                "report_directory": str(report)}

    monkeypatch.setattr(inputs, "load_inputs", load_inputs)
    monkeypatch.setattr(equilibrium, "solve_stationary_ge", solve_stationary_ge)


def _publish_seed(root: Path):
    complete = reserve_case(root)
    seed = _fixture_result()
    seed.closure = {"name": "seed"}
    seed.report_directory = "seed-report"
    save_case(seed, complete, metadata={"closure": seed.closure,
                                        "report_directory": seed.report_directory})
    publish_latest(complete, root)
    return complete.resolve()


def test_roundtrip_restores_namespaces_and_exact_recursive_values():
    with TemporaryDirectory() as directory:
        case = reserve_case(Path(directory))
        original = _fixture_result()
        save_case(original, case, metadata={})
        loaded, _ = load_case(case)
        assert hasattr(loaded.solution, "V")
        assert hasattr(loaded.P, "H_own")
        assert np.array_equal(loaded.solution.V, original.solution.V)
        assert np.array_equal(loaded.solution.type_values, original.solution.type_values,
                              equal_nan=True)
        assert loaded.solution.nested["values"][0].dtype == np.int16
        assert np.isnan(loaded.solution.nested["values"][1])


@pytest.mark.parametrize("failure_point", ["save", "plots", "publish"])
def test_failed_attempt_preserves_previous_latest(monkeypatch, failure_point):
    _fake_core(monkeypatch)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        previous = _publish_seed(root)
        monkeypatch.setattr(workflow, "_cached_plots", lambda case: None)

        def fail(*args, **kwargs):
            raise RuntimeError(f"injected {failure_point} failure")

        if failure_point == "save":
            monkeypatch.setattr(workflow, "save_case", fail)
        elif failure_point == "plots":
            monkeypatch.setattr(workflow, "_cached_plots", fail)
        else:
            monkeypatch.setattr(workflow, "publish_latest", fail)

        with pytest.raises(RuntimeError, match="injected"):
            workflow.run_stationary({}, {}, price_guess=0.7, output_root=root)
        _, latest_case = load_latest(root)
        assert latest_case == previous
        attempt = next(path for path in (root / "cases").iterdir()
                       if path.resolve() != previous)
        failure = json.loads((attempt / "failure.json").read_text())
        assert failure["exception_type"] == "RuntimeError"
        assert "injected" in failure["message"]


def test_workflow_publishes_explorer_case_after_roundtrip(monkeypatch):
    _fake_core(monkeypatch)
    monkeypatch.setattr(workflow, "_cached_plots", lambda case: None)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        fiscal_inputs = {"w_hat": [1.01], "income_age_profile": [1.0], "tau_pay": 0.08}
        result, case = workflow.run_stationary({}, fiscal_inputs, price_guess=0.7,
                                               output_root=root)
        loaded, latest_case = load_latest(root)
        assert latest_case == case.resolve()
        assert hasattr(loaded.solution, "V")
        assert result.closure["status"] == "converged"
        assert (case / "target_fit.csv").is_file()
        assert (case / "parameters.csv").is_file()
        assert (case / "explorer_arrays.npz").is_file()
        assert len(list((case / "standard_diagnostics").glob("*.png"))) == 17
        contract = json.loads((case / "input_contract.json").read_text())
        assert contract["price_guess"] == 0.7
        assert contract["closure"] == "fixed_h0"
        assert contract["budget_seconds"] == 1800
        assert contract["fiscal_mapping"]["edited_primitives"] == [
            "income_age_profile", "tau_pay", "w_hat"]
        assert "fixed-payroll balanced-pension" in contract["fiscal_mapping"]["status"]
        assert "reference entrant-wealth mapping retained" in contract["entry_law"]
        metadata = json.loads((case / "metadata.json").read_text())
        assert metadata["price_guess"] == 0.7
        assert metadata["input_contract_sha256"]
        config = json.loads((case / "explorer_cases.json").read_text())
        assert config["cases"][0]["timing"] == "inherited_only"
        tools_path = str(Path(__file__).resolve().parents[1] / "tools")
        if tools_path not in sys.path:
            sys.path.insert(0, tools_path)
        from economics_explorer import SavedCase
        saved = SavedCase(config["cases"][0], config["common"])
        assert saved.shape == result.solution.V.shape


def test_playbutton_gross_wage_edit_derives_disposable_income_and_pension():
    runner = runpy.run_path(Path(__file__).resolve().parents[1] / "run_model.py")
    from production.parameter_files import load_parameter_file
    config = load_parameter_file(runner["PARAMETER_FILE"])
    parameters = config["parameters"]
    external = config["external_inputs"]
    baseline, baseline_grid = load_inputs(parameters, external)
    edited = dict(external)
    edited["w_hat"] = [external["w_hat"][0] * 1.01]
    changed, changed_grid = load_inputs(parameters, edited)

    assert changed.w_hat[0] == edited["w_hat"][0]
    assert np.array_equal(changed.income_age_profile, baseline.income_age_profile)
    assert not np.array_equal(changed.income, baseline.income)
    assert changed.pension != baseline.pension
    assert np.array_equal(changed_grid, baseline_grid)
    assert np.array_equal(changed.fixed_reference_entry_conditional,
                          baseline.fixed_reference_entry_conditional)
    assert np.array_equal(changed.entry_wealth_ratio_nodes,
                          baseline.entry_wealth_ratio_nodes)
    assert np.array_equal(changed.entry_wealth_ratio_weights,
                          baseline.entry_wealth_ratio_weights)


def test_runner_help_and_unknown_arguments_cannot_start_workflow(monkeypatch, capsys):
    runner = runpy.run_path(Path(__file__).resolve().parents[1] / "run_model.py")
    with TemporaryDirectory() as directory:
        output_root = Path(directory) / "must_not_be_created"
        monkeypatch.setitem(runner, "OUTPUT_ROOT", output_root)
        with pytest.raises(SystemExit) as help_exit:
            runner["main"](["--help"])
        assert help_exit.value.code == 0
        assert "usage:" in capsys.readouterr().out.lower()
        assert not output_root.exists()
        with pytest.raises(SystemExit) as invalid_exit:
            runner["main"](["--unsupported"])
        assert invalid_exit.value.code == 2
        assert not output_root.exists()
