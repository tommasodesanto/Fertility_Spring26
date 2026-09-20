from types import SimpleNamespace
import importlib
import numpy as np
import pytest
import run_e5f_financing_factorial as d


def test_factorial_and_smoke_order_contract():
    full = [(p, l, c) for p in d.PHIS for l in d.LAMBDAS for c in d.CAPS]
    assert len(set(full)) == 8
    assert [(0.8, 0.0, 6.0)] * 2 + [(1.0, 5.0, 10.0)] == [(0.8, 0.0, 6.0), (0.8, 0.0, 6.0), (1.0, 5.0, 10.0)]


def test_dose_design_has_48_grid_arms_and_two_controls():
    full = d.case_order("run", "dose")
    assert len(full) == 50
    assert len(set(full[2:])) == 48
    assert full[:2] == [(0.8, 0.0, 6.0)] * 2


def test_smoke_and_design_defaults():
    for design in ("binary", "dose"):
        assert len(d.case_order("smoke", design)) == 3
        assert d.case_order("smoke", design)[-1] == (1.0, 5.0, 10.0)
    assert len(d.case_order("run")) == 10
    with pytest.raises(ValueError, match="invalid design"):
        d.case_order("run", "invalid")


def test_family_names_are_explicit():
    assert d.FAMILY_NAMES == ("original", "stationary_new_income", "refit_new_income")


def test_population_source_default_is_identity_and_saved_evaluation_is_copy():
    raw = np.zeros((2, 2)); evaluated = raw.copy(); evaluated[0, 0] = 5e-13
    x = {"stationary_g_pre": raw, "evaluation": SimpleNamespace(g_pre=evaluated)}
    assert d.prepare_population(x) is x
    projected = d.prepare_population(x, "saved_evaluation")
    assert projected is not x
    assert np.array_equal(projected["stationary_g_pre"], evaluated)
    assert projected["stationary_g_pre"] is not evaluated
    assert np.array_equal(raw, np.zeros((2, 2)))


def test_population_source_rejects_substantive_difference():
    raw = np.zeros((2, 2)); evaluated = raw.copy(); evaluated[0, 0] = 2e-12
    with pytest.raises(ValueError, match="population mismatch"):
        d.prepare_population({"stationary_g_pre": raw, "evaluation": SimpleNamespace(g_pre=evaluated)}, "saved_evaluation")


def test_population_source_rejects_nonfinite_values():
    raw = np.zeros((2, 2)); evaluated = raw.copy(); evaluated[0, 0] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        d.prepare_population({"stationary_g_pre": raw, "evaluation": SimpleNamespace(g_pre=evaluated)}, "saved_evaluation")


def test_run_case_uses_exact_population_baseline_graph_and_cohort(monkeypatch, tmp_path):
    shape = (1, 2, 1, 1, 1, 1, 1); g = np.zeros(shape); g[0, 0, 0, 0, 0, 0, 0] = .4; g[0, 1, 0, 0, 0, 0, 0] = .6
    policy = SimpleNamespace(**{name: np.ones((1,)) for name in d.POLICY_NAMES})
    ev = SimpleNamespace(g_pre=g.copy(), g_post_fertility=g.copy(), g_current=g.copy(), births=np.array([.1]), policy=policy)
    x = {"parameters": SimpleNamespace(phi=np.array([.8]), lambda_d=0., hR_max=6.), "stationary_g_pre": g.copy(), "evaluation": SimpleNamespace(policy=policy, g_current=g.copy(), births=np.array([.1])), "b_grid": np.array([0., 1.]), "supply_rule": "fixed"}
    monkeypatch.setattr(d, "arm_parameters", lambda *args: SimpleNamespace())
    monkeypatch.setattr(d, "metrics", lambda ev, P: {"birth_flow": .1, "first_birth_flow": .2, "mean_rooms": 2., "ownership": .6, "renter_mass": .4, "owner_mass": .6, "market_residual": 0., "pre_mass": 1., "tenure_probability_sum": 1.})
    calls = []; monkeypatch.setattr(d, "run_cohort", lambda *args: calls.append(args) or {"status": "completed"}, raising=False)
    class Native:
        compare_exact = staticmethod(lambda a, b: (assert_same(a, b)))
    def assert_same(a, b): assert a is policy and b is policy
    class Rental:
        native_solve = staticmethod(lambda x, P: (ev, {"budget_excess_mass": 0., "maximum_occupied_excess": 0.}, object(), object()))
        gates = staticmethod(lambda ev: None)
        standard_graphs = staticmethod(lambda *args: {"status": "completed", "count": 17})
    class Audit: policy_array_audit = staticmethod(lambda *args: {"occupied_negative_steps": 0})
    real = importlib.import_module
    monkeypatch.setattr(d.importlib, "import_module", lambda name: Native if name == "run_e5f_native_financing_diagnostic" else Rental if name == "run_e5f_native_rental_access_diagnostic" else Audit if name == "run_e5f_independent_numerical_audit" else real(name))
    result = d.run_case(x, .8, 0., 6., tmp_path / "case", True)
    assert result["standard_graphs"]["count"] == 17 and len(calls) == 1


def test_run_case_rejects_changed_initial_population(monkeypatch, tmp_path):
    g = np.zeros((1, 2, 1, 1, 1, 1, 1)); bad = g.copy(); bad[0, 0, 0, 0, 0, 0, 0] = 1.
    ev = SimpleNamespace(g_pre=bad, g_post_fertility=bad, g_current=bad, births=np.array([0.]), policy=SimpleNamespace())
    x = {"parameters": SimpleNamespace(), "stationary_g_pre": g, "evaluation": SimpleNamespace(policy=ev.policy, g_current=g, births=np.array([0.])), "b_grid": np.array([0.]), "supply_rule": "fixed"}
    monkeypatch.setattr(d, "arm_parameters", lambda *args: SimpleNamespace())
    class Rental:
        native_solve = staticmethod(lambda *args: (ev, {"budget_excess_mass": 0., "maximum_occupied_excess": 0.}, object(), object()))
        gates = staticmethod(lambda ev: None)
    monkeypatch.setattr(d.importlib, "import_module", lambda name: Rental if name == "run_e5f_native_rental_access_diagnostic" else SimpleNamespace())
    with pytest.raises(ValueError, match="initial population"): d.run_case(x, .8, 0., 6., tmp_path / "bad", False)


def test_parameter_change_contract_only_finance_and_cap(monkeypatch):
    base = SimpleNamespace(phi=np.array([.8]), lambda_d=0., hR_max=6., H_own=np.array([4., 10.]), debt_taper_start_age=80., debt_taper_end_age=84., debt_taper_weights=np.ones(2), debt_caps=np.ones(2))
    from run_e5f_native_financing_diagnostic import changed
    monkeypatch.setattr(d.importlib, "import_module", lambda name: SimpleNamespace(build_debt_caps=lambda P: None,changed=changed))
    altered = d.arm_parameters(base, 1., 5., 10.)
    changed = {k for k in vars(base) | vars(altered) if not np.array_equal(getattr(base, k, None), getattr(altered, k, None))}
    assert changed <= {"phi", "lambda_d", "hR_max", "debt_taper_start_age", "debt_taper_end_age", "debt_taper_weights", "debt_caps"}
