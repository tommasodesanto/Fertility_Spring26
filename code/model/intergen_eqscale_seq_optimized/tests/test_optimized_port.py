"""Regression gates specific to the eqscale/sequential optimized port."""
from __future__ import annotations
import numpy as np
import pytest
from unittest.mock import patch
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters
from intergen_eqscale_seq_optimized.calibration import get_target_system
from intergen_eqscale_seq_optimized.solver import (precompute_shared, run_model_cp_dt,
    solve_markov_income_at_prices, upgrade_fast_markov_solution, refine_one_market_markov_income)
from intergen_eqscale_seq_optimized.calibration import extract_moments

def tiny(**extra):
    base={"J":6,"J_R":5,"A_f_start":1,"A_f_end":4,"Nb":20,"b_min":0.,"b_core_lo":0.,"b_core_hi":3.,"b_mid_hi":6.,"b_max":10.,"n_house":2,"H_own":np.array([2.,4.]),"H0":np.array([4.]),"eta_supply":np.array([1.]),"use_income_types":True,"income_type_transition":"persistent","z_grid":np.array([.8,1.2]),"z_weights":np.array([.5,.5]),"Pi_z":np.array([[.9,.1],[.1,.9]]),"entry_wealth_mode":"scalar","b_entry_fixed":.5,"c_bar_0":.05,"c_bar_n":.02,"h_bar_0":.25,"h_bar_jump":.1,"h_bar_n":.05,"lambda_d":0.,"use_full_kernel":True,"use_tenure_kernel":True,"use_loc_kernel":True,"tenure_choice_kappa":0.,"sequential_births":True,"max_iter_eq":4,"scalar_market_refine_max_expand":4,"tol_eq":1e-3}
    base.update(extra); return base

def test_direct_payload_restores_sequential_side_channel_after_stale_bellman():
    # The direct run establishes the constructed multi-price search condition.
    direct, P, price = run_model_cp_dt(tiny(markov_equilibrium_method="direct_brent"), verbose=False)
    assert direct.timings["unique_fast_price_evaluations"] >= 2
    # Reproduce the critical cache ordering: accepted price A, then rejected B
    # overwrites P._fert2_probs.  Upgrade(A) must restore A before the KFE.
    grid=np.asarray(direct.b_grid); SD=precompute_shared(P,grid)
    accepted=solve_markov_income_at_prices(price,P,grid,fast_stats=True,SD=SD,retain_payload=True)
    rejected=solve_markov_income_at_prices(np.asarray(price)*1.15,P,grid,fast_stats=True,SD=SD,retain_payload=True)
    assert not np.array_equal(accepted._model_payload[-1], rejected._model_payload[-1])
    upgraded=upgrade_fast_markov_solution(accepted,P,grid,SD)
    upgraded_flows = P._second_births_by_age.copy()
    upgraded_attempts = P._second_attempts_by_age.copy()
    fresh=solve_markov_income_at_prices(price,P,grid,fast_stats=False,SD=SD)
    np.testing.assert_array_equal(upgraded_flows, P._second_births_by_age)
    np.testing.assert_array_equal(upgraded_attempts, P._second_attempts_by_age)


def test_integrated_direct_payload_freshness_after_last_price_is_rejected():
    """The accepted payload, not the final rejected Bellman, drives the KFE."""
    from intergen_eqscale_seq_optimized import solver as solver_module
    from intergen_eqscale_seq_optimized.e1_profile import E1_PRICE
    observed_prices: list[float] = []
    original = solver_module.solve_markov_income_at_prices

    def traced(*args, **kwargs):
        if kwargs.get("retain_payload", False):
            observed_prices.append(float(np.asarray(args[0]).reshape(-1)[0]))
        return original(*args, **kwargs)

    overrides = tiny(
        markov_equilibrium_method="direct_brent", p_init_override=np.array([1.30 * E1_PRICE]),
        scalar_market_refine_max_expand=0, max_iter_eq=4, tol_eq=1e-3,
    )
    with patch.object(solver_module, "solve_markov_income_at_prices", side_effect=traced):
        integrated, P, accepted_price = solver_module.run_model_cp_dt(overrides, verbose=False)
    accepted = float(accepted_price[0])
    assert len(observed_prices) >= 3
    assert observed_prices[-1] != accepted
    integrated_births, integrated_attempts = P._second_births_by_age.copy(), P._second_attempts_by_age.copy()
    fresh = original(accepted_price, P, integrated.b_grid, fast_stats=False)
    np.testing.assert_array_equal(integrated.g, fresh.g)
    np.testing.assert_array_equal(integrated_births, P._second_births_by_age)
    np.testing.assert_array_equal(integrated_attempts, P._second_attempts_by_age)

def test_contracts_and_atomic_target_system():
    with pytest.raises(ValueError, match="Unknown parameter"):
        apply_overrides(setup_parameters(), {"definitely_typoed_parameter": 1})
    with pytest.raises(ValueError, match="inconsistent"):
        apply_overrides(setup_parameters(), {"Pi_z":np.array([[1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]),"z_weights":np.array([.2,.3,.5])})
    system=get_target_system("candidate_replacement_income_disciplined_v1")
    assert system.count == 15 and system.fingerprint == get_target_system(system.name).fingerprint


def test_default_nesting_flags_reproduce_production_bitwise():
    """New defaults are opt-in only; these flags restore the inherited nest."""
    from intergen_housing_fertility.solver import run_model_cp_dt as run_production
    overrides = tiny(
        sequential_births=False, fecundity_omega1=0.0,
        preference_spec="stone_geary", markov_equilibrium_method="legacy_damped",
        use_numba_scatter=False, solve_mode="pe", p_fixed=np.array([1.0]),
        w_fixed=np.array([1.0]), entry_shares_fixed=np.array([1.0]),
    )
    optimized, _, _ = run_model_cp_dt(overrides, verbose=False)
    production, _, _ = run_production(overrides, verbose=False)
    for name in ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "g"):
        np.testing.assert_array_equal(getattr(optimized, name), getattr(production, name))


def test_direct_vs_legacy_tiny_e1_root_agreement_and_price_count():
    common = tiny(max_iter_eq=40, tol_eq=1e-3)
    legacy, P_legacy, p_legacy = run_model_cp_dt({**common, "markov_equilibrium_method": "legacy_damped"}, verbose=False)
    direct, P_direct, p_direct = run_model_cp_dt({**common, "markov_equilibrium_method": "direct_brent"}, verbose=False)
    assert legacy.best_max_abs_rel_excess <= P_legacy.tol_eq
    assert direct.best_max_abs_rel_excess <= P_direct.tol_eq
    assert abs(float(p_direct[0]) - float(p_legacy[0])) < 1e-3
    assert direct.timings["unique_fast_price_evaluations"] < legacy.timings["unique_fast_price_evaluations"]
    for name in ("tfr", "childless_rate", "own_rate", "housing_increment_0to1"):
        np.testing.assert_allclose(extract_moments(direct, P_direct)[name], extract_moments(legacy, P_legacy)[name], atol=2e-4, rtol=0.0)


def test_repeat_price_cache_hit_preserves_result():
    overrides = tiny(markov_equilibrium_method="direct_brent", max_iter_eq=4, tol_eq=1e-3)
    initial, P, price = run_model_cp_dt(overrides, verbose=False)
    grid, SD, cache = np.asarray(initial.b_grid), precompute_shared(P, initial.b_grid), {}
    best = solve_markov_income_at_prices(price, P, grid, fast_stats=True, SD=SD, retain_payload=True)
    key = tuple(float(x) for x in np.asarray(price).reshape(-1)); summary = __import__("copy").copy(best); delattr(summary, "_model_payload"); cache[key] = summary
    first = refine_one_market_markov_income(price, best, initial.best_max_abs_rel_excess, P, grid, verbose=False, SD=SD, price_cache=cache)
    second = refine_one_market_markov_income(price, best, initial.best_max_abs_rel_excess, P, grid, verbose=False, SD=SD, price_cache=cache)
    assert second[3]["cache_hits"] > 0
    np.testing.assert_array_equal(first[1], second[1])
    assert first[2] == second[2]
