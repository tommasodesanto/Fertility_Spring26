import math
import os
import sys
import tempfile
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.ces_normalized_shares import adapter
from production.engine import shared


def _point(): return dict(adapter.PARAMETERS)


def test_mapping_bounds_and_no_reference_rent():
    spec=adapter.contract(); assert len(spec['bounds']) == 11 and 'h_P' not in spec['bounds']
    assert spec['bounds']['delta_alpha_jump'] == spec['bounds']['delta_alpha'] == (0., .25)
    family=[r for r in spec['target_fit'] if r['moment']=='family_rooms'][0]
    assert family['role']=='scored' and family['weight']=='280.52808370152104'
    P,g=adapter.load_inputs(_point()); assert not hasattr(P, 'utility_reference_rent') or P.utility_reference_rent != 0
    assert P.child_room_floor is False and P.hbar_first_child_jump == 0


def test_multiplier_all_states_and_sigma_cases():
    for sigma in (1.5,2.,3.):
        P,g=adapter.load_inputs(_point());P.sigma=sigma
        with adapter.install(): sd=shared.precompute_shared(P,g)
        escale=sd.escale_flat.reshape(P.n_parity,P.n_child_states,order='F')
        for n in range(P.n_parity):
            for cs in range(P.n_child_states):
                m=cs if 0 <= cs <= n else 0
                a=.733 if m == 0 else np.clip(.733-P.delta_alpha_jump-P.delta_alpha*m,.05,.95)
                want=((((2+.7*m)/2)**.7)*(a**a*(1-a)**(1-a)))**(sigma-1)
                assert math.isclose(escale[n,cs],want,rel_tol=0,abs_tol=2e-14)
        assert not math.isclose(escale[0,0],1.)


def test_optimized_composite_and_ratio():
    # At equal unit prices, max c^a s^(1-a)/(a^a(1-a)^(1-a)) subject to c+s=E is E at c/s=a/(1-a).
    E=7.;a=.633;c=a*E;s=(1-a)*E
    assert math.isclose(c**a*s**(1-a)/(a**a*(1-a)**(1-a)),E,rel_tol=0,abs_tol=1e-13)
    assert math.isclose(c/s,a/(1-a),rel_tol=0,abs_tol=1e-14)


def test_independent_count_children_and_benefits_are_analytical():
    P, grid = adapter.load_inputs(_point())
    assert str(P.child_state_mode) == "independent_count"
    assert [adapter._children_at_home(P, 3, cs) for cs in (1, 0, 2)] == [1, 0, 2]
    alpha = np.empty((P.n_parity, P.n_child_states)); benefit = np.empty_like(alpha); multiplier = np.empty_like(alpha)
    adapter.normalized_share_callback(P, alpha, benefit, multiplier)
    for cs, m in ((1, 1), (0, 0), (2, 2)):
        assert benefit[3, cs] == (0. if m == 0 else P.psi_child * m ** (1. - P.child_benefit_curvature))
        a = .733 if m == 0 else np.clip(.733-P.delta_alpha_jump-P.delta_alpha*m,.05,.95)
        expected = ((((2. + .7 * m) / 2.) ** .7) * (a ** a * (1. - a) ** (1. - a))) ** (P.sigma - 1.)
        assert math.isclose(multiplier[3, cs], expected, rel_tol=0, abs_tol=2e-14)
    P.child_state_mode = "shared_clock"
    with pytest.raises(ValueError, match="independent_count"):
        adapter._children_at_home(P, 3, 1)


def test_reference_rent_does_not_change_shared_arrays():
    P, grid = adapter.load_inputs(_point()); P.utility_reference_rent = 1.
    Q, _ = adapter.load_inputs(_point()); Q.utility_reference_rent = 999.
    with adapter.install():
        left, right = shared.precompute_shared(P, grid), shared.precompute_shared(Q, grid)
    for key in ("alpha_flat", "psi_flat", "escale_flat"):
        assert np.array_equal(getattr(left, key), getattr(right, key))


def test_clipped_jump_slope_rule_for_matured_and_all_parent_counts():
    P, _ = adapter.load_inputs(_point())
    P.delta_alpha_jump=.12; P.delta_alpha=.20
    alpha=np.empty((P.n_parity,P.n_child_states)); benefit=np.empty_like(alpha); multiplier=np.empty_like(alpha)
    adapter.normalized_share_callback(P,alpha,benefit,multiplier)
    assert alpha[0,0] == .733
    for m in (1,2,3):
        assert alpha[3,m] == max(.05, min(.95, .733-.12-.20*m))


def test_experimental_rescore_promotes_only_family_rooms_and_changes_fingerprint():
    report = adapter.ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge/selected_root'
    residual, fits, _ = adapter.residual_from_report(report)
    spec = adapter.contract()
    family = [r for r in fits if r['moment'] == 'family_rooms'][0]
    assert family['role'] == 'scored' and family['weight'] == '280.52808370152104'
    assert float(family['loss_contribution']) == float(family['weight']) * float(family['gap']) ** 2
    assert residual.shape == (11,)
    assert spec['target_fingerprint'] != spec['baseline_target_fingerprint']
    assert spec['weight_fingerprint'] != spec['baseline_weight_fingerprint']


def test_installation_patches_actual_equilibrium_lookup_and_restores(monkeypatch, tmp_path):
    from production import equilibrium
    P, grid = adapter.load_inputs(_point())
    original = equilibrium.build_context; calls = []
    def fake_context(*args, **kwargs):
        calls.append((args, kwargs)); return {"context": "hooked"}
    monkeypatch.setattr(adapter, "build_reporting_context", fake_context)
    def fake_solver(Q, b, **kwargs):
        assert equilibrium.build_context(Q, b, Path(kwargs["out"]), price_start=1., deadline=1., max_lifecycle=2, closure="population_one") == {"context": "hooked"}
        return {"status": "rejected", "reason": "unit hook only", "lifecycle_solves": 0}
    evaluate = adapter.make_evaluator(tmp_path, "unit", P, grid, deadline=10**12, solver=fake_solver)
    assert evaluate("hook", _point(), 10**12)["status"] == "rejected"
    assert len(calls) == 1
    assert equilibrium.build_context is original


@pytest.mark.actual_context
def test_actual_context_preflight_is_staged_only(tmp_path):
    if os.environ.get("CES_NORMALIZED_SHARES_STAGED_CONTEXT") != "1":
        pytest.skip("actual frozen reporting context requires the isolated staged dependency snapshot")
    with adapter.install(): pass
    receipt=adapter.preflight_contexts(tmp_path)
    assert receipt['actual_lifecycle_solves']==0 and len(receipt['contexts'])==3
    assert all(x['parameter_rows']==31 and x['expected_parameter_rows']==31 for x in receipt['contexts'])
