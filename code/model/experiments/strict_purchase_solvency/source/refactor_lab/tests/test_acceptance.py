"""Consolidated acceptance suite (Torch only; no model solve).

Run through verify_torch.sh PHASE=acceptance, which sets REFACTOR_ROOT,
REFACTOR_SOURCE_ROOT, REFACTOR_EXPORT, REFACTOR_BUNDLE_SHA, REFACTOR_CHECKPOINT.
Sections 1 also covers typed round trips and engine provenance; section 3
adds deterministic boundary cases on actual reference value columns.

The old engine is imported only as a validation oracle. Sections:
  1 input identity      bundle equals the authenticated checkpoint fields
  2 credit fixtures     reference rule equals the oracle; corrected rule,
                        sale boundary and death branch; census never repairs
  3 saving kernels      verbatim source equality; indexed variant bit-identical
"""
from __future__ import annotations

import inspect
import json
import math
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = (Path(os.environ["REFACTOR_ROOT"]) if os.environ.get("REFACTOR_ROOT")
        else Path(__file__).resolve().parents[4])
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]  # oracle + unpickling

from refactor_lab import credit, inputs  # noqa: E402
from refactor_lab.verification import saving  # noqa: E402
from intergen_eqscale_seq_optimized import kernels as oracle_kernels  # noqa: E402
from intergen_eqscale_seq_optimized import parameters as oracle_parameters  # noqa: E402
from intergen_eqscale_seq_optimized import solver as oracle_solver  # noqa: E402

CREDIT_FIELDS = ("J", "use_age_survival", "survival_probs", "debt_caps", "debt_taper_weights",
                 "owner_ltv_multipliers", "phi", "psi", "native_purchase_income",
                 "native_due_stayer_credit", "H_own")


@pytest.fixture(scope="session")
def manifest():
    return inputs.load_manifest(ROOT)


@pytest.fixture(scope="session")
def P(manifest):
    """Credit-relevant reference fields, read directly from the manifest."""
    fields = manifest["actual_serialized_parameters"]
    return SimpleNamespace(**{k: (np.asarray(fields[k]) if isinstance(fields[k], list) else fields[k])
                              for k in CREDIT_FIELDS})


# 1 ---------------------------------------------------------------- inputs
EXPORT = Path(os.environ["REFACTOR_EXPORT"]) if os.environ.get("REFACTOR_EXPORT") else None


def _require_export():
    if EXPORT is None or not os.environ.get("REFACTOR_BUNDLE_SHA"):
        pytest.fail("REFACTOR_EXPORT and REFACTOR_BUNDLE_SHA are required; input identity is not optional")
    return EXPORT / "inputs", os.environ["REFACTOR_BUNDLE_SHA"]


def test_bundle_identity_and_artifact_split(manifest):
    bundle, pin = _require_export()
    normal = inputs.load_inputs(bundle, ROOT, pin)
    full = inputs.load_inputs(bundle, ROOT, pin, load_artifacts=True)
    n_art = sum(inputs.is_artifact(k) for k in manifest["actual_serialized_parameters"])
    assert normal.identity["primitive_fields"] + n_art == 260 and not normal.artifacts
    assert len(full.artifacts) == n_art and not any(inputs.is_artifact(k) for k in vars(normal.parameters))
    assert normal.parameters.psi_child == inputs.PSI_CHILD and normal.parameters.native_due_stayer_credit is True
    with pytest.raises(RuntimeError):
        inputs.load_inputs(bundle, ROOT, "0" * 64)


def test_bundle_equals_checkpoint():
    """Independent re-read of the checkpoint (oracle classes) vs the bundle."""
    import gzip, pickle
    bundle, pin = _require_export()
    ckpt = Path(os.environ["REFACTOR_CHECKPOINT"])
    assert inputs.sha256_file(ckpt) == inputs.CHECKPOINT_SHA256
    with gzip.open(ckpt, "rb") as stream:
        packet = pickle.load(stream)
    loaded = inputs.load_inputs(bundle, ROOT, pin, load_artifacts=True)
    np.testing.assert_array_equal(loaded.b_grid, packet["b_grid"])
    np.testing.assert_array_equal(loaded.reference_price, np.asarray(packet["solution"].p_eq).reshape(-1))
    ref = vars(packet["parameters"])
    for k, v in list(vars(loaded.parameters).items()) + list(loaded.artifacts.items()):
        r = ref[k]
        assert type(v) is type(r), k
        if isinstance(r, np.ndarray):
            assert v.dtype == r.dtype and v.shape == r.shape and np.array_equal(v, r, equal_nan=r.dtype.kind == "f"), k
        else:
            assert inputs.serialized(v) == inputs.serialized(r), k


def test_typed_roundtrip():
    arrays = {}
    value = dict(a=float("nan"), b=float("-inf"), c=np.float64(0.1), d=np.int64(3), e=np.bool_(True),
                 f=(1, 2.5, None), g=[np.zeros(3, dtype=np.float32)], h="x")
    enc = inputs.encode(value, "k", arrays)
    dec = inputs.decode(json.loads(json.dumps(enc)), arrays)
    assert math.isnan(dec["a"]) and dec["b"] == float("-inf")
    for key in ("c", "d", "e"):
        assert type(dec[key]) is type(value[key]) and dec[key] == value[key]
    assert dec["f"] == (1, 2.5, None) and dec["g"][0].dtype == np.float32 and dec["h"] == "x"


def _check_kernels_transformation_chain(engine, mod, text):
    """materialized original kernels -> exact indexed transformation (make_indexed_stage.py) -> live file.

    Only exhaustive_saving_scalar changes: its body is saving.exhaustive_saving_indexed
    (renamed) preceded by four verbatim helpers; every other engine file is unchanged.
    """
    import ast, hashlib
    h = lambda t: hashlib.sha256(t.encode()).hexdigest()
    tr = json.loads((engine / "transform_receipt.json").read_text())
    promo = json.loads((engine / "promotion_receipt.json").read_text())
    k = tr["files"]["kernels.py"]
    assert k["sha256_before"] == mod["output_sha256"]          # starts from the materialized original
    assert h(text) == k["sha256_after"] == promo["chain"]["promoted_kernels_sha256"]
    assert promo["chain"]["tested_stage_kernels_sha256"] == k["sha256_after"]
    saving_path = Path(__file__).resolve().parents[1] / "verification" / "saving.py"
    saving_text = saving_path.read_text()
    assert h(saving_text) == tr["saving_sha256"]
    lines = saving_text.splitlines(keepends=True)
    seg = {}
    for node in ast.parse(saving_text).body:
        if isinstance(node, ast.FunctionDef):
            start = min([node.lineno] + [d.lineno for d in node.decorator_list])
            seg[node.name] = "".join(lines[start - 1:node.end_lineno])
    block = "".join(seg[n] + "\n\n" for n in ("_rank", "_interp_ranked", "_renter_value", "_owner_value")) + \
        seg["exhaustive_saving_indexed"].replace("def exhaustive_saving_indexed(", "def exhaustive_saving_scalar(", 1)
    assert h(block) == tr["new_block_sha256"] and block in text
    assert h(seg["exhaustive_saving_scalar"]) == tr["replaced_segment_sha256"]   # oracle copy = replaced original
    for fname, meta in tr["files"].items():
        if meta.get("unchanged"):
            assert h((engine / fname).read_text()) == meta["sha256"], fname


def test_engine_provenance():
    """Each engine file equals its receipt hash; kept bodies equal the source.

    solver.py was split by apply_split.py: every materialized solver definition
    must appear byte-identically in the stage module named by split_receipt, the
    split source must be the materialized solver, and the facade and stage files
    must match split_receipt hashes and export every definition.
    """
    import ast, hashlib
    h = lambda t: hashlib.sha256(t.encode()).hexdigest()
    engine = Path(__file__).resolve().parents[1] / "engine"
    receipt = json.loads((engine / "materialize_receipt.json").read_text())
    split = json.loads((engine / "split_receipt.json").read_text())
    source_root = Path(os.environ.get("REFACTOR_SOURCE_ROOT", receipt["source_root"]))
    assert split["source_sha256"] == receipt["modules"]["solver"]["output_sha256"]
    for mod, meta in split["modules"].items():
        assert h((engine / f"{mod}.py").read_text()) == meta["sha256"], mod
    dest = {tuple(d["names"]): d for d in split["definitions"]}
    stage_text = {m: (engine / f"{m}.py").read_text() for m in split["modules"]}
    for name, mod in receipt["modules"].items():
        path = source_root / mod["source"]
        if Path(mod["source"]).is_absolute() and not path.exists() and os.environ.get("REFACTOR_OVERLAY_ROOT"):
            path = Path(os.environ["REFACTOR_OVERLAY_ROOT"]) / Path(mod["source"]).name   # reviewed overlay, local copy
        src = path.read_bytes()
        assert hashlib.sha256(src).hexdigest() == mod["source_sha256"], name
        lines = src.decode().splitlines(keepends=True)
        edited = {e["nested_in"] for e in mod["import_edits"] if "nested_in" in e}
        text = None if name == "solver" else (engine / f"{name}.py").read_text()
        if name == "kernels":
            _check_kernels_transformation_chain(engine, mod, text)
        elif text is not None:
            assert h(text) == mod["output_sha256"], name
            ast.parse(text)
        for d in mod["kept"]:
            seg = "".join(lines[d["lines"][0] - 1:d["lines"][1]])
            assert h(seg) == d["sha256"], (name, d["names"])
            if name == "solver":
                target = dest[tuple(d["names"])]
                if d["names"][0] not in edited:
                    assert target["sha256"] == d["sha256"], d["names"]
                assert seg in stage_text[target["module"]] or d["names"][0] in edited, (target["module"], d["names"])
            elif name == "kernels" and d["names"] == ["exhaustive_saving_scalar"]:
                assert seg not in text   # the one intentionally transformed definition
            elif d["names"][0] not in edited:
                assert seg in text, (name, d["names"])
    facade = ast.parse(stage_text["solver"])
    exported = {a.asname or a.name for n in facade.body if isinstance(n, ast.ImportFrom) for a in n.names}
    assert exported == {n for d in split["definitions"] for n in d["names"]}


# 2 ---------------------------------------------------------------- credit
B = np.array([-3.0, -0.2558139535, -1e-12, 0.0, 0.7, 4.0])


def test_reference_renter_floor_matches_oracle(P):
    rule = credit.credit_rule("reference")
    for j in range(int(P.J)):
        expected = oracle_parameters.unsecured_debt_floor(
            B, float(P.debt_taper_weights[j + 1]), float(P.debt_caps[j + 1]))
        np.testing.assert_array_equal(rule.renter_floor(P, B, j), expected)
        np.testing.assert_array_equal(rule.renter_floor(P, B, j), oracle_solver.renter_borrowing_floor(P, B, j))


def test_reference_owner_floors_match_oracle(P):
    rule = credit.credit_rule("reference")
    for j in range(int(P.J)):
        for price in (0.9, 2.3):
            for house in P.H_own:
                base = -float(P.phi[0]) * price * float(house)
                buyer = oracle_solver.owner_borrowing_floor(P, B, base, j)
                np.testing.assert_array_equal(np.full(B.shape, rule.buyer_floor(P, j, price, house)), buyer)
                death = oracle_solver.native_due_death_floor(P, j, price, house)
                stayer = oracle_solver.native_due_owner_floor(
                    B, oracle_solver.effective_owner_collateral_floor(P, base, j), death_floor=death)
                np.testing.assert_array_equal(rule.stayer_floor(P, B, j, price, house), stayer)


def test_corrected_renter_floor_and_death_branch(P):
    for d_bar in (0.0, 0.5):
        rule = credit.credit_rule("corrected", d_bar)
        for j in range(int(P.J)):
            dies = j == int(P.J) - 1 or float(P.survival_probs[j]) < 1.0
            floor = rule.renter_floor(P, B, j)
            assert np.all(floor == (0.0 if dies else -d_bar)), (d_bar, j)
    with pytest.raises(ValueError):
        credit.credit_rule("corrected")          # no implicit d_bar
    with pytest.raises(ValueError):
        credit.credit_rule("corrected", -0.1)
    with pytest.raises(ValueError):
        credit.credit_rule("reference", 0.0)


def test_sale_boundary(P):
    rule = credit.credit_rule("corrected", 0.0)
    price, house = 1.7, 6.0
    edge = -(1.0 - float(P.psi)) * price * house
    allowed = rule.sale_allowed(P, np.array([edge, np.nextafter(edge, -np.inf), edge + 1e-9]), price, house)
    assert allowed.tolist() == [True, False, True]
    assert credit.credit_rule("reference").sale_allowed(P, np.array([edge - 5.0]), price, house).all()


def test_census_reports_without_repair():
    grid = np.array([-0.2558139535, 0.0, 1.0])
    mass = np.array([[2e-6, 0.0], [0.3, 0.2], [0.1, 0.0]])
    value = np.array([[-1e10, -1e10], [-3.0, -2.0], [-1.0, -1e10]])
    before = (mass.copy(), value.copy())
    census = credit.entrant_feasibility_census(mass, value, grid)
    assert not census["feasible"] and census["dead_mass"] == 2e-6 and len(census["cells"]) == 1
    np.testing.assert_array_equal(mass, before[0]); np.testing.assert_array_equal(value, before[1])


# 3 ---------------------------------------------------------------- saving
@pytest.mark.parametrize("name", ["interp_scalar", "renter_wedge_flow", "eval_renter_scalar",
                                  "eval_owner_scalar", "exhaustive_saving_scalar"])
def test_verbatim_extraction(name):
    src = lambda f: inspect.getsource(getattr(f, "py_func", f))
    assert src(getattr(saving, name)) == src(getattr(oracle_kernels, name))


def _grid(bundle_grid=None):
    if bundle_grid is not None:
        return bundle_grid
    return np.concatenate([np.linspace(-12, -5, 13)[:-1], np.linspace(-5, 7, 115), np.geomspace(7.2, 3000, 32)])


def _cases(bg, rng, count):
    ranges = [(bg[0] - 1, bg[-1] + 1), (-1.0, 3.0), (bg[3], bg[3]), (bg[10], bg[40]), (-0.3, 0.2)]
    for k in range(count):
        V = np.cumsum(rng.normal(0.2, 1.0, bg.size)) - 40.0       # nonconcave, sometimes decreasing
        if k % 5 == 0:
            V[: rng.integers(1, 20)] = -1e10                        # dead sentinels
        lo, hi = ranges[k % len(ranges)]
        if k % 7 == 3:
            lo, hi = sorted(rng.uniform(bg[0], 10.0, 2))
        yield dict(lo=float(lo), hi=float(hi), resources=float(rng.uniform(-2, 15)), continuation=V,
                   rent=float(rng.uniform(0.05, 0.6)), hb=float(rng.uniform(0, 1.5)),
                   cb=float(rng.uniform(0, 0.5)), pc=float(rng.normal()), hmax=float(rng.uniform(2, 8)),
                   owner_cost=float(rng.uniform(0, 2)), owner_K=float(rng.uniform(0.2, 2)))


def _call(fn, bg, c, owner, alpha=0.733, oms=-1.0, beta=0.86, es=1.0):
    return fn(c["lo"], c["hi"], c["resources"], c["continuation"], bg, c["rent"], c["hb"], c["cb"], c["pc"],
              c["hmax"], alpha, oms, beta, es, c["owner_cost"] if owner else 0.0,
              c["owner_K"] if owner else 0.0, owner)


@pytest.mark.parametrize("owner", [False, True])
@pytest.mark.parametrize("es", [1.0, 1.37])
def test_indexed_saving_bit_identical(owner, es):
    bundle, pin = _require_export()
    bg = inputs.load_inputs(bundle, ROOT, pin).b_grid
    rng = np.random.default_rng(20260929)
    mismatches = []
    for k, c in enumerate(_cases(bg, rng, 4000)):
        # Put the renter kink exactly on a grid node for a slice of cases.
        if not owner and k % 11 == 0:
            node = bg[int(rng.integers(20, 120))]
            c["resources"] = node + c["cb"] + c["rent"] * c["hb"] + c["rent"] * (c["hmax"] - c["hb"]) / (1 - 0.733)
        ref = _call(oracle_kernels.exhaustive_saving_scalar, bg, c, owner, es=es)
        new = _call(saving.exhaustive_saving_indexed, bg, c, owner, es=es)
        if not (ref[0] == new[0] and (ref[1] == new[1] or (np.isnan(ref[1]) and np.isnan(new[1])))):
            mismatches.append(dict(case=k, ref=ref, new=new))
    assert not mismatches, json.dumps(mismatches[:5], default=str)


def _boundary_cases(bg, V):
    """Deterministic edge cases requested in review."""
    ulp = lambda x, d: np.nextafter(x, d * np.inf)
    base = dict(rent=0.31, hb=0.6, cb=0.2, pc=-0.4, hmax=6.0, owner_cost=0.9, owner_K=1.1)
    for i in (0, 1, 7, len(bg) // 2, len(bg) - 2, len(bg) - 1):
        x = float(bg[i])
        for lo, hi in ((ulp(x, -1), ulp(x, 1)), (x, x), (ulp(x, 1), float(bg[min(i + 3, len(bg) - 1)])),
                       (float(bg[max(i - 3, 0)]), ulp(x, -1)), (x, float(bg[-1]) + 1.0), (float(bg[0]) - 1.0, x)):
            if hi < lo:
                continue
            for res in (lo + 0.5, lo + 5.0, hi + 50.0, lo - 1.0):
                yield dict(base, lo=lo, hi=hi, resources=res, continuation=V)
    # renter kink exactly on, and one ULP around, grid nodes
    for i in (5, 40, 90):
        for d in (-1, 0, 1):
            node = float(bg[i]) if d == 0 else ulp(float(bg[i]), d)
            dc = base["cb"] + base["rent"] * base["hb"]
            cap = base["rent"] * (base["hmax"] - base["hb"]) / (1 - 0.733)
            yield dict(base, lo=float(bg[0]), hi=float(bg[-1]), resources=node + dc + cap, continuation=V)
    # dead-sentinel feasibility boundaries
    for cut in (1, 10, 60):
        W = V.copy(); W[:cut] = -1e10
        yield dict(base, lo=float(bg[cut - 1]), hi=float(bg[cut + 5]), resources=float(bg[cut]) + 2.0, continuation=W)
        yield dict(base, lo=ulp(float(bg[cut]), -1), hi=ulp(float(bg[cut]), 1), resources=float(bg[cut]) + 1.0,
                   continuation=W)


def _reference_columns(limit=40):
    """Actual continuation columns from the saved reference value function."""
    path = EXPORT / "verification" / "reference_solution.npz"
    with np.load(path, allow_pickle=False) as z:
        V = z["V"]
    flat = V.reshape(V.shape[0], -1)
    idx = np.linspace(0, flat.shape[1] - 1, limit).astype(int)
    return [np.ascontiguousarray(flat[:, k], dtype=np.float64) for k in idx]


@pytest.mark.parametrize("owner", [False, True])
def test_indexed_saving_boundaries_and_model_columns(owner):
    bundle, pin = _require_export()
    bg = inputs.load_inputs(bundle, ROOT, pin).b_grid
    bad = []
    for col, V in enumerate(_reference_columns()):
        if not np.isfinite(V).all():
            continue
        for k, c in enumerate(_boundary_cases(bg, V)):
            for es in (1.0, 1.37):
                ref = _call(oracle_kernels.exhaustive_saving_scalar, bg, c, owner, es=es)
                new = _call(saving.exhaustive_saving_indexed, bg, c, owner, es=es)
                if np.float64(ref[0]).tobytes() != np.float64(new[0]).tobytes() or \
                        np.float64(ref[1]).tobytes() != np.float64(new[1]).tobytes():
                    bad.append(dict(column=col, case=k, es=es, ref=ref, new=new))
    assert not bad, json.dumps(bad[:5], default=str)


# 4 ------------------------------------------- corrected credit in the lab engine
# Faithful ports of the reviewed upstream fixed_credit_contract_v1/test_contract.py,
# run against refactor_lab.engine (compiled unless NUMBA_DISABLE_JIT=1).
def _engine():
    from refactor_lab.engine import kernels as k, solver as s
    return k, s


def test_engine_credit_renter_floor():
    _, s = _engine()
    P = SimpleNamespace(J=17, debt_taper_weights=np.linspace(1.0, 0.0, 18), debt_caps=np.linspace(5.0, 0.0, 18),
                        unsecured_credit_limit=None, use_age_survival=False)
    current = np.array([-4.0, -1.0, 2.0])
    for j in (0, 6, 11):
        expected = np.minimum(P.debt_taper_weights[j + 1] * np.minimum(current, 0.0), -P.debt_caps[j + 1])
        np.testing.assert_array_equal(s.renter_borrowing_floor(P, current, j), expected)  # None = legacy identity
    for D in (0.0, 2.5):
        P.unsecured_credit_limit = D
        for j in (0, 6, 11):
            np.testing.assert_array_equal(s.renter_borrowing_floor(P, current, j), np.full_like(current, -D))
    P.unsecured_credit_limit = 2.5
    assert np.all(s.renter_borrowing_floor(P, current, 16) == 0.0)      # terminal estate restriction
    P.use_age_survival = True
    P.survival_probs = np.ones(17); P.survival_probs[6] = 0.999
    assert np.all(s.renter_borrowing_floor(P, current, 6) == 0.0)       # positive death probability


def test_engine_credit_binding_contract():
    for bad in (-0.01, float("inf"), float("nan"), [1.0]):
        with pytest.raises(ValueError):
            credit.bind_engine_credit(SimpleNamespace(), "corrected", bad)
    assert credit.bind_engine_credit(SimpleNamespace(), "corrected", 0.0).unsecured_credit_limit == 0.0
    with pytest.raises(ValueError):
        credit.bind_engine_credit(SimpleNamespace(native_solvency_credit=True), "corrected", 1.0)
    with pytest.raises(ValueError):
        credit.bind_engine_credit(SimpleNamespace(unsecured_credit_limit=0.0), "reference")
    with pytest.raises(ValueError):
        credit.bind_engine_credit(SimpleNamespace(), "corrected", None)
    assert not hasattr(credit.bind_engine_credit(SimpleNamespace(), "reference"), "unsecured_credit_limit")


def test_engine_raw_sale_gate():
    k, _ = _engine()
    b = np.array([-1.0, 0.0, 1.0])
    V = np.zeros((3, 2, 1, 1, 1))
    heq = np.array([[0.0, 0.5]])  # raw owner-sale balances -0.5, 0.5, 1.5
    h = np.zeros((1, 2)); dp = np.zeros((1, 2, 1, 1)); bm = np.full((1, 2, 1, 1), -99.0)
    birth = np.zeros((1, 1, 2, 2), dtype=np.bool_); grant = np.zeros((1, 2, 1, 1))
    _, off = k.tenure_choice_kernel(V, b, heq, h, dp, bm, birth, grant, V, False, False, False)
    _, on = k.tenure_choice_kernel(V, b, heq, h, dp, bm, birth, grant, V, False, False, True)
    assert off[0, 1, 0, 0, 0] == 0 and on[0, 1, 0, 0, 0] != 0
    _, _, prob_on = k.tenure_logit_kernel(V, b, heq, h, dp, bm, birth, grant, 0.1, V, False, True)
    assert prob_on[0, 1, 0, 0, 0, 0] == 0.0
    edge = np.array([[0.0, 1.0]])   # raw balances -1, 0, +1: equality stays feasible
    _, choice_edge = k.tenure_choice_kernel(V, b, edge, h, dp, bm, birth, grant, V, False, False, True)
    _, _, prob_edge = k.tenure_logit_kernel(V, b, edge, h, dp, bm, birth, grant, 0.1, V, False, True)
    assert choice_edge[0, 1, 0, 0, 0] == 0 and prob_edge[0, 1, 0, 0, 0, 0] > 0.0
    if os.environ.get("NUMBA_DISABLE_JIT", "0") != "1":
        assert k.NUMBA_AVAILABLE and k.tenure_choice_kernel.signatures and k.tenure_logit_kernel.signatures


def test_engine_renter_fixed_floor_replaces_legacy_line():
    k, _ = _engine()
    grid = np.array([-2.0, 0.0, 2.0])
    args = (np.array([2.0, 2.0, 2.0]), np.array([2.0, 2.0, 2.0]), np.array([[100.0], [0.0], [-100.0]]), np.zeros((3, 1)), 0, grid,
            np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.5]), np.array([1.0]),
            1.0, 10.0, 1e-6, 0.0, 0.0, 0.5, 0.5, 0.95, 0.0, 0.0, 0.381966, 0.618034, 1e-5)
    legacy = k.full_renter_block_kernel(*args)[1]
    positive = k.full_renter_block_kernel(*args, fixed_renter_floor=-1.0)[1]
    zero = k.full_renter_block_kernel(*args, fixed_renter_floor=0.0)[1]
    assert legacy[1, 0] >= 0.0 and positive[1, 0] < -0.5 and zero[1, 0] >= 0.0
    if os.environ.get("NUMBA_DISABLE_JIT", "0") != "1":
        assert k.full_renter_block_kernel.signatures


def test_engine_owner_kernels_unchanged_by_credit_patch():
    """Buyer/incumbent owner kernel bytes equal the certified pass-2 extraction."""
    lab = Path(__file__).resolve().parents[1]
    now = json.loads((lab / "engine/materialize_receipt.json").read_text())
    before = json.loads((lab / "verification" / "engine_receipt_pass2_reference.json").read_text())
    def seg(receipt, module, name):
        return [d["sha256"] for d in receipt["modules"][module]["kept"] if name in d["names"]]
    for name in ("full_owner_block_kernel", "native_due_owner_floor", "native_due_death_floor",
                 "owner_borrowing_floor", "effective_owner_collateral_floor"):
        module = "kernels" if name.endswith("kernel") else "solver"
        assert seg(now, module, name) == seg(before, module, name) != [], name
    assert set(now["overlay"]) == {"parameters", "solver", "kernels"}


# 5 ------------------------------------ maintained optimized saving kernel
# Default: the live engine's promoted indexed kernel. Optional external stage:
# REFACTOR_INDEXED_STAGE=<make_indexed_stage.py output> and REFACTOR_INDEXED_MODULE.
@pytest.mark.parametrize("owner", [False, True])
def test_transformed_stage_kernel_bit_identical(owner):
    """Maintained engine's optimized kernel (default) or an external stage
    (REFACTOR_INDEXED_STAGE + REFACTOR_INDEXED_MODULE) vs the ORIGINAL scalar oracle."""
    import importlib
    engine = Path(__file__).resolve().parents[1] / "engine"
    stage = Path(os.environ.get("REFACTOR_INDEXED_STAGE", engine))
    receipt = json.loads((stage / "transform_receipt.json").read_text())
    assert receipt["status"] == "experimental_until_113_path_fixed_price_certificate"   # receipt as generated
    assert sum(1 for f in receipt["files"].values() if not f["unchanged"]) == 1   # only kernels.py
    os.environ.setdefault("REFACTOR_INDEXED_MODULE", "refactor_lab.engine.kernels")
    idx = importlib.import_module(os.environ["REFACTOR_INDEXED_MODULE"])
    base = oracle_kernels   # original scalar algorithm (frozen package)
    bundle, pin = _require_export()
    bg = inputs.load_inputs(bundle, ROOT, pin).b_grid
    bad = []
    for col, V in enumerate(_reference_columns()):
        if not np.isfinite(V).all():
            continue
        for k, c in enumerate(_boundary_cases(bg, V)):
            for es in (1.0, 1.37):
                a = _call(base.exhaustive_saving_scalar, bg, c, owner, es=es)
                b = _call(idx.exhaustive_saving_scalar, bg, c, owner, es=es)
                if np.float64(a[0]).tobytes() != np.float64(b[0]).tobytes() or \
                        np.float64(a[1]).tobytes() != np.float64(b[1]).tobytes():
                    bad.append(dict(column=col, case=k, es=es, ref=a, new=b))
    assert not bad, json.dumps(bad[:5], default=str)
