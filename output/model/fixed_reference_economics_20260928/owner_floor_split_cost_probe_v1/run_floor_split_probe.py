"""Owner-ladder floor + split transaction cost probe (experimental; NOT a calibration).

Reference: chain54 quarter point (price 0.6744838540900874, H0 7.288573389887633),
quarter-saving purchase rule, no GE root, no recalibration. One lifecycle solve
per case. Adapted from tenure_barrier_probe_v1/run_barrier_probe.py (verified).

Usage: run_floor_split_probe.py A   -> Part A: 12 owner-ladder solves (orig engine)
       run_floor_split_probe.py B   -> Part B: 4 split-cost + 2 identity solves
                                       (engine_split copy with buyer cost psi_buy)
       run_floor_split_probe.py <label> -> single-label rerun (engine by prefix).

Part A: owner grids A1=(2,4,6,8,10) control, A2=(4,6,8,10), A3=(6,8,10);
  x financed share {0.8, 1.0} x need {need0, need1}. Selling cost unchanged.
Part B (engine_split): psi=0.03 seller + psi_buy=0.03 buyer, need {need0, need1}
  x phi {0.8, 1.0}; plus psi_buy=0 identity checks (need0, both phis) that must
  reproduce A1 bit-for-bit.

Writes ONLY inside this folder. Local one-core execution.
Env threads are pinned by the launcher.
"""
from __future__ import annotations

import ast
import copy
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time
import traceback
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
PACKET = ROOT / "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1"
OLD = ROOT / "output/model/fixed_reference_economics_20260928/utility_floor_round2_v1"
ENGINE_ORIG = PACKET / "engines" / "quarter"  # local chain54 quarter engine
ENGINE_SPLIT = HERE / "engine_split"  # copy + buyer-cost patch (Part B only)
LOCAL_RT = PACKET / "local_runtime"
# TASK-named engine; verified byte-identical (see RECEIPT).
TASK_ENGINE = ROOT / "code/model/experiments/quarter_saving_solvency/source"
CHAIN54 = ROOT / ("output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/"
                 "local_runtime/runs/local10_v1/chain54")

PRICE = 0.6744838540900874
H0FIX = 7.288573389887633
LANE = "floor_s0"
ARM = "floor"
TOTAL_LIMIT_A = 5400.0  # 90 minutes for 12 solves + observers + plots
TOTAL_LIMIT_B = 3600.0  # 60 minutes for 6 solves + observers + plots

NEED0_JUMP = 2.49169824815624
NEED0_ROOMS = 0.0
NEED1_JUMP = 1.49169824815624
NEED1_ROOMS = 1.0

GRID_A1 = (2.0, 4.0, 6.0, 8.0, 10.0)  # control (= reference grid)
GRID_A2 = (4.0, 6.0, 8.0, 10.0)
GRID_A3 = (6.0, 8.0, 10.0)
# Each case: label, hbar_child_rooms (0.0 need0 / 1.0 need1), phi,
# psi_override (None = 0.06 reference), psi_buy (buyer cost rate),
# H_own override (None = reference [2,4,6,8,10]).
CASES_A = (
    {"label": "A1_need0_phi08", "hbar1": 0.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
    {"label": "A1_need0_phi10", "hbar1": 0.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
    {"label": "A2_need0_phi08", "hbar1": 0.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A2},
    {"label": "A2_need0_phi10", "hbar1": 0.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A2},
    {"label": "A3_need0_phi08", "hbar1": 0.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A3},
    {"label": "A3_need0_phi10", "hbar1": 0.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A3},
    {"label": "A1_need1_phi08", "hbar1": 1.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
    {"label": "A1_need1_phi10", "hbar1": 1.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
    {"label": "A2_need1_phi08", "hbar1": 1.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A2},
    {"label": "A2_need1_phi10", "hbar1": 1.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A2},
    {"label": "A3_need1_phi08", "hbar1": 1.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A3},
    {"label": "A3_need1_phi10", "hbar1": 1.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A3},
)
CASES_B = (
    {"label": "S_need0_phi08", "hbar1": 0.0, "phi": 0.8, "psi": 0.03, "psi_buy": 0.03, "H_own": GRID_A1},
    {"label": "S_need0_phi10", "hbar1": 0.0, "phi": 1.0, "psi": 0.03, "psi_buy": 0.03, "H_own": GRID_A1},
    {"label": "S_need1_phi08", "hbar1": 1.0, "phi": 0.8, "psi": 0.03, "psi_buy": 0.03, "H_own": GRID_A1},
    {"label": "S_need1_phi10", "hbar1": 1.0, "phi": 1.0, "psi": 0.03, "psi_buy": 0.03, "H_own": GRID_A1},
    {"label": "ID_need0_phi08", "hbar1": 0.0, "phi": 0.8, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
    {"label": "ID_need0_phi10", "hbar1": 0.0, "phi": 1.0, "psi": None, "psi_buy": 0.0, "H_own": GRID_A1},
)
CASES_ALL = {c["label"]: (c, "A") for c in CASES_A}
CASES_ALL.update({c["label"]: (c, "B") for c in CASES_B})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def write_json(path, val):
    path = Path(path)
    path.write_text(json.dumps(val, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_table(path):
    with Path(path).open(newline="") as s:
        return list(csv.DictReader(s))


def write_table(path, rows):
    with Path(path).open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def jsonable(v):
    import math as _m
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, float) and not _m.isfinite(v):
        return repr(v)
    return v


def install_frozen_overlay():
    """Replica of local_runtime/bootstrap.py read-only overlay.

    Local chains 48-57 (incl. chain54) redirect reads of the two drifted
    code/model/tools files to local_runtime/frozen_sources/. Without this,
    the frozen-integration authenticate step fails on source pins.
    """
    import builtins
    import io
    import importlib.machinery
    digests = {"e5f_exact_policy_cache.py":
               "d51bbd13026288db0194b44b420bb49ff6970588f40a906242978d49038a5d6f",
               "test_e5f_exact_policy_cache.py":
               "50784af81c5e6fd64d345d504c8c7fa81209b05a6a655789a224022485208ae5"}
    mapping = {str(ROOT / "code/model/tools" / n): str(LOCAL_RT / "frozen_sources" / n)
               for n in digests}
    for name, digest in digests.items():
        assert hashlib.sha256((LOCAL_RT / "frozen_sources" / name).read_bytes()
                              ).hexdigest() == digest, name

    def mapped(path):
        try:
            return mapping.get(os.path.abspath(os.fspath(path)), path)
        except TypeError:
            return path

    old_io_open = io.open
    old_open = builtins.open

    def redirected_open(path, mode="r", *args, **kwargs):
        target = mapped(path)
        if target != path and any(k in mode for k in "wax+"):
            raise RuntimeError("Frozen overlay write forbidden")
        return old_io_open(target, mode, *args, **kwargs)

    def redirected_builtin(path, mode="r", *args, **kwargs):
        target = mapped(path)
        if target != path and any(k in mode for k in "wax+"):
            raise RuntimeError("Frozen overlay write forbidden")
        return old_open(target, mode, *args, **kwargs)

    io.open = redirected_open
    builtins.open = redirected_builtin
    Path.open = lambda self, mode="r", *args, **kwargs: redirected_open(
        self, mode, *args, **kwargs)
    old_code = importlib.machinery.SourceFileLoader.get_code

    def source_code(loader, fullname):
        target = mapped(loader.path)
        if target != loader.path:
            return loader.source_to_code(
                old_io_open(target, "rb").read(), loader.path)
        return old_code(loader, fullname)

    importlib.machinery.SourceFileLoader.get_code = source_code


def setup_imports(out, engine_root, check_engine_pins):
    install_frozen_overlay()
    import numpy as _np
    import pathlib as _pl
    import importlib as _il
    sys.modules.setdefault("pathlib._local", _pl)
    sys.modules.setdefault("numpy._core", _np.core)
    sys.modules.setdefault("numpy._core.multiarray",
                           _il.import_module("numpy.core.multiarray"))
    sys.modules.setdefault("numpy._core.numeric",
                           _il.import_module("numpy.core.numeric"))
    sys.path.insert(0, str(PACKET))
    sys.path.insert(0, str(engine_root))
    from small_credit_lab.engine import solver as _executed  # noqa: F401
    if check_engine_pins:
        require("engine_split" not in str(_executed.__file__), _executed.__file__)
    else:
        require("engine_split" in str(_executed.__file__), _executed.__file__)
    from refactor_lab.engine import solver as _checked  # noqa: F401
    sys.path.insert(0, str(OLD))
    import runner as native  # noqa: E402
    native.verify_sources()
    for rel, digest in json.loads((PACKET / "source_pins.json").read_text()).items():
        assert sha(ROOT / rel) == digest, rel
    if check_engine_pins:
        for rel, digest in json.loads((PACKET / "engine_pins.json").read_text()).items():
            assert sha(PACKET / rel) == digest, rel
    else:
        # Split engine: every file except the two patched transaction-cost
        # files must match the executed engine byte-for-byte.
        patched = {"small_credit_lab/engine/household.py",
                   "small_credit_lab/engine/kernels.py"}
        for p in sorted((engine_root / "small_credit_lab").rglob("*.py")):
            rel = p.relative_to(engine_root).as_posix()
            if rel in patched or "__pycache__" in rel:
                continue
            assert sha(p) == sha(ENGINE_ORIG / rel), rel
        for p in sorted((engine_root / "refactor_lab").rglob("*.py")):
            rel = p.relative_to(engine_root).as_posix()
            assert sha(p) == sha(ENGINE_ORIG / rel), rel
    for rel, digest in json.loads((PACKET / "manifest.json").read_text())["sha256"].items():
        assert sha(ROOT / rel) == digest, rel
    # TASK-named engine byte-identity (informational; executed engine is above).
    # shared.py/solver.py are unpatched in the split copy, so these hold there too.
    assert sha(TASK_ENGINE / "small_credit_lab/engine/shared.py") == sha(
        engine_root / "small_credit_lab/engine/shared.py")
    assert sha(TASK_ENGINE / "small_credit_lab/engine/solver.py") == sha(
        engine_root / "small_credit_lab/engine/solver.py")
    return native


def build_base(native, out):
    import_runner_cfg = json.loads((PACKET / "plan.json").read_text())
    chain_completed = json.loads((CHAIN54 / "search/completed.json").read_text())
    sel = chain_completed["selected"]
    require(abs(float(sel["price"]) - PRICE) < 1e-12, "Reference price drift")
    require(abs(float(sel["H0_derived"]) - H0FIX) < 1e-12, "Reference H0 drift")
    point = sel["parameters"]
    if isinstance(point, str):
        point = ast.literal_eval(point)
    require(set(point) == {"beta_annual", "chi", "child_benefit_curvature",
                           "first_birth_fixed_cost", "h_P", "kappa_fert",
                           "kappa_fert_continuation", "psi_child",
                           "tenure_choice_kappa", "theta0"}, "Reference point keys differ")
    seed, bounds, _ = native.inputs.seed_and_bounds(LANE)
    bounds = {k: tuple(x) for k, x in bounds.items()}
    bounds["psi_child"] = tuple(import_runner_cfg["psi_bounds"])
    bounds["h_P"] = (0.1, 2.6)
    native.inputs.LANES[LANE].update(seed=dict(point), bounds=bounds,
                                     free_coordinates=list(point))
    P, grid = native.inputs.proposal(LANE)
    P, entry = native.inputs.entry(P, grid, "nonnegative_mean")
    require(np.asarray(P.phi).shape == (4,), "Expected four financed-share states")
    P.phi = np.full_like(np.asarray(P.phi, dtype=float), 0.8)
    P.experimental_purchase_saving_fraction = 0.25
    from small_credit_lab.engine import solver as _solv
    _ = _solv  # bound after engine import
    Q = native.utility_checks(P, grid, LANE, out)
    P_ref = native.inputs.bind(P, point, bounds, ARM)
    return import_runner_cfg, point, bounds, entry, P_ref, grid, Q


def build_context(native, base, ge, P_ref, grid, point, out):
    native.install_reporter_on_authored(base.authored)
    ctx = base.authored.context_from_bundle(SimpleNamespace(
        bundle=ROOT / "output/model/publication_refactor_20260929/local_export_v1/inputs",
        reference_root=ROOT, out=out))
    base.authored.authenticate_frozen(ctx)
    from small_credit_lab.engine.shared import annual_gross_income_at_state
    from small_credit_lab.engine import solver
    sd = solver.precompute_shared(P_ref, grid)
    require(float(sd.cb_flat[0, 0]) == float(sd.hb_flat[0, 0]) == float(sd.gb_flat[0, 0]) == 0.0,
            "Necessary cash preflight missing childless floors")
    actual_income = np.asarray([annual_gross_income_at_state(P_ref, 0, 0, float(z))
                                for z in P_ref.z_grid])
    np.testing.assert_array_equal(
        actual_income, P_ref.income[0, 0] * P_ref.z_grid / P_ref.period_years / (1 - P_ref.tau_pay))
    cal = ctx["prepared"].rt["primitive"].pf.calendar
    cohort = cal.entrant_cohort(np.asarray([1.0]), P_ref, grid)
    np.testing.assert_allclose(
        cohort.sum(axis=(1, 2, 4, 5)),
        P_ref.fixed_reference_entry_conditional * P_ref.z_weights[None, :], rtol=0, atol=2e-16)
    require(abs(cohort.sum() - 1) < 2e-12, "Entrant mass changed")
    seed, bounds, _ = native.inputs.seed_and_bounds(LANE)
    for row in ctx["manifest"]["full_parameter_table"]:
        if row["parameter"] in bounds:
            row["lower"], row["upper"] = map(str, bounds[row["parameter"]])
    from small_credit_lab import credit
    credit.bind_engine_credit(P_ref, "corrected", 0.0)
    actual = ctx["fp"].actual_parameters(ctx["prepared"], P_ref, grid)
    expected = native.expected_parameters(point, (120, 9), ARM)
    expected.update(H0=H0FIX, financed_share=0.8)
    ge.validate_parameter_estimates(dict(expected_parameters=expected),
                                    ctx["manifest"]["full_parameter_table"], actual)
    live_sd = solver.precompute_shared(P_ref, grid)
    sys.path.insert(0, str(ROOT / "code/model"))
    from refactor_lab.engine import solver as checked_solver
    checked_sd = checked_solver.precompute_shared(P_ref, grid)
    for key in ("h_bar", "c_bar", "g_bar", "alpha_flat", "psi_v", "escale_flat"):
        np.testing.assert_array_equal(getattr(live_sd, key), getattr(checked_sd, key))
    require(float(live_sd.h_bar[1, 1]) == expected["h_P"]
            and float(live_sd.h_bar[1, 0]) == 0.0, "Executed physical floor differs")
    ctx.update(P=P_ref, b_grid=grid, selected_d_bar=0.0,
               reference_psi=float(P_ref.psi_child),
               expected_dimensions={"wealth_grid_nodes": len(grid),
                                    "income_states": len(P_ref.z_grid)},
               free_coordinates=list(point), price_start=PRICE,
               phase_b_max_new_lifecycle=32)
    return ctx


def observe_fixed_price(ge, candidate, live, name, outdir, expected, final_observe=True):
    """Fixed-price, fixed-H0 observer modeled on the pipeline observe_price.

    Identical lifecycle/observer calls; differences: no H0 derivation (P.H0
    stays fixed), no renewal/PAYGO hard gates (measured as diagnostics).
    """
    t0 = time.time()
    fp, prepared, manifest, objective, runtime, reference = ge._observer_context(candidate)
    cal = prepared.rt["primitive"].pf.calendar
    P, grid, sd, sol = (live[k] for k in ("P", "b_grid", "sd", "sol"))
    price = np.asarray(live["price"], dtype=float).reshape(-1)
    ge._require(price.size == 1 and np.isfinite(price[0]) and price[0] > 0,
                "Invalid scalar price")
    ge._require(float(P.unsecured_credit_limit) == float(candidate["selected_d_bar"]),
                "Credit limit drift")
    ge._require(float(P.psi_child) == float(candidate["reference_psi"]),
                "Child benefit drift")
    ge._require(np.array_equal(grid, candidate["b_grid"]), "Wealth grid drift")
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, recon = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(recon["stationary_post_fertility_nesting_l1"], 5e-9,
                             "Stationary cohort reconstruction")
    runtime.require_abs_gate(recon["stationary_feasibility_projection_mass"], 0.0,
                             "Stationary feasibility projection")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]),
        float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
                             supply_rule=supply, supplied_policy=policy)
    renter_floor = ge.audit_realized_renter_floor(P, grid, policy, ev)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=pre, supply_rule=supply,
                  demographic_seed=reference.get("demographic_seed"))
    out = Path(candidate["out"]) / "phase_b_ge" / name
    out.mkdir(parents=True, exist_ok=True)
    try:
        gates = fp.gates(packet, prepared, out, stationary=True)
        gates_status = "passed"
    except Exception as exc:  # fixed-price diagnostic: record, do not stop
        gates = {"status": "failed_fixed_price_diagnostic",
                 "error_type": type(exc).__name__, "error": str(exc)}
        fp.write(out / "gates.json", jsonable(gates))
        gates_status = "failed_fixed_price_diagnostic"
    fiscal = None
    paygo = float("nan")
    try:
        fiscal = gates.get("fiscal_certificate") if isinstance(gates, dict) else None
        if isinstance(fiscal, dict):
            paygo = float(fiscal["actual_accounts"]["scaled_pension_budget_residual"])
    except Exception:
        paygo = float("nan")
    native = prepared.rt["primitive"].pf.transition
    actual_births = native.calendar_topcode_birth_accounting(
        ev.g_pre, ev.g_post_fertility, float(ev.births), P)["topcode_adjusted_birth_children"]
    entry = float(sol.entry_rate)
    demand = float(np.asarray(ev.demand_by_loc).sum())
    physical_supply = float(np.asarray(ev.supply_by_loc).sum())
    require(all(map(np.isfinite, (actual_births, entry, demand, physical_supply))),
            "Nonfinite fixed-price accounting")
    renewal = float(actual_births / (2.1 * entry) - 1)
    closure = dict(price=float(price[0]), d_bar=float(P.unsecured_credit_limit),
        adjusted_births_per_normalized_household=float(actual_births),
        actual_entry_per_normalized_household=entry,
        renewal_residual=renewal, renewal_status="measured_fixed_price_diagnostic",
        population_scale=1.0,
        normalized_housing_demand=demand, physical_housing_supply=physical_supply,
        absolute_housing_demand=demand,
        absolute_housing_residual=demand - physical_supply,
        housing_market_status="measured_fixed_H0_diagnostic",
        actual_paygo_residual=paygo, outside_entry=0.0,
        birth_to_entry_conversion=1 / 2.1,
        housing_supply_elasticity=float(P.xi_supply[0]),
        occupied_renter_floor=jsonable(renter_floor), normalized_population=1.0,
        H0_fixed=float(P.H0[0]), H0_bounds=[0.2, 80.0],
        housing_supply_coefficient_role="fixed at selected-point derived value for diagnostic",
        gates_status=gates_status,
        fixed_price=PRICE, no_renewal_or_housing_root=True)
    fp.write(out / "closure.json", closure)
    fertility = {k: prepared.rt["observe_initial_fertility"](ev, P, age_projection=k)
                 for k in ("uniform_birth_time", "constant_post_cell")}
    housing = prepared.rt["observe_initial_housing_wealth"](
        ev, P, grid, sd, diagnostic_enabled=True,
        age_projection="uniform_within_age_cell", diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = prepared.rt["observe_recent_parent_flow"](
        ev, P, diagnostic_enabled=True, snapshot=prepared.rt["SNAPSHOT"],
        age_projection=prepared.rt["AGE_PROJECTION"], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=name, fixed_price_diagnostic=True))
    tfr = float(prepared.rt["chain"].extract_moments(sol, P)["tfr"])
    fits = runtime.score_targets(objective, fertility, housing, recent["model_value"], tfr)
    total_loss = 0.0
    for row in fits:
        try:
            lc = float(row.get("loss_contribution"))
        except (TypeError, ValueError):
            continue
        if np.isfinite(lc) and str(row.get("role", "scored")) == "scored":
            total_loss += lc
    params = [dict(row) for row in manifest["full_parameter_table"]]
    actual_params, hbar_reporting = actual_parameters_report(fp, prepared, P, grid)
    require(len(fits) == 14 and len(params) == 31, "14 fit/31 parameter rows required")
    for row in params:
        key = row["parameter"]
        row["reference_estimate"], row["estimate"] = (
            row["estimate"], str(float(actual_params[key])))
        if key in candidate.get("free_coordinates", []):
            row["status"] = "reference free coordinate; fixed at selected estimate for diagnostic"
        elif key == "H0":
            row["status"] = "fixed at selected-point derived value for diagnostic"
            row["lower"], row["upper"] = "0.2", "80.0"
            row["near_bound"] = str(min(float(P.H0[0]) - 0.2,
                                       80.0 - float(P.H0[0])) <= 0.01 * 79.8)
        elif key == "financed_share":
            row["status"] = "experimental fixed policy input"
        elif key == "selling_cost" and (candidate.get("barrier_overrides") or {}).get("psi") is not None:
            row["status"] = ("experimental override: selling cost %s" %
                             (candidate["barrier_overrides"]["psi"],))
        elif key == "tenure_choice_kappa" and (candidate.get("barrier_overrides") or {}).get("kappa") is not None:
            row["status"] = "experimental override of reference free coordinate"
        elif key == "h_P":
            row["status"] = ("reference h_P; child-room floor variant applied "
                             "via hbar fields" if float(actual_params[key]) != float(expected.get("h_P", 0))
                             else "reference h_P")
        elif key == "psi_child":
            row["status"] = "fixed reference benefit for diagnostic"
        elif key == "child_benefit_CRRA_coefficient":
            row["status"] = "derived from fixed benefit and proposed curvature"
    fp.table(out / "target_fit.csv", fits)
    fp.table(out / "parameters.csv", params)
    fp.write(out / "observers.json", cal.jsonable(
        dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
    report_packet = dict(packet)
    report_ev = copy.copy(ev)
    report_ev.supply_by_loc = np.asarray(ev.supply_by_loc)
    report_packet["evaluation"] = report_ev
    prepared.rt["audit"].standard_diagnostics(
        report_packet, out, validate_production_young=False)
    plots = sorted(p.name for p in (out / "standard_diagnostics").glob("*.png"))
    closure.update(standard_plot_count=len(plots),
                   target_fit_rows=len(fits), parameter_rows=len(params),
                   completed_fertility_tfr=tfr,
                   total_target_loss=total_loss,
                   hbar_reporting_accommodation=hbar_reporting,
                   observe_seconds=time.time() - t0)
    fp.write(out / "closure.json", closure)
    extra = extra_moments(P, grid, policy, ev, fertility)
    fp.write(out / "extra_moments.json", jsonable(extra))
    return dict(closure=closure, target_fit=fits, parameters=params,
                fertility=fertility, housing=housing, recent=recent,
                extra=extra, gates=gates, total_target_loss=total_loss)


def actual_parameters_report(fp, prepared, P, grid):
    """fp.actual_parameters with a documented need1 reporting accommodation.

    The frozen adapter hard-requires hbar_child_rooms == 0 (calibration
    contract guard). The underlying native mapping does not depend on that
    field except for the h_P row, which the adapter reports as the first-child
    jump. For need1 arms we evaluate the adapter on an identical copy with
    only hbar_child_rooms zeroed, so every other estimate is exact and h_P
    reports the need1 first-child jump. The true hbar pair is recorded in
    case metadata (hbar_first_child_jump / hbar_child_rooms).
    """
    try:
        return fp.actual_parameters(prepared, P, grid), False
    except RuntimeError as exc:
        if "Unexpected later-child floor" not in str(exc):
            raise
        rep = copy.deepcopy(P)
        rep.hbar_child_rooms = 0.0
        return fp.actual_parameters(prepared, rep, grid), True


def extra_moments(P, grid, policy, ev, fertility):
    g = np.asarray(ev.g_current)
    g_pre = np.asarray(ev.g_pre)
    g_post = np.asarray(ev.g_post_fertility)
    hR = np.asarray(policy.hR_pol)
    mass = float(g.sum())
    out = {}
    # Birth flows by order: sums over parity_birth_flows_by_age columns.
    for key in ("uniform_birth_time", "constant_post_cell"):
        acct = fertility[key]["accounting"]
        flows = np.asarray(acct["parity_birth_flows_by_age"], dtype=float)
        out["birth_flows_" + key] = {
            "first": float(flows[:, 0].sum()),
            "second": float(flows[:, 1].sum()),
            "third_bin": float(flows[:, 2].sum()),
            "first_birth_flow": float(acct["first_birth_flow"]),
            "explicit_birth_flow": float(acct["explicit_birth_flow"]),
        }
    # Cross-check from pre/post distributions (parity axis -2, child axis -1).
    by_order = [float(g_post[..., n:, :].sum() - g_pre[..., n:, :].sum())
                for n in range(1, P.n_parity)]
    out["birth_flows_g_arrays"] = {
        "first": by_order[0], "second": by_order[1], "third_bin": by_order[2]}
    # Renter cap shares by children at home (child-state axis -1), renters t=0.
    cap = float(P.hR_max)
    rent_mask = hR >= cap - 1e-9
    ncs = g.shape[-1]
    cap_by_m = {}
    for m in range(ncs):
        gm = g[..., m]
        rm = gm[:, 0].sum()
        tot = gm.sum()
        cap_by_m[str(m)] = {
            "renter_mass": float(rm), "total_mass": float(tot),
            "cap_share_of_renters": float((gm[:, 0] * rent_mask[..., m][:, 0]).sum() / rm) if rm > 0 else float("nan"),
            "cap_share_of_all": float((gm[:, 0] * rent_mask[..., m][:, 0]).sum() / tot) if tot > 0 else float("nan"),
        }
    out["renter_cap_by_children_home"] = cap_by_m
    out["rental_cap"] = cap
    # Matched-entrant cell: youngest age cell j=0 (ages 18-21); entrants fixed.
    j0 = 0
    pre0 = g_pre[:, :, :, j0]
    post0 = g_post[:, :, :, j0]
    cur0 = g[:, :, :, j0]
    entrant_mass = float(pre0.sum())
    entrant_parity_pos = float(pre0[..., 1:, :].sum())
    first0 = float(post0[..., 1:, :].sum() - pre0[..., 1:, :].sum())
    own0 = float(cur0[:, 1:].sum())
    out["entrant_cell_j0"] = {
        "entrant_mass": entrant_mass,
        "entrant_nonzero_parity_mass": entrant_parity_pos,
        "first_birth_probability": first0 / entrant_mass if entrant_mass > 0 else float("nan"),
        "ownership_rate": own0 / float(cur0.sum()) if float(cur0.sum()) > 0 else float("nan"),
    }
    # Ownership / rooms aggregates.
    out["ownership_rate_all"] = float(g[:, 1:].sum()) / mass
    # Ownership of parents by children at home (child-state axis -1).
    own_by_m = {}
    for m in range(ncs):
        gm = g[..., m]
        tot = float(gm.sum())
        own_by_m[str(m)] = {
            "total_mass": tot,
            "ownership_rate": float(gm[:, 1:].sum() / tot) if tot > 0 else float("nan"),
        }
    out["ownership_by_children_home"] = own_by_m
    # Childless (parity 0 = zero children ever born) owners aged 22-29 by rung.
    # Age cells selected by overlap of [age_start + j*da, age_start + (j+1)*da)
    # with [22, 30); tenure axis 1, age axis 3, parity axis -2.
    J = int(P.J)
    a0 = float(P.age_start)
    da = float(P.da)
    require(g.shape[3] == J and g.shape[-2] == int(P.n_parity),
            "g_current axes differ from (b, tenure, loc, age, ..., parity, child)")
    js = [j for j in range(J) if (a0 + (j + 1) * da) > 22.0 and (a0 + j * da) < 30.0]
    sub = np.take(np.asarray(g, dtype=float), js, axis=3)
    tot29 = float(sub.sum())
    par0 = np.take(sub, [0], axis=-2)
    Hown = [float(x) for x in np.asarray(P.H_own, dtype=float)]
    by_rung = {}
    for ti, h in enumerate(Hown, start=1):
        m = float(np.take(par0, [ti], axis=1).sum())
        by_rung[repr(h)] = {
            "rooms": h,
            "mass": m,
            "share_of_22_29": (m / tot29) if tot29 > 0 else float("nan"),
        }
    out["childless_owners_22_29"] = {
        "definition": "parity-0 owners; shares of the age 22-29 cross-section mass",
        "age_cells": js,
        "age_lo": [a0 + j * da for j in js],
        "age_hi": [a0 + (j + 1) * da for j in js],
        "total_mass_22_29": tot29,
        "by_rung_rooms": by_rung,
    }
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "A"
    t_start = time.time()
    if mode in ("A", "B"):
        part = mode
        cases = CASES_A if part == "A" else CASES_B
        labels = [c["label"] for c in cases]
        total_limit = TOTAL_LIMIT_A if part == "A" else TOTAL_LIMIT_B
        done_name = "completed_%s.json" % part
    elif mode in CASES_ALL:
        case0, part = CASES_ALL[mode]
        cases = (case0,)
        labels = [mode]
        total_limit = TOTAL_LIMIT_B if part == "B" else TOTAL_LIMIT_A
        done_name = None
    else:
        raise RuntimeError("unknown mode %r (want A, B, or a case label)" % mode)
    deadline = t_start + total_limit
    engine_root = ENGINE_ORIG if part == "A" else ENGINE_SPLIT
    if done_name is not None:
        require(not (HERE / done_name).exists(),
                "Refusing to overwrite completed part %s" % part)
    native = setup_imports(HERE, engine_root, part == "A")
    out_pre = HERE / ("setup" if done_name is not None else "setup_tmp_%s" % part)
    import shutil as _sh
    if out_pre.exists() and done_name is not None:
        _sh.rmtree(out_pre)
    out_pre.mkdir(parents=True, exist_ok=True)
    CONFIG, point, bounds, entry_report, P_ref, grid, _ = build_base(native, out_pre)
    # Fixed-H0 diagnostic: pin the supply coefficient before validation/solve.
    # H0 is accounting-only (households take price as given).
    P_ref.H0 = np.full_like(np.asarray(P_ref.H0, dtype=float), H0FIX)
    write_json(HERE / "reference_point.json",
               dict(point=point, price=PRICE, H0=H0FIX, entry=entry_report))
    # nk semantics + current floor values.
    from small_credit_lab.engine import solver as eng_solver
    from small_credit_lab.engine import parameters as eng_par
    sd_ref = eng_solver.precompute_shared(P_ref, grid)
    nk_info = dict(
        child_state_mode=str(P_ref.child_state_mode),
        independent_count_active=bool(eng_par.independent_child_maturation_active(P_ref)),
        nk_definition=("children currently at home (child-state axis cs; nk=cs when cs<=parity "
                       "else 0) under independent_count; parity (ever born) otherwise"),
        child_room_floor=bool(getattr(P_ref, "child_room_floor", False)),
        hbar_first_child_jump=float(P_ref.hbar_first_child_jump),
        hbar_child_rooms=float(P_ref.hbar_child_rooms),
        h_bar_parity1_cs1=float(sd_ref.h_bar[1, 1]),
        h_bar_childless=float(sd_ref.h_bar[0, 0]),
    )
    write_json(HERE / "floor_semantics.json", nk_info)
    current_floor_1child = float(P_ref.hbar_first_child_jump) + float(P_ref.hbar_child_rooms) * 1.0
    need1_jump = current_floor_1child - 1.0
    require(abs(float(P_ref.hbar_first_child_jump) - NEED0_JUMP) == 0.0
            and float(P_ref.hbar_child_rooms) == NEED0_ROOMS, "need0 != TASK reference floor")
    require(abs(need1_jump - NEED1_JUMP) == 0.0, "need1 != TASK floor")
    require(abs(float(P_ref.psi) - 0.06) < 1e-12, "Reference selling cost drift")
    # Context (frozen integration).
    sys.path.insert(0, str(native.BASE))
    spec = importlib.util.spec_from_file_location(
        "pilot_matched_workflow", native.BASE / "run_comparison.py")
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    sys.path.insert(0, str(native.HERE))
    import phase_b_pilot as ge  # noqa: E402
    if part == "B":
        # Observer/calendar forward maps live in the production solver, which
        # ignores P.psi_buy. Delegate them to the split-engine implementation
        # (identical signature; inert at psi_buy=0, proven below) so the KFE
        # wealth transition in reconstruction/evaluation matches the solve.
        # Runtime delegation only; no repo source file is edited.
        import inspect as _inspect
        import intergen_housing_fertility_optimized.solver as prod_model
        from small_credit_lab.engine.household import (
            build_forward_tenure_transition_maps as split_forward_maps)
        _sig = ["P", "b_grid", "hc", "he", "phi_choice", "birth_dp", "birth_entry_grant"]
        require(list(_inspect.signature(prod_model.build_forward_tenure_transition_maps).parameters)
                == list(_inspect.signature(split_forward_maps).parameters) == _sig,
                "Forward-map signature drift")
        _orig_forward = prod_model.build_forward_tenure_transition_maps
        _nt = 1 + int(P_ref.n_house)
        _hc = np.zeros((int(P_ref.I), _nt))
        _he = np.zeros((int(P_ref.I), _nt))
        for _i in range(int(P_ref.I)):
            for _ten in range(1, _nt):
                _hs = float(np.asarray(P_ref.H_own, dtype=float)[_ten - 1])
                _hc[_i, _ten] = PRICE * _hs
                _he[_i, _ten] = (1.0 - float(P_ref.psi)) * PRICE * _hs
        _sd0 = eng_solver.precompute_shared(P_ref, grid)
        _o_idx, _o_wt = _orig_forward(
            P_ref, grid, _hc, _he, _sd0.phi_choice, _sd0.birth_dp,
            _sd0.birth_entry_grant)
        _s_idx, _s_wt = split_forward_maps(
            P_ref, grid, _hc, _he, _sd0.phi_choice, _sd0.birth_dp,
            _sd0.birth_entry_grant)
        require(bool(getattr(P_ref, "psi_buy", 0.0) == 0.0), "P_ref must not set psi_buy")
        np.testing.assert_array_equal(_s_idx, _o_idx)
        np.testing.assert_array_equal(_s_wt, _o_wt)
        prod_model.build_forward_tenure_transition_maps = split_forward_maps
        write_json(HERE / "engine_split_hook.json", jsonable(dict(
            hook="intergen_housing_fertility_optimized.solver."
                 "build_forward_tenure_transition_maps -> "
                 "engine_split/small_credit_lab/engine/household."
                 "build_forward_tenure_transition_maps",
            reason="production observer/calendar maps ignore P.psi_buy; delegation "
                   "applies the buyer fee identically in solve and observers",
            inert_at_zero_proven=True,
            split_household_sha256=sha(HERE / "engine_split/small_credit_lab/engine/household.py"),
            split_kernels_sha256=sha(HERE / "engine_split/small_credit_lab/engine/kernels.py"))))
    base_ctx = build_context(native, base, ge, P_ref, grid, point, out_pre)
    results = {}
    if (HERE / "case_results.json").exists():
        results = json.loads((HERE / "case_results.json").read_text())
    for case in cases:
        label, hbar1, phi = case["label"], case["hbar1"], case["phi"]
        if label not in labels:
            continue
        if label in results and results[label].get("status") == "passed" and done_name is not None:
            continue
        case_out = HERE / label
        import shutil
        if case_out.exists():
            shutil.rmtree(case_out)
        case_out.mkdir(parents=True)
        P_case = copy.deepcopy(P_ref)
        if hbar1 == 0.0:
            P_case.hbar_child_rooms = 0.0
            # jump stays at reference value
        else:
            P_case.hbar_child_rooms = float(hbar1)
            P_case.hbar_first_child_jump = float(need1_jump)
        P_case.phi = np.full_like(np.asarray(P_case.phi, dtype=float), float(phi))
        P_case.H0 = np.full_like(np.asarray(P_case.H0, dtype=float), H0FIX)
        if case["psi"] is not None:
            P_case.psi = float(case["psi"])  # selling cost override
        P_case.psi_buy = float(case.get("psi_buy", 0.0))  # buyer cost (split engine)
        require(not bool(getattr(P_case, "joint_nested_choice", False)),
                "Joint-nested path bypasses the buyer-cost kernel args")
        if case["H_own"] is not None:
            P_case.H_own = np.asarray(case["H_own"], dtype=float)
            P_case.n_house = len(case["H_own"])
        from small_credit_lab import credit
        credit.bind_engine_credit(P_case, "corrected", 0.0)
        expected = native.expected_parameters(point, (120, 9), ARM)
        expected.update(H0=H0FIX, financed_share=float(phi))
        if case["psi"] is not None:
            expected.update(selling_cost=float(case["psi"]))
        actual_now, hbar_reporting = actual_parameters_report(
            base_ctx["fp"], base_ctx["prepared"], P_case, grid)
        if float(actual_now.get("h_P", 0.0)) != float(expected.get("h_P", 0.0)):
            expected["h_P"] = float(actual_now["h_P"])
        candidate = dict(base_ctx)
        candidate.update(P=P_case, b_grid=grid, out=case_out,
                         selected_d_bar=0.0, reference_psi=float(P_case.psi_child),
                         expected_parameters=expected,
                         barrier_overrides=dict(psi=case["psi"], H_own=list(case["H_own"])
                                                if case["H_own"] is not None else None,
                                                psi_buy=float(case.get("psi_buy", 0.0))),
                         expected_dimensions={"wealth_grid_nodes": len(grid),
                                              "income_states": len(P_case.z_grid)},
                         free_coordinates=[], fixed_coordinates=list(point),
                         fixed_price=PRICE, deadline_epoch=deadline,
                         price_start=PRICE)
        ge.validate_parameter_estimates(candidate,
                                        candidate["manifest"]["full_parameter_table"],
                                        actual_now)
        budget = base.ArmBudget(case_out, deadline)
        budget.max_lifecycle = 1
        solve_label = label + "_selected"
        stage = case_out / "phase_b_ge" / solve_label / "stage"
        s0 = time.monotonic()
        live = ge.solve_fixed_price(candidate, 0.0, PRICE, budget, solve_label, stage)
        solve_seconds = time.monotonic() - s0
        observed = observe_fixed_price(ge, candidate, live, solve_label, case_out,
                                       expected)
        summary = live["summary"]
        results[label] = dict(
            status="passed", phi=float(phi),
            hbar_first_child_jump=float(P_case.hbar_first_child_jump),
            hbar_child_rooms=float(P_case.hbar_child_rooms),
            selling_cost_psi=float(P_case.psi),
            buyer_cost_psi_buy=float(getattr(P_case, "psi_buy", 0.0)),
            engine="split" if part == "B" else "orig",
            tenure_choice_kappa=float(P_case.tenure_choice_kappa),
            H_own=[float(x) for x in np.asarray(P_case.H_own, dtype=float)],
            n_house=int(P_case.n_house),
            hbar_reporting_accommodation=hbar_reporting,
            price=PRICE, H0=H0FIX,
            purchase_saving_fraction=float(P_case.experimental_purchase_saving_fraction),
            lifecycle_solve_seconds=solve_seconds,
            native_array_count=summary.get("native_array_count"),
            adult_entry_relative_gap=summary.get("adult_entry_relative_gap"),
            closure=jsonable(observed["closure"]),
            extra=jsonable(observed["extra"]),
            target_fit=observed["target_fit"],
            total_target_loss=observed["total_target_loss"],
            report=str(case_out / "phase_b_ge" / solve_label))
        write_json(HERE / "case_results.json", jsonable(results))
        write_json(HERE / "latest_completed.json", jsonable(results[label]))
    if done_name is not None and all(results.get(c["label"], {}).get("status") == "passed"
                                     for c in cases):
        write_json(HERE / done_name, jsonable(
            dict(status="part_%s_solves_complete" % part, cases={c["label"]: results[c["label"]] for c in cases},
                 elapsed_seconds=time.time() - t_start)))
    print(json.dumps({k: v.get("status") for k, v in results.items()}, indent=1))


if __name__ == "__main__":
    main()
