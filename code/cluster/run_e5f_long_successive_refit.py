"""Bounded long-horizon successive permanent-preference refit.

This driver deliberately keeps the original four-vintage household birth queue.
Each historical date solves a new permanent-preference terminal, executes only
the first date of that forecast, and checkpoints the native inherited state
before the next surprise is allowed to start.  The manifest owns all paths,
deadlines, hashes, and Slurm follow-on scripts.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from types import SimpleNamespace as NS
from unittest.mock import patch

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

PSI_SS = 0.1489153145785918
YEARS = (2007, 2011, 2015, 2019)
COUNTS = (104, 103, 102, 101)
TARGETS = (1.974875, 1.861000, 1.755375, 1.645750)
SEEDS = (.12891531457859182, .11696608375682901, .10564290922456478, .09221854783921073)
TOL = .005
NATIVE_TOL = 2e-10


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(driver, path, value):
    driver.save(Path(path), value)


def dump_pickle(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wb") as stream:
        pickle.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)


def exact_gap(left, right):
    a, b = np.asarray(left, float), np.asarray(right, float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("nonfinite replay object")
    return float(np.max(np.abs(a - b))) if a.size else 0.


def manifest_contract(m, path):
    required = {"spec", "output", "numerics_dir", "jacobian_source", "empirical_blocks",
        "shock_years", "year_targets", "psi_seeds", "bounds_relative", "fit_tolerance",
        "max_trials", "max_rounds", "total_deadline_unix", "terminal_year", "scripts"}
    missing = required - set(m)
    if missing:
        raise ValueError("manifest misses: " + ", ".join(sorted(missing)))
    for item, digest in m.get("file_sha256", {}).items():
        if sha(item) != digest:
            raise ValueError("pinned input changed: " + item)
    if not Path(m["output"]).is_absolute():
        raise ValueError("output must be absolute")
    if "absolute_deadline_unix" not in m and "total_deadline_unix" not in m:
        raise ValueError("manifest needs absolute_deadline_unix or total_deadline_unix")
    if sha(m["empirical_blocks"]) != m.get("target_sha256"):
        raise ValueError("empirical fertility block fingerprint differs")
    rows = list(csv.DictReader(Path(m["empirical_blocks"]).open()))
    observed = [(int(r["decision_year"]), float(r["period_tfr_arithmetic_mean"])) for r in rows]
    if observed != list(zip(YEARS, TARGETS)) or list(m["shock_years"]) != list(YEARS):
        raise ValueError("dated fertility target contract differs")
    if [float(m["year_targets"][str(y)]) for y in YEARS] != list(TARGETS):
        raise ValueError("manifest target values differ")
    if [float(x) for x in m["psi_seeds"]] != list(SEEDS):
        raise ValueError("manifest preference seeds differ")
    if [float(x) for x in m["bounds_relative"]] != [-.20, .02] or float(m["fit_tolerance"]) != TOL:
        raise ValueError("manifest fit bounds or tolerance differ")
    if int(m["terminal_year"]) != 2423 or [horizon(y, int(m["terminal_year"])) for y in YEARS] != list(COUNTS):
        raise ValueError("manifest terminal boundary differs")
    if m.get("recovery_sources"):
        parent_path = m["recovery_parent_manifest"]
        if parent_path not in m["file_sha256"]:
            raise ValueError("recovery parent manifest must be pinned")
        parent = read(parent_path)
        for key in ("spec", "target_sha256", "initial_score", "initial_target_fingerprint",
                    "shock_years", "year_targets", "psi_seeds", "bounds_relative",
                    "fit_tolerance", "terminal_year", "contract"):
            if m[key] != parent[key]:
                raise ValueError("recovery changed scientific contract: " + key)
        for source in m["recovery_sources"]:
            for key in ("root_receipt", "terminal_pickle", "terminal_receipt", "native_rows"):
                if source[key] not in m["file_sha256"]:
                    raise ValueError("recovery input must be hash pinned: " + key)
    return sha(path)


def horizon(year, terminal_year=2423):
    """Decision dates before the boundary at ``terminal_year``."""
    if (terminal_year-year) % 4 or terminal_year <= year:
        raise ValueError("terminal year must be a later four-year date")
    return (terminal_year-year)//4


def next_psi(trials, seed, bounds, seen=()):
    """One bounded secant/bracket proposal, without accepting invalid trials.

    ``trials`` is a sequence of records with finite ``psi`` and signed ``gap``.
    The caller records failed equilibrium trials but does not pass them here.
    """
    lo, hi = map(float, bounds)
    if lo >= hi:
        raise ValueError("invalid psi bounds")
    used = [float(x) for x in seen]
    valid = [(float(r["psi"]), float(r["gap"])) for r in trials
             if r.get("gap") is not None and math.isfinite(float(r["gap"]))]
    proposal = float(seed)
    brackets = [(a, b) for i, a in enumerate(valid) for b in valid[i+1:]
                if a[1]*b[1] < 0]
    if brackets:
        x0, y0, x1, y1 = min(((a[0], a[1], b[0], b[1]) for a, b in brackets),
                              key=lambda z: abs(z[0]-z[2]))
        if y1 != y0:
            proposal = x1-y1*(x1-x0)/(y1-y0)
        low, high = sorted((x0, x1)); span = high-low
        proposal = min(high-.15*span, max(low+.15*span, proposal))
    elif len(valid) >= 2:
        x0, y0 = valid[-2]; x1, y1 = valid[-1]
        if y1 != y0: proposal = x1-y1*(x1-x0)/(y1-y0)
    elif valid:
        x, gap = valid[-1]
        proposal = x-.01 if gap > 0 else x+.01
    proposal = min(hi, max(lo, proposal))
    if all(abs(proposal-x) > 1e-10 for x in used):
        return proposal
    # A duplicate secant is not a reason to recycle an endpoint: take a small
    # in-bounds bracket step in the direction implied by the last valid gap.
    direction = -1 if valid and valid[-1][1] > 0 else 1
    for step in (.005, .01, .02, .04, .08):
        candidate = min(hi, max(lo, proposal+direction*step))
        if all(abs(candidate-x) > 1e-10 for x in used):
            return candidate
    return None


def load_frozen_context(m):
    """Import runner only after the manifest has validated frozen inputs."""
    spec = read(m["spec"])
    source = Path(spec["batch"]) / "source"
    sys.path.insert(0, str(source))
    runner = importlib.import_module("run_e5f_original_queue_experiments")
    c = runner.load_context(m["spec"])
    c.spec_path = Path(m["spec"])
    numerics = Path(m["numerics_dir"])
    if not numerics.is_dir():
        raise ValueError("manifest numerics_dir is unavailable")
    sys.path.insert(0, str(numerics))
    scaled = importlib.import_module("e5f_ssj_scaled_step_root")
    toeplitz = importlib.import_module("e5f_ssj_toeplitz_jacobian")
    c.jacobian_receipt = read(m["jacobian_source"])
    if ("lag_profiles" not in c.jacobian_receipt or c.jacobian_receipt.get("status") != "complete"
            or c.jacobian_receipt["baseline_max_abs"] > c.jacobian_receipt["baseline_gate"]):
        raise ValueError("Measured derivative receipt is missing or has not passed")
    if not np.isclose(c.old.parameters.psi_child, PSI_SS, rtol=0, atol=1e-14):
        raise ValueError("Initial calibrated preference changed")
    score = read(m["initial_score"])
    if score.get("contract_sha256") != m["initial_target_fingerprint"]:
        raise ValueError("Initial target-and-weight contract changed")
    c.manifest = m
    c.terminal_max_evaluations = int(m.get("terminal_max_evaluations", 24))
    return c, scaled, toeplitz


def endpoint(c, psi, folder, deadline, start):
    import e5f_original_queue_terminal as terminal
    deadline = min(deadline, time.monotonic() + float(c.manifest["terminal_seconds"]))
    controls = dict(c.controls, max_evaluations=int(getattr(c, "terminal_max_evaluations", 24)))
    kwargs = {} if start is None else {"start": np.asarray(start, float)}
    recovered = next((r for r in c.manifest.get("recovery_sources", [])
                      if float(r["psi"]) == float(psi)), None)
    if recovered:
        with gzip.open(recovered["terminal_pickle"], "rb") as stream:
            result = pickle.load(stream)
        if (float(result.parameters.psi_child) != float(psi)
                or result.receipt != read(recovered["terminal_receipt"])):
            raise ValueError("recovered endpoint does not match its preference/receipt")
        dump_pickle(Path(folder) / "terminal.pkl.gz", result)
        save(c.driver, Path(folder) / "root_receipt.json", result.receipt)
    else:
        result = terminal.solve_terminal(old=c.old, psi=float(psi), audit=c.audit,
            controls=controls, deadline=deadline, folder=Path(folder), **kwargs)
    if not getattr(result, "verified", False):
        raise RuntimeError("candidate terminal endpoint was not verified")
    # A candidate endpoint is only usable after a fresh native one-step audit.
    fresh = terminal._evaluate_trial(old=c.old, psi=float(psi),
        coordinates=np.asarray(result.coordinates, float), audit=c.audit,
        deadline=deadline, trial=1)
    audit = terminal._one_step_audit(fresh, c.old, float(psi))
    if not fresh.mapping_valid or audit.get("status") != "passed" or not all(audit.get("checks", {}).values()):
        raise RuntimeError("fresh terminal audit failed")
    for name in ("asset_price",):
        if exact_gap(getattr(fresh, name), getattr(result, name)) > NATIVE_TOL:
            raise RuntimeError("fresh endpoint differs: " + name)
    if exact_gap(fresh.residual, result.receipt["final"]["residual"]) > NATIVE_TOL:
        raise RuntimeError("fresh endpoint residual differs")
    for name in ("V",):
        if exact_gap(getattr(fresh.policy, name), getattr(result.policy, name)) > NATIVE_TOL:
            raise RuntimeError("fresh endpoint policy differs: " + name)
    for name in ("g_pre", "scheduled_entries", "scheduled_raw_entries"):
        if exact_gap(getattr(fresh.state, name), getattr(result.state, name)) > NATIVE_TOL:
            raise RuntimeError("fresh endpoint state differs: " + name)
    save(c.driver, Path(folder) / "verified_endpoint.json", dict(
        verified=True, psi=float(psi), coordinates=np.asarray(result.coordinates, float), audit=audit,
        endpoint_receipt=getattr(result, "receipt", None)))
    return result


def terminal_distance(path, target, verified=True):
    actual = path.person_tail.terminal_state
    mass = max(float(target.g_pre.sum()), 1e-15)
    return dict(distribution_relative_l1=float(np.abs(actual.g_pre-target.g_pre).sum()/mass),
        population_relative_gap=float(abs(actual.g_pre.sum()/mass-1)),
        queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries)/np.asarray(target.scheduled_entries)-1))),
        raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries)/np.asarray(target.scheduled_raw_entries)-1))),
        terminal_steady_state_verified=bool(verified), horizon_comparison_passed=False,
        production_eligible=False)


def tfr(observation):
    value = observation.get("period_tfr_topcode_adjusted")
    if value is None or not math.isfinite(float(value)):
        raise RuntimeError("missing finite period_tfr_topcode_adjusted")
    return float(value)


def initial_prices(old, endpoint, count, warm=None):
    target = np.asarray(endpoint.coordinates, float)
    if warm is not None:
        candidate = np.asarray(warm, float)
        if candidate.shape == (3, count) and np.isfinite(candidate).all():
            return candidate
    # The initial date uses the inherited calibrated prices; only the far tail
    # is anchored at this candidate's verified terminal.
    start = np.asarray([old.policy.price[0], old.parameters.pension,
        old.parameters.property_tax_lump_sum_transfer], float)
    weight = np.linspace(0., 1., count)
    return (1.-weight)[None, :]*start[:, None] + weight[None, :]*target[:, None]


def scaled_root(scaled):
    """Patch the frozen operator with the explicitly pinned scaled-step root."""
    if not hasattr(scaled, "solve_price_path_scaled"):
        raise AttributeError("scaled module misses solve_price_path_scaled")
    return lambda: scaled.solve_price_path_scaled


@contextmanager
def scaled_root_context(c, scaled):
    solver = scaled_root(scaled)
    with patch.object(c.rebated, "_path_root_solver", solver):
        yield


def controls_for(c, toeplitz, receipt, horizon, round_number):
    controls = dict(c.controls)
    for key in ("automatic_fiscal_polish", "fiscal_tolerance", "fiscal_slope"):
        controls.pop(key, None)
    controls["slope"] = controls.pop("market_slope", 1.63)
    controls["max_evaluations"] = 8
    measured, _ = toeplitz.assemble_from_receipt(c.jacobian_receipt, horizon)
    controls["initial_jacobian"] = measured
    controls["default_jacobian"] = measured
    if receipt:
        jac = receipt.get("final_jacobian")
        if jac is not None:
            controls["initial_jacobian"] = jac
        if receipt.get("final_damping") is not None:
            controls["damping"] = receipt["final_damping"]
    return controls


def run_mapping(c, scaled, toeplitz, *, inherited, endpoint_result, psi, count, year,
                folder, deadline, warm=None, receipt=None, round_number=1, graphs=False,
                max_evaluations=8):
    """One native root round; observer state is cleared on every root mapping."""
    import run_e5f_transition_calibration as fertility
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    observations, snapshot, mapping = [], {}, [0]
    boundary = NS(parameters=endpoint_result.parameters, policy=endpoint_result.policy,
        asset_price=endpoint_result.asset_price)
    def observe(i, evaluation, parameters, grid, shared):
        observations.append(dict(period=int(i), calendar_year=int(year)+4*int(i),
            psi_child=float(parameters.psi_child), **fertility.period_fertility_diagnostics(evaluation, parameters)))
        if i == 0:
            snapshot.update(parameters=parameters, b_grid=grid, evaluation=evaluation,
                shared=shared, supply_rule=c.old.supply_rule)
        if i % 10 == 0:
            save(c.driver, folder / "date_progress.json", dict(period=int(i), dates=count))
    def capture(**kwargs):
        observations.clear(); snapshot.clear()
        path = native(**kwargs); mapping[0] += 1
        here = folder / "mappings" / f"mapping_{mapping[0]:02d}"
        save(c.driver, here / "rows.json", path.rows)
        save(c.driver, here / "fertility.json", observations)
        save(c.driver, here / "terminal_distance.json", terminal_distance(path, endpoint_result.state, endpoint_result.verified))
        save(c.driver, here / "native_queue_gates.json", dict(mass=path.maximum_mass_accounting_error,
            policy=path.maximum_policy_reproduction_error, feasibility=path.maximum_feasibility_projection_mass,
            rows=len(path.rows), values=len(path.values)))
        save(c.driver, folder / "latest_completed_mapping.json", dict(round=round_number,
            mapping=mapping[0], calendar_year=year, path=str(here)))
        plot_mapping(here, path.rows, observations)
        if graphs and mapping[0] == 1 and snapshot:
            from run_e5f_successive_surprises_overnight import standard_graphs
            standard_graphs(snapshot, NS(path=path), here / "graphs")
        return path
    def progress(row):
        row = dict(row, round=round_number, calendar_year=year, psi=float(psi))
        save(c.driver, folder / "latest_completed.json", row)
        if row.get("new_best"):
            save(c.driver, folder / "best_so_far.json", row)
    guess = initial_prices(c.old, endpoint_result, count, warm)
    controls = controls_for(c, toeplitz, receipt, count, round_number)
    controls["max_evaluations"] = int(max_evaluations)
    if count >= 100 and c.manifest.get("reserve_verification_time"):
        controls["max_evaluations"] = budgeted_root_evaluations(
            deadline-time.monotonic(), count, controls["max_evaluations"], c.manifest)
        if controls["max_evaluations"] < 2:
            raise TimeoutError("Insufficient time for initial mapping and verification")
    with c.queue.original_queue_adapter(), scaled_root_context(c, scaled):
        native = c.rebated.evaluate_forecast
        with patch.object(c.rebated, "evaluate_forecast", capture):
            result = c.rebated.solve_rebated_forecast(inherited=inherited, psi=float(psi), old_state=c.old,
                terminal=boundary, demographic_primitives=None, count=count, initial_prices=guess[0],
                initial_pensions=guess[1], initial_transfers=guess[2], audit_controls=c.audit,
                root_controls=controls,
                deadline_monotonic=deadline, callback=progress, observer=observe)
    root = dict(result.root_receipt, calendar_year=year, psi=float(psi), count=count,
        round=round_number, population_closure="original_queue", original_queue=True)
    save(c.driver, folder / "root_receipt.json", root)
    if result.path is not None:
        save(c.driver, folder / "rows.json", result.path.rows)
        save(c.driver, folder / "fertility.json", observations)
        save(c.driver, folder / "terminal_distance.json", terminal_distance(result.path, endpoint_result.state, endpoint_result.verified))
        if root.get("finite_horizon_market_fiscal_converged") and snapshot:
            from run_e5f_successive_surprises_overnight import standard_graphs
            standard_graphs(snapshot, result, folder / "graphs")
            dump_pickle(folder / "first_period_diagnostics.pkl.gz", snapshot)
    return result, observations, root


def plot_mapping(folder, rows, fertility):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    years = [r['calendar_year'] for r in rows]
    for axis, key, label in zip(axes.flat, ['asset_price', 'renter_price', 'adult_population',
            'pension_period_units', 'equal_transfer_period_units'],
            ['House price', 'Rent', 'Household population', 'Pension', 'Equal rebate']):
        axis.plot(years, [r[key] for r in rows]); axis.set_title(label)
    axes.flat[5].plot(years, [r['housing_demand'] for r in rows], label='Demand')
    axes.flat[5].plot(years, [r['housing_supply'] for r in rows], ls='--', label='Supply')
    axes.flat[5].legend(); axes.flat[5].set_title('Housing')
    fig.suptitle('Forecast evaluation — consult convergence receipt')
    fig.savefig(folder / 'path.png', dpi=130); fig.savefig(folder / 'path.pdf'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3.5), constrained_layout=True)
    ax.plot(years, [tfr(r) for r in fertility]); ax.set_ylabel('Period fertility')
    ax.set_xlabel('Decision year'); fig.savefig(folder / 'fertility.png', dpi=130)
    fig.savefig(folder / 'fertility.pdf'); plt.close(fig)


def next_guess(root, count):
    best = root.get("best") or root.get("final") or {}
    prices = best.get("prices")
    if prices is None:
        return None
    value = np.asarray(prices, float)
    return value.reshape(3, count) if value.size == 3*count and np.isfinite(value).all() else None


def accepted(result, observations, target):
    if result.path is None or result.next_state is None:
        return False, None
    if not result.root_receipt.get("finite_horizon_market_fiscal_converged"):
        return False, None
    first = observations[0] if observations else None
    if first is None:
        return False, None
    gap = tfr(first) - target
    return abs(gap) <= TOL, gap


def solve_candidate(c, scaled, toeplitz, *, inherited, psi, year, count, target,
                    folder, deadline, endpoint_start, warm_receipt=None, initial_warm=None):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    terminal = endpoint(c, psi, folder / "terminal", deadline, endpoint_start)
    warm, receipt = (next_guess(warm_receipt, count) if warm_receipt else initial_warm), warm_receipt
    final, completed_rounds = None, 0
    for round_number in range(1, int(c.manifest["max_rounds"])+1):
        maximum = int(c.manifest.get("max_root_evaluations", 8))
        if c.manifest.get("reserve_verification_time"):
            maximum = budgeted_root_evaluations(deadline-time.monotonic(), count,
                maximum, c.manifest)
        if time.monotonic() >= deadline or maximum < 2:
            if final is None:
                raise TimeoutError("Insufficient time for initial mapping and verification")
            break  # Persist the completed candidate instead of losing its root receipt.
        final, obs, root = run_mapping(c, scaled, toeplitz, inherited=inherited,
            endpoint_result=terminal, psi=psi, count=count, year=year,
            folder=folder / f"round_{round_number:02d}", deadline=deadline,
            warm=warm, receipt=receipt, round_number=round_number, graphs=(round_number == 1),
            max_evaluations=maximum)
        completed_rounds += 1
        warm, receipt = next_guess(root, count), root
        if final.path is not None and root.get("finite_horizon_market_fiscal_converged"):
            break
    passed, gap = accepted(final, obs, target)
    record = dict(psi=float(psi), calendar_year=year, target=float(target), model_tfr=(None if not obs else tfr(obs[0])),
        gap=gap, accepted=passed, root_receipt=receipt, endpoint_coordinates=np.asarray(terminal.coordinates, float).tolist(),
        terminal_verified=True, rounds=completed_rounds)
    save(c.driver, folder / "candidate.json", record)
    return final, obs, record, terminal


def budgeted_root_evaluations(remaining, count, maximum, manifest):
    """Count the final replay as a full mapping, with extra time for artifacts."""
    seconds = float(manifest["mapping_seconds_budget_104"]) * count / 104.
    reserve = float(manifest["artifact_reserve_seconds"])
    if not math.isfinite(seconds) or seconds <= 0 or reserve < 0:
        raise ValueError("invalid mapping-time budget")
    return max(0, min(int(maximum), int((remaining-reserve)//seconds)))


def recovery_root(source, year, count):
    """Recover only a matching native price/fiscal point; acceptance is rerun."""
    root = read(source["root_receipt"])
    if (int(root["calendar_year"]) != year or int(root["count"]) != count
            or float(root["psi"]) != float(source["psi"])):
        raise ValueError("recovery root belongs to another candidate or calendar")
    rows = read(source["native_rows"])
    if [r["calendar_year"] for r in rows] != list(range(year, year+4*count, 4)):
        raise ValueError("recovery rows have a different calendar")
    coords = [r[k] for k in ("asset_price", "pension_period_units", "equal_transfer_period_units") for r in rows]
    if (not root["best"]["mapping_valid"] or exact_gap(coords, root["best"]["prices"]) > 0
            or any(float(r["psi_child"]) != float(source["psi"]) for r in rows)):
        raise ValueError("recovery native arrays do not match best coordinates")
    return root


def candidate_values(seed, history):
    candidate = next_psi(history, seed, (PSI_SS-.20, PSI_SS+.02),
        [r["psi"] for r in history if "psi" in r])
    if candidate is not None:
        yield candidate


def checkpoint_state(c, path, *, state, row, terminal, year, psi):
    # The root has already used the adapter's first-period replay.  Persist its
    # native next_state before any later stage can read it.
    if state is None:
        raise RuntimeError("refuse to checkpoint absent next_state")
    dump_pickle(path, state)
    with gzip.open(path, "rb") as stream:
        replay = pickle.load(stream)
    replay_gaps = dict(g_pre=exact_gap(replay.households.g_pre, state.households.g_pre),
        scheduled_entries=exact_gap(replay.households.scheduled_entries, state.households.scheduled_entries),
        scheduled_raw_entries=exact_gap(replay.households.scheduled_raw_entries, state.households.scheduled_raw_entries))
    if replay.year != state.year or max(replay_gaps.values()) > 0.:
        raise RuntimeError("accepted next-state pickle replay failed")
    proof = dict(calendar_year=year+4, psi=float(psi), next_state_pickle=str(path),
        next_state_sha256=sha(path), replay_preserved=True, replay_gaps=replay_gaps, accepted_row=row,
        endpoint_coordinates=np.asarray(terminal.coordinates, float).tolist())
    save(c.driver, Path(path).with_suffix(".json"), proof)
    return proof


def stage(c, scaled, toeplitz, m, stage_index, deadline):
    out = Path(m["output"]); year, count, target = YEARS[stage_index], COUNTS[stage_index], TARGETS[stage_index]
    if count != horizon(year, int(m.get("terminal_year", 2423))):
        raise ValueError("stage horizon does not end at the manifest terminal boundary")
    smoke_receipt = out / "smoke" / "passed.json"
    if not smoke_receipt.exists() or not read(smoke_receipt).get("passed"):
        raise ValueError("long stages require the passed exact-loop smoke")
    folder = out / f"stage_{stage_index}_{year}"
    if folder.exists():
        raise ValueError("refuse duplicate stage output")
    if stage_index == 0:
        inherited = c.rebated.InheritedState(2007, c.old.initial_state)
        stage_warm = None
    else:
        checkpoint = read(out / f"stage_{stage_index-1}_{YEARS[stage_index-1]}" / "accepted.json")["checkpoint"]
        previous = out / f"stage_{stage_index-1}_{YEARS[stage_index-1]}" / "accepted_next_state.pkl.gz"
        if checkpoint["next_state_sha256"] != sha(previous) or int(checkpoint["calendar_year"]) != year:
            raise ValueError("prior accepted checkpoint receipt differs")
        with gzip.open(previous, "rb") as stream:
            inherited = pickle.load(stream)
        if int(inherited.year) != year:
            raise ValueError("prior checkpoint has the wrong surprise date")
        rows = read(out / f"stage_{stage_index-1}_{YEARS[stage_index-1]}" / "accepted_rows.json")
        stage_warm = np.asarray([[r["asset_price"] for r in rows[1:]],
            [r["pension_period_units"] for r in rows[1:]],
            [r["equal_transfer_period_units"] for r in rows[1:]]], float)
        if stage_warm.shape != (3, count): raise ValueError("prior accepted tail has wrong horizon")
    save(c.driver, folder / "contract.json", dict(stage=stage_index, calendar_year=year, count=count,
        target=target, target_measurement="period_tfr_topcode_adjusted", tolerance=TOL,
        psi_bounds=[PSI_SS-.20, PSI_SS+.02], original_queue=True,
        horizon_endpoint_year=2423, production_eligible=False))
    history, best, endpoint_start = [], None, None
    recoveries = m.get("recovery_sources", []) if stage_index == 0 else []
    for trial in range(1, int(m.get("max_trials", 6))+1):
        if time.monotonic() >= deadline: break
        recovered = recoveries[trial-1] if trial <= len(recoveries) else None
        proposed = ([float(recovered["psi"])] if recovered else
                    list(candidate_values(SEEDS[stage_index], history)))
        if not proposed: break
        psi = proposed[0]
        try:
            candidate_deadline = min(deadline, time.monotonic()+float(m.get("candidate_seconds", 36000)))
            warm_receipt = (recovery_root(recovered, year, count) if recovered else
                            (best or {}).get("root_receipt"))
            result, obs, record, terminal = solve_candidate(c, scaled, toeplitz, inherited=inherited,
                psi=psi, year=year, count=count, target=target, folder=folder/f"trial_{trial:02d}",
                deadline=candidate_deadline, endpoint_start=endpoint_start, warm_receipt=warm_receipt,
                initial_warm=stage_warm)
            history.append(record); endpoint_start=terminal.coordinates
            if record.get("gap") is not None and (best is None or abs(record["gap"]) < abs(best.get("gap", float("inf")))):
                best = record
                save(c.driver, folder / "best_so_far.json", best)
            save(c.driver, folder / "latest_completed.json", record)
            if m.get("preserve_unconverged_candidate") and record.get("gap") is None:
                save(c.driver, folder / "failure.json", dict(error="candidate requires continuation",
                    candidate=record, candidates=history, calendar_year=year))
                return
            if record["accepted"]:
                proof = checkpoint_state(c, folder / "accepted_next_state.pkl.gz", state=result.next_state,
                    row=record, terminal=terminal, year=year, psi=psi)
                save(c.driver, folder / "accepted_rows.json", result.path.rows)
                save(c.driver, folder / "accepted_fertility.json", obs)
                save(c.driver, folder / "accepted.json", dict(record, checkpoint=proof, candidates=history))
                dispatch_next(m, stage_index, folder)
                return
        except Exception as exc:
            reject = dict(trial=trial, psi=float(psi), rejected=True, error_type=type(exc).__name__, error=str(exc))
            if m.get("preserve_unconverged_candidate"):
                reject["recoverable_root_receipts"] = [str(p) for p in
                    sorted((folder/f"trial_{trial:02d}").glob("round_*/root_receipt.json"))]
            history.append(reject); save(c.driver, folder / "latest_completed.json", reject)
            if m.get("preserve_unconverged_candidate") and isinstance(exc, TimeoutError):
                save(c.driver, folder / "failure.json", dict(error="candidate requires continuation",
                    candidate=reject, candidates=history, calendar_year=year))
                return
    save(c.driver, folder / "failure.json", dict(error="no accepted candidate", candidates=history,
        valid_best=best, calendar_year=year))


def dispatch_next(m, stage_index, folder):
    scripts = m.get("scripts", m.get("next_scripts", {}))
    key = "policy" if stage_index == 3 else f"stage_{stage_index+1}"
    script = scripts.get(key, scripts.get(key.replace("_", "")))
    receipt = Path(folder) / "dispatch.json"
    if script is None:
        save_path = dict(dispatched=False, reason="no next script in manifest", next_key=key)
    elif receipt.exists():
        return
    else:
        run = subprocess.run(["sbatch", "--parsable", str(script)], check=True, capture_output=True, text=True)
        save_path = dict(dispatched=True, next_key=key, script=str(script), job_id=run.stdout.strip().split(";")[0])
    # c.driver is unavailable here; manifest runs on clusters with atomic outputs
    receipt.parent.mkdir(parents=True, exist_ok=True)
    temporary = receipt.with_suffix(".tmp"); temporary.write_text(json.dumps(save_path, indent=2)+"\n"); temporary.replace(receipt)


@contextmanager
def adopted_tax_validator(c, old, annual_tax):
    """Retain the native supply check while allowing only tax/user-cost changes."""
    import e5f_original_queue_terminal as terminal
    native = terminal._validate_supply
    def validate(candidate):
        if candidate.supply_rule is not c.old.supply_rule:
            raise ValueError("policy changed the fixed supply rule")
        P = candidate.parameters
        if not np.isclose(P.tau_H, P.period_years*annual_tax, rtol=0, atol=1e-15):
            raise ValueError("annual property tax conversion differs")
        if not np.isclose(P.user_cost_rate, P.R_gross+P.delta+P.tau_H-1., rtol=0, atol=1e-15):
            raise ValueError("policy user-cost identity fails")
        for name in ("H0", "r_bar", "xi_supply"):
            np.testing.assert_array_equal(getattr(P, name), getattr(c.old.parameters, name))
        baseline = copy.copy(candidate); baseline.parameters = copy.deepcopy(P)
        baseline.parameters.tau_H = old.parameters.tau_H
        baseline.parameters.user_cost_rate = old.parameters.user_cost_rate
        native(baseline)
    with patch.object(terminal, "_validate_supply", validate):
        yield


def smoke(c, scaled, toeplitz, m, deadline):
    out = Path(m["output"])/"smoke"; out.mkdir(parents=True, exist_ok=False)
    for i, source in enumerate(m.get("recovery_sources", []), 1):
        recovered_root = recovery_root(source, 2007, 104)
        recovered_terminal = endpoint(c, source["psi"], out/f"recovered_terminal_{i}", deadline, None)
        save(c.driver, out/f"recovery_input_{i}.json", dict(passed=True,
            psi=source["psi"], root_source=source["root_receipt"],
            native_coordinates_match=True, terminal_freshly_verified=recovered_terminal.verified,
            old_best_score=recovered_root["best"]["score"], old_root_accepted=False))
    # Exact stationary six-date root smoke.
    psi = float(c.old.parameters.psi_child)
    terminal = endpoint(c, psi, out/"stationary_terminal", deadline,
        [c.old.policy.price[0], c.old.parameters.pension, c.old.parameters.property_tax_lump_sum_transfer])
    stationary = c.rebated.InheritedState(2007, terminal.state)
    result, obs, root = run_mapping(c, scaled, toeplitz, inherited=stationary, endpoint_result=terminal,
        psi=psi, count=6, year=2007, folder=out/"stationary_six_date", deadline=deadline,
        warm=np.repeat(np.asarray(terminal.coordinates)[:, None], 6, axis=1), graphs=False)
    if result.next_state is None or result.path is None or not obs or tfr(obs[0]) <= 0:
        raise RuntimeError("six-date stationary smoke failed first-state/fertility proof")
    proof = checkpoint_state(c, out/"smoke_next_state.pkl.gz", state=result.next_state,
        row=root, terminal=terminal, year=2007, psi=psi)
    drift = {name: exact_gap(getattr(result.next_state.households, name), getattr(terminal.state, name))
             for name in ('g_pre', 'scheduled_entries', 'scheduled_raw_entries')}
    save(c.driver, out/'stationary_drift.json', dict(maximum_absolute_gaps=drift,
        native_stationary_audit=terminal.receipt['one_step_audit']))
    # Two valid root-loop mappings with a fixed-price HOUSEHOLD boundary test
    # the changed preference without imposing the whole terminal capital loss.
    # This boundary is not a stationary general equilibrium.
    shifted = float(m['psi_seeds'][0])
    from e5f_social_security import bind_social_security_income
    pf = c.joined.pf
    Q = copy.deepcopy(c.old.parameters); Q.psi_child = shifted
    bind_social_security_income(Q, pension_period=Q.pension, payroll_tax=.179)
    q = float(c.old.policy.price[0])
    shared = pf.calendar.model.precompute_shared(Q, c.old.b_grid)
    objects = pf.calendar.model.solve_bellman_full_markov_income(
        np.array([Q.user_cost_rate*q]), np.array([q]), Q, c.old.b_grid, shared)
    hh_policy = pf.policy_from_objects(objects, q, Q, c.old.b_grid, shared)
    coordinates = np.array([q, Q.pension, Q.property_tax_lump_sum_transfer])
    shifted_terminal = NS(parameters=Q, policy=hh_policy, asset_price=q,
        coordinates=coordinates, state=c.old.initial_state, verified=False)
    actual = c.rebated.InheritedState(2007, c.old.initial_state)
    _, shifted_obs, shifted_root = run_mapping(c, scaled, toeplitz, inherited=actual, endpoint_result=shifted_terminal,
        psi=shifted, count=2, year=2007, folder=out/"shifted_two_date", deadline=deadline,
        warm=np.repeat(coordinates[:, None], 2, axis=1), graphs=False, max_evaluations=2)
    if not shifted_obs or tfr(shifted_obs[0]) <= 0:
        raise RuntimeError("shifted-psi fertility measurement smoke failed")
    if not shifted_root.get('final', {}).get('mapping_valid'):
        raise RuntimeError('Shifted native root-loop replay is invalid')
    save(c.driver, out/"passed.json", dict(passed=True, stationary_dates=6, shifted_dates=2,
        shifted_mode="two native root mappings with a fixed-price household boundary; no GE claim",
        first_state_checkpoint=proof, changed_preference=shifted,
        measured_derivative_receipt_sha256=sha(m['jacobian_source'])))
    dispatch_next(m, -1, out)


def policy(c, scaled, toeplitz, m, deadline):
    """The policy driver inherits accepted 2023 state and uses a fresh 2% terminal."""
    out = Path(m["output"])/"policy"; out.mkdir(parents=True, exist_ok=False)
    prior = Path(m["output"])/"stage_3_2019"/"accepted_next_state.pkl.gz"
    accepted_row = read(Path(m["output"])/"stage_3_2019"/"accepted.json")
    if (not accepted_row.get('accepted')
            or accepted_row['checkpoint']['next_state_sha256'] != sha(prior)):
        raise ValueError('Policy requires the verified fitted 2023 checkpoint')
    with gzip.open(prior, "rb") as stream: inherited = pickle.load(stream)
    if int(inherited.year) != 2023:
        raise ValueError('Policy inherited year is not 2023')
    inherited_rows = read(Path(m["output"])/"stage_3_2019"/"accepted_rows.json")
    psi = float(accepted_row["psi"])
    save(c.driver, out/'baseline_rows.json', inherited_rows[1:])
    baseline_fertility = read(Path(m['output'])/'stage_3_2019'/'accepted_fertility.json')
    save(c.driver, out/'baseline_fertility.json', baseline_fertility[1:])
    if len(inherited_rows[1:]) != 100 or inherited_rows[1]['calendar_year'] != 2023:
        raise ValueError('Policy and fitted baseline must share 100 dates from 2023')
    # Match the established adopted-tax validator: change only annual tax and
    # its user-cost identity while retaining the calibrated fixed supply rule.
    old = copy.copy(c.old); old.parameters = copy.deepcopy(c.old.parameters)
    old.parameters.tau_H = float(old.parameters.period_years)*.02
    old.parameters.user_cost_rate = float(old.parameters.R_gross)+float(old.parameters.delta)+float(old.parameters.tau_H)-1.
    original = c.old; c.old = old
    try:
        with adopted_tax_validator(c, original, .02):
            terminal = endpoint(c, psi, out/"terminal_tax2", deadline,
                [old.policy.price[0], old.parameters.pension, old.parameters.property_tax_lump_sum_transfer])
            warm = np.asarray([[r["asset_price"] for r in inherited_rows[1:]],
                [r["pension_period_units"] for r in inherited_rows[1:]],
                [r["equal_transfer_period_units"] for r in inherited_rows[1:]]], float)
            result = obs = root = None; receipt = None
            for round_number in range(1, int(m.get("max_rounds", 4))+1):
                result, obs, root = run_mapping(c, scaled, toeplitz, inherited=inherited, endpoint_result=terminal,
                    psi=psi, count=100, year=2023, folder=out/"tax2"/f"round_{round_number:02d}",
                    deadline=deadline, warm=warm, receipt=receipt, round_number=round_number, graphs=(round_number == 1))
                if result.path is not None and root.get("finite_horizon_market_fiscal_converged"): break
                warm, receipt = next_guess(root, 100), root
    finally:
        c.old = original
    save(c.driver, out/"complete.json", dict(completed=result.path is not None,
        psi=psi, annual_property_tax=.02, fixed_supply_curve=True,
        finite_horizon_market_fiscal_converged=root.get("finite_horizon_market_fiscal_converged"),
        terminal_distance_saved=True, horizon_comparison_passed=False))


def heartbeat(c, folder, deadline, stop):
    while not stop.wait(60):
        save(c.driver, Path(folder)/"heartbeat.json", dict(remaining_seconds=max(0., deadline-time.monotonic())))
        if time.monotonic() >= deadline:
            save(c.driver, Path(folder)/"failure.json", dict(error="manifest hard deadline")); os._exit(124)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mode", choices=("smoke", "stage", "policy"), required=True)
    parser.add_argument("--stage", type=int, choices=range(4))
    args = parser.parse_args()
    if (args.mode == "stage") != (args.stage is not None):
        parser.error("--stage is required only in stage mode")
    m = read(args.manifest); manifest_sha = manifest_contract(m, args.manifest)
    hard_deadline = float(m.get("absolute_deadline_unix", m["total_deadline_unix"]))
    budget_key = {"smoke": "smoke_seconds", "stage": "stage_seconds", "policy": "policy_seconds"}[args.mode]
    mode_budget = float(m[budget_key])
    hard_deadline = min(hard_deadline, time.time()+mode_budget)
    remaining = hard_deadline-time.time()
    if remaining < 120: raise TimeoutError("manifest deadline already expired")
    deadline = time.monotonic()+remaining
    c, scaled, toeplitz = load_frozen_context(m)
    out = Path(m["output"]); out.mkdir(parents=True, exist_ok=True)
    run_folder = out/("smoke" if args.mode == "smoke" else ("policy" if args.mode == "policy" else f"stage_{args.stage}_{YEARS[args.stage]}"))
    stop = threading.Event(); threading.Thread(target=heartbeat, args=(c, run_folder, deadline, stop), daemon=True).start()
    try:
        save(c.driver, out/"driver_startup.json", dict(mode=args.mode, stage=args.stage,
            manifest_sha256=manifest_sha, deadline_unix=hard_deadline,
            source="original_queue", terminal_endpoint_per_candidate=True))
        if args.mode == "smoke": smoke(c, scaled, toeplitz, m, deadline)
        elif args.mode == "stage": stage(c, scaled, toeplitz, m, args.stage, deadline)
        else: policy(c, scaled, toeplitz, m, deadline)
    except BaseException as exc:
        save(c.driver, run_folder/"failure.json", dict(error_type=type(exc).__name__, error=str(exc)))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
