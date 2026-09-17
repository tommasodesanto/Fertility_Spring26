"""Experiment-only: re-solve the economy at a larger tenure-smoothing scale and run the 104-date root.

Stages, all in copies of the frozen stack (nothing in production changes):

1. Stationary equilibrium at the frozen scale (control) and at the probe
   scale, both at the original preference, via the frozen
   ``e5f_original_queue_terminal.solve_terminal``.  The two
   ``endpoint_reference`` moment dictionaries form the fit table.
2. Terminal stationary equilibrium at the probe scale and the final
   announced preference, warm-started from the frozen endpoint coordinates.
3. The announced 104-date root through the announced batch's own
   ``run_path`` with: inherited 2007 state = the probe-scale stationary
   state (removes the date-1 re-sorting artifact), terminal = the probe-scale
   endpoint, start guess = a saved 104-date checkpoint, optional measured
   Toeplitz ``initial_jacobian``, optional direction-preserving step and
   trimmed acceptance, and the retained gates unless a tolerance vector is
   requested.  The announced preference path is kept as fitted at 0.005.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from e5f_ssj_scaled_step_root import solve_price_path_scaled
from run_e5f_ssj_announced_rescue import block_tolerance_vector, select_jacobian, warm_coordinates

HORIZON = 104
FIT_KEYS = ("asset_price", "renter_price", "pension_period", "transfer", "owner_rate", "housing_demand",
            "housing_supply", "adjusted_births", "raw_births", "entry_flow", "unit_mass", "population_scale",
            "renewal_ratio")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=_default) + "\n")
    tmp.replace(path)


def _default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return str(value)


def fit_table(control_reference, probe_reference, control_scale, probe_scale, keys=FIT_KEYS):
    """Side-by-side stationary moments at the two smoothing scales (pure)."""
    rows = []
    for key in keys:
        a, b = control_reference.get(key), probe_reference.get(key)
        if a is None or b is None or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            continue
        rel = (b - a) / abs(a) if a != 0 else None
        rows.append(dict(moment=key, control=float(a), probe=float(b), relative_change=rel))
    lines = [f"| moment | scale {control_scale} | scale {probe_scale} | relative change |", "|---|---|---|---|"]
    for r in rows:
        rel = "" if r["relative_change"] is None else f"{100 * r['relative_change']:+.2f}%"
        lines.append(f"| {r['moment']} | {r['control']:.6g} | {r['probe']:.6g} | {rel} |")
    return rows, "\n".join(lines) + "\n"


def with_scale(old, kappa):
    probe = copy.copy(old)
    probe.parameters = copy.deepcopy(old.parameters)
    probe.parameters.tenure_choice_kappa = float(kappa)
    return probe


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    m = read(args.manifest)
    for path, digest in m["file_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Changed pinned input: " + path)
    announced = read(m["announced_manifest"])
    kappa = float(m["tenure_choice_kappa"])
    out = Path(m["output"])
    if (out / "experiment_contract.json").exists():
        raise ValueError("Refusing to overwrite a started experiment")
    sys.path.insert(0, str(Path(m["announced_source_dir"])))
    import run_e5f_announced_original_queue as ann
    spec = read(announced["spec"])
    sys.path.insert(0, str(Path(spec["batch"]) / "source"))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(announced["spec"])
    c.spec_path = Path(announced["spec"])
    import e5f_original_queue_terminal as terminal_solver
    frozen_endpoint, frozen_receipt = ann.endpoint_from_manifest(announced)
    frozen_kappa = float(getattr(c.old.parameters, "tenure_choice_kappa", 0.0))
    psi0 = float(c.old.parameters.psi_child)
    started = time.time()
    end = started + float(m["seconds"])
    deadline = time.monotonic() + (end - time.time())
    stationary_cap = float(m.get("stationary_seconds", 3600))
    save(out / "experiment_contract.json", dict(
        manifest_sha256=sha(args.manifest), tenure_choice_kappa=kappa, frozen_tenure_choice_kappa=frozen_kappa,
        horizon=HORIZON, seconds=m["seconds"], stationary_seconds=stationary_cap, started_unix=started, deadline_unix=end,
        changes_relative_to_retained_root=dict(
            tenure_choice_kappa=f"{frozen_kappa} -> {kappa} on a deep copy of the initial-state parameters",
            initial_2007_state="re-solved stationary state at the probe scale (frozen packet state not used)",
            terminal="re-solved at the probe scale and the final announced preference",
            preference_path="announced levels kept as fitted at the frozen scale",
            start_guess=m.get("warm_start_checkpoint"),
            initial_jacobian=m.get("jacobian_mode", "diagonal"), step_rule=m.get("step_rule", "clipped"),
            mapping_budget=int(m.get("mapping_budget", 8)), trim_count=int(m.get("trim_count", 0)),
            gates=("retained max-abs 2e-4" if m.get("fiscal_gate_scaled") is None
                   else f"housing {m.get('housing_gate', 2e-4)}, fiscal {m['fiscal_gate_scaled']} scaled"),
            per_mapping_plots_skipped=bool(m.get("skip_mapping_plots", True)),
            stationary_start=m.get("stationary_start"), stationary_evaluations=int(m.get("stationary_evaluations", 16)),
            frozen_operator_validation="root_controls carry max_evaluations=8 and market_tolerance=2e-4; the solver copy applies the budget/gate above"),
        production_eligible=False))
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            save(out / "controller_heartbeat.json", dict(remaining_seconds=end - time.time()))
            if time.time() >= end:
                save(out / "controller_failure.json", dict(error="Experiment deadline"))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        controls = dict(c.controls, max_evaluations=16)
        probe_controls = dict(c.controls, max_evaluations=int(m.get("stationary_evaluations", 16)))
        if probe_controls["max_evaluations"] not in (16, 24):
            raise ValueError("stationary_evaluations must be 16 or 24 (frozen terminal solver)")
        probe_start = np.asarray(m["stationary_start"], dtype=float) if m.get("stationary_start") else None
        probe_old = with_scale(c.old, kappa)
        with c.queue.original_queue_adapter(), runner.original_receipts(c), \
                c.cache.policy_cache(c.joined.pf, max_bytes=12 * 1024**3):
            # Stage 1: stationary equilibria at the original preference.
            control = terminal_solver.solve_terminal(old=c.old, psi=psi0, audit=c.audit, controls=controls,
                                                     deadline=min(deadline, time.monotonic() + stationary_cap),
                                                     folder=out / "stationary_control")
            if not control.verified:
                raise RuntimeError("Control stationary solve at the frozen scale did not verify")
            probe = terminal_solver.solve_terminal(old=probe_old, psi=psi0, audit=c.audit, controls=probe_controls,
                                                   deadline=min(deadline, time.monotonic() + stationary_cap),
                                                   folder=out / "stationary_probe", start=probe_start)
            if not probe.verified:
                raise RuntimeError("Probe-scale stationary solve did not verify")
            rows, markdown = fit_table(control.endpoint_reference, probe.endpoint_reference, frozen_kappa, kappa)
            save(out / "fit_table.json", dict(rows=rows, control_coordinates=control.coordinates, probe_coordinates=probe.coordinates,
                 note="Stationary moments only; the full SMM target table needs the calibration collector."))
            (out / "fit_table.md").write_text(markdown)
            save(out / "stage_status.json", dict(stationary="complete", elapsed_seconds=time.time() - started))
            # Stage 2: terminal at the probe scale and final preference.
            endpoint = terminal_solver.solve_terminal(old=probe_old, psi=float(ann.FINAL_PSI), audit=c.audit, controls=probe_controls,
                                                      deadline=min(deadline, time.monotonic() + stationary_cap),
                                                      folder=out / "endpoint_probe",
                                                      start=np.asarray(frozen_endpoint.coordinates, dtype=float))
            if not endpoint.verified:
                raise RuntimeError("Probe-scale terminal solve did not verify")
            save(out / "stage_status.json", dict(stationary="complete", endpoint="complete", elapsed_seconds=time.time() - started))
            # Stage 3: the 104-date announced root at the probe scale.
            probe_old.initial_state = probe.state
            probe_old.policy = probe.policy
            c.old = probe_old
            checkpoint = read(m["warm_start_checkpoint"])
            guess, previous_best = warm_coordinates(dict(best=None), HORIZON, checkpoint)
            toeplitz_receipt = read(m["toeplitz_receipt"]) if m.get("toeplitz_receipt") else None
            jacobian_mode = m.get("jacobian_mode", "diagonal")
            jacobian, jacobian_info = ((None, dict(mode="diagonal")) if jacobian_mode == "diagonal"
                                       else select_jacobian(jacobian_mode, {}, toeplitz_receipt, HORIZON))
            budget = int(m.get("mapping_budget", 8))
            tolerance = (block_tolerance_vector(HORIZON, float(m.get("housing_gate", 2e-4)), float(m["fiscal_gate_scaled"]))
                         if m.get("fiscal_gate_scaled") is not None else None)
            psi_path = np.r_[announced["psi_levels"], np.full(100, ann.FINAL_PSI)]
            evaluate, first = ann.announced_queue_path(c, psi_path)
            native_solve = c.rebated.solve_rebated_forecast
            retained_solver = c.rebated._path_root_solver()

            def solve_with_jacobian(**kwargs):
                controls_ = dict(kwargs["root_controls"])
                if controls_.get("max_evaluations") != 8:
                    raise ValueError("Frozen operator validation expects the retained eight-mapping control")
                if jacobian is not None:
                    controls_["initial_jacobian"] = jacobian
                return native_solve(**dict(kwargs, root_controls=controls_))

            def step_solver(**kwargs):
                if m.get("step_rule", "clipped") != "scaled":
                    return retained_solver(**kwargs)
                kwargs = dict(kwargs, max_evaluations=budget, trim_count=int(m.get("trim_count", 0)))
                if tolerance is not None:
                    kwargs["tolerance_vector"] = tolerance
                return solve_price_path_scaled(**kwargs)
            import run_e5f_successive_surprises_overnight as overnight
            native_plot, native_graphs = ann.plot_packet, overnight.standard_graphs
            skip = bool(m.get("skip_mapping_plots", True))

            def plot_packet(c_, m_, folder_, *a, **k):
                return None if (skip and "mappings" in Path(folder_).parts) else native_plot(c_, m_, folder_, *a, **k)

            def standard_graphs(snapshot, result, folder_):
                return None if (skip and "mappings" in Path(folder_).parts) else native_graphs(snapshot, result, folder_)
            folder = out / "run"
            with patch.object(c.rebated, "evaluate_forecast", evaluate), patch.object(c.rebated, "first_period_state", first), \
                    patch.object(c.rebated, "_path_root_solver", lambda: step_solver), \
                    patch.object(c.rebated, "solve_rebated_forecast", solve_with_jacobian), \
                    patch.object(ann, "plot_packet", plot_packet), patch.object(overnight, "standard_graphs", standard_graphs):
                result = ann.run_path(c, endpoint, announced, folder, deadline, guess=guess)
        history = [dict(evaluation=e["evaluation"], phase=e["phase"], score=e["score"], raw_max_abs=e.get("raw_max_abs", e["score"]),
                        trimmed_score=e.get("trimmed_score", e["score"]), evaluation_seconds=e.get("evaluation_seconds"),
                        safeguard=e.get("safeguard")) for e in result.root_receipt.get("history", []) if "evaluation" in e]
        save(out / "comparison.json", dict(history=history, previous_best_score_at_frozen_scale=previous_best,
             best_score=(result.root_receipt.get("best") or {}).get("score"),
             converged=bool(result.root_receipt.get("finite_horizon_market_fiscal_converged")),
             jacobian=jacobian_info, tenure_choice_kappa=kappa, production_eligible=False))
        save(out / "controller_complete.json", dict(completed=True, tenure_choice_kappa=kappa,
             finite_horizon_market_fiscal_converged=bool(result.root_receipt.get("finite_horizon_market_fiscal_converged")),
             stationary_probe_verified=True, endpoint_probe_verified=True, production_eligible=False,
             elapsed_seconds=time.time() - started))
    except BaseException as exc:
        save(out / "controller_failure.json", dict(error_type=type(exc).__name__, error=str(exc), elapsed_seconds=time.time() - started))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
