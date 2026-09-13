"""Read-only mass trace for one saved e5f successive-surprise forecast.

The driver deliberately does not alter the transition source or its gate.  A
small context pickle is made by the launch wrapper from the saved trial's
``inherited``, ``old_state``, ``demographics``, ``terminal``, ``prices``,
``pensions`` and ``psi`` objects.  This keeps reconstruction explicit and
prevents silently substituting a nearby forecast receipt.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import sys
from pathlib import Path

import numpy as np


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return {"dtype": str(value.dtype), "shape": list(value.shape),
                "sum": float(np.sum(value)), "min": float(np.min(value)),
                "max": float(np.max(value))}
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _mass(x) -> float:
    return float(np.sum(np.asarray(x, dtype=float)))


def _is_expected_gate_error(raised) -> bool:
    return bool(raised and raised.get("type") == "RuntimeError"
                and "mass gate failed" in raised.get("message", ""))


def _profile_factory(output: Path, tolerance: float):
    """Capture transition locals on return; no function arithmetic is patched."""
    state = {"calls": 0, "first_failure": None}
    target_code = None

    def profile(frame, event, arg):
        nonlocal target_code
        if event == "call" and frame.f_code.co_name == "advance_cohort_one_period_markov_income":
            target_code = frame.f_code
        if event != "return" or target_code is None or frame.f_code is not target_code:
            return profile
        state["calls"] += 1
        loc = frame.f_locals
        names = ("gj", "j", "gpl", "gpt", "gps", "g_next", "tenure_probs", "P")
        stage_names = ("gj", "gpl", "gpt", "gps", "g_next")
        masses = {name: _mass(loc[name]) for name in stage_names if name in loc}
        expected = masses.get("gj")
        actual = masses.get("g_next")
        gap = (abs(actual - expected) / expected
               if expected is not None and expected > 0 and actual is not None else None)
        rows = None
        if loc.get("tenure_probs") is not None and "gpl" in loc:
            probs = np.asarray(loc["tenure_probs"])
            # Restrict the stored full policy tensor to this age and income
            # shock.  Its remaining axes match gpl exactly before tenure.
            j = int(loc["j"])
            probs = probs[:, :, :, j, :, :, :, :]
            rowerr = np.sum(probs.astype(float), axis=-1) - 1.0
            rows = {
                "dtype": str(probs.dtype),
                "max_abs_row_sum_minus_one": float(np.max(np.abs(rowerr))),
                "mean_row_sum_minus_one": float(np.mean(rowerr)),
                "mass_weighted_row_sum_minus_one": None,
            }
            gp = np.asarray(loc["gpl"], dtype=float)
            rows["mass_weighted_row_sum_minus_one"] = float(np.sum(gp * rowerr) / max(np.sum(gp), 1e-300))
            rows["gpl_dtype"] = str(gp.dtype)
        record = {"call": state["calls"], "age_index": loc.get("j"), "masses": masses,
                  "relative_gap_gnext_vs_gj": gap, "tenure_rows": rows}
        with (output / "transition_calls.jsonl").open("a") as stream:
            stream.write(json.dumps(_jsonable(record)) + "\n")
        if gap is not None and gap > tolerance and state["first_failure"] is None:
            formal = ("gj", "j", "loc_probs", "tenure_choice", "tenure_probs",
                      "bp_pol", "P", "b_grid", "SD", "lmm_idx", "lmm_wt",
                      "tmx_idx", "tmx_wt", "ust", "Pia", "Pi_z")
            captured = {name: loc[name] for name in formal if name in loc}
            with gzip.open(output / "first_failure_arguments.pkl.gz", "wb", compresslevel=1) as stream:
                pickle.dump(captured, stream, protocol=5)
            np.save(output / 'first_failure_output.npy', loc['g_next'])
            state["first_failure"] = record
        return profile

    return profile, state


def _call_forecast(context, output: Path, tolerance: float):
    import e5f_successive_surprises as surprise

    profile, state = _profile_factory(output, tolerance)
    old_profile = sys.getprofile()
    sys.setprofile(profile)
    try:
        result = surprise.evaluate_forecast(**context)
        raised = None
    except Exception as exc:  # The original gate is expected to raise here.
        result = None
        raised = {"type": type(exc).__name__, "message": str(exc)}
    finally:
        sys.setprofile(old_profile)
    return result, raised, state


def _replay(path: Path, output: Path):
    with gzip.open(path, "rb") as stream:
        captured = pickle.load(stream)
    P = captured["P"]
    original = bool(getattr(P, "use_numba_scatter", False))
    setattr(P, "use_numba_scatter", False)
    import run_e5f_open_population_transition as transition
    try:
        # The wrapper owns the active model through calendar.model; configure
        # it explicitly so this replay cannot accidentally use a stale module.
        transition.configure_sequential_model()
        got = transition.calendar.model.advance_cohort_one_period_markov_income(**captured)
        replay = {"use_numba_scatter": False, "mass": _mass(got),
                  "maximum_abs_change_from_compiled": float(np.max(np.abs(got - np.load(output / 'first_failure_output.npy')))),
                  "input_mass": _mass(captured["gj"]),
                  "relative_gap": abs(_mass(got) - _mass(captured["gj"])) / _mass(captured["gj"])}
        renorm = dict(captured)
        probs = np.asarray(captured["tenure_probs"], dtype=np.float64).copy()
        denom = probs.sum(axis=-1, keepdims=True)
        np.divide(probs, denom, out=probs, where=denom > 0.0)
        renorm["tenure_probs"] = probs
        got64 = transition.calendar.model.advance_cohort_one_period_markov_income(**renorm)
        replay["float64_row_renormalized_tenure_probs"] = {
            "mass": _mass(got64),
            "relative_gap": abs(_mass(got64) - _mass(captured["gj"])) / _mass(captured["gj"]),
            "diagnostic_only": True,
        }
    finally:
        setattr(P, "use_numba_scatter", original)
    (output / "single_cohort_replay.json").write_text(json.dumps(replay, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", type=Path, required=True,
                    help="Explicit pickle of the exact failed trial evaluate_forecast kwargs")
    ap.add_argument("--source-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--tolerance", type=float, default=1.0e-8)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path[:0] = [str(args.source_root / "code/model/tools"), str(args.source_root / "code/model")]
    # Initialize the same sequential runtime before unpickling model objects.
    import run_e5f_open_population_transition as transition
    transition.configure_sequential_model()
    context = pickle.loads(args.context.read_bytes())
    if set(context) != {"inherited", "old_state", "demographics", "prices", "pensions", "psi", "terminal"}:
        raise ValueError("Context keys must be the exact evaluate_forecast kwargs")
    result, raised, state = _call_forecast(context, args.output, args.tolerance)
    receipt = {"context_sha256": _sha256(args.context), "source_root": str(args.source_root),
               "source_sha256": _sha256(args.source_root / "code/model/tools/run_e5f_open_population_transition.py"),
               "raised": raised, "calls": state["calls"], "first_failure": state["first_failure"]}
    (args.output / "trace_receipt.json").write_text(json.dumps(_jsonable(receipt), indent=2) + "\n")
    failure = args.output / "first_failure_arguments.pkl.gz"
    if failure.exists():
        _replay(failure, args.output)
    replay_ok = (args.output / "single_cohort_replay.json").exists()
    reproduced = _is_expected_gate_error(raised)
    return 0 if reproduced and state["first_failure"] is not None and replay_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
