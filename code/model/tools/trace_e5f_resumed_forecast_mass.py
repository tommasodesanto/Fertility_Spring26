#!/usr/bin/env python3
"""Trace one resumed native forecast without changing its model or root controls."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback
from unittest.mock import patch

import numpy as np


class _SingleForecastFinished(BaseException):
    """Private non-Exception stop signal, bypassing driver fallback handling."""


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, Path(path).resolve())
    if spec is None or spec.loader is None:
        raise ImportError("Cannot load pinned module: " + str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def verify(path, expected, label):
    actual = sha256(path)
    if actual != expected.lower():
        raise ValueError(f"Pinned {label} changed: {actual} != {expected.lower()}")


def clean(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def exception_record(exc):
    return dict(type=type(exc).__name__, message=str(exc), repr=repr(exc),
                traceback="".join(traceback.format_exception(type(exc), exc,
                                                              exc.__traceback__)))


def validate_resume_contract(manifest, expected_year):
    resume = manifest.get("resume_history")
    if not isinstance(resume, dict):
        raise ValueError("Trace manifest must resume an accepted history prefix")
    windows = resume.get("windows")
    if (not isinstance(windows, list) or len(windows) != 1
            or windows[0].get("year") != 2007 or expected_year != 2011):
        raise ValueError("Trace requires exactly the accepted 2007 prefix and 2011 restart")
    if manifest.get("reuse_forecast_jacobian", False) is not False:
        raise ValueError("Reset-vintage trace requires reuse_forecast_jacobian=false")


def validate_cache_proof(proof, cache_sha256):
    if (proof.get("status") != "verified"
            or proof.get("cache_sha256") != cache_sha256
            or proof.get("exact_mapping_equal") is not True
            or proof.get("household_and_accounting_mapping_valid") is not True):
        raise ValueError("Exact policy cache lacks its matching native proof")


def manifest_pins(manifest):
    pins = manifest.get("file_sha256")
    if not isinstance(pins, dict) or not pins:
        raise ValueError("Scientific manifest lacks explicit file pins")
    return [(Path(path), digest, "manifest source " + str(path))
            for path, digest in pins.items()]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--driver-sha256", required=True)
    parser.add_argument("--diagnostic", type=Path, required=True)
    parser.add_argument("--diagnostic-sha256", required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--cache-sha256", required=True)
    parser.add_argument("--cache-proof", type=Path, required=True)
    parser.add_argument("--cache-proof-sha256", required=True)
    parser.add_argument("--case", choices=("A0", "A+"), required=True)
    parser.add_argument("--count", type=int, choices=(6, 24, 100), required=True)
    parser.add_argument("--expected-year", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seconds", type=float, required=True)
    args = parser.parse_args(argv)
    if args.expected_year != 2011:
        parser.error("--expected-year must be 2011 for this resumed-prefix trace")
    if not np.isfinite(args.seconds) or not 0 < args.seconds <= 1800:
        parser.error("--seconds must be in (0,1800]")
    pins = ((args.manifest, args.manifest_sha256, "manifest"),
            (args.driver, args.driver_sha256, "history driver"),
            (args.diagnostic, args.diagnostic_sha256, "mass diagnostic"),
            (args.cache, args.cache_sha256, "exact policy cache"),
            (args.cache_proof, args.cache_proof_sha256, "native cache proof"))
    for path, digest, label in pins:
        verify(path, digest, label)
    manifest = read_json(args.manifest)
    validate_resume_contract(manifest, args.expected_year)
    scientific_pins = manifest_pins(manifest)
    for path, digest, label in scientific_pins:
        verify(path, digest, label)
    plan = read_json(manifest["prior_plan"])
    plan_pins = manifest_pins(plan)
    for path, digest, label in plan_pins:
        verify(path, digest, "prior plan " + label)
    proof = read_json(args.cache_proof)
    validate_cache_proof(proof, args.cache_sha256.lower())
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError("Trace output already exists: " + str(output))
    output.mkdir(parents=True)
    trace_output = output / "mass_trace"
    trace_output.mkdir()
    started = time.monotonic()
    source_pins = {str(Path(path).resolve()): digest.lower()
                   for path, digest, _ in pins}
    write_json(output / "trace_contract.json", dict(
        case=args.case, count=args.count, expected_year=args.expected_year,
        maximum_forecasts=1, maximum_seconds=args.seconds,
        initial_jacobian="none_reset_vintage",
        scientific_contract_unchanged=True, policy_reserve_seconds=0,
        cache_proof_status=proof.get("status"),
        source_pins=source_pins, production_eligible=False))

    driver = load_module("run_e5f_final_rebated_history", args.driver)
    diagnostic = load_module("diagnose_e5f_transition_mass", args.diagnostic)
    cache = load_module("e5f_exact_policy_cache", args.cache)
    # This single-forecast diagnostic has no subsequent policy stage to reserve.
    # The resume validator permits these two explicit budget-only differences.
    trace_manifest = dict(manifest, policy_reserve_seconds=0,
                          forecast_seconds=args.seconds)
    trace_manifest_path = output / "trace_manifest.json"
    write_json(trace_manifest_path, trace_manifest)
    resumed_fit = read_json(manifest["resume_history"]["realized_fit"]["path"])
    if len(resumed_fit) != 1 or resumed_fit[0].get("year") != 2007:
        raise ValueError("Resume fit ledger is not exactly the accepted 2007 prefix")
    seed_step = driver.initial_seed_step(plan, manifest)
    if seed_step != -0.01:
        raise ValueError("This trace requires the deterministic center-minus-.01 trial")
    expected_psi = float(resumed_fit[0]["psi"]) + seed_step
    captured = dict(calls=0)
    original = driver.solve_forecast

    def one_forecast(**kwargs):
        captured["calls"] += 1
        if captured["calls"] != 1:
            captured["wrapper_error"] = dict(type="RuntimeError",
                message="Driver attempted more than one forecast")
            raise _SingleForecastFinished()
        captured["folder"] = str(Path(kwargs["folder"]).resolve())
        inherited_year = getattr(kwargs.get("inherited"), "year", None)
        if (inherited_year != args.expected_year or kwargs.get("case") != args.case
                or kwargs.get("count") != args.count):
            captured["wrapper_error"] = dict(type="ValueError",
                message="First intercepted forecast differs from expected year/case/count",
                inherited_year=inherited_year, case=kwargs.get("case"),
                count=kwargs.get("count"))
            raise _SingleForecastFinished()
        if kwargs.get("initial_jacobian") is not None:
            captured["wrapper_error"] = dict(type="ValueError",
                message="Reset-vintage forecast unexpectedly supplied a Jacobian")
            raise _SingleForecastFinished()
        captured["psi"] = float(kwargs["psi"])
        if captured["psi"] != expected_psi:
            captured["wrapper_error"] = dict(type="ValueError",
                message="First resumed trial is not accepted-2007 psi minus .01",
                expected_psi=expected_psi, actual_psi=captured["psi"])
            raise _SingleForecastFinished()
        profile, profile_state = diagnostic._profile_factory(trace_output, 1e-8)
        previous_profile = sys.getprofile()
        cache_gib = 24 if args.count == 100 else 6
        sys.setprofile(profile)
        try:
            import e5f_rebated_surprises as rebated
            _, joined, *_ = rebated._runtime()
            with cache.policy_cache(joined.pf,
                                    max_bytes=cache_gib * 1024**3) as cache_stats:
                try:
                    result, detail = original(**dict(kwargs, initial_jacobian=None))
                    captured["result"] = result
                    captured["detail"] = detail
                    captured["outcome"] = "returned"
                except Exception as exc:
                    captured["outcome"] = "raised"
                    captured["raised"] = exception_record(exc)
                finally:
                    captured["cache"] = cache_stats.snapshot()
        finally:
            sys.setprofile(previous_profile)
            captured["profile_state"] = profile_state
        raise _SingleForecastFinished()

    try:
        with patch.object(driver, "solve_forecast", one_forecast):
            driver.main(["--manifest", str(trace_manifest_path),
                         "--case", args.case, "--count", str(args.count),
                         "--output", str(output), "--seconds", str(args.seconds)])
    except _SingleForecastFinished:
        pass
    except Exception as exc:
        captured["driver_error"] = exception_record(exc)
    elapsed = time.monotonic() - started
    failure_arguments = trace_output / "first_failure_arguments.pkl.gz"
    replay = None
    if failure_arguments.is_file():
        try:
            diagnostic._replay(failure_arguments, trace_output)
            replay = dict(status="completed",
                output=str(trace_output / "single_cohort_replay.json"))
        except Exception as exc:
            replay = dict(status="raised", exception=exception_record(exc))
    else:
        replay = dict(status="not_run", reason="No profiled cohort crossed 1e-8")
    raised = captured.get("raised")
    first_failure = (captured.get("profile_state") or {}).get("first_failure")
    reproduced = diagnostic._is_expected_gate_error(raised)
    if captured.get("driver_error"):
        status = "driver_failed_before_single_forecast_completed"
    elif captured.get("wrapper_error"):
        status = "wrapper_contract_failed"
    elif captured.get("outcome") == "returned":
        status = "forecast_returned_no_native_exception"
    elif reproduced and first_failure is not None:
        status = "native_mass_failure_reproduced_and_profiled"
    elif captured.get("outcome") == "raised":
        status = "forecast_raised_without_profiled_mass_reproduction"
    else:
        status = "no_forecast_intercepted"
    result = captured.get("result")
    root = getattr(result, "root_receipt", None) if result is not None else None
    summary = dict(status=status, case=args.case, count=args.count,
        expected_year=args.expected_year, forecast_calls=captured["calls"],
        forecast_outcome=captured.get("outcome"), psi=captured.get("psi"),
        elapsed_seconds=elapsed, maximum_seconds=args.seconds,
        raised=raised, wrapper_error=captured.get("wrapper_error"),
        driver_error=captured.get("driver_error"),
        transition_calls=(captured.get("profile_state") or {}).get("calls", 0),
        first_failure=first_failure, expected_mass_gate_exception=reproduced,
        replay=replay, cache=captured.get("cache"),
        finite_horizon_market_fiscal_converged=(
            None if root is None else root.get("finite_horizon_market_fiscal_converged")),
        final_reproduction_max_abs=(
            None if root is None else root.get("final_reproduction_max_abs")),
        historical_prefix="accepted_2007_resume_only",
        horizon_verified=False, production_eligible=False)
    write_json(output / "trace_summary.json", summary)
    for path, digest, label in pins:
        verify(path, digest, label)
    for path, digest, label in scientific_pins:
        verify(path, digest, label)
    for path, digest, label in plan_pins:
        verify(path, digest, "prior plan " + label)
    print(json.dumps(clean(summary)), flush=True)
    if (captured.get("calls") != 1 or captured.get("wrapper_error")
            or captured.get("driver_error")):
        return 2
    if replay["status"] == "raised":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
