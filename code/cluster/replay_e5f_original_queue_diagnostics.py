"""Replay one saved original-queue price path with the frozen native kernels.

The worker performs exactly one direct ``queue_path`` evaluation.  It does not
solve a root, update prices, or import the person/headship population law.  The
``--prepare`` front end snapshots the two immutable coordinate files and the
verified stationary reference, then writes an explicit ``scp``/``sbatch``
launch recipe for Torch.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import shlex
import shutil
import subprocess
import sys
import threading
import time
from types import SimpleNamespace as NS

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"

np = None


COUNT = 100
TOL = 2e-10
NUMERICAL_PYTHON = "/share/apps/anaconda3/2025.06/bin/python"
DEFAULT_REMOTE_SPEC = (
    "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/"
    "batches/afternoon_original_queue_20260913a/spec.json"
)
DEFAULT_REMOTE_ENDPOINT = (
    "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/"
    "batches/afternoon_original_queue_20260913a_terminal_restart_v1/endpoint/terminal.pkl.gz"
)
DEFAULT_REMOTE_RECEIPT = (
    "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/"
    "batches/afternoon_original_queue_20260913a_terminal_restart_v1/endpoint/root_receipt.json"
)
DEFAULT_REMOTE_REFERENCE = (
    "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/"
    "batches/afternoon_original_queue_20260913a_terminal_10/transition/stationary_reference.json"
)


def read(path: Path | str):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def jsonable(value):
    if np is not None and isinstance(value, np.ndarray):
        return value.tolist()
    if np is not None and isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def validate_snapshot(snapshot: Path):
    rows = read(snapshot / "rows.json")
    latest = read(snapshot / "latest_completed.json")
    if not isinstance(rows, list) or len(rows) != COUNT:
        raise ValueError("Frozen rows.json must contain exactly 100 rows")
    if not isinstance(latest, dict) or latest.get("evaluation") != 3:
        raise ValueError("Frozen latest_completed.json must be saved iteration3")
    coordinates = latest.get("prices", [])
    if (not isinstance(coordinates, list) or len(coordinates) != 3 * COUNT
            or any(not math.isfinite(float(v)) for v in coordinates)):
        raise ValueError("Saved iteration3 coordinates must contain 300 finite values")
    prices, pensions, transfers = (coordinates[:COUNT], coordinates[COUNT:2 * COUNT],
                                   coordinates[2 * COUNT:])
    if (any(float(v) <= 0 for v in prices)
            or any(float(v) < 0 for v in pensions)
            or any(float(v) < 0 for v in transfers)):
        raise ValueError("Saved path has invalid price or fiscal coordinates")
    expected_years = [2007 + 4 * i for i in range(COUNT)]
    for i, row in enumerate(rows):
        if not isinstance(row, dict) or row.get("calendar_year") != expected_years[i]:
            raise ValueError("Frozen rows are not the 2007--2403 four-year path")
        for key, value in (("asset_price", prices[i]),
                           ("pension_period_units", pensions[i]),
                           ("equal_transfer_period_units", transfers[i])):
            if key not in row or abs(float(row[key]) - float(value)) > TOL:
                raise ValueError(f"Saved coordinate disagrees with frozen row {i}: {key}")
    return rows, latest, prices, pensions, transfers


def compare_rows(actual, expected) -> dict:
    if len(actual) != len(expected):
        raise ValueError(f"Native replay returned {len(actual)} rows; expected {len(expected)}")
    max_gap = 0.0
    failures = []
    for i, (got, want) in enumerate(zip(actual, expected, strict=True)):
        if set(got) != set(want):
            raise ValueError(f"Row {i} keys differ: native={set(got) ^ set(want)}")
        for key in want:
            a, b = got[key], want[key]
            if isinstance(b, (bool, str)) or b is None:
                if a != b:
                    failures.append((i, key, a, b))
                continue
            try:
                gap = abs(float(a) - float(b))
            except (TypeError, ValueError):
                if a != b:
                    failures.append((i, key, a, b))
                continue
            if not np.isfinite(gap) or gap > TOL:
                failures.append((i, key, a, b))
            max_gap = max(max_gap, gap)
    if failures:
        i, key, a, b = failures[0]
        raise ValueError(f"Native row replay differs at row {i}, {key}: {a!r} vs {b!r}")
    return dict(rows=COUNT, numeric_max_abs_gap=max_gap, tolerance=TOL, passed=True)


def endpoint_and_reference(endpoint_path: Path, receipt_path: Path, reference_path: Path,
                           spec: dict):
    receipt = read(receipt_path)
    with gzip.open(endpoint_path, "rb") as stream:
        endpoint = pickle.load(stream)
    if not getattr(endpoint, "verified", False) or not receipt.get("verified", False):
        raise ValueError("Saved endpoint is not verified")
    if getattr(endpoint, "receipt", None) != receipt:
        raise ValueError("Endpoint pickle receipt differs from root_receipt.json")
    if not np.isclose(float(endpoint.parameters.psi_child), float(spec["permanent_psi"]),
                      rtol=0.0, atol=1e-14):
        raise ValueError("Endpoint preference does not match spec permanent_psi")
    reference = read(reference_path)
    terminal = reference.get("terminal", {})
    quantities = terminal.get("quantities", {}) if isinstance(terminal, dict) else {}
    price = float(endpoint.asset_price)
    pension = float(endpoint.parameters.pension)
    rent = float(endpoint.parameters.user_cost_rate) * price
    checks = {"asset_price": (price, quantities.get("asset_price")),
              "renter_price": (rent, quantities.get("renter_price")),
              "pension_period": (pension, quantities.get("pension_period"))}
    for name, (actual, saved) in checks.items():
        if saved is None or not np.isclose(actual, float(saved), rtol=0.0, atol=TOL):
            raise ValueError(f"Verified stationary reference does not match endpoint: {name}")
    if not reference.get("stationary_endpoint_verified", False):
        raise ValueError("Stationary reference is not marked endpoint-verified")
    return endpoint, reference


def distance(actual, target) -> dict:
    actual_pre = np.asarray(actual.g_pre, dtype=float)
    target_pre = np.asarray(target.g_pre, dtype=float)
    target_mass = max(float(target_pre.sum()), 1e-15)
    actual_queue = np.asarray(actual.scheduled_entries, dtype=float)
    target_queue = np.asarray(target.scheduled_entries, dtype=float)
    actual_raw = np.asarray(actual.scheduled_raw_entries, dtype=float)
    target_raw = np.asarray(target.scheduled_raw_entries, dtype=float)
    return dict(
        actual_terminal_population=float(actual_pre.sum()),
        target_terminal_population=float(target_pre.sum()),
        distribution_relative_l1=float(np.abs(actual_pre - target_pre).sum() / target_mass),
        population_relative_gap=float(abs(actual_pre.sum() / target_mass - 1.0)),
        actual_scheduled_entries=actual_queue,
        target_scheduled_entries=target_queue,
        queue_relative_max=float(np.max(np.abs(actual_queue / target_queue - 1.0))),
        actual_scheduled_raw_entries=actual_raw,
        target_scheduled_raw_entries=target_raw,
        raw_queue_relative_max=float(np.max(np.abs(actual_raw / target_raw - 1.0))),
        finite_path_converged=False,
        production_eligible=False,
    )


def worker(args: argparse.Namespace) -> None:
    global np
    import numpy as np

    contract = read(args.input_contract)
    for path, digest in contract['file_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Changed pinned replay input: ' + path)
    spec_path = args.spec.resolve()
    snapshot = args.snapshot.resolve()
    endpoint_path = args.endpoint.resolve()
    receipt_path = args.endpoint_receipt.resolve()
    reference_path = args.stationary_reference.resolve()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing to overwrite a nonempty replay output")
    output.mkdir(parents=True, exist_ok=True)
    frozen_rows, latest, prices, pensions, transfers = validate_snapshot(snapshot)
    spec = read(spec_path)
    sys.path.insert(0, str(Path(spec["batch"]) / "source"))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(spec_path)
    c.spec_path = spec_path
    import run_e5f_transition_calibration as fertility
    smoke = read(c.spec["smoke_summary"])
    if smoke.get("status") != "passed" or smoke.get("spec_sha256") != c.driver.sha(spec_path):
        raise ValueError("Frozen exact-loop smoke prerequisite is missing or mismatched")
    endpoint, reference = endpoint_and_reference(endpoint_path, receipt_path, reference_path, spec)
    remaining = min(float(args.seconds), 3600.0, float(spec['absolute_deadline_unix']) - time.time())
    if remaining < 120:
        raise TimeoutError('Original experiment deadline is too close for a replay')
    deadline = time.monotonic() + remaining
    input_hashes = {str(p): sha(p) for p in
                    (spec_path, snapshot / "rows.json", snapshot / "latest_completed.json",
                     endpoint_path, receipt_path, reference_path)}
    save(output / "experiment_contract.json", dict(
        replay="saved_iteration3_fixed_price_path", count=COUNT,
        years=COUNT * 4, evaluation=latest["evaluation"], coordinates_layout="prices[0:100],pensions[100:200],rebates[200:300]",
        input_sha256=input_hashes, endpoint_verified=True, root_solves=0,
        mappings=1, population_law=spec["population_law"], no_immigration=True,
        production_eligible=False, tolerance=TOL, hard_cap_seconds=3600))
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            save(output / "controller_heartbeat.json", dict(
                status="running", remaining_seconds=max(0.0, deadline - time.monotonic()),
                completed_periods=len(observations), count=COUNT))
            if time.monotonic() >= deadline:
                save(output / "controller_failure.json", dict(error="Replay hard deadline"))
                os._exit(124)
    observations = []
    threading.Thread(target=heartbeat, daemon=True).start()

    def observe(i, evaluation, P, grid, shared):
        from e5f_balanced_terminal import _household_checks
        if (i != len(observations) or float(P.pension) != float(pensions[i])
                or float(P.tau_pay) != 0.179
                or float(P.property_tax_lump_sum_transfer) != float(transfers[i])):
            raise ValueError('Dated fiscal path differs from saved coordinates')
        rents = c.joined.pf.rents_from_asset_prices(prices, endpoint.asset_price, c.old.parameters)
        _, gates = _household_checks(evaluation, P, shared, grid, float(rents[i]), c.primitive, c.audit)
        if not all(gates.values()):
            raise ValueError('Native dated household audit failed')
        diagnostic = fertility.period_fertility_diagnostics(evaluation, P)
        record = dict(period=i, calendar_year=2007 + 4 * i,
                      **jsonable(diagnostic),
                      birth_children_native=float(evaluation.births),
                      birth_children_topcode_adjusted_native=float(
                          np.sum(diagnostic["birth_flow_topcode_adjusted"])),
                      age_mass_native=np.asarray(diagnostic["age_mass"], dtype=float).tolist(),
                      native_g_pre_mass=float(np.sum(evaluation.g_pre)),
                      native_g_current_mass=float(np.sum(evaluation.g_current)))
        observations.append(record)
        c.driver.save(output / "progress" / f"period_{i:03d}.json", record)
        c.driver.save(output / "latest_completed_period.json", record)

    inherited = c.rebated.InheritedState(2007, c.old.initial_state)
    terminal = NS(parameters=endpoint.parameters, policy=endpoint.policy,
                  asset_price=endpoint.asset_price)
    try:
        with c.queue.original_queue_adapter(), c.cache.policy_cache(c.joined.pf, max_bytes=12 * 1024**3):
            result = c.queue.queue_path(inherited=inherited, old_state=c.old,
                prices=prices, pensions=pensions, transfers=transfers,
                psi=float(spec["permanent_psi"]), terminal=terminal,
                observer=observe, demographics=None)
        if len(observations) != COUNT or len(result.values) != COUNT + 1:
            raise ValueError("Native replay did not evaluate exactly 100 dates")
        row_check = compare_rows(result.rows, frozen_rows)
        c.driver.save(output / "rows.json", result.rows)
        c.driver.save(output / "fertility.json", observations)
        shutil.copy2(reference_path, output / "stationary_reference.json")
        c.driver.save(output / "row_reproduction.json", row_check)
        c.driver.save(output / "terminal_distance.json",
                      distance(result.person_tail.terminal_state, endpoint.state))
        c.driver.save(output / "root_receipt.json", dict(
            converged=False, finite_horizon_market_fiscal_converged=False,
            terminal_distance_passed=False, status="Fixed-price-path replay only; no root solve",
            mapping_source="saved iteration3", production_eligible=False))
        c.driver.save(output / 'irf_contract.json', dict(
            label='Permanent preference decline: 100 periods (400 years), iteration 3',
            status_label='UNCONVERGED ROOT; exact replay of saved path.',
            shock_description='Preference falls 38.1% once and remains constant.'))
        subprocess.run([sys.executable, str(args.plotter), '--case-dir', str(output),
                        '--fertility-rate', '--include-pre-shock'], check=True)
        for path, digest in contract['file_sha256'].items():
            if sha(path) != digest:
                raise ValueError('Replay input changed during execution: ' + path)
        save(output / "controller_complete.json", dict(
            status="complete", replay_rows=COUNT, numeric_max_abs_gap=row_check["numeric_max_abs_gap"],
            actual_terminal_population=float(result.person_tail.terminal_state.g_pre.sum()),
            production_eligible=False))
    except BaseException as exc:
        save(output / "controller_failure.json", dict(error_type=type(exc).__name__, error=str(exc)))
        raise
    finally:
        stop.set()


def prepare(args: argparse.Namespace) -> None:
    stage = args.prepare.resolve()
    snapshot = args.snapshot.resolve()
    spec = args.spec.resolve()
    reference = args.stationary_reference.resolve()
    if stage.exists() and any(stage.iterdir()):
        raise ValueError("Refusing to overwrite a nonempty preparation directory")
    validate_snapshot(snapshot)
    if not spec.is_file() or not reference.is_file():
        raise FileNotFoundError("Local spec and stationary reference are required")
    staged_snapshot = stage / "snapshot"
    staged_snapshot.mkdir(parents=True, exist_ok=True)
    for name in ("rows.json", "latest_completed.json"):
        shutil.copy2(snapshot / name, staged_snapshot / name)
    shutil.copy2(reference, stage / "stationary_reference.json")
    worker_source = Path(__file__).resolve()
    remote_input = args.remote_input.rstrip("/")
    remote_output = args.remote_output.rstrip("/")
    remote_script = remote_input + "/replay_e5f_original_queue_diagnostics.py"
    remote_snapshot = remote_input + "/snapshot"
    remote_reference = remote_input + "/stationary_reference.json"
    sbatch = stage / "replay.sbatch"
    worker_argv = [NUMERICAL_PYTHON, remote_script, "--worker",
                   "--spec", args.remote_spec, "--snapshot", remote_snapshot,
                   "--endpoint", args.remote_endpoint, "--endpoint-receipt", args.remote_receipt,
                   "--stationary-reference", remote_reference, "--output", remote_output,
                   "--seconds", "3600", '--input-contract', remote_input + '/worker_contract.json',
                   '--plotter', remote_input + '/build_e5f_stationary_shock_figures.py']
    sbatch.write_text("\n".join([
        "#!/bin/bash", "#SBATCH --job-name=e5f_orig_replay3", "#SBATCH --account=torch_pr_570_general",
        "#SBATCH --cpus-per-task=1", "#SBATCH --mem=32G", "#SBATCH --time=60",
        f"#SBATCH --output={remote_input}/replay_%j.log", "set -euo pipefail",
        "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1",
        shlex.join(worker_argv), ""]), encoding="utf-8")
    launch = stage / "launch.sh"
    plotter = args.plotter.resolve()
    remote_digests = subprocess.run(
        ['ssh', args.remote_host, 'sha256sum', args.remote_endpoint, args.remote_receipt],
        check=True, capture_output=True, text=True)
    hashes = {line.split(maxsplit=1)[1].strip(): line.split(maxsplit=1)[0]
              for line in remote_digests.stdout.splitlines() if line.strip()}
    hashes.update({args.remote_spec: sha(spec), remote_script: sha(worker_source),
        remote_input + '/build_e5f_stationary_shock_figures.py': sha(plotter),
        remote_snapshot + '/rows.json': sha(staged_snapshot / 'rows.json'),
        remote_snapshot + '/latest_completed.json': sha(staged_snapshot / 'latest_completed.json'),
        remote_reference: sha(reference)})
    save(stage / 'worker_contract.json', dict(file_sha256=hashes, mappings=1, root_solves=0))
    launch.write_text("\n".join([
        "#!/bin/bash", "set -euo pipefail",
        shlex.join(["ssh", args.remote_host, "mkdir", "-p", remote_input, remote_snapshot, remote_output]),
        shlex.join(["scp", "-p", str(worker_source), f"{args.remote_host}:{remote_script}"]),
        shlex.join(["scp", "-p", str(staged_snapshot / "rows.json"), str(staged_snapshot / "latest_completed.json"),
                    f"{args.remote_host}:{remote_snapshot}/"]),
        shlex.join(["scp", "-p", str(stage / "stationary_reference.json"), f"{args.remote_host}:{remote_reference}" ]),
        shlex.join(['scp', '-p', str(plotter), str(stage / 'worker_contract.json'),
                    f'{args.remote_host}:{remote_input}/']),
        shlex.join(["scp", "-p", str(sbatch), f"{args.remote_host}:{remote_input}/replay.sbatch"]),
        shlex.join(["ssh", args.remote_host, "sbatch", f"{remote_input}/replay.sbatch"]), ""]), encoding="utf-8")
    os.chmod(launch, 0o755)
    save(stage / "prepare_manifest.json", dict(spec=str(spec), spec_sha256=sha(spec),
        snapshot=str(snapshot), snapshot_sha256={name: sha(snapshot / name) for name in ("rows.json", "latest_completed.json")},
        stationary_reference_sha256=sha(reference), remote_host=args.remote_host,
        remote_input=remote_input, remote_output=remote_output, remote_spec=args.remote_spec,
        remote_endpoint=args.remote_endpoint, remote_receipt=args.remote_receipt,
        launch_script_sha256=sha(launch), sbatch_sha256=sha(sbatch), mappings=1, root_solves=0))
    print(json.dumps(dict(stage=str(stage), launch=str(launch), command=f"bash {shlex.quote(str(launch))}")))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", type=Path)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--stationary-reference", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--endpoint", type=Path, default=Path(DEFAULT_REMOTE_ENDPOINT))
    parser.add_argument("--endpoint-receipt", type=Path, default=Path(DEFAULT_REMOTE_RECEIPT))
    parser.add_argument("--seconds", type=float, default=3600)
    parser.add_argument('--input-contract', type=Path)
    parser.add_argument('--plotter', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'model/tools/build_e5f_stationary_shock_figures.py')
    parser.add_argument("--remote-host", default="torch")
    parser.add_argument("--remote-input", default="/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/replay_e5f_original_queue_iteration3")
    parser.add_argument("--remote-output", default="/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/replay_e5f_original_queue_iteration3/output")
    parser.add_argument("--remote-spec", default=DEFAULT_REMOTE_SPEC)
    parser.add_argument("--remote-endpoint", default=DEFAULT_REMOTE_ENDPOINT)
    parser.add_argument("--remote-receipt", default=DEFAULT_REMOTE_RECEIPT)
    args = parser.parse_args()
    if bool(args.prepare) == bool(args.worker):
        parser.error("Choose exactly one of --prepare or --worker")
    if args.stationary_reference is None:
        if args.prepare:
            parser.error("--stationary-reference is required with --prepare")
        args.stationary_reference = Path(DEFAULT_REMOTE_REFERENCE)
    if args.snapshot is None or (args.worker and args.output is None):
        parser.error("--snapshot is required, and --output is required with --worker")
    if args.worker and args.input_contract is None:
        parser.error('--input-contract is required with --worker')
    return args


if __name__ == "__main__":
    arguments = parse_args()
    if arguments.prepare:
        prepare(arguments)
    else:
        worker(arguments)
