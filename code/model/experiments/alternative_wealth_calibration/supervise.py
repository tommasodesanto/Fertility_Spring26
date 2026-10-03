"""Persistent, memory-aware local launch queue for ten independent chains.

The supervisor never retries a chain or changes the calibration contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
MIN_SECONDS_FOR_START = 2400
MIN_DISK_BYTES = 20 * 1024**3
POLL_SECONDS = 30


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--starts-file", required=True, type=Path)
    ap.add_argument("--starts-sha256", required=True)
    ap.add_argument("--hard-deadline-epoch", required=True, type=float)
    ap.add_argument("--max-live", type=int, default=6)
    ap.add_argument("--per-chain-wall-seconds", type=float, default=9900)
    ap.add_argument("--min-available-gib", type=float, default=7.0)
    ap.add_argument("--max-started", type=int, default=10)
    args = ap.parse_args()
    if args.out.exists():
        raise RuntimeError("Refusing existing supervisor output")
    if args.max_started != 10 or not 1 <= args.max_live <= 5:
        raise RuntimeError("Ten-chain count or memory cap drift")
    if not 2400 <= args.per_chain_wall_seconds <= 21600:
        raise RuntimeError("Invalid bounded per-chain wall budget")
    starts = args.starts_file.resolve()
    if sha(starts) != args.starts_sha256:
        raise RuntimeError("Start table hash drift")
    plan = json.loads(starts.read_text())
    if len(plan["starts"]) != 10:
        raise RuntimeError("Start table is not ten distinct starts")
    if len({tuple(sorted(x.items())) for x in plan["starts"]}) != 10:
        raise RuntimeError("Duplicate starts")
    driver, entry = HERE / "calibrate.py", HERE / "local_entry.py"
    source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in (driver, entry, Path(__file__))}
    args.out.mkdir(parents=True)
    out = args.out.resolve()
    env = os.environ.copy()
    for key in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
                "NUMEXPR_MAX_THREADS"):
        env[key] = "1"
    records = {str(i): dict(chain=i, status="queued") for i in range(10)}
    live = {}
    launch = dict(status="running", started_epoch=time.time(),
                  hard_deadline_epoch=args.hard_deadline_epoch,
                  max_live=args.max_live, min_available_gib=args.min_available_gib,
                  maximum_objective_calls_per_chain=250,
                  maximum_wall_seconds_per_chain=args.per_chain_wall_seconds,
                  final_native_reserve_seconds=1800,
                  max_started=args.max_started,
                  starts_file=str(starts), starts_sha256=args.starts_sha256,
                  original_target_fingerprint=plan["original_target_fingerprint"],
                  experimental_target_fingerprint=plan["target_fingerprint"],
                  experimental_weight_fingerprint=plan["weight_fingerprint"],
                  numeric_wealth_weight_unchanged=True,
                  alternative_interest_timing_only=True,
                  source_sha256=source_hashes,
                  pid=os.getpid(), no_automatic_restart=True)
    save(out / "launch_receipt.json", launch)
    stop = False

    def request_stop(*_):
        nonlocal stop
        stop = True

    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    try:
        while True:
            now = time.time()
            for idx, (proc, stream) in list(live.items()):
                rc = proc.poll()
                if rc is None:
                    heartbeat = out / f"chain_{idx:02d}" / "heartbeat.json"
                    if heartbeat.exists():
                        age = now - heartbeat.stat().st_mtime
                        records[str(idx)]["heartbeat_age_seconds"] = round(age, 1)
                        records[str(idx)]["health"] = "stale_30_minutes" if age >= 1800 else "heartbeat_recent"
                    continue
                stream.close()
                receipt = out / f"chain_{idx:02d}" / "completed.json"
                records[str(idx)].update(status="completed" if rc == 0 and receipt.exists() else "failed",
                                         returncode=rc, ended_epoch=now,
                                         completion_receipt_present=receipt.exists())
                live.pop(idx)
            remaining = [i for i in range(10) if records[str(i)]["status"] == "queued"]
            if now >= args.hard_deadline_epoch:
                for i in remaining:
                    records[str(i)].update(status="not_started_deadline", ended_epoch=now)
                for idx, (proc, _) in live.items():
                    try:
                        prior = records[str(idx)].get("hard_deadline_signal_epoch")
                        if prior is None:
                            os.killpg(proc.pid, signal.SIGTERM)
                            records[str(idx)]["hard_deadline_signal_epoch"] = now
                        elif now - prior >= 60:
                            os.killpg(proc.pid, signal.SIGKILL)
                            records[str(idx)]["hard_deadline_kill_epoch"] = now
                    except ProcessLookupError:
                        pass
            elif stop:
                for i in remaining:
                    records[str(i)].update(status="not_started_supervisor_stop", ended_epoch=now)
                for idx, (proc, _) in live.items():
                    try:
                        os.killpg(proc.pid, signal.SIGTERM)
                        records[str(idx)]["stop_signal_epoch"] = now
                    except ProcessLookupError:
                        pass
            else:
                current_hashes = {str(p.relative_to(ROOT)): sha(p) for p in (driver, entry, Path(__file__))}
                source_ok = current_hashes == source_hashes and sha(starts) == args.starts_sha256
                disk_ok = os.statvfs(out).f_bavail * os.statvfs(out).f_frsize >= MIN_DISK_BYTES
                available = psutil.virtual_memory().available
                free_ok = available >= args.min_available_gib * 1024**3
                launch["source_pin_current"] = source_ok
                launch["disk_budget_current"] = disk_ok
                launch["available_memory_gib"] = round(available / 1024**3, 3)
                if source_ok and disk_ok and free_ok and len(live) < args.max_live and remaining:
                    i = remaining[0]
                    deadline = min(now + args.per_chain_wall_seconds,
                                   args.hard_deadline_epoch)
                    if deadline - now < MIN_SECONDS_FOR_START:
                        records[str(i)].update(status="not_started_insufficient_time", ended_epoch=now)
                    else:
                        chain_out = out / f"chain_{i:02d}"
                        log = (out / f"chain_{i:02d}.log").open("w")
                        cmd = [sys.executable, str(entry), str(driver),
                               "--arm", "alternative", "--chain", str(i),
                               "--out", str(chain_out), "--deadline-epoch", str(deadline),
                               "--starts-file", str(starts),
                               "--starts-file-sha256", args.starts_sha256]
                        proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log,
                                                stderr=subprocess.STDOUT,
                                                start_new_session=True)
                        live[i] = (proc, log)
                        records[str(i)].update(status="running", pid=proc.pid,
                                               started_epoch=now, deadline_epoch=deadline,
                                               log=str(out / f"chain_{i:02d}.log"))
            launch["updated_epoch"] = time.time()
            launch["running_count"] = len(live)
            launch["queued_count"] = sum(r["status"] == "queued" for r in records.values())
            launch["completed_count"] = sum(r["status"] == "completed" for r in records.values())
            launch["failed_count"] = sum(r["status"] == "failed" for r in records.values())
            save(out / "supervisor_heartbeat.json", launch)
            save(out / "chains.json", records)
            if not live and not any(r["status"] == "queued" for r in records.values()):
                break
            time.sleep(POLL_SECONDS)
    finally:
        launch["status"] = "terminal"
        launch["updated_epoch"] = time.time()
        save(out / "supervisor_heartbeat.json", launch)
        save(out / "chains.json", records)


if __name__ == "__main__":
    main()
