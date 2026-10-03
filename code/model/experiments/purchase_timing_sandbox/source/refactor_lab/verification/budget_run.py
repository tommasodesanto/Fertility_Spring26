#!/usr/bin/env python3
"""Run ONE command under a wall-clock and resident-memory budget (macOS/Linux, stdlib).

    python3 budget_run.py --seconds 1200 --max-rss-gib 12 --report R.json -- CMD ...

The child runs in its own session; every 2 s the RSS of the whole process
tree (via `ps`) is summed. On wall or memory excess the session gets SIGTERM,
then SIGKILL after 20 s. Exit codes: child's own, 124 wall, 125 memory.
No RLIMIT_AS (it breaks Numba/LLVM JIT mappings on macOS). Peak RSS is recorded.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time


def tree_rss_kib(root: int) -> int:
    out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True).stdout
    kids, rss = {}, {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 3:
            pid, ppid, r = map(int, parts)
            kids.setdefault(ppid, []).append(pid)
            rss[pid] = r
    total, stack = 0, [root]
    while stack:
        pid = stack.pop()
        total += rss.get(pid, 0)
        stack.extend(kids.get(pid, []))
    return total


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seconds", type=float, required=True)
    ap.add_argument("--max-rss-gib", type=float, required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    start = time.monotonic()
    child = subprocess.Popen(cmd, start_new_session=True)
    peak, reason = 0, None
    while child.poll() is None:
        time.sleep(2)
        peak = max(peak, tree_rss_kib(child.pid))
        if time.monotonic() - start > a.seconds:
            reason = "wall"
        elif peak > a.max_rss_gib * 1024 ** 2:
            reason = "memory"
        if reason:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            break
    rc = {"wall": 124, "memory": 125}.get(reason, child.returncode)
    with open(a.report, "w") as f:
        json.dump(dict(cmd=cmd, rc=rc, stop_reason=reason, wall_seconds=time.monotonic() - start,
                       peak_sampled_tree_rss_gib=peak / 1024 ** 2, sampling_interval_seconds=2,
                       memory_cap_kind="sampled every 2 s (not an instantaneous hard ceiling)",
                       cap_seconds=a.seconds, cap_rss_gib=a.max_rss_gib), f, indent=1)
    return rc


if __name__ == "__main__":
    sys.exit(main())
