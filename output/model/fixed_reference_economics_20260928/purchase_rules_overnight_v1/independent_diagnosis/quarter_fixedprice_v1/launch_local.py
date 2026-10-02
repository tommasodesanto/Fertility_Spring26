"""One-core local fallback with a 15-minute wall and 24-GiB address cap."""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path
import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
FACTORIAL = "--price-factorial" in sys.argv
SHAPLEY = "--price-shapley" in sys.argv
if FACTORIAL and SHAPLEY:
    raise ValueError("Choose one branch")
OUT = HERE / ("local_run/shapley" if SHAPLEY else
              "local_run/factorial" if FACTORIAL else "local_run/production")
OUT.mkdir(parents=True, exist_ok=True)

env = os.environ.copy()
for name in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    env[name] = "1"
env["PYTHONDONTWRITEBYTECODE"] = "1"
env["MPLCONFIGDIR"] = "/private/tmp/quarter_fixedprice_mpl"
cmd = [str(ROOT / "code/model/.venv/bin/python"), str(HERE / "run.py"),
       "--out", str(OUT)]
if FACTORIAL:
    cmd.append("--price-factorial")
if SHAPLEY:
    cmd.append("--price-shapley")
with (OUT / "run.log").open("w") as log:
    proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    (OUT / "pid.txt").write_text(str(proc.pid) + "\n")
    started = time.monotonic()
    high_water = 0
    while proc.poll() is None:
        if time.monotonic() - started > 900:
            proc.kill()
            proc.wait()
            (OUT / "launch_status.txt").write_text("failed: 900-second wall limit\n")
            raise TimeoutError("900-second wall limit")
        try:
            process = psutil.Process(proc.pid)
            rss = process.memory_info().rss
            high_water = max(high_water, rss)
        except psutil.NoSuchProcess:
            rss = 0
        if rss > 24 * 1024**3:
            proc.kill()
            proc.wait()
            (OUT / "launch_status.txt").write_text("failed: 24-GiB resident memory limit\n")
            raise MemoryError("24-GiB resident memory limit")
        time.sleep(2)
    if proc.returncode != 0:
        (OUT / "launch_status.txt").write_text(f"failed: exit {proc.returncode}\n")
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    (OUT / "launch_status.txt").write_text(
        f"completed: exit 0; peak observed RSS {high_water / 1024**3:.3f} GiB\n")
