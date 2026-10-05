#!/usr/bin/env python3
"""Submit the two separately pinned overnight arrays exactly once on Torch."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ARMS = {
    "a": ("estate_birth_overnight_20261004_v1", "ESTATE_RUN_MODE"),
    "b": ("soft_timing_overnight_20261004_v1", "CONTINUE_RUN_MODE"),
}


def main(arm):
    name, variable = ARMS[arm]
    root = Path("/scratch/td2248/projects") / name
    control = root / "control"
    receipt = control / "production_submission.json"
    with (control / "production_submission.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if receipt.exists():
            raise SystemExit(f"Already submitted: {receipt.read_text()}")
        completed = json.loads(next((root / "results").glob("mock_*/run/completed.json")).read_text())
        terminal = json.loads(next((root / "results").glob("mock_*/launcher_terminal.json")).read_text())
        if not (completed.get("status") == "mock_loop_passed_zero_solves"
                and completed.get("objective_calls") == 2
                and completed.get("lifecycle_solves") == 0
                and terminal.get("exit_code") == 0):
            raise SystemExit("Exact-loop two-case zero-solve mock failed")
        subprocess.run([sys.executable, str(root / "verify_stage.py"), "--host"], check=True,
                       stdout=subprocess.DEVNULL)
        jobs = subprocess.check_output(["squeue", "-u", os.environ["USER"], "-h", "-o", "%j"], text=True).splitlines()
        if any(name in job for job in jobs):
            raise SystemExit(f"Possible duplicate job for {name}")
        command = ["sbatch", "--parsable", "--begin=now+5minutes", "--array=0-11%12",
                   "--time=12:00:00", "--export=ALL," + variable + "=production",
                   str(root / "launch_torch.sh")]
        output = subprocess.check_output(command, text=True).strip()
        jobid = output.split(";")[0]
        data = {"arm": arm, "job_id": jobid, "submitted_utc": datetime.now(timezone.utc).isoformat(),
                "command": command, "mock": {"objective_calls": 2, "lifecycle_solves": 0},
                "no_auto_retry": True}
        tmp = receipt.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2) + "\n")
        tmp.replace(receipt)
        print(json.dumps(data))


if __name__ == "__main__":
    main(sys.argv[1])
