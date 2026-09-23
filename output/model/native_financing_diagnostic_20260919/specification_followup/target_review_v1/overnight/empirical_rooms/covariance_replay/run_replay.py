"""One bounded exact Stata replay; never retries or changes its inputs."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parent
RECEIPT = ROOT / "execution.json"
SCRIPT = ROOT / "replay_v1.do"
CAP_SECONDS = 1200


def write_receipt(value):
    temporary = RECEIPT.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(RECEIPT)


def main():
    if RECEIPT.exists():
        raise RuntimeError("Replay already attempted; refusing duplicate launch.")
    digest = hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
    review = json.loads((ROOT / "lead_preflight.json").read_text())
    assert digest == review["staged_sha256"]
    start = time.monotonic()
    now = lambda: dt.datetime.now(dt.timezone.utc).isoformat()
    receipt = dict(status="starting", started_utc=now(), supervisor_pid=os.getpid(),
                   script_sha256=digest, cap_seconds=CAP_SECONDS, processors=4,
                   automatic_retry=False)
    write_receipt(receipt)
    with (ROOT / "stata_batch_console.log").open("wb") as log:
        process = subprocess.Popen(
            ["/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp", "-q", "-b", "do", str(SCRIPT)],
            cwd=ROOT / "output", stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        receipt.update(status="running", pid=process.pid)
        while process.poll() is None:
            elapsed = time.monotonic() - start
            receipt.update(heartbeat_utc=now(), elapsed_seconds=round(elapsed, 2))
            write_receipt(receipt)
            if elapsed >= CAP_SECONDS:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                receipt["status"] = "timed_out"
                break
            time.sleep(min(30, CAP_SECONDS - elapsed))
        if receipt["status"] != "timed_out":
            receipt["status"] = "process_completed_pending_scientific_review" if process.returncode == 0 else "process_failed"
        receipt.update(finished_utc=now(), elapsed_seconds=round(time.monotonic()-start, 2),
                       exit_code=process.returncode)
        write_receipt(receipt)


if __name__ == "__main__":
    main()
