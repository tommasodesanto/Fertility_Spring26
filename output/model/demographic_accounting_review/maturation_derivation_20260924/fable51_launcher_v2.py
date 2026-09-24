"""One fresh, bounded Fable 5.1 first-party Claude Max review."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import signal
import subprocess
import sys
import time

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
OUT = ROOT / "output/model/demographic_accounting_review/maturation_derivation_20260924"
PROMPT = OUT / "fable51_handoff.md"
CLI = "/Users/tommasodesanto/.local/bin/claude"
MODEL = "claude-fable-5-1"
LIMIT_SECONDS = 1800


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(receipt):
    tmp = OUT / "fable51_receipt_v2.tmp"
    tmp.write_text(json.dumps(receipt, indent=2) + "\n")
    tmp.replace(OUT / "fable51_receipt_v2.json")


def main():
    if not PROMPT.is_file():
        raise SystemExit("Exact handoff prompt is missing")
    for name in ("fable51_receipt_v2.json", "fable51_stream_v2.jsonl", "fable51_final_v2.md"):
        if (OUT / name).exists():
            raise SystemExit(f"Refusing a second attempt: {name} already exists")

    env = os.environ.copy()
    for key in ("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_BASE_URL"):
        env.pop(key, None)
    auth = subprocess.run([CLI, "auth", "status"], env=env, capture_output=True,
                          text=True, check=True)
    info = json.loads(auth.stdout)
    if not all((info.get("loggedIn"), info.get("authMethod") == "claude.ai",
                info.get("apiProvider") == "firstParty",
                info.get("subscriptionType") == "max")):
        raise SystemExit("Claude Max first-party authentication required; no fallback")

    prompt_bytes = PROMPT.read_bytes()
    args = [CLI, "--model", MODEL, "--effort", "max", "--print",
            "--output-format", "stream-json", "--verbose", "--max-turns", "60",
            "--tools", "Read,Glob,Grep", "--allowedTools", "Read,Glob,Grep",
            "--permission-mode", "dontAsk", "--mcp-config", "{\"mcpServers\":{}}",
            "--strict-mcp-config"]
    receipt = {
        "status": "running", "started_utc": now(), "model_requested": MODEL,
        "route": "Claude Code CLI; claude.ai Max; firstParty",
        "prompt_path": str(PROMPT), "prompt_sha256": hashlib.sha256(prompt_bytes).hexdigest(),
        "time_limit_seconds": LIMIT_SECONDS, "max_turns": 60,
        "tools": ["Read", "Glob", "Grep"], "mcp_servers": "disabled",
        "fresh_session": True, "supervisor_pid": os.getpid(),
        "command_redacted": [CLI, "--model", MODEL, "--effort", "max", "--print",
                              "--output-format", "stream-json", "--verbose",
                              "--max-turns", "60", "--tools", "Read,Glob,Grep",
                              "--allowedTools", "Read,Glob,Grep", "--permission-mode",
                              "dontAsk", "--mcp-config", "{\"mcpServers\":{}}", "--strict-mcp-config"]
    }
    save(receipt)
    start = time.monotonic()
    stream_path, stderr_path = OUT / "fable51_stream_v2.jsonl", OUT / "fable51_stderr_v2.txt"
    with stream_path.open("w") as stream, stderr_path.open("w") as err:
        child = subprocess.Popen(args, cwd=ROOT, env=env, stdin=subprocess.PIPE,
                                 stdout=stream, stderr=err, text=True, start_new_session=True)
        receipt["pid"] = child.pid
        save(receipt)
        child.stdin.write(prompt_bytes.decode("utf-8"))
        child.stdin.close()
        while child.poll() is None and time.monotonic() - start < LIMIT_SECONDS:
            receipt.update(heartbeat_utc=now(), elapsed_seconds=round(time.monotonic() - start, 1))
            save(receipt)
            time.sleep(10)
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            receipt["status"] = "timed_out"
    final = None
    for line in stream_path.read_text().splitlines():
        try:
            item = json.loads(line)
        except json.JSONDecodeError:
            continue
        if item.get("type") == "system" and item.get("subtype") == "init":
            receipt.update(session_id=item.get("session_id"), model_actual=item.get("model"))
        if item.get("type") == "result":
            final = item.get("result")
            receipt.update(result_subtype=item.get("subtype"), is_error=item.get("is_error"),
                           usage=item.get("usage"))
    if final:
        (OUT / "fable51_final_v2.md").write_text(final + "\n")
    receipt.update(exit_code=child.returncode, finished_utc=now(),
                   elapsed_seconds=round(time.monotonic() - start, 1))
    if receipt.get("is_error") or receipt.get("model_actual") != MODEL:
        receipt["status"] = "failed"
    elif receipt.get("status") == "running":
        receipt["status"] = "completed" if child.returncode == 0 else "failed"
    save(receipt)
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"] in ("running", "completed") else 1


if __name__ == "__main__":
    sys.exit(main())
