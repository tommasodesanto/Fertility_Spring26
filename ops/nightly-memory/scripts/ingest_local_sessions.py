#!/usr/bin/env python3

import argparse
import json
import re
import shutil
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Iterable


@dataclass
class SessionArtifact:
    provider: str
    session_id: str
    title: str
    source_file: Path
    raw_copy: Path
    normalized_copy: Path
    first_timestamp: str | None
    last_timestamp: str | None
    message_count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect local Codex and Claude sessions for a repo/date.")
    parser.add_argument("--date", default=date.today().isoformat(), help="Target date in YYYY-MM-DD.")
    parser.add_argument("--repo-root", required=True, help="Absolute repo root path.")
    parser.add_argument("--output-dir", help="Destination transcript directory. Defaults to repo memory/transcripts/YYYY-MM-DD.")
    parser.add_argument("--codex-home", default=str(Path.home() / ".codex"))
    parser.add_argument("--claude-home", default=str(Path.home() / ".claude"))
    return parser.parse_args()


def read_jsonl(path: Path) -> Iterable[dict]:
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def norm_repo(repo_root: Path) -> str:
    return str(repo_root.resolve())


def claude_project_slug(repo_root: Path) -> str:
    normalized = re.sub(r"[^A-Za-z0-9]+", "-", repo_root.resolve().as_posix().lstrip("/")).strip("-")
    return "-" + normalized


def extract_codex_text(obj: dict) -> tuple[str | None, str | None]:
    payload = obj.get("payload", {})
    if payload.get("type") != "message":
        return None, None
    role = payload.get("role")
    if role not in {"user", "assistant"}:
        return None, None
    parts: list[str] = []
    for item in payload.get("content", []):
        item_type = item.get("type")
        if item_type in {"input_text", "output_text"}:
            text = item.get("text", "").strip()
            if text:
                parts.append(text)
    text = "\n\n".join(parts).strip()
    if not text:
        return None, None
    phase = payload.get("phase")
    if phase and role == "assistant":
        role = f"{role}:{phase}"
    return role, text


def extract_claude_text(obj: dict) -> tuple[str | None, str | None]:
    if obj.get("type") not in {"user", "assistant"}:
        return None, None
    message = obj.get("message", {})
    role = message.get("role")
    if role not in {"user", "assistant"}:
        return None, None
    parts: list[str] = []
    for item in message.get("content", []):
        if isinstance(item, str):
            text = item.strip()
            if text:
                parts.append(text)
            continue
        if not isinstance(item, dict):
            continue
        item_type = item.get("type")
        if item_type == "text":
            text = item.get("text", "").strip()
            if text:
                parts.append(text)
        elif item_type == "tool_use":
            name = item.get("name", "tool")
            tool_input = json.dumps(item.get("input", {}), ensure_ascii=True)
            parts.append(f"[tool_use] {name} {tool_input}")
    text = "\n\n".join(parts).strip()
    if not text:
        return None, None
    return role, text


def write_normalized_markdown(header: str, messages: list[tuple[str, str, str | None]], path: Path) -> None:
    lines = [header, ""]
    if not messages:
        lines.append("_No normalized user/assistant messages found._")
    else:
        for ts, role, text in messages:
            label = role.upper()
            if ts:
                lines.append(f"## {label} [{ts}]")
            else:
                lines.append(f"## {label}")
            lines.append("")
            lines.append(text)
            lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def infer_codex_title(messages: list[tuple[str, str, str | None]]) -> str:
    for _, role, text in messages:
        if role != "user":
            continue
        stripped = text.strip()
        if not stripped:
            continue
        if stripped.startswith("<environment_context>"):
            continue
        return stripped.splitlines()[0][:120]
    return "Codex session"


def collect_codex_sessions(target_date: str, repo_root: Path, codex_home: Path, output_dir: Path) -> list[SessionArtifact]:
    day_dir = codex_home / "sessions" / target_date[0:4] / target_date[5:7] / target_date[8:10]
    artifacts: list[SessionArtifact] = []
    if not day_dir.exists():
        return artifacts

    raw_dir = output_dir / "raw" / "codex"
    norm_dir = output_dir / "normalized" / "codex"
    raw_dir.mkdir(parents=True, exist_ok=True)
    norm_dir.mkdir(parents=True, exist_ok=True)

    for path in sorted(day_dir.glob("*.jsonl")):
        rows = list(read_jsonl(path))
        if not rows:
            continue
        session_meta = rows[0]
        if session_meta.get("type") != "session_meta":
            continue
        cwd = session_meta.get("payload", {}).get("cwd")
        if cwd != norm_repo(repo_root):
            continue

        session_id = session_meta.get("payload", {}).get("id", path.stem)
        messages: list[tuple[str, str, str | None]] = []
        first_ts = None
        last_ts = None
        for row in rows:
            ts = row.get("timestamp")
            role, text = extract_codex_text(row)
            if role and text:
                messages.append((ts, role, text))
                first_ts = first_ts or ts
                last_ts = ts

        raw_copy = raw_dir / path.name
        shutil.copy2(path, raw_copy)
        normalized_copy = norm_dir / f"{path.stem}.md"
        title = infer_codex_title(messages)
        write_normalized_markdown(
            f"# Codex Session `{session_id}`",
            messages,
            normalized_copy,
        )
        artifacts.append(
            SessionArtifact(
                provider="codex",
                session_id=session_id,
                title=title,
                source_file=path,
                raw_copy=raw_copy,
                normalized_copy=normalized_copy,
                first_timestamp=first_ts,
                last_timestamp=last_ts,
                message_count=len(messages),
            )
        )
    return artifacts


def collect_claude_sessions(target_date: str, repo_root: Path, claude_home: Path, output_dir: Path) -> list[SessionArtifact]:
    project_dir = claude_home / "projects" / claude_project_slug(repo_root)
    artifacts: list[SessionArtifact] = []
    if not project_dir.exists():
        return artifacts

    raw_dir = output_dir / "raw" / "claude"
    norm_dir = output_dir / "normalized" / "claude"
    raw_dir.mkdir(parents=True, exist_ok=True)
    norm_dir.mkdir(parents=True, exist_ok=True)

    for path in sorted(project_dir.glob("*.jsonl")):
        rows = list(read_jsonl(path))
        if not rows:
            continue

        title = "Claude session"
        messages: list[tuple[str, str, str | None]] = []
        first_ts = None
        last_ts = None
        matched = False

        for row in rows:
            if row.get("type") == "ai-title":
                title = row.get("aiTitle") or title
            row_cwd = row.get("cwd")
            ts = row.get("timestamp")
            if row_cwd != norm_repo(repo_root):
                continue
            if not ts or not ts.startswith(target_date):
                continue
            role, text = extract_claude_text(row)
            if role and text:
                matched = True
                messages.append((ts, role, text))
                first_ts = first_ts or ts
                last_ts = ts

        if not matched:
            continue

        raw_copy = raw_dir / path.name
        shutil.copy2(path, raw_copy)
        normalized_copy = norm_dir / f"{path.stem}.md"
        write_normalized_markdown(
            f"# Claude Session `{path.stem}`",
            messages,
            normalized_copy,
        )
        artifacts.append(
            SessionArtifact(
                provider="claude",
                session_id=path.stem,
                title=title,
                source_file=path,
                raw_copy=raw_copy,
                normalized_copy=normalized_copy,
                first_timestamp=first_ts,
                last_timestamp=last_ts,
                message_count=len(messages),
            )
        )
    return artifacts


def build_combined_summary(path: Path, artifacts: list[SessionArtifact]) -> None:
    lines = [f"# Combined User/Assistant Transcript Extract", ""]
    if not artifacts:
        lines.append("_No local Codex or Claude sessions matched this repo and date._")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    for artifact in sorted(artifacts, key=lambda item: (item.first_timestamp or "", item.provider, item.session_id)):
        lines.append(f"## {artifact.provider.title()} | {artifact.title}")
        lines.append("")
        lines.append(f"- Session ID: `{artifact.session_id}`")
        if artifact.first_timestamp:
            lines.append(f"- First timestamp: `{artifact.first_timestamp}`")
        if artifact.last_timestamp:
            lines.append(f"- Last timestamp: `{artifact.last_timestamp}`")
        lines.append(f"- Message count: {artifact.message_count}")
        lines.append(f"- Normalized transcript: `{artifact.normalized_copy}`")
        lines.append("")
        snippet = artifact.normalized_copy.read_text(encoding="utf-8")
        body = "\n".join(snippet.splitlines()[2:]).strip()
        lines.append(body if body else "_No extracted text._")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    try:
        datetime.strptime(args.date, "%Y-%m-%d")
    except ValueError as exc:
        raise SystemExit(f"Invalid --date value: {args.date}") from exc

    output_dir = Path(args.output_dir).resolve() if args.output_dir else (repo_root / "memory" / "transcripts" / args.date)
    output_dir.mkdir(parents=True, exist_ok=True)

    codex_home = Path(args.codex_home).expanduser()
    claude_home = Path(args.claude_home).expanduser()

    artifacts = []
    artifacts.extend(collect_codex_sessions(args.date, repo_root, codex_home, output_dir))
    artifacts.extend(collect_claude_sessions(args.date, repo_root, claude_home, output_dir))

    manifest = {
        "date": args.date,
        "repo_root": str(repo_root),
        "generated_at": datetime.now().astimezone().isoformat(),
        "sessions": [
            {
                "provider": item.provider,
                "session_id": item.session_id,
                "title": item.title,
                "source_file": str(item.source_file),
                "raw_copy": str(item.raw_copy),
                "normalized_copy": str(item.normalized_copy),
                "first_timestamp": item.first_timestamp,
                "last_timestamp": item.last_timestamp,
                "message_count": item.message_count,
            }
            for item in sorted(artifacts, key=lambda item: (item.first_timestamp or "", item.provider, item.session_id))
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    build_combined_summary(output_dir / "combined_user_assistant.md", artifacts)

    print(f"Collected {len(artifacts)} sessions into {output_dir}")


if __name__ == "__main__":
    main()
