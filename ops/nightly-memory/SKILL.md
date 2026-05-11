---
name: nightly-memory-refresh
description: Use when the user wants to roll a day's work into durable repo memory. Supports automated nightly ingestion of local Codex and Claude JSONL sessions, then refreshes memory/daily/YYYY-MM-DD.md and memory/AGENT_MEMORY.md.
---

# Nightly Memory Refresh

Use this workflow at the end of a work session or when the user asks to rebuild project memory from saved chats.

## Inputs

- Core project context: `CLAUDE.md`, `CALIBRATION_STATUS.md`, `SESSION_DIARY.md`
- Today's daily note: `memory/daily/YYYY-MM-DD.md`
- Automated transcript capture: `ops/nightly-memory/scripts/ingest_local_sessions.py`
- Nightly runner: `ops/nightly-memory/scripts/run_nightly_memory.sh`
- Optional raw transcript files or copied captures: `memory/transcripts/YYYY-MM-DD/`
- Repo diff and status from `ops/nightly-memory/scripts/update_memory_scaffold.sh`
- Detached macOS automation runtime: `~/Library/Application Support/FertilityNightlyMemory/<project-slug>/`

## Procedure

1. Prefer `ops/nightly-memory/scripts/run_nightly_memory.sh --summarize-with-codex`.
2. For fully automated macOS runs, use `ops/nightly-memory/scripts/install_launchd.sh`. It deploys a detached runtime outside `Desktop/...` and ingests both `~/.codex/sessions/...` and `~/.claude/projects/...` every night.
3. If you want repo-local visibility while keeping detached automation reliable, run `ops/nightly-memory/scripts/install_launchd.sh --link-repo-memory`. This replaces repo `memory/` with a symlink to detached memory.
4. If only doing a manual refresh, run `ops/nightly-memory/scripts/update_memory_scaffold.sh`.
5. Read the generated note for today in `memory/daily/YYYY-MM-DD.md` or the detached memory root.
6. If transcript files exist, extract only high-signal facts:
   - decisions made
   - parameter changes
   - diagnostics run
   - unresolved issues
   - next actions
7. Update `memory/daily/YYYY-MM-DD.md`:
   - complete `Summary`, `Decisions`, `Open Questions`, and `Next Actions`
   - keep the generated block intact
8. Update `memory/AGENT_MEMORY.md`:
   - overwrite stale status instead of appending duplicate history
   - keep it short enough to read at the start of a session
9. Update `SESSION_DIARY.md` only when the day produced durable solver, calibration, or economics findings.
10. Update `CALIBRATION_STATUS.md` only when the canonical calibration picture changed.

## Rules

- Do not paste full raw chats into tracked memory files.
- Keep raw transcripts under `memory/transcripts/`.
- Use exact dates.
- Prefer replacing stale bullets over appending near-duplicates.
- Treat generated logs and outputs as evidence, not as canonical project memory.
- Claude is ingested from local JSONL under `~/.claude/projects/...`.
- Codex is ingested from local JSONL under `~/.codex/sessions/...`.
- The detached launchd runtime writes nightly memory outside the repo to avoid macOS privacy restrictions on background jobs reading `Desktop/...`.
