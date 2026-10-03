#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") [--date YYYY-MM-DD] [--repo-root PATH] [--memory-root PATH] [--detached-memory] [--summarize-with-codex]" >&2
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root_default="$(cd "${script_dir}/../../.." 2>/dev/null && pwd || true)"
repo_root="${AGENT_REPO_ROOT:-${repo_root_default}}"
memory_root="${AGENT_MEMORY_ROOT:-}"
note_date="$(date +%F)"
summarize_with_codex=0
detached_memory=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --date)
            [[ $# -ge 2 ]] || usage
            note_date="$2"
            shift 2
            ;;
        --repo-root)
            [[ $# -ge 2 ]] || usage
            repo_root="$2"
            shift 2
            ;;
        --memory-root)
            [[ $# -ge 2 ]] || usage
            memory_root="$2"
            shift 2
            ;;
        --detached-memory)
            detached_memory=1
            shift
            ;;
        --summarize-with-codex)
            summarize_with_codex=1
            shift
            ;;
        *)
            usage
            ;;
    esac
done

if [[ -z "${repo_root}" ]]; then
    echo "repo root is required" >&2
    exit 1
fi

if [[ -z "${memory_root}" ]]; then
    memory_root="${repo_root}/memory"
fi

daily_dir="${memory_root}/daily"
transcript_dir="${memory_root}/transcripts/${note_date}"
automation_dir="${memory_root}/automation"
logs_dir="${automation_dir}/logs"
agent_memory_file="${memory_root}/AGENT_MEMORY.md"
daily_file="${daily_dir}/${note_date}.md"
summary_log="${logs_dir}/summary-${note_date}.log"
summary_status_file="${logs_dir}/summary-${note_date}.status"

mkdir -p "${logs_dir}" "${daily_dir}" "${transcript_dir}"

if [[ ! -f "${agent_memory_file}" ]]; then
    cat > "${agent_memory_file}" <<EOF
# Agent Memory

## Project

- Repo root: \`${repo_root}\`
- Automated memory root: \`${memory_root}\`

## Current State

- Fill this file with durable cross-session facts, not raw chat copies.

## Next Session Start Here

- Read the latest note in \`daily/\`
- Use the combined transcript for the day as evidence, not as canonical memory
EOF
fi

python_bin="${PYTHON_BIN:-$(command -v python3 || true)}"
if [[ -z "${python_bin}" ]]; then
    echo "python3 not found; cannot ingest local sessions" >&2
    exit 1
fi

"${python_bin}" "${script_dir}/ingest_local_sessions.py" \
    --date "${note_date}" \
    --repo-root "${repo_root}" \
    --output-dir "${transcript_dir}" \
    > "${logs_dir}/ingest-${note_date}.log" 2>&1

/bin/bash "${script_dir}/update_memory_scaffold.sh" \
    --date "${note_date}" \
    --repo-root "${repo_root}" \
    --memory-root "${memory_root}" \
    > "${logs_dir}/scaffold-${note_date}.log" 2>&1

if [[ "${summarize_with_codex}" -eq 1 ]]; then
    codex_bin="${CODEX_BIN:-$(command -v codex || true)}"
    if [[ -z "${codex_bin}" ]]; then
        echo "codex CLI not found; skipping AI summary step" >> "${summary_log}"
        printf 'SKIPPED: codex CLI not found\n' > "${summary_status_file}"
        exit 0
    fi

    prompt_file="$(mktemp)"
    trap 'rm -f "${prompt_file}"' EXIT
    cat > "${prompt_file}" <<EOF
Refresh project memory for ${note_date} in this repository.

Read:
- ${agent_memory_file}
- ${daily_file}
- ${transcript_dir}/manifest.json
- ${transcript_dir}/combined_user_assistant.md
EOF

    if [[ "${detached_memory}" -eq 0 ]]; then
        cat >> "${prompt_file}" <<EOF
- ${repo_root}/CLAUDE.md
- ${repo_root}/CALIBRATION_STATUS.md
- ${repo_root}/SESSION_DIARY.md
EOF
    fi

    cat >> "${prompt_file}" <<EOF

Update only these files when justified:
- ${daily_file}
- ${agent_memory_file}
EOF

    if [[ "${detached_memory}" -eq 0 ]]; then
        cat >> "${prompt_file}" <<EOF
- ${repo_root}/SESSION_DIARY.md
- ${repo_root}/CALIBRATION_STATUS.md
EOF
    fi

    cat >> "${prompt_file}" <<EOF

Rules:
- Prefer concise durable facts over long transcript copies.
- Do not create new files.
- Keep the generated block in ${daily_file} intact.
EOF

    if [[ "${detached_memory}" -eq 0 ]]; then
        cat >> "${prompt_file}" <<EOF
- Only update SESSION_DIARY.md if the day changed the durable technical story.
- Only update CALIBRATION_STATUS.md if the canonical calibration picture changed.
EOF
    else
        cat >> "${prompt_file}" <<EOF
- This is detached automation outside the repo. Do not assume repo files are readable.
- Focus on durable workflow state and actionable next-session memory.
EOF
    fi

    workdir="${repo_root}"
    if [[ "${detached_memory}" -eq 1 ]]; then
        workdir="${memory_root}"
    fi

    codex_args=(
        exec
        --cd "${workdir}"
        --full-auto
        --color never
        -o "${logs_dir}/codex-last-message-${note_date}.txt"
    )
    if [[ "${detached_memory}" -eq 1 ]]; then
        codex_args+=(--skip-git-repo-check)
    fi
    codex_args+=(-)

    set +e
    "${codex_bin}" "${codex_args[@]}" \
        < "${prompt_file}" \
        > "${summary_log}" 2>&1
    codex_status=$?
    set -e

    if [[ "${codex_status}" -ne 0 ]]; then
        {
            echo
            echo "[nightly-memory] Codex summary failed with exit code ${codex_status}."
            echo "[nightly-memory] Ingestion and scaffold may still have succeeded."
            echo "[nightly-memory] Check Codex auth. If the log shows token expiry or refresh-token reuse, run:"
            echo "[nightly-memory]   codex logout"
            echo "[nightly-memory]   codex login"
            echo "[nightly-memory] Then rerun:"
            echo "[nightly-memory]   ops/nightly-memory/scripts/install_launchd.sh --link-repo-memory"
        } >> "${summary_log}"
        printf 'FAILED: codex summary exit %s\n' "${codex_status}" > "${summary_status_file}"
        exit "${codex_status}"
    fi

    printf 'OK\n' > "${summary_status_file}"
fi
