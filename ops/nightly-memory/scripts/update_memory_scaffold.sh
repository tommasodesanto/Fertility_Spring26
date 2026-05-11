#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") [--date YYYY-MM-DD] [--repo-root PATH] [--memory-root PATH]" >&2
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root_default="$(cd "${script_dir}/../../.." 2>/dev/null && pwd || true)"
repo_root="${AGENT_REPO_ROOT:-${repo_root_default}}"
memory_root="${AGENT_MEMORY_ROOT:-}"
note_date="$(date +%F)"

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
        *)
            usage
            ;;
    esac
done

if [[ -z "${memory_root}" ]]; then
    [[ -n "${repo_root}" ]] || usage
    memory_root="${repo_root}/memory"
fi

daily_dir="${memory_root}/daily"
transcript_root="${memory_root}/transcripts"
transcript_dir="${transcript_root}/${note_date}"
daily_file="${daily_dir}/${note_date}.md"

mkdir -p "${daily_dir}" "${transcript_dir}"

branch="unknown"
commit="unknown"
generated_at="$(date '+%Y-%m-%d %H:%M:%S %Z')"
git_status="(repo snapshot unavailable)"
status_display_file="$(mktemp)"

if [[ -n "${repo_root}" && -d "${repo_root}" ]]; then
    branch="$(git -C "${repo_root}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")"
    commit="$(git -C "${repo_root}" rev-parse --short HEAD 2>/dev/null || echo "unknown")"
    git_status="$(git -C "${repo_root}" status --short 2>/dev/null || true)"
    if [[ -z "${git_status}" ]]; then
        git_status="(clean working tree)"
    fi
fi

transcript_lines_file="$(mktemp)"
payload_file="$(mktemp)"
tmp_file="$(mktemp)"
trap 'rm -f "${transcript_lines_file}" "${payload_file}" "${tmp_file}" "${status_display_file}"' EXIT

if find "${transcript_dir}" -mindepth 1 -maxdepth 1 | grep -q .; then
    while IFS= read -r transcript_file; do
        entry_name="$(basename "${transcript_file}")"
        if [[ -d "${transcript_file}" ]]; then
            entry_name="${entry_name}/"
        fi
        printf -- "- %s\n" "${entry_name}" >> "${transcript_lines_file}"
    done < <(find "${transcript_dir}" -mindepth 1 -maxdepth 1 | LC_ALL=C sort)
else
    printf -- "- None yet. The nightly ingester will create files here when it finds local sessions.\n" >> "${transcript_lines_file}"
fi

if [[ "${git_status}" == "(clean working tree)" || "${git_status}" == "(repo snapshot unavailable)" ]]; then
    status_count=0
    printf '%s\n' "${git_status}" > "${status_display_file}"
else
    status_count="$(printf '%s\n' "${git_status}" | awk 'NF {count++} END {print count+0}')"
    printf '%s\n' "${git_status}" | sed -n '1,120p' > "${status_display_file}"
    if [[ "${status_count}" -gt 120 ]]; then
        printf '... (%s total status entries; truncated to first 120)\n' "${status_count}" >> "${status_display_file}"
    fi
fi

{
    echo "<!-- BEGIN GENERATED -->"
    echo "## Generated Snapshot"
    echo "- Generated at: ${generated_at}"
    echo "- Git branch: \`${branch}\`"
    echo "- HEAD commit: \`${commit}\`"
    echo "- Git status entries: ${status_count}"
    echo "- Transcript folder: \`${transcript_dir}\`"
    echo
    echo "### Raw Transcript Files"
    cat "${transcript_lines_file}"
    echo
    echo "### Git Status Snapshot"
    echo '```text'
    cat "${status_display_file}"
    echo '```'
    echo "<!-- END GENERATED -->"
} > "${payload_file}"

if [[ ! -f "${daily_file}" ]]; then
    cat > "${daily_file}" <<EOF
# Nightly Memory ${note_date}

## Summary

- Fill this in after reviewing today's work, diffs, and any saved chat transcripts.

## Decisions

- 

## Open Questions

- 

## Next Actions

- 

EOF
fi

if grep -q "<!-- BEGIN GENERATED -->" "${daily_file}"; then
    awk -v payload_file="${payload_file}" '
        BEGIN {
            while ((getline line < payload_file) > 0) {
                payload = payload line ORS
            }
            close(payload_file)
            inside = 0
        }
        /<!-- BEGIN GENERATED -->/ {
            printf "%s", payload
            inside = 1
            next
        }
        /<!-- END GENERATED -->/ {
            inside = 0
            next
        }
        !inside { print }
    ' "${daily_file}" > "${tmp_file}"
    mv "${tmp_file}" "${daily_file}"
else
    printf '\n' >> "${daily_file}"
    cat "${payload_file}" >> "${daily_file}"
    printf '\n' >> "${daily_file}"
fi

echo "Updated ${daily_file}"
echo "Transcript folder: ${transcript_dir}"
