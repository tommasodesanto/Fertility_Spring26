#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") [--repo-root PATH] [--detached-memory-root PATH]

Creates a repo-local memory/ symlink that points to the detached nightly-memory
store under ~/Library/Application Support/FertilityNightlyMemory/...
EOF
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root_default="$(cd "${script_dir}/../../.." 2>/dev/null && pwd || true)"
repo_root="${AGENT_REPO_ROOT:-${repo_root_default}}"
runtime_base_dir="${HOME}/Library/Application Support/FertilityNightlyMemory"
detached_memory_root=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root)
            [[ $# -ge 2 ]] || usage
            repo_root="$2"
            shift 2
            ;;
        --detached-memory-root)
            [[ $# -ge 2 ]] || usage
            detached_memory_root="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

if [[ -z "${repo_root}" || ! -d "${repo_root}" ]]; then
    echo "repo root not found: ${repo_root}" >&2
    exit 1
fi

if [[ -z "${detached_memory_root}" ]]; then
    project_slug="$(printf '%s' "${repo_root}" | tr '/[:space:]' '-' | tr -cd 'A-Za-z0-9._-')"
    detached_memory_root="${runtime_base_dir}/${project_slug}/memory"
fi

repo_memory_path="${repo_root}/memory"
mkdir -p "${detached_memory_root}"

if [[ -L "${repo_memory_path}" ]]; then
    current_target="$(readlink "${repo_memory_path}")"
    if [[ "${current_target}" == "${detached_memory_root}" ]]; then
        echo "Repo memory already linked: ${repo_memory_path} -> ${detached_memory_root}"
        exit 0
    fi
    echo "repo memory is already a symlink to a different target: ${current_target}" >&2
    exit 1
fi

if [[ -e "${repo_memory_path}" ]]; then
    if [[ -d "${repo_memory_path}" ]]; then
        if command -v rsync >/dev/null 2>&1; then
            rsync -a "${repo_memory_path}/" "${detached_memory_root}/"
        else
            cp -R "${repo_memory_path}/." "${detached_memory_root}/"
        fi
    fi

    backup_root="$(dirname "${detached_memory_root}")/repo-memory-backups"
    timestamp="$(date +%Y%m%d_%H%M%S)"
    backup_path="${backup_root}/memory-${timestamp}"
    mkdir -p "${backup_root}"
    mv "${repo_memory_path}" "${backup_path}"
    echo "Moved existing repo memory to ${backup_path}"
fi

ln -s "${detached_memory_root}" "${repo_memory_path}"
echo "Linked repo memory: ${repo_memory_path} -> ${detached_memory_root}"
