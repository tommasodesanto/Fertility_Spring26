#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") [--link-repo-memory]

Options:
  --link-repo-memory   Replace repo memory/ with a symlink to detached memory.
EOF
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
launch_agents_dir="${HOME}/Library/LaunchAgents"
runtime_base_dir="${HOME}/Library/Application Support/FertilityNightlyMemory"
project_slug="$(printf '%s' "${repo_root}" | tr '/[:space:]' '-' | tr -cd 'A-Za-z0-9._-')"
runtime_dir="${runtime_base_dir}/${project_slug}"
runtime_scripts_dir="${runtime_dir}/scripts"
memory_root="${runtime_dir}/memory"
plist_path="${launch_agents_dir}/com.tommasodesanto.fertility-nightly-memory.plist"
label="com.tommasodesanto.fertility-nightly-memory"
codex_bin="${CODEX_BIN:-$(command -v codex || true)}"
python_bin="${PYTHON_BIN:-$(command -v python3 || true)}"
link_repo_memory=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --link-repo-memory)
            link_repo_memory=1
            shift
            ;;
        *)
            usage
            ;;
    esac
done

if [[ -z "${codex_bin}" ]]; then
    echo "codex CLI not found in PATH; set CODEX_BIN and retry." >&2
    exit 1
fi

if [[ -z "${python_bin}" ]]; then
    echo "python3 not found in PATH; set PYTHON_BIN and retry." >&2
    exit 1
fi

mkdir -p "${launch_agents_dir}"
mkdir -p "${runtime_scripts_dir}"
mkdir -p "${memory_root}/automation/logs"

cp "${script_dir}/run_nightly_memory.sh" "${runtime_scripts_dir}/run_nightly_memory.sh"
cp "${script_dir}/update_memory_scaffold.sh" "${runtime_scripts_dir}/update_memory_scaffold.sh"
cp "${script_dir}/ingest_local_sessions.py" "${runtime_scripts_dir}/ingest_local_sessions.py"
cp "${script_dir}/link_repo_memory.sh" "${runtime_scripts_dir}/link_repo_memory.sh"
chmod 755 \
    "${runtime_scripts_dir}/run_nightly_memory.sh" \
    "${runtime_scripts_dir}/update_memory_scaffold.sh" \
    "${runtime_scripts_dir}/ingest_local_sessions.py" \
    "${runtime_scripts_dir}/link_repo_memory.sh"

cat > "${plist_path}" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key>
  <string>${label}</string>
  <key>EnvironmentVariables</key>
  <dict>
    <key>CODEX_BIN</key>
    <string>${codex_bin}</string>
    <key>PYTHON_BIN</key>
    <string>${python_bin}</string>
  </dict>
  <key>ProgramArguments</key>
  <array>
    <string>/bin/bash</string>
    <string>${runtime_scripts_dir}/run_nightly_memory.sh</string>
    <string>--repo-root</string>
    <string>${repo_root}</string>
    <string>--memory-root</string>
    <string>${memory_root}</string>
    <string>--detached-memory</string>
    <string>--summarize-with-codex</string>
  </array>
  <key>WorkingDirectory</key>
  <string>${runtime_dir}</string>
  <key>StartCalendarInterval</key>
  <dict>
    <key>Hour</key>
    <integer>23</integer>
    <key>Minute</key>
    <integer>50</integer>
  </dict>
  <key>RunAtLoad</key>
  <true/>
  <key>StandardOutPath</key>
  <string>${memory_root}/automation/logs/launchd-stdout.log</string>
  <key>StandardErrorPath</key>
  <string>${memory_root}/automation/logs/launchd-stderr.log</string>
</dict>
</plist>
EOF

launchctl unload "${plist_path}" >/dev/null 2>&1 || true
launchctl load "${plist_path}"
echo "Installed launchd job at ${plist_path}"
echo "Detached runtime dir: ${runtime_dir}"
echo "Automated memory root: ${memory_root}"

if [[ "${link_repo_memory}" -eq 1 ]]; then
    /bin/bash "${script_dir}/link_repo_memory.sh" \
        --repo-root "${repo_root}" \
        --detached-memory-root "${memory_root}"
fi
