#!/usr/bin/env bash
# Safe, non-interactive Codex worker launcher.  It defaults to read-only.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "$script_dir/../../.." && pwd -P)"
config_file="$repo_root/ops/codex-workers/config/models.env"

usage() {
    cat <<'EOF'
Usage:
  codex-worker.sh --profile NAME --prompt-file PATH [options]

Run a named Codex worker from a prompt file. The worker is read-only by
default, runs ephemerally, and receives access only to the selected repository
workdir.

Required:
  --profile NAME          explorer_fast, mechanic_fast, worker_fast, or reviewer_strong
  --prompt-file PATH      Markdown or text file containing the task prompt

Options:
  --write                 Allow workspace writes for this invocation
  --workdir PATH          Working directory inside this repository (default: repository root)
  --final-message PATH    Write Codex's final message to a file inside this repository
  --time-limit-min N      Stop the worker after N minutes (default: profile limit)
  --full-context          Override a minimal profile and require full project startup
  --dry-run               Print the resolved route without invoking Codex
  -h, --help              Show this help text

Profiles:
  explorer_fast   Terra, low reasoning, read-only by default
  mechanic_fast   Codex Spark, low reasoning, read-only by default
  worker_fast     Terra, medium reasoning, read-only by default
  reviewer_strong Terra, high reasoning, read-only by default

The wrapper always uses --ephemeral, passes an explicit model and
model_reasoning_effort, and instructs the worker never to commit or push.
EOF
}

fail() {
    printf 'codex-worker.sh: %s\n' "$*" >&2
    exit 2
}

path_is_within_repo() {
    local path="$1"
    [[ "$path" == "$repo_root" || "$path" == "$repo_root/"* ]]
}

profile=""
prompt_file=""
workdir="$repo_root"
final_message=""
sandbox="read-only"
dry_run=0
force_full_context=0
time_limit_min=""

while (($#)); do
    case "$1" in
        --profile)
            (($# >= 2)) || fail '--profile requires a value'
            profile="$2"
            shift 2
            ;;
        --prompt-file)
            (($# >= 2)) || fail '--prompt-file requires a path'
            prompt_file="$2"
            shift 2
            ;;
        --workdir)
            (($# >= 2)) || fail '--workdir requires a path'
            workdir="$2"
            shift 2
            ;;
        --final-message)
            (($# >= 2)) || fail '--final-message requires a path'
            final_message="$2"
            shift 2
            ;;
        --write)
            sandbox="workspace-write"
            shift
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        --full-context)
            force_full_context=1
            shift
            ;;
        --time-limit-min)
            (($# >= 2)) || fail '--time-limit-min requires a positive integer'
            time_limit_min="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "unknown option: $1"
            ;;
    esac
done

[[ -n "$profile" ]] || fail '--profile is required'
[[ -n "$prompt_file" ]] || fail '--prompt-file is required'
[[ -f "$prompt_file" && -r "$prompt_file" ]] || fail "prompt file is not readable: $prompt_file"
[[ -f "$config_file" ]] || fail "missing profile configuration: $config_file"

if ! workdir="$(cd "$workdir" 2>/dev/null && pwd -P)"; then
    fail "workdir is not an accessible directory: $workdir"
fi
path_is_within_repo "$workdir" || fail "workdir must be inside repository root: $repo_root"

if [[ -n "$final_message" ]]; then
    final_parent="$(dirname "$final_message")"
    if ! final_parent="$(cd "$final_parent" 2>/dev/null && pwd -P)"; then
        fail "final-message parent is not an accessible directory: $(dirname "$final_message")"
    fi
    path_is_within_repo "$final_parent" || fail "final-message must be inside repository root: $repo_root"
    final_message="$final_parent/$(basename "$final_message")"
fi

# shellcheck source=/dev/null
source "$config_file"

case "$profile" in
    explorer_fast)
        model="$EXPLORER_FAST_MODEL"
        reasoning_effort="$EXPLORER_FAST_REASONING_EFFORT"
        context_mode="$EXPLORER_FAST_CONTEXT_MODE"
        profile_time_limit_min="$EXPLORER_FAST_TIME_LIMIT_MIN"
        ;;
    mechanic_fast)
        model="$MECHANIC_FAST_MODEL"
        reasoning_effort="$MECHANIC_FAST_REASONING_EFFORT"
        context_mode="$MECHANIC_FAST_CONTEXT_MODE"
        profile_time_limit_min="$MECHANIC_FAST_TIME_LIMIT_MIN"
        ;;
    worker_fast)
        model="$WORKER_FAST_MODEL"
        reasoning_effort="$WORKER_FAST_REASONING_EFFORT"
        context_mode="$WORKER_FAST_CONTEXT_MODE"
        profile_time_limit_min="$WORKER_FAST_TIME_LIMIT_MIN"
        ;;
    reviewer_strong)
        model="$REVIEWER_STRONG_MODEL"
        reasoning_effort="$REVIEWER_STRONG_REASONING_EFFORT"
        context_mode="$REVIEWER_STRONG_CONTEXT_MODE"
        profile_time_limit_min="$REVIEWER_STRONG_TIME_LIMIT_MIN"
        ;;
    *)
        fail "unknown profile: $profile"
        ;;
esac

if [[ -z "$time_limit_min" ]]; then
    time_limit_min="$profile_time_limit_min"
fi
[[ "$time_limit_min" =~ ^[1-9][0-9]*$ ]] || fail '--time-limit-min must be a positive integer'

if [[ "$sandbox" == "workspace-write" && "$profile" != "worker_fast" && "$profile" != "mechanic_fast" ]]; then
    fail "--write is allowed only for worker_fast or mechanic_fast; $profile is read-only"
fi

((force_full_context)) && context_mode="full"

launch_mode="repository-full"
launch_dir="$workdir"
if [[ "$context_mode" == "minimal" ]]; then
    launch_mode="detached-minimal"
    if ((dry_run)); then
        launch_dir="<temporary-directory>"
    else
        launch_dir="$(mktemp -d "${TMPDIR:-/tmp}/codex-worker.XXXXXX")"
        cleanup() {
            rm -rf -- "$launch_dir"
        }
        trap cleanup EXIT
    fi
fi

codex_args=(
    exec
    --ephemeral
    --sandbox "$sandbox"
    --model "$model"
    --config "model_reasoning_effort=\"$reasoning_effort\""
    --color never
)

if [[ "$context_mode" == "minimal" ]]; then
    codex_args+=(--skip-git-repo-check --cd "$launch_dir" --add-dir "$workdir")
else
    codex_args+=(--cd "$workdir")
fi

if [[ -n "$final_message" ]]; then
    codex_args+=(--output-last-message "$final_message")
fi

if ((dry_run)); then
    printf 'profile=%s\nmodel=%s\nreasoning_effort=%s\ncontext_mode=%s\nlaunch_mode=%s\nsandbox=%s\ntime_limit_min=%s\nworkdir=%s\n' \
        "$profile" "$model" "$reasoning_effort" "$context_mode" "$launch_mode" "$sandbox" "$time_limit_min" "$workdir"
    [[ -z "$final_message" ]] || printf 'final_message=%s\n' "$final_message"
    exit 0
fi

command -v codex >/dev/null 2>&1 || fail 'codex executable not found on PATH'

{
    cat -- "$prompt_file"
    printf '\n\n## Wrapper safety constraints\n'
    printf '%s\n' 'Never run git commit, git push, git reset, or any command that publishes or changes Git history.'
    printf '%s\n' 'Stay within the task scope. Stop and report if the task requires an action outside the stated scope.'
    printf 'The repository root is `%s`; the authorized task workdir is `%s`.\n' "$repo_root" "$workdir"
    if [[ "$context_mode" == "minimal" ]]; then
        printf '%s\n' 'Context mode is minimal under the repository worker exception. The process starts in an empty temporary directory so project instructions are not loaded automatically. Run repository commands from the authorized task workdir named above.'
        printf '%s\n' 'Do not load AGENTS.md, CLAUDE.md, memory files, daily notes, CALIBRATION_STATUS.md, or unrelated project context. If the task requires economics, calibration, targets, model results, or model-code interpretation, stop and request a full-context route.'
    else
        printf '%s\n' 'Context mode is full: follow the repository Mandatory Startup instructions before substantive work.'
    fi
} | perl -e '
    my $limit = shift @ARGV;
    my $pid = fork();
    die "codex-worker.sh: cannot fork for time limit: $!\\n" unless defined $pid;
    if ($pid == 0) {
        exec @ARGV;
        die "codex-worker.sh: cannot launch worker: $!\\n";
    }
    $SIG{ALRM} = sub {
        warn "codex-worker.sh: time limit of ${limit} minutes reached; terminating worker\\n";
        kill "TERM", $pid;
        sleep 5;
        kill "KILL", $pid;
        waitpid $pid, 0;
        exit 124;
    };
    alarm($limit * 60);
    waitpid $pid, 0;
    my $status = $?;
    alarm 0;
    exit($status >> 8);
' "$time_limit_min" codex "${codex_args[@]}" -
