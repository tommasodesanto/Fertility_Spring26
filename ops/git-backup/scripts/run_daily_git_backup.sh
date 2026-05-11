#!/usr/bin/env bash
set -euo pipefail

REPO="${FERTILITY_REPO:-/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26}"
REMOTE="${FERTILITY_GIT_BACKUP_REMOTE:-origin}"
REMOTE_BRANCH="${FERTILITY_GIT_BACKUP_BRANCH:-clean-main-2026-05-11}"
LOCAL_BRANCH="${FERTILITY_GIT_BACKUP_LOCAL_BRANCH:-main}"
LOG_DIR="${FERTILITY_GIT_BACKUP_LOG_DIR:-$REPO/logs/git-backup}"
MAX_FILE_BYTES="${FERTILITY_GIT_BACKUP_MAX_FILE_BYTES:-20971520}"

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/git-backup-$(date +%F).log"
exec >>"$LOG_FILE" 2>&1

echo
echo "=== Git backup start: $(date '+%Y-%m-%d %H:%M:%S %Z') ==="
cd "$REPO"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: $REPO is not a Git repository"
    exit 2
fi

current_branch="$(git branch --show-current)"
if [[ "$current_branch" != "$LOCAL_BRANCH" ]]; then
    echo "ERROR: expected branch $LOCAL_BRANCH, found $current_branch"
    git status -sb
    exit 3
fi

if ! git remote get-url "$REMOTE" >/dev/null 2>&1; then
    echo "ERROR: remote $REMOTE is not configured"
    exit 4
fi

echo "Remote: $(git remote get-url "$REMOTE")"
echo "Target: $REMOTE/$REMOTE_BRANCH"
echo "Status before staging:"
git status -sb

# Stage tracked edits/deletions and non-ignored active source/documentation files.
# Gitignore remains the main guardrail against generated outputs and raw data.
git add -u
git add \
    .gitignore \
    AGENTS.md \
    CLAUDE.md \
    CALIBRATION_STATUS.md \
    SESSION_DIARY.md \
    README.md \
    fertility_dt_latex.code-workspace \
    code \
    docs \
    latex \
    ops

blocked_paths="$(
    git diff --cached --name-only |
    grep -E '(\.dta|\.rds|\.dat|\.dat\.gz|\.zip|\.xlsx|\.pkg|\.aux|\.log|\.fls|\.fdb_latexmk|\.synctex|\.blg|\.bbl)$' || true
)"
if [[ -n "$blocked_paths" ]]; then
    echo "ERROR: refusing to commit blocked generated/raw paths:"
    echo "$blocked_paths"
    exit 5
fi

large_paths="$(
    git diff --cached --name-only |
    while IFS= read -r path; do
        if [[ -f "$path" ]]; then
            size="$(stat -f %z "$path")"
            if [[ "$size" -gt "$MAX_FILE_BYTES" ]]; then
                printf '%s (%s bytes)\n' "$path" "$size"
            fi
        fi
    done
)"
if [[ -n "$large_paths" ]]; then
    echo "ERROR: refusing to commit staged files larger than $MAX_FILE_BYTES bytes:"
    echo "$large_paths"
    exit 6
fi

if git diff --cached --quiet; then
    echo "No staged changes to back up."
    untracked="$(git status --porcelain --untracked-files=all | grep '^??' || true)"
    if [[ -n "$untracked" ]]; then
        echo "Untracked non-ignored paths remain:"
        echo "$untracked"
    fi
    echo "=== Git backup end: $(date '+%Y-%m-%d %H:%M:%S %Z') ==="
    exit 0
fi

echo "Staged change summary:"
git diff --cached --stat

commit_message="Daily project backup $(date '+%Y-%m-%d %H:%M %Z')"
git commit -m "$commit_message"
git push "$REMOTE" "HEAD:$REMOTE_BRANCH"

echo "Status after push:"
git status -sb
echo "=== Git backup end: $(date '+%Y-%m-%d %H:%M:%S %Z') ==="
