#!/usr/bin/env bash
set -euo pipefail

LABEL="com.tommasodesanto.fertility-git-backup"
REPO="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
PLIST_SRC="$REPO/ops/git-backup/$LABEL.plist"
PLIST_DEST="$HOME/Library/LaunchAgents/$LABEL.plist"
LOG_DIR="$REPO/logs/git-backup"

mkdir -p "$HOME/Library/LaunchAgents" "$LOG_DIR"
cp "$PLIST_SRC" "$PLIST_DEST"

launchctl bootout "gui/$(id -u)" "$PLIST_DEST" >/dev/null 2>&1 || true
launchctl unload "$PLIST_DEST" >/dev/null 2>&1 || true

if ! launchctl bootstrap "gui/$(id -u)" "$PLIST_DEST" >/dev/null 2>&1; then
    launchctl load "$PLIST_DEST"
fi

launchctl enable "gui/$(id -u)/$LABEL" >/dev/null 2>&1 || true

echo "Installed $LABEL"
echo "Schedule: daily at 23:40 local time"
echo "Log directory: $LOG_DIR"
