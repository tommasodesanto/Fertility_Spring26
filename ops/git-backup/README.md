# Daily GitHub Backup

This folder contains the local `launchd` job that backs up the active research
repository to GitHub once per day.

The active GitHub backup target is:

```text
tommasodesanto/Fertility_Spring26:clean-main-2026-05-11
```

The job runs `scripts/run_daily_git_backup.sh` at 23:40 local time. It stages
tracked changes plus non-ignored source/documentation files under the active
project folders, creates a dated backup commit if anything changed, and pushes
to the backup branch.

It deliberately refuses to commit obvious raw/generated data extensions and any
single staged file larger than 20MB. Logs are written under:

```text
logs/git-backup/
```

Install or refresh the launch agent with:

```bash
ops/git-backup/scripts/install_launchd.sh
```

Run a manual backup check with:

```bash
ops/git-backup/scripts/run_daily_git_backup.sh
```
