# Git Backup Routine

Current GitHub backup branch:

```text
tommasodesanto/Fertility_Spring26:clean-main-2026-05-11
```

The local branch `main` tracks that GitHub branch. The old GitHub `main` branch
has not been overwritten.

An automatic local `launchd` job is installed from:

```text
ops/git-backup/
```

It runs daily at 23:40 local time. The job stages tracked changes and
non-ignored active source/documentation files, creates a dated backup commit if
anything changed, and pushes to `clean-main-2026-05-11`. It refuses to commit
obvious raw/generated file extensions and any staged file larger than 20MB.

Logs are written to:

```text
logs/git-backup/
```

Manual backup check:

```bash
ops/git-backup/scripts/run_daily_git_backup.sh
```

Install or refresh the daily job:

```bash
ops/git-backup/scripts/install_launchd.sh
```

Even with the automatic job, use this routine at the end of important working
sessions:

```bash
git status -sb
git diff --stat
```

If the changed files are the intended source/documentation changes:

```bash
git add <explicit files>
git commit -m "Short description"
git push
```

Avoid `git add .` when generated outputs or raw data may have appeared. The
`.gitignore` excludes the main generated folders, but explicit staging is still
the safer default.

If a session only produced generated outputs that should not be versioned, leave
the repository clean or document the run in the appropriate status file before
committing.

Do not force-push to the old GitHub `main` branch unless we deliberately decide
to replace the old remote history.
