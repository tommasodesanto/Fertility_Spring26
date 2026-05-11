# Git Backup Routine

Current GitHub backup branch:

```text
tommasodesanto/Fertility_Spring26:main
```

The local branch `main` tracks GitHub `main`. The cleaned repository history was
put on GitHub `main` on 2026-05-11 after the user requested that the current
work live on `main`.

An automatic local `launchd` job is installed from:

```text
ops/git-backup/
```

It runs daily at 23:40 local time. The job stages tracked changes and
non-ignored active source/documentation files, creates a dated backup commit if
anything changed, and pushes to `main`. It refuses to commit
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

Do not force-push GitHub `main` again unless the user explicitly requests a
history rewrite.
