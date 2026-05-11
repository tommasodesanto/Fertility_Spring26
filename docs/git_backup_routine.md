# Git Backup Routine

Current GitHub backup branch:

```text
tommasodesanto/Fertility_Spring26:clean-main-2026-05-11
```

The local branch `main` tracks that GitHub branch. The old GitHub `main` branch
has not been overwritten.

Use this routine at the end of each working session:

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
