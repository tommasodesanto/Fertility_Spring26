# Git Alignment, 2026-05-11

This note records the repository alignment after the May cleanup and model-writeup
work.

The old Git metadata was very large and the working tree had thousands of
generated or archived files. To make future edits recoverable, the repository was
restarted locally with a curated initial commit. The old Git metadata was not
deleted; it was moved aside under `.git_legacy_2026-05-11_70506af`.

The cleaned history is now the active GitHub `main` branch for
`tommasodesanto/Fertility_Spring26`. During the alignment, the cleaned history
was first pushed to `clean-main-2026-05-11` as a non-destructive backup branch;
after user confirmation, `main` was updated to this cleaned history.

The new Git baseline is intended to track:

- project-control files in the root;
- active Python model code under `code/model/`;
- active empirical scripts under `code/empirical/` and selected data scripts;
- active cluster launch and collection scripts;
- active LaTeX sources, bibliography files, and current PDFs under `latex/`;
- compact documentation under `docs/`.

Generated outputs, raw data, large calibration archives, LaTeX build products,
and old archive folders are kept on disk but ignored by default.
