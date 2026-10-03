# Portable review bundle — paused, not verified for delivery

Tommaso paused this work on October 3, 2026, to resume at home. No package has
been sent. Do not present an existing ZIP or staging folder here as ready.

The scoped worker created `code/model/tools/export_review_bundle.py` and bundle
staging under this directory. Its first clean-unzip native check failed during
dependency authentication: the exported tree omitted
`tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/production_contract.json`.
The worker added a dependency change and began rebuilding; it and its subprocess
were stopped on the author's pause. That change and rebuild are unverified.

Resume from `work/worker_task.md`, `work/worker.log`, `build.log`, and
`verification_fresh_ge.log`; `verification_unzip_root.txt` names the retained
temporary extraction. Inspect the exporter diff before running it. Verify a
complete rebuild, clean extraction, no external links or required home-directory
reads, actual one-core GE, full fit/parameter/array comparison, plots and explorer.
Only then mark the ZIP suitable to send. Keep production sources unchanged.

The intended package includes the runnable model, derived inputs, editable
parameters, clear setup instructions and plotting/explorer tools, with current
net-estate birth-menu experiments identified separately from the older default.
Exclude credentials, raw survey microdata, environment directories and unrelated
research history. The worker's numerical-check budget was one GE, 1,200 seconds,
32 lifecycle solves, on one local core; the failed import did not establish a GE.
