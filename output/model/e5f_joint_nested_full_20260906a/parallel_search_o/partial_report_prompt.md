# Codex worker task

## Goal
Extend the experimental PDF builder to honestly report PARTIAL policy results, preserving fail-closed validation. Lead will review all changes and generate/deliver PDF.

## Scope
Exclusive write ownership ONLY tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_review.py and sibling test_build_e5f_joint_nested_review.py. There are existing lead changes adding policy_overview() and supplemental chart: PRESERVE and adapt them. Working directory is repository root; experimental files are in tmp worktree. Read exact file before edits.

## Context
Full project startup required. New confirmed fact at 12:27 UTC September7: job17106283 finished calibration verification but overall FAILED because property-tax-2pct-no-rebate could not clear market in2051, residual2.791e-4>2e-4. Baseline,supply-plus-20,dependent-child-ltv95 complete11dates2023–2063. Tax has7 validated dates2023–2047 ONLY. Original overall receipt status partial_policy_failures with three complete cases and tax in failures. Lead preserves source/targets/gates and may diagnose separately. Local actual evidence: output/model/e5f_joint_nested_full_20260906a/parallel_search_o/output/model/joint_nested_overnight/equilibrium_path/. Original selected summary local support_repair_m/smoke/smoke_histories/task_004/ within same experiment output. Final controller search/final_verification.json exists after fresh collection; don't modify.

## Do not touch
No numerical/model/production code, no source receipts/CSVs/graphs or signatures, no narratives/status/other sources, no cluster jobs, no git actions, no protected latex/JMP_DS_draft. Do not certify a completed44date path. Do not waive validation.

## Required output
Add an explicit opt-in --allow-partial-policies flag (default refuses incomplete as before). Under opt-in validate original partial_policy_failures receipt, exact selectedhash/bundle/target/inheritedstate, complete cases and all existing dated budget/probability/value/17graph gates. For failed tax require matching failure.json and strict consecutive prefix dates2023–2047 from policy_path_progress.csv; every reported point passes its dated gates. Ensure cases and failures partition four expected policies with baseline complete. Full policy_effects.csv currently includes only four rows for completed nonbaseline cases; validate those exactly and independently derive tax impact2023 / lastvalid2047 for reporting with explicit incomplete/fails2051 wording. NEVER call2047 a2063 endpoint. Chart may include tax only through2047, no interpolation/extrapolation beyond that. Adjust policy_overview per-policy years to handle prefixes. Full path behavior unchanged; all source checks stay. Report status must explicitly say33 certified full-branch dates plus7 valid tax-prefix dates,40total,4unavailable; no overall pass assertion. Receipt and output metadata retain incomplete status.

## Verification
Existing6 report tests with bundledPython /Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 and --fixture-root /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/exhaustive_smoke_c. Add meaningful partial-policy validation coverage for missing/incorrect failedcase, mismatchedhash/gates, incorrectdateprefix; use copied fixtures and do not alter original evidence. Report concise changed paths, validation results, limitations. The lead installed matplotlib in tmp/pdfs/plot_dependencies because bundled runtime lacks it; use PYTHONPATH pointing there for chart checks if needed. Do NOT render final PDF.

## Stop and report if
Stop after passing focused tests and complete diff, or at wrapper30minute limit; no broad audit or additional delegation. Unresolved scientific ambiguity goes to lead rather than changing the contract.
