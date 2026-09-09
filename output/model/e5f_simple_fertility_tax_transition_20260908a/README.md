# Rebated property-tax paths through 2063

Purpose: establish persistence of the new calibration's tax-policy fertility and housing responses, with special attention to young parents. Annual 1% versus 2%, equal rebates in both paths. No calibration or production promotion.

Both paths begin with the same verified 2023 inherited population, selected parameters, historical queue and dated supply rule (elasticity 0.63). Preferences stay at the selected 2023 values. Closed-national diagnostic: no outside entries, full retention, adjusted births converted to household entrants using 1/2.1 with the maintained twenty-year lag. No future Census reweighting, person-demography replacement, supply reanchoring or fiscal changes. This is a path of household-unit temporary equilibria, not a resident population forecast or perfect-foresight transition.

Outstanding assumptions remain visible: household formation/headship and the transition from historical empirical age weights to endogenous entrants are unresolved. Record entry cohorts and queue explicitly. Policy-induced household-mass differences must be zero before 2043.

## Run design and stop rules

Two parallel smoke paths through 2023 and 2027 use the exact intended solve/replay/audit/advance loop. Each2023 endpoint must reproduce the independently verified tax equilibrium. Only after BOTH smokes pass, run the two full paths at eleven dates (2023–2063). Full paths repeat their own 2023/2027 smoke checkpoints and quantities. A final collector requires both full paths.

26 coupled equilibrium solves plus26 fresh fixed-price replays, including four smoke dates and22 full dates. Observed tax endpoint time:2m28s–2m45s, peak memory6.49–6.79GiB. Expected elapsed compute roughly35–45minutes with the two paths parallel, allowing extra checkpoint/advance costs; queue time additional. Each path reserves two CPUs and48GiB, with one numerical thread and memory headroom for harder roots. Each date has30minutes; each stage has4hours and a4h05 Slurm cap. Failed numerical gates stop that path and dependent jobs cancel; no retries, threshold relaxation or parameter changes.

Every date restores the selected fiscal transfer after the cached root, independently checks the ledger, performs a fresh full-policy replay, checks distribution/budget/feasibility/probabilities/value monotonicity, advances the existing cohort law and verifies mass/queue accounting. Save the complete dated state, next state, lifecycle/family tables, heartbeat every30seconds, latest completed date and selected-calibration summary. No new figures or monitoring automation.

The morning report compares births per household AND total births, household mass, entrants, rooms, ownership, prices and rebates, including young/dependent-child groups. If a path fails, retain and label the completed prefix. Complete calibration fit and parameter tables remain in `../e5f_simple_fertility_overnight_20260907a/morning_review/MORNING_REVIEW.md`.

Frozen contract SHA256: `ce80b6ad241bec9556d2f3d3cdccb9f89e0ccfdd5b2e80ddcd167ed5156b5b68`. Model/source/closure and both verified impact receipts are pinned in contract.json.

## Submitted jobs

Smoke array **17250629** runs both tax paths through2027. Full array **17250630** depends on both smoke tasks succeeding; collector **17250631** depends on both full paths. Invalid dependencies cancel automatically. Source `e483254a` on isolated `codex/fertility-nest-computation`. Compilation, CLI import and nine pure queue/policy/failure checks pass; lead reviewed the numerical loop and independent review checked cohort routing. All source/contract/endpoint artifact hashes and twenty-year lag verified on Torch. Actual two-date numerical smoke pending at submission.

Login-node `/tmp` was full during a harmless font-cache attempt; scratch has ample capacity. Job temporary files and Matplotlib cache explicitly use scratch. No user data was deleted.
