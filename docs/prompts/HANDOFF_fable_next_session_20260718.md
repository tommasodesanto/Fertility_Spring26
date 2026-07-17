# HANDOFF: next Fable session (written 2026-07-17 evening)

Paste this into a fresh session. Mandatory startup still applies (read
`memory/AGENT_MEMORY.md`, latest `memory/daily/`, `CALIBRATION_STATUS.md`
top sections). This file adds the working queue and the decisions already
made — do not re-litigate them.

## Where the project stands

- **M5 is the working calibration** (loss 9.044, established-12 block 3.166,
  best to date; kappa_T estimated interior 0.0100 for the first time; both
  late-life wealth targets hit; the single disclosed miss is old-age
  ownership 0.954 vs 0.764). Canonical artifacts:
  `output/model/intergen_income_disciplined_recalibration_20260716/report/`
  (+ `identification/`: Jacobian rank 9/14 at 1e-2, condition 2926, weak
  direction jointly mixes h_bar_0/theta1/theta0/kappa_T).
- **The calibration-section rewrite** (`latex/quantification_rewrite_review.tex`)
  was circulated to advisors (Guido and Corina) on 2026-07-17 with all M5
  numbers verified, the income-variance provisional note, and the honest
  old-age-ownership discussion. The paper draft's page-1 status, intro,
  policy prose, and Table 5 now carry the M5 policy deltas (grant: births
  +2.1%, price +0.90%; grant+property-tax: +2.7%, −18.4%); the two model
  figures in `latex/figures/` are regenerated from M5 (backups in
  `latex/figures/backup_pre_m5_20260717/`).
- **Income risk remains the flagged provisional item**: SS-range sigma
  (0.12–0.20) is infeasible under the subsistence minima with no safety net
  (verified by probes; census: lowest-persistent-state households cannot
  afford c_bar_0 + rent(h_bar_0) ≈ 32% of mean income). Full solution space
  and falsifications:
  `docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md`
  and the Desktop advice folder `~/Desktop/income_risk_advice_20260717/`
  (files 1–10, incl. the ChatGPT-Pro advice and both probe rounds).
- **Decided, do not reopen without new evidence**: theta_n = 0 external
  (statistically insignificant PSID gap; child-blind literature convention);
  Nb=120 is the project standard — never propose Nb=240 verification;
  deterministic replays/figures/policy solves run LOCALLY (torch is for
  searches); the income note ships as written (feasibility defense, no BPP
  overclaim); estate p90/p50 and the family-size estate gap stay untargeted
  diagnostics.

## Priority queue (in order)

1. **Fold in advisor feedback.** Guido and Corina received the note + the
   week's summary (entry-wealth fix, flat-target swap, old-age ownership
   restored with kappa_T estimated, income-risk impasse with the three
   candidate fixes). Their replies may reorder everything below.
2. **Merge the calibration rewrite into the paper draft.** Tommaso's
   paragraph-approval flow (see the comment header of
   `quantification_rewrite_review.tex`) governs; once he approves the
   remaining paragraphs, replace the draft's old calibration section and
   verify every number and cross-reference survives the merge. The draft's
   status bullet already says calibration work is ongoing.
3. **The income-type robustness (pre-designed, awaiting go).** The
   prevalent-persistent income-type configuration has teeth: probe 3
   (`~/Desktop/income_risk_advice_20260717/10_topstate_probe_results.txt`,
   corrected cells) moved estate p90/50 from 1.75 to 2.74 at FIXED
   parameters with TFR, young liquid wealth, and aggregate ownership all
   moving toward targets. ChatGPT-Pro's advice (spot-verified): run ONE
   externally disciplined recalibration as robustness. Design constraints:
   first-stage estimate the type objects from PSID (education-group mass,
   age–income profiles, within-type rho/sigma — the HSZ-by-education
   template; our pipeline in
   `code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R`
   extends naturally); keep 15 moments / 14 parameters with the type process
   frozen; prefer a 2-type-by-residual grid over one top atom (a 12% atom
   puts p90 inside it — report p85/p90/p95 and masses either side of 10%);
   NEVER identify the type parameters off the estate ratio (circular);
   predeclare which headline policy effects must survive (sign + order of
   magnitude) BEFORE the run. Chain budgets: 8 chains × 3:55 was right for
   M5; reuse. This is a data day + one overnight.
4. **Feasibility appendix** (cheap, ~an hour): the (rho, sigma) frontier —
   dead-node mass and first-failure age — from the existing probe machinery
   (`scratchpad probes are gone; rebuild from the census/probe scripts
   described in the decision memo). It hardens the draft's income footnote
   and is the referee-proofing item the advice memo ranked first.
5. **Old-age exit margin (design memo only, no code)**: the 0.954-vs-0.764
   miss needs a mechanism (health/LTC expense risk à la De Nardi–French–
   Jones is the leading candidate — it would also eventually unlock income
   risk via the same floor infrastructure). Write the spec memo with
   identification bookkeeping; decide after advisor input.
6. **Parked, conditional**: equivalence-scale preference redesign (open only
   if the measured support schedule rejects c_bar_0 or the type robustness
   breaks the headline results); transfer-floor model version (open if the
   paper later makes welfare claims).

## Operational gotchas (this arc's lessons)

- **Torch ssh rides the user's login** (ControlMaster socket,
  `~/.ssh/sockets/`, 10-min persist): check `ssh -o BatchMode=yes torch
  'echo OK'` before promising cluster work; if dead, ask Tommaso to open
  `ssh torch` and keep the window up. The remote repo copy is NOT a git
  clone — deploy changed files individually (rsync -c or ssh cat).
- **Codex CLI is not installed** (`npm install -g @openai/codex` pending) —
  coding tasks route to built-in subagents; objective-critical diffs get
  line-by-line lead verification regardless.
- **Cross-platform reproducibility band**: local solves reproduce Torch
  records only to ~0.5% loss-level (BLAS/SIMD; documented in
  `output/model/intergen_m5_draft_refresh_20260717/gate0/PLATFORM_NOTE.md`).
  Torch records are canonical for reported numbers; local work must compute
  benchmark and counterfactual in the SAME run so the offset cancels.
- **Document-as-you-go discipline** (Tommaso's standing order): commit
  protectively and update `CALIBRATION_STATUS.md` at every launch/collect;
  never leave run machinery untracked; run contracts written before launch
  with all gates machine-encoded (the M4 audit is the cautionary tale:
  `docs/model/intergen_m4_calibration_audit_20260716.md`).
- **Style**: the word "parity" is banned in all writing; concise responses;
  present plans and WAIT for confirmation before new runs or spends; the
  full free/fixed parameter table goes in front of Tommaso before any
  launch — no restriction travels silently inside a phrase.
