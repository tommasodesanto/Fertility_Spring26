# Prompt: bequest/exit-margin battery, real calibration, and combined note
# (v3, 2026-07-14 late — launch corrections after source/unit verification)

Use this entire document as the task. Do not ask the user to restate context.
Repository: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

## Role

Run the decisive bounded experiment on the bequest block of the one-market
intergenerational housing--fertility model, then write the combined
specification-and-results note. You are continuing work by other agents; do
not relitigate what is settled — verify against the cited artifacts.

## Mandatory startup reading (in order)

1. `memory/AGENT_MEMORY.md`, latest `memory/daily/*.md`, `CALIBRATION_STATUS.md`
   (top July-14 section = the completed bounded diagnostic).
2. `output/model/full_calibration_audit_20260713_claude/HANDOFF_20260714_v2.md`.
3. `output/model/full_calibration_audit_20260713_claude/round3/REPORT.md` and
   `round3/BEQUEST_LIT_MEMO.md`.
4. `docs/model/bequest_specification_memo_20260714.tex` (+ PDF) — to extend.
5. `output/model/intergen_bequest_calibration_exercise_20260714/README.md` —
   the July-14 opt-in switches and results.

## Settled facts (verify, then use)

- Anchor: the clean frontier, tight loss `6.964712360220218`, theta at
  `/scratch/td2248/projects/Fertility_Spring26_20260711_feasibility_recal/tmp_claude_audit_20260713/round3/seeds/sweep_th00.0_thn0.json`
  — chain_a coordinates with `h_bar_0=0.3774141`, `theta0=theta_n=
  tenure_choice_kappa=0` exactly. Verified bit-identical twice. All
  comparisons anchor here.
- Evaluator: search loose (`max_iter_eq=10, tol_eq=1e-4`), reserve the final
  five minutes of every chain, then solve its winner twice and report only
  TIGHT records (`max_iter_eq=40, tol_eq=2.5e-5`); Nb=120, J=17. Never compare
  loose losses across thetas (~0.08 cross-environment spread).
- Old bequest form rejected (n-scaling ~7x penalty; theta1=0.01 necessity
  bequests; no exit margin -> 74+ ownership 0.98 at theta0=0.3).
- July-14 diagnostic: parent-gated normalized glow alone loss-neutral
  (6.9774); linear age-66-to-82 LTV-to-zero rule reached 4.3718 only via an
  ownership cliff (0.860@70, 0.561@74, ~0@78) parked outside the measured
  65-75 window, with no liquid decumulation (74+/62-74 ratio 1.016).
  Schedule rejected; the ~2.6-point upside behind a CORRECT old-age
  borrowing rule motivates this battery.
- Terminology (use in code comments and the note): the schedule below is a
  REDUCED-FORM AGE-DEPENDENT BORROWING RESTRICTION on owners, not a genuine
  amortizing mortgage contract.
- kappa (tenure smoothing) = 0 exactly, everywhere. phi = 0.80 fixed at
  working ages. The h_bar_0 >= 1 decision (D-D) is SEPARATE: do not impose
  that restriction here. Instead, h_bar_0 remains one of the 11 core
  coordinates re-estimated jointly on the declared relaxed domain in every
  arm; the anchor value is only the chain starting point.
- Identification state (do not ignore): the 15-moment system has ~11
  effective directions at the anchor family; the old-wealth median is
  grid-frozen (zero local information); the decumulation ratio is
  diagnostic-only until pinned. Consequence: (theta0, theta1) CANNOT both be
  treated as identified by simply letting both vary. The battery design
  below reflects this.

## Authorized code changes (opt-in switches only; production defaults inert)

Same pattern as the July-14 exercise; each with unit tests + an inert-diff
regression:
1. `theta1` as an overridable parameter of the parent-gated normalized glow
   v(b,n) = 1{n>=1} * theta0 * ((b+theta1)^(1-sigma) - theta1^(1-sigma))/(1-sigma).
2. Equal-division variant v(b,n) = 1{n>=1} * theta0 * n * u_norm(theta1 + b/n)
   behind its own flag.
3. The age-dependent owner borrowing restriction: collateral share
   phi_j = 0.80 for age < 66, linear in age to a TERMINAL FLOOR Lbar at 82.
   Lbar is an EXTERNAL input, never searched against the loss.
4. Expose ownership 74+ and the 74+/62-74 mean liquid-wealth ratio as solver
   statistics if not already exposed (a driver-level version exists in
   `.../full_calibration_audit_20260713_claude/round3/wealth_profile.py`).
5. For any arm estimating theta0: theta0 must be searched on a domain that
   CONTAINS ZERO EXACTLY (softzero transform on [0, 1.5], as used for
   c_bar_0 in the wide domain), NOT a log domain. This is required for the
   "does theta0 return to zero" question to be answerable.

Run the full unittest suite before launching anything.

## Step 0 — external pins (before any calibration)

0a. **Lbar (required).** Pin from an LTV-or-mortgage-balance-BY-AGE profile
    for elderly owners (SCF or HRS; public summary tables acceptable;
    document source, year, and the mapping in the run README). Mortgage
    INCIDENCE (share with any mortgage) is NOT an acceptable pin — it is a
    different object. If no defensible balance/LTV-by-age source is
    attainable tonight, run Lbar in {0.2, 0.4, 0.6} as three EXTERNAL
    variants, report all three, and state that the schedule is unpinned. Do
    not select Lbar by scalar loss.
0b. **theta1 (required for headline arms).** The published
    Nakajima--Telyukova (2017) calibration reports an annual bequest shifter
    of $7,619 in 2000 dollars, not $19,600; their earlier related model had
    $45,714. There is no source-verified one-to-one conversion to this model's
    annual-gross-income normalization. Therefore use theta1=0.25 as the
    predeclared main external scale and 0.50 as its x2 sensitivity, label both
    as external scale variants rather than literature estimates, and never
    select between them by loss. theta1-free arms are diagnostic only.
0c. **Acceptance band (required).** Build the empirical ownership-by-age
    profile for ages 62-82 from the existing ACS/MMS pipeline. Use
    `code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv`:
    household heads in DUE housing, MMS center+periphery, 2012--2023, with a
    predeclared +/-7.5pp model-tolerance envelope in every bin. This envelope
    is an economic/model acceptance tolerance, not a sampling confidence
    interval. Do not use the old provisional band.
0d. **75+ moments (attempt).** Extract ownership 75+ and the 75+/62-74
    wealth ratio with provenance + bootstrap SEs from the existing ACS/PSID
    pipelines if feasible; otherwise keep them as reported diagnostics. Do
    NOT add hard targets without provenance.

## The battery

All arms: NM polish chains from the anchor seed; per-chain budget 230
minutes AND `max_evals=2000` (state both); checkpoint per eval; Torch
cpu_short, account `torch_pr_570_general` (never project 571); smoke each
arm through its full transformed initial simplex first (at least 12--14
exact-loop evaluations, depending on dimension); deterministic tenure;
Nb=120.
**Chain diversification:** chain 1 starts at the anchor unit vector; chain 2
uses `start-mix=0.08` (a documented, materially distinct starting simplex).
Note: RNG seeds in this codebase do jitter the initial simplex steps, but
start-mix is the required diversification device.

| arm | free parameters | fixed by design | question |
|---|---|---|---|
| A0 | the 11 anchor coordinates | theta0=theta1-block off, no schedule | control re-anchor (~6.9647 tight) |
| A1 | 11 + theta0 (softzero, includes 0) | theta1 pinned (0b), NO schedule | do bequests pay before the exit margin? (expect no) |
| A2 | the 11 | no bequests; schedule at pinned Lbar | the borrowing restriction alone |
| A3 (headline) | 11 + theta0 | theta1 pinned; schedule at pinned Lbar | does theta0 come back positive? |
| A4 | 11 + theta0 | as A3, equal-division variant | robustness of the child link |
| A5 (diagnostic only) | 11 + theta0 + theta1 | schedule at pinned Lbar | labeled UNIDENTIFIED; ridge exploration, never a candidate |

Because Lbar remains a 3-variant external profile, run A2/A3 at each variant,
A4 and A5 at the middle one, and the theta1=0.50 A3 sensitivity at the middle
one. This is 11 specifications x 2 chains = 22 main chains. At the observed
~31 seconds per tight solve, the 230-minute time cap implies roughly 440
evaluations per chain; `max_evals=2000` will not bind. The main array is at
most 84.3 CPU-hours and about four hours of concurrent wall time. The
nuisance-reoptimized theta0 profile adds 10 chains, at most 38.3 CPU-hours
and another four hours of wall time. State this 122.6 CPU-hour, two-stage
upper bound before launch; do not exceed it without saying so.

## Identification protocol (mandatory, after the search)

1. **theta0 profile** at the primary A3 configuration: fix theta0 in
   {0, 0.05, 0.1, 0.2, 0.4}, keep theta1=0.25 and Lbar=0.4 externally fixed,
   and RE-OPTIMIZE all eleven nuisance coordinates at every grid value with
   two starts. A slice holding nuisance coordinates at the A3 winner is not a
   profile and is not an inference object. Report the reoptimized curve tight.
2. **Fresh Jacobian** at the A3 winner (pattern:
   `.../tmp_claude_audit_20260713/audit_jacobian.py`, tight evaluator);
   report rank at rel. tol 1e-2/1e-3 and whether theta0's column is above
   the noise floor.
3. State parameter/moment counts per arm explicitly; A5 is presented as
   underidentified by construction.

## Acceptance rules (hard)

- Reject any configuration, regardless of loss, whose ownership-by-age
  profile violates the empirical band (0c) in ANY 4-year bin 62-82, or
  whose ADJACENT-BIN drop exceeds 15 percentage points anywhere in 62-82
  (this generalizes the July-14 cliff test, which was gameable at 70-74).
- Report the 74+/62-74 decumulation ratio for every finalist.
- Tight-evaluator numbers only in tables; determinism repeat (x2) on each
  arm winner; boundary checks on the A3 winner.
- No production promotion. The battery informs the pending decisions
  (exit margin, bequest swap, D-D); it does not take them.

## Deliverable: the combined note

Extend `docs/model/bequest_specification_memo_20260714.tex` into
`docs/model/bequest_exit_note_20260715.tex` (leave the memo untouched):
1. Carry over sections 1-3 (question; evidence with citations and the
   verification footnote; model evidence). Trim, do not re-derive.
2. "Candidate solutions": one subsection per mechanism with 3-5 sentences of
   economics and its literature line — parent gate (Hurd 1989;
   Kopczuk-Lupton 2007); luxury shifter (De Nardi 2004; Lockwood 2018;
   Nakajima-Telyukova 2017); the reduced-form age-dependent borrowing
   restriction (motivated by the July-14 cliff failure; retention evidence);
   equal-division variant (McGarry; Wilhelm; the n^(eps+sigma-1) algebra
   already in the memo); inter-vivos transfers as the child-intensity
   channel (Kvaerner 2023; ledger D7 — future work).
3. "Overnight results": per-arm tight tables (loss, full 15-moment fit,
   parameters/bounds), the acceptance-band outcomes, the theta0 profile
   curve, the Jacobian rank statement, and the verdict on the falsifiable
   prediction: if theta0 still runs to zero in A3 under a sane ownership
   profile, the data genuinely reject bequests in this model.
4. "Recommendation": which configuration to carry toward production, what
   remains unpinned (Lbar source, 75+ moments, theta1 estimation), and the
   explicit list of model changes then needing the author's approval.
Style: `docs/style/econ_writing_style_guide.md`; the word "parity" is
banned; compile standalone; clean aux files.

## Reporting

End with: job IDs and exit states; the per-arm tight table; acceptance
outcomes; theta0-profile curve; note path + compile status; a one-paragraph
append for `CALIBRATION_STATUS.md` (append, never rewrite history). Preserve
raw records on scratch; mirror summaries under
`output/model/intergen_bequest_exit_battery_20260714/`.
