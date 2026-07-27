# Handoff: E-strand calibration hardening (author away ~10 days)

You are taking over the E-strand (equivalence-scale, sequential-fertility)
calibration of the housing-fertility lifecycle model in this repository,
working autonomously for about a week while the author is unreachable. This
document is your contract. Where it conflicts with your defaults, it wins;
where the repository's `AGENTS.md` is stricter, `AGENTS.md` wins.

## 0. Mandatory startup

Read, in this order, before doing anything:

1. `AGENTS.md` (root) — house rules, verification, delegation, cluster safety.
2. `memory/AGENT_MEMORY.md` and the latest `memory/daily/`.
3. `CALIBRATION_STATUS.md`.
4. `docs/model/e5_target_review_20260724.md` — the target-system ledger.
   Every decision of the last week is registered there, including the E5b
   certification, the probe verdicts, and four open items.
5. `docs/model/e5b_calibration_autopsy_20260726.md` — the diagnosis you are
   acting on.
6. `docs/model/entry_postponement_mechanism_memo_20260725.md` — the
   mechanism menu and its evidence.

Scope boundary: you own `code/model/intergen_eqscale_seq_optimized/` and its
outputs. You must NOT touch `intergen_housing_fertility/` or
`intergen_housing_fertility_optimized/` (the M strand, owned by another
agent), their calibration contracts, or `latex/` (the author edits the note
and paper personally). `code/data/cps_fertility/` is M-owned; read only.

## 1. Where things stand

- Operative winner: E5b, certified 2026-07-26, loss 385.1 on the signed
  twelve-row / ten-parameter system
  (`output/model/eqscale_seq_e5b_recalibration_20260725/report/`). Eight of
  twelve rows fit within ~1.3 SE. Winner: annual beta 0.9911, psi 0.585,
  kappa_E 3.552, kappa_C 1.254, chi 1.038, H0 8.162, theta0 0.168, theta1
  0.050, delta_jump 0.0239, delta_a 0.0122. All interior.
- The four misses and their diagnosis (autopsy, verified at machine zero
  against the certified winner):
  a. Childlessness 0.080 vs 0.188 — the model has no permanent
     heterogeneity; childlessness is flat across income states while the
     CPS gradient runs 0.068 (<HS) to 0.245 (advanced degree)
     (`code/data/cps_childlessness_gradient/output/`).
  b. Wealth p90/p50 (76-84) 2.07 vs 3.45 and aggregate wealth/earnings
     5.45 vs 6.87 — model cross-sectional income p90/p50 is 1.61 vs 3.77
     measured in the PSID
     (`code/data/psid_followup_mar2026/output/psid_income_between_within.csv`).
     The model's AR(1) matches the measured WITHIN component (1.94); the
     missing dispersion is the BETWEEN layer. Caveat registered: the naive
     70 percent between share is contaminated by slow mean reversion
     (within deviations show 2-year autocorrelation 0.04 where the FL
     component implies 0.84); the fixed-vs-persistent split needs a proper
     minimum-distance decomposition.
  c. Share of first births 30+ 0.332 vs 0.270 — the excess is ALL at ages
     38+ (9.7 vs 3.0 percent of first births), priced by the unfitted
     fecundity tail (pi = 0.50 at 42 vs ~97 percent permanently infertile
     by 45 in Trussell-Wilson/Sommer). Separately, the data show a
     readiness hump at 25-30 the single-logit entry cannot produce.
  d. Price pass-through to births is ~zero (10 percent rent move = half a
     percent of the child surplus; GE H0 +-7 percent moves tfr by 0.0006).
     This is the paper's central mechanism and it is currently absent. Any
     housing-gate design must raise the price share of the entry decision
     by roughly two orders of magnitude.
- Loss forensics: childlessness (200) + wealth dispersion (108) carry 80
  percent of the loss under every reweighting tried. The weights are not
  the binding problem.

## 2. The mission

In one week, deliver a calibration the author can stand behind, with a
decision package. Priorities in order:

1. Fix what is fixable without changing the model (search, seeds, fecundity
   tail constants — see 3.1-3.2).
2. Where a model change is unavoidable, implement it as a GATED VARIANT
   (default-off, bitwise-identical default path, suite green), refit, and
   quantify exactly what it buys on which rows and what it costs elsewhere.
3. End state: a ranked "what buys what" table across all variants tried,
   a recommended configuration, and full certified tables for it. The
   author decides adoption; you do not promote any variant to "the model."

Explicitly hoped-for outcome: no big model changes. Explicitly acceptable
outcome: a small set of changes, each with its price tag.

## 3. Work program

Treat this as the plan of record; deviate when evidence justifies it and
register every deviation in the ledger.

### 3.1 Fecundity tail (no model change; cheapest first)

The schedule is `1 - pi_a = w1 * exp(w2 (a-18))`, zero from 45, currently
(0.02, 0.134), hand-set (open item #37 in the ledger). Fit the tail to the
demographic evidence in `docs/model/fecundity_schedule_note_20260720.tex`
(Trussell-Wilson shares permanently infertile; Menken et al.; Leridon
four-year conception probabilities at 30/35/40 = 0.91/0.84/0.64). The job
is specifically to kill the 38+ excess while keeping the 30/35/40 levels
close to Leridon. Note the object: pi_a is a FOUR-YEAR conception
probability for a < 45 and zero after. A steeper tail may need the
functional form loosened (e.g., free terminal decay) — that is a parameter
change, not an architecture change; gate it anyway. Refit (E6a) and record
the deltas on all twelve rows.

### 3.2 Search hygiene on every refit

The E5 winner was dominated by its own restriction because all eight chains
sat in one basin. For every refit: seed from the current certified winner
AND from at least two dispersed seeds (start_mix or perturbed seeds);
compare basins before certifying. Strict twice-repeated tight winners only;
the collector is the only reportable source.

### 3.3 Income fixed-effect layer (the one likely-necessary model change)

a. First measure properly: STY-style minimum-distance decomposition of PSID
   household log earnings into fixed effect + AR(1) + transitory
   (extend `build_psid_income_between_within.R`; the raw-data conventions
   are in `audit_aggregate_wealth_earnings_ratio.R`). Deliver the fixed-
   effect variance with a defensible split from the persistent component.
b. Implement permanent income levels as a gated variant (E6b): 2-3
   fixed-effect states at entry, spacing/weights from (a), multiplying the
   existing z-chain. This is the standard lifecycle structure (Huggett
   1996; Storesletten-Telmer-Yaron 2004; Kaplan-Violante 2010) — cite it
   as such, do not present it as a novel type model. State-space cost:
   solve time roughly times the number of levels; budget cluster time
   accordingly.
c. Identification: the new dispersion parameter(s) are pinned by the
   measured fixed-effect variance (external), NOT freely estimated —
   the free-parameter count stays ten unless you register otherwise with
   named identifying moments.
d. Refit; report what it buys on childlessness (level AND gradient — the
   model gradient across permanent levels vs the CPS education gradient),
   wealth p90/p50, aggregate wealth, and what it costs on the housing rows.

### 3.4 Timing shape (only if 3.1 leaves the 25-30 hump missing)

Readiness arrival: a binary "settled" state with an age-dependent hazard
gating entry, shrinking the entry logit's role. Gated variant E6c. The two
timing rows identify the hazard's location and scale. This is an
architecture change — implement only if the fecundity tail fix plus income
levels leave share30+ badly missed, and say so in the decision package.

### 3.5 What NOT to do

- No housing-gate mechanism (chi_O(n), rooms requirements, Stone-Geary
  re-entry). The author explicitly reserved this decision. Your job is to
  make everything else right so that decision is clean.
- No changes to the signed target system or weights as the reportable
  baseline. Diagnostic reweighting runs are fine if labeled; the certified
  baseline stays the signed system.
- No touching phi = 0.80, tenure determinism (kappa_T = 0), theta_n = 0,
  or the beta_annual >= 0.94 external restriction. Never hard-target
  median moments. The word "parity" must not appear in anything written
  for the author.
- No edits to latex/ or the M strand. No force-push. No policy-experiment
  claims in any author-facing text (policy numbers remain diagnostics).

## 4. Operating rules

- Cluster: anything over ~30 minutes goes to Torch via the patterns in
  `code/cluster/submit_intergen_e5b.sh` (array + afterok collector).
  Partition cpu_short caps at 3:55; account torch_pr_570_general. The
  scratch copy at `/scratch/td2248/projects/Fertility_Spring26` is an
  rsync mirror, NOT a git clone — rsync exact files, never git-pull there.
  Smoke-test the exact loop before every array (the `--smoke` flag on
  `run_e1_chain.py`). Budgets, checkpoints, and best-so-far artifacts per
  `AGENTS.md`; a run without a heartbeat for 30 minutes is unhealthy.
- Do not submit cluster jobs with `sbatch --wrap` (it runs /bin/sh and the
  module preamble silently fails — this bit us on 2026-07-26); use real
  submit scripts copied from the working ones.
- Verification: every figure or diagnostic solve at a winner must
  reproduce the certified twelve moments to 1e-6 before its output is
  used (`run_e5b_autopsy.py` shows the pattern). The full test suite
  (`code/model/.venv/bin/python -m pytest code/model/intergen_eqscale_seq_optimized/tests -q`,
  123 tests) must pass with all gates off after every code change.
- Known traps: theta stores the FOUR-YEAR discount factor (annual =
  value^0.25) and the probe collector's `beta_annual` column actually
  holds the period value; `pi_a` is a four-year conception probability;
  the completed-fertility moment uses the CPS capped top-group mean
  3.602359422009 on both model and data sides; utility-side hooks must go
  through Bellman AND KFE (see `bellman-only-hooks-gotcha` in memory).
- Agent economy: delegate mechanical work (scripts, extractions, table
  builds) to cheap workers; keep your own context for economics,
  identification, and verification. Never accept a delegated calibration
  conclusion without checking the full target-fit table yourself.
- Reporting discipline: any calibration readout includes the full
  twelve-row target-fit table (target, model, gap, weight, contribution)
  and every parameter with its bound status. Loose mid-search losses are
  never reportable. Never compare losses across different target systems
  or weight schemes.

## 5. Record keeping and the final package

- Register every decision, launch, and verdict in
  `docs/model/e5_target_review_20260724.md` (append dated sections, same
  style as the existing entries). Commit and push after each coherent unit
  of work; the daily 23:40 backup job also pushes.
- Keep a daily progress note in `memory/daily/YYYY-MM-DD.md`: what ran,
  what it showed, what is queued, next decision.
- Final deliverable, in `docs/model/e6_decision_package_<date>.md`:
  1. The ranked what-buys-what table: for each variant (E6a fecundity
     tail, E6b income levels, E6c readiness, and combinations tried),
     the certified loss, the twelve-row fit deltas vs E5b, parameters
     with bound flags, solve-cost change, and the one-line economics of
     what the variant is doing.
  2. A recommended configuration with the complete certified tables.
  3. The open questions that are the author's to decide (housing gate;
     any target-system change; adoption of any variant as the spec).
  4. An honest list of what was attempted and failed.
- The standing baseline for all comparisons is E5b. If nothing beats it
  defensibly, say so — that is an acceptable outcome.
