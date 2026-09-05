# Simplified OLG: development plan and overnight record

Authorized by Tommaso on 2026-09-04: “yeah, let's do this.. so you can work on
this overnight while i sleep hopefulyl.” This begins development after the
discussion pass. It does not convert provisional economic choices into final
author decisions. The issue tree remains in `ACTIVE_DECISION_LEDGER.md`.

## Scope, deadline, and deliverables

Work through the night of September 4–5, with a morning checkpoint by
**2026-09-05 08:00 America/New_York (12:00 UTC)**. Stop earlier if the packet
is complete. Scheduled continuation resumes this same task once per hour;
each turn advances one bounded unfinished phase for at most 45 minutes.
At the deadline, finish the checkpoint, report remaining work honestly, and
pause the schedule. If resumed after the deadline, wrap up rather than start
another development phase. Do not repeat completed checks without a new reason.

Produce one consolidated, readable LaTeX/PDF amendment proposal, supporting
reproducible checks, and a short morning readout in this file. For each claim,
distinguish **proved under stated assumptions**, **numerically illustrated**,
**candidate**, **not established**, and **requires author decision**. A failed
conjecture is a valid research outcome; do not replace a missing theorem with
assertions. No particular theorem is guaranteed by the deadline.

Proposed manuscript: `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`.
Readable output: `output/pdf/simplified_olg_amendment_proposal.pdf`.
Supporting evidence: `output/model/simplified_olg_amendments/`.
Use existing theory drivers where suitable; one focused verification driver
may own the new analytical checks. Preserve the independent review as historical
evidence. No files in `latex/JMP_DS_draft/` may be written, including build files.

## Accepted directions and retained decisions

- W0/W1: establish an allocation gain at fixed fertility and cohort masses;
  borrowing limits may be relaxed. Compare housing willingness to pay in
  consumption units, and respect real resources, creditors, and estates.
- W3/W4 are **provisional**: household-specific compensation and a nonnegative
  real reallocation cost. Show the actual transactions and a simpler-transfer
  comparison. Return these concrete choices to the author before presenting
  the welfare benchmark as final. Keep marginal costs and fixed moving costs
  distinct. Seller compensation is a transfer, not resource destruction.
- I0/Q0/W2: capital-gains taxation is deferred. Develop the core at zero gains
  tax. Property tax is separate and is not silently removed. A zero-property-tax
  specialization can be informative if explicitly labeled as such.
- H0/H1: one young completed-fertility choice and one old period. Seek an
  explicit sufficient condition, preferably with no prices or multipliers;
  distinguish household partial effects from equilibrium policy effects.
- D0/D1/D2: compare policy and baseline transitions from the same initial
  state, then total population at their endpoints. Seek the broadest valid
  transition result. Do not assume small reforms in advance or equate finite
  terminal closure with an infinite-horizon theorem. Fixed housing stock is
  the benchmark; assess a simple stationary supply extension.
- U0: repair inconsistent units by an explicit change of variables; keep
  adult households, children, and resident persons distinct. Do not change
  the quantitative model's demographic contract.
- P0/P1: presentation is secondary. Preserve all 18 issue IDs and 14 review
  repairs, including parked alternatives and the conditions for revisiting them.

## Work sequence and acceptance checks

| Phase | Economic task | Required evidence | Status |
|---|---|---|---|
| 1 | Construct the compensated old-to-young reallocation at fixed fertility. | Dated goods, title, debt, sale, and estate ledger; admissibility/slackness conditions; costs and welfare coverage; simpler-transfer alternative. | CHECKED local proposition; compensation/cost remain provisional |
| 2 | Derive and audit the fertility threshold. | Exact sign calculation; a condition without multipliers; examine a genuinely primitive condition; negative-response example; capped-old-renter and common-cap correction. | CHECKED local threshold and primitive stationary specialization; aggregate policy sign open |
| 3 | Establish the transition result and its limits. | Separate existence, convergence, uniqueness, and policy signs; attempt broader continuation under explicit bounds; identify failures or missing hypotheses; check the zero-shock system when using a local theorem. | NEXT: broad transition argument and tax-free diagnostic |
| 4 | Reconcile population units and stationary supply. | Exact rescaling, transition population identity, common initial state, terminal household/person definitions; supplier-income/financing qualifications. | Algebra drafted; literal-person convention still requires author review |
| 5 | Consolidate amendments and deliver the PDF. | Reconcile proposal against section, appendix, builder, and claim ledger; relevant small checks; compile twice, inspect rendered pages; preserve frozen evidence and unrelated edits. | Seven-page proposal compiled and inspected; full source reconciliation pending |

Prioritize phases 1–2, then 3–4, then presentation. If a broad transition
theorem cannot be justified, deliver the strongest valid result and the precise
remaining obstacle. Do not silently weaken the requested objective.

## Execution and file ownership

The lead owns economic judgment, derivations, synthesis, and final verification.
Use at most one independent worker at a time for a bounded adversarial check
through the repository wrapper. Default reviewer limit: 20 minutes. Inspect its
actual output; do not accept its conclusion without checking the argument.
No duplicate worker if an earlier one is still active. No silent retries after
timeout. Record any worker PID/session, output, scope, and completion below.

Only small illustrative theory computations are in scope. No calibration,
production policy runs, target changes, or quantitative sweeps are authorized
by this theory request. Before any numerical experiment, record solve count,
estimated time, and stopping rule; smoke-test any loop before expansion.
Keep individual local experiments under 10 minutes and total local numerical
experiments under 30 minutes. Anything requiring a longer run needs a new
explicit plan under the project's cluster procedure, not an unattended launch.

The checkout contains substantial unrelated work, including another overnight
quantitative task. Never stage all changes. Prefer the new suggestion/check
files and targeted discussion-map edits; inspect diffs before changing existing
theory sources. Source reconciliation can document a proposed replacement
without promoting provisional assumptions into the active manuscript. Keep
`CALIBRATION_STATUS.md` and quantitative code unchanged. Commit and push only
coherent, verified changes owned by this task. Do not force push.

## Checkpoints and morning readout

### Initial checkpoint — 2026-09-04, approximately 22:10 EDT

- Startup context and decision tree loaded; active theory section inspected.
- Work plan recorded. Hourly same-task continuation is active as `simplified-olg-overnight-theory`;
  it must pause at completion or the stated morning deadline.
- New fertility lead: on the strictly binding down-payment branch, the
  stationary housing expenditure at the cap can be independent of the house
  price. This may yield a primitive sufficient condition when lifetime
  transfers are exogenous (in particular in a zero-property-tax specialization).
  It has not yet passed an independent check and is not a transition theorem.
- Next: write and verify the candidate, and construct the dated compensation
  ledger while one bounded reviewer checks the fertility argument.

Append substantive checkpoints here; update the phase table rather than
creating a separate status file each hour. The morning readout should lead with
what is established, then what failed or remains conditional, then the small
set of concrete choices for Tommaso. Include source/PDF/check links.

### First analytical verification budget

- Worker: `reviewer_strong`, read-only, wrapper session 39561, 20-minute limit;
  output `output/model/simplified_olg_amendments/independent_review.md`.
  Scope: the two new welfare/fertility arguments only; no transition audit.
- Local verification: 6,075 candidate reduced-household configurations, with
  only binding cases solved and three scalar roots per retained case; smoke
  the same loop on one preference triple first. Add 108 direct welfare-ledger
  evaluations and one common-cap finite-difference check. Expected under one
  minute; stop at ten minutes or any failed assertion. No equilibrium solve.
- Source proposal created. It gives an exact seller compensation function,
  an alternative small uniform premium, explicit debt repayment, the primitive
  stationary fertility restriction, and exact population accounting. Claims
  remain subject to checks and the independent review below.

### Verified checkpoint — 2026-09-04 22:13 EDT

- **Phase 1:** explicit exact compensating transfer, bond reserve for property
  tax, seller's estate reserve, buyer repayment out of old consumption, and
  unchanged government receipts and next-date housing supply. A small uniform
  seller premium is a simpler local alternative to exact household-specific
  compensation. These are checked conditional propositions, not final W3/W4
  instrument choices. Fixed moving costs require the separate finite-trade test.
- **Phase 2:** a sufficient threshold with no multiplier, plus a stationary
  exogenous-resource specialization with no price either. The primitive
  inequality is `(1-q)b/[(1-phi)(y+b)] >= alpha/(rho+alpha)` when property tax
  and its rebate are zero and old choices are interior. It is conditional on
  strict down-payment binding and maintained branches; it is not a GE sign
  theorem. The report includes a complete conditional-owner example satisfying
  saving, retention, and estate inequalities.
- **Verification:** 4,280 binding-cap configurations;
  796 satisfy the primitive sufficient condition
  and all have positive local responses. 548 other
  configurations have negative responses. Maximum derivative/finite-difference
  discrepancy is 6.69e-08. The 108 constructed
  welfare cases satisfy the budget/repayment checks; exact compensation keeps
  the seller indifferent and a uniform premium makes both parties better off.
  The common-cap derivative correction also agrees with finite differences.
  These are constructed household checks, not calibrated or GE results.
- **Reviewer:** session 39561 completed successfully, within its 20-minute
  limit. Both arguments passed the local, fixed-price/branch scope. Its two
  accounting clarity corrections (explicit tax-reserve asset and explicit
  seller title receipt) were adopted and checked by the lead. Review output:
  `output/model/simplified_olg_amendments/independent_review.md`.
- **PDF:** seven pages, compiled twice; no undefined references or overfull
  boxes. All pages visually inspected, with the transaction table inspected
  at full size. The analytical diagnostic figure was also inspected. An early
  summary-output attempt encountered a NumPy integer JSON serialization error;
  it was corrected, and the complete verification driver now passes in under
  two seconds. No model result was accepted from the failed output attempt.
- **Phase 4:** exact transition-product and terminal-scale identities, complete
  child-unit rescaling, and stationary supply accounting drafted. A resident-
  person counting convention remains an author choice. No quantitative units
  or targets were changed.
- **Next wakeup:** phase 3. Read the existing transition residual and local
  descriptor code; assess a broad continuation result, then run only the
  smallest justified tax-free diagnostic with explicit overrides and outputs
  under this workflow. Existing driver defaults use positive gains tax and
  overwrite older figures: do not invoke its main unchanged. Preserve those
  historical outputs. No worker is still active for this theory task.
- **Still unfinished:** broad transition existence/convergence/uniqueness;
  aggregate policy fertility and terminal-population signs; reconciliation of
  all 14 original repairs against the final source proposal; final morning
  summary. Do not pause the automation merely because the first PDF exists.
