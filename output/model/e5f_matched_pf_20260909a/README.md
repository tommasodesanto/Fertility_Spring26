# Matched perfect-foresight implementation

Author instruction, September 9: approximate initial steady state; households
learn the entire transition immediately at its onset; calibrate observations
along the transition; pursue sequential and simultaneous fertility nests in
parallel under the presentation time constraint.

## Delivery order

1. Verify dated household choices and forward accounting in both arms at common
   parameters. Existing estimates initialize these tests; they are not PF fits.
2. Integrate historical PF evaluation with the retained historical age margins
   and target measurement. Keep 2023 as an observation on the path.
3. Append a consistent post-2023 demographic and equilibrium continuation;
   solve the baseline before estimating the historical PF objective.
4. Recalibrate both arms under the same resolved empirical target contract;
   prioritize a complete sequential delivery if only one arm can be finished.
5. Compare the same rebated property-tax experiment, using coherent person/head
   accounting and separately verified baseline/reform transitions in each arm.

The scientific comparison concerns choice information and correlation, not only
shock timing. Retain housing dispersion 0.005 and housing-supply elasticity 0.63
for initial tests. Re-estimation of fertility dispersions remains part of the
eventual calibration. Do not choose the baseline for its preferred policy sign.

## Source and isolation

Worktree: `tmp/e5f_matched_pf`, branch `codex/matched-perfect-foresight`, starting
at clean nested-source commit `d122eb52`. The working project and frozen prior
results are preserved. Source changes are not promoted to production. Implementation commits
`6d04d342`, `a77abebd`, `9d4aa48b`, `eb3152eb`, and `649da484` are backed up
on `origin/codex/matched-perfect-foresight`. Large state arrays remain in the
hash-pinned cluster run folders; this local folder retains the lightweight
contracts, launch scripts, stage receipts and independent diagnosis.

## Test scope and limits

Runtime planning: the measured household pilot takes approximately three to
three-and-a-half minutes per arm; the sequential twelve-call joined pilot takes
just over three minutes. These timings are for supplied price paths and cached
terminal values. A market-clearing price iteration, a terminal root, and an
estimation search require repeated paths; do not quote these smoke times as
full-equilibrium or calibration times. Each long follow-up needs its own measured
solve count, checkpoint plan and stopping budget.

The first paired test uses the same selected no-policy checkpoint in both arms,
the original 120-node wealth grid, one fixed-price stationary household solve,
one constant-continuation reproduction and a two-date backward/forward replay:
six Bellman calls per arm, twelve total. Two one-core jobs, each capped at
15 minutes (14-minute internal limit), measure this updated computation.
Historical full-grid timing is not a forecast for this pilot.

Checks cover choice-array reproduction, conception and distribution accounting,
actual dated rental costs, preserved policy objects, and checkpoint output.
Fixed stationary entry is used only to test household-operator invariance.
These are household primitives: no claim of market equilibrium, person/head
closure, fiscal balance, terminal equilibrium, historical fit, or full equilibrium
precision. No full calibration may be launched on the strength of this pilot
alone. Each stage writes timing and the latest completed-stage receipt; failures
stop the job. Saved arrays permit independent review. No figures are requested.

## Historical integration contract

The existing historical preference path moves linearly from the old intercept
to its 2023 value over four model periods. Immediate full information means
knowing that path from 2007. The initial conditional household states come from
an old stationary solution, then are reweighted to observed 2007 ages. The
2011/2015/2019/2023 age margins remain empirically imposed in the first diagnostic.
That is conditional historical calibration, not an endogenous population fit.

Retain the existing cohort timing ledger and 2019-to-2023 birth-housing branch;
the latter must not receive the aggregate age reweighting. Retain the full
twelve-target objective for diagnostics, without claiming the unresolved
ownership-family groups and date window have been aligned.

The post-2023 continuation must be specified and solved before a historical PF
objective is certified. A supplied terminal value in a short plumbing test is
an explicit boundary condition, not evidence of terminal convergence. Historical
2007 person/head inputs cannot be manufactured from the existing 2023 builder.

## Current execution

- Primitive pilot **17276750** (snapshot `...matched_pf_20260909e`): sequential
  passed six full-grid Bellman calls in172.17seconds. Values and all choice
  arrays reproduce exactly. Maximum mass residual4.44e-16; stationary current
  distribution L1 gap7.32e-15, next-state invariant gap6.98e-10; no feasibility
  projection or budget excess mass. Markets are prescribed, not solved.
- Nested in the same pilot: stationary solve succeeds and constant dated values
  and native joint arrays reproduce; effective tenure probabilities differ by
  up to4.963e-5 on1374cells. Independent diagnostic17276871
  traced the discrepancy to roundoff in the conditioning population: probabilities
  on the same population reproduce exactly; product-mass error8.50e-17, current
  distribution L1 error2.24e-14. Corrected nested pilot17276953 **passed
  all six calls in206.76seconds** on fresh snapshot `...matched_pf_20260909f`
  with unchanged tolerances.
- Joined historical/person diagnostic **17276868**: **passed in188.09seconds**,
  six dates2007--2027, twelve
  Bellman calls, original120nodes, one core/32GB/15minutes. The passed sequential
  primitive packet is the seed. Constant prices/preference and zero transfers
  are supplied; observed historical age bridges and the frozen2023 person/head
  mapping exercise the actual composition. All12 accounting/reproduction gates
  passed, zero feasibility projection, all dated budgets passed. Maximum market
  residual15.83% confirms this prescribed-price path is not an equilibrium. Synthetic current-preference seed,
  not old2.1 normalization; supplied terminal value, not endpoint convergence.
  Inputs and source hashes: `history_contract.json`; launch: `submit_history_smoke.sh`;
  receipt: `history_submission.json`. Remote snapshot:
  `/scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a`.

- Matched nested joined-history test **17277058** **passed in246.68seconds**
  from snapshot F, using its passed primitive seed and identical diagnostic
  design. All12 gates and dated budgets passed; exact value reproduction,
  zero projection mass. Prescribed-price market residual15.95% is not gated. Contract and
  launch are `history_nested_contract.json` and `submit_history_nested_smoke.sh`.

## Numerical diagnosis and verification

Pilot17276031 stopped before solving because a new guard checked legacy supply
fields rather than the authoritative saved supply rule. The corrected guard
retains elasticity0.63. Reduced12-node pilot17276095 then failed stationary
feasibility in both arms; neither reached a PF date. No gate was weakened.

Original-grid pilot17276230 passed stationary feasibility but failed strict
constant-value reproduction. Diagnostic17276528 established that rent arithmetic
changed constant rent by1.39e-17 and values on six unoccupied states by3.73e-9.
Writing the identical user-cost relation as
`rent = user_cost * current_price + (current_price - next_price)` restores
bitwise stationary rents/values in both arms. Native nested arrays also reproduce;
the remaining tenure check concerns distribution conditioning and is separate.
The new formula checks consistency with the underlying carrying-cost parameters.

Thirty-one affected regression tests passed on Torch, job17276658, including
existing PF numerical checks. Eight joined-composition/CLI tests passed in the
new job before numerical startup. Tests cover tail-to-history anticipation,
2023 counted once, historical age bridges, dated choice ownership, rents,
checkpoint/source guards and household/person accounting. No completed historical
PF objective or matched production policy effect is claimed.

## Remaining delivery gates

1. Household primitive tests: passed in both arms.
2. Joined historical/person calculation: passed in both arms.
3. Restore the normalized old initial state, extend and converge the terminal
   tail, and solve market/fiscal paths under the explicit baseline contract.
4. Attach unchanged historical moment measurement, reconcile outstanding
   ownership-group/date definitions, and reproduce a complete target ledger.
5. Only then launch matched estimation and the rebated property-tax comparison.

The earlier verified sequential PF pair remains a separate fallback artifact;
current temporary-expectations results are not relabeled as perfect foresight.

## Final verification state

All numerical jobs listed above have finished. Both arms pass the original-grid
household primitives and the six-date conditional historical/person composition.
Downloaded rows, observations and initialization records match the SHA-256
hashes in their successful receipts. The isolated branch is clean and pushed.
No estimation, converged baseline, production policy run, or figure refresh was
launched in this implementation pass. The remaining delivery gates above are
substantive work, not a claim that the desired research result is already done.
