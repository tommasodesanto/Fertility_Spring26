# Matched perfect-foresight implementation

Latest overnight plan: [OVERNIGHT_PLAN.md](OVERNIGHT_PLAN.md). The100-date anchor has completed; the new price panel is queued. Read this before the older dated checkpoints below.

Author instruction, September 9: approximate initial steady state; households
learn the entire transition immediately at its onset; calibrate observations
along the transition; pursue sequential and simultaneous fertility nests in
parallel under the presentation time constraint.

## Read first: what is established and what remains

Both the sequential model and the simultaneous fertility-nest model pass the
full historical perfect-foresight calculation at supplied prices. The old
steady-state normalization, stationary terminal equilibrium, twelve target
measurements, and household/person accounting have also passed. This establishes
that the matched calculation works; it does not yet establish a calibrated
historical equilibrium. All eleven inherited parameter values and all target
weights are unchanged. Households know the transition from 2007. Historical
age composition is imposed from the observed age profiles through 2023; the
person-demographic model advances the population afterward. Thus this is a
conditional historical fit, not yet a fit to the historical population path.

The 12-date paths in both arms and the 28-date sequential path now clear
markets and reproduce exactly. The longer sequential objective is94.5223,
versus129.2431 on the short path, at identical parameters. The change is
material; horizon stability remains unverified. See
[HORIZON_COMPARISON.md](HORIZON_COMPARISON.md) for all targets and restrictions.
A single100-date prescribed-price diagnostic17300115 is submitted; this does
not yet clear its markets or re-estimate parameters.

The main remaining empirical decision concerns two ACS family-group definitions
and the comparison of four pooled ACS rows with a 2023 model cross-section.
The new date comparison preserves the original samples and weights. It finds
material changes in the two family gaps, smaller changes in aggregate ownership
and rooms. It does not change the calibration contract or reopen the
Sun–Abraham childbirth regression.

The historical baseline has a 1% annual property tax and no rebate. Holding the
2023 fertility-preference intercept constant afterward is an explicit diagnostic
continuation. Neither this continuation nor a new matched policy comparison
has been promoted as a final benchmark.

## Meeting work, September 9

Normalized old states passed in both arms (jobs17277586/17277587). Terminal
price/person roots also passed (sequential17277907, nested17277683): eight fresh
endpoint evaluations per arm, including exact final replay. Prices0.46420547 and
0.46391036, signed housing residuals4.619e-6 and−1.991e-5, all population gates
passed. These retain annual1%tax, zero transfer and the inherited2007 supply
normalization. The terminal preference holds normalized2023psi constant; the
post2100 demographic primitives are frozen. These are explicit diagnostic
assumptions, not an estimated continuation or a complete historical equilibrium.

See `meeting_receipts/READOUT_VERIFIED.md` and `verification.json` for small
receipts and independent checks. The next stage measures all12 historical
targets along a normalized PF path. The first horizon has12dates2007–2051,
24 Bellman solves, with terminal distribution2055. Expected6–8minutes per
prescribed-price arm from prior measured timings; watchdog28minutes/onecore32GB.
This short horizon is an exact-loop smoke and finite-boundary diagnostic,
not a horizon-converged historical objective. Conditional price probes will
follow only after it passes; each writes complete fit/parameter tables,
all dated markets and terminal-distance checks. No target or parameter search
has changed. Initial startup17278252/17278253 stopped before model work on a
stale test-message assertion; source and failed logs remain in snapshot
`...matched_pf_paths_20260909a`. The corrected retry also records signed excess
demand for the forthcoming price Jacobian.

The corrected paired paths17278316/17278317 passed in420.68/631.59seconds,
including all12target rows and every inherited coordinate/bound. Receipts:
`meeting_receipts/path_anchor_02/{sequential,nested}/`. Table hashes and summed
loss contributions were verified after collection. These prescribed-price
losses121.223/139.517 are **not equilibrium calibration comparisons**: maximum
market gaps are71.08%/70.93%, and the2055 population is far from the endpoint.

Released after these exact-loop passes:

| Arm | 12-coordinate price panel | 28-date horizon check | Dependent collector |
|---|---:|---:|---:|
| Sequential |17278556|17278557|17278711|
| Nested |17278629|17278633|17278712|

Each panel changes exactly one log asset price by0.01 with fixed parameters,
initial state, supply, demographic inputs and terminal endpoint. Per arm:
12×24=288 Bellman solves; expected7.0/10.5minutes per case, with independent
cases parallelized. Each horizon check has56 Bellman solves, expected16.4/24.6
minutes. Every job has a28-minute watchdog and writes progress during the
forward pass plus a15-second health heartbeat. Queues can add wall time.
The original sequential jobs reserved32GB; measured peak was3.1GB. A pending
memory-update attempt raced with scheduling and was rejected once the jobs
started; no running job was interrupted. Newly submitted nested probes reserve
16GB (measured peak about7.4GB), and the longer nested path20GB.

The collector requires all coordinates and verifies signed excess demand from
hash-pinned market CSVs; it rejects mixed source/target/input fingerprints. It
cannot launch calibration. The subsequent bounded Broyden root uses a reserved
fresh final replay and maintains separate market, terminal-distance and horizon
extension classifications. Source commit595690e9 is backed up on the isolated
branch. All83 focused tests pass on Torch job17278618;49 root/collector tests
also passed locally with `/usr/bin/python3` (plain `python` lacks NumPy here).

Sequential panel17278556 and collector17278711 completed successfully. Its
log-price market Jacobian has condition number3.681 and negative own-price
entries at all12dates; the reviewed first Newton step is clipped at0.10 in
absolute log-price change. Bounded sequential historical root17279004 uses
at most six full paths (144 Bellman solves including its reserved final replay),
expected42minutes at measured anchor timing, with a60-minute watchdog.
It writes every trial's complete fit/parameters/markets separately and saves
latest and best residuals. Contract: `historical_root_sequential_contract.json`;
validated matrix: `meeting_receipts/jacobian_sequential.json`.

Nested panel17278629 and collector17278712 also passed: condition number3.634,
with a closely similar local price response. Nested root17279843 permits five
full paths (120 Bellman solves), estimated53minutes with a60-minute watchdog;
it was submitted at18:51UTC and initially queued for memory. Its saved matrix
is `meeting_receipts/jacobian_nested.json`.

The28-date sequential anchor17278557 passed in898.78seconds (15minutes),
with measured peak memory3.41GiB. It remains off-equilibrium (maximum gap39.63%)
and its2119person total is40.67%above the stationary endpoint. Its conditional
loss87.555 has complete fit/parameter tables in
`meeting_receipts/path_long_anchor_01/sequential/`. The longer numerical seed
also changes provisional post2023prices, so this is **not a controlled horizon
comparison** and the loss change is not attributed solely to horizon length.
Longer sequential price panel17279762 has28cases×56=1568Bellman solves, initially at
most16concurrent one-core/8GB jobs, expected15minutes each (two waves plus
queue overhead), and28-minute case watchdogs. Collector17279797 waits for
all28cases. No longer-horizon market root or parameter calibration was launched
yet. The short roots and longer panel use distinct saved output paths.

The four ACS pooled-versus2023 diagnostics are now complete in
`meeting_receipts/acs_date_diagnostic/`. A single110MBcache read in7.4seconds
reproduced the authoritative pooled room counts, weights and points before
extracting2023means. Original family/sample definitions and all calibration
weights are unchanged. This resolves missing descriptive evidence, not the
author's choice about model/empirical groups or calendar alignment.

## Conditional next calibration design—not launched

Once a baseline price path is accepted, write equilibrium residuals as
F(x,theta)=0 with x=logprices. The price panel supplies F_x; the same saved
moment tables can supply M_x. At fixed **absolute solved prices**, eleven
parameter probes would rebuild the old2.1normalization, old cohort prehistory,
2007supply anchor, mechanically implied demographic alignment and consistent
terminal endpoint. Then solve F_x Z=F_theta and use G=M_theta−M_x Z for the
local equilibrium moment derivative. This saves complete price roots for
actual candidate tests rather than every derivative column. A rent projection
must not silently move prices inside a parameter derivative. The current
price panel is off-equilibrium and is a root preconditioner; it requires a
locality check or refresh before use as an equilibrium derivative.

A future round must inspect rank/conditioning of the weighted12×11moment
matrix, use bounded steps in the existing parameter coordinates and verify
actual complete equilibrium candidates before accepting any improvement.
Estimated post-baseline runtime at current short horizon: up to19minutes for
parallel normalized parameter probes and about59minutes per candidate allowing
six full price mappings, with three step lengths parallel if resources permit.
This is a recoverable proposal, not a launched calibration, identified global
result, production target approval, or horizon-convergence certificate.

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

## Earlier implementation-pass verification state

All numerical jobs listed above have finished. Both arms pass the original-grid
household primitives and the six-date conditional historical/person composition.
Downloaded rows, observations and initialization records match the SHA-256
hashes in their successful receipts. The isolated branch is clean and pushed.
No estimation, converged baseline, production policy run, or figure refresh was
launched in this implementation pass. The remaining delivery gates above are
substantive work, not a claim that the desired research result is already done.

## Meeting checkpoint, 19:15 UTC

Sequential short root17279004 has five valid evaluations: maximum market gap
0.710809,0.291869,0.071806,0.016467,0.003107. The reserved sixth path is an
uncached replay. This is still above the unchanged0.0002 acceptance tolerance;
a bounded restart can retain the verified prices and approximate Jacobian.
Nested root17279843 began at19:00UTC; its first mapping exactly reproduces the
original supplied-price residual. Both are finite-horizon diagnostics.

The updated root driver supports a verified completed-root restart and an
optional pre-choice2023 stock checkpoint. It changes no economic kernel or
numerical tolerance. Source883b49f9 is pushed, and all95 focused checks passed
on Torch17280998. The checkpoint distinguishes the inherited stock from2023
choices and remains uncertified until its parent baseline passes. New work
uses the immutable source snapshot `...matched_pf_path_roots_20260909b`;
running roots retain snapshot A.

Both longer anchors passed. The nested28-date anchor took1270.84seconds; its
complete twelve-target and eleven-parameter tables are saved beside the
sequential long anchor. Long sequential panel17279762 has completed its first
sixteen cases; the remaining twelve are running. The concurrency limit was
raised to28 after measured memory use supported it; no numerical case was
interrupted. All per-case budgets, thread limits and scientific inputs remain
unchanged.

Endpoint interpretation: the28-date final population gap does not by itself
measure historical moment bias. The active household has17four-year age cells
and direct utility from its own bequest, without dynastic continuation. At
fixed dated prices/primitives, a2023age18household's last decision is2087 and
the last rental price also involves the2091asset price. A2119terminal value
therefore cannot directly enter its value. Equilibrium prices can nevertheless
transmit population/tail changes backward through overlapping lifetimes. A
controlled horizon extension must preserve all shared prices and primitives
first, then re-clear the extended path and compare all historical moments.
Population closeness and historical horizon stability remain separate checks;
no terminal gate is waived. The first minimal extension after a solved28-date
path can append one four-year period, with broader extensions if sensitivity
remains. Source review: optimized solver age recursion and direct bequest
utility; dated PF rental identity.

## Verified continuation launches, 19:23 UTC

Short sequential root17279004 finished its six-call budget at maximum signed
market gap0.0031072113, above tolerance0.0002. Evaluations5and6 have exactly
identical full target, parameter, market and measurement files. Complete tables
and verification are in `meeting_receipts/historical_root_01/sequential/`.
The conditional loss129.480313 is not a converged calibration result. The
last state arrays remain remote; their hash was not recomputed locally.

Continuation17281783 uses the verified best prices and approximate Jacobian
from that completed run. It has four fresh paths including another final
replay, expected26minutes,30-minute watchdog,35-minute allocation,onecore8GB.
The original scientific inputs,12dates,tolerances and target system are
unchanged. Its three restart receipts are hash-pinned; source D adds only
checkpoint saving, restart validation and budget controls. Contract:
`historical_root_restart_sequential_contract.json`. Output: source D
`output/historical_root_restart_01/sequential/`.

Long panel17279762 and collector17279797 completed all28cases successfully.
The28×28log-price Jacobian has condition number2.256; all own-price derivatives
are negative (−2.076to−1.742). Long sequential root17281784 now solves28dates
using that verified matrix: six paths including final replay,336Bellman
solves,expected90minutes,110-minute watchdog,two-hour allocation,onecore8GB.
Contract: `historical_root_long_sequential_contract.json`; output in source D
`output/historical_root_long_01/sequential/`. This work can continue beyond
the meeting. Both new jobs preserve the source snapshot and write individual
case tables,latest/best summaries and15-second heartbeats. If a mapping fails,
the run stops; no target, parameter or tolerance is relaxed.

Nested short root17279843 remains in its original snapshot C; its second
valid mapping reduced the maximum gap from0.709334to0.293138. It is separate
from these sequential continuations. No parameter calibration, new production
policy, or benchmark promotion has been launched.

The first D restart mapping reproduces the C best market residual exactly.
All twelve target rows, eleven parameter rows and measurement records are
byte-identical across versions. Every dated economic value also agrees; the
only four path-file differences are source-location labels in historical age
audit records, with identical input hashes. The optional inherited2023 stock
checkpoint was written and reloaded successfully (13.5MB, kept on Torch).
Receipt: `meeting_receipts/restart_driver_reproduction.json`. This checks the
updated driver in the real full-grid loop, not just unit tests.

## Latest handoff, 19:52 UTC

The short sequential continuation17281783 finished four calls with exact full
fit/parameter/market/measurement replay. Its gap0.00031196 remains above the
0.0002gate. Further continuation17282612 has four fresh paths, a30-minute
watchdog and35-minute allocation. Nested root17279843 finished five calls
with the same complete replay checks; its gap is0.01576076. Continuation
17282613 permits six fresh paths, expected about58minutes,70-minute watchdog
and75-minute allocation. Both run from hash-pinned completed receipts and
unchanged source D. The longer sequential root17281784 has completed two
valid mappings, reducing its maximum gap from0.396347to0.195309.

The corrected fixed-policy tail-sizing diagnostic17282572 passed two actual
stationary updates (maximum distribution/person difference below1e-8), two
anchor updates, and an exactly matching start of the full forecast. It required
zero Bellman solves. Starting from the supplied-price2119anchor, the existing
terminal level/distribution thresholds first pass in2387 after67additional
four-year updates (95dates from2007), taking109.55seconds. This is conditional
on fixed terminal policy/prices, not a prediction of the horizon needed by a
re-cleared equilibrium. It confirms that terminal population adjustment can
be much slower than household lifetimes; historical price/moment stability
still needs its separate controlled check. Receipts and verification:
`meeting_receipts/tail_sizing_02/`.

The first tail-sizing job17282386 is invalid and must not be used: its extra
diagnostic wrapper omitted activation of the sequential calendar routines.
The corrected wrapper matches the maintained forward block line by line and
adds the stationary check. Main market-root drivers already performed the
required activation and are unaffected. The first invalid script/output stays
on Torch under tail_sizing_01; the local exclusion receipt records the cause.
Before accepting a reused smoke launcher, verify both smoke receipts contain
exactly two completed updates; bounded timeout exit alone is insufficient.

No calibrated PF benchmark or matched production policy result is claimed.
Current decisions remain empirical family-group/calendar alignment and the
future preference continuation; no target or weight has changed.

September10,04:09UTC: source96a41873 passes58cluster tests. Full100-date market root17304265 is queued after the price collector. See OVERNIGHT_PLAN.md for its bounded budget and source reconciliation.

## Two-page algorithm and calibration primer

The author requested a short teaching note on the implemented solution and calibration architecture. Read [the two-page PDF](../../pdf/model_solution_and_calibration_sketch.pdf). Editable LaTeX: [model_solution_and_calibration_sketch.tex](../../../docs/model/model_solution_and_calibration_sketch.tex). Revised at the author’s request: page 1 gives a literal stationary price-guess/household-solve/cohort-aggregation/Brent-update algorithm, including conditional consumption and housing, exhaustive saving optimization, and sequential fertility/outcome/housing timing. Page 2 gives the complete-price-path/backward-values/forward-population/Broyden-update algorithm and the intended outer parameter-search loop, alongside observation dates and the 12 moment groups. It explicitly distinguishes completed equilibrium work from pending horizon certification and re-estimation. No new numerical result is reported.

Delegated draft and independent lead technical/visual review are recorded in solution_primer_qa.json. Exactly two pages, compiled twice, no overfull boxes; both pages inspected. Rebuild with pdflatex from the source, directing all output to tmp/pdfs/model_solution_and_calibration_sketch, then copy the checked PDF to output/pdf. No slide or manuscript file changed.


September10,04:59UTC: both100-date smoke probes passed and were independently collected/validated (horizon100_smoke_verification.json). All98remaining price probes now RUNNING after the same array concurrency limit increased32to98; collector17303126 and root17304265 remain dependent. Read the00:59EDT scheduling update in OVERNIGHT_PLAN.md for measured runtime, resource budget and limits. No new equilibrium or calibration.


September10,05:42UTC: original case20hit its60-minute watchdog at date90/100. Single reviewed recovery17309087uses identical source/contract/gates oncs693and a separate output, preserving the failure. Existing collector must be superseded by explicit validated-recovery selection before existing root17304265can proceed; see OVERNIGHT_PLAN.md.35completion receipts and64running original cases at this check.
