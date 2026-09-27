# September 27 daytime continuation and diagnostics

Latest search status: Torch job18645479 was cancelled unstarted at12:23:35 EDT after the12:20:29 cutoff. It never received CPUs:0 of24 requested workers,0 new cases. No late replacement was launched; see `search/completion.json`. Local numerical diagnostics and native replays/smoke are complete; full equilibrium transition remains outstanding.

Author authorization: continue the calibration search; do policy, stationary
population, Jacobian and three-household-life diagnostics locally; compute a
same-parameter benchmark without the borrowing constraint. This supersedes the
completed overnight schedule's launch stop; the existing heartbeat now supervises this daytime work. The reference is overnight continuation de_0093, with the exact source,
parameters and targets recorded in ../supervised_calibration_20260927/primary_final/.

Separate work and ownership:

- `search/`: Torch continuation, at most 24 workers and 192 search cases,
  two exact-loop smokes and two final repeats. Absolute 3.5-hour budget leaves
  90 minutes for verification/export. Economic model, targets, weights, bounds
  and grids unchanged; new starting point/seed and deadline are search changes.
- `jacobian/`: local central differences for all nine searched coordinates,
  18 objectives and at most two workers, 90-minute total/15-minute case caps.
  Each objective retains the endogenous child-benefit normalization. Weighted
  local sensitivity is not proof of global identification. Step definitions
  and fresh immutable fingerprints will be saved with the diagnostic.
- `households/`: saved-policy and full stationary-distribution checks,
  supplemental zoom plots and three reproducible illustrative household lives.
  No new equilibrium solution or change to the standard 17 diagnostics.
- `credit_benchmark/run_v2/`: completed same-parameter experiment replacing
  artificial credit limits by net-estate solvency at every possible death and
  feasible continuation. One worker, 30-minute total/10-minute case caps;
  baseline replay, benchmark and exact repeat. All preference parameters,
  including the child-benefit level, remain fixed. No fertility renormalization.
  Conservative feasible-grid boundary and native value cutoff are explicit
  limitations; no continuum frictionless-credit certification. All three cases
  passed, exact benchmark repeat; see credit_benchmark/RESULTS.md. Birth renewal
  and estate funding are reported, not repaired. See plan.json for source pins.
- `lifecycle_dashboard/`: first-pass population age profiles versus actual ACS,
  PSID and CPS data, with complete CSVs and provenance. Supplemental diagnostic,
  not new calibration targets; empirical sample differences are documented.

Torch authentication succeeds at the start of this session. Local diagnostics
have a shared cap of three simultaneous model evaluators (two Jacobian plus
one benchmark); saved-checkpoint work uses one additional process without
solves. Every new experiment writes to a fresh folder; overnight contracts,
checkpoints and final results remain immutable. All numerical failures are
preserved. Full fits and parameters accompany every new reported calibration.

Half-hour supervision is active on the existing task heartbeat, with these new
budgets. It will pause again after the authorized jobs finish and are reported.

Credit source review: the authenticated purchase-income adapter enforces
`x+y/R >= -phi*pH` and owner saving `b_next >= -phi*pH`; changing phi to one
relaxes both coherently but retains renter debt limits. The two tests INSIDE the purchase routine are algebraically equivalent on
transaction-grid-supported states. This does not identify the purchase screen
with the separate end-of-period saving constraint, and must not be projected
back onto the September14 pre-income purchase restriction.
The estate audit explicitly leaves negative-estate creditor treatment unresolved.
A sure-repayment benchmark must respect liquidation solvency at possible death
dates; it cannot simply borrow to the numerical grid minimum. Full review was
read-only and did not change any economic mechanism.

## Production transition work authorized subsequently

The author subsequently requested a closed demographic endpoint, a transition initialized at the exact current calibrated steady state, speed testing and integration of ALL accepted recent changes into production code. The September 14 reference remains immutable; its timing/closure discipline is retained while later accepted primitives replace its old numbers.

- `credit_closed_endpoint/README.md`: closed same-parameter endpoint complete and repeated; complete fit/parameter tables and 17 plots retained.
- `credit_transition/preparation/PRODUCTION_INTEGRATION.md`: full object-by-object integration and historical reconciliation. Native entry/purchase/allocation, split entry queue and dated estate ledger implemented; native numerical replay pending before adoption.
- `credit_transition/smoke_v3/`: completed two-date diagnostic operator/cache test using the frozen reference runtime. Both arms are exact cache on/off; constant-path speedups are not evidence of a solved transition or guaranteed changing-path speedups. Earlier failed smoke attempts remain preserved.
- `jacobian/README.md`: completed 18/18 local sensitivity evaluations, with all targets and parameters and finite-step nonlinearity checks.

Production integration is not yet certified. It must pass native reference reproduction, benchmark reproduction, dated population/fiscal/estate checks and an actual price/pension transition solve with terminal-distance diagnostics. No old deadline is extended automatically by this new work.

## September 27 — credit-constraint clarification and author decision

Payment-to-income is deferred by the author, retained as an item to revisit;
no threshold, payment convention or model change has been adopted.

The lead previously conflated two duplicate transaction-screen tests with the
separate purchase and saving constraints shown in the talk. That statement was
too broad. The maintained September14 presentation source, lines173–175, shows
a pre-income down-payment test and an end-of-period collateral limit. The
September15 timing review in docs/model/POST_PRESENTATION_ISSUES.md, M06,
explicitly records trade before earnings: these are distinct restrictions.
Example: house100, financed share.8, initial cash10, gross return1.1, income40,
consumption5 and holding costs5 imply post-trade wealth−90 and final saving−69.
The final collateral limit−80 passes, but the cash down payment20 fails.

Current native code differs: define Q=pH, S=net sale proceeds (zero for a
renter buyer), x=b+S−Q, y=current income and R=gross period return. It tests
b+S >= (1−phi)Q−y/R and x >= −phi Q−y/R, the SAME purchase test, apart from
separate numerical grid support. Sources: solver.py3372–3378 and
kernels.py757–780 in intergen_eqscale_seq_optimized. Income is not added to x;
it enters the conditional budget c+b_next+o=R*x+y once (solver.py2794–2795),
where o denotes nonnegative owner holding costs. Final saving must separately
satisfy b_next >= −phi Q (native owner floor; solver.py165–168).

Under current positive-return, nonnegative-cost baseline primitives, feasible
consumption and final collateral imply x+y/R=(c+b_next+o)/R >= −phi Q/R
>= −phi Q. Thus the current income-aware purchase screen is implied by the
conditional household budget and final collateral in the continuous accounting.
It is not the same algebraic inequality: the reverse implication need not hold.
This is not a numerical regression test authorizing deletion of grid screening
or of the final collateral floor. No kernel or pinned run is changed.

The slide's indicator only on the down-payment RHS also leaves a condition on
nonmoving owners; the code bypasses purchase screens when tenure/product is
unchanged. Presentation should state the purchase condition only for buyers.
Historical sources and author text are preserved; reconciliation remains open.

### Primary-source lending comparison (September27; proposed, not adopted)

- DUE, section2.1.7, equations2.2–2.3: postpurchase liquid wealth b>=−phi*pH;
  between trades prohibit further borrowing if at/below ceiling, without forcing
  repayment after price declines. Primary project text: tmp/september_slides_review/
  greaney_reference.txt:440–460; https://www.nber.org/papers/w33512.
- Sommer–Sullivan–Verbrugge2013, equation4/footnote14: LTV applies to higher
  mortgage balances or changed housing; no mandatory principal reduction for
  nonmoving owners following price declines. https://kamilasommer.net/RentPriceRatio.pdf
- Boar–Gorea–Midrigan, Liquidity Constraints in the U.S. Housing Market, pp10–12,
  equations1–7: LTV and PTI at origination/refinancing; existing borrowers owe
  contractual payments, not new collateral tests.
  https://www.virgiliumidrigan.com/uploads/1/3/9/8/13982648/paper_bgm_v1.pdf
- Kaplan–Mitman–Violante2020, pp3294–3295: mortgage LTV/PTI at origination,
  contractual payments thereafter; separate one-period HELOC collateral limits
  do apply each period.
  https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf

Author concern is the distinction between origination and ongoing collateral.
Current native baseline marks the final debt floor to current house prices for
all owners, so it is stricter than DUE on inherited above-limit debt after price
falls. Lead recommendation to adopt a DUE-like distinction is only a proposal.
No author adoption and no code changes. PTI remains explicitly deferred.
