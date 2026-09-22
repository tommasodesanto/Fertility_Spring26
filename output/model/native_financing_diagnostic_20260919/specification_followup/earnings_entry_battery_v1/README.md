# Earnings and entry-wealth diagnostic battery

Author-authorized September 22 evening comparison; reviewed opinion scheduled for around 21:00 America/New_York, not earlier. This is a finite diagnostic, not model adoption or a converged calibration. Original submission at 19:00 EDT (superseded below): smokes A18288614, B18288617, C18288619, D18288635; dependent ten-worker production arrays A18288616, B18288618, C18288629, D18288688. All 244 remote plans and four zero-solve native preflights passed. See [submission](submission.json) and [lead launch review](launch_review.json). Each worker has one CPU, 32 GiB and one hour.

| Cell | Earnings risk, plus the common age profile | Entry wealth |
|---|---|---|
| A | One persistent AR(1), seven states | Zero assets at age 18 |
| B | One persistent AR(1), seven states | Inherited heterogeneous wealth marginal |
| C | Persistent AR(1) plus iid, seven by three states | Zero assets at age 18 |
| D | Persistent AR(1) plus iid, seven by three states | Inherited heterogeneous wealth marginal |

Both income processes have no fixed permanent type. All four retain the same age profile, current-income purchase accounting, parenthood housing floor, power equivalence scale, target system, nine structural search coordinates, bounds, and fiscal/fertility normalization. The earlier floor-versus-child-dependent-share experiment is deferred, rather than crossed with these four cells.

## What is recalculated

The single process matches the saved four-year log-earnings variance and first covariance: persistence 0.7345934905942886 and innovation standard deviation 0.4838308245314463. Its second covariance is 5.35% below the saved data moment. The two-component process retains the direct four-year estimates: persistence 0.7761446698586146, persistent innovation standard deviation 0.43743545312066157, and iid standard deviation 0.1649906105257306. These are period coefficients, not annual coefficients. Existing source data and bootstrap draws are reused; [external-estimate receipt](single_process_external_estimate.json) records transformed uncertainty and its limitations.

Each new process requires its own transition matrix, stationary weights, mean-one income grid and entry-income coupling. The deterministic age profile is unchanged. Every scored parameter point resolves the existing stationary equilibrium and fertility normalization. Paired parameter draws permit common-parameter contrasts as well as within-cell finite-search comparisons.

## Explicit assumptions and numerical limitations

Zero entry assets is an experimental restriction only. The heterogeneous alternative holds the inherited model wealth marginal fixed and remaps old total-income ranks onto the new persistent-income ranks; iid draws are independent conditional on persistent rank. This is a diagnostic coupling, not an estimated joint distribution.

The inherited five-node distribution comes from childless renter reference persons aged 18–24, using nonhousing net worth divided by annual family income. Multiplying these ratios by annualized gross model labor income retains a denominator mismatch. It is therefore an inherited empirical proxy, not a newly harmonized wealth estimate. The original checkpoint had no entry mass altered by grid clipping or value-based frontier censoring; each new solution still requires its own check. See [entry audit](entry_reference_audit.json).

The wealth grid stays at 160 nodes with upper endpoint 3000 in all cells, including the original 120 knots. The income grids are deliberately coarser than V5's 45 states. They reproduce targeted log moments but understate some continuous level covariances by about 10–11%. This separate exploratory plan allows a maximum 15% relative level-covariance approximation error; the older V5 5% contract is preserved. Household accounting, feasibility, market-clearing and value-quality gates are unchanged. Neither the coarse income grid nor the enlarged wealth domain is certified as converged. A finer-income-grid comparison and exact repeats remain necessary before adoption. Finite infeasibility values in continuation interpolation remain unresolved.

## Bounded execution

Four common-parameter smoke cases, followed only on each cell's own success by ten single-thread workers per cell (40 workers total). Each production worker receives at most six predeclared parameter vectors and one hour; it stops if its remaining budget cannot accommodate another case based on observed runtime. Maximum 244 full objectives including smoke, each with at most eight nested stationary solves (1,952 theoretical maximum; runtime makes the realized count much smaller). Native/wrapper/controller case caps are 3,100/3,150/3,200 seconds. Every case writes progress and latest/best receipts. Failed smoke cells block their dependent arrays; there is no fallback entry law or relaxed gate.

The final comparison will show the complete target and actual-bound parameter tables, standard figures, valid/rejected/incomplete counts, and common-parameter contrasts. Raw scalar loss alone is insufficient to choose the specification. Exact repetition and richer-grid work will be labeled completed or unrun from actual receipts.


Resource-only correction at 19:07 EDT: the original CPU-short jobs never started because that partition limits each user to 120 GiB and unrelated jobs already reserved 96 GiB. Scheduler rejected an in-place partition update. Only the eight pending battery records were cancelled, and the identical immutable bundle/budgets were submitted to the allowed standard CPU partition `cs`: A smoke18288966/array18288967; B18288969/18288970; C18288971/18288972; D18288973/18288974. These are the live IDs. Original receipts remain preserved; see [resource plan](resource_requeue_plan.json), [new submission](submission_cs.json), and exact scripts in `launch/`. No numerical model was restarted or changed.


At 19:24 EDT, four bounded income-resolution checks were submitted: A18290336, B18290337, C18290338, D18290340. They repeat the common smoke parameters at 15 states for the single AR1 and 45 states for AR1+iid, with the same wealth grid and all economic parameters, targets and household gates. A/B start after their reviewed completed smoke receipts; C/D use scheduler afterok dependencies on their coarse smokes. Each receives one hour, one CPU and 32 GiB, with at most eight stationary solves. Actual level-covariance approximation errors are 5.13% and 4.61%; this is a resolution sensitivity check, not a numerical-convergence certificate. The four zero-solve preflights passed. See [resolution plan](resolution_plan.json), [submission](resolution_submission.json) and saved launch scripts. The first submission attempt issued no jobs because completed A/B parent IDs had already left the scheduler; its receipt is preserved.


The 19:32 collection is an interim snapshot, not the scheduled author review. All four common-parameter smokes completed and scored; all four production arrays are released. The verified snapshot is [readout_1932](readout_1932/README.md), with complete tables and standard figures. Independent checkpoint audits `smoke_A_checkpoint_audit.json` through `smoke_D_checkpoint_audit.json` record zero entry censoring and planned-versus-realized entry wealth marginal L1 gaps below 2.6e-16. These statements concern the saved smoke checkpoints; they do not establish entry-law suitability or freedom from continuation-interpolation effects on choices. Final selection and numerical-resolution comparisons remain pending.

The reusable workflow uses `code/model/tools/prepare_e5f_earnings_entry_battery.py`, `run_e5f_earnings_entry_battery.py`, `collect_e5f_earnings_entry_battery.py`, `audit_e5f_earnings_entry_checkpoint.py` and `build_e5f_earnings_entry_review_pdf.py`; exact cluster launch scripts are saved in `launch/`. The collector requires a fresh output directory and verifies pinned input contracts, complete tables and original diagnostic hashes before selection. The PDF is generated from that verified readout and a separately reviewed narrative. Focused constructor/builder/controller checks: 22 passed; launcher shell syntax and reporting-script compilation passed.

The separate low-value audit finds zero saved beginning-of-period mass with V<=-1e6 in A/B/C, but 2.06e-11 in D (5.18e-15 at V<=-1e9). This tiny mass is reported explicitly; passing the scored gates does not eliminate the continuation-sentinel concern. See [interim review](interim_review.json).

The D richer-income check (45 states with inherited wealth, job18290340) failed after one started and zero completed stationary solves: at age30, dead-node mass2.50e-12 exceeded the unchanged1e-12 gate. The coarse21-state smoke passed, so the result is sensitive to income support/resolution. This tiny failing mass is neither silently discarded nor treated as proof of economic impossibility. Failure evidence is preserved at `resolution_failure_D/`; no retry or gate change, and no dependents require cancellation. Other checks and the separately authorized coarse diagnostic continue.


At the 19:42 resolution collection, A (single AR1 with zero entry wealth) is verified at 15 states; B/C remain running and D is the preserved failure. [Complete A finer-grid tables and original figures](resolution_readout_1942/README.md) and [all 13 common-parameter moment changes](resolution_interim_review.json) are retained. Relative to seven states, mean first-birth age rises by 0.523 years and childlessness by 1.85 percentage points, while mean rooms changes by only -0.006. Thus the simple process can solve at richer resolution, but seven-state numerical adequacy is not established. No substantive author review has been delivered early.

For the final reader PDF, pass both the fresh main readout and `--resolution-readout` with the fresh resolution readout to `build_e5f_earnings_entry_review_pdf.py`. The optional comparison page checks identical objectives and structural parameters before showing all13 finer-minus-coarse moment changes. The new page was rendered and visually checked; the final complete PDF still requires its final visual QA. Preview files remain under `tmp/pdfs/` and are not reader deliverables.
