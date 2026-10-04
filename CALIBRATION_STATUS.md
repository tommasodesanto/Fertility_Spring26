# Calibration and model status

**Reconciled October 3, 2026.** This is the consolidated current-state note. The
complete previous status, previous memory, reconstruction evidence and preservation
manifest are in [the context archive](calibration_archive/context_refresh_20261002/README.md).
Read historical chronology there, in dated daily notes or in the named experiment
packets when needed.

**Working continuation convention, author-adopted October 3:** net housing
transactions occur after interest on beginning net financial wealth,
\(b'=Rb+S-Q+y-c-K\). The verified soft-constraint, old-target post-interest
chain 13 (loss **13.771131463467462**) is the working continuation anchor,
not a certified paper baseline, global optimum, or grid/transition validation.
The earlier original-timing soft selected point (loss 23.078309) and chain 15
(loss 18.445305) remain comparison evidence. The matched timing search is
terminal. Verified
as of **October 3, 2026, 06:00 New York**, arrays **19086987 and 19087556**
have 48 terminal chains: 46 passed fresh native selected-point checks and two
(chain 23 in each arm) found no admissible candidate. The lowest verified
original-timing loss is **18.44530519407432** (chain 15); the post-interest
loss is **13.771131463467462** (chain 13). These searches do not certify
optimizer convergence.
An isolated post-interest-timing recalibration with a narrower experimental
PSID wealth numerator has also finished: ten of ten chains passed native
verification; its lowest new-contract loss is **48.170707377609034** (chain 2).
The target contracts differ, so these native losses cannot be ranked directly.

## Reference identities and navigation

| Object | Identity and role | Authoritative evidence |
|---|---|---|
| Frozen paper reference | `paper-baseline-2026-09-14`; retained checkout `tmp/paper_baseline_sep14/`. Preserve source and results. | `PAPER_BASELINE.md` in that checkout; [baseline checker](code/model/tools/check_paper_baseline.py) |
| September 28 fixed-economics reference | Older equilibrium/normalization and utility objects. A refactor oracle, not interchangeable with the current soft calibration. | [refactor report](output/model/publication_refactor_20260929/REPORT.md), [refactor runtime](code/model/refactor_lab/README.md) |
| Earlier soft selected point | Original timing; chain 16 / case 0046; verified loss 23.078309. Historical comparison, not the continuation anchor. | [selection](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json), [verification](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_verification.json) |
| Working revised-timing anchor and matched comparison | Soft constraint; post-interest chain 13 with the retained 6.926584 PSID wealth target, verified loss 13.771131. Original-timing chain 15 loss 18.445305 is comparison evidence. Both arms had 24 matched starts, 48 terminal chains, 46 verified. | [48-chain receipt](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json), [complete fit and parameter readout](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/RESULTS.md), [driver plan](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/driver_plan.json) |
| Experimental wealth-numerator comparison | Post-interest timing with one new PSID aggregate wealth/earnings target; ten verified chains, lowest new-contract loss 48.170707. No adoption. | [ten-chain verification](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/collection/verification.json), [three-arm comparison and full tables](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/comparison/COMPARISON.md) |
| Historical purchase-rule comparison | Fresh hard and quarter results: 88.588403 and 48.319938. Older policy exercises use earlier points with losses 97.011220 and 51.556036. | [fresh calibration readout](output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/README.md) |
| Historical transition initializer | Normalized-v1 chain 20 / case 0028_nm, loss 30.371888. Its one-shock transition attempt failed acceptance. | [transition deployment/readout](output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/v2/README.md) |

The canonical local stationary-GE implementation and editable runner are indexed in [`code/model/README.md`](code/model/README.md) and [`code/model/production/README.md`](code/model/production/README.md); the current Python interface and commands are in [MODEL_PLAYGROUND.md](code/model/tools/MODEL_PLAYGROUND.md). The October 3 local deployment passed scoped verification: fresh same-host original-reference comparisons for unchanged chain 13 and annual \(\beta-0.001\) matched all 91 arrays, every field in the 14-row target table, all 31 numeric parameter values and bounds, and all 17 standard-plot hashes; only 13 enumerated descriptive role/status metadata differences remain ([unchanged](output/model/production_deployment_20261003/compare_unchanged/comparison.json), [beta](output/model/production_deployment_20261003/compare_beta/comparison.json)). The default cached case produced the standard 17, policy 8, aggregate 7, and explorer artifacts. The \(\sigma=2.01\) GE input experiment and algebraic scale exercise are diagnostics, not recalibration or adoption ([input receipt](output/model/production_deployment_20261003/external_input_propagation.json), [scale receipt](output/model/production_deployment_20261003/scale_exercises.json)). Four transition modules import, but no dated transition was replayed. No new calibration or baseline was adopted; chain 13 remains the working continuation anchor, and its empirical identification and grid adequacy remain open. Authenticated observer and oracle bundles remain read-only dependencies; only superseded frontends were archived. See the production guide for exact limits and commands.

## Working economic and accounting contract

**Isolated birth-menu and Estate-A comparison, verified October 3, 2026,
17:32 New York.** The author selected Estate A for this isolated test:
\(W=b'+(1-\psi)Ph'\), with selling cost \(\psi=0.06\) and no additional interest on \(b'\), applied consistently
in utility, native death-flow accounting and the empirical bequest observer.
There is no adult estate-recipient mapping. Both one-intended-birth and
up-to-three-intended-birth arms retain all ten current chain-13 coordinates,
fixed inputs, timing and grid. The only empirical target change for common
scoring is aggregate wealth/earnings **4.45838713455674**; every other target and
weight, including the provisional SCF bequest target, is retained. This is a
fixed-parameter experiment, not recalibration or a promoted default.

| Intended-birth cap | No A, common new-target loss | Estate A, common new-target loss | Ownership at age 82, no A | Ownership at age 82, A |
|---|---:|---:|---:|---:|
| One | 53.064 | 57.595 | 98.739% | 96.225% |
| Three | 1714.289 | 1728.298 | 97.897% | 95.010% |

The A runs' old-wealth-target diagnostic losses are **16.895** and **1685.273**;
these are distinct scoring contracts. The [complete four-arm fit and parameter
readout](output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md)
contains every 14-row fit and 31-row parameter table, plus raw/net estates and
age-82 financial saving. [Source identities and exact
results](output/model/experiments/birth_count_choice/estate_a_v1/comparison.json)
and [estate diagnostics](output/model/experiments/birth_count_choice/estate_a_v1/paired_estate_diagnostics.csv)
are retained. A reduces terminal ownership but does not eliminate its high
level or settle SCF wealth-scope, recipient or creditor comparability. The
one-birth A price and ownership match the independently saved Claude-A solution;
its utility-only gross-flow report is replaced by a net-flow calculation in
[the separate comparison](output/model/experiments/birth_count_choice/estate_a_v1/claude_a_comparison.csv).

Both A GEs passed native acceptance and exact repeats, using nine and eight
lifecycle solves, respectively, and each wrote all 32 standard/policy/aggregate
plots. The final [35-test log](output/model/experiments/birth_count_choice/estate_a_v1/tests_verified.log)
passed; [36 production source hashes](output/model/experiments/birth_count_choice/estate_a_v1/source_snapshot.json)
remained unchanged. All-zero negligible dead menus retain the existing 1e-12
mass gate. Descriptive count hazards use pre-birth exposure; active empirical
target observers remain separate. The original no-A count experiment and its
full fit remain in [current_params_v1](output/model/experiments/birth_count_choice/current_params_v1/RESULTS.md).

**Matched Estate-A recalibration: reviewed v3 smoke passed; production array
19127370 submitted October 3, 19:21:31 New York; all ten tasks RUNNING at
the latest scheduler check.**
Smoke array **19125188**, tasks **0 and 5**, has a 1.5-hour limit, two objective
calls plus fresh selected-point GE verification per arm; both smokes must pass
before production release. Both tasks completed with exit code 0, in 15m40s
and 15m52s. Remote `verify_smoke_gate.py` passed with
`matched_both_arm_smoke_passed` for inventory SHA
`d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe`,
target fingerprint `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`,
and weight fingerprint `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`.
The original `collect_torch.py` failed on duplicate `arm` and `chain` keys.
A collector-only fix asserts the receipt identity and merges its fields once;
its SHA-256 is `2ed523f66299ed6542ff01adde80bb563137d7f6bcc29aa27f257506d6529967`.
The reviewed copy ran outside the immutable stage, passed both actual receipts,
identity-negative tests and remote strict collection. The
[reviewed smoke collection](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/smoke_collection_reviewed.json)
contains the full 14-target fits and 31 parameter records for both arms.
Independent review checked every target row, parameter and ten free bounds,
gaps, weighted losses and sums; the selected root and fresh repeat matched on
14 new-target rows, 31 parameter rows and 17 PNG hashes. Native fresh-child,
exact-repeat and search-full gates passed. These are two-call smoke diagnostics,
not optimized calibration results. The earlier failure is retained in the
[historical blocker receipt](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/monitor_blocked.json).
The [smoke submission receipt](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/smoke_submission_receipt.json)
and [actual-context preflight](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/context_preflight_receipt.json)
pin the 412-file v3 package (inventory SHA
`d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe`).
Both arms built their actual reporting contexts with zero solves. The reporting
adapter now installs count corrections once per actual module and rebinds the
facade for each context; 35 tests pass, including three successive count-menu
contexts. The retained first-call local A GE results above remain valid.

The [start plan](output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json)
uses five identical starts per arm, ten tasks total, population-one GE and beta
range \([0.93,0.99]\). Starts include the preserved provisional new-wealth
candidate (chain 6, `0060_nm`, search loss **22.141841386410267**), current
anchor and three perturbations. The provisional search value awaits fresh native
verification; it is distinct from the earlier fully verified new-target loss
**48.170707377609034**. [Complete paused fit, bounds and provenance](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/PAUSED_BEST.md)
are retained. Previous old calibration searches are stopped with no automatic
restart. Smoke submissions and zero-solve checks are validation, not completed
calibration results; earlier smoke chronology is in the October 3 daily note.

The duplicate-guarded `submit_torch.sh` passed its 412-source inventory pins,
both smoke gates and storage reserve, then submitted production array
**19127370** exactly once (`0-9%10`; five matched starts per cap). Each task has
one CPU, 24 GiB, six hours, at most 500 objective calls and a 1,800-second
final native reserve. All ten tasks were RUNNING at about 19:22:40 New York,
elapsed 1m09s each, on `cs602`, `cs604`, `cs617` and `cs629`. See the
[production receipt](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/production_submission_receipt.json)
and [release review](output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/production_release_review.json).
This launch is experimental; it does not adopt a new baseline or certify
optimization convergence.

**Overnight Estate-A continuation queued October 3, 23:34 New York.**
Controller **19136605** is pending `afterok:19127370`; the original ten-task
array remains running. The separate [continuation stage](code/cluster/estate_birth_calibration/continuation/README.md)
and [submission receipt](output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/deployment/controller_submission.json)
pin the same economic inputs, ten free-parameter bounds, target and weight
fingerprints, and numerical gates. If every parent task exits successfully with
a fresh native exact repeat, each new chain will begin from its **own** verified
endpoint with a fresh optimizer simplex. The one-core controller first runs two
native optimizer calls and fresh selected-point verification in each arm; only
if both pass does it release up to ten one-core production tasks. Search retains
the 500-call cap per task, reserves 1,800 seconds for native verification and
stops by October 4 **10:00 New York** (epoch `1791122400`). The production job
ID is not yet assigned. All parent-final and continuation-native gates remain
pending; this is an experimental search, not an adopted calibration or exact
optimizer-state resume. No failed chain will be restarted automatically.

**Five additional count-three starts queued October 4, about 00:32 New York.**
The isolated [count-three expansion](code/cluster/estate_birth_calibration/count3_expansion/README.md)
controller **19139361** is PENDING on `afterok:19136605`; both controllers
wait behind parent array **19127370**, whose ten tasks remain running. The
new controller will select the best verified binary and count-three parent
endpoints and construct five bounded, nonduplicate search starts (binary best,
the midpoint, and three deterministic count-three perturbations). All ten
parameters, bounds, target values, weights, economics and native gates remain
unchanged. One fresh two-call count-three native smoke must pass before a
five-task array can be released. Combined with the existing continuation,
this caps Estate-A production at **15 one-core tasks**. The extra tasks retain
500 calls each, a 1,800-second fresh-native reserve and the same October 4
**10:00 New York** absolute stop. A 7,800-GiB free-space gate covers the
7,500-GiB combined retained-case planning maximum; the actual staging check
found about 434,526 GiB free. The new smoke, final parent verification and
five-task production remain pending. No production job ID, adopted result,
automatic retry or exact optimizer-state resume is claimed.

The hourly heartbeat **Monitor matched estate and birth calibrations**
(automation `finish-and-monitor-soft-timing-calibration`) in chat
**Compare interest timing and credit** (`01a0ff53-5843-73a2-aaf3-0e1f1313e91f`)
is ACTIVE after the reviewed release. No automatic retry
is authorized.

These are the retained objects behind the soft-constraint comparison. Numerical
examples explicitly sourced to the earlier 23.078309 point below remain
historical; the working revised-timing anchor is chain 13. Neither resolves
every empirical or publication issue.

**Units and household choices.** One period is four years. Households choose
consumption, housing, tenure, saving and fertility over the lifecycle. Children
ever born, `n`, and children currently at home, `m`, are different states.
The highest birth-count state represents three or more children. Housing is
measured in rooms; owner products have 2, 4, 6, 8 or 10 rooms, and rental
housing has the retained six-room cap. The active solution has 120 asset nodes and
nine income states. Legacy labels containing “B15” do not establish 15 executed
income states.

**Utility.** Retained risk aversion is \(\sigma=2\), the consumption share is
\(\alpha=0.733\), and there is no adopted child-dependent consumption-share
change. The equivalence scale is
\[
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]
The parenthood housing floor is a physical-room object applied before the housing
taste parameter \(\chi\). The child benefit is \(\psi m^{1-\gamma}\), with zero
benefit at \(m=0\); the first-birth cost and bequest motive are separate objects.
Historical lower-\(A(m)\), renter-borrowing and other utility probes are
diagnostics, not silently adopted components of this reference.

**Normalized CES-limit share test; four-chain array 19133352 running; verified October 3 at 22:32:52 New York.** The
experimental composite denominator applies in all family states:
\(Q=c^{\alpha(m)}s^{1-\alpha(m)}/[\alpha(m)^{\alpha(m)}(1-\alpha(m))^{1-\alpha(m)}]\).
There is no childless-share numerator or reference-rent correction. The
author-authorized rule sets \(\alpha(0)=.733\) and
\(\alpha(m)=\operatorname{clip}(.733-\delta_{\rm jump}-\delta_{\rm slope}m,.05,.95)\)
for parents; both parameters are free on \([0,.25]\), and \(h_P=0\).
The experimental target contract has 11 free coordinates and 11 scored moments
among 14 rows. It promotes `family_rooms` (target 0.38509964969278165, weight
280.52808370152104); national uncertainty is unavailable, so this is the
inherited 42-metro bootstrap weight, and the model-dependent-child observer
remains a proxy.

The post-interest chain-13 starting reference retains old wealth target
6.92658379107299, timing, earnings, birth architecture and other economic
inputs. There is no `rstar` or `alpha0` numerator, added birth shock, or utility
cost rescaling. The original starting guess did not bracket a stationary root.
The four starts now use `first_birth_fixed_cost=1.9`, still estimated on
\([0,8]\); this changes a starting guess, not a restriction or utility formula.
V5 includes the missing native reference artifacts directly in its immutable
3,008-file source. Inventory SHA-256 is
`39fe1c3316919fb3d8672030afda6fb5fcc29aef0969181d206c9fc0757ca127`.
Three actual contexts passed with zero lifecycle solves
([receipt](output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt5/preflight/completed.json)).
V4 smoke **19132298** passed two distinct numerical cases, fresh selected
native postcheck and exact repeat. Its [full starting-fit tables](output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt4/smoke/RESULTS.md)
are a feasibility check, not a calibrated result. Final V5 smoke **19132940**
passed its exact-loop gate and full native collector
([complete tables and parameters](output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt5/smoke/RESULTS.md)).
The search requires 900 seconds of GE time after reserving 1,800 seconds for
final verification. Array **19133352**, `0-3%4`, runs four six-hour chains on
one CPU and 24 GiB each, up to 500 objective calls and 32 lifecycle solves per GE
([submission receipt](output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt5/submission_receipt.json), [four active first-GE heartbeats](output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt5/launch_health.json)).
The reviewed launcher SHA is
`e092a10197e3269dce0acfee125a7dc408b74b68026738f1cbe9abfc530ec05f`;
it repairs only EXIT-receipt bookkeeping, which caused the earlier Slurm
failure after successful numerical verification. Its four success/failure
receipt fixtures passed locally and on Torch. Model and search source remain
the immutable v5 package. All 17 final smoke plots match the inspected v4 packet.
Authenticated unbracketed candidates receive a numerical penalty; budget exhaustion
stops search and verifies the best completed result. Other errors are terminal.
Local checks are 11 passed/1 skipped; all four mock loops passed.
No overnight calibration result or adoption is yet established.
The native postcheck requires 11 coordinates, full 14-row experimental target
CSV, 31 parameter records, 17 plots, and exact repeat.
See [experiment specification](code/model/experiments/ces_normalized_shares/README.md),
[deployment workflow](code/cluster/ces_normalized_shares_calibration/README.md), and
[complete proposed targets, bounds and starts](output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json).

Historical-vector diagnostic: September 27 de_0093 supplies all eleven transferable coordinates. The current numerical price caps rejected the point. A separate wider-cap diagnostic found root 7.670590398 with exact repeat and identical economic-input fingerprint, but derived housing supply coefficient 0.198580 is below the retained 0.2 lower bound. Diagnostic loss 4417.699 is not an accepted calibration candidate. Fertility fits closely while housing collapses; [full 14-row fit, 31 parameters and plots](output/model/experiments/ces_normalized_shares/overnight_v1/historical_start_diagnostic/README.md). No production bound, utility, target or running search was changed.

**Earnings and entry.** Earnings use a deterministic age profile and a single
persistent four-year process: \(\rho=0.7345934906\), innovation standard deviation
\(0.4838308245\). There is no separately adopted permanent-type plus transitory
decomposition. See the [external earnings estimate](output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/single_process_external_estimate.json).

The selected input uses the provisional `nonnegative_mean` entry-wealth mapping:
negative five-bin mean wealth ratios are floored at zero after binning, then
positive means are rescaled to preserve the retained mean. Effective ratios are
\([0,0,0.037866,0.127705,1.127294]\). Mean entry wealth is \(0.186520\);
mean annual entrant income is \(0.720554\); the zero-wealth share is \(0.666948\).
This is not a recoding of every raw survey observation. Tommaso provisionally
preferred this option; production paper adoption and donor/recipient accounting
remain separate decisions. Earlier negative-entry-cell infeasibility is not a
current blocker at this point. The executed entry contract is in
[the selected input contract](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/input_contract.json).

**Finance and timing.** Annual real interest is 0.02, giving
\(R=1.02^4=1.08243216\). The financed share is \(\phi=0.8\), so a down payment
uses \(1-\phi\). Renter debt capacity is zero. Selling cost is 0.06.
Annual depreciation \(0.0141614372\) is compounded to a four-year rate
\(0.0554537908\); annual property tax \(0.0105983608\) is multiplied by four to
\(0.0423934431\). These conventions are explicit, not interchangeable
annualizations.

Let \(b\) be beginning net financial wealth, \(S\) net sale proceeds, \(Q\) the
purchase cost, \(y\) current income, \(c\) consumption and \(K\) other costs.
The historical original-timing convention is
\[
b'=R(b+S-Q)+y-c-K.
\]
The author-adopted working convention is
\[
b'=Rb+S-Q+y-c-K.
\]
The revised convention moves current net transaction financing outside the
interest factor.
At held choices the difference is \((R-1)(Q-S)\). Sale, forward-budget and
solvency maps must use the chosen convention consistently. Existing debt is
already in \(b\); do not subtract it again. The applicable owner ending-debt
floor remains \(b'\geq-\phi Q\). Tommaso adopted the revised timing for
continuation on October 3 after reviewing both recalibrated arms. Hard versus
soft purchase constraints remain a distinct specification choice.

**Demography and closure.** The earlier original-timing soft equilibrium fixes
household scale
\(N=1\), solves price using stationary birth renewal, and derives the housing
supply scale \(H_0=6.757074\) from market clearing. \(\psi\) is jointly free
in the ten-parameter calibration; it is not independently normalized back to
2.1 at each proposal. The unscored normalization row checks adjusted births
divided by entry.

The top-count correction adds \((w-3)\) times the flow into the three-or-more
state to recorded births, where \(w\) is its representative count. Potential
entry is adjusted births divided by 2.1, once. A half-age-16 / half-age-20
entry queue is implemented, approximating mean entry age 18. Child departure
uses the retained \(2/9\) rule, without a newborn exemption; parental death
removes dependency without an additional entry flow. See
[adult-entry accounting](code/model/intergen_eqscale_seq_optimized/adult_entry.py).
This is a working household demography approximation, not explicit tracking of
individual offspring genealogies or certified resident-person accounting.
The October 2 check found no extra factor of two in this path; different
three-plus weighting and age windows still require care when comparing fertility
statistics.

At that earlier selected point, entry is \(0.0617334562\), adjusted births
\(0.1296402579\), the birth-renewal residual is approximately
\(-7.08\times10^{-10}\), the housing residual is zero, and the PAYGO residual
is \(2.71\times10^{-14}\). Payroll tax is \(0.0802807096\) and the period
pension is \(0.9177840475\). The accepted CPS ASEC pension/gross-earnings ratio
\(0.2294460119\) underlies the retained PAYGO normalization; a transition
holds the adopted payroll tax and balances period pensions, rather than adding
unrecorded government spending. See
[verified closure](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/selected_postcheck/phase_b_ge/selected_repeat_final/closure.json).

Housing supply elasticity 0.63 is externally fixed provisionally.
The current equilibrium is a closed single housing market, with no adopted
outside-entry/geographic migration valve. For a policy endpoint under the
retained supply function, the reported population scale uses
\(N=H_{0,\mathrm{base}}/H_{0,\mathrm{derived}}\). Holding \(H_0\) fixes the supply
schedule, not the quantity at every price. This population-scale effect must
not be described as a total fertility effect.

## Earlier original-timing selected fit and parameters

The earlier local selected-point repeat passed, reproducing loss **23.07830929416065**.
The largest saved-moment discrepancy was \(5.33\times10^{-15}\). The packet has
14 target rows, 31 parameter records and the established 17 plots. A first
cross-platform exact-double comparison was rejected for roundoff; the accepted
local repeat did not relax economic acceptance gates. This is selected-point
verification, not proof of optimizer convergence.

The gap below is model minus target. Blank weights identify the separate
normalization; zero weights identify validation rows. Rounding is for display;
the linked CSVs retain full precision.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Stationary birth-renewal normalization | 2.100 | 2.100 | -1.486e-09 | — | — |
| Childless women, ages 40–44 | 0.198 | 0.202 | 0.004 | 35532.304 | 0.553 |
| One child among mothers, ages 40–44 | 0.214 | 0.213 | -0.001 | 26952.821 | 0.028 |
| Mean age at first birth | 25.976 | 26.032 | 0.056 | 139.828 | 0.436 |
| First births at ages 30+ (validation) | 0.249 | 0.234 | -0.015 | 0 | 0 |
| Net wealth / labor earnings | 6.927 | 6.551 | -0.376 | 7.595 | 1.073 |
| Child-directed bequest flow / wealth | 0.007 | 0.007 | -4.449e-04 | 5165289.256 | 1.023 |
| Older-household wealth dispersion (validation) | 3.516 | 2.968 | -0.548 | 0 | 0 |
| Mean occupied rooms | 5.729 | 5.957 | 0.228 | 128.021 | 6.639 |
| Ownership, heads ages 30–55 | 0.676 | 0.665 | -0.011 | 2339.362 | 0.282 |
| First-birth room response | 1.465 | 1.324 | -0.141 | 137.565 | 2.726 |
| Rooms: 3+ versus 1–2 children (validation) | 0.385 | 0.304 | -0.081 | 0 | 0 |
| Recent-parent ownership gap | 0.128 | 0.118 | -0.009 | 27055.823 | 2.371 |
| Children ever born, capped at 3, age 25 | 0.810 | 0.528 | -0.282 | 100.000 | 7.947 |

Source: [complete target-fit CSV](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_target_fit.csv).
Early fertility and mean rooms account for approximately 63% of loss. Exactly
one child is conditional on being a mother, not an unconditional share of women.
The ten positive-weight moments match the ten free parameters in count;
informative rank and weak identification at this point have not been certified.

| Free parameter | Estimate | Lower | Upper | Near bound |
|---|---:|---:|---:|---|
| `beta_annual` | 0.967 | 0.940 | 0.990 | no |
| `chi` | 1.097 | 0.100 | 5.000 | no |
| `first_birth_fixed_cost` | 0.352 | 0 | 8.000 | no |
| `kappa_fert` | 0.124 | 0.020 | 50.000 | yes |
| `kappa_fert_continuation` | 0.363 | 0.020 | 50.000 | yes |
| `theta0` | 0.106 | 0 | 8.000 | no |
| `child_benefit_curvature` | 0.101 | 0 | 0.800 | no |
| `tenure_choice_kappa` | 0.014 | 0.001 | 0.100 | no |
| `psi_child` | 0.178 | 0.010 | 0.500 | no |
| `h_P` | 2.504 | 0.100 | 2.600 | no |

Source: [all 31 parameter records](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_parameters.csv).
“Near” follows the saved screen: within 1% of the search interval, not actual
endpoint contact. The derived \(H_0\) has admissibility range \([0.2,80]\);
\(\theta_1=0.0081930841\) remains fixed. Fixed objects and the executed input
contract are linked above.

## Empirical target provenance and measurement limits

The row-by-row [provenance review](calibration_archive/context_refresh_20261002/target_provenance_review.json)
preserves builders, source records, estimators, samples, controls/clustering,
uncertainty where available, current weights, observer definitions and warnings.
It is a reconstruction of the executed contract, not a new target registry or
a certification of statistical design. Governing reviewed sources are
[the September 27 measurement review](docs/model/e5f_target_measurement_review_20260927.md)
and [accepted-input reconciliation](docs/model/accepted_input_reconciliation_20260926.md).

- **CPS fertility stocks:** pooled June 2004/2006 women ages 40–44, valid
  `FREVER` 0–20 and positive supplement weights. Childlessness is unconditional;
  exactly one child is conditional on mothers. Official annual generalized
  variance approximations exist, but pooled covariance is not design-certified.
  The model uses uniform within-period birth timing and a reproductive household
  member proxy, not a demonstrated female-exposure reconstruction.
- **Early fertility:** the same supplements at exact age 25, with children ever
  born capped at three. The recorded person bootstrap standard error is 0.028035,
  stratified by year; it is not CPS design-consistent. Weight 100 is a working
  calibration choice, not inverse survey variance. The model's uniform
  age interpolation is an approximation.
- **NCHS timing:** first births in 2003–2006, ages 12–49, mapped to retained
  four-year age-cell representative values. Annual dispersion, including
  0.084567 for mean age, is not a sampling standard error. The age-30-plus
  share is validation only.
- **PSID wealth:** 2005/2007 weighted net-wealth to gross labor-earnings ratio;
  wealth ages 18–85 and earnings ages 18–65 under the retained sample definitions.
  Reference-person bootstrap standard error is 0.417310. Wealth includes home
  equity. Older-household wealth dispersion is validation only; the retained
  PSID income filter/age support and model pension denominator do not coincide.
- **SCF bequests:** the 2007 mortality-weighted child-directed annual estate flow
  divided by net wealth is 0.007291. The empirical calculation is complete;
  recipient eligibility, timing and model entry mapping remain open. No survey
  standard error was established; a synthetic working uncertainty is not
  empirical precision.
- **AHS rooms:** 2007 occupied households, head ages 18–85, positive weights and
  valid rooms; 37,793 observations, literal public-use room topcode 21. Fay BRR
  standard error is 0.008934. The retained loss weight 128.020702 is not the
  inverse of this variance, and model room exposure remains an approximation.
- **ACS housing/ownership:** accepted national 2005/2006 targets. Ownership and
  recent-parent rows retain a DUE structure restriction absent from the model;
  family rooms uses resident own children, a minor-child screen and rooms capped
  at nine. The old 42-metro bootstrap weights were retained; they do not
  establish national sampling variance. The recent-parent model observer uses
  births into previously empty-dependent homes versus currently empty homes,
  including former parents. It is not the exact ACS oldest-child-age estimator.
- **PSID first-birth rooms:** the authoritative A2h Sun–Abraham contrast is
  calendar years \(+3/+4\) versus \(-3/-2\), baseline reference persons/spouses,
  biological first births and confirmed-childless controls. It uses person/year
  fixed effects, age/education controls and individual clustering: 117,853 fitted
  observations, 9,310 clusters, standard error 0.050270. Tommaso selected the
  rounded 1.465 target for this working calibration. The model destination-period
  contrast does not replicate the empirical panel estimator or cohort selection.
  Historical 0.600/0.720 targets are not current. Contemporaneous-income-controlled
  OLS is a robustness exercise: 1.300 versus approximately 1.461 without income on
  the matched sample. It has not replaced the target.

These limits remain visible. Removing, reweighting or replacing a target needs an
identification argument and an explicit contract change; this refresh makes none.

## Active timing comparison: execution state and acceptance

The deployment receipt identifies remote root
`/scratch/td2248/projects/soft_timing_calibration_20261002_v2`.
Attempt-one smoke 19085105 failed repeat verification in a reused interpreter.
The repaired driver uses a fresh child interpreter within the same absolute
deadline. Attempt-two smoke 19086529 passed both arms, including the full
14-row fit, 31-parameter record, 17-plot packet and selected-point repeat. No
economic change or gate relaxation was introduced by this process repair.
The production array 19086987 was subsequently submitted.

The alternative timing at the common starting coordinate has loss
**57.186333** and price **0.730719**, versus **23.078309** and **0.726639** for
the original reference. This is a held-coordinate comparison after equilibrium
solution, not a ranking of separately recalibrated specifications. Both complete
native smoke tables and parameter records are linked by
[the smoke collection readout](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/smoke_collection/readout.json).

The author subsequently expanded the design to 48 chains, 24 matched starts per
arm. Array 19087556 adds 40 chains (indices 4–23 per arm) in remote v3. The
[expanded plan](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/expanded_start_plan.json)
retains bounds, targets, gates and the passed v2 economic implementation. Starts
combine the original four, seven other historical soft candidates, eight nearby
variations and five broader starts. The 240-chain proposal was superseded.

**Final matched search, verified as of October 3, 2026, 06:00 New York.** Both
arrays are terminal, with zero active Slurm tasks at collection. Of 48 chains,
46 passed the fresh native selected-point and exact-repeat gates. Original and
alternative chain 23 ended with no admissible selected point; neither enters
the winner comparison. The lowest verified original-timing loss is
**18.44530519407432** (v3 chain 15), and the lowest verified
post-interest loss is **13.771131463467462** (v3 chain 13). Both winner packets
have 14 target-fit rows, 31 parameter rows, 17 standard plots in each native
root and exact repeat, and exact repeated tables and plot hashes. See the
[complete 48-chain collection receipt](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json),
[original fit](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/original_target_fit.csv),
[alternative fit](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/alternative_target_fit.csv),
[original parameters](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/original_parameters.csv), and
[alternative parameters](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/alternative_parameters.csv).
The selected original estimate of \(h_P=2.6\) contacts its upper bound; the
alternative estimate is \(2.59376\), near the same bound. Early fertility is
0.533079 and 0.533805, respectively, against a target of 0.809528. The fit
improvement therefore does not resolve this target, and neither optimizer has
a convergence certificate. The common target fingerprint is
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`,
the weight fingerprint is
`2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`,
and the selected soft checkpoint SHA-256 is
`b5e21a8584fa536a3740039b3e54480f3c318a63b96266bd31406712c5da991c`.
The revised timing was author-adopted October 3; chain 13 is the working
continuation anchor under the retained target system. Its numerical optimum,
grid adequacy and dated transition remain uncertified.

Each chain had a six-hour budget, at most 250 objective evaluations,
a retained 1,800-second finalization reserve, one CPU and 24 GiB; case lifecycle
limits and checkpoints are defined in the driver plan. The combined maximum
budget was 288 core-hours, with at most 48 concurrent single-core chains.
No failed chain was silently retried, and the selected points remain experimental.

The input contract pins target SHA-256
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`
and weight SHA-256
`2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`.
The v2 upload archive SHA-256 is
`f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b`;
the manifest verifies the pinned source set and enumerates the driver repair.
The v3 expansion archive SHA-256 is
`604096fe0bc9d8d7f7f3de7b52fbf3a7cad411dd4478f4405e4add68f54371a0`;
its [stage and submission evidence](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/expanded/status.json)
pins the parent v2 production/archive and new start table. The parent v2
receipt above identifies the passed two-arm native smoke. The terminal deployment state and
submission receipt supersede older preparation prose/flags.
Use the owning “Compare interest timing and credit” chat and its existing
monitor; do not create a duplicate monitor from this documentation refresh.
Older “prepared/not submitted” prose in design notes is superseded by the
structured deployment receipt.

## Experimental narrower wealth target: completed local comparison

**Verified October 3, 2026, 08:07 New York.** The isolated local post-interest
transaction-timing search ended with ten terminal chains, all ten passing fresh
native selected-point and exact-repeat checks. Chain 2 has the lowest verified
loss **48.170707377609034** under its *experimental* target contract. The
new pooled 2005/2007 PSID aggregate net-wealth/annual gross-labor-earnings
target is **4.45838713455674**, excluding business/farm equity, other real
estate and vehicles while retaining catch-all other assets; the model moment
at chain 2 is **6.561732864006831**. The old numerical weight
**7.595098472533724** was retained for a controlled sensitivity, and no new
standard error has been estimated. Entry wealth and income distributions, the
bequest target and all other targets, model observers and economic inputs remain
unchanged relative to the alternative-timing arm. The model bequest denominator
and the new PSID wealth numerator have not been reconciled.

The [complete new-target collection](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/collection/RESULTS.md)
includes the authoritative 14-row fit and 31 parameter records; its native
`target_fit.csv` uses the **old** wealth target only as a diagnostic, while
`target_fit_new_contract.csv` is authoritative for this experimental arm. The
[three-arm comparison](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/comparison/COMPARISON.md)
links every target, model value, gap, weight, contribution, parameter bound and
source hash. Arithmetic rescoring of the same saved moments yields old-contract
scores **18.445305**, **13.771131**, and **15.580542**, and new-contract scores
**51.765179**, **53.064444**, and **48.170707**, for original timing,
alternative timing and new-wealth chain 2 respectively. These scores are not
new solves. The age-25 children-ever-born stock is **0.535240** at chain 2
against **0.809528**, essentially unresolved across the three points. The ten
positive-weight moments equal the ten free coordinates in count, but informative
rank and optimizer convergence are uncertified. The unequal 48-chain and
10-chain budgets do not establish that the new wealth target is unreachable.
Neither the target nor the chain-2 parameter vector is adopted; the working
continuation anchor remains the old-target post-interest chain 13. A separate
ten-start old-target search under this timing, with
\(\beta_{\mathrm{annual}}\in[0.94,0.99]\), was submitted on October 3 as
Torch array **19112020** (`0-9%10`) under
`/scratch/td2248/projects/soft_timing_continuation_20261003_v1`.
Its ten starts are the ten lowest-loss distinct numerically verified
post-interest endpoints from the completed comparison, beginning with chain 13.
It retains all original targets, weights, economic inputs, and ten parameter
bounds. Each chain has one core, 24 GiB, six hours, at most 500 objective calls,
and 1,800 seconds reserved for fresh native verification. The
[start plan](output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/start_plan.json),
[stage and submission receipt](output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/deployment/status.json),
and [chain-0 smoke receipt](output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/smoke_receipt/gate.json)
pin the source, objective and start identities. The exact two-evaluation smoke
passed the 14-row target fit, 31-parameter record, 17 standard plots and
selected-point repeat. The production search had no fallback or automatic
retry. At Tommaso's request, all ten tasks were stopped on October 3 at
approximately 17:00 New York time, before final optimizer completion or native
verification. The per-chain best/latest/start checkpoints saved before the
stop are preserved in the [old-target pause receipt](output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/pause_20261003.json).
They contain no resumable optimizer state; restarting would be a new search.
There is no final verification or new adoption from this array.

The isolated **new-wealth-target Torch continuation** was submitted as
array **19111687**, tasks `0-9%10`, in
`/scratch/td2248/projects/alternative_wealth_calibration_20261003_v1`.
It searches ten distinct starts with
\(\beta_{\mathrm{annual}}\in[0.93,0.99]\), retaining the new wealth target
4.45838713455674, its numerical weight 7.595098472533724, post-interest
timing, and all other economic inputs. Each task has one core, 24 GiB, six
hours, a 500-objective-call cap and 1,800 seconds reserved for final native
verification. There is no fallback, automatic retry, or adoption. The pinned
[source archive](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/deployment/stage.tar.gz)
and [production receipt](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/deployment/production_submission_receipt.json)
identify the submitted source and complete target/weight fingerprints. A
two-call Torch smoke at the new lower bound, job **19111315** task 9, passed
the exact full-equilibrium repeat, new-contract fast/full comparison, 14-row
target table, 31-row parameter table and 17 standard plots; its
[collection receipt](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/smoke_collection/smoke_collection.json),
[new-contract target table](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/deployment/smoke_target_fit.csv),
and [parameter table](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/deployment/smoke_parameters.csv)
is diagnostic, not a calibrated result. At Tommaso's request, all ten
production tasks were stopped on October 3 at approximately 17:00 New York
time, before final optimizer completion or native verification. The per-chain
best/latest/start checkpoints saved before the stop are preserved in the
[new-wealth pause receipt](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/pause_20261003.json).
They contain no resumable optimizer state; restarting would be a new search.
There is no final verification or new adoption from this array.

## Interactive inspection and numerical readiness

**Verified October 3, 2026.** The October 3 stationary-GE deployment
commit is `94c0a6e3`; it is distinct from the chain-13 scientific
anchor and from either live search array below. The editable [`run_model.py`](code/model/run_model.py)
now runs the adopted post-interest soft chain-13 snapshot as a stationary GE
under `fixed_h0`. Its validated saved case is [`output/model/local_solution/latest`](output/model/local_solution/latest/SUMMARY.md),
with a complete 14-row target table, 31-row parameter table, 17 standard plots,
8 policy plots and 7 aggregate plots. The deployment's fresh same-host reference
checks matched 91 arrays, the full target table, numeric parameter estimates and
bounds, and all 17 standard-plot hashes for both unchanged inputs and a beta
lowered by 0.001; only 13 documented descriptive role/status fields differed.
See the [deployment report](output/model/production_deployment_20261003/README.md)
and [source map](code/model/README.md). This verifies the local deployment in
that scope; it establishes no new calibration, grid convergence or dated transition.

The runner now selects an external input file with `PARAMETER_FILE`:
[`best_params.py`](code/model/parameters/best_params.py) contains that same working
anchor, while [`toy_params.py`](code/model/parameters/toy_params.py) is an independent
editable example with annual beta lower by 0.001. Toy outputs have their own
`output/model/experiments/toy_params/latest` pointer. Plotters and the explorer
follow the selected file, with `--params` available for a one-off selection.
The [parameter guide](code/model/parameters/README.md) records the controls and
output routing. These files isolate parameter experiments and their outputs;
they do not isolate edits to shared solver source code.

The parameter-file workflow passed 40 focused checks and a fresh toy stationary
solve, with all 32 figures and the full tables saved separately; production
`latest` and its solution hash stayed unchanged. Both cached plotters and the
toy explorer were checked without additional solves
([verification receipt](output/model/production_deployment_20261003/parameter_file_workflow.json)).
The canonical calibration driver also passed initialization with zero solves.
It now makes the authenticated zero unsecured-credit binding explicit and
exports run-local `best_params.py` only after final native verification, exact
repeat and matching full input identity. No optimization was launched for this
workflow change, and existing cluster deployments were not rewritten. Older
receipts without this input fingerprint are not eligible for automatic export.

Cleanup retained 24 frozen source-authentication dependencies: the two archived
packages have original-path compatibility symlinks, and the historical runner
and five helpers retain their exact pinned bytes. These are not alternative
production entry points; see the [archive record](calibration_archive/model_legacy_20261003/README.md).

The runner holds the physical housing coefficient fixed while finding the
birth-renewal price root. It reports implied H0 and population scale N
algebraically from the same solved policies and price; this adds no solve and
depends on conditional scale independence under the fixed-payroll mapping.
The `population_one` calibration normalization instead fixes household scale at
one and derives H0. The production guide documents this distinction and the
ordinary input controls. Supported primitive edits through `model.P` are carried
into the next solve; entry-distribution and structural grid edits are explicitly
unsupported. The active grid has 120 wealth nodes and 9 income states. Grid
adequacy and optimizer convergence remain uncertified.

For a quick fixed-price household inspection, use
[`MODEL_PLAYGROUND.md`](code/model/tools/MODEL_PLAYGROUND.md) and launch
[`start_model_playground.command`](code/model/tools/start_model_playground.command).
Its `model.params` interface edits the ten displayed calibration coordinates;
ordinary supported primitive edits are available through `model.P`. A call to
`model.solve()` solves the lifecycle at a fixed price. It does not solve the
stationary GE root, clear housing markets, or fit the targets. The cached
plotters ([policies](code/model/plot_model_policies.py),
[aggregates](code/model/plot_model_aggregates.py)) read `latest` without solving.
For the saved-case browser, launch
[`start_model_explorer.command`](code/model/tools/start_model_explorer.command)
and use the local URL and case it reports; the server reads saved artifacts and
does not solve the model. Do not open the HTML directly or assume a remembered
port or selected case is still current.

The adopted working continuation anchor remains post-interest chain 13 under
the old wealth target (6.926584), with loss 13.771131. Old-target array
**19112020** and experimental new-wealth-target array **19111687** were both
stopped at the author's request on October 3 at approximately 17:00 New York
time. A follow-up scheduler check at approximately 17:01 New York time found
the queue empty and all 20 tasks `CANCELLED`. The pause receipts preserve each chain's best/latest/start
checkpoints: [old target](output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/pause_20261003.json)
and [new wealth target](output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/pause_20261003.json).
Neither array reached final native verification, and neither produced an
adopted result. Their optimizer states cannot be resumed exactly; no restart
or replacement job has been launched. The monitor automation
`finish-and-monitor-soft-timing-calibration` is paused. Do not rank their
losses across target systems.

The [asset-grid diagnosis](output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md)
is historical fixed-price evidence under unchanged prices and entry rules. Its
tail-extension and occupied-grid refinements do not certify the current
calibration, target fit, stationary equilibrium, dated transition, or a
coarse-to-fine production handoff. Historical frontend details remain in the
[archived exact model-README copy](calibration_archive/model_frontend_20261003/README_code_model_before_map.md).
The earlier [overnight-default replay receipt](output/model/fixed_reference_economics_20260928/model_control_scripts_v1/overnight_default.json)
and [workflow verification receipt](output/model/fixed_reference_economics_20260928/model_control_scripts_v1/verification.json)
are historical checks, not current runner or runtime claims.

## Historical policy and transition results to retain

**Hard/quarter purchase rules.** All 24 fresh-calibration selected postchecks
completed; best hard loss is 88.588403, quarter loss 48.319938. The hard housing
floor reaches its upper bound 2.6; the quarter floor is 2.572098. Full fit,
bounds and diagnostics are in the fresh-calibration readout linked above.
These are accepted selected solutions, not proof of optimizer convergence or
author adoption.

The earlier 97.011220/51.556036 points produced accepted temporary tax-shock
paths at 48/64 periods. First-period birth-flow changes were approximately
−0.542%/−0.539% for hard and −0.431%/−0.429% for quarter, respectively.
These are first-four-year birth-flow changes at those historical calibrations.
Permanent terminal steady states passed, implying population-scale changes
−2.324% hard and −2.105% quarter. All four permanent dated transitions failed
root or terminal gates, including the recovered hard-64 run. There is no
accepted permanent transition impact. See
[mechanism deployment](output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism_deployment/README.md).
The older quarter fixed-price birth increase, approximately +0.184%, and its
negative equilibrium response support a rent-offset mechanism at that point;
they are not a demonstrated mechanism for the current soft solution.

**Historical successive-surprise contract.** The authorized historical design
fits four successive surprises, each believed permanent until the next arrives,
or one permanent 2007 shock to the final 2020–2023 window. Carry the complete
household state and both entry queues between surprises; do not substitute an
announced four-shock path. The block0506 attempts 18801439/18801451 failed before
fitting shocks. A later pension-correction diagnostic is not an estimated history.
The [historical estimator packet](output/model/fixed_reference_transition_20260928/four_shock_v1/README.md)
and [bounded diagnostic](output/model/fixed_reference_transition_20260928/four_shock_v1/budget_diagnostic_v2/README.md)
preserve the request and launch history; their old draft-plan disabled flags do
not negate subsequent author launch authorization. No accepted four-surprise fit
is established by the reviewed receipts.

**Normalized-v1 one-shock transition.** Historical job 18995772 hit its time
cap after two trials. Its best shock benefit was \(\psi=0.129489443\), within
bounds \([0.00171989,0.34397799]\); the 2020–2023 model statistic 1.703658
missed target 1.645750 by 0.057908. Four-window target fit and all shock
parameters are linked in the transition readout above. Exact replay and
intermediate housing roots passed, but strict terminal-state/horizon acceptance
failed and final diagnostic plots were not produced. Monitors were paused;
there was no approved automatic budget extension. It is not a certified policy
initializer for the current soft point.

**Paper-facing artifacts.** The [Corina source map](latex/corina_progress_20260930/source_map.md)
and [slides](latex/corina_progress_20260930/corina_progress.tex) use the earlier
hard/quarter calibration and policy packets. They do not display the latest
88.588/48.320 results or the revised-timing soft 13.771 anchor. The coordinated draft,
slides and mock sources are listed in [latex/README.md](latex/README.md).
Experimental physical-floor/nonlinear-benefit objects and the newly adopted
timing are not fully synchronized with the author manuscript or continuing
slides. The specific slide timing discrepancy is recorded in
[latex/README.md](latex/README.md). This status update does not change
manuscript wording or certify a new paper calibration.

## Outstanding decisions and checks

1. **Timing and purchase constraint:** propagate the adopted post-interest
   transaction budget consistently through active implementation, saved-input
   defaults and paper representations before claiming full synchronization.
   The author authorized resuming the current one-birth Estate-A fertility-shock
   estimation on October 3 and clarified on October 4 that different shock
   guesses must run simultaneously, with at most 48 nodes. The current design
   uses twelve independently refining scalar-fit starts, each with the same
   saved baseline, target contract and numerical gates. Local checks passed;
   Torch smoke 19139732 is running and array release is tracked in
   [the current transition readout](output/model/transition_readiness_v1/current_baseline_20261003/README.md).
   This authorization does not certify a transition result. The hard-versus-soft
   economic choice remains distinct.
2. **Identification and weights:** count is ten moments for ten free parameters.
   Check informative rank/substitution and the two near-bound fertility/continuation
   noise parameters. National ACS uncertainty and the early-fertility weight remain
   working choices requiring their own decision.
3. **Empirical observers:** resolve or explicitly maintain female versus household
   exposure, Sun–Abraham versus stationary first-birth observation, recent-parent
   residence/age definitions, and old-wealth income/age support. Do not call current
   proxies exact measurement equivalence.
4. **Entry wealth, estates and creditors:** current entry mapping is provisional.
   A retained gate reports positive entrant funding 0.011515 per period against
   provisional available net estates 0.120274, but this is not full donor-recipient,
   spouse, negative-estate/creditor or physical-resource settlement. The SCF
   calculation being complete does not close these mappings.
5. **Demographic and population interpretation:** explicitly state household,
   reproductive-member and resident-person units; retain the implemented queue and
   three-plus adjustment. Do not infer literal genealogical tracking, a factor-two
   correction, total fertility effects or a geographic migration closure.
6. **Numerical readiness:** complete relevant grid/target/GE convergence checks
   before production policy claims. Fixed-price grid tests and refactor parity
   do not certify a dated transition or coarse-to-fine handoff.
7. **Policy acceptance:** historical permanent paths and the one-shock normalized-v1
   path failed required gates. A current soft-policy result needs its own reconciled
   closure, reference bridge and complete visual packet.
8. **Publication contract:** distinguish an exploratory working calibration from
   author-adopted economics, preserve the frozen September 14 reference, and
   reconcile accepted changes with the authorized paper representations.

## How to maintain this note

Replace superseded state within its section. Record the source identity and
verification time for mutable claims; do not paste another chronological
transcript here. Put detailed chronology in daily notes or the experiment README,
durable preferences/gotchas in `memory/AGENT_MEMORY.md`, and old snapshots in
`calibration_archive/`. Advisory size budgets are a review signal, not permission
to truncate unresolved items, discard evidence or silently adopt experiments.
