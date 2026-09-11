# Overnight proposal: utility, pensions, calibration and policies

**Status: Fable 5.1 review completed; revised by the lead; awaiting Tommaso's approval. No model implementation or cluster job has been launched in this preparation.** The author requests a broad review first, then approval and launch. Maximum proposed computation window: 12 hours from approved launch; a morning readout is delivered even if the bounded computation window extends beyond it. The prior September 9–10 execution record is retained below as history, not a live queue.


### What the author is being asked to approve

Approve implementation and a **gated 12-hour run**, using up to **96 concurrent one-CPU jobs** across nodes. Retain the agreed parenthood-only utility, fixed 17.9% tax with actual-budget pensions, the observed-age 2007 initialization, and all proposed housing identification moments. Approve the **new early-target design described below**, conditional on its builders, observers, weights and population identities being certified; it does not silently inherit the old late-target fingerprint. The initial model fertility normalization remains 2.1, with its empirical approximation disclosed.

The core morning deliverable is the verified utility/fiscal implementation, balanced initial/terminal accounts where solved, the strongest completed initial fit or explicitly provisional diagnostics, and measured runtimes. A fully fitted historical path is the next priority. Certified policies are contingent on that success. An optional short policy screening can be produced from a numerically valid finite-horizon diagnostic baseline and clearly labeled as such; it is not a production policy result. The primary policy remains the equal-rebate property-tax comparison.

Use a single supply elasticity of 0.63 in the new initial and dated model. At the mapped seed, rebase its level to preserve the old price/quantity point; thereafter estimate the scale and keep the selected supply law fixed across baseline and policies. This is an explicit proposed closure revision, not an assertion that changing elasticity leaves every recalibrated equilibrium unchanged.

Fable agrees on the utility algebra and pension repair and warns correctly that parallel core-hours do not remove the sequential critical path. Its review is **plan-only**, not source verification. The lead rejected unsupported recommendations to drop the family-size rooms target, remove the observed-age initialization, declare a stationary endpoint impossible, treat unchanged continuous controls as proof of a cache bug, or certify 28 dates from a household's finite life alone. Full adjudication and unedited review: `FABLE_OVERNIGHT_REVIEW.md`.

No guarantee is made that calibration, horizon verification and policies will all finish by morning. If a prerequisite fails, the remaining budget goes to the smallest decisive repair/diagnostic and the failure is reported. Nothing here authorizes weaker gates or an unreported change in economic assumptions.

## 1. Objective and the decisions already made

Deliver the strongest verified quantitative benchmark feasible in the approved window: the simplified utility, balanced Social Security, a defensible initial calibration and an announced historical fertility-preference path. Conditional on that benchmark passing, solve matched housing policies. A complete successful calibration and policy package by morning is an objective, not a guaranteed outcome.

The experiment is about housing-policy scenarios conditional on a fertility-preference decline, not identifying the cause of the U.S. fertility decline. Approximate the pre-announcement economy by a stationary benchmark with the author-selected fertility normalization 2.1. Households learn the full preference path in 2007. Fit historical observations over 2007–2023; hold the preference intercept constant at its 2023 value thereafter. Population, prices, pensions and fertility remain endogenous after 2023. Keep the sequential fertility architecture; this request does not reopen the choice-shock nests.

Use the isolated checkout:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf`, branch `codex/balanced-social-security`, current repair commit `a654219c`.
Main checkout contains unrelated edits. `latex/JMP_DS_draft/` remains strictly read-only. No manuscript or slide rewrite is part of this run.

The exact user utility note is retained at:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/utility_change_author_specification.md`.

## 2. Utility specification and what must be tested

Let c be nonhousing consumption, s housing services, and m the number of children currently at home, not lifetime parity. Maintain

\[
e(m)=((2+0.7m)/2)^{0.7},\quad Q=c^{\alpha_0}[s-h_P\mathbf1\{m>0\}]^{1-\alpha_0},\quad
u=-e(m)/Q+\psi_t m,\quad\sigma=2.
\]

Maintain the fixed scale coefficients, constant alpha, zero nonhousing floor, linear direct child benefit and every other lifecycle/fertility component. Do not add a household multiplier, a new curvature parameter or a child reward merely to offset the scale's concavity. The author explicitly accepts economies of scale among additional children. Housing and the scale are coherent distinct cost components; no double-counting theorem has been established.

Replace two housing-floor coordinates with one. At the unchanged inherited parameter point set h_P=hbar_first_child_jump+hbar_child_rooms, preserving the first-child requirement; set the old slope to exactly zero and exclude it from the searched vector. The retained code already permits this floor shape. Prefer a clear parameter adapter and contract over a new household solver. Every launch/reload must enforce the zero slope; old calibration bounds must not silently reintroduce it. Provisional h_P search bounds are [0.1,2.3] rooms, the sum of the two inherited bound intervals, subject to the existing feasibility checks. These are proposed numerical bounds, not empirical restrictions; review before approval.

Verification must check the actual solver's m=0,1,2,3 type arrays, effective CRRA multiplier, parenthood floor disappearing when all children leave, and annual-to-four-year beta conversion. Include psi=0 as well as nonzero psi so the equivalence scale cannot be lost when parent types have identical housing floors and no direct child reward. The inherited beta is annual 0.9952765791 and enters as beta^4; do not change it outside calibration. For an interior renter at fixed expenditure X and rent r, verify s*=(1-alpha)X/r+alpha*h_P for every m>0. Test constrained renters and owners separately; that formula is not their universal policy rule. Test that the first-child floor is preserved in the initial mapping, without requiring equilibrium moments to remain unchanged.

A four-cell conditional comparison separates sources of change: (A) old utility/old pensions; (B) old utility/balanced pensions; (C) new utility/old pensions; (D) new utility/balanced pensions. A and C are explicitly invalid-fiscal diagnostics, not selectable benchmarks. Run a frozen-psi comparison first to isolate the direct utility/pension changes; report any separately renormalized 2.1 comparison as such. Do not allow renormalizing psi to conceal the first comparison.

## 3. Social Security: restore the actual budget condition

Tax tau=0.179 remains externally fixed. Pension b is determined internally:

\[
\tau B_t=b_t E_t,\quad
B_t=\Delta\sum_{i,j<J_R,z}w_i a_j z\,m_{tijz},\quad
E_t=\sum_{i,j\ge J_R,z}[1+s_z(z-1)]m_{tijz}.
\]

Here m_t is actual household-head mass summed over wealth, housing, parity and dependent-child states; Delta is the model flow-period scale, and s_z is the inherited retirement-income dispersion coefficient (zero in the common-benefit case). Resident nonheads and property-tax rebates do not enter this budget. Do not multiply period pensions by four twice.

Rebuild the pre-announcement stationary endpoint, the reweighted announcement-state accounts and the terminal endpoint. Use the same dated pension path in household backward expectations, forward choices and fiscal accounting. Housing prices and fiscal consistency must both clear. The previously implemented generic joint root still needs explicit production adapters and actual endpoints; do not call it a completed fiscal repair.

Pure tests pass, but the actual compiled smoke is not certified: its universal occupied-saving-response assertions failed. A follow-up showed future-income effects on values and conditional controls at unoccupied states. Replace the invalid universal assertion with a meaningful predeclared anticipation test, preserve the original failure, and rerun the complete compiled loop; do not change the household optimizer just to force a response. The independent evidence is in `social_security_repair/README.md` and `anticipation_diagnostic.json`.

Existing housing residual gate remains 2e-4. Proposed additional fiscal gate: abs(revenue-outlays)/max(revenue,outlays,1e-12) <= 1e-6 for each positive-budget state/date, with absolute residual also saved. The zero-budget case requires zero outlays and a defined handling rule. Keep household budgets, feasibility, probabilities, mass/entry queues, population identities and fresh replay checks. Document all fiscal bounds as numerical bounds and flag binding bounds. Existing terminal/horizon gates are not waived.

## 4. Calibration design and empirical prerequisites

First calibrate common parameters in the approximate pre-2007 stationary economy; then freeze them while estimating the historical preference path. This differs from the old joint late-moment calibration and requires approval of a named new target contract. The proposed early data and original builders already exist. Their observers, remaining uncertainty and weights are not all certified. We must finish these before calling any early search calibrated SMM.

With jump-only housing, the initial searched structural vector has **nine** coordinates: beta_annual, kappa_fert, kappa_fert_continuation, chi, H0, theta0, theta1, first_birth_fixed_cost, h_P. The initial psi_0 is separately normalized to 2.1; the historical preference parameters are fitted later. All nine structural coordinates vary jointly in the main search. The equivalence scale, alpha, sigma, tenure dispersion, income process and survival remain external. The new proposal uses housing-supply elasticity 0.63 in both initial and dated economies, replacing the inherited initial 1.75 explicitly. Rebase the seed supply level so the old initial price/quantity point lies on the new curve, then re-estimate the scale; hold the selected initial supply function fixed through history and policy. Verify the code's supply units before calculating the rebase. This proposed revision requires the author's approval with this plan.

The proposed initial target system retains 13 restrictions (including the 2.1 normalization), hence 12 scored moments for nine searched structural parameters. Count is not proof of identification. Removing the housing slope does **not** remove the family-size rooms moment.

| Restriction or moment | Early source/value | Status for approval |
|---|---|---|
| Initial model fertility | Author normalization 2.1 | Keep; explicitly distinguish its model units from empirical female period TFR and old completed cohort fertility |
| Childlessness at 40–44 | CPS 2004/06, 0.198279 | Point available; uncertainty and exact model observer to verify |
| Mean age at first birth | NCHS 2003–06, 25.976264 | Period-flow observer and uncertainty/scale to finish |
| First births at 30+ | NCHS same period, 0.249278 | Same issue |
| Exactly one child among mothers 40–44 | CPS 2004/06, 0.213655 | Point available; uncertainty and observer to verify |
| Wealth / gross earnings | PSID 2003/05, 6.145861, SE 0.362855 | Verified early builder/uncertainty; confirm model denominator |
| Annual bequests / wealth | Inherited restriction 0.0088 | Retain with external status and inherited documented weight |
| Old wealth/income p90/median | PSID 2003/05, 3.515935, SE 0.306911 | Preserve ages 76–84 and exact ratio definition |
| Mean rooms | ACS 2005/06, 5.561097, SE 0.088381 | min(rooms,9) observer; same 42 metro codes/new footprint |
| Ownership at 30–55 | ACS 2005/06, 0.648334, SE 0.020675 | Match household-head sample and structures |
| First-birth room response | Reviewed PSID Sun–Abraham, 0.720246, SE 0.085260 | Retain exact source/contrast; early application assumes stability |
| Rooms: 3+ versus 1–2 resident children | ACS 2005/06, 0.347067, SE 0.059705 | Retain; reconcile observed children and model dependents |
| Recent-parent ownership gap | ACS 2005/06, 0.162896, SE 0.006080 | Retain; exact recent-parent/comparison groups must match |

Young ownership 25–34 (0.431158), old median wealth/income (7.285793) and old completed fertility (1.856608 capped) remain visible validation rows in this proposed contract. The author has not approved treating the old completed fertility observation as literally equal to the 2.1 normalization. The current late-target contract remains preserved and losses are not compared across contracts.

The source report is `design_research/README.md`; candidate table and limitations are in `design_research/DECISION_REPORT.md`. That older report's 2.0605 anchor and two housing parameters are superseded by the author decisions above. Pin each authoritative builder, sample, measurement date, geography, uncertainty, weight and source hash. No substitute empirical scale is invented just to start a search. Use verified inverse-variance diagonal weights where defensible; inherited synthetic scales remain explicitly synthetic. Complete missing CPS/NCHS uncertainty from the original samples/definitions. If that cannot be completed, early diagnostic panels may run under explicitly provisional weights, but cannot become the calibrated benchmark or launch production policies.

The 2.1 normalization also appears in birth-to-entry conversion, initial queues and an entry gate. Preserve or explicitly justify each population identity; never substitute female TFR into a household conversion automatically. The inherited 2007 age reweight is conditional demographic initialization, not an endogenous population fit. Four-year birth dating, female exposure and maternal age/parity mapping must be certified before fitting a fertility-rate path. National fertility/PSID plus metro housing remains an explicit geographic approximation.

## 5. Historical fit, continuation and policy specification

The historical fitting dates are 2007–2023, not 2403. Subsequent dates close expectations. The old 100-date horizon was motivated by slow demographic convergence, not a demonstration that all 100 dates are needed for historical moment accuracy. Historical rate observations need their own female exposure. Do not relabel the existing birth-count decline as TFR. Under the current bridge, a decision at t produces births in t+1,...,t+4; the 2020–2023 block therefore maps to the 2019 decision.

Proposed first shock fit: one amplitude with a linear 2007–2023 decline; evaluate its entire available historical profile, not only one endpoint. In parallel, prepare a parsimonious two-parameter alternative (amplitude plus one timing/curvature parameter), activated only if the linear profile demonstrably misses intermediate windows and enough independent observations identify the extension. No arbitrary annual shock series, stochastic preference process, or search over future post-2023 trends. The complete path is announced in 2007 and flat after 2023 in both cases.

Use 28- and 48-date continuations as screening comparisons, with the 100-date horizon as an inherited reference, not a certified truth. Re-clear each compared path with the same parameter and fiscal contract and compare every scored historical moment and historical price/pension. A proposed horizon-accuracy criterion for reviewer assessment is changes <=0.05 empirical standard errors per scored moment, with explicit absolute tolerances for synthetic/no-SE rows and <=0.1% for historical prices/pensions. This is an additional proposed accuracy test, not permission to waive existing terminal gates. Retain long-horizon verification for selected candidates; if no shorter horizon is certified, label shorter searches approximate and re-evaluate their finalists before promotion.

Do not recompute identical empirical inputs, initial states or derivative panels unnecessarily. Reuse a numerical value only when parameter/input/source hashes match; otherwise reuse it as a starting guess and re-solve. For the initial stationary economy, first exploit the analytic survival/earnings marginal recursion if its exogeneity is verified against the actual distribution; the source uses exogenous entry-income weights and structural age survival, so a separate pension root may be unnecessary there. On the dated path, a damped forward pension update followed by a new backward/forward replay is the first numerical method to test. It is not a proof of convergence or permission to skip the anticipation-consistency residual. Retain the bounded joint solver as a fallback; never assume exogeneity throughout the post-2023 demographic transition. Warm-start housing/pension paths and their numerical derivative matrices, with fresh checks after each changed candidate.

After a selected benchmark passes, the primary policy is the property tax 1% to 2% with equal household rebates **in both baseline and reform**. Pensions still have their separate payroll budget. Compare common 2023 inherited states, use anticipated policy paths consistently, and recompute policy-specific endpoints. Do not use zero starting transfers in the equal-rebate solver to masquerade as a no-rebate case. Supply +20% and dependent-child LTV95% are next priorities; any unrebated tax decomposition is subordinate and must name the disposition of retained revenue. Post-2023 preference alternatives are separately labeled scenarios. No production policy launches from an uncertified calibration/closure.

## 6. Proposed parallel jobs, solve counts and budgets

**Resource ceiling proposed for approval:** up to 96 simultaneous one-CPU jobs across cluster nodes, 16 GB/job for full-grid model paths (about 1.5 TiB if all occupied); use the measured memory requirement for smaller stationary jobs. Account `torch_pr_570_general`. Set numerical libraries to one thread per allocated CPU. Do not request 96 full exclusive nodes or pretend an individual solver is an MPI program. Increase concurrency only within this envelope and cluster/account limits. At preparation, `torch.sh status` succeeds and the user's queue is empty; this is not a reservation of capacity.

Time budgets begin at author-approved launch and are shared, not independently stackable:

| Phase | Concrete work | Budget/parallelism | Progression rule |
|---|---|---|---|
| A, first 0–2 hours | Utility adapter, fiscal endpoint/path adapters, empirical contract; independent scoped reviews; compiled exact-loop smoke | Distinct file ownership; up to 8 small test jobs after approval | Source/math checks and genuine loop outputs pass; new target contract complete before SMM |
| B, diagnostic wave | Four utility/fiscal cells; new-spec stationary normalized baseline; 19-case +/- local panel for 9 coordinates | At most 23 normalized cases plus conditional cells; 19 independent derivative cases | Complete all fits, bounds, budget checks and rank/sensitivity; no target dropped for poor fit |
| C, stationary search | 96 space-filling or bounded local perturbation candidates across all nine coordinates, then up to 6 promising multistart refinements | First wave 96 cases; refinements use batched local derivatives/trial steps, max 48 evaluations/chain and 3 hours/chain; total cap 384 normalized evaluations, excluding B and final repeats | Budget derivative evaluations explicitly; reuse/update local matrices with safeguards; at most 3 finalists, each repeated twice |
| D, historical fitting/horizon work | Prioritize one initial finalist: up to 4 amplitude candidates, then up to 2 adaptive refinement candidates; 3 matched 28/48/100-date horizon checks; optional no-shock demographic control and at most 2 backup/timing cases | Global cap 12 PF roots, including horizon/control roots; each <=8 full mappings including fresh replay; shared 8-hour stage ceiling and remaining global deadline, whichever binds | Begin with 28-date screening; launch the selected long reference as soon as a credible candidate exists. No claim that 28 dates are adequate before horizon checks; unfinished long references keep results provisional |
| E, policies | Rebated-tax pair first, then supply/LTV if time permits | At most 4 paths including baseline, each <=8 mappings; same 96-job ceiling and global deadline | Production only after certified common benchmark. Optional 28-date screens require a valid matched finite-horizon baseline and explicit diagnostic label; no screens with uncertified fiscal/measurement inputs |
| F, final verification/readout | Fresh repeats, complete target/bound tables, stationary and dated fiscal ledgers, stable diagnostic graph packet | Reserve final hour; bounded jobs may finish within the approved 12-hour window | Explicit complete/partial/failure classification, no claims based on queued jobs |

Phases overlap only when their scientific dependencies allow it. Empirical/adapter work can run in parallel. Calibration refinements run on independent starts. Transition/horizon experiments on provisional parameters are labeled diagnostic; final policies wait for the selected calibrated benchmark. Distinct update iterations inside one equilibrium remain sequential.

Observed **old-code**, unbalanced-pension times: about 5 minutes for one normalized initial candidate (four stationary GE calls in the observed case), 41–55 minutes for a 100-date mapping, and about 7.5 hours for the last 9-call price root, after derivative preparation. These are not promises for the revised model. At 5 minutes/candidate, C's 384-evaluation cap is 32 core-hours and the 96-candidate first wave is 8 core-hours. Each 48-evaluation refinement is 4 core-hours but has a 3-hour wall limit, so batching is necessary to use its entire evaluation budget at that timing. A central nine-coordinate Jacobian costs 19 candidate evaluations including its baseline, not one solve; for six unrelated refinement centers that is 114 evaluations, and further trial/rebuild calls consume the same 288-refinement cap. Do not promise a fixed number of nonlinear iterations. B's 19 normalized cases are 1.58 core-hours. If each initial candidate still requires four GE calls, C implies up to 1,536 stationary GE calls; the new fiscal method may change that count. Measure and resize in the exact-loop smoke before full dispatch.

For D's worst-case 100-date horizon, 12 roots x 8 mappings x 100 dates x 2 household passes = 19,200 dated Bellman calls and about 66–88 core-hours at inherited timing, excluding new endpoints, derivative construction and compilation. E's corresponding cap is 4 x 8 x 100 x 2 = 6,400 calls, about 22–29 core-hours. These are worst-case solve caps, not additions automatically scheduled; the 12-hour global deadline can bind far earlier. Eight 100-date mappings take about 5.5–7.3 hours per root even when other roots run in parallel, so launching full-length policies only after a late-night baseline finishes is not a credible morning promise. A fresh 200-column joint price/pension finite-difference matrix could itself cost 200 full mappings; do not launch that by default. First test block structure/warm-started safeguarded updates; a new large panel requires explicit sizing inside the overall approved budget. Do not use simultaneous colored/comb perturbations unless the required Jacobian sparsity or an error-controlled approximation is actually established. Shorter verified horizons reduce these costs, but that reduction is not assumed in the launch promise.

Actual cluster queue, allocation speed, utility effects and fiscal convergence may prevent D/E completion by morning. The plan's fallback is more completed repair/calibration evidence, not invalid policy output. At 90 minutes without a viable fiscal exact loop or certified empirical contract, report the specific barrier and restrict subsequent jobs to clearly labeled diagnostics until it is resolved. Stop launching a stage when its reserved completion/replay cannot fit the remaining budget. No blind restarts or threshold relaxation.

Every case writes source/target/parameter contracts, progress at least every case or five minutes, latest completed and best-so-far records, full target fit and all parameter bounds, and explicit invalid-case reasons. Investigate 30 minutes without progress. Checkpoint finished cases and collect them while later cases run. Monitor only meaningful changes/failures/completion. This preparation does not reactivate the old paused monitor; activate the approved run's monitor only after launch. Cluster jobs survive a laptop disconnect, but local monitoring does not acquire that property automatically.

## 7. Review questions and approval deliverable

Fable should assess the complete sequence, not assume the stage budgets guarantee feasibility:

1. Does the new utility match the author's note exactly, and which family-size/housing moments now test its strongest restrictions?
2. Is the fixed-tax pension rule carried through initial normalization, observed-age reweighting, terminal demographics and both PF passes? What is the smallest decisive missing test?
3. Is the proposed initial target system coherent and sufficiently identified? Which unresolved observer/weight/demographic issue must block SMM versus remain a declared approximation?
4. Is two-stage calibration plus a parsimonious historical shock fit defensible for this policy question? What is the sharpest achievable alternative if the early contract cannot be certified tonight?
5. Are the solve counts, parallelism, dependency graph, horizon strategy and stop rules credible? Identify wasted computation and changes that materially shorten time without changing economics.
6. What should be cut or reordered for the morning deadline? Return a concrete revised plan, major objections and exact approval choices; do not launch anything or edit model/empirical/paper files.

The review has returned and is retained in `FABLE_OVERNIGHT_REVIEW.md`, together with the lead's adjudication. The executed review used only a self-contained plan/utility message, with all built-in tools and MCP access disabled; Fable did not read these evidence files. Lead source checks, rather than reviewer prose, support the corrections recorded there. No new numerical model result was produced by review.

## Evidence for the reviewer

All relative paths below are under `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

- Canonical live state: `CALIBRATION_STATUS.md` (top September 11 fiscal/history corrections supersede older blocks).
- Fiscal implementation/evidence: `output/model/e5f_matched_pf_20260909a/social_security_repair/README.md`, `budget_audit.json`, `historical_reconciliation.json`, `anticipation_diagnostic.json`.
- Actual isolated sources: `tmp/e5f_matched_pf/code/model/tools/e5f_social_security.py`, `e5f_social_security_root.py`, `test_e5f_social_security_compiled.py`, `e5f_matched_pf_initial_state.py`, `run_e5f_matched_pf_baseline.py`, `run_e5f_matched_pf_history.py`; package `intergen_eqscale_seq_optimized/parameters.py` and `solver.py` type construction around 2220–2305.
- Early data, every proposed target and caveat: `output/model/e5f_matched_pf_20260909a/design_research/README.md`, `DECISION_REPORT.md`, `computation.md`.
- Complete inherited target-fit table: `output/model/e5f_matched_pf_20260909a/design_research/computation/final_replay/evaluation_003/target_fit.csv`; all inherited parameters/bounds: adjacent `parameters.csv`. These are pre-repair inherited estimates, not new calibration results.
- Historical timing evidence: `output/model/e5f_matched_pf_20260909a/path_pilot_20260910/README.md` and `RESULTS.md`; old 100-date prices clear but fiscal and horizon certification fail.

---
# Historical record: September 9–10 (superseded launch instructions)

# Overnight work: September 9–10

Author instruction: work across the quantitative dimensions overnight and provide a rigorous morning readout. Target delivery08:30 EDT September10; no guarantee of a completed calibration or journal readiness. The previous deleted generic routine stays deleted; the existing paused overnight heartbeat has been updated to this bounded assignment.

## What is being fitted

The current observer measures its12-row vector at2023, with reconstructed cohort histories and a2019–2023 birth/control branch. Initial2007 conditions approximate a preceding stationary regime normalized to completed fertility2.1. Households learn the full transition immediately in2007. This is not a fit to an annual historical time series.2023 is an observation along a PF path, not the terminal SS. Data combine2024 CPS completed fertility,1979–84 NCHS cohorts, reviewed PSID event-time evidence, four pooled2012–2023 ACS observations, approved older PSID wealth vintages and an external bequest normalization.

Read overnight_target_mapping_review.md for all12 date/operator mappings. Four ACS pooled-date comparisons and two family-group mismatches remain outstanding. Preserve the current objective; any diagnostic measurement change needs a separate operator identifier and stays outside that loss. No empirical row is silently replaced, dropped or reweighted.

## Current numerical evidence

Anchor17300115 completed47m02; model time2805.147s,200Bellman calls,100four-year dates2007–2403, final2407.55cluster startup tests pass. Collector load_case independently accepts all mapping gates and recomputes market residuals from the CSV. Four local artifact hashes match; large state/checkpoint binaries remain remote. Maximum market residual0.1819903189: not equilibrium. Terminal person/head gaps0.00822214/0.00806724 pass1%; normalized household distribution L1.0501169 still exceeds.02. No horizon or calibration certification. Off-equilibrium loss is not a horizon-improvement comparison; complete tables are in meeting_receipts/horizon100_anchor_01/sequential/.

## Submitted overnight chain

Immutable source E: /scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a.

-17303124: first/last-coordinate price probes0,99. Two independent one-core16GB cases,200Bellman calls each,60min internal and65minSlurm limits.
-17303125: remaining98coordinates, dependent afterok17303124, maximum32simultaneous one-core16GB cases.
-17303126: collector afterok both arrays, five-minute limit. Requires all100distinct probes, identical scientific/target contracts, valid source pins and mapping gates; saves output/jacobian_horizon100_sequential.json.

All100contracts differ from the successful anchor only by probe_coordinate. Each changes one logprice by.01 with parameters/targets/weights/gates fixed. Full batch20,000Bellman calls, about77.9CPU hours from measured runtime. With32slots plus the two-case smoke, roughly4hours before collection excluding queue delays; actual memory and priority may lengthen this. Existing28date exact perturbation loop and full100date anchor passed; fresh first/last probes gate the enlarged panel. Each case retains heartbeat/latest_completed/best_so_far, full fit and final state. Failures preserve their outputs; no blanket retry. Investigate any running case without heartbeat for30minutes.

## Active independent reviews

Existing in-thread agents, separate20minute bounds, no source edits or cluster launches:
- policy_compatibility: overnight_target_mapping_review.md (written; lead checked observer chronology and authoritative target ledger).
- pf_smoke: overnight_numerical_review.md (readiness and efficient full-root/horizon plan).
- tax_transition_contract: overnight_policy_readiness.md (paired PF tax/rebate closure and exact adaptation route).

Lead owns economic judgment, objective, provenance and final verification. Finish reviewing pending files before applying their recommendations.

## Next decisions and execution

1. Verify smoke, full price panel and collector; inspect conditioning and signs. Do not call a probe an equilibrium.
2. Prepare and test an isolated runtime allowance for100date roots: current root cap7200seconds is too small for six evaluations at47minutes each. Use a new immutable source snapshot and honest source-map reconciliation; never silently loosen provenance. A six-evaluation root needs about4h42plus overhead. Consider a smaller bounded round with fresh initial/final and verified continuation, but do not reserve too few calls to take an actual step. Every market and reproduction gate remains unchanged.
3. Solve full PF markets, then evaluate terminal distance AND historical price/moment stability under an appropriate extension. The full household distribution still fails at the supplied anchor; market clearing may change that. No automatic claim that100dates is sufficient.
4. In parallel, prepare date/group diagnostic aggregation and local identification calculations without altering the active objective. Reuse required equilibrium replays for supplemental measurements where source provenance can be explicitly maintained. Exact recent-parent/adult-child residence cannot be inferred from count-only states.
5. With a valid stable equilibrium evaluator, re-normalize all affected initial/terminal objects for each parameter probe, estimate the equilibrium-adjusted moment Jacobian and assess rank/tradeoffs before a bounded joint parameter step. Current inherited11coordinates are not re-estimates. Existing target count12does not by itself establish informative identification.
6. Advance matched nested PF comparison and policy adapters where independent progress is useful. R1/R2 equal-rebate tax comparisons must start from the SAME arm-specific pre-choice2023 historical state. Old person-policy CLI must not silently rebuild a temporary-equilibrium state, re-anchor supply or change headship. Rebuild tax/rebate-specific stationary endpoints. Preserve original verified PF policy as fallback evidence; do not attach it to a new calibration.

## Morning deliverable and persistence

One advisor readout with empirical/model date map, scientific specification, complete target/parameter tables, numerical/horizon qualification, identification diagnosis, actual policy findings and decisions needing the author. A consolidated PDF can use the earlier explicit request when useful; no decorative illustrations or unrequested changes to the standard graph set. No unsupported significance or causal/policy claims.

Heartbeat overnight-rebated-tax-results-and-graphs is ACTIVE every15minutes in this task until08:30 EDT; stay quiet on unchanged state, advance useful authorized work on each wakeup, and pause after morning readout. Do not count a paused/queued job as scientific progress. Mac detached caffeinate PID17077 asserts idle/system sleep prevention for32400seconds; both assertions were verified. Cluster jobs continue independently of the laptop. Protected author draft and unrelated dirty files stay untouched; back up only owned work.


## 00:09 EDT execution update: full root queued

The first and last price probes17303124 are running with fresh15-second heartbeats; the remaining panel and collector wait on success. No new completed equilibrium is available.

Root17304265 is now queued after collector17303126. Snapshot F is `/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a`, isolated source commit96a41873, pushed. The driver and shared contract loader now permit an explicit six-hour ceiling; the actual root budget is21000seconds, Slurm21600seconds, at most six complete paths including fresh initial/final evaluation. Expected compute is about4h41m at the observed2805seconds/path, excluding queue delays; it may still be running at the morning cutoff. Early convergence saves unused evaluations. Market, fiscal/accounting, feasibility and replay gates are unchanged.

The three-file diff changes only runtime bounds and explicit source-provenance reconciliation, plus regression tests. The prior panel snapshot's root and history hashes were checked unchanged after atomic staging. Source changes require exact old/new pins and a declared restricted scope; altered solver or other economic files still fail. The shared history function's backward/forward evaluation is untouched.49pure local tests and58tests on supported cluster Python pass, including accepted long budgets, rejected excess budgets and rejected missing/forged/broadened source reviews. The launcher repeats58tests, verifies source/input hashes, reconstructs the Jacobian from all100receipts and only then finalizes and hashes its root contract. Launch/provenance: submit_horizon100_root.sh and horizon100_root_source_review.json.

The two remaining review agents did not produce their promised files within the20-minute window and were still pending_init at the next wakeup. They were stopped, rather than silently extended. Their preliminary messages are advice, not completed independent reviews. Lead independently verified the runtime finding against both loaders and the actual58-test suite.

Lead also checked the policy-interface finding in run_e5f_perfect_foresight_person_demography_policy.py: its legacy main reconstructs2023, rebuilds demographic primitives and reanchors supply to a rebated-tax baseline. Those operations do not implement the present matched historical contract. The callable solve_person_funded_path already accepts explicit initial_state, primitives, supply_rule and terminal; a matched adapter should use those inputs from the saved historical state. Baseline endpoints presently represent1%annual tax and zero transfers. Rebate and higher-tax endpoints must be rebuilt. For a policy announced in2023, all branches within an arm must start from that same pre-choice historical state.

Let N1,N2 be outcomes under1%/2%tax without rebates, and R1,R2 the corresponding outcomes with equal household rebates. The exact level identity R2−R1=(N2−N1)+[(R2−N2)−(R1−N1)] separates the tax contrast without rebates from the change in rebate effects; it is not a capitalization Shapley decomposition. Percentage reporting must use a common stated denominator for the identity. Introducing R1 from historical N1 is itself a fiscal-policy intervention. No policy has been submitted or certified under this new matched contract yet.


## 00:59 EDT — exact-loop probes pass; scheduling adjustment

Both first/last probes17303124 completed in56m26/56m28. Independent local collector validation accepts all gates, shared anchor contract and eight remote receipt hashes. Both own-price signed excess-demand derivatives are negative. Evidence: horizon100_smoke_verification.json and full per-case fit/parameter tables under meeting_receipts/horizon100_prices_01/sequential/{0,99}. These remain supplied-price probes, not equilibria.

The remaining98-case panel released automatically;32cases were healthy with all heartbeat ages at most15seconds. At the slower measured56-minute mapping time, retaining32slots would make a six-path root unlikely to finish by08:30. The author explicitly authorized extensive cluster parallelism. Increase only ArrayTaskThrottle from32to98 on the existing array: same98cases, same single-core16GB requests, same immutable source/contracts, same60/65-minute case limits and success dependencies; no additional solves or retries. Maximum reserved panel resources become98cores/1568GiB across nodes. Observed cpu_short snapshot has5152idleCPUs; actual scheduler/QOS availability controls allocation. Expected panel work remains about92CPUhours for98cases; more concurrency reduces elapsed time, not total solves. If all remaining cases start promptly, completion near02:00EDT and a five-to-six-hour root could finish near08:00; no guarantee. Original21000-second root ceiling remains sufficient for six measured3374-second paths with about754seconds overhead margin; investigate if path timings rise.

Scheduler update returned an unspecified-error line, but subsequent scontrol confirmed ArrayTaskThrottle=98 and squeue independently confirmed all98cases RUNNING. The mutation succeeded; no retry or duplicate array was submitted. Current collector/root remain dependent.


## 01:34 EDT — partial collection and runtime watch

Twenty of100receipts validate through the unchanged collector.load_case and full shared-contract comparison, including the two original smokes; no failures at this check. Partial receipt/provenance packet: horizon100_panel_partial_validation.json. Most first-wave cases finish48–52minutes. Case20oncs754 reached65/100forward dates after3123seconds, compared with16at2162seconds; extrapolated completion exceeds3600seconds. It is making progress, not stalled. No premature rerun or running-source/contract mutation. Check actual failure and resource evidence before any recovery. A one-time allocation on a demonstrably faster node could test a hardware-runtime explanation while preserving the exact source,60-minute cap and economic contract, but it has NOT been submitted and needs a complete explicit collector/recovery path preserving original failure artifacts. Temporary monitor interval5minutes until this near-timeout panel phase is resolved, then15minutes again.


## 01:41 EDT — reviewed single-case timeout recovery

Case20failed after60m16 with exit124, at forward date90/100, error wall-time budget exhausted; CPU59m45 and5.36GBRSS. It made continuous progress and did not report a numerical-gate failure. Other31first-wave cases completed under the same contract, including case16oncs693 in48m19. The specific revised hypothesis is node/allocation-dependent mapping speed, not a changed numerical problem. Submit exactly one fresh case20 tocs693, which has4unallocatedCPUcores and sufficient scheduler memory at review. Keep the original sourceE, exact input contract/hash, price perturbation, single core,16GB,3600-second internal/65-minute Slurm caps and every numerical gate. New cache/output directory; original failed output remains untouched. Expected48–56minutes based on that node; one attempt only, no automatic retry on another timeout. Added budget at most65CPUminutes. This execution placement is an operational test, not proof of the cause.

Recovery output: E/output/horizon100_price_recovery_01/sequential/20. Do not use the old collector blindly: its original case20directory has a failure. After the full original panel ends, create an explicit reviewed collection list replacing only any individually validated recovery directories, retain all100coordinates and source contracts, queue that collector, then repoint the existing pending root17304265 to its success. Cancel the obsolete pending collector only after the root dependency is updated. No root duplication or weakened validation.

Recovery17309087submitted; remote launch-script SHA256cc11688a40705f784fdb8e0f02c93280fa0765d13d0ddb796343cd8aa130bbed matches local. Failure JSON/contract/last completed row collected under meeting_receipts/horizon100_prices_01/sequential/20. Root and original collector remain pending; dependency rerouting is outstanding.


## 02:00 EDT — original panel verified; collection dependency repair

Alloriginalprobesended:99successfulandonlycase20timedout. Unchangedremoteandlocalcollector.load_case validates all99sharedscientificcontracts and numericalgates;297remoteprovenancehashes and396localavailableartifactpins match. Full CSV/JSONreceipts collected; binariesremainremote. Updated horizon100_panel_partial_validation.json records everycase.

Explicit100-coordinatecollector17309766queuedafterokrecovery17309087. Recipe substitutesonly20andpreservesallvalidator/source/outputpins; scriptSHA318c8bf448c8ea9378b74ec8d1a81307f5a1c60da41297df64beddbcbb54d662. It will reject anymissingorinvalidrecovery; noconditioning/rootclaimuntilcomplete.

Two attempts to update pendingroot17304265dependency returned unspecifiederror and independent scontrol/squeue both confirmed the olddependencyremained. Do not clear its dependency and risk an ungated start. Instead cancel the never-started root and obsoletecollector17303126, then queue one replacementroot with the unchanged exactlauncher and afterok17309766. Oldroot/collectorusedzeroruntime;this changes schedulingonly,addsnomodelsolveandcreatesnoduplicate. Rootlauncherlocal/remoteSHA829f629a98ee51e8592b1da71d7dcb51ee2b7409ffdd69388e2d68cb38d4bd41. All21000second/sixevaluation/sourceF96a41873/gaterequirementsunchanged.

Oldcollector17303126androot17304265confirmedCANCELLEDwith00:00:00elapsed. Replacementroot17309773submittedafterokcollector17309766usingunchangedlauncher/source/budget. Recovery17309087stillrunning. Newchain:17309087→17309766→17309773.


## Policy handoff details checked while the recovery runs

The saved `initial_2023.pkl.gz` contains parameters, wealth grid, `initial_state`, `demographic_primitives`, inherited `supply_rule`, and the originating contract hash. Its state is explicitly before 2023 choices, after the observed historical age bridge. `run_e5f_matched_pf_baseline.py` saves this object from the final historical state and the 2023 person/head state, then verifies the checkpoint reload. It intentionally writes `equilibrium_certified=False` and `policy_announcement_included=False` even when the same evaluator is called inside a price root. A policy handoff must therefore link the chosen final/reproduced root evaluation and its external market/horizon certificates; it must not flip those flags or regard any supplied-price checkpoint as an equilibrium.

`solve_person_funded_path` is specifically the equal-rebate solver. At each trial it sets target transfers to property-tax revenue divided by household heads and drives both market and fiscal residuals to zero. Passing zero initial transfers does not turn it into a no-rebate policy: it will update them toward full rebates. Use this callable for R1/R2 only. N1/N2 diagnostic comparisons require the fixed-zero-transfer price evaluator and an explicitly declared use of retained revenue. The existing person evaluator reports `government_budget_residual = property_tax_revenue - equal_transfer_outlays`; that residual cannot be called a balanced no-rebate fiscal ledger without recording the corresponding non-rebated spending/use. Do not silently route N1/N2 through the funded solver or weaken its fiscal gate. The author's main comparison remains R2 minus R1, with equal rebates in both paths; the optional N1/N2 mechanism decomposition is additional work.

These are source-verified adapter requirements, not implemented or simulated policies. Relevant source anchors: baseline checkpoint saving at lines 308–331, funded-path target transfers at policy-driver lines 290–310, and person-path fiscal ledger at lines 908–911. No model source changed. Recovery was healthy at 25 minutes and forward date 15/100 at 02:07 EDT.


## 02:34 EDT — full price matrix verified; root running

Recovery17309087completed45m54, collector17309766completed5seconds. Root17309773RUNNINGcs616 sinceabout02:28EDT,58testsPASSandfirstbackward/forwardpathhascurrentheartbeat. Full local reconstruction of100by100J exactlyequalsremotepacketSHA750f83b982bec4a4d401cc1187f4e2b7ed434b4209647a206565d6f9f6f86657. Condition3.1395,rank100,singularvalues.9023to2.8329,all100ownderivativesnegative; unconstrainedlinearNewtonmaxlogstep.1997exceedsroot.10stepcap,whichremainsunchanged. Goodlocalpriceconditioningisnotparameteridentification,globalconvergenceorhorizoncertification. Originalfailed20ispreserved;newrecoverytimeisconsistentwithallocation-dependentperformance,butdoesnotisolatehardwarecausally.

All1212anchor/probefitrowsandweightedlossesrecomputed; target/weightidentityverified. Saved12by100moment-priceforwarddifferences reuseexistingcasesandrequirednoadditionalsolves. Theyareatthenonequilibriumanchor; do not silently use them as derivatives at a newequilibrium. ExistingMtheta-minus-Mx-Fx-inverse-Ftheta identificationplanstillrequiresparameterprobesandappropriatebaseline. Monitoringreturns15minutesduringthelongroot.


## 03:40 EDT — first full root evaluation reproduces the anchor

Evaluation 1 completed in 3310 seconds. All numerical mapping gates pass; the residual vector and loss match the original supplied-price anchor exactly. Target-fit, parameter and measurement files are byte-identical. The transition CSV differs only in the ACS audit source path in four historical rows, from immutable snapshot E to F, with identical data hashes and all remaining nested audit fields equal. See horizon100_initial_root_reproduction.json. This verifies the runtime/source-wrapper adaptation on an actual full evaluation; it does not establish equilibrium. Evaluation 2, the first changed-price trial, is now running.

Monitoring detail: the shared progress dictionary retains completed_dates=100 and current_year=2403 from the prior trial during the next trial's backward phase. Interpret those fields only when phase is historical_forward and associate them with the current evaluation. A fresh heartbeat in historical_backward_and_forward does not mean the new trial already completed 100 dates. No running source is changed for this reporting convention.


## 04:29 EDT — first price step improves market clearing

Trial2passedallmappinggates in3256seconds and reducedmaxmarketgap .1819903→.0698461 (61.62%reduction), stillabove2e-4. Lead recomputed100datedresiduals, verifiedartifactpins/all12weightedlossrows and unchangedtarget/weight/parameter/sourcecontracts for bothcompletedtrials. Bothterminaldistancesremainnot_converged. Trial3isrunning; preserveallgatesandsix-pathbudget. No newparameterestimatesorequilibrium/policyclaim. Evidence:horizon100_root_progress_review.json.


## 05:22 EDT — third trial and reusable provisional readout

Trial3passes mapping gates, reducesmaxmarketgap to.0191497, and leaves onlyunitrentoutside terminal-distance tolerances(.0125439vs.01). HouseholddistributionL1.0175702nowpasses.02. Noneofthiscertifieshistoricalhorizonstabilityorfinalmarketclearing. Trial4running. review_horizon100_root_progress.py validatesallcollectedtrials andbuilds HORIZON100_PROGRESS.md withcompletefit/parameter/restrictionandtailtables, selectingbysmallestmarketresidual. Actual3trialverificationpassed; labelscheckedagainstauthoritativestatistic definitions,includingp90/p50wealth/incomeratio. This is a morningreadoutbuildingblock,notnewcalibration.


## 06:18 EDT — fourth trial and remaining horizon work

Trial4passed mapping checks and lowered marketgap to.0022353236;trial5running,sixthcallreservedforfreshreplay. Fulltrialreview/tablebuilderpassesfourcases. Unit-rentterminaldistance.0108614stillfails.01whileallotherrecordedtailchecks pass. Keepthatgate. Longer-horizoncheckcannotbesilentlycalledcertifiedfromnearstationarypersonsorhouseholdL1alone. Currentrun_history_probe caps100dates and itsCLI3600seconds; tostudy128dates, prepareanisolatedsourcechangeextendingonlyoperationalbounds and a hash-pinned explicit pricepath through the existingprices_override hook. KeepF/runningrootimmutable. Aconstant-terminalprice extensionis initiallya prescribed-price diagnostic; itestablisheshorizonmarketvalidityonlyifalldatedresidualspass,andotherwise needsitsownpriceroot. No H128codechangeorjobwasmadeinthisturn.


## 07:12 EDT — bounded continuation after the reserved replay

The fifth price trial passes mapping checks and reduces the maximum market residual to 0.000452719, still above 0.0002. The sixth evaluation is the already reserved fresh replay, not another improvement step. Four successive price updates improve the residual; there is evidence for one more numerical step rather than a change of model or tolerances.

Prepare one continuation of at most three complete paths: fresh evaluation of the reproduced best, one Broyden step using its saved updated Jacobian, and fresh final reproduction. This uses the already implemented and tested load_restart contract; the earlier 28-date continuation exercised the same exact structure. No source changes. Expected 2.3–2.8 hours on the observed nodes, including 600 Bellman calls; hard internal budget 10800 seconds, Slurm 11100 seconds, one core and 16GB on previously faster cs693. This is a new bounded numerical round after review, not a replay of a failed candidate. No automatic further continuation is authorized by this recipe.

The launcher depends on parent17309773 completing successfully, repeats58tests, verifies the exact raw parent contract and identical source map, pins all three completed parent receipts, and requires successful independent parent reproduction through load_restart. It skips computation if the parent already converged, and refuses duplicate contracts or output directories. Same targets, parameters, weights, fiscal/demographic closure, numerical gates and source F96a41873. It writes a separate continuation directory. The original six-call budget is untouched. This continuation will likely run beyond the08:30readout; report that honestly and still deliver the scheduled partial assessment.

Submitted continuation17319026 at07:15EDT with afterok:17309773. Local shell and inline-Python syntax checks pass; remote launcher SHA256a0c04fee368c33e2900af773e1f5a6ed0f200e77ad8d96679ded4c8108591a4e matches exactly. The queued allocation repeats the existing58-test exact-loop startup before checking completed parent receipts. No new economic source was staged.


## 08:06 EDT — exact final replay and continuation underway

Parent17309773COMPLETED5h24m28, status evaluation_budget, marketgap.0004527188458867619 and finalreproduction0. Lead checked all six trial mapping gates/artifact pins/12fit rows and independent100-date residuals. Best5/replay6 prices/residuals match exactly; fit/parameters/measurement/transitionCSV bytes identical. observed_dates differs only in100 elapsed_seconds fields; every other field matches. Completedroot contract/summary/history hashes match remote and continuation restart pins. Reusable report driver now verifies this replay and writes horizon100_final_reproduction_review.json, without claiming equilibrium.

Continuation17319026RUNNINGcs693,58testsPASS and existing exact-restart guards passed; first backward/forward evaluation healthy. Final continuationcontractSHA00b08089bab8d1d61dc385726a5f009ae1c7f688d90c07fc2dcfbd29ca1500e1collected and remotely matched. Still only3paths/10800sec, no further rounds. No modelsourcechange,horizoncertificate,re-estimationorpolicy.


## Morning handoff packet

Three-page partial PDF output/pdf/matched_pf_morning_readout_20260910.pdf is prepared. All pages visually inspected;60target and37parameter/bound numeric cells match the alreadyverifiedCSVtables. Full12moments/15parameters/9terminalmetrics retained, no newgraphs or modelwork. PDFSHA2b7db84d1ea7cdd58ae9f873cfa4638fca72d71dcd28b516160b1232563b4b49. Status snapshot08:23EDT:continuation17319026firstpath44/100, freshheartbeat. Scheduled overnightmonitorpausesafterthemorningreadout; the alreadylaunchedboundedjobcontinuesindependently. Before any further action inspectits receipts, validateunchangedcontracts andcomparealltargetrows. Do notduplicateit or claimcompletedcalibration/policy.
