# Handoff: discuss the calibration plan and the theory before September 14

Prepared September 5, 2026. This is a discussion prompt and a dated orientation snapshot, not a replacement for the live decision records.

## Your role and the request

You are an independent senior quantitative macroeconomist and a critical discussion partner for Tommaso De Santo. You have access to the local filesystem but not to the conversations that produced this work. Read the evidence below, form your own assessment, and discuss it directly with Tommaso. Do not assume that the previous agent's interpretation or proposed plan is correct.

The project studies housing costs, fertility, tenure choice, and population dynamics. Tommaso presents on **September 14, 2026**. He wants to assess both the current quantitative calibration plan and the overall simplified theoretical argument, then bring your discussion back to the original agent. The objective is a clear, defensible economic contribution and credible quantitative evidence. Stronger fertility effects or a particular welfare theorem are not required outcomes.

Working project:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

This request is for **reading, independent assessment, and conversation**. Do not change code, targets, calibration, manuscript text, decision records, or cluster jobs. Do not launch another agent or an overnight search. A separate Astra task owns code correctness and simplification in an isolated checkout; the original task owns calibration assessment. Your role is economic judgment and prioritization. Read existing verification before proposing additional computation. If Tommaso later asks for implementation, establish its specific scope then.

The entire `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/` subtree is author-controlled and strictly read only, including compilation artifacts. Suggestions elsewhere are proposals, not automatically adopted manuscript text. Other tasks are active in the repository: dirty files are not yours to fix, revert, stage, or commit.

## Read this evidence in order

Start with the project instructions, then the mandatory context sequence:

- [Project instructions](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/AGENTS.md).
- [Working memory](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/memory/AGENT_MEMORY.md): focus on the September 4–5 entries and the final focused-code-repair entry. This file also contains much older historical material.
- [September 5 daily note](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/memory/daily/2026-09-05.md); read a newer daily note instead if one exists.
- [Live calibration status](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md): start with the September 5 and September 4 blocks. Read older sections only to resolve a specific provenance question.
- [Live theory decision ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/ACTIVE_DECISION_LEDGER.md): first read “Simplified theory discussion map,” especially the latest committed-transfer result and proposed route, then the branch table. Later historical statements can be superseded by the dated opening entries.

Memory and earlier reviews can be stale. The current calibration status governs quantitative results; the current theory ledger governs author decisions and the status of theoretical claims. If an apparent contradiction remains after checking dates and source evidence, identify it explicitly instead of choosing silently. Past execution permissions in a work log do not turn this discussion request into permission to resume that work.

### Core calibration reading

1. [Bounded refinement review, with both complete fit tables and every parameter/bound](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_bounded_calibration_refinement_review.md). [PDF with the standard diagnostic graphs](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/e5f_bounded_calibration_refinement_review.pdf).
2. [Concise overnight quantitative assessment](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_overnight_morning_review.md).
3. [Full independent quantitative audit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_independent_quantitative_audit.md), particularly Sections 16–18 for measurement alignment, numerical checks, and rental-opportunity sensitivity.
4. [Current code reading map](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/README.md) and [sequential package entry map](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/README.md). Read the relevant source when a claim depends on implementation; do not turn this into a second whole-codebase audit.

The refinement review's paragraph saying calendar policy storage still needs repair is superseded by the September 5 repair entry in `CALIBRATION_STATUS.md` and [the repair/replay evidence](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_policy_cleanup_verification_20260905a/README.md). Preserve this distinction: a diagnosed defect, a verified repair, and an unchanged historical calibration are different findings.

### Core theory reading

1. [Consolidated amendment proposal](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex), [readable PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/simplified_olg_amendment_proposal.pdf). This supplies the environment, household problems, compensated allocation, conditional fertility prediction, and transition/population arguments. Some presentation conventions are awaiting restoration to the author's chosen notation.
2. [Latest constrained-efficiency construction](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_constrained_efficiency.tex), [PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/simplified_olg_constrained_efficiency.pdf). This is the latest positive complete-path result under stated extra transfer/ownership permissions. Read its assumptions and complete resource account, not just the opening claim.
3. [Direct-allocation and constrained-planner comparisons](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_planner_benchmarks.tex), [PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/simplified_olg_planner_benchmarks.pdf). This explains the distinct feasible sets. Its older statement that no complete-path constrained result exists is superseded by item 2; its direct-allocation analysis remains relevant.
4. [OLG efficiency reading guide](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_efficiency_reading_guide.tex), [PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/simplified_olg_efficiency_reading_guide.pdf). The literature discussion is supporting context, not proof of the model-specific claims.
5. [Verification and review index](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/README.md), then [the complete-path review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/constrained_full_path_review.md) and [the alternative-instruments review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/constrained_instruments_review.md) where relevant.

Use these only as needed to resolve a particular issue:

- [Original independent theory review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_independent_review.tex), [PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/simplified_olg_independent_review.pdf): a dated critique, not a statement that every original objection remains unresolved.
- [Theory development and reconciliation record](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/simplified_olg_overnight_work.md): completed checks, qualifications, and resolutions of the original review.
- [Author's current analytical section](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/sections/02_analytical_model.tex) and [derivations](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/appendices/A_analytical_derivations.tex), both read only.
- [Pre-amendment theory section](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/simplified_olg_paper_theory_section.tex) and [appendix](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/simplified_olg_paper_theory_appendix.tex): comparison sources for the author's notation, separate flow utilities, explicit value functions, and variable names. Do not treat their older economics as automatically superseding reviewed amendments.
- [Author's quantitative-model section](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/sections/03_quantitative_model.tex) and [calibration section](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/sections/04_calibration.tex), both read only. Check their claims against live status; do not presume manuscript results were refreshed.

## Quantitative position as of this handoff

The retained production calibration is September 4 `task_010`, with loss **30.4829667**. The independently reproduced review candidate is `round2/task_009`, loss **23.1534447**, a **24.04% improvement in the identical twelve-moment objective**. It has not been promoted. There are eleven estimated parameters. The dated supply elasticity is externally fixed at 0.63 and tenure taste dispersion at 0.005. No targets, weights, search bounds, grids, or numerical acceptance gates changed in that comparison. The appendix below reproduces the full fit and parameter tables.

The candidate fits completed fertility and childlessness closely, but produces a first-birth housing response of **0.4660 rooms against 0.7202**, while mean occupied housing is **6.3286 rooms against 5.7800**. First births at ages 30+ and prime-age ownership also remain low. The joint problem is not simply insufficient housing in aggregate. Establish whether the remaining tension comes from empirical-to-model measurement, search limitations, or the economic restrictions on housing demand and selection. A first-birth housing event-study response is not automatically a causal estimate of fertility's response to a housing intervention.

The completed work already includes a 23-case small coordinate panel, 12 joint moves, exact historical repetitions, and independent numerical checks. Raising the first-child room intercept while reducing the per-child floor helped the fit; repeating that idea alone is not a new research design. The smaller-step Jacobian has full algebraic rank, but its weakest direction is overwhelmingly the bequest shifter and the supply-intercept one-sided responses disagree. Full rank does not establish precise identification or justify blindly inverting this Jacobian.

Young ownership remains a validation moment; its optional calibration row is default off. Earlier ACS and model percentages used different age conventions, which need explicit alignment. The first-birth rooms mapping, existing age-window definitions, late-life ownership, and rental opportunities remain substantive discussion items. The inherited old equilibrium uses supply elasticity 1.75 before the dated rule uses 0.63; the economic interpretation of that initialization is explicitly outstanding.

All existing headline policy responses belong to **retained production**, not the better-fitting candidate. The evidence suggests that housing supply has a materially larger fertility effect than dependent-child mortgage relief, even though mortgage relief changes ownership. The numerical grid/global-saving checks broadly preserve that contrast. A larger-rental-opportunity diagnostic materially changes ownership responses, and changes the already-small credit fertility coefficient. This is an uncalibrated economic sensitivity, not a promoted alternative. Full policy paths, fiscal closures, population-entry conventions, and their limitations are documented in the audit. These are not established welfare rankings or a solved perfect-foresight transition.

The calendar saved-policy defect has been repaired and independently replayed. Exact historical reproduction does not establish that every original economic assumption is correct, nor that every off-support policy is globally optimal. The new code-cleanup task must not silently determine the model specification for this discussion.

### The original agent's proposed plan — critique it, do not merely endorse it

1. Diagnose the remaining first-birth housing response and excessive mean rooms, beginning with measurement comparability and then the economic mechanism.
2. Use the completed sensitivity evidence to decide whether one additional bounded, same-target calibration round is likely to improve the important misses. Do not remove or demote a difficult target as a shortcut. Any proposed new model, target, or weight system needs a separate identifying argument and cannot inherit the old loss comparison.
3. Recompute matched policy comparisons at the selected calibration before combining its fit table with policy magnitudes in the presentation.
4. Aim to freeze the quantitative package by **September 11**, leaving preparation time before September 14. September 11 is the original agent's proposed planning date, not an author-approved deadline.

Assess whether that order is right. In particular, should the policy comparison be brought forward to determine whether the better fit changes the economic conclusion? Which measurement issue could invalidate further calibration work? Which outstanding policy/population closure matters for the claims actually intended for the talk? What minimum evidence supports using the new candidate, retaining the old benchmark, or reporting both as sensitivity? Recommend a small set of useful decisions and calculations, not an unrestricted new research program.

## Theory position and open choices

The simplified theory is a two-period overlapping-generations model with young housing/fertility decisions, tenure and financing restrictions, old-age housing and estates. Its fertility choice is one-shot completed fertility; the quantitative model has sequential births. A theory derivative cannot be read as a sequential first- or second-birth hazard effect without a mapping. Symbols such as chi and kappa can denote different primitives across the two models; define each object before making a comparison.

The desired argument has three distinct parts:

1. **Housing allocation at fixed individual fertility.** Establish whether financial restrictions leave housing with households who value it less, and whether an admissible reallocation can make someone better off without hurting anyone. The author wants both a direct-allocation comparison and a constrained comparison where transfers are followed by market clearing. A practical policy implementation is not required for an efficiency benchmark, but feasible instruments, resources, initial titles, debts, estates, and all affected households must be specified. A compensated local improvement is not a characterization of the full first best.
2. **A conditional fertility prediction.** More housing can raise fertility, but the resources paid for that housing matter. The amendment derives a local condition using the net marginal payment, plus a price-free sufficient restriction in a stationary, zero-property-tax, interior-old specialization with strict borrowing constraints and positive saving. These are conditional household results. They do not prove a uniformly positive aggregate equilibrium response to credit or supply policy. A negative response outside the condition and the effect of changing tenure composition are retained in the evidence.
3. **Transition and population implications.** The author wants fertility along a transition and the population level at its endpoint. A positive stationary population satisfies replacement under its demographic law; a permanently higher stationary fertility rate is not the intended conclusion. A finite-reform continuation argument is proved under explicit regularity, feasibility, and uniqueness hypotheses, but those hypotheses have not been verified over the full model policy interval. Numerical credit examples are illustrations, not general existence/convergence/sign theorems. The largest illustrative credit reform initially lowers fertility while raising terminal household population. Keep households, births, children, and resident persons distinct.

### The latest welfare result, and what has not been chosen

The latest note establishes a **local complete-path Pareto improvement under a specific committed-transfer class**, in a particular stationary group structure with strict branch/tenure margins and a contraction condition. It uses household-specific payments when young and old, enforceable future lump-sum taxes, and a passive initial title owner with outside bond access and participation in the fiscal account. Young housing increases at each finite date; future housing markets clear; every initial claimant is counted; the allocation returns to the original steady state. Individual fertility is held fixed. The proof and checks also give a strict initial-young gain variant.

Those are substantive financial and ownership permissions. Funding a young grant through a committed old-age tax provides public financing across ages even if the private mortgage share is unchanged. The author has **not** accepted those permissions as the final welfare benchmark. Do not label this a pure price externality or a result under one-time cash gifts. The same verified regime has a local obstruction for one-time gifts; that is not a global impossibility theorem. The expanded-instrument fallback is not established as the uniquely minimal intervention.

Direct allocation remains a distinct valid route where its feasibility is established. Real reallocation costs remain tentative: the cost base, incidence, functional form, and magnitude are not selected. Tenure segmentation must be respected, but its exact planner-feasibility interpretation is still part of the discussion. Capital-gains taxation is deferred; it is distinct from property taxation, which remains in the maintained environment except in explicitly labeled special cases.

The immediate theoretical decision in the live ledger concerns **the financial powers of the compensation planner** (IDs W1/W4/I1). Full Pareto coverage includes the initial old and all future cohorts, but does not require changing every cohort's allocation. A one-time intervention is an instrument restriction, not the definition of the OLG welfare test. Diamond-style overaccumulation and welfare comparisons with different populations are separate questions, currently parked. The reading guide points to the relevant literature; verify any new attribution from the primary work instead of invoking a citation as proof of our result.

Tommaso wants a simple illustrative theory, explained plainly, with short transparent proofs. He specifically requested preserving his notation, separate flow utilities, explicit value functions, and variable names. The desired two conceptual figures concern housing misallocation and adjustment between steady states; a simulated numerical path is not automatically the second theory figure. Do not allow a search for a stronger theorem to obscure what the model already establishes.

## Questions your assessment should resolve

- What is the strongest accurate economic claim we can currently make, separately for allocation, constrained efficiency, fertility, and population? Which statements are proved conditionally, illustrated only, unproved, or awaiting an author choice?
- Is the committed-transfer benchmark a useful way to expose the housing/credit distortion, or do its added financial powers make a simpler direct-allocation result more appropriate for this paper? Identify the decisive assumption and its economic role before discussing technical extensions.
- What genuinely links the theory to the quantitative mechanism? Does a weak fertility effect of mortgage relief contradict the conditional theory, or identify a missing premise, selection effect, or weak empirical channel? Do not presume an answer.
- Can the aggregate housing level and first-birth housing change be jointly disciplined under the current measurement and model? What existing evidence separates a measurement failure, a search problem, and a model restriction?
- Which next steps materially improve the September 14 presentation? Which are essential now, useful robustness, or post-presentation research? Include explicit stopping conditions and the smallest decisive calculation for any proposed run.

## How to conduct the conversation and hand it back

First give Tommaso a concise independent assessment: your central diagnosis, the two or three decisions that matter most, and any material correction to this handoff. State what you actually read and what remains unverified. Then begin with **one** substantive question, explain why it matters, and discuss that question before moving on. Do not present an eighteen-item questionnaire or force a choice among unexplained labels. You may recommend a different starting point from the original agent if the evidence warrants it.

Use plain economics, define technical objects, and distinguish quantities from labels. When stating a calibration result, retain the complete fit and parameter/bound tables or point to the exact complete artifact. Do not infer causality from a fit, identification from a parameter count, or a theorem from a numerical example. When Tommaso explores an alternative, do not record it as accepted unless he explicitly adopts it.

At the end, provide a compact **return memo in the conversation for Tommaso to bring back**. Include:

- Your agreed central interpretation and any disagreement with the original plan.
- Explicit author decisions, separately from your recommendations and unresolved questions; retain W1/W4/I1 and other ledger IDs where applicable.
- Any error or source discrepancy, with its exact absolute path and equation/section or line reference.
- An ordered pre-presentation plan with essential tasks, stopping rules, and deferred work.
- The strongest wording currently supported for the theory and quantitative findings, with its qualifications.

Do not update the live ledger or code yourself during this discussion. The original agent will reconcile your return memo with Tommaso before implementation.

## Additional exact paths for checking evidence

- Selected calibration outputs: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/README.md`.
- Complete selected-case fit: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/round2/task_009/target_fit_long.csv`.
- Complete selected-case parameters: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/round2/task_009/parameter_table.csv`.
- Standard candidate graphs: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/`.
- Numerical robustness receipts: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/README.md`.
- Dated calibration/measurement driver: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py`.
- Policy driver: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_post2023_policy_mechanisms.py`.
- Population drivers: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py`.
- Theory checking code, for inspection: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_amendments.py` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_planner_benchmarks.py`. These drivers write outputs when executed; read saved receipts first.

## Appendix: complete calibration tables at the handoff date

The following is copied from the bounded refinement review named above. Consult the live status and that review if newer work supersedes this snapshot. Gaps are model minus target; each contribution is its existing weight times squared gap. Synthetic scales must not be described as empirical standard errors, and the scalar loss is not a specification-test statistic.

### Complete target fit: verified candidate

Gaps equal model minus target. Contributions are weight times squared gap.
The weights retain their existing empirical or synthetic interpretations;
this loss is not a specification-test statistic.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.921698 | +0.003698 | 1425.74 | 0.0195 |
| Childlessness | 0.188000 | 0.189859 | +0.001859 | 17180.7 | 0.0594 |
| Mean age at first birth | 26.044627 | 26.251703 | +0.207076 | 44.4444 | 1.9058 |
| First births at age 30+ | 0.260327 | 0.237163 | -0.023165 | 10000 | 5.3660 |
| Rooms response to first birth | 0.720246 | 0.466044 | -0.254202 | 137.565 | 8.8893 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.379198 | +0.011499 | 2958.51 | 0.3912 |
| Family ownership gap | 0.167662 | 0.169578 | +0.001917 | 14229.6 | 0.0523 |
| Ownership, ages 30-55 | 0.575472 | 0.546618 | -0.028854 | 1207.85 | 1.0056 |
| Mean occupied rooms | 5.779970 | 6.328594 | +0.548624 | 11.9732 | 3.6038 |
| Wealth / annual earnings | 6.873100 | 7.031196 | +0.158096 | 6.28767 | 0.1572 |
| Bequest flow / wealth | 0.008800 | 0.008485 | -0.000315 | 5.16529e+06 | 0.5139 |
| Older wealth/income p90 / p50 | 3.448111 | 3.303595 | -0.144516 | 56.9598 | 1.1896 |

### Complete target fit: retained production

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.922907 | +0.004907 | 1425.74 | 0.0343 |
| Childlessness | 0.188000 | 0.189407 | +0.001407 | 17180.7 | 0.0340 |
| Mean age at first birth | 26.044627 | 26.256032 | +0.211404 | 44.4444 | 1.9863 |
| First births at age 30+ | 0.260327 | 0.237579 | -0.022748 | 10000 | 5.1748 |
| Rooms response to first birth | 0.720246 | 0.436418 | -0.283828 | 137.565 | 11.0821 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.403441 | +0.035741 | 2958.51 | 3.7793 |
| Family ownership gap | 0.167662 | 0.161699 | -0.005963 | 14229.6 | 0.5060 |
| Ownership, ages 30-55 | 0.575472 | 0.544488 | -0.030984 | 1207.85 | 1.1596 |
| Mean occupied rooms | 5.779970 | 6.317291 | +0.537321 | 11.9732 | 3.4568 |
| Wealth / annual earnings | 6.873100 | 6.932652 | +0.059552 | 6.28767 | 0.0223 |
| Bequest flow / wealth | 0.008800 | 0.008433 | -0.000367 | 5.16529e+06 | 0.6971 |
| Older wealth/income p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.9598 | 2.5503 |


### Complete parameter comparison

| Parameter | Production | New candidate | Bounds / restriction | Near bound? |
|---|---:|---:|---|---|
| Annual patience $\beta$ | 0.995117 | 0.995739 | [0.94, 0.9995] | No |
| First-birth dispersion $\kappa_1$ | 2.168173 | 2.168173 | [0.02, 50] | No |
| Later-birth dispersion $\kappa_C$ | 1.736471 | 1.736471 | [0.02, 50] | No |
| Owner-service premium $\chi$ | 1.043472 | 1.043472 | [0.1, 5] | No |
| Supply intercept $H_0$ | 14.562959 | 14.562959 | [0.2, 80] | No |
| Bequest strength $\theta_0$ | 0.528428 | 0.549189 | [0, 8] | No |
| Bequest shifter $\theta_1$ | 0.107249 | 0.107249 | [0.02, 16] | Raw scale only |
| Per-child room floor | 0.282210 | 0.248877 | [0.1, 1.8] | No |
| First-birth fixed cost | 4.559138 | 4.559138 | [0, 8] | No |
| First-child room intercept | 0.364931 | 0.464931 | [0, 0.5] | No |
| 2007-2023 child-value change | -0.328714 | -0.328714 | [-1.5, 0.2] | No |
| Tenure dispersion $\kappa$ | 0.005000 | 0.005000 | Externally fixed | Not estimated |
| Dated supply elasticity $\eta$ | 0.630000 | 0.630000 | Externally fixed | Not estimated |

The old fertility-preference intercept is re-normalized to replacement for
every parameter vector. New candidate values are
$\psi_{2007}=0.28947879$ and
$\psi_{2023}=-0.03923512$.
The stored near-bound flag uses raw parameter units. For the bequest shifter,
interpret it alongside its position in the logarithmic search interval; a
raw-scale flag alone does not establish a binding search constraint.

The supply restriction 0.63 applies to dated markets. The inherited initialization
retains elasticity 1.75 before the dated supply intercept is normalized. Both
values and that order are unchanged across this comparison; whether this is the
intended economic contract remains a discussion item in the separate code audit.

