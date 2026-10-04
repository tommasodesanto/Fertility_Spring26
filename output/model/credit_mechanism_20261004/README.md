# Credit, family space and fertility: bounded October 4 investigation

**Start with [the decision memo](DECISION_MEMO.md), approximately four pages.**
The adopted model is unchanged. The evidence supports a conditional housing-cost
result alongside weak mortgage-fertility effects; it does not establish that
fixing the income gradient will restore a mortgage channel.

## Verified facts and reference

- Author request: investigate overnight with Fable 5.1 CLI, Astra max debate and Sol diagnostics; concise decision memo by October 4 around 10:00 America/New_York.
- Reference: October 3 post-interest, soft-financing chain 13, old empirical target contract; canonical `code/model/parameters/best_params.py`. No baseline, target, model or manuscript change is authorized as adoption.
- Stationary root: adjusted births/(2.1 * household entry) - 1. Holding H0, household population scale clears housing. This endpoint does not supply a transition.
- Coded supply uses gross rent `(R_gross - 1 + delta + tau_H) * price`; supplier-tax accounting is a separate issue from the credit mechanism.
- Original scratch scripts and results have been recovered and hash checked (`evidence/RECOVERY.md`). Their arithmetic is reproducible; their complete executed fixed-price inputs and acceptance packets were not retained. The strong-space native equilibrium attempt failed its dated-budget gates.
- Existing saved phi=0.8 versus 1 cases are fixed-price diagnostics; phi=1 aggregate-policy reporter has retained infeasible-value failures, while distribution-only checks passed. Exact common-state accounting gives raw births −0.874%, comprising +0.166% changed policies and −1.040% changed exposure. This is accounting, not a causal debt decomposition.
- The 2024 CPS comparison reproduces, and its young-age gradient differs sharply from the model under several grouping checks. It compares family money income with current Markov earnings and a different calibration vintage; it does not identify the cause of weak credit effects (`measurement/AUDIT.md`).
- The recovered benefit-world experiment adds an unfinanced, child-dependent cash floor and renormalizes child utility. It cannot isolate a correction to income sorting.
- The theoretical purchase condition uses cash before current income; the quantitative purchase screen includes four years of current income. Equal borrowing/saving rates alone do not eliminate the shadow value of a credit constraint.
- Consumption and ending assets are jointly optimized within each housing alternative. A screen-passing purchaser can still hit the ending-debt floor. The cached age-22 renter envelope checks 249 susceptible states and every positive-probability owner alternative, covering 21.5% of total first-birth susceptibility. The average child-relative current borrowing effect is negative in this slice, with heterogeneous signs by income. No before-income specification was substituted.
- The fresh, fixed-price 80%→95% financing diagnostic changes explicit births −0.250% and young ownership +10.112 percentage points. A matched 10% imposed price/rent increase changes explicit births −6.434%. Neither is a market-clearing policy equilibrium, transition or welfare calculation.
- The complete working-baseline fit is retained: loss 13.771, all 14 target/diagnostic rows, and all 31 parameter/restriction rows. The early-fertility stock miss and unavailable matched prebirth-resource validation remain material limits.

## Scope and budgets

Few questions: (1) is the weak fertility-credit link a numerical/accounting failure; (2) is wrong birth-margin sorting by income a cause or a correlate; (3) can the unchanged model make a defensible paper point, and if not what minimal evidence-grounded change is justified?

First-pass workers: 45 minutes each, exclusive output subfolders, stop on missing premise or completed bounded evidence. Fable review: at most 45 minutes per round, up to three substantive rounds. Micro solves only after lead approval of exact paired design, serial single-core locally, one solve at a time, explicit case/total budgets and standard diagnostics. No sweeps, production-policy certification or automatic changes to other jobs.

Final artifacts planned: `DECISION_MEMO.md` (roughly 2–4 pages), a compact evidence index, bounded debate rounds, reproducible diagnostic driver and essential figures/tables. Full target and parameter tables are linked/copied in supporting evidence when reporting any baseline calibration fit.

## Progress

2026-10-04 01:18 New York: startup and CLI routing; baseline source remains unchanged. Existing worktree has unrelated changes; this investigation owns only its output directory and any explicitly scoped diagnostic driver.

2026-10-04 01:45 New York: Fable 5.1 CLI round one completed on the requested model; Astra max is challenging its universal structural-null claim using the exact decomposition. Sol recovery, measurement and primary-source literature passes are complete. One additional phi=.95 diagnostic is being prepared with a matched reached-runtime source contract; the initial launch stopped before solving because macOS rejected an address-space resource limit. Its replacement uses an external resident-memory supervisor. A separate zero-solve financing census is underway. No economic specification has changed.

2026-10-04 02:00 New York: all three Fable CLI rounds and Astra replies complete. Fresh phi=.95 and price +10% cases saved full arrays and 17 standard graphs each. Explicit births change −0.250% and −6.434%, respectively. Native distribution/probability and aggregate-policy checks pass; full dated-budget/purchase auditing is unexecuted because frozen ancestry authentication reaches a preexisting deleted archive initializer (exact path in `audit/summary.json`). The author explicitly requests continued overnight work with economical token use and emphasizes simultaneous consumption, saving and purchase choices. An active goal now tracks the task; the next check verifies that joint budget and its actual binding margins. Purchase-screen counts alone are not evidence of unconstrained finance.

2026-10-04 02:35 New York: the joint-budget illustration and broader age-22 renter envelope are complete. The lead checked the implemented budget, continuation, interpolation, floor multiplier and fertility-logit inversion against source, independently reaggregated the 249-node CSV, and verified all 35 source/reference pins and five baseline copies unchanged. The compact memo incorporates this evidence. No additional solves, sweeps or specification changes were launched. Account usage reads 44% used / 56% remaining; the initial check read 42% used.

## Evidence and reproduction

| Read when needed | Contents |
|---|---|
| [Decision memo](DECISION_MEMO.md) | Main judgment, worked budget example, two fresh diagnostics and next empirical requirement |
| [Diagnostics](diagnostics/README.md) | Exact birth accounting, standard graph sets, paired controls, source pins and runtime limits |
| [Joint budget](diagnostics/joint_budget/README.md) | Consumption/saving checks, collateral shadow values, income heterogeneity and coverage |
| [Recovered experiments](evidence/RECOVERY.md) | Original scratch scripts/results, tax arithmetic, bundled changes and failed acceptance gates |
| [Full fit](evidence/baseline/target_fit.csv), [parameters](evidence/baseline/parameters.csv) | Every target row and estimated/fixed parameter, restrictions and bounds |
| [Measurement](measurement/AUDIT.md), [prebirth feasibility](measurement/prebirth/README.md) | CPS reproduction and robust comparisons; precise missing PSID resource inputs |
| [Purchase census](census/README.md) | Origination screens, realized debt-floor incidence and common-state room use |
| [Literature](literature/EVIDENCE.md) | Eight primary sources and the limits of their relevance |
| [Full audit limitation](audit/README.md) | Unexecuted canonical dated-budget/purchase checks and exact preexisting source blocker |
| [Fable final reply](reviews/fable_round3.md), [Astra final review](reviews/astra_round3.md) | Three-round debate; earlier inputs/replies retained alongside them |

Fable was run through the CLI on `claude-fable-5-1`, effort max; every round's
launch receipt and returned model identity were checked. The lead rejects its
unproved attribution of the positive policy component to one specific borrowing
channel and its binary test equating down-payment funding sources with validation
of a four-year budget screen. Universal credit irrelevance, welfare inference
from birth counts and sorting-only causation were withdrawn during the debate.
The memo states the retained conclusions rather than adopting any review wholesale.

Each diagnostic README provides the exact single-core reproduction command;
the joint-budget README also supplies a function to regenerate the canonical
17-graph packet from a saved case without solving. Fresh cases retain complete
arrays and executed inputs locally under `diagnostics/phi_095_run1/` and
`diagnostics/price_110_run1/`. These large files, raw CLI streams and private
transcript excerpts remain outside Git. Compact public inputs, source/hash
receipts, essential tables, scripts and the memo are preserved in Git. Recovery
scripts use pinned local case paths; reproducing on another machine requires
the retained numerical cases or an explicitly authenticated fresh solve.
