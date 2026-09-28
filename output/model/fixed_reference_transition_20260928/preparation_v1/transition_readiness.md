# Transition readiness: block0506

September 28, 2026. Initial preparation budget: 45 active minutes, excluding
external queue time. This note prepares experiments; it does not report a
completed credit policy or historical equilibrium transition.

## Judgment and exact reference

The dated household/population machinery is reusable. The existing production
entrypoint is not ready to run the requested experiments unchanged: it pins an
older checkpoint, forces an approximate natural-credit option incompatible
with block0506's owner-debt flag, and inherits an elastic supply curve. A
small output-only bridge now imports the **unchanged frozen September 28
runtime**, authenticates the complete saved parameters, and bypasses all
fertility normalization. Shared code and measurements remain unchanged.

Use **2007 stationary reference — block0506, September 28 verified export**.
Its common-primary loss is 19.581; this is a reference identity, not a new fit.
The current reference has 14 reported targets, 10 searched parameters, and a
separately normalized child-benefit parameter. Counterfactuals freeze every
preference, including \(\psi_{child}=0.136\) (rounded here; exact saved value retained). The full fit and
31-parameter tables are linked in the parent README and reproduced in the PDF.
The estimated first- and later-birth taste scales carry near-lower-bound flags.

The source manifest is `07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d`;
the identification contract is `68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf`;
the objective is `04528ece6513e4da9436bb8b35b54cd7e9cecd3a6695d050a0a7aea5a71d5e60`.
The Torch project is `/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project`.
The reference checkpoint remains there. Constructor defaults and an improved
future calibration cannot replace this state.

## What the estimated historical transition was

The September 12 user messages chose successive unexpected preference shocks:
at each date households believed the newly observed preference would persist.
Only the first realized period of each forecast was carried into the next
surprise. Four scalar preference levels were fitted to national period-fertility
windows, keeping the then-calibrated 2007 structural parameters fixed.

| Surprise date | Observation window | Old preference | Fertility target | Old short-forecast result |
|---|---|---:|---:|---:|
| 2007 | 2008-2011 | 0.129 | 1.975 | 1.973 |
| 2011 | 2012-2015 | 0.117 | 1.861 | 1.861 |
| 2015 | 2016-2019 | 0.106 | 1.755 | 1.755 |
| 2019 | 2020-2023 | 0.092 | 1.646 | 1.646 |

Full precision is retained in
`output/model/e5f_final_night_20260913/history_A0_6/realized_fit.json`.
The source receipt is `corrected_history_source_v2_receipt.json`, SHA256
`087a8a34f3b2f909d7cb624cf457ad9a667b7815f0ee8363b42877ae9286aeac`;
the fitted-series SHA256 recorded in the original-queue specification is
`e970dee333af959fc2f377e5fb11f55df4ddd7800e577550f5d90bde0f28e9fd`.
The initial calibration came from `corrected_initial_replay_v6/summary.json`.
Full source and target paths are in `historical_source_identity.json`.
These are retrieved identities, not newly recomputed historical hashes.

Two separate experiments reused those estimates: a single permanent full
decline, and the September 13 author-requested announced four-step path. The
latter revealed all four changes at time zero and used 104 decision dates.
Neither should be confused with the intended successive-surprise history.

The inherited population law used adjusted births divided by 2.1 and a legacy
entry queue, with no immigration or historical age reweighting. Its payroll
tax was 0.179, property tax was 1% annually with equal rebates, and housing
supply elasticity was 0.630. The phrase 'fixed supply curve' in those notes
does not mean fixed physical stock. These fiscal, entry, earnings and utility
objects are not block0506's contract.

The later 104-date successive-surprise refit ended September 16 with **zero
accepted fitted shocks**. One finite market/fiscal root passed but missed its
fertility target; its terminal household mass still differed by 1.742% and its
distribution by 9.446% in relative L1 distance. Another candidate failed the
unchanged root gate. Thus a short fitted history exists; a certified long
historical transition does not. See the final assessment in
`output/model/e5f_original_queue_20260913a/long_successive_refit/README.md`.

The measured block-Toeplitz Jacobian was a useful numerical improvement. A
10-date continuation cleared its finite market/fiscal system, but terminal
population was still far away. It was a measured approximate Jacobian, not a
complete fake-news construction. Recompute derivatives for block0506 and its
two-block price/pension system; do not reuse the old three-block Jacobian,
which included a rebate equation. See `docs/model/e5f_sequence_space_prototype.md`.

## Closures and their status

| Object | Classification | Contract for this preparation |
|---|---|---|
| Preferences | Estimated or calibration-normalized, then fixed | All saved preferences; no post-shock normalization. Historical preference shocks are a separate, pending specification. |
| Earnings | Externally estimated, fixed | Approved B15 15-state Markov process and saved age profile; use completed measurement audit, not constructor defaults. |
| Initial population | Empirically normalized reference | Exact saved pre-choice distribution, mass one. No age reweighting, reset or rescaling along a path. |
| Adult entry | Author-fixed reduced-form normalization | Half of each birth vintage enters after 16 years and half after 20; adjusted children map to households at 1/2.1. Preserve both raw and adjusted queues. This is not a new estimate of child survival. |
| Entry wealth/income | Empirically normalized; coupling approximation retained | Full saved conditional entrant distribution and grid. No zero-wealth fallback or frontier censoring of inherited population. |
| Survival | Externally fixed | Exact saved age survival and terminal death; no mortality change. |
| Geography and outside entry | National target scope; closed-path assumption explicit | One pooled market; no spatial migration in I=1. No new immigration, retention parameter, age bridge or old quota defaults. National calibration does not itself validate a population forecast. |
| PAYGO pension | Empirical baseline normalization; fixed tax closure | Hold saved payroll tax at 8.028%; solve each date's equal pension from actual worker/retiree masses. Baseline benefit 0.918 is a starting value, not a fixed transition benefit. |
| Property tax | Externally fixed | Saved annual rate 1.060%, period rate 4.239%, zero household rebate. Do not import old equal rebates. |
| Estates/entry funding | Outstanding substantive settlement; provisional reference ledger retained | Net positive estates fund actual next-period positive entrant assets; residual sink; funding shortage and negative estates fail. Donor utility stays unchanged. Lender counterparties, recipient law and physical settlement remain unresolved. |
| Housing, credit experiment | Estimated intercept, fixed external elasticity | Retain saved absolute supply curve (elasticity 0.630). Do not rescale it by population. |
| Housing, fixed-stock experiment | Author-fixed counterfactual | H equals actual reference supply at the reference equilibrium price, not H0. Prices/rents and individual housing/tenure choices may change. Constant gross stock entails replacement of depreciation; it is not zero gross construction. |
| Credit experiment | Author-adopted comparison; implementation outstanding | Replace artificial purchaser and incumbent debt restrictions by lifetime no-default solvency; keep interest, repayment, prices and death settlement. No arbitrary debt-floor substitute. |
| Terminal population | Endogenous equilibrium object; outstanding | Solve renewal, PAYGO and housing jointly at fixed preferences. A normalized stationary price root is not enough. |

The reference has a small measured renewal discrepancy: entry exceeds potential
birth-derived entry by \(4.889\times10^{-8}\) per reference household. Seed
prehistory from actual saved entry, then let actual births enter the queue.
Report resulting no-shock drift rather than adjusting child benefit or queues.

## Natural solvency and the new steady state

A natural borrowing limit is the most debt a household can repay under every
modeled future event with positive probability. Construct feasible sets backward
at each existing decision node, preserving the order of fertility realizations
and subsequent housing choices. Feasibility must be a separate Boolean object;
large negative utility is not itself proof of infeasibility.

For each chosen saving/tenure branch, require the resulting successor state to
be feasible for every reachable income and family-state outcome. Whenever death
has positive probability, require the inherited timing's net liquidation estate
\(b'+(1-\psi_{sell})q_t h\geq0\). The terminal condition is the same repayment
condition. With positive death risk and no default or life insurance, this
condition can itself exclude unsecured debt even when future earnings are
positive; removing an artificial limit does not authorize unpaid death-state
liabilities. The current model values this estate at today's price after saving;
changing its timing would be an additional economic change.

The current `native_solvency_credit` uses value cutoffs, first feasible grid
nodes and the grid's lower end; it is a prototype, not yet a certified natural
limit. It also conflicts with saved `native_due_stayer_credit=True`. Replacing
the incumbent-owner restriction is part of the authorized credit change, but
must be explicit. Setting financed share to one only removes a down payment;
it does not remove the collateral-based debt limit. Verify the new feasible
sets against terminal/worst-income cases and an expanded/refined debt grid,
without altering the reference entry distribution or relaxing occupied gates.

For a closed positive stationary population, let \(B(q,p)\), \(E(q,p)\) and
\(d(q,p)\) denote adjusted births, new household entry and housing demand per
normalized household, with pension \(p\). At fixed preferences the endpoint
must satisfy
\[
B(q,p)/(2.1E(q,p))=1,\qquad
\tau_{pay}Y_W(q,p)=p N_R(q,p),\qquad
N d(q,p)=H^S(q).
\]
Here \(Y_W\) is gross worker earnings and \(N_R\) is retiree exposure in the
normalized distribution. After solving renewal and PAYGO, housing determines
the level \(N=H^S(q)/d(q,p)\); for fixed stock replace \(H^S(q)\) by \(\bar H\).
One-step native distribution and queue reproduction must then pass. If no
positive root is found, report that failure and the search range; do not force
replacement fertility, rescale entrants or claim nonexistence from a timeout.

## Supply scaling: a candidate endpoint, not a transition

For a pure 10% increase in the housing-supply intercept, a candidate endpoint
keeps prices, pension, household policies and the normalized distribution
unchanged and multiplies population, both entry queues, births, deaths, housing
demand and fiscal flows by 1.1. This follows from the level-linear forward
population operator and the per-household earnings and entry laws; both sides
of PAYGO and the provisional estate-funding ledger scale by the same factor.
The absolute supply schedule also scales by 1.1. The reference's small renewal
residual scales in levels and is unchanged in relative terms.

This algebra applies to the present single-market, exogenous-earnings closure
with no outside entry and no fixed aggregate transfer. A fixed immigration
flow, aggregate fiscal grant, population-dependent wages or amenities, a
non-scaling estate allocation, or normalization of population at each date
would break it. The reviewed estate ledger is homogeneous in its economic
flows; its absolute numerical tolerance is not an economic transfer. Native
one-step and scaling checks are still required, and no uniqueness, stability
or transition claim follows. The economic-analysis chat owns this supply
comparison; no additional supply solve was launched here.

## What each calculation can establish

| Calculation | What is solved | What it cannot establish |
|---|---|---|
| Prescribed-price impact | One dated household response at the exact inherited state with explicit continuation assumptions | Market clearing or a permanent-policy equilibrium |
| Heuristic partial path | A stated short price/pension guess, with households and population carried forward | An equilibrium price path |
| Cleared finite path | Backward household choices and forward population jointly clear housing/PAYGO at all included dates | Correct endpoint or sufficient horizon |
| Genuine equilibrium transition | Cleared finite path plus terminal distribution/queue approach and stable early responses under horizon extension | Uniqueness without further evidence |
| Demographic steady state | Joint renewal, fiscal and housing conditions above | Reachability from the inherited initial state |

The no-credit-limit experiment starts at the frozen 2007 distribution. It is
not the old 2007-2023 preference transition. The fixed-stock historical arm
requires the same agreed preference shocks in its elastic-supply control.
Old absolute preference levels cannot silently transfer across utility changes.
The author has been asked whether to leave historical shocks pending or use
their proportional declines only as an explicitly illustrative diagnostic.

## Finite budgets, gates and next steps

The isolated bridge is preserved in successive immutable scripts. The checks
are serial on Torch with the reference project mounted read-only.

| Attempt | Scope and result |
|---|---|
| 18737340 / version 1 | Stopped before any solve: an old serializer expanded large arrays instead of using the manifest's array descriptors. |
| 18737715 / zero-solve diagnosis | Located that representation mismatch; no change in the 260 saved fields. Its large raw JSON remains remote. |
| 18738048 / version 2 | All 260 fields authenticated; one stationary solve completed, then diagnostic output hit the read-only reference directory. No comparison certificate was produced. |
| 18738593 / version 3 | Imported the other chat's independently verified control: 113 exact arrays, 14 fit rows, 31 parameters, and 17 identical PNGs. The two-date operator timed out at 180 seconds before a dated audit. Six-date and fixed-stock stages were not reached. |
| 18739319 / version 4 | PASS: one date, two Bellman calls, 111.860 seconds for the map (142.013 seconds including authentication). Household/accounting and aggregate gates pass; initial/final distribution distance is 4.010e-14, housing residual 1.760e-9, PAYGO residual 3.172e-13, and projection/mass-accounting errors zero. The adjusted queue changes by 2.445e-8 as expected from the retained renewal discrepancy. |

The attempted diagnostic write came from
`_require_exact_inherited_distribution`: it records any positive mass at the
infeasible-value cutoff, including mass below the retained-tail tolerance.
A write attempt alone did not establish a failed feasibility gate. Version 3
and later authenticate the original parameters first, then redirect only
`native_inherited_distribution_evidence_dir` to the new output directory.
This is an output-only change. No reference checkpoint or scientific source
was writable. The earlier scripts and all failure receipts remain preserved.

The one-date narrowing addresses an unknown dated-solve cost; it does not
relax an economic or numerical acceptance condition. The full six-date test
that crosses both entry lags remains required before a policy path.

Every case checks source/checkpoint identity and all 260 serialized fields.
Replay must match every saved numerical array, all 14 fit/31 parameter rows
and all 17 PNG hashes. Dated checks retain household budget/debt/transaction,
occupied-value, probability, estate-funding, population-mass and zero-projection
gates. No-shock housing tolerance is \(2\times10^{-4}\), PAYGO
\(10^{-6}\), mass accounting \(2\times10^{-8}\), backward/forward reproduction
\(10^{-10}\); small drift remains reported. Dates and complete-case receipts
are written as progress, with latest and furthest-passed summaries.

Proposed subsequent finite stages, **not launched by this note**:

1. Natural-solvency unit and occupied-state verification: at most 12 constructed
   boundary cases and two fixed-price household solves, 20 minutes, one worker.
   Stop on repayment/support disagreement; retain complete controls.
2. One credit endpoint: at most 16 fixed-price stationary evaluations (including
   two repeat/one-step checks), one 45-minute Torch allocation. Bracket the
   renewal root; do not search preferences. The authenticated fixed-price control
   took 65.589 s, suggesting about 18 minutes of solving before setup/audits.
   If measured cost cannot fit this budget, return the plan rather than launch.
3. One six-date diagnostic path: at most 8 root mappings plus one exact final
   replay, 108 dated Bellman calls, a proposed two-hour allocation. The measured
   52-second policy cost implies about 94 minutes before setup/audit overhead.
   First smoke-test one exact six-date mapping, estimated at 11 minutes with
   a 15-minute cap. A measured Jacobian requires a separate counted budget;
   do not hide its evaluations. The
   six-date endpoint gap must remain visible, so this is not horizon-certified.
   Use measured no-shock mapping time to approve or shorten the budget first.
4. Historical fixed-stock/elastic comparison and longer horizon certification
   require their own agreed shock contract and finite plan. Support for a
   supply-intercept multiplier belongs in the generic housing-rule interface,
   but the other chat owns and launches its +10% comparison.

All subsequent solution reports retain the full fit/parameter tables and the
stable 17 plots. Supplemental price, rent, demand/supply, pension, births,
entry, deaths, population and boundary-distance paths are additions. No result
is promoted because a job finished or a smaller subset of gates passed.
