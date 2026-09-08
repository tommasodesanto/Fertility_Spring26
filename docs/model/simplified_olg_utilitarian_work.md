# Utilitarian housing note: work record

## Scope and time budget

September 8, 2026, starting 15:26 UTC. A focused two-hour pass, with research
and consolidation targeted by 17:26 UTC. Author requests a separate note:
direct utilitarian allocation, transfers-only constrained welfare, then
fertility along a transition separated from convergence-dependent conclusions.
Existing Pareto results remain supporting appendix material. Existing notes
and the protected author manuscript are not overwritten.

The author additionally requests a ChatGPT Pro review. The exact prompt is
`docs/prompts/oracle_simplified_olg_utilitarian.md`; its narrow source scope is
the current conventional-finance LaTeX note and the economics writing guide.
No calibration data, logs, credentials, browser history or personal records
are included. Use the visible Pro model at highest available reasoning level.

## Ownership and stop conditions

Three 25-minute Astra/max mathematical passes have distinct deliverables in
`output/model/simplified_olg_amendments/`: `utilitarian_direct_review.md`,
`utilitarian_transfers_review.md`, and `utilitarian_fertility_path_review.md`.
The lead owns the new note, welfare conventions, integration, verification,
and delivery. Each pass stops with a proved proposition or an explicit
obstruction; no numerical search or automatic restarts. Subsequent hostile
review will be scoped to the actual proposed statements.

## Decisions and questions

- Preserve current income available at purchase, old income, ordinary mortgage
  repayment, old resizing, heterogeneous endowments, ownership choice, positive
  child costs, and both physical caps in the current proposal.
- Initially fix fertility and cohort masses for welfare. Household fertility
  responses are a separate result; the planner does not value creating people.
- Equal remaining-lifetime utility and equal birth-normalized lifetime utility
  imply different initial-old weights. Keep the convention explicit and do not
  select it merely to obtain a sign. Full future welfare needs an explicit
  summability or finite-change convention.
- Direct allocation can waive household financing. Transfers-only must retain
  household optimization and market clearing; affordability alone is not
  implementation. Government borrowing/commitment and fiscal timing must be
  stated, not silently enlarged.
- The previous finite transition proofs belong to a different household
  specification. No inherited-state, convergence, or general-equilibrium policy
  sign is imported without a new proof.

## Initial analytical lead (historical)

At a stationary equilibrium, match a young type to the same type in the old
cohort. With slack young/old caps and old estate floor, the housing marginal
utility difference under equal remaining-utility weights is
\((\beta/q-1)pm+\mu L\), where \(m=1/c^2\), \(p\) is housing service
cost, \(L\) the cash requirement, and \(\mu\) the credit multiplier.
Thus \(\beta\ge q\) and strict finance appear sufficient for a direct
old-to-young utilitarian housing improvement. At \(\beta=q\) only the
credit term remains. The multiplier is a proof device; the intended theorem
should state a parameter and household-regime condition, with primitive
coverage supplied separately. The direction of the global housing optimum,
transfer implementation and fertility implications remain under investigation.

## Checked findings at 16:30 UTC

The separate source is `latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex`.
Its environment and all four household problems are copied verbatim from the
conventional-finance note; the small checker verifies that preservation.
The existing Pareto note and author-controlled manuscript remain unchanged.

1. **Direct allocation:** beta >= q is sufficient under equal weights on
   remaining utility, with the existing primitive conditions supplying a
   positive strictly constrained owner group and slack relevant caps/floors.
   The exact gap is (beta/q-1)pm+mu L. At beta=q only the mortgage term remains.
   The settlement holds consumption, estates, later real allocations and net
   external assets fixed. The conditional global young-housing ordering needs
   additional cap/floor conditions and is not asserted generally.
2. **Original-choice constrained utilitarian benchmark:** an unanticipated
   one-time balanced transfer from richer to poorer old owners in their common
   uncapped regime strictly raises welfare and preserves every aggregate and
   all current young/future choices. This is ordinary within-old redistribution,
   not the desired housing shift toward young. A primitive specialization uses
   phi=q and separated old-income groups among covered constrained owners.
3. **Young-directed transfers with committed fertility:** the two-date grant
   construction is verified, including private mortgages, estate/title sales,
   rich capped old funders, the complete government account and external
   resource settlement. Its sufficient preference restriction is
   alpha(1+omega_B) <= gamma(1-phi+q tau)/(1-q+q tau), with phi>=q and beta>=q.
   This is an actual conditional equilibrium; it is not an unrestricted-n
   equilibrium. The note states that additional choice-timing restriction.
4. **Fertility:** a current cash gift raises conditional fertility with every
   cap regime retained. Exact finite bundle and tenure-decomposition tests
   apply on existing dated paths without convergence. Endpoint population
   comparisons require positive stationary limits. Free fertility changes
   housing demand and future cohort masses, so the fixed-n transfer theorem
   cannot be imported unchanged. The original nominal payments remain funded;
   the larger grants considered separately implement a different allocation.

Independent first derivations and hostile checks are saved in the named
`utilitarian_*_review.md` files in the amendments output folder. Final assembled
transfer and fertility reviews required three small corrections, all adopted:
remove the unsupported sign of housing under mortgage relaxation; state tenure
fertility ordering for each affected endowment type; restrict the repayment
formula to uncapped young housing. A 32-identity symbolic check passed before
final formatting, without a model solve or numerical reference equilibrium.

## External Pro review

Oracle prepared and copied the original two-file context packet. Browser
submission did not complete. The author instead pasted a self-contained prompt
from the chat and reported it running at
https://chatgpt.com/c/6aa02e1e-89dc-83e9-91da-029d96f71d57 .
The manually submitted text is not asserted to be identical to the earlier
Oracle bundle. Initially, the in-app tab was absent from browser-control
inventory. It became accessible later; the completed 6 Pro response was
retrieved at approximately 16:46 UTC, with a visible duration of 29m 37s.
The response is saved in `oracle_utilitarian_response_capture.txt`; the
lead's assessment is `oracle_utilitarian_review.md` in the amendments output
folder. No recurring monitor was created and no automatic PDF preview opened.

Two independently verified improvements were incorporated: the common old
housing cap already implies the needed donor/funder marginal-utility ordering,
so a redundant assumption and funder-bound term were removed; and the original
pair of current/future grants strictly raises the treated owner's fertility
after reoptimization at fixed prices and rebates. The latter does not preserve
the original market-clearing allocation once fertility changes. Pro reviewed
the submitted packet, not the subsequent four-group appendix.

## Final bounded extension: reviewed and included in Appendix E

At 16:24 UTC, a materially different 20-minute proof attempt began: phi=q,
freely chosen individual fertility, a young uncapped grant recipient, a young
capped taxpayer, an uncapped old taxpayer, and a capped old rebate recipient.
Current young taxes offset the grant recipient's additional births exactly,
so the next cohort size can stay fixed without fixing individual fertility.
The old tax clears housing; its funding surplus goes to capped old owners.
All young old-age resources then remain unchanged. A strict ownership-value
margin can preserve voluntary tenure choice under identity-specific cash
payments. The complete finite policy, original budget and estate accounting,
information assumptions, welfare inequality, and a nonempty positive-tax
primitive family passed an independent hostile review. The assembled appendix
was also reviewed, and the warning about taste-dependent targeting and the
ordinary logit formula was incorporated.

This establishes a special-case equilibrium utilitarian improvement toward
young housing with individual fertility and tenure freely chosen. It requires
phi=q, strong individual targeting, and offsetting fertility responses. It
does not establish an increase in aggregate births or population, and it is
not a claim about the already calibrated quantitative model. The result stays
in the appendix because its construction is more involved than the main
illustration.

## Delivery and open decisions

The separate note is 14 pages: seven pages for the environment, four main
propositions and transition discussion; seven pages of proofs and extensions.
The source is `latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex` and
the deliverable is `output/pdf/simplified_olg_utilitarian.pdf`. The amendments
README indexes all reports, the Pro capture, and the verification receipt.

Forty-seven exact symbolic/source-preservation checks passed. Two successive
final LaTeX builds completed; the last log has no warnings or overfull boxes.
Final page inspection and PDF/source hashes are recorded in
`utilitarian_checks.json`. These checks support the analytical proofs; no
equilibrium simulation or numerical neighborhood certificate was used.

The author still needs to decide the welfare weights and the transfer
authority's information and timing. Equal weights on each living household's
remaining utility are explicitly proposed, not treated as an adopted change.
The result that a policy raises fertility at every date in a market-clearing
demographic transition, and consequently raises limiting population, remains
open. The finite bundle test and population accounting are valid conditional
comparisons on existing paths; they do not supply that missing policy result.
