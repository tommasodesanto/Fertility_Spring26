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

## September 8 discussion and theory-slide revision

The author requests parallel development of the discussion draft and the
existing theory slides. The deadline is **Thursday, September 10**, confirmed
explicitly. The active seminar source remains
`latex/september_14_presentation.tex`; only its seven main theory frames and
contiguous theory appendix are in scope. Quantitative sections are preserved.
The author-controlled manuscript remains read-only.

Two presentation clarifications are immediate: \(T_t\) in the household
budgets is the common property-tax rebate, not the later policy grants;
\(\mathbb E_t h\) denotes a cross-sectional cohort average over endowments
and tenure, not a forecast. One explanatory sentence was added to the
discussion note's equilibrium definition. Its notation and household equations
are preserved.

The old-estate condition remains an open substantive issue. With old housing
uncapped, put \(A=p+qP\), \(K=1+\gamma+\omega_B\) and \(m=K/z\).
When \(e=Ph^o\), optimal old housing is
\[
h^o=\frac{(\gamma+\omega_B)z}{KA},
\qquad
\frac\alpha s-\frac\gamma{h^o}
=\left(\frac\beta q-1\right)pm+\mu L
-m\frac{q\gamma P-\omega_Bp}{\gamma+\omega_B}.
\]
Thus \(\beta\ge q\) and binding young finance alone do not sign the gap.
Zero financial estate is a regular parameter regime, not an exceptional
boundary. An independent Astra/max check verified both old regimes and a
matched-household limiting counterexample; no full-equilibrium counterexample
or new primitive theorem is asserted. Weak bequest tastes reduce old housing
within the binding-estate regime; strong tastes instead reduce housing within
the financial-estate regime as saving shifts to financial assets.

Liquidation costs were discussed but **not adopted**. A hypothetical loss of
fraction \(\delta\) at death changes the financial-estate threshold to
\(\omega_B/\gamma>(1-\delta)/(R_f-1+\tau^p+\delta)\).
That is not the same specification as a cost of downsizing while alive.
Death-only costs also make young and old service costs differ, and real losses
must enter planner feasibility; the previous welfare theorem cannot simply
reuse the easier threshold. The independent reviewer checked this distinction.

The updated diagrams use raw housing utility for the allocation comparison
and two conditional fertility/population panels for demographic adjustment.
The latter assumes the policy fertility ordering and applies cohort accounting;
it does not import the earlier financing model's equilibrium transition or
assert new stationary equilibria. Earlier figure files and supporting notes
are retained. The existing build driver now generates these explicit
illustrations without running an equilibrium solver.

The completed first slide revision contains seven main theory frames (PDF
pages 6–12) and nine supporting appendix frames (pages 44–52) in the existing
71-page deck. The seven-frame review extract is
`output/pdf/simplified_olg_theory_slides.pdf`; the full reader copy is
`output/pdf/september_14_presentation.pdf`. Source outside the two authorized
theory spans is byte-for-byte unchanged. All internal links resolve. Both
decks were compiled twice after final edits; all changed frames were rendered
and inspected, with no layout overflow. The full deck retains its existing
non-visible appendix-bookmark warning. The slide receipt records final hashes
and the scope of the schematic figures.

The discussion-note checker again passed all 47 exact identities and source
checks. Its clarification changes only page 3; that page was rendered and
inspected after two compilation passes. The remaining pages have identical
extracted text and layout to the previously inspected version. No PDF preview
was opened. The estate-regime issue, welfare weights, transfer timing and
information, and the equilibrium policy transition remain discussion items;
this slide revision does not close them.

The author requests that subsequent theory changes be made in both the
discussion note and the theory slides. Pure LaTeX changes go to a fast
subagent; the lead checks the scope and equations. Cross-sectional averages
use bars, with integrals where an explicit aggregation operator is needed.
Expectation notation is reserved for uncertainty. This convention applies
throughout the current note and the deck's theory appendix, as well as the
main equilibrium slide.

## September 8: full planner and stationary comparison

The author stopped further note/slide revisions to settle the planner itself.
The intended direct planner chooses the full allocation, including consumption;
the housing-only variation belongs in a proof and does not define its choice
set. The author favors holding fertility fixed for the first comparison of
competitive and planner stationary allocations, and is open to two welfare
criteria informed by dynastic OLG literature.

The handoff recommends freezing each type's competitive fertility and the same
cohort mass N. These precise restrictions are proposals for the benchmark,
not a newly implemented household problem. Stationarity requires average
fertility 1/nu but does not determine N. A stationary cohort lifetime objective
weights old utility by beta; the old note's current-living remaining-utility
sum is a different object. Neither automatically represents welfare over an
attainable transition from a common inherited state. The model has warm-glow
estate utility, not descendants' continuation utility; importing dynastic
altruism would require new household/inheritance equations.

An independent Astra/max audit verified the permanent matched-owner variation:
young housing rises by epsilon, old housing falls by epsilon, young next-period
net wealth falls by P epsilon, and old financial saving rises by q P epsilon.
Age transfers +p epsilon to young and -p epsilon to old settle the budgets at
reference prices; estates and consumption remain fixed. Entering old total
resources and aggregate bond investment are unchanged. With slack relevant
caps and a positive financial estate, the stationary lifetime-welfare
derivative is
\[
\frac{\alpha}{s}-\beta\frac{\gamma}{h^2}
=\beta(1/q-1)pm+\mu L.
\]
This is positive even with a slack mortgage when q<1. It therefore does not
isolate a mortgage-induced inefficiency, nor characterize the direction of
housing in the full optimum. No stronger theorem is promoted.

Full resource accounting remains unresolved: estate recipients, the source of
entrant wealth, rental-intermediary/title ownership, and the real definition
of bequests for a direct planner. The individual estate floor is a financial
restriction, not a physical housing cap. Free choice of endowed aggregate
external wealth would invalidate a welfare comparison. The Pro packet asks
for the smallest explicit completion, with every new assumption distinguished
from the maintained household model.

The focused prompt is `docs/prompts/oracle_simplified_olg_stationary_planner.md`;
the self-contained packet is
`output/model/simplified_olg_amendments/oracle_stationary_planner_bundle.md`.
It includes the exact current note but supersedes its welfare agenda for this
review. Primary source pointers include Becker–Barro, Golosov–Jones–Tertilt
Section 3.4, and Farhi–Werning. Transfers, fertility welfare and demographic
transition proofs are deferred. The author will paste the packet manually;
preparation does not mean a Pro run has been launched. No note, slide, PDF,
household equation or protected author draft was edited in this step.

A bounded independent Astra/max handoff review found no wrong equations or
scope contradictions. Its two accounting corrections were incorporated:
state both renter and owner estate definitions, and qualify internal property
tax transfers by domestic ownership/consolidation when rental title ownership
is still unresolved. This review validates the handoff's stated distinctions,
not the existence or solution of a completed planner problem.

## September 8: Pro response captured; first assessment only

The Pro chat at
`https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354` completed after
20m32s. Its full response is preserved as
`output/model/simplified_olg_amendments/oracle_stationary_planner_response.md`.
The browser's completed-response controls confirmed completion. Capture checks
preserved all 31 numbered equations and 128 mathematical expressions; display
formatting was reconstructed from visible text and DOM equation labels.

Pro recommends stationary cohort lifetime welfare, with a discounted sequence
of cohort lifetimes as the comparator. Its proposed accounting treats entrant
wealth as outside remittances and estates as real goods paid to outside
recipients, with domestic zero-equity rental intermediaries. It fixes aggregate
external wealth before estate remittances. These are explicit new accounting
assumptions for discussion, not decisions already adopted by the author.

Under that completion, Pro derives the full goods allocation and claims a
general stationary utilitarian welfare gap for q<1, including with slack
financial constraints. Its consumption comparison follows from
\[
\frac{1}{c_i-\chi n_i}-\frac{\beta}{c_i^2}
=\beta(q^{-1}-1)\frac{1}{c_i^2}+\mu_i>0.
\]
The lead checked this immediate first-order-condition subtraction and the
stationary accounting substitution; a full independent audit is outstanding.

For the full optimum to have strictly greater total young housing, its main
proposition additionally requires the uncapped planner candidate to respect
the actual tenure caps, all competitive old housing caps to be slack, and
\[
(\gamma+\omega_B)(1-q+q\tau^p)>q\gamma(1+q\tau^p).
\]
This has no upper bound on beta, but the cap restrictions are substantive.
The lead checked how the inequality signs the old-owner marginal housing
comparison in both estate regimes. No claim of mildness or calibration
compatibility is made. Pro also supplies a proposed heterogeneous analytical
family with less young housing under the planner when old financial estates
are zero; its equilibrium construction remains to be checked. This reinforces
the need to separate the welfare criterion's redistribution from an effect
specifically caused by young mortgage constraints.

The `watch-pro-planner-review` heartbeat was paused after capture. No follow-up
was sent to Pro and no note, slide, PDF, household equation or protected draft
was changed. Next discussion should settle the proposed accounting and welfare
criterion before adopting or extending a theorem.

## September 8: author redirects to one dated allocation

The author rejected the stationary-cohort/permanent-transfer argument as the
answer to the housing question. The active exercise is now a one-date
utilitarian comparison of current young and old households, with fertility
fixed and both current consumption and housing free. Pro's stationary
criterion and outside-estate accounting have NOT been adopted. Do not resume
that branch as though the author accepted it.

The dated benchmark proposed in discussion fixes current tenure and preserves
young continuation resources and old estate payments. Equal weights apply to
current utility. Housing and current goods are separately conserved. Individual
financial positions may adjust to preserve those future real entitlements;
individual financing restrictions can be relaxed. For an owner housing change,
use changes in young next-period wealth and old financial saving of
\(-P_{t+1}\Delta h^y\) and \(-qP_{t+1}\Delta h^o\), respectively.
Current transfers are \(\Delta c+u_t\Delta h\), which sum to zero.
This is a full current-consumption/current-housing allocation problem
conditional on the fixed objects, not a full dynamic planner.

### A proved static housing ordering with binding caps allowed

Normalize the same type measure Q to one in each age group, each of mass N.
Only a common tenure/cap distribution is needed for the planner argument;
income heterogeneity is unrestricted. Let \(\bar h_i\) be the type's cap,
identical across its two ages, and \(\kappa n_i<\bar h_i\) its child-space
requirement. Assume \(n_i>0\), \(\alpha\ge\gamma>0\), and
\[
\kappa N\int n_i\,dQ<\bar H<2N\int\bar h_i\,dQ.
\]
The upper inequality says the stock does not exhaust every household's cap.
Some individual caps may bind in either allocation.

With fixed fertility, utility is additively separable between goods and space,
so optimizing consumption jointly does not change the housing solution. For
the housing-resource multiplier \(\lambda>0\), that solution is
\[
h_i^{y*}=\min\{\bar h_i,\kappa n_i+\alpha/\lambda\},\qquad
h_i^{o*}=\min\{\bar h_i,\gamma/\lambda\}.
\]
The common multiplier clears the housing market. Therefore
\(h_i^{y*}\ge h_i^{o*}\), strictly wherever the old household is below its cap.
If all old households were at their caps, all young households would also be
at their caps, contradicting the strict aggregate capacity inequality. Hence
\[
H_Y^*>\bar H/2.
\]
Consequently the following is a sufficient, interpretable cross-sectional
condition for the FULL dated optimum to increase total young housing:
\[
H_Y^{\rm eq}\le H_O^{\rm eq}
\quad\Longrightarrow\quad H_Y^*>H_Y^{\rm eq}.
\]
It compares total housing with equal cohort masses, or average home sizes.
It is a theorem conditional on the observed equilibrium ordering, not yet a
generic implication of mortgage constraints. The assumption alpha>=gamma is
explicit: young adults value housing at least as strongly as old adults,
before accounting for their children's space. No high-estate-taste condition,
bound on beta, or universally slack housing caps is needed for this planner
ordering. It establishes utilitarian improvement, not a Pareto improvement.

### Primitive sufficient conditions that deliver the market ordering

A separate Astra/max derivation, checked by the lead against both old estate
regimes, gives a conservative primitive route in a stationary competitive
reference. Stationarity is used only to match the current young/old type and
tenure distributions; welfare remains entirely one-date.

Let \(d_p=1-q+q\tau^p\), \(a_p=1+q\tau^p\),
\(d_L=1-\phi+q\tau^p\), \(K=1+\gamma+\omega_B\), and
\(J=\gamma+\omega_B\). Old housing is exactly
\[
h_o^d=\min\{h_d^{\max},g_d z/P\},\qquad
g_R=\frac{\gamma}{K d_p},\qquad
g_O=\frac1K\min\left\{\frac{\gamma}{d_p},\frac{J}{a_p}\right\}.
\]
The owner coefficient covers both financial-estate regimes. Clipping the
uncapped optimum at the physical cap is valid by concavity even if that
changes the estate regime at the capped choice.

For owners, the original cash constraint gives \(Ph_y<w/d_L\), and mortgage
repayment gives \(z\ge v+(1-\phi/q)Ph_y\). Since \(g_O<1\), the coefficient
\(g_O^{-1}-1+\phi/q\) is positive. It follows that
\[
\frac vw\ge A_O\equiv\frac{g_O^{-1}-1+\phi/q}{d_L}
\quad\Longrightarrow\quad g_Oz/P>h_y.
\]
Clipping old housing at the shared cap preserves \(h_o\ge h_y\).
For renters, \(Ph_y<w/d_p\), \(z\ge v\), and the corresponding sufficient
threshold is \(A_R=K/\gamma\). Thus requiring \(v_i/w_i\ge A\), where
\(A=\max\{A_R,A_O\}>1\), for all endowment types gives the market ordering
in both tenures. Logistic ownership tastes are fully retained.

In the explicitly identified zero-tax subcase, prices and rebates disappear:
\[
\frac{y_i^o}{y_i^y+b_i}\ge
\max\left\{\frac K\gamma,
\frac{\max\{K(1-q)/\gamma,K/J\}-1+\phi/q}{1-\phi}\right\}.
\]
This is a sufficient restriction on primitive income/wealth and parameters,
not a numerical-neighborhood proof. It can demand high old income relative
to young cash; it has not been shown mild or quantitatively relevant.

For retained positive property tax, the current note's existing bound
\[
\bar T=\frac{\tau^p\int(y^y+b+qy^o)\,dF}{(1-q)(2-\tau^p)}
\]
is valid for \(0\le\tau^p<2\). A fully primitive sufficient restriction is
\[
y_i^o\ge A(y_i^y+b_i)+(A-1)\bar T.
\]
Indeed it implies \((y_i^o+T)/(y_i^y+b_i+T)\ge A\) for every
\(T\le\bar T\). The lead checked the bound from the lifetime budget and
rebate identity. It can be stringent or even vacuous at large taxes; it must
not be presented as necessary. These restrictions do not independently
establish a binding young multiplier or isolate its causal contribution.

### Verification limit that must survive future rewrites

The cap audit also disproves the inference that paired marginal housing
utility gaps always imply MORE aggregate young housing in the full optimum.
As a purely allocation-level counterexample, take two types of unit mass at
each age, alpha=gamma=1, common cap4, child-space needs (3,1), young housing
(3.9,3.5), old housing (1,2.6). Each paired young marginal housing utility is
larger. With total housing11, the static optimum is young(4,3), old(2,2), so
young aggregate housing falls from7.4 to7. This example is not asserted to be
a competitive equilibrium of the household model. It demonstrates why a
local improving variation is insufficient to establish the full aggregate
direction when heterogeneous households meet housing caps.

Two independent bounded Astra/max checks supplied the cap proof and the
competitive bounds. The lead verified the formulas, strictness and financial
settlement. No model solve, paper/slide revision or new welfare convention was
implemented. The open task is to sharpen the conservative equilibrium
conditions and decide whether this static statement captures the author's
intended claim; do not replace it with the rejected permanent-transfer result.
