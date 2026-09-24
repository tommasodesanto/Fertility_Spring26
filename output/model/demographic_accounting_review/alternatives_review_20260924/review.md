**Recommendation for discussion, not adoption.** Prefer **A: one departure-to-entry rule, including release when the parental household dies**, while retaining the proposed independent (2/9) departure probability. This makes the smallest change to household dependence dynamics and avoids inventing guardianship. Do not add released dependents to the existing birth queue. However, A is defensible only as an **economic lifecycle approximation**: it cannot deliver literal chronological age-18 entry. Nor does it settle the separate conversion or top-bin problems. If exact chronological membership and entry are indispensable, the current constraints cannot all hold simultaneously.

Reference completed fertility of 2.1 remains accepted. It neither selects households per departing person nor proves replacement. The decision now is which approximation to accept, before implementation; no new calibration, model run, or specification change underlies this review.

**Objects and identities.** A person is one living individual. A household (H_t) is the accepted two-adult decision unit. (n) records children ever born, capped in the state representation; (m) records currently economically dependent children, including potentially adult offspring. It is not automatically a count of minors. Adult offspring means chronological adults descended from earlier households; (n-m) does not identify their current survival, residence, or wealth. An entrant is a newly created decision household, not a newborn or a person.

Take stocks immediately before births. Let (C_t=sum m g_t) count represented dependents; (B_t) count births in identical units; (X_t) count exiting households; (M_t) count ordinary departures among surviving households; (R_t) count dependents released from exiting households; and (D_t) count genuine dependent deaths. Set (F_t=M_t+R_t). Parent death occurs before survivor maturation, so (R_t) and (M_t) are disjoint. Terminal households belong in (X_t,R_t).

Let (c) be households per retained departing person, (
ho) the retained fraction, (I_t) outside entrant households, and (J^H_t,J^C_t) any additional household/dependent adjustments from an imposed empirical distribution. Then

\[
C_{t+1}=C_t+B_t-F_t-D_t+J^C_t,\qquad
H_{t+1}=H_t-X_t+c\rho F_t+I_t+J^H_t.
\]

If every resident belongs exactly once to either category, (P_t=2H_t+C_t), hence

\[
\Delta P_t=B_t-2X_t-D_t-(1-\rho)F_t
 +(2c-1)\rho F_t+2I_t+2J^H_t+J^C_t.
\]

Literal conservation requires (c=1/2), or an explicit disposition of the residual ((2c-1)\rho F_t). With (c=1/2.1), (1/21) of retained departing persons lacks a destination: calling this normalization does not identify mortality, emigration, or continued non-household residence. A residual resident pool would need resources and subsequent transitions. **This derivation does not authorize replacing the retained conversion.**

In a closed stationary population without dependent deaths, these identities require (B=F), (E=X=cB). Thus births per entering household must satisfy (B/E=1/c). With consistent literal fertility units, no reproductive mortality, and 2.1 births per couple, (c=1/2) gives 1.05 descendant couples per entering couple, not exact replacement. Age-dependent fertility and prices may change renewal away from that reference. A normalized stationary household distribution is useful conditional on its entry rule; it is not proof of a closed, constant-size population. A growing age distribution is also distinct from a stationary population.

**What the source currently establishes.** The stationary solver's *potential maturation-entry counter* uses literal departures times 0.5; its actual entry/scale closure must be distinguished from this counter. Transition drivers instead use adjusted births times (1/2.1), subsequently applying any retention/outside-flow/observed-age bridge. The independent queue does not remove offspring from (m) when they enter. A person can therefore remain represented as a dependent after queued entry; conversely an early departure leaves (m) before queued entry. Parent death removes (m), while the queue survives. Neither household mass checks nor matching an empirical age profile resolves these overlaps and gaps. [S1–S4]

| Alternative | What becomes coherent | Approximation or cost | Assessment |
|---|---|---|---|
| **A. Ordinary departure plus death release; one entry stream** | Every represented dependency exit has one destination; no parent-death disappearance if offspring survive. Same accounting in stationary and dated calculations. | Early geometric exits and late/death exits all start the same adult lifecycle; conversion/top-bin units still need decisions. Dependents are assumed to survive until release. | Smallest change to the retained household problem; preferred with explicit economic-age interpretation. |
| **B. One birth-age queue; (m) is only a support obligation** | An aggregate pre-entry stock (Q), with (Q'=Q+B-q) and (E=cq), can supply a separate demographic ledger. | Must count (2H+Q), never also (m). Parent housing/cost obligations can continue after the offspring has an independent household, or disappear before entry. Requires a new interpretation of costs, residence, and any transfers. | Does not answer the author's membership objection merely by renaming (m); reject as the default. Existing four waiting slots imply a 20-year effect lag, not exact 18. |
| **C. Constant departure plus parent-age cap, one entry stream** | Eliminates dependents at parental death under a compatible mortality schedule. | Changes continuation values and makes later births have shorter dependence. A lower hazard can restore a reference average, not every child's duration or chronological entry age. | Additional substantive restriction, not needed merely to close counts. Consider only if the tail itself is rejected. |

No fourth alternative is materially better without adding a missing-person, guardian, age, or history object. An aggregate residual ledger can expose a discrepancy; it cannot itself supply economic behavior. No child-age state expansion is proposed.

**What an (m)-only state can honestly implement.** Two households at the same parent age with identical (n,m) can contain children born at different dates. Their admissible age-18 departures differ. A common transition from that state cannot implement both exactly. Parent-age caps do not remove this information loss. Nor can an average entrant-age reassignment recover actual child ages under arbitrary policy-induced fertility timing.

Under the retained constant law, newborns face the first draw at the end of their birth period; there is no newborn exemption. Absent death, four-year counted dependence has mean (4/\mu=18). With a cap allowing (K(a)) periods after birth at parent age (a),

\[
d(a,\mu)=4\sum_{k=0}^{K(a)-1}(1-\mu)^k.
\]

A reference restriction (sum_a w_a d(a,\mu)=18) must use **all-birth** weights. If weights are model-generated, they and the hazard must be solved consistently; empirical fixed weights are another restriction. Freeze the selected hazard for counterfactuals. Do not repeatedly renormalize it to erase timing effects. A newborn exemption would change this formula and is a separate choice. [S1,S7]

**Example 1 — exact counts through births, maturation, and death.** Begin with two couples, no children: (H=2,C=0,P=4). Each has one birth: (B=2,P=6). On the next transition one parental couple exits with its dependent, while the other household's child departs normally: (X=1,R=1,M=1). Under literal pairing, the two released people form one couple: (H'=2,C'=0,P'=4=6-2). No child dies. With (c=1/2.1), new household mass is (20/21), so (P'=82/21=3.905): the missing (2/21) people require a destination. Leaving both births in a queue to generate households later would duplicate their entry.

**Example 2 — the same count is not the same age.** Consider a child born at parent node 42. The first ordinary draw, 42→46, can remove the child at chronological age four and feed an entrant labeled 18. If instead it remains dependent until parent node 66, it is already 24. Death release during 66→70 produces an entrant at the next date when the child would be 28, again labeled 18. The label resets earnings, fertility opportunity, and remaining lifetime; this is economically consequential. A certain draw evaluated at node 62 removes the child on 62→66, at age 24. To clear it on arrival at 62, the certain draw belongs at 58. These are distinct deadlines, neither adopted. These examples condition on stated paths; they are not cohort loss probabilities.

**Death, resources, and omitted heterogeneity.** In the historical measured checkpoint, births end at 42 and mortality first applies at node 66. Thus the remaining dependents exposed to death are chronologically adult under its timing, not orphans requiring a care model. The historical loss is \(5.972\times10^{-3}\) dependent units per four-year step and unit household mass, or 1.200% of the dependent stock; 48.985% occurs at terminal exit. Its 5.181% ratio to the birth flow is not a cohort probability. It is not the latest checkpoint. The newer frozen-parameter receipt confirms retirement survival but does not independently certify every maturation/birth setting. Verify the actual selected parameter object and adapter overrides before implementation. [S5,S8]

Household death is the accepted estate event. Counting it as two adult deaths is a joint-exit abstraction, not a derived two-spouse mortality model; an individual life table does not automatically validate it. A assumes zero offspring mortality while dependent. These approximations require disclosure, not automatic spouse or orphan states.

Inheritance is separate from entry. Warm-glow utility is a preference, not a transfer. Source supports a pooled estate payment to ages 45–65, while entrants draw exogenous wealth. Neither establishes that released offspring receive their own parents' estates. If a new inheritance allocation is ever added, enforce

\[
\text{net estates}=\text{paid to existing households}+\text{paid to entrants}
+\text{taxes/external disposition},
\]

with entry endowments' financing separately identified. Never pay the same estate through the pool and again to released children, or add it to entry wealth already interpreted as inheritance. Money conservation does not require genealogical matching; matching is an optional extension. [S6]

The top state is literally three for decisions/dependence but has reporting weight 3.602; queue births add 0.602 units at entry into that bin. Those units currently lack corresponding dependent resource/membership states. A blanket multiplier on top-bin dependents also reweights earlier births when (n) changes. An age-free extra-child pool cannot enforce parent-age caps or parent-death release. The smallest honest ledger uses literal state births/dependents and reports the top-bin measurement mapping separately; that would change current renewal and the relationship between the accepted 2.1 observer and literal reproduction. Retaining expanded births instead needs a justified auxiliary membership/resource mapping. Neither choice is adopted here. [S2,S3,S7]

**Author choices and implementation boundary.** Resolve: economic versus chronological entry age; the conversion and its population interpretation; literal versus expanded demographic units; and, only if desired, a cap/mean-duration restriction. Core implementation then has one birth/departure/death ledger shared by stationary and both transition drivers, terminal release, one entrant insertion, explicit migration/bridge residuals, and consistent initialization. Historical age conditioning may remain an explicit bridge; a zero-migration forecast must not silently inherit it or outside entrants.

After approval, the smallest checks are deterministic birth/death examples (including terminal/top-bin states), an actual-parameter contract check, and fixed-policy no-shock propagation against the chosen stationary identity. Compare ordinary/death-release flows and entry timing before any equilibrium run. A changes demographic feedback and therefore can change equilibrium prices, pensions, policies, and calibration even if its household dependence kernel is unchanged; C directly changes that kernel too. Holding units and conversion fixed, A increases departure-only stationary renewal by \(cR\); with no dependent deaths its stationary flow becomes \(cB\). Transition responses can begin at the next four-year date instead of waiting twenty years, with a dispersed lag thereafter. Their equilibrium sign and size are not established. Existing equilibria are comparison evidence, not automatically valid initial states. Only then scope the required fresh solve and standard diagnostics.

Evidence identifiers S1–S8, exact lines, hashes, established facts, and limitations are in the adjacent `receipt.json`. No raw rejected worker arithmetic is used.
