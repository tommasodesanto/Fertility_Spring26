# Resolve the illustrative housing–fertility theory

You are reviewing and completing the analytical theory section of a paper about housing, ownership, and fertility. Work as an independent economic theorist. We need a short, transparent argument supported by correct proofs, not a long catalogue of sufficient conditions. The intended contribution is simple: young households need space for children, face mortgage limits and limited rental sizes, and may occupy too little housing; improving housing access may increase fertility during a demographic transition and change the eventual population level.

Several rounds of work have produced useful results but repeatedly answered slightly different welfare questions. The definitions below now take priority. Do not change them to obtain a desired sign. The author's latest additional request is a **separate general-preferences investigation**, starting from an unspecified utility function, with no logarithms or linear child-needs shifts imposed at the outset. Complete the existing model first where possible, and use the separate investigation to identify what is structural and what depends on the functional form.

Use the maximum substantive mathematical effort available in this run. The author can leave this running overnight. Do not stop after discovering the first sufficient condition: try to weaken it, test its boundaries, and complete the proof and economic interpretation. This is a request for sustained work, not for a claimed number of hours. Do not claim to have run tools, consulted sources, or proved results unless you actually have. Give checkable proofs, counterexamples and derivations; there is no need to narrate private deliberation. The outcome should be a usable theory section or an exact account of the remaining obstruction, with the strongest valid result delivered in either case.

## 1. Authority and scope

Read the attached files before answering. Their order of authority is:

1. This prompt: the author's current question, welfare definition, and restrictions.
2. `oracle_consolidated_theory_context.md`: the exact maintained household environment, a curated decision history, and the status of checked arguments. Its original-source excerpts are explicitly identified.
3. `simplified_olg_general_preferences.md`: the corrected, separate exploration. Its claims are still to be checked independently.
4. The two earlier Pro responses: advisory mathematical material, including counterexamples. They are not adopted specifications or certified proofs. Their qualifications and notation do not override this prompt.
5. The writing and presentation guides: conventions for the final exposition, not mathematical assumptions.

The earlier chronological work record, quantitative calibration, and manuscript are not being supplied wholesale because they include different questions and superseded models. The context digest preserves the relevant decisions instead. Do not infer that the quantitative model is literally this two-period model. The former has a richer lifecycle, income, and housing architecture; this exercise should explain a mechanism and guide interpretation, not duplicate its numerical work.

Preserve income–wealth heterogeneity, the two-age structure, endogenous tenure when young and its persistence into old age, warm-glow estates, mortgage limits, and the two physical housing-size limits. Current income **is available for the down payment**. Old households can sell or resize their houses, and do not face a cap equal to their inherited house. They cannot take out a new loan in the competitive economy. Distinguish that restriction from the planner's powers below. Do not silently add transaction costs, altruism toward children, endogenous inherited entry wealth, another age group, a new rental technology, or a new policy instrument.

Use small superscripts (y,o) for age in the proposed exposition. (q=1/R_f) is a bond price, not a housing parameter. The discount factor and return must refer to the same model period. Use the existing symbols where possible. No expectation operator for deterministic future prices or merely to denote a cross-sectional average. Define every newly introduced object before using it. Do not create an alphabet of auxiliary bounds unless they actually shorten the final proposition.

## 2. The welfare question is fixed

Start from a **positive stationary competitive equilibrium**, so that each current age group has mass (N) and the same distribution (Q) of matched endowment types and retained tenure. Stationarity supplies a convenient reference allocation. The welfare comparison is **at one date**, not a ranking of stationary cohort lifetimes that omits an initial old generation.

Give each currently living household weight one on its remaining utility, on the maintained cardinal scales. The young household's own continuation utility carries its private discount factor \(\beta\). It is constant in the first dated comparison because its future real opportunities are held fixed. The old household's fixed net estate is also a constant in the objective. Ownership tastes are fixed when tenure is fixed. There is no additional social utility assigned to unborn people.

First fix each young household's fertility (n_i), tenure, current cohort masses, incumbent future real continuation opportunities, and old net estates (e_i). Let the planner choose **all current consumption and housing** subject to total current goods (C^{eq}), housing stock \(\bar H\), and the retained physical caps (H_R<H_O). It may redistribute resources and relax both the young financial constraint and the old nonnegative-financial-saving restriction, while honoring existing claims. This changes financial positions as needed; it does not create external resources or forgive inherited obligations. Audit the supplied settlement equations carefully, including rental intermediaries.

For the maintained log model the nonconstant objective is
\[
\max N\int\left[\log(c_i^y-\chi n_i)+\alpha\log(h_i^y-\kappa n_i)
+\log c_i^o+\gamma\log h_i^o\right]dQ,
\]
subject to
\[
N\int(c_i^y+c_i^o)dQ=C^{eq},\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,\qquad
h_i^y,h_i^o\le H_{d_i},
\]
and positive utility arguments. This is the full current allocation problem. A small housing transfer holding consumption fixed may help prove a result but is not a replacement for this optimization problem.

Then consider the same dated planner also choosing fertility through the utility of current parents. Child goods and space must be counted once. This is a static parental-welfare comparison; changing births cannot simultaneously leave every future population unchanged. Keep that limit explicit when moving from this comparison to a demographic equilibrium.

An equally weighted utility sum need not be maximized by a frictionless equilibrium with unequal endowments. Therefore distinguish redistribution inherent in this welfare criterion from an additional effect of financial constraints. Do not call every utilitarian difference a Pareto inefficiency, a failure of the first welfare theorem, or constrained inefficiency. Identify precisely which claim the final theorem warrants. We want the young-housing direction and its economic source, not just the truism that equal-weight utilitarian redistribution can raise welfare.

## 3. Main task A: a usable housing-allocation theorem

Seek the weakest transparent conditions under which the full dated planner gives **more aggregate housing to the young** than the competitive equilibrium. Give an individual characterization too, without claiming that every young household gains unless proved. Keep separate:

- a local improving direction;
- a comparison of a given household with its planner allocation;
- a comparison of aggregate young housing at the full optimum.

The preferred final statement has the form: “Under [a short set of economically interpretable restrictions on preferences, endowments, finance and housing capacity], any positive stationary competitive equilibrium has [specified misallocation], and the planner increases young housing.” The proof must derive the relevant constrained group and equilibrium comparison, not assume that young marginal housing utility is already larger. Multiplier conditions and an exact allocation test can be useful intermediate results, but they are not the entire answer.

Do not impose either \(\beta R_f\ge1\) or \(\beta R_f<1\) by default. The author is concerned about obtaining the result only on one side of that threshold. A natural combined restriction involving patience, borrowing limits, old-age income, housing preferences and estate motives may be necessary; derive it and explain the economics. First test whether the desired claim can hold generally. If not, give a simple analytical counterexample and then a genuine sufficient condition. Do not use a numerical reference point plus continuity as the main “open set” proof. Do not remove ex ante heterogeneity to make the proof easy.

We already have an explicit planner solution in the log case and partial conditions using \(\beta\Gamma/q\), as described in the digest. Earlier results required old total housing to exceed young total housing. The weaker comparison involving young **adult space** may suffice; investigate this first. The actual age profile of total housing remains an interesting separate implication. Do not claim it follows from age alone in this two-period model.

Housing caps are physical maxima for rental and owner homes. Identify separately which competitive choices and which planner choices must be below their caps. A fully capped young household cannot receive strictly more housing without relaxing its cap, which is not allowed in the baseline. Avoid requiring caps to be so large that renting becomes irrelevant without acknowledging that cost. Test whether useful aggregate results survive binding caps and genuine tenure heterogeneity.

Assess the old estate restriction and both of its regimes. Explain how estate preferences can affect the age comparison even though old households can sell their houses. Do not describe all old households as mechanically locked into their inherited houses. If an estate-floor assumption or a no-borrowing restriction defeats the desired conclusion, locate the exact term and assess whether the weakest useful result can avoid it.

For any final primitive conditions, check that they are jointly compatible with a positive equilibrium, a positive mass that **strictly** wants more borrowing, and the stated cap regime. A borrowing equality with a zero multiplier is not enough for strictness. Prefer an analytical nonempty parameter region and economic comparative statics to six layers of conservative min/max bounds.

## 4. Main task B: fertility at the dated comparison

Develop three distinct statements, in this order:

1. Given a changed consumption–housing bundle, when does privately chosen fertility increase? Derive the marginal and, where possible, finite-change condition. More housing may cost consumption, so the full tradeoff matters.
2. Under the allocation change actually proved in task A, what can be said about private fertility when households rechoose it? Do not substitute an arbitrary increase in (c,h) for the planner's or policy's actual changes. Explain any difference between a fixed-fertility allocation followed by private adjustment and a joint optimum.
3. When the dated planner jointly chooses consumption, housing and fertility through parents' utility only, does average fertility increase? Derive the comparison with heterogeneous reference households, retaining the exact cap assumptions. Validate the aggregation step rather than replacing heterogeneous households by their means without proof.

The existing responses contain an interior fertility derivative and a joint-planner scalar equation. These are useful ingredients, not substitutes for the complete comparison. Verify the corrected Cauchy–Schwarz step in the digest. Aim for a proposition whose qualifications can be explained in a few sentences.

## 5. Main task C: a constrained reform and the demographic transition

The theory must connect two ideas. First, an exogenous decline in the taste for fertility starts a transition from an old steady state toward another. We do **not** seek to explain that original decline. Second, a housing reform introduced along this transition can alter fertility along the path and the eventual population level. A pair of unrelated steady-state comparisons or a simulated fertility path imposed by hand does not establish this mechanism.

Use a permanent increase in the existing property-tax rate, with revenue rebated equally to all young and old household decision units, as the first market-based candidate to analyze. It is a candidate whose signs need proof, not an adopted claim of success. Begin both paths with the same inherited populations, asset positions, housing titles and obligations at an intervention date (t_p) along the baseline transition. For a definite exercise, treat the reform as unexpected at (t_p), then deterministic and known thereafter. If a different announcement convention changes the result, show it separately instead of comparing inconsistent histories.

The competitive households retain their original borrowing restrictions. They reoptimize. Prices, rebates, tenure, future saving, and estates can change. The policy authority does not acquire the direct planner's credit powers, and the dated planner's fixed-price settlement does not prove implementation by a tax. Specify the authority, feasible instrument and welfare criterion before any constrained-welfare claim. For a welfare assessment at the intervention date, keep equal weights on the living households' remaining utilities, include the young continuation consequences and old warm-glow estate consequences, and distinguish these from the separately reported outcomes of later cohorts. Do not silently freeze future prices or add future-cohort welfare weights.

Establish as much as possible analytically: the housing and fertility response, fiscal balance and market clearing, and the terminal population effect. If the tax cannot deliver the claim generally, give the precise obstruction and a clean sufficient condition. A different instrument may be proposed as an explicitly unadopted alternative, but do not replace the candidate or the model without saying so. A result for unrestricted transfers is not automatically a result for this tax-and-rebate rule.

Demography is
\[
Y_{t+1}=\nu\bar n_tY_t,\qquad O_{t+1}=Y_t.
\]
For common initial young mass,
\[
\frac{Y_T^{P}}{Y_T^{B}}=\prod_{t=t_p}^{T-1}
\frac{\bar n_t^{P}}{\bar n_t^{B}}.
\]
If both paths converge to finite positive stationary cohort masses, both terminal fertilities equal (1/\nu). Their population levels can differ because of cumulative fertility differences along the path. Prove the required sign from the equilibrium mechanism, or state exactly what remains conditional. A positive cumulative log gap need not require higher fertility at every date. Distinguish permanent reforms from temporary reforms under a unique terminal attractor.

These variables count adult households. Report (Y+O) and stationary (2N) accordingly. Do not label them total resident persons without specifying the conversion, including children. Transition existence and convergence are mathematical obligations, not implications of the accounting identity. If a full global theorem is out of reach, prioritize a rigorous local transition result under explicit stability and nondegeneracy conditions, derive the relevant derivatives, and identify any assumptions still involving equilibrium objects. Do not call a simulation or an assumed stability condition a proof of global convergence.

## 6. Side task, kept separate: general preferences without imposed shifts

The author explicitly rejected starting the general branch with
\(f(c-\chi n)+\alpha g(h-\kappa n)+v(n)\).
That already imposes the child-needs structure we are trying to investigate. Start with unspecified gross-bundle utility
\[
U^y(c,h,n),\qquad U^o(c^o,h^o,e),
\]
and minimal, clearly stated assumptions. A useful intermediate class the author suggested is
\[
U^y(c,h,n)=a(c,n)+b(h,n)+v(n).
\]
The separate function names distinguish the goods and housing components; they are not a substantive notation decision for the paper. Specify monotonicity in (c,h), concavity, domain, fertility tradeoffs and any cross-partial restrictions needed. Do not require total utility to be globally increasing in (n) at a fixed gross bundle while omitting every fertility cost or upper bound: that need not give an attained interior optimum. Child costs may enter the unspecified interactions; do not silently reinstall linear offsets through a domain restriction.

Determine: (i) which competitive wedges and planner conditions survive; (ii) what additional age-comparison and complementarity restrictions suffice for the relevant housing and fertility directions; (iii) which results need separability, common cardinal housing preferences or aggregation restrictions; and (iv) whether the existing shifted-log specification satisfies these conditions and where it fails. Provide a counterexample when mere increasing concavity is insufficient. Map back to Stone–Geary-type utility only after the general analysis. Preserve the maintained budgets and planner powers throughout this branch.

This is an exploration, not authority to rewrite the main model. If it reveals a much simpler theorem, explain exactly which assumptions produce that simplification and what it says about the existing form. Neither a planner-only allocation ordering nor a cross-partial fertility derivative by itself proves the equilibrium-to-planner comparison.

## 7. Existence, uniqueness, literature, and proposed extensions

Independently check the supplied stationary existence certificate and its compatibility with your proposed inequalities. Distinguish unique household choices conditional on tenure, a unique planner allocation, a unique stationary equilibrium within a specified regime, global equilibrium uniqueness, and uniqueness/convergence of a deterministic transition. Do not import the infinite-lived precautionary-saving restriction on \(\beta R_f\) into this finite-lived warm-glow model without deriving its relevance.

Use primary literature where it resolves a real conceptual issue. Relevant starting points are Coven, Golder, Gupta and Ndiaye on property taxes and intergenerational housing; van Doornik, Fazio, Ramadorai and Skrastins on housing and fertility; and the fertility/OLG work of Doepke and coauthors. Links and a conservative account of what has been checked are in the context. Verify actual statements before attributing an efficiency theorem. The goal is to use familiar ingredients, not invent a welfare criterion or claim another paper proves ours.

If two periods, frictionless resizing, or our utility form prevents an important result, demonstrate the obstruction before recommending a change. Transaction/sale costs, a third age, and larger rental homes are possible extensions only. Real moving or conversion costs must enter physical feasibility, including the planner's problem. Rank any proposed amendment by how much of the original setup and proven work it preserves. Do not turn this illustrative exercise into a new quantitative model.

## 8. Deliverable and completion standard

Return a finished analytical assessment and a proposed **theory section**, not a claim to have completed the empirical or quantitative paper. Organize the answer so the author can read a short conclusion first and inspect the proofs afterward:

1. A short verdict: the strongest housing, fertility, and transition statements actually established; the indispensable qualifications; and the best path forward. Distinguish proved, conditional, and unresolved claims.
2. A coherent LaTeX-ready main text, approximately 5–8 pages if typeset, with environment, household problems, equilibrium, the dated planner, and a small number of propositions in the requested order. If the full chain cannot be proved, write the strongest honest main text and state the missing step exactly. Do not bury it in the appendix.
3. A proof appendix with complete derivations, boundary/cap/estate cases, analytical counterexamples where needed, and compatibility checks for the conditions. Auxiliary results may be longer than the main text, but remove repetitive variants.
4. A separate concise general-preferences result and mapping back to the existing form.
5. A 5–7-slide outline focused entirely on theory. Keep two illustrations: a housing-misallocation movement, and a two-panel transition showing the initial fertility decline, the baseline destination, and a policy introduced along that same path leading to another destination. Specify axes, curves, timing and which arrows follow from proved equilibrium results. Do not relabel a chosen illustrative path as a solved transition.
6. A brief unresolved-decisions list containing only choices that actually require the author. If a desired theorem fails, supply the obstruction and the smallest promising alternative rather than a menu of unrelated models.

Before finalizing, conduct three distinct checks: reconstruct the proposed derivations from the budgets; challenge the result at equality, borrowing, cap and estate boundaries and under heterogeneity; then check economic interpretation, resource accounting and transition timing. Correct problems found in these checks. Numerical experiments may diagnose a claim but cannot replace the requested analytical proof. If response limits bind, preserve the principal theorem, its proof, and the exact unresolved transition issue rather than spending space repeating the setup.

Write simply and precisely, as in a clear economics problem set that could support a paper. The author values the direct exposition of Guido Menzio and Raquel Fernández: environment first, defined primitives, equations with a short economic explanation. Avoid unnecessary technical-sounding labels, claims of grandeur, and sentences explaining what the model does not contain unless an actual ambiguity requires it. We need a simple result that is true, with the qualifications visible.

### File: output/model/simplified_olg_amendments/oracle_consolidated_theory_context.md
````md
# Consolidated context: current model, author decisions, and checked work

Prepared for a new Pro review, September 8, 2026. This digest is research context, not a claim that the theory is complete. The accompanying consolidated prompt has priority. The exact source excerpts at the end describe the maintained household model and an existing sufficient existence argument; they do not reinstate any superseded planner.

## A. Decisions that must survive the handoff

The author wants a simple illustrative theory of housing allocation, ownership and fertility. The intended message is that finance and limited rental sizes can prevent young families from occupying the housing a planner would assign them. The final exposition should resemble a transparent problem set with a few extensions. It should retain the notation and substance the author chose, and avoid technical language that adds no content.

Earlier work pursued compensated/Pareto reallocations and constrained inefficiency. Those results are now background or possible appendix material. The present main route is equally weighted utilitarian welfare. A stationary lifetime comparison that improves later cohorts while sacrificing initial old households was explicitly rejected as answering a different question. The current comparison is a full allocation at one date, optimizing current consumption as well as housing. Do not switch back to a housing-only objective or assume future-cohort compensation solves the current question.

The planner's accepted direct powers are broad enough to relax private finance, including the old financial-saving floor, while honoring obligations and preserving fixed continuation opportunities and net estates. Old competitive households themselves still cannot borrow. They may sell and resize their houses. The author was repeatedly confused by descriptions that made old households mechanically unable to sell; that is not the model. The retained-home size is not a physical constraint.

Tenure is fixed for the dated comparison and persists across the two ages. “Same tenure distribution” means the matched stationary cohorts have the same distribution of owner/renter status, not that young and old occupy the same amount of housing. A statement about a young individual, a matched pair, a positive-mass group, and the entire young cohort are different statements. We would like both an individual characterization and a useful aggregate theorem.

The author initially expected old households to occupy larger total homes. That remains a prediction to investigate. It is acceptable to use a weaker comparison with the young household's adult space if that is enough for the housing-allocation result. Do not assume the desired equilibrium ordering or a marginal-utility gap as the entire theorem. Derive conditions from preferences, endowments, credit and physical constraints. Multiplier conditions and numerical reference-point/continuity arguments were rejected as the main final proposition.

The author is concerned about a standalone restriction on \(\beta R_f\) in either direction. We have not adopted one. Recent inequalities using \(\beta\ge q\) are partial results to improve, not instructions to impose that condition. High annual discount factors in the richer quantitative model do not automatically map into the two-period parameter without matching horizons.

The demographic illustration must start with an exogenous fertility-taste decline that initiates a transition from one steady state to another. Housing policy is introduced along this transition and may lead to a different eventual population level. We are not trying to explain the initial fertility shock. Both new positive stationary endpoints have replacement fertility; the difference in population comes from fertility along the path. Do not replace this with a static policy comparison or two arbitrarily selected fertility paths.

Transaction costs, larger rental homes, another lifecycle stage, and an alternative policy instrument are not adopted. They may be useful if a precise obstruction is established. The current baseline already includes old-age income and current earnings available for the down payment. Entry wealth is exogenous and the estate is warm glow; do not accidentally impose an endogenous dynastic inheritance law.

The separate general-preferences branch has just been corrected at the author's request. Begin with gross \(U^y(c,h,n)\), optionally \(a(c,n)+b(h,n)+v(n)\), without linear offsets. The shifted-log model is a specialization to test afterward. Both the main model and this exploratory branch must remain clearly identifiable.

## B. Useful formulas and their status

All identities below should be independently verified against the original budgets. They summarize existing progress, not an instruction to accept every earlier auxiliary bound.

### Stationary reduction

Write \(d_p=1-q+q\tau^p\), \(p=d_pP\), \(L_R=p\), and \(L_O=(1-\phi+q\tau^p)P\). Current cash is \(w=y^y+b+T\); old income including its rebate is \(v=y^o+T\). Conditional on tenure, use total resources on entering old age \(z\). The young restrictions reduce to
\[
c+ph+qz=w+qv,\qquad c+L_dh\le w.
\]
Let \(K=1+\gamma+\omega_B\), \(E=1+\alpha+\vartheta\), and \(D=E+\beta K\). Define
\[
\Gamma_R=\gamma,\qquad
\Gamma_O=\Gamma=\min\left\{\gamma,
\frac{(\gamma+\omega_B)(1-q+q\tau^p)}{1+q\tau^p}\right\}.
\]
When old housing is uncapped, \(c^o=z/K\) and \(ph^o=\Gamma_dc^o\), including either old-owner estate regime. A sufficient primitive restriction making the financial-estate floor slack is
\[
\omega_B(1-q+q\tau^p)>q\gamma.
\]
Its economic role must be explained; it is not a condition for the old to be able to sell.

### Full dated planner

Let \(x_i=c_i^y-\chi n_i\), \(s_i=h_i^y-\kappa n_i\). An overbar is a mean under the common stationary law \(Q\), not an expectation about aggregate uncertainty. At fixed fertility the planner's adult consumption is common:
\[
x^F=c^{o,F}=(\bar x+\bar c^o)/2.
\]
Its housing choices are
\[
h_i^{y,F}=\min\{H_{d_i},\kappa n_i+\alpha/\lambda\},
\qquad h_i^{o,F}=\min\{H_{d_i},\gamma/\lambda\},
\]
where \(\lambda\) clears housing. This solves the full consumption–housing problem because the fixed-fertility log objective separates. That separation need not hold for arbitrary utility.

For owners the proposed settlement uses
\[
\Delta a'_i=-P_{t+1}\Delta h_i^y,\qquad
\Delta a_j^e=-qP_{t+1}\Delta h_j^o,
\qquad t_k=\Delta c_k+u_t\Delta h_k.
\]
It preserves young future net resources and old net estates. Rental allocations require the corresponding intermediary accounting. With total goods and housing fixed, transfers sum to zero. These are dated settlement identities; they are not an equilibrium tax-policy implementation or permission for competitive old households to borrow.

### Partial equilibrium-to-planner comparisons

The latest attached Pro response used the strong condition
\[
\alpha\ge\gamma,\qquad \beta\Gamma/q\ge\alpha+\vartheta
\]
to put old total housing above young total housing and derive further results. We then investigated weaker adult-space comparisons. Under appropriate slack caps, \(\beta\Gamma/q>\gamma\) is sufficient for the resource comparisons used by the aggregate housing and joint-planner fertility argument. With \(\alpha\ge\gamma\), binding competitive caps can be handled in some comparisons; binding **planner** caps remain an important unresolved generalization. Do not collapse these different cap requirements.

If the old financial-estate floor is everywhere slack, \(\Gamma=\gamma\). At \(\beta=q\), the uncapped-old first-order conditions imply
\[
\alpha h_i^o-\gamma s_i
=(\mu_iL_{d_i}+\eta_i^y)s_i h_i^o\ge0,
\qquad c_i^o=x_i(1+\mu_i/\Lambda_i).
\]
Here \(\Lambda_i\) is the young lifetime-budget multiplier, \(\mu_i\) the young financing multiplier, and \(\eta_i^y\) the young housing-cap multiplier. With \(\alpha\ge\gamma\) and the relevant planner caps slack, positive mass with \(\mu_i>0\) makes the aggregate housing and joint-fertility comparisons strict. These are useful partial results; the goal is to avoid using patience relative to \(q\) as a standalone final assumption if possible.

In the explicit zero-tax, \(\phi=q\), market-cap-slack, positive-financial-estate subcase,
\[
x_i=\min\{w_i/E,(w_i+qv_i)/D\},\qquad
p=(\nu\vartheta\bar x-\chi)/\kappa,
\qquad n_i=\vartheta x_i/(\chi+p\kappa).
\]
At \(\beta=q\), \(F\{v>(K/E)w\}>0\) gives positive strictly constrained mass. This demonstrates a tractable primitive income–cash condition, but \(\phi=q\) is a special benchmark, not an adopted restriction on the main model. If all financial and physical constraints are slack at \(\beta=q\), aggregate housing and average fertility can be unchanged at the planner optimum even with heterogeneous endowments. Individual redistribution may nevertheless improve utilitarian welfare.

The earlier dated Pro response contains a counterfamily with binding young finance but a reversed aggregate housing direction when the old estate regime differs. Preserve and check it before asserting that young borrowing constraints alone suffice. It also contains conservative primitive inequalities; the author found those too elaborate for the main statement. They may inform a proof without becoming the presentation.

### Fertility identities and the corrected aggregation step

For the maintained log utility, a private interior fertility optimum at a given gross bundle satisfies
\[
\frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s}.
\]
Its conditional response is
\[
dn=\frac{\chi\,dc/x^2+\alpha\kappa\,dh/s^2}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
\]
For housing rising by \(\Delta h\ge0\) and consumption falling by \(\delta\ge0\), the finite bundle comparison, evaluated at the initial fertility and positive adult bundle, gives
\[
n_1\ge n_0\quad\Longleftrightarrow\quad
\delta\le\frac{\alpha\kappa x^2\Delta h}
{\chi s(s+\Delta h)+\alpha\kappa x\Delta h},
\]
with the usual feasibility/interiority conditions. This is not a funded policy theorem.

If the joint dated planner is uncapped, it chooses common fertility \(n^J\) satisfying
\[
\frac{\vartheta}{n^J}=
\frac{2\chi}{\bar c-\chi n^J}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa n^J},
\qquad \bar c=C^{eq}/N,\quad \bar h=\bar H/N.
\]
The heterogeneity step uses Cauchy–Schwarz after multiplying the reference household fertility equation by **\(n_i^2\)**:
\[
\vartheta\bar n=
\chi\int n_i^2/x_i\,dQ+\alpha\kappa\int n_i^2/s_i\,dQ
\ge\bar n^2(\chi/\bar x+\alpha\kappa/\bar s).
\]
The earlier response said multiply by \(n_i\); the squared version is the corrected argument. Additional comparisons of old versus young mean resources are still needed to sign the joint optimum's fertility gain. The relevant cap assumptions must be checked for the **joint** optimum, not just the reference or fixed-fertility optimum.

### General-preferences identity

For the same stationary reduced budgets and regular interior goods choices, write \(m_i=V_d'(z_i)=U_c^o\), and put \(\rho_i\ge0\) on the old-owner estate floor \(e-Ph^o\ge0\), zero for renters. Then
\[
U_h^y-U_h^o=
(\beta/q-1)p m_i+L_d\mu_i+\eta_i^y-P\rho_i-\eta_i^o.
\]
No logarithm is required. The identity does not prove the full planner's aggregate direction, especially with nonseparable consumption–housing utility. It also shows why eliminating logarithms alone does not eliminate every appearance of \(\beta/q\).

### Existence and remaining transition work

The existing sufficient stationary existence proof is reproduced below. Its low-price fertility condition contains \(D=E+\beta K\), so it has a finite upper bound on \(\beta\) at fixed other primitives. This is a conservative certificate, not a necessary existence condition and not the infinite-lived precautionary-saving restriction.

Strict concavity gives conditional household real-choice uniqueness and a unique dated planner real allocation where an optimum exists. The explicit \(\phi=q\), zero-tax, cap-slack benchmark gives one positive stationary price within that regime when \(\nu\vartheta\bar x>\chi\). This does not exclude additional equilibria involving capped households. General stationary uniqueness and deterministic transition existence, uniqueness and convergence have not been established. None of the recent static results proves a funded general-equilibrium property-tax transition. Older transition arguments used earlier household specifications and should not be silently imported.

## C. Literature starting points and visual structure

These links were checked in the preceding work; verify the actual statements if using them in the new result.

- Coven, Golder, Gupta and Ndiaye, August 1, 2026 version: <https://abdouecon.github.io/research/papers/Property_Tax.pdf>. The simple two-period section establishes capitalization and intergenerational redistribution with consumption utility. Its quantitative lifecycle model includes housing, finance, transaction costs and estates. It does not supply a theorem that our equally weighted planner gives more housing to young households. Welfare and allocation effects of the reform need separate assessment.
- van Doornik, Fazio, Ramadorai and Skrastins, Housing and Fertility, December 2024 version: <https://bcb.gov.br/content/publicacoes/WorkingPaperSeries/WP612.pdf>. Goods and space are useful fertility ingredients. The simple model does not provide our endogenous tenure/mortgage architecture. Verify exact formulations before borrowing a result.
- Aiyagari (1994), especially pp. 668–670: <https://www.liuyanecon.com/wp-content/uploads/Aiyagari-1994.pdf>. Its stationary precautionary-saving condition concerns an infinite-lived income-risk environment. It should not be applied as a theorem about this finite-lived model with exogenous entry wealth.

The author wants two familiar illustrations and only 5–7 theory slides. First show housing moving from an old household to a young household, with labeled market and planner allocations and housing marginal-utility curves that match the welfare criterion. An older compensated-Pareto figure is a visual reference, not the current welfare theorem; do not carry over compensation claims automatically.

The second illustration combines a fertility decline and an intervention during the resulting transition. It should show an initial steady state \(S_-\), a baseline destination \(S_0\), and a reform introduced at a common inherited state that leads toward \(S_1\). Two panels should make the fertility and population/equilibrium movements readable together. Curves, axes and arrows must correspond to the proven mechanism; do not represent all transition points as lying on a stationary schedule. The existing deck's recent curves were acknowledged to assume fertility paths. The author explicitly rejected treating those as a solved equilibrium transition. Numerical illustrations belong mainly in the quantitative model, so an analytical or clearly labeled schematic construction is preferable here.

## D. Maintained household model: source excerpt

The following is extracted from `latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex`, from Environment through the stationary household reduction. Only old-age superscripts have been harmonized from the stale `2` to the author's requested `o`. This is the household model, **not** the new full dated planner definition, which is supplied by the consolidated prompt above. Allowing the existing property-tax parameter to be date-dependent for a permanent reform is the explicit policy comparison requested in that prompt.

```latex
\section{Environment}

\textbf{Households.} At date $t$ there are $Y_t$ young and $O_t$ old households.
A young household has liquid wealth $b_i>0$, income $y_i^y>0$ when young, and
income $y_i^o\ge0$ when old. The triple $(y_i^y,b_i,y_i^o)$ has distribution $F$
in each entering cohort. Future income is known. Households choose completed
fertility once and retain their tenure when old.

\textbf{Preferences.} Each child uses $\chi>0$ units of goods and $\kappa>0$
units of space. Total nondurable expenditure is $c$ and housing is $h$; the
bundles left for adults are:
\begin{equation}
x=c-\chi n>0,\qquad s=h-\kappa n>0.
\end{equation}
Young and old utility are:
\begin{align}
u_t^y(c,h,n)&=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta_t\log n,
\label{eq:new_uy}\\
u^o(c^o,h^o,e)&=\log c^o+\gamma\log h^o+\omega_B\log e.
\label{eq:new_uo}
\end{align}
All preference weights are positive. Households discount future utility by
$\beta>0$. Fertility $n>0$ is continuous; the superscript $o$ labels old-age
quantities. The estate $e$ consists of financial assets and the proceeds from
selling housing at death. It enters the parent's utility but does not determine
an entering household's wealth $b_i$.

A household draws an ownership preference $\xi_i$ before making its choices.
The draw is independent of its endowments and logistic with location $\bar\xi$
and scale $\sigma_\xi>0$.

\textbf{Housing and financial markets.} The housing stock is $\bar H>0$.
Rental units satisfy $h\le h_R^{\max}$ and owner units satisfy
$h\le h_O^{\max}$, where $0<h_R^{\max}<h_O^{\max}$. Goods and bonds trade with
the rest of the world. The bond price is $q\in(0,1)$ and its gross return is
$R_f=1/q$. House prices are $P_t>0$. Owners pay property tax at rate $\tau^p$.
Competition among rental intermediaries gives:
\begin{equation}
u_t\equiv qr_t=(1+q\tau^p)P_t-qP_{t+1}>0.
\label{eq:new_usercost}
\end{equation}
Here $r_t$ is rent paid at the end of the period, and $u_t$ is its value at
the beginning. Property-tax revenue is rebated equally to young and old
households; $T_t$ denotes the rebate at the beginning of the period.

\textbf{Home finance.} Current income and liquid wealth are available at
purchase. A mortgage finances at most the share $\phi_t\in(0,1)$ of the house
price. Principal and accumulated interest are repaid on entering old age.
Unsecured borrowing is unavailable. Old owners can buy or sell housing within
the owner size limit, using their income and wealth without new borrowing.
Thus an old owner's choice is not bounded by the size of its previous home.

\textbf{Demography.} A child produces $\nu>0$ young households next period,
after survival and household formation. Cohort masses satisfy:
\begin{equation}
Y_{t+1}=\nu\bar n_tY_t,\qquad O_{t+1}=Y_t,
\label{eq:new_demography}
\end{equation}
where $\bar n_t$ is average fertility among the young. Population here counts
adult households.

\section{Household choices}

\textbf{Young renters.} Financial wealth on entering old age is $a'$.
Conditional on renting, the young household solves:
\begin{align}
W_t^R(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^R(a';i),\nonumber\\
&c+qa'+u_th=y_i^y+b_i+T_t,\qquad a'\ge0,\quad h\le h_R^{\max}.
\label{eq:new_young_renter}
\end{align}
The renter pays for current housing services and saves for old age.

\textbf{Young owners.} For an owner, $a'$ is financial wealth net of mortgage
repayment. Its problem is:
\begin{align}
W_t^O(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^O(a',h;i),\nonumber\\
&c+qa'+(1+q\tau^p)P_th=y_i^y+b_i+T_t,\nonumber\\
&qa'+\phi_tP_th\ge0,\qquad h\le h_O^{\max}.
\label{eq:new_young_owner}
\end{align}
To see the mortgage directly, let $d$ be principal borrowed at purchase and
$k$ the amount invested in bonds. Then $0\le d\le\phi_tP_th$, $k\ge0$,
$a'=(k-d)/q$, and the budget is
$c+k+(1+q\tau^p)P_th=y_i^y+b_i+T_t+d$.
Eliminating $k,d$ gives \eqref{eq:new_young_owner}.

\textbf{Old households.} An old owner has net financial wealth $a$ and title
to $H$ units of housing. Let $a^e\ge0$ be financial saving during old age.
The estate is:
\begin{equation}
e=\begin{cases}
q^{-1}a^e,&\text{renter},\\
q^{-1}a^e+P_{t+1}h^o,&\text{owner}.
\end{cases}
\label{eq:new_estate}
\end{equation}
The retained house is sold at death, at the end of old age. The two old-age
problems are:
\begin{align}
V_t^R(a;i)=\max_{c^o,h^o,e}\;&u^o(c^o,h^o,e),\nonumber\\
&c^o+qe+u_th^o=a+y_i^o+T_t,\quad h^o\le h_R^{\max},
\label{eq:new_old_renter}\\
V_t^O(a,H;i)=\max_{c^o,h^o,e}\;&u^o(c^o,h^o,e),\nonumber\\
&c^o+qe+u_th^o=a+P_tH+y_i^o+T_t,\nonumber\\
&h^o\le h_O^{\max},\qquad e\ge P_{t+1}h^o.
\label{eq:new_old_owner}
\end{align}
The owner's estate restriction is equivalent to $a^e\ge0$. Its cash budget
before substituting for the estate is
$c^o+a^e+(1+q\tau^p)P_th^o=a+P_tH+y_i^o+T_t$.
Income, liquid wealth, and sale proceeds therefore enter the same budget.

\textbf{Tenure.} The household owns if $W_t^O(i)+\xi_i\ge W_t^R(i)$.
Its ownership probability is:
\begin{equation}
\pi_t^O(i)=
\frac{\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}
{\exp\{W_t^R(i)/\sigma_\xi\}+\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}.
\label{eq:new_tenure}
\end{equation}
The taste draw affects tenure but not choices conditional on tenure.

\Needspace{10\baselineskip}
\section{Equilibrium and the welfare comparison}

\begin{definition}
An equilibrium consists of household choices, tenure probabilities, prices,
rebates, and cohort masses satisfying the household problems, rental pricing,
and demographic equations above. Housing and the property-tax budget clear:
\begin{equation}
Y_t\bar h_t^y+O_t\bar h_t^o=\bar H,
\qquad (Y_t+O_t)T_t=q\tau^pP_t\bar H.
\label{eq:clearing}
\end{equation}
Here $\bar h_t^y$ and $\bar h_t^o$ are average housing occupied by young and old households.
The old distribution is generated by the preceding cohort's choices.
\end{definition}
At a positive stationary equilibrium, $Y=O=N$, $\bar n=1/\nu$, and
$N=\bar H/(\bar h^y+\bar h^o)$.

For the stationary results assume $0\le\tau^p<2$. At stationarity, let
$w=y^y+b+T$ be current cash and $v=y^o+T$ be old income
including the rebate. Write:
\begin{equation}
\begin{gathered}
p=(1-q+q\tau^p)P,\qquad L=(1-\phi+q\tau^p)P,\\
z=a'+Ph+v,\qquad K=1+\gamma+\omega_B,\qquad
D=1+\alpha+\vartheta+\beta K.
\end{gathered}
\label{eq:reduction_objects}
\end{equation}
Here $p$ is the cost of housing services, $L$ is the cash required per unit
of owner housing, and $z$ is the owner's resources on entering old age.
The young owner's two financial restrictions become:
\begin{equation}
c+ph+qz=w+qv,\qquad c+Lh\le w.
\label{eq:reduced}
\end{equation}
The second inequality is the original mortgage limit. If old housing and the
estate restriction are slack, $V(z)=K\log z+C$, with
$c^o=z/K$, $h^o=\gamma z/(Kp)$, and $e=\omega_Bz/(Kq)$.
```

## E. Existing sufficient stationary existence certificate: source excerpt

This excerpt is a proof to independently audit. Its reference to the original numbered Proposition does not adopt that older proposition. No claim of general uniqueness or transition stability accompanies it.

```latex
The following sufficient conditions produce the equilibrium and constrained
owner group in Proposition~\ref{prop:direct}. They are conservative analytical
bounds, with no equilibrium price or multiplier as an input. Assume bounded
endowments, $w_{0i}=y_i^y+b_i\ge\underline w>0$, $v_{0i}=y_i^o\ge0$, and
$0\le\tau^p<2$. Define:
\begin{equation}
\begin{gathered}
d_p=1-q+q\tau^p,\quad d_L=1-\phi+q\tau^p,\quad
d_{\max}=\max\{d_p,d_L\},\\
\bar M_0=\int(w_0+qv_0)\,\mathrm dF,\quad
\bar T=\frac{\tau^p\bar M_0}{(1-q)(2-\tau^p)},\quad
\bar M=\bar M_0+(1+q)\bar T.
\end{gathered}
\end{equation}
Require enough potential fertility when housing is inexpensive:
\begin{equation}
h_R^{\max}>\kappa/\nu,\qquad
\vartheta\nu>\chi D/\underline w+
\frac{\alpha\kappa}{h_R^{\max}-\kappa/\nu}.
\label{eq:existence_conditions}
\end{equation}
Then a positive stationary equilibrium exists, and every such equilibrium has
$0\le T\le\bar T$ and $P_-<P<P_+$, where:
\begin{equation}
P_-\equiv\frac{\alpha\underline w}{D d_{\max}h_R^{\max}},
\qquad P_+\equiv\frac{\nu\bar M}{d_p\kappa}.
\label{eq:price_bounds}
\end{equation}

To verify existence, combine the original budgets in either tenure:
\begin{equation}
x+ps+(\chi+p\kappa)n+qc^2+qp h^2+q^2e=w+qv.
\label{eq:lifetime}
\end{equation}
Scaling the first-order conditions gives
$\lambda(w+qv)+\mu w\le D$, while $1/x=\lambda+\mu$; hence
$x\ge w/D\ge\underline w/D$. Below $P_-$, the housing condition
$\alpha x/s\le d_{\max}P$ when uncapped, or the binding cap itself,
implies $h\ge h_R^{\max}$. Equation \eqref{eq:new_fertility} and
\eqref{eq:existence_conditions} then give fertility above replacement for
all conditional choices. Above $P_+$, \eqref{eq:lifetime} gives mean
fertility below replacement for every $T\le\bar T$.

The rebate map, writing $\bar h=\bar h^y+\bar h^o$, satisfies:
\begin{equation}
\mathcal T(P,T)=q\tau^pP\bar h/2
\le\frac{\tau^p}{2d_p}\{\bar M_0+(1+q)T\}\le\bar T.
\end{equation}
The final inequality uses
$2d_p-\tau^p(1+q)=(1-q)(2-\tau^p)>0$.
Conditional allocations and tenure probabilities are continuous. On the
compact rectangle of prices and rebates, move the price in the direction of
$\nu\bar n-1$, clipping it to $[P_-,P_+]$, and update the rebate by
$\mathcal T$. A Brouwer fixed point has an interior price by the strict
boundary signs, and therefore replacement fertility. Setting
$N=\bar H/\bar h$ completes the stationary equilibrium. This argument
establishes existence, not uniqueness or transition stability.
```
````

### File: docs/model/simplified_olg_general_preferences.md
```md
# General preferences: housing allocation and fertility

Working memo, September 8, 2026. This separate branch starts from unspecified
preferences over **gross** goods, housing, and fertility. It changes neither
the household model nor the dated planner agreed in
`simplified_olg_utilitarian_work.md`. Linear child needs enter only
the final special case.

## 1. Unspecified utility and the maintained benchmark

Let young utility be \(U^y(c,h,n)\). Consumption \(c\) and housing \(h\) are
the gross resources in the existing budgets. Specify its domain as a primitive;
do not initially impose \(c>\chi n\) or \(h>\kappa n\). Assume twice continuous
differentiability, joint concavity, and \(U^y_c,U^y_h>0\). Old utility
\(U^o(c_o,h_o,e)\) is increasing and jointly concave. All comparisons retain
the specified cardinal utility scales.

**Do not assume \(U^y_n>0\) everywhere at fixed \(c,h\).** Without a separate
fertility resource cost or an upper domain bound, that assumption rules out
a finite interior fertility optimum. Net marginal utility of children may
become negative because of costs or crowding represented inside \(U^y\).
Concavity alone does not guarantee existence on an unbounded fertility domain.
The derivative results below concern an existing regular interior optimum;
active domain restrictions would add their own multipliers.

Each current age group has mass \(N\) and the same stationary probability law
\(Q\) of paired types and retained tenures. The physical cap is
\(H_i\in\{H_R,H_O\}\), with \(H_R<H_O\). First fix individual fertility and
tenure. The planner chooses **all current consumption and housing**, fixes
future real continuation opportunities and old net estates, and relaxes
individual financing restrictions while honoring obligations. It maximizes
\[
N\int[U^y(c_i^y,h_i^y,n_i)+U^o(c_i^o,h_i^o,e_i)]\,dQ
\]
subject to
\[
N\int(c_i^y+c_i^o)dQ=C,\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,\qquad h_i^a\le H_i.
\]
Here \(C,\bar H\) are the reference resource totals. Write
\(H_Y=N\int h_i^y\,dQ\), \(H_O=N\int h_i^o\,dQ\); superscripts \(eq,*\)
denote the competitive reference and planner allocation.
Continuation and fixed-tenure taste terms are constant. The owner's financial
estate floor is relaxed; the promised net estate is unchanged.

With resource multipliers \(\lambda_C,\lambda_H\) and planner cap multipliers
\(\eta_i^{a*}\ge0\), interior-domain optimality requires
\[
U^y_c=U^o_c=\lambda_C,\qquad
U^y_h=\lambda_H+\eta_i^{y*},\qquad
U^o_h=\lambda_H+\eta_i^{o*}.
\]
Together with feasibility and complementary slackness these characterize an
optimum under concavity, if one exists. Strict concavity ensures uniqueness.
Without separability, consumption changes housing marginal utilities.

## 2. The competitive wedge is fully general

At a positive stationary equilibrium, let \(q=1/R_f\),
\(p=(1-q+q\tau^p)P\), \(L_R=p\), and
\(L_O=(1-\phi+q\tau^p)P\). Here \(p\) is housing's service cost and \(L_d\)
its coefficient in the current financing constraint. Current cash is
\(w_i=y_i^y+b_i+T\), old income including rebate is \(v_i^o=y_i^o+T\), and
\(z_i\) is total old-age resources after repayment of the young mortgage.
The reduced budgets remain
\[
c+ph+qz=w_i+qv_i^o,\qquad c+L_dh\le w_i.
\]
The young objective is \(U^y(c,h,n)+\beta V_d(z)\). The maintained continuation
value has **no direct dependence on \(n\)**; \(\beta\) discounts the parent's
own old age. For budget and financing multipliers \(\Lambda_i,\mu_i\), and
competitive young cap multiplier \(\eta_i^y\),
\[
q\Lambda_i=\beta m_i,\qquad
U^y_c=\Lambda_i+\mu_i,\qquad
U^y_h=p\Lambda_i+L_d\mu_i+\eta_i^y,
\quad m_i=V_d'(z_i)>0.
\]
For the old budget \(c_o+ph_o+qe=z\), put \(\rho_i\ge0\) on the owner's
floor \(e-Ph_o\ge0\), and set \(\rho_i=0\) for renters. Let \(\eta_i^o\)
be the competitive old housing-cap multiplier. The envelope and old
first-order conditions give
\[
U^o_c=m_i,\qquad U^o_e=qm_i-\rho_i,\qquad
U^o_h=pm_i+P\rho_i+\eta_i^o.
\]
Stationarity matches the young household's future old allocation with its
current old counterpart. Consequently
\[
\boxed{U^y_h-U^o_h
=\left(\frac{\beta}{q}-1\right)pm_i
+L_d\mu_i+\eta_i^y-P\rho_i-\eta_i^o.}
\]
This uses neither logarithms nor linear child costs. The old estate floor
raises old direct housing marginal utility. No restriction on \(\beta R_f\)
has been imposed. A bound used to sign this identity is an additional
sufficient restriction, not a consequence of logarithmic utility.

A positive gap permits a small improving transfer to an **uncapped** young
recipient, holding consumption fixed. It does not establish the full
optimum's aggregate direction. Young finance need not dominate old finance
or the age-weight term.

Cardinal comparison matters independently of curvature: replacing \(U^o\)
by \(A U^o\) and \(\beta\) by \(\beta/A\) leaves all competitive choices
unchanged but multiplies the planner's weight on old utility by \(A>0\).
Thus concavity alone cannot establish a universal age direction. This is a
family of different social comparisons, not a normalization of one fixed
criterion.

## 3. Fertility and a useful intermediate class

Given a gross bundle and fixed continuation opportunities, an interior
fertility choice satisfies \(U^y_n=0\). If \(U^y_{nn}<0\), then
\[
dn=-\frac{U^y_{nc}\,dc+U^y_{nh}\,dh}{U^y_{nn}}.
\]
Positive consumption and housing responses therefore require positive
cross-partials with fertility. Joint concavity does not sign them.

For a counterexample with no linear resource offsets, take
\[
U^y(c,h,n)=\sqrt c+\sqrt{h+n}+A\sqrt n-kn,\qquad A,k>0.
\]
This is strictly concave and increasing in \(c,h\). At each positive bundle
there is a unique interior fertility optimum: \(U_n\) decreases from
\(+\infty\) to \(-k\). Yet
\(U_{nh}=-[4(h+n)^{3/2}]^{-1}<0\), so more housing lowers fertility.

An intermediate specification separates the two interactions:
\[
\boxed{U^y(c,h,n)=a(c,n)+b(h,n)+v(n).}
\]
Assume \(a_c,b_h,v'>0\); \(a_n,b_n\) may be negative due to goods costs and
crowding. Require joint concavity and an existing interior optimum. Define
\(D=a_{nn}+b_{nn}+v''<0\). Its fertility condition and response are
\[
a_n+b_n+v'=0,\qquad
dn=-\frac{a_{cn}\,dc+b_{hn}\,dh}{D}.
\]
Thus \(a_{cn}>0\) and \(b_{hn}>0\) suffice for both resource effects to be
positive, without a shift specification. They say that additional resources
raise the marginal utility of children. If consumption falls while housing
rises, this numerator gives the local tradeoff.

At fixed fertility this class also separates the planner's consumption and
housing problems. A useful **additional age-comparison restriction** is
\[
U^o=f_o(c_o)+b(h_o,0)+B(e),\qquad b_{hh}<0,\quad b_{hn}>0.
\]
The old housing component is cardinally matched to the childless young
component. Suppose the cross-partial restriction holds between \(0\) and
each \(n_i>0\), on a common feasible housing domain, and positive solutions
exist. Then \(b_h(h,n_i)>b_h(h,0)\). A common planner housing multiplier
therefore yields
\[
h_i^{y*}\ge h_i^{o*},
\]
strictly wherever the old counterpart is uncapped. If
\(\bar H<2N\int H_i\,dQ\), not all old households can be capped, so
\(H_Y^*>\bar H/2\). Hence \(H_Y^{eq}\le H_O^{eq}\) suffices for a strict
aggregate young housing gain. This is a theorem under explicit preference
and reference-allocation restrictions, not a generic implication of finance.
Neither this result nor \(b_{hn}>0\) proves that every young household gains.

## 4. Joint fertility choice by the dated planner

If the planner also chooses \(n_i\), the gross-variable conditions are simply
\[
U^y_c=\lambda_C,\qquad
U^y_h=\lambda_H+\eta_i^{y*},\qquad U^y_n=0
\]
at an interior choice. There is no additional fertility resource term in
the maintained gross budgets. The planner values children through currently
living parents; it attaches no new welfare weight to future people.
Conditional fertility responses do not establish greater average fertility
at this joint optimum. Changing births also changes future entrant masses,
so this is not a completed dynamic allocation holding every future
population fixed.

Separate explicit child-resource constraints would constitute another model
architecture. They require their own definitions and are not adopted here.
The same child goods or space already included in gross \(c,h\) must not be
charged again as an additional aggregate resource requirement.

## 5. Restricted mapping: the existing shifted specification

Only now impose the existing form
\[
a(c,n)=\log(c-\chi n),\qquad
b(h,n)=\alpha\log(h-\kappa n),\qquad
v(n)=\vartheta\log n.
\]
Its domain \(c>\chi n,\ h>\kappa n,\ n>0\) belongs to **this special case**.
Old utility is \(\log c_o+\gamma\log h_o+\omega_B\log e\).
Here \(a_{cn}=\chi/(c-\chi n)^2>0\) and
\(b_{hn}=\alpha\kappa/(h-\kappa n)^2>0\), so the intermediate fertility result
reproduces the existing positive conditional responses. The gross condition
\(U^y_n=0\) becomes exactly
\[
\frac{\vartheta}{n}
=\frac{\chi}{c-\chi n}+\frac{\alpha\kappa}{h-\kappa n}.
\]

A still-restricted nonlogarithmic extension replaces the two logarithms by
\(f(c-\chi n)\) and \(\alpha g(h-\kappa n)\), with increasing, strictly
concave functions. With additively separable old housing utility
\(\gamma g(h_o)\), fixed-fertility planner housing is
\[
h_i^{y*}=\min\{H_i,\kappa n_i+r_y\},\quad
h_i^{o*}=\min\{H_i,r_o\},\quad
r_y=(g')^{-1}(\lambda_H/\alpha),\quad
r_o=(g')^{-1}(\lambda_H/\gamma).
\]
Assume positive solutions and the appropriate derivative range, extending
the inverse to infinity when desired housing exceeds every finite cap.
Let \(\bar n=\int n_i\,dQ\). For \(\alpha=\gamma\) and slack planner caps,
\[
r_y=r_o=\frac{\bar H/N-\kappa\bar n}{2},\qquad
H_Y^*=\frac{\bar H+\kappa N\bar n}{2}.
\]
This equalizes **adult space**, a concept specific to the shifted class.
Young household \(i\) gains housing exactly when
\(h_i^{y,eq}-\kappa n_i<r_y\). Aggregate young housing increases exactly when
old mean reference housing exceeds young mean reference adult space.
For \(\alpha\ge\gamma\), the paired planner housing ordering survives binding
caps. These formulas are not conclusions for unspecified \(U^y(c,h,n)\).

The fully general competitive age ordering, equilibrium existence, aggregate
fertility comparison, and funded policy transition remain unproved.
```

### File: output/model/simplified_olg_amendments/oracle_housing_fertility_response.md
```md
# Pro response — housing restrictions and fertility through parents

Source: https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354

Third response, following attachment Pasted text(20260908-235319).txt. Completed in 22 minutes 19 seconds. Captured September 8 Eastern (September 9 UTC), from visible text and mathematical source labels. Display delimiters reconstructed for readability. This is an external review, not an adopted revision to the model, note or slides.

# Housing allocation and fertility valued through parents

**The half-stock lemma is correct. A savings-based restriction can generate the required market housing ordering without the earlier lower-price certificate. Fertility, however, requires a separate argument: neither more housing nor higher parental welfare alone establishes its direction.**

Throughout, retain the specified one-date benchmark: existing adults, retained tenure, total current goods and housing, existing obligations, incumbent continuation opportunities, and old estate payments. The planner relaxes both young financing and old nonnegative financial saving; the competitive economy does not. Pasted text

Write \(H_d=h_d^{\max}\), \(d\in\{R,O\}\), and normalize resources by the reference cohort mass:

\[
\bar c=\frac{C^{eq}}N,\qquad
\bar h=\frac{\bar H}N.
\]
These are totals across both ages **per cohort mass**, not per living household. Reference young quantities are

\[
x_i=c_i^{y,eq}-\chi n_i,\qquad
s_i=h_i^{y,eq}-\kappa n_i,\qquad
\bar n=\frac1\nu.
\]
Superscripts \(F,S,J\) below denote the fixed-fertility planner, subsequent conditional fertility choice, and joint planner.

## 1. Housing: a theorem allowing binding physical caps

Define

\[
\begin{gathered}
E=1+\alpha+\vartheta,\qquad K=1+\gamma+\omega_B,\\
d_p=1-q+q\tau^p,\qquad a_p=1+q\tau^p,\\
\Gamma_R=\gamma,\qquad
\Gamma_O=\min\left\{\gamma,\frac{(\gamma+\omega_B)d_p}{a_p}\right\},
\qquad \Gamma=\Gamma_O.
\end{gathered}
\]
A transparent sufficient preference restriction is

\[
\boxed{\alpha\ge\gamma,\qquad
       \frac{\beta\Gamma}{q}\ge\alpha+\vartheta.}
\tag{H}
\]

### Housing proposition

Take any positive stationary competitive equilibrium satisfying (H) and

\[
\bar H<2N\int H_{d_i}\,dQ.
\]
**Individual comparison.** Every young household in

\[
\boxed{\mathcal I=\{i:h_i^{y,eq}<H_{d_i},\ s_i\le\bar s\}}
\tag{1}
\]
receives strictly more housing at the full fixed-fertility optimum. A positive-mass conclusion requires \(Q(\mathcal I)>0\); an explicit primitive recipient group is supplied below.

**Aggregate comparison.** Equilibrium generates

\[
h_i^{o,eq}\ge h_i^{y,eq}
\quad\text{for every paired type,}
\]
and the full optimum satisfies

\[
\boxed{H_y^F>\frac{\bar H}{2}\ge H_y^{eq}.}
\tag{2}
\]
No income ratio, slack old estate restriction, or absence of binding housing caps is assumed in this proposition. It is conditional on a positive equilibrium, not a general existence theorem.

### Proof

The full fixed-fertility solution is

\[
\begin{aligned}
c_i^{y,F}&=\chi n_i+x^F,&c_i^{o,F}&=x^F,
&x^F&=\frac{\bar x+\overline{c^{o,eq}}}{2},\\
h_i^{y,F}&=\min\{H_{d_i},\kappa n_i+\alpha/\lambda_F\},&
h_i^{o,F}&=\min\{H_{d_i},\gamma/\lambda_F\},
\end{aligned}
\tag{3}
\]
where \(\lambda_F>0\) clears housing.

To establish the market ordering, let

\[
p=d_pP,\qquad L_R=p,\qquad L_O=(1-\phi+q\tau^p)P.
\]
Let \(\Lambda_i,\mu_i,\eta_i\) be the lifetime-budget, young-finance, and young-cap multipliers. Household optimality gives

\[
\begin{gathered}
\Lambda_i=\frac{\beta}{q c_i^{o,eq}},\qquad
\frac1{x_i}=\Lambda_i+\mu_i,\\
s_i=\frac{\alpha}{p\Lambda_i+L_{d_i}\mu_i+\eta_i},\\
n_i=\frac{\vartheta}
{\chi(\Lambda_i+\mu_i)+\kappa(p\Lambda_i+L_{d_i}\mu_i+\eta_i)}.
\end{gathered}
\]
Consequently,

\[
h_i^{y,eq}
\le
\frac{\displaystyle\alpha+
          \frac{\vartheta p\kappa}{\chi+p\kappa}}
     {p\Lambda_i}
<
\frac{\alpha+\vartheta}{p\Lambda_i}.
\tag{4}
\]
If old housing is uncapped, the two estate regimes give

\[
h_i^{o,eq}
=\frac{\Gamma_{d_i}c_i^{o,eq}}p
=\frac{\beta\Gamma_{d_i}}{q p\Lambda_i}.
\]
Condition (H) therefore implies the ordering. If old housing is capped, it equals the young household’s retained cap and the ordering follows immediately.

Inside the planner allocation, \(\alpha\ge\gamma\) and \(n_i>0\) imply \(h_i^{y,F}\ge h_i^{o,F}\), strictly wherever old housing is uncapped. If every old household were capped, every young household would also be capped, contradicting the stock restriction. This proves (2).

Finally, put

\[
S=\bar h-\kappa\bar n=\bar s+\overline{h^{o,eq}}.
\]
Market ordering gives \(S\ge2\bar s+\kappa\bar n\), while (3) gives

\[
S\le\frac{\alpha+\gamma}{\lambda_F}.
\]
Hence \(\alpha/\lambda_F>\bar s\). Both arguments of the minimum in (3) exceed reference housing for every \(i\in\mathcal I\). ∎

### Interpretation and sharper restrictions

The old-income and home-sale proceeds are already included in old resources. Nothing in this proof prevents old owners from selling or downsizing. The smaller coefficient \(\Gamma_O\) captures their inability to borrow against the estate **after** selling and choosing their new home. Pasted text

The restriction \(\alpha\ge\gamma\) compares adult-space weights relative to the common consumption-log coefficient. Children weaken the restrictions needed for particular conclusions:

\[
\frac{\alpha}{s_i}>\frac{\gamma}{h_i^{o,eq}}
\quad\Longleftrightarrow\quad
\frac{\alpha}{\gamma}>
\frac{h_i^{y,eq}-\kappa n_i}{h_i^{o,eq}}.
\]
Thus, with \(h_i^{o,eq}\ge h_i^{y,eq}\), the sufficient threshold is
\(\alpha/\gamma>1-\kappa n_i/h_i^{y,eq}\), not necessarily one.

Separately, when the fixed-fertility planner is uncapped,

\[
H_y^F>\bar H/2
\quad\Longleftrightarrow\quad
\boxed{\frac{\alpha}{\gamma}>
1-\frac{2\kappa\bar n}{\bar h}.}
\]
With binding caps, the exact minimum formulas must be used; a sufficient typewise comparison is \(\alpha+\lambda_F\kappa n_i\ge\gamma\), with strictness on a positive-mass uncapped-old set. These weakenings do not automatically preserve every individual conclusion above.

Condition (H) can be demanding. When the estate restriction binds and \(d_p\) is small, \(\Gamma_O\) can be much smaller than \(\gamma\). The required two-age patience can then be high. No calibration compatibility follows without a mapping from the quantitative model.

Saving and old income can substitute for one another. Define

\[
A(p)=\alpha+\frac{\vartheta p\kappa}{\chi+p\kappa}.
\]
For renters, and for owners in the explicit mortgage subcase \(\phi=q\), the sharper condition

\[
\boxed{
\Gamma_d\max\left\{\frac{\beta}{q},
                  \frac{E v}{K w}\right\}\ge A(p)
}
\tag{5}
\]
suffices for \(h^{o,eq}\ge h^{y,eq}\), including caps. With both caps slack, it is exact. The patience part follows from (4). For the income part, current expenditure \(B=c+ph\le w\) implies

\[
h^y\le\min\{H_d,A(p)w/(Ep)\},\qquad z\ge v,
\]
which yields the result. This is considerably more informative than bounding housing by the entire cash endowment.

### An explicit primitive recipient group

For a solved illustration, take **\(\phi=q,\ \tau^p=0\)** and require the following computed market demands to lie below their finite tenure caps. Define

\[
\begin{gathered}
w_{0i}=y_i^y+b_i,\qquad v_{0i}=y_i^o,\qquad D=E+\beta K,\\
x_i=\min\left\{\frac{w_{0i}}E,\frac{w_{0i}+qv_{0i}}D\right\},\\
c_i^{o,eq}=\frac{w_{0i}+qv_{0i}-Ex_i}{qK},\qquad
p=\frac{\nu\vartheta\bar x-\chi}{\kappa}>0,\\
n_i=\frac{\vartheta x_i}{\chi+p\kappa},\qquad
h_i^{y,eq}=\frac{A(p)x_i}{p},\qquad
h_i^{o,eq}=\frac{\Gamma_{d_i}c_i^{o,eq}}p.
\end{gathered}
\tag{6}
\]
These equations solve the stationary subcase, with \(N\) determined by housing clearing.

Conditional young choices coincide across tenures. The old value difference is constant across endowments, so the original logistic taste gives a constant owner share \(\pi\in(0,1)\). Explicitly,

\[
\pi=\operatorname{logit}^{-1}
\left(\frac{\bar\xi+\beta\Delta}{\sigma_\xi}\right),
\]
where \(\Delta=0\) on the unrestricted-estate branch; otherwise, writing \(J=\gamma+\omega_B\),

\[
\Delta=J\log J+\gamma\log\frac{1-q}{\gamma}
                    +\omega_B\log\frac q{\omega_B}.
\]
Let \(\bar\Gamma=(1-\pi)\gamma+\pi\Gamma_O\). The **aggregate**, rather than all-type, market-ordering restriction is exactly

\[
\boxed{\bar\Gamma\,\overline{c^{o,eq}}\ge A(p)\bar x.}
\tag{7}
\]
Under (7) and \(\alpha\ge\gamma\), every \(x_i\le\bar x\) gains both consumption and housing in the fixed-fertility planner, even if planner caps bind. A primitive positive-mass group with **strictly binding young finance** is

\[
\boxed{
F\left\{w_{0i}<E\bar x,\quad
qv_{0i}>\frac{\beta K}{E}w_{0i}\right\}>0.
}
\tag{8}
\]
Thus the relevant group is cash-poor relative to the cross-section but sufficiently future-resource-rich. Heterogeneity and both tenures remain.

These are redistributive utilitarian results. Under (H), patience itself contributes to the age gap. Only in the diagnostic case \(\beta=q,\Gamma_d=\gamma\), with slack physical caps, does the paired marginal housing gap reduce to \(\mu_iL_d\).

Finally, your correction to the previous certificate is right: it implied

\[
\beta<\frac{\nu\vartheta\bar w_0/\chi-E}{K}.
\]
Calling that certificate “beta-unrestricted” was incorrect. The new general ordering proof does not use it; the solved subcase has its own explicit \(p>0\) feasibility requirement. Pasted text

## 2. Fertility chosen privately within the assigned bundle

Let \(\mathfrak n(c,h)\) denote the parent’s optimal fertility when the current bundle and continuation opportunities are fixed. It is the unique zero, on \(0<n<\min\{c/\chi,h/\kappa\}\), of

\[
f(n;c,h)=\frac{\vartheta}{n}
-\frac{\chi}{c-\chi n}
-\frac{\alpha\kappa}{h-\kappa n}.
\]
This is exactly the maintained parental fertility condition, not a valuation of unborn utility. Pasted text

Since \(f_n<0\), differentiation confirms

\[
\boxed{
dn=
\frac{(\chi/x^2)\,dc+(\alpha\kappa/s^2)\,dh}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
}
\tag{9}
\]
For a finite change accommodating baseline fertility \(n_0\),

\[
\boxed{
n_1\ge n_0
\Longleftrightarrow
\frac{\chi}{c_1-\chi n_0}
+\frac{\alpha\kappa}{h_1-\kappa n_0}
\le\frac{\vartheta}{n_0}.
}
\tag{10}
\]
If the new bundle cannot accommodate \(n_0\), fertility falls.

In particular, let housing rise by \(\Delta h>0\), while consumption falls by \(\delta\ge0\). Evaluating \(x,s\) at the original bundle,

\[
\boxed{
n_1\ge n_0
\Longleftrightarrow
\delta\le
\frac{\alpha\kappa x^2\Delta h}
{\chi s(s+\Delta h)+\alpha\kappa x\Delta h}.
}
\tag{11}
\]
Strict inequality gives strictly higher fertility.

Locally, the allowable consumption loss per unit of housing is

\[
-\frac{dc}{dh}<\frac{\alpha\kappa x^2}{\chi s^2}.
\]
For the illustrative exchange \(dc=-p\,dh\), at an uncapped, financially unrestricted reference bundle \(s=\alpha x/p\), this becomes

\[
\boxed{p\kappa>\alpha\chi.}
\]
The space component of child costs must be sufficiently important. This is a bundle-composition test, not an assertion that a particular tax or credit policy implements the exchange.

### Applying the test to the actual fixed-fertility optimum

Put \(s_i^F=h_i^{y,F}-\kappa n_i\). Then

\[
\boxed{
n_i^S>n_i
\Longleftrightarrow
\frac{\chi}{x^F}+\frac{\alpha\kappa}{s_i^F}
<\frac{\vartheta}{n_i}.
}
\tag{12}
\]
The low-resource recipients identified in (6)–(8) gain both components and therefore have strictly higher conditional fertility. Other young households may face a consumption–space tradeoff.

Crucially, \(x^F,s_i^F\) evaluate the new bundle **at original fertility**. After fertility responds, adult goods and space are

\[
x_i^S=x^F-\chi(n_i^S-n_i),\qquad
s_i^S=s_i^F-\kappa(n_i^S-n_i).
\]
An aggregate result needs to include losers. One cap-valid sufficient test is

\[
\boxed{
\int\mathfrak n(x^F,s_i^F)\,dQ
>\frac{1+\alpha}{E}\bar n
\quad\Longrightarrow\quad
\bar n^S>\bar n.
}
\tag{13}
\]
To prove it, \(\mathfrak n\) is homogeneous and concave: its upper contour sets are convex by (10), and homogeneity turns this into concavity. Hence it is superadditive, and

\[
n_i^S
=\mathfrak n(x^F+\chi n_i,s_i^F+\kappa n_i)
\ge\mathfrak n(x^F,s_i^F)+\frac{\vartheta}{E}n_i.
\]
Integrating proves (13). This is a sufficient allocation test, not a primitive equilibrium restriction.

No market-policy sign has been claimed. The direct assignments are financed by the maintained balanced transfers

\[
t_k=\Delta c_k+u_t\Delta h_k,
\]
with the previously specified financial adjustments preserving incumbent future resources and estates. A cash grant or credit reform instead requires household reoptimization, fiscal funding, price clearing, and—when tenure is free—the ownership-selection contribution to average fertility.

## 3. The same planner chooses fertility, valuing parents only

Use \(x_i=c_i^y-\chi n_i\) and \(s_i=h_i^y-\kappa n_i\). The normalized problem is

\[
\max\int\left[
\log x_i+\alpha\log s_i+\vartheta\log n_i
+\log c_i^o+\gamma\log h_i^o
\right]dQ
\]
subject to

\[
\begin{aligned}
\int(x_i+\chi n_i+c_i^o)dQ&=\bar c,\\
\int(s_i+\kappa n_i+h_i^o)dQ&=\bar h,\\
s_i+\kappa n_i&\le H_{d_i},\qquad h_i^o\le H_{d_i},
\end{aligned}
\tag{14}
\]
and positive log arguments. This merely rewrites total \(c,h\); it does not add child costs a second time. No replacement-fertility constraint is imposed. Pasted text

The objective is strictly concave in these variables and the constraints are affine. Averaging within retained tenure preserves feasibility and raises utility, reducing existence and characterization to a finite-dimensional problem. The positive reference supplies feasibility.

Let \(\lambda_C>0,\lambda_H>0\) be the goods and housing resource multipliers. Let

\[
r_d=\lambda_H+\eta_d,\qquad \eta_d\ge0
\]
include the young tenure-cap multiplier. The unique optimum satisfies

\[
\boxed{
\begin{aligned}
x_d^J=c_d^{o,J}&=\frac1{\lambda_C},\\
n_d^J&=\frac{\vartheta}{\chi\lambda_C+\kappa r_d},\\
h_d^{y,J}&=\frac{\alpha}{r_d}+\kappa n_d^J,\\
h_d^{o,J}&=\min\{H_d,\gamma/\lambda_H\}.
\end{aligned}}
\tag{15}
\]
For an uncapped young tenure, \(r_d=\lambda_H\). Otherwise \(r_d>\lambda_H\) is the unique solution of

\[
\frac{\alpha}{r_d}
+\frac{\kappa\vartheta}{\chi\lambda_C+\kappa r_d}=H_d.
\]
Writing \(\pi_d=Q(d)\), the resource equations determining the multipliers are

\[
\frac2{\lambda_C}+\chi\sum_d\pi_dn_d^J=\bar c,\qquad
\sum_d\pi_d(h_d^{y,J}+h_d^{o,J})=\bar h.
\tag{16}
\]
Thus individual fertility rises exactly for \(n_i<n_{d_i}^J\); individual housing rises exactly for \(h_i^{y,eq}<h_{d_i}^{y,J}\). The fixed-fertility recipient set need not remain unchanged.

### Useful cap-aware sufficient inequalities

The optimality conditions imply

\[
\lambda_C\le U_C:=\frac{2+\vartheta}{\bar c},\qquad
r_d\le U_d:=
\max\left\{\frac{\alpha+\gamma+\vartheta}{\bar h},
           \frac{\alpha+\vartheta}{H_d}\right\}.
\]
Indeed, \(\chi\lambda_C n_d^J<\vartheta\); summing the goods and housing optimality identities gives the resource bounds, while a binding young cap implies \(r_d<(\alpha+\vartheta)/H_d\).

Consequently,

\[
\underline n_d=
\frac{\vartheta}{\chi U_C+\kappa U_d},\qquad
\underline h_d=\frac{\alpha}{U_d}+\kappa\underline n_d
\]
are lower bounds on joint-planner fertility and young housing. Therefore

\[
\boxed{
n_i<\underline n_{d_i}\Rightarrow n_i^J>n_i,\qquad
h_i^{y,eq}<\underline h_{d_i}\Rightarrow h_i^{y,J}>h_i^{y,eq},
}
\tag{17}
\]
and, separately,

\[
\boxed{\sum_d\pi_d\underline n_d>\frac1\nu
\Rightarrow\bar n^J>\bar n.}
\tag{18}
\]
These are conservative **resource restrictions**, not primitive equilibrium inequalities. In the solved benchmark (6), the resource totals are explicit functions of primitives.

Importantly, under \(\alpha\ge\gamma\), (15) also gives

\[
H_y^J>\bar H/2.
\]
If young housing is capped it is at least old housing; otherwise \(r_d=\lambda_H\) and positive child space makes it larger. Thus the housing theorem’s aggregate gain survives free fertility **even with caps**, independently of whether average fertility rises.

### A sharper fertility theorem when the joint optimum is uncapped

Then fertility is common across all young households and is the unique solution of

\[
\boxed{
\frac{\vartheta}{n^J}
=\frac{2\chi}{\bar c-\chi n^J}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa n^J}.
}
\tag{19}
\]
The exact aggregate comparison is

\[
\boxed{
n^J>\bar n
\Longleftrightarrow
\frac{\vartheta}{\bar n}>
\frac{2\chi}{\bar c-\chi\bar n}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa\bar n}.
}
\tag{20}
\]
**Under (H), this inequality follows from household optimality**, provided the solution of (19) respects the finite caps. The Euler condition gives \(\overline{c^{o,eq}}>\bar x\), and the housing theorem gives

\[
\overline{h^{o,eq}}\ge\bar s+\kappa\bar n>
\frac{\gamma}{\alpha}\bar s.
\]
Moreover, multiplying the reference fertility condition by \(n_i\), integrating, and applying Cauchy–Schwarz yields

\[
\frac{\vartheta}{\bar n}
\ge\frac{\chi}{\bar x}+\frac{\alpha\kappa}{\bar s}.
\]
The right side of (20) is strictly smaller than this last expression. Hence \(n^J>\bar n\), and every initially below-average-fertility parent gains fertility.

If **both** planner solutions are uncapped,

\[
H_y^J-H_y^F
=\frac{N\gamma\kappa}{\alpha+\gamma}(n^J-\bar n)>0.
\tag{21}
\]
Concavity of \(\mathfrak n\), applied to the sequential bundles, also gives \(\bar n^S<n^J\) in this regime. There is no corresponding unrestricted componentwise ordering between \(J\) and \(S\). Generally, only their welfare ordering is automatic:

\[
\mathcal W^J\ge\mathcal W^S\ge\mathcal W^F.
\]
Free fertility need not rise without the restrictions. In the uncapped, young-finance-slack, unrestricted-estate subcase,
\(c_i^{o,eq}=(\beta/q)x_i\) and \(h_i^{o,eq}=\gamma c_i^{o,eq}/p\). Setting \(r=\beta/q\), the fertility test has the sign of

\[
(r-1)\left[\frac{\chi}{1+r}
+\frac{p\kappa\gamma}{\alpha+\gamma r}\right].
\]
Thus \(r<1\) gives \(n^J<\bar n\), including with heterogeneous incomes. This is an age-weighting effect, not a value assigned to unborn people.

## Bridge to the transition

The conditional tests (9)–(12) apply at each date using the same exogenous \(\vartheta_t\) on both paths; they require no convergence. The stationary housing-pairing shortcut does not automatically apply when age distributions differ. A policy path additionally requires fiscal funding, prices, financial obligations, housing clearing, and tenure selection. From common initial cohorts, date-by-date fertility ordering implies cohort ordering through

\[
\frac{Y_t^{policy}}{Y_t^{baseline}}
=\prod_{s=t_0}^{t-1}
\frac{\bar n_s^{policy}}{\bar n_s^{baseline}},
\qquad O_{t+1}=Y_t.
\]
It does not prove that a policy produces that fertility ordering. Positive stationary endpoints both have \(\bar n=1/\nu\); population-level differences can reflect fertility accumulated along the path. Today’s parent-only exercise neither imposes replacement fertility after reallocation nor fixes the future allocations of an unchanged set of unborn households.
```

### File: output/model/simplified_olg_amendments/oracle_dated_allocation_response.md
```md
# Pro response — one-date allocation, individual and aggregate comparisons

Source: https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354

Completed response: “Worked for 23m 41s”. Captured 2026-09-08 from visible browser text and mathematical source labels. Display delimiters reconstructed for readability. This is an external review, not an accepted change to the model or paper.

## Result

**A full-optimum individual result is possible, but “the young household is borrowing-constrained” is not sufficient.** The additional restriction must compare its current cash with the resources available across the cross-section. Below, that comparison is expressed entirely through parameters and the endowment distribution—not equilibrium housing choices or assumed marginal-utility gaps.

The aggregate conclusion requires a separate condition. In fact, there is an analytical family in which **every young household is financially constrained, yet the full dated planner reduces aggregate young housing**.

## 1. The dated planner and financial feasibility

Write \(H_d=h_d^{\max}\), for tenure \(d\in\{R,O\}\). Let \(Q\) be the common stationary distribution of complete types and retained tenure, and let \(N\) be each age-group’s reference mass. Throughout, individual fertility, tenure, future real allocations, and estate payments remain fixed. This is the one-date, equally weighted remaining-utility comparison specified in the follow-up. Pasted text

Let \(C^{\mathrm{eq}}\) be total reference current consumption expenditure. The planner solves

\[
\max N\int\left[
\log(c_i^y-\chi n_i)+\alpha\log(h_i^y-\kappa n_i)
+\log c_i^o+\gamma\log h_i^o
\right]dQ
\]
subject to

\[
N\int(c_i^y+c_i^o)dQ=C^{\mathrm{eq}},
\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,
\]
the positive log arguments, and

\[
h_i^y,h_i^o\le H_{d_i}.
\]
Consumption and housing separate. Consequently,

\[
c_i^{y,SP}=\chi n_i+x^*,\qquad
c_i^{o,SP}=x^*,\qquad
x^*=\frac12\left(\frac{C^{\mathrm{eq}}}{N}-\frac{\chi}{\nu}\right),
\]
and

\[
\boxed{
h_i^{y,SP}=\min\{H_{d_i},\kappa n_i+\alpha/\lambda\},
\qquad
h_i^{o,SP}=\min\{H_{d_i},\gamma/\lambda\}.
}
\tag{1}
\]
Here \(\lambda>0\) clears housing when the stock is below total retained capacity. These are the **full dated optimum’s** housing choices, not a consumption-fixed approximation.

### Financial implementation

For any such reallocation, hold the price path fixed and use the proposed owner adjustments

\[
\Delta a_i'=-P_{t+1}\Delta h_i^y,
\qquad
\Delta a_j^e=-qP_{t+1}\Delta h_j^o.
\]
They preserve, respectively,

\[
\Delta(a_i'+P_{t+1}h_i^y)=0,
\qquad
\Delta(q^{-1}a_j^e+P_{t+1}h_j^o)=0.
\]
The required current transfer to each household is exactly

\[
t_k=\Delta c_k+u_t\Delta h_k.
\]
Thus \(\int t_k\,dk=0\). These identities follow from the original purchase, mortgage-repayment, and estate budgets. Pasted text

Including rental-intermediary financing, the change in aggregate net financial payoffs is

\[
-P_{t+1}\Delta H^{\mathrm{own}}
-P_{t+1}\Delta H^{\mathrm{rent}}=0.
\]
Existing obligations need not be written down; offsetting financial positions implement the changes.

**No choice between internal and external estate recipients is needed:** each estate payment is unchanged. The old owner’s retained physical cap is \(H_O\), not \(\min\{H_O,e/P_{t+1}\}\). The latter would incorrectly reimpose nonnegative old financial saving, which this planner is authorized to relax. Pasted text

## 2. One proposition with primitive sufficient conditions

The following conditions are conservative certificates, but they allow binding physical caps and both old-owner estate regimes.

Define

\[
\begin{gathered}
w_0=y^y+b,\qquad v_0=y^o,\qquad
a=\alpha+\vartheta,\qquad E=1+a,\\
K=1+\gamma+\omega_B,\qquad D=E+\beta K,\\
d_p=1-q+q\tau^p,\qquad d_L=1-\phi+q\tau^p,\\
\ell_R=1,\qquad \ell_O=d_L/d_p,\\
\Gamma=\min\left\{\gamma,\frac{(\gamma+\omega_B)d_p}{1+q\tau^p}\right\}.
\end{gathered}
\tag{2}
\]
Assume bounded endowments, \(w_0\ge\underline w>0\), \(v_0\ge0\), and

\[
\boxed{0\le\tau^p<2,\qquad \phi\ge q.}
\tag{3}
\]
The second restriction gives \(0<\ell_d\le1\). There is **no restriction on \(\beta>0\)** and no requirement that \(\alpha\ge\gamma\).

The following are primitive upper bounds on the rebate and service price:

\[
\bar M_0=\int(w_0+qv_0)dF,\qquad
\bar T=\frac{\tau^p\bar M_0}{(1-q)(2-\tau^p)},\qquad
\bar p=\frac{\nu[\bar M_0+(1+q)\bar T]}{\kappa}.
\tag{4}
\]
Choose a primitive lower price certificate \(\underline p>0\) as described in the short appendix below.

Define the **cross-sectional adult-space resource bound**

\[
\boxed{
\mathcal R_-=
\int\left[
\alpha\min\left\{\frac{w_0}{D},\frac{\underline p H_R}{a}\right\}
+
\min\left\{\underline p H_R,
\frac{\beta\Gamma(w_0+qv_0)}{qD}\right\}
\right]dF.
}
\tag{5}
\]
This uses clipped moments of current and lifetime resources; it does not require high old income for every type.

Finally, put

\[
B_d=
1+\alpha\ell_d+
\vartheta\frac{\chi+\ell_d\bar p\kappa}{\chi+\bar p\kappa},
\qquad
k_d=\frac{D}{B_d}-1.
\tag{6}
\]

### Proposition

Consider any positive stationary competitive equilibrium of the maintained model.

**Individual conclusion.** For tenure \(d\), define the endowment set \(\mathcal S_d\) by the two inequalities

\[
\boxed{
\begin{aligned}
w_0+\bar T
&<
E\ell_d
\min\left\{
\frac{\mathcal R_-}{\alpha+\gamma},
\frac{\underline p H_d}{a}
\right\},\\
qv_0
&>
k_dw_0+(k_d-q)_+\bar T.
\end{aligned}
}
\tag{I}
\]
Every current young household whose endowments belong to \(\mathcal S_d\) and whose equilibrium tenure is \(d\):

- is strictly below its market housing cap;

- has a strictly positive multiplier on its young financing restriction;

- receives strictly more housing at the full dated optimum:

\[
\boxed{h_i^{y,SP}>h_i^{y,\mathrm{eq}}.}
\]
If \(F(\mathcal S_d)>0\) for either tenure, these households have positive \(Q\)-mass. Logistic tastes give each feasible tenure strictly positive probability at every endowment type.

**Separate aggregate conclusion.** Let \(\widehat w=w_0+\bar T\), and define

\[
G_d(w_0)=
\min\left\{
\left[\underline p H_d-\frac{a\widehat w}{E\ell_d}\right]_+,\;
\alpha\left[
\frac{\mathcal R_-}{\alpha+\gamma}
-\frac{\widehat w}{E\ell_d}
\right]
\right\}.
\tag{7}
\]
If the additional distributional inequality

\[
\boxed{
\mathcal G\equiv
\int\min\{G_R(w_0),G_O(w_0)\}\,dF>0,
}
\tag{A}
\]
holds, then

\[
\boxed{
H_Y^{SP}-H_Y^{\mathrm{eq}}
\ge\frac{N}{p}\mathcal G>0.
}
\tag{8}
\]
Condition (A) allows some young households to lose housing: their potentially negative contributions are explicitly subtracted. It does not infer aggregation from the individual conclusion.

The conclusions hold in **every positive reference equilibrium** satisfying these primitive restrictions. This is not an equilibrium-existence or uniqueness theorem.

## 3. Proof

### Market bounds, including binding caps

Fix an equilibrium and a conditional tenure. Write

\[
w=w_0+T,\qquad v=v_0+T,\qquad M=w+qv,
\]
and let \(\Lambda,\mu,\eta_y\) be the lifetime-budget, young-finance, and young-cap multipliers. With \(x=c-\chi n\) and \(s=h-\kappa n\),

\[
\frac1x=\Lambda+\mu,\qquad
\frac{\alpha}{s}=\Lambda p+\mu\ell_dp+\eta_y,\qquad
\frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s}.
\tag{9}
\]
Multiplying the household first-order conditions by quantities and summing gives

\[
\Lambda M+\mu w\le D.
\tag{10}
\]
The omitted terms are nonnegative physical-cap contributions; the homogeneous estate restriction contributes zero. In particular,

\[
x\ge w/D,\qquad \Lambda\le D/M.
\]
Because \(\ell_d\le1\), the effective current housing price

\[
\rho=\frac{\alpha x}{s}
\]
satisfies \(\rho\ge\ell_dp\). Combining this with the fertility condition and the cash constraint gives

\[
\boxed{
ps\le\frac{\alpha w}{E\ell_d},
\qquad
ph<\frac{aw}{E\ell_d}.
}
\tag{11}
\]
For example, the fertility condition implies

\[
ax=\rho h+\chi n,
\]
so \(\ell_dph<ac\), which proves the second inequality. For the first, \(x\ge\ell_dps/\alpha\) and

\[
n\ge\frac{\vartheta\ell_dps}
{\alpha(\chi+\ell_dp\kappa)};
\]
substitution into \(c+\ell_dph\le w\) proves the claim.

These bounds also establish a lower resource bound. If young housing is uncapped, \(\rho\le p\), hence \(ps\ge\alpha w/D\). If it is capped, the fertility condition gives \(s>\alpha H_d/a\). Therefore

\[
ps\ge\alpha\min\{w/D,pH_d/a\}.
\tag{12}
\]
For old housing, define \(\Gamma_R=\gamma\) and \(\Gamma_O=\Gamma\). In either uncapped old-estate regime,

\[
c^o=z/K,\qquad ph^o=\Gamma_d z/K.
\]
The minimum defining \(\Gamma\) explicitly includes zero old financial saving; it is not an assumption of a slack estate restriction. Pasted text

Let \(m=1/c^o\). Since \(\Lambda=\beta m/q\), (10) gives

\[
m\le\frac{qD}{\beta M}.
\]
If old housing is uncapped, \(ph^o=\Gamma_d/m\); if capped, \(ph^o=pH_d\). Consequently,

\[
ph^o\ge
\min\left\{pH_d,\frac{\beta\Gamma_dM}{qD}\right\}.
\tag{13}
\]
Using \(p\ge\underline p\), \(w\ge w_0\), \(M\ge w_0+qv_0\), and taking the lower envelope across tenures in (12)–(13), we obtain

\[
\boxed{
p\left(\frac{\bar H}{N}-\frac{\kappa}{\nu}\right)
=p\int(s_i+h_i^o)dQ\ge\mathcal R_-.
}
\tag{14}
\]

### Comparing the same household with the full optimum

From the capped planner formula (1),

\[
\frac{\bar H}{N}-\frac{\kappa}{\nu}
\le\frac{\alpha+\gamma}{\lambda}.
\]
Thus

\[
\frac{p\alpha}{\lambda}
\ge\frac{\alpha\mathcal R_-}{\alpha+\gamma}.
\tag{15}
\]
The first inequality in (I), together with (11), proves both

\[
h_i^{y,\mathrm{eq}}<H_d,
\qquad
ps_i^{\mathrm{eq}}
<
\frac{\alpha\mathcal R_-}{\alpha+\gamma}
\le\frac{p\alpha}{\lambda}.
\]
Both arguments of the minimum defining \(h_i^{y,SP}\) therefore exceed \(h_i^{y,\mathrm{eq}}\). This proves the individual full-optimum comparison.

### Establishing binding finance

Suppose instead that \(\mu=0\) for one of these households. Its young cap is already known to be slack. Equation (9) then implies

\[
s=\frac{\alpha x}{p},\qquad
n=\frac{\vartheta x}{\chi+p\kappa},
\]
and its cash expenditure is

\[
c+\ell_dph
=
x\left[
1+\alpha\ell_d+
\vartheta\frac{\chi+\ell_dp\kappa}{\chi+p\kappa}
\right].
\]
The bracket decreases with \(p\), so it is at least \(B_d\). Moreover, (10) with \(\mu=0\) gives \(x\ge M/D\). Cash feasibility therefore requires

\[
w\ge B_dM/D,
\quad\text{equivalently}\quad qv\le k_dw.
\]
But the second inequality in (I) guarantees \(qv>k_dw\) for every \(T\in[0,\bar T]\). Contradiction. Hence \(\mu>0\).

Notice that this argument **does not require the old housing cap or old estate restriction to be slack**.

### Aggregation with losses included

Exactly,

\[
p(h_i^{y,SP}-h_i^{y,\mathrm{eq}})
=
\min\left\{
p(H_d-h_i^{y,\mathrm{eq}}),\;
\frac{p\alpha}{\lambda}-ps_i^{\mathrm{eq}}
\right\}.
\]
Equations (11) and (15) bound this below by \(G_d(w_{0i})\). Integrating, and using \(G_{d_i}\ge\min\{G_R,G_O\}\), proves (8). ∎

## 4. Interpretation and a precise obstruction

The first inequality in (I) identifies **cash-poor households relative to the economy’s available adult space**, while also guaranteeing room below their retained tenure cap. The second establishes that future resources are sufficiently large relative to current cash to make finance strictly restrictive.

For renters, \(B_R=E\), so the finance test without taxes is simply

\[
\frac{v_0}{w_0}>\frac{\beta K}{qE}.
\]
At \(\beta=q\), this is \(K/E\), rather than the conservative \(K/\gamma\) needed to force old housing above young housing type by type. **These tests establish different facts:** the lower threshold establishes binding finance; the cash and distributional conditions establish the full housing comparison.

Binding finance alone cannot replace those additional conditions. Here is an analytical counterfamily.

Take the explicitly tax-free subcase

\[
\phi=q,\qquad \beta=q,\qquad
\omega_B(1-q)<q\gamma,
\]
with finite caps large enough to be inactive. Then

\[
\Gamma=(\gamma+\omega_B)(1-q)<\gamma.
\]
Conditional tenure-value differences are constant across endowments, so logistic tastes can generate any constant owner share \(\pi\in(0,1)\). Put

\[
\bar\Gamma=(1-\pi)\gamma+\pi\Gamma,
\]
and choose

\[
\boxed{\frac KE<r<\frac{\gamma K}{E\bar\Gamma}.}
\tag{16}
\]
Let current cash be genuinely heterogeneous and let \(v_0=rw_0\).

Every young household has strictly binding finance, because \(r>K/E\). Nevertheless,

\[
ps_i=\frac{\alpha w_{0i}}E,
\qquad
p\bar h^o=\frac{\bar\Gamma r}{K}\,\overline{w_0}.
\]
The uncapped full dated optimum consequently satisfies

\[
\boxed{
H_Y^{SP}-H_Y^{\mathrm{eq}}
=
\frac{N\alpha\overline{w_0}}{p(\alpha+\gamma)}
\left(\frac{\bar\Gamma r}{K}-\frac{\gamma}{E}\right)<0.
}
\tag{17}
\]
This is compatible with \(\alpha\ge\gamma\). It is also a stationary family: choose positive

\[
\chi<\frac{\vartheta\nu\overline{w_0}}E,
\qquad
p=\frac{\vartheta\nu\overline{w_0}/E-\chi}{\kappa},
\]
and determine \(N\) from housing clearing. No numerical example or neighborhood argument is involved.

The obstruction is substantive: **old owners’ inability to borrow against their estate changes their housing demand sufficiently that relaxing finance need not move housing toward the young.**

Conversely, equal-weight redistribution can produce individual housing gains without any financing distortion. In the uncapped, slack-finance, positive-financial-estate benchmark with \(\beta=q\),

\[
h_i^{y,SP}-h_i^{y,\mathrm{eq}}
=\frac{\alpha}{p}(\bar x-x_i),
\qquad
H_Y^{SP}-H_Y^{\mathrm{eq}}=0.
\]
Thus the proposition identifies households that **both are constrained and gain housing**. It does not attribute their entire gain exclusively to the borrowing constraint.

## Technical appendix: the primitive price certificate

The bounds \(T\le\bar T\) and \(p<\bar p\) follow by aggregating the original lifetime budget and using replacement fertility; the packet provides the corresponding rebate and price bounds. Pasted text

A sharper lower bound than one based only on minimum wealth can use the distribution. For each \(w_0\), let \(r_c(w_0)\in(0,H_R/\kappa)\) solve the elementary quadratic equation

\[
\frac{\vartheta}{r_c}
=
\frac{\chi D}{w_0}
+\frac{\alpha\kappa}{H_R-\kappa r_c}.
\]
Choose any \(\underline p>0\) satisfying

\[
\boxed{
\int
\min\left\{
\frac{\vartheta w_0}{D(\chi+\underline p\kappa)},
r_c(w_0)
\right\}dF
>\frac1\nu.
}
\tag{18}
\]
To verify it, uncapped young households satisfy

\[
n\ge\frac{\vartheta w_0}{D(\chi+p\kappa)};
\]
capped young households satisfy \(n\ge r_c(w_0)\), because \(x\ge w_0/D\) and their cap is at least \(H_R\). If \(p\le\underline p\), (18) would imply fertility above replacement, contradicting stationarity.

All quantities in (I) and (A) are therefore functions of primitives and the chosen scalar certificate. Failure of this sufficient certificate—or of (A)—makes these bounds inconclusive; it does not establish the opposite allocation.

**The resulting scope is precise:** the individual theorem identifies a positive-mass group without assuming its constraint pattern; the aggregate theorem separately accounts for possible young losers; and both retain the original mortgage, positive child costs, persistent tenure, and physical caps.
```

### File: docs/style/econ_writing_style_guide.md
```md
# Economics Writing Style Guide

Tool-agnostic version of the `writing-econ-papers` Claude skill. Paste or attach this
file in any LLM prompt (ChatGPT, Codex, Claude) that drafts or revises paper-facing
text: model sections, abstracts, theory notes, slide decks.

Audience for this guide: an LLM (or human) drafting or revising paper-facing text in
quantitative macro, housing, spatial, urban, or household finance. Distilled from real
revision rounds on theory drafts and decks; every rule below exists because its violation
was produced by a strong model and rejected by the author.

## 1. Section architecture (the Menzio rule)

Model sections follow this order, and the boundary between "environment" and "problems"
is strict:

1. **Environment** — labeled primitive paragraphs, in bold lead style:
   - **Agents.** Who exists, masses, types and their distributions.
   - **Preferences.** Utility functions for every agent, stated as primitives.
     Interpret parameters in words. If a structural identity makes the economics vivid
     (e.g. substituting residuals into a budget to show "the full price of a child"),
     show it here.
   - **Endowments / Technology / Markets.** What is traded, at what prices, what clears
     against what. Choose the simplest market structure that carries the point — one
     clearing condition if possible. Richer structure (multiple stocks, construction,
     segment-specific supply) is a footnoted extension, not the baseline.
   - **Institutions.** Financing limits, tax rules, assessment rules — as rules of the
     game, not as constraints inside someone's problem yet.
   - **Entry / demography** if relevant.
   - NO maximization problems, Bellman equations, FOCs, or value functions here.
2. **One section per agent's problem.** "The young household's problem" contains the
   problem, taking prices as given, and nothing that belongs in the environment.
3. **Equilibrium** — a numbered Definition listing the objects first, then the
   conditions. Lead the section with one paragraph saying what the equilibrium object is
   and what it does NOT pin down (e.g. "clearing pins aggregates, not assignment").
4. **Planner / efficiency** — open with the question the planner answers.
5. **Results, then policy.** Comparative statics before policy formulas; honest scoping
   sentences on what is local/fixed-price/conditional.

If the paper has a toy model and a full model (Coven et al. structure), the toy comes
first under its own heading, then "Setup", then the framework.

## 2. Primitives before derived objects

Preserve the author's existing notation, utility definitions, value-function
presentation, and variable names when amending a draft. Permission to develop
economic results or simplify prose does not authorize changing those choices.
If a correction requires a convention to change, explain the specific reason
and propose the change separately; otherwise work within the existing convention.

A derived object (a multiplier, gap, wedge, sufficient statistic) may not do work in the
text before it has been expressed in primitives and tied to observables.

**Bad (rejected in a real round):**

> Define the goods-equivalent value of relaxing the rental size limit by
> $\zeta_i^R=\theta_i^R/\lambda_i^R\ge0$.

**Good (accepted):**

> Converting the cap multiplier into consumption goods, $\zeta_i^R=\theta_i^R/\lambda_i^R$,
> the renter's marginal value of housing is $\alpha c_i^R/s_i^R = q+\zeta_i^R$.
> A renter against the size cap values one more unit of space above the rent it would
> command, and the gap is not an abstract multiplier: it is computable from the observed
> bundle as $\zeta_i^R=\alpha c_i^R/s_i^R-q$. An unconstrained renter has $\zeta_i^R=0$.

Corollaries of the rule:
- Convert every multiplier to goods units before interpreting it.
- Decompose composite gaps by the primitive constraint that generates each term
  ("the collateral margin... the debt-service margin... the menu").
- Name objects by their economics: "incumbent tax discount", "subsidy to staying put",
  "effective family-housing cap" — never decorative or purely technical labels.

## 3. Words around math

- Every display gets a one-line lead-in ending in a colon or verb, and a one-sentence
  interpretation immediately after. No orphan displays.
- Mechanisms are stated in numbered plain prose before any notation:
  "The mechanism has two parts. First, ... Second, ..."
- Facts first. Open sections and abstracts with the economic fact or question, not with
  notation or with the literature.
- Honest scoping is part of the result: say "holding prices, tenure, and active sets
  fixed", say what reactivates the omitted margin, say where the general case is handled.
- Do not overclaim: a conditional marginal comparison is not a welfare theorem; a local
  formula is not a GE counterfactual.

## 4. Remarks, footnotes, and apparatus

- Remark environments are for standalone formal caveats only. If it reads like
  commentary, fold it into prose; if it is bookkeeping, make it a footnote.
- Status notes, normalization conventions, and "we will extend this" content belong in
  footnotes, written as economics ("Two supply extensions are deliberately left out of
  the baseline and developed later: ...").
- Keep equation labels stable across revisions. Avoid notation collisions (x vs X, y vs
  Y); if a local construction needs symbols that collide, scope them explicitly in a
  footnote.
- Use bars for cross-sectional averages and integrals for aggregation. Reserve
  expectation notation for uncertainty; do not use it merely to average households.
- Assumptions stated where used; iff results get both directions discussed; proofs
  include the nontrivial computation (no "it is easy to see").

## 5. Citation discipline

- Never describe another paper's model from memory. Open the PDF and read the model
  section, the equilibrium definition, and the institutional details before writing
  "as in X (year)".
- "As in X" must name exactly which feature is borrowed ("all floorspace clears in a
  single market, as in X"), and a footnote states where you differ ("X impose a minimum
  owned size; we cap rental size instead").
- Cite the version you actually read; update years and working-paper numbers.
- If a mechanism you use is one that a cited paper explicitly abstracts from, say so —
  that is positioning gold, not a problem.

## 6. Slides

For Beamer decks and presentation-facing text, read and follow
`docs/style/econ_presentation_style_guide.md`. Its defaults are intentionally conventional:
one idea per frame, noun-phrase titles, full-width exposition, sparse figures, and standard
production-paper tables. In particular, a calibration slide defaults to
`Moment | Target | Model`; diagnostic columns appear only when the author requests them.

## 7. Workflow checks (LaTeX)

- Compile twice; require zero errors, zero undefined references, zero overfull boxes
  above ~2pt. Render pages to images and look at them before declaring victory —
  overlapping TikZ labels and clipped boxes do not show up in logs.
- After structural edits, grep for orphaned `\ref`/`\eqref` of deleted labels and for
  leftover symbols from the old structure.
- Commit a checkpoint before and after any restructure.
- Before delivery, reread the entire main text against Sections 2 and 9.
  Check against the author's reference draft as well as the previous version;
  a recent assistant rewrite is not evidence that its conventions were approved.
  Correct avoidable notation and prose drift before handing the draft back.

## 8. Red flags → corrections

| Thought | Reality |
|---|---|
| "I remember what that paper does" | Read the model section first; quote the clearing condition. |
| "This object is standard, no need to unpack" | Express it in primitives and observables first. |
| "A remark environment is safer" | Fold into prose or a footnote. |
| "The title can carry the claim" | Conventional title; the claim goes in the text. |
| "More structure shows rigor" | Simplest market structure that carries the point; extensions in footnotes. |
| "The math speaks for itself" | One sentence of economics before and after every display. |
| "Creative presentation helps the reader" | More normal than you think. One idea per frame. |
| "The preferences belong with the problem" | Preferences are primitives; they go in the environment. |

## 9. Plain explanations and simple proofs

Tommaso reaffirmed on September 4 that simple theory should read simply.
Earlier instructions in the `Reconstruct theory and draft history` conversation
were to “make simple things simple,” make proofs more verbal and clear, and
avoid announcing a complex proof strategy for an elementary exercise.

- Use the prose of **Guido Menzio and Raquel Fernández** as writing references,
  as Tommaso requested. Revisit relevant passages when drafting or revising;
  take cues from the actual prose, not just the authors' names.
- Present this theory as a **simple, illustrative exercise**. Keep the claims,
  explanations, notation, and proofs proportionate to that purpose. Do not
  make it sound more ambitious than it is or add complexity to signal rigor.
- Start with what the result says about households, housing, or children.
- Use a technical term only when it saves a necessary distinction. Explain
  “the same constraints continue to bind” before using “fixed active set.”
- Give the short argument directly. Avoid announcing a framework, strategy,
  architecture, or mechanism map for a few lines of algebra.
- State the scope once where it matters. Repeated lists of what a result does
  not prove make the argument harder to follow.
- Keep workflow, verification receipts, issue IDs, and implementation status
  in the work record. They do not belong in proposed paper prose.
- Prefer a compact appendix that explains the nontrivial steps over a catalogue
  of near-identical cases. Preserve every condition needed for correctness.

Writing references consulted for the illustrative theory:

- Guido Menzio, [*A Theory of Partially Directed Search*](https://web-facstaff.sas.upenn.edu/~gmenzio/linkies/PDS.pdf),
  JPE 2007, Section II, pp. 751–754: agents and timing precede the equilibrium
  definition, whose conditions are then explained in words.
- Raquel Fernández and Richard Rogerson,
  [*Income Distribution, Communities, and the Quality of Public Education*](https://drive.google.com/file/d/1QQA7trmlvZeW5a0HFZ2OMH9vTsPm8g8T/view),
  QJE 1996, pp. 137–138 and 140–141: the question and chosen simplifications
  are stated before the formal analysis.

Use these as references for exposition. Their economic assumptions and the
complexity of their results are not a template for this illustrative exercise.
```

### File: docs/style/econ_presentation_style_guide.md
```md
# Economics Presentation Style Guide

Use this guide for paper talks, seminar decks, advisor updates, and quantitative-result
slides in this project. The intended reader is technically sophisticated but does not know
the codebase or the sequence of internal experiments.

## 1. Default look

- Follow the visual language of a conventional economics job-market or production-paper
  presentation: restrained, flat, and easy to scan.
- Use a single full-width composition for explanatory text and tables. Two columns are
  appropriate for two graphs or a genuine side-by-side comparison, not for prose plus a
  second wall of prose.
- Use conventional noun-phrase titles: `Toy economy`, `Competitive equilibrium`,
  `2023 calibration`, `Population dynamics`. Put the argument in the frame body.
- Use academic language even in advisor updates. Avoid deictic titles such as `Today`,
  `Now`, `Current status`, `What we established`, or `Where this leaves us` unless the
  date or sequence is itself substantive.
- Use `booktabs` tables without vertical rules, boxes, cards, dashboard elements, or
  decorative labels.
- Keep the main deck readable at presentation distance. Move robustness, derivations,
  implementation details, and alternative diagrams to the appendix.

## 2. Frame architecture

- One idea per frame.
- Keep bullets to one line whenever possible. If a bullet needs a paragraph, shorten it or
  split the argument across frames.
- A standard frame contains at most: one short framing sentence, one display/table/figure,
  and one interpretation sentence.
- Define each parameter or variable in plain English where it first appears. Do not use a
  generic state such as `$x$` as a substitute for naming the economically relevant states.
- Describe what the model does. State limitations only when they change the interpretation
  of the result; do not fill slides with lists of what the exercise does not do.
- Never expose workflow language such as `current reading`, `what we tried`, `the model
  shown to X`, `implementation status`, or internal experiment names. Translate it into the
  economic object an audience needs.
- Keep issue-ledger language out of outward-facing slides: `provisional`, `pending`,
  `not final`, `needs review`, run readiness, and implementation status belong in project
  notes. Include them on a visible slide only when the author explicitly requests them.

## 3. Calibration and quantitative tables

- The default calibration table is exactly:

  `Moment | Target | Model`

- Do not add `Gap`, `Weight`, `Loss`, `Loss contribution`, standardized residuals, parameter
  counts, or objective diagnostics unless the author explicitly asks for them.
- Do not replace a conventional calibration table with a target-fit bar chart.
- If two or three misses deserve attention, bold those rows and name them in one sentence
  below the table. Do not add a diagnostic side panel.
- A parameter table defaults to `Parameter | Economic interpretation | Value`. Bounds,
  transforms, Jacobians, and near-bound diagnostics belong in the appendix unless requested.
- A policy table may show `Baseline | Reform | Change` when the change is the result being
  presented. Use percent changes for levels and percentage points for rates.
- Preserve units in row labels or column headers. Round for reading, while retaining enough
  precision to distinguish the compared objects.

## 4. Model and equilibrium slides

- Introduce the simplified economy before presenting its results: agents, timing, choices,
  prices, and market clearing.
- Define a steady state as a constant allocation, population/cohort structure, and price
  satisfying household optimality, demographic reproduction, and market clearing.
- When presenting a shock, separate three objects when they matter:
  1. the initial steady state;
  2. the impact or partial-equilibrium response;
  3. the demographic adjustment and new equilibrium.
- Use equations selectively. Every display gets one sentence before it and one economic
  interpretation after it. Derivations and proofs go to the appendix.

## 5. Figures and dynamics

- Use a figure only when it carries a mechanism or result more clearly than a table or one
  equation. Avoid diagrams of workflow, pipelines, or internal architecture.
- Label axes and equilibrium points economically. The primitive shock and the direction of
  movement should be visible without narration.
- For dynamics, show the impact point, any population/price trough, and the long-run point
  when those are the economic content. Choose illustrative parameters that make distinct
  stages visually distinct, without changing the qualitative economics.
- Keep the preferred mechanism graph in the main deck. Put alternative visualizations in
  the appendix and link to them only when useful during discussion.

## 6. Quantitative claims and caveats

- Distinguish targets, imposed inputs, untargeted validation series, and model outcomes.
- A path fixed to observed totals or age shares is an input, not a model fit. Say this once,
  plainly, where the path is introduced.
- Do not call a diagnostic closure a forecast, a finite transition a new steady state, or an
  impact calculation a complete policy experiment.
- Keep caveats short and economically substantive. Source notes belong in a small footer;
  code hashes, run identifiers, residual tolerances, and solver details do not belong in the
  visible main deck.
- Source footers should be standard citations or dataset names, not explanatory prose. If a
  measurement mismatch changes the economics, explain it once in the frame body.

## 7. Final check

Before delivering a deck:

1. Check every main frame against this guide.
2. Compile twice and require no errors, undefined references, or material overfull boxes.
3. Render the changed frames at full size and inspect them visually.
4. Confirm tables reproduce their source data and use consistent rounding.
5. Confirm the main deck contains economic content rather than a technical progress report.
```
