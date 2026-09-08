# Full planner with fixed fertility: settle the benchmark before the theorem

You are reviewing a simple, illustrative two-age OLG model of housing and
fertility. We need a precise analytical foundation for a paper and a few theory
slides. Previous iterations have produced conditional housing variations but
have not specified the full planner satisfactorily. Resolve that problem first.
Do not simply improve the wording of the attached proposition.

The author wants to compare a competitive steady state with a planner's
stationary allocation, initially holding fertility fixed. The planner chooses
consumption as well as housing and the other components of the allocation.
The hoped-for economic conclusion is that borrowing limits and differences in
the housing available to renters and owners can leave too little housing with
young households. This is a hypothesis to establish or qualify, not an answer
you must obtain. We have spent too many iterations selecting assumptions to
rescue a predetermined conclusion.

Please give one recommended formulation and at most one substantive welfare
alternative. The author is willing to compare welfare weights and generally
prefers conventions established in dynastic OLG. Explain what that preference
does and does not imply for this particular household model. Preserve the model
unless a change is necessary; identify every proposed change explicitly.

## 1. Author instructions and scope

- First compare stationary allocations with the SAME positive masses
  \(Y=O=N\), housing stock, income distribution, and aggregate resources.
  Fix the fertility of each household type at its competitive choice, not just
  average fertility. A type includes endowments and the ownership taste draw.
  If the planner changes tenure, its fertility remains the frozen individual
  choice. This is a conditional welfare exercise, not a claim that fertility
  would still be privately optimal after reallocation.
- Stationarity requires \(\nu\bar n=1\), but this does not determine \(N\).
  Take \(N\) from the reference competitive steady state. Do not let the
  planner get more housing per household by choosing a smaller population.
- The direct planner chooses the FULL feasible allocation: young and old
  consumption, housing, permitted tenure assignments, saving and estates (or
  their correctly defined real counterparts). It can redistribute resources
  and relax individual financing restrictions. Physical housing limits and
  tenure persistence across a household's two ages remain. Assume full
  information for this direct benchmark, including realized ownership tastes.
  Fixing tenure or consumption in a proof variation is permissible; defining
  the entire planner by that restricted variation is not.
- Welfare weights must be justified independently of the desired housing
  direction. Distinguish failure to maximize the chosen utilitarian objective
  from Pareto inefficiency and from constrained inefficiency.
- The two-period structure, heterogeneous young wealth and young/old income,
  housing size limits by tenure, current income available for purchase, standard
  mortgage repayment, positive child goods and space costs, and warm-glow
  estate preferences are the maintained starting point. Retain the author's
  notation. Do not remove heterogeneity, assume a very low household discount
  factor, add liquidation costs, or replace the mortgage to manufacture a sign.
- THIS RUN DOES NOT need a transfer-policy theorem, endogenous-fertility
  welfare, a demographic transition proof, numerical examples, calibration,
  plots, or slides. Those come after this benchmark is settled. Explain how a
  stationary comparison differs from welfare along an attainable transition,
  but do not spend the run solving the transition.

## 2. What the attachment establishes, and what it does not

The attached `simplified_olg_utilitarian.tex` is the exact current discussion
note. Its environment and household equations are the maintained specification.
Its welfare section and propositions use an earlier, dated remaining-lifetime
criterion; they are material to reassess, not an authoritative definition of
the new full planner. Later transfer and fertility results are background only.

Each household lives young and old. Lifetime utility, including the tenure
taste once, is
\[
u^y(c,h,n)+\xi\mathbf1\{\mathrm{owner}\}+\beta u^o(c^2,h^2,e),
\]
\[
u^y=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^2+\gamma\log h^2+\omega_B\log e.
\]
The note's superscript \(2\) denotes old-age quantities; \(y,o\) label utilities.
The endowment triple \((y^y,b,y^o)\) has distribution \(F\), and the independent
ownership taste \(\xi\) is logistic. Future income is known. Use integrals or
bars for aggregation, not expectations suggesting aggregate uncertainty.

The estate enters parental utility directly. There is NO continuation utility
of children and NO equation making entrant wealth \(b\) inherit the estate.
Thus this is currently warm-glow OLG, not a Barro–Becker dynasty. If genuine
dynastic altruism is essential to a suggested benchmark, show the extra
preference and inheritance equations and label this a different model. Do not
silently replace \(\omega_B\log e\) by descendants' utility.

Goods and bonds trade with the outside world at bond price \(q\in(0,1)\).
At stationarity let
\[
p=(1-q+q\tau^p)P,\quad L=(1-\phi+q\tau^p)P,
\quad w=y^y+b+T,\quad v=y^o+T,
\]
\[
z=a'+Ph+v,\qquad K=1+\gamma+\omega_B.
\]
Here \(p\) is the housing service cost, \(L\) the cash requirement per owner
unit, and \(z\) resources on entering old age. The young owner's constraints
reduce exactly to
\[
c+ph+qz=w+qv,\qquad c+Lh\le w.
\]
Thus current income CAN finance the down payment. The original constraints,
renter problems, estate equation, rebates and housing clearing are attached.

## 3. Resolve welfare weights first

Use primary literature to distinguish the objects below; do not call all of
them interchangeable versions of utilitarian welfare.

1. A stationary cohort's lifetime welfare, integrating
   \(u^y+\xi\mathbf1\{\mathrm{owner}\}+\beta u^o\) across the cohort.
2. Welfare over an infinite sequence of cohort lifetimes, with explicit social
   generation weights, a social discount factor distinct from \(\beta\), and
   an explicit treatment of the initial old. Derive the relative age weights
   in its dated first-order conditions. Explain why maximizing welfare only
   over stationary allocations need not give the steady state approached by
   this dynamic planner.
3. Welfare of initial dynasties when parental utility recursively includes
   descendants. Explain separately the social weights across dynasties and
   the private altruism weights within each dynasty.

Recommend ONE primary criterion for our immediate stationary comparison and
ONE useful comparator if needed. Tell us whether the recommendation preserves
the household model. Do not introduce a new dynasty merely because the author
asked for literature-consistent weighting. With fixed population, state when
total and average welfare differ only by a constant normalization.

The current note gives weight one to each living household's remaining utility:
young lifetime utility plus current old utility. At stationarity its old flow
terms consequently have combined coefficient \(1+\beta\), whereas a single
cohort lifetime criterion has coefficient \(\beta\). Neither criterion should
be conflated with equal weights on every cohort along a whole transition.

Useful PRIMARY starting points, to verify rather than invoke as authority:

- Becker and Barro, *A Reformulation of the Economic Theory of Fertility*,
  QJE (1988); working-paper version: https://www.nber.org/papers/w1793 and
  https://www.nber.org/system/files/working_papers/w1793/w1793.pdf.
  Parental altruism generates a dynastic utility formulation.
- Golosov, Jones and Tertilt, *Efficiency with Endogenous Population Growth*,
  Econometrica (2007):
  https://tertilt.vwl.uni-mannheim.de/research/optimality_Econometrica.pdf.
  Section 3.4, especially Results 1–2 around printed pp. 1055–1056, distinguishes
  planning problems with weights on potential people versus initial agents.
  These results do not select a unique equal-weight social objective for us.
- Farhi and Werning, *Inequality and Social Discounting*, JPE (2007):
  https://www.journals.uchicago.edu/doi/10.1086/518741;
  author manuscript: https://web.mit.edu/iwerning/Public/inequality_social_screen_old.pdf.
  Check its separation of parental altruism and social concern for descendants.

Give precise equation/section references for the convention you actually use.
Report a lack of a unique convention honestly. We need a small, relevant
literature foundation, not a literature review.

## 4. Specify full feasibility instead of assuming it

An independent accounting check found that the current note's household
budgets alone do not uniquely specify the price-free planner. Please resolve:

- Who receives estates and who supplies the exogenous entrant wealth \(b\)?
  Internal bequests are transfers; they are not lost goods. Outside recipients
  imply external outflows. Identify which closure preserves the note.
- Who owns rental intermediaries and housing titles outside living owners?
  Their pricing condition alone does not specify ownership or their resource
  contribution. Preserve the aggregate housing stock; do not count purchases
  of existing housing as newly produced goods.
- What is the estate object in a direct allocation? The owner's definition
  \(e=q^{-1}a^e+P_{t+1}h^2\) is wealth including house resale value; the
  renter's estate is \(e=q^{-1}a^e\). A planner
  cannot choose arbitrary nominal house prices to create utility/resources.
  Explain a consistent real delivery/valuation rule. If reference-price
  valuation is retained, say exactly what that means for the benchmark.
- Which individual financial constraints does the planner relax? The owner's
  \(e\ge P_{t+1}h^2\) is equivalent to nonnegative financial saving when old;
  it is not itself a physical housing limit. Do not retain it as physical
  feasibility without explanation.
- What fixes aggregate external wealth? Do not give the planner an arbitrary
  asset endowment. Individual saving and wealth distribution may change while
  the aggregate endowment is held fixed. Distinguish that restriction from
  requiring an attainable path from a common inherited state.

For orientation ONLY, if all domestic claims can be consolidated, a candidate
dated ledger is
\[
C_t+qB_{t+1}^{\mathrm{ext}}+X_t
=Y_t^g+B_t^{\mathrm{ext}}+I_t,
\]
where \(C_t\) is total nondurable expenditure (already including child goods),
\(Y_t^g\) goods income, \(B_t^{\mathrm{ext}}\) external bond payoffs available
at \(t\), and \(I_t,X_t\) actual outside inflows and outflows. This is NOT a
closed specification until recipients, ownership and timing are fixed. At
stationarity it implies
\(C+X=Y^g+I+(1-q)B^{\mathrm{ext}}\).
Freely choosing endowed \(B^{\mathrm{ext}}\) can make welfare unbounded.
With domestic taxpayers and rebate recipients, fully rebated property taxes
are internal transfers, not resource losses. Any foreign titleholders' tax
payments must be recorded consistently with outside rental/title cash flows.

Recommend the smallest explicit completion consistent with the household
model, label its added assumptions, and then write the full planner with all
its controls, objective, feasibility constraints, and fixed objects. If no
such completion preserves the claimed economy, identify the precise conflict
and the smallest proposed revision. Do not hide it in a theorem's hypotheses.

## 5. Analytical result and a diagnostic that needs interpretation

On the uncapped old-owner branch with positive financial estate,
\[
c^2=z/K,\quad h^2=\gamma z/(Kp),\quad e=\omega_Bz/(Kq).
\]
This branch requires
\[
\omega_B(1-q+q\tau^p)>q\gamma.
\]
It is a substantive restriction, not a harmless regularity condition.
With \(m=1/c^2\), \(s=h-\kappa n\), an uncapped young owner satisfies
\[
\alpha/s=(\beta/q)pm+\mu L,\qquad \gamma/h^2=pm,
\]
where \(\mu\ge0\) is the multiplier on \(c+Lh\le w\).

For the stationary cohort lifetime criterion, a candidate permanent housing
variation across otherwise identical young and old owner types gives
\[
\Delta h=\varepsilon,\quad\Delta h^2=-\varepsilon,
\quad\Delta a'=-P\varepsilon,\quad\Delta a^e=qP\varepsilon,
\]
with internal age transfers \(+p\varepsilon\) to young and
\(-p\varepsilon\) to old. Each entering old household then inherits an extra
\(\varepsilon\) housing and \(-P\varepsilon\) financial wealth. At reference
prices, entering-old total resources, estates and both ages' consumption stay
unchanged; \(q\Delta a'+\Delta a^e=0\). Repeated cohort by cohort, this is
a stationary candidate variation. It relaxes private financing when needed.

The resulting derivative appears to be
\[
\frac{\alpha}{s}-\beta\frac{\gamma}{h^2}
=\beta(1/q-1)pm+\mu L.
\]
Independently verify the algebra AND feasibility under your proposed complete
closure. It is positive for \(q<1\) even if \(\mu=0\). Interpret that fact:
it cannot by itself attribute the welfare difference to a binding mortgage.
Nor does it account for the initial old's transition welfare.

We need more than another sufficient condition saying that marginal utilities
are unequal. Seek a simple statement linking the model's parameters and
competitive choices to failure of the chosen FULL planner optimality
conditions. Wherever possible give interpretable primitive inequalities;
the attached note has existing income/capacity restrictions to inspect, not
automatically reuse under a new benchmark. Check the regular regime with zero
financial estate as well. If a positive-estate or cap restriction is necessary,
say how much of the result it limits and why.

Most importantly distinguish:
(a) an improving feasible housing variation; (b) the competitive allocation
not solving the full planner problem; (c) the full optimum allocating MORE
aggregate housing to young households after consumption, estates and tenure
are jointly optimized. Prove (c) if possible, but do not infer it from (a).
Separate ordinary redistribution under the welfare weights from the extra
effect of financing and tenure restrictions. If equal utility weights already
generate a gap without those restrictions, say so clearly. A frictionless
comparison with the SAME objective and resources is a useful diagnostic, not
an invitation to choose convenient weights after observing the answer.

## 6. Required output

Give a short verdict, the recommended welfare criterion and one comparator,
the fully specified planner, and one main analytical proposition with proof.
Use a compact appendix only for essential algebra or accounting. End with the
few author decisions that actually remain, separating normative choices from
missing accounting assumptions and changes to household behavior. Keep the
main answer roughly within 2,500 words, with essential proof detail beyond
that if necessary. Use plain economics prose and explicit equations. Take the
reasoning time needed to check the accounting and theorem; verbosity is not
the objective. No numerical-point-plus-open-neighborhood proof, no claims of
universal inefficiency without proof, and no promise of an unproved fertility
effect. A clear limitation is preferable to another superficially complete
formulation that requires restarting tomorrow.
