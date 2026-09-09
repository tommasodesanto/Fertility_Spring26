# Independent second review of the dated housing and fertility results

September 9, 2026 UTC. Bounded analytical review of
`overnight_housing_independent.md` and Section 6 of
`overnight_fertility_independent.md`. No numerical run, browser review, model
change, or build was used. This review owns only this file.

## Verdicts

| Statement | Verdict | Scope or clarification |
|---|---|---|
| Primitive cap-valid fixed-fertility housing theorem | **Pass** | For the stated dated planner, with both young finance and old nonnegative financial saving relaxed, while current goods, total housing, incumbent continuation opportunities, old net estates, tenure, and physical caps are retained. |
| Analytical compatibility and stationary existence construction | **Pass** | It proves a nonempty parameter family and existence of a positive stationary price, not uniqueness. The intermediate young Euler bound is valid; its derivation is recorded below. |
| Construction with binding planner rental caps | **Pass** | The owner-emulates-renter comparison requires both the stated `phi >= q` condition and the stated strict estate condition. They suffice even when renter caps bind. |
| Exact capped equilibrium counterexample | **Pass, with a scope clarification** | It disproves sufficiency of paired marginal-utility ordering plus a **positive mass** of strictly constrained young households. High-type renters have slack finance, so it is not a counterexample imposing strictly binding finance on every young household. |
| Joint-fertility theorem under the tenure-conditional adult-consumption cushion | **Pass** | Conditional tenure weights and physical caps must agree across the two incumbent cohorts, as stated. This is the joint dated optimum, not private fertility rechoice or a demographic-transition theorem. |
| Conservative primitive certificate for that cushion | **Pass** | Requires a finite upper bound on current resources, as used explicitly in the certificate. It is substantially stronger than necessary. |

I found no substantive algebraic or logical gap in these stated results. Before
incorporation, state the planner's relaxation of **old financial saving as well
as young finance explicitly**. Keeping the old restriction `e >= P h` would
change the planner's housing solution and would not be the problem reviewed
here. State finite aggregate current goods and well-defined utility integrals
as routine regularity conditions if they are not already maintained.

The financial settlement checks after including the rental intermediary.
At the same zero-tax stationary price, fixed young continuation resources and
fixed old estates require, for owners,
\[
\Delta a_i'=-P\Delta h_i^y,\qquad
\Delta a_i^e=-qP\Delta h_i^o.
\]
For renters those two financial changes are zero: their continuation wealth
and estate contain no housing title. Household current net bond purchases
therefore change by \(-qP\Delta H_O\), where \(\Delta H_O\) is the change
in total occupied owner housing across the two ages. The rental intermediary
changes its financial position by \(-qP\Delta H_R\), financing the change
in its housing stock against the corresponding future resale value. Adding
that account gives
\[
\Delta B_{\mathrm{all}}=-qP(\Delta H_O+\Delta H_R)=0.
\]
Each household's remaining current cash adjustment, in either tenure, is
\(\Delta c_i+p\Delta h_i\). Those adjustments sum to zero because
aggregate consumption and occupied housing are fixed. This is a feasibility
check for the direct planner, not a proof that households voluntarily choose
the allocation under the original financial restrictions. The rental account
is essential; applying the owner financial-wealth adjustment to renters would
be incorrect.

## 1. Independent household bounds and the housing theorem

Use the housing report's notation. Let `lambda = beta/(q c^o)` be the young
lifetime-budget multiplier, `mu >= 0` its cash multiplier, and `eta_y >= 0`
its housing-cap multiplier. The goods and housing conditions give

\[
\lambda+\mu=1/x,\qquad
\rho\equiv\alpha x/s=x(\lambda p+\mu L_d+\eta_y)\ge mp.
\]

The fertility condition and the identity `h = s + kappa n` imply

\[
\vartheta x=n(\chi+\kappa\rho),\qquad
ax=\rho h+\chi n,\qquad
ac=\rho h+E\chi n>mph.
\]

Combining this with `c + L_d h <= w` and `L_d >= mp` proves

\[
ph<\frac{aw}{Em},\qquad
qz=w+qv-c-ph\ge qv-(1-m)ph,\qquad
z\ge v-\delta w.
\]

Thus these bounds allow both young cap and finance regimes. Their strict
housing inequality uses the maintained positive goods cost and positive
fertility.

For old owners without the physical cap, the estate-slack and estate-binding
housing choices are respectively `gamma z/(Kp)` and
`(gamma + omega_B) z/(K P)`. The active branch is exactly their minimum.
Optimizing consumption and estate conditional on housing gives a strictly
concave housing objective, so adding its upper cap clips that unconstrained
optimum. This proves the stated formula with `Gamma`.

The old Euler identity is particularly useful. If `eta_o` is the housing-cap
multiplier and `m_o=1/c^o` the budget multiplier, then

\[
m_o z=K-\eta_o h^o\le K.
\]

The estate-floor term vanishes by complementarity. Therefore `c^o >= z/K`
continues to hold when either old constraint binds. Under (P), an uncapped old
house is strictly larger than its paired young house; a capped old house is
weakly larger because both retain the same tenure cap. If young finance were
slack, `x=q c^o/beta>w`, contradicting current cash feasibility. No inequality
between `beta` and `q` enters this argument.

The full fixed-fertility planner separately allocates adult goods and housing.
Its housing formulas imply `h_i^{y,F} >= h_i^{o,F}`, strictly for each uncapped
old household, because `alpha >= gamma` and `kappa n_i > 0`. If all old were
capped, all paired young would also be capped. Strict unused total capacity
excludes this, giving

\[
H_Y^F>\bar H/2\ge H_Y^{eq}.
\]

The individual statement also checks: if `t=alpha/lambda_H`, then the strict
aggregate gain implies `t > bar s`, since each planner adult-space allocation
is at most `t`. A reference young household below its cap with
`s_i <= bar s` therefore gains housing. A paired local housing transfer has a
strictly positive derivative for any reference young household below its cap,
since `alpha/s_i > gamma/h_i^o` under the theorem's total-house ordering.

## 2. Compatibility and existence: the missing intermediate derivation

The bound asserted in the construction follows from the young Euler identity:

\[
E=\lambda(c+ph)+\mu(c+L_dh)+\eta_yh
  =(\lambda+\mu)w+\lambda qv-\lambda qz+\eta_yh.
\]

Since `lambda qz = beta z/c^o <= beta K`, nonnegative old income gives

\[
\frac wx\le E+\beta K,\qquad x\ge\frac{w}{E+\beta K}.
\]

For bounded endowments and a positive lower bound on `w`, the goods and
finance multipliers are uniformly bounded. As price tends to zero, an
uncapped young housing condition would have a right-hand side tending to
zero, while `alpha/s >= alpha/H_O > 0`. Hence all conditional young houses
eventually reach their caps. Inequality (4), the adult-goods lower bound,
and `H_d >= H_R` then force fertility above replacement.

At high prices, `n < h/kappa` and the proved current-cash housing bound force
fertility uniformly below replacement. Conditional choices and values are
continuous because the conditional optimization is strictly concave; logistic
tenure probabilities are continuous too. Thus average fertility crosses
replacement. Housing clearing sets positive cohort mass. If every reference
young household were capped, (4) would force mean fertility strictly above
replacement, so the reference equilibrium has unused aggregate capacity.

The bound on the stationary price follows by averaging the same housing
inequality and using replacement. Increasing bounded old incomes can then
make every old household capped without violating (4).

For the owner-emulates-renter comparison, write the renter's next financial
wealth as `a_R'` and set `a_O'=a_R'-P h`. At `phi >= q`, the owner's finance
restriction becomes `q a_R' + (phi-q)P h >= 0`. Both young budgets and old
total resources coincide. At the renter optimum,

\[
e/c^o=\omega_B/q,\qquad P h^o/c^o\le\gamma/(1-q),
\]

so the stated estate condition makes that old allocation feasible for an
owner, including capped renters. Ownership value is consequently at least
rental value. Taste location zero implies owner probability at least one
half, and the housing-stock contradiction used to force planner young
renters to their cap follows exactly as written.

## 3. Exact capped counterexample

All displayed low-type and high-owner choices satisfy the original budgets,
fertility equations, slack caps, estate floors, and multipliers. In particular
the high owner's multiplier is `2/5 - 5/13 = 1/65`.

For the high renter, substituting its lifetime budget and both binding caps
gives `x_R=421/100-(2/5)n_R`. Its fertility equation reduces exactly to

\[
3600n_R^2-20900n_R+12209=0,
\]

whose admissible root is the reported expression. Its stated bounds imply
slack young finance and strictly positive young and old housing-cap
multipliers. The value-difference expression includes all young and discounted
old utility terms. Owner feasibility of the renter optimum and the unique
owner optimum's violations of the rental caps establish its strict positivity.

The specified logistic parameters yield exactly the asserted tenure shares.
The choices of `nu` and `bar H` clear replacement and housing with cohort mass
two. Both tenures and both endowment types have positive masses.

With both renter ages capped and both owner ages uncapped, housing clearing
gives exactly

\[
t\equiv1/\lambda_H=91/40-(53/40)\varepsilon.
\]

For `0 < epsilon < 1/2`, `t > H_R` and `t + 5/4 < H_O`, which verify the
claimed planner regime over the full interval. The total young-housing change
is exactly `(-5+7 epsilon)/40 < 0`. Every paired young housing marginal utility
exceeds its old counterpart. The exception in young finance is essential to
report accurately; it does not invalidate the intended counterexample to the
weaker positive-mass formulation.

## 4. Joint fertility with the consumption cushion

The fertility map is the smaller root of

\[
\chi n^2-[\chi t+(\alpha+\vartheta)X]n+\vartheta Xt=0,
\qquad t=h/\kappa.
\]

The stated formula is correct. Its square-root term is the Euclidean norm of
the nonsingular linear map
`(t,X) -> (chi t + (alpha-vartheta)X, 2 sqrt(alpha vartheta)X)`.
It follows that the map is jointly concave and homogeneous of degree one.
Implicit differentiation of the fertility condition proves that both first
derivatives are strictly positive. At fixed positive `X`, its second housing
derivative is strictly negative.

The joint planner equalizes adult goods at `X^J`. The goods identity is exact:

\[
2X^J+\chi m^J=\bar x+\bar c^o+\chi m.
\]

Under a hypothetical `m^J <= m`, condition (C) gives `X^J >= bar x_d` for
each retained tenure. Conditional Jensen is therefore applicable exactly as
written. It is crucial to use tenure-conditional goods means here; an overall
mean comparison alone does not establish this step.

At the joint optimum all uncapped young houses have common size `h_*`;
capped ones have size `H_d`. This allocation also maximizes the integral of
the concave fertility map at its own total young housing: uncapped housing
marginals coincide, and capped marginals are weakly larger. The all-capped
case is trivially the only feasible allocation at that young housing total.

The joint old-young housing comparison and unused capacity give strictly
more total young housing than the reference. From reference tenure-mean
housing, one can add housing componentwise to attain that larger total while
respecting the retained caps. Some positive-mass group receives a strictly
positive addition. Strict monotonicity of the fertility map gives the strict
inequality needed for the contradiction. Thus the proof concerns the actual
joint optimum and establishes `m^J > m`.

Finally, the conservative primitive certificate gives
`bar c^o >= (bar v - delta_0 bar w)/K >= 2 w_max`, while every reference
`x_i < w_i <= w_max`. This implies (C), with slack. The certificate can be
satisfied by raising bounded old incomes in the construction because these
means average over all retained tenures and hence preserve the underlying
endowment distribution.

The result leaves private fertility rechoice, the welfare of additional
future people, and a funded market-policy transition outside its conclusion.
