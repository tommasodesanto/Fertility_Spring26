# Analytical efficiency result: assessment

7 September 2026. Pro conversation: https://chatgpt.com/c/6a9f2beb-c8b0-83e9-b649-3304f9801ad3

The response gives a valid analytical route to the requested equilibrium result, in the specialization with identical entrant income and wealth and zero property tax. The household and equilibrium arguments received separate max-reasoning reviews; the lead checked the welfare comparison and settlement. The manuscript has not been changed. This is a research assessment, not proposed paper prose.

## What the result establishes

The model does not imply inefficiency merely because financial and rental-size restrictions are present. Pro constructs a class of stationary equilibria where they are slack and the marginal housing values of young and old coincide. Concavity rules out any Pareto improvement in the specified comparison: consumption and housing can be reallocated at one date, with fertility, tenure, estates and future real allocations fixed. This is not a claim of efficiency under arbitrary intergenerational reallocations or endogenous population comparisons.

Pro then provides explicit inequalities in primitives under which a positive stationary equilibrium exists and **every positive stationary equilibrium** has a housing reallocation that benefits young owners and compensates old owners. The inequalities establish the household constraint pattern; the theorem does not assume that pattern or a marginal-value gap. Its primitive region is shown analytically to be nonempty.

The restrictions are sufficient, not necessary. The result does not cover arbitrary heterogeneity, positive property-tax rebates, arbitrary transitions or a planner restricted to transfers followed by market clearing. It keeps positive child goods and space costs, finite physical housing caps, endogenous household fertility and mixed tenure with finite logistic taste dispersion.

## The simplest useful version

For exposition, retain the branch where old owners leave positive financial assets, and suppose the financed share satisfies $\phi\ge q$. Write

\[
M=1+\beta(1+\gamma+\omega_B),\qquad
\lambda=\frac{(1-q)b}{(1-\phi)(y+b)}.
\]

Here $\lambda$ is the user cost of the largest home the down payment permits, relative to the young household's income plus initial wealth. The central condition simplifies to

\[
\frac{\beta\gamma}{qM+\beta\gamma}
<\lambda<
\frac{\alpha}{M+\alpha+\vartheta},
\qquad (1-q)\omega_B>q\gamma.
\]

The upper bound makes young owners want more housing. The lower bound ensures that old owners sell part of the house. The estate restriction makes old owners' marginal housing value equal user cost. Positive gross saving follows from $\phi\ge q$ and the household saving condition.

These inequalities are **not the complete stationary-equilibrium theorem**. Condition (C) in the response also restricts replacement fertility $1/\nu$ between explicit endpoint fertilities. Its price endpoints are functions of primitives, not unknown equilibrium prices. It ensures that every stationary equilibrium has a binding young rental cap and a slack owner physical cap, and lies in the region where the owner result applies. Pro proves that this demographic interval is nonempty. The larger theorem also permits old owners with zero financial estates, using a stricter upper liquidity bound.

The stationary cohort mass is endogenous: after the fertility equation determines a stationary price, housing clearing determines population. No convergence or reachability from a prescribed initial state is proved.

## Verification and exposition repairs

1. **Owner choice.** In the relaxed problem let $z=a'+Ph$ be total wealth entering old age and $A=1+\gamma+\omega_B$. Joint homogeneity gives $zV_z+hV_h=A$. The saving condition and $V_h\ge0$ imply $qz\le\beta Ax$, including binding retention. This validates the lower bound on young consumption and the proof that owners choose $\min\{h_O^{\max},b/[(1-\phi)P]\}$. Pro restores the omitted nonnegative gross-bond constraint using condition (B); it does not confuse gross bonds with net $a'$.

2. **Old estate and fertility formulas.** With retention slack, the old-owner formulas and marginal housing value $vP$ are correct, where $v=\max\{1-q,\gamma/(\gamma+\omega_B)\}$. The explicit fertility root solves the stated quadratic and increases strictly in both available resources and housing. The renter's old-age cap is slack inside the bracket, so its reduced coefficient is $M$.

3. **Renter demand wording.** Equation (15) is the demand when both rental caps are removed. It is not globally the solution obtained by removing only the young cap, since at low prices the old cap can bind. Under the theorem's restrictions it is valid for the comparison inside the stated price bracket. A rewrite should state this domain explicitly. The household proof does not fail.

4. **All stationary equilibria.** Condition (A) implies $\alpha>\beta\gamma/q$ and $t_-<\ell<\alpha w/D<t_+<w$. At $P_0=B/\rho$, inside the bracket, the current owner and renter allocations coincide, even in the old estate-binding branch. This proves the demographic interval is nonempty. Each tenure's fertility falls globally with price; the logistic mixture need not. Endpoint bounds on each tenure separately exclude all stationary roots outside the bracket.

5. **Regime changes.** The response's short assertion about global regime changes can be expanded without adding assumptions. For fixed current housing $H$, fertility on each branch is $\mathcal N_m(w-aPH,H)$, so both fertility and $x/P$ decrease with $P$. Let $r=z/(PH)$. If $(1-q)\omega_B\ge q\gamma$, retention changes status at $r=(1-q)A/\gamma$. If $(1-q)\omega_B<q\gamma$, as $r$ falls the regimes are retention only, both old constraints, then estate only, with thresholds $1+q/\omega_B$ and $A/(\gamma+\omega_B)$. The policies join continuously. These observations justify the global monotonicity used in the proof.

6. **Welfare.** The exact old-age consumption compensation has derivative $MV^o=vP$, strictly below the young owner's housing value. A sufficiently small reallocation raises young utility and preserves old utility. The financial settlement preserves each estate, future real allocations and existing repayments. It works also when the old donor initially has no financial saving: the transaction replaces estate housing with bonds. Physical housing caps remain respected.

No numerical equilibrium, parameter sweep or computer-assisted neighborhood was used to establish these claims. All 122 rendered TeX labels were preserved in the readable copy of Pro's response. The response and original model packet are retained alongside this assessment.

## Next decision

Discuss whether the homogeneous, zero-property-tax specialization is an acceptable main illustrative theorem. If so, present the simpler positive-financial-estate case first, with the demographic restriction stated honestly and derived in the proof. Keep the old estate-binding extension separate. Preserve the author's notation when drafting. The existing transition work remains separate and has not been replaced.
