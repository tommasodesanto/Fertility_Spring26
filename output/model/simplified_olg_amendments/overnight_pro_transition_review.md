# Adversarial review of the Pro transition and welfare arguments

September 9, 2026 UTC. Reviewed the transition and living-household welfare
appendices of `oracle_consolidated_theory_source.tex`, together with its compact
Proposition T. Cross-checked against the frozen household budgets and the
earlier independently verified sequence-space system. No browser, model run,
numerical root calculation, or build was used. Only this file was written.

## Verdict

**Pass, with two proof/scope clarifications below.** I found no fatal algebraic
or economic-accounting error in the proposed transition region. Its analytic
root argument can replace the earlier conservative small-ownership inverse
bound for this more restricted region. It permits a positive owner share
without requiring that share itself to be small, while imposing explicit
restrictions on child-cost composition, relative old housing, and the old
housing expenditure share.

| Claim | Result |
|---|---|
| Exact owner/renter transition equations and endogenous tenure | Pass |
| Cubic root bound: two stable roots and one simple unstable root | Pass |
| Two additional zero roots and projection onto inherited states | Pass |
| Smooth local nonlinear forward map despite zero roots | Pass |
| Geometric local transition after permanent shocks | Pass |
| Uniqueness among all nearby bounded paths | Valid, but add the unweighted contraction sentence below |
| Tax impact fertility, price overshooting, and initial price recovery | Pass |
| Fertility-preference impact and endpoint signs | Pass |
| Nonempty analytic primitive family with finite caps | Pass |
| General local living-household welfare envelopes | Pass |
| Closed-form welfare threshold | Pass at the stationary intervention reference |

Two clarifications should be incorporated:

1. A contraction on a geometrically decaying sequence space alone establishes
   uniqueness within that space. Add that the same stable forward and unstable
   backward kernels are absolutely summable on unweighted
   \(\ell^\infty\). Making the nonlinear derivative sufficiently small also
   makes the Lyapunov–Perron map a contraction there. Every nearby bounded
   trajectory satisfies that system, so it coincides with the geometric
   solution. This completes the stated bounded-path uniqueness under the
   existing assumptions.
2. The equality \(N^{-1}d\mathcal W/d\tau=M A_W-B_W\) uses a stationary
   reference: equal current cohort masses, reference old title holdings, and
   a date-independent first-order rebate \(T'_0\). At a nearby nonstationary
   policy reset, use the general envelopes with actual masses, claims, and
   continuation derivatives. A **strict** welfare sign at the stationary
   reference extends to sufficiently nearby reset states by continuity; the
   same closed-form equality is not exact at those nonstationary states.

## 1. Exact map and equivalence to the earlier budget-based linear system

The binding owner condition with \(\phi=q\) implies
\(a'=-P_t h\). The future resource equation is therefore exactly
\[
z_{i,t+1}=v_i+T_{t+1}+(P_{t+1}-P_t)h_i^O.
\]
Its current budget uses \(L_t\), whereas a constrained renter with zero
saving uses \(u_t=L_t-q\Delta P_t\). The owner first-order conditions
are those obtained by differentiating these original budgets. The next old
indirect-utility term \(-\beta\gamma\log u_{t+1}\) is common to both
tenures and cancels from their value difference. Consequently its omission
from the current choice map does not omit an additional forward price variable.
That term remains in the welfare envelope.

The old aggregate equation is exact in the maintained regime because old
housing is linear in total resources. At the initial date,
\[
Z_0=\bar a_0+P_0\bar H_0+V+T_0.
\]
For a reference cohort generated under binding finance,
\(\bar a_0=-P_{-1}B_{-1}\) and \(\bar H_0=B_{-1}\). Holding those
claims and titles fixed retains the initial capital gain or loss. The
normalization \(z_{-1}=0\) fixes the inherited debt-issue price; it does not
revalue debt at the new price or impose anticipation before the surprise.

To map the Pro notation to the earlier checked system, let
\(d=1-q\), let \(F_t=\delta P_t/p\) denote the earlier normalized
price, and let \(B_f\) be the earlier heterogeneous owner-response moment.
Then
\[
z_t=dF_t,\qquad \epsilon=\eta,\qquad
R=\pi B_f/q,\qquad
\bar h^y=\frac{WA}{Ep},\qquad
\frac{B}{A}=\frac{a+\eta c}{a+c}
\]
in the earlier quantity notation. Substitution maps the two derivative
systems exactly. In particular,
\(\mathsf f=(q\eta-\pi d_n)/d\), and the housing coefficient and
initial-old title coefficient likewise coincide. There is no discrepancy in
the owner cash term, heterogeneous moment, or lagged capital revaluation.

The endogenous tenure derivatives cancel from first-order real aggregates
because the conditional reference quantities coincide. The derivative of
\(B_{t-1}\) is multiplied by the zero reference price change in old resources.
Neither cancellation freezes tenure or title choices in the nonlinear system.

## 2. Independent root verification

The stated restrictions imply
\[
E\epsilon-A=(1+\alpha)\epsilon-\alpha\ge0,
\qquad
EB-A^2=\alpha+\alpha\vartheta(1-\epsilon)^2+
\vartheta\epsilon^2>0.
\]
Together with \(0<R<\pi<1\), these give
\[
0\le\mathsf f\le q\epsilon/d,\qquad
\mathsf a\ge q\zeta/d>0.
\]
The restriction \(qK\zeta>\pi\gamma\) gives
\(\mathsf a>\mathsf c\).

The signs of \(\mathcal P(1)\) and \(\mathcal P(-1)\) are correct.
The positive cubic leading coefficient therefore gives a real root
\(r_u>1\). After factoring that root, the quadratic factor is positive at
both \(1\) and \(-1\). Its constant coefficient is exactly
\[
c_q=\frac{\zeta(\epsilon+\mathsf f)-\mathsf c}
          {\mathsf a r_u}<1,
\]
because
\(\zeta(\epsilon+\mathsf f)\le\zeta\epsilon/d<q\zeta/d\le\mathsf a\).
The two positive endpoint values imply
\(c_q>-1\) and \(|b_q|<1+c_q\). These are the strict quadratic Jury
inequalities, so both remaining roots are inside the unit disk. This also
proves that the unstable root is unique and simple. No numerical root or
small-owner-share continuity argument enters the proof.

## 3. State dimension, projection, and nonlinear existence

Solving the homogeneous housing row for the next price yields a four-state
population-price matrix. In the coordinate order used by the source, its rows
satisfy
\[
\mathrm{row}_1=\mathrm{row}_2+
\mathsf f\,\mathrm{row}_4-(\epsilon+\mathsf f)\mathrm{row}_3.
\]
This supplies one zero root. The separate title state has a zero column in
the first-order feedback block and supplies another zero root. The total map
therefore has four stable dimensions and one unstable dimension, including
any possible multiplicity of zero among the cubic's stable roots.

For a left unstable eigenvector whose current-price coefficient were zero,
the first two coordinate equations indeed give
\[
r_u-1=-\frac{\mathsf f}{\mathsf a}
       \left(1+\frac\zeta{r_u}\right),
\]
which is impossible. The coefficient on the independent title state is zero.
Thus the stable tangent space projects invertibly onto the four inherited
coordinates. A zero stable root does not obstruct that projection and does
not require inversion of the full forward map.

The finite-dimensional nonlinear map exists locally. At zero tax, the
implicit next-rebate/demographic block has identity derivative in the next
population coordinate. After eliminating it, the housing derivative with
respect to the next price is proportional to \(\mathsf a>0\). Both
properties persist locally. Conditional demands and their endogenous tenure
weights are smooth under the compact-support and uniform regime assumptions.

The stated stable forward and unstable backward sums establish a local stable
manifold for this possibly noninvertible map. They give geometric convergence
in a weighted norm. Adding the unweighted-norm argument stated above establishes
uniqueness among all nearby bounded paths. This is a finite-dimensional
counterpart of the earlier sequence-space inverse and does not conflict with
it.

The reduced inherited state is sufficient for **aggregate clearing** in the
uniformly slack old regime. Actual individual inherited claims and titles must
still satisfy those regime margins, and they are needed for welfare. The
common-state experiments in the source generate such states from a nearby
baseline path; an arbitrary aggregate title total alone is not a certificate
of individual feasibility.

## 4. Impact formulas and signs

Normalizing the current-price component of the left unstable eigenvector to
one gives its other three components as
\[
-\frac{B_u}{D_u},\qquad
-\frac{\zeta(r_u-1)}{r_uD_u},\qquad
\frac{\mathsf c(r_u-1)}{r_uD_u}.
\]
Projecting the initial deviation from the new stationary level onto its stable
hyperplane gives exactly both forms of the source's impact-price formula.
The date-zero housing row gives its price-recovery formula. Direct substitution
into initial demographic growth, using the reported root identity, gives the
stated \(L_N\) and \(L_P\). Both signs are correct.

For the tax perturbation, \(V\ge W\) makes \(\ell_*>0\), and the upper
bound on \(\zeta\) makes \(z_*<0\). Therefore
\[
\ell_1>0,\qquad z_0<z_*<0,\qquad z_1-z_0>0.
\]
The proof does not imply monotonicity of the subsequent price path or a
positive impact response of young mean housing. The separate housing-impact
inequality in the source is needed for the latter conclusion.

For a fertility-preference increase, both direct forcing coefficients are
positive. Independent simplification reproduces
\[
\mathsf b-\frac{\epsilon C_{\vartheta}}{g_{\vartheta}}
=\zeta+\frac{\alpha(1+\alpha+\vartheta\epsilon)}
                 {(1+\alpha)A}>1.
\]
The reported ratio of \(L_P\) to \(L_N/(1+\zeta)\) is correct and is
strictly below one under \(\pi\gamma/K<d^2\). This proves the initial
fertility sign. The endpoint derivatives and the stationary tax formulas also
check directly.

Strict signs and the stable-manifold projection extend to sufficiently nearby
inherited states. At a policy reset one must preserve actual claims, titles,
and cohort masses, rather than preserving old resources evaluated at each
scenario's new house price. The baseline preference and tax sequences are
permanent after their respective shocks, so geometric tail convergence and
absolute summability of the tail log fertility differences follow. This is
stronger than the mere convergence conclusion for arbitrary convergent forcing
in the earlier sequence-space result, and uses the permanent-shock restriction.

## 5. Welfare envelopes and the proportional-income threshold

The old envelope includes the inherited-title gain
\(m_tH\,dP_t\) and the estate-floor term
\(-\rho_t h_t^o\,dP_{t+1}\), with the correct signs. The young envelope
follows from the reduced lifetime budget and cash restriction. In particular,
the future title effect is already contained in the current user-cost term;
adding it separately would double count it. The future old user-cost and
estate-floor terms must remain, as they do in the source.

Integrating maximized young utility including its ownership taste introduces
no extra switching term: the two utilities agree at the switching boundary.
This does not imply cancellation of tenure changes from fertility away from
the symmetric reference.

For \(v_i=a_vw_i\), write \(M=W\mathbb E(1/w_i)\). At a stationary
intervention reference, the integrated old contribution is
\[
M\frac{T'_0}{W}\frac K{a_v}
-\gamma U_0+\frac{\pi\gamma}{\zeta d}z_0.
\]
The integrated young contribution is
\[
M\frac{T'_0}{W}\left(E+\frac{\beta K}{a_v}\right)
-A U_0-\frac qd\pi A(1-r)\Delta z_0-\beta\gamma U_1.
\]
Their sum is exactly \(M A_W-B_W\). In particular, the negative initial
price change creates an old title loss in this formula; it has not been
discarded or compensated by assumption.

The heterogeneity construction is valid. Holding \(W\) fixed and taking
two equally weighted positive endowments with one approaching zero makes
\(M\) arbitrarily large. Proportional old income keeps all aggregate
linear coefficients and price derivatives fixed. A single pair of finite
caps can cover the entire family by placing them above the limiting upper
endowment demands. Strictly speaking, the individual maximum-demand bound
changes as the upper endowment approaches \(2W\); the same cap inequalities
can remain satisfied, rather than their numerical bounds being unchanged.

For every selected economy the lower endowment is positive, so its local
welfare derivative and implicit-function neighborhood are well defined. No
uniform positive tax magnitude is established as that lower endowment tends
to zero. Select a finite distribution satisfying the strict threshold first,
then choose a sufficiently small tax change. The criterion concerns the
specified sum of living households' remaining utility; it is not a Pareto
claim, a full future-population welfare criterion, or implementation of the
dated planner.

## 7. Primitive region and remaining scope

The proposed interval for \(\zeta\) is nonempty exactly because
\(E\epsilon>A\). Increasing \(K\) can satisfy the estate and capital-
feedback restrictions while proportional old income maintains the selected
\(\zeta\). The financing ratio is then
\(r_i=\beta\gamma/(qA\zeta)<1\). In the displayed nonempty family,
\(\alpha\ge\gamma\) implies \(A\zeta>\gamma\), so the allowed
\(\beta\) interval contains values on both sides of \(q\).

The theorem retains \(\phi=q\), inactive competitive physical caps, and a
slack old financial-estate floor. It does not extend the transition proof to
the separate active-rental-cap dated-planner construction. Its intermediary
capital-loss treatment remains the maintained external-financier convention;
a different domestic landlord balance sheet would need its own accounting
and welfare analysis. Policy and initial-state neighborhoods remain local and
are not given explicit numerical sizes.
