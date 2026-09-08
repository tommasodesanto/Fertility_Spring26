# Assembled utilitarian note: final transfer review

September 8, 2026. Read-only review of the complete 811-line
`latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex`, with a subsequent
targeted read confirming the corrected opening of “Housing allocation.”
No builds, simulations, or source edits were performed.

**Final verdict: PASS within the stated scope. No outstanding required
correction was found in the transfer theorem, primitive coverage, financial
accounting, welfare normalization, or separation from free fertility.**

## Required correction found and resolved

The previous opening said that the relevant constrained young owners “would
use more housing if [the mortgage restriction] were relaxed.” The stated
primitive conditions did not ensure that comparative static. It would need
the earlier additional condition \(\beta KL\ge p-L\), or the corresponding
household-specific condition. Neither the new direct-utilitarian proof nor
the transfer implementation requires this claim.

The lead removed it. I verified that the current source asserts only the
strictly binding mortgage and slack young housing cap, old housing cap, and
old estate floor. Those properties are supplied by the stated primitive
appendix. No additional finance-relaxation restriction is needed.

## Transfer theorem and optimization — PASS

Proposition “A transfer improvement with committed fertility” inherits
\(\beta\ge q\) from the preceding direct proposition and adds
\(\phi\ge q\), \(\alpha(1+\omega_B)\le\gamma L/p\), and the capped
low-marginal-utility funder group. Its interpretation as an unexpected
intervention after fertility and tenure commitments is explicit.

The grants remain predetermined household-specific amounts rather than
subsidies contingent on subsequent housing. Setting
\(\Delta a'=-\phi P\varepsilon/q\) verifies the original young cash budget,
the original mortgage covenant, and
\(\Delta a'+P\varepsilon+J=0\). Thus inclusive old resources remain exactly
unchanged. The full positive-log objective is strictly concave on an affine
feasible set when fertility and tenure are fixed; the displayed first-order
conditions and persistent strict margins establish global conditional
optimality. The original old estate restriction and nonnegative gross
financial saving are retained.

The residual \(R\) and its derivative bound are correct. The welfare
derivative correctly includes the matching old donor's utility loss and the
funders' marginal utility of cash. Equal collection across a funder group
produces its mean marginal utility in the formula; integrating over selected
recipients handles differing group masses. The common finite step follows
from positive-mass subsets with uniform strict margins, not an assumption of
identical endowments.

## Primitive coverage and nonemptiness — PASS

The price and rebate bounds are evaluated on the full endowment distribution.
The eligible-income condition makes finance strictly restrictive despite the
common rebate. The cap and estate conditions supply the precise branches used
in the proof. Stationarity supplies matching current old owner measures, and
logistic tastes provide positive owner mass.

The funder bound correctly uses the covenant lower bound
\(z_f\ge y_f^o-(\phi/q-1)P_+h_O^{\max}\), the old cap threshold
\(Kp h_O^{\max}/\gamma\), and the recipient bound
\(m_i\ge qK/M_{\mathcal S}\). The maintained estate condition makes the
funder estate floor slack whenever its housing cap is strictly binding.

The nonemptiness argument is self-consistent: choose \(\phi>q\) near enough
to \(q\) that large \(\vartheta\) gives \(k<q\), choose sufficiently small
\(\chi>0\) for stationary existence, and give funders income \(V\) with
mass \(c_*/(V-V_S)\). Their contribution keeps
\(\bar M_0=W+qV_S+qc_*\) fixed as \(V\) increases. All bounds on the
right side of the funder inequality therefore remain fixed, while its left
side increases. A sufficiently large finite \(V\) works with positive mass
and bounded endowments. This is analytical compatibility, not a numerical
reference-equilibrium argument.

## Dated accounting and welfare weights — PASS

Net current program revenue is \(qJ\), whose bond return funds \(J\) exactly
next date. The larger inherited titles of current young offset the smaller
death sales of current old. The affected young choose their original old
allocations, and the inherited household state is restored from date two.
All new mortgage and government transactions settle finitely. Fixed occupied
stock and cohort masses preserve property-tax revenue and the original
common rebate.

The explicit identity \(\Delta C_0+q\Delta E_1=0\) includes the decrease in
initial old estates. There is no omitted external financing requirement or
terminal refinancing assumption. Estates do not enter new households' wealth
in this model, which is essential to the finite-tail result.

The welfare objective counts individual households and gives current young
weight one on lifetime utility, including their private \(\beta\), and
current old weight one on remaining utility. The note correctly distinguishes
this from weighting each cohort's utility at birth, which would put
\(\rho=\beta\) on initial old utility. The direct gap with general old
weight \(\rho\) is correct. No Pareto claim is attached to uncompensated old
losses.

## Free fertility and shortened prose — PASS

The fixed-fertility transfer is not presented as a simultaneous-choice
equilibrium. Restored fertility uses its own household condition, corrected
grants, and the additional funding term
\([\chi+\kappa(p-1/g)]\Delta n\). The note explicitly says the previous
conditions neither complete this funding/welfare argument nor keep the next
cohort mass unchanged. It therefore does not recycle the finite-tail theorem
as an endogenous-fertility transition result.

The knife-edge example remains valid conditional on committed fertility and
tenure; “sufficiently small” is justified by the surrounding strict-margin
argument. The dated fertility and cohort identities are explicitly
conditional on existing equilibrium paths, with no unsupported convergence
or population-welfare conclusion. No remaining shortening overstates the
verified transfer result.
