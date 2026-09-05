# Independent review: local OLG transition

September 5, 2026. Read-only reviewer; original household specification retained.
Scope: recurrence, household conditions, local path, and initial proof gaps. The later closure review checks their resolution.

The all-owner local argument is mathematically sound, conditional on a local branch interpretation. The mixed-tenure extension is plausible but not yet a theorem: it needs an infinite-dimensional implicit-function setup with an explicitly verified bounded inverse and smoothness of the initial-old boundary map.

1. All-owner derivation

From the owner budget in [simplified_olg_amendment_proposal.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:110), with \(\chi=\tau^p=0\), the binding down-payment constraint gives
\[
h_t=\frac{d}{P_t},\qquad d=\frac{b}{1-\phi}.
\]
The fertility FOC [eq. `amend_fertility_foc`] gives
\[
n_t=\frac{\vartheta}{\kappa(\alpha+\vartheta)}h_t=\frac12 h_t.
\]
Since \(\nu=2\),
\[
Y_{t+1}=Y_t h_t,\qquad P_t=d\frac{Y_t}{Y_{t+1}}.
\]

The young owner’s intertemporal FOC, correctly retaining the anticipated resale/rental value, is
\[
\rho c_t=w-d+qd\frac{P_{t+1}}{P_t},\qquad
\rho=1+\beta A,\quad A=1+\gamma+\omega_B.
\]
For an old household at \(t\), substituting the preceding young choice yields
\[
qR_t=\frac{\beta A}{\rho}
 \left[w-d+qd\frac{Y_t^2}{Y_{t-1}Y_{t+1}}\right].
\]
Its interior housing choice is
\[
h_t^2=\frac{\gamma R_t}{q r_t A},\qquad
qr_t=P_t-qP_{t+1}.
\]
Hence housing clearing is exactly
\[
H=Y_{t+1}
+D\frac{\frac{r-1}{q}Y_{t-1}+Y_t^2/Y_{t+1}}
 {Y_t/Y_{t+1}-qY_{t+1}/Y_{t+2}},
\quad D=\frac{\beta\gamma}{\rho},\quad r=\frac wd.
\]
So the stated recurrence is correct. The \(P_{t+1}/P_t\) term in the young saving FOC is essential; dropping it would give the wrong coefficient \(D\).

The inherited old cohort has title \(1\) and \(a_0=-301/928\). Its date-0 resource is \(a_0+d/Y_1\), while \(qr_0=d/Y_1-qdY_1/Y_2\). Thus the stated boundary equation
\[
F_0=Y_1+\frac1{11}
 \frac{-301/928+d/Y_1}{d/Y_1-qdY_1/Y_2}-H=0
\]
is also correct. Direct differentiation gives exactly
\[
(F_{Y_1},F_{Y_2},F_d,F_H)
=\left(\frac{33783}{10208},-\frac{285}{232},
\frac{1505}{10208},-1\right).
\]

The stationary expressions check:
\[
h^{2*}=\frac{75}{232}\left[\frac wd-(1-q)\right],\qquad
Y^*=\frac{H}{1+h^{2*}},
\]
and at \(d=1\),
\[
\frac{dY^*}{dd}=\frac{345}{1213}>0.
\]
The reported baseline constraints are strict where required:
\[
a'=\frac{627}{1160}>0,\quad h^2=\frac{285}{928}<1,\quad
e=\frac{475}{928}>P h^2,\quad MV^Y=\frac{19}{87}>\frac15.
\]
Old financial wealth may indeed be negative: [eq. `amend_young_owner`] constrains \(a'\ge0\), not the post-mortgage old \(a\).

Linearizing gives coefficients
\[
\left(\frac{45}{928},-\frac{945}{928},
\frac{3253}{928},-\frac{285}{232}\right),
\]
hence
\[
1140z^3-3253z^2+945z-45=0.
\]
The supplied rational brackets are valid: the polynomial has signs
\[
(-,+)\text{ on }(.059,.060),\quad (+,-)\text{ on }(.261,.262),
\quad(-,+)\text{ on }(2.53,2.54).
\]
Because these are three disjoint intervals for a cubic, they establish exactly three positive real roots.

The boundary/stable-mode system is transverse. Its solution has
\[
c_1\simeq0.87792408,\qquad c_2\simeq-1.16234288,
\]
and
\[
y_1-y_0\simeq .03265445,\qquad y_2-y_1\simeq .17533322.
\]
Thereafter
\[
y_{t+1}-y_t=
c_1(\ell_1-1)\ell_1^t+c_2(\ell_2-1)\ell_2^t>0:
\]
the negative term decays at \(\ell_1^t\), the positive term at \(\ell_2^t\), and positivity at \(t=0\) suffices. Thus first-order fertility is positive at every date. The stock derivative has both transient coefficients negative, so it is likewise monotone.

Initial-boundary nonsingularity, \(F_{Y_{t+2}}\ne0\), and stable-boundary transversality do prove a unique **local nonlinear path converging to the new steady state among paths remaining in the chosen neighborhood**. They do not prove global uniqueness of all convergent nonlinear paths.

2. Mixed tenure and changed dynamic order

The concern is real. Under [eq. `amend_young_renter`] and [eq. `amend_old_renter`], a renter’s young value contains the old rental cost
\[
qr_{t+1}=(1+q\tau^p)P_{t+1}-qP_{t+2}.
\]
Therefore \(W_t^R\), and thus
\[
\pi_t^O=\frac1{1+\epsilon
 \exp\{(W_t^R-W_t^O)/\sigma_\xi\}},
\qquad \epsilon=e^{-\bar\xi/\sigma_\xi},
\]
depends on \(P_{t+2}\). This adds a genuine lead to clearing; it cannot be treated as the same scalar recurrence merely with perturbed coefficients.

A sequence-space IFT can handle this. Extend the equations smoothly to
\[
\mu=(\epsilon,\chi,\tau^p)=(0,0,0),
\]
fix \(\sigma_\xi>0\), and work in an exponentially weighted space
\[
\|v\|_\eta=\sup_{t\ge0}|v_t|/\eta^t,\qquad \ell_2<\eta<1.
\]
The all-owner recurrence plus the initial boundary has a bounded inverse on this space: the unstable root is excluded by boundedness, while the two stable modes and the transverse boundary determine the solution. Finite forward/backward shifts are bounded on this space. If conditional renter policies and values are uniformly \(C^1\) on a neighborhood, their added lead terms are \(O(\epsilon)\) bounded shift operators; a Neumann-series/IFT argument then works for sufficiently small \(\mu\).

But this requires one additional substantive condition beyond “conditional renter choices are well defined”:

\[
\boxed{\text{The full initial-old map }G_0(\mu)=K(\eta_{-1}^{\rm pre}(\mu))
\text{ from [eq. `amend_finite_memory`] must be }C^1
\text{ into the boundary operator, uniformly on the weighted neighborhood.}}
\]

For \(\epsilon>0\), the pre-reform old cohort includes an \(O(\epsilon)\) renter mass; retaining the exactly all-owner inherited cohort is not the stationary original-model boundary. Strict renter caps help establish this smoothness, but one must also verify interiority/strictness of every remaining renter and owner constraint and uniform boundedness of the relevant value derivatives. Without that condition, ordinary finite-dimensional continuity does not justify the extension.

3. What may be claimed

Retain the following short theorem:

> In the all-owner \((\chi,\tau^p)=(0,0)\) limiting economy, there is a unique local convergent equilibrium path on the specified smooth branch after a sufficiently small permanent increase in \(d\) or \(H\). A rise in \(d\) raises terminal \(Y^*\), and the derivative of fertility is strictly positive at every finite date.

Do not yet claim positive fertility at every date for a finite reform. Pointwise continuity gives any fixed finite prefix only. Since the first-order effect tends to zero as \(t\to\infty\), a uniform \(O(\delta^2)\) remainder is insufficient. One needs either an invariant cone for the nonlinear transition or a tail asymptotic result showing that the \(\ell_2^t\) stable coefficient keeps its strict sign under the reform, combined with a finite-prefix continuity argument.