# Independent review: local OLG transition

September 5, 2026. Read-only reviewer; original household specification retained.
Scope: revised proof sections 5–7, including the mixed-model infinite-path argument and uniform limiting tail. The point-mass entrant qualification below has been adopted. Section 8 was outside this review.

## Judgment

Both newly supplied steps are valid under their stated local-branch conditions.

Section 5 now closes the prior initial-old and changed-order objections. Section 6 supplies the previously missing uniform tail argument, so the all-owner limiting economy supports an all-date finite-reform fertility sign. Neither result establishes all-date fertility signs in nearby mixed-tenure economies.

## Checks performed

- Read sections 5–7 against the dated original household, clearing, demographic, tenure, and inherited-state equations in [simplified_olg_amendment_proposal.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:55).
- Checked the earlier scoped review and the verification script/JSON as arithmetic support only.
- Derived the weighted convolution bounds and the nonlinear tail comparison independently.

### 1. Initial old and stationary endpoint

Valid.

For the constructed point-mass entrant specialization, the inherited old cohort is completely represented by
\[
(Y^{\rm pre},O^{\rm pre},\pi^{O,\rm pre},
a^{O,\rm pre},H^{O,\rm pre},a^{R,\rm pre}).
\]
At date zero, these determine old demand through owner/renter masses, financial assets, and the owners’ inherited title. This is exactly what the original equilibrium requires after an unexpected reform: \(G_0=K(\eta_{-1}^{\rm pre})\), not a reoptimized post-reform old cohort.

The normalized atomic distribution is indeed not TV-differentiable when its atom locations move. That is irrelevant here: the boundary operator uses the finite moments/states above, which are smooth because the conditional policies are smooth on the stated fixed-active-set neighborhood.

The renter checks are sufficient:

\[
(x^R,h^R,n^R,a_R')=(53/110,1/4,1/8,34/55),
\]
with positive saving and all log arguments; both caps bind with strict multipliers,
\[
MV_Y^R=848/825>1/5,\qquad MV_O^R=159/550>1/5.
\]
Thus there is no hidden renter kink at the \(\epsilon=0\) limit.

The stationary system is correctly reduced to replacement fertility and the rebate equation. At \((\epsilon,\chi,\tau^p,P,T)=(0,0,0,1,0)\), its \((P,T)\) Jacobian is
\[
\operatorname{diag}(-1,1):
\]
at \(\chi=0\), \(\nu n=d/P\), while the tax term vanishes to first order in \(P\) when \(\tau^p=0\).

One scope qualification should be explicit: the finite vector proves the theorem for the proof’s homogeneous/point-mass entrant specialization, not automatically for arbitrary heterogeneous \(F\) in the broader model description. The latter would require the corresponding smooth aggregate integrals.

### 2. Weighted infinite-path inverse and mixed-tenure lead

Valid.

With \(\eta=.4\),
\[
\lambda_1\simeq .0596<\lambda_2\simeq .2616<\eta<\Lambda\simeq2.5323.
\]
For arbitrary \(g\in X_\eta\), the stable convolution is bounded because
\[
\sum_{j<t}\lambda_2^{t-1-j}\eta^{j+1}
\leq C\eta^t/(\eta-\lambda_2),
\]
and the unique \(X_\eta\)-bounded unstable solution is bounded because
\[
\sum_{j\ge t}\Lambda^{t-1-j}\eta^{j+1}
\leq C\eta^t/(\Lambda-\eta).
\]
The fixed initial young-population condition plus the date-zero inherited-old clearing condition determine the two stable coordinates. The determinant in (T5) is nonzero (\(\simeq .58887\)), so this is a bounded isomorphism including arbitrary dated and boundary residuals.

Linearized demography recovers prices through a finite shift/difference (at the limit, the log-linear relation is of the form \(p_t=y_t-y_{t+1}\), plus the reform term), hence boundedly on \(X_\eta\).

The original renter value does contain the genuine lead
\[
qr_{t+1}=(1+q\tau^p)P_{t+1}-qP_{t+2}.
\]
Section 5 correctly does not retain the old scalar dynamic-map claim. In \(X_\eta\), the two-date lead is a bounded shift; it enters aggregate clearing multiplied by the renter mass, \(O(\epsilon)\). Uniform conditional \(C^1\) policy/value maps and the strictly preserved active set make this a small bounded perturbation. Thus the Neumann-series/IFT conclusion is justified.

### 3. Uniform finite-reform fertility sign

Valid in the all-owner limit only.

The tail algebra in section 6 is correct. For
\[
v_{t+1}=\operatorname{diag}(\lambda_1,\lambda_2)v_t+R(v_t),
\quad |v_t|\le C\delta\eta^t,
\]
the first coordinate satisfies
\[
v_{1t}=v_{10}\lambda_1^t+O(\delta^2\eta^{2t}).
\]
Defining
\[
B_2=v_{20}+\sum_{j\ge0}\lambda_2^{-1-j}R_2(v_j)
\]
is legitimate because \(\eta^2<\lambda_2\), and gives
\[
v_{2t}=B_2\lambda_2^t+O(\delta^2\eta^{2t}).
\]
Transversality yields the limiting coefficients \(c_1,c_2\). Continuity preserves the positive two-mode increment bound, at least \((m/2)\delta\lambda_2^t\). Since
\[
\frac{\delta^2\eta^{2t}}{\delta\lambda_2^t}
\leq \delta\sup_t(\eta^2/\lambda_2)^t\to0,
\]
the remainder is uniformly dominated over every date. Hence \(Y_{t+1}>Y_t\), and the demographic identity gives \(n_t>1/2\) at every date for sufficiently small \(\delta>0\).

### 4. Welfare derivative

Correct:
\[
x^*=95/232,\qquad
-\frac{1-q}{x^*}-\frac{\beta\gamma}{d}
=-\frac{232}{475}-\frac{57}{475}
=-\frac{289}{475}.
\]
It is properly separated from the population result and does not support a Pareto claim.

## Minimal corrected theorem statement

For the homogeneous-entrant mixed-tenure specialization, on the stated strict conditional-policy branch, there is a neighborhood of \((\epsilon,\chi,\tau^p)=(0,0,0)\) with \(\epsilon,\chi,\tau^p>0\) and of \(d=1\) in which a sufficiently small permanent increase in \(d\) has a unique local exponentially convergent equilibrium path. It raises initial fertility and terminal household population. In the all-owner limiting economy, the same reform raises fertility at every date and population monotonically along the transition.