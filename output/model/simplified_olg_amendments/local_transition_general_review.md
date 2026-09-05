# Independent review: general convergence and fertility conditions

September 5, 2026. Scope: proof sections 10.1–10.3, the r>1 convergence result, initial-condition accounting, exact/sufficient fertility tests, and counterexample. The requested choice of decay weight and inherited-state/operator qualifications have been made explicit. The later broader condition in 10.4 was outside this scoped report; it is checked separately in `local_transition_sharp_stability_review.md` and by the lead’s second algebraic proof and symbolic/root receipts.

Classification: the all-owner convergence theorem and both fertility conditions are valid, conditional on the already stated strict-branch and sequence-space assumptions. The mixed-model extension is incomplete only in wording: it must explicitly retain Section 5/9’s \(C^1\) initial-old boundary and bounded-operator conditions, with a weight \(\eta\) chosen above the particular stable spectral radius.

What I checked:

- The general normalization is correct. With
  \[
  K=\frac{\nu\vartheta}{\kappa(\alpha+\vartheta)},\quad
  P_t=K d\,\frac{Y_t}{Y_{t+1}},
  \]
  multiplying clearing by \(K\) gives exactly (T1), with \(K\bar H\) on its left-hand side, and its stationary scale is
  \[
  Y^*=\frac{K\bar H}{1+C}.
  \]
  This follows directly from the dated owner and old-age budgets in [simplified_olg_amendment_proposal.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:110).

- (T11) is correct for a proportional rise in \(d=E[b_i]/(1-\phi)\), with both initial cohort masses normalized to one:
  \[
  F_1=1+\frac{C(1+q)-L}{\ell},\quad
  F_2=-\frac{Cq}{\ell},\quad
  F_d=\frac{L}{\ell}-C.
  \]
  The inherited old financial assets and titles are held fixed. In particular, the \(F_d\) term comes from changing current prices relative to those predetermined assets—not from reoptimizing them.

- The root argument for (T10) is valid for every \(q\in(0,1)\). Under \(r>1\), \(0<D<\ell C\), hence \(\mathcal P(0)<0\). Together with
  \[
  \mathcal P(1)=-\ell(1+C)<0,\qquad
  \mathcal P'(1)=-\ell(2+C)<0,
  \]
  there is one simple root \(\Lambda>1\). The displayed identity for \(\mathcal P(R)\) is correct and gives \(0<R/\Lambda<1\). The \(\mathcal P(-1)\) argument correctly rules out negative stable roots outside the unit circle:
  \[
  \mathcal P(-1)=4D-\ell-(3+q)C<0.
  \]
  This covers real, complex, and repeated stable roots. Thus precisely two roots lie strictly inside the unit disk.

- Boundary transversality is also valid. For the stable-root sum \(S\),
  \[
  F_1+F_2S
  =1+\frac{C(1+q-qS)-L}{\ell}
  >1-\frac L\ell>0.
  \]
  For a complex pair, \(S\) remains real and lies in \((-2,2)\); for a repeated stable root, the \(\lambda^t,t\lambda^{t-1}\) basis gives the same determinant. Nearby price recovery is legitimate because \(P_t=KdY_t/Y_{t+1}>0\) locally.

- (T12) is correct. The identity
  \[
  (1-\lambda_1)(1-\lambda_2)
  =\frac{\ell(1+C)}{Cq(\Lambda-1)}
  \]
  yields the stated formula for \(m\). Therefore:
  \[
  C\ge L/\ell \implies m>0,
  \]
  and, if \(C<L/\ell\),
  \[
  m>0\iff \Lambda<Z
  \iff \mathcal P(Z)>0,\qquad
  Z=1+\frac{Dr}{q(L-\ell C)}.
  \]
  This is genuinely primitive: it evaluates a known cubic at an explicit point and does not conceal a dynamic-root or equilibrium-price calculation.

- The \(q\ge1/2\) sufficient condition (T13) is valid. From \(q\ge\ell\), \(C>D/\ell\), and \(R/\Lambda< R\), one obtains \(S>0\), hence
  \[
  \Lambda<\frac{\ell-D+C(1+q)}{Cq}.
  \]
  Substitution produces exactly (T13). The quoted baseline gap checks:
  \[
  \ell C[\ell+(1+q)C]-L(\ell-D+C)
  =\frac{229247}{47365120}>0.
  \]

- The counterexample is valid. Its primitives give
  \[
  (C,D,L)=(.05,.02,.1),\quad J=\frac{2}{21}>0,
  \]
  with the claimed stationary owner allocation
  \[
  (x,h,n,a')=(.48,1,.5,.52),\quad
  (c^2,h^2,e)=(.1,.05,.425).
  \]
  The purchase cap binds; saving, retention, estate, and owner-cap inequalities are strict; \(MV^Y=.768>.6\). The polynomial has \(\Lambda\simeq32.485\) and a stable complex pair, and (T12) gives
  \[
  m\simeq-0.11952631.
  \]
  Thus a small credit reform raises terminal \(Y^*\) but lowers initial fertility. It does not establish a negative fertility response at every date, monotone population adjustment, or any welfare conclusion.

Minimal correction: replace “the weighted operator perturbation proof in section 5 applies” with “the Section 5 argument applies after choosing \(\eta\in(\rho_s,1)\), where \(\rho_s\) is this configuration’s stable spectral radius, and retaining its \(C^1\) actual-pre-reform-old boundary and uniformly bounded operator assumptions.” The theorem itself need not be weakened.