# Independent review: sharp primitive stability condition

September 5, 2026. Scope: section 10.4, including the Cayley proof, zero/repeated stable roots, actual initial-old boundary, and carryover of the fertility tests. The standalone decay-weight and operator-scope qualification has been incorporated.

**Valid.** Section 10.4 correctly relaxes \(r>1\) to \(r>\ell\) under the maintained strict branch assumptions.

I checked the displayed algebra and the original dated initial-old setup.

- The Cayley transform is correct:
  \[
  (1-v)^3\mathcal P\!\left(\frac{1+v}{1-v}\right)
  =A_3v^3+A_2v^2+\ell(C-1)v-\ell(C+1),
  \]
  with
  \[
  A_3=\ell+(3+q)C-4D,\quad
  A_2=\ell+4D+(7q-3)C.
  \]
  Under \(A_3>0\), Descartes gives exactly one positive \(v\)-root; the transformed polynomial is negative at \(v=0\) and \(8Cq>0\) at \(v=1\), so that root is simple and lies in \((0,1)\), hence gives one simple \(\Lambda>1\). The remaining roots have positive product and negative sum (using \(v_+>0\) and the negative pairwise-sum coefficient), hence negative real parts and therefore map to \(|z|<1\). The excluded point \(v=1\) is not a root.

- Necessity is also correct. Since
  \[
  \mathcal P(-1)=4D-\ell-(3+q)C=-A_3,
  \]
  equality gives a unit root at \(-1\); reversal gives a root below \(-1\), while \(\mathcal P(1)<0\) and positive leading coefficient still give a root above one. Thus precisely two strict stable roots fail.

- The inherited-boundary factor remains positive without \(r>1\):
  \[
  F_1+F_2S
  =1+\frac{C(1+q-qS)-L}{\ell}>1-\frac L\ell>0,
  \]
  because strict stability gives \(S<2\) and \(L<\ell\). At \(r=1\), \(D=\ell C\), so zero is a stable root; when additionally \(q=1/2\), it is double. The generalized sequence is correctly understood as
  \[
  \left.\partial_\lambda\lambda^t\right|_{\lambda=0},
  \]
  i.e. \((0,1,0,\ldots)\), not as an undefined literal \(t0^{t-1}\). The generalized boundary determinant reduces to \(F_1+F_2S\), so no transversality loss occurs.

- There is no new original-recurrence or price-recovery issue. At the stationary point the dated numerator is proportional to \(r-\ell>0\), the user-cost denominator is \(\ell>0\), and locally \(P_t=Kd\,Y_t/Y_{t+1}>0\). The date-zero condition remains the actual pre-reform old cohort condition, rather than recurrence (T1) incorrectly applied at \(t=0\).

- Consequently, (T12), its primitive polynomial test, and \(J>0\) carry over under (T14). T13 properly remains only a sufficient result under its explicit additional assumptions \(q\ge1/2\) and \(r>1\).

Exact correction: none to the mathematics. If §10.4 is meant to stand alone, add this scope phrase to its final sentence: “with a decay weight \(\eta\in(\rho_s,1)\), retaining the Section 5/9 \(C^1\) actual-pre-reform-old boundary and uniformly bounded-operator assumptions, including dated tenure leads.” This preserves local uniqueness only among exponentially convergent nearby sequences; it does not imply global uniqueness or a finite positive-parameter bound.