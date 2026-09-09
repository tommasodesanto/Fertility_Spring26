# Independent second review: all-share transition and k = 6 welfare

September 9, 2026 UTC. Bounded analytical review of
`overnight_binding_ownership_extension.md` and
`overnight_binding_welfare_independent.md`. No numerical roots, model runs,
browser use, or builds. Only this review file was written.

## Verdict and scope

**Pass: no algebra failure or fatal inverse gap found.** I independently
expanded the cubic and its unit-pair expression, reconstructed the generating-
function boundary equation, and derived both all-share welfare coefficient
polynomials. The stated lower bounds are valid.

| Result | Checked scope |
|---|---|
| Cubic has two stable roots and one unstable root | Every \(k\ge4\), \(0\le\pi\le1\) |
| Unstable-root and boundary-determinant bounds | Same interval |
| Bounded inverse on \(\ell^\infty\) and \(c\) | Same fixed parameter family; actual logistic calibration uses \(0<\pi<1\) |
| Tax impact fertility increase and purchase-price decrease | \(4\le k\le10\), every interior share |
| Positive impact service-price derivative | Same signed-tax family |
| Preference impact fertility sign | Every \(k\ge4\), every interior share |
| Both current cohort welfare averages increase | The specified \(k=6\) family, every interior share |

This is a working extension in a specialized parameter family. It is distinct
from ProT: here \(\epsilon=q=\alpha/(1+\alpha)=1/2\), so ProT's strict
region does not apply. The extension's independent root and boundary proof
is therefore necessary.

The price paths also differ. Here the tax lowers the impact purchase price
but raises its terminal level, while the impact service rent rises. ProT's
selected region instead has a lower terminal purchase price. These diagram
interpretations should not be combined.

## 1. Household and coefficient checks

The reference choices satisfy the original budgets and all stated regimes.
A young owner purchases \(h_i^y=w_i/2\) at price two, borrows principal
\(w_i/2\), and owes face value \(w_i\) on entering old age. The unchanged
reference sale price then offsets that face debt. The finance multiplier is
\[
\mu_i=\frac{3-4/k}{w_i}>0.
\]
Old housing exceeds young housing for \(k\ge4\), so the stated renter
cap bound also covers both conditional young choices. The old financial-
estate floor is strictly slack.

Direct substitution into the earlier budget-derived response moments gives
\[
B_f=\frac{2}{3k},\quad d_n=\frac14,\quad
d_h=\frac54-\frac{2}{3k},\quad d_o=\frac34.
\]
No cross-sectional variation is suppressed: proportional old income makes
the required ratio moment constant while the income distribution remains
heterogeneous.

At \(k=6\), the terminal estate is \(6w_i\), of which \(3w_i\) is the
terminal financial component and \(3w_i\) is housing liquidation proceeds.
Current old financial saving is \(q\cdot3w_i=3w_i/2\). Thus the note's
phrase “financial estates” is consistent if it refers to the terminal
component, rather than current saving.

## 2. Independent cubic and root-count verification

Expanding the earlier characteristic equation and multiplying by \(24k\)
gives exactly the four coefficients in equation (6). The leading coefficient
is uniformly negative. Direct substitution verifies
\[
C(1)=9k(k+2),\qquad
\min_{\pi\in[0,1]}C(-1)=117k^2-30k+64.
\]

A nonreal unit pair requires
\(a_3a_1-a_0a_2-a_3^2+a_0^2=0\). Independent expansion gives the three
reported coefficients of \(3kG\). For example, its quadratic coefficient
in \(\pi\) factors as
\[
9k(k+2)(9k^2-12k+16)
=3k(27k^3+18k^2-24k+96).
\]
The reported expression for \(G(1,k)\), after writing \(k=4+u\), is
correct and strictly negative. Since its first two coefficients are positive,
\(G\) has no zero on the share interval.

At zero ownership the factorization in the source is exact, and its quadratic
has strict Jury inequalities. No root can subsequently cross the unit circle,
and the degree cannot drop. This establishes the root count over the entire
share interval. The unstable root is necessarily real and simple: a nonreal
unstable root would require its conjugate, and a multiple unstable root would
also contradict the count with multiplicity.

The evaluations at two and at \(5/2\) are correct. The coefficient of
\(\pi\) in \(C(5/2)\) is positive for \(k\ge4\), so its maximum is
at one. The resulting quadratic is strictly negative with the source's sign.
Hence \(2\le r_u<5/2\), and \(2/5<s\le1/2\).

## 3. Boundary inverse reconstructed directly

Let \(Y(z)=\sum_{t\ge0}y_tz^t\) with \(y_0=0\). Summing the two
linear rows gives
\[
(1-z)Y+(mz-n)\widehat F+nF_0=z\widehat f,
\]
\[
-z(3+bz)Y+(-D+Jz+\pi d_oz^2)\widehat F+DF_0=z\widehat g.
\]
Eliminating \(Y\) reproduces equation (12), including both the power of
\(z\) on the fertility forcing and the boundary polynomial
\[
B(z)=(1-z)D+nz(3+bz).
\]
The determinant is precisely \(z^3C(1/z)/(24k)\). Its sole zero inside
the unit disk is \(s\), so analyticity of a bounded sequence's generating
function supplies the unique boundary value (13).

The bounds on \(D\), \(B(s)\), and \(F_0\) are correct. In particular,
\[
D\ge D_{\min}=\frac{3k}{4}+\frac{2}{3k}\ge\frac{19}{6},
\qquad B(s)\ge\frac{19}{12}.
\]
The boundary does not become singular at an ordinary positive share.

After this boundary value is inserted, the determinant numerator vanishes
at \(s\). The coefficients of division by \(z-s\) are exactly the
discounted forward sum stated in the note. The remaining quadratic has both
zeros outside the unit disk and a nonzero constant coefficient, so its causal
inverse has absolutely summable coefficients. Repeated roots are harmless.

Cramer's rule for \(Y\) has the same determinant. Its numerator also
vanishes at \(s\), because \(1-s\ne0\) and the boundary compatibility
condition makes the two rows consistent there. There is no uncancelled
unit-root pole from the first row. These facts establish an inverse for
arbitrary bounded forcing, not only constant policy forcing. Every operation
also preserves convergent sequences.

Thus the original fiscal and price rows can be eliminated and reinstated for
arbitrary residuals, and the full derivative is invertible. The Banach
implicit-function argument applies at each fixed interior logistic share.
The endpoint shares in the polynomial proof are limits, not finite-logistic
calibrations. Shock and inherited-state neighborhood sizes remain local.

The coefficient on the initial price includes the inherited old owners'
capital revaluation. Nearby actual inherited claims and titles enter finite
boundary forcing terms. Individual old regime margins must still be
preserved; a nearby aggregate title total by itself does not establish them.

## 4. Impact signs and stationary comparisons

The expression for \(C_*=mD-nJ\) is exact and nonnegative. Tax forcing
is exactly \(f=(k-2)/8\), \(g=-(15k+22)/16\). The polynomial \(Q\)
is increasing in \(s\), and its upper bound at \(s=1/2\) is negative
throughout \([4,10]\). Consequently \(F_0<0\); all terms in the initial
fertility expression then have the required signs.

For the rent sign, \(2D-J=\pi(d_o-d_h)\le0\), so its product with
the negative \(F_0\) contributes nonnegatively. The remaining lower bound
is positive because
\[
32D_{\min}-(15k+22)=9k-22+\frac{64}{3k}>0
\quad(k\ge4).
\]
This establishes the claimed rent sign on the signed-tax interval; it should
not be presented as a proof of the tax price signs for all larger \(k\).

For the preference derivative, the positive forcing and
\(C_*\le mD\le D/2\) give exactly the lower bound in (20), whose
numerator is positive for \(k\ge4\).

The stationary tax, price, and population formulas check. The reported
log-population derivatives are correct. In particular, the positive
terminal purchase-price derivative is \((k-2)/2\), whereas its impact
derivative is negative in the signed family. Both limiting fertility rates
are replacement fertility. Common-state comparisons require the same actual
inherited obligations and the same subsequent preference sequence.

## 5. Independent k = 6 welfare bounds

Set
\[
Q=207-41\pi,\quad F=P'_0,\quad
P'_1=\frac{(414-68\pi)F+252}{Q}.
\]
The uniform bound \(-3/4<F<0\) is valid. The polynomial used to prove
the lower bound transforms exactly into
\[
\frac54-34u+449u^2-108u^3
\ge\frac54-34u+395u^2>0,
\quad 0\le u\le1/2.
\]
The final discriminant is \(-819\), so its positivity is strict.

Independent elimination of the initial housing and demographic rows gives
\[
\ell_0=\frac{81-41\pi-7\pi F}{Q},\qquad
y_1=\frac{333-167\pi}{2Q}
 -\frac{\pi(193-27\pi)}{4Q}F,
\]
\[
\ell_1=
\frac{1681\pi^2-18081\pi+25758}{Q^2}
-\frac{\pi(8298-1394\pi)}{Q^2}F.
\]
These identities reproduce the source's next-period rent expression.

The old and young cohort envelopes, at the stationary reference, are
\[
\frac{W'_o}{N}=\frac23B_w-\ell_0+\frac\pi3F,
\]
\[
\frac{W'_y}{N}=\frac{10}{3}B_w-\frac32\ell_0
 -\frac12\ell_1-\frac{7\pi}{12}(P'_1-F).
\]
Substitution at \(B_w=1\) gives the exact old rational expression and
the exact young constant and \(F\)-coefficient polynomials in the note.

For old welfare, \(F>-3/4\) bounds the numerator below by
\(171-130\pi+(123/4)\pi^2\), which decreases to \(287/4\).
Its denominator is at most \(621\). Hence
\[
W'_o/N\ge287/2484>0.
\]
For young welfare, the coefficient of \(F\) is nonpositive because
\(2583\pi^2-28334\pi+74691\) decreases to the positive value
\(48940\). Dropping its nonnegative contribution leaves the source's
constant numerator, which decreases to \(316874\), over a denominator
at most \(257094\). Thus
\[
W'_y/N\ge158437/128547>0.
\]
Heterogeneity adds exactly the stated multiples of \(B_w-1\). Both
cohort-average lower bounds therefore hold for every interior share in this
family. They are not individual-household or Pareto bounds.

The earlier zero-share fractions, general Jensen certificate, and conservative
positive-share welfare bound also check. The direct dated-planner allocation
and joint-fertility allocation quoted for this family satisfy their resource
identities and fit the stated caps.

## 6. Presentation limits

The same finite caps work on the specified bounded support; the local path
must remain in the strict borrowing, estate, and cap regimes. The exact
all-share result is specialized to the displayed primitive family. Small
primitive perturbations can preserve the strict conclusions on compact
interior share intervals, but no quantitative neighborhood width is supplied.

The welfare rational formulas use a stationary intervention reference. The
note correctly says not to reuse that pairing literally at a nonstationary
reset; actual cohort measures and inherited titles must be used in the general
envelopes. Strict reference signs persist sufficiently locally.

The reported outside rental-financier loss uses the entire inherited rental
stock, including both ages, and is negative because \(F_0<0\). It is an
external monetary incidence account, not a welfare term with an assigned
utility weight. The result does not establish a Pareto gain or identify
mortgages as the necessary source of the welfare or population effects.
