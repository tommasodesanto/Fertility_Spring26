# Binding transition with arbitrary positive ownership shares

September 9, 2026. Exact analytical extension of
`overnight_binding_transition_independent.md`. The earlier memo is unchanged.
Only exact rational polynomial arithmetic was used; no numerical roots,
model runs, or live Pro output were used. This is a local nonlinear result
in a parameterized subcase of the maintained household model.

## Result

The small-ownership restriction can be removed in the family below. For
every logistic owner share \(0<\pi<1\), there is a unique nearby bounded
equilibrium transition, and it converges. A sufficiently small permanent
rebated property tax raises initial fertility and limiting household
population. A sufficiently small permanent fertility-taste decline lowers
both. Initial capital prices fall after the tax, while initial service rents
rise. These are exact local conclusions for all shares, not a numerical
continuation from the renter limit.

The proof has two distinct parts: a cubic calculation establishes two stable
roots and one unstable root; an explicit initial-state determinant then
selects a unique bounded jump. The latter includes initial-old title
revaluation with inherited face obligations held fixed.

## 1. An admissible heterogeneous family

Take the exact log-offset household utility and set
\[
q=\beta=\phi=\frac12,\qquad
\alpha=\gamma=\vartheta=\kappa=\nu=\chi=1,
\qquad \omega_B=2,\qquad \tau=0.
\tag{1}
\]
Let young available wealth \(w_i=y_i^y+b_i\) have any compact positive,
nondegenerate distribution with mean six. Let old income be
\(v_i=k w_i\). The root and boundary proof works for \(k\ge4\); the
signed-impact theorem below uses the useful closed interval
\[
4\le k\le10.                                                     \tag{2}
\]
Both incomes remain heterogeneous. Positive components \(y_i^y,b_i\)
can be chosen to sum to \(w_i\). Logistic taste scale is any positive
number, and its finite location selects any \(0<\pi<1\).

Here \(E=1+\alpha+\vartheta=3\), \(K=1+\gamma+\omega_B=4\). The stationary
service price is one, the capital price is two, and individual choices are
\[
x_i=w_i/3,\quad n_i=w_i/6,\quad h_i^y=w_i/2,
\qquad c_i^o=h_i^o=k w_i/4.
\tag{3}
\]
The cash-constraint condition is \(k>4/3\), uniformly strict in (2).
Every owner borrows the maximum principal \(w_i/2\) and repays \(w_i\)
when old. The old financial-estate floor has strict margin because
\(\omega_B(1-q)-q\gamma=1/2\). Both caps remain finite; it suffices to take
\[
h_R^{\max}>k w_{\max}/4,\qquad h_O^{\max}>h_R^{\max}.
\tag{4}
\]
Strict gaps preserve the regime on nearby paths. The size restrictions are
economically inactive in this diagnostic family. The finance restrictions
are strictly active for every young household, including every owner.

Mean young and old housing are \(3\) and \(3k/2\), respectively, so
\[
H_*=3+3k/2,\qquad N=\bar H/H_*.
\]
At the reference, conditional tenure choices coincide in real quantities.
Their future transition quantities and logistic probabilities are still
computed from the exact household problems; they are not held constant.

## 2. The cubic has the desired root count for every share

Use the first memo's normalized derivatives
\(y_t=dY_t/N\), \(F_t=dP_t/p=dP_t\). In this family,
\[
d=(5+3k)/2,\quad h_y=3,\quad b=3k/2,\quad
\eta=1/2,\quad d_n=1/4,\quad d_h=5/4-2/(3k),\quad d_o=3/4.
\tag{5}
\]
The heterogeneity moment is exactly \(B_f=2/(3k)\), without suppressing
the distribution of \(w_i\). Multiply the cubic in the first memo by
\(24k\) and write it as \(C(r)=a_3r^3+a_2r^2+a_1r+a_0\), where
\[
\begin{aligned}
a_3&=-18k^2-30k+\pi(30k-16),\\
a_2&=54k^2+72k+\pi(-60k+32),\\
a_1&=-45k^2-24k+\pi(9k^2+48k-16),\\
a_0&=18k^2-\pi(9k^2+18k).
\end{aligned}                                                     \tag{6}
\]
For \(k\ge4\), \(0\le\pi\le1\), its degree never drops:
\(a_3\le-18k^2-16<0\). A real root cannot cross the unit circle because
\[
C(1)=9k(k+2)>0,
\]
\[
C(-1)\ge117k^2-30k+64>0.
\tag{7}
\]
A nonreal unit-conjugate pair would require
\(a_3a_1-a_0a_2-a_3^2+a_0^2=0\). Exact expansion gives this expression
as \(3kG(\pi,k)\), with
\[
\begin{aligned}
G={}&(27k^3+18k^2-24k+96)\pi^2\\
 &+(216k^2+264k-32)\pi-54k^3-198k^2-60k.
\end{aligned}
\]
The first two coefficients are positive for \(k\ge4\), so
\(G(\pi,k)\le G(1,k)\). Writing \(u=k-4\ge0\),
\[
G(1,k)=-27u^3-288u^2-828u-368<0.                                 \tag{8}
\]
Thus no unit-conjugate crossing is possible. At \(\pi=0\),
\[
C(r)=-3k(r-2)
 \{2(3k+5)r^2-2(3k+2)r+3k\}.
\]
The quadratic satisfies the strict Jury inequalities. Its two roots are
inside the unit circle and the third root is two. Equations (7)–(8) and
the nonzero leading coefficient preserve this root count throughout
\(0\le\pi\le1\). This is an analytic root-count argument over the whole
share interval, not a check at one numerical point.

The unique unstable root \(r_u\) is real and simple. In fact
\[
C(2)=\pi(9k^2+78k-32)\ge0,
\]
\[
C(5/2)\le C(5/2)|_{\pi=1}
 =-\tfrac94(11k^2-52k+40)<0.
\]
The last quadratic is \(11u^2+36u+8>0\). Hence
\[
2\le r_u<5/2,\qquad 2/5<s:=1/r_u\le1/2.                         \tag{9}
\]
The equality at the upper bound for \(s\) occurs only at zero ownership.

## 3. Initial-state selection and a bounded inverse

For arbitrary reduced forcing sequences \(f_t,g_t\), the linearized system
is
\[
y_{t+1}-y_t+mF_t-nF_{t+1}=f_t,
\]
\[
-DF_{t+1}+JF_t+\pi d_oF_{t-1}-3y_t-by_{t-1}=g_t,                \tag{10}
\]
where
\[
m=(2-\pi)/4,\quad n=(1-\pi)/4,\quad
D=d/2-\pi d_h>0,\quad J=d-\pi(d_h+d_o)>0.
\tag{11}
\]
The inherited state fixes \(y_{-1}=y_0=0\) and \(F_{-1}=0\) in the
homogeneous state-variation problem. The date-zero term \(-\pi d_oF_0\)
already contained in \(JF_0\) is exactly the initial-old capital
revaluation. Thus this boundary condition does not erase the initial old's
title exposure. Nearby changes in inherited populations, face claims, or
titles become finite-date forcing terms; they do not change the homogeneous
boundary operator.

Let \(\widehat f(z)=\sum_{t\ge0}f_tz^t\), and similarly for the other
sequences. Direct generating-function elimination gives
\[
\mathcal D(z)\widehat F(z)+B(z)F_0
 =z(1-z)\widehat g(z)+z^2(3+bz)\widehat f(z),
\tag{12}
\]
\[
\mathcal D(z)=\frac{z^3C(1/z)}{24k},\qquad
B(z)=(1-z)D+nz(3+bz).
\]
The only zero of \(\mathcal D\) in the unit disk is \(s=1/r_u\).
Boundedness therefore requires the unique jump
\[
\boxed{
F_0=\frac{s(1-s)\widehat g(s)+s^2(3+bs)\widehat f(s)}{B(s)}.}
\tag{13}
\]
This is the initial-condition determinant. It is uniformly separated from
zero in the stated family:
\[
B(s)\ge(1-s)D\ge\tfrac12D_{\min}\ge19/12,
\quad D_{\min}=3k/4+2/(3k)\ge19/6.
\tag{14}
\]
No additional initial-state degeneracy occurs at ordinary owner shares.

For completeness, (13) proves more than uniqueness of a constant forcing
solution. For any bounded \(f,g\), it obeys
\[
|F_0|\le
\frac{\|g\|_\infty+(3+3k/4)\|f\|_\infty}{D_{\min}}.
\]
After inserting it into (12), the numerator \(A(z)\) vanishes at \(s\).
Division by \(z-s\) has coefficients
\[
\left[\frac{A(z)}{z-s}\right]_t
 =\sum_{j\ge0}s^j A_{t+j+1}.
\]
This is a bounded discounted forward sum. The remaining quadratic factor
of \(\mathcal D\) has both zeros outside the unit disk, so its causal
inverse has absolutely summable coefficients. Cramer's rule gives the same
factorization for \(\widehat y\); its numerator also vanishes at \(s\).
These formulas construct a bounded inverse on \(\ell^\infty\). They
preserve convergent forcing sequences, and so construct the inverse on
the Banach space \(c\) as well. Repeated stable roots cause no problem:
their polynomial-times-geometric coefficients remain summable.

The full fiscal/price derivative reduces to (10), as in the first memo;
arbitrary fiscal and price residuals just supply bounded \(f,g\).
The strict household regime makes the exact four-equation residual map
continuously differentiable on both sequence spaces. Banach IFT therefore
gives a convergent local nonlinear solution, uniquely among nearby bounded
solutions, for every fixed \(0<\pi<1\). This includes sufficiently small
convergent preference paths and admissible nearby inherited states. The
argument does not assert global uniqueness or exclude distant paths.

## 4. Exact impact formulas and signs

For constant forcing, (13) becomes
\[
F_0=\frac{gs+fs^2(3+bs)/(1-s)}{B(s)},\qquad
F_1=\frac{JF_0-g}{D}.
\tag{15}
\]
Define a nonnegative quantity
\[
C_*:=mD-nJ
 =\frac{\pi\{9k^2+9k(1-\pi)+8\}}{48k}\ge0.
\]
Then initial fertility, equivalently \(y_1\) at this normalization, is
\[
y_1=f-\frac{ng}{D}-\frac{C_*F_0}{D}.                            \tag{16}
\]
Equations (6), (9), and (15) are exact algebraic expressions for the jumps;
they do not require selecting approximate numerical roots.

### Permanent rebated property tax

For a unit tax-rate derivative, the initial rebate derivative is
\(T'=3(k+2)/4\), and
\[
f=(k-2)/8>0,\qquad g=-(15k+22)/16<0.
\tag{17}
\]
The numerator of \(F_0\), multiplied by \(16(1-s)/s>0\), is
\[
Q(s)=-(15k+22)+(21k+10)s+3k(k-2)s^2.
\]
It increases in \(s\). For \(4\le k\le10\),
\[
Q(s)\le Q(1/2)=\tfrac14(3k^2-24k-68)<0.
\tag{18}
\]
Therefore \(P'_0=F_0<0\) for every owner share. Equations (16)–(17)
then give \(y'_1>0\), since all three displayed contributions are
nonnegative and \(f>0\).

The initial service-price derivative is positive:
\[
\ell_0=1+F_0-\tfrac12F_1
=1+\left(1-\frac J{2D}\right)F_0+\frac g{2D}>0.                 \tag{19}
\]
Indeed, \(2D-J=\pi(d_o-d_h)\le0\), and
\[
1+g/(2D)\ge1-\frac{15k+22}{32D_{\min}}>0\quad(k\ge4).
\]
Thus the immediate capital-price decline is not a decline in the service
rent. Both prices adjust endogenously in the proof.

### Permanent fertility-taste change

For a unit derivative of \(\vartheta\) at one, the forcing is
\(f=2/3\), \(g=0\). Hence \(F_0>0\). Also \(C_*\le mD\le D/2\),
so (9), (14)–(16) give the explicit uniform lower bound
\[
\frac{y'_1}{2/3}
\ge1-\frac{3+3k/4}{2D_{\min}}
=\frac{9k^2-36k+16}{18k^2+16}>0\quad(k\ge4).
\tag{20}
\]
A small permanent taste decline therefore lowers initial fertility and
initial capital prices. This statement concerns that specified permanent
perturbation; arbitrary time profiles are not assigned an unconditional
fertility sign.

## 5. Endpoints, a common-state policy comparison, and limits

At a stationary tax, let \(H_*=3+3k/2\). The exact endpoint satisfies
\[
T(\tau)=\frac{12H_*\tau}{24+13\tau},\quad
p(\tau)=1+T(\tau)/3,\quad
N(\tau)=\frac{\bar H[1+T(\tau)/3]}{H_*+(11/12)T(\tau)}.
\tag{21}
\]
In particular,
\[
\left.\partial_\tau\log N\right|_0=\frac{6k+1}{24}>0,
\qquad
\left.\partial_\vartheta\log N\right|_{\tau=0,\vartheta=1}
=\frac{4(5+3k)}{9(k+2)}>0.
\tag{22}
\]
The stationary capital-price derivative is \((k-2)/2>0\), while the
stationary service-price derivative is \((k+2)/4>0\). The negative
initial capital-price derivative thus reverses later. No prices, rebates,
estates, or young future resources are fixed along these paths.

Start from the reference inherited state and introduce a small permanent
taste decline. The theorem gives an actual convergent baseline transition,
with lower initial fertility and a lower endpoint. At a policy date on a
nearby such path, compare zero tax with a small permanent tax, using exactly
the same inherited state and subsequent taste sequence. The local theorem
and its strict derivative signs persist in a sufficiently small parameter
neighborhood. The policy path has higher initial fertility and a larger
limiting household population than that common-state baseline. The policy
is unexpectedly announced at that date and is subsequently anticipated
perfectly. Different pre-policy anticipation would be a different initial-
state comparison.

Both limiting fertility rates are \(1/\nu=1\). The different population
levels come from the accumulated fertility gap along actual transitions:
\[
\log\frac{N^{\rm pol}}{N^{\rm base}}
=\lim_{T\to\infty}\sum_{t=0}^{T-1}
\log\frac{\bar n_t^{\rm pol}}{\bar n_t^{\rm base}}>0.
\]
There is no claim of higher fertility at every date or restoration to the
pre-decline population. Population counts adult households; its stationary
total is \(2N\).

This family has strict mortgages for a positive, unrestricted owner share.
Owner capital gains, collateral finance, endogenous tenure selection, and
initial-old capital losses enter the exact local system. The population
sign also holds in the renter limit, so mortgage finance is not established
as its sole cause. The old/young income ratio is restricted, and both
tenure-specific size caps are slack. The special parameter choice has
\(\beta R_f=1\); it is a diagnostic construction, not a general restriction
on the model. All established inequalities are strict, so at any compact
interior share interval there is a neighborhood of these primitives with
the same conclusions, allowing either side of this discount-return equality.
No quantitative width for that additional primitive neighborhood is claimed.

The ownership extension is not itself a welfare theorem. The separate
`overnight_binding_welfare_independent.md` uses the \(k=6\) jump formulas
for living-household welfare. External rental financiers bear their own
capital losses under the maintained outside-financing convention. Neither
tax rebates nor those losses are silently assigned to extra domestic
resources. Initial household face obligations are preserved; subsequent
estates and financial positions remain endogenous.

### Verification

Exact symbolic checks verified the cubic coefficients, its unit-pair
expression, its reciprocal generating-function identity, the initial
boundary coefficient, and the closed jump formulas. The root count, strict
sign bounds, and bounded-inverse construction above provide the proof;
there were no numerical root evaluations or model simulations.
