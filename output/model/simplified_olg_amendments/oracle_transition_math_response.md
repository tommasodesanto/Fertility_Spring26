I found no substantive error in the household reduction or the actual-inherited-old boundary. The scalar recurrence in Section 4 cannot simply be reused: positive child goods costs break proportionality between fertility and housing, and endogenous tenure introduces future rental costs into current choices. That is an obstruction to the **scalar proof method**, not to the proposed transition.

**An explicit finite extension is possible at the genuine mixed-tenure example.** Below is a computer-assisted proposition with an infinite-horizon proof. I constructed and ran a standalone exact-arithmetic verifier; it does not require the missing repository helpers.

[Full proof note](sandbox:/mnt/data/finite_mixed_certificate/FINITE_MIXED_TRANSITION.md) · [Proof, executable verifier, coefficients, and certificate](sandbox:/mnt/data/finite_mixed_transition_proof.zip)

## 1. A finite mixed-tenure proposition

Take exactly the original TeX’s example, with taste scale \(\sigma_\xi=4\). In particular,

$$
\vartheta_0=\frac{141}{400},\qquad
\phi_0=\frac45,\qquad
\chi=\frac3{20},\qquad
\tau^p=\frac{467}{9250}.
$$

Its initial stationary equilibrium has \(P=Y=O=1\), owner share \(11/21\), and the conditional allocations stated in the packet. Set the taste location at that reference equilibrium and then **hold both taste parameters fixed**. All remaining primitives are unchanged.

**Proposition.** For every

$$
\boxed{
0<\epsilon\le\frac1{20000},\qquad
0<\delta\le\frac1{20000},\qquad
t_p\ge1,
}
$$

let the unexpected permanent decline to \(\vartheta_1=\vartheta_0-\epsilon\) occur at date zero. At \(t_p\), compare continuation of that baseline with an unexpected permanent increase to \(\phi_1=\phi_0+\delta\), starting from exactly the same actual inherited state.

Both infinite equilibrium paths exist, are unique within the neighborhood specified below for their respective inherited data, and converge exponentially to positive stationary endpoints. Every original household constraint is satisfied. Moreover,

$$
\boxed{
\bar n_0^B<\frac12,\qquad
N_B^*<2,\qquad
\bar n_{t_p}^{1}>\bar n_{t_p}^{B},\qquad
N_1^*>N_B^*,
}
$$

where \(N^*=Y^*+O^*=2Y^*\).

Throughout these paths,

$$
0.523<\pi_t^O<0.525.
$$

Thus **more than 47.5% rent**, while child goods costs and rebated property taxation remain strictly positive.

The limitation is quantitative scope: the largest certified credit change is from **80% to 80.005% financing**. This is an explicit finite rectangle, not an economically large reform or a broad primitive-inequality theorem.

## 2. Use an infinite boundary-value operator

Let

$$
X_t=(P_t,Y_{t+1}),\qquad t\ge0,
$$

with \(Y_0,O_0\) inherited and \(O_t=Y_{t-1}\) thereafter. Define

$$
u_t=(1+q\tau^p)P_t-qP_{t+1},\qquad
T_t=\frac{q\tau^pP_t\bar H}{Y_t+O_t},
$$

and obtain conditional choices and tenure probabilities from the original household problems.

For clarity, those conditional equations remain

$$
h_t^O=\frac{b}{(1-\phi)P_t},\qquad h_t^R=a,
$$

$$
\rho_Ox_t^O+\chi n_t^O=w_t-u_th_t^O,\qquad
\rho_Rx_t^R+\chi n_t^R=w_t-au_t-qa u_{t+1},
$$

$$
\frac{\vartheta}{n_t^m}
=\frac{\chi}{x_t^m}
+\frac{\alpha\kappa}{h_t^m-\kappa n_t^m},
$$

where \(a=h_R^{\max}\), \(w_t=y+b+T_t+qT_{t+1}\),
\(\rho_O=1+\beta(1+\gamma+\omega_B)\), and
\(\rho_R=1+\beta(1+\omega_B)\). The probabilities use the original \(W^O,W^R\), not replacement preferences.

Write the equilibrium residual as

$$
\mathcal F_t(X)=
\begin{pmatrix}
Y_t\bar h_t^Y+M_t/u_t+R_t-\bar H\\
Y_{t+1}-\nu Y_t\bar n_t
\end{pmatrix},
$$

with, for \(t\ge1\),

$$
M_t=\frac{\beta\gamma}{q}Y_{t-1}\pi_{t-1}^Ox_{t-1}^O,
\qquad
R_t=aY_{t-1}(1-\pi_{t-1}^O).
$$

### The initial old are not regenerated

Use the inherited vector

$$
I=(Y_0,O_0,p^-,a_O^-,H_O^-,a_R^-),
$$

where the last four entries describe the old tenure share, per-owner net financial assets and purchased title, and per-renter financial assets. At the initial date,

$$
M_0=
\frac{\gamma}{1+\gamma+\omega_B}
\,O_0p^-\bigl(a_O^-+P_0H_O^-+T_0\bigr),
\qquad
R_0=aO_0(1-p^-).
$$

In particular, at the later reform,

$$
a_O^-=
\frac{a_{t_p-1}^{\prime,B}
-\phi_0P_{t_p-1}^{B}h_{t_p-1}^{O,B}}{q}.
$$

Saving on the right was chosen under the **baseline forecast**. It is not recalculated using the reform continuation. The new price revalues \(H_O^-\), but the new financed share does not rewrite the mortgage. This is precisely the packet’s actual-old boundary.

## 3. The finite certificate

Use the supremum norm on bounded two-component sequences. Define

$$
R=\frac1{2000},\qquad
R_B=\frac1{5000},\qquad
R_P=\frac1{4000}.
$$

The verification covers the entire box

$$
\|X-\mathbf1\|_\infty\le R,\qquad
\vartheta\in[\vartheta_0-1/20000,\vartheta_0],\qquad
\phi\in[\phi_0,\phi_0+1/20000],
$$

and the following inherited-state box. Its decimal endpoints are exact rational numbers.

| Inherited component |              Interval |
| ------------------- | --------------------: |
| \(Y_0,O_0\), each   |   \([0.9998,1.0002]\) |
| \(p^-\)             | \([0.52349,0.52413]\) |
| \(a_O^-\)           | \([0.34569,0.34904]\) |
| \(H_O^-\)           | \([0.99980,1.00021]\) |
| \(a_R^-\)           | \([1.61670,1.61740]\) |

The original stationary inherited data lie in this box. More importantly, **every inherited state generated by a baseline inside radius \(R_B\) lies in it**, using the original budgets and pre-reform forecasts.

Let

$$
\mathcal L=D_X\mathcal F
$$

at the reference stationary economy, including its actual-old boundary. The accompanying coefficient file specifies an exact dyadic, finite-band operator \(A\). Its first 24 block rows are individually specified; subsequent rows are translates of a fixed stencil with offsets \(-24,\ldots,24\).

Exact rational calculations establish the half-line operator bounds

$$
\|I-A\mathcal L\|<2\times10^{-9},\qquad
\|I-\mathcal LA\|<3\times10^{-9},\qquad
\|A\|<4.01.
$$

The second inequality makes \(A\) injective.

Throughout the **whole nonlinear box**,

$$
\boxed{\|I-A D_X\mathcal F\|<\frac1{20}},
$$

and the particular row determining the first subsequent young cohort satisfies

$$
\|(I-A D_X\mathcal F)_{(0,Y),\cdot}\|_1<\frac2{125}.
$$

For \(c_\lambda=-A\mathcal F_\lambda\), with inherited data fixed, the uniform bounds are:

| Quantity                        | \(\lambda=\vartheta\) | \(\lambda=\phi\) |
| ------------------------------- | --------------------: | ---------------: |
| \(\|c_\lambda\|_\infty\)        |             \(<10/3\) |     \(<451/100\) |
| First \(Y\) component           |          \(>839/500\) |      \(>29/250\) |
| Stationary-tail \(Y\) component |          \(>237/100\) |     \(>109/100\) |

These are interval enclosures, not evaluations at sampled economies. The checker uses rational endpoints, outward rounding, rational Taylor remainders for logarithms and exponentials, and integer-square-root enclosures. A numerical inverse was used only to **propose** \(A\); its accepted properties are verified exactly on the infinite half-line.

## 4. Proof of existence, convergence, and signs

Define

$$
\mathcal T(X)=X-A\mathcal F(X).
$$

It contracts with constant below \(1/20\). Because \(A\) is injective, its fixed points are exactly the zeros of the original equilibrium residual.

**Baseline.** At the original constant sequence, the preconditioned residual after the preference decline is bounded by \((10/3)\epsilon\). Since

$$
\frac{10}{3}\frac1{20000}+\frac1{20}R_B<R_B,
$$

\(\mathcal T\) maps the radius-\(R_B\) ball into itself. Banach’s theorem gives the entire infinite baseline.

**Later reform.** Fix any actual baseline date and its inherited state, and center the new problem on the baseline tail. At \(\delta=0\), its residual is exactly zero. At \(\delta>0\), its preconditioned residual is at most \((451/100)\delta\). Now

$$
\frac{451}{100}\frac1{20000}+\frac1{20}R_P<R_P,
\qquad
R_B+R_P<R.
$$

Thus the same contraction gives the reform continuation. Nothing in these inequalities depends on \(t_p\).

**Infinite tail.** The operator has a finite, time-homogeneous tail and therefore preserves convergent sequences. Starting iteration at the constant reference sequence, then at the convergent baseline tail, proves convergence of both fixed points because convergent sequences form a closed subspace of the supremum-norm space.

There is also an explicit exponential estimate. Beyond date 26, an operator row depends only on coordinates between \(t-27\) and \(t+26\). For

$$
E_j=\sup_{t\ge j}\|X_t-X^*\|_\infty,
$$

the contraction gives

$$
E_j\le\frac1{20}E_{j-27}\quad(j\ge27),
$$

hence

$$
\boxed{
E_j\le\frac1{1000}
\left(\frac1{20}\right)^{\lfloor j/27\rfloor}.
}
$$

No terminal date or terminal allocation is imposed. This contraction is a boundary-value proof device—not a claim that the six-dimensional forward map contracts.

**Comparative statics.** Differentiation at fixed inherited data gives

$$
w_\lambda=K_Xw_\lambda+c_\lambda,\qquad
K_X=I-A D_X\mathcal F,
$$

so

$$
\|w_\lambda\|_\infty
\le\frac{\|c_\lambda\|_\infty}{1-1/20}.
$$

The first-row and stationary-tail bounds imply, uniformly,

$$
\boxed{
\partial_\vartheta Y_1>\frac85,\qquad
\partial_\phi Y_1>\frac1{25},\qquad
\partial_\vartheta Y^*>2,\qquad
\partial_\phi Y^*>\frac45.
}
$$

For example, the credit-impact lower bound is explicitly

$$
\frac{29}{250}
-\frac2{125}\frac{451/100}{19/20}
=\frac{951}{23750}>\frac1{25}.
$$

Taking constant sequences gives the corresponding invertibility estimate for the stationary system, so these endpoint derivatives do not rely on an unjustified interchange of differentiation and a time limit.

Integrate the preference derivatives with the original inherited state fixed, and the credit derivatives with the actual intervention-date state fixed. Since

$$
\bar n_0=\frac{Y_1}{\nu Y_0},
\qquad N^*=2Y^*,
$$

all four inequalities in the proposition follow.

## 5. Household feasibility and interpretation

The interval verification establishes the original branches throughout, rather than assuming them. Uniform lower margins include owner and renter saving above \(0.972\) and \(0.808\); young-owner housing value minus service cost above \(0.111\); the young-renter cap gap above \(3.36\); and the old-renter cap gap above \(0.757\). Old-owner retention and estate slacks exceed \(0.541\) and \(0.180\), respectively, **including the actual initial old**.

All log arguments and old resources remain positive. The owner physical cap is slack, the down payment holds with equality, and its restriction is strict. Generated estates use the original timing

$$
e_{t+1}=\frac{\beta\omega_Bx_t}{q^2}.
$$

Conditional concavity and the verified first-order-condition margins establish household optimality; the unchanged logistic rule establishes tenure choice. Homogeneous entrants make these checks individual-household checks, not merely checks of aggregate means.

Both endpoints have replacement fertility. Exponential convergence also yields

$$
\sum_{t=t_p}^{\infty}
\log\frac{\bar n_t^1}{\bar n_t^B}
=
\log\frac{N_1^*}{N_B^*}>0.
$$

Neither a fertility ordering at every future date nor a welfare improvement follows.

**Recommendation:** retain the short mixed-tenure proposition in the illustrative model and place this finite certificate in an appendix. A broad theorem under simple primitive inequalities, economically large shocks, and transitions across binding-constraint changes should remain separate extensions.
