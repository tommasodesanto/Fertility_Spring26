# Independent overnight check: property taxes and a common-state transition

Read-only mathematical investigation except for this report. September 9, 2026.
No model runs, numerical searches, or browser conversations were used. No
household equation was changed. The restrictions below define an explicit
diagnostic subcase of the consolidated prompt's maintained model.

## Verdict and scope

A nonlinear local transition theorem is established below in an explicit
subcase. A sufficiently small permanent, equally rebated property-tax increase
has a unique nearby bounded perfect-foresight path from the common inherited
state. That path converges, and its limiting household population is higher.
This is a local theorem, not global convergence or the general finance result.

The subcase has uniformly slack young finance, old financial-estate floors,
and both tenure caps. It retains heterogeneous endowments, positive owner and
renter shares, ordinary mortgages, and the original tax/rebate mechanism.
It therefore diagnoses the transition mechanism without proving that a
binding young mortgage is its cause. No restriction \(\phi=q\) is needed
when both financial inequalities are inactive.

The second bounded pass closes the nonlinear gap by constructing a bounded
linear inverse on both convergent and bounded sequence spaces. Initial-old
capital revaluation and bounded-price selection are included explicitly.

## 1. Exact equations in the specified subcase

Write \(p_t=u_t\) for prepaid housing services and let the post-decline
fertility taste be constant \(\vartheta>0\). Put
\[
K=1+\gamma+\omega_B,\qquad D=1+\alpha+\vartheta+\beta K,
\qquad M_0=\overline{y^y+b}+q\overline{y^o}.
\]
For each entering type, the exact unconstrained choices are
\[
x_{it}=\frac{y_i^y+b_i+T_t+q(y_i^o+T_{t+1})}{D},\quad
s_{it}=\frac{\alpha x_{it}}{p_t},\quad
n_{it}=\frac{\vartheta x_{it}}{\chi+\kappa p_t},\quad
z_{i,t+1}=\frac{\beta Kx_{it}}q.
\tag{1}
Old households choose \(c^o=z/K\), \(h^o=\gamma z/(Kp_t)\),
and \(e=\omega_Bz/(Kq)\). The same old value function applies to each
tenure; hence ownership depends only on the logistic taste, with a constant
\(\pi\in(0,1)\). Heterogeneity need not disappear.

Let \(X_t=\overline{x_{it}}\). For dates after the unexpected policy jump,
the preceding young cohort's planned old resources are realized. Thus
\[
\begin{aligned}
X_t&=[M_0+T_t+qT_{t+1}]/D,\\
Y_{t+1}&=\frac{\nu\vartheta X_t}{\chi+\kappa p_t}Y_t,\\
\bar H&=Y_t\left[\frac{\alpha X_t}{p_t}
       +\frac{\kappa\vartheta X_t}{\chi+\kappa p_t}\right]
       +Y_{t-1}\frac{\beta\gamma X_{t-1}}{qp_t},\\
T_t&=\frac{q\tau_tP_t\bar H}{Y_t+Y_{t-1}},\\
p_t&=(1+q\tau_t)P_t-qP_{t+1}.
\end{aligned}
\tag{2}
At the intervention date, the old term must instead use their actual
inherited net resources
\[
\bar Z_0=\bar a_0+P_0\bar H_0+\overline{y^o}+T_0.
\tag{3}
Here \(\bar H_0\) is inherited owned housing per old household, zero for
renters before averaging. Replace the old term in (2) by
\(O_0\gamma\bar Z_0/(Kp_0)\). Omitting (3) erases the initial old's
capital gain or loss and fails the common-inherited-state comparison.

Neither price nor rebate is held fixed. Old estates change endogenously;
creditor obligations inherited at announcement do not. The private financial
restrictions must be checked along the resulting path, even though this
subcase presumes strict margins. A stationary \(\phi=q\) simplification
cannot replace (2): outside stationarity \(L_t-p_t=q(P_{t+1}-P_t)\)
when \(\phi=q\).

## 2. Exact stationary endpoint and its tax derivative

Define
\[
A=\alpha+\beta\gamma/q,\quad B=A+\vartheta,
\quad d_p=1-q+q\tau,\quad
\Psi(\tau)=\frac{(1+q)q\tau}{2d_p}.
\]
The stationary solution in this regime is
\[
X(\tau)=\frac{M_0-\Psi(\tau)\chi/\nu}{D-\Psi(\tau)B},\quad
p(\tau)=\frac{\nu\vartheta X(\tau)-\chi}{\kappa},\quad
P(\tau)=p(\tau)/d_p,
\tag{4}
\]
\[
h_{pair}(\tau)=\frac\kappa\nu+\frac{AX(\tau)}{p(\tau)},\qquad
N(\tau)=\bar H/h_{pair}(\tau).
\tag{5}
These follow from replacement fertility, housing clearing, and the fiscal
identity, not from fixed external wealth. They require positive numerator,
denominator and service price and satisfaction of the maintained inequalities.

At a zero-tax reference \(X=M_0/D\), assume \(\nu\vartheta X>\chi\).
Then
\[
X_\tau(0)=\frac{(1+q)q}{2(1-q)D}
                  \left[BX-\frac\chi\nu\right]>0,
\qquad
\frac{N_\tau(0)}N=
\frac{A\chi}{\kappa p^2h_{pair}}X_\tau(0)>0.
\tag{6}
The positive child-goods cost matters: at replacement fertility, more
resources raise service prices and reduce housing per adult-household pair.
The fixed housing stock then accommodates more households. This effect does
not, by itself, establish an improvement in welfare or identify young finance
as the mechanism. Capital prices need not move in the same direction as rents.

## 3. A verified linearized transition result

At that zero-tax stationary reference, define the following positive housing
quantities and child-cost share:
\[
a=\alpha X/p,\quad b=\beta\gamma X/(qp),\quad c=\kappa/\nu,\quad
\eta=\frac{\kappa p}{\chi+\kappa p}\in(0,1),\quad
d=a+b+\eta c.
\]
Let \(y_t=dY_t/N\), \(\ell_t=dp_t/p\), and \(\xi=dX_t/X\).
For a small permanent tax introduced at date zero, the first-order rebate is
constant across dates: at \(\tau=0\), changes in its price and population
base multiply zero. Consequently
\[
dT_t=\frac{qP\bar H}{2N}\,d\tau,
\qquad \xi=\frac{1+q}{DX}\,dT_t>0.
\]
For \(t\ge1\), linearized housing clearing and demography imply
\[
\ell_t=\frac{(a+c)y_t+by_{t-1}+(a+b+c)\xi}{d},
\]
\[
\boxed{y_{t+1}=\mathcal A y_t-ky_{t-1}+f\xi},\qquad
\mathcal A=1-\frac{\eta(a+c)}d,\quad
k=\frac{\eta b}d,\quad f=\frac{(1-\eta)(a+b)}d.
\tag{7}
Here \(0<\mathcal A<1\), \(0<k<1\). The polynomial
\(r^2-\mathcal A r+k\) satisfies
\(1-\mathcal A+k>0\), \(1+\mathcal A+k>0\), and \(1-k>0\).
Both population roots are therefore strictly inside the unit circle. The
unique forced limit is
\[
y_\infty=\frac{(1-\eta)(a+b)}{\eta(a+b+c)}\xi>0,
\tag{8}
\]
which agrees with (6). This allows damped oscillations; no date-by-date
positive fertility response has been proved.

The initial jump is not supplied by (7). Under the bounded-price condition,
the linear asset-price equation is
\[
dP_0=\sum_{t\ge0}q^t\{p\ell_t-qP\,d\tau\}.
\tag{9}
Let \(G=1-q\mathcal A+q^2k>1-q\). In the homogeneous perturbation problem with
inherited populations fixed, \(y_0=0\), \(y_1=-\eta\ell_0\), and
(7) gives \(dP_0=p(1-q)\ell_0/G\). Equation (3) then leaves the scalar
initial-price coefficient
\[
\mathcal J=d-\frac{\gamma\bar H_0(1-q)}{KG}.
\tag{10}
Thus \(\mathcal J\ne0\) gives a unique bounded linearized jump and,
together with (7), a unique bounded linearized path. The sufficient condition
\(\gamma\bar H_0<Kd\) guarantees \(\mathcal J>0\).
At the stationary reference \(\bar H_0=\pi(a+c)\); a sufficiently small
positive owner share can satisfy it without eliminating either tenure.

These linear calculations are sufficient for the nonlinear local result,
once the following inverse construction is supplied.

## 4. Nonlinear local theorem and inverse construction

**Theorem.** Suppose the positive zero-tax stationary reference has bounded
endowment support, uniform strict margins on all private financial and
housing constraints, and \(\mathcal J\ne0\). For sufficiently small permanent
\(\tau\ge0\), convergent preference paths \(\vartheta_t\) uniformly close
to the reference constant, and admissible nearby inherited states, there is
a unique nearby bounded solution of (1)–(3). It converges to the stationary
solution (4)–(5) evaluated at \(\vartheta_\infty\). Its real allocations,
net financial positions, prices and rebates are locally unique. Against the
zero-tax path from the same state and same preference path, a sufficiently
small positive tax gives a strictly larger limiting \(N\).

Here uniform margins apply separately to both tenures. The unbounded logistic
taste causes no difficulty: real conditional choices are independent of that
draw, and its owner probability stays constant in this regime.

**Spaces and residual map.** Let \(\mathscr C\) be the Banach space of
convergent real sequences with the supremum norm, and let
\(\mathscr B=\ell^\infty\). The five unknown sequences are
\[
(Y_{t+1},p_t,P_t,T_t,X_t)_{t\ge0}.
\]
Fix \(Y_0,O_0\) and initial financial/title endowments. Use the five
residuals in (2), replacing the date-zero housing row by (3), and set
\(D_t=1+\alpha+\vartheta_t+\beta K\). This defines a continuously
differentiable map \(\mathscr C^5\to\mathscr C^5\) near the constant
reference; the same holds on \(\mathscr B^5\). Shifts and the exceptional
date-zero row are bounded operators. Positive denominators are uniformly
bounded away from zero. Nearby inherited distributions can be parameterized
in the essential-supremum norm; only their means enter these aggregate rows.

**Surjectivity and inverse bound.** Take any five residual sequences \(r\).
At \(\tau=0\), the fiscal linearization first gives \(dT=r_T\), and
\(D\,dX_t=dT_t+q\,dT_{t+1}+r_{X,t}\). Hence
\(\|dT\|+\|dX\|\le C\|r\|\).
Eliminating post-zero housing prices from housing and demography yields
\[
y_{t+1}=\mathcal A y_t-k y_{t-1}+u_t,\quad t\ge1,
\qquad \|u\|\le C\|r\|,
\]
where \(u\) converges when \(r\) does. For an initially unspecified
\(y_1\), the companion matrix
\[
M=\begin{pmatrix}\mathcal A&-k\\1&0\end{pmatrix}
\]
has \(\|M^j\|\le C_0r_0^j\) for some \(r_0<1\), by (7).
The variation-of-constants formula therefore gives
\[
\|y\|\le C_1(|y_1|+\|r\|).
\]
For convergent forcing it converges to the forced fixed point, because
the geometric convolution of a sequence tending to zero tends to zero.
Housing then determines \(\ell_{t\ge1}\) with the same bound.

The date-zero demographic row gives \(y_1=-\eta\ell_0+v(r)\).
The price row has the unique bounded inverse
\[
dP_t=\sum_{j\ge0}q^j[p\ell_{t+j}-r_{P,t+j}],\qquad
\|dP\|\le(p\|\ell\|+\|r_P\|)/(1-q).
\]
Its date-zero dependence on \(\ell_0\) is exactly (9)–(10); all remaining
terms are bounded linear functions of \(r\). The initial housing row is
therefore \(\mathcal J\ell_0=L(r)\), with
\(|L(r)|\le C_2\|r\|\). Since \(\mathcal J\ne0\), this uniquely
determines \(\ell_0\), and all five unknowns satisfy
\[
\|du\|\le C_3(1+|\mathcal J|^{-1})\|r\|.
\]
This proves a bounded bijective derivative on BOTH sequence spaces.

**Implicit-function step.** The Banach implicit-function theorem on
\(\mathscr C\) supplies an actual convergent solution, differentiable in
the tax, preference sequence and inherited state. The theorem on
\(\mathscr B\) makes that solution unique among nearby bounded paths;
convergence is therefore not merely assumed of an otherwise unselected
local bounded equilibrium. Uniform margins preserve the original household
optima throughout. Taking limits in the exact residuals gives (4)–(5).
Their positive tax derivative (6), continuously maintained near the reference,
proves the strict endpoint comparison. The algebra may be smoothly extended
to small signed taxes for the theorem; only \(\tau\ge0\) is reported.

## 5. Exogenous fertility decline and intervention along its path

At zero tax let \(D=D_0+\vartheta\),
\(D_0=1+\alpha+\beta K\). Then
\[
\frac Xp=\frac{\kappa M_0}{\nu\vartheta M_0-\chi D},\qquad
\frac{d(X/p)}{d\vartheta}
=-\frac{\kappa M_0(\nu M_0-\chi)}
       {(\nu\vartheta M_0-\chi D)^2}<0.
\tag{11}
\]
Price positivity implies \(\nu M_0>\chi\), and \(A\) in (4) is
independent of \(\vartheta\); hence \(N_\vartheta>0\). A small permanent
preference decline has a lower stationary population endpoint.

An initial fertility decline itself also has a precise sufficient condition.
Put \(j=\vartheta/D\in(0,1)\), \(\zeta=d\vartheta/\vartheta\),
\(\xi=-j\zeta\), and \(g=(1-j)\zeta\). At zero inherited owner exposure,
\[
\ell_0=(a\xi+cg)/d,\qquad
y_1=\frac{(a+b)(1-j)+\eta a j}{d}\,\zeta.
\tag{12}
\]
Thus a negative preference shock immediately reduces fertility. For positive
\(\pi\), set \(F_\vartheta=(a+b)[1-(1-\eta)j]/d\) and
\(L_\vartheta=(a+b)(1-j)+\eta a j\). The same discounted-price calculation
including initial capital revaluation gives
\[
\frac{y_1}{\zeta}=
\frac{L_\vartheta-\pi\gamma(a+c)K^{-1}
 [(1-j)/(1-q)-qF_\vartheta/G]}{\mathcal J}.
\]
The explicit sufficient restriction
\[
\boxed{\pi\gamma(a+c)<K(1-q)(a+b)}
\tag{13}
\]
makes both numerator and \(\mathcal J\) strictly positive. This follows
by dropping the negative \(-qF_\vartheta/G\) term when bounding the
bracket above. It admits strictly positive owner and renter shares.
Differentiability then gives the same initial sign for sufficiently small
finite preference declines; this is an analytic bound, not a numerical example.

A small ongoing preference transition is a \(\mathscr C\)-parameter, with
nearby inherited states generated by its baseline history. Resetting time at
an intervention date gives baseline and tax paths from exactly that state
and the same remaining preference sequence. The local theorem proves that
they converge to their respective endpoints, and (6) orders the endpoints.
It does not require fertility to be higher at every policy date.

## 6. Nonemptiness and remaining scope

The reference subcase is nonempty: take bounded heterogeneous positive
current resources, sufficiently small positive old incomes to preserve
financial saving, \(\phi\ge q\),
\(\omega_B(1-q)>q\gamma\), sufficiently small \(\chi>0\), and finite
caps above the computed demands. Choose an interior logistic owner share
satisfying (10). Every \(\beta>0\) can be considered, but the income and
price-positivity conditions depend on it. No standalone patience restriction
is inferred.

A large transition, an already positive baseline tax, or active household
constraints needs additional work. The theorem uses an explicit bounded-price
selection condition, not a claim that the original equations exclude every
explosive asset-price solution. Gross saving/borrowing decompositions can also
remain nonunique when the same net position implements the real allocation.

The result is not a welfare theorem. Under the maintained outside rental-
financing convention, unexpected revaluation of inherited rental titles,
\(\Delta P_0 H^{rent}_{-1}\), belongs to the outside residual financier's
account, not the living-household utility sum. Equal tax rebates are transfers,
not newly produced resources. The no-arbitrage equation alone does not prove
solvency of a different interpretation with zero-equity intermediaries and
fixed face-value debt: that interpretation needs an explicit loss-bearing
account before asserting all creditor payments are honored. No domestic
landlord wealth or bailout is introduced here. The local equilibrium theorem
inherits, rather than separately derives, the existing outside settlement.

For an actual existing pair of policy and baseline paths,
\[
\log\frac{Y_T^P}{Y_T^B}=
\sum_{t=0}^{T-1}\log\frac{\bar n_t^P}{\bar n_t^B}.
\]
Finite positive limiting cohort masses turn its limit into
\(\log(N^P/N^B)\). Total adult households then converge to \(2N\);
these are not resident-person counts. The theorem supplies the local path;
the identity accounts for its population effect rather than proving its sign.

**One precise follow-up for Pro:** Audit the inverse construction above, then
identify whether a corresponding locally convergent permanent-tax result
survives a positive mass with binding young finance. Keep original old
financial constraints, initial capital revaluation, future rebates and bounded
price selection explicit. Do not replace the common-state path by two
stationary allocations or attribute this unrestricted-finance result to a
mortgage mechanism.
