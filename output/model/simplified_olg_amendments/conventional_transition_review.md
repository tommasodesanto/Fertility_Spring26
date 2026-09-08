# Transitions under conventional mortgage finance

Bounded analytical review, September 8, 2026. This review concerns the proposed
ordinary mortgage covenant, known old-age income, and old owners' ability to
resize within the owner cap. It does not amend the proposal or the protected
author manuscript. No model simulation, numerical anchor, or numerical
eigenvalue calculation is used.

The bounds below maintain \(q\in(0,1)\), \(\phi_t\in(0,1)\),
\(\tau^p\ge0\), and the proposal's strictly positive preference weights,
goods requirement \(\chi\), and space requirement \(\kappa\).

**Result.** The dated equilibrium system, several finite population and price
bounds, and stationary endpoint comparative-static identities can be derived
exactly while retaining heterogeneous entrants and logistic tenure choice.
They do not yet give a primitive theorem establishing a convergent transition,
a positive fertility effect at every future date, or a higher final population
after a credit reform. Stationary existence and the fixed-price fertility
theorems leave identifiable, substantive gaps; those gaps are listed below.

## 1. The policy experiment and inherited state

Write a type as \(i=(w_i^0,v_i^0)\), where
\(w_i^0=y_i^y+b_i>0\) is income plus liquid wealth available when young and
\(v_i^0=y_i^o\ge0\) is known old-age income. Keep their fixed, potentially
nondegenerate joint distribution \(F\). Separating \(y_i^y\) from \(b_i\)
does not change the proposed household feasible set because both are available
at purchase. It can still matter for empirical interpretation.

A permanent preference decline is \(\vartheta_t=\vartheta_-\) before date
zero and \(\vartheta_t=\vartheta_+<\vartheta_-\) from date zero onward.
The initial state \((Y_0,O_0,G_0)\) is inherited from the original stationary
equilibrium. A later unanticipated credit reform at date \(k\) starts from
the actual state \((Y_k,O_k,G_k)\) on the no-reform preference-shock path.
It changes \(\phi_t\) for the newly young from its announced implementation
date onward. It does not refinance the old cohort's contracted \(a\) or
change its inherited title \(H\).

The normalized old distribution must retain \((a,H,v^0,m)\), where \(m\)
is lifetime tenure. Omitting \(v^0\), replacing actual titles by a stationary
title distribution, or resetting old net debt after the reform changes this
experiment. Future entrant wealth remains exogenous; estates do not replenish
the entrant distribution.

## 2. Exact dated household equations

Define the beginning-of-date housing-service price and owner cash requirement
by
\[
p_t=(1+q\tau^p)P_t-qP_{t+1}>0,\qquad
L_t=(1-\phi_t+q\tau^p)P_t>0.
\]
Let \(w_{it}=w_i^0+T_t\), \(v_{i,t+1}=v_i^0+T_{t+1}\), and
\(M_{it}=w_{it}+qv_{i,t+1}\). The old value functions at date \(t\) are
\[
\begin{aligned}
\mathcal V_t^R(z)
 &=\max\{\log c^2+\gamma\log h^2+\omega_B\log e:
 c^2+p_th^2+qe=z,\quad h^2\le H_R\},\\
\mathcal V_t^O(z)
 &=\max\{\log c^2+\gamma\log h^2+\omega_B\log e:
 c^2+p_th^2+qe=z,\quad h^2\le H_O,\quad e\ge P_{t+1}h^2\}.
\end{aligned}                                                    \tag{1}
\]
All log arguments are positive, and \(H_R<H_O\) are the physical caps.
An inherited renter has \(z=a+v^0+T_t\); an inherited owner has
\(z=a+P_tH+v^0+T_t\).

Young owners solve
\[
\begin{aligned}
W^O_{it}=\max\;&\log(c-\chi n)+\alpha\log(h-\kappa n)
       +\vartheta_t\log n+\beta\mathcal V_{t+1}^O(z),\\
&c+p_th+qz=M_{it},\qquad c+L_th\le w_{it},\qquad h\le H_O.
\end{aligned}                                                    \tag{2}
\]
Recover their next-date net assets as
\[
a_{it}^{\prime O}=z-P_{t+1}h-v_{i,t+1}.
                                                                    \tag{3}
\]
Young renters solve (2) with \(\mathcal V^R\), cap \(H_R\), and financial
restriction \(z\ge v_{i,t+1}\); their next-date assets are
\(a_{it}^{\prime R}=z-v_{i,t+1}\). Equivalently their cash inequality is
\(c+p_th\le w_{it}\).

Tenure retains the stated logistic probability,
\[
\pi_{it}=\Lambda\!\left(
 \frac{W^O_{it}-W^R_{it}+\bar\xi}{\sigma_\xi}\right),\qquad
\Lambda(u)=\frac1{1+e^{-u}}.
                                                                    \tag{4}
\]
Thus a young allocation generally depends on
\((P_t,P_{t+1},P_{t+2},T_t,T_{t+1},\vartheta_t,\phi_t)\). The second
future price enters the old-age estate and service-price problem; it must not
be silently dropped when setting a terminal condition.

## 3. Exact cohort, housing, and fiscal equations

For \(j\in\{h,n\}\), set
\[
\bar j_t^Y=\int[(1-\pi_{it})j_{it}^R+\pi_{it}j_{it}^O]\,dF(i).
\]
For the inherited old cohort,
\[
\bar h_t^O=\int h_t^m(a+\mathbf 1_{m=O}P_tH+v^0+T_t)\,dG_t.
\]
Equilibrium requires
\[
\begin{gathered}
Y_t\bar h_t^Y+O_t\bar h_t^O=\bar H,\qquad
T_t=\frac{q\tau^pP_t\bar H}{Y_t+O_t},\\
Y_{t+1}=\nu\bar n_t^YY_t,\qquad O_{t+1}=Y_t.                         \tag{5}
\end{gathered}
\]
The old-state law is the push-forward of \(F\) through the dated young
policies: weight the point
\((a_{it}^{\prime R},0,v_i^0,R)\) by \(1-\pi_{it}\), and the point
\((a_{it}^{\prime O},h_{it}^O,v_i^0,O)\) by \(\pi_{it}\). These weights
sum to one, whereas the old cohort's mass is separately \(O_{t+1}=Y_t\).
The taste distribution need not be carried into old age because old tenure
is retained and the taste was additive to the young tenure value.

After an unanticipated reform, only \(G_k\) is predetermined. Reoptimizing
old housing through (1), including the current title value \(P_kH\), is
essential. Treating \(G_k\) as a price-dependent stationary distribution
would eliminate precisely the inherited-state wealth effect at issue.

## 4. Finite bounds that actually follow from the model

The following are necessary conditions for an equilibrium, not a proof that
an equilibrium exists inside the resulting bounds.

**Price growth and solvency.** Positive service prices imply
\[
0<\frac{P_{t+1}}{P_t}<\frac{1+q\tau^p}{q}.                           \tag{6}
\]
There is no corresponding positive lower growth bound from rental entry.
For inherited owners, old-age feasibility requires
\[
a+P_tH+v^0+T_t>0\quad G_t\text{-almost surely}.                       \tag{7}
\]
Inherited renters analogously require \(a+v^0+T_t>0\).
For every occupied inherited type with \(H>0\), this is the exact restriction
\(P_t>(-a-v^0-T_t)/H\). Uniform strict slack requires a positive margin,
rather than merely weak inequality at the essential supremum of these
thresholds. A convenient sufficient condition following from the preceding
young mortgage covenant is
\[
qP_t\ge\phi_{t-1}P_{t-1},
\]
together with strict positivity of either remaining net resources or
\(v^0+T_t\). This is sufficient, not necessary: old income may repay
underwater debt. Following the announcement, anticipated household choices
themselves enforce future solvency. The difficult inherited-state check is
the surprise date, when the old loans were chosen under a different path.

The household condition \(L_t\ge p_t\) is exactly
\(P_{t+1}/P_t\ge\phi_t/q\). Its stationary version is \(\phi\le q\).
Consequently that particular housing-and-fertility theorem cannot cover a
transition converging to a stationary equilibrium with \(\phi>q\).
The alternative fertility theorem for \(L\le p\) does not require this
lower price-growth restriction.

**Population.** Put \(g_{\max}=\nu H_O/\kappa\). Since every young choice
satisfies \(h>\kappa n\) and \(h\le H_O\),
\[
Y_{t+1}<\frac{\nu\bar H}{\kappa},\qquad
Y_{t+1}<g_{\max}Y_t,\qquad
Y_t+O_t\ge\frac{\bar H}{H_O}.                                      \tag{8}
\]
If date \(t+1\) also clears, then
\[
\boxed{\quad
Y_t>\frac{\bar H}{H_O(1+g_{\max})}.
\quad}                                                           \tag{9}
\]
Indeed, next-date capacity requires
\(\bar H\le H_O(Y_t+Y_{t+1})<H_O(1+g_{\max})Y_t\). Thus every infinite
fully occupied equilibrium has a primitive positive lower cohort bound and,
after the initial date, a primitive finite upper cohort bound. These bounds
also expose a possible failure of existence: a predetermined young cohort
below (9) cannot produce enough next-date households to occupy the stock,
even at the maximum feasible fertility. Vacancy is absent from this model.

**Rebates.** Housing capacity gives
\[
0\le T_t/P_t\le q\tau^pH_O.
                                                                    \tag{10}
\]
For \(t\ge2\), (8) also gives
\(T_t/P_t>q\tau^p\kappa/(2\nu)\) when \(\tau^p>0\). These bounds
control the rebate relative to the price, not the price level or the rebate
level. A finite stationary rebate bound cannot automatically be applied to
a transition with inherited title revaluations and net debt.

## 5. An exact finite-dimensional reduction with heterogeneous entrants

Let \(K=1+\gamma+\omega_B\), \(m_t=(1+q\tau^p)P_t\), and
\(r_t^e=qP_{t+1}/m_t\in(0,1)\). This superscript distinguishes the estate
value share from the original end-of-period rent notation. Direct solution
of (1) gives old housing, even when its size cap binds:
\[
h_t^R(z)=\min\{H_R,g_{Rt}z\},\qquad
h_t^O(z)=\min\{H_O,g_{Ot}z\},                                      \tag{11}
\]
where
\[
g_{Rt}=\frac{\gamma}{Kp_t},\qquad
g_{Ot}=\begin{cases}
\gamma/(Kp_t),&r_t^e\le\omega_B/(\gamma+\omega_B),\\
(\gamma+\omega_B)/(Km_t),&r_t^e>\omega_B/(\gamma+\omega_B).
\end{cases}                                                       \tag{12}
\]
The cap clips the unconstrained housing optimum; its presence may change
whether the estate floor binds at that clipped allocation. Formula (11)
concerns housing alone, not the capped value function.

If old size caps are slack for all relevant types, define tenure masses
\(M_t^m\), tenure-weighted asset means \(A_t^m\), tenure-weighted old-income
means \(V_t^m\), and the inherited owner-title mean \(I_t^O\). All means
are normalized by the whole old cohort, not conditional on tenure. Then
\[
\begin{aligned}
\bar h_t^O={}&g_{Rt}(A_t^R+V_t^R+T_tM_t^R)\\
&+g_{Ot}(A_t^O+P_tI_t^O+V_t^O+T_tM_t^O).                          \tag{13}
\end{aligned}
\]
For example,
\[
\begin{aligned}
M_{t+1}^O&=\int\pi_{it}\,dF,&
A_{t+1}^O&=\int\pi_{it}a_{it}^{\prime O}\,dF,\\
I_{t+1}^O&=\int\pi_{it}h_{it}^O\,dF,&
V_{t+1}^O&=\int\pi_{it}v_i^0\,dF,
\end{aligned}                                                     \tag{14}
\]
with corresponding renter equations. Since \(M^R+M^O=1\) and
\(V^R+V^O=\int v_i^0dF\) for cohorts generated by this entrant law, only
five old-state statistics are independent. This is exact aggregation of
arbitrary \(F\), not a representative-household approximation. If old caps
bind, the distribution inside the clipped expressions in (11) generally
matters and these means alone are insufficient.

With slack old caps, \(\mathcal V_t^m(z)=K\log z+C_{mt}\). Conditional
young real choices therefore depend on \(P_t,P_{t+1},T_t,T_{t+1}\), while
\(P_{t+2}\) affects tenure through the old constants. Their difference is
\[
C_{Ot}-C_{Rt}=
\begin{cases}
0,&r_t^e\le\omega_B/(\gamma+\omega_B),\\
c_*+\gamma\log(1-r_t^e)+\omega_B\log r_t^e,
 &r_t^e>\omega_B/(\gamma+\omega_B),
\end{cases}                                                       \tag{15}
\]
where
\(c_*=(\gamma+\omega_B)\log(\gamma+\omega_B)
-\gamma\log\gamma-\omega_B\log\omega_B\).
The second branch is nonpositive and decreasing in \(r_t^e\); it joins the
first with value and derivative zero. This supplies the explicit forward
price channel through logistic tenure. If both old estate floors and old
size caps are slack, this extra lead cancels from the tenure comparison.

## 6. Preference shocks: a stronger fixed-price sign

At a fixed price and rebate path, fertility is nondecreasing in
\(\vartheta\) for each conditional tenure by revealed preference. This
does not need a cap regime or a restriction on \(\chi/\kappa\).
For differentiable conditional policies, the envelope theorem gives
\[
\pi_{i,\vartheta}
=\frac{\pi_i(1-\pi_i)}{\sigma_\xi}
   \log(n_i^O/n_i^R),
\]
and consequently
\[
\begin{aligned}
\bar n_{\vartheta}=
\int\bigg[&\pi_i n_{i,\vartheta}^O+(1-\pi_i)n_{i,\vartheta}^R\\
&+\frac{\pi_i(1-\pi_i)}{\sigma_\xi}
 (n_i^O-n_i^R)\log(n_i^O/n_i^R)\bigg]\,dF\ge0.                   \tag{16}
\end{aligned}
\]
The selection term has the correct sign regardless of which tenure has
higher fertility. A finite-change proof fixes each household's taste draw
and applies revealed preference to its complete tenure-and-allocation menu:
\((\vartheta_2-\vartheta_1)(\log n_2-\log n_1)\ge0\). Integration
then preserves the ordering, including constraint and tenure changes.

Thus a permanent preference decline has a nonpositive direct fertility effect
at each affected date. It still need not lower equilibrium fertility at
every date, because prices and rebates respond. Current young utility has
no direct dependence on future \(\vartheta\); anticipation of the permanent
shock works through the price and transfer path.

For credit policy, the extra selection term is instead
\(\pi_{i,\phi}(n_i^O-n_i^R)\), and its sign needs the conditions in
`conventional_fertility_review.md`. Equation (16) cannot be recycled as a
credit-policy proof.

## 7. Stationary endpoint equations and exact derivatives

Let \(f(P,T;\vartheta,\phi)\) be stationary average young fertility, and
let \(d(P,T;\vartheta,\phi)\) be young plus old housing per young cohort.
The old distribution in these functions is generated by the same stationary
young policies, including tenure selection, debt and old income. Both
functions retain \(F\); they are not moments at a frozen old distribution.

Every positive stationary equilibrium satisfies
\[
\boxed{\quad f=1/\nu,\qquad
T=kPd,\qquad Y^*=O^*=\bar H/d,\qquad k=q\tau^p/2.\quad}             \tag{17}
\]
The separate primitive stationary-existence result supplies an intersection
of the first two equations under its stated restrictions. It supplies
neither uniqueness nor a sign for their Jacobian.

For a scalar parameter \(j\), such as \(\phi\) or \(\vartheta\),
differentiate (17) along a regular stationary branch:
\[
\begin{bmatrix}
f_P&f_T\\
-k(d+Pd_P)&1-kPd_T
\end{bmatrix}
\begin{bmatrix}P_j\\T_j\end{bmatrix}
=\begin{bmatrix}-f_j\\kPd_j\end{bmatrix}.                          \tag{18}
\]
Let
\(\Delta=f_P(1-kPd_T)+kf_T(d+Pd_P)\ne0\). Then
\[
\begin{aligned}
P_j&=\frac{-f_j(1-kPd_T)-kf_TPd_j}{\Delta},\\
T_j&=\frac{kPf_Pd_j-k(d+Pd_P)f_j}{\Delta},\\
\frac{Y_j^*}{Y^*}&=-\frac{d_PP_j+d_TT_j+d_j}{d}.
\end{aligned}                                                     \tag{19}
\]
With positive tax, the last identity is equivalently
\(Y_j^*/Y^*=P_j/P-T_j/T\). With zero tax, \(T=0\), and
\[
P_j=-f_j/f_P,\qquad
\frac{Y_j^*}{Y^*}=\frac{d_Pf_j/f_P-d_j}{d}.                         \tag{20}
\]
Even \(f_\phi\ge0\), \(f_P<0\), and \(d_P<0\) do not determine the
population sign: one must compare the direct housing response \(d_\phi\)
with the housing reduction caused by the equilibrium price adjustment.

There is an especially transparent exact finite endpoint criterion. Since
stationary fertility equals replacement,
\[
d=\kappa/\nu+\overline{s^Y}+\overline{h^O},\qquad
s^Y=h^Y-\kappa n^Y.                                               \tag{21}
\]
For two stationary equilibria, final adult-household population is higher
if and only if average adult young space plus average old housing is lower.
The replacement children's space per cohort, \(\kappa/\nu\), is unchanged.
This criterion is an accounting identity, not a reason to presume that a
credit reform satisfies it. Both endpoint fertility levels are exactly
\(1/\nu\), even when their population levels differ.

## 8. Finite stationary ordering needs additional monotonicity

A usable analytical route first solves the fiscal equation for
\(T=\mathcal T(P,j)\). On a finite domain, a sufficient regularity condition
is \(b=1-kPd_T>0\), together with fiscal boundary signs ensuring a root.
Then
\[
\mathcal T_P=\frac{k(d+Pd_P)}b,\qquad
\mathcal T_j=\frac{kPd_j}b.
\]
Define the resulting reduced functions \(\tilde f(P,j)\) and
\(\tilde d(P,j)\). If \(\tilde f_P<0\) throughout the domain and
\(\tilde f_j\) has a stated sign, the stationary price comparison follows
for finite changes staying in that domain. Population additionally requires
the sign of
\[
\tilde d_j-\tilde d_P\tilde f_j/\tilde f_P.                         \tag{22}
\]
Uniform signs in (22) can be integrated along a finite regular branch.
These are transparent missing monotonicity restrictions, not established
primitive facts of the current proposal. In particular, the direct
fixed-price theorem does not already sign \(\tilde f_j\), because the
fiscal solution changes the rebate as well.

## 9. Why market monotonicity cannot be presumed

Old owners can have negative net financial wealth. For an uncapped old
owner whose estate floor binds,
\[
h^2=\frac{\gamma+\omega_B}{K(1+q\tau^p)}
 \left[H+\frac{a+v^0+T_t}{P_t}\right].
\]
Holding the future price, transfer and inherited state fixed within that
regime,
\[
\frac{\partial h^2}{\partial P_t}
=-\frac{(\gamma+\omega_B)(a+v^0+T_t)}
        {K(1+q\tau^p)P_t^2}.                                     \tag{23}
\]
This derivative is positive when debt exceeds old income plus the rebate,
even though total old resources can remain positive due to the inherited
house. A price rise then increases old housing demand through the title
wealth effect. With a slack estate floor, the derivative is
\[
\frac{\gamma[-qP_{t+1}H-(1+q\tau^p)(a+v^0+T_t)]}{Kp_t^2},
\]
which can also be positive. Nonnegative \(a+v^0+T_t\) is a sufficient
condition for these old-demand slopes to be nonpositive. It is stronger
than solvency and is not part of the proposal.

There is a useful positive result for the young block. Suppose old size caps
are slack and hold future prices and both relevant rebates fixed. Then each
conditional young housing demand is nonincreasing in current \(P_t\).
For a strictly constrained, uncapped young owner, write
\(B=\beta K\), \(x=c-\chi n\), \(s=h-\kappa n\),
\(S=qz\), \(V=q(v^0+T_{t+1})\ge0\), and \(\Delta_h=L-p\). Define
\[
\begin{aligned}
a_h&=L^2/x^2+\alpha/s^2+B\Delta_h^2/S^2,\\
b_h&=-L\chi/x^2+\alpha\kappa/s^2,\\
d_n&=\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2,\\
J_h&=L+b_h\chi/d_n
=\frac{L\vartheta/n^2+\alpha\kappa(\chi+L\kappa)/s^2}{d_n}>0.
\end{aligned}
\]
Strict concavity gives \(a_hd_n-b_h^2>0\). Since
\(L_{P_t}=\ell=1-\phi+q\tau^p>0\) and
\((\Delta_h)_{P_t}=-\phi\), direct differentiation yields
\[
\boxed{\quad
\frac{\partial h^O}{\partial P_t}
=-\frac{d_n}{a_hd_n-b_h^2}
\left[\frac{\ell(x+hJ_h)}{x^2}+\frac{B\phi V}{S^2}\right]<0.
\quad}                                                          \tag{23a}
\]
The derivative is also negative in the uncapped unconstrained allocation;
a binding young cap makes it zero. Renter formulas give the same weak sign.
Continuity permits the corresponding finite comparison across financial and
young-cap changes, while old caps remain slack.

Conditional fertility need not be inferred merely from that housing sign.
In the constrained regime,
\[
n_{P_t}=\frac{b_hh_{P_t}-\chi\ell h/x^2}{d_n}.                    \tag{23b}
\]
Thus \(b_h\ge0\) suffices for \(n_{P_t}<0\). A sufficient common-price
restriction is
\(\min\{L,p\}^2\ge\alpha L\chi/\kappa\): the effective marginal
space price \(\rho=\alpha x/s\) lies between \(L\) and \(p\), so this
restriction gives \(b_h=(\kappa/x^2)(\rho^2/\alpha-L\chi/\kappa)\ge0\).
These are partial derivatives holding the future price fixed; they are not
the stationary derivative \(f_P\), which moves all future stationary prices.

Endogenous tenure still enters aggregate young demand. If \(\mu_O\) is
the multiplier on the owner's cash constraint, the envelope theorem gives
\[
\partial_{P_t}(W^O-W^R)
=-(1+q\tau^p)\left(\frac{h^O}{x^O}-\frac{h^R}{x^R}\right)
 +\phi\mu_Oh^O.                                                  \tag{23c}
\]
Its sign, and hence the sign of
\(\pi_{P_t}(h^O-h^R)\), is not fixed by (23a). In addition, the fiscal
rebate rises with \(P_t\) at a fixed cohort state. Conditional concavity and
the own-price lemma therefore do not establish aggregate market monotonicity
or a unique clearing price.

## 10. An analytical finite-transition existence route

The following route avoids a numerical reference equilibrium. It is a
conditional proof program; the required uniform market inequalities have
not been established for the proposed primitives.

1. Choose finite positive price bounds, finite rebate bounds, a compact
   support for \(F\), and a uniform positive user-cost margin. Verify initial
   old solvency over the proposed price-rebate region. The household problems
   then have continuous unique conditional real allocations, and logistic
   mixing preserves continuity. Future old states are the exact push-forwards.
2. For any finite horizon, compute household policies, cohorts and old states
   from a candidate price/rebate path, retaining the two future prices needed
   by (1). Prove uniformly that housing excess demand is positive on each
   lower-price face and negative on each upper-price face, and that the fiscal
   update maps the rebate box into itself. These bounds must cover all
   candidate paths in the box, not just equilibrium paths. The necessary
   population bounds (8)-(9) alone do not establish this invariance.
3. A clipped price adjustment and the fiscal update then form a continuous
   self-map of a compact box. Brouwer gives a fixed point, and the strict face
   signs rule out clipped prices, producing a finite-horizon equilibrium.
4. If all these bounds are uniform in the horizon, diagonal subsequences
   produce an infinite equilibrium: each fixed dated equation eventually lies
   away from the terminal boundary and passes to the limit. This establishes
   existence of an infinite path, **not convergence to the chosen stationary
   endpoint**.
5. Convergence needs an additional uniform stability or tail bound, such as
   a proved global contraction in a finite domain or a proved order-preserving
   trapping region whose only invariant limit is the stationary point. A
   finite-horizon terminal price pinned to that point does not itself supply
   such a theorem.

A possible purely analytical implementation uses the finite-dimensional
state in (13)-(14), bounds the derivatives of the explicit household maps
over a declared finite primitive domain, and proves a uniform contraction
or a suitable comparison principle. For example, strict row diagonal
dominance with positive diagonals and nonpositive off-diagonals in an
appropriately ordered residual system supplies a positive inverse and a
comparison theorem. The residual forcing from policy and the inherited
boundary must then also have the required signs. None of these matrix sign
conditions is automatic here. They must not be substituted by a numerical
eigenvalue calculation at a conveniently chosen point.

There is a still narrower exact testbed: with zero tax, slack caps, and
\((w_i^0,v_i^0)=s_i(w_0,v_0)\) for a nondegenerate positive distribution
of \(s_i\), every conditional real allocation scales with \(s_i\), and
conditional values share the additive term
\((1+\alpha+\vartheta+\beta K)\log s_i\). Logistic tenure is therefore
independent of \(s_i\), and aggregation uses \(E[s_i]\) exactly. This
retains wealth heterogeneity and nondegenerate tenure shares, but restricts
all income timing ratios and removes active size segmentation. It could be
an explicitly labeled analytical subcase; it does not establish the general
heterogeneous model. Positive uniform tax rebates break this scale reduction.

## 11. What an all-future fertility claim would have to prove

For two equilibria starting from the same state at policy date \(k\), the
cohort comparison is exactly
\[
\frac{Y_t^{\mathrm{reform}}}{Y_t^{\mathrm{baseline}}}
=\prod_{s=k}^{t-1}
 \frac{\bar n_s^{\mathrm{reform}}}{\bar n_s^{\mathrm{baseline}}}.
                                                                    \tag{24}
\]
Equivalently, log cohort differences are cumulative log fertility
differences. If both paths converge to positive stationary equilibria,
the partial sums in (24) converge to the log endpoint cohort ratio. A
higher endpoint population therefore constrains the cumulative effect;
it does not require positive effects at every date.

In a differentiable path family, the contemporaneous fertility response is
\[
\frac{d\bar n_t}{dj}
=f_{\phi,t}\frac{d\phi_t}{dj}
 +f_{\vartheta,t}\frac{d\vartheta_t}{dj}
 +\sum_{r=0}^{2}f_{P_{t+r},t}\frac{dP_{t+r}}{dj}
 +\sum_{r=0}^{1}f_{T_{t+r},t}\frac{dT_{t+r}}{dj}.                   \tag{25}
\]
There is no omitted direct old-distribution argument in (25): with entrant
\(F\) fixed, the inherited old state affects young fertility through
equilibrium prices and rebates. In the uncapped-old regime, its second
future-price derivative acts only through tenure as in (15).

An all-future sign requires controlling every term in (25) along the actual
inherited-state transition, including its terminal cancellation: for a
permanent policy change, stationary fertility remains \(1/\nu\), so its
total long-run derivative is zero. A uniform strictly positive direct-effect
bound that ignores this cancellation cannot prove the desired asymptotic
comparison. A Pareto-improving feasible reallocation, conditional owner
fertility growth, higher ownership, and a higher stationary population are
four different claims; none automatically implies this dated fertility sign.

## 12. Claim ledger and recommended next step

| Claim | Analytical status |
|---|---|
| Correct dated mortgage, old-income, rebate and cohort accounting | Derived in (1)-(5) |
| Finite necessary population and price-growth bounds | Proved in (6)-(10) |
| Exact moment reduction retaining arbitrary \(F\) and logistic tenure | Proved when old size caps are slack, (11)-(15) |
| Conditional young housing falls with current price, future prices/rebates fixed | Proved with slack old size caps, (23a); aggregate tenure and old-debt effects remain |
| Preference decline weakly lowers fertility at fixed prices, including tenure changes | Proved for finite changes, (16) |
| Stationary endpoint fertility equals replacement | Exact, (17) |
| Higher final population after credit relaxation | Unproved; requires (19), (21) or (22) to have the appropriate sign |
| Existence of a perfect-foresight transition from the inherited state | Unproved; initial solvency and uniform market/fiscal bounds remain |
| Convergence to the new stationary equilibrium | Unproved; stationary existence does not establish stability or tail selection |
| Credit increases fertility at every future date | Unproved; full path response (25) remains unsigned |

The strongest next analytical step is to use the exact moment system, state
finite primitive intervals, and establish uniform market and fiscal signs
there. This targets the unresolved general-equilibrium channel directly.
Until such inequalities are proved, present the household fertility theorem
and stationary existence theorem with their actual scopes, and leave the
transition and population comparison as outstanding results. The old
separate-liquid-wealth down-payment transition formulas do not apply to the
proposed ordinary mortgage covenant.

Verification: direct symbolic differentiation gives zero residuals for the
two equations in the stationary Jacobian solution (19), the value and first
derivative joins in (15), and the housing own-price numerator, cross-score,
and positive-\(J_h\) identities underlying (23a)-(23b). These are algebraic
checks; no equilibrium was computed and no transition result is inferred
from them. Review closed within its 25-minute analytical budget.
