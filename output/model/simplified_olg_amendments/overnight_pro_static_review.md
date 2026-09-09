# Independent review of the Pro static results

September 9, 2026. **PASS, with minor boundary and wording corrections.**
Reviewed the local oracle_consolidated_theory_source.tex: housing Proposition
“retained capacity” and Appendix B; the combined primitive corollary and
Appendix C; joint fertility Proposition and Appendix D; and the sequential
fertility counterfamily. No browser, numerical model, or build was used.
This review does not certify the transition or welfare appendices.

## 1. Fixed-fertility capacity condition CF: pass

Write \(A_i=H_{d_i}-\kappa n_i\), \(r=\gamma/\alpha\), and
\[
F_y(t)=\int\min\{A_i,t\}\,dQ,\qquad
F_o(t)=\int\min\{H_{d_i},rt\}\,dQ.
\]
The actual full planner satisfies
\(F_y(t_F)+F_o(t_F)=\bar s+\bar h^{o,eq}\).
If \(\bar s<\int A_i\,dQ\), there is a unique \(t_y\) with
\(F_y(t_y)=\bar s\). Monotonicity proves the exact test
\[
H_Y^F>H_Y^{eq}\iff \bar h^{o,eq}>F_o(t_y).
\]
Under CF, \(A_i\ge\bar s\) almost everywhere and \(Q(A_i>\bar s)>0\);
hence \(t_y=\bar s\), and
\(F_o(\bar s)\le r\bar s<\bar h^{o,eq}\).
Therefore \(t_F>\bar s\) and \(F_y(t_F)>\bar s\).
The same strict residual-capacity mass ensures total capacity is not exhausted.

The individual difference is exactly
\[
h_i^{y,F}-h_i^{y,eq}
=\min\{H_{d_i}-h_i^{y,eq},\,t_F-s_i\}.
\]
No \(\alpha\ge\gamma\) assumption is needed for this proposition.
CF is a real restriction on the residual menu. It excludes the earlier
counterexample in which small residual capacities diverted housing toward
the old; it does not invalidate that counterexample.

## 2. Joint fertility condition CJ: pass

CJ requires
\[
\bar c^o\ge\bar x,\qquad
\bar h^o\ge(\gamma/\alpha)\bar s,\qquad
H_R>\bar h^{y,eq},
\]
with at least one resource inequality strict. It also needs no
\(\alpha\ge\gamma\) restriction.

Let \(m=\bar n^{eq}\). The corrected \(n_i^2\) Cauchy–Schwarz calculation gives
\(\vartheta/m\ge\chi/\bar x+\alpha\kappa/\bar s\).
Suppose joint mean fertility does not exceed \(m\). Goods clearing then
gives \(x^J\ge\bar x\). Some tenure has \(n_d^J\le m\), so its fertility
condition implies \(s_d^J\le\bar s\). Its gross house is consequently below
\(\bar s+\kappa m=\bar h^{y,eq}<H_R\), and it is uncapped.

Its fertility is the common unconstrained choice \(n_u\). Capped tenures
have \(r_d\ge\lambda_H\), so they choose \(n_d^J\le n_u\le m\).
Repeating the preceding argument makes every young tenure uncapped.
Joint total housing is then at most
\[
(1+\gamma/\alpha)\bar s+\kappa m.
\]
A strict old-housing gap contradicts clearing. If only the old-consumption
gap is strict, \(x^J>\bar x\) implies \(s_u<\bar s\), again making total joint
housing strictly smaller than the fixed stock. Thus \(\bar n^J>m\).

The young-housing conclusion also checks. If a young cap binds, every young
home is at least \(H_R>\bar h^{y,eq}\). If none binds, old housing is bounded
by \((\gamma/\alpha)s_u\), and clearing plus \(\bar n^J>m\) gives the strict
young aggregate gain.

This is a different sufficient route from our tenure-consumption-cushion
theorem. CJ weakens its age-resource comparisons, but adds a minimum-cap
condition that our theorem did not require. Neither theorem generally
contains the other.

## 3. Mixed-finance benchmark and both estate regimes: pass

For \(\phi=q,\tau^p=0\), both tenures have the same current cash coefficient
\(p\). Write \(E=1+\alpha+\vartheta\) and \(K=1+\gamma+\omega_B\).
With competitive caps slack,
\[
x_i=\min\{w_i/E,(w_i+qv_i)/(E+\beta K)\},\qquad
\frac{c_i^o}{x_i}=\max\{\beta/q,Ev_i/(Kw_i)\}.
\]
The strict borrowing condition is exactly \(qEv_i>\beta Kw_i\).
At equality the multiplier is zero.

For owners, optimizing the old estate gives
\[
\Gamma_O=\min\{\gamma,(\gamma+\omega_B)(1-q)\},\qquad
c^o=z/K,\qquad ph^o=\Gamma_Oc^o.
\]
The two formulas coincide at \(\omega_B(1-q)=q\gamma\). In the strictly
binding-floor regime the owner-minus-renter old value constant is
\[
\Delta=(\gamma+\omega_B)\log(\gamma+\omega_B)
+\gamma\log[(1-q)/\gamma]+\omega_B\log(q/\omega_B).
\]
Thus \(\pi=\operatorname{logit}^{-1}[(\bar\xi+\beta\Delta)/\sigma_\xi]\)
is constant across income types and prices. The weighted primitive inequality
is exactly \(\bar c^o/\bar x>\gamma/\bar\Gamma\), where
\(\bar\Gamma=(1-\pi)\gamma+\pi\Gamma_O\). It gives both required strict
mean resource gaps. Neither side of \(\beta/q=1\) is imposed.

Replacement fixes the displayed positive price; household choices and the
constant tenure probability then fix housing per paired cohort and \(N\).
This establishes uniqueness within the verified uncapped regime only.
The proportional-income nonemptiness construction is valid in either estate
regime because \(\bar\Gamma\) is independent of that income ratio.

## 4. Appendix C active-rental-cap construction: pass

On \(0<n_H<H_R/\kappa\), the high-renter interval is exactly the pair of
required inequalities:
\[
n_H<\frac{\vartheta H_R}{\kappa(\alpha+\vartheta)}
\iff x_H>0,
\]
\[
n_H>\frac{\vartheta pH_R}
{\alpha\chi+p\kappa(\alpha+\vartheta)}
\iff \alpha/s_H>p/x_H.
\]
It is nonempty for \(\chi>0\). The defined \(w_H\) satisfies the current
budget and fertility condition, and the last inequality gives a strictly
positive housing-cap multiplier. Strict concavity certifies the choice.

The uncapped owner solution at that cash endowment must have
\(h_H^O>H_R\). The static Euler identity yields
\[
\alpha/s=pE/w+\eta(1-pH_R/w)>pE/w
\]
at the high renter cap, because \(w-pH_R=c>0\). Hence \(s_H<s_H^O\).
Choosing \(H_O>h_H^O\) makes its owner residual
\(H_O-\kappa n_H^O>s_H^O>s_H\). Low renter residuals exceed \(s_H\) because
\(n_L<n_H\); low owner residuals are larger still.

For any resulting tenure probabilities, high-type space and total housing
are bounded above by \(s_H^O,h_H^O\), respectively. Consequently the two
displayed bounds on its mass \(\delta\) imply
\[
\bar s\le(1-\delta)s_L+\delta s_H^O<s_H,\qquad
\bar h^y\le(1-\delta)h_L+\delta h_H^O<H_R.
\]
These bounds are uniform over the actual endogenous logistic probabilities;
no arbitrary tenure weights were substituted.

The old-income certificate makes current spending strictly optimal:
its marginal utility exceeds \(1/w_{\max}\), whereas the marginal continuation
benefit is at most \(\beta K/(qv_{\min})\).
Thus \(z=v_i\) in both tenures. The further bound
\(\Gamma_Ov_i/(Kp)>H_O\) caps every old household, including renters because
\(\Gamma_R\ge\Gamma_O\). At an owner cap \(H_d\), the remaining estate check
is explicit:
\[
(c^o,e)=
\begin{cases}
((v-pH_d)/(1+\omega_B),\,\omega_Bc^o/q),
&v\ge PH_d(1+q/\omega_B),\\
(v-PH_d,\,PH_d),&v<PH_d(1+q/\omega_B).
\end{cases}
\]
The strict cap test ensures positive consumption; equality at the estate
threshold joins the two formulas. Thus active physical caps do not require
silently imposing one estate regime. The old consumption bound
\(c^o\ge v_i/K>w_i>x_i\)
and \(\alpha\ge\gamma\) give both strict resource comparisons.
Owner mortgage principal is \(qPh^y\), with \(a'=-Ph^y\); the repayment and
resale terms cancel on entering old age. Renter saving is zero.
A finite taste location and positive logistic scale give both tenures
positive probability.

Finally choosing \(\nu=1/\bar n\) using those actual probabilities and
\(N=\bar H/(\bar h^y+\bar h^o)\) closes a positive stationary equilibrium.
The construction legitimately fixes \(\nu\) as part of its primitive family;
it is not a result holding every preassigned \(\nu\) fixed.

## 5. Corrections and comparison with earlier work

1. In the joint-optimum appendix, “a young tenure is capped” implies
   \(r_d\ge\lambda_H\), not necessarily \(r_d>\lambda_H\). Strict inequality
   requires a positive cap multiplier. The proof already uses the weak form.
2. The primitive proof writes \(H_R>\sup_i h_i^y\) from strict cap slackness.
   State uniform slackness on the support, or use essential suprema with
   weak inequalities and retain strictness at the mean. Pointwise slackness
   can approach a cap. CF and CJ still follow: proportional \(s_i,n_i,h_i\)
   give residuals at least \(\operatorname{ess\,sup}s_i\), and the relevant
   expectation is strictly below its cap.
3. In the main discussion, repeat the qualifier \(\beta=q\) in the statement
   that slack finance gives no aggregate age-housing gap. The appendix states
   this correctly; without that qualifier the claim is false.
4. Describe an estate floor with \(\Gamma_O<\gamma\) as strictly binding.
   Equality of the floor with a zero multiplier is already handled correctly
   in the appendix.

The sequential-fertility counterfamily also passes. Its use of homogeneous,
strictly concave-off-rays private fertility gives \(G<m\). With
\(0<\delta<2(m/G-1)\), the actual planner bundles imply
\(\bar n^S\le(1+\delta/2)G<m\), while joint fertility rises. This is a cleaner,
more general bound than our earlier explicit quadratic example, with the
same economic distinction between sequential and joint allocation.
