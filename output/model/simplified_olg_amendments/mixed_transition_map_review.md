Verdict: the six-variable representation is exact on the stipulated uniform branch. It correctly retains positive \(\chi\), positive \(\tau^p\), logistic tenure choice, and heterogeneous \(F\). I find no additional endogenous aggregate state after date \(1\). The important qualification is that \((B_0,A_0,H_0)\) suffices for date-0 *market clearing*, but not to reconstruct every initial old household’s consumption and estate; the original initial distribution \(G_0\) must remain given data.

Let \(a=h_R^{\max}\) and \(A=1+\gamma+\omega_B\). On the old-owner interior branch, the dated old FOCs and the young saving FOC imply
\[
c^2_{t+1,i}=\frac{\beta x^O_{t,i}}{q},
\qquad
h^2_{t+1,i}=\frac{\gamma c^2_{t+1,i}}{u_{t+1}}
=\frac{\beta\gamma}{q\,u_{t+1}}x^O_{t,i}.
\]
On the old-renter cap branch, \(h^{2,R}_{t+1,i}=a\). Hence, exactly,
\[
\bar H_t^O=\frac{M_t}{u_t}+R_t,
\]
with the proposed
\[
M_t=\frac{\beta\gamma}{q}Y_{t-1}\!\int\!\pi^O_{t-1}(i)x^O_{t-1}(i)\,dF,
\quad
R_t=aY_{t-1}\!\int\![1-\pi^O_{t-1}(i)]\,dF.
\]
Thus \(M_t\) is not financial wealth: it is the sufficient aggregate coefficient on \(1/u_t\) in old-owner housing demand.

Given \((P_t,u_t,Y_t,O_t,M_t,R_t)\),
\[
P_{t+1}=\frac{(1+q\tau^p)P_t-u_t}{q},\qquad
T_t=\frac{q\tau^pP_t\bar H}{Y_t+O_t}.
\]
For candidate \((u_{t+1},Y_{t+1})\),
\[
T_{t+1}=\frac{q\tau^pP_{t+1}\bar H}{Y_{t+1}+Y_t}.
\]
Evaluate the original conditional owner and renter problems at
\((P_t,u_t,u_{t+1},T_t,T_{t+1})\), and use their logistic probability \(\pi_i\). The two equations are precisely
\[
0=Y_t\!\int [\pi_i h_i^O+(1-\pi_i)a]\,dF+\frac{M_t}{u_t}+R_t-\bar H,
\]
\[
Y_{t+1}=\nu Y_t\!\int[\pi_i n_i^O+(1-\pi_i)n_i^R]\,dF.
\]
The stated updates for \(O,M,R\) then give the exact map.

The tenure lead is handled correctly. Conditional policies do not depend on the ownership taste draw because it enters additively; its only effect is the logistic \(\pi_i\). That same \(\pi_i\) enters current young demand, fertility, and the next old-owner/renter aggregates. No lagged ownership-share state is missing.

The useful sign check is also correct. Holding \(T_{t+1}\) fixed,
\[
\frac{\partial(W_i^O-W_i^R)}{\partial u_{t+1}}
=-\frac{\beta\gamma}{u_{t+1}}+\frac{qa}{x_i^R}.
\]
The old renter’s cap strictly binding means its unconstrained demand exceeds \(a\):
\[
\frac{\beta\gamma x_i^R}{q u_{t+1}}>a.
\]
This is exactly the condition making the displayed derivative negative. At \(\tau^p=0\), current conditional housing levels do not themselves depend on \(u_{t+1}\); with \(h_i^O>a\), logistic ownership falls strictly in \(u_{t+1}\), so young housing demand falls strictly. This supplies the claimed local inversion without an all-owner or small-\(\chi\) argument.

The date-0 correction is essential and correctly stated. For actual initial owners,
\[
h^2_{0,i}=\frac{\gamma}{A\,u_0}
(a_{0,i}+P_0H_{0,i}+T_0).
\]
Therefore aggregate initial owner housing is \(M_0/u_0\), where
\[
M_0=\frac{\gamma}{A}\{A_0+P_0H_0+T(P_0,Y_0,O_0)B_0\},
\]
and \(R_0=a(O_0-B_0)\). Thus \(M_0\) must vary with the surprise \(P_0\) and rebate. Treating it as fixed would revalue neither title nor rebate correctly.

One qualification: \(A_0,H_0,B_0\) do not by themselves reconstruct each initial old household’s \(c_0^2,h_0^2,e_0\), nor initial renters’ consumption/estate. Keep \(G_0\), or the type-level list/distribution of \((a_{0,i},H_{0,i},m_i)\), as exogenous initial data. This is not an additional *dynamic aggregate state*: under the cap/slack branch it does not affect prices beyond the four aggregates, but it is required for a full original-equation equilibrium and to verify individual slackness. These claims must also be feasible pre-reform claims and remain uniformly slack after the small surprise.

Precise regularity conditions:

- \(F\) is fixed; a compact support plus uniform branch margins is sufficient, or otherwise dominated first derivatives permitting differentiation under the integral.
- \(\sigma_\xi>0\), \(P,u,Y,O,Y+O>0\), and adult consumption, adult space, fertility, and estate are uniformly positive.
- Purchase limits bind for every young owner; owner physical caps are slack; young and old rental caps bind strictly; saving is positive in both tenures; old owner retention and estate-composition constraints are slack uniformly.
- The autonomous steady-state claim also requires time-invariant \(\vartheta_t\) and all other primitives after the permanent reform.
- The two-equation Jacobian in \((u_{t+1},Y_{t+1})\) is nonsingular in a neighborhood, not only at the fixed point.

Finally, the proposed local theorem is sufficient if stated as a sequence-IFT/exponential-dichotomy result. Let the four boundary restrictions have rank four and let their derivative restricted to the four-dimensional stable generalized eigenspace be nonsingular. With four roots strictly inside the unit circle and two strictly outside, this is exactly the required transversality. Zero stable roots cause no obstruction: they are strictly stable. One should avoid relying solely on the diffeomorphism version of the stable-manifold theorem when the map is noninvertible; the one-sided sequence IFT works because the unstable block is invertible and the stable block is propagated forward. Under these conditions, the boundary intersects the local stable set uniquely and yields exponential convergence for sufficiently small permanent reforms.