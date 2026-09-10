Continue our theory discussion. Independently verify the attached all-tenure simplification and resolve its precisely stated fertility distinction. Return one short useful allocation proposition with its proof.

### File: docs/prompts/oracle_all_tenure_simplification.md
```md
# A shorter allocation theorem, allowing every young tenure choice

Continue this existing discussion. The previous common-retirement-income theorem passed our independent mathematical checks, but the author found its many restrictions unsuitable for a simple illustrative exercise. His presentation is in four days. We have now derived a more general argument. Please verify it independently and improve its economic statement. Do not reopen the lifecycle design or return another broad menu. We need one short, honest proposition and an indispensable proof, not a declaration that the paper is solved.

The main comparison remains the FULL equally weighted dated planner choosing consumption AND housing for all currently living households. First fix competitive fertility. Then distinguish (a) a compensated local reallocation followed by private fertility choice, (b) private fertility after receiving the fixed-fertility full-planner bundles, and (c) the joint parents-only planner choosing fertility too. These are different experiments. No independent welfare weight is attached to unborn people. No transition or tax-policy appendix is requested in this round.

## Freeze this specification

There are two periods. Young resources have a positive heterogeneous distribution w (two types wL<wH with probabilities f and 1-f suffice); everyone receives the same positive retirement income y. All resources already received are usable. Keep

u^y = log(c-chi*n) + alpha*log(h-kappa*n) + vartheta*log(n),
u^o = log(c^o) + alpha*log(h^o) + omega_B*log(e),
U = u^y + beta*V^o(z).

All preference and child-cost coefficients are positive. The serviced-interest owner rule is

c+(1+tau)*P*h+a=w, z=y+a/q+P*h, a>=-q*phi*P*h.

Renters satisfy c+p*h+a=w, a>=0, h<=r. There is a common divisible housing stock, no owner size minimum, upper bound, taste term or adjustment cost. Take 0<q,phi<1 and tau>0. Define

d=1+tau-q, p=d*P, Gamma=q*(1-phi)/d>0,
K=1+alpha+omega_B, J=1+alpha+vartheta.

The owner consolidated constraint is c+p*h+q*z=w+q*y, q*(z-y)>=Gamma*p*h. Old households use the same rule with current resources z, subsequent estate e and no further income. Under omega_B>=alpha*Gamma, unrestricted old choices are feasible at every z>0:

c_o=z/K, h_o=alpha*z/(K*p), q*e=omega_B*z/K,
V^o(z)=K*log(z)+constant.

The equality case may have a binding constraint with zero multiplier. External finance and estates independent of entrants' endowments remain the proposed closure. Positive property tax finances public services, fixed in the dated comparison, rather than rebates. Do not claim these proposals have already been adopted into the manuscript.

The dated planner preserves individual young continuation resources z, individual old estates e, and continuation prices. It preserves current goods and housing totals. It can reassign tenure and relax private finance; 100%-balance finance suffices because young z>=y. Transfers satisfy T=Delta c+p*Delta h, sum T=0; owned-title adjustments Delta a=-q*P*Delta h_owned preserve continuation targets, with matching intermediary positions. This compares one allocation at one date, including the initial old. It is not a stationary lifetime-welfare comparison.

For stationary equilibrium retain your previous demographic closure: nu*n entrants per parent, fixed type mix independent of parent type, equal endogenous young/old mass N at replacement. At a stationary equilibrium all inherited old resources are >=y. Away from stationarity, old inherited titles must be revalued at current P; fixed old purchasing power is not an interchangeable equilibrium closure.

## 1. Exact affected-income interval, without preselecting tenure

Let delta=beta*K, Y=q*y,
A(p)=alpha+vartheta*kappa*p/(chi+kappa*p), xr=p*r/A(p).

The uncapped rental optimum first meets the cap at

wr(p)=J*xr+max{delta*xr-Y,0}.

The unrestricted lifetime plan first becomes financeable through ownership at

wO(p)=Y*[J+Gamma*A(p)]/[delta-Gamma*A(p)] if delta>Gamma*A(p),
wO(p)=infinity otherwise.

Claim: EVERY globally optimal tenure choice has housing MRS>=p; strict inequality holds precisely for wr(p)<w<wO(p). The interval is nonempty iff delta*p*r<A(p)*[Y+Gamma*p*r]. When the interval is empty there is no strict distortion. A tenure tie within the interval has strictly distorted choices on both branches.

Reason: ownership is a subset of the hypothetical uncapped-rental menu, which still imposes z>=y. Below wr that menu's optimum is actually feasible as a renter. Above wO the fully unrestricted lifetime optimum is ownership-feasible. Between those thresholds both candidate optima have strictly positive housing wedges. This avoids solving the rent-own value crossing. We separately found a single crossing within this interval, but it is unnecessary to the main theorem. Household choices with constrained ownership need a scalar cubic; do not falsely claim universal elementary closed forms.

A simpler sufficient test for a young household is
w<J*q*y/(beta*K) and
(w/J)*[alpha/p+vartheta*kappa/(chi+kappa*p)]>r.
It wishes to bring retirement resources forward and wants more than the largest rental. It may optimally rent OR own. Very poor households can be undistorted because a rental below r suffices.

## 2. A spending identity replaces the strong consumption-gain condition

For each competitive young choice define adult goods x_i=c_i-chi*n_i, adult-space value v_i=p*(h_i-kappa*n_i)/alpha, and t_i=x_i/v_i>=1. These are proof objects, not additional primitives. Use a different symbol from estate e for the following index:

xhat_i=(c_i+p*h_i)/J=[w_i-q*(z_i-y)]/J.

Let C=chi+kappa*p, lambda=chi/(kappa*p), and g(x,v)=(chi/x+kappa*p/v)^(-1). Private fertility satisfies n_i=vartheta*g(x_i,v_i). Direct algebra gives

xhat_i-v_i = v_i*(t_i-1)/J * [1+vartheta*lambda/(lambda+t_i)],
x_i-xhat_i = v_i*(t_i-1)/J * [alpha+vartheta*t_i/(lambda+t_i)],
xhat_i/C-g(x_i,v_i) = v_i*(t_i-1)*(t_i-alpha*lambda)/[J*C*(lambda+t_i)],
C*g(x_i,v_i)-v_i = lambda*v_i*(t_i-1)/(lambda+t_i).

Therefore, under kappa*p>=alpha*chi,

v_i <= C*n_i/vartheta <= xhat_i <= x_i,

with strict relevant inequalities for t_i>1, including at equality in the cost bound. Economically, children's housing-to-goods spending ratio kappa*p/chi is at least the unrestricted adult ratio alpha. This matches the author's intuition that children's needs may be more housing intensive. It is not the assumption that the young must receive both more goods and more housing.

Write B=mean old consumption, mu=N_o/N_y>0. At fixed competitive fertility the full planner has common adult goods and space

X=(mean x+mu*B)/(1+mu),
V=(mean v+mu*B)/(1+mu), S=alpha*V/p.

Young i gets (X+chi*n_i,S+kappa*n_i); old gets (X,S). The sharper resource condition is B>=mean xhat. Positive mass in the affected-income interval then implies B>mean v, and hence

mean(h_y^F)-mean(h_y^E)=alpha*mu*(B-mean v)/[p*(1+mu)]>0.

This part does not require the child-cost ratio restriction. At B=mean xhat, mean young consumption actually falls because mean x>mean xhat.

A simple fully primitive sufficient income restriction is

y/K >= mean(w)/J.

Indeed xhat_i<=w_i/J, and B>=y/K at a stationary equilibrium. This avoids any beta/q ordering and works across all young tenure regimes. It is conservative: it ignores old accumulated wealth and young saving. It is NOT uniformly weaker than your former regime-specific income bound; do not claim nested parameter regions or empirical mildness.

## 3. Joint parents-only fertility

Normalize current total resources per young household as Ctotal and Htotal. After optimizing c,h, the joint planner's reduced objective is

(1+mu)*log[X(n)]+(1+mu)*alpha*log[S(n)]+vartheta*log(n)+constants,
X(n)=(Ctotal-chi*n)/(1+mu), S(n)=(Htotal-kappa*n)/(1+mu).

Its derivative D(n)=vartheta/n-chi/X(n)-alpha*kappa/S(n) is strictly decreasing. The function g is concave:

d2g=-2*chi*kappa*p*(v*dx-x*dv)^2/(chi*v+kappa*p*x)^3.

Under B>=mean xhat, kappa*p>=alpha*chi and positive housing distortion, B/C>mean n/vartheta. Jensen gives

g(X,V)>mean n/vartheta,

where X,V are evaluated at competitive mean fertility. Thus D(mean n)>0 and the JOINT planner chooses nF>mean n. Young aggregate housing rises further when joint fertility rises. Two independent analytical cross-reviews passed the combined proof, including arbitrary mu and weak cost inequality, but please examine it independently.

## 4. Private fertility is a separate question

For a compensated small housing transfer to any distorted young household, old exact compensation costs p*epsilon+O(epsilon^2) goods. The recipient's private fertility increases locally iff kappa*p*t_i^2>alpha*chi. Therefore the weak cost bound above suffices for a strictly distorted recipient. This is a local compensated reallocation, not the full utilitarian optimum.

After receiving the fixed-n full-planner bundle, parent i increases private fertility iff n_i < vartheta/(chi/X+alpha*kappa/S). The joint proof puts this cutoff above competitive MEAN fertility; it does not by itself establish a rise in mean PRIVATE fertility after all parents rechoose at their assigned total bundles. Please settle whether the stronger primitive floor y/K>=mean(w)/J proves that mean-private result or find an analytical counterexample. Also establish the individual low-income scope if a short result is possible without restoring the forced-tenure regime.

Our independent check has just produced a counterfamily even under the STRONGER primitive floor. Please verify this rather than chase a mean-private theorem if it is false. Set alpha=vartheta=chi=kappa=1; p=P=2; q=1/2, beta=1/8, tau=1/2, phi=1/2 (Gamma=1/4); omega_B=2, K=4; y=8; equal young incomes 3 and 9; r=5/2-epsilon for 0<epsilon<=1/10000. Both types rent and save zero; the richer young hit the rental cap. Old resources are all 8, so B=2=mean xhat=y/K=mean(w)/J. At epsilon=0, competitive fertility is 1/3 and 1. The fixed-n planner assigns common adult goods X=2 and adult space S=1. Private rechoice gives fertility (11-sqrt(37))/9 and (5-sqrt(7))/3, whose mean is below 2/3 by (sqrt(37)+3*sqrt(7)-14)/18>1/1000. The agent derives an explicit rechoice perturbation bound 19*epsilon/8, so the entire positive-epsilon family retains the mean loss. Owner deviations are excluded globally by your concavity-support inequality xi<Gamma*p*zeta. Per-unit-cohort housing clearing uses Hstock=8/3-epsilon/2; choose nu=1/competitive mean n for stationary replacement. This analytical counterexample is meant to DISPROVE a universal claim, not to prove model existence through a numerical reference point. The full joint planner fertility still rises. The poorer individual gains, the cap-distorted richer individual loses under the full utilitarian redistribution. That is a useful warning about claiming every constrained recipient gains. Verify the bounds or give a simpler exact argument, and make the surviving fertility statement precise.

## 5. Analytical nonvacuity and the requested finished answer

We verified compatibility with an ACTUAL stationary equilibrium using your last response's analytical construction. In your notation choose omega_B>alpha*Gamma, beta*K>Gamma*(alpha+vartheta); wL in your prior strict lower-income interval so b<y/K and wL<J*y/K; and wH above your finance threshold and J*y/K. Choose f sufficiently close to one to satisfy both your old resource bound and f*wL+(1-f)*wH<J*y/K. Impose the explicit positive bound chi<kappa*wL/[r*(1+alpha+k)], where k=b/(wL/J). Then at p0=alpha*chi/kappa, M(p0)>k*p0, so p0 lies below your lower price cutoff. Every price in the old verified regime satisfies the child-cost condition. Finally choose replacement 1/nu strictly between F at its two cutoffs, and set N=Hstock/D(p). All inherited old wealth is generated by the previous young choices. This is an analytical open family, not a numerical-reference-point proof.

Qualification: this compatibility witness retains the old stronger resource bound and does not exhibit young consumption losses. If it is simple to exhibit an analytical actual-equilibrium region where the new result works with consumption losses, do so. Do not add its construction inequalities to the main proposition's list merely because they prove nonvacuity.

Please deliver: (1) your verdict on the combined result; (2) ONE concise proposition with visible assumptions and the appropriate aggregate/individual scope; (3) a transparent proof and short economic assessment of the resource and child-intensity bounds; (4) a clear verdict on private versus joint fertility. If a materially better simplification exists within this frozen model, develop it. Otherwise tell us which three or four conditions carry the economics and put pure regime/existence algebra after the main argument. Keep the finished answer to a few readable pages. Avoid new timing conventions, unrelated literature digressions, a transition appendix, or claims stronger than the proof.
```
