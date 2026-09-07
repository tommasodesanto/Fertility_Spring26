# Pro response: analytical equilibrium inefficiency

Source: https://chatgpt.com/c/6a9f2beb-c8b0-83e9-b649-3304f9801ad3

Retrieved 7 September 2026. External model response, preserved for independent review; this is not an adopted manuscript result. Equations extracted from the rendered math labels. Attachment citation buttons omitted.

## Verdict

**The unrestricted claim is false.** The model admits stationary equilibria that are efficient for the packet’s fixed-fertility, fixed-tenure, datewise welfare comparison.

However, a stronger result is available than existence of one inefficient equilibrium: **under explicit primitive inequalities, every stationary equilibrium is inefficient.** The result below retains endogenous fertility, positive child goods and space costs, finite housing caps, and finite-scale logistic tenure choice. It uses homogeneous income and wealth and zero property tax. It also covers old owners whose estate constraint binds—not only those leaving positive financial assets.

---

## 1. Setup and an analytical efficient class

Set $\tau^p=0$ , let all entrants have the same $(y,b)$ , and keep all preference weights and child costs strictly positive. Write

$$
\begin{gathered} w=y+b,\qquad d=1-q,\qquad u=qr=dP,\\ A=1+\gamma+\omega_B,\qquad M=1+\beta A,\qquad D=M+\alpha+\vartheta,\\ \rho=h_R^{\max},\qquad \bar h=h_O^{\max},\qquad B=\frac{b}{1-\phi}. \end{gathered}
$$

Thus $Ph\le B$ is the down-payment constraint. Importantly, the owner’s gross bond holdings are $qa'+\phi Ph$ , not $qa'$ .

I use **stationarity in levels** : $Y=O=N>0$ and

$$
\nu\bar n=1.
$$

Individual fertility remains a household choice; this is an equilibrium consistency condition. The stationary population $N$ is part of the stationary state, with its scale determined by housing clearing.

### A class of efficient stationary equilibria

Define, entirely from primitives,

$$
u_*=\frac{\nu\vartheta w/D-\chi}{\kappa}, \qquad P_*=\frac{u_*}{d}, \qquad x_*=\frac wD, \qquad n_*=\frac1\nu,
$$

and

$$
\begin{aligned} h_*^Y&=\frac{\alpha x_*}{u_*}+\frac{\kappa}{\nu}, & z_*&=\frac{\beta A x_*}{q},\\ c_*^2&=\frac{\beta x_*}{q}, & h_*^2&=\frac{\beta\gamma x_*}{q u_*}, & e_*&=\frac{\beta\omega_B x_*}{q^2}. \end{aligned}
$$

Suppose the following explicit inequalities hold:

$$
\boxed{ \begin{gathered} u_*>0,\qquad d\omega_B>q\gamma,\\ h_*^2<h_*^Y<\rho<\bar h,\\ P_*h_*^Y<B,\qquad qz_*+(\phi-q)P_*h_*^Y>0. \end{gathered}} \tag{E}
$$

These formulas solve the relaxed household problems, and (E) verifies that the omitted constraints are slack. Both tenures therefore choose the displayed real allocations, with

$$
a_R'=z_*, \qquad a_O'=z_*-P_*h_*^Y.
$$

Their conditional values coincide, so finite logistic tastes produce strictly positive shares of both tenures. Setting

$$
N=\frac{\bar H}{h_*^Y+h_*^2}
$$

completes the stationary equilibrium.

This allocation is **globally efficient within the specified datewise feasible set** , not merely immune to one old-to-young transfer. At fixed fertility and estates,

$$
\frac{\alpha x_*}{h_*^Y-\kappa n_*} = \frac{\gamma c_*^2}{h_*^2} =u_*.
$$

Concavity gives, for each young household,

$$
x_*\Delta U_i^Y\le \Delta c_i+u_*\Delta h_i,
$$

and, for each old household,

$$
c_*^2\Delta U_j^O\le \Delta c_j^2+u_*\Delta h_j^2.
$$

Summing over households makes the right-hand side zero for any consumption- and housing-preserving reallocation. A Pareto improvement on a positive-measure set is therefore impossible.

Thus the **presence** of credit and rental-size limits does not establish inefficiency. The next result makes both young-household limits bind with strictly positive multipliers.

---

## 2. Primitive sufficient conditions for inefficiency of every stationary equilibrium

Define

$$
v=\max\left\{d,\frac{\gamma}{\gamma+\omega_B}\right\}, \qquad \ell=dB,\qquad E=w-\ell,
$$

and

$$
\underline\ell=\frac{\beta\gamma w}{qM+\beta\gamma}.
$$

The central restrictions are

$$
\boxed{ \underline\ell<\ell< \frac{\alpha d\,w}{\alpha d+v(M+\vartheta)} } \tag{A}
$$

and

$$
\boxed{ \frac{\beta A E}{M+\vartheta}>(q-\phi)_+B. } \tag{B}
$$

Condition (B) is automatic when $\phi\ge q$ , but that restriction on $\phi$ is not required.

We next give a demographic restriction using only explicit functions of primitives.

### An explicit fertility function

For $m,E,H>0$ , define

$$
\mathcal N_m(E,H)= \frac{2\vartheta EH} {K_m(E,H)+ \sqrt{K_m(E,H)^2- 4\chi\kappa(m+\alpha+\vartheta)\vartheta EH}},
$$

where

$$
K_m(E,H) =\kappa E(\alpha+\vartheta)+\chi H(m+\vartheta).
$$

It is the unique feasible solution of

$$
mx+\chi n=E, \qquad \frac{\vartheta}{n} =\frac{\chi}{x}+\frac{\alpha\kappa}{H-\kappa n}. \tag{1}
$$

It is strictly increasing in both $E$ and $H$ .

### Primitive price endpoints

Set

$$
t_-=\max\left\{\underline\ell,\frac{\ell\rho}{\bar h}\right\}, \qquad g=\frac{\chi\rho}{\kappa},
$$

and

$$
t_+= \frac{(\alpha+\vartheta)w-Dg+ \sqrt{\big[(\alpha+\vartheta)w-Dg\big]^2+4D\alpha wg}} {2D}.
$$

Define

$$
P_-=\frac{t_-}{d\rho}, \qquad P_+=\frac{t_+}{d\rho},
$$

and the four endpoint fertilities

$$
n_O^\pm=\mathcal N_M\!\left(E,\frac{\ell\rho}{t_\pm}\right), \qquad n_R^\pm=\mathcal N_M(w-t_\pm,\rho).
$$

The final restriction is

$$
\boxed{ \max\{n_O^+,n_R^+\} < \frac1\nu < \min\{n_O^-,n_R^-\}. } \tag{C}
$$

This is an **analytical root-bracketing condition** , rather than a closed-form equilibrium price. Its endpoints and inequalities contain no unknown equilibrium quantities.

### Proposition

Under (A)–(C), with $0<\rho<\bar h<\infty$ , $\bar H>0$ , finite $\bar\xi$ , and $0<\sigma_\xi<\infty$ :

**A stationary state exists. Every stationary equilibrium has $P\in(P_-,P_+)$ , and every such equilibrium is Pareto inefficient for the packet’s datewise comparison.**

In every such equilibrium:

- Young owners have a strictly binding down-payment constraint, a slack physical cap, and positive gross bond holdings.
- Young renters have a strictly binding rental-size cap.
- Old owners sell part of their inherited housing; old renters’ housing caps are slack.
- Both tenures have positive measure.

The old owner’s financial estate is positive when $d\omega_B>q\gamma$ , zero with a strictly binding estate constraint when $d\omega_B<q\gamma$ , and zero with a zero estate-constraint multiplier at equality.

---

## 3. Proof: establish the constraints rather than assume them

### Old owners and the relevant housing value

Let

$$
z=a'+Ph
$$

be the young owner’s total wealth entering old age, including the inherited house. Its lifetime budget becomes

$$
c+qz+uh=w. \tag{2}
$$

The old problem has resources $z$ , housing-retention limit $h^2\le h$ , and estate bound $e\ge Ph^2$ , exactly as obtained from the packet’s budgets.

When retention is slack, solving the old problem gives

$$
c^2=\frac zA, \qquad h^2=\frac{\gamma z}{AvP}, \tag{3}
$$

and

$$
e= \begin{cases} \displaystyle\frac{\omega_B z}{Aq}, &d\omega_B\ge q\gamma,\\[5pt] \displaystyle Ph^2, &d\omega_B<q\gamma. \end{cases} \tag{4}
$$

Consequently,

$$
MV^O=\frac{\gamma c^2}{h^2}=vP. \tag{5}
$$

Thus the donor’s fixed-estate housing value equals user cost only in the financially slack branch. When the estate bound binds strictly, $v>d$ : compensation must cover a higher housing value.

### The owner chooses the largest privately admissible house

Temporarily drop the nonnegative-gross-bond constraint, retaining all other constraints. Write $K=Ph\le B$ .

The fertility first-order condition implies

$$
\chi n<\vartheta x. \tag{6}
$$

Homogeneity of the old problem, together with the saving first-order condition, gives

$$
qz\le\beta A x,
$$

with equality when old housing retention is slack. Combining this with (2),

$$
w-dK=x+\chi n+qz<(M+\vartheta)x.
$$

Therefore

$$
x>\frac{w-dK}{M+\vartheta}. \tag{7}
$$

Since $s<h$ , condition (A) now implies

$$
\begin{aligned} MV^Y =\frac{\alpha x}{s} &> \frac{\alpha(w-dK)P}{(M+\vartheta)K}\\ &\ge \frac{\alpha E P}{(M+\vartheta)B} >vP\ge u. \end{aligned} \tag{8}
$$

Increasing current housing also weakly relaxes the future retention limit. Hence the household’s optimized value is strictly increasing in $h$ up to its current limits:

$$
\boxed{h_O^Y(P)=\min\{\bar h,B/P\}.} \tag{9}
$$

The omitted gross-bond constraint is indeed slack. If retention is slack,

$$
qa'+\phi Ph =qz+(\phi-q)K =\beta A x+(\phi-q)K>0
$$

by (7) and (B). If retention binds, $h^2=h$ , the old budget and $e\ge Ph$ imply

$$
z=c^2+qe+uh>Ph=K,
$$

so $a'=z-K>0$ , which also gives positive gross bonds. Thus dropping that constraint changed nothing.

### Policies inside the bracket

For $P\ge B/\bar h$ , consider

$$
h_O^Y=B/P,\qquad n_O=\mathcal N_M(E,B/P),\qquad x_O=\frac{E-\chi n_O}{M}. \tag{10}
$$

The lower inequality in (A) verifies old retention slack:

$$
\frac{h_O^2}{h_O^Y} =\frac{\beta\gamma x_O}{qvB} \le\frac{\beta\gamma x_O}{q\ell} < \frac{\beta\gamma E}{qM\ell}<1. \tag{11}
$$

Hence (10) is the actual owner solution, with

$$
z_O=\frac{\beta A x_O}{q}, \qquad c_O^2=\frac{\beta x_O}{q}, \qquad h_O^2=\frac{\beta\gamma x_O}{qvP}. \tag{12}
$$

For renters inside the proposed bracket,

$$
h_R^Y=\rho,\qquad n_R=\mathcal N_M(w-dP\rho,\rho),\qquad x_R=\frac{w-dP\rho-\chi n_R}{M}. \tag{13}
$$

Their old housing is

$$
h_R^2=\frac{\beta\gamma x_R}{qdP}<\rho, \tag{14}
$$

because $dP\rho\ge t_-\ge\underline\ell$ .

Finally, without the young rental cap, renter housing demand is

$$
h_{R,\mathrm{free}}^Y(P) =\frac wD\left[ \frac{\alpha}{dP} +\frac{\kappa\vartheta}{\chi+\kappa dP} \right]. \tag{15}
$$

The equation $h_{R,\mathrm{free}}^Y=\rho$ has the explicit solution $dP\rho=t_+$ . Thus the young rental cap binds with positive multiplier whenever $P<P_+$ .

At an interior-bracket price, the owner’s physical cap is slack, its gross bonds are positive, and old retention is slack. Its down-payment multiplier $\mu$ therefore satisfies

$$
MV_O^Y-u=(1-\phi)P\mu x_O.
$$

Equation (8) proves $\mu>0$ . **No assumption that the future estate bound is slack enters this envelope calculation.**

---

## 4. Why the bracket captures every stationary equilibrium

The crucial global fact is:

$$
\boxed{\text{Conditional owner and renter fertility each decrease strictly with }P.}
$$

Their logistic mixture need not be monotone.

For owners above $B/\bar h$ , this follows from (10): $E$ is constant and $B/P$ decreases. Below that price, (9) gives $h=\bar h$ . Depending on the old constraints, fertility is one of

$$
\begin{aligned} &\mathcal N_M(w-dP\bar h,\bar h),\\ &\mathcal N_{M_0}(w-(1+q)dP\bar h,\bar h), \qquad M_0=1+\beta(1+\omega_B),\\ &\mathcal N_{1+\beta}(w-P\bar h,\bar h). \end{aligned} \tag{16}
$$

These correspond respectively to slack retention; binding retention with slack estate bound; and both old constraints binding. Each decreases strictly on its applicable range. The branches join continuously. The relevant wealth-to-house-value ratios decrease along each branch, so the old regimes cannot switch back and forth.

For renters, (A) implies $\alpha>\beta\gamma/q$ . An old rental cap therefore cannot bind while the young rental cap is slack. The successive fertility schedules are

$$
\mathcal N_{M_0}(w-(1+q)dP\rho,\rho),\qquad \mathcal N_M(w-dP\rho,\rho),\qquad \frac{\vartheta w}{D(\chi+\kappa dP)}. \tag{17}
$$

Again, each decreases strictly and the branches join continuously.

Now observe

$$
t_-<\ell<\frac{\alpha w}{D}<t_+<w.
$$

At the intermediate price

$$
P_0=\frac B\rho,
$$

both tenures choose the same current housing, consumption, and fertility:

$$
n_O(P_0)=n_R(P_0)=n_0:=\mathcal N_M(E,\rho).
$$

Consequently,

$$
\max\{n_O^+,n_R^+\} <n_0< \min\{n_O^-,n_R^-\}. \tag{18}
$$

This also proves that the demographic interval in (C) is nonempty.

Let $\pi(P)$ be the logistic owner share and set

$$
\bar n(P)=\pi(P)n_O(P)+(1-\pi(P))n_R(P).
$$

It is continuous. Condition (C) gives

$$
\bar n(P_-)>\nu^{-1}>\bar n(P_+),
$$

so a stationary fertility root exists. Global conditional monotonicity also gives

$$
P\le P_-\implies\bar n(P)>\nu^{-1}, \qquad P\ge P_+\implies\bar n(P)<\nu^{-1}.
$$

**Every stationary root is therefore inside the bracket.**

For any root, let $\pi=\pi(P)$ and set

$$
N= \frac{\bar H} {\pi(h_O^Y+h_O^2)+(1-\pi)(\rho+h_R^2)}. \tag{19}
$$

Together with the induced asset distribution, this satisfies housing clearing and stationary demography. Finite logistic tastes give $0<\pi<1$ .

---

## 5. The compensated Pareto improvement

At every covered equilibrium,

$$
MV_O^Y>vP=MV_O^O. \tag{20}
$$

Select a common positive submass of young and old owners. This is a proof device, not an assumption about matching or cohort shares.

Give each selected young owner $\epsilon$ additional housing and take it from a selected old owner. Keep fertility, tenure, and the old estate fixed. The donor’s exact consumption compensation is

$$
D(\epsilon) =c^2\left[ \left(\frac{h^2}{h^2-\epsilon}\right)^\gamma-1 \right].
$$

It preserves the donor’s utility exactly, and

$$
D'(0)=\frac{\gamma c^2}{h^2}=vP.
$$

The young owner’s gain is

$$
\Delta U^Y = \log\!\left(1-\frac{D(\epsilon)}x\right) +\alpha\log\!\left(1+\frac{\epsilon}{s}\right),
$$

so

$$
\left.\frac{d\Delta U^Y}{d\epsilon}\right|_{\epsilon=0} =\frac{MV_O^Y-vP}{x}>0.
$$

A sufficiently small positive transfer therefore benefits the young strictly while respecting their physical cap and keeping all consumption positive.

The packet’s settlement applies unchanged. Set

$$
L(\epsilon)=D(\epsilon)-u\epsilon.
$$

The old owner receives $P\epsilon+L(\epsilon)$ , consumes $D(\epsilon)$ , and buys a bond costing $qP\epsilon$ , replacing the title removed from its estate. The young owner’s new borrowing proceeds are also $qP\epsilon$ ; next period it sells the additional title for $P\epsilon$ and repays that new debt. Existing estates, future real allocations, and promised repayments on existing loans remain unchanged.

This also works when the old owner initially has **zero** financial saving: the settlement creates positive financial saving to replace the transferred estate housing.

---

## 6. Economic meaning and scope

The lower bound on $\ell$ prevents the inherited owner house from being so small that the old household wants to retain all of it. The upper bound makes liquid wealth sufficiently scarce relative to income and young housing preferences that young owners value additional housing more than eligible old donors do. Condition (B) ensures positive **gross** saving; it does not require positive net $a'$ .

The factor $v$ is important. When the estate constraint binds, old housing is worth more than ordinary user cost at fixed estate. The upper liquidity bound automatically becomes tighter rather than incorrectly treating such donors as having marginal value $u$ .

The region is demonstrably nonempty. For arbitrary positive old-age preference weights, choose

$$
\alpha> \frac{\beta\gamma v(M+\vartheta)}{dqM}.
$$

Then the interval in (A) is nonempty. Choosing $\phi\ge q$ makes (B) automatic; for any admissible $\ell$ and $w>0$ , set

$$
b=\frac{(1-\phi)\ell}{d}, \qquad y=w-b>0.
$$

Positive $\chi,\kappa$ and any finite ordered caps are admissible, and (18) supplies a nonempty interval for $\nu$ .

For **existence alone** , (C) can be weakened to the taste-dependent endpoint condition

$$
\bar n(P_+)<\nu^{-1}<\bar n(P_-),
$$

using the explicit endpoint policies in the logistic rule. That weaker condition guarantees an inefficient stationary state but does not exclude other stationary equilibria.

The proved result is for homogeneous entrants and zero property tax. It does not establish price uniqueness, convergence from arbitrary initial states, or inefficiency outside the displayed region. If the stationary cohort mass $N$ is prescribed independently, housing-clearing compatibility in (19) must additionally hold. No unresolved step remains inside the proposition; heterogeneous entrants and positive, price-dependent tax rebates require a separate equilibrium argument.
