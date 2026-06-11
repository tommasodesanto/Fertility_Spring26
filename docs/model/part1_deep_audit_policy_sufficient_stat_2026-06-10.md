# Part I Deep Audit And Policy Sufficient Statistic

Date: 2026-06-10

Scope: standalone Part I analytical model in `latex/intergenerational_housing_fertility_part1.tex`. This note is a working audit, not a replacement draft.

## 1. Deep Audit

The current Part I is on the right path. It now starts from an economy: young potential parents, old incumbent owners, rental and owner housing segments, construction, aggregate market clearing, a planner, and then marginal comparisons. That is the right macro-theory order. It no longer reads primarily as a list of wedges.

The strongest part is the young household block. Tenure is now primitive, and the child-input structure is economically coherent: children require adult residual consumption and adult residual housing services through
\[
c_i=C_i-\chi_i n_i,\qquad s_i=h_i-\kappa n_i.
\]
The branch representation with residual consumption \(c_i\) is the right way to write the FOCs. The text must keep reminding the reader that \(c_i\) is adult residual consumption, not total nonhousing consumption.

The Section 9 derivative is now doing useful work. Under equal scaling,
\[
\chi_i=\frac{\kappa q^m}{\alpha},
\]
one child uses consumption and housing in the same expenditure proportions as the adult residual Cobb-Douglas bundle. Then any relaxation of a strictly binding effective family-housing cap raises fertility until the cap stops binding. This is the cleanest fertility statement in the current draft.

The main remaining economic issue is the distinction between three objects that are still too easy to blur:

1. Real size menus: \(h\le \bar h^R\) and \(h\le \bar h^O\). These are real constraints unless the planner has construction, conversion, or assignment instruments that change the menu.
2. Financial implementation constraints: down-payment and payment-to-income constraints. These can create a CE/planner gap if the planner has feasible instruments that relax them at resource cost.
3. Aggregate stock scarcity: \(H^m=\bar H^m+I^m\), cleared by prices and construction.

The current text now handles this better, but the paper should keep these objects separate in every efficiency and policy statement. A binding rental size cap is not by itself an inefficiency if the planner faces the same cap. It is a shadow value of relaxing a real constraint. A binding down-payment constraint is a private implementation gap only if the planner can relax it through an explicit instrument.

There is also a subtle construction/menu issue. Current construction \(I^m\) expands aggregate segment services and affects \(q^m\). It does not, by itself, raise the individual size cap \(\bar h^m\). If the policy is "more family-sized rental units" or "conversion of owner units into family-capable rentals," the model needs an instrument that changes the effective cap or menu, not merely an increase in aggregate stock. Otherwise the policy changes prices but not the family-size access constraint that drives the clean Section 9 result.

The planner section is improved but still should be sold cautiously. It is a conditional marginal comparison, not a full welfare theorem. It is conditional on entry, tenure assignments, interior construction, and being away from real menu corners. That is acceptable, but the prose should keep saying "conditional marginal comparison" rather than "the CE is inefficient" without qualifications.

The old-side property-tax mechanism is conceptually clean. \(L_j^\tau\) is not old utility and not a real resource saving. It lowers the old owner's private carrying cost below the social opportunity cost. The old-side channel can amplify young access problems, but the young-access/fertility mechanism does not depend on property tax.

Notation cleanup still matters before this becomes paper-facing. Superscript \(O\) means young owner tenure and old household variables. That is readable locally but fragile. A later pass should use \(B\), \(\mathrm{own}\), or \(b\) for young buyers and reserve "old" notation for old households. Also, \(A_i(x)\) should be described as pledgeable liquid resources for purchase constraints, not full wealth, unless it is put into the lifetime budget.

## 2. Sufficient Statistic For Policy Effectiveness

The simplest sufficient-statistic object is local and branch-specific. It should not pretend to solve the full equilibrium. Fix a policy parameter \(z\). Fix a type \(i\), tenure branch \(m\), and active set. Let the policy change the household's effective family-housing cap
\[
H_i^m(z).
\]
For renters, \(H_i^R=\bar h^R\). For owners,
\[
H_i^O=
\min\left\{
\bar h^O,
\frac{A_i(x)\rho_p}{(1-\phi)q^O},
\frac{\psi y_i}{q^O}
\right\}.
\]

For a type whose cap binds, define the branch fertility response to a cap relaxation:
\[
\Psi_i^m
=
\frac{
\frac{\alpha\kappa}{s_i^2}
-
\frac{\chi_i q^m}{c_i^2}
}{
\frac{\beta}{n_i^2}
+
\frac{\chi_i^2}{c_i^2}
+
\frac{\alpha\kappa^2}{s_i^2}
}.
\]
This is \(\partial n_i^m/\partial H_i^m\) holding the branch, prices, income, child costs, and preferences fixed. Under equal scaling and a strictly binding cap, \(\Psi_i^m>0\).

Define policy pass-through to the effective cap:
\[
P_i^m(z)=\frac{\partial H_i^m}{\partial z}.
\]
Examples:

- Rental cap relaxation: if \(z=\bar h^R\), then \(P_i^R=1\) for capped renters.
- Down-payment liquidity grant: if the down-payment component uniquely binds and \(z=A_i\), then
\[
P_i^O=\frac{\rho_p}{(1-\phi)q^O}.
\]
- Payment-to-income relaxation: if \(z=\psi\), then
\[
P_i^O=\frac{y_i}{q^O}.
\]
- A policy that makes the wrong units cheaper has \(P_i^m=0\) for fertility if it does not relax the binding family-space cap.

With fixed entry, tenure, and prices, the local aggregate fertility effect is
\[
\frac{\dd N^{cond}}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\pi_i^E\,
\Psi_i^m\,
P_i^m(z)
\,
\dd G_Y(i),
\]
where \(\mathcal B(z)\) is the set of entered types in the active branch whose effective family-housing cap binds.

This is the clean fertility sufficient statistic:
\[
\text{fertility effect}
=
\text{mass exposed}
\times
\text{policy pass-through to effective family space}
\times
\text{fertility response to family space}.
\]

If entry is allowed to respond but tenure and prices are still held fixed, the logit entry block adds an entry term. Since
\[
\frac{\partial \pi_i^E}{\partial W_i^Y}
=
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\]
and the utility value of relaxing the cap is \(\lambda_i^m\zeta_i^m\), the local effect becomes
\[
\frac{\dd N}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\left[
\pi_i^E\Psi_i^m
+
n_i^m
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\lambda_i^m\zeta_i^m
\right]
P_i^m(z)
\,
\dd G_Y(i).
\]
The first term is the intensive fertility effect among entrants. The second is the entry effect: policies that relax binding family-space constraints raise inside value and can draw in additional households.

This statistic is sufficient for a local, fixed-price, fixed-branch policy experiment. It is not sufficient for the full general equilibrium response. A full GE policy effect also needs price pass-through, construction responses, tenure switching, changes in old retention, and the resource cost of instruments.

For welfare, the analogous local statistic is even simpler:
\[
\frac{1}{\Lambda}\frac{\dd W}{\dd z}
=
\bar M
\int
\pi_i^E
\zeta_i^{relaxed}
P_i^m(z)
\dd G_Y(i)
-
\frac{\dd \Xi^\Pi}{\dd z}
\]
plus any old-side term if the policy changes old retained housing:
\[
\int L_j^\tau
\left(-\frac{\partial h_j^O}{\partial z}\right)
\dd F_O(j).
\]
Thus the birth statistic uses \(\Psi_i^m\), while the welfare statistic uses \(\zeta_i\). They are related but not the same object.

## 3. Simple Paper Statement

The clean statement is:

> A policy affects fertility only to the extent that it relaxes the effective family-housing cap of households for whom that cap binds. Under equal scaling of child consumption and housing requirements, the fertility effect of a local cap relaxation is positive for every strictly constrained household. Aggregate effectiveness is the mass of constrained households times the policy pass-through to their effective cap times their branch fertility response, plus any entry response.

That sentence avoids two common overclaims. It does not say every housing policy raises fertility. It says the policy must relax the binding family-space constraint. It also does not say the planner always raises fertility in GE. It says the local fertility response is signed when the policy relaxes the effective cap under the stated primitive condition.

## 4. Two Figures To Use

Figure 1 should show the branch mechanism. Put the effective cap \(H\) on the horizontal axis and fertility \(n(H)\) on the vertical axis. Under equal scaling, the curve rises over the binding region \(H<h^u\) and flattens once the cap no longer binds. Label the shaded binding region "family-space constrained" and annotate the policy arrow \(\dd H>0\), \(\dd n>0\). This picture should be the visual anchor for Section 9.

Figure 2 should show the sufficient-statistic decomposition. The picture should not be a supply-demand graph. It should be a pipeline:
\[
\text{policy}
\to
\text{pass-through to } H_i^m
\to
\text{fertility response } \Psi_i^m
\to
\text{aggregate births}
\]
with a separate entry branch. The bottom of the figure should display the formula
\[
\dd N
\approx
\bar M\int
\left[
\pi_i^E\Psi_i^m
+
n_i^m
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\lambda_i^m\zeta_i^m
\right]
\dd H_i^m
\dd G_Y(i).
\]
The visual point is that the policy works only through \(\dd H_i^m\), the relaxation of the effective family-space cap. If \(\dd H_i^m=0\), the policy may affect ownership labels or prices but not this fertility mechanism.

