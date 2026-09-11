## Utility specification: the issue, discussion, and agreed solution

We reviewed how children enter household utility because the existing specification combined an equivalence scale with child-dependent minimum housing requirements. The objective was to understand the economic restrictions and choose a simpler specification to carry into the pension-corrected recalibration.

### 1. The previous specification

Let \(c\) denote nonhousing consumption, \(s\) housing services, and \(m\) the number of currently dependent children.

The consumption–housing composite was

\[
Q(c,s;m)
=
c^{\alpha_0}
\left[s-\bar h(m)\right]^{1-\alpha_0},
\]

with housing requirement

\[
\bar h(m)=h_1\mathbf 1\{m>0\}+h_m m.
\]

Thus, becoming a parent increased required housing by \(h_1+h_m\), and each subsequent dependent child added another \(h_m\).

Flow utility was

\[
u_t(c,s;m)
=
\frac{[Q(c,s;m)/e(m)]^{1-\sigma}}{1-\sigma}
+\psi_t m,
\qquad \sigma=2,
\]

where

\[
e(m)
=
\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]

The scale divides the entire consumption–housing composite. There is no nonhousing consumption floor and no additional household multiplier outside CRRA.

### 2. What initially caused concern

Children affected material utility through two channels:

- The equivalence scale increased the bundle needed to maintain a given effective living standard.
- The housing minimum increased both at entry into parenthood and with every additional child.

These channels are mathematically coherent. We did not establish that their coexistence necessarily constitutes double counting.

However, borrowing an equivalence scale does not independently validate the total expenditure requirement after adding a separate housing minimum. The scale applies to the bundle remaining after required housing; its original empirical interpretation cannot automatically be transferred to that combined specification.

For an interior renter, expenditure required to attain an effective material living standard \(q=Q/e(m)\) can be written as

\[
E(m,q,r)
=
r\bar h(m)+P_Q(r)e(m)q,
\]

where \(r\) is the rental price of housing services and \(P_Q(r)\) is the price of the discretionary consumption–housing composite.

This makes the two components explicit: an additive housing requirement and a scaled discretionary bundle.

### 3. Calibration partly addresses cost levels

Tommaso correctly emphasized that calibration can partly accommodate the overall cost of children.

If the imposed scale raises material costs, stronger child preferences or smaller housing requirements can compensate, subject to fitting the other moments. Therefore, an apparently large combined cost is not, by itself, proof of misspecification.

The more substantive question concerns the **functional form**: how costs vary across family sizes and household resources.

A linear child benefit, \(\psi_t m\), can change the overall attractiveness of children. It cannot independently change curvature across family sizes. Housing parameters can affect that curvature, but their adjustment is constrained by housing evidence.

### 4. What the scale’s concavity means

The scale is increasing and concave:

| Dependent children | Scale | Additional increment |
|---|---:|---:|
| 0 | 1.0000 | — |
| 1 | 1.2338 | 0.2338 |
| 2 | 1.4498 | 0.2160 |
| 3 | 1.6528 | 0.2030 |

Successive children require progressively smaller additional resources to maintain the same material living standard. This represents economies of scale from sharing housing and household goods.

**Tommaso considers declining marginal child costs economically reasonable. They are not a defect that needs to be counteracted mechanically.**

This concavity is distinct from CRRA’s diminishing marginal utility of consumption. The former governs economies of scale across children; the latter governs how painful consumption sacrifices are at different resource levels.

In the previous specification, the increasing housing floor could make successive children increasingly costly by reducing the resources left for discretionary consumption. That force opposed the scale’s economies of scale. Its presence nevertheless requires an economic justification beyond producing convenient fertility curvature.

### 5. The agreed simplification

Retain the current equivalence scale, but replace the housing jump plus slope with a single housing requirement for parenthood:

\[
\boxed{
\bar h(m)=\bar h_P\mathbf 1\{m>0\}.
}
\]

The interpretation is:

> Starting a family requires suitable housing. Additional children share that housing and other household resources, with their additional needs captured by the equivalence scale.

The proposed composite is therefore

\[
Q(c,s;m)
=
c^{\alpha_0}
\left[
s-\bar h_P\mathbf 1\{m>0\}
\right]^{1-\alpha_0},
\]

and utility remains

\[
\boxed{
u_t(c,s;m)
=
\frac{[Q(c,s;m)/e(m)]^{1-\sigma}}{1-\sigma}
+\psi_t m,
\qquad \sigma=2.
}
\]

Equivalently, at the maintained curvature,

\[
u_t(c,s;m)=-\frac{e(m)}{Q(c,s;m)}+\psi_t m.
\]

The scale coefficients remain externally fixed. One estimated parenthood housing requirement replaces the two housing-requirement parameters. Constant \(\alpha_0\), zero nonhousing floor, and the existing direct benefit from children are retained. Other lifecycle and fertility components are outside this change.

For an initial comparison preserving the previous first-child requirement, use \(\bar h_P=h_1+h_m\). Its final value would be re-estimated.

### 6. What this solution implies

Among parents, there is no further unavoidable increment in minimum housing for another child. Additional children still increase material needs through the scale.

For an interior renter with current expenditure \(X=c+rs\), optimal housing among parents satisfies

\[
s^*
=
\frac{(1-\alpha_0)X}{r}
+\alpha_0\bar h_P.
\]

Therefore, at the same current expenditure and rental price, additional children do not independently shift spending toward housing after the first child.

This does not imply that larger families occupy the same-sized homes in the full model. They can choose different total expenditure, saving, tenure and locations. Housing costs also continue to affect fertility because housing remains part of the scaled bundle.

At fixed current expenditure, the material utility cost of additional children declines among parents under the maintained curvature. This is an intended implication of sharing. It does not determine the full dynamic family-size distribution on its own.

### 7. Why this is a defensible choice

The additional housing slope should be retained only if there is a convincing reason that every child requires an unavoidable increment of housing services. A first-birth housing response does not independently establish that requirement.

The jump-only specification makes a simpler economic distinction between entering parenthood and expanding an existing family. Its implications should be assessed against housing differences among parents, family-size distributions, and consumption and saving behavior.

The literature supports the ingredients, although not an exact replication of our combined model:

- Scholz et al. provide the borrowed scale shape. Their full utility aggregation differs and should not be attributed to our model.
- Dustmann et al. combine household scaling with a housing minimum.
- De la Croix and Pommeret provide an endogenous-fertility precedent for inside consumption scaling with a separate child reward.

The agreed solution is therefore to **retain the concave scale and use a parenthood-only housing requirement**. The pension correction provides the occasion to recalibrate the revised specification, but it is a separate issue: the housing change is an economic simplification, whereas the pension change repairs fiscal accounting.