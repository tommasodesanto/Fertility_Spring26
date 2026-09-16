---
title: "How the literature models the benefit and the cost of children, and what it implies for P1"
subtitle: "Fertility–housing project, decision note"
date: "September 16, 2026"
---

# Purpose

The consolidated review (P1) found that the model's benefit of children is
linear in children at home while its cost is concave, so that a household
without taste shocks wants zero children or the biological maximum, and the
one-child and three-plus shares are carried by the continuation taste scale.
This note does the reading that choice deserves. Every functional form below
was read from the paper's own PDF by a verifier and quoted with a page; papers
not on disk are marked and cited only for a headline claim. The July 18
readings memo (`equivalence_scales_sequential_fertility_readings_20260718.tex`)
already covers the equivalence-scale literature (Scholz–Seshadri–Khitatrakun,
Fernández-Villaverde–Krueger, Attanasio et al., Comin–Lashkari–Mestieri) and is
not repeated.

Notation: \(m\) children at home, \(n\) children ever born, \(X=c^{\alpha}(s-\bar h)^{1-\alpha}\)
the composite, \(e(m)=((2+0.7m)/2)^{0.7}\) the scale, \(\sigma=2\). The model's
flow utility is \(u=-e(m)/X+\psi m\).

# 1. What the model does, in one line each

Benefit: \(\psi m\), a constant flow \(\psi\) per child per period while the
child is at home. Cost: the scale \(e(m)\) deflating the composite, so at
\(\sigma=2\) the cost of the \(k\)-th child in units of \(1/X\) is
\(e(k)-e(k-1)=0.234,\,0.216,\,0.203\), falling in \(k\); plus a parenthood-only
space floor \(\bar h(m)=h_P\mathbf 1\{m>0\}\), which costs only at the first
child; plus a one-time utility cost \(\xi\) at the first birth. Hence the net
flow benefit of the \(k\)-th child rises in \(k\). No income gradient is built
in: the scale is homothetic and the floor is an absolute quantity, so relative
to income the cost of children falls as income rises.

# 2. The benefit side: three families

**Concave in the number ever born.** Barro and Becker (1989, p. 2):
\(U_i=v(c_i)+a(n_i)\,n_i\,U_{i+1}\) with \(a(n)=a\,n^{-\varepsilon}\), so the
weight on children is \(n^{1-\varepsilon}\), concave. Sommer (2016, p. 34):
\(U=c^{1-\gamma}/(1-\gamma)+\zeta(nq)^{1-\kappa}/(1-\kappa)\), concave in
\(n\) at given quality, entering every period the children are at home. De la
Croix and Doepke (2003, p. 2): \(\ln c+\beta\ln d'+\gamma\ln(nh')\). Baudin,
de la Croix and Gobbi (2015, p. 1861): \(\ln c+\ln(n+\nu)\), with \(\nu>0\) so
that childlessness is defined. Kim, Tertilt and Yum (2024, p. 1587):
\(\phi(n)\log(h'-\chi\tilde h')\) with \(\phi(n)\) nonparametric for
\(n=1,2,3\). Daruich and Kozlowski (2020, p. 226): \(b(n)u(c_k)\), \(b\)
increasing and concave. The handbook exposition (Doepke, Hannusch, Kindermann
and Tertilt 2023, p. 161) uses \(\log c+\delta\log(nh)\). This is the
majority form: diminishing marginal benefit gives an interior family size
from preferences alone.

**Linear, one-time at birth.** Doepke and Kindermann (2019, p. 3284):
\(u=c-d+v\,b\), utility from a child only in the period of birth, with
\(v\) heterogeneous across couples and a veto rule. Here the intensive margin
comes from heterogeneity in \(v\) and from per-child costs, not from
curvature. Moreno-Maldonado and Santamaría (2024, p. 17) likewise put an
idiosyncratic multiplicative child taste \(\eta^i\) on utility and choose
among young, postponed and never.

**Linear flow while at home.** This is the model's form. I found no paper on
disk that uses exactly it. The closest are the two linear one-time forms
above, and both pair linearity with per-child costs that are linear or
convex in the number of children.

The lesson is not that linear is wrong; it is that linear benefits are always
paired in the literature with either heterogeneity in the taste (Doepke–
Kindermann, Moreno-Maldonado–Santamaría) or costs that do not fall in
\(k\). The model has the taste heterogeneity (the fertility shocks) but
pairs it with falling costs.

# 3. The cost side: four devices and their curvature in \(k\)

**Time cost proportional to the wage.** De la Croix and Doepke (2003, eq. 1):
income \(w h(1-\phi n)\), \(\phi=0.075\) per child. Baudin et al. (2015):
a fixed time cost \(\eta\) of becoming a parent plus \(\phi\) per child, no
goods cost. Kim, Tertilt and Yum (2024, eq. 7): \(\lambda n\) of parental time
plus education spending \(xn\). Doepke and Kindermann (2019): the mother's
foregone wage or the childcare price while a child is under three, plus a
fixed goods cost \(\phi_c\) and a fixed utility cost \(\phi_u\) per child.
Moreno-Maldonado and Santamaría embed an estimated child earnings penalty in
the income process. Curvature: linear in \(n\) (constant marginal cost), and
proportional to the wage, which is what produces the negative income–fertility
gradient (Jones, Schoonbroodt and Tertilt 2010, not on disk, as cited by
Baudin et al. fn. 16: with a log form and no non-labor income the gradient
needs a wage-proportional cost).

**Goods cost per child.** Barro–Becker: \(P_i\) per child. Couillard (2025,
eq. 4.1): \(c+ni+nrh=y\), a per-child goods investment and a per-child
housing cost \(rh\) per child. Curvature: linear in \(n\).

**Equivalence scale.** Two versions. Deflating consumption only, \(U(c/e)\):
Kim, Tertilt and Yum (p. 1588, OECD-modified weights), Borella–De Nardi–Yang
(divisor \((j+0.7f)^{0.7}\)), and this model. Multiplying utility as well,
\(e\,U(c/e)\): Scholz, Seshadri and Khitatrakun (2006, p. 615,
\(E\sum\beta^{j-S}n_jU(c_j/n_j)\)). Curvature at \(\sigma=2\) with the
Citro–Michael scale: deflating only gives cost increments \(0.234, 0.216,
0.203\), falling; multiplying gives \(-e^2/X\), increments \(0.522, 0.580,
0.630\), rising. The sign of the curvature flips with the outer weighting.
Neither SSK nor Borella et al. choose fertility, so neither paper had to
confront this; Kim, Tertilt and Yum do choose fertility and pair the deflating
scale with a concave \(\phi(n)\) and linear time and education costs.

**Minimum quality or space.** Sommer (2016, p. 5): \(q\ge\bar q\) if \(n>0\),
with quality produced from time and goods, so the floor costs \(n\bar q\)
worth of inputs, linear in \(n\). Couillard: \(rh\) per child, linear. This
model: \(h_P\) once, at the first child, with the per-child slope restricted
to zero in the current vintage. Curvature: the literature's floors are per
child; ours is a step.

Summary of curvature. Every paper that chooses fertility has a marginal cost
that is constant or rising in the number of children. Ours is the only one
with a falling marginal cost, and it comes from the one device (a deflating
scale under \(\sigma>1\) with exponent below one) that the fertility papers
either pair with curvature elsewhere or do not use.

# 4. The income gradient

Cross-sectional fertility falls with income in the United States. In the
models above this comes from a cost proportional to the wage (time cost or
earnings penalty), or from quantity–quality with a high income elasticity of
quality (Becker–Lewis, not on disk; the handbook's eq. 2). A homothetic scale
gives a flat gradient; an absolute floor gives a rising one. The model has
neither device, and the E6b diagnostic on July 27 found the childlessness
gradient across permanent types reversed relative to the CPS. This is the same
symptom as P3's \(\xi\): a first-birth cost in utility units is income-neutral.

# 5. Three candidate specifications

Each is stated as primitives, with what it adds or removes and which moments
identify it. None is implemented by this note.

**S1. Concave benefit, everything else unchanged.** Replace \(\psi m\) by
\(\psi\,v(m)\), \(v(m)=\log(1+m)\) or \((1+m)^{1-\varepsilon}\) with
\(\varepsilon\) fixed externally. Net flow benefit of the \(k\)-th child then
falls if \(v\) is concave enough to beat the falling scale cost; with
\(v=\log(1+m)\) the benefit increments are \(0.693, 0.405, 0.288\), against
cost increments \(0.234, 0.216, 0.203\), so the net benefit falls in \(k\)
and an interior family size exists. Adds no parameter if \(\varepsilon\) is
fixed. Identified by the one-child and three-plus shares, which then stop
being \(\kappa_C\)'s job. Does not fix the income gradient.

**S2. The SSK weighting.** Keep \(\psi m\), replace \(U(X/e)\) by
\(e\,U(X/e)\), so flow utility is \(-e(m)^2/X+\psi m\). Costs now rise in
\(k\); an interior size exists with a linear benefit. Adds no parameter and
matches the cited paper's objective. Two costs. The first child's utility cost
more than doubles, so \(\psi,\xi,\kappa_1,\kappa_C\) all move and the
calibration restarts. And it puts the negative level of CRRA utility to work:
multiplying a negative number by family size penalizes children through the
level normalization, the same objection ChatGPT raised against child-scaled
bequests; Daruich and Kozlowski (fn. 11) and Jones–Schoonbroodt–Tertilt note
that Barro–Becker weighting needs positive utility to be meaningful. S2 is
defensible only if the paper states this and shows the \(\sigma=2\) level is
not driving the result. Does not fix the income gradient.

**S3. A wage-proportional child cost.** Keep the current benefit and scale;
add a time or earnings cost while \(m>0\), \(y^d(1-\tau_c(m))\), with
\(\tau_c\) taken from the child-penalty literature (Kleven et al., not on
disk; Moreno-Maldonado–Santamaría estimate it on U.S. data), and make the
space floor per child as in Couillard. Marginal cost becomes constant in
\(k\) and proportional to the wage, so the income gradient appears and the
intensive margin is set by cost against a linear benefit plus shocks, the
Doepke–Kindermann configuration. Adds one external schedule \(\tau_c(m)\) and
possibly retires \(\xi\) (P3). Identified by fertility and childlessness by
income or education (CPS), which the model currently does not target, and by
the three-plus versus one-to-two rooms gap for the per-child floor.

**What I would do.** S3 first, because it is the only one that answers both
objections (curvature and gradient) with an externally measured object, and
it is the configuration the two quantitative fertility references use. Then
test whether S1's curvature is still needed for the one-child share. S2 is a
one-line change worth running as a diagnostic to show the author what the
cited objective implies, not a baseline. Whatever is chosen, the paper's
preference paragraph must state the benefit form, the four cost devices in
play, and why the marginal cost does not fall in the number of children.

# 6. What each reference actually assumes, for the paper's citations

| Paper | Benefit of children | Cost of children | Scale | Children as |
|---|---|---|---|---|
| Barro–Becker 1989 | \(n^{1-\varepsilon}U_{i+1}\), concave | goods \(P_i\) per child | none | ever born, once |
| Sommer 2016 | \(\zeta(nq)^{1-\kappa}\), concave, per period | time and goods into quality; floor \(\bar q\) per child | none | at home, stock |
| Doepke–Kindermann 2019 | \(v\,b\), linear, at birth, \(v\) heterogeneous | fixed goods \(\phi_c\), fixed utility \(\phi_u\), mother's wage under 3 | none | birth flow; stock for care cost |
| de la Croix–Doepke 2003 | \(\gamma\ln(nh')\) | time \(\phi n\) at own wage; education at teacher wage | none | ever born, once |
| Baudin et al. 2015 | \(\ln(n+\nu)\) | time: fixed \(\eta\) for first child, \(\phi\) per child | none | ever born |
| Kim–Tertilt–Yum 2024 | \(\phi(n)\log(\cdot)\), nonparametric | time \(\lambda n\), education \(xn\) | deflates \(c\) only, OECD | ever born, once |
| Daruich–Kozlowski 2020 | \(b(n)u(c_k)\), concave | \(C(h,n)\) rising in wage | none | ever born, at fixed age |
| Couillard 2025 | general \(U(n,q,c,h)\) | goods \(ni\), housing \(nrh\) per child | none | number of children |
| Moreno-Maldonado–Santamaría 2024 | multiplicative taste \(\eta^i\) | earnings penalty by timing | none | timing choice |
| Scholz et al. 2006 | none (no fertility choice) | scale only | multiplies and deflates | at home |
| Borella et al. 2023 | none | scale only | deflates | at home |
| This model | \(\psi m\), linear, per period at home | deflating scale; step floor \(h_P\); one-time \(\xi\) | deflates only | \(n\) and \(m\) |

Not on disk and cited from memory only: Becker and Lewis (1973), Greenwood,
Seshadri and Vandenbroucke (2005), Adda, Dustmann and Stevens (2017), Guner,
Kaygusuz and Ventura (2020), Jones, Schoonbroodt and Tertilt (2010), Kleven et
al. (2019). Verify before any of them appears in the paper.
