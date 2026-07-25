# Why the model does not postpone births, and what would make it

2026-07-25. Decision memo for the mechanism discussion. Evidence base: the
kappa_E-profile probe and wealth slice of 2026-07-25 (ledger, "OVERNIGHT
PROBE VERDICTS") and the certified E5 run. Nothing here is implemented;
every option needs the author's sign-off.

## 1. The evidence

With entry noise switched almost off (kappa_E <= 1), every household tries
for a first birth at age 18: mean age at first birth 18.2, childlessness
0.000. As noise rises the mean jumps to ~26 and stays there, with the
30-plus share stuck near 0.36. Nowhere on the frontier does the model
approach the data pair (mean 25.3, share 0.270), and childlessness never
exceeds 0.118 against 0.188. Conclusion: in the deterministic model,
immediate entry dominates for essentially everyone; all postponement the
calibrated model displays is taste noise, and prices play no role in the
timing of births.

## 2. Why immediate entry dominates

Consider the deterministic trade-off at age 18 between trying now and
waiting one period.

Gains from trying now: the per-child utility flow psi starts earlier and
runs longer (children stay ~18 years regardless of when they arrive);
fecundity is highest now and only falls (0.98 at 18 to 0.50 at 42), so an
attempt now is strictly more likely to succeed; with three children
possible sequentially, starting early preserves the option of a large
family.

Costs of trying now: the equivalence scale divides the consumption
composite, and with sigma = 2 the utility cost of a child is proportional
to e(n)-growth over current resources -- children are more painful in
utility terms when the household is poor, which is the one force pushing
toward delay. The housing tilt raises desired rooms, which is expensive for
a young renter.

The probe says the psi-plus-biology side wins everywhere: the
poverty-of-the-young channel is too weak at plausible psi to delay anyone.
Two things the model simply does not have: children cost no time and no
earnings (income is an exogenous profile, so there is no career cost), and
nothing about a birth requires housing the household cannot rent (renters
can carry any family size up to the 6-room cap). The two empirical
first-order reasons to wait -- career and housing readiness -- have no
teeth in the current architecture.

## 3. What the data require

A first-birth age distribution with mean 25.3 and only 27 percent of first
births after 30: substantial delay, but a tighter right tail than logit
noise produces. And 18.8 percent never entering. So the mechanism must
deliver (i) a reason to wait that fades with age, (ii) a persistent group
that never enters, and (iii) less reliance on the noise scale, so that
prices can matter.

## 4. Candidate mechanisms

### A. Career cost of children (child penalty)

Births reduce the earnings path (level or growth) while children are home
or permanently. Standard empirically (Kleven-Landais-Sogaard child-penalty
estimates; Adda-Dustmann-Stevens 2017 JPE structurally). Buys: genuine
postponement (waiting until the flat part of the profile is cheaper),
price-relevant timing, and some childlessness among steep-profile states.
Cost: earnings must depend on fertility history -- either a cheap version
(a proportional penalty while children are at home, one parameter, no new
state) or the honest version (permanent penalty from first birth, which
adds birth-timing dependence to income). Identification: the timing pair
plus the child-penalty magnitude imported externally. The cheap version is
one parameter and reuses existing states.

### B. Permanent heterogeneity in psi (taste types)

Two or three permanent types drawn at entry, one with psi at or below
zero. Baudin-de la Croix-Gobbi (2015 AER) style: childlessness becomes
partly chosen by type rather than manufactured by noise. Buys: the
childlessness level directly (the 0.118 ceiling breaks), and lets kappa_E
fall because noise no longer has to produce non-entrants -- which tightens
the 30-plus tail toward the data. Does not by itself delay the entrants.
Cost: one type dimension, solve time roughly times the number of types; no
new dynamics. Identification: childlessness pins the low-type share;
completed fertility pins psi of the parent type.

### C. Readiness arrival replacing the entry logit

Entry possible only after a stochastic "settled" event (partnership, first
stable job) with an age-dependent arrival hazard; the entry logit shrinks
toward deterministic. This is the structural version of what kappa_E
currently proxies ("partnership timing and career events"). Buys: the
timing distribution takes the hazard's shape, so the mean and tail can be
matched jointly -- the one option that directly fixes the pair. Cost: one
binary state with an exogenous hazard (two or three parameters); fits the
existing sequential machinery naturally. Identification: the two timing
rows pin the hazard location and spread; childlessness still needs B or A.

### D. Make housing gate fertility (the thesis mechanism)

Currently a renter can raise any family within the rental cap, so the
down-payment margin never binds on the entry decision. Options: a
family-size-dependent ownership premium chi_O(n), or a rooms requirement
for children that the rental market cannot supply beyond a threshold.
Buys: births become price-sensitive -- the paper's central channel starts
to exist -- and saving-for-downpayment becomes a real postponement motive.
Cost: this changes the economics of the tenure block and would need the
family-housing moments re-examined; highest risk, highest relevance.
Identification: the family ownership gap and the rooms moments, which the
system already carries.

## 5. Recommendation

B plus C is the standard combination and the cheapest route to a system
that fits: types deliver childlessness, the arrival process delivers
timing, and both let the noise scales shrink until prices matter. A (cheap
version) is a defensible substitute for C with better economics but an
imported parameter. D is not an alternative to these -- it is the thesis
mechanism and should be evaluated at the E5b winner via the policy packet
first: if prices move births once noise is moderate, D may be partially
present already through wealth effects; if they do not, D is what the
paper needs, with B/C supplying the demographic margins.

Separate note: no candidate above touches the late-life p90/p50 failure
(2.2 vs 3.45 at every visited point); that requires an earnings or returns
tail and is logged separately in the ledger.

Citation hygiene: the three references named here are from memory anchors
used in prior sessions; verify each against the PDF before any of this
enters paper prose.
