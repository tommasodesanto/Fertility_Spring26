# Population closure: what the literature actually does, and what that implies

Synthesis of a two-scout literature survey (2026-08-06) on population closure
for the housing–fertility model. Verification discipline: items below marked
"read directly" were verified against the paper's own text this session;
"publisher-verified" means journal/year/abstract confirmed; working papers
labeled as such. Companion document:
`docs/prompts/PROMPT_population_closure_fresh_take_20260806.md` (give the
fresh chats that prompt FIRST, this memo only afterward, to avoid anchoring).

## 1. The headline: every branch of the literature refuses this question

The population claim we want — endogenous below-replacement fertility moving
a stationary population in general equilibrium — is made by no published
paper we could find, because each adjacent literature closes it off:

- **Quantitative fertility-policy models are partial equilibrium.** Doepke
  and Kindermann (AER 2019; read directly): wages are exogenous individual
  draws, no market clearing anywhere; policy experiments hold all prices
  fixed. Sommer (JME 2016; read directly): exogenous wage process, interest
  rate a calibrated constant (1.04), no markets at all. The Handbook chapter
  (Doepke, Hannusch, Kindermann, Tertilt 2023; read directly) keeps
  "household fertility choice" and "macro consequences of demography" in
  separate sections and calls the latter "in the early stages."
- **Structural housing–fertility papers fix the population count.**
  Moreno-Maldonado and Santamaría (CEPR DP 19002, working paper; read
  directly): the equilibrium definition takes the demographic composition
  \(L_t(a,z)\) as an exogenous *input*; fertility timing is a household
  choice, but births never become future households. Greaney, Parkhomenko,
  Van Nieuwerburgh (NBER WP 33512, working paper; read directly): unit mass
  of households fixed forever, exogenous turnover, newborn households enter
  with zero wealth, bequests leave the economy — architecturally our model's
  closest relative, with turnover exogenous.
- **Demographic-GE (pension/aging) models make fertility exogenous.** The
  whole Auerbach–Kotlikoff lineage runs transitions on projected vital
  rates; when a terminal steady state is needed it is manufactured by
  assumption — Kitao (JEDC 2015; read directly) literally assumes Japanese
  fertility recovers so cohort growth hits zero by 2150. Immigration, where
  present, is an exogenous policy-set flow: Storesletten (JPE 2000; read
  directly) — government picks an age/skill inflow rule; nobody models entry
  as a value-based choice against an outside option.
- **Urban/spatial models fix the national total and reallocate it.**
  Hsieh–Moretti (AEJ:Macro 2019; read directly): aggregate labor normalized
  to one, each city faces infinitely elastic supply at the common utility.
  Roback (JPE 1982; read directly): the "open city" takes reservation
  utility as given *because the system of cities pins that utility*; it is a
  partial-equilibrium device for one location inside a closed system.
- **Dynasty models (Barro–Becker, Econometrica 1989)** close population
  through fertility choice itself — but that requires dynastic altruism, not
  a warm-glow lifecycle architecture, and admits no immigration margin.

Two consequences. First, the "suspect" feeling about our construction is
explained: there is no convention to inherit — anything we do here is a
modeling contribution and must be defended on its own economics. Second, the
author's requirement that population respond to policy already exceeds the
published frontier (the nearest competitor holds the count fixed while
letting timing move); the right response to novelty is maximal simplicity
and transparency in the closure.

## 2. Two sharpened arguments the survey delivers

**(a) The outside-VALUE layer is a category error for a national model.**
The fixed-outside-utility device is the open-*city* closure, and in its own
tradition (Roback) it is licensed only because a surrounding system of
locations determines the reservation utility. A one-region national model
has no system to pin \(\bar V^{out}\); importing the open-city device there
answers "which city grows," not "does the country have more children." Its
honest home is the future spatial paper, where the reservation utility is an
equilibrium object across locations — exactly Roback's original setup.

**(b) Scarcity under population change has a standard device, and we
already have it.** The urban literature keeps land scarcity meaningful while
population moves via *non-proportional* local supply — the Saiz (QJE 2010)
geography/regulation elasticities that parameterize Hsieh–Moretti. Our
levels supply curve \(H^S = H_0 (r/\bar r)^\xi\) is that device: when the
stationary population rises 3 percent, demand climbs a fixed schedule and
prices rise. A balanced-growth path would force per-capita (proportional)
supply and give exactly this up — which settles the BGP question on
literature grounds, not just taste: the growth-rate object also has no
reporting precedent as a policy outcome, while the scarcity margin has a
canonical empirical anchor.

## 3. Recommended closure, restated with its supports

**Stationary equilibrium; entrants = matured city-born children (two
children form one household) + an exogenous immigration inflow \(M\); \(M\)
calibrated so outside-origin entrants are 16.9 percent of the flow (the ACS
across-CBSA measure already adopted in May); counterfactuals hold \(M\)
fixed; population then responds to policy only through births, with
multiplier \((1-s^{out})/s^{out} \approx 4.9\).** The outside value, taste
scale, and entry logit are retired to default-off code.

Supports, one line each: the inflow-as-exogenous-policy-object has the
Storesletten precedent; the "stationary despite below-replacement fertility,
because immigration" statement is the actual US demographic regime, not a
trick; the turnover architecture is DUE's, with our single change being that
turnover is generated by the model's own births; the fixed-population
decomposition rows remain reportable as the Moreno-Maldonado–Santamaría-
comparable object; and the novelty — fertility policy moving stationary
population against a scarce housing supply — is exactly the gap both scouts
confirmed. What it gives up: local/city experiments (which belong to the
spatial model anyway) and any migration *response* margin (unidentified in
any case; the survey found no credible identifying moment for it at the
national level).

## 4. One trap to preempt, and two options we still owe fresh eyes

- **Welfare across different population sizes is a known minefield**
  (Handbook Section 7: efficiency concepts change when the set of people
  alive differs across allocations). Keep welfare statements per-household
  for fixed cohorts; never total-utilitarian across counterfactuals with
  different populations.
- The fresh-take prompt asks two questions this memo does not settle: are
  the symmetric-outside and exogenous-inflow closures truly equivalent in
  every reportable object, and is there an off-menu closure we have not
  considered? If either fresh chat surfaces something new, it should be
  weighed against Section 1's finding that simplicity is the binding
  criterion.

## 5. Sources (as verified this session)

Read directly: Roback (JPE 1982); Hsieh–Moretti (AEJ:Macro 2019);
Doepke–Kindermann (AER 2019); Sommer (JME 2016); Doepke–Hannusch–
Kindermann–Tertilt (Handbook 2023, §7–8); Storesletten (JPE 2000); Kitao
(JEDC 2015); Greaney–Parkhomenko–Van Nieuwerburgh (NBER WP 33512, 2025,
working paper); Moreno-Maldonado–Santamaría (CEPR DP 19002, working paper).
Publisher-verified: Saiz (QJE 2010); Barro–Becker (Econometrica 1989);
De Nardi–İmrohoroğlu–Sargent (RED 1999); Attanasio–Kitao–Violante (JME
2007); Blanchard (JPE 1985). Field-journal leads, unread: Okamoto (J. Japanese
Int. Economies 2021 — endogenous fertility + immigration lever, transition-
based); Hagiwara (Japanese Economic Review 2025 — family-policy GE-OLG,
closure unverified).
