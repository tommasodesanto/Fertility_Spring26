# Calibration reconciliation: the July-22 provisional ledger and the E-series architecture

Date: 2026-07-22 (evening). Author: lead agent session, independent assessment requested by Tommaso.
Status: conceptual reconciliation only. No model code changed, no calibration launched.
The running Torch diagnostic `14579021_[1-8]` (collector `14579038`) was not touched.

REVISED 2026-07-22 (late), after author review. Four corrections: (i) the p90
"impossibility" argument is reframed as a model-analogue order constraint —
the data values are not mutually contradictory; (ii) the housing-response
analogue on the ACTIVE Markov path is a horizon-0 within-birth-period
diff-in-diff, not the 12-year window an earlier draft described (that block
is legacy-path-only); (iii) the parity-adjusted child-cost number is removed
— it must come from the A2 CEX rerun, not a back-of-envelope; (iv) a
fertility-units decision section is added: the sequential architecture
inherits `tfr = 2 x mean parity` while treating parity increments as literal
single births, and this must be resolved before the Part-3 table is frozen.
**The 12/12 table is a proposal, not a settled contract**, pending the A2/A3
reruns and the units decision.

Sources verified for this note (not taken on trust): `new_moment_profile.py`,
`m5_profile.py`, `calibration_search.py`, `solver.py` moment implementations
(lines cited below), `run_e1_chain.py` DOMAIN, the E-series preference/fertility
kernels, `latex/calibration_strategy_provisional.tex`,
`docs/model/eqscale_seq_specification_note_20260720.tex`,
`docs/model/fecundity_schedule_note_20260720.tex`,
`docs/model/sequential_fertility_design_note_20260718.tex`,
`output/model/eqscale_seq_optimized_recalibration_20260719/report/` (collected
winner), the CEX/PSID empirical drivers, and De Nardi and Yang (2014), NBER WP
20058, Table 2 and surrounding text.

---

## 0. Source-of-truth discrepancies found during verification

These are reported per the AGENTS.md rule before any substantive proposal.

1. **CALIBRATION_STATUS.md is missing the E2 collection.** The July 19 section
   says "E2 overnight"; the July 20 spec note reports the collected result
   (loss 5.49) and `output/model/eqscale_seq_optimized_recalibration_20260719/report/`
   contains the strict winner (`rank_loss = 5.485557807431926`, 8/8 chains
   eligible). The canonical status file was never updated. Until it is, "E2"
   facts live only in the spec note and the output directory.
2. **The `c_bar_n` bound claim is stale.** CALIBRATION_STATUS says the CEX
   target 0.043655 "is just below the current `c_bar_n` search lower bound
   0.05". The live domain used by both `--profile m5` and `--profile
   new-moments` is `("c_bar_n", 0.0, 3.0, "softzero")`
   (`calibration_search.py:54`, enforced against `promotion_contract.py:18`).
   The 0.05 bound exists only in the historical `production_profile.py`
   bounds table, which is not wired into the running search. The CEX target is
   NOT clipped by the live domain.
3. **`e1_profile.py` constants are stale as a "winner".** `E1_THETA` /
   `E1_LOSS = 12.608` are the parent-fork E1 point retained as a numerical
   parity benchmark for the optimized port. The collected production winner is
   loss 5.4856 with a very different theta (e.g. `kappa_fert` 1.80 vs 10.77).
   Anyone reading `e1_profile.py` as the current E-series answer will be wrong.
4. **`intergen_eqscale_seq_optimized/IMPLEMENTATION_STATUS.md` is inherited
   boilerplate** and states fertility "is not a sequential parity hazard",
   which is false for this package.
5. **The four brand-new ledger moments have no `TARGET_MOMENT_OBJECTS`
   documentation entries** in `calibration.py` (the model/data/status ledger
   used for every other target).

---

## 1. Audit of the provisional 14-parameter / 14-moment ledger (current M5 architecture)

Verified implementation facts used throughout:

- The four CEX LES quantities are **free parameters with identity moment rows**
  (`calibration.py:1373-1376` returns the live parameter values as "moments";
  all four appear in the search DOMAIN). They are soft penalty rows (weight
  `1/target^2`), not first-stage restrictions. Only `theta_n = 0` is external.
- The tenure object is `E[p(1-p)]` over the full beginning-of-period state
  distribution, ages 25–55, one 4-year transition, where `p` is the summed
  owner-option probability of the multinomial tenure logit with scale
  `tenure_choice_kappa` (`solver.py:5439-5473`).
- `tfr = 2.0 * mean_parity` pooled over all post-fertility-window ages;
  parity states are {0, 1, 2+}; `childless_rate = parity_dist[0]` on the same
  pooled population.
- `housing_increment_0to1` on the ACTIVE Markov-income path is a horizon-0
  difference-in-differences: birth-cohort rooms within the birth period minus
  a no-birth control (`forward_distribution_markov_income`,
  `solver.py:4284-4299`; `production_profile.py:22` pins
  `PRODUCTION_HOUSING_EVENT_HORIZON = 0`, inherited by `m5_overrides`).
  Cohort = all parity-0 exits including direct jumps to parity 2+
  (`solver.py:4348`). The 12-year `event_horizon = 3` block exists only in
  the legacy non-Markov `forward_distribution` (`solver.py:3782`, `3909`) and
  must not be quoted for the running system; the `_eventstudy_t3` stat-name
  suffix is a legacy misnomer on the active path.
- `aggregate_wealth_to_annual_after_tax_earnings`: wealth = `b + p*H_own`
  (gross housing) over all ages; earnings = working-age after-tax labor income
  only, pensions and lump-sum transfers excluded (the July-22 fix, verified at
  `solver.py:5476-5536`).
- `annual_bequest_flow_to_aggregate_wealth`: genuinely decedent-weighted —
  age-specific death probabilities (survival < 1 only from 66), forced
  bequest of the full surviving mass at the terminal age (82), nonnegative
  estates, flow annualized by /4.
- `old_total_estate_wealth_to_annual_income_p90_7684`: p90 of `(b + p*H)/own
  annual gross income` over the **living** population aged 76–84
  (`solver.py:5078-5160`). No death weighting. Own-income denominators use the
  retirement-income branch and include the lump-sum transfer.

### Block A: saving and bequests — (beta, theta0, theta1) vs (6.90, 0.0088, 4.53)

The block design is soundly precedented: De Nardi and Yang (2014, Table 2) use
exactly these three moments to calibrate (beta, phi1, phi2), and their stated
definitions are: wealth-to-after-tax-earnings from Hendricks (2007), the
bequest-flow-to-wealth ratio from Gale and Scholz (1994), and the 90th
percentile of the **bequest distribution at death, single decedents**,
normalized by income, from Hurd and Smith (2002).

**beta <-> 6.90.** Identifying variation: patience shifts saving at every age,
so the aggregate wealth stock relative to the earnings flow is the natural
level moment; it replaces young liquid wealth, which honest income risk breaks
(E2: 0.52 vs 0.18). Definitional mapping in code is clean (earnings-only
denominator, pension exclusion verified). Two genuine concerns:
(i) *Reachability.* DNY reach 6.90 in an economy with wide permanent earnings
heterogeneity, intergenerational ability persistence, and no housing. The M5
income process carries roughly a third of the PSID innovation s.d.
(0.0645), and E2's untargeted `wealth_to_income` is 2.40. The model object may
sit far below 6.90 everywhere in the box, in which case beta will pin at its
0.9995 annual bound and the row becomes a constant torque on the whole system
rather than a moment. Check the running collector for `beta_annual` at bound.
(ii) *Cross-loading through the price level.* Wealth includes `p*H_own`; any
parameter moving the equilibrium price (H0, chi, alpha block) mechanically
moves the numerator. Joint-discipline language is required, as in the M5
identification audit.

**theta0 <-> 0.0088.** Sound as the level-of-bequests moment; the flow is
correctly decedent-weighted, and with `theta_n = 0` external the warm glow is
child-blind, matching the DNY concept. Three honest caveats to disclose:
mortality exists only from age 66, so the model flow has no young-death
component; the forced terminal bequest at 82 truncates late-life dissaving
that the data flow embodies; and part of the model flow is accidental
(incomplete annuitization), so theta0 is identified by the flow *conditional*
on the survival schedule — fine, but the 0.0088 is then partly demographic
bookkeeping. None of these is disqualifying for a diagnostic.

**theta1 <-> 4.53. This row is broken as implemented, three ways.**

1. *Population concept.* The model moment is the p90 of living 76–84
   wealth-to-own-income. The 4.53 is the p90 of **estates at death of single
   decedents** normalized by income. These are different objects in the data
   by nearly an order of magnitude: the project's own July-16 decomposition
   measured the PSID living 76–84 estate-wealth p90 at ~36 median-income
   units, while Hurd–Smith single-decedent estates give 4.53. Within the
   model the two populations coincide (mortality is age-only), which is
   precisely why the *data* concept must be chosen deliberately — the model
   cannot bridge them.
2. *Model-analogue order constraint (framing corrected in the revision).*
   As data, 4.53 and the PSID median 6.50 are NOT mutually contradictory:
   they come from different samples and concepts (single-decedent estates at
   death vs living 76–84 households) under different normalizations. The
   bite is entirely model-side: both targets are assigned to one and the
   same model statistic family — quantiles of the living-76–84 own-income
   wealth ratio — and within that family p90 >= p50 holds identically. So
   the 4.53 target is unreachable by construction whenever the median
   analogue stays anywhere near its previously targeted (and E2-delivered)
   6.5, and chasing it means crushing the whole model estate distribution
   against the bequest-flow and beta rows. Predicted symptoms in the
   collector: theta1 at a bound and/or large opposing residuals across rows
   1–3. The root cause is the population/normalization mismatch of item 1 —
   the model analogue erases a distinction the two data numbers depend on —
   not any inconsistency in the data themselves.
3. *Denominator convention.* The model divides by own annual gross income
   (retirement branch, transfer-inclusive) at 76–84; DNY normalize by economy
   income. Retiree own-income denominators inflate the model ratio relative
   to the DNY convention even at identical wealth.

*Recommendation:* keep theta1 on a **dispersion ratio, not a level**: the
project-measured PSID living 76–84 estate p90/p50 = 3.448 (bootstrap SE known
from the July-15 work; the model already computes
`old_total_estate_wealth_to_annual_income_p90_p50_7684`, E2 untargeted value
3.92). The ratio isolates the curvature/luxury margin theta1 actually governs,
is scale-free (no double-counting of the level already assigned to theta0),
and is measured on this project's own data with a defensible model analogue.
Retain 4.53 only as a cross-reference after a genuine estates-at-death
reconstruction, if ever.

### Block B: the four CEX LES rows — first stage in name, penalty in practice

As implemented these are prior-style penalty rows: the parameter is free, the
"moment" is the parameter itself, the weight is `1/target^2`. Two consequences:

- The "14 moments for 14 parameters" count is bookkeeping. The four identity
  rows carry no model information; the effective system is 10 model-generated
  moments for 10 structural parameters plus four penalized parameters that
  also move all 10 moments. That is a perfectly respectable penalized-SMM
  design, but the strategy note's "one-to-one mapping" language and the
  implementation should be reconciled explicitly: either fix the four
  parameters (true two-step, 10/10 system, with a bootstrap/delta-method
  first-stage SE feeding a sensitivity appendix), or keep the penalty rows and
  say so, with weights justified by the first-stage SEs rather than
  `1/target^2`.
- The good news: pinning `h_bar_0` (and `alpha_cons`) attacks exactly the
  weakest directions of the M5 Jacobian audit (the rank-9/14 weak direction
  loaded on `h_bar_0`, `theta1`, `theta0`, `tenure_choice_kappa`). This is the
  single best identification improvement in the ledger.

Row-specific issues:

- **alpha_cons = 0.733, c_bar_0 = 0.397, h_bar_0 = 2.630.** The LES auxiliary
  regression (`rent ~ allocatable expenditure + cross-fitted price per room`)
  and the mapping (alpha = 1 - b, h_bar_0 = d/(1-b), c_bar_0 = -a/b) are
  internally consistent with the model's interior-renter LES demand. The
  quiet assumption is that CEX childless renters are interior Stone–Geary
  consumers; in the model, renters face `hR_max = 6` and tenure selection.
  For a paper target the honest analogue is the same regression on simulated
  interior renters; for the diagnostic the identity rows are fine.
  Note the distances from the M5 winner: alpha 0.591 -> 0.733, c_bar_0 1.259
  -> 0.397, h_bar_0 0.392 -> 2.630. This is a *large* reparameterization of
  the demand system; expect substantial re-equilibration. A welcome side
  effect: c_bar_0 = 0.397 four-year units is ~0.10 annual — the exact "low
  bundle" the July-18 feasibility frontier showed relaxes the income-risk cap
  to sigma ~ 0.18 even without a floor.
- **c_bar_n = 0.0437 is the weakest row in the ledger.** Three compounding
  problems. (i) *Reduced-form vs structural:* the CEX estimate is
  d E[nondurable spend | income, rooms, age]/d n; even inside the model this
  is not `c_bar_n` — conditional on total expenditure the LES split loads
  only (1-alpha) of a floor shift onto consumption, and conditional on income
  the response is smoothed intertemporally. Equating coefficient and
  parameter is a category error the ledger already flags as a placeholder.
  (ii) *Units:* the CEX coefficient is per dependent child; the model
  parameter multiplies the parity index, and model parity 1 is the empirical
  "1–2 children" bin (the rooms-gap row and the 2x TFR bridge both encode
  this). At ~1.5–1.7 children per parity-1 family the per-parity target
  should be ~0.065–0.075, not 0.0437. (iii) *Consumption concept:* CEX
  nonhousing nondurables vs model `c` = all nonhousing consumption. All three
  are fixable by re-specifying the moment as a parity-binned conditional
  increment replicated identically on model data (see Part 3/4). Meanwhile
  the gap to the M5 winner (0.393, i.e. 9x) makes this the row most likely to
  be sacrificed by the running search.

### Block C: children and housing needs

- **h_bar_jump <-> first-child rooms response 0.664.** CORRECTED in the
  revision: the active analogue is a **horizon-0 within-birth-period
  difference-in-differences**, not a 12-year window. The Markov-income path
  that M5 and the running diagnostic use reads
  `event_horizon = getattr(P, "housing_event_horizon", 0)`
  (`solver.py:4285`) and the production/M5 overrides pin it to 0
  (`production_profile.py:22`), with an explicit code comment that the PSID
  ~3-year post-birth controlled response lands inside the 4-year birth
  period, so horizon 0 is the deliberate match and the no-birth control
  reduces to the contemporaneous committed-rooms jump. The hardcoded
  `event_horizon = 3` block an earlier draft of this note described exists
  only in the legacy non-Markov `forward_distribution` and is inert for the
  running system; the `_eventstudy_t3` suffix is a legacy misnomer worth
  renaming. Remaining flags: the birth cohort pools all parity-0 exits
  including direct 2+ entrants (`solver.py:4348`; the one-child-only variant
  is a separate stat), and the mapping is triangular, not one-to-one: the
  parity-1 committed-rooms increment is `h_bar_jump + h_bar_n`, so this row
  identifies the sum, the next row identifies `h_bar_n`, and `h_bar_jump` is
  the residual. The strategy note should say "triangular block", not two
  independent one-to-one rows.
- **h_bar_n <-> 3+ vs 1–2 rooms gap 0.368.** Direct and persuasive given the
  documented bin bridge (model parity 2+ maps to empirical 3+; parity 1 to
  1–2). Cross-loads with tenure composition (owners sit on discrete rungs),
  but the primary assignment is credible.

### Block D: fertility — (psi_child, kappa_fert) <-> (TFR-equivalent 1.918, childlessness 0.188)

The pair is the standard level/dispersion split and is fine as a primary
mapping. Three flags. (i) The 2x parity-to-births bridge is a strong units
convention; it should be stated once, prominently, in the strategy note — all
per-"child" parameters are actually per parity state. (ii) Both data moments
are cohort stocks; the model pools all post-window ages of the stationary
population — consistent under stationarity, worth one sentence. (iii) The
one-shot completed-parity logit makes kappa_fert both the dispersion and the
*only* smoothing device; in the E-series winner the analogous configuration
drifts to psi < 0 with fertility carried by logit noise. A
fertility-gradient validation moment (fertility by income or by housing-cost
exposure) should police "fertility as noise" in both architectures.

### Block E: tenure and supply

- **chi <-> prime-age ownership 0.575.** Clean. Cross-loads: phi, entrant
  wealth, kappa_T, price level. No objection.
- **kappa_T <-> E[p(1-p)] vs cross-fitted Brier 0.1171.** The structural
  analogue is a *lower bound* for the data object at any given kappa_T, for
  two reasons pushing the same direction: (a) the model conditions on the
  complete state, and E[Var(O|full state)] <= E[Var(O|X)] for the coarser
  PSID X by the law of total variance; (b) the data Brier adds the logit's
  approximation error E[(p_hat - p_X)^2] on top of E[Var(O|X)]. Hence
  calibrating E[p(1-p)] to 0.1171 **overstates kappa_T** — and may not even
  be reachable: at the M5 winner kappa_T = 0.0100 the tenure choice is nearly
  deterministic and E[p(1-p)] is tiny, while the domain caps kappa_T at 0.12.
  Predicted symptom in the running collector: kappa_T pinned at 0.12 with a
  large residual on this row. If confirmed, the correct response per the
  identification discipline is not to drop the row but to (i) implement the
  identical simulated auxiliary exercise, and (ii) re-run the *data* side
  with model-feasible covariates only — the PSID spec conditions on `married`
  and `year`, which do not exist in the model; dropping them raises the data
  Brier somewhat and closes part of the definitional gap from the other side.
  Under the E-series, honest income risk (sigma = 0.20) generates real
  state-driven tenure churn, so far less taste noise is needed — the row and
  the architecture fit together much better there.
- **H0 <-> aggregate occupied rooms 5.78.** Clean supply-scale mapping;
  eta external. The "18_85" label is cosmetic (terminal model age is 82).

### Ledger-level points

- **Weights.** `1/target^2` equal relative gaps are fine for a diagnostic and
  deliberately transparent. For a paper system, project-measured rows have
  bootstrap SEs (tenure 0.0021; estate ratios; CEX first-stage) and should be
  weighted by them; borrowed rows need declared synthetic SEs. Keep the
  promised joint-covariance step on the roadmap.
- **Demotions are properly replaced** (beta: young wealth -> aggregate wealth;
  theta0/theta1: estate median/nonhousing share -> flow/p90; kappa_T: old-age
  ownership -> tenure dispersion), so the AGENTS.md replacement rule is
  satisfied in form. But the demoted set contains the paper's central
  mechanism moment (family ownership gap) and the old-age
  portfolio-composition failure (the July-16 diagnosis). Both must survive as
  *named validation moments with acceptance bands*, not silently vanish from
  reports. The upside is real: keeping the mechanism moment out of the loss
  turns it into an honest overidentification test.

**Bottom line for Part 1.** The architecture of the ledger — a DNY saving
block, a CEX-disciplined demand block, direct physical child-housing moments,
a fertility level/dispersion pair, one tenure-dispersion moment, one supply
moment — is the right shape and a clear improvement over the M5 target
system's tangle. Three rows are not usable beyond the diagnostic as
implemented: the p90 row (wrong population, internally impossible target,
wrong denominator convention — replace with the PSID p90/p50 ratio), the
c_bar_n row (reduced-form-as-structural, wrong units, wrong consumption
concept — re-specify parity-binned), and the tenure row (known analogue gap,
acknowledged; add the covariate-contract fix). The four identity rows should
either become true first-stage restrictions or keep their penalty
interpretation explicitly, with SE-based weights.

---

## 2. What the E-series replaces, and what survives

### Architecture deltas (verified in code)

| Object | M5 (production) | E-series (`intergen_eqscale_seq_optimized`) |
|---|---|---|
| Child costs | Floors `c_bar(n,s)`, `h_bar(n,s)`: 5 SMM params | Share tilt `alpha(n,s) = clip(alpha0 - (d_jump + d_a*n)*1{raising}, .05, .95)`, scale `e = 1 + gamma_e*n`: 3 SMM params. Children have **zero budget-side cost**; preferences only, active only while a child is at home |
| Feasibility | Committed bundle must be affordable; gates/censoring | Every y > 0 feasible; no floor apparatus |
| Income risk | sigma_z = 0.0645 (feasibility-capped) | sigma_z = 0.20, rho = 0.9602, external |
| Fertility | One-shot completed-parity logit (choice set {0,1,2+}, childless state only) | Sequential: childless {wait, try}; one-child {stop, try-second} while the child is at home; success prob pi_j external (Leridon-anchored); second birth restarts the single child clock |
| Free params / moments | 14 / 15 | 12 / 15 (E1/E2 contract) |

E2 collected winner (strict, repeated; 8/8 chains eligible): loss 5.4856 vs
M5's 9.0444 on the identical 15 moments. Winner theta (all interior):
beta_annual 0.9846, alpha0 0.6595, d_jump 0.0388, d_a 0.0216, gamma_e 0.2025,
psi -0.3000, kappa_f 1.798, chi 1.0760, H0 7.155, theta0 0.1480, theta1
0.2851, kappa_T 0.0231. Resolved relative to M5: old-age ownership (0.742 vs
0.764 target), estate median (6.55 vs 6.50), estate p90/p50 untargeted 3.92 vs
data 3.45, old nonhousing share (0.605 vs 0.608). Remaining misses: young
renter liquid wealth 0.521 vs 0.179 (the honest-risk cost), **family
ownership gap 0.042 vs 0.168 (the paper's mechanism moment, off 4x)**,
childless renter rooms 4.29 vs 3.81, childlessness 0.229 vs 0.188, 3+/1–2
rooms gap 0.536 vs 0.368.

### Disposition of the 14 ledger rows under the E-series

**Survive unchanged (7):** wealth/after-tax-earnings (beta), bequest
flow/wealth (theta0), TFR-equivalent, childlessness, prime-age ownership
(chi), tenure dispersion (kappa_T), aggregate rooms (H0). None of these
depends on the Stone–Geary or one-shot structure.

**Redefined (4):**
- p90 estate row -> p90/p50 living 76–84 ratio (theta1), for the Part-1
  reasons, identical in both architectures.
- `alpha_cons_parameter` -> the LES *slope* evidence still identifies alpha0
  (for childless interior renters the model's housing expenditure share is
  exactly 1 - alpha0), but the identity-row shortcut should become either a
  hard first-stage fix or a simulated auxiliary slope.
- The two rooms rows survive as data moments but now identify
  (d_jump + d_a) and d_a — the share-tilt block — instead of
  (h_bar_jump + h_bar_n) and h_bar_n. Same triangular structure.

**Dead (3):** `c_bar_0_parameter` and `h_bar_0_parameter` have no
counterpart — there are no intercepts. They must not be force-mapped onto
anything. The CEX LES intercept and price coefficient become the sharpest
available *specification test*: the E-series predicts a = 0 and d = 0 for
childless renters; the data reject that (a < 0, d > 0). Report as an honest
overidentification contrast — it is the empirical case that committed bundles
exist, i.e. the cost of the E-series' feasibility gain.
`c_bar_n_parameter` likewise dies as a parameter row.

**Re-purposed (1):** the CEX child *consumption-increment* evidence is
exactly what identifies `gamma_e`. Economics: conditional on total
expenditure, children tilt spending toward housing (d_a lowers the
consumption share); the only force that makes measured consumption spending
*rise* with children conditional on income — as the CEX says it does
(+$1,289/child/year) — is the scale channel: `e(n) > 1` raises the marginal
utility of the whole composite, pulling expenditure into child-present
periods (cutting saving). So within the 3-parameter block:
rooms 0->1 response and 3+/1–2 gap discipline (d_jump, d_a) — the allocation
tilt; the conditional consumption increment disciplines gamma_e — the level
scale. This replaces the old, wrong reading of the same evidence as a
committed-consumption floor. The moment must be re-specified parity-binned
(0 vs 1–2 children) and replicated identically on model data.

**New fertility objects.** The sequential architecture makes the previously
banned timing margins measurable and already computes them (attempt and
realized first- and second-birth hazards by age, first-birth age
distribution, mean age at first birth 28.18 at the E2 winner, share of first
births 30+ = 0.484, flow parity progression 1->2 = 0.285, chosen/clock
childlessness split 0.132/0.097). `parity_progression_1to2` is now a
legitimately defined flow object for the first time — the AGENTS.md
prohibition applies to the one-shot architecture only. Fecundity
(omega1, omega2, close at 45) stays external: biology identifies it, and
pinning it is what keeps kappa_f identified (both smooth hazards).

---

## 3. Proposed reconciled calibration system for the E-series (proposal, not settled — see the revision note)

### 3.1 Externally estimated or fixed (never in the SMM)

| Object | Value / source |
|---|---|
| Timing, ages, child duration | 4-year period; 18/66/82; 18-year maturation |
| sigma (CRRA) | 2 |
| Survival | SSA 2023 from 66 |
| Income process | rho = 0.9602, sigma_z = 0.20 (SS range; the E-series' point) |
| Fecundity | (omega1, omega2) least-squares fit to Leridon four-year anchors (pending half-day), close at 45 |
| Housing finance | phi = 0.80, sale cost 6%, depreciation 1.1%, property tax 1.0% |
| Real return | 2.0% annual |
| Housing menus | H_own = {2,4,6,8,10}, hR_max = 6 (owner-ladder density caveat stands) |
| Supply elasticity | eta = 1.75 (Saiz) |
| Entrant wealth | PSID 18–24 childless renters, frontier-censored |
| Bequest child shifter | theta_n = 0 |

### 3.2 Internally calibrated (12) and their hard moments (12) — exactly identified

| # | Parameter | Hard moment | Target | Status |
|---|---|---|---|---|
| 1 | beta | Aggregate wealth / annual after-tax labor earnings | 6.90 (borrowed DNY; remeasure from FoF/NIPA before paper) | port moment to E-package |
| 2 | theta0 | Annual bequest flow / aggregate wealth | 0.0088 (borrowed, Gale–Scholz; label as borrowed) | port |
| 3 | theta1 | Estate p90/p50, living 76–84, own-income normalized | 3.448 (PSID, project-measured, SE known) | stat exists in E-package |
| 4 | alpha0 | Childless-renter LES expenditure slope (1 - b) | b = 0.267 -> alpha0 = 0.733 (CEX) | first-stage fix or simulated auxiliary slope |
| 5 | d_jump | First-child rooms response, horizon-0 diff-in-diff (block with #6) | 0.664 (PSID) | exists; active path matches the PSID window by construction |
| 6 | d_a | Rooms gap, top family-size bin vs parity-1 bin, ages 30–55 (bin labels pending the Section-3.5 units decision) | 0.368 (data) | exists |
| 7 | gamma_e | Conditional consumption increment on the model's family-size bins, income-controlled, identical model replication | **TBD — must be re-estimated from the CEX on the chosen bin convention (A2); no number is adopted here** | new data + model moment |
| 8 | psi | Completed-fertility equivalent (2x mean parity) | 1.918 | exists |
| 9 | kappa_f | Completed childlessness | 0.188 | exists |
| 10 | chi | Ownership rate, ages 30–55 | 0.5755 | exists |
| 11 | kappa_T | Cross-fitted 4-year tenure Brier, model-feasible covariates, identical simulated auxiliary regression | re-estimate data side without `married`/`year` (current all-covariate value 0.1171, SE 0.0021) | new machinery |
| 12 | H0 | Aggregate occupied rooms per household | 5.780 | exists |

Block structure: {1,2,3} x (beta, theta0, theta1); {5,6,7} x (d_jump, d_a,
gamma_e) — triangular in the same way as the floors block; {8,9} x (psi,
kappa_f); rows 4, 10, 11, 12 one-to-one. Nothing is underidentified; the
count is 12/12 with every external restriction named.

### 3.3 Hard-moment design decisions (and why)

- **CEX evidence.** Slope -> alpha0. Intercept and price coefficient ->
  specification test, never a target. Child consumption increment ->
  gamma_e, on the model's family-size bins, identically replicated. This is
  the proposed answer to "how should the CEX child-cost evidence identify
  the equivalence-scale parameters": the *allocation* evidence (rooms)
  identifies the share tilt, the *level* evidence (conditional consumption
  increment) identifies the scale, and the *nonhomotheticity* evidence
  (intercept, price coefficient) is the test the specification is allowed to
  fail. Stated honestly, the rooms->deltas / increment->gamma_e assignment
  is a **conjectured primary loading inside a jointly identified 3x3
  block**, not an established one-to-one map: the increment coefficient also
  absorbs the share tilt, the saving response, housing adjustment, and
  fertility selection into parenthood. The assignment is adopted only
  conditional on (i) the identical simulated auxiliary regression (B2) and
  (ii) the T5 Jacobian confirming the block's loading pattern; if the
  Jacobian shows the increment loading mainly on the deltas, the block is
  re-blocked rather than relabeled.
- **Fertility moments for a sequential architecture.** Keep the two stocks
  (TFR-equivalent, childlessness) as the hard pair: the sequential model
  still determines them, and psi/kappa_f remain the level/dispersion
  parameters. Do **not** promote timing moments to hard targets in the
  baseline: every candidate timing parameter is already externally pinned
  (fecundity by biology) or already disciplined (kappa_f by the stocks), so a
  hard timing row would either double-count or silently re-identify biology.
  The timing battery — mean age at first birth, share of first births 30+,
  first-birth hazard profile, flow PP 1->2, chosen/clock split — is the
  overidentification test of the postponement mechanism, which is the paper's
  actual claim. Costed option if the parity composition misfits after
  recalibration: split psi into per-parity (psi1, psi2) and promote the flow
  PP 1->2 to a hard third fertility row (3 params / 3 moments). Not
  recommended now: E2's parity composition is roughly right, and the E2 tell
  (psi < 0, fertility carried by logit dispersion) argues for *fewer* free
  utility knobs plus a gradient validation moment, not more.
- **Weights.** SE-based where measured (CEX first stage, PSID tenure and
  estate rows, rooms moments); declared conservative relative SEs for the two
  borrowed DNY rows; document as the interim step to the promised joint
  covariance.

### 3.4 Overidentifying validation battery (reported with acceptance bands, never in the loss)

Mechanism: family ownership gap 0.168 (E2 currently 0.042 — see risk below);
fertility-income gradient (to be measured; polices noise-driven fertility).
Demand system: LES intercept and price coefficient contrast; childless renter
mean rooms 3.805; owner-renter rooms gap 2.419; childless owner 6+ share
0.596. Lifecycle/tenure: own_rate_2534 0.341; old-age ownership 0.764; tenure
switch rate. Saving: young renter liquid wealth 0.179 (disclosed E-series
miss); estate median 6.50; old nonhousing >= 1x income share 0.608. Timing:
the full battery above vs NCHS/NSFG/PSID counterparts (data pass pending).

### 3.5 Residual risks, stated plainly

1. **The mechanism moment collapsed under eqscale.** Floors created a
   parent-specific housing *need*; share tilts only re-slice a budget. E2's
   family ownership gap is 0.042 vs 0.168. If no (d_jump, d_a, gamma_e)
   region restores a substantial gap while holding the rest, the equivalence-
   scale form fails the paper's central mechanism and needs a committed-
   housing device (e.g. a minimum-rooms constraint for parents — a
   constraint, not a utility intercept, preserving feasibility). This must be
   probed *before* full recalibration (Part 4, test T4).
2. **psi < 0 at the E2 winner.** Children are flow-utility bads and fertility
   is carried by logit dispersion. The stocks pair cannot exclude this
   configuration; only a gradient moment can. Measure the fertility-income
   (or housing-cost-exposure) gradient and gate acceptance on its sign.
3. **gamma_e vs psi collinearity.** Both scale the utility consequence of
   children. The consumption-increment row is what separates them (gamma_e
   has a budget-behavior signature; psi does not). If the Jacobian still
   shows a joint weak direction, profile gamma_e.
4. **Borrowed saving targets under honest risk.** 6.90 may be unreachable
   (E2 untargeted wealth/income 2.40); watch beta at bound in any run and
   remeasure the target properly before paper use.
5. **Fertility units are unresolved and gate rows 5–9 of the table.**
   (Expanded in the revision; this is a decision, not a footnote.)

   *The fact pattern.* The E-package inherits `"tfr": 2.0 *
   mean_completed_fertility` unchanged (`calibration.py:1276-1287`,
   `n_parity = 3`), and E2's collected tfr = 1.981 = 2 x 0.9905 confirms it
   is live. Under the one-shot model the x2 could be read as a bin bridge
   (parity 1 ~ the empirical 1–2-children family, parity 2+ ~ 3+), because a
   childless household could jump directly to the top bin. Under the
   sequential architecture, parity increments are **literal single births by
   construction** — the first attempt lands at parity 1, the second attempt
   moves 1 -> 2 — and the design notes speak of literal first and second
   births. The code's flow bookkeeping is literal too (births are counted
   parity-weighted, no doubling), while the headline fertility moment
   doubles. These two conventions cannot both be right.

   *The arithmetic that forces a decision.* Under the literal reading with
   the current cap (parity <= 2), dropping the x2 makes the completed-
   fertility target unreachable: with childlessness 0.188, max mean parity =
   2 x 0.812 = 1.624 < 1.918. So there are exactly two coherent options:

   **Option B (bins, no code change).** Declare parity 1 = empirical 1–2
   children, parity 2 = 3+. Then tfr = 2p1 + 4p2 is an approximate
   births-per-household bridge (E2: 1.98, target 1.918); first-birth timing
   moments stay valid (entering parenthood is a literal first birth under
   both readings), and the rooms-gap row keeps its current empirical bins.
   The cost: the model's "second-birth" margin is empirically the
   progression from the 1–2 bin to 3+, i.e. a THIRD child. The flow PP
   "1->2" (E2: 0.285) must then be validated against the data 1–2 -> 3+
   progression (roughly 0.3 in completed CPS-type distributions), NOT the
   literal second-birth progression (roughly 0.8) — coincidentally, 0.285
   sits near the former and would be a catastrophic miss against the
   latter, which is itself weak evidence the model currently behaves like
   the bin reading. All per-parity parameters (gamma_e, delta_a, psi) are
   per bin, and the spec note's "gamma_e = 0.203 reproduces the OECD child
   weight" claim is overstated: per actual child it is ~0.12–0.13 against
   the OECD-modified 0.3.

   **Option L4 (literal, extend parity to 0/1/2/3+).** Drop the x2; target
   completed fertility = mean parity directly (feasible with a 3+ state);
   the second-birth margin, PP 1->2 (~0.8), spacing, and the tempo/quantum
   narrative all become literal and usable as advertised in the design
   notes. Cost: one more parity state (~+33% on that state dimension), a
   third attempt margin with a fert3 side channel, KFE splits, and
   re-measured child-binned moments. Unlike the costed v2 redesign this
   does NOT require per-child age tracking — parity is already a state
   dimension and the youngest-child clock approximation carries over.

   *Recommendation.* Option B for the next diagnostic round (zero code
   change, but relabel the second-attempt objects and remap PP validation
   to the 3+ progression); Option L4 for the paper if second-birth timing
   is to carry the tempo/quantum story the sequential design notes promise.
   The literal-cap-2-without-doubling reading is arithmetically excluded.
   Until this is decided, the bin definitions in rows 5–9 and every timing
   moment's data counterpart are provisional.

---

## 4. Implementation plan (no launch)

### Phase A — measurement and contracts (no model code)

A1. Freeze this reconciled contract as the E-series target-system spec:
    names, definitions, parity-bin conventions, the 2x bridge, weight rule.
A2. CEX: rebuild the child-cost moment on the family-size bins chosen in A8
    (both binned reads: first bin vs childless, top bin vs first bin),
    income/age/rooms/tenure-controlled, in model units with the documented
    income normalization; bootstrap SEs; keep the per-child linear slope as
    robustness. No target number exists until this runs. Driver:
    `code/data/cex_child_cost/build_child_cost_target.R` (extend, do not
    overwrite).
A3. PSID tenure: re-run `build_tenure_residual_variance.R` with model-
    feasible covariates only (drop `married`, `factor(year)`); freeze fold
    seeds and the exact spec in a versioned contract file consumed by both
    the data and (future) simulation sides; report both covariate sets.
A4. Estate row: adopt p90/p50 = 3.448 with its July-15 bootstrap SE; retire
    the 4.53 target from any model-facing role pending a real estates-at-
    death reconstruction.
A5. Fecundity: run the pending least-squares fit of (omega1, omega2) to the
    Leridon anchors; re-solve the E2 winner at the fitted schedule as the
    robustness line (the spec note's half-day item).
A6. Timing battery data pass: NCHS/NSFG mean age at first birth and 30+
    share; PSID flow PP 1->2 in 4-year windows; SEs; validation-only.
A7. Saving-block provenance: either keep DNY values explicitly labeled
    borrowed, or remeasure wealth/after-tax-earnings (FoF B.101 net worth /
    NIPA compensation net of taxes) and bequest-flow (literature survey) —
    decide before the paper system, not before the next diagnostic.
A8. Fertility-units decision memo (gates A2, A6, and rows 5–9): choose
    Option B (bins) vs Option L4 (literal 0/1/2/3+) per Section 3.5 risk 5;
    if L4, scope the `n_parity = 4` extension (state-dimension growth, third
    attempt margin and fert3 side channel, KFE splits, re-measured binned
    moments) before any model edit.

### Phase B — model code (E-package only; production untouched)

B1. Port from `intergen_housing_fertility_optimized` into
    `intergen_eqscale_seq_optimized`: `add_aggregate_wealth_bequest_flow_moments`
    (with the pension-exclusion earnings denominator), the tenure
    residual-variance function, and their unit tests. The p90/p50 stat
    already exists in the E-package.
B2. Implement the parity-binned conditional consumption-increment analogue
    (deterministic, from the stationary distribution: within income-state x
    age cells, parity-bin consumption differences aggregated with fixed
    weights mirroring the data spec). No simulation needed.
B3. Implement the childless-renter LES slope analogue (deterministic:
    interior-renter housing expenditure share). Decide alpha0 treatment
    (first-stage fix vs penalty row) and encode it.
B4. Simulated-panel auxiliary machinery for the exact cross-fitted tenure
    Brier (seeded, deterministic). Required before any *paper* calibration;
    the deterministic E[p(1-p)] may stand in for diagnostics only, with the
    bias direction disclosed.
B5. New profile module (e.g. `reconciled_profile.py`): 12 targets, SE-based
    weights, `require_identified(12)`, complete `TARGET_MOMENT_OBJECTS`
    entries for every row (also backfill the four missing entries in the M5
    package).
B6. Housekeeping (same edits, no behavior change): add the missing July-20
    E2 section to CALIBRATION_STATUS.md; correct the stale c_bar_n-bound
    sentence; replace the inherited `IMPLEMENTATION_STATUS.md` in the
    E-package; comment `e1_profile.py` constants as parity-benchmark-only.

### Phase C — tests gating any Torch search

T1. Unit tests per new moment (hand-computed small configs); deaths-
    accounting identity for the bequest flow; parity-bin consistency between
    rooms and consumption rows; auxiliary-regression determinism (bitwise
    repeat under fixed seed).
T2. Nesting: with the new moments wired as diagnostics-only, the E2
    15-moment evaluation must be bit-identical to the collected winner.
T3. Reachability precheck: one strict solve at the E2 winner reporting all
    12 reconciled rows plus the validation battery; any row >3x off triggers
    a target-or-analogue review before search (this is where 6.90 and the
    Brier get their first honest look).
T4. Mechanism probe: bounded profile over (d_jump, d_a) (other coordinates
    fixed or nuisance-reoptimized, ~25 cells, the July-13 pattern) asking
    whether own_family_gap >= ~0.10 is reachable at all. If not, stop and
    redesign the child-housing device before any recalibration.
T5. Identification smoke: 12x12 SMM-weighted finite-difference Jacobian at
    the E2 winner under the reconciled system (the July-15 audit pattern);
    require rank 12 at 1e-3, report condition number and the weakest
    direction; pre-register the response if deficient (externalize or
    re-block the implicated parameter).
T6. Only after T1–T5: exact-loop smoke, then an 8-chain Torch search with
    strict-repeat reserves, checkpoints, budgets, and a collector enforcing
    the full target-fit table — the standing Long-Run Search Safety rules.

### Post-run checklist for the live diagnostic `14579021` (do not touch the run)

When the collector lands, audit against the predictions here: (1) is
`beta_annual` at 0.9995? (2) is `kappa_T` at 0.12 with a large tenure-row
residual? (3) does the p90 row dominate the loss and drag theta1 to a bound
while the flow row resists? (4) how far do the four LES parameters drift from
their CEX targets (penalty-row trade-off)? (5) where does `c_bar_n` settle
relative to 0.0437 vs the M5 winner's 0.393? Each confirmed prediction is
evidence for the corresponding Part-1 fix; each refuted one weakens it.
