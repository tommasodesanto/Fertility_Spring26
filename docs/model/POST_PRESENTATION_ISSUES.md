# Post-presentation issues ledger

## Author decision — September 12, 2026

For the September presentation, use the historical May empirical plots and
original regression specification. Tommaso has not had time to review the
subsequent measurement and design revisions. This is a presentation-version
decision, not certification that the historical measurement is correct.
Suggested disclosure: **Original May specification; measurement revisions
under review.** Do not describe that specification as timing-corrected.

The May rooms graph is reproducible: all 18 saved points match the recognized
original specification to numerical precision. The rooms survey-year assignment
problem is verified against the assembled source panel, rather than merely
conjectured. Preserve both findings when discussing the historical results.

The new annual and binned comparisons remain diagnostics. Do not replace the
calibration target or its weight, relabel the historical graph, or promote a
new specification without reconciling the empirical and model definitions.
Empirical follow-up is deferred; the completed data jobs and their collection
automation are stopped/paused. No new empirical regressions are launched by
this decision. Resume the items below after the presentation review.

## Numbers that must not be conflated

| Object | Rooms | Interpretation |
|---|---:|---|
| May graph coefficient at +3 | 0.796858555 | Original normalization omits both −2 and −6. |
| May graph +3 minus estimated −1 | 0.740737457 | Four-year contrast calculated from that same historical curve. |
| Current pinned calibration target | 0.720246262 | Later August specification; +3 minus −1, different sample/weighting/control construction. |
| Historical slide scalar | 0.66 / 0.664 | Does not equal the saved May graph's +3 coefficient; source remains to be reconciled. |

The comparable-horizon scalars 0.740737457 and 0.720246262 differ by
0.020491195 rooms (about 2.8% of the May contrast). This numerical proximity
does not establish equivalent estimands or validate either specification.
Tommaso's preferred pre-birth reference remains −2; subtracting −1 answers a
different question. A −2-to-+3 comparison spans five calendar years, whereas
the currently pinned model mapping spans four years.

## Deferred empirical work

All rows remain **open / deferred**, except the reproduction result recorded
above. Close each with a reproducible result and an explicit specification
decision; smoothness or statistical significance alone is not a criterion.

| Priority | Issue and established evidence | Required resolution |
|---|---|---|
| 1 | Rooms answers are stored one interview early in the merged panel. Every assigned value in the 352,250-row prepared common sample matches its survey-year source after alignment. | Reconstruct rooms directly from the source-year crosswalk; recover the 52,317 source-observed values not restored by shifting; validate year-specific response codes. Preserve May reproduction separately. |
| 1 | Annual event-time support alternates after PSID becomes biennial. Some cohorts lack the intended −2 reference; May also omits −6. | Choose exact years versus two-year windows, retain a reference before the final pre-birth year, and document cohort support and aggregation weights. The tested baseline −3/−2 is a changed object, not exactly −2. |
| 1 | The first binned diagnostic is much smaller: 149,402 fitted rows versus 345,751 in the prior annual common fit. | Review the restrictions individually: exclude 2019+ dates, exclude no-recorded-birth-year observations, require all six displayed windows, then estimator exclusions. Do not attribute the full difference to binning or data transfer. |
| 1 | The binned driver excludes 101,848 rows without a recorded first-birth year after its date restriction. These may include childless people and unknown histories. | Reconstruct history status before deciding which observations can be controls. Compare admissible last-treated and confirmed-childless designs separately; report their populations and identifying assumptions. |
| 1 | A last-treated 2019 cohort cannot remain an untreated comparison group in 2019 and later. | Either restrict admissible dates/cohort comparisons or choose and justify another control group. Preserve this change separately from outcome timing. |
| 2 | Ludovica code constructs rounded frequency weights and uses them in csdid blocks; its separate Sun–Abraham command reproducing May has no weight argument. August uses direct IW probability weights. | Compare unweighted and correctly defined survey-weighted estimates on identical samples; state the population represented. Do not claim the original analysis generally forgot weights. |
| 2 | August restricts to women who are reference persons/spouses, selects one per household-year, and excludes multi-family-unit dwellings. | Decide the observational unit and population; quantify each restriction, repeated household outcomes, and appropriate weighting/clustering. These are not interchangeable with mechanical corrections. |
| 2 | August defines first biological birth across 20 child records and all available histories, unlike the original first-child field. | Compare birth dates and exclusions person by person; distinguish measurement corrections from a changed biological-parent population. |
| 2 | HOMEOWN has a recovered 41-wave source-year mapping; no evidence justifies applying the rooms shift to it. Original regressions exclude tenure 'neither owns nor rents'. | Numerically validate the mapping, decide the ownership denominator, then rerun ownership with the chosen event design. For ownership transitions, use the previous observed interview rather than calendar-year L. |
| 2 | Moving variables lack a recovered upstream year crosswalk. In 2019 the question covers moves since January 2017 and dates the most recent move; recall conventions vary by vintage. | Recover/validate source mapping, move dates and recall intervals; distinguish most recent move from any move, first-mentioned reason from all reasons, non-movers from missing responses. Bin only after establishing the outcome clock. |
| 2 | Reason code 3 includes expansion/better housing; code 6 includes neighborhood, schools and proximity to friends/relatives. | Check vintage-specific codes and use accurate outcome names. Revisit move gating, missing-to-zero errors and reason-response denominators. |
| 3 | Presentation numbers, empirical contrasts and calibration mapping differ. | Reconcile 0.66/0.664, 0.797, 0.741 and 0.720; choose the intended horizon and population; regenerate the target estimate, covariance-based uncertainty, weight and provenance together before recalibration. |

## Model and slide issues — September 13, 2026

Items raised while editing the September 14 deck. Numbering continues the mock-feedback ledger (M01–M18 in `MOCK_PRESENTATION_FEEDBACK.md`). None is authorized for implementation by the slide task.

### M19 — Tenure smoothing $\kappa_H$

`tenure_choice_kappa = 0.005` is the lower search bound set on June 28, not an interior estimate. Searches returned the bound; the June sweep showed that 0.05 changes the economics materially (old-age ownership near 0.90). The only dedicated moment ever proposed was the PSID four-year-ahead ownership Brier score (0.117113, SE 0.002102); the matching simulated-history exercise was never implemented. Removed from all September 14 slides.

**Decision needed:** treat as external numerical smoothing with a stated value, or implement the auxiliary prediction target and estimate it.

### M20 — Parameter-table vintage

Deck tables report $\psi_0 = 0.160$ and the earlier parameter vintage (e.g. $\kappa_1 = 0.279$, $\xi = 0.127$). The September 13 corrected initial packet gives $\psi_0 = 0.149$, $\kappa_1 = 0.338$, $\xi = 0.266$. Refresh all initial tables together or not at all.

### M21 — Recent-parent ownership gap

Removed from the Identification slide and both target tables on September 13. It remains a weighted target in the live run (13 rows; deck shows 12). If it is dropped from the calibration, name the replacement discipline for the housing block; if kept, the paper table must show it.

### M22 — Owner housing grid in exposition

The draft budget-constraint slide states continuous sizes $h \in [0,\bar h]$; the household-problem and state slides still write $h_{t+1} \in \{0,H_1,\dots,H_K\}$. Choose one exposition (discreteness as computational detail) and harmonize.

### M23 — Equilibrium definition loose ends

Sequence definition adopted September 13. Open: the set $\{V_t,g_t,G_t,N_t,P_t,r_t,T_t,\varpi_t\}$ is not named element by element; $N_t$ is first defined inside the definition after the Population and Households slide was cut; the rebate condition $T_t\int dG_t=\tau_t^p P_t\int h_{t+1}dG_t$ should read $T_t=0$ in the no-rebate baseline; $\mathcal A_t,\mathcal F_t$ defined in words only; four appendix person-accounting frames are now unlinked. Supersedes the exposition part of M14/M16.

### M24 — Preferences cross-partial

Slide states $\partial^2 u_t/\partial s_t\partial m_t>0$ on the utility, not on the aggregator, because $\mathcal C_{sm}$ has ambiguous sign in the implemented form (housing requirement raises it, equivalence scale lowers it). Confirm the sign at $\sigma=2$ or drop the bullet.

### M25 — Identification mapping unverified

Slide assigns $\xi \to$ childlessness, $\kappa_1 \to$ first-birth timing, $\kappa_C \to$ one-child families, $h_P \to$ both room moments. This is a reading of the parameter table, not a Jacobian check. Verify against the local weighted Jacobian before the paper.

### M26 — Population law in slides versus code

The fitted 2007–2023 history imposes observed household age masses (births/2.1 entry queue rescaled to data); the post-2023 forecast converts annual persons to heads with fixed 2023 ACS headship. The deck now says neither. Decide how much to state on the calibration slide.

### M27 — Underwater-debt rollover

No unsecured credit line ($\lambda_d=0$). Debt below the collateral floor arises only after a price fall or a sale with shortfall and rolls over at share $\lambda_{a+1}$ (1 before age 42, linear to 0 at 62). Removed from the deck as niche; belongs in the paper appendix.

### M28 — Schematic frames overfull

Initial Steady State, Impact, and Demographic Adjustment each overflow by 18pt (minipage heights). Cosmetic.

### M29 — Supply elasticity source

The deck cites Baum-Snow and Han (2024) for $\eta = 0.63$. Their headline is an average urban floor-space supply elasticity near 0.5; the independent quantitative audit found no primary receipt deriving 0.63. Establish the derivation or change the value/citation before the paper.

### M30 — First-birth fixed cost $\xi$

The author was not aware the model carries a one-time utility cost at the first birth (`first_birth_fixed_cost`, default zero, estimated at 0.127 in the deck vintage and 0.266 in the September 13 corrected initial). It exists because the first-child housing jump $h_P$ is pinned by rooms moments and childlessness needed its own lever. Decide whether to keep it as a fixed cost of parenthood, replace it with a per-period time/goods cost of children, or test whether childlessness can be matched by $h_P$ alone. Requires recalibration; not for the September 14 deck.

### M31 — Earnings process documentation

Correction (September 14): the live calibration layer (`code/model/intergen_eqscale_seq_optimized/local_panel.py`) builds the earnings state as a five-point Rouwenhorst discretization of an AR(1) crossed with three permanent income types; the three-point grid with persistence 0.85 in `parameters.py` is a module default that the calibration overrides. The slides say only "discretized AR(1)" and omit the permanent types. Open: (i) document the annual persistence/innovation parameters and their four-year conversion in one place with a source; (ii) verify that the PSID-based permanent-type variance and the literature-based persistent process do not double count dispersion; (iii) decide whether the three permanent types stay in the paper exposition.

### M32 — Sources for the selling cost and the rental size cap

The deck cites Greaney, Parkhomenko and Van Nieuwerburgh (2025) for both the 6% selling cost and the 6-room rental cap. The selling cost is a standard transaction-cost value in that literature; the rental cap is a maintained restriction whose level is contentious (earlier sensitivity work found it non-monotone and load-bearing for the ownership fit). Establish a proper empirical basis for the cap, for example the share of 6+ room units that are renter-occupied in the AHS/ACS, or reframe it as a calibrated object.

## Evidence and reproduction

- Recognized original code: `/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/code_per tommi_addingcontrolsandfixingthings.do`.
- May rooms graph: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Graphs/rooms_f_c_y_all.png`.
- May saved coefficient table: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/rooms_f_c_y_all_estimates.dta`.
- Consolidated audit: [first-birth correction review](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/README.md).
- Binned diagnostic: [results and sample flow](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/summary.csv), [sample exclusions](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/sample_flow.json).
- Source-year ownership construction: `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/Construction_Files/Code/01 Collect housing variables.do`, lines 95–118.
- Moving definitions: [PSID 2019 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM2019ER_codebook.pdf), pp. 54–56; [1984 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM1984_codebook.pdf), p. 161.

Presentation asset changes are owned by the existing September presentation
task. This data task communicated the author's choice and the distinction
between the historical graph and the pinned calibration scalar. The ledger
does not certify restoration of other historical figures whose source assets
have not yet been verified.
