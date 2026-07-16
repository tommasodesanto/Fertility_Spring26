# Review: intergenerational bequest calibration and the old-age balance sheet

Date: 2026-07-16. Reviewer: Fable (independent referee pass per
`docs/prompts/intergen_bequest_balance_sheet_fable_review_20260716.md`).
Scope: review and decision memo only; no code edits, no runs launched.

All numbers below were verified against the working tree on 2026-07-16:
the active calibration/moment code, the four evidence packets, and the PSID
construction chain (including the external PSID-SHELF `.do` files). Line
references are to the current uncommitted tree.

---

## 1. Executive verdict

- **The p90/p50 failure is primarily an asset-composition / missing-mechanism
  failure, not a taste-curvature failure.** The model's 76–84 households hold
  zero positive nonhousing wealth through the 90th percentile (PSID:
  2.37/10.34/26.27 median-income units at p50/p75/p90); the top estate decile
  is 99.9% housing (PSID: 18.6%). No value of \(\theta_1\) repairs this
  (frontier verified: dispersion reaches 3.45 only when the median collapses
  to ~2.2–2.5), and retirement-income dispersion alone does not either
  (verified: p90/p50 ≤ 2.37 at \(s_R=2\)). A nontrivial part of the empirical
  tail is also *non-comparable*: PSID nonhousing net worth includes business,
  IRA/annuity, vehicle, and other-real-estate wealth with no model
  counterpart, while DB pension and Social Security wealth are excluded from
  the numerator though their income flows sit in the denominator.
- **But the universal negative-\(b\) old age is a property of the M3 winner
  point, not of the model.** The M1 (mortality-only) winner delivered a
  positive old nonhousing median of **2.00 vs target 2.23** at ages 65–75
  while fitting the 12 shared non-estate moments at a summed loss of **3.84
  versus 16.92 at M3** (4.4× worse). Chasing the unreachable dispersion
  target — 139.93 of the 166.65 loss, i.e. 84% — is what pushed the
  calibration into the all-housing, high-debt regime and flipped the young
  renter liquid-wealth moment to the wrong sign (−0.183 vs +0.179).
- **The accounting is coherent; there is no measurement bug on the model
  side.** The estate object is \(W=b+pH\) with housing gross and debt netted
  one-for-one, deliberately and documented (`solver.py:4796-4802`); the
  utility and the moment use the same \(W\); the diagnostic reproduces all
  three calibrated moments to \(10^{-10}\) from an independently coded
  quantile routine. The one real asymmetry — no \(\psi\) wedge at death versus
  \((1-\psi)\) on living sales — is an economic incentive to die in place,
  not an accounting error.
- **Keep exactly one estate moment as a hard target:** the 76–84 median
  estate/income (disciplines \(\theta_0\), weakly). **Demote p90/p50 and the
  2+-minus-1-child gap to diagnostics.** The gap target is statistically zero
  (0.101, bootstrap SE 0.563, t ≈ 0.18), is a verified marital-composition
  cancellation (−0.886 married, +0.884 nonmarried) in a model with no
  marital state, and its target-scaled Jacobian row is so inflated by the
  0.101 denominator that *removing* it raises the identifiable rank from 7 to
  11 at the 1e-2 threshold.
- **Restore composition discipline:** reinstate the old nonhousing median at
  65–75 (2.23, dropped in the M3 redesign at inert legacy weight 0.8) as a
  hard target with a freshly bootstrapped inverse-variance weight. It is the
  moment that would have caught this failure, and M1 proves it is reachable.
  Resulting system: **14 moments, 12 free parameters** (11 clean-frontier
  + \(\theta_0\)); \(\theta_1=0.25\) external (0.50 sensitivity),
  \(\theta_n=0\) external.
- **Next action: one bounded M4 refit** of that system with warm starts that
  include the exact \(\theta_0=0\) nested seed (enforcing nested-model
  dominance, per the July-15 lesson), with ex ante pass/fail criteria below.
  No model-mechanism change and no new estate-tail target until it reports.

---

## 2. Verification table

| # | Claim (source) | Verdict | Evidence |
|---|---|---|---|
| 1 | Live system has 15 targets, 14 internally estimated parameters; `tenure_choice_kappa=0` external | **Verified / qualified** | 14-moment set `candidate_replacement_bequest_internal_v1` (`calibration.py:206-220`) + rooms moment appended by the M3 driver (`run_intergen_bequest_exit_chain.py:153-159`); domain of 11+3 at `:80-95, 169-171`. Qualification: 12 (not 11) non-estate moments are retained, and many other objects are pinned by convention (\(q=1.02^4-1\), \(\delta\), \(\phi=0.80\), \(\psi=0.06\), \(\sigma=2\), \(\eta_{supply}=1.75\), \(\lambda_d=0\), SSA survival, income process, entry-wealth distribution). |
| 2 | Estate targets 6.5013 / 3.4481 / 0.1011 with bootstrap SEs 0.2320 / 0.1325 / 0.5630; weights = 1/SE² | **Verified** | Values verbatim at `calibration.py:217-219`; hard-coded weights 18.5858 / 56.9794 / 3.1548 at `:583-585` reproduce 1/SE² to every printed digit. Off-diagonal bootstrap covariance exists on disk but is unused. Person-cluster bootstrap (499 reps, seed 20260715) verified in `audit_intergen_bequest_family_size_targets.R:339-354`. |
| 3 | Model counterpart is \(W=b+pH\), gross housing, no sale wedge, household-level annual gross income denominator | **Verified, deliberate** | `add_annual_gross_estate_wealth_moments`, `solver.py:4789-4871`; docstring at `:4796-4802` states the gross-housing convention and why. Denominator is `annual_gross_income_at_state` per state; ages 76–84 map to model start-ages {78, 82} (j=15,16) and 65–75 to {66, 70, 74} (j=12–14), all retired (\(J_R=12\), age 66). |
| 4 | Strict M3 winner: loss 166.65, residual 1.22e-6, \(\theta_0=0.310,\ \theta_1=0.536,\ \theta_n=0.710\), all 14 interior | **Verified** | `report/target_fit_full.csv` sums to 166.6536757 exactly; `parameter_table_full.csv` — no parameter within 2% of its range width of a bound (that is the `near_bound` definition, `run_intergen_bequest_exit_chain.py:375`). |
| 5 | Estate fit: median 6.267 (contrib. 1.02), p90/p50 1.881 (139.93), family gap 1.770 (8.79) | **Verified** | Same CSV. p90/p50 alone is 84.0% of the total loss. |
| 6 | Jacobian: SMM-weighted rank 9/14 at 1e-2, 12/14 at 1e-3, condition 5,499; target-scaled weaker (7/11, cond 29,944); weakest direction mixes \(\theta_n\), `h_bar_0`, \(\theta_1\), `h_bar_n` | **Verified, plus a new finding** | `identification/summary.json`, `weakest_direction.csv` (loadings −0.798, −0.435, −0.303, −0.253). New: the \(\theta_1\) and \(\theta_n\) target-scaled columns are numerically near-identical (median loadings 3.5867 vs 3.5861; p90/p50 loadings 1.9566 vs 1.9565) — the two parameters are locally indistinguishable in this moment set. Leave-one-out: dropping the *family-gap* moment **improves** target-scaled conditioning (rank 7→11 at 1e-2, condition 29,944→5,615) because its 0.101 target inflates the scaled row; dropping either 76–84 moment worsens conditioning. |
| 7 | Reachability frontier: no cell matches median + gap within 1 SE; \(\theta_1=1.2\): median 6.278, p90/p50 1.739, gap 1.448; \(\theta_1=8/16\): dispersion 3.57/3.95 but median 2.51/2.19 | **Verified** | `reachability_frontier.csv`; `admissible_1se_cells: 0`. Cells 7 (\(\theta_1=0.8\)) and 11 (\(\theta_1=8\)) are non-strict, honestly excluded. Additional reading: across all cells with \(\theta_1\le 4\), the median stays in [6.19, 6.28] while the re-optimized \(\theta_0\) ranges 0.10–1.31 — the median is nearly flat in the whole bequest block given the M3 housing coordinates. |
| 8 | Retirement-income diagnostic: best strict cell \(s_R=0\); no cell within 1 SE on all three; \(s_R=2\): p90/p50 2.371, median 5.505, \(\theta_0=0\), non-strict | **Verified** | `retirement_income_dispersion_frontier.csv`; non-strict cells {5,6,7}; \(\theta_0\) at its zero bound in cell 7. |
| 9 | Decomposition table (PSID 6.867/17.101/36.335 vs model 6.267/8.866/11.789; nonhousing 2.368/10.336/26.272 vs 0/0/0; correlations 0.393 vs 0; top-decile housing share 0.186 vs 0.999; raw p90/p50 5.291 vs 1.881; model b −6.44/−5.14/−3.84; pH 11.97/15.96/19.95) | **Verified bit-level** | `headline_comparison.csv`; three calibrated moments reproduced to better than 1e-10 by an independently coded quantile routine (`diagnose_intergen_bequest_distribution.py:329-331`). Model "correlation 0" is a 0/0 float artifact of degenerate pension income (1.6e-16; the asinh variant correctly reports NaN). |
| 10 | Diagnostic component mapping (owner equity \(=pH+\min\{b,0\}\), owner nonhousing \(=\max\{b,0\}\), renter \(b\)) is dictated by the model | **Qualified: it is a convention, but the only sensible one** | The model has a single net asset state; `kernels.py:930-932` confirms "b already contains secured mortgage debt." Any split is a diagnostic overlay (`diagnose_intergen_bequest_distribution.py:154-165`); raw \(b\) and \(pH\) are saved alongside, so it is auditable. |
| 11 | Old-age model income degenerate | **Verified, intended restriction** | `retirement_income_z_scale=0.0` default (`parameters.py:138`, comment says diagnostic-only); scalar pension \(=0.4296\) annual = the normalizing median (`parameters.py:542-557`). Not a lost state (z keeps mixing, payoff-irrelevant) and not an aggregation error. |
| 12 | Old-age borrowing: unsecured capacity 0 after 62; owner floor \(-\phi p H\) at every retired age; no amortization | **Verified** | Taper `parameters.py:435-488` (s=0, D=0 for age ≥ 62); owner floor `solver.py:118-133`, `2529-2534`; tenure-staying transition re-tests nothing (`solver.py:3047-3048`). Old renters cannot carry \(b<0\); all old debt is collateralized owner debt. |
| 13 | "The model's age-82 state forces terminal liquidation" (July-15 status/notes) | **Rejected for the production/M3 spec** | Terminal continuation is `Vnr = Vbq` with tenure-indexed gross housing (`solver.py:2189-2190, 2155-2162`) — owners die in place, no forced sale, no wedge. The forced-liquidation language is accurate only for the rejected owner-LTV-taper battery arms (A2–A5), where the tightened floor can make selling the only feasible branch. Status prose should be corrected. |
| 14 | Bequest utility matches \(B(W,n)=\theta_0\max\{1+\theta_n n,0\}\frac{(\theta_1+W)^{1-\sigma}-\theta_1^{1-\sigma}}{1-\sigma}\) | **Verified for M3, with a latent-default trap** | `bequest_utility_vec` (`solver.py:6013-6047`) implements scale, curvature, and — when `normalize_bequest_utility=True` — the \(-\theta_1^{1-\sigma}/(1-\sigma)\) subtraction exactly; the M3/M-arm chain overrides normalization on (`run_intergen_bequest_exit_chain.py:191,203`), so the estimated spec matches the intended formula. Trap: the `parameters.py:67` default is `False`, and the un-normalized level shift scales with \(1+\theta_n n\) (a hidden child utility penalty at low estates) — any future driver that forgets the override changes the economics silently. \(W\) is clamped at zero inside the utility; the measured moment does not clamp (immaterial at old ages: owner estates are bounded below by \((1-\phi)pH>0\) and old renters have \(b\ge 0\)). |
| 15 | Family-gap target is a marital-composition cancellation; model has no marital state | **Verified** | PSID 65–75 median gaps: −0.886 (married), +0.884 (nonmarried), aggregate +0.101 (`psid_marital_composition.csv`, recomputed exactly); `summary.json: model_has_marital_state: false`. |
| 16 | PSID estate = `NETWORTHR` = `NETWORTH2R` + `HOMEEQUITYR`; denominator `INCFAMR` | **Verified, with comparability caveats** | Additivity holds to rounding on 16,923 triples (`definition_consistency.csv`). Nonhousing includes business/farm, stocks + IRA/annuity balances, vehicles, other real estate, net of non-mortgage debt; DB pension and Social Security *wealth* are absent from the 28-component list while their income flows are (per standard PSID documentation) in `INCFAMR`; income refers to the prior tax year, wealth to the interview date; raw-response top-coding compresses the tail, so 3.448 understates true dispersion if anything. |
| 17 | Timing of measured distribution | **Verified, one edge case** | Moments use beginning-of-period \(b\) conditioned on current-period tenure (`stats.wealth_moment_timing`, `solver.py:4417`) — a survivor cross-section, which is the right counterpart to PSID living households (both "estate" labels are misnomers in the same direction, so no asymmetry). Edge case: a household transacting this period contributes pre-transaction \(b\) with post-transaction tenure; second-order at old ages with little churn. |
| 18 | Numerics create the negative-\(b\) old distribution | **Rejected as primary cause** | The b-grid spans [−12, 30] with an interior floor: at p=0.8569 collateral floors are −1.37 (H=2) to −6.86 (H=10) in model units, all far from the grid edge. Quantiles snap to grid nodes (reported p50 −2.7674 is the node just below the H=4 floor −2.7421), and the H≥8 floors sit in the coarse 0.7-spacing segment, so quantile *values* carry ±(one node) noise — but the sign and scale of old debt are economics, not discretization. `report_clamp_hits` is a dead switch; no boundary-mass diagnostic exists. |

---

## 3. Economic and numerical diagnosis

**The causal chain producing the M3 old-age balance sheet.** After age 62 the
unsecured line is zero and the only borrowing is collateralized owner debt
with the *constant* floor \(b \ge -\phi p H\), \(\phi=0.80\), never amortized
(§2, rows 12–13). Returns are symmetric at \(R=1.02^4=1.0824\) per period.
At the winner \(\beta=0.9866^4=0.9475\), so \(\beta R \approx 1.026\), but
with SSA survival the effective weight \(\beta s_j R\) falls from 0.963 at 66
to 0.851 at 78 — retired households optimally decumulate. The pension is a
scalar (0.4296 annual, identically the normalizing median), so there is no
income risk, no expense risk, and hence *no precautionary motive of any kind
after 66*. The only reason to hold wealth into death is the warm glow, and
the warm glow is cheapest to satisfy with housing: the estate values housing
gross (\(W=b'+pH\), no \(\psi\)), a living sale nets only \((1-\psi)pH\)
(dying in place beats sell-then-die by \(\psi pH \approx 0.7\) annual incomes
at H=6), and the house pays service flow while waiting. Optimal behavior is
therefore: keep the house, run \(b\) down toward (not necessarily to) the
collateral floor, die with home equity as the entire estate. Every household
does this, so nonhousing wealth is identically zero through p90 and the
estate support collapses to the few owner rungs {6, 8, 10} net of debt —
p90/p50 of 1.88 with a hard ceiling near 2 while the pension denominator is
common. That is the composition failure, and it is fully consistent with the
model's constraint set; nothing is broken numerically.

**Why this is a point property and not a model property.** The M1
mortality-only winner — same model, same constraints, \(\theta_0=0\),
\(\theta_1=0.25\), \(\theta_n=0\) external — delivered an old (65–75)
nonhousing *median* of +2.00 against the 2.23 target
(`intergen_mortality_recalibration_20260715/report/target_fit_full.csv`),
with the young renter liquid-wealth moment nearly exact (+0.190 vs +0.179).
Positive late-life liquid saving is attainable in this model class. The M3
optimizer abandoned it because the p90/p50 target, at inverse-bootstrap
weight 56.98, offered ~140 loss points that no feasible point could collect;
the search wandered on that plateau and paid for it with the established
block: summed loss on the 12 shared moments 3.84 (M1) → 16.92 (M3), TFR
2.01→2.20, ownership 0.607→0.713, young liquid wealth +0.19→−0.18 (wrong
sign), childless renter rooms 3.83→4.76. The three estate rows also dominate
every column of the SMM-weighted Jacobian (loadings ~100/~50 on nearly every
parameter), which is why the 14-parameter system shows rank 9/14: the
weighting collapsed the effective moment space onto the estate block.

**What the frontier really rules out.** Holding M3's 11 non-bequest
coordinates fixed, the 76–84 median estate is *flat* in the entire bequest
block (6.19–6.28 for \(\theta_1\le4\) while \(\theta_0\) ranges 0.10–1.31)
— the estate level is set by the housing block (\(\chi\), \(H_0\),
\(\bar h\)-parameters), not by tastes over estates. Dispersion becomes
reachable only at \(\theta_1=8\)–16, where the near-linear warm glow makes
most households stop holding housing to death and the median collapses. So
the frontier rules out *any* repair through \((\theta_0,\theta_1,\theta_n)\)
at this equilibrium; it does not rule out bequest motives, and it also warns
that \(\theta_0\)'s identification off the median alone is weak.

**The 65–75 window makes the profile failure visible too.** Model estate
medians at 65–75 are 11.97 (one child) and 13.74 (2+) versus PSID 4.80 and
4.90 — the model is ~2.7× too wealthy at 66–74 and then decumulates to 6.27
by 78–82, while the PSID profile is flat-to-rising (survivorship included).
The family-gap miss (+1.77 vs +0.10) is therefore *not* evidence about
\(\theta_n\) alone: with \(\theta_n=0.71\) the 2+ scale is
\(1+0.71\cdot 2 \approx 2.4\), and families with more children also carry
larger houses through the \(\bar h_n\) requirements and the \(\psi\)-induced
tenure hysteresis, so the whole 2+ distribution shifts up mechanically. In
the data the aggregate gap is a knife-edge cancellation across marital
strata that the model cannot represent — and at t ≈ 0.18 it is statistically
zero anyway.

**Measurement verdict.** On the model side, no accounting bug: utility and
moment use the same \(W\); the diagnostic's split is a disclosed convention
forced by the single-asset state; timing is survivor-stock vs survivor-stock.
On the data side, the *level* target (median) is clean, but the *tail*
targets embed wealth categories the model does not own (business equity,
IRA/annuity balances, other real estate — dominant in the PSID top decile,
where housing is only 18.6% of wealth) and an annuitization asymmetry (DB/SS
wealth excluded from the numerator, its income in the denominator). Matching
3.448 exactly was never a well-posed demand on this model; the moment can
return as a hard target only after (i) a mechanism exists that can move it
and (ii) its empirical counterpart is restated on a model-consistent wealth
concept.

**Answers to the eight standing questions.**
1. *Taste-curvature vs composition:* composition, with a measurement wedge on
   top. Separating evidence: frontier (no \(\theta_1\) works), dispersion
   diagnostic (income heterogeneity insufficient), decomposition (nonhousing
   0/0/0 vs 2.4/10.3/26.3; top-decile housing share 0.999 vs 0.186), and raw
   p90/p50 (5.29 data vs 1.88 model — the tail is missing in levels, not
   created by the denominator).
2. *Why \(b<0\) through p90 at 76–84:* zero unsecured capacity after 62 plus
   a constant \(-\phi pH\) collateral floor with no amortization, symmetric
   2% returns, \(\beta s_j R<1\), no late-life risk, and a warm glow that is
   cheapest to satisfy with gross-valued housing. Economically intended in
   the code; empirically indefensible (most US 76–84 owners are
   mortgage-free) and, at the M3 point, an artifact of tail-chasing — M1 did
   not look like this.
3. *Terminal coherence:* yes — debt netted one-for-one, clamp at zero
   (limited liability at death, estate never negative for owners since
   \(W\ge(1-\phi)pH\)), no estate recycling (dying mass vanishes; acceptable
   for calibration, worth remembering for welfare statements). The only
   asymmetry is the missing \(\psi\) at death, which is a deliberate
   modeling choice that subsidizes dying in place.
4. *p90/p50 as hard target:* no. Demote to diagnostic until the balance
   sheet has a mechanism that can move it. It currently identifies nothing
   (\(\theta_1\)'s column is collinear with \(\theta_n\)'s) and corrupts
   everything else through its weight.
5. *Family gap as hard target:* no. Statistically zero, marital-composition
   confounded, mechanically contaminated by housing hysteresis, and its
   target-scaled Jacobian row actively damages the conditioning. \(\theta_n\)
   should be externally restricted to 0 (consistent with the July-8 and
   July-13 verdicts that \(\theta_n\) is null/unidentified); no credible
   replacement moment exists without a marital state.
6. *Median as the single bequest target:* yes, with eyes open — it is also
   (mainly) a housing-block moment, so \(\theta_0\)'s practical
   identification rests on the nested \(\theta_0=0\) comparison. If the
   refit sends \(\theta_0\) to zero, report the zero envelope honestly
   rather than forcing positivity; that is a finding about survivor-wealth
   moments in this model class, not proof that bequest motives are absent.
7. *Existing moments with retirement-saving/leverage information:* exactly
   one — young childless renter liquid wealth 25–35 (currently wrong-signed
   at M3, nearly exact at M1). The redesign dropped the two old-age
   composition moments; that is precisely how the model matched an estate
   median with the wrong balance sheet without paying any loss for it.
8. *Minimal defensible system:* below.

---

## 4. Minimal calibration recipe

| Moment | Status now | Recommendation | Parameter it disciplines | Replacement / restriction |
|---|---|---|---|---|
| Median total estate/income, 76–84 (6.5013, SE 0.2320) | hard | **Keep hard** (weight 1/SE²) | \(\theta_0\) (weakly; also housing block) | — |
| p90/p50 total estate/income, 76–84 (3.4481, SE 0.1325) | hard | **Demote to diagnostic** | was \(\theta_1\) | \(\theta_1\) external at the Nakajima–Telyukova (2017, JF) anchor: shifter \(\zeta\approx\$19{,}600\approx\) a quarter of average annual household income \(\Rightarrow \theta_1\approx0.25\) in model units (July-14 spec memo anchor; NT numbers tagged [E] there — verify against the paper before wiring). Sensitivity cell at a De Nardi (2004)-scale shifter (several times average income; exact value verified from the paper), the regime where the frontier's dispersion appears. Moment returns as a hard target only with a tail mechanism *and* a model-consistent empirical counterpart. |
| 2+ minus 1-child median estate/income gap, 65–75 (0.1011, SE 0.5630) | hard | **Demote to diagnostic** | was \(\theta_n\) | \(\theta_n=0\) external. No replacement moment exists absent a marital state; the alternative is modeling marriage, which this gap alone does not justify. |
| Old nonhousing (net liquid) median/income, 65–75 (2.2305) | dropped in M3 redesign | **Reinstate as hard target** with fresh person-bootstrap SE and 1/SE² weight | late-life portfolio composition (loads on \(\beta\), \(\alpha_{cons}\), housing block; gates against the all-housing regime) | its person-bootstrap SE must be computed first via the existing R pipeline (minutes of work; same 499-rep design). |
| 12 established non-estate moments | hard | keep hard, unchanged weights | as before | — |
| Old-age ownership level (0.7643) | dropped (July-15 identification audit) | stays a reported validation moment | — | — |
| Full decomposition packet (nonhousing/equity quantiles at 76–84, top-decile housing share, estate–income correlation, 65–75 profile) | new | standing diagnostic in every calibration readout | — | — |

Count: **14 hard moments, 12 free parameters** (11 clean-frontier +
\(\theta_0\)); externally restricted: \(\theta_1=0.25\), \(\theta_n=0\),
`tenure_choice_kappa` \(=0\), plus the standing pinned conventions. The
system is counted-identified with two wealth moments (one level, one
composition) against one free wealth-preference parameter; conditioning
should improve materially since both pathological rows (the weight-57
dispersion row and the scale-inflated gap row) leave the objective.

## 5. Parameter recommendation

- \(\theta_0\): **keep internally estimated** against the estate median, but
  enforce nested-model dominance: inject the exact \(\theta_0=0\) M4-nested
  seed into every chain and fail the collector if the free search ends above
  it (the July-15 A3 failure must not recur). Accept a zero outcome as a
  result, not a defeat.
- \(\theta_1\): **external restriction at the verified Nakajima–Telyukova
  (2017, JF) anchor** (\(\gamma=0.43\), \(\zeta\approx\$19{,}600\approx\) a
  quarter of average annual household income \(\Rightarrow\theta_1\approx
  0.25\) with the model's average-income normalization of 1.0), with a
  De Nardi (2004)-scale sensitivity cell (several times average income —
  verify \(\phi_2\) and units from the paper). Provenance: this is the
  July-14 spec memo's own anchor (`bequest_specification_memo_20260714.tex`,
  §(d) and the Option-B1 discussion), which also stated that \(\theta_1\)
  "must remain an explicitly external scale until a model-specific estate
  moment is added." M3 was that attempt; the reachability frontier and the
  \(\theta_1\)/\(\theta_n\) Jacobian collinearity show it failed, so the
  external anchor stands. NT identify \(\zeta\) from variation this model
  lacks (medical risk, downsizing, reverse mortgages) — internal estimation
  here is noise-fitting, not identification. State plainly in any writeup
  that \(\theta_1\) is externally restricted.
- \(\theta_n\): **external restriction at 0**, per the child-blind
  convention of the structural literature (De Nardi 2004; Lockwood 2018;
  Nakajima–Telyukova 2017; De Nardi–French–Jones–McGee 2025, where children
  are not even a conditioning variable; Kværner 2023 models the motive
  child-blind), the classic Hurd (1989, Econometrica) null, Kopczuk–Lupton
  (2007, REStud — motive *incidence*, not scale), and our own PSID gap of
  0.101 (SE 0.563), a within-sample replication of that null. If the child
  margin of bequests becomes a research question, it needs a marital state
  and a within-stratum moment, not this aggregate.
- Housekeeping: the chain driver already sets
  `normalize_bequest_utility=True`, but the `parameters.py` default is
  `False`; flip the default (or assert it in the spec) so no future driver
  silently reintroduces the un-normalized level term, which scales with
  \(1+\theta_n n\) and functions as a hidden child utility penalty
  (§2 row 14).

## 6. Next action: one bounded experiment (M4 refit)

**Design.** New arm M4 in the existing chain driver: 14 moments (12
established + estate median 76–84 + reinstated nonhousing median 65–75 with
fresh bootstrap weight), 12 free parameters (11 clean-frontier +
\(\theta_0\)), externals \(\theta_1=0.25\), \(\theta_n=0\),
`tenure_choice_kappa=0`, SSA survival on (M1 convention). Warm starts: M1
winner, M3 winner (mapped), and the exact \(\theta_0=0\) nested seed; 4–6
chains, search evaluator \((10,10^{-4})\), two tight \((40, 2.5\times10^{-5})\)
winner repeats, strict collector, ≤ 15 CPU-hours, per-case checkpoints,
Torch smoke of the exact loop first — all per the standing contract
conventions. Precondition (data, minutes): extend
`audit_intergen_bequest_family_size_targets.R` to bootstrap the 65–75
nonhousing median SE with the identical person-cluster design.

Second precondition (literature verification, half a day): open Nakajima–
Telyukova (2017) and De Nardi (2004) and verify \(\gamma\), \(\zeta\), and
\(\phi_2\) with their exact units from the papers (the July-14 memo tags the
NT numbers [E], extracted but not re-checked); convert through the model's
average-income normalization (the pension formula hardcodes mean worker
income \(=1.0\); confirm against the solved profile) and wire the two
verified values as the \(\theta_1\) baseline and sensitivity cells.

**Pass (all required):**
1. Strict, exactly repeated tight winner (bit-identical two solves).
2. Estate median within 1 bootstrap SE (|gap| ≤ 0.232).
3. Nonhousing median within 2 of its fresh bootstrap SEs, and strictly
   positive at the winner.
4. Summed loss on the 12 established moments ≤ 1.15 × the M1 winner's 3.84
   (i.e. ≤ 4.42) — the estate block must no longer tax the established fit.
5. Free-\(\theta_0\) winner weakly below the injected \(\theta_0=0\) seed.

**Fail / stop conditions.** If (2) and (3) cannot hold jointly with all 12
coordinates free, the level and composition moments genuinely conflict at
this model structure — that is new information and the answer becomes "no
further calibration until a late-life mechanism is added," with the
mechanism decision escalated (leading candidate: out-of-pocket
medical/LTC expense risk after retirement, \(m_j(z,\eta)\) entering the
budget as \(c + m_j \le Rb + y_j - \dots\) with the \(m\)-process calibrated
*externally* from HRS/MEPS age-profiles à la De Nardi–French–Jones (2010),
so it consumes no SMM moments; second candidate: annuitized-vs-liquid
retirement wealth heterogeneity, external from HRS/SCF annuitization
shares). Diagnostics (p90/p50, family gap, decomposition packet) are
reported either way but decide nothing here.

**What this run decides.** Whether the minimal target redesign restores a
single calibration that (i) keeps the established block at M1 quality,
(ii) matches late-life wealth level *and* sign-correct composition, and
(iii) says whether the estate median lends any support to \(\theta_0>0\)
once the unreachable rows stop paying. It does not — and is not meant to —
produce the PSID tail.

## 7. Open risks

- **\(\theta_0\) may be practically unidentified even in the reduced
  system** (the frontier shows the median nearly flat in the bequest block
  at fixed housing coordinates). The nested-zero profile is the honest
  instrument; a zero envelope would leave the paper's bequest channel
  externally parameterized, which referees will accept only if stated
  plainly. This is the classic bequest-vs-precautionary separation problem
  (survivor wealth levels cannot distinguish them without expense risk or
  end-of-life data).
- **The old-age leverage counterfactual remains in the model** even after
  the redesign: nothing amortizes owner debt, so any policy experiment
  operating through late-life home-equity margins (reverse-mortgage-like
  capacity, old-age property-tax changes) is unreliable until either an
  amortization/LTV-age schedule or a liquid-savings motive exists. The
  rejected owner-LTV taper battery showed crude schedules break the
  ownership path; a mechanism, not a constraint patch, is the likely fix.
- **The dying mass vanishes** (no estate recycling to entrants). Fine for
  moment-matching; not fine for any welfare or intergenerational-transfer
  statement built on this block later.
- **Tail non-comparability is unresolved by design:** the model cannot and
  should not chase PSID's business/IRA-driven top decile. If upper-tail
  wealth ever becomes substantive for the housing–fertility question, that
  is a different model (earnings tail or entrepreneurial asset à la
  De Nardi 2004 / Castañeda et al. 2003), and the target should be restated
  on a housing+liquid concept first.
- **Numerical hygiene, second order:** estate quantiles snap to b-grid nodes
  (±1 node), the H≥8 collateral floors sit in the coarse 0.7-spacing grid
  segment, `report_clamp_hits` is a dead switch, and the golden-section
  ±2-unit search window near the floors is a heuristic. None of this drives
  the diagnosis; all of it belongs on the list before any future run leans
  on old-age wealth quantiles at Nb=120 (verify at 240 per convention).
- **Status-file correction needed:** the "age-82 forced terminal
  liquidation" line in the July-15 status/notes is wrong for the production
  and M3 specifications (it described the rejected taper arms); it should
  not propagate into future ownership-path acceptance rules.

---

## Appendix A. Full 15-moment target fit at the strict M3 winner (loss 166.653676, residual 1.22e-6)

| moment | target | model | gap | weight | loss contribution |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 2.20285 | +0.28485 | 20.0 | 1.6228 |
| childless_rate | 0.188 | 0.16672 | −0.02128 | 20.0 | 0.0091 |
| own_rate | 0.57547 | 0.71325 | +0.13777 | 100.0 | 1.8982 |
| own_family_gap | 0.16766 | 0.09396 | −0.07370 | 45.0 | 0.2445 |
| housing_increment_0to1 | 0.66443 | 0.42378 | −0.24066 | 14.0 | 0.8108 |
| prime30_55_childless_renter_mean_rooms | 3.80529 | 4.76457 | +0.95928 | 6.0 | 5.5213 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.59613 | 0.89666 | +0.30053 | 25.0 | 2.2580 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.17923 | −0.18296 | −0.36219 | 12.0 | 1.5742 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.41876 | 2.09249 | −0.32627 | 12.0 | 1.2774 |
| own_rate_2534 | 0.34117 | 0.33206 | −0.00911 | 80.0 | 0.0066 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36770 | 0.27258 | −0.09512 | 8.0 | 0.0724 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.50132 | 6.26738 | −0.23394 | 18.5858 | 1.0171 |
| old_total_estate_wealth_to_annual_income_p90_p50_7684 | 3.44811 | 1.88101 | −1.56710 | 56.9794 | 139.9301 |
| old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575 | 0.10111 | 1.76985 | +1.66874 | 3.1548 | 8.7850 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997 | 6.30057 | +0.52060 | 6.0 | 1.6262 |

## Appendix B. Full parameter table at the strict M3 winner

All 14 estimated parameters are interior (none within 2% of its range width
of a bound). Bounds are the audited relaxed discovery domain, not the
narrower production box; note the domain's `beta_annual` floor of 0.80 sits
below the standing external restriction \(\beta_{annual}\ge 0.94\) (the
estimate, 0.9866, satisfies it).

| parameter | estimate | bounds | transform | role |
|---|---:|---|---|---|
| beta_annual | 0.98662 | [0.8, 0.9995] | discount | estimated |
| alpha_cons | 0.54857 | [0.02, 0.98] | logit | estimated |
| c_bar_0 | 0.96666 | [0.0, 2.0] | softzero | estimated |
| c_bar_n | 0.30589 | [0.0, 3.0] | softzero | estimated |
| h_bar_0 | 0.35206 | [0.05, 5.8] | log | estimated |
| h_bar_jump | 1.83102 | [0.0, 8.0] | softzero | estimated |
| h_bar_n | 0.22525 | [0.0, 5.0] | softzero | estimated |
| psi_child | −0.03693 | [−3.0, 3.0] | asinh | estimated |
| kappa_fert | 1.30727 | [0.02, 50.0] | log | estimated |
| chi | 1.13098 | [0.1, 5.0] | log | estimated |
| H0 | 7.76407 | [0.2, 80.0] | log | estimated |
| theta0 | 0.31025 | [0.0, 1.5] | softzero | estimated |
| theta1 | 0.53614 | [0.1, 2.0] | log | estimated |
| theta_n | 0.71008 | [0.0, 1.5] | softzero | estimated |
| tenure_choice_kappa | 0.0 | — | external restriction | fixed |

Other pinned objects (not searched): \(q=1.02^4-1\),
\(\delta=1-(1-0.011)^4\), \(\phi=0.80\), \(\psi=0.06\), \(\sigma=2\),
\(\tau_{pay}=0.179\), \(\eta_{supply}=1.75\), \(\lambda_d=0\) with the
42–62 unsecured taper, SSA-2023 post-retirement survival, the matched
annual income process (\(\rho=0.9602\), \(\sigma_\epsilon=0.0645\)), the
PSID entry-wealth ratio distribution, `H_own = [2,4,6,8,10]`,
`hR_max = 6.0`, and `retirement_income_z_scale = 0`.

## Appendix C. M1 vs M3 on the 12 shared non-estate moments (loss contributions)

| moment | M1 (loss 6.860 total) | M3 (loss 166.654 total) |
|---|---:|---:|
| tfr | 0.1811 | 1.6228 |
| childless_rate | 0.0003 | 0.0091 |
| own_rate | 0.1002 | 1.8982 |
| own_family_gap | 0.0073 | 0.2445 |
| housing_increment_0to1 | 0.0000 | 0.8108 |
| prime30_55_childless_renter_mean_rooms | 0.0047 | 5.5213 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.8829 | 2.2580 |
| young_childless_renter_liquid_wealth (2535) | 0.0013 | 1.5742 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 0.3780 | 1.2774 |
| own_rate_2534 | 0.9577 | 0.0066 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.0000 | 0.0724 |
| aggregate_mean_occupied_rooms_18_85 | 1.3278 | 1.6262 |
| **sum** | **3.841** | **16.922** |

M1 also delivered `old_nonhousing_wealth_to_income_median_6575` = 2.0029 vs
target 2.2305 and young liquid wealth +0.1898 vs +0.1792 — positive
late-life and young liquid balances are reachable in this model at a point
that fits the established block.
