# Mock presentation feedback

Running record of the September presentation mock: slide fixes, substantive questions, and decisions.

First-pass review task: **Mock presentation: concise conceptual review** (`01a09818-86af-75a1-9ff4-3a396e394c4d`), max reasoning. Replies are advisory; only individually authorized fixes may be implemented.

## Open issues

### M02 — Adult aging and death

Do people die randomly and age randomly at any age? Clarify the actual timing.

**Status:** open; concise first-pass review requested. No changes authorized.

### M03 — Child versus adult aging

Children mature stochastically; do adults age deterministically?

**Status:** resolved for author-facing wording. The concern is preserved for review: the stochastic event is the exit from dependency, while adult household age advances deterministically. No implementation verification is claimed.

### M04 — Number of children: m versus n

Do we still need the distinction between children at home (m) and children ever born (n), now that bequests do not depend on n?

**Status:** notation consolidation discussed and deferred; keep children ever born $n$ distinct from dependents at home $m$.

### M05 — Utility specification and literature

The utility function is broadly unjustified and hard to parse. Is it common in the literature and quantitative work? Explain the role of every component: why divide by e(m) outside; why retain the Stone–Geary housing floor; why combine Cobb–Douglas and CRRA; why is the preference for children linear and outside the consumption/housing aggregator? Saverio suggests showing only a generic u and its shape properties in the model section, then specifying functional forms in Quantification.

**Status:** exposition edit authorized and completed. The generic model utility and exact quantitative functional form are now separated; justification and literature review remain open.

### M06 — Budget-constraint slide

The slide is messy; the constraints need to be made clear.

**Status:** timing checked; slide simplification explicitly deferred. No budget-slide or model edits made.

**Timing finding:** households trade housing out of beginning-of-period liquid wealth, then earn/pay the gross bond return on the resulting liquid position. Current earnings/pensions and rebates enter afterward; consumption, rent, maintenance, and property tax are period flows. The ordinary initial down-payment test uses wealth available before current earnings. These conventions support the budget slide's gross-return factors on housing purchases and sale proceeds.

The transition rental-pricing identity is consistent with this convention: current rent plus next-period asset value, less current maintenance and property tax, equals the gross bond return on the current purchase price. This verifies timing, not the economic justification of frictionless rental pricing for constrained owners with transaction costs (M10).

**Evidence:** read-only inspection and source-extracted checks of the optimized solver, tenure kernel, and perfect-foresight rent function; nine interior tenure alternatives, the pre-income down-payment restriction, the dated rental identity, and constant-price nesting passed. The quantitative task independently confirmed that all four source hashes match the fitted-patch terminal contract and stationary-policy preparation. No model solve was run. Grid-clipping incidence and the estate price-date convention were not audited here. Scratch check: `tmp/september_slides_review/budget_timing_check.py`.

### M07 — Fertility and child-aging slide

The slide is messy and confusing. Some material concerns calibration and should not be in the model exposition.

**Status:** event-order clarification checked; slide cleanup and any model change remain unapproved.

**Timing finding:** the household chooses whether to attempt a birth; success raises both children ever born and current dependents. Housing, consumption, and current utility use this post-birth child count. Maturation then determines next period's dependents. Under independent-count maturation, conditional on the parent's survival, $m_{t+1}\sim\operatorname{Binomial}(m_t+d_t,1-\mu)$. The newborn is included in this first maturation draw: there is no minimum childhood duration. With a four-year period, a child can therefore exit dependency by the next model date after birth. This is a substantive approximation to flag, not evidence of a discrepancy between the household problem and forward transition. No change implemented; parental-death accounting remains M13.

**Evidence:** pinned optimized `solver.py` uses post-birth values at lines 2770–2821, current child-count utility adjustments at 2241–2271, and child aging in continuation values at 2521 and 7056–7079. `parameters.py:896` builds the independent binomial transition. The perfect-foresight population path calls `run_e5f_open_population_transition.py:828`, which advances the post-fertility survivor distribution using that same child transition. Inspection only; no numerical solve.

### M08 — Age earnings profile

Is e_a exogenously fixed, and is it really needed?

**Status:** open; concise first-pass review requested. No changes authorized.

### M09 — Housing-supply timing

Is housing supply instantaneous? Confirm.

**Status:** open; concise first-pass review requested. No changes authorized.

### M10 — User-cost equation

Clarify the user-cost equation very deeply, for Tommaso as well as the audience.

**Status:** open; concise first-pass review requested. No changes authorized.

### M11 — Gumbel shocks, expectations, and constraints

If a household is constrained, it still receives the taste shock: is that appropriate? Clarify the Gumbel shocks and expectations.

**Status:** open; concise first-pass review requested. No changes authorized.

### M12 — Household-problem slides

The slides are messy, the notation is not transparent, and it is unclear whether m or n is needed. The preceding introductory slide may be unnecessary.

**Status:** open; concise first-pass review requested. No changes authorized.

### M13 — Parental death and children

Very important: how do survival s_a and the number of children interact? If parents die, do their children mature? If so, is market clearing/population accounting consistent?

**Status:** substantive accounting gap flagged after authorized read-only investigation. No model or slide changes authorized for this issue.

**Measured size, initial stationary state:** using checkpoint `120ffc45...` (the recovered no-rebate history's initial equilibrium), parent-household exit removes 0.005972242 literal dependents per unit household mass every four years: **1.19965% of current dependents**, or **5.18089% of births**. All loss is at parent ages 66+, and 48.9851% comes from the forced final-age exit. These are flows, not a measured unassigned-child stock; the 2023 measurement is reported below. Full age table, source hash, reproduction script and checks: `output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/child_accounting/README.md`.

**Author decision:** retain stochastic maturation to avoid tracking child birth cohorts. Fixed or bounded dependency durations are set aside. The quantitative task has been asked to recommend an accounting repair that preserves stochastic maturation, explains care/housing needs after parental death, and avoids person/entrant loss or double counting. Recommendation only; no implementation authorized. This task remains slides-only.

**Quantitative task recommendation, awaiting author choice:** dependents of exiting households enter an aggregate publicly financed care pool and retain stochastic maturation. Specify their consumption and housing resources and balance funding, potentially before rebating property-tax revenue. Surviving households' children-ever-born state remains unchanged. Estates follow the existing allocation unless an explicit care-fund claim is chosen; do not allocate assets twice. Parent death neither destroys child persons nor itself creates adult entrants. Pool maturation must be reconciled with the separate adult-entry register without adding a second entrant stream. A scalar pool addresses missing care demand but does not by itself close the existing dependent/person-stock bridge. Care costs, funding, and the entry reconciliation remain design decisions; no implementation approved.

**2023 measurement from the quantitative task:** recovered-history current literal dependents 0.525204862; births 0.093827175; parent-exit dependent loss 0.005310226 per four-year transition. Loss is **1.011077% of dependents** and **5.659582% of births**; 36.787057% comes from forced terminal-age exit, and all losses occur at parent ages 66+. Source: `recovered_sequence/source/readout_2023/dependent_loss_summary.json` within the same overnight result folder. This is a flow, not the untracked pool's stock. Arrays were obtained during the task's separately authorized fixed-coordinate observer replay; no extra solve for M13.

**Alternative considered and set aside:** bounded child dependency would address the measured source of the gap. Last birth is at model age 42; a maximum 20-year dependency duration ends before parental mortality begins in this specification. A 50/50 departure at 16 or 20 years illustrates how to preserve the current 18-year mean, while removing immediate maturation and the unbounded tail. This would need child-age/cohort information and renewed fit assessment. If younger-age parent mortality were introduced, explicit guardian/care assignment would still be necessary.

$s_a$ is the parent household's survival probability to the next age cell. The household transition multiplies the entire post-birth distribution by $s_a$ before maturation. Dependents attached to dead households therefore leave the household distribution; they are neither reassigned nor included in that transition's reported mature-child flow. At the terminal parent age, the entire household cohort exits.

This does **not** imply that the fitted transition deletes their future adult entry: the separate recorded-birth queue continues independently of parental survival. The newer person-demography operator instead advances child persons using their own age/sex survival and determines household-head totals from headship rates. Neither mechanism, in the inspected implementation, explicitly places surviving children of dead parents into another household or reconciles those children's ongoing housing needs with household dependent states. Thus aggregate demographic/household-head identities do not certify that dependent-child accounting link. Equilibrium and welfare effects of a repair have not been measured.

**Evidence:** pinned `run_e5f_open_population_transition.py:847–895` (survivor-only dependent transition and maturation; terminal exit); `run_e5f_perfect_foresight_transition.py:650–694` and `advance_birth_vintage_queue` at `run_e5f_open_population_transition.py:351` (separate future-entry queue); `run_e5f_perfect_foresight_person_demography.py:807–838` and `demographic_transition/household_person_coupling.py:33–99` (separate person law and age-head-mass reconciliation). The fitted-patch receipt retains the historical head-age bridge and frozen 2023 person anchor; do not represent it as the newer person-cohort operator. The quantitative task independently confirmed the gap: separate person and head accounting gates pass, but no invariant reconciles household-dependent children with demographic children, and there is no explicit orphan reassignment. No solve or implementation change.

### M14 — Equilibrium and population accounting

The within-period equilibrium detail is unclear and probably unnecessary. Give a general equilibrium definition, presumably including the population law of motion. The preceding population slide may be unnecessary. Very important: clarify for Tommaso and the audience the distinction between households, agents, children, and the other population objects.

**Status:** open; concise first-pass review requested. No changes authorized.

### M15 — Speed of convergence

Why does convergence to a new equilibrium take so long? Why not one full generation, or two?

**Status:** open; concise first-pass review requested. No changes authorized.

### M16 — Equilibrium exposition as a whole

The sequence of slides around the equilibrium concept is confusing and needs substantial clarification.

**Status:** open; concise first-pass review requested. No changes authorized.

### M17 — Interest rate and calibration

Is i too low? More generally, the whole calibration needs a much clearer explanation and must be clearer in Tommaso’s own understanding.

**Status:** open; concise first-pass review requested. No changes authorized.

### M18 — Fertility measurement

Are the slides showing TFR or completed fertility?

**Status:** open; concise first-pass review requested. No changes authorized.


Each issue will have a stable number and contain:

- **Feedback / question:** preserve the concern as Tommaso states it.
- **Slide or topic:** identify the affected material.
- **Next step and owner:** Tommaso's hand edit, slides edit, or a specific question for another task.
- **Status:** open, investigating, awaiting an answer, or resolved.
- **Answer / decision:** record the resolution and supporting source beside the original issue.

**Current authorization:** record the issues and request one separate task, at max reasoning, to give concise first-pass answers to all 18 points. That task may flag problems and uncertainties only; no code inspection, implementation, slide edits, or model changes. Further investigation and fixes require individual authorization; M01 has been authorized and completed. Keep each eventual answer beside its original issue, without treating an unanswered concern as resolved.

## Resolved issues

### M01 — Earnings terminology

Replace “heterogeneous income” with “idiosyncratic earnings risk.”

**Status:** resolved; Tommaso authorized this edit with “do it.”

**Resolution:** On “This Paper,” replaced “Income heterogeneity, earnings risk, and bequests” with “Idiosyncratic earnings risk and bequests.” The model section still identifies permanent income groups separately.

## Related material

- [Working presentation](../../output/pdf/september_14_presentation.pdf)
- [Slides collaboration handoff](../prompts/HANDOFF_september14_slides.md)
- [Existing empirical follow-up ledger](POST_PRESENTATION_ISSUES.md) — earlier empirical issues, not yet classified as mock feedback.
