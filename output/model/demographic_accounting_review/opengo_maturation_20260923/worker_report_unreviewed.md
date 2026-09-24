 

 

 

 

 

 

 

 

# Population accounting review for the September 24 author decisions

**Source reviewed:** live working tree `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26` (main checkout, files read directly on 2026-09-23; commit hash not verifiable because shell use was prohibited). Prior packet reused, not redone: `/Users/tommasodesanto/.codex/worktrees/d606/Fertility_Spring26/output/model/demographic_accounting_review/` (`final.md`, `population_conversion_literature.md`, `lead_review.md`). No serialized checkpoint values were re-verified; per the standing gotcha, source defaults are not run configuration, so run values below are cited as source defaults/profile values.

## Findings (six)

1. **Two disconnected demographic clocks coexist.** In the stationary solver, new households are generated from *measured household departures*: each child leaving home credits `entrant_conversion_factor = 0.5` households to the entrant pool (`solver.py:4952-4957`, `5653-5658`; entry equation `solver.py:606-610`; factor set in `e5_profile.py:182`). In the dated transition, entrants come from a *birth-vintage queue*: top-code-adjusted births × `1/2.1`, delayed four slots (20 years), times a retention factor (`run_e5f_open_population_transition.py:351-370`, `1228-1244`). The transition computes measured maturation (`run_e5f_open_population_transition.py:898-906`, returned at `1214`) and does not use it for entry.
2. **The two clocks also use different child units and different conversion factors.** Stationary entry counts literal departures (at-home count capped at 3) × 1/2. The queue counts top-code-adjusted births (explicit births + 0.602 per family entering the 3+ bin, `run_e5f_open_population_transition.py:768`; weight 3.602359422009 in `e5_profile.py:181`) × 1/2.1. The extra 0.602 children per 3+ family exist only in aggregate renewal accounting, never in any household state.
3. **Children of dying households vanish.** Survival multiplies the whole household cohort before aging (`solver.py:4871-4872`, `5503-5504`; transition `run_e5f_open_population_transition.py:868-871`, terminal node `:908`). Dependents in the dying share are not matured, not entered, not orphaned; a dedicated read-only tool measures this loss from saved states (`measure_e5f_dependent_parent_death.py:13-36`). With no maturation cutoff, my hand arithmetic gives 4.0 to 14.2 per 100 births lost this way, depending on birth timing (below).
4. **Maturation is a household departure draw, not biological aging.** Each child at home leaves independently with probability μ = 1/4.5 per four-year period (`parameters.py:75`, `1212-1243`), giving a geometric duration with mean 18 years and no maximum; a child can still be "at home" when the parent is 82. Whether a departure equals the child reaching adulthood is an interpretation, not something the code checks.
5. **Normalization does mask nonconservation in the normalized closures.** Under `normalize_population_mass`, total mass is rescaled to `N_target` every forward pass (`solver.py:4989-4996`, `213-216`), absorbing any entry/death imbalance into scale. The live `outside_option_benchmark_normalized` closure instead solves scale from the entry identity (README "Population Closure"); the transition imposes household totals and age shares externally through the Census/ACS bridge (`run_e5f_open_population_transition.py:1251-1265`). Conservation is imposed, not generated, in all three places.
6. **No source-proven bug found; the established objects are inconsistencies between two documented laws, plus the dependent-death gap.** The default-off parent-age maturation machinery (`parameters.py:82-85`, `1094-1113`, `1193-1205`) already implements a ramp to certain exit by age 62 with a newborn exemption — the cutoff option tomorrow does not require new code, only activation and re-normalization.

## What each object counts (Q1)

- **n (children ever born):** integer state 0/1/2/3+, monotone increasing; a realized birth moves the household (n, m) → (n+1, m+1) (`solver.py:5239-5245` first birth, `5252-5281` continuation; destination `solver.py:301-310`). Fertile nodes are ages 18–42 (`parameters.py:45-46`, `solver.py:5217`).
- **m (children at home):** integer 0..min(n,3) (`parameters.py:1044-1056`). Drives child costs, the child-room floor, and the children-at-home moments. Top-coded: a 3+ household's m never exceeds 3.
- **Maturation:** per-child binomial thinning, μ = 1/4.5 (`parameters.py:1229`); a fall m→m′ credits (m−m′) × 0.5 to `entrants_mature_total` (`solver.py:4952-4957`).
- **Adult entrants:** enter at age 18 with the exogenous entry-wealth distribution. Stationary: `(outside_flow + 1.0 × mature_flow) × city_probs` (`solver.py:606-610`). Transition: `outside_flow × outside_shares + retention × queue_B` (`run_e5f_open_population_transition.py:1240-1247`), queue conversion `1/2.1` (`:302-311`), retention anchored at the old state (`:336-348`).
- **Household death:** post-retirement only; survival 1 until 66, then 0.939/0.918/0.885/0.830 at 66/70/74/78, certain exit at 82 (`run_e1_chain.py:22-27`, `364-370`, activated `:383,392`). Deaths drive the warm-glow bequest flow (`solver.py:1148-1175`, `6851-6893`).

## Worked examples (Q2)

**Example A — stationary law, 100 child units born in one period.** Each period each child departs with μ = 2/9; (7/9)^k remain after k periods; expected counted duration 4.5 periods = 18 years. If no parent ever died: all 100 eventually depart, crediting 50 entrant households (×1/2), and with completed fertility 2.1 the next generation's births from 50 households… note the asymmetry: 50 households at 2.1 adjusted births give 105 adjusted child units, but only ~100 × (literal/adjusted ratio) literal departures next generation — the 3+ adjustment (0.602 per top-bin family) never re-enters as a dependent, so literal-unit replacement is not exact even before deaths.

**Example B — same cohort with the actual survival schedule, worst case (all born at parent 42).** Fraction at home at parent ages 66/70/74/78/82: 0.221/0.172/0.134/0.104/0.081. Losses to parental death: 0.221×0.061 + 0.172×0.082 + 0.134×0.115 + 0.104×0.170 + 0.081×1.0 = 1.3+1.4+1.5+1.8+8.1 ≈ **14.2 of 100 omitted**; 85.8 mature → 42.9 entrant households. Best case (all born at parent 22): same calculation gives ≈ **4.0 of 100 omitted**, 48.0 entrant households. The queue would deliver 100/2.1 = 47.6 households regardless of timing. So the death gap is material for late births and is dominated by the terminal-node tail — the geometric tail, not working-age mortality, does the damage.

**Example C — transition law, same 100 (adjusted) child units.** They ride in household m as dependents (costs, rooms), and independently the queue schedules 47.6 households due in 20 years. Each child is counted once as a dependent; entry is tied to the birth event, not the departure event, so no child is double-counted in a given period, but the departure margin and the entry margin can drift apart arbitrarily over a path (that is the renewal-residual object in `audit_closed_reproductive_closure.py`).

## Roles of the factors (Q3)

- **1/2.1** is an external lineage normalization for the transition queue only: at completed fertility 2.1 it makes a generation replace itself in household units (`run_e5f_open_population_transition.py:302-311`). It is not derived from the two-adult abstraction; the gap between 1/2 and 1/2.1 is an unmodeled 4.8% loss between birth and household formation (prior packet, `population_conversion_literature.md`).
- **1/2** (`entrant_conversion_factor`, `e5_profile.py:182`) is the literal two-adult conversion applied to measured departures in stationary accounting. The author has kept these distinct; no change adopted.
- **Survival** is whole-household death; there is no child survival state and no spousal survivor state.
- **Normalization** (`N_target`, `solver.py:4989-4996`) is a scale device for normalized stationary closures; stationary *distributions* are scale-free shares, while transition *counts* are levels driven by the queue plus the external bridge. Yes: under normalized closures the rescale hides any entry/death imbalance; under the live benchmark closure scale is solved from the entry identity instead.

## Decision table (Q4)

| # | Current rule (source) | Economic reading | Established issue | Minimal alternative | Steady-state / transition consequence | Recalculate after choice |
|---|---|---|---|---|---|---|
| 1 | Constant per-child exit hazard μ=1/4.5, no cutoff (`parameters.py:75,1212-1243`) | 18-year average at home, geometric tail | Dependents survive to parent 82; death gap 4–14 per 100 births (Ex. B) | Forced exit on the 62→66 draw; re-set μ externally to keep birth-weighted 18 years (~0.185 illustrative, prior packet) | Removes death gap identically; shortens late-birth durations; every child-cost and room-floor path shifts | Pi_child, all policies, stationary distribution, children-at-home moments, 2.1 normalization (psi), full refit |
| 2 | Dependents of dying households vanish (`solver.py:4871-4872`) | Orphans/guardianship unmodeled | Violates dependent identity C′=C+B−M−L with L>0; bequests never reach dependents | Option 1 makes L≡0; else explicit orphan-entry or measured-loss disclosure | If L=0 via cutoff: estates unchanged in structure; entry flow rises slightly | Dependent-loss receipt (`measure_e5f_dependent_parent_death.py`), estate-flow target, entry accounting |
| 3 | Stationary entry = 1/2 × measured literal departures (`solver.py:606-610,4952-4957`) | Two adults per formed household | Units (literal ≤3) differ from queue units (adjusted) | Keep; or unify units before unifying factors | Scale and age-profile anchors move with any unit change | Old-state renewal anchor (retention, outside flow), population scale, housing levels |
| 4 | Transition entry = 1/2.1 × adjusted births, 20-yr queue (`run_e5f_open_population_transition.py:351-370,1228-1244`) | External replacement normalization | Disconnected from household maturation; measured mature flow discarded | Common clock: entry = conversion × measured matured flow (prior packet recommendation) | Transition timing of entry changes; bridge still pins totals | Full transition re-solve; historical validation; policy paths |
| 5 | 3+ children only in aggregate renewal (+0.602/top-bin entry, `:768`) | Top-bin mean 3.602 vs state cap 3 | Represented children never incur costs or mature | Aggregate pool with its own maturation, or 1.201 weight on top-bin dependents | Renewal flow and scale only; household block untouched if pool stays aggregate | Adjusted-birth series, queue, scale; no policy re-solve if household state unchanged |
| 6 | Entry factors 1/2 vs 1/2.1 kept distinct (author, Sept 23) | Formation loss vs couple pairing | Replacement exactness holds only for the queue law | Defer (author's low priority) | — | None now; document which law each result uses |

No alternative is recommended here on numerical convenience; rows 1–2 are tomorrow's agenda, rows 3–6 are recording/deferral decisions.

## Missing evidence

- Serialized run values in the frozen September 14 checkpoint and the selected E5F candidate (survival schedule, ecf, top-bin weight as actually run) — not re-verified; no execution permitted.
- The measured dependent-loss L for the current selected state (tool exists; was not run).
- The production birth-age distribution needed to pin the re-normalized μ (~0.185 is illustrative from the F3 note, not the run's distribution).
- Whether `run_e5f_transition_calibration.py`'s old-state parameters carried `use_age_survival=True` in every certified run (it inherits `old_parameters.survival_probs`, `:2189-2190`).

## Files inspected (main checkout unless noted)

`memory/AGENT_MEMORY.md`; `memory/daily/2026-09-23.md`; `CALIBRATION_STATUS.md`; `code/model/README.md`; `code/model/intergen_eqscale_seq_optimized/parameters.py` (45-46, 60-139, 1044-1270); `.../e5_profile.py` (150-191); `.../e5_maturation_repair.py`; `.../run_e1_chain.py` (15-64, 360-404); `.../solver.py` (200-330, 590-649, 1134-1175, 4621-4746, 4860-5009, 5133-5289, 5500-5689, 6668-6893, 7206-7213); `code/model/tools/run_e5f_open_population_transition.py` (295-404, 740-939, 1210-1269); `code/model/tools/measure_e5f_dependent_parent_death.py`; d606 packet `final.md`, `population_conversion_literature.md`, `README.md`.

## Verification checklist

**Source-established:** the two clocks and their line citations; μ=1/4.5 and the binomial form; ecf=0.5 and 1/2.1 queue conversion; survival schedule and its whole-household application; dependents vanishing at death; N_target rescaling; the 0.602 top-bin adjustment; the default-off parent-age machinery. **Hand arithmetic (verified twice):** geometric shares (7/9)^k; Example B loss sums 14.2/4.0 per 100; queue 100/2.1=47.62. **Hypotheses/illustrations, not established:** the ~0.185 re-normalized hazard; the birth-weighted average loss (~5/100) — bounded by 4.0–14.2 but not measured; that the stationary and transition laws agree at the reference state beyond the constructed old-state anchor. **Not done (prohibited):** any execution, checkpoint inspection, or git identity check; source identity is the working tree as read.
