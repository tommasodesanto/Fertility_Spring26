# New Chat Handoff: Bequests Must Stay; Identify Them Correctly

Date: 2026-07-15

You are taking over the quantitative housing-fertility project after a sequence
of old-age ownership, mortality, and bequest calibration exercises. The author
wants a clean conceptual reset and a concise advisor-level discussion before
any more long runs.

## Author's Non-Negotiable

The model should retain a bequest motive. Do **not** infer from a near-zero
estimate under the current target system that bequests should be deleted.
Bequests matter economically for late-life saving, and the author wants them in
the model. The question is:

> What is the most parsimonious bequest specification, and which wealth or
> estate moments identify each of its parameters separately from patience,
> housing demand, and tenure choice?

The author is rightly frustrated by a loop in which each failed fit produces a
new mechanism and then a recommendation to discard it. Do not continue that
loop. Separate model content, empirical identification, and goodness of fit.

## Mandatory Startup

Before giving substantive advice, read in this order:

1. `memory/AGENT_MEMORY.md`
2. the latest `memory/daily/YYYY-MM-DD.md`
3. `CALIBRATION_STATUS.md`
4. `code/model/README.md`

Then read the specific files below. Treat `CALIBRATION_STATUS.md` and the
listed result packets as canonical when older notes disagree.

## Active Model And Bequest Code

The active strand is the one-market lifecycle model under:

- `code/model/intergen_housing_fertility/`

Read:

- `code/model/intergen_housing_fertility/solver.py`, especially
  `bequest_utility_vec` and the construction of `Vbq`.
- `code/model/intergen_housing_fertility/parameters.py`
- `code/model/intergen_housing_fertility/calibration.py`
- `code/model/intergen_housing_fertility/tests/test_bequest_normalization.py`

Current bequest utility is evaluated on terminal estate resources

\[
W=b+pH,
\]

where liquid wealth `b` may include mortgage debt, so this is net estate wealth
under the model's balance-sheet convention. The normalized child-blind form is

\[
B(W)=\theta_B\left[
\frac{(\bar b+W)^{1-\sigma}-\bar b^{1-\sigma}}{1-\sigma}
\right].
\]

In code this is `bequest_spec="linear_child_scale"`,
`normalize_bequest_utility=True`, and `theta_n=0`. The luxury shift is
`theta1=\bar b`; overall strength is `theta0=\theta_B`.

Other implemented opt-in forms are:

- `parent_gated_luxury`: zero bequest utility for childless households.
- `equal_division_luxury`: estate divided equally across children, with total
  joy of giving summed across children.

Do not add dynastic altruism. It is computationally expensive and outside the
paper's intended scope.

## Settled Numerical Evidence

### Mortality

The clean mortality exercise re-estimated all 11 core parameters against all
15 current moments:

- M0, no mortality/no bequest: loss `6.973326`.
- M1, externally pinned SSA survival/no bequest: loss `6.860325`.

Mortality improves the population and housing-mass accounting but does not
materially change ownership conditional on survival. Mortality is an externally
measured transition, not a parameter to be identified by the SMM objective.

Packet:

- `output/model/intergen_mortality_recalibration_20260715/report/`

### Clean child-blind bequest arm

M2 combined the same mortality schedule with the normalized child-blind warm
glow, no post-retirement owner-LTV taper, and a fixed luxury shift
`theta1=0.25`. It re-estimated the 11 core parameters plus `theta0`:

- M1 loss: `6.8603247`.
- M2 loss: `6.8367871`.
- M2 `theta0=0.00030099`, at the soft-zero boundary.
- Old ownership and old housing mass barely move.
- The strict winner reproduces exactly.

Packet:

- `output/model/intergen_mortality_bequest_recalibration_20260715/report/`

This result says the **current objective** does not reward a positive
child-blind bequest coefficient. It does not establish that bequests are
economically irrelevant or correctly identified.

### Earlier no-taper and equal-division evidence

The July battery also found:

- A0, no bequest/no owner-LTV taper: loss `6.7881`.
- A1, parent-gated bequest/no owner-LTV taper: loss `6.7430`,
  `theta0=0.0211`.
- A4, equal division under the now-rejected borrowing schedule:
  `theta0=0.0106`.

The age-dependent owner borrowing schedule is rejected because it creates a
late-life ownership cliff. Do not use those schedule arms to judge bequests.

Packet:

- `output/model/intergen_bequest_exit_battery_20260714/report/`
- `docs/model/bequest_exit_note_20260715.tex`

## The Central Identification Problem

The active objective already contains two late-life wealth moments:

1. `old_nonhousing_wealth_to_income_median_6575`
   - target `2.23046078`
   - weight `0.8`
2. `old_parent_childless_nonhousing_wealth_to_income_gap_6575`
   - target `1.00744952`
   - weight `2.0`

But `calibration.py` currently labels both moments `needs-fix`:

- The model and data use inconsistent annual-versus-period income
  denominators.
- The parent-childless object mixes a model ratio-of-sums with a data
  mean-of-ratios.
- Both are **nonhousing wealth** moments, while bequest utility is defined over
  the total estate, including housing equity.
- Their weights are tiny relative to old ownership (`160`) and aggregate/young
  ownership (`100`/`80`).

Therefore the existing calibration is not a clean test of whether wealth data
identify a bequest motive. Audit and repair these empirical objects before
interpreting `theta0≈0`.

The model already computes diagnostic total-wealth objects, including:

- `old_total_wealth_to_income_6575`
- `old_total_wealth_to_income_median_6575`
- `old_parent_childless_total_wealth_to_income_gap_6575`
- median parent-childless total-wealth gaps

Verify their exact definitions and whether matching data objects can be built.
Do not promote a diagnostic simply because it exists in code.

## Identification Audit Already Completed

A strict 15-by-11 Jacobian at the mortality-only M1 winner shows that old-age
ownership is locally redundant for the current 11-parameter block:

| Scaling | Full rank at `1e-2` / `1e-3` | Full condition | Drop old ownership rank | Drop old ownership condition |
|---|---:|---:|---:|---:|
| Target-scaled | 9 / 11 | 470.13 | 9 / 11 | 487.02 |
| SMM-weighted | 10 / 11 | 191.43 | 10 / 11 | 191.65 |

Old ownership loads strongly on `chi`, but aggregate and young ownership load
on `chi` at least as strongly. The remaining weak direction combines
`h_bar_0`, `h_bar_n`, `beta_annual`, and `kappa_fert`; old ownership does not
resolve it.

Packet:

- `output/model/intergen_mortality_identification_20260715/README.md`
- `output/model/intergen_mortality_identification_20260715/main/leave_one_moment_out.csv`
- `output/model/intergen_mortality_identification_20260715/main/parameter_column_summary.csv`

Important limitation: this Jacobian fixes bequests off. The earlier 12-column
bequest Jacobian was evaluated at a rejected, non-optimal borrowing-schedule
point. There is not yet a valid bequest-block Jacobian under a repaired,
wealth-based target system.

## Literature Position To Use

The main specification should rest on published lifecycle/bequest work, not on
the unpublished Scholz--Seshadri paper.

- A De Nardi-style shifted warm glow over total estate wealth is the clean
  published benchmark.
- Fertility plus equal division is a useful robustness specification when the
  paper wants the number of children to dilute per-child inheritances.
- Scholz--Seshadri is only a feasibility precedent for combining endogenous
  fertility with joy-of-giving; it is not the foundation.
- Do not import health states, medical expenditures, or a full dynasty.

Refresh the literature with recent published top-five and top-field work. The
author has explicitly rejected a stale review consisting mainly of 1990s
papers. Distinguish published benchmarks from working papers.

## Your Task

Start with a concise advisor-level recap. Then deliver a disciplined bequest
identification strategy. Do not launch a calibration or edit production model
code before the author agrees on the moments.

### 1. Audit what is already targeted

For each existing old-age wealth moment, report:

- exact data definition and source;
- exact model definition;
- annual/period denominator convention;
- mean, median, ratio-of-sums, or mean-of-ratios convention;
- nonhousing wealth versus total net worth/estate wealth;
- which bequest or saving parameter it can plausibly identify;
- whether the object is currently clean enough to remain a hard target.

### 2. Map parameters to identifying variation

At minimum distinguish:

- `theta0`: overall strength of the bequest motive;
- `theta1`: luxury/threshold shift;
- `theta_n` or equal division: how family size changes the marginal value of an
  estate;
- `beta`: ordinary lifecycle patience;
- `chi` and tenure/housing parameters: desire to retain housing.

Do not claim that a single old ownership level identifies this block.

The likely empirical logic is:

- overall late-life wealth retention or decumulation identifies `theta0`
  relative to `beta`;
- the wealth gradient or upper-tail estate pattern identifies `theta1`;
- parent-childless and completed-fertility wealth/estate gradients identify the
  family-size component;
- housing equity and nonhousing wealth separately help distinguish bequest
  saving from housing retention.

Verify this logic rather than accepting it automatically.

### 3. Propose the smallest credible target package

Give a table with columns:

| Candidate target | Data source | Model analogue | Parameter direction | Measurement status | Keep/fix/drop |

Consider, without automatically requiring all of them:

- total net worth or estate wealth at ages 65--75 and 75+;
- late-life decumulation between age bins;
- wealth or estate percentiles, especially the upper tail;
- housing equity and nonhousing wealth separately;
- parent-childless wealth gaps;
- wealth/estate gradients by completed fertility or number of heirs;
- direct bequest/estate-at-death data if a credible source is available.

The final package must have enough independent variation for every free
bequest parameter. If `theta1` or the family-size shifter lacks a direct moment,
fix it externally or profile prespecified values; do not let the loss choose an
unidentified parameter.

### 4. Compare minimal specifications

Discuss these as choices, not as an invitation to expand the model:

- **B1:** child-blind De Nardi warm glow over total estate;
- **B2:** parent-gated warm glow;
- **B3:** equal division across children.

Recommend one benchmark and one robustness specification. Explain treatment of
childless households, computational cost, literature precedent, and exact
identifying moments.

### 5. End with an executable next step

Only after the target definitions are repaired and approved, propose a compact
calibration battery that:

- keeps mortality and removes the rejected owner-LTV taper;
- re-estimates all 11 core parameters jointly;
- adds only genuinely identified bequest parameters;
- keeps the full target-fit and parameter/bound reporting contract;
- computes a fresh Jacobian at the valid winner;
- runs on Torch with a smoke test, explicit solve count, budget, checkpoints,
  and strict winner repetitions.

Do not use `theta0≈0` as the starting conclusion. The first conclusion to test
is whether the current wealth targets are the wrong objects, inconsistently
measured, or too weakly weighted to identify the bequest block.

## Tone And Deliverable

Write for an advisor who understands the paper but has not followed every
calibration run. Be concise, explicit, and decisive. Lead with:

1. what the current calibration actually established;
2. what it did **not** establish;
3. which wealth moments are already present but defective;
4. the smallest identification-preserving repair;
5. the exact decision the author must make before another run.

Do not respond with another long list of rejected mechanisms. The destination
is a parsimonious model **with** bequests and a defensible empirical strategy
for identifying them.
