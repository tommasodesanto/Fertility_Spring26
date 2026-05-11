# Twins-IV Extension Plan (PSID)

## Objective

Estimate a surprise-child effect (LATE) using twins, and compare it to the general childbirth effect from the existing event-study.

## Baseline to Compare Against

Current baseline in this project:

- Sun-Abraham event-study around first childbirth (`f_c_y`) for outcomes like `own`, `rooms`, `mv_s`, `mv_n`.
- This identifies the average dynamic association around childbirth (not necessarily causal for fertility surprises).

## Causal Target

Use twins as an instrument for an unexpected increase in child quantity.

Two comparisons are needed:

1. **General effect** of childbirth timing (existing SA event-study).
2. **Twins-IV effect** (surprise variation in children) on housing/tenure outcomes.

## Data Build Required

Construct an event-level panel with:

- Fertility timing: first birth year, parity path, spacing.
- Twin indicator at first birth (or at target birth margin).
- Outcome panel: `own`, `change_own`, `moved_to_own`, `rooms`, `move` reasons.
- Baseline controls: age, education, sex, pre-birth income/wealth.

## Estimation Roadmap

### Step 1: Replicate baseline SA moments

- Keep existing SA specification for key outcomes.
- Freeze this as the "general effect" benchmark.

### Step 2: First-stage diagnostics (twins relevance)

At the person/event level:

- `additional_child_i = pi0 + pi1 * twins_i + X_i'pi + u_i`

Key diagnostics:

- First-stage coefficient sign and size.
- F-statistic / weak-IV diagnostics.
- Balance checks around twins indicator.

### Step 3: Main IV outcome models

For post-birth windows (for example +0 to +3 years):

- `y_i = beta0 + beta1 * additional_child_i + X_i'beta + e_i`
- Instrument `additional_child_i` with `twins_i`.

Run for:

- `own_post3`
- `moved_to_own_post3`
- `rooms_change_post3`

### Step 4: Dynamic IV version (optional extension)

If feasible, estimate dynamic responses by instrumenting event-time treatment interactions:

- `y_it = sum_k beta_k * D_{ik} + FE + e_it`
- Instrument each `D_{ik}` with `twins_i * 1[event_time=k]`.

This is heavier and should come after the static-window IV checks.

## Comparison Table to Produce

One compact table per outcome:

- Column A: SA general-effect estimate (event-study summary statistic or post-window average).
- Column B: OLS post-window effect of additional child.
- Column C: IV post-window effect (twins instrument).
- Row notes: sample restrictions, years, controls, clustering.

## Practical Guardrails

- Keep all new work in separate scripts/files under `code/data/psid_followup_mar2026`.
- Do not overwrite legacy outputs under `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs`.
- Start with years where fertility and wealth coverage are reliable for the chosen design.

## Immediate Next Script (when resuming)

Create a new standalone file:

- `twins_iv_v1.do`

with sections:

1. Sample construction.
2. First stage.
3. OLS vs IV window regressions.
4. Output tables and diagnostics.
