# One-Market Intergenerational Housing Fertility Audit

Date: 2026-06-08

## Scope

This audit covers the diagnostic one-market implementation in
`code/model/intergen_housing_fertility/`, especially the 2026-06-07 Torch run
`intergen_old_nonlocation_age_tfr_overnight_20260607`.

The run was a diagnostic random screen, not a production calibration. It used
`64` tasks, `48` cases per task, `J=16`, `Nb=70`, `n_house=6`, and target set
`old_nonlocation`.

## Run Status

- All `64` tasks finished.
- All `3,072` cases produced valid JSON records.
- `2,868` cases cleared the one-market housing residual to
  \(\le 5\times 10^{-3}\).
- `2,650` cases cleared to \(\le 10^{-4}\).

The scalar-best converged case is not economically usable:

| Moment | Model | Target |
|---|---:|---:|
| `tfr` | `3.548` | `1.700` |
| `childless_rate` | `0.003` | `0.150` |
| `mean_age_first_birth` | `30.483` | `26.000` |
| aggregate ownership | `0.937` | diagnostic |
| prime-age ownership | `1.000` | `0.575` |
| old-age ownership | `1.000` | `0.764` |
| family ownership gap | `0.000` | `0.168` |
| old parent-childless ownership gap | `0.000` | `0.070` |

## Coding Findings

### Fixed: ownership target mapped to the wrong statistic

The old workhorse objective maps the `own_rate` target to
`sol.own_rate_3055`, i.e. prime-age ownership. The one-market diagnostic
extractor was mistakenly mapping `own_rate` to aggregate ownership
`sol.own_rate`.

This has now been patched:

```python
"own_rate": float(getattr(sol, "own_rate_3055", np.nan)),
"aggregate_own_rate": float(getattr(sol, "own_rate", np.nan)),
```

The saved Torch run was re-ranked under the corrected mapping. The scalar-best
candidate remains the same high-fertility, high-ownership corner, so this bug
is real but not the sole cause of the bad result.

### Verified: TFR convention

The model stores household/unitary-agent mean parity as
`mean_completed_fertility`. The paper-facing TFR target uses

\[
\text{tfr}=2\times \text{mean_completed_fertility}.
\]

The extractor and logs now follow this convention.

### Verified: child maturity timing

The code now uses 4-year periods with

\[
\text{stage\_durations}=[A_m/\text{period\_years}]=[18/4]=[4.5].
\]

Thus the single dependent-child state has expected duration 18 years. This is
still a simplification: child age is not deterministic or age-resolved.

## Main Economic/Numerical Failure

The first-birth-age target dominates the raw squared-loss objective. In the
current one-shot fertility architecture, the model has no independent
age-specific fertility shifter. Fertility timing is generated only by the
lifecycle budget, housing, and continuation-value tradeoff. The search
therefore lowers the age-at-first-birth loss by pushing fertility early, which
also drives childlessness toward zero and TFR too high.

For the scalar-best case, the largest loss components are:

| Component | Loss Contribution |
|---|---:|
| `mean_age_first_birth` | `241.149` |
| `tfr` | `40.979` |
| all remaining targeted moments combined | about `4.323` |

This explains why the optimizer accepts an absurd corner: once age at first
birth is hard to move, it trades off the rest of the objective poorly.

## Re-Rank Diagnostics

Under the corrected ownership mapping, the best converged scalar case still has
loss `287.041` and remains a corner:

| Moment | Model |
|---|---:|
| `tfr` | `3.548` |
| `childless_rate` | `0.003` |
| `mean_age_first_birth` | `30.483` |
| `own_rate_3055` | `1.000` |
| `old_age_own_rate` | `1.000` |

If only the first-birth-age moment is removed from the audit loss, the best
candidate is much more informative:

| Moment | Model | Target |
|---|---:|---:|
| `tfr` | `1.739` | `1.700` |
| `childless_rate` | `0.294` | `0.150` |
| `mean_age_first_birth` | `33.119` | `26.000` |
| `own_rate_3055` | `0.617` | `0.575` |
| aggregate ownership | `0.655` | diagnostic |
| `old_age_own_rate` | `0.836` | `0.764` |
| `housing_increment_0to1` | `0.702` | `0.664` |
| `housing_increment_1to2` | `-0.059` | `0.566` |
| `young_liquid_wealth_to_income` | `0.511` | `0.600` |
| `old_age_parent_childless_gap` | `0.036` | `0.070` |

This says the implementation can get some non-age moments into the right
region, but the current architecture/search cannot jointly match early
first births, childlessness, and completed fertility.

## Search Design Problem

`calibrate-small` is a random screen. It is not an optimizer:

- no incumbent seeding;
- no local refinement;
- no staged objective;
- no parameter transformation tuned to economically admissible regions;
- no automatic rejection of high-fertility/high-ownership corners except through
  the raw moment loss.

Running longer random screens is therefore low value until the objective and
timing problem are fixed.

## Current Model Limitations Relevant To The Failure

- Fertility is one-shot completed-family choice for childless fertile
  households, not sequential parity hazards.
- There is no age-specific fertility utility or age-specific biological/time
  cost shifter.
- `mean_age_first_birth` is computed over 4-year period starts
  \(26,30,34,38,42\). Using period midpoints would make the measured model age
  even later, not earlier.
- The single dependent-child state is stochastic with mean 18 years, not
  age-resolved.
- The diagnostic search varies `c_bar_n` but leaves `c_bar_0` fixed, so it has
  limited control over the extensive fertility margin separately from intensive
  parity.

## Required Next Fixes Before Another Overnight Run

1. Keep the corrected `own_rate -> own_rate_3055` mapping.
2. Produce diagnostics for the corrected no-age-audit best and the scalar-best
   corner before any new search.
3. Add an explicit fertility-timing margin if `mean_age_first_birth` remains a
   hard target. Minimal options:
   - an age-specific fertility utility shifter;
   - an age-specific child cost schedule;
   - or a staged objective that first fits TFR/childlessness/ownership and then
     tunes timing.
4. Replace blind random search with a staged local search around economically
   admissible candidates.
5. Rename reported ownership fields in CLI/diagnostics so aggregate ownership
   and prime-age target ownership cannot be confused again.

