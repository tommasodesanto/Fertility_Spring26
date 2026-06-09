# One-Market Policy Proof Of Concept

Date: 2026-06-09

## Scope

This note documents fixed-theta policy comparisons from the final global-DE
toy best of the one-market intergenerational housing-fertility model:

- run directory:
  `results_intergen_housing_fertility_intergen_candidate_no_timing_v0_globalde_3g_20260609`
- task `27`, case `236`, label `de_g008_i011`
- baseline rank loss under `candidate_no_timing_v0`: `11.503191936648555`

The exercise does not recalibrate. It asks whether the current toy model can
produce inspectable proof-of-concept policy responses under one common code
path and one common diagnostic packet.

The reproducible driver is:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
code/model/.venv/bin/python code/model/tools/run_intergen_policy_poc.py
```

Primary output:

```text
output/model/intergen_policy_poc_20260609/
```

Comparison panels:

```text
output/model/intergen_policy_poc_20260609/figures/
```

## Cases

| Case | Policy object |
|---|---|
| `baseline` | Final global-DE toy calibration with zero estate tax |
| `new_parent_ltv95` | Raises financed share to 95 percent only for owner choices by new parents in the birth child-state |
| `property_tax_up_1pp` | Raises annual property tax from 1 percent to 2 percent and re-clears the house price under \(q=(r+\delta+\tau^p)P\) |
| `estate_tax_30pct` | Applies a 30 percent tax to positive terminal resources before bequest utility |

The bequest-tax case is intentionally narrow. It changes terminal bequest
utility only. The model still has no estate-tax revenue rebate, no inheritance
transfer kernel, and no bequest-principal adding-up. Since the baseline has
zero estate tax, bequest-tax removal is read as `baseline` relative to
`estate_tax_30pct`.

## Results

| Case | Loss | Price | User cost/rent | TFR | Childless | Prime own | Old own | \(H_{0\to1}\) | \(H_{1\to2}\) | Young liquid wealth/income |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `baseline` | `11.503` | `0.635` | `0.182` | `1.660` | `0.238` | `0.222` | `0.783` | `0.782` | `0.501` | `0.464` |
| `new_parent_ltv95` | `11.459` | `0.635` | `0.182` | `1.677` | `0.236` | `0.227` | `0.793` | `0.791` | `0.453` | `0.466` |
| `property_tax_up_1pp` | `11.435` | `0.557` | `0.182` | `1.663` | `0.239` | `0.233` | `0.783` | `0.799` | `0.476` | `0.468` |
| `estate_tax_30pct` | `11.429` | `0.634` | `0.182` | `1.661` | `0.238` | `0.224` | `0.789` | `0.779` | `0.495` | `0.465` |

All four cases clear the aggregate housing-services residual below
`4.1e-05`.

## Interpretation

The proof-of-concept policy code is usable, but the current toy calibration
does not yet provide credible quantitative policy estimates.

First, the parent-targeted credit relief moves outcomes in the expected
direction, but weakly. TFR rises from `1.660` to `1.677`, childlessness falls
from `0.238` to `0.236`, and prime-age ownership rises from `0.222` to
`0.227`. The response is small because the baseline lifecycle ownership margin
is already distorted.

Second, the property-tax experiment shows the capitalization channel directly:
the asset price falls from `0.635` to `0.557`, while the flow user cost/rent is
essentially unchanged at `0.182`. This is exactly the channel intended by the
stationary user-cost closure, but it does not repair the lifecycle ownership
profile.

Third, the terminal bequest-tax wedge has small aggregate effects at this
theta. Relative to the taxed case, tax removal lowers old ownership from
`0.789` to `0.783`, lowers aggregate ownership slightly, and barely moves
fertility. Because the wedge is only inside terminal utility, this should not
be interpreted as a complete fiscal counterfactual.

## Plot Packet

The policy runner writes one standard comparison packet:

| Panel | File |
|---|---|
| Prices and user costs/rents | `figures/01_prices_user_costs.png` |
| Headline target moments | `figures/02_headline_moments.png` |
| Lifecycle ownership and fertility | `figures/03_lifecycle_profiles.png` |
| Wealth-distribution diagnostics | `figures/04_distributions.png` |
| Owner-entry policy at age 30, \(z=1\) | `figures/05_owner_entry_policy_age30_z1.png` |
| Owner-entry policy at age 42, \(z=1\) | `figures/06_owner_entry_policy_age42_z1.png` |
| Owner rung and tenure service use | `figures/07_tenure_owner_rungs.png` |

Each case also writes the full existing diagnostic packet under
`<case>/diagnostics/`.

## Caveats

The same pathologies documented in
`docs/model/intergen_housing_fertility_pathology_audit_20260609.md` remain
visible:

1. Prime-age ownership is far below target while old ownership is near target.
2. Fertility remains too late in the lifecycle, although first-birth age is not
   a target in `candidate_no_timing_v0`.
3. Owner-entry policies retain large nonmonotone spikes and drops over liquid
   wealth.
4. Owner services remain concentrated at the temporary `5.2` owner rung.

The policy-comparison machinery is now tidy enough for proof-of-concept
experiments. The economics still require a better calibrated baseline, a
deliberate owner ladder, and a cleaner ownership-entry policy before policy
results can be presented as quantitative evidence.
