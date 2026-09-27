# First-pass lifecycle comparison

`lifecycle_dashboard.png` compares the verified overnight de_0093 population with actual survey microdata/sufficient statistics. It is a supplemental diagnostic, not a new target system or a replacement for the standard 17 plots. No model solve or new stochastic simulation was needed.

## Read the four panels

- Ownership and rooms: national ACS 2005–2006, all structures. This differs deliberately and visibly from the DUE-structure restriction on the ownership calibration target. Physical rooms are capped at 9 in both data and model; this ACS diagnostic is distinct from the uncapped AHS 2007 calibrated rooms mean.
- Net worth: PSID 2005/2007 reference-person family-years, exactly the accepted aggregate wealth/earnings sample. Model net worth is measured before housing transactions using g_post_fertility, matching the calibrated observer. Each age-cell mean net worth is divided by the overall working-age mean annual household gross earnings in its own population. These are comparable normalized profiles, not dollar predictions. PSID mean wealth is noisy, especially in the oldest cells; uncertainty bands remain future work.
- Children: CPS June 2004/2006 women, literal children-ever-born counts capped at 3. Model values average its pre-birth and post-birth distributions within each four-year cell, the existing uniform birth-time approximation. Model households are not a literal female survey sample. Only full CPS four-year cells, ages 18–41, are included; available ages 42–44 do not form the model 42–45 cell. The model terminal capped mean is not the 2.1 fertility normalization, which uses a different final-bin weight.

Each point is a four-year age cell, placed at its midpoint. Lines only connect these observations. These are cross-sectional age profiles, not observed or simulated lifetime trajectories. No lifecycle comparisons have been added as calibration targets.

## Files and reproducibility

- `lifecycle_comparison.csv`: all 125 plotted observations, source labels, age intervals, counts and weight sums.
- `model_profiles.csv`: model population means, mass, pre/post birth counts.
- `psid_four_year_profiles.csv`: empirical wealth levels, scaling denominator, sample counts and accepted aggregate replay.
- `provenance.json`, `source_receipt.json`: exact definitions, references, source identities and authenticated saved checkpoint.
- `report_qa.json`: first-pass checks.

From the repository root:

```
code/model/.venv/bin/python code/model/tools/build_e5f_lifecycle_data_dashboard.py --output output/model/daytime_calibration_20260927/lifecycle_dashboard
```

The builder sets all numerical libraries to one thread. It loads the authenticated checkpoint, reuses small ACS sufficient statistics, reads exactly two hash-verified June CPS partitions from the existing compressed file, and performs one selected-column PSID read with a 240-second cap. It does not run the model. The emitted R extraction code is retained for transparency.

## First-pass findings

The model has lower ownership through most working ages but continues rising after retirement: ages 82–85 are 0.954 versus ACS 0.749. Housing is larger than the ACS comparison through much of working life (ages 42–45: 6.547 versus 5.923 rooms), then smaller in old age. Mean net worth decumulates much more sharply after retirement (ages 82–85: 2.243 versus 6.209 times mean working-age earnings). Births are delayed (ages 30–33 capped children: 1.050 versus 1.446). These are descriptive flags, not identified causal mechanisms, and the different empirical samples matter.
