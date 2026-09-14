# Matched policy transition: preliminary slide

Figure: `policy_transition.pdf` / `.png`. Slide:
`output/pdf/property_tax_transition.pdf`; source:
`latex/appendix_property_tax_transition.tex`.

The author requests original model units in the graphs, with any national
population conversion discussed in words only. Accordingly the two panels
show period total fertility and household mass, with no conversion to people.
Common 2023 household mass is approximately one. The policy raises the annual
property tax from 1% to 2%, with equal rebates in both paths, after the single
permanent preference decline. Both arms inherit the same 2023 state; this is
not the announced four-shock history. The original birth-entry queues, no
immigration, housing supply and fixed-payroll PAYGO contract are retained.

`frozen_policy_paths.json` preserves native rows, fertility observations,
gates, terminal distances, best-iterate receipts, remote paths and input hashes.
The baseline is round 1 mapping 8 (exact replay of best mapping 6); the policy
is round 2 mapping 2. `comparison.csv` contains the full 100-date series; only
2023–2063 is plotted. The builder verifies all three price/fiscal paths against
the saved best coordinates and the plotted values against the native arrays.
Mass, household-budget reproduction and feasibility checks pass. Market and
fiscal convergence is incomplete; horizon adequacy has not passed. The saved
paths and the figures are preliminary, not certified equilibrium results.
`verification.json` distinguishes artifact checks from equilibrium status.

At 2063, the 2% arm has 0.11358% more households and 1.03101% more births;
TFR is 1.78412 versus 1.77070. Birth counts and TFR are distinct because TFR
aggregates age-specific rates, whereas total births depend on age composition.
The separately verified stationary endpoints have 1.80506% more households;
that is not the 2063 effect or a proof of convergence to those endpoints.

For an illustration in words only: the 2023 ACS total U.S. population is
334,914,896 ([Census B01003](https://data.census.gov/table/ACSDT1Y2023.B01003?t=Population+Total)).
The project's national ACS head count for ages 18–85 is 128,086,487, from
`code/data/Spatial_aggregate_withmicrodata/output/national_householder_housing_path/national_householder_housing_path.csv`.
Holding heads per total resident at that ratio means multiplying relative
household mass by the initial population. This yields roughly 329,000 extra
people-equivalent in 2063 and 2.11 million at the computed stationary endpoint.
The latter applies 1.805% to the much smaller baseline stationary scale, not
to today's population. The headship factor cancels in this proportional
conversion. It neither counts dependent children separately nor permits
changing household composition; it is not a resident-population forecast.
The mixed national/42-metro geographic scope is a very urgent open ledger item.

Regenerate without model solves:

```sh
python code/model/tools/build_e5f_policy_transition_slide.py
```

Compile the generated TeX twice from `latex/`, directing outputs to this
folder's `build/`, copy the PDF to the path above, and inspect a Poppler render.
