# Matched policy transition: preliminary slide

**Author's final display choice:** use `policy_fertility.pdf` / `.png`, showing
only fertility without headline effect numbers. The main deck's Policy Results
frame and standalone slide now use that figure. The earlier two-panel graph is
retained as a diagnostic. The paths are from the one-permanent-shock history,
not the announced four-shock history. The first fertility observation is after
the unexpected 2023 tax increase, so the two impact decisions can differ despite
identical inherited states.

The author subsequently requests axes and legend only in the slide, with no
status or experiment annotations. Numerical limitations remain in this note
and the verification receipt; removing slide annotations does not certify the
equilibrium.

## Comparison with Coven et al. (2025 version)

Their June 19, 2025 paper (local `docs/reference/coven2025_property_tax.txt`,
sections 5.1–5.5) compares stationary equilibria after California's property
tax rises from 0.8% to 2%, with equal rebates. Their
[July 2025 NBER presentation, slide 28](https://conference.nber.org/conf_papers/f222903/f222903.slides.pdf#page=28)
reports California prices -11.2%, aggregate ownership +6 percentage points
(61% to 67%), and ownership at ages 25–44 +8 points (35% to 43%).
Our verified stationary 1%-to-2% comparison instead gives prices +1.1553%
and aggregate ownership -0.6265 points. Our preliminary 2063 transition gives
prices -6.2005% and ownership +0.8188 points, but that date is not a stationary
comparison. We have not established a matched age-25–44 ownership response.
The models also differ in geography, migration, demographic adjustment and
baseline conditions. These differences preclude treating our experiment as a
replication or mechanically scaling their numbers by the tax change. The
present fertility result is not evidence that their housing response has been
quantitatively reproduced. Falling aggregate housing alone would not refute
their mechanism: their paper also has a smaller stock and a tenure reallocation.

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
