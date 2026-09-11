# Utility and pension visual review

Open `index.html` for the transition overview, direct old/new utility lifecycle
comparison, policy functions and all eight original 17-graph sets. The concise
sixteen-page PDF is `../../../pdf/e5f_utility_pension_visual_review.pdf`.

Regenerate with one command from the repository root:

```sh
python3 code/model/tools/build_e5f_utility_review_packet.py
```

Requires matplotlib and reportlab (available in the local system Python).
No solver is imported or run. The measured build stage takes about six seconds;
first-use Python imports/font discovery can take longer. Original PNG bytes and
filenames stay unchanged; all136 hashes are checked before packaging. Missing
data or changed hashes stop the build. `plotted_series.json` preserves every
supplemental series, and `review_qa.json` records source-column and render review.

The transition is the converged six-date diagnostic17370186, not a fitted or
horizon-certified historical equilibrium. Its final2023 policies are distinct
from the stationary pre-2007 candidate17370427. The four-cell utility comparison
17358470 holds preference and structural parameters fixed. Its lifecycle filenames
contain2023 for historical plotting-interface reasons; the experiment is pre-2007.
The two old-pension cells are fiscally invalid diagnostics, clearly labeled in
the gallery. The new late-preference stationary comparison is separately scoped
in `../utility_fiscal_decomposition/static_2023/README.md`.

Original plots contain dense income-state legends and visible policy reversals.
These are preserved for inspection, not smoothed or presented as globally verified.

The completed late-preference comparison appears in the appendix and in the gallery, with complete target/parameter links. Both arms reproduce exactly;74 lightweight artifacts were hash-verified. These cases were not recalibrated.


## The fixed tester

The first PDF page and the top of `index.html` answer the same questions every
iteration. `scorecard.json` records the verdicts; `assessment.json` records complete
tables, loss checks, data comparisons and source hashes. Numerical passes apply
only to the named solved experiment. Fit, global policy accuracy, matched2023
recalibration, fitted history and horizon certification are separate checks.
The 12+1 early rows and 12 old/new rows are shown in full; shaded contributions
above4 flag a gap larger than2 objective scales for review, not statistical
rejection. Some scales are synthetic. No targets, weights or solver gates change.

`transition_vs_data.png` uses verified NCHS blocks and Census HH-3 inputs from
`../path_pilot_20260910/fertility_data/`. Decision2007 maps to2008-11 births;
2023 maps to2024-27 and has no observed block in this dataset. Both birth series
are indexed to their first block, comparing shape, not absolute levels or TFR.
Household counts are conditioned on historical information, not independent fit.

`current_policy_extraction/` preserves the read-only cluster recipe and extracted
current2023 policies. Its result receipt pins643 sources and the checkpoint;
17374949 took11seconds and ran no Bellman solves. Probability-weighted housing
is rental probability times conditional rental housing plus each owner-product
probability times its physical housing size. Conditional renter consumption is
not aggregate consumption. The supplemental slice is selected by mass-weighted
ownership declines, not assumed representative. Original graphs remain intact.

The next scientific experiment is a matched2023 recalibration; the present
fixed-parameter late-preference comparison is explicitly not that experiment.
