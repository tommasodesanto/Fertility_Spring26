# Utility and pension visual review

Open `index.html` for the transition overview, direct old/new utility lifecycle
comparison, policy functions and all eight original 17-graph sets. The concise
nine-page PDF is `../../../pdf/e5f_utility_pension_visual_review.pdf`.

Regenerate with one command from the repository root:

```sh
python3 code/model/tools/build_e5f_utility_review_packet.py
```

Requires matplotlib and reportlab (available in the local system Python).
No solver is imported or run. The measured build stage takes about four seconds;
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

The completed late-preference comparison now appears on PDF page9 and in the gallery, with complete target/parameter links. Both arms reproduce exactly;74 lightweight artifacts were hash-verified. These cases were not recalibrated.
