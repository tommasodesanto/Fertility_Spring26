# Family housing requirement diagnostic for mechanism slides

Saved-input, zero-solve comparison of the post-interest chain-13 fixed-price
lifecycle cases in `../../price_decomposition/cases/baseline_resolve` and
`../../price_decomposition/cases/arm_floor_only`. The executed public inputs
differ only in `hbar_first_child_jump` (the estimated (h_P)): 2.593759507364224
versus 2.7814454323543467 rooms, a 7.23605733134638% increase. The
`hbar_child_rooms` multiplier in the arm design operates on zero and changes no
executed value. House price 0.7760569760205563, rent, the 120-node wealth
grid, earnings, entry wealth, financing, preferences, birth timing and all other
economic inputs are identical. Both cases use the no-Estate-A lifecycle model.
The comparison is partial equilibrium; it does not clear a new housing market
or refit the calibration.

`main_floor_state.pdf` is the slide figure (and `.png` its raster copy). It
compares age-22 initially childless renters in the fifth of nine earnings states
at common occupied starting wealth nodes. The 12 shown nonnegative nodes cover
99.66758178490113% of this state's baseline pre-birth exposure. Its three
panels show first-birth probability, expected physical rooms conditional on a
successful first birth, and expected ending net financial wealth $b'$ on that
branch. The latter two average over each arm's tenure choice using the engine's
own transaction map and selected saving controls, including renter stayer
controls. The figure divides beginning and ending wealth by the same baseline
household's annual gross earnings, 0.5602745665736877, calculated as four-year
aftertax income 2.0611813071495764 divided by $4(1-\tau_{pay})$. The CSV keeps
raw wealth and plotted ratios. `main_floor_state.csv` contains all 28
positive-exposure nodes; `main_floor_state_receipt.json` records support,
source hashes and checks. The
baseline branch extraction agrees with the independent household figure table
at all 28 nodes to within 3e-14. Maximum positive-option budget residual is
5.33e-15; minimum consumption and housing surplus are positive. Rooms can
exceed the six-room **rental** cap where owner options receive weight.

At beginning wealth $b=0$, first-birth probability falls from 17.207% to
11.081%, successful-birth rooms rise from 5.383 to 5.395, and expected $b'$
rises from 0.236 to 0.302 (0.422 to 0.540 relative to annual gross earnings).
This is a matched state response, not an aggregate
fertility effect. The higher floor raises the minimum housing service requirement
from 2.594 to 2.781 rooms; it is itself well below the six-room rental cap.
Chosen rooms for this median-earner group approach that cap. The source
`crowding_cost_fable` cases combine changes in $h_P$, per-child rooms and,
in other worlds, consumption share or the rental cap, so they cannot isolate
this effect.

`first_birth_floor_common_states.pdf` and `.png` retain broader appendix
evidence at ages 22 and 30 by fixed earnings and wealth bins. The 42-row
`first_birth_by_age_earnings_wealth.csv` uses baseline pre-birth exposure weights
and records each group's support. The exact total first-birth flow changes from
0.05007322691325753 to 0.048718040658776045 (−2.70640886961232%). Its
baseline-exposure policy term is −0.004369135214085933; changed exposure adds
+0.0030139489596044604. The identity residual is 1.34e-17. `receipt.json`
contains the source hashes and exact inversion checks.

Regenerate from saved arrays, without a model solve:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg NUMBA_NUM_THREADS=1 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
code/model/.venv/bin/python output/model/credit_mechanism_20261004/mechanism_slides/floor/extract_floor.py
PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg NUMBA_NUM_THREADS=1 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
code/model/.venv/bin/python output/model/credit_mechanism_20261004/mechanism_slides/floor/build_main_figure.py
```
