# Purchase access and rental space at common prebirth states

Verified October 4, 2026. This zero-solve diagnostic reads the retained
post-interest soft chain-13 fixed-price cases, with diagnostic financed shares
\(\phi=.8\) and \(1\), without Estate A. It adopts no model change.

The purchase screen already admits a parent-suitable product for essentially all
locally responsive first-birth states. Current period income explains much of
this screen passage. This necessary condition does not reserve income for a
purchase: consumption and ending assets still share the household budget and
are optimized conditional on each housing choice. Rental space is tight in the first-parent menu, but relaxing credit
changes common-state housing only slightly. These are descriptive mechanism
counts, not an identified causal effect of debt or a test of another timing model.

## Definitions and implemented timing

The exposure is the baseline distribution immediately before birth, reconstructed
with the checked helper in `../diagnostics/extract_common_states.py`, at fertile
ages 18–42 and childless states \((n,m)=(0,0)\). Signed roundoff is retained.
Its total is 0.206400 household mass.

Let \(a\) be first-birth attempt probability, \(\pi_j\) conception probability,
and \(\kappa=0.116522\) the first-birth logit scale. The susceptibility weight
\(W=G^{pre}\pi_j^2a(1-a)/\kappa\) is the derivative of first-birth flow with
respect to successful-birth-minus-wait cardinal utility, holding exposure fixed.
It totals 0.138189; it is not a population probability. This follows the
expected-try logit in `household.py:972`; it differs from a derivative with
respect to the expected try-minus-wait gap.

Price is 0.776057 per room, gross four-year return is 1.082432, and owner
products are 2, 4, 6, 8 and 10 rooms. Budget income is after-tax four-year
income \(y=P.income_{ij}z\), not annual gross earnings. Purchase access is
\(Rb+y+S\ge(1-\phi)Q\), \(Q=Ph\), with net sale
\(S=(1-.06)Ph^{old}\) for owners and zero for renters. `household.py:894`
transforms the threshold by \(R\); `kernels.py:390` applies it at transaction
wealth \(b+(S-Q)/R\). Its transformed borrowing-floor screen is algebraically
equivalent. Same-product owners stay without this new-purchase screen; the
census excludes that product and separately records current-home suitability.

The active first-parent floor is 2.593760 rooms. Strict owner feasibility requires
\(h-\bar h>0\) before the service premium \(\chi=1.050076\) is applied
(`kernels.py:854`). Thus 2 rooms suit childless households but fail the parent
floor; 4+ rooms pass. Screen passage does not imply unconstrained credit. Ending-debt floors, lifetime feasibility and preferences remain active.

## Decisive counts

Percentages divide by baseline exposure or susceptibility of the stated group.

| Soft screen | Household mass | Susceptibility |
|---|---:|---:|
| Parent-suitable purchase already accessible at .8 | 96.253% | 100.000% |
| Only unsuitable small purchase accessible at .8 | 3.525% | \(5.060\times10^{-10}\)% |
| Newly parent-suitable at .95 or 1 | 3.747% | \(5.060\times10^{-10}\)% |
| .8 access depends on including current income | 78.981% | 71.575% |

For beginning renters, .95 newly opens parent-suitable purchases on 4.185% of
mass but only \(5.114\times10^{-23}\) of susceptibility. Beginning owners'
corresponding figures are 0.259% and \(3.803\times10^{-11}\).

The diagnostic screen \(Rb+S\ge(1-\phi)Q\), removing current income while
holding everything else fixed, yields parent-suitable access of 17.272%,
46.266%, 99.948% of mass at .8/.95/1; susceptibility access is
28.425%, 62.299%, 99.980%. This is **not** an exact alternative-theory timing
experiment; it changes only a diagnostic inequality and solves no choices.

Susceptibility shares passing each purchaser screen (same-product stay excluded):

| Rooms | .8 | .95 or 1 | Newly passing |
|---|---:|---:|---:|
| 4 | 93.754% | 93.754% | 0 |
| 6 | 97.270% | 97.270% | \(1.810\times10^{-9}\)% |
| 8 | 99.273% | 99.273% | 0.000513% |
| 10 | 99.782% | 99.792% | 0.009933% |

At age 22, the median childless income node has \(z=.778\), so current
period budget income is 2.061. A four-room down payment is only 0.621: a
renter at \(b=0\) passes the soft screen but fails the before-income diagnostic.
That comparison does not establish that the desired home, consumption and
saving plan are affordable together.
The corresponding asset thresholds are −1.331 and 0.574.

## Housing at common states

The CSV also compares saved childless and one-child-at-home housing menus,
weighted by the same baseline \(W\), including tenure probabilities. Expected
physical rooms move from 4.635 to 4.652 in the childless menu and from 6.285 to
6.287 in the parent menu at .8 versus 1. Among the menu's renter choice mass,
baseline cap incidence is 5.970% versus 66.588%; parent conditional renter rooms
are 5.842 against a six-room cap. Conditional renter housing is **not** chosen
housing. Owner sales use \(b+S/R\), and cap indicators are interpolated over
the native forward-map nodes. The executed inputs use constant child maturation,
so newborn-exempt reoptimization is inactive; these saved parent policies also
represent the successful-first-birth branch. The exact check is in
`../diagnostics/joint_budget/README.md`.

Common live menu support excludes dead parent menus: baseline dead exposure is
0.000924, relaxed dead exposure 0.000874. Their susceptibility is zero. Native forward probability normalization is applied; live sums agree within \(3.331\times10^{-16}\). Realized baseline owner ending-debt-floor incidence is 4.974% overall and
11.536% among fertile childless owners (tolerance \(10^{-6}\), using distinct
new-purchaser and stayer saving policies). Hence origination access cannot
establish that credit is irrelevant. Source pins and checks are in `receipt.json`.

Reproduce from the project root:

```sh
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python output/model/credit_mechanism_20261004/census/census.py
```
