# Empirical Roundup Assessment

## Bottom Line

The earlier MSA-based audit overstated what the current PSID workflow can do.

The repo already has a coherent empirical package, but it is split across two
data systems:

- `mms_center_periphery` carries the **spatial** evidence.
- `psid_followup_mar2026` carries the **birth, housing, ownership, and wealth**
  evidence.

So the right assessment is not “most key predictions are untested.” It is:

- the **sorting** side is already well covered in ACS
- the **housing-step / ownership / moves-for-space** side is already well
  covered in PSID
- the truly missing pieces are:
  - parent-vs-childless housing elasticities
  - a clean slack/nonlinearity test
  - a direct wealth gradient of **spatial** family flight

## Corrected Audit Table

| Prediction | Theory source | Status | Assessment | Best evidence |
|---|---|---|---|---|
| P1a. Parent share falls with local rent, conditional on income | `simple_spatial.tex` | Partially tested | The repo strongly supports the broader sorting fact that parents are less central, but I do **not** see the exact conditional-on-income rent-slope regression on disk. | `../code/data/mms_center_periphery/output_middle_center/mms_location_by_parent_summary.csv` |
| P1b. First birth triggers relocation toward lower-price areas | `simple_spatial.tex` | Partially tested | There is strong indirect support: ACS movers with children choose peripheral destinations within metros, and PSID shows birth raises moves for space. But there is no direct PSID “birth -> cheaper metro/area” event study. | `mms_origin_transition_summary_allparent.csv`, `mms_move_regressions_allparent.csv`, legacy `mv_s_f_c_y_all.png` |
| P2a. Parents have lower housing income elasticity | `simple_spatial.tex` | Untested but feasible | This is not on disk yet, but the current PSID data do support a rooms/income or ownership/income specification by parental status. | Potentially `PSIDSHELF_MOBILITY.dta` via `INCFAMR`, `ACTUALROOMS_`, `HOMEOWN` |
| P2b. Parents have lower housing price elasticity | `simple_spatial.tex` | Untested and currently not feasible in the existing workflow | No live PSID metro panel or cleaned price merge is in the current codebase. | None in current pipeline |
| P3a. Family flight is decreasing in pre-birth wealth | `simple_spatial.tex` | Untested spatially | Wealth clearly matters for ownership transition and fertility, but not yet for spatial relocation at birth. | `within_wealth_gap_v9.csv`, `dual_rates_quintile_v9.csv` |
| P3b. Wealth threshold above which parents stay central | `simple_spatial.tex` | Untested | No direct spatial wealth-threshold exercise is on disk. | None |
| L3a. First child causes a discrete housing-demand jump | `claude_slides.tex`, broader theory block | Tested, slide-ready | This is one of the cleanest PSID facts in the repo. | `rooms_f_c_y_all_summary.csv`, `rooms_first_birth_one_vs_two_horizon_v1.csv` |
| L3b. Low slack makes the first-child housing cost sharply nonlinear | `claude_slides.tex` | Untested | The theory is explicit, but there is no finished empirical slack-hazard or low-slack spline design on disk. | None |
| L3c. Relaxing the ownership / family-space constraint raises fertility | main deck, wealth scripts | Strongly supported, but mostly descriptive / mediation rather than final causal proof | This is already the strongest reduced-form constraint-margin fact in the repo. | `within_wealth_gap_v9.csv`, `fertility_wealth_v9.log`, `childless_by_transition_v7_bar.png` |

## What Is Already Solid

### 1. ACS spatial sorting facts are already there

The current center/periphery build already supports the main spatial claims used
in the deck.

From `../code/data/mms_center_periphery/output_middle_center/`:

- `mms_location_by_parent_summary.csv`
  - non-parents center share: `0.494`
  - new parents center share: `0.416`
- `mms_move_regressions_allparent.csv`
  - parent indicator lowers within-CBSA center destination by about `11.4 pp`
  - parent indicator lowers across-CBSA moving by about `5.2 pp`
- `mms_origin_transition_summary_allparent.csv`
  - among center-origin movers, center destination share:
    - non-parents `0.670`
    - parents `0.580`
  - among periphery-origin movers, center destination share:
    - non-parents `0.374`
    - parents `0.263`

This is already enough to support:

- parents are less central
- parents re-sort mainly within metros
- within-metro movers with children are more peripheral in destination choice

### 2. PSID first-birth housing-demand facts are already there

From `../code/data/psid_followup_mar2026/output/`:

- `sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv`
  - `+0.664` rooms at `k = +3`
  - `+0.843` rooms at `k = +5`
- `rooms_first_birth_one_vs_two_horizon_v1/rooms_first_birth_one_vs_two_horizon_v1.csv`
  - one kid by +2: `+0.259`
  - two-plus by +2: `+1.368`
- `sa_rooms_first_birth_grouped_v1/rooms_f_c_y_one_kid_by3_summary.csv`
  - one-kid path by +3: `+0.289`
- `sa_rooms_first_birth_grouped_v1/rooms_f_c_y_two_plus_by3_summary.csv`
  - two-plus path by +3: `+1.031`
- `sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all_summary.csv`
  - second birth: `+0.566` rooms at `k = +3`

So the discrete family-space step is already very well established.

### 3. The ownership-constraint result is already central and usable

From `../code/data/psid_followup_mar2026/output/fertility_wealth_v1/within_wealth_gap_v9.csv`
and the corresponding script/log:

- among renters at 25-30, those who become owners by 35 are much less likely to
  remain childless by 45
- the gap is present in every wealth quintile
- the v9 mediation output shows:
  - total wealth effect on fertility: small / insignificant (`p = 0.619`)
  - wealth effect conditional on ownership: essentially zero (`p = 0.917`)
  - ownership coefficient in fertility regression: about `0.321`, highly
    significant

That is why the main deck says the effect runs through ownership. On the current
evidence, that is a fair characterization of the repo’s reduced-form story.

### 4. The hump-shaped wealth-fertility relationship is already documented

From:

- `lpoly_parenthood_v12.png`
- `parenthood_by_quintile_v12.csv`

This is already slide-ready and is important because it blocks an overly simple
“more wealth -> more fertility” interpretation.

## What Exists But Is Not Yet Clean

### 1. First-birth hazard by wealth

These files currently collapse to zeros:

- `first_birth_hazard_by_tercile.csv`
- `first_birth_hazard_by_quintile.csv`

So they should **not** be treated as evidence until debugged.

### 2. Wealth-transition-after-birth outputs

`wealth_v1` is interesting, but I would treat it as exploratory:

- `wealth_tercile_outcomes_v1.csv`
- `wealth_transition_regressions_v1.tex`

The regression table is informative, but the raw tercile summary is not yet
something I would elevate into the core narrative without another pass.

## What The Earlier Audit Got Wrong

### 1. It treated missing PSID metro identifiers as if they were available

That is the biggest error. The current PSID shelf workflow is not set up for a
direct “move to lower-rent MSA” design.

### 2. It understated how much is already done

The repo already has:

- ACS spatial sorting
- ACS within-metro mover resorting
- PSID first-birth rooms jump
- PSID second-birth rooms jump
- PSID ownership-constraint / within-wealth gap
- PSID hump-shaped wealth-fertility
- PSID moves-for-size in the legacy event-study stack

That is already a substantial empirical package.

### 3. It blurred “exists on disk” with “cleanly exportable”

Some things are truly slide-ready.
Some are only exploratory.
Some are broken.
Those should be separated.

## Best Next Tests Given The Real Data

### Best next test: wealth heterogeneity in first-birth housing adjustment

This is the highest-return extension that actually matches the current data.

Run first-birth event studies interacted with pre-birth wealth for:

- ownership
- rooms
- moves for size

This would test whether lower-wealth households make larger housing adjustments
at birth even if you cannot observe metro-level “flight.”

### Second-best test: descriptive slack gradient

Use pre-birth rooms or renter slack bins to estimate a descriptive first-birth
hazard gradient.

This is not yet a full structural validation of the power-law mechanism, but it
is the closest empirical step toward Lemma 3 with the current data.

### Third-best test: income elasticity attenuation

Use `INCFAMR` and `ACTUALROOMS_` in PSID to estimate a parental-status
interaction.

That is a real feasible extension.

### Not recommended right now: price elasticity or cheaper-MSA flight

Those both require geography / price objects the current PSID workflow does not
actually use.

## Practical Conclusion

If I had to summarize the repo in one sentence:

The paper already has good empirical support for **family-space activation** and
for **spatial sorting in ACS**, but it does **not** yet have a direct reduced-form
test of price-elasticity attenuation, slack nonlinearity, or wealth-moderated
spatial flight.

## Files Added In This Turn

- [EMPIRICAL_ROUNDUP_PROMPT.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/EMPIRICAL_ROUNDUP_PROMPT.md)
- [EMPIRICAL_ROUNDUP_ASSESSMENT.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/EMPIRICAL_ROUNDUP_ASSESSMENT.md)
