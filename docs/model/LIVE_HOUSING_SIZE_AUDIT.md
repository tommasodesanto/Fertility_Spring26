# Live Housing-Size Audit

Date: `2026-04-18`

## 2026-04-18 Addendum

This note originally documented the older April 13 housing audit. The live
room-fit reference has now changed.

Current room-audit convention:

- sample: **MMS large metros only**, using the existing
  [puma_mms_lookup_2020.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv)
  join on `statefip`, `puma`, `met2013`
- geography reconciliation: canonical `middle -> center`
- room bins: **split `5` and `6`**, not pooled `5-6`

Current live comparison object:

- prime-age renters `25-45`, all:
  - ACS MMS: `<=4 = 65.1%`, `5 = 16.0%`, `6 = 10.2%`, `7-8 = 6.5%`
  - model: `<=4 = 1.9%`, `5 = 98.1%`, `6 = 0.0%`, `7-8 = 0.0%`
- prime-age owners `25-45`, all:
  - ACS MMS: `5 = 17.1%`, `6 = 20.0%`, `7-8 = 27.6%`, `9-10 = 12.6%`, `11+ = 6.5%`
  - model: `5 = 6.5%`, `6 = 24.6%`, `7-8 = 68.9%`, `9-10 = 0.0%`, `11+ = 0.0%`

Interpretation:

- restricting to the actual large-metro sample makes the comparison more
  appropriate for the project
- splitting `5` from `6` makes the renter-cap bunching much clearer
- the owner-side fit improved somewhat under the current `H0` experiment
- the renter-side compression remains the main housing miss

Everything below this addendum is older context from the April 13 audit and
should not be treated as the current room-fit reference by default.

This note is for the **live April discrete-time benchmark**, not the older March planning files.

## 0. Executive Summary

This is the shortest accurate read of the housing block tonight.

### Main conclusion

The current live benchmark should **not** be trusted for housing-policy interpretation yet.

Why:

- it is too high in **physical rooms** relative to the ACS, not just for families but also for childless households
- the weak starter-home policy effects therefore look **calibration-driven at first order**
- the birth-related housing moments `H01` and `H12` were not enough to validate the physical housing distribution

### What is clearly not the main issue

- the existence of an `11`-room top rung is **not** the first-order problem
- ACS owners do have a real `11+` tail
- the bigger issue is that the benchmark puts far too much mass in `8+` rooms and too little mass on the lower and middle support

### What is clearly wrong in the live benchmark

- all-household owner median rooms:
  - ACS `6`
  - model `8.2`
- all-household renter median rooms:
  - ACS `4`
  - model `7.89`
- prime-age owners age `25-45`, child bins `0/1/2+`:
  - ACS `6 / 6 / 7`
  - model `9.6 / 9.6 / 9.6`
- prime-age renters age `25-45`, child bins `0/1/2+`:
  - ACS `4 / 4 / 5`
  - model `7.86 / 8.00 / 8.31`

### What this means for the parameter interpretation

- the family-space schedule is **too aggressive relative to the ladder**
- but that is **not the whole problem**
- the benchmark is already allocating too much housing before the family thresholds bite

So the current miss is best read as:

- a broadly too-high physical housing allocation
- amplified by a too-steep family-space schedule

### What I would do

I would **freeze the current benchmark as policy-invalid**, keep it as a diagnostic anchor, and move to a disciplined acceptance rule:

1. Require the benchmark to match the ACS room distribution at least coarsely.
2. Fix the non-family housing allocation first.
3. Then re-evaluate the family-space schedule.
4. Only then read starter-home policy.

## 1. Live Objects To Audit

Separate the housing block into three different objects:

- `H_own`: physical owner size ladder
- `chi`: owner-specific quality / service premium
- `h_bar_0`, `h_bar_jump`, `h_bar_n`: family-space need schedule

Do not treat these as interchangeable.

- Physical owner size is `H_own`
- Effective owner services are `chi * H_own`
- Child housing need is `h_bar`
- The physical owner size required to clear the threshold is `h_bar / chi`

In the live code, the benchmark setup still carries:

- `H_own = linspace(4.0, 11.0, 6)` in [build_calibration_setup.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/build_calibration_setup.m:89)
- `hR_max = 8.0` in [build_calibration_setup.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/build_calibration_setup.m:93)
- `h_bar_0 = 4.0` in [build_calibration_setup.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/build_calibration_setup.m:72)

Owner services are implemented as:

- `hsrv(i,ten) = P.chi * hs` in [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:439)

Family housing need is implemented as:

- `h_bar = h_bar_0 + h_bar_jump + h_bar_n * nk` in active parent states in [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:333)

## 2. Live Benchmark Mapping

Using the live benchmark values you flagged:

- `H_own = [4.0, 5.4, 6.8, 8.2, 9.6, 11.0]`
- `hR_max = 8.0`
- `h_bar_0 = 4.0`
- `h_bar_jump = 2.3`
- `h_bar_n = 1.0`
- `chi = 1.09`

The implied thresholds are:

| State | Service threshold `h_bar` | Physical owner threshold `h_bar/chi` | Smallest owner rung that clears |
| --- | ---: | ---: | --- |
| childless | 4.00 | 3.67 | `H1 = 4.0` |
| one child | 7.30 | 6.70 | `H3 = 6.8` |
| two children | 8.30 | 7.61 | `H4 = 8.2` |
| three children | 9.30 | 8.53 | `H5 = 9.6` |

That means:

- the first child effectively wipes out `H1` and `H2`
- the second child effectively wipes out `H1`, `H2`, and `H3`
- bottom-rung owner policy experiments are mechanically weak for families

This is the central live concern.

## 3. What Is Already Data-Disciplined

### PSID housing response moments

The repo already contains direct housing-response evidence:

- first birth baseline event-study: `+0.664` rooms at `+3`, `+0.843` at `+5` in [rooms_f_c_y_all_summary.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv:1)
- second birth baseline stay-at-one control event-study: `+0.566` rooms at `+3`, `+0.576` at `+5` in [rooms_s_c_y_all_summary.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all_summary.csv:1)
- second birth last-treated robustness: `+0.704` at `+3` in [rooms_s_c_y_all_summary.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_second_birth_lasttreated_controls_v1/rooms_s_c_y_all_summary.csv:1)
- second birth tighter robustness: `+0.845` at `+3` in [rooms_s_c_y_no_third_by3_gap5_summary.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_no_third_by3_gap5_summary.csv:1)

These moments discipline:

- the total first-child housing response
- the marginal additional-child housing response

They do **not** separately identify:

- the raw owner ladder
- the owner premium `chi`
- the decomposition of first-birth response into `h_bar_jump` versus `h_bar_n`

### Tenure moments

The live calibration setup also uses:

- `own_rate`
- `own_gradient`
- `own_family_gap`

in [build_calibration_setup.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/build_calibration_setup.m:17)

These partially discipline `chi`, but only after the physical menu and rental cap are judged coherent.

## 4. What Is Still Effectively Ad Hoc

### `H_own`

Still largely a physical menu choice, not directly estimated from a dated ACS target file.

### `hR_max`

Still a segmentation device. The older planning note explicitly described the rental cap as state-space engineering in [Calibration_Plan_Merged.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/Calibration_Plan_Merged.tex:126).

### `h_bar_jump`

Still not separately locked. The repo already says this in [first_birth_housing_target_april.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/archive/april_guides_2026-05-07/first_birth_housing_target_april.md:37):

- the data are informative about first-birth housing response
- but the mapping from that response to `h_bar_jump` is not one-to-one

### `h_bar_0`

This is functioning like a normalization, but with the live scale at `4.0` it still needs a childless-size check against ACS/PUMS.

## 5. Important Repo Hygiene Point

Do not rely on the older merged calibration plan as the live benchmark description.

That file still describes:

- `h_bar_0 = 0.50`
- `h_bar_jump = 0.50`
- `H_own = [4, 7.5, 11]`
- `hR_max = 4.5`

in [Calibration_Plan_Merged.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/Calibration_Plan_Merged.tex:120).

Those are not the current live benchmark objects.

## 6. Tonight's Audit Workflow

### A. Audit `H_own`

Use 2024 ACS:

- `B25020` Tenure by Rooms
- `B25021` Median Number of Rooms by Tenure
- `B25042` Tenure by Bedrooms

Build:

- owner room distribution
- owner bedroom distribution
- owner medians and upper quartiles
- the same by households with `0`, `1`, and `2+` children if using PUMS

Goal:

- decide whether the physical ladder itself is plausible
- do not let `chi` fix a bad ladder

### B. Audit `hR_max`

Use 2024 ACS/PUMS renter upper-tail size distribution.

Build:

- renter p90 and p95 rooms
- renter p90 and p95 bedrooms
- share of renters with children above candidate caps

Goal:

- justify the rental cap as a real tail restriction
- or conclude it is mostly a stylized segmentation wedge

### C. Audit `chi`

Use ACS tenure facts conditional on size.

Build:

- `P(own | rooms bin, bedrooms bin, children bin, age 25-45)`

Goal:

- keep `chi` as a residual tenure/services wedge
- not as a substitute for missing family-sized physical units

### D. Audit `h_bar_0`, `h_bar_jump`, `h_bar_n`

Use:

- ACS/PUMS childless size distribution for `h_bar_0`
- PSID first-birth housing response for `h_bar_jump + h_bar_n`
- PSID second-birth housing response for `h_bar_n`

Goal:

- back out whether the live thresholds are too high relative to the physical ladder

## 7. Minimum Output For Tonight

Produce four objects:

- one calibration table listing live values, implied thresholds, and smallest clearing owner rung
- one empirical table listing ACS tenure-size moments and PSID birth-housing moments
- one diagnostic plot of ACS size distributions with ladder and cap overlaid
- one diagnostic plot of threshold schedule versus owner ladder

## 8. Working Prior Before Running ACS

The live setup probably overstates the first-child family-space threshold relative to the ladder.

Why:

- the first-child threshold jump in owner physical units is about `3.03` rooms
- the second-child threshold jump in owner physical units is about `0.92` rooms
- the PSID realized housing responses are much smaller than that raw first threshold jump
- the live threshold map forces one child to `H3` and two children to `H4+`

That does not prove the ladder is wrong.
It points more directly to the need schedule being too aggressive relative to the ladder.

## 8A. What The ACS Distribution Audit Changed

The direct ACS room-distribution comparison sharpened the diagnosis:

- the live benchmark misses the **overall physical room distribution**, not just the birth-related housing response moments
- prime-age childless owners and prime-age childless renters are already too large in physical rooms relative to the ACS
- so the family-space schedule is part of the problem, but it is not the whole problem

That means:

- `h_bar_jump` and `h_bar_n` still matter for family sorting
- but the room-distribution miss cannot be blamed only on the child-need schedule
- the benchmark is broadly placing too much mass on high housing states

## 9. ACS References For Tonight

Use official Census sources:

- ACS 2024 detailed table groups: <https://api.census.gov/data/2024/acs/acs1/groups.html>
- ACS 2024 PUMS overview: <https://api.census.gov/data/2024/acs/acs1/pums.html>
- ACS 2024 PUMS variables: <https://api.census.gov/data/2024/acs/acs1/pums/variables.html>
- ACS 2024 `B25115` tenure by household type / own children: <https://api.census.gov/data/2024/acs/acs1/groups/B25115.html>

For PUMS, the key variables are:

- `TEN`
- `RMSP`
- `BDSP`
- `WGTP`
- `NRC`
- `HUPAOC`
- `PAOC`

## 10. Local Data Constraint For Tonight

The local bridge extract at:

- [processed_acs_data_slim.rds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/data/processed_acs_data_slim.rds)

is useful for some background ACS work, but it is **not** the clean same-night source for this audit because:

- it runs through `2023`, not `2024`
- it contains `bedrooms`, but not `rooms`
- it is not already packaged as a clean one-record-per-housing-unit tenure-size target file for this purpose

So for the housing-size audit itself:

- use official 2024 ACS detailed tables for the room objects
- use 2024 ACS PUMS only if you need the tenure × children × bedroom cross-tabs
- do not confuse the local slim bridge extract with the final audit source

## 11. ACS Room-Distribution Audit

The direct ACS room-distribution comparison is now complete using the full ACS household extract on disk:

- [fertility_microdata_2023.rds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/processed_data/yearly_rds_v7/fertility_microdata_2023.rds)

Sample:

- household-head records only: `pernum = 1`, `relate = 1`
- occupied units only: `gq = 1`, `ownershp in {1,2}`
- household weights: `hhwt`
- prime-age comparison: householder age `25-45`
- child bins: `0`, `1`, `2+` using household-head `nchild`

Raw outputs:

- ACS summary:
  [acs_2023_room_distribution_audit.txt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/acs_2023_room_distribution_audit.txt)
- model summary:
  [model_room_distribution_audit.txt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/model_room_distribution_audit.txt)
- side-by-side comparison:
  [room_distribution_compare_report.txt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/room_distribution_compare_report.txt)
- prime-age figure:
  [room_distribution_compare_25_45.png](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/room_distribution_compare_25_45.png)

Headline result:

- the live benchmark does **not** target room-distribution moments directly
- the implied room distribution is much too large relative to the ACS

Top-line all-household medians:

- owners: ACS median `6` rooms, model median `8.2`
- renters: ACS median `4` rooms, model median `7.89`

Prime-age (`25-45`) medians by tenure and current-child bin:

- owners, `0/1/2+` children:
  - ACS: `6 / 6 / 7`
  - model: `9.6 / 9.6 / 9.6`
- renters, `0/1/2+` children:
  - ACS: `4 / 4 / 5`
  - model: `7.86 / 8.00 / 8.31`

So the live benchmark is too high in physical rooms broadly, not just at the first-birth and second-birth margins.

## 12. What The Model Actually Uses

The new support-usage diagnostic shows where the benchmark puts its mass:

- [model_support_usage.txt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/model_support_usage.txt)

Prime-age (`25-45`) benchmark usage:

- childless owners:
  - `H1 = 0%`
  - `H2 = 0%`
  - `H3 = 9.0%`
  - `H4 = 23.6%`
  - `H5 = 65.0%`
  - `H6 = 2.4%`
- one-child owners:
  - `H1-H3 = 0%`
  - `H4 = 46.2%`
  - `H5 = 38.3%`
  - `H6 = 15.6%`
- two-plus-child owners:
  - `H1-H3 = 0%`
  - `H4 = 12.3%`
  - `H5 = 87.7%`
  - `H6 = 0%`

Prime-age renters:

- childless renter cap share: `34.2%`
- one-child renter cap share: `100%`

This is the crucial distinction:

- the **support** is not the first-order problem
- the **mass placement on the support** is

In particular:

- an `11+` owner tail is not itself absurd in the ACS
- ACS owners overall have about `6.3%` in `11+` rooms
- the model has about `7.0%` there overall

So `H6 = 11` existing is not the main issue.
The main issue is that the benchmark places far too little mass on the lower and middle parts of the support and far too much mass in `8+` rooms.

## 13. What Is Actually Going Wrong

There are no obvious algebraic or coding mistakes in the housing mapping itself.

The actual mistakes are calibration / validation mistakes:

- we let the benchmark pass without checking the ACS room distribution
- we asked `H01` and `H12` to do too much
- we interpreted policy nulls from a benchmark where the relevant support barely carries mass

More concretely:

- `H01` and `H12` are dynamic response moments
- they do **not** tell you whether the benchmark sits at the right physical housing levels
- the live benchmark is too large in physical rooms even in childless states
- then the family-space schedule pushes family states even higher

So the current weak starter-home policy effects should be read first as a benchmark-quality problem, not as deep structural evidence.

## 14. What We Should Be Doing Differently

The practical change is straightforward:

1. Treat ACS room-distribution fit as a required benchmark screen.
2. Separate:
   - matching **levels / distributions** of physical rooms
   - matching **birth-related changes** in rooms
3. Reject any benchmark that has:
   - prime-age renters piled at the cap
   - prime-age owners barely using lower / middle rungs
   - childless households already far too large in rooms
4. Only after that should starter-home and family-rung policies be interpreted.

This is not “recalibrate everything.”
It is a disciplined benchmark acceptance rule.

## 15. Object-Level Read After The ACS Audit

### `H_own`

- the existence of `H1-H3 = 4.0, 5.4, 6.8` means the ladder does include smaller units
- the problem is that the benchmark barely uses them in prime ages
- so the first-order issue is not “the ladder lacks small houses”
- it is “the benchmark allocates too little mass to the small and middle rungs”

### `hR_max`

- `hR_max = 8` is not innocuous in the live benchmark
- it is a live pile-up point
- that is especially clear because prime-age one-child renters are at the cap `100%` of the time

### `h_bar_0`

- childless renters and childless owners are already too large relative to the ACS
- so the miss is not only a first-child or second-child threshold problem

### `h_bar_jump`, `h_bar_n`

- they still matter because they push family states even further into `8+` rooms
- but they are amplifying an already-high physical housing allocation
- they are not the sole source of the room-distribution miss

### `chi`

- `chi` remains a plausible contributor because the owner side is too large even before children
- but this audit alone does not isolate `chi` cleanly from the broader tenure / housing-demand mapping

## 16. Current Bottom Line

The cleanest conclusion now is:

- the current live benchmark is too high in physical housing across the board
- the family-space schedule then pushes family mass even higher
- the benchmark therefore makes lower owner rungs and moderate rentals too irrelevant
- weak starter-home policy effects look calibration-driven at first order

So if the question is:

- “Are we making actual mistakes?”

the answer is:

- **yes, mainly calibration-screening mistakes**

and if the question is:

- “Is the current benchmark physically too high relative to the real ACS room distribution?”

the answer is:

- **yes**

## 17. What I Would Do

If the goal is to get to a benchmark you can actually interpret, I would do the following in this order.

### Step 1. Mark the current benchmark as failed on physical housing levels

Do **not** discard it.
Keep it as the reference point that taught us the failure mode.

But do not treat it as policy-ready, because:

- it fails the ACS room distribution badly
- it piles renters at the cap
- it leaves too little mass on the lower and middle owner rungs

### Step 2. Add a hard ACS room-distribution screen

At minimum, require the benchmark to be roughly right on:

- all-household owner median rooms
- all-household renter median rooms
- prime-age owner median rooms by child bin `0/1/2+`
- prime-age renter median rooms by child bin `0/1/2+`
- prime-age owner and renter shares in coarse bins:
  - `<=4`
  - `5-6`
  - `7-8`
  - `9-10`
  - `11+`

If a candidate misses these badly, do not read policy on it.

### Step 2A. Change the live calibration architecture only where needed

Do **not** add a huge new block of housing moments to the loss.

The disciplined Stage A move is:

- keep the current core behavioral targets
- keep `H01` hard
- keep `H12` only as a **soft** target for now
- add exactly two ACS housing-level moments:
  - prime-age childless renter median rooms = `4`
  - prime-age childless owner median rooms = `6`
- free `h_bar_0` in a narrow band
- keep `H_own` fixed on this first pass

The reason is simple:

- the room-distribution miss is already visible for childless households
- so the first repair has to discipline **baseline physical housing allocation**
- not just the family-jump margins

On the current live benchmark, these two new moments are badly missed already:

- prime-age childless renter median rooms: model `7.86`, ACS target `4`
- prime-age childless owner median rooms: model `9.6`, ACS target `6`

### Step 3. Fix the non-family housing allocation before touching family policy interpretation

The benchmark is already too large for childless households.

So I would **not** start by saying:

- “just lower `h_bar_jump`”

That would be incomplete.

I would first ask:

- why are childless renters already near `8` rooms?
- why are prime-age childless owners already concentrated at `H4-H5`?

This points first to the **baseline housing allocation**, not just the family schedule.

### Step 4. Only after Step 3, revisit the family-space schedule

Once childless and prime-age non-family room distributions are in range, then look again at:

- `h_bar_jump`
- `h_bar_n`

At that point the question becomes clean:

- conditional on a reasonable physical distribution, are children still pushing households too high too fast?

Right now that question is contaminated by the benchmark already being too high before children.

### Step 5. Re-run the starter-home policy only after the benchmark passes the screen

Right now the nulls are not very informative because:

- the benchmark barely uses the lower owner rungs
- the benchmark already pushes renters to the top of the rental support

So I would read the current policy nulls as:

- evidence about a **bad benchmark**

not:

- evidence that the policy is structurally weak in the model

## 18. What I Would Not Do First

I would **not** start by:

- deleting `H6 = 11`
- recalibrating everything at once
- reading the current starter-home nulls as structural
- adjusting only `h_bar_jump` while ignoring the childless room miss

## 19. Practical Decision Rule

A housing benchmark should be treated as usable only if all three are true:

1. It is roughly right on the ACS room distribution.
2. Lower and middle owner rungs carry meaningful prime-age mass.
3. `H01` and `H12` are matched without pushing almost everyone into `8+` rooms.

The current live benchmark fails that test.

## 20. Stage A Directional Search

After wiring the two ACS childless room medians into the live objective, I ran fast saved-anchor probes around the live benchmark region.

The useful baseline fast anchor is:

- `h_bar_0 = 4.0`, `chi = 1.09`
- childless renter median rooms: `7.66`
- childless owner median rooms: `8.50`
- `H01 = 0.333`
- `H12 = 0.664`
- total loss on the Stage A fast surface: `252.6`

### What lowering `h_bar_0` does

Holding `chi = 1.09`:

- `h_bar_0 = 3.0`
  - renter median falls to `7.07`
  - owner median stays at `8.50`
  - `H01 = 0.677`
  - `H12 = 0.245`
  - total loss falls to `92.3`
- `h_bar_0 = 2.75`
  - renter median falls to `6.33`
  - owner median stays at `8.50`
  - `H01 = 0.798`
  - `H12 = 0.161`
  - total loss falls to `72.8`
- `h_bar_0 = 2.5`
  - renter median falls to `6.13`
  - owner median stays at `8.50`
  - `H01 = 0.898`
  - `H12 = 0.137`
  - total loss is `77.9`

Interpretation:

- lowering `h_bar_0` clearly moves the renter-side baseline room miss in the right direction
- it also strongly moves `H01`
- but it does **not** fix the owner-side baseline-size problem by itself
- and if pushed too low, it overshoots `H01` and crushes `H12`

### What lowering `chi` does

At `h_bar_0 = 3.0`:

- `chi = 1.07`
  - renter median `7.04`
  - owner median `6.00`
  - `H01 = 0.688`
  - `H12 = 0.240`
  - own rate `0.442`
  - total loss `93.0`
- `chi = 1.05`
  - renter median `7.01`
  - owner median `6.00`
  - `H01 = 0.676`
  - `H12 = 0.234`
  - own rate `0.398`
  - total loss `96.3`

At lower `h_bar_0` values, the owner median does **not** flip:

- for `h_bar_0 = 2.75`, the owner median stays at `8.50` across `chi = 1.03, 1.05, 1.07, 1.09`
- for `h_bar_0 = 2.5`, the owner median stays at `8.50` across `chi = 1.03, 1.05, 1.07, 1.09`

Interpretation:

- the owner-side baseline-size problem is **nonlinear and discrete**
- it is not a smooth “lower `chi` a bit and the owner median drifts down” margin
- instead there is a regime switch near the `h_bar_0 = 3.0` row where lowering `chi` flips the owner median from `8.5` to `6.0`
- but that owner-side fix is too blunt by itself because it drives ownership too low and still leaves `H12` weak

### What this means

The Stage A directional search clarified the mechanism:

- renter-side baseline housing is mainly tied to `h_bar_0`
- owner-side baseline housing is not just `h_bar_0`; it is tied to a discrete owner-attractiveness margin
- the old benchmark was hiding these two problems inside a single housing block

So the next real search should not restart from the generic Stage A `x0`.
It should start from the live benchmark region with:

- `h_bar_0` lower than `4.0`
- `chi` explored carefully around the owner-side regime switch
- the understanding that fixing the owner median and fixing ownership are not the same thing

### Full benchmark confirmations run on 2026-04-13 evening

I re-ran the most informative Stage A overrides on the full six-house benchmark setup, using the saved live benchmark anchor rather than the lighter `fast` surface.

Results:

- `h_bar_0 = 2.75`, `chi = 1.07`
  - total loss `49.583`
  - `H01 = 0.638`
  - `H12 = 0.126`
  - renter childless median rooms `6.761`
  - owner childless median rooms `8.200`
  - own rate `0.299`
  - old-age parent-childless gap `0.0026`
- `h_bar_0 = 3.00`, `chi = 1.07`
  - total loss `51.174`
  - `H01 = 0.546`
  - `H12 = 0.143`
  - renter childless median rooms `6.933`
  - owner childless median rooms `8.200`
  - own rate `0.306`
  - old-age parent-childless gap `0.0058`
- `h_bar_0 = 2.75`, `chi = 1.09`
  - total loss `55.104`
  - `H01 = 0.613`
  - `H12 = 0.193`
  - renter childless median rooms `6.780`
  - owner childless median rooms `8.200`
  - own rate `0.371`
  - old-age parent-childless gap `-0.0083`

Interpretation:

- the Stage A direction survives the full benchmark surface
- lowering `h_bar_0` and modestly lowering `chi` is not a fake `fast`-mode result
- among these simple two-parameter overrides, `h_bar_0 = 2.75`, `chi = 1.07` is the best benchmark confirmation
- but the repair is still incomplete:
  - renter baseline rooms remain far above the ACS target `4`
  - owner baseline rooms remain far above the ACS target `6`
  - ownership collapses too much
  - `H12` remains far too low

So the current conclusion is sharper than before:

- `h_bar_0` really does need to come down
- `chi` likely needs to come down somewhat too
- but a pure two-parameter Stage A override is not enough to make the housing block benchmark-valid
- the owner-side baseline-size problem is still not solved

### Final Torch scout read on 2026-04-13 evening

The two-hour Torch Stage A scout finished. It did not produce a clean no-penalty benchmark candidate, but the checkpoint ranking is still informative.

Best checkpoint:

- loss `10107.875`
- `beta = 0.9389`
- `b_entry_fixed = 0.2567`
- `psi_child = 0.1152`
- `h_bar_jump = 2.3097`
- `h_bar_n = 0.8751`
- `c_bar_n = 0.1259`
- `kappa_fert = 2.5897`
- `chi = 1.0595`
- `kappa_loc = 1.9303`
- `mu_move = 0.0499`
- `theta0 = 0.5409`
- `theta_n = 0.2357`
- `h_bar_0 = 2.6826`

Interpretation:

- the full benchmark still prefers `h_bar_0` well below `4.0`
- it also prefers `chi` below `1.09`
- the best full-objective region is not the old live benchmark neighborhood
- but the scout still sits on a hard guardrail, so this is not yet a valid replacement benchmark

### Torch-best local benchmark frontier

I also took the best Torch checkpoint region and re-ran a small local benchmark frontier directly off the saved live anchor. This is not a PSO; it is a focused full-benchmark comparison around the Torch point.

Completed cases:

- Torch-style candidate with `chi = 1.07`, `h_bar_n = 0.95`
  - total loss `118.622`
  - `H01 = 0.664`
  - `H12 = 0.560`
  - renter childless median rooms `6.615`
  - owner childless median rooms `8.200`
  - own rate `0.442`
  - old-age parent-childless gap `-0.0885`
- Torch-style candidate with `chi = 1.0595`, `h_bar_n = 0.95`
  - total loss `113.688`
  - `H01 = 0.835`
  - `H12 = 0.598`
  - renter childless median rooms `6.691`
  - owner childless median rooms `5.400`
  - own rate `0.446`
  - old-age parent-childless gap `-0.0652`
- Torch-style candidate with `chi = 1.07`, `h_bar_n = 0.8751`
  - total loss `104.232`
  - `H01 = 0.706`
  - `H12 = 0.150`
  - renter childless median rooms `6.443`
  - owner childless median rooms `5.400`
  - own rate `0.442`
  - old-age parent-childless gap `-0.0133`

Interpretation:

- the Torch-centered region can get very close to the dynamic housing moments
- but it does so by breaking the old-age parent-childless gap and still leaving childless room levels far too high
- the simple Stage A benchmark override `h_bar_0 = 2.75`, `chi = 1.07` remains better on total benchmark loss than these Torch-style variants
- raising `h_bar_n` restores `H12`, but the cost is a more negative old-age gap and worse geography

### Five-rung boundary-shift probe on 2026-04-14 morning

The next structural test was the one implied by the overnight segmentation work: do not rely on the rental cap alone, and instead shift the owner boundary upward while dropping the `11`-room rung from the active search structure.

Tested owner ladder:

- `H_own = [4.5, 6.0, 7.2, 8.8, 10.5]`

Tested caps and seed parameters:

- `hR_max = 6.8` or `7.2`
- `h_bar_0 = 2.75`
- `h_bar_jump = 1.90`
- `h_bar_n = 0.90`
- `chi = 1.07`, plus one lower-`chi` check at `1.05`

Benchmark results:

- `five68_seed`
  - total loss `76.046`
  - own rate `0.590`
  - childless renter median rooms `6.711`
  - childless owner median rooms `8.800`
  - `H01 = 0.330`
  - `H12 = 0.960`
  - starter-margin shares from renter birth-state parents:
    - `H1_H2 = 0.0558`
    - `H3+ = 0.3812`
- `five72_seed`
  - total loss `60.158`
  - own rate `0.467`
  - childless renter median rooms `6.659`
  - childless owner median rooms `8.800`
  - `H01 = 0.355`
  - `H12 = 0.767`
  - starter-margin shares:
    - `H1_H2 = 0.0328`
    - `H3+ = 0.2392`
- `five68_chi105`
  - total loss `78.611`
  - own rate `0.539`
  - childless renter median rooms `6.711`
  - childless owner median rooms `8.800`
  - `H01 = 0.363`
  - `H12 = 0.918`
  - starter-margin shares:
    - `H1_H2 = 0.0318`
    - `H3+ = 0.2748`

Mechanical-bite interpretation:

- This is the first structure that creates a nontrivial `H2` wedge on the relevant parent margin without an extreme rental cap.
- In the center, the first owner switch at age `32` is now to `H2`, with:
  - `age32_first_owner_switch_b = 12.85`
  - `age32_first_H1_H2_choice_b = 12.85`
  - `age32_first_H3plus_choice_b = 17.97`
- In the periphery, the first owner switch still goes straight to `H3+`.
- The policy still does not work through financing feasibility:
  - the birth-state renter mass is already above the baseline `H1/H2` down-payment thresholds
  - the new wedge comes from boundary support placement, not from relaxing the down-payment constraint

Bottom line:

- Lowering the rental cap alone was not enough to make the starter-owner margin relevant.
- Shifting the owner ladder near the boundary does make the tenure angle show up mechanically.
- But this five-rung structure is not yet a usable benchmark:
  - childless room levels remain far too high relative to ACS
  - `H01` is too low
  - `H12` is too high
  - childless owner median rooms remain stuck at `8.8`

So the structural lesson is sharper now:

- the current six-rung `[4.0, ..., 11.0]` ladder places the entry-family boundary badly
- stronger tenure segmentation requires boundary-support changes, not just a lower rental cap
- but the first five-rung boundary-shift tested here still needs a proper recalibration before it can be treated as an economically credible benchmark

### Boundary-overlap probe centered closer to ACS `5-6` rooms

I then tested a more data-centered five-rung ladder:

- `H_own = [4.5, 5.8, 7.0, 8.5, 10.0]`

with moderate rental caps and milder housing-need parameters:

- `hR_max = 6.4` and `6.0`
- `h_bar_0 = 2.75`
- `chi = 1.05`
- `h_bar_jump = 1.80`
- `h_bar_n = 0.80`

Completed benchmark cases:

- `bo64_seed`
  - total loss `62.287`
  - own rate `0.635`
  - childless renter median rooms `6.400`
  - childless owner median rooms `8.500`
  - `H01 = 0.341`
  - `H12 = 0.601`
- `bo60_seed`
  - total loss `68.057`
  - own rate `0.805`
  - childless renter median rooms `6.000`
  - childless owner median rooms `8.500`
  - `H01 = 0.364`
  - `H12 = 0.575`

Interpretation:

- This support is more in line with the ACS room boundary, and it improves the dynamic housing moments relative to some earlier structural probes.
- But it kills the starter-owner margin again:
  - `H1_H2 = 0` in both completed cases
  - the first owner switch at age `32` goes directly to `H3+`
  - there is no `H2` entry region left

So the lesson is now more precise:

- moving the boundary down toward the ACS room bins helps fit
- moving it down too far collapses the tenure-entry margin we were trying to recover
- the relevant structural frontier is therefore **between** the earlier `H = [4.5, 6.0, 7.2, 8.8, 10.5]` probe and this more ACS-centered `[4.5, 5.8, 7.0, 8.5, 10.0]` probe

### Wide boundary-margin search on 2026-04-14 evening

I then widened the search around the five-rung frontier and ran a two-hour cluster search over:

- `lambda in [0.70, 0.95]`
- `hR_max in [6.45, 6.90]`
- `chi in [1.050, 1.075]`
- `h_bar_0 in [2.50, 2.85]`
- `h_bar_jump in [1.95, 2.25]`
- `h_bar_n in [0.60, 0.82]`

The objective was also adjusted to reflect what the audit had learned:

- `H01` received more weight
- `H12` received less weight
- candidates were penalized if the `H2` starter-family margin collapsed

Best wide-run candidates:

- `task_10`
  - `lambda = 0.756`
  - `hR_max = 6.859`
  - `chi = 1.058`
  - `h_bar_0 = 2.508`
  - `h_bar_jump = 2.075`
  - `h_bar_n = 0.645`
  - `own = 0.454`
  - `mR0 = 6.567`
  - `mO0 = 8.727`
  - `H01 = 0.378`
  - `H12 = 0.683`
  - `H1_H2 = 0.037`
- `task_14`
  - `lambda = 0.894`
  - `hR_max = 6.869`
  - `chi = 1.070`
  - `h_bar_0 = 2.506`
  - `h_bar_jump = 2.067`
  - `h_bar_n = 0.691`
  - `own = 0.508`
  - `mR0 = 6.541`
  - `mO0 = 8.768`
  - `H01 = 0.361`
  - `H12 = 0.776`
  - `H1_H2 = 0.064`

Interpretation:

- The wider search confirms a stable neighborhood:
  - `h_bar_0` low, around `2.5-2.6`
  - `chi` around `1.06`
  - `h_bar_jump` around `2.0-2.1`
  - `h_bar_n` around `0.65-0.70`
  - `hR_max` around `6.7-6.9`
- But the main tradeoff remains:
  - better `H01` and `H12` push toward a weaker `H2` starter margin
  - keeping the `H2` margin alive still leaves owner childless rooms far too high

This means the search is no longer drifting. It has localized the structural region, but it has **not** found a policy-ready benchmark.

### Current interpretation of the starter margin

The key additional diagnostics were run on the lead wide-run candidates.

First, the ownership-entry report shows that the relevant new-parent renter mass is **not** down-payment constrained in the current best candidates:

- in both `task_10` and `task_14`, essentially all screened renter new-parent mass can already afford some owner rung
- large shares still rent anyway
- for `task_10`:
  - periphery: `feasible_any_owner_share = 1.000`, `feasible_any_owner_but_rent_share = 0.652`
  - center: `feasible_any_owner_share = 1.000`, `feasible_any_owner_but_rent_share = 0.926`
- for `task_14`:
  - periphery: `1.000` feasible, `0.606` still rent
  - center: `1.000` feasible, `0.875` still rent

So the current wedge is **not** “wants to buy but cannot make the down payment.”

Second, the size-need report shows that the relevant starter-family rung is `H2`, not `H1`:

- in `task_10`, the birth-state physical thresholds are:
  - parity-1 parent: `4.94`
  - parity-2 parent: `5.55`
  - parity-3 parent: `6.16`
- `H1 = 4.5` is therefore physically too small for the birth-state parent margin
- `H2 = 5.95` is physically adequate for essentially all of the screened mass:
  - periphery: `H2_physically_adequate_share = 0.990`
  - center: `H2_physically_adequate_share = 0.998`
- but most of that `H2`-adequate mass still rents:
  - periphery: `H2_adequate_but_rent_share = 0.651`
  - center: `H2_adequate_but_rent_share = 0.924`

So the clean current interpretation is:

- `H1` is not a family starter house in the current calibrated region
- `H2` is the relevant starter-family rung
- the current barrier is crossing into ownership at `H2`, not down-payment feasibility and not lack of physical space at `H2`

This is a real ownership-related wedge, but it is **not** a strong finance-binding wedge.

A direct value-gap diagnostic on the wide-run leader `task_10` sharpens this further:

- In the **periphery**, among states where `H2` is feasible and physically adequate:
  - `Rent - H2 = +1.15` on average
  - `H3 - H2 = +0.59` on average
  - so `H2` is not the relevant owner package there; conditional on owning, the model tends to prefer the larger rung
- In the **center**, among states where `H2` is feasible and physically adequate:
  - `Rent - H2 = +0.44` on average
  - `H2 - H3 = +1.48` on average
  - so `H2` is the relevant starter-family owner package there, but some renter new-parent mass still prefers renting to crossing into ownership at `H2`

This is exactly the distinction needed for interpretation:

- In the periphery, the issue is partly a **package mismatch** (`H2` versus `H3+`)
- In the center, the issue is a genuine **ownership wedge at `H2`**

So the model is no longer saying only “starter houses are too small” or “down payments block entry.”
It is saying:

- `H1` is too small
- `H2` is the relevant starter-family rung
- and, at least in the center, there is a real ownership-package wedge even when `H2` is feasible and large enough

### Presentation-ready bottom line

For presentation purposes, the cleanest current message is:

- The ACS supports a real but moderate tenure-size margin.
- The relevant empirical boundary is around `5-6` rooms, not `8`.
- The live benchmark is invalid because it puts both renters and owners too high in rooms.
- Cap-only segmentation does not work.
- A shifted five-rung boundary can create a real `H2` starter-family margin.
- In the best current candidates, the wedge is ownership-related, but not down-payment-related:
  - `H1` is too small
  - `H2` is large enough
  - many new-parent renters can buy `H2` and still choose to rent

So the current model evidence supports a **tenure-size / ownership-package** mechanism much more than a **strict down-payment-constraint** mechanism.

### Policy Battery From The Mechanism-Preserving Wide Candidate

I then ran a policy battery from the current mechanism-preserving wide candidate (`best_task_14`).

Baseline for that candidate:
- `tfr = 2.321`
- `mean_age_first_birth = 26.80`
- `own = 0.508`
- `H01 = 0.361`
- `H12 = 0.776`
- renter new-parent starter share `H1/H2 = 0.064`

Completed policy cases so far:
- `10% parents`
- `0% parents`
- `birth_grant`
- `10% H2 at birth`
- `0% H2 at birth`
- `0% center H2 at birth`
- `10% H2-H3 at birth`
- `0% H2-H3 at birth`
- `0% center H2-H3 at birth`

#### What the battery says

The generic parent-only down-payment waivers do **not** preserve the starter-family ownership margin:

- `10% parents`
  - `tfr = 2.439`
  - `mean_age_first_birth = 26.23`
  - `own = 0.290`
  - `H01 = 0.739`
  - `H12 = 0.069`
  - `H1/H2 share = 0.000`
- `0% parents`
  - `tfr = 2.454`
  - `mean_age_first_birth = 26.19`
  - `own = 0.300`
  - `H01 = 0.823`
  - `H12 = 0.057`
  - `H1/H2 share = 0.000`

So those broad finance relaxations do raise fertility and shift birth timing earlier, but they do so while collapsing ownership and eliminating the `H2` starter-family margin.

The narrow `H2` and `H2-H3` zero-down policies are also not producing a clean starter-home channel:

- `0% H2 at birth`
  - `tfr = 2.425`
  - `own = 0.312`
  - `H01 = 0.746`
  - `H12 = 0.066`
  - `H1/H2 share = 0.000`
- `0% H2-H3 at birth`
  - `tfr = 2.425`
  - `own = 0.312`
  - `H01 = 0.746`
  - `H12 = 0.066`
  - `H1/H2 share = 0.000`

The center-only variants are essentially the same. So the current model is **not** saying “targeted zero-down on the relevant starter rungs unlocks family ownership.”

The one finance policy that creates a large direct owner-entry response is the `birth_grant` case:

- `birth_grant`
  - `tfr = 2.298`
  - `mean_age_first_birth = 25.97`
  - `own = 0.260`
  - `H01 = 0.305`
  - `H12 = 0.315`
  - renter new-parent `H1/H2 share = 0.342`
  - renter new-parent rent share falls from `0.752` to `0.346`
  - center age-32 first `H1/H2` choice falls from `12.85` to `0.89`

This is a real finance response, but it is also a **very strong** intervention. In the current code, `birth_grant` effectively bypasses the down-payment feasibility test for birth-state entrants into ownership, across owner rungs, rather than just making `H2` slightly easier to buy.

#### Policy interpretation

- The model can produce strong finance-driven action, but only under the much more aggressive `birth_grant` design.
- The targeted down-payment-waiver policies are not currently producing the clean starter-ownership story.
- So in the present calibration region, the ownership wedge is real, but it is **not** well summarized as a simple starter-rung down-payment barrier.
- The honest policy reading at this stage is:
  - the data support a moderate tenure-size boundary,
  - the model supports an ownership-related wedge,
  - and strong finance action appears only under the much more aggressive birth-grant intervention.

The two supply-style cases are still pending on the cluster and are not yet incorporated here.
