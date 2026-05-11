# Housing Moment Mapping

Date: `2026-04-13`

This note clarifies the exact mapping between:

- the empirical housing moments `H01` and `H12`
- the model objects `h_bar_0`, `h_bar_jump`, `h_bar_n`
- the owner quality premium `chi`
- the physical housing menus `H_own` and `hR`

## 1. What The Empirical Target Is

The data moment is a **realized average housing change in physical units**:

- `H01_data = E[h_post - h_pre | first birth treated cohort]`
- `H12_data = E[h_post - h_pre | second birth treated cohort / control design]`

In the current live target sheet:

- `H01 = 0.66443467`
- `H12 = 0.56581378`

These are rooms.

They are **not** direct estimates of `h_bar_jump` or `h_bar_n`.

## 2. What The Model Moment Is

The model computes the event-style housing moments using **actual physical housing**, not services:

- renters contribute actual rental housing `hR`
- owners contribute actual physical owner size `H_own`

This is explicit in [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:1363):

- renter contribution: `hR_pol`
- owner contribution: `P.H_own(ten-1)`

So the model-side `H01` / `H12` are physical-housing moments, directly comparable to rooms in the data.

## 3. Where `h_bar` Enters The Model

`h_bar` is **not** chosen housing.

It is the Stone-Geary housing subsistence floor:

- active-parent housing need is set in [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:333)
- owner services are `chi * H_own` in [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:439)

So the within-period utility problem is built around:

- renters: effective housing surplus `hR - h_bar`
- owners: effective housing surplus `chi * H_own - h_bar`

That is why:

- `h_bar` is a **need threshold**
- `H01` / `H12` are **realized choice responses**

These are different objects.

## 4. Exact Renter Mapping

For an unconstrained interior renter, the code implies:

- total resources after savings choice: `Rv - bp`
- subsistence outlay: `c_bar + r * h_bar`
- remaining surplus is split in Cobb-Douglas shares

From [run_model_cp_dt.m](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/root_scripts/run_model_cp_dt.m:570):

- `ct = alpha * surplus`
- `ht = (1-alpha)/r * surplus`
- actual housing is `h = h_bar + ht`

Substituting gives:

\[
h^R_* = h_{\text{bar}} + \frac{1-\alpha}{r}\left(Rv - bp - c_{\text{bar}} - r h_{\text{bar}}\right)
\]

which simplifies to:

\[
h^R_* = \alpha h_{\text{bar}} + \frac{1-\alpha}{r}(Rv - bp - c_{\text{bar}})
\]

So, holding savings fixed and staying in the unconstrained renter interior:

\[
\frac{\partial h^R_*}{\partial h_{\text{bar}}} = \alpha
\]

With the live calibration:

- `alpha = 0.70`

So a `+1.0` increase in `h_bar` does **not** raise realized renter housing by `+1.0`.
It raises it by about `+0.70`, before any extra savings response.

This is already enough to show why a `0.6` room data moment does **not** map one-for-one into a `0.6` change in `h_bar`.

## 5. Exact Owner Mapping

For owners, physical housing is discrete:

- actual physical housing is the chosen rung `H_k`
- effective services are `chi * H_k`

At a fixed owner rung:

- actual physical housing does not change when `h_bar` changes
- only utility changes through the surplus term `chi * H_k - h_bar`

So on the owner side, `h_bar` affects physical housing only indirectly:

- by changing whether renting or owning is preferred
- by changing which owner rung becomes worth paying for

To keep the same effective housing surplus at a fixed utility margin, the owner-side physical translation is:

\[
\Delta H_{\text{phys, equiv}} = \frac{\Delta h_{\text{bar}}}{\chi}
\]

This is a **threshold-equivalent physical change**, not an average realized moment.

## 6. Why A `+0.6` Data Moment Can Coexist With A Much Larger `h_bar` Shift

Because the empirical moment is:

- an average over many heterogeneous households
- after endogenous savings, location, and tenure adjustment
- measured in physical housing

while `h_bar` is:

- an internal subsistence shift in utility
- partly absorbed through consumption and savings
- partly absorbed through no adjustment
- partly absorbed through discrete tenure jumps

In the model, a large `h_bar` shift can generate a modest average physical response if:

- many households cannot or do not adjust
- some households hit the rental cap
- some households stay on the same owner rung
- only a subset switch from rent to own or move up the ladder

So it is normal that:

- `H01_data` is much smaller than `h_bar_jump + h_bar_n`

The key question is whether it is **too much smaller**.

## 7. What The Live Benchmark Implies

Under the current canonical local benchmark:

- `h_bar_jump = 2.3`
- `h_bar_n = 1.0`
- `chi = 1.09`
- `alpha = 0.70`

The first-child increase in `h_bar` is:

\[
\Delta h_{\text{bar},01} = h_{\text{bar jump}} + h_{\text{bar n}} = 3.3
\]

Three different translations matter:

### A. Raw subsistence shift

- `3.3` service units

### B. Owner threshold-equivalent physical shift

\[
3.3 / 1.09 = 3.03
\]

This is where the “extra three rooms” number comes from.

It is **not** the model housing moment.
It is the physical owner-size increase needed to keep the same surplus at the threshold.

### C. Interior renter choice translation

\[
0.70 \times 3.3 = 2.31
\]

This is the approximate increase in actual renter housing if the same renter stays interior and savings do not change.

Again, this is **not** the final event-study moment.
It is the within-regime partial-equilibrium translation.

## 8. Why The Current Benchmark Still Looks Too Aggressive

The live benchmark gives:

- model `H01 = 0.622`
- model `H12 = 0.420`

while the benchmark internal need shift implies:

- owner threshold-equivalent first-child shift `= 3.03`
- interior-renter first-child partial mapping `= 2.31`

So the current benchmark is using a very large internal first-child need shock to generate a realized average response of about `0.62`.

That is the core coherence concern.

In plain language:

- the model can mechanically compress a big need shock into a smaller observed room response
- but the current first-child shock looks very large relative to the resulting average physical adjustment

## 9. What The Live `H01` Is Actually Coming From

Using the canonical local benchmark decomposition in
[h01_birth_cohort_decomposition.txt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/h01_birth_cohort_decomposition.txt):

- saved `H01 = 0.622`
- recomputed `H01 = 0.626`

So the decomposition is matching the production event-study moment.

The key benchmark facts are:

- `68.5%` of the first-birth treated cohort is renting before the birth
- that renter group contributes `0.609` of the total `0.626` `H01`
- only `3.3%` of the treated cohort is below the one-child threshold before the birth
- `96.7%` is already above that threshold before the birth

Post-birth at `t+3`, the treated cohort is:

- `23.7%` capped renters at `hR = 8.0`
- `27.4%` owners at `H4 = 8.2`
- `35.2%` owners at `H5 = 9.6`
- `13.7%` owners at `H6 = 11.0`
- essentially `0%` interior renters

This clarifies the compression:

- the model is **not** producing `H01` through a smooth interior-renter response close to `alpha * Δh_bar = 2.31`
- instead, it is producing `H01` mostly through pre-birth renters whose average physical housing rises by only `0.889`
- by `t+3`, those households are sitting on discrete margins: rental cap or owner rungs `H4-H6`

So the benchmark `H01` does **not** strongly identify the one-child threshold level itself.
Most of the treated cohort is already above the threshold in effective services before the birth, and the realized response is being generated on flat or discrete post-birth margins.

That is the core mapping result:

- `Δh_bar` is large internally
- realized `H01` stays modest because the treated cohort is mostly operating on capped/discrete housing margins rather than a smooth renter interior

## 9. Practical Rule Going Forward

Do **not** read:

- `H01 = 0.66`

as implying:

- `h_bar_jump + h_bar_n = 0.66`

That mapping is false.

The correct rule is:

- choose `h_bar_jump`, `h_bar_n`, `chi`, `H_own`, and `hR_max`
- solve the full model
- compute model-side event-style physical housing moments
- compare those to the empirical `H01` and `H12`

Then separately check whether the implied threshold mapping:

- `h_bar / chi`

is interpretable relative to the physical ladder.

## 10. Bottom Line

The “extra three rooms” is not a bug in the moment definition.

It is the owner-side threshold-equivalent translation of the current first-child subsistence shift:

\[
(h_{\text{bar jump}} + h_{\text{bar n}})/\chi
\]

That object is conceptually different from the empirical `+0.6` average realized room response.

The reason this matters is not that those two numbers should be identical.
It is that, in the current live benchmark, the internal first-child need shift may be too large relative to the moderate realized room response it produces.
