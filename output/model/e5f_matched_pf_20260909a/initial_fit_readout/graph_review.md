# Focused review of the four stationary diagnostic graph packets

Read-only review, September 11, 2026. Scope: source definitions, six original
PNGs (age-30 policies, age-by-income ownership and age-30 wealth distributions
for `old_old` and `new_balanced`), and existing small summaries/CSVs for all four
cases. No model solve, cluster submission, checkpoint download or graph change.

**Finding:** the saved population-weighted ownership gradient across permanent
income groups is positive in every case. The age-30 owner-entry oscillations
are already visible in `old_old`; they are not introduced by the new utility
plus balanced pension. Whether the oscillations materially affect occupied
decision states is still unresolved. Existing numerical receipts do not settle
that question.

## What the two graphs measure

1. `ownership_by_age_income_state.png` is an actual-distribution statistic.
   For age \(j\) and combined earnings state \(z\), it divides owner mass in
   `evaluation.g_current` by total mass at that same age and state, summing
   over wealth, location, parity and child state. `g_current` is **after**
   fertility and tenure/location choices. It is not an unweighted policy curve
   or a controlled income comparative static. Groups with mass at most
   \(10^{-14}\) are omitted. The age/type cells need checkpoint extraction
   for exact numeric comparisons; the collected lifecycle CSV pools all types.
2. Each \(z\) label is the product of a permanent level and a changing
   five-state earnings component. The source uses three permanent levels
   `0.2774924094, 0.8219540209, 2.4346915070` with weights
   `1/6, 2/3, 1/6`. Within each permanent group it lists five changing states;
   permanent types never switch (`kron(I_3, Pi_base)`). Therefore the 15 legend
   entries are not globally income-sorted: `0.584903` is followed by `0.328467`,
   and `1.732527` by `0.972944`. The changing component is a Markov state,
   not an independent one-period transitory shock. Moreover the 15 plotted
   lines cycle through only ten colors: low/high permanent groups reuse blue,
   orange, green, red and purple. An upper blue line cannot be identified as
   the first, low-income blue legend entry solely from color.
3. `policy_childless_renter_age30.png` conditions on age 30, inherited renter
   tenure, location 0, parity 0 and the readiness-settled child state. It plots
   **every** wealth node satisfying `V > -1e9`, without any mass filter. Its
   owner-entry line is the sum of saved tenure probabilities for owner products.
   In this sequential model tenure is chosen after fertility, so the relevant
   exposure for the childless tenure decision is the matching slice of
   `g_post_fertility`, not simply `g_pre` or `g_current`. The policy graph is
   conditional on the childless outcome, not an ex ante birth-weighted owner
   probability. The saved wealth-distribution graph instead uses `g_current`
   and normalizes each type separately; it cannot by itself measure how much
   pre-tenure mass encounters an owner-entry reversal.
4. The housing panel uses `tenure_choice` to select either a discrete owner
   product or the conditional renter housing policy. It is not expected housing
   weighted by tenure probabilities. Its consumption panel directly uses the
   renter slice of `c_pol`. These plotted panels should not be read as one
   common realized consumption/housing bundle. The reducer exports both
   conditional renter housing and the graph's selected-tenure housing values.
5. Large owner-entry reversals across wealth are visually present in both
   `old_old` and `new_balanced`, especially over roughly 0–15 model wealth
   units. This is an observation, not a diagnosis of either coding error or
   valid discrete-product economics. All four receipts report zero occupied
   adjacent-wealth value drops at the \(10^{-7}\) gate, zero household budget
   violations, and probabilities within [0,1]. Those checks do not require
   ownership probabilities to be monotone or establish global policy
   optimality. Discrete owner products are a candidate mechanism to test,
   not a demonstrated explanation.

## Verified existing numbers

Population-weighted ownership at ages 30–55, percent, by permanent income
level; directly from each repetition summary's
`permanent_income_own_rate_3055_by_level`:

| Case | Low permanent income | Middle | High |
|---|---:|---:|---:|
| old_old | 9.954666 | 65.426921 | 95.658722 |
| old_balanced | 4.869347 | 62.783881 | 95.497067 |
| new_old | 10.557487 | 65.598430 | 95.611206 |
| new_balanced | 5.196397 | 62.851672 | 95.462838 |

Ownership rises with the changing earnings state within each permanent group
in the **all-age** saved summary for every case. This does not assert an
identical ordering in every age cell. Across all 15 combined states sorted by
current earnings, a single monotone ordering need not hold because the permanent
type and earnings history also differ.

Selected all-type age ownership rates, percent, taken from `lifecycle_2023.csv`
(the inherited filename does not make these stationary initial cases a 2023
historical observation):

| Case | Age 30 | Age 42 | Age 62 | All ages |
|---|---:|---:|---:|---:|
| old_old | 38.098455 | 64.277034 | 85.472167 | 64.370972 |
| old_balanced | 36.598272 | 61.644835 | 80.970131 | 62.057774 |
| new_old | 38.354294 | 64.483734 | 85.532334 | 64.546638 |
| new_balanced | 36.774314 | 61.703037 | 80.998777 | 62.145477 |

These four cases hold the preference intercept and structural coordinates fixed
and each solve its own stationary housing equilibrium. They are a utility/pension
decomposition, not re-estimated early calibrations. The two old-pension cases are
explicitly fiscally unbalanced diagnostics. No target fit or loss is inferred here.

## One bounded next check, prepared but not submitted

`reduce_graph_checkpoints.py` in this directory is a read-only reducer for the
two exact `old_old` and `new_balanced` checkpoint hashes in the four-cell receipt.
It validates the pinned contract/source files and both checkpoint SHA256 values
before unpickling. It configures the sequential runtime only; it calls no
household solver, equilibrium routine, population operator or plotting routine.

Lead can run it in the immutable c6dd3508 snapshot with one CPU, 16 GB and an
external five-minute cap, writing a new output directory. The script accepts
only `--output`; the snapshot, cases and hashes are frozen inside it. It loads
the two checkpoints sequentially and writes six CSVs and a hash receipt:

- All wealth nodes, both ages 30/42 and all 15 income states: the exact graph
  mask, values, consumption, conditional renter housing, selected tenure product,
  graph housing, owner probability and every destination probability; matching
  pre-fertility, post-fertility/pre-tenure and current/post-tenure node masses.
- Adjacent valid-node owner-probability drops, both unweighted and evaluated at
  positive lower-node mass, for each of the three mass timings. Thresholds are
  explicit: probability drop \(10^{-7}\), occupancy \(10^{-12}\). A screen flag
  is diagnostic and does not change any acceptance gate.
- Exact ownership numerator, denominator and ratio by every age and combined
  earnings state, plus genuinely population-weighted aggregation by permanent
  income group. These are the data behind the age/type interpretation.

The primary comparison should use post-fertility/pre-tenure mass for tenure
policy exposure; pre-fertility mass remains the relevant timing for the existing
value screen. Current mass is retained to reproduce the original ownership
and wealth-distribution graphs. If occupied probability reversals remain,
inspect their saved destination probabilities and product switches before
choosing a further Bellman comparison. No such new solve is proposed here.

Validation of the reducer: Python AST parsing passed locally. The two large
checkpoints remain on Torch, so runtime extraction and its numeric readout are
pending. No cluster submission was made by this review.

## Source evidence

All paths below are relative to the project root unless noted.

- Original four-case evidence: `output/model/e5f_matched_pf_20260909a/utility_fiscal_decomposition/`,
  specifically `contract.json`, `case_summary.csv`, and each
  `collected/<case>/repetition_01/{summary.json,policy_array_summary.json,lifecycle_2023.csv,standard_diagnostics/summary.json}`.
- Exact plot definitions: `tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/diagnostics.py`,
  lines 147–190 and 291–374. SHA256
  `3d77fc4ebddad9635407c65cd73f1add3cea8fd2f4a8250d125329c9113a654b`
  matches the frozen four-cell contract.
- Exact plotting adapter: `tmp/e5f_matched_pf/code/model/tools/run_e5f_independent_numerical_audit.py`,
  lines 107–164 and 248–276. SHA256
  `656d766e0197e854c73e6b6fb7a3e97cd4fef57e75ad4c7a924d8b1ca5212528`
  matches the frozen contract.
- Income-group construction: `tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/e6b_profile.py`,
  lines 18–66. Distribution ordering:
  `tmp/e5f_matched_pf/code/model/tools/run_dynamic_population_transition.py`,
  `evaluate_period`, lines 469–527.
- The current probe driver has changed since c6dd3508 during the parent's active
  implementation. This review relies on the contract-matching plotting source
  and saved outputs, not an assumption that the current driver is the old snapshot.

