# Final-night quantitative work

**New author cutoff:12:30EDT (16:30UTC).** Access restored. Numerical jobs stop
at12:05EDT, collectors by12:25EDT; then slides only. These deadlines supersede
the older18UTC references below. Four conditional24-period policy checks are
queued behind exact six-period baseline/tax replays. They retain the same2023
households and preference; they do not refit the history. Contracts and job IDs
are under `cutoff_horizon/` and `jobs.json`.

Current evidence: September 13, 13:51 UTC. Fixed deadline: September 13, 18:00 UTC.
Cluster computation and collection are independent of the laptop. Latest account
check,13:52UTC:38%weekly remaining; keep the author’s20%floor.

**Access update: September 13, 14:33 UTC.** Local SSH now requires login renewal. Already
submitted jobs and collectors are independent; later progress has not been
verified. Do not read the snapshot below as a current queue observation.

## Discussion packet

The September 13 [terminal-method review](../../../docs/model/transition_terminal_method_review.md)
compares the current boundary with a transition between steady states, reviews
the official sequence-space package, and records measured timing and the
prerequisites for a clean single-shock experiment. No new model run was launched
for that review.

**`verified_history_readout.pdf` is the current 41-page report.** The first page
summarizes the decision points. Pages 2–8 contain the complete initial target
fit, historical fit, untargeted 2023 validation, paired tax tables, and parameter
restrictions. Pages 9–35 preserve all 17 native diagnostics for each of the
initial equilibrium, baseline policy, and higher-tax policy. Pages 36–41 show
the five historical/2023 comparison figures and supplemental policy comparison.

All 41 draft pages were rendered and visually inspected. The final version only
updates its actual timestamp. All 225 displayed source numbers checked;
`verified_history_readout_qa.json` binds the report to its inputs. Full-precision
CSV/JSON files remain the numerical source. `verified_initial_readout.pdf` is
an older initial-only readout.

## Established results

- Corrected initial equilibrium: two exact numerical repetitions, loss
  179.2984242480. This is a verified candidate from an incomplete search, not
  a converged optimizer. The probability normalization repair changes its
  loss by less than one millionth. Complete 13-row target and 17-row parameter
  tables: `corrected_initial/`. Nine free structural coordinates; annual beta
  capped at 0.99; targets and weights unchanged.
- All four ordinary/seeded A0/A+ six-period histories fit all four observed
  fertility windows. Ordinary A0 model/data: 1.973124/1.974875 (2008–2011),
  1.861042/1.861000 (2012–2015), 1.755338/1.755375 (2016–2019), and
  1.645559/1.645750 (2020–2023). Every accepted root passes the finite housing,
  PAYGO, equal-rebate and exact-replay gates.
- Native A0 2023 readout passes: maximum aggregate discrepancy 3.55e-15,
  including the dated first-birth room response. Full 13-row untargeted table:
  `history_A0_6/validation/validation_2023.csv`. Actual empirical vintages and
  measurement qualifications remain explicit.
- The A0 baseline and doubled-property-tax policy forecasts both converge from
  exactly the same inherited 2023 households, grid and supply curve. Separate
  pension and property-tax budgets balance. Native proof:
  `history_A0_6/policy_state_verification_v2.json` (job 17691079).

## Paired policy comparison

Both 1% and 2% annual property taxes return revenue equally per current head.
The comparison uses each policy's own 2023 forecast; it does not substitute the
2023 row from the earlier 2019 forecast. Six periods span 24 years; the final
2043 decision produces the 2044–2047 fertility flow.

| Change from 1% to 2%, both rebated | Initial decision 2023 | Final decision 2043 |
|---|---:|---:|
| Birth flow | −0.2478% | +0.4718% |
| Period TFR | −0.2212% | +0.6223% |
| Resident persons | 0% | +0.03251% |
| Physical rooms per head | −4.2367% | −4.8879% |
| Ownership | −1.4867 percentage points | +1.4263 percentage points |
| House price per room | −6.6489% | −7.6500% |
| Rent per room | +20.2015% | +14.1891% |

These are dated flow/level effects, not cumulative births or long-run stationary
comparisons. Supply is static-elastic with elasticity 0.63 in both policies;
its exact price/supply log ratio verifies 0.63. Lower prices therefore reduce
housing supply. No fixed-stock sensitivity has been launched or substituted.
All 66 dated comparison rows and 132 original level cells were checked against
the native outputs; see `history_A0_6/policy_comparison/qa.json`. The standard
policy graphs are under `history_A0_6/policies/*/graphs/standard_diagnostics/`.

A separate read-only comparison also verifies the ordinary A+ and seeded A0
pairs (`policy_sensitivity/summary.json`). A+ changes the initial birth flow by
−0.2631% and final flow by+0.5196%; the seeded A0 reproduces the corresponding
ordinary A0 effects within0.000002percentage points for births. These are
additional provisional checks; only A0 is in the main report.

## Outstanding work and limits

The six-period results remain **provisional**: horizon adequacy is unverified,
and production eligibility remains false. Longer 24- and 100-period fits are
still running. The 24-period tracks have not yet accepted their first fitted
window. The new 100-period tracks have passed their first two complete mappings;
these are valid evaluations, not converged equilibria. A0's second mapping
requires a safeguarded step; A+ improves. No tolerance was loosened.

Economic fit remains weak, especially housing. Untargeted 2023 completed
fertility is 1.69387 versus CPS2024 1.91842; capped mean rooms 6.25308 versus
ACS2023 5.58372; ownership ages 30–55 is 49.4514% versus 58.7409%; the first-birth
room response is 1.01368 versus retained 0.720246. See the complete table for all
moments, not just these examples. Physical rooms in policy graphs are distinct
from the capped empirical room measure.

A0 removes all post-2023 migration; A+ retains the supplied migration sensitivity.
Historical head-age conditioning through 2023 remains an imposed bridge. The
older orphan-care omission remains the explicitly approved fallback. B0/B+ are
not launched: converting child/person mass into new household heads remains
unresolved, and B+ additionally needs signed migrant-state allocation. The
0.5647956 diagnostic coefficient is not a production default. The optional
age-profile initial pilot did not improve its age component and was not promoted.

## Active calculations and numerical repair

| Job array | Work | Last established state |
|---|---|---|
| 17658836 | Ordinary A0/A+, six periods | Both histories and both paired policies complete |
| 17661737 | Price-seeded A0/A+, six periods | Both histories complete; A0 pair complete; A+ tax finishes |
| 17663940 | Ordinary A0/A+, 24 periods | First-window preference search |
| 17664449 | Price-seeded A0/A+, 24 periods | First-window preference search |
| 17686968 | Corrected boundary initialization, A0/A+, 100 periods | Two valid mappings each; root solving |
| 17676958 | Earlier 32-GiB 100-period starts | Independent old-helper attempt; see live receipts |
| 17680316 | Earlier extended 100-period starts | A0 failed premature audit; A+ independent attempt |

The long-boundary initializer previously audited inherited 2007 households at
terminal conditions before carrying the population. The isolated helper now
constructs lifetime household values without that substitute-population audit.
Every actual dated and carried-endpoint population still passes the same checks.
Native verification 17686318 reproduces policies, fiscal inputs, full accepted
six-period paths, terminal households and ledgers bitwise. Source:
`history_source_policy_seed_v1`; model kernels are unchanged. Immutable first-
and second-mapping evidence: `policy_seed_first100_mapping/` and
`policy_seed_second100_mapping/`. The new A0 first mapping matches all 303 old
coordinates and residuals exactly.

The 48-GiB duplicate array 17663986 was retired after actual 32-GiB mappings and
memory use were verified. Numerical cache results remain exact. One 100-period
mapping currently takes about 30 minutes, so completion by the deadline is
uncertain. Each track retains checkpoints, latest/best summaries, five-minute
heartbeats, bounded root/preference attempts and independent failure handling.
Automatic native readout collectors cover every active horizon; see `jobs.json`
and `readout_verification/`. Failed and superseded results are preserved.

## Scientific contract and regeneration

Preferences unexpectedly change at decisions 2007, 2011, 2015 and 2019. At each
vintage households expect that preference to persist; the solver computes the
corresponding price/pension/rebate path, matches the four-year fertility window,
and carries only the first realized period to the next surprise. Preferences
stay fixed after the final shock. The boundary values remaining lifetimes under
constant conditions, clearing markets and budgets on the actual carried endpoint
population. It does not reset to a stationary population or prove feasibility
forever beyond the horizon.

The fixed-preference migration and 6-versus-24-period comparisons remain separate
diagnostics under `fixed_preference_migration_comparison/` and
`fixed_preference_horizon_comparison/`. The latter changes first-window TFR by
−0.01424, which motivates the longer-horizon refit; it is not a horizon certificate.

Regenerate the supplemental comparison without a model solve:

```sh
code/model/.venv/bin/python code/model/tools/build_e5f_final_policy_readout.py --case-dir output/model/e5f_final_night_20260913/history_A0_6 --out output/model/e5f_final_night_20260913/history_A0_6/policy_comparison --state-verification output/model/e5f_final_night_20260913/history_A0_6/policy_state_verification_v2.json
```

The full report builder is `code/model/tools/build_e5f_final_night_report.py`.
Use the bundled Python runtime with reportlab and pypdf, supplying `--packet`,
`--history-case`, `--policy-dir`, `--policy-state-verification`, `--output`, and
an actual `--as-of` timestamp. The retained five historical figures use
`build_e5f_final_history_plots.py`; all 13 validation rows use
`build_e5f_final_history_validation.py`. Existing full-precision inputs are
sufficient; none of these readers runs the model.

Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/`.
The author-approved plan is `docs/model/e5f_two_closure_overnight_plan.md`.
`CALIBRATION_STATUS.md` is canonical; `jobs.json` records job IDs and provenance.
