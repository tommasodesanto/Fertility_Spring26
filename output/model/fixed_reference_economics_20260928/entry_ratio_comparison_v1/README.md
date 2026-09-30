# Five-ratio entrant diagnostic — completed, small changes at fixed parameters

## Verified result — September 30

Torch **18888956 completed successfully in 18m28s**, with seven lifecycle
solves in each arm, including exact repeats. The five-ratio/current-income
initialization changes the main fitted moments only slightly at the inherited
parameter vector. Weighted loss rises from **29.477 to 30.030 (+1.876%)**.
This is a fixed-parameter sensitivity result, not a re-estimated calibration.

| Equilibrium quantity | Retained entry | Five-ratio entry | Change |
|---|---:|---:|---:|
| Childlessness | 20.203% | 20.215% | +0.012 percentage points |
| Ownership, ages 30–55 | 61.890% | 61.874% | −0.017 percentage points |
| Mean age at first birth | 25.806 | 25.797 | −0.009 years (about 3 days) |
| Wealth / earnings | 6.150 | 6.145 | −0.005 |
| Mean rooms | 5.764 | 5.762 | −0.003 |
| Old-age wealth p90/p50 (validation) | 3.068 | 3.125 | +0.057 |

Equilibrium price rises **0.089%** and household population **0.100%**.
Completed fertility is approximately 2.1 in both arms because price clears birth
renewal; that equality is not independent validation. At the same prescribed
reference price, completed fertility rises by **0.001 children**, ownership
ages 30–55 rises **0.013 percentage points**, and childlessness falls
**0.018 percentage points**. All six checked seed policy arrays (value, saving,
consumption, rental housing, tenure probabilities, fertility probabilities)
are exactly equal: only initial probability weights differ at that price.

The test supports the economic transparency of the five-ratio rule without
evidence of a major disturbance to the main fitted moments. It does not prove
that re-estimated parameters would remain unchanged, or certify grid convergence,
transitions, the income-proxy assumption, or the empirical independence assumption.
The older-wealth dispersion diagnostic moves more than the main fitted moments
(about 1.9%); it remains below its 3.516 target. Neither initialization is being
selected by which gives the lower loss. No candidate or credit rule is adopted.

- [Full equilibrium comparison: all 14 targets, weights, roles, gaps and loss contributions](collected/full/comparison_target_fit.csv).
- [Full same-price comparison: all 14 rows](collected/full/comparison_reference_price_target_fit.csv).
- [All 31 parameter values, restrictions and bound flags](collected/full/candidate_five_ratios/phase_b_ge/selected_root/parameters.csv). Every row is identical across arms; the inherited fertility shock-scale parameters retain their near-bound flags. No parameters were estimated in this test.
- [Initial distribution comparison](collected/full/initial_entry_distribution.png), [control standard diagnostics](collected/full/control_fixed_reference/phase_b_ge/selected_root/standard_diagnostics/), and [candidate standard diagnostics](collected/full/candidate_five_ratios/phase_b_ge/selected_root/standard_diagnostics/). All 17 standard plots per arm are retained, unchanged in definition. Regeneration uses the runner's inherited `ge.observe_price(..., final=True)` with the saved native solution and authenticated context; `launch_torch.sh` regenerates the complete two-arm run and packets.
- [Lead verification](collected/verification.json), [137 downloaded artifact hashes](collected/remote_hash_receipt.json), [native completion](collected/full/completed.json), and [exact standard-plot repeat hashes](collected/full/standard_plot_repeat_hashes.json).

The lead independently checked all 137 downloaded hashes, complete 14-row fit
and 31-row parameter equality within each exact repeat, all 17 selected/repeat
PNG hashes per arm, zero actual calendar-entry mapping error, and zero forward
feasibility projection. The control reproduces the prior D=0.53 target table
exactly. Ownership and fertility lifecycle plots and candidate policy plots
were visually inspected. Existing estate-counterparty limitations remain as
recorded in native gate receipts; this test does not resolve them.

## Specification and preparation record

The lead submitted **18888956** after independent zero-solve smoke, formula and
calendar-entry review, five-file archive authentication, and remote SHA check.
The exact submitted archive is retained as `deployed_stage.tar.gz`, SHA
`4a0381f6d0f66d52f46c05124cd61442bc499ca734cf847fb2a86bccf57c3853`.
`launch.json` records the run; the verified result is reported above.

Provenance clarification: the worker prepared a later six-file archive
(`c319cea586dab734f4400dd9dbfd4a27439ce5977732cdf4d69a1b90c3f69766`)
as the earlier archive was being deployed. It was **not** submitted. Numerical
driver and launcher bytes are identical. The later archive adds explicit
objective/config pins and a plan; those checks already existed in the submitted
workflow: frozen observer authentication hashes the objective and target-weight
fingerprint before solving, and the launcher authenticates the config in the
immutable 81-file source inventory. Use the deployed archive for exact replay.

Both arms use the matched, authenticated 160 wealth × 15 income-state indexed engine, fixed D = 0.53, corrected full-sale/no-taper credit rules, fixed preferences and fiscal/supply inputs. Price clears actual birth renewal, and household population clears physical housing supply. This is an isolated diagnostic, not recalibration, adoption, or transition validation.

Control uses the exact retained `fixed_reference_entry_conditional`. Candidate draws the retained **five** empirical ratio nodes independently of current income, sets raw wealth to ratio times **current annual gross entrant income**, and calls the native linear grid projection separately at every income node. It installs this verified conditional matrix with `native_fixed_reference_entry=True` so the frozen calendar observer uses the same law. Only that matrix changes. The full empirical observation distribution is not an alternative arm.

The native finite-support approximation clips raw wealth outside the retained grid before projection: approximately 0.0183% of draw probability. Raw candidate mean wealth is 0.18651967924681828, projected mean is approximately 0.18663361710229673. This changes the very small out-of-support tail; the test does not preserve raw debt support exactly. No extra forward feasibility relocation is allowed. Native mass, feasibility, estate, budget, PAYGO, renewal, parameter, and reporting gates remain binding.

`run_comparison.py` imports the immutable matched workflow and unchanged `phase_b_grid.py` under `output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner/`. Each arm authenticates the frozen checkpoint **before** replacing entry. The lead's completed D53 control could otherwise be reused, but it lacks a complete 14-row prescribed-reference-price fit table. The new pair computes that useful comparison from each necessary price-seed solution using zero extra lifecycle calls; those seed results are explicitly not GE. This avoids a saved-state reconstruction adapter.

Budget: one Torch CPU, 24 GiB, 2400 seconds globally from launcher entry, 300 seconds per lifecycle/observer case, 20 lifecycle calls per arm (40 total), including the seed and exact final repeat. No automatic retry or extension. Prior 160×15 workflow used seven calls and about 469 seconds, so two comparable arms plausibly take about 16 minutes, conditional on convergence; the deadline remains binding.

Local exact-loop smoke is `smoke_pinned/completed.json`: both fresh-process arms traverse the original mocked price loop with zero lifecycle solves. Input checks authenticate unchanged sources, verify exactly five ratio nodes, all non-entry fields, identical grids/income, native lookup agreement, probability normalization and necessary current entrant cash feasibility. Two earlier preparation failures are retained under `smoke_v1` (pin-generation path issue) and `smoke_v2` (Numba cache created arm directory); both were repaired without any model solves. This smoke does not authenticate the heavyweight frozen observer or certify native equilibrium.

Full outputs retain all 14 target rows, 31 parameter rows, 17 standard plots for each arm and exact repeat, direct policy differences on common feasible states, price/population closures, the initial distribution plot, and prescribed-reference-price fit comparison. `calendar_entry_verification.json` checks the actual frozen calendar joint entry against C × income weights before any lifecycle solve. All 31 parameter rows must agree across arms; both final plot sets must equal their repeat hashes.

Stage only the compact overlay archive into the new remote root `/scratch/td2248/projects/entry_ratio_comparison_v1`. `launch_torch.sh` verifies the existing immutable 81-file remote stage at `/scratch/td2248/projects/grid_resolution_credit053_v2` and binds those files read-only alongside the new overlay. The original bundle/checkpoint is reused in place; no large source or checkpoint copying is required. Lead review and submission are separate.

```bash
sbatch /scratch/td2248/projects/entry_ratio_comparison_v1/launch_torch.sh
```

Local zero-solve invocation (fresh output path required):

```bash
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1/run_comparison.py preflight --out output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1/new_preflight
```
