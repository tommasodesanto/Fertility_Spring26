# Five-ratio entrant diagnostic — Torch 18888956 running

The lead submitted **18888956** after independent zero-solve smoke, formula and
calendar-entry review, five-file archive authentication, and remote SHA check.
The exact submitted archive is retained as `deployed_stage.tar.gz`, SHA
`4a0381f6d0f66d52f46c05124cd61442bc499ca734cf847fb2a86bccf57c3853`.
`launch.json` records the run. No result is certified yet.

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
