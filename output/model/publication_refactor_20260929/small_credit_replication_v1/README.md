# Matched corrected-credit scalar/indexed replication

Completed matched replication. Job **18876666** compared the sole saving
optimization under the **same corrected D=.14 contract and full birth-renewal,
population-scaled housing closure**. It is not a comparison with
pre-correction economics or the normalized-population lab GE.
No earnings, entry, grids, preferences, mortgage, fiscal, supply, measurement,
or gate changes. Both arms copy the completed experiment's driver and sources.
Scalar substitutes certified `scalar_src_gridfix` kernels SHA19dceb70...;
indexed retains SHAfe7d43af.... All other files are identical. AST verification
finds only `exhaustive_saving_scalar` differs, plus four indexed helper additions.

`prepare.py` uses only stdlib, creates immutable arm trees, verifies this
one-file difference and emits `preparation.json`/`source.sha256`. Preparation
ran successfully. `launch_torch.sh` passed shell syntax verification. These
statements describe preparation before the launch records below.

New independent budget: 2400 seconds from launcher entry; one CPU24GiB; all
numerical thread limits1. Two zero-lifecycle authenticated loop preflights
within300seconds, then scalar and indexed sequentially on the same node,
separate fresh production Numba caches, at most6 lifecycle evaluations each,
300seconds per evaluation, with existing root/repeat reserves unchanged.
Six per arm matches the completed price path including exact root repeat.
Insufficient convergence or time stops; no automatic extension/retry.
Observed indexed job was467seconds total; scalar overhead unknown. Rough
paired estimate15–20minutes, remaining margin within40minutes. No grid sweep.

Stage this whole directory read-only at
`/scratch/td2248/projects/small_credit_replication_v2/source`; create sibling
logs/results/cache directories. Refuse overwrite. Authenticate `source.sha256`
before submission. The bundle and frozen project mount are the same as18869900.
Submit `sbatch source/launch_torch.sh` only after lead review. Production output
is kept in remote results/{scalar,indexed}; never overwrite existing small-credit
results. The copied driver checks pinned bundle/checkpoint/frozen observer at runtime.

`compare_pair.py` compares every price's closure JSON; exact finite saved
solution/shared arrays; both selected/repeated full14-fit/31-parameter tables;
and actual17PNG hashes in both final reports. It fails on missing paths/arrays
or unequal results. Outputs one compact comparison receipt. Original driver
already checks each arm's exact root repeat. Historical18869900 cross-comparison
should use its remote arrays plus compact collected tables/closures/PNG evidence;
those arrays are not locally collected. The comparer additionally requires exact historical-indexed identity against
`/scratch/td2248/projects/small_credit_v1/results/full` using the same comparisons.
It records per-arm wall time, six-case solve sum, residual overhead, and speedup;
it makes no historical hardware timing claim. The lead owns collection.

Important scope: scalar is the pre-optimization implementation in the extracted
engine with the accepted credit correction. There is no independent replay of
an unmodified legacy pre-credit engine; such a replay would change economics.

## Completed replication

The first submission, job **18876610**, stopped after 31 seconds because its
smoke observer omitted the `phase_b_ge/` output directory. It completed zero
preflights and zero lifecycle solves. The narrow retry repair creates that
mock output directory and separates the external preflight timeout from the
driver's existing 700-second reserve; numerical model sources were unchanged.

Retry job **18876666** completed in 17m32s under a new 40-minute, one-CPU,
24-GiB budget. It passed both authenticated zero-solve preflights, completed
12 lifecycle evaluations and passed the matched comparison. The staged remote
`source.sha256` manifest was checked and remained unchanged. Full receipts are
in [`collected/`](collected/) and [`launch.json`](launch.json).

| Measure | Scalar saving | Indexed saving |
|---|---:|---:|
| Workflow elapsed | 549 s | 417 s |
| Sum of lifecycle solve time | 416.4315 s | 283.3924 s |
| Other workflow time | 132.5685 s | 133.6076 s |

The indexed version reduced total workflow time by **24.04%** and measured
lifecycle solve time by **31.95%**. The comparison passed all eight closure
paths; two saved solution-array files each had 87 paths; both arms matched on
the 14-row fit, 31-row parameter table and 17 final PNGs. The historical
indexed output from job 18869900 also matched exactly in the same comparisons.
The complete fit and parameter tables are linked at
[`indexed/phase_b_ge/selected_repeat_final/target_fit.csv`](collected/indexed/phase_b_ge/selected_repeat_final/target_fit.csv)
and
[`indexed/phase_b_ge/selected_repeat_final/parameters.csv`](collected/indexed/phase_b_ge/selected_repeat_final/parameters.csv).
The same 17 diagnostic plots are already retained in
[`small_credit_v1 collected output`](../../fixed_reference_economics_20260928/credit_no_taper_v1/small_credit_v1/collected_v1/full/phase_b_ge/selected_repeat_final/standard_diagnostics/);
no duplicate plot set was added.

For context, the earlier 886.5269-second credit GE used 513.9394 seconds in
solves and 372.5875 seconds elsewhere; it used a 262-point wealth grid, versus
160 points here. Job 18869900 took 467 seconds total, 314.199 seconds in solves
and 152.801 seconds elsewhere, also on the 160-point grid. These are not
matched timing comparisons: the earlier GE had a different credit contract,
and the 18869900 workflow differs in driver/run context. The residual workflow
time includes unprofiled setup and reporting overhead, not all report
generation. The new scalar/indexed comparison holds D=.14 and the full closure
fixed and isolates the saving implementation.

## Source lineage and scope

The refactor's authenticated source and checkpoint lineage is the frozen
**2007 stationary reference — block0506, September 28 verified export**,
checkpoint `repeat_0212` (SHA prefix `b15ba92d`). See the
[calibration reference record](../../fertility_identification_20260928/README.md)
and the [refactor source and verification record](../REPORT.md). The later
overnight `one_birth_024_gn1_0` E01 candidate and experimental
`two_birth_024_gn1_0` E02 candidate are distinct, unadopted calibration results
([overnight readout](../../fertility_identification_20260928/two_stream_overnight_v1/morning_readout_v1/RESULTS.md));
they are not the refactor input. Equivalence evidence therefore applies to the
pinned block0506 reference and tested stationary calculations, not all
calibration candidates or parameter regimes.

The refactor verification does not test the original calibration's outer
normalization/search loop, which adjusts child-benefit preference to completed
fertility 2.1. It also does not test dated transition dynamics. The separate
[four-shock transition record](../../fixed_reference_transition_20260928/four_shock_v1/README.md)
uses block0506 through the frozen legacy source tree; its reported shock fits
failed before estimating shocks.

## Narrow preflight repair

Job 18876610 failed in 31 seconds before lifecycle evaluation because the mock
observer omitted creation of `phase_b_ge/`, while the real observer creates it.
The smoke-only wrapper now creates that directory after driver output creation.
Production numerical files are unchanged. The external preflight timeout is
300 seconds; the driver receives the global 2400-second deadline so the
unchanged mock price loop can check its 700-second selected/repeat reserve.
The mock selected price factor is 1.02; qref is reused, followed by
lower_085, upper_115, root_01 and selected_repeat mock solve calls. The normal
six lifecycle evaluations per arm and 12 total cap remain. The failed v1
receipt and successful retry under the immutable v2 remote root are retained.

The staged `source.sha256` authenticates the immutable Torch source snapshot;
this local README and `launch.json` record post-launch outcomes, and the staged
manifest was intentionally left unchanged. The 138 collected-result files were
separately hash-verified in `collected/collection_verification.json`.
