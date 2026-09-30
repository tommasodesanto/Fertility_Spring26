# Stationary refactor: result and limits

The isolated test package is [code/model/refactor_lab](../../../code/model/refactor_lab/README.md).
It uses **2007 stationary reference — block0506, September 28 verified export**,
including its complete saved parameters, wealth grid and entrant distribution.
There is no recalibration, change in targets or normalization of fertility.
The active model, manuscripts and September 14 reference are preserved.

## What changed

The normal calculation has one entry point and separates shared state
construction, household decisions, distribution evolution and equilibrium.
The engine has 18 Python files (including its initializer), 8,872 lines;
the old directory has 53 files and 22,871 lines, including calibration and
other paths outside this stationary calculation. This is a scope comparison,
not replacement of every old feature. The former 8,722-line solver is split
into four modules; 99 moved definitions retain their source and normalized
bytecode. Verification tools are in one separate folder and tests in one suite.

The package contains the separately reviewed borrowing correction described
below. In reference mode it reproduces the old borrowing behavior for an exact
comparison. The only performance transformation replaces repeated sorting and
interpolation searches inside saving maximization with indexed interval access.
It preserves all candidates, formulas, order and strict-improvement tie-breaking,
without assuming that the continuation value is globally concave. The 18 final
engine files are byte-identical to the tested indexed source stage.

Static renter allocation was already analytical. Saving already evaluates the
analytical interior first-order-condition candidates on every piecewise-linear
continuation interval, alongside its endpoints. The one-market pension also
already uses an analytical age/income recursion. Missing elementary derivatives
were not the source of the improvement found here.

## Fixed-population baseline timing

A matched local pair used one core on an Apple M5 Pro, Python 3.13.15,
NumPy 2.2.6 and Numba 0.61.2. Both used fresh compilation caches, identical
profiling, all the same parameters and tolerances, and an initial price 5%
above the saved equilibrium price. Both performed four household solves,
five forward-distribution passes and one final output upgrade that reused
household policies.

| Measured stage | Original | Refactored |
|---|---:|---:|
| Full stationary GE, including compilation and identical profiling | 152.65 s | 99.55 s |
| Whole benchmark process, including input/output | 156.00 s | 103.10 s |
| Household work, included within GE time | 122.20 s | 68.11 s |
| Distribution work, included within GE time | 30.42 s | 31.40 s |
| Sampled peak process-tree memory | 1.865 GiB | 1.909 GiB |

Observed GE time fell **34.79%**, a **1.53× speedup**. This is one matched
cold pair at the named reference and start, not a guarantee for every
calibration or initial price. It excludes the historical calendar observers,
full tables and 17-plot certificate. All 90 native arrays, their complete key
sets and effective public parameters matched exactly; market and pension
checks passed. [Timing and comparison receipts](native_local_pair_v1/summary.json).

The earlier 886.527-second credit GE had a different credit contract and six
lifecycle calls, including an exact repeat. It is not a matched before/after
benchmark. The improvement here is not described as reducing that run from
15 minutes to 100 seconds. A household inner-loop time is never substituted
for full GE time.

## Corrected-credit scalar/indexed replication

A separate Torch comparison (job **18876666**) held the corrected unsecured
credit limit (D=0.14), birth renewal, and population-scaled housing closure
fixed while comparing the original scalar saving search with indexed interval
access. Both arms completed six lifecycle evaluations. Workflow time was
549 s for scalar and 417 s for indexed, a **24.04% reduction**; measured solve
time was 416.4315 s versus 283.3924 s, a **31.9474% reduction (31.95%)**. Remaining
workflow time was 132.5685 s versus 133.6076 s. The exact comparison passed
eight closure paths, 87 solution-array paths in each of two saved NPZ files,
all 14 fit rows, all 31 parameter rows, and all 17 final plot files. Indexed
outputs also match the earlier corrected-credit run, job 18869900. The full
receipt and tables are in the [replication packet](small_credit_replication_v1/README.md).

This matched pair isolates the saving-search implementation under the same
corrected-credit contract. For context, the earlier 886.5269 s credit GE used
513.9394 s in solves and 372.5875 s elsewhere on a 262-point wealth grid;
the current pair uses 160 points. Job 18869900 used 467 s total, 314.199 s in
solves and 152.801 s elsewhere, also on 160 points. Those earlier runs differ
in credit contract and run context, so they do not measure a saving-code speedup.
The residual workflow time includes unprofiled setup and reporting overhead,
not all report generation. The older 15-minute result is not a matched
comparison with this scalar/indexed pair.

## Numerical verification

- **Independent fixed-price replay:** scalar and indexed packages, two fresh
  repetitions each. Every repetition passed 113 exact, finite nested array
  paths, all 14 target rows, all 31 parameter rows, and identical hashes for
  the 17 standard plots. These include shared arrays computed by the lab.
  [Collected evidence](fp_pair_gridfix_18850069/COLLECTION.md).
- **Final package:** 27 component tests pass, including borrowing boundaries,
  typed input identity, provenance and the maintained compiled saving kernel
  against the original on saved model continuation columns. Driver smoke
  checks cover sequencing, fail-fast, timeouts and solve budgets.
- **Local full GE:** the strict symmetric 90-array comparison and complete
  effective parameter comparison pass. Household solve counts and accepted
  prices agree. [Strict array comparison](native_local_pair_v1/comparison.json).
- **Full Torch GE/reporting pair:** job 18851943 passed. The strict comparison
  matched all 87 saved solution/shared arrays with identical complete key sets,
  all 14 fit rows and 31 parameter rows were byte-identical, and the actual
  original/refactored files for all 17 plots had identical hashes. Both engines
  passed the unchanged budget, purchase, estate, fiscal, policy and distribution
  checks, plus the GE market check. [Collected certificate](ge_pair_18851943/COLLECTION.md).
  The x86 Torch solve stages took 275.95 s original and 204.87 s refactored
  (25.8% less); the two-engine validation batch took 631 s including reporting.

The complete unchanged-reference [target-fit table](fp_pair_gridfix_18850069/indexed/oracle/rep1/target_fit.csv)
contains every target, model value, gap, weight and loss contribution. The
[parameter table](fp_pair_gridfix_18850069/indexed/oracle/rep1/parameters.csv)
contains all estimates, bounds/restrictions and bound flags. These are inherited
reference values, not new estimates. The [17 standard plots](fp_pair_gridfix_18850069/indexed/oracle/rep1/standard_diagnostics)
are retained; representative market, policy and fertility plots were inspected.
No ad hoc replacement graph set was introduced. The fresh GE
[target-fit table](ge_pair_18851943/verify/lab_certificate/certificate/target_fit.csv),
[parameter table](ge_pair_18851943/verify/lab_certificate/certificate/parameters.csv),
and [standard plots](ge_pair_18851943/verify/lab_certificate/certificate/standard_diagnostics)
are also retained. Plot equality in the GE comparison is old versus new at
their common solved price; it does not require the slightly different frozen
checkpoint price to generate identical pictures.

The GE renewal gap is **1.7028122e-6**, against **7.9199670e-7** in the frozen
reference. Their difference is below the lead-selected diagnostic threshold
of 1e-6. This is **not** a claim that the absolute renewal gap passes a
1e-6 gate. The reference manifest supplies no absolute renewal tolerance;
no preference or demographic closure was changed to improve that number.

## Borrowing correction and required author decision

The reviewed contract makes the unsecured limit an explicit finite scalar
\(D\geq0\). Renters require \(b'\geq-D\), with the existing zero floor at
terminal or positive-mortality ages. An owner may sell into renting only if
raw net sale wealth \(b+(1-\psi)pH\geq0\), before interpolation or clipping.
Buyer and incumbent-owner rules are preserved. The lab implements this rule;
reference mode remains available solely for controlled reproduction.

At the selected \(D=0\), two retained age-18 entrant cells have no feasible
choice. They represent about **0.00802% of entrant mass**. At wealth
\(-0.2558139535\), their cash before rent, consumption and saving is already
**−0.1340497885** and **−0.0677697562**. Buying resources are negative too;
positive housing prices cannot fix those cells under the unchanged contract.
The other chat's authenticated compiled test independently confirms this.

The code reports the cells and stops. No asset truncation, deleted mass, debt
forgiveness, new transfer, chosen positive credit or recalibration was used.
A corrected-credit GE therefore awaits the author's decision about inherited
entrant debt. The timings above are reference-mode timings, not a solved
zero-credit economy.

## Review and publication boundary

Claude Opus 5.5 implemented the package and addressed review findings. Sol
reviewed source, validation coverage, the actual receipts and final wiring;
the lead reviewed the borrowing mathematics and saving transformation and
verified the tested/final source identity. Luna handled bounded discovery,
receipt collection and mechanical fixes. See [evidence review](sol_result_evidence_review.md)
and [package review](sol_final_package_review.md).

This is a reviewed stationary test package, not a standalone publication
release of the full project. Transition code, dated observers and full
reporting remain outside this extraction. Local historical certification is
intentionally blocked by a changed pinned helper in the working tree; the
normal model runs locally, and full historical certificates use the intact
Torch snapshot. No source-authentication check was weakened.

Individual one-core runs, including overnight, are now explicitly permitted
in both project instruction files. Torch remains the route for long batches
and parallel work. All numerical checks and final source review are complete. The launcher also
passes a relocated-Slurm-wrapper check, including rejection of a missing source
path. Final source files, compact evidence and this report are backed up together;
large checkpoints, arrays, caches and environments stay outside the commit.

Codex account usage was 1% at entry and 7% at the last closeout check, an
account-wide increase of 6 percentage points shared with other chats. This
is below the requested 20-point alert threshold. Claude Max has a separate
allowance not exposed by the Codex usage tool; no combined usage claim is made.

## Verified production grid inventory

[Exact grid inventory](production_grid_inventory.txt) records every wealth node and all household state axes from the pinned September 28 block0506 manifest. The dense value array has shape `(160, 6, 1, 17, 15, 4, 4)`. The active child-state mode is `independent_count`; the shared-clock fallback does not apply. Both later one-birth and two-birth candidate parameter tables also record 160 wealth nodes. This is a read-only extraction, not a new numerical accuracy test.
