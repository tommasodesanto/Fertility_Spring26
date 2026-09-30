# Native validation packet — prepared, not submitted

This is the bounded compiled validation for the isolated renter-taper removal
from **2007 stationary reference — block0506, September 28 verified export**.
It runs exactly two fresh-process lifecycle solves: first a flag-off exact
control, then `renter_no_taper_estate_bound=True` with `lambda_d=0` and rebuilt
debt caps on an authenticated deep copy of the reference parameters.  The asset
grid remains the original 160 nodes.  Positive unsecured credit, natural credit,
birth renormalization, preferences, income, entry, fiscal rules, prices/rents,
targets and timing are not changed.

The driver registers the existing hash-pinned overlay before it imports any
runtime package.  It authenticates the checkpoint/source/target/parameter pins,
requires the control to reproduce all 113 arrays, 14 fit rows, 31 parameter rows
and the standard 17 diagnostics, and then retains the parent household, mass,
budget, probability, occupied monotonicity, transaction, estate and fiscal
gates.  A renter-by-age table uses realised current renter masses and actual
branch policies; it reports negative current assets, negative policy saving,
old/new lower-bound binding and violations without assuming every `loc_probs`
entry is an active choice.

`launch_native_validation.sh` is deliberately **not submitted**. Before the two
numerical cases it runs zero-solve controller mock tests for success, failed
control, timeout, duplicate output, and changed pin.  It requests one CPU,
24 GiB and 20 minutes; the controller enforces 360 seconds per case and 1,200
seconds total, with no retry. Results are written only on Torch under
`/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929`.

Before lead submission, stage this four-file packet plus the already reviewed
`overlay_v1/parameters.py` into the immutable remote source directory
`/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929/source_native_validation_v1`.
Verify its hashes against `plan.json` and retain no mutable files below that
source directory. Exact prepared command after that review:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v1/launch_native_validation.sh
```

This is conditional fixed-price validation, not a cleared equilibrium,
demographic stationary endpoint, transition, recalibration, or adoption. Any
maintained fiscal-gate failure is written and stops the run; the pension is never
altered to clear it.
