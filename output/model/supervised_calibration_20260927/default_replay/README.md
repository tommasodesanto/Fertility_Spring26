# Optional-off native regression

One bounded Torch job, `18624270`, was submitted with one CPU, 16 GB RAM and
a 30-minute wall-clock limit. Its current numerical result is **pending
collection**. Subsequent SSH access failed when the shared connection closed;
submission itself returned the Slurm job ID successfully.

The purpose is exact reproduction of the September 25 selected **old-floor**
economy using the September 27 v2 source with the new preference and warm-price
options disabled. This is distinct from accepting or fitting the new no-floor
economy. Original payroll tax, target system, wealth/earnings inputs, normalizer
start and default step remain unchanged.

Driver: `code/cluster/validate_e5f_default_replay.py`. It validates the original
pair lock and the complete new-source inventory, runs a zero-solve preflight
comparing recomputed shared arrays with the retained checkpoint, then uses the
original complete objective evaluator in a fresh process. The retained point
previously required six stationary solves and 599.629 seconds of solve time.

Remote output:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/supervised_calibration_20260927/default_replay`.
Inspect `slurm-18624270.log`, `preflight/preflight.json`, and, on completion,
`replay/comparison.json`, `replay/target_comparison.csv`,
`replay/parameter_comparison.csv`, plus `replay/case/target_fit.csv` and
`replay/case/parameters.csv`. The case requests all 17 standard diagnostic plots.

Acceptance requires exact equality of every cell in all 13 inherited target
rows and all inherited parameter rows, all direct solution/shared/evaluation
array fields and the wealth grid/pre-choice population, price, loss, normalized
benefit, completed fertility, market residual and solve count. Gzip hashes are
authenticated for the reference but are not used as a numerical-equality test,
because gzip timestamps and additional metadata can differ.

The new review-only contract SHA256 is
`399abb6e9eab0d447dca627f94e3de4e6a8006d8920d02241fadac21f6c6ebae`.
The original pair lock SHA256 is
`6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1`.
Both immutable source bundles and original selected results remain untouched.
