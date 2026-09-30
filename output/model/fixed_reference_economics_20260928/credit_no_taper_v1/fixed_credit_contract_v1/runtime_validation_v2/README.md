# Fixed-credit runtime validation v2

Reference: **2007 stationary reference — block0506, September 28 verified export**.

## Reviewed launch, September 29 evening

**Terminal submission failure:** job 18845938 failed in two seconds, before
Python or checkpoint authentication. Slurm copied the submitted wrapper to its
spool directory, so its `$0`-relative sibling-launcher lookup failed. Zero
lifecycle evaluations and zero numerical budget were consumed. Preserve this
source version and receipt; a narrowly corrected launcher will use a new
immutable `runtime_validation_v3/` packet. Do not resubmit v2.

Torch zero-lifecycle smoke **18845938** was submitted once. Do not duplicate it.
Immutable remote source:
`/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v2/source`.
Results and scheduler logs are under the sibling `results/` directory.
Driver SHA256: `f9940add31491a6132cff52867d469f91b6134b1e5572b42310aea5a80caf4c1`.
The smoke tests original checkpoint authentication, actual overlay module
origins, compiled boundary fixtures and the strict-zero affordability proof.
It permits zero lifecycle evaluations and has a five-minute timeout.

The authorized successor is **one** exact baseline control, only after the lead
reviews a passing smoke receipt and creates `results/SMOKE_LEAD_REVIEWED`.
Submit `source/launch_control.sh` once and record its job ID here. Its fixed clocks
begin at launcher entry: 300 seconds for the complete control case, 900 seconds
total, one thread and 24 GiB. No retry, strict-zero lifecycle evaluation,
recalibration, positive-credit choice, price search or GE attempt is allowed.
Read `launch.json` for immutable start/deadline timestamps. All existing reference
17 diagnostic plots remain retained; this replay is not a new model solution.
The hourly heartbeat follows this launch and may dispatch that reviewed successor.

This packet authenticates the immutable block0506 reference with the original fixed-price authenticator before importing the approved three-file fixed-credit overlay as a separate namespace package. Dependencies not in the overlay resolve from the frozen original package path. It never mutates the canonical package or replaces the original ancestry guard.

`smoke` performs zero lifecycle solves. It authenticates the checkpoint, verifies all overlay hashes and module origins, executes the compiled renter-saving and owner-sale fixtures, and derives the two strict-\(D=0\) entrant cash failures from the actual reference parameters and overlay functions. The documented negative cash cells make a strict lifecycle solve mathematically infeasible under unchanged entry, so this packet does not run one.

`control` is allowed only after a successful smoke receipt has been reviewed. It leaves the scalar absent/`None`, deep-copies the authenticated parameter object without rebuilding it, performs one lifecycle solve at the retained reference price, and requires exact equality for every saved solution array. It also runs the established gates and verifies the complete 14-row fit and 31-parameter table. Any mismatch stops with shapes and maximum absolute differences in the failure text.

On Torch, stage once with `bash stage_source.sh`. Then submit the smoke with `sbatch launch_smoke.sh`. After lead review of `results/smoke/receipt.json`, create `/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v2/results/SMOKE_LEAD_REVIEWED` and submit the control with `sbatch launch_control.sh`. No command in this packet submits a job automatically.
