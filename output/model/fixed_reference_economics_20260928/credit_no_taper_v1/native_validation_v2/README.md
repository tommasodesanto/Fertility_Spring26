# Native validation v2 — prepared only

This immutable, two-case Torch harness validates the isolated renter taper
removal against **2007 stationary reference — block0506, September 28 verified
export**. It is not submitted or executed. The first fresh process is the
flag-off exact control; the second enables only
`renter_no_taper_estate_bound=True`, holds `lambda_d=0`, and rebuilds caps.
The original 160-node grid, prices/rents, psi, fiscal rules, owner rules and
zero-estate requirement remain fixed. No new unsecured borrowing, endogenous
prices, borrowing experiment, renormalization, retry or relaxed scientific gate
is present.

The plan pins the v2 driver (`93606aa1a23a0f6cc18d3e8c5b66374f70cdd1383805415b83ad66bf59d19765`),
the overlay (`9c6def300f76b2d5ac55c392e8a595fca881b1d78b15c1cba96016a30a3b83b9`),
the frozen native parameters/solver/kernels (pinned in the driver), the
reference manifest, and every helper. The launcher sets one start/deadline at
entry: 1 CPU, 24 GiB, 360 seconds per entire case and 1,200 seconds total.
The controller rejects any reset and re-authenticates all pins after both cases.

Before numerical cases, the same controller runs only zero-solve Torch tests:
tiny subprocess children cover success, failed control, timeout, changed pin and
duplicate output; the real package bootstrap verifies that the native solver's
builder resolves to the overlay; and a real seven-axis branch fixture exercises
realised renter incidence. These tests gate the numerical controller and are
not claimed as run here. Successful cases require the complete 14-row fit table,
31-row parameter table, all 17 standard plots, output hashes, fixed-rule fiscal
readout, all parent gates, and operative (not descriptive old-taper) renter
floor compliance. This remains fixed-price conditional validation, not market
clearing, a demographic endpoint, transition, recalibration, or adoption.

Single immutable staging/verification recipe for the lead, after reviewing the
four v2 files: copy `run_native_validation.py`, `launch_native_validation.sh`,
`plan.json`, `README.md`, and the already authenticated `overlay_v1/` directory
into a new read-only
`/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929/source_native_validation_v2/`.
Do not alter any staged byte. Then submit exactly:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v2/launch_native_validation.sh
```

Inside the container, `verify_plan` checks every container-visible source pin;
the original project mount is read-only and the cache/result mounts are writable.
