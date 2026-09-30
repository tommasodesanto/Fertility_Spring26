# Native validation v3 — prepared only

This immutable two-case Torch harness tests the isolated renter-taper removal
against **2007 stationary reference — block0506, September 28 verified export**.
It is not submitted or executed.  The first fresh process is the exact flag-off
control; the second enables only `renter_no_taper_estate_bound=True`, retains
`lambda_d=0`, and rebuilds the zero debt caps.  The original 160-node grid,
fixed prices/rents, fixed psi and fiscal inputs, owner rules, and zero-estate
gate remain unchanged.  This is fixed-price conditional validation, not a
cleared equilibrium, demographic endpoint, transition, calibration, or new
credit design.

The driver contains the only runtime exception: before parent authentication it
installs an import finder for precisely
`intergen_eqscale_seq_optimized.parameters` at the plan-pinned overlay SHA and
an in-memory extension of the native identity guard.  It never changes an
original source file or any module `__file__`.  The extension requires the
original current-source solver path and SHA, rejects every other external
`intergen_eqscale_seq_optimized*` module, and falls back to the original guard
when no overlay is present.  After authentication it verifies that the actual
model resolved both `setup_parameters` and `unsecured_debt_floor` from the
registered overlay.  Each case receipt discloses this sole authenticated-module
exception; all existing source/contract/runtime authentication and numerical
gates still run.

The plan pins the v3 driver
(`1493e2a591118671b7fffd900e3955ea8c49d1a457fa8fccb98e5da98881708d`),
overlay (`9c6def300f76b2d5ac55c392e8a595fca881b1d78b15c1cba96016a30a3b83b9`),
native solver/parameters/kernels in the driver, current identity-runtime hash,
reference manifest, and each staged helper.  The launcher preserves the
entry-time 1,200-second total budget, two 360-second fresh cases, 1 CPU and
24 GiB.  Successful cases still require all 14 fit rows, 31 parameter rows,
17 standard plots, fixed-rule fiscal readout, original gates, and operative
renter lower-bound compliance.

Run a Torch-only zero-solve preflight first (it uses imports, actual controller
fixtures, shaped arrays, and tiny child processes only; it reads no checkpoint,
lifecycle, KFE, or GE):

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v3/launch_native_validation.sh --self-test-only
```

The preflight tests include a rejected external package module, owner-stayer
mass consistency, a permitted zero terminal renter saving, a rejected negative
terminal renter saving, and actual controller timeout termination.  Production
is refused until the caller pins the SHA256 of that successful receipt and the
launcher confirms the receipt's plan, driver, and all scientific plan fields
match exactly:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v3/launch_native_validation.sh --preflight-receipt-sha RECEIPT_SHA256
```

Immutable staging recipe: copy exactly `run_native_validation.py`,
`launch_native_validation.sh`, `plan.json`, `README.md`, and the already
authenticated `overlay_v1/` directory into the new read-only
`/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929/source_native_validation_v3/`.
Verify the plan's `pinned_files` and driver SHA against those staged bytes
before submitting.  Do not alter any staged byte.  The launcher refuses
duplicate preflight or production output directories.
