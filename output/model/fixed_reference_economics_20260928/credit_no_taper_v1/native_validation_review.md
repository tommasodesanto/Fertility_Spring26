# Lead review of the prepared native validation harness

Do not submit this draft until these concrete issues are addressed. These are
validation-harness bugs found before any numerical launch, not new baseline bugs.

1. `register_overlay` imports the package before replacing parameters. Its
   `__init__.py` eagerly imports parameters AND solver, so the assertion that
   parameters is absent always fails. Bootstrap a proper package/module spec,
   register overlay parameters BEFORE executing package __init__, and verify
   the actual runtime solver's imported builder resolves to the overlay.
2. `renter_incidence` uses `np` without a local/global numpy import. With branch
   arrays [b,h,I,j,n,m,z], renter slice [:,0,:,j,:,:,:] has 5 dimensions. Floor
   broadcasts need grid[:,None,None,None,None], not six dimensions. Boolean
   indexing does not broadcast: sum negative current mass with mass[grid<0].
   Test the actual native-shaped synthetic incidence fixture on Torch.
3. Its `old_floor` is evaluated using changed P, so flag-on results falsely call
   the NEW floor old. Compute the original floor from authenticated reference P
   explicitly. Gate violations against the OPERATIVE floor only: flag-on debt
   may legitimately violate the old taper. Keep old-floor exposure descriptive.
   Verify owner-stayer branches have no renter mass and do not mix inactive
   branch policies into realised renter incidence.
4. `changed` requires exact equality with four expected names even though
   debt_caps and mean_labor_income_by_age will probably be unchanged. Require
   the flag and debt weights change, allow only the specified set, and certify
   any rebuilt bookkeeping field is identical if it existed. No economic
   change may hide in an allowed bookkeeping field. Confirm caps stay zero and
   owner_ltv_multipliers unchanged.
5. Launcher host preflight cannot resolve container-root pinned paths. Verify
   those inside the container, or translate paths explicitly. Bind cache path
   as writable inside container and keep original source mount read-only.
6. Total budget starts at launcher entry, before mock tests/authentication.
   Set/export started/deadline once there and preserve in controller launch.
   Reject resetting them. Case cap360, total1200, exactly2 lifecycle solves.
7. `mock_tests` recreates a toy loop rather than testing the real controller.
   Test the actual controller with tiny subprocess child fixtures for success,
   failed control, timeout, changed pin and duplicate output; no numerical
   solves. Never introduce an exposed production bypass for these tests.
8. Do not write misleading gates: cohort at unchanged prices is not clearing
   equilibrium. Keep all existing scientific/fiscal gates; report any failure.
   Add completion hashes of actual outputs, full tables and plot set, and
   authenticate reference/effective source identity again after both cases.

Prepare all mock tests and compiled-registration/incidence zero-solve tests on
Torch before numerical execution. New tested code is immutable before launch.
No numerical submission by worker. No new borrowing access, grid extension,
source mutation, relaxed gates or automatic retries.
