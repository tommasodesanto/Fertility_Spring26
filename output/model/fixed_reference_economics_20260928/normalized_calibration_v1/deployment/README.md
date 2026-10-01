# Normalized floor calibration deployment

The user authorized a 24-chain experimental calibration with the economics of the preceding free-\(\psi\) floor search,
the 120-by-9 wealth and income grids, the existing
14-row target contract and original base weights, and free \(\psi_{child}\).
The housing supply coefficient \(H_0\) is derived internally to clear housing
at normalized population \(N_0=1\); it is checked against its [.2, 80] bound
only at an accepted root. No price anchor or grid change is used. The run is
not an adopted calibration.

## Submission and resource contract

- Native incumbent gate: Slurm job `18979627_0`.
- Search array: `18979628_[0-23]`, submitted with `afterok:18979627`.
- Each task has 1 CPU, 24 GiB, one numerical thread, a 7,200-second budget
  from actual launcher start, at most 100 objective calls, and a 900-second
  reserve for selected native verification. Original base weights only; no
  automatic retries or budget extension.
- Stage archive SHA-256: `2f4af6da9c289d677f9437d9b5986105bbde269ff550a04d6411baf0cb054529`.
  It contains 229 source files, including all 223 final source pins. The
  submission receipt is `submission_receipt.json`.

## Native incumbent gate

Job `18979627_0` completed with exit code zero in 2:16. The full native
incumbent root and selected repeat each ran once. All 14 target rows and the 30 non-`H0` estimates in the 31-row
parameter report match the verified center exactly. The internally derived
`H0` changes to express the benchmark at population one. ROOT and REPEAT
match each other exactly for all 14 target rows and all 31 parameter rows. The repeat
receipt reports `exact_full_ge_repeat_passed` and includes hashes for all 17
standard diagnostic plots. The independently reviewed root closure has
`H0=6.838655072054291`, physical housing supply and demand both
`5.988286068783477`, zero housing residual, and renewal residual
`6.1768e-7`.

The full compact gate receipt and ROOT/REPEAT tables, parameters and closure
files are in `monitor_snapshot/first_array_startup/`. Startup-only receipts
and normalization source review are in `monitor_snapshot/smoke_startup/`.

## First array startup snapshot

The recorded matrix shows all 24 exact initializers passed, all 24 search
processes entered their first full-GE objective case, and no startup failure
receipts. Every chain records the same target fingerprint
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1` and weight
fingerprint `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`,
with `N0=1` and derived `H0`. Chain 0 had completed its first objective case
and entered the next case at the snapshot; the other 23 were in case `0000_nm`.
See `monitor_snapshot/first_array_startup/chain_startup_matrix.json` and the
compact scheduler snapshots beside it. Monitoring stopped after these
startup and first-case checks; the array continues under its submitted budget.

The four previous continuous-v2 floor calibration tasks (`18973391_[0-3]`)
were collected and cancelled before this launch. Transition job `18974228`
was identified separately and left running.
