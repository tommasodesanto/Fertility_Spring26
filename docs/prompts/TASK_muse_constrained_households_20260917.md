# Task for Muse (OpenCode): who is constrained, and does scarcer space change it?

Repository root is the current directory. Read `CLAUDE.md` and
`code/model/sandbox/README.md` first. Python: `code/model/.venv/bin/python`.
Permitted: read anything; write only under `code/model/sandbox/` (one new
script `diagnostics_constrained.py` and new spec files under `specs/`) and
`output/model/sandbox/constrained/`. Not permitted: editing the package
`code/model/intergen_eqscale_seq_optimized/`, `code/model/tools/`, targets,
git, cluster jobs. Time limit: 3 hours. Never use the word "parity" in prose
or comments you write.

## Part A. Who is constrained, from the two saved steady states

Inputs: `output/model/sandbox/baseline_code/baseline_psi_fixed/` (down payment
share 1-phi = 0.20) and `output/model/sandbox/frictionless/frictionless_nodp_psi_fixed/`
(phi = 1, no down payment), both solved at the same parameters and the same
psi. If the saved folders do not contain the solution objects (they hold only
summary.md, moments.csv, parameters.csv, graphs.pdf), re-solve both specs with
`run_ss.py` in-process (about 5 minutes each) and keep the solution objects in
memory; see how `diagnostics_dependents_by_age.py` reuses `run_ss.py`'s solve.

Produce one table by model age (18, 22, ..., 62) and by children at home
(0, 1, 2, 3+), for the baseline: (1) share of households renting; (2) share of
renters exactly at the rental cap (6 rooms); (3) share of owners in the
owner-only sizes (8 rooms and above); (4) share of households whose housing
choice differs between the two solutions at the same state (this is the
share for whom the down payment binds); (5) median liquid wealth of renters
relative to the down payment on the smallest family-sized owner unit (8 rooms
at the solved price); (6) the attempt probability for a first birth (n = 0)
and the share of first births coming from renters at the cap. Also report the
same rows for the no-down-payment solution where they change.

## Part B. Does scarcer space change the answer?

Mean occupied rooms is 5.72 in the sandbox baseline against a data target of
5.56 (and 6.42 in the retained calibration against 5.56). Make space scarce at
fixed parameters by lowering the supply scale `H0` (a spec override). Find, by
two or three trial solves at fixed psi (`psi_mode: fixed`, spec files
`scarce_h0_*.yaml`), the value of `H0` at which mean occupied rooms is within
0.05 of 5.56; then solve two more steady states at that `H0`: down payment kept
(phi default) and no down payment (`phi: 1.0`). Report for both, side by side
with the two originals: completed fertility, childless share, mean first-birth
age, first births at 30+, ownership 30–55, mean rooms, first-birth rooms
response, the three-plus versus one-to-two rooms gap, price, and the Part A
row for "share whose housing choice changes when the down payment is removed"
at ages 26–38 with children at home 0 and 1.

## Deliverable

`output/model/sandbox/constrained/README.md` with both tables and a
ten-line reading: does the down payment bind for family-forming households at
the baseline, and does it start to once space is scarce? Plus the CSVs behind
the tables and the spec files used. Final message: the receipt in the format
Outcome / Verification / Artifacts / Unresolved / Reported cost. If a solve
fails (the package raises an infeasibility error), record the exact error and
the spec, and continue with the remaining solves.
