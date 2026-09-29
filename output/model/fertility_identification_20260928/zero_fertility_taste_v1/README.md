# E04/E05 exact-zero fertility taste preparation

This isolated packet prepares two separate experiments against the original
one-birth E01 selected case `one_birth_024_gn1_0` (loss
`7.826226594410982`, first-birth fixed utility cost
`0.35270914196085973`, normalized `psi_child=0.12184551693359474`). It
does not alter the shared model or promote a candidate. **No model run has been
launched.**

`first_scale_zero` fixes `kappa_fert=0`; `later_scale_zero` fixes
`kappa_fert_continuation=0`. The other taste scale stays at E01 for the
initial comparison and is free in the subsequent local refit. Each lane also
frees the first-birth fixed utility cost within `[0,8]` during refitting. The
zero restriction is external to the native `[0.02,50]` search bound and must
be reported as such. Conception risk, the one-birth opportunity, birth timing,
earnings, entry distributions, transfers and floors, housing and credit
primitives remain unchanged. The same ten scored targets, three validation
rows, completed-fertility normalization to `2.1`, and demographic renewal
gate apply. Successful evaluations retain 14 fit rows, 31 parameter rows,
and all 17 standard plots.

`patch_solver.py` changes only the active sequential first/later fertility
choice branches in an authenticated source copy. At scale zero, the value is
the maximum of the unchanged wait and conception-risk-adjusted attempt values.
The probability is one for the maximizing action, zero for the other, and
zero for both in dead states. Exact ties choose wait by numerical convention;
the positive-scale softmax limit splits ties. Positive-scale expressions keep
their original arithmetic order. The overlay is grafted into loaded canonical
and byte-identical private solver modules. The isolated binder validates the
nine free native bounds, uses `0.02` only as an internal construction surrogate,
then overwrites the actual `P` scale to zero before evaluation; first-birth
`eps_fert` is set to zero with `kappa_fert`. Native receipts and saved parameter
rows must contain the actual zero, not the surrogate.

The two serial controllers each allow at most 16 cases: anchor replay,
zero-scale center, nine one-sided probes, three damped Gauss–Newton proposals,
and two exact final repeats. Each lane requests 1 CPU and 24 GiB on Torch,
with a 4-hour controller budget, at most 8 stationary solves and 35 minutes
per case, and 70 minutes reserved for repeats. The controller writes heartbeat,
latest completed case and best-so-far files. Missing derivatives are never
imputed. A Jacobian is a local numerical diagnostic; exact-zero choices can
make the map nonsmooth, and full rank alone cannot certify statistical
identification or global feasibility.

Required Torch sequence, after lead review:

1. Run `run.sh tests` under Slurm. This runs the two-lane synthetic controller
   and exact-zero menu tests with no model solves and writes `TESTS.json`.
2. Run `run.sh prepare` under Slurm. This authenticates the E01 anchor,
   source, objective and tests, and creates a fresh expiring `config.json`.
3. Pin its SHA256 in `EXPECTED_ZERO_TASTE_CONFIG_SHA256`; submit each lane as
   `run.sh search first_scale_zero` or `run.sh search later_scale_zero` under
   Slurm. Search requires a fresh unexpired config and exclusive output path.

There is no automatic launch, adoption, transition or fallback. Full model
smokes, positive-scale replay against E01, KFE/cache checks and source review
remain required before either search submission.
