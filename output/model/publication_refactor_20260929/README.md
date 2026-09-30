# Publication refactor laboratory

**Final status:** the stationary test package and reference-mode verification
are complete. The zero-credit entrant decision remains open. See the
[final report](REPORT.md) and [verification record](final_verification.json).

Author request, September 29, 2026: simplify the executed model, investigate the
cost of a complete equilibrium, and incorporate the borrowing correction in an
isolated test folder. No recalibration or other economic changes. This packet
records the specification, review and verification; it is not a new baseline.

## Fixed reference and permitted difference

Use **2007 stationary reference — block0506, September 28 verified export**.
The authoritative identity is
`../fertility_identification_20260928/fixed_reference_manifest.json`.
Its `checkpoint` names repeat_0212_primary, SHA256
`b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
The original, seed and other repeat checkpoint hashes have different roles.
Load the actual serialized parameters, grids and entrant distribution, including
`psi_child=0.1355551166583114`. Do not use constructor defaults or a newer fit.

The other chat, **Analyze fixed calibration mechanisms**, owns the borrowing
correction, whose small compiled checks passed on September 29. Its author instruction is an
explicit nonnegative unsecured-credit limit `d_bar`, with renter saving
`b_next >= -d_bar`, and mortgage repayment from financial assets and net sale
proceeds before becoming a renter. Start the separately labeled correction at
zero. Preserve the adopted incumbent-owner and mortality rules; any ambiguity
requires lead/author resolution, not an implementation guess. Positive credit
has not been selected. The old age taper and inherited renter principal rollover
must not reappear as hidden defaults in the corrected mode.

Known economic blocker: at zero credit, two retained age-18 entrant states are
infeasible. Preserve and report them. No asset truncation, transfers, deletion,
renormalization, positive-credit fallback or relaxed feasibility gate is allowed.
A corrected baseline GE cannot be certified until this is resolved by the author.
Independent code extraction and reference-equivalence checks may continue.

## Implementation scope and ownership

Opus 5.5 owns `code/model/refactor_lab/` exclusively. Create a genuinely readable,
self-contained model package for the maintained specification: explicit inputs,
economic primitives/credit, household choices, distribution, equilibrium and a
small driver. Prefer a small number of cohesive modules and one consolidated
acceptance test suite. Parameters/data may be external; the normal execution
path must not depend on historical experiment scripts, source-string rewriting,
runtime monkey-patching or a facade that merely imports the old solver. Retain
the original engine unchanged as the validation oracle. Do not delete existing
code, move archives or edit manuscripts, shared model files or other chats' work.

Do not chase an arbitrary line count by removing required checks or economic
branches. Clearly identify any supported-scope limits. Prefer exact extraction
first, then separately review optimizations. Preserve outputs needed by the
maintained stationary calculation; do not claim transition release readiness
without corresponding verification.

Sol independently reviews the diff and tests, then Opus addresses concrete
findings. The lead reviews all borrowing/Bellman/distribution mathematics and
the acceptance evidence. Luna's bounded read-only scopes provide pointers,
not accepted scientific conclusions. First implementation/review limits are
30/20 minutes; another pass needs a specific changed scope or finding.

## Performance and acceptance

Historical completed credit GE: `../fixed_reference_economics_20260928/credit_ge_v1/solve_v1/completed.json`
records 886.527 seconds, six lifecycle calculations including a final exact
repeat, and 107 exactly repeated numeric arrays. This is a different credit
experiment, not the fixed-reference benchmark or an estimate for every GE.

Benchmark time to a verified stationary equilibrium, separately recording cold
import/JIT, household policies, distributions, equilibrium iterations, checks,
serialization and the existing 17 diagnostic plots. No calibration-normalizer
loop. No runtime claim may substitute an inner-step timing for the total.
Keep parameters, grids, tolerances and thread count fixed. Profile before
choosing a major algorithm change. Static renter allocation and saving FOC
candidates already exist; test redundant interpolation/segment bookkeeping
first. Do not replace global saving maximization with an unverified local root.

The author's latest September29 instruction permits individual model runs on
one local core, including overnight. Torch is for long batches and parallel
work. Use one numerical CPU and explicit memory and wall-clock caps; same-machine
old/new comparisons are required when moving from Torch x86 to Apple Silicon.
Before a numerical batch, smoke its exact controller,
pin source/input identities, estimate solve count/time, and retain progress and
latest/best summaries. Do not silently extend budgets or relax gates.

Acceptance requires a complete dependency inventory; effective parameter/input
identity; reference-mode comparisons of values, policies, probabilities,
distributions, moments and prices; unchanged budget, market, fiscal, population,
estate and feasibility checks; borrower/sale boundary fixtures; full 14-row
target-fit and 31-row parameter tables when reporting a computed result; and
the standard 17 plots. Record exact equivalence versus within-gate numerical
changes explicitly. Failed checks remain failures. No production promotion.

Usage starts at 1% of the account's weekly Codex allowance (Luna reserve 0%).
The intended ceiling is another 20 percentage points, with an alert rather than
an automatic stop. Account use is shared across chats; Claude Max has a separate
allowance that this Codex usage tool cannot measure.
