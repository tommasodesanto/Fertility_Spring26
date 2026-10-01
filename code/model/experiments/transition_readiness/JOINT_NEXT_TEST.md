# Bounded joint six-date follow-up

This follows the fixed-price test, whose input and corrected-pension housing
errors were 0.0024688895 and 0.0024695655 against 0.0002. Its corrected fiscal
error was approximately 3.5e-16. Those are failures of a fixed-price candidate;
they do not establish failure of the joint price/pension solver. The collected
endpoint at psi_child=0.1356906717749697 passed its original root, fresh repeat
and native one-step gates at price 0.7910731369532722. Its checkpoint SHA is
3d4d08e323be40f6a6c19267aeea167e79ad1133f6e78a55ecec79cc6fc00f06.

The new `joint_six_date.py` reuses that exact endpoint and unchanged reference
state. No endpoint or stationary solve is repeated. It uses the original
`e5f_four_shock_acceleration.solve_joint_with_acceleration` and its unchanged
physical residual scaling, Newton/Broyden step, projection and reproduction
gates. Initial paths follow the original `_path` linear reference-to-endpoint
guess. Every price remains endogenous; the final date is never overwritten to
force the endpoint price. The small psi increase is experimental; all other
historical earnings, entry, credit, housing, fiscal and timing objects remain
unchanged. This retains the historical repayment limitation and does not test
the current floor economy or adopt any calibration.

Required seed: the authenticated original-reference **12-date measured**
two-block Jacobian receipt and its matrix. The original estimator's scientific
provenance checks are reused. The saved matrix must equal the original exact lag
reconstruction at 12 dates; the same helper constructs the six-date approximate
seed. Its scaled condition number must pass the original threshold, so a
guessed-slope first step cannot replace the measurement. Reuse at psi ×1.001 is
a numerical warm start; native nonlinear mappings and exact replay still decide
certification. A missing seed or unreviewed source mismatch stops before model setup.

The budget is 1950 seconds including setup, mappings and plotting, one numerical
thread and 64-GiB cache. Three maps imply at most 36 conservative backward/forward
policy calls: input, at most one joint step, then fresh replay only when a
candidate clears both physical gates. Each native map is capped at 600 seconds,
with remaining-time checks reserving future maps and 120 seconds for the stable
17-plot packet. The original generic root would otherwise replay an unqualified
best candidate: this driver's evaluator refuses that third call, without
changing the root mathematics. No additional numerical step, automatic restart,
104/128-date run or fit is authorized by preparation.

Output separates `root_certified` from `terminal_pass`. All original terminal
distribution, population, price and birth-queue tolerances remain 1e-3, including
the raw queue. Terminal metrics are recorded for every completed candidate,
but a terminal failure becomes interpretable as a horizon failure only after the
joint root clears its physical and fresh-replay gates. Three dates (0,3,5) are
captured and the same 17 diagnostics are rendered from the last completed map.
On later failure, its checkpoint and captures remain and plotting is attempted
only within the deadline. Heartbeat is every 30 seconds; latest/best/start/root,
checkpoint and complete/failure receipts persist. Neither a six-date root pass
nor a terminal pass substitutes for a full horizon comparison.

Config fields are exact pinned endpoint/seed receipt records, `schema` equal to
`historical_joint_six_date_v1`, `runner_sha256`, `total_seconds=1950`, `horizon=6`
and `maximum_maps=3`. The endpoint's own checkpoint pin and seed's own matrix pin
must authenticate. Relative paths are resolved inside the unchanged original
project root, which allows reuse of the real remote receipt without rewriting
its content.

```sh
python code/model/experiments/transition_readiness/joint_six_date.py \
  --config output/model/transition_readiness_v1/joint_preparation/config.json \
  --config-sha256 VERIFIED_CONFIG_SHA256 --preflight
python code/model/experiments/transition_readiness/joint_six_date.py \
  --config output/model/transition_readiness_v1/joint_preparation/config.json \
  --config-sha256 VERIFIED_CONFIG_SHA256 \
  --output output/model/transition_readiness_v1/joint_six_date_run
```

The concrete config pins the real remote seed receipt (SHA d6dfbf6eedbb23ebefea681dccb2b098ea2514dcc407e0ab8cc9c89f6237974d) and matrix (SHA b48bd10fa7ac2806a4db34b5eb373b8b22748b0160b1b08fbccf9e8a1d176860). No synthetic seed or production manifest was created. Focused tests use the actual
unchanged joint root with synthetic linear residuals and stubbed native callbacks
to validate orchestration only: one step plus exact qualifying replay, refusal
of an unqualified replay, time reserves and three-map cap. They run no lifecycle
solve and carry no scientific validation claim. The original 18 frozen files,
execution pins and failed-case packet remain intact.


The lead reviewed and approved three source provenance differences explicitly:
cache d51bbd…→1e8a02…, estimator e39873…→bd9382…, and shock-fit 50a91b…→495b61….
The six other measured-map numerical modules match exactly. The cache proof
checks that current bytes equal the archived seed-era bytes after precisely two
`except TimeoutError: raise` insertions; successful cache/mapping computations
are unchanged. The estimator proof compares ASTs for the actually used
`draft_plan`, `_reconstruct_seed` and `initial_jacobian`; initializer bookkeeping
adds only `self.stage_deadline=self.deadline`. Original archived nine-source
provenance validation is invoked with a separately labeled seed validation plan.
Runtime source authentication still checks the frozen current package. The
shock-fit source is archived/pinned but its fit functions are never called;
endpoint/path/fit estimator routines are also never called by this follow-up.
No seed receipt byte, measured derivative or runtime cache was replaced.
`seed_compatibility.json` records every differing pin rather than claiming all
nine runtime sources match. The actual zero-solve reconstructed matrix matches
exactly; six-date measured lag coverage is complete and scaled condition number
is 213.2244369653 against the original 1e8 ceiling.

`joint_preparation/run_joint.sh` is a concrete 4-CPU allocation, 96-GiB, 33-minute launcher,
with the numerical driver capped at 1950 seconds. It uses one numerical thread; four allocated CPUs satisfy the scheduler memory rule. Five focused local control-flow/provenance tests passed; the launcher runs authenticated metadata preflight before native setup and maps. It uses a new `source_joint/transition_readiness` snapshot and preserves the earlier execution package. The preflight explicitly reports `native_setup_verified=false`;
loading the actual reference/endpoint state is still a separate runtime gate.
Preparation did not submit this launcher or run any model. Additional native
execution remains subject to the human's requested budget extension.

## Executed result and unlaunched numerical continuation

The first joint package was subsequently approved and executed as job 18926856.
Two native maps completed in 11:02. Housing passed after one damped update;
fiscal error 1.44704e-5 still failed the 1e-6 gate, so no qualifying fresh replay
ran. See output/model/transition_readiness_v1/README.md for verified accounts.

The evolving driver also supports the explicit full_step_1200 preset. It starts
from the authenticated previous best price/pension paths and the actual Broyden
Jacobian. A fresh map must first reproduce the saved residuals; one full Newton
step follows, with a fresh replay only if the updated candidate clears both
gates. The numerical proposal uses damping 1, a 1200-second total cap, at most
three maps and 360 seconds per map. Original gates, endpoint, inherited
population, fiscal/housing closure and other economic primitives remain fixed.
The original default retains damping 0.7, a 1950-second total cap and 600 seconds
per map. Archived executed sources retain the prior source and configuration.

The exact proposed config and launcher are in
output/model/transition_readiness_v1/joint_full_step_preparation/. Eight
zero-lifecycle tests and the lead's actual Broyden reconstruction check passed.
Native reproduction remains unverified. No full-step job was launched; the
latest sleep/coordination instruction permits no further budget extension.
