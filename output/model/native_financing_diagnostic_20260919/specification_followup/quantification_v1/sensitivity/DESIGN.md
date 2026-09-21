# Controlled calibration sensitivity design

Status: specified by lead; implementation and launch gates pending. No job
submission is implied by this design.

Question: at the selected new-income fit, which structural parameters move
the housing, fertility and wealth residuals together or against one another,
and are the local responses stable when perturbations are halved? This is a
controlled diagnostic of the calibration map. It does not validate the income
measurement contract, prove identification or determine a globally attainable
fit. It is distinct from the earlier coordinate poll around the worse pilot
and from the adaptive multivariate proposal search.

Anchor: case 60 of the completed 96-proposal search, selected objective
353.6588729140903; original initialization psi 0.23949950404168222. Verified
reference is `income_overnight_v1/production/selected_verification/evaluation`
on Torch under the existing native-financing root. Use full-precision parameter
values from its score/summary, not rounded prose. Both fresh anchor repetitions
must reproduce the selected target/parameter tables, objective and native
arrays before dependent production.

Keep the frozen economic source, gross-earnings candidate (15 joint states),
entry rule, full target definitions and weights, structural bounds, numerical
settings, stationary equilibrium procedure and separate fertility normalization
at 2.1 unchanged. Each point re-solves equilibrium and derives its own child
preference scale through that same normalization; prices and the stationary
distribution are endogenous here. There are no policy arms.

Full increments: annual beta 0.001; first-child housing requirement 0.05;
first-birth fixed cost 0.01; each of the other six coordinates 5% of its anchor.
Use both signs if inside the actual search bounds. Beta and the housing
requirement are at their upper bounds, so use inward one-sided perturbations
and omit the infeasible outward directions without clipping. This gives 16
full-step points. Repeat at half step for beta, the housing requirement,
ownership preference chi, housing supply level H0 and first-birth taste scale
kappa_fert: eight additional points. Total: 24 probe evaluations.

Smoke the same controller/evaluator loop with two fresh anchor evaluations and
the full inward-beta probe. Production reuses that probe, evaluates the other
23 points, freezes the lowest original-objective valid candidate including the
anchor, and verifies it in two exact fresh repetitions. Maximum: 28 new full
objective evaluations, including smoke and final repetitions. Each full
objective contains multiple stationary solves; record their actual count and
the adapter's enforced maximum separately in the launch manifest.

Observed case-60 runtime is 530 seconds; its earlier two-repeat verification
was 1,043 seconds. With eight single-threaded workers, the 23 production probes
need three waves; allow 900 seconds per full evaluation and 1,800 seconds for
the two final repetitions. Planning estimate: roughly 50–75 minutes of running
work including smoke, plus queue waiting; reserve one hour for the smoke and
three hours for production as hard Slurm caps. Absolute numerical deadline
07:30 America/New_York September 21; no new batch after 02:00.

All source/input/target checks must pass before solving. Stop on unknown,
source, target, budget, accounting or validation failures. Only narrowly
recognized numerical infeasibility/nonconvergence may be retained as a rejected
probe; then the affected local sensitivity is missing, never imputed. A failed
smoke blocks production. No automatic retries or relaxed gates.

Keep 30-second active-case receipts, latest completed and best-so-far summaries,
all full target/parameter tables and the unchanged 17 plots for each completed
case. The final report may compute finite-difference weighted moment responses
with explicit step/coordinate scaling, one-sided columns and half-step checks.
Rank and condition numbers describe this local finite-difference map only;
unstable steps, missing points and normalization-induced responses remain
visible. No proposed target deletion or parameter fixing follows from rank
alone. The lead reviews the frozen controller/launcher before submission.
