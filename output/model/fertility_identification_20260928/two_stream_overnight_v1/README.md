# Bounded two-stream overnight launcher

This is an experimental calibration and identification diagnostic, not a
production/adopted calibration.  It keeps the approved original fourteen-row
target system, weights, and bounds.  The reference is **2007 stationary
reference — block0506, September 28 verified export**; its actual normalized
`psi_child` is distinct from both stream seed values.

Each one-CPU, 24GB stream starts with its authenticated anchor, then (time
permitting) runs the anchor replay, two rounds of ten local Jacobian probes and
three ridge-Gauss--Newton trials, seven deterministic explorations, and two
final repeats.  Each case is capped at 2,100 seconds and eight stationary
solves.  Each stream has 25,200 seconds total, including a 4,200-second final
repeat reserve, and at most 36 objectives.  The hard end is 2026-09-29
11:30 UTC; time and scientific failures take precedence over completing the
case count.  Two consecutive admissibility/resource censors stop the chain;
all other unknown failures are fatal.  There is no resume or automatic restart.

Run-size estimate: at most 72 candidate evaluations across the two streams,
with an absolute ceiling of 576 stationary solves. Recent observed stationary
solves took roughly 170–242 seconds, before export overhead. One-solve warm
evaluations took 272–318 seconds and a three-solve evaluation took 580 seconds;
36 such evaluations would take roughly three to six hours per stream. More
normalization iterations can exhaust the seven-hour budget sooner, so the case
ceiling is not a completion promise. Every case is independently supervised;
heartbeat is written every 15 seconds and latest/best summaries after each case.

The starting checkpoints are the original-model numerical-pair retained-start
candidate and the normalized two-birth v2 diagnostic, respectively. These are
experimental search seeds; neither replaces the named stationary reference.
Preparation records their exact checkpoint, parameter, target and source
identities before dispatch.

The overnight follow-up is read-only monitoring for major failures. It may
inspect Slurm, heartbeat, latest/best and failure receipts, but must not repair,
retry, restart, extend budgets, change specifications, or promote a candidate.
Older completed search automations remain paused.

The one-birth stream uses the verified original model.  The two-birth stream
alone installs the tested immutable `two_births_optimized_v2` overlay in its
isolated worker: after a successful birth it permits one extra attempt, so at
most two births can occur in one four-year cell; it uses the existing
later-birth Gumbel scale/inclusive value and an independent age-specific second
conception draw.  Common-event linear within-period age projection remains;
separate birth dates are not modeled.  Earnings, entry distributions,
transfers/floors, housing/credit primitives, targets and weights are unchanged.
In both lanes `psi_child` is separately normalized to completed fertility 2.1
and the existing demographic-renewal gate remains binding.

Every accepted case must export the full 14-row target fit, frozen 31-name
parameter table, and all 17 standard plots.  The controller cross-checks the
worker success receipt, native `case/receipt.json`, scientific identity, and
selected checkpoint hash.  The actual selected checkpoint is hashed once on
Torch before final repeats.  Jacobian rank is reported only as local numerical
sensitivity: a full numerical rank does not establish statistical
identification, and a deficient rank is reported as provisional/underidentified
locally without dropping moments or parameters.

Before a `search` invocation, the parent must create a source-pinned
`launch_approval.json` with schema `two_stream_launch_approval_v1`.  It must
pin the exact configuration, both lane source fingerprints, one typed synthetic
receipt, and two typed lane-specific integration receipts (each containing two
real dispatcher evaluations).  This packet intentionally never creates that
approval.

## Preflight evidence

The initial controller test job 18764691 passed 15 of 16 tests; immediate
post-SIGKILL inspection of a descendant raced its observed death. The revised
test allows up to two seconds to observe missing/zombie state and still fails
on a live descendant. Controller behavior was unchanged. The first preparation
18764700 and its source/configuration are preserved on Torch under
`failed_preflight_v1/`; they cannot authorize this search.

Fresh controller tests 18765038 pass all 16 tests in 88.857 seconds. Fresh
preparation 18765039 passes, with zero model solves and configuration SHA
`419b46d7cffdf51d95760f2999aec67a66b10b08bb9198164ac3d6f7d96323d3`.
The exact config contains every source/input pin, both seed checkpoints and
parameter vectors, the original bounds and complete target fingerprint.
`synthetic_receipt.json` records the lead's verification of the completed
Torch test job. Integration array 18765327 must also pass before launch.

Small launcher edits used the repository Terra medium worker profile and an
explicit Terra medium subagent; the lead reviewed their diffs and ran Torch
verification. No local model import, numerical test or rendering was used.

The time reserve protects repeat capacity but does not guarantee two successful
repeats at the maximum case duration: controller bookkeeping also uses time.
If the remaining budget cannot accommodate a full case, the run reports an
incomplete repeat packet. It must not be described as a verified final result.

## Launch and scheduled monitoring

Integration array **18765327** passed both repeated evaluations in both streams:
610.979 seconds for the original model and 524.593 seconds for the two-birth
model, one stationary solve per evaluation. The original-model maximum weighted
residual difference from its anchor is 7.362e-05; the two-birth difference is zero.
Both are within the unchanged 0.01 screen, with loss differences below 0.05.
Full first-evaluation fits and all parameters are collected under
`integration_smoke_v1/{one_birth,two_birth}/first_case/`; all four native 17-plot
packets and checkpoints remain on Torch. This is a reproducibility smoke, not
new calibration evidence or an adopted specification.

The source-pinned readiness approval is `launch_approval.json`, SHA256
`77be2cf5e988f754ef5650f5d1a4ae2626f9b77165937cc6387fd5ff5c91a32a`.
Both searches launched as Slurm array **18766206** at 23:24:50 EDT September 28:
task 0 original model, task 1 experimental two-birth model. Both were confirmed
RUNNING, with fresh heartbeats and the initial replay dispatched. Their current
seven-hour clocks end around 06:25 EDT September 29; the separate hard end
remains 07:30 EDT. Exact clocks live under `run_v1/<lane>/clock.json` on Torch.

The remote packet is
`/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fertility_identification_20260928/two_stream_overnight_v1`.
App heartbeat `monitor-two-stream-fertility-calibration` checks every 30 minutes,
with one read-only inspection per wake-up and new-major-failure reporting only.
It pauses after both jobs terminate or by 12:30 UTC September 29. It does not
repair, restart, retry, cancel, extend budgets, change targets, promote candidates
or message other chats. Old completed search monitors remain paused.

The existing search command `run.sh search` invokes the source-pinned
controller only with its exact configuration and live readiness approval; it
must not be rerun to resume this job. Each owned worker invocation writes the
full diagnostic packet through the native `graphs=True` evaluator. The saved
request and source/configuration pins recover the command for a separately
authorized future replay. No new replay is authorized by this documentation.
