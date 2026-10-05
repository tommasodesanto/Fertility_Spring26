# Two isolated overnight calibration continuations

Submitted October 4, 2026, 10:49 PM New York. Each arm has 12 distinct
one-core Torch tasks (`0-11%12`), a 12-hour wall limit, at most 500 objective
calls per task, a 1,800-second reserve for unchanged fresh native selected-point
and exact-repeat checks, and no automatic restart. The exact two-case launcher
loop passed with two objective calls, zero lifecycle solves and exit 0 in each
arm. This execution mock does not validate a numerical result. Search-best
losses remain provisional until the native gates pass.

| Arm | Economic contract | Target/weight SHA-256 | Verified anchor | Job |
| --- | --- | --- | ---: | ---: |
| A | At most one intended birth, Estate A, new wealth target 4.45838713455674, post-interest soft financing | `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70` / `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4` | 21.275413361071312 | 19194495 |
| B | Original fertility architecture and at-death mapping, no Estate A, old wealth target 6.92658379107299, post-interest soft financing | `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1` / `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0` | 13.771131463467462 | 19194496 |

The two losses are across different contracts and should not be ranked against
each other. Within each arm, earnings, entry distributions, transfers, floors,
preferences, ten free coordinates and bounds, 14 target rows and weights, and
native numerical gates are retained from its named reference. The new code
only repairs orchestration: recognized native `uncomputed_price_unbracketed`
and typed infeasible-point exceptions become rejected cases without a finite
valid loss, allowing other starts to continue. A penalty used inside an
optimizer cannot be selected as a calibrated solution. All other runtime or
contract failures remain fatal. No price cap or economic mechanism changed.

`build.py` regenerates each stage by overlaying the narrow driver and start-plan
changes on the pinned prior deployment archives. A uses ten earlier feasible
continuation starts plus two canceled provisional full-GE points; B uses ten
verified alternative endpoints plus two verified prior-chain endpoints. Each
arm's own SHA-checked receipts and bounds authenticate these twelve distinct
starts. `submit_once.py` checked the exact mock, source manifest and queue,
then wrote a durable exactly-once submission receipt. It uses a delayed Slurm
start so the receipt is on disk before production tasks begin.

Remote stage roots are `/scratch/td2248/projects/estate_birth_overnight_20261004_v1`
and `/scratch/td2248/projects/soft_timing_overnight_20261004_v1`. Local source
plans, manifests and copies of submission receipts are in
`output/model/overnight_dual_20261004_v1/{a,b}/deployment/` (generated output).
Manifest SHA-256: A `d86616474ca232071305fc7d877df33c24917d49898a461ab0c4dc1a5802fe23`,
B `cef5869f2e18079a08134a32f12cb706a7b1613396b8f6b7c3f77ee651108171`.
The updated shared heartbeat `monitor-one-birth-estate-a-calibration` watches
both arrays and stays quiet while healthy; it never retries or submits.
