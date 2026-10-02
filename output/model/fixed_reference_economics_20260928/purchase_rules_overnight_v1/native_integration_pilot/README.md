# Provisional native integration pilot

`prepare.py` freezes the lowest loss currently postchecked checkpoint per arm
from the read-only collector into this folder's `manifest.json` and
`selected_{arm}.json`. These are **pilot inputs**, not final calibrated
baselines or production selections. The frozen October 2 pilot selected hard
Torch chain 9 and quarter local chain 55.

`launch_torch.sh` submitted the first one-date controls under job IDs 19017655
and 19017678. Both failed before a model solve because the then-current
constructor tried to read full CSV reports from the deliberately compact saved
`selected_repeat` folder. Their complete failure receipts remain under
`/scratch/td2248/projects/purchase_native_integration_pilot_v1/results/`.

`launch_torch_retry.sh` uses the lead-reviewed compact-repeat authentication
source as a hash-pinned, read-only file override. It otherwise uses the same
staged engine, model sources, selected checkpoints, one-date control,
one-core/32 GiB/20-minute limits, and 64-policy-call cap. The retry is job
19018245, array tasks 0 (hard) and 1 (quarter), with separate results under
`/scratch/td2248/projects/purchase_native_integration_pilot_v2/results/`.
No final policy or calibration result should be inferred from a pilot smoke.

The v2 retry failed before solving because the dated integration imported the
legacy model before native authentication. The lead approved deferring only
that import. `launch_torch_retry3.sh` binds hash-pinned copies of the corrected
`run_case.py` and `selected_runtime.py` into a third isolated pilot source;
the prior v1/v2 sources and results remain intact. Hard job 19020019 and
quarter job 19020025 use `/scratch/td2248/projects/purchase_native_integration_pilot_v3/`.

The v3 hard pilot reproduced the selected stationary solve twice but then
stopped at a missing import in the copied dated path. The v3 quarter pilot
passed the numerical reference tables and closure, then stopped because its
saved plots were rendered with Matplotlib 3.7.1 while the fresh Torch plots
used 3.10.0. The reviewed v4 source authenticates all 95 saved solution and
observer arrays to the original numerical tolerance, keeps exact same-renderer
plot checks, and reads renderer identity from each PNG. It also applies the
existing typed stationary audit gate. `launch_torch_retry4.sh` freezes these
three corrected source files in
`/scratch/td2248/projects/purchase_native_integration_pilot_v4/`; hard job
19022024 and quarter job 19022025 are isolated one-date control checks on the
same provisional points. Their results do not select final baselines.
Both v4 checks passed with three native policy calls each. Both fresh replays
matched all 95 saved fields; maximum absolute difference was zero for hard and
`6.37e-12` for quarter. The dated one-period control mapping passed for both.
