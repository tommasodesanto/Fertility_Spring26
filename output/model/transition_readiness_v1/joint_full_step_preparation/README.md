# Full-step joint continuation preparation

The previous joint attempt completed two fresh maps in 658.54 seconds. Its
second map cleared the housing gate (0.00016565 <= 0.0002) but failed the fiscal
gate (0.00001447 > 0.000001), so it correctly refused an unqualified replay.

This preparation adds two closed numerical presets to the evolving driver.
The existing default remains damping 0.7, total 1950 seconds and mapping cap
600 seconds. The explicit full-step continuation uses damping 1, total 1200
seconds and mapping cap 360 seconds, reserving 120 seconds for plots. Economic
primitives, endpoint, derivative source compatibility, original root/scaling,
physical gates, terminal tolerances and both queues remain unchanged.

The proposed new test uses the authenticated previous best coordinates and
its genuine final Broyden matrix. The matrix reconstructs from the unchanged
measured seed and two actual native histories with maximum difference
1.82e-14. Source/config/reference/preflight/root/latest artifacts are pinned;
first fresh input must reproduce both residual blocks within 1e-10. Only then
can the original root take one full Newton step and a qualifying fresh replay,
with at most three fresh native maps. There is no cached initial evaluation,
additional iteration, guessed derivative, population rescaling or endpoint-price
overwrite. Root and terminal certification stay separate. Neither is full
horizon certification or validation of the current floor economy.

Eight focused zero-lifecycle tests pass, including default behavior, full-step
linear root/replay, actual saved Broyden reconstruction and provenance failures.
The launcher uses a new source_joint_full_step snapshot, 20 minutes, four
allocated CPUs but one numerical thread, 96 GiB memory and 64 GiB policy cache.
No job or restart was performed during preparation. Native behavior and timing
of this proposed continuation remain unverified, and execution awaits explicit
lead/user authorization. Earlier source snapshots and the original 18 frozen
files are preserved.
