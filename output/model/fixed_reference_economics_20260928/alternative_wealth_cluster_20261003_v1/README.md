# Alternative timing, model-matched wealth, widened-beta cluster input

This packet prepares ten one-core Torch chains; it does **not** record a
production cluster result. The base economics use soft financing and housing
transactions after interest. The experimental target replaces only the
`wealth_earnings` PSID 2005/07 pooled ratio, from 6.92658379107299 to
4.45838713455674 (net worth excluding business/farm equity, other real estate,
and vehicles, over gross head-plus-spouse earnings). The original numerical
weight 7.595098472533724 is retained as a controlled sensitivity weight,
not a new precision estimate. The old target remains a native diagnostic;
`target_fit_new_contract.csv` is authoritative for experimental scores. The
bequest target and all other economic, observer, entry, earnings, and floor
objects are unchanged. The bequest/wealth denominator compatibility remains
unresolved; informative rank for 10 scored moments and 10 free parameters is
unverified.

Relative to the verified local new-wealth search, this variant changes only
the annual $\beta$ search bound to $[0.930,0.990]$, the objective-call cap to
500, and the ten starting points. The upper bound and other nine parameter
bounds are unchanged. Each chain has a six-hour wall cap including a
1800-second native verification reserve, and uses one BLAS/OpenMP/Numba
thread. No solver acceptance threshold is relaxed. The 500-call cap is an
upper bound; at recent solve speeds the wall cap is likely to bind first.

`start_plan.json` has SHA-256
`752d1f02b414745cb391a499d5d519811f897ea56cee4c28ef1547472a7e47e1`.
It pins verified new-wealth local chain 02 as the anchor. The ten annual-$\beta$
seeds are 0.964494931582842, 0.970, 0.966, 0.960, 0.955, 0.950, 0.945,
0.940, 0.935, and 0.930; the other nine seed coordinates equal that anchor.
The target fingerprint is
`c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`;
the target-and-weight fingerprint is
`f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`.
The driver SHA-256 at preflight is
`b08ad1f1436f6ff3c498ed13527c59cc36a1b16340d18ea9c557a5dc12e049c5`.

Local validation used the exact driver loop. `mock_low_beta/` confirms the
bound-constrained optimizer path at seed $\beta=0.930$ without model solves.
`preflight_beta930_child/` confirms the fresh native child accepts the same
bound and target contract without model solves. `native_smoke_chain0/` contains
two completed full-GE objective calls and a fresh selected/repeat native
postcheck: 14 target rows, 31 parameter rows, 17 standard PNGs, identical
experimental target-fit tables on repeat, and exact search/native new-target
rows. The native `parameters.csv` reports annual-$\beta$ bounds
$[0.930,0.990]$. The smoke status is `selected_numerically_verified`; it does
not establish optimizer convergence or rank identification.
