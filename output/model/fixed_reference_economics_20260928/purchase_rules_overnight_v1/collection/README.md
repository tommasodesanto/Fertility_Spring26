# Postchecked calibration selection

`collect.py` reads the original Torch chains, the isolated Torch restart
controller, the approved 16-slot `broader_regions_v1` Torch array, and the ten local original chains. It also reads local restart
results from `local_runtime/restart_v2/runs/chainN` when present, using the
local worker-terminal and PID-registry receipts rather than the Torch launcher
receipt. A restart contributes a
candidate only when its summary names `winner=restart`, binds the unchanged
parent receipts by hash, respects the parent's deadline and remaining 250-call
budget, and has a fresh successful selected-point postcheck. Original parent
postchecks stay in their original directories and remain eligible.
The local controller writes `restart_search_finished` while its postcheck is
running; the collector records it as pending. A failed or deadline-unverified
postcheck is visible but cannot be selected.

The broader-region array has physical paths `purchase_broader_regions_v1/results/slot_N` while its native runner keeps original chain IDs 0–7 (hard) and 24–31 (quarter). Selection records `origin=torch_regions`, `region_slot=N`, and that original chain ID separately. The collector authenticates the launched source pins and design, submission and per-stage contracts, 80-call and 75-minute limits, 05:15 cutoff, original ten bounds, targets, weights and 80% rule before applying the same fresh native 14-target/31-parameter/17-plot postcheck gates. Running and terminal-pending slots stay nonselectable. The immutable policy snapshot uses `selected_postchecks/chain_N` only as a canonical **alias** for the authenticated physical `slot_N`; it does not replace the original chain data. No result from this array is selected merely because its search has a low provisional loss.

`selected_hard.json` and `selected_quarter.json` record the winning physical
`remote_root`, `source_run`, `origin`, chain ID, full 14-row target table and
31-row parameter table, source report hashes, selected coordinates, and
target/weight fingerprints. For a restart winner they also record the parent
root, restart contract/summary hashes, and optimizer source hash. The next
stage must derive its postcheck path from `remote_root` and preserve these
provenance fields. It must not reconstruct a path from chain number alone or overwrite the parent
postcheck.

Run `python collection/test_routes.py` for a synthetic provenance and budget
check, then `python collection/test_regions.py` and `python mechanism_deployment/test_prepare_selection.py` for actual-schema slot and snapshot-routing fixtures. These tests do not run the model or connect to Torch. Run `collect.py` for an
actual read-only snapshot; use `--fetch-reports` or `--fetch-arrays` only when
the selected native files are needed locally.
