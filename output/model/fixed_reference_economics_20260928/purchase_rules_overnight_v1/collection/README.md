# Postchecked calibration selection

`collect.py` reads the original Torch chains, the isolated Torch restart
controller, and the ten local original chains. It also reads local restart
results from `local_runtime/restart_v2/runs/chainN` when present, using the
local worker-terminal and PID-registry receipts rather than the Torch launcher
receipt. A restart contributes a
candidate only when its summary names `winner=restart`, binds the unchanged
parent receipts by hash, respects the parent's deadline and remaining 250-call
budget, and has a fresh successful selected-point postcheck. Original parent
postchecks stay in their original directories and remain eligible.

`selected_hard.json` and `selected_quarter.json` record the winning physical
`remote_root`, `source_run`, `origin`, chain ID, full 14-row target table and
31-row parameter table, source report hashes, selected coordinates, and
target/weight fingerprints. For a restart winner they also record the parent
root, restart contract/summary hashes, and optimizer source hash. The next
stage must derive its postcheck path from `remote_root` and preserve these
provenance fields. It must not reconstruct a path from chain number alone or overwrite the parent
postcheck.

Run `python collection/test_routes.py` for a synthetic provenance and budget
check. It does not run the model or connect to Torch. Run `collect.py` for an
actual read-only snapshot; use `--fetch-reports` or `--fetch-arrays` only when
the selected native files are needed locally.
