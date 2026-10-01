# Independent Fable calibration audit

Author-requested read-only review of the floor economy, joint psi calibration,
weights, entry wealth, childbirth timing, identification and fast search strategy.
The exact scope is in `prompt.md`. Existing first-party Claude Max authentication
is required. Fable is requested with max effort; only Read, Grep and Glob are
enabled. Safe/restricted mode and strict MCP configuration prevent hooks,
command execution and external tools. The supervisor has a hard 30-minute
deadline and no automatic retry. No model runs or edits are authorized.

`launch.json` records the actual initialized model, session and PIDs;
`progress.json` tracks response progress. `completion.json` or `failure.json`
records termination. `final.md` holds the final audit, or labels a partial result
at timeout. `stream.jsonl` and `stderr.txt` preserve execution evidence.
