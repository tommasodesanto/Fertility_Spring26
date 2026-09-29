# Launcher contract fixes

- `worker.py` reads the exact `reference_parameters.csv` pin key written by `prepare.py`.
- Synthetic worker successes now include a checkpoint and matching native/scientific receipts, so the production `authenticate_best` checkpoint path is exercised.
- The launch-approval fixture supplies both lane source fingerprints required by `verify_launch`.
- The owned-timeout test waits at most two seconds for descendant reaping after the unchanged process-group `SIGKILL`; it still accepts only an absent or zombie descendant.
