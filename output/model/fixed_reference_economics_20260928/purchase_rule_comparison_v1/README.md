# Purchase-rule comparison

`build_readout.py` checks the saved strict-80 and quarter-saving-80 native reports and renders the complete 14-moment and 31-parameter tables. It also accepts the isolated hard-100 and quarter-saving-100 retry completion receipts through `--hard100` and `--quarter100`. Missing reports remain visibly blank, and failed first attempts are disclosed. The script reads saved artifacts only and makes no model solve.

Run `code/model/.venv/bin/python build_readout.py` here after collection. The default paths point to the two isolated solvency retries; use the CLI overrides if a retry is collected elsewhere. The generated `readout/verification.json` records source hashes, gate states, and arithmetic checks.
