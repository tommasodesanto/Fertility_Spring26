# Five-date PAYGO checkpoint comparison — prepared, not launched

This read-only packet reconstructs the 2007 distribution from the original verified initial checkpoint using the completed historical run's exact age factors. It reuses the actual age reweighting function, fiscal accounting, and five-date prefix helper. No GE, Bellman, empirical reload, 2027 extrapolation, or root modification occurs.

Frozen source: `527ab397218e18611ba6930e62283eefc4225688`; **643 complete code/model files, 539 Python files**, hashed from committed archive bytes. The contract additionally pins this runner, the original checkpoint, historical root, preflight, age CSV, and historical contract. Source hashes do not use the changing working tree.

Stage the frozen commit's complete `code/model/` under `/scratch/td2248/projects/Fertility_Spring26_paygo_prefix_527ab397` and copy this entire folder to `/scratch/td2248/projects/Fertility_Spring26_paygo_prefix_527ab397/prefix_comparison/`. The four small inputs are bundled under `inputs/`; the initial checkpoint remains at its original pinned cluster path. The launch script embeds contract SHA256 `698b77087e08964c04d4331dad5c99f9fdec063d04624b5d9ee3ea37a380ddc1`. Any edit requires refreshing that hash; runner edits also require refreshing its contract hash first.

Prepared Slurm budget: **1 CPU, 8 GB, 5 minutes**. Lead must review and launch `submit.sh`; no job has been submitted. The runner fails before unpickling unless every source/input pin matches. Existing output directories are refused. Outputs: `prediction.json`, `comparison.json`, `predicted_marginals.npz`, `summary.json`, and `heartbeat.json` in a job-specific output directory.

Numerical gates: five-date fiscal and factor comparison at 1e-6 or tighter; full 2007 preflight factor reproduction and absolute age-mass checks at 1e-12. Actual trial pension balance is reported separately from the proposed pension's balance on actual ledger factors. The completed historical trial remains uncertified regardless of this comparison. Full 2011–2023 age/income marginals were not saved and cannot be certified here. The reconstructed full 2007 pre-choice distribution is compared with preflight; root fiscal ledgers reflect actual dated post-choice heads.

Preparation checks: Python AST parse, shell syntax, 85 unique date/age CSV rows, original checkpoint linkage, and frozen source inventory count. Full-checkpoint execution awaits the cluster; no local model run was performed.
