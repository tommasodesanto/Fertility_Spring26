# CES normalized-share overnight calibration

This is a prepared, non-adopted four-chain bounded Nelder--Mead search for the author-authorized normalized CES-share experiment. It starts from post-interest soft chain 13 and retains the old wealth target (6.92658379107299), timing, earnings, and all other inherited economics. The experiment has 14 target rows, 11 scored moments, and 11 free coordinates: the nine inherited coordinates plus `delta_alpha_jump` and `delta_alpha`, each in `[0,.25]`. The added `family_rooms` target is 0.38509964969278165 with weight 280.52808370152104.

The share rule is `alpha(m)=.733` when childless and `clip(.733-delta_alpha_jump-delta_alpha*m,.05,.95)` for parents. The normalized material denominator applies in all states and `h_P=0`. The existing birth menu, raw utility costs, no-estateA restriction, grids, fixed inputs, timing, and earnings remain unchanged. There is no `r*` correction or `alpha0` numerator, and no added birth shock or cost rescaling. The `family_rooms` weight is the inherited 42-metro bootstrap weight because national uncertainty is unavailable; the model-dependent-child observer remains a proxy.

The v1 and v2 zero-solve preflights failed because historical packaging dependencies were missing. V3 attempt 3 was built, but Torch SSH authentication expired during transfer; it requires NYU/Duo reauthentication before staging can resume. No Slurm job was submitted. No calibration or adoption has occurred. Stage and submit only after the v3 zero-solve context preflight and exact-loop smoke are recorded as passed. The native selected-point postcheck is configured to require all 11 coordinates, the full 14-row experimental target CSV, 31 parameter records, 17 plot hashes, and an exact repeat. The prior v1/v2 packages and receipts remain historical evidence.

Per chain: CPU 1, 24 GiB, six hours, no more than 500 objective calls and 32 lifecycle solves per GE, with 1,800 seconds reserved for a fresh-child selected native postcheck. The nominal maximum is 2,000 GE / 64,000 lifecycle solves, although wall time will bind earlier. Chain-0 smoke runs exactly two distinct objective points followed by the fresh native selected postcheck and exact repeat; it is capped at 1.5 hours. Production is guarded by `verify_smoke_gate.py`, then submits `0-3%4`. Errors are terminal: no retry, fallback, adoption, or overwrite.

Lead commands, only after the adapter is generated and reviewed:

```bash
python3 code/cluster/ces_normalized_shares_calibration/prepare_plan.py
PYTHON=output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python
for chain in 0 1 2 3; do "$PYTHON" code/model/experiments/ces_normalized_shares/calibrate.py --chain "$chain" --out "/private/tmp/ces_mock_${chain}" --deadline-epoch $(( $(date +%s)+60 )) --mock-smoke; done
bash code/cluster/ces_normalized_shares_calibration/stage_torch.sh
ssh torch 'cd /scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3 && sbatch --parsable --array=0 --time=01:30:00 --export=ALL,CES_RUN_MODE=smoke launch_torch.sh'
ssh torch 'cd /scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3 && /share/apps/anaconda3/2025.06/bin/python verify_smoke_gate.py'
ssh torch 'cd /scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3 && ./submit_torch.sh'
```

Each selected native postcheck runs one fresh GE in a child process. The native solver itself performs the exact selected repeat, including arrays; the child independently requires the root and internal-repeat 14-row target CSVs and 31-row parameter CSVs to match exactly, all 17 diagnostic PNG hashes to match, normalized-utility contract pins to match, and the eleven selected coordinates to equal the root parameter table within their approved bounds. The smoke gate and collector reject any mismatch in target, weight, start-plan, selected-source, source-checkpoint, or stage-inventory identity. `collect_torch.py` exits nonzero for every incomplete or mismatched chain and writes `RESULTS.md` only as a readable non-adoption receipt.

The immutable v3 stage retains its bundled collector. Use the separately hashed `followup_tools/collect_torch.py` overlay for host-side collection: it maps container result paths to the chain launch folder and checks all 17 diagnostic hashes. This tool-only correction changes no solver, utility, target or launch contract.
