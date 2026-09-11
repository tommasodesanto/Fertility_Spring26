# Stationary old/new utility test at the retained 2023 preference

Job17373666: two parallel arms, old utility and new parenthood-only utility.
Both balance the payroll/pension budget, fix payroll tax0.179, and retain the same
inherited structural parameters, source snapshot, initial supply curve and grid.
The original selected2023 preference is-0.03200017062648686; no2.1 normalization
is applied. The stationary age distribution is endogenous to the stationary
operator; it is not the observed-age2023 historical distribution.

This is a controlled stationary counterpart, not a re-estimated2023 calibration
or a reproduction of the old transition. Each arm repeats twice,4GE total,
roughly240 seconds per arm anticipated,720-second internal loop limit and
840-second Slurm limit. Each case writes progress during GE, latest/best results,
full original12-row target fits with unchanged weights, parameter restrictions,
checkpoints and17 standard diagnostics. There is no automatic follow-on search.

`run_comparison.py` validates all641 frozen source hashes and the original2023
summary/target fingerprint before invoking the unchanged solver driver. It uses
the existing TargetSystem loss function and retains both repetition tables.
`manifest.json`, `contract.json` and `submit.sh` completely specify the run.
The original early-target objective is not used. The legacy observer/sample
approximations remain those of the old target system; no empirical rows change.

Remote bundle:
`/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/static_2023_utility_comparison/`.
Results appear in `results/old_balanced/` and `results/new_balanced/` there.
No previous output or model source is modified.

## Completed result

Both arms completed in about151 seconds and reproduce exactly. All74 lightweight
result hashes match the remote collection manifest; all34 PNG hashes also pass
the review builder. Relative pension residuals are below2e-9. Full12-row fits and
parameter restrictions are local in `results/old_balanced/` and
`results/new_balanced/`. The old/new losses1335.71971 and1361.32977 are fixed-point
diagnostic scores, not estimates. The stationary fertility stocks1.39999/1.40925
show why the retained historical2023 preference cannot be mistaken for a2023
stationary calibration. No automatic further search runs.
