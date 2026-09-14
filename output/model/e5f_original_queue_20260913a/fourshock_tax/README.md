# Four-shock inherited-2023 tax comparison

Update: recovery 17723529 reproduced the five native rows and saved the 2023
snapshot, then failed while serializing NumPy arrays in `fertility.json`.
The array/scalar serializer and atomic replacement were tested. Recovery retry
17734560 is submitted in sibling immutable batch `fourshock_tax_20260914b`;
it repeats the native recovery checks and then dispatches the 2% policy.
The existing baseline transition is reused. Retry receipts are saved here.

Recovery job 17723529 waits for baseline job 17711519 (`afterany`). It freezes
a completed 104-date native baseline evaluation, preferentially the saved
best, and replays the household problem to recover its own 2023 distribution
and original birth queues. A fallback to the latest completed evaluation is
explicitly labeled. No new 1% equilibrium path is solved.

After exact row reproduction and native household checks pass, the recovery
job automatically submits a 100-period 2% annual property-tax continuation.
All receipts are equally rebated; PAYGO balances, the supply curve is fixed,
and there is no immigration or population rescaling. The shared 2% terminal
equilibrium is freshly verified before reuse. Both histories have the same
structural parameters and final fertility preference.

The submitted recovery has not yet run. A completed but unconverged baseline
can support only a provisional policy comparison; market convergence, terminal
distance and horizon robustness remain separate checks.

Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/fourshock_tax_20260913a`

Recovery outputs: `output/recovery/`. Automatic policy submission receipt:
`policy/dispatch.json`; reform outputs: `policy/tax2/` within the remote batch.
Pinned manifest and initial submission receipt are saved beside this note.

Sources: `code/cluster/recover_e5f_fourshock_tax_state.py` and
`code/cluster/run_e5f_inherited_2023_tax_long.py`. Local syntax checks and the
104-date evaluation-selection fixture passed. Native recovery verification is
pending the baseline dependency.
