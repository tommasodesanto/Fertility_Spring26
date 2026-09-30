# Torch collection: fixed-price pair/grid-fix job 18850069

Collected compact evidence from `/scratch/td2248/projects/publication_refactor_20260929/results/fp_pair_gridfix_18850069/`; no model work was run during collection. Large NPZ/pickle arrays and Numba caches remain on Torch.

Completed receipt: `rc=0`, phase `all_passed`, elapsed 529 s; maximum four lifecycle solves. Both scalar and indexed branches report all five steps exit 0.

## Tests and numerical checks

The compiled component check reports 2 passed and 25 deselected in 9.14 s. In each scalar/indexed × rep1/rep2 frozen-oracle receipt: 113/113 nested numeric arrays are exact and finite; 14 fit rows and 31 parameter rows are exact; all 17 standard plot hashes were checked. The one collected standard plot packet is indexed oracle rep1; its 17 files (1.7 MB) match the receipt hashes.

The lab-to-checkpoint comparisons report 67 common arrays exact, zero reference arrays missing, and 20 lab-only shared arrays; these receipts predate the `strict_paths` field. The rep2-to-rep1 comparison reports 87 arrays exact with no missing, extra, or differing keys.

## Time and resources

Scalar phase elapsed 292 s; lab lifecycle timings are 86.07 s and 72.24 s. Indexed phase elapsed 223 s; lab lifecycle timings are 56.60 s and 46.23 s. Inputs, precompute, and serialization timings remain in each lab receipt.

The saved plan sets one thread, 12 GiB, a 2,400 s total cap, a 300 s component cap, and 900 s per fixed-price phase. Slurm accounting:

```text
18850069|refactor_fp_pair|COMPLETED|0:0|00:08:56|1|12G||
18850069.batch|batch|COMPLETED|0:0|00:08:56|1||6796984K|0
18850069.extern|extern|COMPLETED|0:0|00:08:57|1|||
```

## Collected files

Both branches include plan/summary/steps/cache state, lab and oracle summaries, lab receipt, all three comparison JSONs, both oracle repetition receipts/control arrays/target-fit and parameter tables, and oracle diagnostic summaries. Root completion, plan, component XML/log, driver logs, and `slurm_resources.tsv` are also included. No requested compact file was missing.
