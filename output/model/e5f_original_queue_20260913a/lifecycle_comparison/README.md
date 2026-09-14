# Lifecycle profiles: 2007 and 2023

Author-requested comparison built September 14, 2026, with no Bellman,
equilibrium, calibration or transition solve.

- `lifecycle_2007.pdf` / `.png`: the initial stationary model versus ACS 2007.
- `lifecycle_2007_2023.pdf` / `.png`: the same panels and axes in both years.
- `lifecycle_comparison.csv`: all plotted age cells, model/data values and gaps.
- `verification.json`: input hashes, plotted arrays, and exact reproduction of
  the 2023 lifecycle figure's verified arrays.

The 2007 model is extracted from the same frozen initial checkpoint used by
the original-household-queue exercises. The existing 2023 age observer is
reused, with exact checks against the initial prices/distribution and native
aggregate housing demand. `model_2007.json` records the frozen specification
hash and observer hash. Extraction sources are mirrored in the Torch batch
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/lifecycle_comparison_20260914a/source/`.

The data reuse `code/model/tools/extract_e5f_patch_validation.py`, setting its
source-year selector to 2007. The household sample remains the same 42 matched
metros, GQ 1/2, PERNUM=1, RELATE=1, positive HHWT, ages18–85, OWNERSHP1/2 and
positive rooms. Ownership is OWNERSHP1, rooms are capped at9, and resident
minor children use NCHILD>0 and YNGCH<18. The extraction reads only the 2007
row range and finds 516,656 eligible household heads; no room codes>=99 occur.
The household metadata correct the legacy builder's literal2023 labels;
the extraction filters themselves already use the selected year.

The 2023 model/data inputs are the existing, verified one-permanent-shock
lifecycle cross-section used in the presentation. This figure does not use
the announced four-shock path. The 2023 transition is not converged. The
initial2007 economy is the calibration benchmark, but the full age curves
are descriptive comparisons, not separately targeted calibration moments.
In particular, model dependents and empirical resident own children under18
are distinct measurement objects. These figures make that mismatch visible.

Regenerate the figures locally from saved inputs:

```sh
python code/model/tools/build_e5f_lifecycle_comparison.py --mode render
```

Re-extract the 2007 empirical profiles locally with `--mode data`. To extract
the model on Torch, run the same script with `--mode model`, `--spec` pointing
to the frozen `afternoon_original_queue_20260913a/spec.json`, `--profile-helper`
pointing to the mirrored `recover_e5f_permanent_2023_profile.py`, and `--outdir`
pointing to the diagnostic batch's `output` directory. Copy `model_2007.json`
back into this folder before rendering. No numerical model run is needed.
