# Torch staging preparation — not submitted

This pass writes only a compact local source archive and launcher. No SSH/scp/sbatch or model run occurred. `build_stage.py` verifies every runner/preparation source pin before archiving80 files (1,341,730 uncompressed bytes), including complete local refactor/credit packages, numerical proposal arrays and all required pinned configurations. It omits the original bundle/checkpoint and reuses their existing remote read-only locations.

Archive: `grid_resolution_stage.tar.gz`, SHA256 `fa243bde59afe3b49e2975fd52a061a74a5a1dcefc19bb967b09cf384b7dd483`, 330,695 bytes. `inventory.json` records every payload path/SHA; `mounts.txt` gives directory overlays and individual pinned-file mounts. The exact source tree is mounted over the immutable frozen-root bind. No frozen-root file is changed. The indexed package in this compact archive is byte-pinned locally; the older remote matched-pair directory is checked for existence but does not supply unverified replacement files.

Remote destination must be new: `/scratch/td2248/projects/grid_resolution_120x9_v1`. Existing inputs: `/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs`; frozen project: `/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project`. Apptainer maps the frozen project to `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`; read-only source overlays retain the identical repository-relative paths expected by the runner. All output/caches go through `/work/results`, backed by the new destination's results directory.

Lead staging sequence, after reviewing the archive and launcher (replace `TORCH_HOST` with the established login alias):

```bash
ssh TORCH_HOST 'mkdir /scratch/td2248/projects/grid_resolution_120x9_v1'
scp grid_resolution_stage.tar.gz TORCH_HOST:/scratch/td2248/projects/grid_resolution_120x9_v1/
ssh TORCH_HOST 'cd /scratch/td2248/projects/grid_resolution_120x9_v1 && echo "fa243bde59afe3b49e2975fd52a061a74a5a1dcefc19bb967b09cf384b7dd483  grid_resolution_stage.tar.gz" | sha256sum -c - && tar -xzf grid_resolution_stage.tar.gz && mkdir logs'
ssh TORCH_HOST 'sbatch /scratch/td2248/projects/grid_resolution_120x9_v1/launch_torch.sh'
```

Exact submission command: `sbatch /scratch/td2248/projects/grid_resolution_120x9_v1/launch_torch.sh`.

The launcher requests one CPU,24GB,40min, partition `cs`, account `torch_pr_570_general`. Its single2400s clock starts before module loading/source verification; a zero-solve preflight uses at most300s from that start, followed by full production only within remaining job time. It invokes `run_comparison.py` directly with the original deadline, because the standalone shell launcher deliberately creates a new clock. GNU `timeout --signal=KILL` bounds each container by its remaining external deadline; Slurm's40min cap independently applies to the entire job. Full production arms use separate initially empty Numba cache directories, fixedD=.14, six lifecycle calls each/twelve total. Input/source/gate failures stop without retry. The final job timing includes verification, preflight, compilation, authentication, solves, observers, plots and comparisons.

Local checks: shell syntax passes; archive extraction verified80 individual hashes and all mount sources; no checkpoint was bundled. The archive is deterministic (fixed tar/gzip metadata) and rebuilt only from unchanged pinned files. Container bind creation on Torch and full frozen-observer compatibility with120x9 remain untested. Lead must inspect the zero-solve job receipt before interpreting production results; no native numerical result is certified by these preparation checks.
