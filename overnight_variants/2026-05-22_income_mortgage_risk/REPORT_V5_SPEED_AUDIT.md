# V5 HANK-z Speed Audit

This audit runs the benchmark-normalized outside-option V5 closure for a fixed
small number of GE iterations. These are timing probes, not accepted
calibrations. The goal is to separate Bellman time, fast forward-distribution
time, and residual final-statistics / overhead time as \(N_z\) and \(N_b\)
change.

## Fixed-Iteration Audit

| \(N_b\) | \(N_z\) | GE iters | elapsed sec | Bellman sec | forward sec | residual sec | Bellman / iter | forward / iter |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 3 | 4 | 22.84 | 4.85 | 3.01 | 14.98 | 1.21 | 0.75 |
| 30 | 5 | 4 | 73.32 | 10.91 | 5.71 | 56.70 | 2.73 | 1.43 |
| 30 | 7 | 4 | 142.36 | 23.06 | 10.71 | 108.59 | 5.76 | 2.68 |
| 20 | 7 | 4 | 108.57 | 14.97 | 6.68 | 86.92 | 3.74 | 1.67 |
| 30 | 7 | 8 | 164.44 | 38.95 | 18.06 | 107.43 | 4.87 | 2.26 |

## Accepted Full V5 References

Cold standalone accepted run:

- elapsed seconds: `620.13`
- iterations completed: `17`
- Bellman seconds: `257.46`
- forward-distribution seconds: `103.82`
- residual final-stat / overhead seconds: `258.86`
- Bellman seconds per GE iteration: `15.14`
- forward seconds per GE iteration: `6.11`

Same-process warm accepted run after a `30x7x4` warm-up:

- elapsed seconds: `200.82`
- iterations completed: `17`
- convergence reason: `strict_tol`
- Bellman seconds: `79.81`
- forward-distribution seconds: `36.38`
- residual final-stat / overhead seconds: `84.63`
- Bellman seconds per GE iteration: `4.69`
- forward seconds per GE iteration: `2.14`

Standalone `30x5` accepted plotting run:

- elapsed seconds for solve plus figure packet: `83.10`
- iterations completed: `12`
- convergence reason: `strict_tol`
- Bellman seconds: `28.69`
- forward-distribution seconds: `14.25`
- residual final-stat / plotting / overhead seconds: `40.17`
- Bellman seconds per GE iteration: `2.39`
- forward seconds per GE iteration: `1.19`

## Read

The runtime problem is not the benchmark normalization arithmetic itself.
There are two separate problems:

1. Cold-start / compilation overhead is large. The accepted `30x7` solve falls
   from `620.13` seconds cold to `200.82` seconds when the same process has
   already warmed the relevant HANK-\(z\) kernels.
2. The recurring warm solve is still too slow for calibration. At `200.82`
   seconds per accepted solve, `500` evaluations would take about `27.9` hours
   on one process before failed candidates, diagnostics, or search overhead.

The repeated-cost bottleneck is split across the structural HANK-\(z\) Bellman,
the HANK-\(z\) forward distribution, and the final full-statistics pass. The
fixed-iteration rows show the state-space scaling clearly: at \(N_b=30\), moving
from \(N_z=3\) to \(N_z=7\) raises Bellman time per iteration from `1.21` to
`5.76` seconds and forward time per iteration from `0.75` to `2.68` seconds.
The final-statistics residual also rises from about `15` seconds to about
`108` seconds.

## Implications

- Use \(N_z=5\), not \(N_z=7\), for any exploratory calibration unless the
  point of the run is explicitly to validate the final \(N_z=7\) economics.
- Add a cheap calibration/reporting split. Calibration should avoid the full
  Python final-statistics pass where possible; accepted/reporting runs can keep
  the full moment table.
- Add a HANK-\(z\) Howard/evaluation path analogous to the baseline solver. The
  current HANK-\(z\) GE loop resolves the full Bellman every price iteration.
- Warm each calibration worker before timing or comparing solve costs. The
  cold first solve is not representative of a worker that evaluates many
  candidates, but it is representative of one-off local experiments.
