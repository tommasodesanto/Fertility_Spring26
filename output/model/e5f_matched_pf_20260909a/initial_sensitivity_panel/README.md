# Initial sensitivity panel

Source 7e872053; all Python sources and inherited seed pinned. Diagnostic only, no objective, weights, empirical target activation or SMM promotion. All nine structural coordinates are bound explicitly. Each case separately normalizes fertility to 2.1 with new parenthood utility, balanced initial pension, fixed payroll tax .179, exhaustive saving and common supply elasticity .63.

Design: baseline plus two directions for each of nine coordinates, 19 cases total. Eight coordinates move ±2% in levels. Beta moves by ±2% of its annual log discount rate, preserving admissible bounds and a local discount perturbation. Cases record actual level differences for derivatives. No bound clips.

First run two fresh exact loops at baseline, including early measurement and the stable17 graph packet. Require success before launching19 cases. Each case: at most8 GE solves,30 minutes internal,32 minutes Slurm, one CPU16GB and one numerical thread. Expected approximately5 minutes per perturbed case, based on4 solves at65–75seconds; 19 concurrent jobs, well below96 ceiling. Save every completed equilibrium, each final checkpoint, parameter table, observer outputs and standard diagnostics. Failed cases remain failures; no loosened gates.

CPS within-cell age projection is saved in two versions. Family-size housing uses an explicit dependent-count proxy; exact recent-parent ownership remains unavailable. The panel diagnoses moment sensitivity without pretending these measurement decisions or weights are certified.

## Prepared diagnostic analysis

The analysis script refuses to construct a full Jacobian until all19 planned directions pass. It uses log coordinates; beta uses the log annual discount rate. It reports centered derivatives, one-sided discrepancies, singular values and numerical rank. Numerical rank does not certify empirical identification.

Before looking at the completed panel, two diagnostic scaling choices are declared: synthetic5% target scales for all11 observed scored restrictions, and a reference profile using the audited CPS correlation-one approximation, NCHS annual temporal SD (not a sampling SE), other listed standard errors and the inherited synthetic bequest scale. Actual SMM weights remain unset. The twowithin-cell CPS projections are examined separately and jointly. The unavailable recent-parent ownership row is retained as missing and is not replaced.

Bounded ridge least squares proposes joint parameter directions with three penalty strengths and three step fractions. These are linear predictions, not calibration results. Nonlinear verification, full target/parameter tables, robust comparisons across measurement/scaling choices, and a separately documented launch budget are required before expanding the run. No proposals are generated from the incomplete panel.
