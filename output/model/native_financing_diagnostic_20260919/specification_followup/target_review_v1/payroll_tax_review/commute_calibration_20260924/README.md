# Experimental commute calibration — September 24

The detached Torch chain was submitted as smoke `18471839`, dependent eight-worker
array `18471840`, and dependent selected-export job `18471841`. The remote run
root is `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/commute_calibration_20260924_v1/results/run_001`.
The shared deadline, measured from first smoke start and covering search and
export, is **2026-09-24 20:56:08 EDT**. Smoke must complete its full loop and
17 standard figures before the worker array can run. Production allows at most
five proposals per worker, so the entire chain has at most 41 objectives
including smoke. Latest completed and best-so-far summaries are written under
the remote run root. This launch record does not report a completed fit.

The objective contract is `objective.json` (SHA-256
`c28e5d620d463dc3a592c038ca2771b47bf37ed1a9cf5ba31f7dc793b58b61db`):
eight free structural coordinates, twelve positive-weight moments with every
previous working weight retained, and a separate completed-fertility
normalization at $2.1$. The selected B-floor checkpoint SHA-256 is
`83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee`.
The frozen native source inventory SHA-256 is
`76406fcc10206d9e30bcc29d4e18accdf7fffdd9219e50e6e11c1c3360f01336`.
The manifest pins the new driver, Slurm scripts, objective, selected plan,
receipts, and reviewed PAYGO comparison driver. `preflight_receipt.json` records
the Torch zero-solve checks. Source and numerical gates are otherwise unchanged.

Relative to the selected September 23 configuration, the experiment applies
the author-accepted national housing targets, PSID wealth/earnings ratio
$6.92658379107299$, SCF bequest-flow/wealth ratio $0.007291023472616158$,
first-birth rooms response $1.465$, externally fixed $\theta_1=
0.008193084126995582$, annual depreciation $0.01416143718381309$, and annual
property tax $0.010598360773872594$. The PAYGO payroll tax
$0.08751017424959717$ is **experimental, not adopted**. The floor utility,
B15 income process, inherited heterogeneous entry wealth, supply elasticity
$0.63$, other unchanged targets and weights, and stationary entry law remain.
The adopted 16/20 birth-entry queue has not been implemented, so this is not a
new demographic steady state. The model's stationary housing observer is not
the PSID panel estimator; its bequest observer counts positive estates across
all child states while the accepted SCF target is child-directed. These two
measurement mismatches remain explicit approximations; no model observer or
recipient rule was changed. Do not compare this objective's loss with old
target systems.
