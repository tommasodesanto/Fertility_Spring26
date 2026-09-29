# Morning two-stream readout (September 29)

Experimental only; no adoption or promotion. Slurm `18766206_0` and `_1` completed (`0:0`). Reference: **2007 stationary reference — block0506, September 28 verified export**.

| Lane | selected case | attempts / successes | loss | repeats |
|---|---|---:|---:|---|
| Original one-birth | `one_birth_024_gn1_0` | 34 / 32 | 7.826 | 2/2 passed |
| Experimental two-birth | `two_birth_024_gn1_0` | 32 / 31 | 7.842 | 2/2 passed |

The two-birth experiment alone permits one optional extra attempt after first success, capped at two births per four-year cell; it uses the later-birth scale/inclusive value and independent conception draw, retaining common-event age projection. Original earnings, entry, transfers/floors, preference functional forms, targets, and bounds are retained subject to recalibration and `psi_child` normalization. Seed losses are not a controlled causal comparison.

## Full target fits and parameters

Full precision: [one-birth fit](run_v1/one_birth/one_birth_024_gn1_0/case/target_fit.csv), [two-birth fit](run_v1/two_birth/two_birth_024_gn1_0/case/target_fit.csv), [one-birth parameters](run_v1/one_birth/one_birth_024_gn1_0/case/parameters.csv), [two-birth parameters](run_v1/two_birth/two_birth_024_gn1_0/case/parameters.csv), [one diagnostics](run_v1/one_birth/one_birth_024_gn1_0/case/standard_diagnostics/), [two diagnostics](run_v1/two_birth/two_birth_024_gn1_0/case/standard_diagnostics/). Display values are three decimals.

### Normalization

| Moment | target | one model / gap | two model / gap |
|---|---:|---:|---:|
| Completed fertility | 2.100 | 2.100 / -0.000 | 2.100 / 0.000 |

### Scored moments

| Moment | target | one model/gap/weight/loss | two model/gap/weight/loss |
|---|---:|---|---|
| CPS childlessness | 0.198 | 0.198 / -0.000 / 35532.304 / 0.000 | 0.200 / 0.002 / 35532.304 / 0.152 |
| CPS one child among mothers, ages 40–44 | 0.214 | 0.214 / 0.000 / 26952.821 / 0.004 | 0.219 / 0.006 / 26952.821 / 0.875 |
| Mean first-birth age | 25.976 | 25.965 / -0.011 / 139.828 / 0.017 | 25.961 / -0.016 / 139.828 / 0.034 |
| Wealth/earnings | 6.927 | 6.896 / -0.031 / 7.595 / 0.007 | 6.830 / -0.097 / 7.595 / 0.071 |
| Bequest wealth | 0.007 | 0.007 / -1.049e-05 / 5165289.256 / 0.001 | 0.007 / -0.000 / 5165289.256 / 0.060 |
| Mean rooms | 5.729 | 5.730 / 0.001 / 128.021 / 3.471e-05 | 5.733 / 0.004 / 128.021 / 0.002 |
| Ownership, 30--55 | 0.676 | 0.676 / -0.000 / 2339.362 / 7.850e-06 | 0.673 / -0.003 / 2339.362 / 0.019 |
| First-birth rooms | 1.465 | 1.464 / -0.001 / 137.565 / 0.000 | 1.450 / -0.015 / 137.565 / 0.029 |
| Recent-parent ownership | 0.128 | 0.127 / -0.001 / 27055.823 / 0.009 | 0.118 / -0.010 / 27055.823 / 2.459 |
| Early fertility | 0.810 | 0.530 / -0.279 / 100.000 / 7.789 | 0.606 / -0.203 / 100.000 / 4.140 |

### Validation moments (untargeted)

| Literal identifier | target | one model/gap/weight/loss | two model/gap/weight/loss |
|---|---:|---|---|
| First births at age 30 or later | 0.249 | 0.226 / -0.024 / 0 / 0 | 0.230 / -0.020 / 0 / 0 |
| Old dispersion | 3.516 | 2.919 / -0.596 / 0 / 0 | 2.935 / -0.581 / 0 / 0 |
| Family rooms | 0.385 | 0.473 / 0.088 / 0 / 0 | 0.289 / -0.096 / 0 / 0 |

### Ten free coordinates and normalization

| parameter | one | two | shared bounds | one/two near bound |
|---|---:|---:|---|---|
| H0 | 6.104 | 6.108 | [0.200, 80.000] | no/no |
| beta_annual | 0.969 | 0.969 | [0.940, 0.990] | no/no |
| chi | 1.098 | 1.098 | [0.100, 5.000] | no/no |
| first_birth_fixed_cost | 0.353 | 0.364 | [0.000, 8.000] | no/no |
| kappa_fert | 0.109 | 0.094 | [0.020, 50.000] | yes/yes |
| kappa_fert_continuation | 0.222 | 0.405 | [0.020, 50.000] | yes/yes |
| theta0 | 0.156 | 0.135 | [0.000, 8.000] | no/no |
| delta_alpha_jump | 0.115 | 0.112 | [0.000, 0.250] | no/no |
| child_benefit_curvature | 0.066 | 0.057 | [0.000, 0.800] | no/no |
| tenure_choice_kappa | 0.013 | 0.011 | [0.001, 0.100] | no/no |
| psi_child (normalized) | 0.122 | 0.032 | completed-fertility normalization | n/a |

Near-bound flags use 1% of the broad inherited intervals; neither fertility scale equals its lower bound.

## Recorded local numerical diagnostics

Jacobians were evaluated at round centers before the selected Gauss--Newton step, not freshly at the selected candidate. Rank 9 uses relative cutoff `1e-6`; it is a local numerical warning, not proof of infeasibility, global identification, or statistical identification.

| lane | round | rank | condition | singular values |
|---|---:|---:|---:|---|
| one-birth | 0 | 10 | 11195.834 | 106.796, 44.832, 33.712, 32.912, 6.855, 4.130, 2.618, 1.937, 0.199, 0.010 |
| one-birth | 1 | 10 | 20147.537 | 93.390, 46.705, 33.675, 30.344, 6.583, 3.946, 2.722, 2.022, 0.197, 0.005 |
| two-birth | 0 | 10 | 7160.226 | 107.633, 39.946, 33.424, 27.359, 8.301, 3.982, 2.613, 1.893, 0.153, 0.015 |
| two-birth | 1 | 9 | 1.144e+06 | 97.975, 42.528, 33.704, 27.740, 6.593, 3.788, 2.524, 1.902, 0.102, 8.564e-05 |

Censored authenticated eight-solve-resource-cap cases: `one_birth_029_explore2`, `one_birth_031_explore4`, and `two_birth_029_explore2`; each reports `Old-steady-state fertility normalization missed tolerance: 23-solve cap`. No other failed records.

## Lead conclusion

Original early fertility contributes 99.52% of selected loss. After recalibration, two-birth closes about 27.1% of the original early-fertility gap, while overall loss is tied by losses elsewhere. No adoption.

## Provenance and collection limitation

Remote packet: `/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fertility_identification_20260928/two_stream_overnight_v1`; config SHA-256 `419b46d7cffdf51d95760f2999aec67a66b10b08bb9198164ac3d6f7d96323d3`.

Selected SUCCESS hashes matched FINAL; listed small artifacts passed hashes (31 one-birth, 33 two-birth, 17 PNGs/lane); selected case checkpoint identity, lane/source/target identity, 14 target rows, 31 parameter rows, ten free coordinates, identical bounds, and repeats passed. Large checkpoints were not newly hashed; identity is controller-recorded/audited only.

Collection mistake: initial recursive rsync copied 1,394 small files (113,908 KiB), with zero checkpoints. Before the instruction to stop deletion arrived, only unselected directories in this owned readout were removed. Retained readout: 140 files / 11,768 KiB. No further copies or deletions were made.
