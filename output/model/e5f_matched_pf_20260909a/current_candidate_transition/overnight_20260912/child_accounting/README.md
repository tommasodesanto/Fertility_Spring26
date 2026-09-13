# Dependent children and parental-household death

Read-only diagnostic for M13. This is the initial stationary checkpoint used for the recovered no-rebate historical sequence, **not its 2023 distribution** and not the later selected initial calibration. Checkpoint SHA-256: `120ffc45c0fb8756f4182f999c96b7c0236adf315cb938190ec31cd2068c87c2`. Full source path and unrounded results are in `initial_parent_death.json`.

Per four-year transition and one unit of household mass:

- Current dependents: 0.497832619.
- Dependents attached to exiting parent households: 0.005972242, or **1.1996% of current dependents**.
- Births: 0.115274548; dependent loss equals **5.1809% of births**.
- Survivor-household maturation: 0.109302306; loss equals 5.4640% of that flow.
- 48.985% of dependent loss occurs at the forced final household exit, from age cell 82–85.
- Parent survival is one through the 62–65 cell. All affected children are attached to parents aged 66 or above.

These are transition **flows**, not a measured stock of unassigned children. The model has no such stock. Literal dependent counts are used consistently: the 3+ state counts as three here, without demographic top-bin expansion. The flow does not measure a general-equilibrium or welfare effect of a repair.

For post-birth dependent stock C_a and parent survival s_a, loss is D_a=(1-s_a)C_a, with final-age survival set to zero exactly as in the forward operator. Maturation is computed from the saved child transition, conditional on parent survival. The identity C_a=D_a+M_a+remaining_a holds to 0; post-choice and post-birth child stocks agree within 8.74e-15. In this stationary state, births minus maturation minus parent-exit dependent loss is -9.01e-11.

| Parent age cell | Current dependents | Household exit probability (%) | Dependent loss |
|---|---:|---:|---:|
| 18–21 | 0.015804885 | 0.000 | 0.000000000 |
| 22–25 | 0.030416503 | 0.000 | 0.000000000 |
| 26–29 | 0.043226471 | 0.000 | 0.000000000 |
| 30–33 | 0.052939999 | 0.000 | 0.000000000 |
| 34–37 | 0.059061892 | 0.000 | 0.000000000 |
| 38–41 | 0.060817183 | 0.000 | 0.000000000 |
| 42–45 | 0.056993008 | 0.000 | 0.000000000 |
| 46–49 | 0.044327895 | 0.000 | 0.000000000 |
| 50–53 | 0.034477252 | 0.000 | 0.000000000 |
| 54–57 | 0.026815640 | 0.000 | 0.000000000 |
| 58–61 | 0.020856609 | 0.000 | 0.000000000 |
| 62–65 | 0.016221807 | 0.000 | 0.000000000 |
| 66–69 | 0.012616961 | 6.087 | 0.000768041 |
| 70–73 | 0.009215827 | 8.150 | 0.000751112 |
| 74–77 | 0.006583667 | 11.505 | 0.000757436 |
| 78–81 | 0.004531513 | 16.995 | 0.000770145 |
| 82–85 | 0.002925508 | 100.000 | 0.002925508 |

## Candidate repair, not implemented

The measured problem comes from the unbounded tail of memoryless dependency: the current departure probability is 0.222222222 per four-year period (18-year mean), while the last possible birth is at model age 42. A child-aging scheme with a minimum duration and a maximum of 20 years would make every child independent by parent age 62, before parental mortality begins in this saved specification. As a concrete illustration, equally likely departure after 16 or 20 years preserves the current unconditional 18-year mean duration. This would require child-age/cohort information and a new solution/fit assessment; matching the mean does not preserve incentives or housing demand automatically. It also eliminates first-period maturation. This is a candidate to discuss, not an approved specification.

If mortality at younger parent ages is introduced, bounded child duration alone will not suffice: surviving children need an explicit guardian/care assignment that preserves their consumption and housing needs. Immediate adulthood on parental death is not the proposed repair.

## Reproduction

`code/model/tools/measure_e5f_dependent_parent_death.py CHECKPOINT.pkl.gz` reads and aggregates the saved checkpoint without importing the solver or solving any model. Run with NumPy available, redirecting stdout to JSON. The source checkpoint is on `ssh torch`; the script can be passed on stdin to `/share/apps/anaconda3/2025.06/bin/python - CHECKPOINT`. The source-reduction formulas received an independent read-only review.
