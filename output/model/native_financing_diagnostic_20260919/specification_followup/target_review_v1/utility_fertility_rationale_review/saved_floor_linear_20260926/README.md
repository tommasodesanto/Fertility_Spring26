# Saved floor/linear fertility incentives — September 26

No new equilibrium, household solve, recalibration or economic change. This
diagnostic reads the verified selected floor/linear point from the September 26
comparison. The original search stopped in its initial population, so this is
an experimental selected point, not an optimized or adopted final baseline.
The original ACS mean-rooms target and unmatched PSID/estate observers remain.

The recovered gap compares trying for a birth now with waiting, before the
current taste draw, but includes future shock-inclusive continuation values.
It does not compare parenthood with remaining childless forever. Negative
utility levels, negative current flow differences and negative action-value
gaps are different objects. No zero-shock equilibrium is calculated.

## Findings

| Share of expected births | First birth | Later births |
|---|---:|---:|
| Positive finite attempt-minus-wait gap | 54.654% | 29.763% |
| Negative finite attempt-minus-wait gap | 45.346% | 70.237% |
| Probability boundary | 0% | 0% |

The interior at-risk distribution's median normalized gap is -1.031 for first
birth and -0.224 for later births. These are decision-population weights, not
birth weights. At-risk states with zero attempt probability account for 4.872%
of the first-birth pool and 0.005% of the later-birth pool, and are excluded from
finite-gap quantiles. Unavailable choice mass is separately disclosed below.

The universal negative-gap claim is false. Favorable current taste draws
nevertheless overturn a systematic preference to wait for many births,
especially later births. That fact alone does not establish that the model is
misspecified or that economic incentives have little influence. Later-birth
gaps tend to become positive near the end of the fertility window; the age
table and figure make this timing pattern visible.

## Method and verification

The exact comparison collector authenticates contract, source, checkpoint and
all original target and parameter rows before loading. Both stored choice
probabilities are used: Delta/kappa = log(p_try) - log(p_wait), avoiding loss of
precision from subtracting a near-one probability. Only positive probability
pairs are inverted. At zero endpoints, the audit does not distinguish numeric
underflow from an infeasible alternative and assigns no finite gap.

Weights are `stationary_g_pre`, immediately before fertility choices. Each
age/family risk pool is snapshotted before births. Expected births use the
attempt probability times the verified conception probability. A manual
application of these flows matches the saved post-fertility distribution with
L1 error 3.204e-15 and the unadjusted birth total with error 4.163e-16. Later
births include entry into the 3+ group; this statistic does not apply its
separate 3.602 demographic weight.

Torch job 18602545 stopped because the native solver sets both probabilities
to zero at dead states. The revised job 18602589 explicitly classified that
mass (4.943e-15 across risk pools), retained the original 1e-12 gate, and passed
in 14 seconds. Both jobs performed zero model solves. The source and checkpoint
hashes and actual configuration are in `receipt.json`. Figure visually checked.

## Files and reproduction

- `summary.csv`: all risk and birth shares, censoring and interior quantiles.
- `by_age_birth_number.csv`: the same objects by age and birth number.
- `incentives.png`: supplemental age profile, not a replacement diagnostic set.
- `target_fit.csv` and `parameters.csv`: complete original 13-row fit and
  28-row parameter/restriction tables, unchanged from the selected point.
- `receipt.json`: source identity, actual settings, checks and interpretation.

On Torch, load `anaconda3/2025.06`, use one thread, and set `PYTHONPATH` to the
frozen comparison's `tools` directory and `/scratch/td2248/commute_pdf_qa_deps`.
Then execute the repository helper with a fresh output directory:

```bash
export EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256=c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a
python3 code/model/tools/diagnose_saved_utility_birth_incentives.py \
  --contract /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/launch_v1/contract.json \
  --run-root /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/results/run_001 \
  --output /path/to/fresh/output
```

Retained remote run directory:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/saved_utility_incentives_20260926/`.
