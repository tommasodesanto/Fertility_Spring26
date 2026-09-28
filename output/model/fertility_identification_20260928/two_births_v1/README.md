# Two births in one four-year period: fixed-policy diagnostic

Author-authorized diagnostic against **2007 stationary reference — block0506,
September 28 verified export**. No specification is adopted by this experiment.

This first test replays the saved household choices, allowing one additional
birth opportunity immediately after a successful birth. The extra opportunity
uses the reference's subsequent-birth attempt probability at the new child
state and the same age-specific conception probability. A household can have
at most two new births in the period; a failed attempt or a decision to wait
does not open another opportunity. This is an accounting diagnostic with fixed
policy functions, not a solved two-birth preference specification or a model of
twins. It adds no Gumbel draw or inclusive value to household optimization.

**Only experimental economic change:** the extra within-period birth
opportunity. Earnings, entry, housing/saving/location choices, prices,
preferences, targets and parameter bounds remain at the reference. The saved
child benefit is fixed. Consequently completed fertility can change: this
diagnostic does not enforce replacement or claim to be an equilibrium.
Any later calibrated experiment must restore the author's 2.1 normalization
and demographic renewal check. No transition or search is launched.

Budget: zero Bellman/market/normalization solves; two full cohort replays,
one control and one diagnostic, on one Torch CPU with 24 GB and a 15-minute
Slurm cap (14-minute internal cap). Existing fixed-policy propagation is a
small part of the approximately 971 seconds used by the reference's six-solve
normalization; the cap is a limit, not a promised runtime. No retries or search.
Progress prints each age and saves after each completed replay.

Validation before interpreting the diagnostic: source/checkpoint/30 artifact
hash authentication; toy mass and birth accounting including no more than two
new births; control replay against every saved age's full pre/post distribution;
reference birth flow, age-25 count and first-birth mean reproduction; unchanged
first-birth cell flows; nonnegative mass, cohort mass and dead-state checks.
The saved DUE stayer savings policy is used explicitly.

Age-window counts keep the existing linear interpolation of pre/post stocks.
This does not specify the ordered dates of two births within a four-year cell;
age-cell endpoints and birth flows are retained separately in `result.json`.
The original 17 plots and full fit/parameter tables remain in
`../resume_v1/selected_export/primary/`; this replay is not scored as a new
calibration and does not replace that standard packet.

Run on Torch using `sbatch run.sh` in a fresh versioned output location. The
completed version-one artifacts below must not be overwritten.

## Verified result

Torch **18744518** completed successfully in 45 seconds (40 seconds within the
Python driver), with zero model solves. The control reproduces every saved
pre/post distribution to maximum absolute differences of 3.096e-17 / 4.224e-17;
all three birth-order flows, age-25 count and first-birth mean reproduce.
Maximum survival-adjusted cohort mass error is 1.774e-14. Maximum occupied
dead-state exposure is 9.368e-15, below the unchanged 1e-12 gate. Original
source ancestry, checkpoint and 30 export files including all 17 plots pass
authentication. A separate read-only reviewer checked routing, array axes and
native renewal accounting before launch; all numerical tests ran on Torch.

**Age 25 — the existing targeted count, plus its diagnostic decomposition**

| Object | Data | Reference | Two-birth replay |
|---|---:|---:|---:|
| Children per woman, capped at 3 (targeted) | 0.810 | 0.535 | 0.698 |
| Mothers, percent (diagnostic) | 45.725 | 45.011 | 45.011 |
| Children among mothers (diagnostic) | 1.770 | 1.190 | 1.551 |

The increase is 0.163 children per woman, closing 59.409% of the reference's
0.274 age-25 shortfall. This is a direct mechanical effect in the stipulated
fixed-policy replay, not an estimate of how much survives reoptimization and
replacement normalization. Mapped first-birth age remains 25.933 against the
25.976 target; first-birth cell flows and motherhood are unchanged by design.

**Untargeted lifecycle windows — children ever born capped at 3**

| Age window | Data | Reference | Two-birth replay |
|---|---:|---:|---:|
| 20–24 | 0.475 | 0.308 | 0.410 |
| 25–29 | 0.990 | 0.697 | 0.890 |
| 30–34 | 1.476 | 1.096 | 1.336 |
| 35–39 | 1.709 | 1.453 | 1.692 |
| 40–44 | 1.718 | 1.731 | 1.939 |

The replay improves earlier counts but overshoots at 40–44. It both advances
births and raises lifetime counts. Terminal capped fertility rises from 1.870
to 2.053; using the saved top-bin weight, the separate normalization object
rises from **2.100 to 2.362**. The implied potential household entrants exceed
retained entry by **12.494%**. No renewal pass is claimed for this experiment.
The reference itself reproduces its renewal restriction.

This establishes that the one-birth-per-period restriction can materially
depress early fertility at the reference's saved choices. It does not establish
that a two-birth model can jointly match early fertility, first-birth timing,
all other targets and replacement fertility. That next question requires a
fully specified household choice problem, coherent birth/housing observers,
and a normalized equilibrium solve; it is not answered by this replay.

Artifacts: [full results and identities](result.json),
[all model lifecycle windows](lifecycle.csv), [Torch log](replay_18744518.log),
and [reproduction source](replay.py). Empirical lifecycle data are the unchanged
`../measurement_audit_v1/fertility_lifecycle_matched_windows.csv`.
The complete reference target/parameter tables, with targeted and untargeted
rows separated, are in [the measurement audit](../measurement_audit_v1/README.md#complete-reference-fit-and-parameter-restrictions).
The replay keeps every parameter at that estimate and does not estimate new
parameters, change bounds, score a loss or produce a replacement solution.
