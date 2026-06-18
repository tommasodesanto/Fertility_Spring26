# Intergen Candidate Replacement V1 Panel

Date: 2026-06-18

This note records a local diagnostic panel for the one-market
intergenerational model using the candidate replacement target set
`candidate_replacement_v1`. It is a frontier check, not a production SMM
calibration.

## Command

```bash
PYTHONPATH=code/model code/model/.venv/bin/python -m intergen_housing_fertility.cli local-panel \
  --cases 72 --seed 20260618 --J 16 --Nb 60 --income-states 5 --n-house 5 \
  --max-iter-eq 3 --workers 6 --minutes 20 --diagnostic-best 2 \
  --target-set candidate_replacement_v1 \
  --outdir output/model/intergen_candidate_replacement_v1_local_panel_20260618
```

The run completed all 72 submitted cases in 562.6 seconds. The best finite
record is `draw_0055` with rank loss `34.009`, old nonlocation loss `994.965`,
and market residual `3.24e-05`. The baseline under this target system has rank
loss `48.883`, so the diagnostic search improves the scalar objective but does
not yet deliver a substantively good fit.

## Target Set

`candidate_replacement_v1` keeps the core fertility, tenure, ownership,
housing-response, and young-liquid-wealth moments, and adds replacement
moments for:

- old nonhousing wealth to income, ages 65--75;
- parent-minus-childless old nonhousing wealth to income, ages 65--75;
- childless renter and owner mean rooms, ages 30--55;
- childless renter and owner shares with at least six rooms, ages 30--55.

The target count remains 13 for 13 varied internal parameters. This preserves
formal count identification, but the previous Jacobian audit showed that count
identification is not enough: the chosen moments must also be informative and
economically mapped to the intended parameter blocks.

## Best Case

Best case parameter vector:

| Parameter | Value |
|---|---:|
| `beta` | 0.9361 |
| `alpha_cons` | 0.6977 |
| `b_entry_fixed` | -0.2718 |
| `c_bar_0` | 0.3405 |
| `c_bar_n` | 0.4384 |
| `h_bar_0` | 3.2761 |
| `h_bar_jump` | 1.2460 |
| `h_bar_n` | 0.9286 |
| `psi_child` | 0.0528 |
| `kappa_fert` | 6.8353 |
| `chi` | 1.2382 |
| `theta0` | 0.8548 |
| `theta_n` | 0.1949 |

Moment fit:

| Moment | Target | Model | Gap |
|---|---:|---:|---:|
| `tfr` | 1.7000 | 1.9516 | +0.2516 |
| `childless_rate` | 0.1500 | 0.2303 | +0.0803 |
| `own_rate` | 0.5755 | 0.7501 | +0.1747 |
| `own_family_gap` | 0.1677 | 0.1416 | -0.0261 |
| `housing_increment_0to1` | 0.6644 | 1.2260 | +0.5616 |
| `housing_increment_1to2` | 0.4880 | 0.4153 | -0.0727 |
| `young_liquid_wealth_to_income` | 0.1792 | 0.1888 | +0.0096 |
| `old_nonhousing_wealth_to_income_6575` | 6.4185 | 2.2138 | -4.2047 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.0074 | 0.2683 | -0.7391 |
| `prime30_55_childless_renter_mean_rooms` | 3.8053 | 4.6205 | +0.8152 |
| `prime30_55_childless_owner_mean_rooms` | 6.2240 | 5.2725 | -0.9515 |
| `prime30_55_childless_renter_share_rooms_ge6` | 0.1377 | 0.0085 | -0.1291 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.5961 | 0.5960 | -0.0001 |

## Frontier Read

The run reveals a sharper tradeoff than the previous target bundles. The best
ranked candidate matches young liquid wealth and the childless owner large-home
share, but it overstates ownership, overstates first-birth housing growth,
keeps old nonhousing wealth far below the PSID mean target, and compresses the
childless owner-renter room separation.

High-old-wealth candidates exist, but they move into a different basin. The top
old-nonhousing-wealth records reach `5.77`, `5.24`, and `5.12` wealth-to-income,
but their ownership rates are about `0.140`, near zero, and `0.077`,
respectively, with very large renter rooms. Conversely, records with very large
owner mean rooms generally have near-zero ownership and renter rooms that are
also too large.

The current mechanism therefore appears to struggle with the joint pattern:

1. moderate prime-age ownership;
2. high old nonhousing wealth;
3. childless owners in substantially larger homes than childless renters;
4. not putting childless renters into too much space.

## Interpretation

This is a partial "yes" only in the narrow sense that the replacement target
system is runnable and the panel improves rank loss from `48.883` to `34.009`.
It is not a good enough match to treat as a calibrated benchmark.

The leading issue is not that we should drop the old-age wealth moments. Under
the SMM identification rule, any demotion or replacement must identify the same
parameter block or fix the affected parameter externally. Before changing the
objective, re-audit the PSID old-wealth object:

- mean versus median nonhousing wealth to income;
- nonhousing versus total wealth;
- denominator income definition;
- age window and household head definition;
- treatment of bequests, home equity, and business wealth;
- whether the model statistic should be an unconditional old-age ratio or a
  ratio by parental state.

The old-age mean target `6.419` is much harder than the median target `2.230`;
which object is appropriate materially changes the calibration problem.

## Next Diagnostic

The next useful run should not be a blind longer search on this same objective.
Run a targeted two-basin diagnostic instead:

1. around `draw_0055`, to see how much ownership and room separation can be
   improved without losing young liquid wealth;
2. around high-old-wealth candidates such as `draw_0064` and `draw_0070`, to
   identify exactly why ownership collapses as old nonhousing wealth rises;
3. with side-by-side ranking under old nonhousing mean, old nonhousing median,
   and old total wealth targets, after the data-object audit.

The replacement target system remains a diagnostic candidate until that audit is
complete.
