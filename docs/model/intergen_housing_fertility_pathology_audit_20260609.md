# One-Market Final-Best Pathology Audit

Date: 2026-06-09

## Scope

This note documents fixed-theta local diagnostics for the final saved
global-DE best from the Torch shutdown run:

- run directory:
  `results_intergen_housing_fertility_intergen_candidate_no_timing_v0_globalde_3g_20260609`
- task `27`, case `236`, label `de_g008_i011`
- rank loss under `candidate_no_timing_v0`: `11.503191936648555`

The audit does not recalibrate. It asks whether the visibly broken policies and
wealth distributions are caused by the asset grid, the temporary owner menu, or
the economics/objective itself.

The reproducible audit driver is:

```bash
code/model/.venv/bin/python code/model/tools/audit_intergen_final_best_pathologies.py
```

Primary output:

```text
output/model/intergen_final_best_pathology_audit_20260609/
```

## Tests

The same final-best parameter vector was re-solved under four variants:

| Variant | Change relative to final saved run |
|---|---|
| `baseline_nb60_cluster_ladder` | Exact final saved settings: `Nb=60`, `H_own=linspace(2,10,6)` |
| `fine_nb120_cluster_ladder` | Same owner ladder, but `Nb=120` |
| `baseline_nb60_room_ladder` | Same `Nb=60`, deliberate owner ladder `[2,4,6,8,9.5,11]` |
| `fine_nb120_room_ladder` | `Nb=120` and deliberate owner ladder `[2,4,6,8,9.5,11]` |

All four variants clear the one-market residual tightly, so the pathology is
not a market-clearing failure.

## Results

| Variant | Loss | Prime own | Old own | Housing share | TFR | Childless | Owner median rooms |
|---|---:|---:|---:|---:|---:|---:|---:|
| `baseline_nb60_cluster_ladder` | `11.503` | `0.222` | `0.783` | `0.327` | `1.660` | `0.238` | `5.2` |
| `fine_nb120_cluster_ladder` | `11.942` | `0.184` | `0.758` | `0.326` | `1.700` | `0.226` | `5.2` |
| `baseline_nb60_room_ladder` | `45.546` | `0.225` | `0.851` | `0.327` | `1.689` | `0.234` | `4.0` |
| `fine_nb120_room_ladder` | `45.947` | `0.174` | `0.747` | `0.324` | `1.711` | `0.225` | `4.0` |

### Asset Grid

Moving from `Nb=60` to `Nb=120` with the same owner ladder does not rescue the
economic fit. Prime-age ownership falls from `0.222` to `0.184`, old-age
ownership remains near target, and the old-minus-prime ownership gap remains
about `0.57`.

The finer grid does smooth some age-specific wealth distributions. For the
cluster ladder:

| Age | Metric | `Nb=60` | `Nb=120` |
|---:|---|---:|---:|
| `30` | max asset-bin share | `0.257` | `0.311` |
| `30` | effective asset bins | `6.28` | `6.67` |
| `42` | max asset-bin share | `0.112` | `0.074` |
| `42` | effective asset bins | `14.95` | `28.62` |
| `46` | max asset-bin share | `0.103` | `0.059` |
| `46` | effective asset bins | `16.19` | `32.66` |

The `max asset-bin share = 1.0` reported in the automated summary is partly
mechanical: entrants start with deterministic liquid wealth at age 22. The
more informative reading is that midlife wealth distributions improve with
`Nb=120`, but the early-life distribution remains concentrated and the
lifecycle ownership path remains wrong.

### Owner Menu

Replacing the temporary cluster ladder

\[
H^{own}=(2,3.6,5.2,6.8,8.4,10)
\]

with the deliberate room ladder

\[
H^{own}=(2,4,6,8,9.5,11)
\]

raises the loss from about `11.5` to about `45.5` at the same theta. The
largest visible reason is that prime childless owner median rooms falls from
`5.2` to `4.0`, against the target `6.0`.

Interpretation: the final best is not robust to the owner menu. It exploited
the temporary `5.2` rung to partially fit the owner-room moment. A production
calibration should not continue with `linspace(2,10,6)` if owner-room levels
are targeted.

### Owner-Entry Policy

The childless-renter owner-entry policy is not merely kinked. Under the
baseline variant, `93.8%` of checked age-income pairs have at least one
material downward drop in owner-entry probability as liquid wealth increases.
With `Nb=120`, that share is still `81.3%`.

The worst baseline downward drop is `-0.583`. In the plotted policy slices,
owner-entry probabilities spike sharply near low liquid wealth and then drop.
Some threshold behavior is economically expected from down-payment and PTI
screens, but this shape is too jagged to treat as a credible smooth lifecycle
margin.

The deliberate room ladder reduces the material-drop share, especially with
`Nb=120`, but it does not repair the lifecycle ownership profile and it
destroys the scalar fit. This points to a joint menu/finance/policy-function
problem, not a small plotting artifact.

### Lifecycle Ownership

The central economic failure survives all variants:

| Variant | Prime ownership | Old ownership | Old minus prime |
|---|---:|---:|---:|
| `baseline_nb60_cluster_ladder` | `0.222` | `0.783` | `0.561` |
| `fine_nb120_cluster_ladder` | `0.184` | `0.758` | `0.574` |
| `baseline_nb60_room_ladder` | `0.225` | `0.851` | `0.626` |
| `fine_nb120_room_ladder` | `0.174` | `0.747` | `0.573` |

This is not an acceptable lifecycle tenure profile. The current objective
allows the optimizer to match old ownership and some fertility/housing moments
while suppressing prime-age ownership.

## Interpretation

The pathologies are mixed, but the main problem is not only numerical:

1. `Nb=120` improves some midlife distribution smoothness, so coarse asset-grid
   interpolation contributes to the visual problem.
2. The temporary owner ladder materially contaminates the fit. The owner-room
   moment is not robust to using a production-style ladder.
3. The owner-entry policy remains highly nonmonotone even with a finer asset
   grid.
4. The lifecycle ownership path remains economically broken under every
   variant.

Therefore, the next calibration should not simply continue optimizing the
current `candidate_no_timing_v0` objective. Before another long run, the model
needs:

- a deliberate owner ladder that includes the empirical room supports being
  targeted;
- lifecycle ownership-shape restrictions or moments, especially young and
  prime-age ownership;
- a policy-function audit of down-payment/PTI feasibility sets and owner-entry
  probability monotonicity;
- distribution diagnostics that distinguish deterministic entry mass from
  endogenous lifecycle concentration;
- possibly a finer asset grid for production diagnostics, even if the search
  itself uses a coarser grid.

