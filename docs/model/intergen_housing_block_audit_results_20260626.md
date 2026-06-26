# Intergen Housing-Block Audit Results

Date: 2026-06-26

This note records the first housing-first audit run after the June 26 calibration
readout. The point is to diagnose the housing and wealth block before doing
more fertility calibration.

## Source

- Baseline cache:
  `output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/solution_cache.pkl`
- Baseline candidate:
  `de_g044_i022`
- Target set:
  `candidate_replacement_young_old_roomgap_v1`
- Baseline grid:
  `J=16`, `Nb=60`, `income_states=5`, `H_own=[2,4,6,8,10]`,
  `hR_max=6.0`, `interp_method=linear`

## New Reproducible Packet Command

The read-only packet builder is:

```bash
code/model/.venv/bin/python code/model/tools/build_intergen_housing_block_packet.py \
  --cache output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/solution_cache.pkl \
  --outdir output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/housing_block_audit
```

It evaluates policies using the active post-tenure convention: for each current
state, it integrates over the tenure logit probabilities, maps the household to
branch liquid wealth after the tenure transaction, then evaluates consumption,
savings, and housing on the target-tenure policy branch.

This is the object to use going forward, not the earlier ad hoc line plot that
read `c_pol` at the current-tenure index.

## Output Locations

- Baseline housing-block packet:
  `output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/housing_block_audit/`
- Grid and owner-ladder audit root:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/`
- Grid comparison plot:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/current_liquid_consumption_grid_comparison.png`
- Owner-ladder comparison plot:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/current_liquid_consumption_ladder_comparison.png`
- Summary table:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/summary.md`
- Curve pathology table:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/curve_pathology_summary.md`

## Main Findings

### 1. The grid is not innocuous

The default liquid-wealth grid is much wider than the occupied support.
`grid_interp_lab.py --diagnose` on the baseline cache reports:

- grid range: `[-35, 100]` with `Nb=60`;
- occupied support at mass `>1e-9`: `[-6.56, 16.89]`;
- 99.9 percent wealth quantile: `7.56`;
- mass at `b>=6`: `0.017`;
- median spacing above `b=10`: `4.611`;
- suggested dense core: roughly `[-8, 11]`, sparse buffer to about `25`.

Fixed-price dense-grid runs changed the solution materially even without
changing theta or equilibrium price:

| case | loss | young own | own | old own | TFR |
|---|---:|---:|---:|---:|---:|
| current grid `Nb=60`, GE | `17.0710` | `0.124` | `0.684` | `0.978` | `1.867` |
| dense core `Nb=120`, fixed price | `15.0458` | `0.174` | `0.690` | `0.976` | `1.887` |
| dense core `Nb=240`, fixed price | `15.3824` | `0.194` | `0.698` | `0.976` | `1.891` |
| dense core `Nb=120`, GE re-solve | `15.0154` | `0.172` | `0.689` | `0.974` | `1.886` |

The GE dense-grid price moved only from `0.7946` to `0.7960` and cleared tightly
with residual `9.39e-07`. So the improvement is not just a fixed-price
disequilibrium artifact.

Interpretation: do not treat the current `Nb=60`, `[-35,100]` grid as
production-ready. Calibration and policy results are grid-sensitive.

### 2. Grid refinement does not remove the policy wiggles

The active-policy consumption function remains non-smooth after the dense-grid
audit.

Largest significant adjacent drops in \(E[c\mid b]\), restricting to exact
wealth points with mass at least `1e-4`:

| case | largest drop | drop at `b` | mass at point |
|---|---:|---:|---:|
| current grid `Nb=60`, GE | `-0.5039` | `-2.6667` | `0.0017` |
| dense core `Nb=120`, GE | `-0.3426` | `0.1465` | `0.1025` |
| dense core `Nb=240`, fixed price | `-0.6068` | `4.2487` | `0.0150` |

The current-grid high-mass zero kink is not an artifact of the earlier plot
convention. Under the active-policy convention, the baseline still has a large
drop around the near-zero liquid-wealth mass point:

- nearest-zero grid node: `b=0.146514`;
- mass at that node: `0.137367`;
- adjacent consumption drop there: about `-0.432`.

On the dense `Nb=120` GE re-solve, mass is split across exact `b=0` and
`b=0.146514`, but a large drop remains:

- mass at `b=0`: `0.0470`;
- mass at `b=0.146514`: `0.1025`;
- largest significant drop: `-0.3426` at `b=0.146514`.

Interpretation: grid redesign is necessary, but it is not sufficient. The
remaining wiggles are tied to tenure/product threshold behavior.

### 3. Total wealth does not make the problem disappear

The packet defines post-tenure total wealth as
\[
W' = b^{branch} + (1-\psi)pH
\]
for owner outcomes, and \(W'=b^{branch}\) for renter outcomes. This is the
right convention for active consumption policies because consumption is
evaluated after tenure choice.

The exact \(E[c\mid W']\) line is still very jagged. Some of this is mechanical:
post-tenure total wealth interleaves renter outcomes and several discrete owner
rungs. But this is itself informative. The current model does not have a smooth
housing-equity state; it has discrete tenure/product jumps plus liquid wealth.

Interpretation: total wealth plots should be kept, but exact total-wealth lines
must be read together with tenure/rung composition. Binned or tenure-split
versions can be added for exposition, but the exact line exposes a real
discreteness problem.

### 4. The owner ladder is a structural object, not a harmless grid

Fixed-price owner-ladder variants at dense `Nb=120` show that the owner menu
changes behavior sharply:

| owner menu | loss | young own | own | old own | TFR |
|---|---:|---:|---:|---:|---:|
| `[2,4,6,8,10]` dense `Nb=120`, fixed price | `15.0458` | `0.174` | `0.690` | `0.976` | `1.887` |
| `[3,4,5,6,8,10]`, fixed price | `125.6887` | `0.335` | `0.753` | `0.997` | `1.896` |
| `[4,5,6,7,8,10]`, fixed price | `17.3794` | `0.111` | `0.687` | `0.981` | `1.863` |
| `[4,5,6,7,8,9,10]`, fixed price | `17.5256` | `0.111` | `0.687` | `0.981` | `1.863` |

The `[3,4,5,6,8,10]` menu gets young ownership near the target, but at the
fixed old theta/price it drives old ownership essentially to one and blows up
the objective. Removing the `H=2` rung and starting at `H=4` avoids the massive
loss explosion but does not smooth the consumption function.

Interpretation: the owner menu is doing real economic work. The bad `H=2` rung
is not simply removable without rethinking owner entry, old retention, and
price/utility calibration. A better owner-product design may still be needed,
but it cannot be treated as a cosmetic grid refinement.

### 5. Old ownership is mechanically near absorbing

The housing-block packet computes old-owner transition shares, ages 65-75, using
tenure logit probabilities:

| transition | share |
|---|---:|
| same rung | `0.9727` |
| downsize | `0.0107` |
| upsize | `0.0150` |
| sell to rent | `0.0016` |

Interpretation: the old-age ownership miss is a housing-block mechanism
failure. The model has almost no old-owner exit margin. This is separate from
fertility and should be fixed before interpreting fertility policy
counterfactuals.

## Bottom Line

The audit supports three conclusions.

1. The wealth grid must be redesigned before another serious calibration. The
   current grid is too wide and too sparse in relevant regions, and dense-grid
   re-solves materially change moments.
2. The weird consumption/housing shapes are not only grid coarseness. They
   persist under denser grids and are tied to discrete tenure/product threshold
   behavior.
3. The housing block needs a mechanism audit before fertility work: owner menu,
   branch wealth accounting, young owner value gaps, and old-owner exit must be
   separated.

## Recommended Next Step

Do not recalibrate yet. The next safe step is a targeted mechanism packet on
the dense `Nb=120` GE re-solve:

1. value gaps between renting and each owner rung around `b=0`;
2. exact purchase cash needs \((1-\phi)pH\);
3. branch liquid wealth and post-tenure total wealth around tenure switches;
4. old-owner value gaps for same-rung, downsize, sell-to-rent;
5. a binned total-wealth plot split by target tenure/rung.

Only after that should we decide whether the fix is mainly a better wealth
grid, a redesigned owner ladder, an old-age exit mechanism, or a deeper
housing-equity/mortgage-state change.
