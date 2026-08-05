# Intergenerational Housing Fertility Model

This package is the new one-market quantitative implementation for the
intergenerational housing allocation and fertility project. It starts from the
active workhorse lifecycle code under `code/model/dt_cp_model/`, but keeps that
code path unchanged.

For the simplest current-model run, open and run:

```bash
code/model/run_intergen_model.py
```

That file lists the current source theta, target set, grid, and output folder
at the top.

In Spyder, open `code/model/run_intergen_model.py` and press Run. The file runs
the model through `code/model/.venv/bin/python` without replacing the Spyder
kernel, then loads `solution_summary`, `moments`, `target_fit`, `age_profiles`,
`room_bin_fit`, `first_look_density_path`, `first_look_total_wealth_path`,
`first_look_total_wealth_density_path`, `solution_cache_path`,
`first_look_policy_lines`, `first_look_market_summary`, and the total-wealth
density/support tables for inspection in the Variable Explorer.
The runner also prints macOS `CPU_Speed_Limit` before and after the solve. If a
slow run has a low speed limit, the same model path is being throttled by the
laptop rather than taking a different equilibrium path.
After one successful solve, set `REFRESH_PLOTS_FROM_SAVED_SOLUTION = True` in
the runner to rebuild the graph packet from `solution_cache.pkl` without
solving the model again. Leave `FAST_REFRESH_FROM_SAVED_SOLUTION = True` for
quick first-look plot iteration; it skips the larger standard-diagnostic,
contact-sheet, age-profile, room-bin, rung, and threshold-plot rebuilds while
still refreshing the first-look plots and CSVs.

Run from `code/model` when using the package CLI directly:

```bash
.venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
.venv/bin/python -m intergen_housing_fertility.cli solve --max-iter-eq 20 --quiet
.venv/bin/python -m intergen_housing_fertility.cli diagnostics --fixed-prices --outdir ../../output/model/intergen_housing_fertility_smoke_fixed --quiet
.venv/bin/python tools/build_intergen_mechanics_packet.py --max-iter-eq 10
```

The current pass uses one aggregate housing-services market, a 4-year decision
period, one dependent-child stage, a persistent Markov income state, lifecycle
income by age, owner housing rungs, continuous renter housing, down-payment
constraints, and transaction/sale wedges. The live diagnostic convention caps
renter housing at `hR_max=6.0` while the owner ladder remains
`H_own=[2,4,6,8,10]`; this keeps strong tenure/product segmentation while
allowing a thin first large-rental margin and reserving the upper 8/10-room
rungs for owners. A payment-to-income screen exists as an optional diagnostic
switch, but the default model leaves it off and uses the collateral/down-payment
constraint as the active finance restriction. It is not calibrated.

`IMPLEMENTATION_STATUS.md` is the live implementation record. Any future
simplification or deferred object should be added there in the same edit as the
code change.

`tools/build_intergen_mechanics_packet.py` is a non-production inspection
driver for the June 2026 one-market strand. It reads a saved theta, re-solves
the active intergen model, and writes standard diagnostics, target-fit tables,
room-bin/rung shares, simple and full first-look policies/market panels, the
same first-look policy panels against total wealth after the chosen tenure
branch, standalone aggregate liquid-wealth and total-wealth density plots with
childless/parent and renter/owner splits, wealth-support coverage diagnostics,
owner-rung diagnostics, owner-entry thresholds, a local solved-object cache, and
optional policy
proof-of-concept cases under
`output/model/intergen_mechanics_packet_YYYYMMDD/`.

## E-strand hardening

Author decision, August 5: the certified pre-E6 E5b specification remains the
maintained model. E6a, E6b, E6c, and their refits are retained below only as
default-off experimental and reproducibility paths; none is adopted.

The autonomous E-strand hardening program starts from the certified E5b
twelve-row system. Its first default-off extension is E6a, an externally
fitted late-age fecundity tail with no added free parameter. Reproduce the
strict fixed-winner frontier with
`intergen_eqscale_seq_optimized/run_e6a_fecundity_tail_frontier.py`; the
cluster refit uses `code/cluster/submit_intergen_e6a.sh` followed by its
strict collector.

E6b adds an externally measured permanent income level at entry. The PSID
long-lag minimum-distance decomposition fixes log-level variance at `0.3931`;
`e6b_profile.py` maps it to a mean-one, three-node Gauss--Hermite rule with
weights `(1/6, 2/3, 1/6)` and takes the Cartesian product with the existing
five-state FL--HSV Rouwenhorst chain. The resulting 15-state transition
matrix is block diagonal, so permanent levels never change. It adds no free
calibration parameter and is default-off.

`run_e6b_income_level_diagnostic.py` strictly compares the E5b control, E6b
alone, and E6a+E6b at the certified E5b parameter vector. It reports all
twelve target rows plus childlessness, completed fertility, and ownership by
permanent level. `code/cluster/submit_intergen_e6b_fixed_diagnostic.sh` runs
that comparison on Torch.

E6c is a default-off binary readiness state. An irreversible logistic
age-arrival moves a childless household from unsettled to settled, and
first-child entry is available only when settled. `e6c_profile.py` registers
the location and spread domains; those two free coordinates are identified by
the two signed timing moments, making the combined E6a+E6b+E6c refit a
twelve-parameter / twelve-moment system. Before refitting,
`run_e6c_timing_jacobian.py` checks the strict local timing Jacobian. The
production launchers are `code/cluster/submit_intergen_e6abc.sh` and its
strict collector.

After strict collectors finish, `build_e6_decision_tables.py` verifies that
every report uses the same twelve targets and weights, then writes a ranked
loss summary, the complete target-fit deltas against E5b, and every free
parameter with bounds and external-restriction checks. It reads certified
reports only and does not run the model.

`build_e6_winner_diagnostics.py` reconstructs an E6 winner through the exact
calibration-chain runtime, requires all twelve certified moments to reproduce
within `1e-6`, and then writes the existing standard graph set. It also writes
a clearly labeled supplemental model-versus-NCHS first-birth age-bin table,
plot, and preregistered timing-shape flag without changing the standard set.
The Torch launcher is
`code/cluster/submit_intergen_e6_winner_diagnostic.sh`.

`build_e6_plain_comparison.py` compares a certified E6AB refit under the
block-equal proportional-gap objective with the canonical E6AB rescue. It
requires the same twelve-target contract, evaluates both estimates under both
objectives, reports every free parameter and bound, and writes the standard
signed-gap comparison plot. It is a reporting driver and does not solve the
model.

The follow-up L1 diagnostic uses
`--weight-scheme target_relative_block_equal_l1`. Its loss is the sum, across
the three economic blocks, of each block's mean absolute proportional target
gap. Unlike the squared version, it does not give a disproportionately large
penalty to the worst proportional miss. The production Torch launcher is
`code/cluster/submit_intergen_e6ab_l1.sh`; its strict collector is the matching
`submit_intergen_e6ab_l1_collector.sh`.

`build_e6_l1_comparison.py` evaluates the canonical, squared-proportional,
and absolute-proportional certified winners under all three criteria, raw
absolute loss, and mean absolute percentage loss. It writes the complete
twelve-row target comparison, all ten parameters and bounds, and a
three-estimate signed-gap plot. The standard winner diagnostic reads the
certified objective scheme, so its loss and row contributions reproduce either
quadratic or absolute-gap reports exactly.
