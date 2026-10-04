# Bounded diagnosis of task 6 inherited-state rejection

Verified 2026-10-04 by transition_contract_review. Read-only inspection of current source and remote task 6 records; no model solve, numerical gate change, or source edit.

## Established evidence

Local receipt: `deployment_v6/failure_review/task_6.json` in this handoff folder. Start psi was 0.15745023418116633 (0.88 times the selected reference). Failure type was `InheritedDistributionInfeasible`, with total rejected mass 2.916407549662773e-11.

Remote evidence root: `/scratch/td2248/projects/current_estate_transition_20261003_v6/results/panel_guess_6/run/candidate_0001/horizon_024/`.

- `driver.log` identifies `e5f_ssj_scaled_step_root.py:128`, `current = sample(projected(p0), phase='initial')`: this was the INITIAL trial of the first candidate's 24-date path root, after a stationary endpoint had been obtained.
- `map_001/dated_audits.json` contains only completed period 0, with every recorded gate true. `heartbeat.json` reports `phase=forward, period=0`. The next guard rejection therefore occurred at period 1, year 2011 under the unchanged four-year calendar.
- `map_001/inherited_state_evidence/inherited_infeasible_1791093666341315314.json` reports three affected cells, all renters (`tenure=0`), zero wealth (`wealth_index=46`), age 30 (`age_index=3`), income index 0, location 0. Their children-ever-born/children-at-home indices are (1,1), (2,1), and (2,2). Masses are 2.8970890659652093e-11, 1.2171913891416726e-13, and 7.146569806146772e-14. Values are -9999999999.72225. The price is 0.7723857081350192. Origin mass is 1.000948933926429; projection mass is zero and `distribution_modified=false`.
- The endpoint's `latest_completed.json` reports price 0.690643719043182. Initial path construction is linear (`one_shock_floor.py:806`), from reference price 0.776101253093739 to this endpoint. The rejected price equals period 1 of that linear path, with price step approximately -0.00371554495872. `warm_start.json` labels this `fresh_measured_seed_initialization`.

## Relevant implemented conditions

`code/model/tools/run_dynamic_population_transition.py:416-453` checks the exact pre-fertility distribution against its dated policy. Positive mass with value at or below -1e9 (or nonfinite value) is recorded; rejection occurs above the existing feasibility mass tolerance, here 1e-12. It preserves the distribution and applies zero projection.

`code/model/experiments/birth_count_choice/model/engine/household.py:34-48` gives renters a zero borrowing floor when the fixed unsecured credit limit is zero. `engine/kernels.py:469-474` returns the infeasible sentinel for renter surplus

\[
R_t b+y_t-\bar c-r_t\bar h-b'\leq 10^{-10}.
\]

The renter block sets the saving lower bound from that unsecured floor and the grid minimum (`kernels.py:690-708`). Child-state-specific minimum consumption and housing enter through `cb_v` and `hb_v` (`household.py:222-223`). Values at the sentinel can also depend on a dead continuation; the rejection record alone does not distinguish these causes.

`code/model/tools/run_e5f_perfect_foresight_transition.py:358-375` implements

\[
r_t=(R+\delta+\tau_H)p_t-p_{t+1}.
\]

A falling trial price path can therefore raise current rent through the capital-loss term even when current house prices are lower. The affected states are renters, so a direct owner collateral/default rejection is ruled out for these cells.

The forward evaluator supplies current dated parameters/current price and continuation value `values[period+1]` (`run_e5f_perfect_foresight_transition.py:690-713`), advances the actual current evaluation (`:737-742`), and uses that resulting next distribution (`:872-873`). The inherited state observer captures copies after evaluation (`floor_runtime.py:597-602`); it does not replace the advancing state. No timing or state-capture defect was demonstrated by this bounded inspection.

## Interpretation and smallest next check

The evidence establishes rejection of an initial numerical trial path. It does **not** establish that the equilibrium for this shock is economically infeasible or unreachable. Current rental-budget infeasibility, dead continuation, or grid interpolation into a later dead state remain possible explanations; the census does not prove which mechanism generated the small mass.

The smallest useful next check is a three-cell audit using the existing date-0 diagnostic packet and dated path inputs: report each cell's period-1 renter resources minus minimum consumption/housing costs at the saving lower bound; identify its incoming date-0 contribution, including saving interpolation and income transitions; and inspect the corresponding continuation dead mask. No full equilibrium solve is needed if the dated values are available. No mass projection, relaxed gate, or modification to tonight's running pinned source is proposed.
