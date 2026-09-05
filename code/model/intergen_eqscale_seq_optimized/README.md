# Sequential housing and fertility model

Read [`../../../CALIBRATION_STATUS.md`](../../../CALIBRATION_STATUS.md) for the
maintained calibration, outstanding issues, and approved policy contract.
This package contains the sequential household solver; a bare package CLI does
not reproduce the dated calibration's parameter, target, and population contract.

## Where to start

| Purpose | File / saved artifact |
|---|---|
| Understand the household state, budgets and choices | `parameters.py`, `solver.py` |
| Understand the dated calibration and target measurement | `../tools/run_e5f_transition_calibration.py` |
| Follow the population and housing-market calculation | `../tools/run_e5f_open_population_transition.py`, `../tools/run_dynamic_population_transition.py` |
| Inspect the calibration already computed, without a new solve | [Complete calibration review PDF](../../../output/pdf/e5f_bounded_calibration_refinement_review.pdf), including all target rows, parameters and the standard diagnostic plots |
| Read the code audit and proposed cleanup order | [Code review](../../../docs/model/e5f_full_code_correctness_efficiency_review_20260905.md) |

For exact reruns, use the source snapshot and fingerprinted commands recorded
with the selected result in `CALIBRATION_STATUS.md`. Later code repairs have
different source fingerprints; do not relax a production fingerprint gate to
make an old command run with new code.

`../run_intergen_model.py` is the **historical June one-shot diagnostic runner**.
Its source, target system and grid differ from the maintained sequential dated
calibration. It remains available to reproduce that historical strand.
`IMPLEMENTATION_STATUS.md` is likewise a June historical record, not live status.

## Historical E-strand experiments

The following records describe their named experiments at the time they ran.
Current inclusion of any mechanism is determined by `CALIBRATION_STATUS.md`.

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
