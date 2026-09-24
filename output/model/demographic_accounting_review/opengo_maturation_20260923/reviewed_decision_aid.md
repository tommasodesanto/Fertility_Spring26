# Maturation and household accounting: decisions for September 24

OpenCode Go / Kimi K3 completed the requested source review in six minutes.
The lead checked the cited source rules and corrected material errors below.
This is a source review, not a checkpoint audit or a new numerical result.

## What is established in the inspected source

- The independent-child process removes each child from the at-home count with a constant per-period probability. The 18-year mean and four-year period imply 2/9 per draw. This is a geometric departure clock, not tracking each child's biological age. Sources: `parameters.py:1212` and `e5_maturation_repair.py:16` under `code/model/intergen_eqscale_seq_optimized/`.
- The stationary forward accounting credits measured departures times the entrant conversion factor; the E5 profile sets that factor to 1/2 (`e5_profile.py:182`, `solver.py:5653`). The dated transition driver instead schedules top-code-adjusted births times 1/2.1 in a queue (`code/model/tools/run_e5f_open_population_transition.py:302`, `:1228`). Queue length and timing must be taken from the run configuration; its metadata reports effect lag as waiting slots plus one, so four waiting slots at four-year periods mean 20 years.
- Household survival acts before advancing surviving dependents. Terminal households exit. No child reassignment is present in the inspected advancement path (`solver.py:4871`, `run_e5f_open_population_transition.py:868–908`). The transition's separate birth queue can still generate future entrants: disappearing from the household dependent count is not the same as disappearing from that queue.
- The household child count is capped, while aggregate renewal expands the top birth category. These are different units. Neither normalizing household mass nor matching an external age distribution establishes conservation of people across the two systems.

## Worked examples

Two couples each have one child: initially two households and six people. When the two children form one new couple, while parents remain alive, there are three households and still six people. This is the accepted abstraction; it motivates separating household counts from person counts.

For 100 child units, literal pairing gives 50 new couples. The separate 1/2.1 convention gives 47.619 before retention. The 4.762% difference is an arithmetic normalization difference, not evidence that this share actually dies or fails to form households.

The worker's parental-death example overstated loss because it omitted survival through previous parental ages. For an illustration with births at parent age42, six maturation draws before66, per-draw retention7/9 and survival0.939/0.918/0.885/0.830 at66/70/74/78, start with C66=100(7/9)^6. At each age, loss=(1-s)C; next C=s(7/9)C; remove all remaining dependents at82. This gives10.479 lost dependent units, rather than the worker's14.2. Births at22 give2.983 under the same conventions. These are hand examples, not measured population losses, bounds, or verified checkpoint settings. The exact birth-stage timing matters.

## Tomorrow's decisions

| Order | Choice to settle | Consequence and next operation |
|---|---|---|
| 1. Maturation | Keep the geometric18-year average, add a parent-age cutoff, or track actual child ages? | A cutoff changes the duration of child costs and housing needs, particularly for late births. Existing ramp code is a candidate, not an adopted or validated solution. Specify the intended clock before recalculating policies. |
| 2. Parental death | What happens to dependent children when the household exits? | First establish actual exposure in a selected checkpoint. A cutoff before mortality can remove exposure under that survival schedule but forces earlier departures; reassignment requires an explicit accounting rule. No automatic model addition is recommended. |
| 3. Adult entry | Should entries follow measured departures or a separately specified birth-cohort clock? | Specify one coherent person/household identity, including retention and top-bin units, before revising stationary or transition accounting. Accepted two-adult abstraction does not by itself select a clock. |
| 4. Formation factor | Keep the documented1/2.1 convention distinct from literal1/2. | No change was authorized; empirical/stochastic formation-factor robustness remains low priority. |

## Corrections and limits

Do not use the unreviewed worker report's4–14% loss range, its assertion that no transition double-counting occurs, or its implication that normalizations prove conservation. Do not reuse the illustrative0.185 hazard as a calibrated alternative. Extra top-bin accounting can affect equilibrium through entry and population scale even if household state dimensions stay fixed; no promise of avoiding policy recomputation is warranted.

The actual saved checkpoint's maturation mode, survival, conversion factor and dependent-loss mass were not reverified. No run, source change, calibration choice or paper edit was made. Source hashes, execution evidence and lead corrections are in `lead_review.json`; the original report is preserved as `worker_report_unreviewed.md`.
