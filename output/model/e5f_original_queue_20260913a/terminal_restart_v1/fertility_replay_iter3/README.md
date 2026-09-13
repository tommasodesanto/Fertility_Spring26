# Native fertility diagnostics for the saved 100-period path

## Requested historical overlay

`output/fertility_overlay.png` and `.pdf` overlay the native iteration3 path
with the exact prior slide's historical fertility, constant-preference
continuation, and four-year data averages over2007–2039. Both use the prior
slide's start-of-window dating:2007labels the first window, whose empirical
average covers2008–2011. Source and artist values are checked without changing
any model output. CSV and verification JSON are alongside the figures.

This is a visual comparison, not a controlled timing experiment: the prior
recovered_sequence run uses no rebate, demographic conditioning and six-period
forecasts; the current run uses equal rebates, the original closed household
queue and100periods. The current root is unconverged and the prior horizon
remains unverified. The pre-shock2.1point belongs to the current IRF only.
The existing presentation and standard six-panel figure remain unchanged.

Regenerate with Python containing NumPy and Matplotlib:
`python code/model/tools/build_e5f_fertility_path_overlay.py`.
The lead visually inspected the complete one-panel output.

**Completed and collected September13 at20:24UTC.** All100 iteration3 rows
reproduce exactly. `output/lead_verification.json` records independent
age-specific fertility reconstruction and visual inspection. Deliverables:
`output/irf_fertility.png` and `.pdf`; full native data and terminal gaps are
alongside them. Fertility2.10 ->1.68186 on impact ->2.08277 at the final period;
year400 household mass0.387755 versus stationary0.349454. The path remains
unconverged; diagnostic replay success does not change its equilibrium status.

Job 17702691 is confirmed running on Torch (cs609), September 13 at 19:44 UTC.
It performs one exact replay of iteration 3, not a price-path search. Runtime
estimate is 28.3 minutes for one native mapping, with a 60-minute Slurm cap
and the original global deadline. One CPU, 32 GiB, numerical threads one.

All 300 saved coordinates, structural parameters, fiscal rules and original
household entry queues are frozen. Inputs and code are SHA-256 pinned. Native
household gates are reapplied; every saved dated row must reproduce to 2e-10.
Only then are fertility.json, row_reproduction.json, terminal_distance.json and
the new fertility-rate figure written. The figure includes the initial steady
state and impact jump. Every fertility rate must reproduce its native
age-specific birth flows and household masses. A successful replay is still
an unconverged equilibrium iterate, explicitly labeled as such.

Remote directory:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a_fertility_replay_iter3`

Collect its `output/` here, inspect row_reproduction.json and
irf_fertility_qa.json, visually inspect irf_fertility.png, then deliver the graph.
The updated half-hour thread monitor includes this job. It requires no laptop
connectivity after dispatch. Submission, immutable inputs, worker_contract.json,
launch.sh and replay.sbatch are retained alongside this README.

Lead checks: plotting validated using completed ten-period native observations;
300-coordinate validation and changed-row rejection passed; exact-loop smoke
receipt is required on the node. Lead repaired the delegated launcher's log
placement, NumPy import, input pin checks, terminal-value count (100 dated
values plus the terminal value), and post-replay figure generation before launch.
