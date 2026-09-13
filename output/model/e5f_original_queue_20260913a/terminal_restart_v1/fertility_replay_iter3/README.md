# Native fertility diagnostics for the saved 100-period path

The requested simple view is `output/fertility_data_model.png` and `.pdf`:
orange permanent-shock model path and empirical data only. The old fitted
sequence is excluded. Regenerate with
`python code/model/tools/build_e5f_fertility_path_overlay.py --data-model-only`.
Exact artist/data checks and visual inspection pass; the model remains
unconverged iteration3. CSV and verification JSON accompany the figure.

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

## Mock transition panels

`output/macro_transition_mock.png/pdf` shows2003–2103: period fertility and
US four-year-average data, reconstructed completed fertility, housing services
demanded/supplied and household mass, and housing services per household.
`output/macro_transition_full_horizon.png/pdf` uses the full saved horizon.
These are the unchanged one-permanent-shock iteration3, **not** the new four-shock
announced path or a converged historical fit. All labels retain start-of-window
dating:2007 data are2008–2011 and2019 data are2020–2023. No model solve.

Completed fertility is reconstructed from saved age-specific birth rates along
cohorts, using the stationary prehistory. Age-only survival cancels in cohort
means; entry starts childless; non-birth transitions preserve children ever born.
The mean stock recursion is C[t,j]=C[t-1,j-1]+births[t,j]/mass[t,j]. Output uses
the model's final fertile age cell42–45 after its births and its fixed3+ bin
weight. It is not a directly saved-distribution measurement or an empirical
ages40–44 match. No CPS data overlay is asserted in that panel.

Both figure canvases were visually inspected. Sidecar
`output/macro_transition_verification.json` records source hashes, the native
flow/rate and independent cohort-diagonal checks, definitions and selected
values. `output/macro_transition_mock.csv` supplies all100observations.
Rebuild: `python code/model/tools/build_e5f_fertility_path_overlay.py --macro-transition`.

At2063, housing used is11.61% below initial, housing supplied9.81% below,
household mass13.59% below and services per household2.29% above. The
market gap remains1.99%, so neither housing curve is a cleared allocation.
At2103, housing used is22.15% below and household mass31.45% below.
