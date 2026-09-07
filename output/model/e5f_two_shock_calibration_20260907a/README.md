# Two-shock calibration test — active

The author clarified that the priority is actual calibration fit and parameter
re-estimation, not additional fixed-price diagnostics. The first full historical
objective evaluations are Torch array17125770 (cases1–2), pending priority at
submission on September7. Each case runs the normalized old steady state and
five historical dates; no policy counterfactuals. These are two identical retained
anchors for an exact-loop/reproducibility smoke before the bounded search.

Source: experimental branch codex/two-contemporaneous-shocks,25ab08c4.
Remote snapshot: /scratch/td2248/projects/Fertility_Spring26_two_shock_calibration_20260907a.
Plan SHA256:19d96150d356fde8ee5b8532fee1aecfac3a2e8f5c7d61b4ec7ec42720df7a7e.

The preflight verifies all twelve target values/weights and all eleven original
parameter names and bounds against the retained parameter table. First/later
fertility scales are estimated, housing dispersion remains externally.005 and
supply elasticity.63. First-child jump upper remains.5, not the overnight GEV
expansion to2.0. No GEV lambda. All source, target, domain and numerical gates
are pinned; independent operator and calibration-contract tests pass.

Budget: two cases with one core and24GB each, six-hour hard cap per case.
The prior nested-GEV histories took roughly25–30minutes; the measured additive
Bellman is slower, so allow roughly2–4hours per history pending actual cluster
measurement. Do not claim an objective value until the historical gates pass.
Every candidate writes a heartbeat, full fit/parameter tables and terminal
checkpoint. No figures under the author's explicit instruction.

Once both histories pass and exactly reproduce, compare all twelve fit rows
with retained task010; use observed runtime to budget and launch the bounded
recalibration of all eleven parameters. A restored10-minute thread heartbeat
monitors this work. No production changes or policy runs.
