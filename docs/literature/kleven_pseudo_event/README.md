# Kleven pseudo-event study source packet

Downloaded September 17, 2026 from the author's [research page](https://www.henrikkleven.com/research).

- [Paper, November 2025](Kleven_2025_pseudo_event_study.pdf).
- [Original replication code ZIP](replication_code.zip): 50 R files; no code executed or modified.
- [Replication README](replication_README.pdf).
- [Replication details](replication_details_2025.pdf).
- [Download URLs, sizes and SHA-256 hashes](manifest.json).

The full Dropbox archive advertises 4,986,377,562 bytes. Its partial download was discarded; only the small Code subfolder and the two guides were retained. No raw microdata were downloaded. The guides describe separate access requirements for PSID and GSS geography.

Key files inside the ZIP: `MASTER.R` (entire pipeline), `setup.R` (parameters), `clean_acs.R`, `clean_psid.R`, `matching.R` (ACS/CPS matching), `matching_panel.R` (panel validation), `functions.R`, `fig_validation_eventstudy.R`, and `table_validation_assumption_conditional_independence.R`.

Selected source checks: `matching.R:35–75` implements exact matching with replacement and all ties; `matching.R:106–109` applies parent eligibility; `matching.R:132–149` anchors matching on oldest child age zero and shifts age/year. `matching_panel.R:287` sets reference event time −2. `setup.R:99–100` sets first-birth ages 25–45. These are checked excerpts, not a full replication audit.

The [revised project plan](../../model/acs_fertility_pseudopanel_feasibility.md) follows this method. Preserve these downloaded originals; put any later adaptation under `code/empirical/acs/`.
