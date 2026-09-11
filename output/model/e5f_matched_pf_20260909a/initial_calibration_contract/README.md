# Complete early working minimum-distance objective

The lead decision and working_weights.csv retain all13 restrictions: twelve
scored rows and the separately normalized2.1. working_contract.json freezes the
observation/provenance/weight choices. It is a working diagonal criterion, not
an optimal covariance estimate or a certified production benchmark.

Completed and checked after the September11 morning restart:

- score_initial.py: pure scorer,18 tests passed.
- complete_panel_analysis.py: all19 saved sensitivity cases,247 fit rows,
  323 parameter rows,108 derivatives; original eleven-row and recent-parent
  derivatives independently reproduced. Input hashes in the output JSON.
- score_saved_cases.py: full starting-point and reproduced-joint scores under
  the frozen contract, supported by prior verified cluster collections and
  fresh checks of saved observation hashes. Remote bytes were not re-read during
  the current authentication outage. Both complete fits/parameter tables and
  scoring receipts are in saved_case_scores/.

Baseline loss1499.851825; joint loss1428.171606. The recent-parent row accounts
for84.27% and97.03%, respectively. No final parameter estimate is selected.
All source observations retain their original diagnostic flags. The explicit
maintained recent-parent proxy is certified for this working comparison, not
claimed to reconstruct exact annual ACS interviews.

**Unfinished after credit interruption:** run_scored_candidate.py has not had
its requested independent review/tests or an exact cluster-loop smoke. Do not
use it to launch a search yet. The isolated historical-continuation driver also
contains an unverified partial edit; do not infer a committed/deployed change.

Current working-contract canonical SHA256:
c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2.
Byte SHA256:
e0bd8316a19bb197ab0fe9adaf25cb3173ff3b4072c14515cdd4e34b516ff43c.
