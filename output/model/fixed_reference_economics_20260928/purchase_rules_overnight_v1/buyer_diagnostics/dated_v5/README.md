# Accepted hard T48 date-zero financial access

Saved-policy buyer diagnostics jobs **19027866** (80% control) and **19027868** (100% temporary) both completed with exit `0:0` in 14 seconds. The two accepted source packets, exact input hashes, output hashes, and selected-fit identity are in `../v5_submission_receipt.json`. All six collected output files match the remote result hashes. No Bellman, equilibrium, or dated root solve was run by these readouts.

| Dated packet | Observed financed share | First-birth flow from renter origins | Among that flow, financially ineligible at 80% but eligible at 100% |
| --- | ---: | ---: | ---: |
| Hard control, date zero | 80% | 0.0462729133 | 61.8292887% |
| Hard temporary, date zero | 100% | 0.0460044154 | 61.6134213% |

These matched-state financial-access shares use each path's own date-zero distribution. They do not estimate desired ownership or the causal first-birth response. The [control summary](case_00_hard_control_h48_date_000/dated_access.json) and [temporary summary](case_01_hard_temporary_h48_date_000/dated_access.json) contain the full native access and observed-choice audits. The mechanism's accepted dated-path report is the source for the first-birth response.
