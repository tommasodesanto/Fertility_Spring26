clear all
set more off
version 17.0
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/covariance_replay/output"
matrix b_interact = (1.2345678901234567, -.5 \ 2, 3)
matrix rownames b_interact = 1970 1975
matrix colnames b_interact = L3event F1event
matrix V_interact = (.012345678901234567, .25 \ .16, .36)
matrix rownames V_interact = 1970 1975
matrix colnames V_interact = L3event F1event
matrix b_full = (1.2345678901234567, 2, -.5, 3)
matrix colnames b_full = c1 c2 c3 c4
matrix V_full = (.012345678901234567, .01, .02, .03 \ .01, .16, .04, .05 \ .02, .04, .25, .06 \ .03, .05, .06, .36)
matrix rownames V_full = c1 c2 c3 c4
matrix colnames V_full = c1 c2 c3 c4
local cohortnames : rownames b_interact
local eventnames : colnames b_interact
local rawnames : colnames b_full
local ncohort = rowsof(b_interact)
local nevent = colsof(b_interact)
local n_interactions = `ncohort' * `nevent'
file open bfh using "`outdir'/smoke_interaction_coefficients.csv", write replace
file write bfh "matrix_index,cohort_index,cohort,event_index,event,stata_interaction_name,estimate,marginal_variance" _n
forvalues i = 1/`n_interactions' {
    local event_index = ceil(`i' / `ncohort')
    local cohort_index = mod(`i' - 1, `ncohort') + 1
    local cohort : word `cohort_index' of `cohortnames'
    local event : word `event_index' of `eventnames'
    local rawname : word `i' of `rawnames'
    local beta = b_interact[`cohort_index', `event_index']
    local beta_full = b_full[1, `i']
    local vari = V_interact[`cohort_index', `event_index']
    local vari_full = V_full[`i', `i']
    assert abs(`beta' - `beta_full') <= 1e-14
    assert abs(`vari' - `vari_full') <= 1e-14
    local beta_export : display %24.17g `beta'
    local vari_export : display %24.17g `vari'
    file write bfh `"`i',`cohort_index',`cohort',`event_index',`event',`rawname',`beta_export',`vari_export'"' _n
}
file close bfh
file open vfh using "`outdir'/smoke_interaction_covariance.csv", write replace
file write vfh "row_index,column_index,covariance" _n
forvalues i = 1/`n_interactions' {
    forvalues j = 1/`n_interactions' {
        local cov = V_full[`i', `j']
        local cov_export : display %24.17g `cov'
        file write vfh `"`i',`j',`cov_export'"' _n
    }
}
file close vfh
import delimited using "`outdir'/smoke_interaction_coefficients.csv", clear asdouble
assert _N == 4
assert cohort[1] == 1970 & event[1] == "L3event"
assert abs(estimate[1] - 1.2345678901234567) < 5e-16
assert abs(marginal_variance[1] - .012345678901234567) < 5e-16
assert abs(estimate[2] - 2) < 1e-14
assert abs(estimate[3] + .5) < 1e-6
assert stata_interaction_name[4] == "c4"
import delimited using "`outdir'/smoke_interaction_covariance.csv", clear asdouble
assert _N == 16
assert row_index[2] == 1 & column_index[2] == 2
assert row_index[7] == 2 & column_index[7] == 3
assert abs(covariance[3] - .02) < 5e-16
assert abs(covariance[9] - .02) < 5e-16
di "COVARIANCE_EXPORT_SMOKE_PASS interaction_coefficients=4 full_covariance_cells=16"
