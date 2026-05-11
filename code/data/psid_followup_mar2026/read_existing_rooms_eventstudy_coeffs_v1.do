clear all
set more off
version 17.0

local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/read_existing_rooms_eventstudy_coeffs_v1"
cap mkdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/read_existing_rooms_eventstudy_coeffs_v1.log", replace text

tempname posth
postfile `posth' str20 series double coef_p3 coef_p5 using "`outdir'/existing_rooms_eventstudy_coeffs_v1.dta", replace

use "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/rooms_f_c_y_all_estimates.dta", clear
quietly summarize b if relative_time == 3
local f_p3 = r(mean)
quietly summarize b if relative_time == 5
local f_p5 = r(mean)
post `posth' ("first_birth_all") (`f_p3') (`f_p5')

use "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_second_birth_v1/rooms_s_c_y_all_repl_estimates.dta", clear
quietly summarize b if relative_time == 3
local s_p3 = r(mean)
quietly summarize b if relative_time == 5
local s_p5 = r(mean)
post `posth' ("second_birth_all") (`s_p3') (`s_p5')

postclose `posth'

use "`outdir'/existing_rooms_eventstudy_coeffs_v1.dta", clear
export delimited using "`outdir'/existing_rooms_eventstudy_coeffs_v1.csv", replace
list, clean noobs

log close _all
