clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_second_birth_v1/rooms_s_c_y_all_repl_estimates.dta"

use "`dta'", clear
sort relative_time
list relative_time b se ci_lo ci_hi pre_event_mean, noobs sep(0)
