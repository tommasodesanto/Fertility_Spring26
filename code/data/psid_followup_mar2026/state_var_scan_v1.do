clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/state_scan_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/state_var_scan_v1.log", replace text

use "`dta'", clear

di as text "=== lookfor state ==="
lookfor state

di as text "=== lookfor stfips ==="
lookfor stfips

di as text "=== lookfor region ==="
lookfor region

di as text "=== lookfor move/mig ==="
lookfor mig
lookfor move

log close _all
