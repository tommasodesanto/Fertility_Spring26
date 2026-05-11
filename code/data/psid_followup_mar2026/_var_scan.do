clear all
set more off
use "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", clear
capture log close _all
log using "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/var_scan.log", replace text

lookfor twin
lookfor multiple
lookfor child
lookfor birth
lookfor sex
lookfor os71
lookfor cah
lookfor relchi

capture noisily ds RELCHI*
capture noisily ds *OS71*
capture noisily ds *CAH*
capture noisily ds *SEX*
capture noisily ds *BIRTH*

log close _all
