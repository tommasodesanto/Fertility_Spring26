clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/ado/"
adopath ++ "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/ado/"
log using "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/full/sa_rooms_first_birth_household_aligned_v1.log", replace text
use "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/analysis_sample.dta", clear
do "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/reference_check.do" full "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/full"
