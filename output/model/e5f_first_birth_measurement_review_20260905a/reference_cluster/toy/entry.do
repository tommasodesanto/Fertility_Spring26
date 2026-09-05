clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/ado/"
adopath ++ "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/ado/"
log using "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/toy/sa_rooms_first_birth_household_aligned_v1.log", replace text
do "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/reference_check.do" toy "/scratch/td2248/projects/e5f_first_birth_reference_20260905a/toy"
