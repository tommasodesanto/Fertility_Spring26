version 17.0
clear all
set more off

* -------------------- USER PATHS --------------------
local ACS   "/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta"
local SAIZ  "/Users/tommasodesanto/Desktop/ACS_microdata/extract_10/saiz_stataformat.dta"
local METRO "main_metro.dta"
local OUT   "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs"
local TBL   "`OUT'/Tables_ACS"
local GFX   "`OUT'/Graphs_ACS"

cap mkdir "`OUT'"
cap mkdir "`TBL'"
cap mkdir "`GFX'"

* -------------------- LOAD ACS --------------------
use "`ACS'", clear

* Tenure DV (unambiguous)
cap drop ownhome
gen byte ownhome = .
replace ownhome = 1 if ownershp==1
replace ownhome = 0 if ownershp==2
label var ownhome "Homeowner (1=own, 0=rent)"

* New parent indicator (eldest child <4 vs none)
cap drop littlechild_var2
gen byte littlechild_var2 = .
replace littlechild_var2 = 1 if eldch<4 & eldch!=.
replace littlechild_var2 = 0 if nchild==0
label var littlechild_var2 "New Parent (eldest child <4 vs none)"

* Preferred income control
local INC inctot
capture confirm variable hhincome
if _rc==0 local INC hhincome

* MSA id (prefer met2013, else cbsa, else metarea)
cap drop msa
gen long msa = .
capture confirm variable met2013
if _rc==0 replace msa = met2013 if met2013<.
capture confirm variable cbsa
if _rc==0 replace msa = cbsa if cbsa<. & missing(msa)
capture confirm variable metarea
if _rc==0 replace msa = metarea if metarea<. & missing(msa)
drop if missing(msa)
label var msa "MSA id (CBSA if available, else metarea)"

* PUMA density (if present)
capture confirm variable density
if _rc==0 {
    cap drop lndensity_puma
    gen double lndensity_puma = log(density) if density>0
    label var lndensity_puma "Log density (PUMA)"
}

* -------------------- SAIZ ELASTICITIES (by CBSA) --------------------
tempfile SAIZTMP
preserve
capture confirm file "`SAIZ'"
if _rc==0 {
    use "`SAIZ'", clear
    local geocand GEOID geoid GEOID10 cbsa cbsa_code met2013
    local elascand saiz_eta elas elasticity eta supply_elasticity
    local GEO ""
    foreach v of local geocand {
        capture confirm variable `v'
        if _rc==0 local GEO "`v'"
    }
    local ELAS ""
    foreach v of local elascand {
        capture confirm variable `v'
        if _rc==0 local ELAS "`v'"
    }
    if "`GEO'"!="" & "`ELAS'"!="" {
        cap drop cbsa5
        cap confirm numeric variable `GEO'
        if _rc==0 gen long cbsa5 = `GEO'
        else {
            tostring `GEO', gen(_gstr) force
            replace _gstr = subinstr(_gstr,"US","",.)
            replace _gstr = subinstr(_gstr,"CBSA","",.)
            destring _gstr, gen(cbsa5) force
            drop _gstr
        }
        rename `ELAS' saiz_eta
        keep if cbsa5<. & saiz_eta<.
        keep cbsa5 saiz_eta
        save `SAIZTMP', replace
    }
}
restore

* Attach CBSA code in ACS for Saiz merge
cap drop cbsa5
gen long cbsa5 = .
capture confirm variable met2013
if _rc==0 replace cbsa5 = met2013 if met2013<.
capture confirm variable cbsa
if _rc==0 replace cbsa5 = cbsa if cbsa<. & missing(cbsa5)

capture confirm file `"`SAIZTMP'"'
if _rc==0 {
    merge m:1 cbsa5 using `SAIZTMP', nogen keep(master match)
    cap drop saiz_inveta
    gen double saiz_inveta = 1/saiz_eta if saiz_eta>0
    label var saiz_inveta "Inverse Saiz elasticity (tightness ↑)"
}

* -------------------- MSA DENSITY from main_metro.dta --------------------
tempfile MDEN
capture confirm file "`METRO'"
if _rc==0 {
    preserve
    use "`METRO'", clear
    local msa_id ""
    capture confirm variable met2013
    if _rc==0 local msa_id "met2013"
    if "`msa_id'"=="" {
        capture confirm variable metarea
        if _rc==0 local msa_id "metarea"
    }
    if "`msa_id'"!="" {
        capture confirm variable densityw2010
        if _rc==0 {
            keep `msa_id' densityw2010
            keep if densityw2010<.
            collapse (mean) densityw2010, by(`msa_id')
            gen long msa = `msa_id'
            gen double lndensity2010 = log(densityw2010) if densityw2010>0
            keep msa lndensity2010
            save `MDEN', replace
        }
    }
    restore
}
capture confirm file `"`MDEN'"'
if _rc==0 merge m:1 msa using `MDEN', nogen keep(master match)

* -------------------- SAMPLE --------------------
cap drop S
gen byte S = inrange(age,22,45) & (nchild==0 | (eldch<4 & eldch!=.)) & !missing(ownhome, msa, year, perwt)

* -------------------- Z-SCORES ON SAMPLE --------------------
cap drop z_inveta z_lndens_msa z_lndens_puma
egen z_inveta      = std(saiz_inveta)    if S & saiz_inveta<.
egen z_lndens_msa  = std(lndensity2010)  if S & lndensity2010<.
egen z_lndens_puma = std(lndensity_puma) if S & lndensity_puma<.

* -------------------- REGRESSIONS --------------------
local ctrls "c.`INC' i.educd i.empstatd i.sex c.age"

eststo clear

* (A) MSA Saiz: absorb CBSA FE (main level of Saiz is time-invariant -> dropped)
quietly count if S & !missing(z_inveta, cbsa5)
if r(N)>1000 {
    eststo saiz_msa: areg ownhome c.littlechild_var2##c.z_inveta `ctrls' i.year if S & !missing(z_inveta,cbsa5) [pw=perwt], absorb(cbsa5) vce(cluster cbsa5)
}

* (B) MSA Density: absorb MSA FE (main level of density is time-invariant -> dropped)
quietly count if S & !missing(z_lndens_msa)
if r(N)>1000 {
    eststo dens_msa: areg ownhome c.littlechild_var2##c.z_lndens_msa `ctrls' i.year if S & !missing(z_lndens_msa) [pw=perwt], absorb(msa) vce(cluster msa)
}

* Export Table 1: MSA Saiz + MSA Density (only NP main + interaction)
capture noisily esttab saiz_msa dens_msa using "`TBL'/tenure_newparent_MSA_saiz_density.tex", replace keep(littlechild_var2 *#*) b(3) se(3) label stats(N r2, fmt(0 3) labels("Obs." "R^2")) title("New parent × tightness: MSA (Saiz 1/elasticity and MSA density), FE absorbed")

* (C) PUMA Density: keep main + interaction (absorb MSA FE; PUMA varies within MSA)
quietly count if S & !missing(z_lndens_puma)
if r(N)>1000 {
    eststo dens_puma: areg ownhome c.littlechild_var2##c.z_lndens_puma `ctrls' i.year if S & !missing(z_lndens_puma) [pw=perwt], absorb(msa) vce(cluster msa)
    capture noisily esttab dens_puma using "`TBL'/tenure_newparent_PUMA_density.tex", replace keep(littlechild_var2 z_lndens_puma *#*) b(3) se(3) label stats(N r2, fmt(0 3) labels("Obs." "R^2")) title("New parent × tightness: PUMA density (with MSA and year FE)")
}

* Optional sanity lines (quiet)
quietly sum z_inveta if e(sample)
quietly corr z_inveta saiz_inveta if e(sample)
quietly tab msa cbsa5 if e(sample) & msa!=cbsa5 & !missing(msa,cbsa5)

di as result "Exported: `TBL'/tenure_newparent_MSA_saiz_density.tex"
di as result "Exported: `TBL'/tenure_newparent_PUMA_density.tex"
