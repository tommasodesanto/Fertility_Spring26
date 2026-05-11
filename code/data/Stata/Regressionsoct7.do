*******************************************************
* ACS: Fertility × Housing × Location — final one-shot
*******************************************************
version 17.0
clear all
set more off

* -------------------- PATHS --------------------
local ACS   "/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta"
local SAIZ  "/Users/tommasodesanto/Desktop/ACS_microdata/extract_10/saiz_stataformat.dta"
local METRO "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/main_metro.dta"
local OUT   "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs"
local TBL   "`OUT'/Tables_ACS"
local GFX   "`OUT'/Graphs_ACS"

cap mkdir "`OUT'"
cap mkdir "`TBL'"
cap mkdir "`GFX'"

* -------------------- LOAD --------------------
use "`ACS'", clear

* Tenure DV
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

* -------------------- SAIZ (by CBSA) --------------------
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
        if _rc==0 & "`GEO'"=="" local GEO "`v'"
    }
    local ELAS ""
    foreach v of local elascand {
        capture confirm variable `v'
        if _rc==0 & "`ELAS'"=="" local ELAS "`v'"
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

* -------------------- DENSITY (build lndensity2010 robustly) --------------------
tempfile _MDEN
preserve
capture confirm file "`METRO'"
if _rc==0 {
    use "`METRO'", clear
    local msa_id ""
    foreach v in met2013 cbsa metarea msa {
        capture confirm variable `v'
        if _rc==0 & "`msa_id'"=="" local msa_id "`v'"
    }
    local dens ""
    foreach d in lndensity2010 densityw2010 density2010 dens2010 density dens {
        capture confirm variable `d'
        if _rc==0 & "`dens'"=="" local dens "`d'"
    }
    if "`msa_id'"!="" & "`dens'"!="" {
        tempvar _ld
        gen double `_ld' = `dens'
        quietly summarize `_ld'
        capture drop lndensity2010
        if r(min)>0 & r(max)>20 {
            gen double lndensity2010 = log(`_ld')
        }
        else rename `dens' lndensity2010
        keep `msa_id' lndensity2010
        rename `msa_id' msa
        collapse (mean) lndensity2010, by(msa)
        save "`_MDEN'", replace
    }
}
restore
capture confirm file "`_MDEN'"
if _rc==0 merge m:1 msa using "`_MDEN'", nogen keep(master match)

capture confirm variable lndensity2010
if _rc {
    local trydens densityw2010 density Density dens
    local found 0
    foreach d of local trydens {
        capture confirm variable `d'
        if _rc==0 & `found'==0 {
            tempvar _ld
            gen double `_ld' = .
            replace `_ld' = log(`d') if `d'>0
            bys msa: egen lndensity2010 = mean(`_ld')
            replace lndensity2010 = . if lndensity2010==. & missing(`_ld')
            local found 1
        }
    }
}
label var lndensity2010 "Log density (MSA)"

* -------------------- SAMPLE --------------------
cap drop S
gen byte S = inrange(age,22,45) & (nchild==0 | (eldch<4 & eldch!=.)) & !missing(ownhome, msa, year, perwt)

* Z-scores
cap drop z_inveta z_lndens_msa
egen z_inveta     = std(saiz_inveta)   if S & saiz_inveta<.
egen z_lndens_msa = std(lndensity2010) if S & lndensity2010<.

* -------------------- PRINCIPAL-CITY STATUS + FLOWS --------------------
cap drop pc_status
gen byte pc_status = .
replace pc_status = 1 if inlist(metro,1,3)
replace pc_status = 0 if inlist(metro,2,4)
label var pc_status "Lives in"

cap drop metro1y
gen byte metro1y = .
replace metro1y = 1 if inlist(migtype1,2,3)
replace metro1y = 0 if inlist(migtype1,1,4,5)

cap drop went_metro left_metro
gen byte went_metro = (pc_status==1 & metro1y==0)
replace went_metro = . if metro1y==. | pc_status==.
gen byte left_metro = (pc_status==0 & metro1y==1)
replace left_metro = . if metro1y==. | pc_status==.
label var went_metro "Moved To"
label var left_metro "Left From"

* -------------------- (A) COMBINED PRINCIPAL-CITY TABLE --------------------
preserve
keep if !missing(perwt, pc_status)
eststo clear
estpost summarize pc_status went_metro left_metro [aweight=perwt] if nchild==0
eststo np
estpost summarize pc_status went_metro left_metro [aweight=perwt] if littlechild_var2==1
eststo newp
estpost summarize pc_status went_metro left_metro [aweight=perwt]
eststo full
esttab np newp full using "`TBL'/principalcity_means_flows_ALL.tex", replace booktabs nonumbers ///
    title("Baseline Weighted Means by Population Group") ///
    collabels("Non-Parents" "New Parents" "Full Sample") ///
    cells(mean(fmt(%6.3f))) stats(N, fmt(%12.0f) labels("Observations"))
restore

* -------------------- (B) BASELINE HOMEOWNERSHIP TABLE --------------------
preserve
keep if S
eststo clear
estpost summarize ownhome [aweight=perwt] if nchild==0
eststo non_par
estpost summarize ownhome [aweight=perwt] if littlechild_var2==1
eststo new_par
estpost summarize ownhome [aweight=perwt]
eststo full_samp
esttab non_par new_par full_samp using "`TBL'/tenure_baseline_means.tex", replace booktabs nonumbers ///
    title("Baseline Homeownership (weighted means)") ///
    collabels("Non-Parents" "New Parents" "Full (22--45, NP or New Parent)") ///
    cells(mean(fmt(%6.3f))) stats(N, fmt(%12.0f) labels("Observations"))
restore

* -------------------- (C) TIGHTNESS: CLEAN TABLES (NO MOVERS) --------------------
capture noisily eststo clear
areg ownhome c.littlechild_var2##c.z_inveta ///
    c.`INC' i.educd i.empstatd i.sex c.age i.year if S & !missing(z_inveta,cbsa5) [pw=perwt], absorb(cbsa5) vce(cluster cbsa5)
eststo t_saiz

areg ownhome c.littlechild_var2##c.z_lndens_msa ///
    c.`INC' i.educd i.empstatd i.sex c.age i.year if S & !missing(z_lndens_msa,msa) [pw=perwt], absorb(msa) vce(cluster msa)
eststo t_dens

esttab t_saiz t_dens using "`TBL'/tenure_newparent_tightness_clean.tex", replace ///
    booktabs label b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nocons ///
    keep(littlechild_var2 c.littlechild_var2#c.z_inveta c.littlechild_var2#c.z_lndens_msa) ///
    coeflabel( littlechild_var2                    "New Parent" ///
               c.littlechild_var2#c.z_inveta      "New Parent $\times$ inverse \(h^s\) selasticity" ///
               c.littlechild_var2#c.z_lndens_msa  "New Parent $\times$ density") ///
    stats(N r2, fmt(0 3) labels("Obs." "R^2")) ///
    mgroups("Homeowner", pattern(1 1) span) mtitles("(1)" "(2)") ///
    title("New parent on MSA (inverse \(h^s\) selasticity and MSA density), FE absorbed") ///
    addnotes("Standard errors in parentheses", "* p < 0.10, ** p < 0.05, *** p < 0.01")



* -------------------- (D) TIGHTNESS × MOVER: TRIPLE-INTERACTION CLEAN TABLE --------------------
* Ensure mover flag
capture confirm variable movedthisyear
if _rc {
    cap drop movedthisyear
    gen byte movedthisyear = .
    replace movedthisyear = 0 if migrate1d==10
    replace movedthisyear = 1 if migrate1d!=10 & migrate1d<.
}
cap drop mover1y
gen byte mover1y = movedthisyear
label var mover1y "Mover"

* Estimate (Saiz + Density)
capture noisily areg ownhome c.littlechild_var2##c.z_inveta##i.mover1y ///
    c.`INC' i.educd i.empstatd i.sex c.age i.year if S & !missing(z_inveta,cbsa5,mover1y) [pw=perwt], absorb(cbsa5) vce(cluster cbsa5)
eststo saiz_triple

capture noisily areg ownhome c.littlechild_var2##c.z_lndens_msa##i.mover1y ///
    c.`INC' i.educd i.empstatd i.sex c.age i.year if S & !missing(z_lndens_msa,msa,mover1y) [pw=perwt], absorb(msa) vce(cluster msa)
eststo dens_triple

esttab saiz_triple dens_triple using "`TBL'/tenure_newparent_tightness_mover_triple_clean.tex", replace ///
    booktabs label b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nocons ///
    keep(littlechild_var2 1.mover1y ///
         c.littlechild_var2#c.z_inveta 1.mover1y#c.littlechild_var2 1.mover1y#c.z_inveta 1.mover1y#c.littlechild_var2#c.z_inveta ///
         c.littlechild_var2#c.z_lndens_msa 1.mover1y#c.z_lndens_msa 1.mover1y#c.littlechild_var2#c.z_lndens_msa) ///
    coeflabel( littlechild_var2                               "New Parent" ///
               1.mover1y                                      "Mover" ///
               c.littlechild_var2#c.z_inveta                  "New Parent $\times$ inverse \(h^s\) selasticity" ///
               1.mover1y#c.littlechild_var2                   "New Parent $\times$ Mover" ///
               1.mover1y#c.z_inveta                           "Mover $\times$ inverse \(h^s\) selasticity" ///
               1.mover1y#c.littlechild_var2#c.z_inveta        "New Parent $\times$ inverse \(h^s\) selasticity $\times$ Mover" ///
               c.littlechild_var2#c.z_lndens_msa              "New Parent $\times$ density" ///
               1.mover1y#c.z_lndens_msa                       "Mover $\times$ density" ///
               1.mover1y#c.littlechild_var2#c.z_lndens_msa    "New Parent $\times$ density $\times$ Mover") ///
    stats(N r2, fmt(0 3) labels("Obs." "R^2")) ///
    mgroups("Homeowner", pattern(1 1) span) mtitles("(Saiz)" "(Density)") ///
    title("Homeownership: New Parent $\times$ inverse \(h^s\) selasticity $\times$ Mover (and density), FE absorbed") ///
    addnotes("Standard errors in parentheses", "* p < 0.10, ** p < 0.05, *** p < 0.01")



display as result "WROTE:"
display as text "`TBL'/principalcity_means_flows_ALL.tex"
display as text "`TBL'/tenure_baseline_means.tex"
display as text "`TBL'/tenure_newparent_tightness_clean.tex"
display as text "`TBL'/tenure_newparent_tightness_mover_triple_clean.tex"
*******************************************************
* ===== REBUILD USING THE ORIGINAL DEFINITIONS =====
preserve

* 0) Original status logic and sample prep
drop if missing(metro)
drop if metro==4                         // exclude non-metro, as before

cap drop metroyn metro1yyn went_metro left_metro
gen byte metroyn = .                     // "Lives in" = Not Principal City
replace metroyn = 1 if metro==2
replace metroyn = 0 if inlist(metro,1,3)
label var metroyn "Lives in (Not Principal City)"

gen byte metro1yyn = .                   // Metro last year (migtype1)
replace metro1yyn = 1 if inlist(migtype1,2,3)
replace metro1yyn = 0 if inlist(migtype1,1,4,5)

gen byte went_metro = (metroyn==1 & metro1yyn==0)
replace went_metro = . if metro1yyn==. | metroyn==.
label var went_metro "Moved To"

gen byte left_metro = (metroyn==0 & metro1yyn==1)
replace left_metro = . if metro1yyn==. | metroyn==.
label var left_metro "Left From"

* 1) LEGACY TABLE (2 columns): NP+New Parents vs Full
eststo clear
keep metroyn went_metro left_metro perwt littlechild_var2 nchild

estpost summarize metroyn went_metro left_metro [aweight=perwt] ///
    if littlechild_var2==1 | nchild==0
eststo legacy_specific

estpost summarize metroyn went_metro left_metro [aweight=perwt]
eststo legacy_overall

esttab legacy_specific legacy_overall using ///
"/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables_ACS/principalcity_means_flows_LEGACY.tex", ///
replace booktabs label nonumbers ///
mtitles("Non-Parents and New Parents" "Full Sample") ///
cells(mean(fmt(%6.3f))) ///
stats(N, fmt(%12.0f) labels("Observations")) ///
title("Baseline Weighted Means by Population Group (legacy)")

* 2) NEW SPLIT TABLE (3 columns) from the SAME logic
eststo clear

estpost summarize metroyn went_metro left_metro [aweight=perwt] if nchild==0
eststo NP

estpost summarize metroyn went_metro left_metro [aweight=perwt] if littlechild_var2==1
eststo NEWP

estpost summarize metroyn went_metro left_metro [aweight=perwt]
eststo FULL

esttab NP NEWP FULL using ///
"/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables_ACS/principalcity_means_flows_ALL.tex", ///
replace booktabs label nonumbers ///
mtitles("Non-Parents" "New Parents" "Full Sample") ///
cells(mean(fmt(%6.3f))) ///
stats(N, fmt(%12.0f) labels("Observations")) ///
title("Baseline Weighted Means by Population Group")

restore
*******************************************************
* ===== LIFECYCLE OWNERSHIP RATES FOR MODEL COMPARISON =====
*******************************************************
preserve

* Keep relevant sample
keep if !missing(ownhome, age, perwt)
keep if age >= 25 & age <= 84

* Create 5-year age bins
gen age_bin = 5 * floor(age/5)

* Collapse to get ownership rates by age bin
collapse (mean) own_rate=ownhome [pw=perwt], by(age_bin)

* Convert to percentage
gen own_pct = own_rate * 100
format own_pct %5.1f

* Display in console
list age_bin own_pct, clean noobs

* Export - use EXPLICIT path
export delimited age_bin own_pct using "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/ownership_by_age.csv", replace

* Graph
twoway (connected own_pct age_bin, lcolor(navy) mcolor(navy) msymbol(circle)), ///
    xlabel(25(5)80) ylabel(0(20)100, grid) ///
    xtitle("Age") ytitle("Homeownership Rate (%)") ///
    title("Lifecycle Homeownership (ACS)") ///
    graphregion(color(white))
graph export "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Graphs_ACS/lifecycle_ownership_ACS.pdf", replace

restore
