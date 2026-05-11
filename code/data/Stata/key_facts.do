********************************************************************************
* KEY FACTS: Two analyses to guide model design
*
* Analysis 1: MSA-level fertility vs ownership/price-to-income
*   → Does fertility vary with housing affordability across MSAs?
*
* Analysis 3: Within-metro vs across-metro moves for new parents
*   → Are parents moving city→suburb or leaving metros entirely?
********************************************************************************

clear all
set more off

global path "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs"
global tables "$path/Tables_ACS"
global graphs "$path/Graphs_ACS"

* Output directory for new results
cap mkdir "$path/KeyFacts"
global out "$path/KeyFacts"

********************************************************************************
* Load data
********************************************************************************

use "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta", clear

* Basic sample restrictions
drop if age < 22
drop if age > 45   // Focus on fertile ages for fertility analysis
drop if gq != 1 & gq != 2  // Households only (not group quarters)

* Key variables
gen owner = (ownershp == 1) if ownershp != 0
gen renter = (ownershp == 2) if ownershp != 0

* Fertility measures
gen has_youngchild = (nchlt5 > 0 & nchlt5 != .)
gen newparent = (eldch < 4 & eldch != 99 & nchild > 0)

* Price to income ratio (owners only, exclude top/bottom coded)
gen pti = valueh / hhincome if owner == 1 & valueh > 0 & valueh < 9999999 & hhincome > 0 & hhincome < 9999999

* Migration variables (1-year)
gen moved1y = (migrate1 >= 2 & migrate1 != .) if migrate1 != .

* Metro status categories
* metro: 1=not metro, 2=principal city, 3=metro not principal, 4=mixed
gen in_principal = (metro == 2) if metro >= 1 & metro <= 3
gen in_metro_notprincipal = (metro == 3) if metro >= 1 & metro <= 3
gen in_metro = (metro == 2 | metro == 3) if metro >= 1 & metro <= 3

* Movement type classification using 2013 delineations
* migmet131 = MSA 1 year ago, met2013 = current MSA
gen same_msa = (met2013 == migmet131) if moved1y == 1 & met2013 != 0 & migmet131 != 0
gen diff_msa = (met2013 != migmet131) if moved1y == 1 & met2013 != 0 & migmet131 != 0


********************************************************************************
* ANALYSIS 1: MSA-level fertility × ownership × price-to-income
********************************************************************************

di as text ""
di as text "============================================="
di as text "ANALYSIS 1: MSA-level fertility vs affordability"
di as text "============================================="

preserve

* Restrict to identifiable MSAs
keep if met2013 > 0

* Collapse to MSA-year level
collapse (mean) fert_rate=has_youngchild own_rate=owner ///
         (median) med_pti=pti ///
         (mean) mean_hhinc=hhincome mean_rent=rent mean_valueh=valueh ///
         (sum) pop=perwt ///
         [aw=perwt], by(met2013 year)

* Drop tiny MSAs
drop if pop < 50000

* Average across years for more stable MSA-level estimates
collapse (mean) fert_rate own_rate med_pti mean_hhinc mean_rent mean_valueh ///
         (mean) pop, by(met2013)

* Compute rent-to-income
gen rti = mean_rent * 12 / mean_hhinc

* Label
label var fert_rate "Share with child < 5"
label var own_rate "Ownership rate"
label var med_pti "Median price-to-income"
label var rti "Rent-to-income ratio"

* Summary stats
di as text ""
di as text "--- MSA-level summary statistics ---"
sum fert_rate own_rate med_pti rti, detail

* Correlations
di as text ""
di as text "--- Correlations ---"
pwcorr fert_rate own_rate med_pti rti, sig

* Scatter 1: Fertility vs Ownership
twoway (scatter fert_rate own_rate [aw=pop], msymbol(Oh) mcolor(navy%50)) ///
       (lfit fert_rate own_rate [aw=pop], lcolor(cranberry) lwidth(medthick)), ///
       ytitle("Fertility rate (share with child < 5)") ///
       xtitle("Homeownership rate") ///
       title("Fertility vs Ownership across MSAs") ///
       subtitle("Ages 22-45, ACS") ///
       legend(off) ///
       note("Bubble size proportional to population. Weighted OLS fit line.")
graph export "$out/scatter_fertility_ownership.pdf", replace
graph export "$out/scatter_fertility_ownership.png", replace width(1200)

* Scatter 2: Fertility vs Price-to-income
twoway (scatter fert_rate med_pti [aw=pop], msymbol(Oh) mcolor(navy%50)) ///
       (lfit fert_rate med_pti [aw=pop], lcolor(cranberry) lwidth(medthick)), ///
       ytitle("Fertility rate (share with child < 5)") ///
       xtitle("Median price-to-income ratio") ///
       title("Fertility vs Housing Costs across MSAs") ///
       subtitle("Ages 22-45, ACS") ///
       legend(off) ///
       note("Bubble size proportional to population. Weighted OLS fit line.")
graph export "$out/scatter_fertility_pricetoincome.pdf", replace
graph export "$out/scatter_fertility_pricetoincome.png", replace width(1200)

* Scatter 3: Fertility vs Rent-to-income
twoway (scatter fert_rate rti [aw=pop], msymbol(Oh) mcolor(navy%50)) ///
       (lfit fert_rate rti [aw=pop], lcolor(cranberry) lwidth(medthick)), ///
       ytitle("Fertility rate (share with child < 5)") ///
       xtitle("Rent-to-income ratio") ///
       title("Fertility vs Rental Costs across MSAs") ///
       subtitle("Ages 22-45, ACS") ///
       legend(off) ///
       note("Bubble size proportional to population. Weighted OLS fit line.")
graph export "$out/scatter_fertility_renttoincome.pdf", replace
graph export "$out/scatter_fertility_renttoincome.png", replace width(1200)

* Scatter 4: Ownership vs Price-to-income (sanity check)
twoway (scatter own_rate med_pti [aw=pop], msymbol(Oh) mcolor(navy%50)) ///
       (lfit own_rate med_pti [aw=pop], lcolor(cranberry) lwidth(medthick)), ///
       ytitle("Homeownership rate") ///
       xtitle("Median price-to-income ratio") ///
       title("Ownership vs Housing Costs across MSAs") ///
       subtitle("Ages 22-45, ACS") ///
       legend(off)
graph export "$out/scatter_ownership_pricetoincome.pdf", replace
graph export "$out/scatter_ownership_pricetoincome.png", replace width(1200)

* Regression: fertility on ownership and price-to-income
di as text ""
di as text "--- Regressions: MSA-level fertility ---"
reg fert_rate own_rate [aw=pop], r
reg fert_rate med_pti [aw=pop], r
reg fert_rate own_rate med_pti [aw=pop], r
reg fert_rate rti [aw=pop], r

restore


********************************************************************************
* ANALYSIS 3: Within-metro vs across-metro moves
********************************************************************************

di as text ""
di as text "============================================="
di as text "ANALYSIS 3: Within vs across metro moves"
di as text "============================================="

preserve

* Focus on movers only
keep if moved1y == 1

* We need non-missing metro codes for classification
keep if met2013 > 0 & migmet131 > 0

* Classify move type
gen move_type = .
replace move_type = 1 if same_msa == 1   // Within-MSA move
replace move_type = 2 if diff_msa == 1   // Across-MSA move
label define move_type_lbl 1 "Within MSA" 2 "Across MSA"
label values move_type move_type_lbl

* Overall: what share of moves are within vs across MSA?
di as text ""
di as text "--- All movers: within vs across MSA ---"
tab move_type [aw=perwt]

* By parenthood status
di as text ""
di as text "--- By new parent status ---"
tab move_type newparent [aw=perwt], row

* Statistical test
di as text ""
di as text "--- Regression: P(across-MSA move) on new parent ---"
reg diff_msa newparent i.age i.year i.sex [pw=perwt], r
eststo reg_across: reg diff_msa newparent i.age i.year i.sex i.race [pw=perwt], r

* Now the key question: among within-MSA movers who are new parents,
* what's happening with principal city status?
di as text ""
di as text "--- Within-MSA movers: principal city status ---"
di as text "Among within-MSA new parents vs non-parents:"
tab in_principal newparent [aw=perwt] if same_msa == 1, col

* For within-MSA movers: are new parents less likely to be in principal city?
di as text ""
di as text "--- Regression: P(in principal city) among within-MSA movers ---"
reg in_principal newparent i.age i.year i.sex [pw=perwt] if same_msa == 1, r
eststo reg_principal: reg in_principal newparent i.age i.year i.sex i.race [pw=perwt] if same_msa == 1, r

* Summary table
di as text ""
di as text "--- Summary: Move patterns by parent status ---"
di as text "(weighted shares)"

* Non-parents
sum same_msa [aw=perwt] if newparent == 0
local np_within = r(mean)
sum diff_msa [aw=perwt] if newparent == 0
local np_across = r(mean)

* New parents
sum same_msa [aw=perwt] if newparent == 1
local p_within = r(mean)
sum diff_msa [aw=perwt] if newparent == 1
local p_across = r(mean)

di as text ""
di as text "                    Non-Parents    New Parents"
di as text "Within MSA:     " %8.3f `np_within' "       " %8.3f `p_within'
di as text "Across MSA:     " %8.3f `np_across' "       " %8.3f `p_across'

* Bar chart
graph bar (mean) same_msa diff_msa [aw=perwt], ///
      over(newparent, relabel(1 "Non-parents" 2 "New parents")) ///
      bar(1, color(navy)) bar(2, color(cranberry)) ///
      legend(label(1 "Within MSA") label(2 "Across MSA")) ///
      ytitle("Share of movers") ///
      title("Move Type by Parent Status") ///
      subtitle("ACS movers, ages 22-45") ///
      note("Movers with identifiable current and previous MSA")
graph export "$out/bar_move_type_by_parent.pdf", replace
graph export "$out/bar_move_type_by_parent.png", replace width(1200)

* Export regression table
esttab reg_across reg_principal using "$out/reg_move_type.tex", ///
       keep(newparent) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
       label replace ///
       mtitle("P(Across MSA)" "P(In Principal City | Within MSA)") ///
       title("Move Patterns and Parenthood")

restore


********************************************************************************
* BONUS: Age-split for Analysis 1
* Is the fertility-ownership gradient steeper for young households?
********************************************************************************

di as text ""
di as text "============================================="
di as text "BONUS: Age-split fertility-ownership gradient"
di as text "============================================="

preserve

keep if met2013 > 0

gen young = (age <= 32)

* Collapse to MSA × age-group
collapse (mean) fert_rate=has_youngchild own_rate=owner ///
         (median) med_pti=pti ///
         (sum) pop=perwt ///
         [aw=perwt], by(met2013 young)

drop if pop < 20000

label define young_lbl 0 "Age 33-45" 1 "Age 22-32"
label values young young_lbl

* Scatter by age group
twoway (scatter fert_rate own_rate [aw=pop] if young == 1, msymbol(Oh) mcolor(navy%50)) ///
       (lfit fert_rate own_rate [aw=pop] if young == 1, lcolor(navy) lwidth(medthick)) ///
       (scatter fert_rate own_rate [aw=pop] if young == 0, msymbol(Th) mcolor(cranberry%50)) ///
       (lfit fert_rate own_rate [aw=pop] if young == 0, lcolor(cranberry) lwidth(medthick)), ///
       ytitle("Fertility rate") ///
       xtitle("Homeownership rate") ///
       title("Fertility vs Ownership by Age Group") ///
       legend(order(2 "Age 22-32" 4 "Age 33-45")) ///
       note("Bubble size proportional to population")
graph export "$out/scatter_fertility_ownership_byage.pdf", replace
graph export "$out/scatter_fertility_ownership_byage.png", replace width(1200)

* Regressions by age group
di as text ""
di as text "--- Young (22-32) ---"
reg fert_rate own_rate [aw=pop] if young == 1, r
di as text ""
di as text "--- Older (33-45) ---"
reg fert_rate own_rate [aw=pop] if young == 0, r

restore

di as text ""
di as text "============================================="
di as text "DONE"
di as text "============================================="
