** Try to generalize the cleaner to make it apt for all files

clear

set more off



**setting up overleaf
global path "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs"
global tables "$path/Tables_ACS"
global graphs "$path/Graphs_ACS"
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg" ) lcolor("$bg" )) plotregion(fcolor("$bg" ) lcolor("$bg" ))


** big faile is main
use /Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta, replace


replace rent=. if rent==0
gen logrent=log(rent)
gen logrthhi=log(rent/hhincome)
gen logrtpi=log(rent/inctot)
gen loghhi=log(hhi)
gen logpinc=log(inctot)
label variable logrent "log(rent)"

drop if age<18

*replace Population2010=log(Population2010)
*replace densityw2010=log(densityw2010)

 drop if hhincome==9999999
drop if inctot==9999999

 label define agelab  0 "Under 40" 1 "40 to 60" 2 "Over 60"
gen agegroup=.
replace agegroup= 0 if age<40
replace agegroup= 1 if age>=40 & age<=60 
replace agegroup =2 if age>60
label values agegroup agelab
 tabulate agegroup, summarize()

  *now with extended categories
  gen agegroup_decades=.
replace agegroup_decades= 0 if age<30
replace agegroup_decades= 1 if age>=30 & age<=40 
replace agegroup_decades= 2 if age>=40 & age<=50 
replace agegroup_decades= 3 if age>=50 & age<=60 
replace agegroup_decades= 4 if age>=60 & age<=70 
replace agegroup_decades= 5 if age>70

 label define agegroup_decades_lab  0 "Under 30" 1 "30 to 40" 2 "40 to 50" 3 "50 to 60" 4 "60 to 70" 5 "Older than 70"
 label values agegroup_decades agegroup_decades_lab
 
 
 replace famsize=10 if famsize >10 
  
 
   ** export delimited using "/Users/tommasodesanto/Desktop/ACS_microdata/extract_10/export_main.csv", replace
   global y_choices rent logrent logrthhi logrtpi
    global sizevar  Population2010 densityw2010
	global incvar loghhi logpinc
	global demos_noinc sex educd empstatd 
	
	
	
	
	gen hadfirstchild=( eldch==0 ) 
	**review the above
	gen movedthisyear=.
	replace movedthisyear=0 if migrate1d==10
replace movedthisyear=1 if (migrate1d!=10 & migrate1d!=.)



	gen movedinthelast5=.
	replace movedinthelast5=0 if migrate1d==10
replace movedinthelast5=1 if (migrate5d!=10 & migrate5d!=.)


gen littlechild=(eldch<6)

label variable movedthisyear "Moved this year"
label variable movedinthelast5 "Moved in last 5 years"



**just a test
gen inc_weighted = inctot * perwt
gen weight_total = perwt

bysort metarea year: gen sum_inc_weighted = sum(inc_weighted)
bysort metarea year: gen sum_weight_total = sum(weight_total)

gen avg_income_metro = sum_inc_weighted / sum_weight_total

label variable hadfirstchild "Had First Child"
label variable littlechild "Child less than 5"

*preserve
*drop if eldch>15

**for 1 year migration


local controls "No"
eststo reg1: qui reg movedthisyear hadfirstchild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg movedthisyear hadfirstchild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
 esttab reg* using	"$tables/moved_andfirstchild.tex", drop(*.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) 		 stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")		
estimates drop _all

local controls "No"
eststo reg1: qui reg movedthisyear littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg movedthisyear littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg movedinthelast5 littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg movedinthelast5 littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
 esttab reg* using	"$tables/moved_andfirstchild_robust.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")	 			
estimates drop _all

local controls "No"
eststo reg1: qui reg rooms hadfirstchild##movedthisyear i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg rooms hadfirstchild##movedthisyear i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
 esttab reg* using	"$tables/rooms_firstchildandmoves.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) stats(Controls, labels("Controls"))  addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment status, sex, and county.")	 					
estimates drop _all

gen hadfirstchild_var= 1 if eldch==0
replace hadfirstchild_var=0 if nchild==0
gen littlechild_var= 1 if eldch<6
replace littlechild_var=0 if nchild==0

gen littlechild_var2= 1 if eldch<4 & eldch!=.
replace littlechild_var2=0 if nchild==0

label var hadfirstchild_var "Had First Child"
local controls "No"
eststo reg1: qui reg rooms hadfirstchild_var##movedthisyear i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg rooms hadfirstchild_var##movedthisyear i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
 esttab reg* using	"$tables/rooms_firstchildandmoves_var.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) stats(Controls, labels("Controls"))    addnotes("Only non-parents and new parents. Controls: income, education, employment, sex, and county.")	 					
estimates drop _all


local controls "No"
eststo reg1: qui reg rooms littlechild##movedthisyear i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg rooms littlechild##movedthisyear i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
 esttab reg* using	"$tables/rooms_littlechildandmoves.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")	 		 			
estimates drop _all


eststo reg1: qui reg rooms hadfirstchild##movedinthelast5 i.age i.year [aweight=perwt],r
eststo reg2: qui reg rooms hadfirstchild##movedinthelast5 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r

eststo reg3: qui reg rooms littlechild##movedinthelast5 i.age i.year [aweight=perwt],r
eststo reg4: qui reg rooms littlechild##movedinthelast5 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
 esttab reg* using	"$tables/rooms_littlechildandmoveslast5.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) 			
estimates drop _all

**try plot moved by age
* Step 1: Create the age group variable
gen agegrp_eldch = .
replace agegrp_eldch = eldch if eldch < 10
replace agegrp_eldch = 10 if eldch >= 10
replace agegrp_eldch = -1 if eldch==99


* Step 2: Label the age group categories
label define agegrp_eldch_lbl -1 "No Child" 0 "0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5"  6 "6" 7 "7" 8 "8" 9 "9" 10 "10+"
label values agegrp_eldch agegrp_eldch_lbl
	local graphregion  $graphregion_slides
* Step 3: Plot the proportion of people who moved this year by age group
graph bar (mean) movedthisyear [aweight=perwt], over(agegrp_eldch) ascategory b1title("Age of eldest child") ytitle("Proportion Moved")  ylabel(0(.1)0.3, format(%2.1f))  legend(off) `graphregion'  
graph export "$graphs/moved_andageofeldestchild.pdf",replace

* Step 3: Plot the proportion of people who moved this year by age group
graph bar (mean) rooms [aweight=perwt], over(agegrp_eldch) ascategory  ytitle("Percent Moved") ylabel(0(10)100) blabel(bar, format(%2.0f)) `graphregion' 

* Create a conditional version with controls
* Create positive-valued version of agegrp_eldch
gen agegrp_eldch_pos = agegrp_eldch + 2  // This shifts everything up by 2, so -1 becomes 1, 0 becomes 2, etc.

* Create labels for the new variable
label define agegrp_eldch_pos_lbl 1 "No Child" 2 "0" 3 "1" 4 "2" 5 "3" 6 "4" 7 "5" 8 "6" 9 "7" 10 "8" 11 "9" 12 "10+"
label values agegrp_eldch_pos agegrp_eldch_pos_lbl

* Run regression with the new positive-valued variable
quietly reg movedthisyear i.agegrp_eldch_pos i.age i.year inctot i.educd i.empstatd i.sex [aweight=perwt], r

* Use margins to get adjusted predictions for each age group category
margins agegrp_eldch_pos, atmeans

* Create the bar graph
marginsplot, noci recast(bar) xtitle("Age of eldest child") ytitle("Adjusted probability of moving") ///
  title("Probability of moving by age of eldest child", size(medium)) ///
  subtitle("Adjusted for demographics, income, education, and location", size(small)) ///
  ylabel(0(.05).25, format(%3.2f)) `graphregion' plot1opts(barw(0.8)) ///
  xlabel(1 "No Child" 2 "0" 3 "1" 4 "2" 5 "3" 6 "4" 7 "5" 8 "6" 9 "7" 10 "8" 11 "9" 12 "10+") ///
  xsize(8) ysize(5) ///
  saving("$graphs/adjusted_margins.gph", replace)

* Add bar labels
quietly margins agegrp_eldch_pos, atmeans
matrix b = r(b)'
local n = rowsof(b)
* Create a temporary dataset from the matrix
preserve
clear
svmat b
gen group = _n
label define group_lbl 1 "No Child" 2 "0" 3 "1" 4 "2" 5 "3" 6 "4" 7 "5" 8 "6" 9 "7" 10 "8" 11 "9" 12 "10+"
label values group group_lbl

* Create the bar graph with the background color settings
global bg "250 250 250"
graph bar b1, over(group, label(angle(45))) ///
    ytitle("Adjusted probability of moving") ///
    title("Predicted Probability of moving by age of eldest child", size(medium)) ///
    subtitle("Controlling for demographics, income and education", size(small)) ///
    ylabel(0(.05).25, format(%3.2f)) ///
    graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg")) ///
    blabel(bar, format(%4.3f))

* Save the graph
graph export "$graphs/moved_andageofeldestchild_adjusted.pdf", replace
restore



label define moved_lbl 0 "Did not Move" 1 "Moved"
label values movedthisyear moved_lbl
graph bar (mean) rooms [aweight=perwt], over(agegrp_eldch, label(angle(45))) over(movedthisyear)   ascategory ytitle("Mean Number of Rooms") legend(title("Mover Status")) bargap(15)  `graphregion' 
graph export "$graphs/numberofrooms_moversornot.pdf",replace


preserve
collapse (mean) mean_rooms=rooms [aweight=perwt], by(agegrp_eldch movedthisyear)
twoway (line mean_rooms agegrp_eldch if movedthisyear==0, sort lcolor(blue) lpattern(solid))   (line mean_rooms agegrp_eldch if movedthisyear==1, sort lcolor(red) lpattern(dash)),    legend(order(1 "Non-Movers" 2 "Movers") title("Mover Status"))  ytitle("Mean Number of Rooms")  xtitle("Age Group of Eldest Child")   xlabel(0(1)10, angle(45))

restore


* Strict definitions
gen byte moved_excl_puma  = inlist(migrate1d,24,31,32,40) if !missing(migrate1d)  // 24=between PUMAs (within state)
gen byte moved_excl_state = inlist(migrate1d,31,32,40)     if !missing(migrate1d)  // across states or abroad

* (Optional lenient variant that counts "unknown within-state" as mover)
gen byte moved_excl_puma_lenient = inlist(migrate1d,24,25,31,32,40) if !missing(migrate1d)
* Strict definitions
gen byte moved_excl_puma  = inlist(migrate1d,24,31,32,40) if !missing(migrate1d)  // 24=between PUMAs (within state)
gen byte moved_excl_state = inlist(migrate1d,31,32,40)     if !missing(migrate1d)  // across states or abroad

* (Optional lenient variant that counts "unknown within-state" as mover)
gen byte moved_excl_puma_lenient = inlist(migrate1d,24,25,31,32,40) if !missing(migrate1d)


* --- Predicted prob: excl. within-PUMA ---
qui reg moved_excl_puma i.agegrp_eldch_pos c.age i.year c.hhincome i.educd i.empstatd i.sex [pw=perwt], r
margins agegrp_eldch_pos, atmeans
marginsplot, noci recast(bar) `graphregion'
graph export "$graphs/moved_age_excl_puma.pdf", replace

* --- Predicted prob: excl. within-state (also excludes within-PUMA) ---
qui reg moved_excl_state i.agegrp_eldch_pos c.age i.year c.hhincome i.educd i.empstatd i.sex [pw=perwt], r
margins agegrp_eldch_pos, atmeans
marginsplot, noci recast(bar) `graphregion'
graph export "$graphs/moved_age_excl_state.pdf", replace

	



preserve
collapse (mean) mean_rooms=rooms [aweight=perwt], by(agegrp_eldch migrate1d)

* Step 1: Calculate base_mean_rooms for each migrate1d where agegrp_eldch == -1
bysort migrate1d: egen base_mean_rooms_temp = mean(mean_rooms) if agegrp_eldch == -1

* Step 2: Propagate base_mean_rooms to all observations within the same migrate1d
 bysort migrate1d: egen base_mean_rooms= max(base_mean_rooms_temp)

* Step 3: Drop observations where base_mean_rooms is missing (no base value for that migrate1d)
drop if missing(base_mean_rooms)

* Step 4: Calculate normalized mean rooms
gen normalized_mean_rooms = mean_rooms / base_mean_rooms



**still have to fix this but right direction
twoway (line normalized_mean_rooms agegrp_eldch if migrate1d==10, sort lcolor(blue) lpattern(solid)) (line normalized_mean_rooms agegrp_eldch if migrate1d==23, sort lcolor(red) lpattern(dash)) (line normalized_mean_rooms agegrp_eldch if migrate1d==24, sort lcolor(green) lpattern(dot)) (line normalized_mean_rooms agegrp_eldch if migrate1d==25, sort lcolor(orange) lpattern(longdash)) (line normalized_mean_rooms agegrp_eldch if migrate1d==31, sort lcolor(purple) lpattern(shortdash)) (line normalized_mean_rooms agegrp_eldch if migrate1d==32, sort lcolor(brown) lpattern(dash_dot)) (line normalized_mean_rooms agegrp_eldch if migrate1d==40, sort lcolor(black) lpattern(solid)), legend(order(1 "Same house" 2 "Different house, moved within state, wi" 3 "Different house, moved within state, be" 4 "Different house, unknown within state" 5 "Moved between contiguous states" 6 "Moved between non-contiguous states" 7 "Abroad one year ago") position(6) ring(1) col(1) width(8) size(small)) title("Migration Status") ytitle("Normalized Mean Number of Rooms") xtitle("Age Group of Eldest Child") xlabel(0(1)10, angle(45)) ylabel(, angle(0) format(%9.2f)) graphregion(color(gs14)) plotregion(color(white))

restore

graph export "$graphs/rooms_byeldestchildage.pdf",replace

preserve
drop if metro==4

gen metroyn=.
replace metroyn=1 if metro== 2 | metro==4
replace metroyn=0 if metro==1| metro==3
gen metro1yyn=.
replace metro1yyn=1 if migtype1== 2 | migtype1==3
replace metro1yyn=0 if migtype1==1| migtype1==4 |migtype1==5

gen went_metro=(metroyn==1 & metro1yyn==0)
replace went_metro=. if (metroyn==. |metro1yyn==. )
gen left_metro= (metroyn==0 & metro1yyn==1)
replace left_metro=. if (metroyn==. |metro1yyn==. )



gen metroyn_variant=.
replace metroyn_variant=1 if metro==2
replace metroyn=0 if metro==1 | metro==3

gen metro1yyn_variant=.
replace metro1yyn_variant=1  if migtype1==3
replace metro1yyn_variant=0 if migtype1==1| migtype1==4 |migtype1==5

gen left_metro_variant= (metroyn_variant==0 & metro1yyn_variant==1)
replace left_metro_variant=. if (metroyn_variant==. |metro1yyn_variant==. )

gen went_metro_variant=(metroyn_variant==1 & metro1yyn_variant==0)
replace left_metro_variant=. if (metroyn_variant==. |metroyn_variant==. )




gen metro5yyn=.
replace metro5yyn=1 if migtype5== 2 | migtype5==3
replace metro5yyn=0 if migtype5==1| migtype5==4 |migtype5==5

gen went_metro_inlast5=(metroyn==1 & metro5yyn==0)
replace went_metro_inlast5=. if (metroyn==. |metro5yyn==. )
gen left_metro_inlast5= (metroyn==0 & metro5yyn==1)
replace left_metro_inlast5=. if (metroyn==. |metro5yyn==. )


label variable metroyn "Lives in"
label variable left_metro "Left"
label variable went_metro "Moved to"

label variable metroyn_variant "Lives in"
label variable left_metro_variant "Left"
label variable went_metro_variant "Moved to"



**for 1 year migration
local controls "No"
eststo reg1: qui reg metroyn hadfirstchild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg metroyn hadfirstchild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro hadfirstchild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro hadfirstchild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro hadfirstchild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro hadfirstchild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"


 esttab reg* using	"$tables/metro_status_andchildren.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip)   scalars("sum Sum" "k Metro Area Status \\ \hline")   replace label se(3) 	stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")			
estimates drop _all




local controls "No"
eststo reg1: qui reg metroyn littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"

local controls "Yes" 
eststo reg2: qui reg metroyn littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"


 esttab reg* using	"$tables/metro_status_andchildren_littlechild.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3)	stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")	
estimates drop _all



local controls "No"
eststo reg3: qui reg went_metro_inlast5 littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro_inlast5 littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro_inlast5 littlechild i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes"
eststo reg6: qui reg left_metro_inlast5 littlechild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"

 esttab reg* using	"$tables/metro_status_andchildren_littlechild5y.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3) 	stats(Controls, labels("Controls"))   addnotes("Controls include income, education, employment status, sex, and county FIPS.")		
estimates drop _all


**for 5 year migration


**variant
**for 1 year migration
local controls "No"
eststo reg1: qui reg metroyn hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg metroyn hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"

 esttab reg2 reg4 reg6 using	"$tables/metro_status_andchildren_var.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip)   scalars("sum Sum" "k Metro Area Status \\ \hline")   replace label se(3) 	stats(Controls, labels("Controls"))     addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment, sex, and county.")
 
 
 
estimates drop _all








*variant definition of metro
**variant
**for 1 year migration
local controls "No"
eststo reg1: qui reg metroyn_variant hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg metroyn_variant hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro_variant hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro_variant hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro_variant hadfirstchild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro_variant hadfirstchild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"




 esttab reg2 reg4 reg6 using	"$tables/metro_status_var_andchildren_var.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip)   scalars("sum Sum" "k Metro Area Status \\ \hline")   replace label se(3) 	stats(Controls, labels("Controls"))     addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment, sex, and county.")
 


local controls "Yes" 
eststo reg2: qui reg metroyn littlechild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro littlechild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro littlechild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro littlechild_var i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro littlechild_var i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"


 esttab reg2 reg4 reg6 using	"$tables/metro_status_andchildren_littlechild_var.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3)	stats(Controls, labels("Controls"))      addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment status, sex, and county.")
estimates drop _all

label var littlechild_var2 "New Parent"

local controls "Yes" 
eststo reg2: qui reg metroyn littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro littlechild_var2 i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro littlechild_var2 i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"

			 
			 
 esttab reg2 reg4 reg6 using	"$tables/metro_status_andchildren_littlechild_var2.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip)   scalars("sum Sum" "k Metro Area Status \\ \hline")   replace label se(3) 	stats(Controls, labels("Controls"))     addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment, sex, and county.")
 
 
estimates drop _all


*variant definition of metro
**variant
**for 1 year migration
local controls "No"
eststo reg1: qui reg metroyn_variant littlechild_var2 i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg2: qui reg metroyn_variant littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg3: qui reg went_metro_variant littlechild_var2 i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg4: qui reg went_metro_variant littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"
local controls "No"
eststo reg5: qui reg left_metro_variant littlechild_var2 i.age i.year [aweight=perwt],r
estadd local Controls "`controls'"
local controls "Yes" 
eststo reg6: qui reg left_metro_variant littlechild_var2 i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
estadd local Controls "`controls'"




 esttab reg2 reg4 reg6 using	"$tables/metro_status_var_andchildren_var_littlechildvar2.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip)   scalars("sum Sum" "k Metro Area Status \\ \hline")   replace label se(3) 	stats(Controls, labels("Controls"))     addnotes("Only non-parents and new parents" ///
             "Controls include income, education, employment, sex, and county.")


*different definition of metro




gen haschild=1 if nchild>0
replace haschild=0 if nchild==0


**new export for notes
preserve
**select sample
drop if age<22 
replace nchild=6 if nchild>6
label define nchild_lbl 0 "0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "more than 5"
label values nchild nchild_lbl
eststo reg1: qui reg metroyn haschild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r
eststo reg2: qui reg metroyn i.nchild i.age i.year inctot i.educd i.empstatd i.sex i.countyfip [aweight=perwt],r

 esttab reg* using	"$tables/metrostat_nchildren.tex",nobaselevels drop(inctot *.sex *.age *.educd *.empstatd  *.year   *.countyfip) replace label se(3)
 estimates drop _all

 
 
 
 
 
 **baseline means for migration analysis
 * Ensure variable labels are set (as in your original code)

* --- Calculate Means ---
preserve
keep metroyn went_metro left_metro perwt littlechild_var2 nchild
* 1. Specific Sample (Non-parents and New Parents)
*    Adjust the 'if' condition below if 'hadfirstchild_var' is not the right identifier
*    or if the sample is defined differently.
estpost summarize metroyn went_metro left_metro [aweight=perwt] if littlechild_var2==1 | nchild==0
eststo means_specific

* 2. Overall Sample (All observations after initial cleaning)
estpost summarize metroyn went_metro left_metro [aweight=perwt]
eststo means_overall

* --- Generate LaTeX Table ---



* Now run esttab
esttab means_specific means_overall using "$tables/baseline_means_simple.tex", ///
    cells(mean(fmt(%9.3f)))  /* Show mean with 3 decimals */ ///
    stats(N sum_w, fmt(%15.0f) labels("Observations" "Sum of Weights")) /* Explicit format for N and sum_w */ ///
    label /* Use variable labels */ ///
    booktabs /* Use booktabs rules */ ///
    replace /* Overwrite file */ ///
    nonumbers /* Remove (1), (2) column numbers */ ///
    nogaps /* Remove extra space between rows */ ///
    title("Baseline Weighted Means by Population Group") /// /* Add a caption title */ ///
    collabels("Specific Sample" "Overall Sample") /// /* Use simple, single-line column labels */ ///
    addnotes("Means weighted using perwt. Specific Sample = Non-parents and new parents (hadfirstchild\_var non-missing).") // Simplified note with escaped underscoreAdd note, ESCAPED ampersand and underscore

* --- Clean up estimates store ---
estimates drop means_specific means_overall

di "Corrected LaTeX table saved to $tables/baseline_means.tex"

* --- Clean up estimates store ---
estimates drop means_specific means_overall

di "LaTeX table saved to $tables/baseline_means.tex"
