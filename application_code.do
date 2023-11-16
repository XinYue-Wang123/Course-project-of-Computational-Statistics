**********  Stata code for application of Lasso   **********



use "Z:\Desktop\code.dta", clear
* keep hhid region rural pop a2000a E_indirect_energy E_direct_energy total_income_w edu_mean health_mean age_mean house_area

**********  the LASSO regression  **********
ssc install asdoc
ssc install lassopack
ssc install estout

* number of people in househould *
drop pop
gen pop = a2000a + 1

* dependent variable*
gen per_direct    = E_direct_energy / pop
gen ln_per_direct = log(per_direct)

* independent variables*
gen per_income     = total_income_w / pop
gen ln_per_income  = log(per_income)
gen ln_pop         = log(pop)
gen ln_edu_mean    = log(edu_mean)
gen ln_health_mean = log(health_mean)
gen ln_age_mean    = log(age_mean)
gen ln_house_area  = log(house_area)

* descriptive statistics *
asdoc sum ln_per_direct ln_per_income ln_pop ln_edu_mean ln_health_mean ln_age_mean ln_house_area 

* plot direct carbon emissions pre person for Urban and Rural *
bysort region rural:egen direct_emission = mean(per_direct)
graph bar direct_emission if region==1 , over(rural) scheme(s2mono)

* standardisation *
egen direct_s = std(per_direct)
egen income_s = std(per_income)
egen pop_s    = std(pop)
egen edu_s    = std(edu_mean)
egen health_s = std(health_mean)
egen age_s    = std(age_mean)
egen area_s   = std(house_area)

*** LASSO on direct carbon emissions  ***
* Urban *
lasso linear direct_s income_s pop_s edu_s health_s age_s area_s if rural==0, selection(cv, alllambdas) stop(0) rseed(12345) nolog
estimates store cv
cvplot
lassocoef, display(coef,penalized) sort(coef,penalized)

lasso2  direct_s income_s pop_s edu_s health_s age_s area_s  if rural==0, plotpath(lambda)
cvlasso direct_s income_s pop_s edu_s health_s age_s area_s  if rural==0, lopt seed(123) plotcv
esttab, star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) se(%6.4f) mtitle(`models') r2 sca(r2_w)

* Rural *
lasso linear direct_s income_s pop_s edu_s health_s age_s area_s if rural==1, selection(cv, alllambdas) stop(0) rseed(12345) nolog
estimates store cv
cvplot
lassocoef, display(coef,penalized) sort(coef,penalized)

lasso2  direct_s income_s pop_s edu_s health_s age_s area_s  if rural==1, plotpath(lambda)
cvlasso direct_s income_s pop_s edu_s health_s age_s area_s  if rural==1 , lopt seed(123) plotcv 
esttab, star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) se(%6.4f) mtitle(`models') r2 sca(r2_w)
