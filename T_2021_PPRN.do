********************************************************************************
**  Regression Discontinuity Designs 
** Prepared for the Public Policy Research Network Methods Masterclass
**  March 12, 2021
** Author: Rocio Titiunik, Matias Cattaneo, Sebastian Calonico
********************************************************************************
** RDROBUST:  		net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
** RDLOCRAND: 	net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
** RDDENSITY: 		net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
** RDPOWER:   		net install rdpower, from(https://raw.githubusercontent.com/rdpackages/rdpower/master/stata) replace
** LPDENSITY :  	net install lpdensity, from(https://sites.google.com/site/nppackages/lpdensity/stata) replace
********************************************************************************

********************************************************************************
clear all
set more off
set linesize 200

******************************************************************************
**  U.S. Senate incumbency advantage --- Cattaneo, Frandsen, Titiunik (2015, Journal of Causal Inference)
********************************************************************************
use senate, clear
global x demmv
global y demvoteshfor2
global c = 0
gen T=.
replace T=0 if $x<0 & $x!=.
replace T=1 if $x>=0 & $x!=.
label var T "Democratic Win at t"
order $x $y T

********************************************************************************
** Summary Stats & Diff-in-means
********************************************************************************
sum $y $x demvoteshlag1 demvoteshlag2 dopen
ttest $y, by(T)

********************************************************************************
** Part 2:  Graphical Illustration of RD Models
********************************************************************************

** Quick RD plots (default: mimicking variance)
rdplot $y $x

** Evenly-spaced
rdplot $y $x, c($c) binselect(es)

** Quantile-spaced
rdplot $y $x, c($c) binselect(qs)

** Global 2nd order polynomial
rdplot $y $x, c($c) p(2) 

** Select manually number of bins
rdplot $y $x, c($c) nbins(10)  

** Add confidence intervals
rdplot $y $x, c($c) nbins(10) ci(95) 

** Add confidence intervals w/shade
rdplot $y $x, c($c) nbins(10) ci(95) shade

** Generate variables
rdplot $y $x, c($c) p(2) ci(95) shade genvars
drop rdplot_*

** Use support option
rdplot $y $x, nbins(20) support(-100 100)

** Stored output
ereturn list			    


********************************************************************************
** Part 3:  RD Estimation and Inference
********************************************************************************	    
			    
********************************************************************************
** Part 3a:  Continuity-Based RD Estimation and Inference
********************************************************************************

** Equivalence between least-squares and local polynomial estimation
* Output from rdrobust, uniform kernel
rdrobust $y $x, h(10) kernel(uniform)

* Output from regressions on both sides
reg $y $x if $x<0 & $x>=-10
matrix coef_left=e(b)
local intercept_left=coef_left[1,2]
reg $y $x if $x>=0 & $x<=10
matrix coef_right=e(b)
local intercept_right=coef_right[1,2]
local difference=`intercept_right'-`intercept_left'
display "The RD estimator is `difference'"

* Output from joint regression
gen T_x=$x*T
reg $y $x T T_x if abs($x)<=10

* To get the same standard errors, use vce() option
reg $y $x T T_x if abs($x)<=10, robust
rdrobust $y $x, h(10) kernel(uniform) vce(hc1)

** Local polynomial estimation and inference choosing bandwidth by hand (not recommended)
rdrobust $y $x,  h(10)  

** Choosing MSE-optimal bandwidth (using rdrobust default options)
rdbwselect $y $x, kernel(triangular) p(1) bwselect(mserd) 

** Allowing different bandwidths
rdbwselect $y $x, kernel(triangular) p(1) bwselect(msetwo) 

*  Local polynomial estimation and robust inference using MSE-optimal bandwidth 
rdrobust $y $x
ereturn list

* showing the associated rdplot)
rdrobust $y $x, p(1) kernel(triangular) bwselect(mserd)
local bandwidth=e(h_l)
rdplot $y $x if abs($x)<=`bandwidth', p(1) h(`bandwidth') kernel(triangular)

* CER bandwidth 
rdrobust $y $x, kernel(triangular) p(1) bwselect(cerrd)  

* Using rdrobust default options and showing all the output)
rdrobust $y $x, kernel(triangular) p(1) bwselect(mserd) all

* rdrobust with covariates 
global covariates "presdemvoteshlag1 demvoteshlag1 demvoteshlag2 demwinprv1 demwinprv2 dmidterm dpresdem dopen"
rdbwselect $y $x, covs($covariates) p(1) kernel(triangular) bwselect(mserd) scaleregul(1)
		    
********************************************************************************
** Part 3b: Local Randomization RD Estimation and Inference
********************************************************************************
* Step 1: select window
*  Select the window with covariates: Using rdwinselect with the wobs Option
global covs "presdemvoteshlag1 demvoteshlag1 demvoteshlag2 demwinprv1 demwinprv2 dmidterm dpresdem dopen"
rdwinselect $x $covs, seed(50) wobs(2) nwindows(15) 

* Using rdwinselect with the wobs, nwindows and plot Options
global covs "presdemvoteshlag1 demvoteshlag1 demvoteshlag2 demwinprv1 demwinprv2 dmidterm dpresdem dopen"
rdwinselect $x $covs, seed(50) wobs(2) nwindows(3) plot 

*  Using rdwinselect with the wstep and nwindows Options
global covs "presdemvoteshlag1 demvoteshlag1 demvoteshlag2 demwinprv1 demwinprv2 dmidterm dpresdem dopen"
rdwinselect $x $covs, seed(50) wstep(0.1) nwindows(3) 

* Step 2: Estimation and inference in chosen window
global window = 0.75

* Using rdrandinf in the Ad-hoc Window
rdrandinf $y $x, wl(-$window) wr($window) seed(50)

* Using rdrandinf in the Ad-hoc Window with the Bernoulli Option
gen bern_prob=1/2 if abs($x)<=0.75
rdrandinf $y $x, wl(-$window) wr($window) seed(50) bernoulli(bern_prob)

*  Using rdrandinf in the Ad-hoc Window with the Kolmogorov-Smirnov Statistic
rdrandinf $y $x, wl(-$window) wr($window) seed(50) statistic(ksmirnov)

*  Using rdrandinf in window with the Confidence Interval Option --- use large increment because it takes long-
rdrandinf $y $x, wl(-$window) wr($window)  seed(50) ci(0.05 -10(0.5)10)

* Using rdrandinf with the Alternative Power Option
global covs "presdemvoteshlag1 demvoteshlag1 demvoteshlag2 demwinprv1 demwinprv2 dmidterm dpresdem dopen"
rdrandinf $y $x, covariates($covs) seed(50) d(7.416) 
			    
********************************************************************************
** Part 4: RD Falsification Analysis
********************************************************************************
* -------- Density Tests ---------------------*
* Continuity-based : local polynomial density estimation
rddensity $x, plot

* Local randomization: binomial test
rdwinselect $x, wmin($window) nwindows(1) 
bitesti 39 15 1/2	    

* ---------------- Covariate Balance and Placebo Outcomes Tests -----------------------------------*
* Continuity-based:  Using rdrobust on demvoteshlag1
rdrobust demvoteshlag1 $x
* Local randomization: Using demvoteshlag1
ttest demvoteshlag1 if abs($x)<=$window, by(T)
rdrandinf  demvoteshlag1 $x, wl(-$window) wr($window)
