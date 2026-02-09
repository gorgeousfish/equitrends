*! example_empirical_application.do
*! Replication of Dette & Schumann (2024), Section 7: Empirical Illustration
*!
*! Data: Di Tella & Schargrodsky (2004) Buenos Aires crime panel
*!   - 876 city blocks observed monthly (April-December 1994)
*!   - Treatment: 37 blocks received police protection after July terrorist attack
*!   - Outcome: Monthly car theft counts per block
*!
*! This file demonstrates the full workflow of the equitrends package:
*!   1. Data preparation
*!   2. Table 9 replication (minimum equivalence thresholds)
*!   3. Hypothesis testing with a specified threshold
*!   4. Coefficient plots via equivtest_plot
*!   5. Variance estimation options
*!   5b. Control variables (conditional parallel trends)
*!   6. Monte Carlo simulation via equivsim

version 16.0
clear all
set more off


* =============================================================================
* 1. DATA PREPARATION
* =============================================================================

* Load bundled dataset (after installation: equitrends_data, clear)
capture noisily equitrends_data, clear
if _rc != 0 {
    capture use "`c(pwd)'/data/MonthlyPanel.dta", clear
    if _rc != 0 {
        di as error "Cannot find MonthlyPanel.dta. Install the package first."
        exit 601
    }
}

* Keep pre-treatment months: April(4)-July(7). July = base period.
drop if mes == 72 | mes == 73
drop if mes > 7

rename observ ID
rename totrob Y
rename mes period

gen G = (distanci == 0)
label variable G "Treatment (block with Jewish institution)"

xtset ID period


* =============================================================================
* 2. TABLE 9 REPLICATION — Minimum Equivalence Thresholds
* =============================================================================
*
* delta*/tau*/zeta* = smallest threshold at which H0 is rejected (alpha=0.05)

* --- T=1: June only (pretreatment 6 7, base 7) ---

equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(6 7) baseperiod(7) vce(cluster) cluster(ID)
local delta_iu_T1 = e(min_threshold)

equivtest Y, type(max) method(boot) id(ID) group(G) time(period) ///
    pretreatment(6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_boot_T1 = e(min_threshold)

equivtest Y, type(max) method(wild) id(ID) group(G) time(period) ///
    pretreatment(6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_wild_T1 = e(min_threshold)

equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(6 7) baseperiod(7) vce(cluster) cluster(ID)
local tau_T1 = e(min_threshold)

equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(6 7) baseperiod(7) seed(2024)
local zeta_T1 = e(min_threshold)

* --- T=2: May & June (pretreatment 5 6 7, base 7) ---

equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(5 6 7) baseperiod(7) vce(cluster) cluster(ID)
local delta_iu_T2 = e(min_threshold)

equivtest Y, type(max) method(boot) id(ID) group(G) time(period) ///
    pretreatment(5 6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_boot_T2 = e(min_threshold)

equivtest Y, type(max) method(wild) id(ID) group(G) time(period) ///
    pretreatment(5 6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_wild_T2 = e(min_threshold)

equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(5 6 7) baseperiod(7) vce(cluster) cluster(ID)
local tau_T2 = e(min_threshold)

equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(5 6 7) baseperiod(7) seed(2024)
local zeta_T2 = e(min_threshold)

* --- T=3: April-June (pretreatment 4 5 6 7, base 7) ---

equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)
local delta_iu_T3 = e(min_threshold)

equivtest Y, type(max) method(boot) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_boot_T3 = e(min_threshold)

equivtest Y, type(max) method(wild) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(2024) nodots
local delta_wild_T3 = e(min_threshold)

equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)
local tau_T3 = e(min_threshold)

equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) seed(2024)
local zeta_T3 = e(min_threshold)

* --- Display Table 9 ---
display as text ""
display as text "Table 9: Minimum equivalence thresholds (alpha = 0.05)"
display as text "  (cf. Dette & Schumann, 2024, Table 9)"
display as text "{hline 60}"
display as text %20s "Test" %13s "T=1" %13s "T=2" %13s "T=3"
display as text "{hline 60}"
display as text %20s "delta*_IU"    %13.4f `delta_iu_T1'   %13.4f `delta_iu_T2'   %13.4f `delta_iu_T3'
display as text %20s "delta*_Boot"  %13.4f `delta_boot_T1'  %13.4f `delta_boot_T2'  %13.4f `delta_boot_T3'
display as text %20s "delta*_c.Boot" %13.4f `delta_wild_T1' %13.4f `delta_wild_T2' %13.4f `delta_wild_T3'
display as text %20s "tau*"         %13.4f `tau_T1'         %13.4f `tau_T2'         %13.4f `tau_T3'
display as text %20s "zeta*"        %13.4f `zeta_T1'        %13.4f `zeta_T2'        %13.4f `zeta_T3'
display as text "{hline 60}"

* =============================================================================
* 3. TESTING WITH A SPECIFIED THRESHOLD
* =============================================================================
*
* Threshold must exceed the minimum threshold for rejection.
* From Part 2: delta*=0.147, tau*=0.093, zeta*=0.191

* Max IU: delta=0.15 > delta*=0.147 → should reject
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.15) ///
    vce(cluster) cluster(ID)

* Mean: tau=0.10 > tau*=0.093 → should reject
equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.10) ///
    vce(cluster) cluster(ID)

* RMS: zeta=0.20 > zeta*=0.191 → should reject
equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.20) seed(2024)


* =============================================================================
* 4. VISUALIZATION
* =============================================================================
*
* equivtest_plot reads e() results and auto-generates subtitle + note.
* ±threshold lines are EXACT for max test, CONSERVATIVE for mean/rms.

* Fig 1: Max IU + threshold + CI
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.15) ///
    vce(cluster) cluster(ID)
equivtest_plot, ci ///
    title("Di Tella & Schargrodsky (2004): Max IU Test") ///
    xtitle("Month (relative to July)") ///
    name(fig1_max_iu, replace)

* Fig 2: Mean test + connected points
equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.10) ///
    vce(cluster) cluster(ID)
equivtest_plot, connect ///
    title("Di Tella & Schargrodsky (2004): Mean Test") ///
    xtitle("Month (relative to July)") ///
    name(fig2_mean, replace)

* Fig 3: RMS with minimum threshold (auto from e())
equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) seed(2024)
equivtest_plot, connect ///
    title("Di Tella & Schargrodsky (2004): RMS Test") ///
    xtitle("Month (relative to July)") ///
    name(fig3_rms, replace)

* Fig 4: Max IU minimum threshold + CI
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)
equivtest_plot, ci ///
    title("Di Tella & Schargrodsky (2004): Min. Threshold") ///
    xtitle("Month (relative to July)") ///
    name(fig4_max_min, replace)

* Fig 5: Custom styling
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) threshold(0.15) ///
    vce(cluster) cluster(ID)
equivtest_plot, ci ///
    msymbol(D) msize(medlarge) mcolor(dkgreen) ///
    cilcolor(dkgreen) ///
    threshlcolor(cranberry) threshlwidth(medthick) threshlpattern(shortdash) ///
    title("Custom Styled Plot") ///
    xtitle("Month (relative to July)") ///
    ytitle("Placebo coefficient estimate") ///
    name(fig5_custom, replace)

* Fig 6: Override auto-generated subtitle and note
equivtest_plot, ci ///
    title("Di Tella & Schargrodsky (2004)") ///
    subtitle("IU test, cluster-robust SE, delta = 0.15") ///
    note("Source: Di Tella & Schargrodsky (2004). Buenos Aires crime panel.") ///
    xtitle("Month (relative to July)") ///
    name(fig6_override, replace)


* =============================================================================
* 5. VARIANCE ESTIMATION OPTIONS
* =============================================================================
*
* Max (IU) and Mean tests support multiple VCE types.
* RMS uses self-normalized subsampling (no VCE needed).

foreach vce in ols robust hc2 hc3 hac {
    equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
        pretreatment(4 5 6 7) baseperiod(7) vce(`vce')
    local delta_`vce' = e(min_threshold)
}

* Cluster-robust variants require cluster()
foreach vce in cluster cr1 hc1_cluster {
    equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
        pretreatment(4 5 6 7) baseperiod(7) vce(`vce') cluster(ID)
    local delta_`vce' = e(min_threshold)
}


* =============================================================================
* 5b. CONTROL VARIABLES (Conditional Parallel Trends, Section 5)
* =============================================================================
*
* The x() option includes additional regressors in the TWFE placebo regression.
* Time-invariant covariates are absorbed by individual fixed effects after
* double demeaning. For the full conditional PTA (Eq. 5.5), construct
* covariate-by-time interactions in your data and pass them via x().

equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) x(edpub estserv banco) ///
    vce(cluster) cluster(ID)

equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) x(edpub estserv banco) ///
    vce(cluster) cluster(ID)


* =============================================================================
* 6. MONTE CARLO SIMULATION
* =============================================================================
*
* equivsim generates panel data following the DGP in Dette & Schumann (2024).

* Under H0: non-zero placebo coefficients
equivsim, n(200) preperiods(3) beta(0.1 0.1 0.1) dgp(ar3) het(1) seed(42) clear
xtset id period
equivtest Y, type(max) method(iu) id(id) group(G) time(period) threshold(0.3)
equivtest_plot, ci title("Simulated: beta = (0.1, 0.1, 0.1)") name(fig7_sim_h0, replace)

* Under H1: perfect parallel trends
equivsim, n(200) preperiods(3) beta(0 0 0) dgp(ar3) het(1) seed(42) clear
xtset id period
equivtest Y, type(max) method(iu) id(id) group(G) time(period) threshold(0.3)
equivtest_plot, ci title("Simulated: beta = (0, 0, 0)") name(fig8_sim_h1, replace)

* Other DGP options
equivsim, n(100) preperiods(4) dgp(spherical) seed(123) clear
equivsim, n(100) preperiods(4) dgp(ar1) het(1) seed(123) clear
equivsim, n(100) preperiods(4) dgp(ar3) phi(0.4 0.2 0.05) het(0.5) dip(-0.3) seed(123) clear
equivsim, n(100) preperiods(4) piatt(0.5) trend(0.02) seed(123) clear


* =============================================================================
graph close _all
