{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{vieweralsosee "equivsim" "help equivsim"}{...}
{vieweralsosee "equivtest_plot" "help equivtest_plot"}{...}
{viewerjumpto "Syntax" "meanequivtest##syntax"}{...}
{viewerjumpto "Description" "meanequivtest##description"}{...}
{viewerjumpto "Options" "meanequivtest##options"}{...}
{viewerjumpto "Stored results" "meanequivtest##results"}{...}
{viewerjumpto "Examples" "meanequivtest##examples"}{...}
{viewerjumpto "References" "meanequivtest##references"}{...}
{viewerjumpto "Authors" "meanequivtest##authors"}{...}
{viewerjumpto "Also see" "meanequivtest##alsosee"}{...}
{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{cmd:meanequivtest} {hline 2}}Mean equivalence test for pre-trends in difference-in-differences estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:meanequivtest}
{depvar}
{cmd:,}
{opt id(varname)}
{opt g:roup(varname)}
{opt t:ime(varname)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt id(varname)}}panel individual identifier{p_end}
{synopt:{opt g:roup(varname)}}treatment group indicator (0/1); alias: {cmd:g()}{p_end}
{synopt:{opt t:ime(varname)}}time period variable; alias: {cmd:period()}{p_end}

{syntab:Optional}
{synopt:{opt x(varlist)}}control variables{p_end}
{synopt:{opt pre:treatment(numlist)}}pre-treatment periods to include{p_end}
{synopt:{opt base:period(#)}}base period; default is max of pretreatment{p_end}
{synopt:{opt thresh:old(#)}}equivalence threshold τ{p_end}
{synopt:{opt alpha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt vce(vcetype)}}variance-covariance estimator type{p_end}
{synopt:{opt cluster(varname)}}cluster variable for clustered standard errors{p_end}
{synopt:{opt nodis:play}}suppress output display{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:vcetype} may be {opt ols}, {opt hc1} (or {opt robust}), {opt hc2}, {opt hc3}, {opt hac}, 
{opt cluster}, {opt cr0}, {opt cr1}, or {opt hc1_cluster}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:meanequivtest} implements the Mean Equivalence Test for pre-trends in 
Difference-in-Differences (DiD) estimation.

{pstd}
The test examines whether the mean of placebo coefficients is within an 
equivalence bound, testing:

{pmore}
H0: |β̄| >= τ  vs  H1: |β̄| < τ

{pstd}
where β̄ = (1/T) Σ β_t is the mean of placebo coefficients from a two-way 
fixed effects regression.

{pstd}
If {opt threshold()} is specified, the command performs a hypothesis test 
and reports the critical value, p-value, and rejection decision.

{pstd}
If {opt threshold()} is not specified, the command computes the minimum 
equivalence threshold τ* such that the null hypothesis would be rejected 
at the specified significance level.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt id(varname)} specifies the variable that identifies the panel units 
(individuals, firms, etc.).

{phang}
{opt group(varname)} specifies the treatment group indicator variable. This 
should be a binary variable equal to 1 for treated units and 0 for control 
units. Abbreviation {cmd:g()} is accepted.

{phang}
{opt time(varname)} specifies the time period variable. 
Alias {cmd:period()} is also accepted for backward compatibility.

{dlgtab:Optional}

{phang}
{opt x(varlist)} specifies control variables to include in the regression.

{phang}
{opt pretreatment(numlist)} specifies which pre-treatment periods to include
in the analysis. If omitted, all periods before the maximum period are used.
The base period will be excluded from the placebo coefficients.

{phang}
{opt baseperiod(#)} specifies the base (reference) period for placebo construction.
The default is the maximum value in {opt pretreatment()} or the last pre-treatment
period if pretreatment is not specified.

{phang}
{opt threshold(#)} specifies the equivalence threshold τ. Must be 
strictly positive. If not specified, the minimum threshold is computed.

{phang}
{opt alpha(#)} specifies the significance level. The default is 0.05.

{phang}
{opt vce(vcetype)} specifies the type of variance-covariance matrix 
estimator. The following types are available:

{phang2}
{opt ols} uses the standard OLS variance estimator (homoskedastic). This 
corresponds to the default in R's {cmd:plm} package.

{phang2}
{opt hc1} or {opt robust} uses the HC1 heteroskedasticity-robust estimator 
(White 1980). This corresponds to {cmd:vcovHC(type="HC1", method="white1")} 
in R's {cmd:plm}.

{phang2}
{opt hc2} uses the HC2 leverage-adjusted heteroskedasticity-robust estimator 
(MacKinnon & White 1985). Better finite-sample properties than HC1. Uses 
Stata's native {cmd:vce(hc2)}.

{phang2}
{opt hc3} uses the HC3 more conservative leverage-adjusted estimator 
(Davidson & MacKinnon 1993). Recommended for small samples. Uses Stata's 
native {cmd:vce(hc3)}.

{phang2}
{opt hac} uses the Arellano HAC estimator. This corresponds to 
{cmd:vcovHC(type="HC3", method="arellano")} in R's {cmd:plm}.

{phang2}
{opt cluster} uses the CR0 cluster-robust estimator without small-sample 
adjustment (custom Mata implementation). This corresponds to 
{cmd:vcovCR(type="CR0")} in R's {cmd:clubSandwich} package. 
Requires {opt cluster()} to be specified.

{phang2}
{opt cr0} uses the CR0 cluster-robust estimator without small-sample 
adjustment (native Stata implementation). Requires {opt cluster()} to be 
specified.

{phang2}
{opt cr1} uses the CR1 cluster-robust estimator with G/(G-1) small-sample 
adjustment (Stata default). More conservative than CR0. Uses Stata's native 
{cmd:vce(cluster)}. Requires {opt cluster()} to be specified.

{phang2}
{opt hc1_cluster} uses the HC1 cluster-robust estimator with small-sample 
adjustment. Requires {opt cluster()} to be specified.

{phang}
{opt cluster(varname)} specifies the variable for cluster-robust standard 
errors. Required when {cmd:vce(cluster)}, {cmd:vce(cr0)}, {cmd:vce(cr1)}, 
or {cmd:vce(hc1_cluster)} is specified.

{phang}
{opt nodisplay} suppresses the output display. Estimation results are still 
stored in {cmd:e()}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:meanequivtest} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of panel units{p_end}
{synopt:{cmd:e(T)}}number of time periods{p_end}
{synopt:{cmd:e(no_placebos)}}number of placebo coefficients{p_end}
{synopt:{cmd:e(abs_mean_placebo)}}absolute mean of placebo coefficients |β̄|{p_end}
{synopt:{cmd:e(var_mean_placebo)}}variance of mean placebo σ̂²{p_end}
{synopt:{cmd:e(se_mean_placebo)}}standard error of mean placebo σ̂{p_end}
{synopt:{cmd:e(alpha)}}significance level{p_end}
{synopt:{cmd:e(base_period)}}base period{p_end}
{synopt:{cmd:e(is_balanced)}}1 if balanced panel, 0 otherwise{p_end}
{synopt:{cmd:e(T_min)}}minimum periods per individual (unbalanced panels only){p_end}
{synopt:{cmd:e(T_max)}}maximum periods per individual (unbalanced panels only){p_end}
{synopt:{cmd:e(threshold)}}equivalence threshold (if specified){p_end}
{synopt:{cmd:e(critical_value)}}critical value (if threshold specified){p_end}
{synopt:{cmd:e(p_value)}}p-value (if threshold specified){p_end}
{synopt:{cmd:e(reject)}}1 if H0 rejected, 0 otherwise (if threshold specified){p_end}
{synopt:{cmd:e(min_threshold)}}minimum threshold (if threshold not specified){p_end}
{synopt:{cmd:e(threshold_specified)}}1 if threshold specified, 0 otherwise{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:meanequivtest}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(idvar)}}name of panel ID variable{p_end}
{synopt:{cmd:e(groupvar)}}name of treatment group variable{p_end}
{synopt:{cmd:e(timevar)}}name of time variable{p_end}
{synopt:{cmd:e(vce)}}VCE type{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable (if used){p_end}
{synopt:{cmd:e(placebo_names)}}names of placebo coefficients{p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt:{cmd:e(b_placebo)}}placebo coefficient vector (no_placebos × 1){p_end}
{synopt:{cmd:e(V_placebo)}}placebo VCE matrix (no_placebos × no_placebos){p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup: Create panel data{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. gen id = ceil(_n/5)}{p_end}
{phang2}{cmd:. bysort id: gen t = _n}{p_end}
{phang2}{cmd:. gen treat = (id > 50)}{p_end}
{phang2}{cmd:. gen y = rnormal()}{p_end}

{pstd}Basic usage without threshold (compute minimum threshold){p_end}
{phang2}{cmd:. meanequivtest y, id(id) group(treat) time(t) pretreatment(1 2 3) baseperiod(3)}{p_end}

{pstd}Test with specified threshold{p_end}
{phang2}{cmd:. meanequivtest y, id(id) group(treat) time(t) pretreatment(1 2 3) baseperiod(3) threshold(0.1)}{p_end}

{pstd}With heteroskedasticity-robust standard errors{p_end}
{phang2}{cmd:. meanequivtest y, id(id) group(treat) time(t) pretreatment(1 2 3) baseperiod(3) vce(hc1)}{p_end}

{pstd}With cluster-robust standard errors{p_end}
{phang2}{cmd:. meanequivtest y, id(id) group(treat) time(t) pretreatment(1 2 3) baseperiod(3) cluster(id) vce(cr0)}{p_end}

{pstd}Access stored results{p_end}
{phang2}{cmd:. display "Mean placebo effect: " e(abs_mean_placebo)}{p_end}
{phang2}{cmd:. display "Standard error: " e(se_mean_placebo)}{p_end}
{phang2}{cmd:. matrix list e(b_placebo)}{p_end}

{pstd}Using short option names (g and period are also accepted){p_end}
{phang2}{cmd:. meanequivtest y, id(id) g(treat) period(t) pretreatment(1 2 3) baseperiod(3)}{p_end}


{marker references}{...}
{title:References}

{phang}
Dette, H. and M. Schumann. 2024.
Testing for equivalence of pre-trends in Difference-in-Differences estimation.
{it:Journal of Business & Economic Statistics} 42(4): 1289-1301.
{browse "https://doi.org/10.1080/07350015.2024.2308121"}
{p_end}


{marker authors}{...}
{title:Authors}

{pstd}
Cai Xuanyu ({browse "mailto:xuanyuCAI@outlook.com":xuanyuCAI@outlook.com})
{p_end}

{pstd}
Xu Wenli ({browse "mailto:wlxu@cityu.edu.mo":wlxu@cityu.edu.mo})
{p_end}


{marker alsosee}{...}
{title:Also see}

{psee}
Online: {helpb equivtest}, {helpb maxequivtest}, {helpb rmsequivtest}, 
{helpb equivsim}, {helpb equivtest_plot}
{p_end}
