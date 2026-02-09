{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{vieweralsosee "equivsim" "help equivsim"}{...}
{vieweralsosee "equivtest_plot" "help equivtest_plot"}{...}
{viewerjumpto "Syntax" "equivtest##syntax"}{...}
{viewerjumpto "Description" "equivtest##description"}{...}
{viewerjumpto "Options" "equivtest##options"}{...}
{viewerjumpto "Examples" "equivtest##examples"}{...}
{viewerjumpto "Stored results" "equivtest##results"}{...}
{viewerjumpto "References" "equivtest##references"}{...}
{viewerjumpto "Also see" "equivtest##alsosee"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:equivtest} {hline 2}}Equivalence testing for pre-trends in difference-in-differences estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:equivtest}
{depvar}
{cmd:,}
{opt type(string)}
{opt id(varname)}
{opt g:roup(varname)}
{opt t:ime(varname)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt type(string)}}test type: {cmd:max}, {cmd:mean}, or {cmd:rms}{p_end}
{synopt:{opt id(varname)}}panel identifier variable{p_end}
{synopt:{opt g:roup(varname)}}treatment group indicator (0/1); alias: {cmd:g()}{p_end}
{synopt:{opt t:ime(varname)}}time period variable; alias: {cmd:period()}{p_end}

{syntab:Model}
{synopt:{opt x(varlist)}}control variables{p_end}
{synopt:{opt pre:treatment(numlist)}}pre-treatment periods to include{p_end}
{synopt:{opt base:period(#)}}base period for placebo construction{p_end}

{syntab:Test}
{synopt:{opt thresh:old(#)}}equivalence threshold; if omitted, computes minimum threshold{p_end}
{synopt:{opt alpha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt method(string)}}method for {cmd:type(max)}: {cmd:iu}, {cmd:boot}, or {cmd:wild}; default is {cmd:iu}{p_end}
{synopt:{opt nboot(#)}}bootstrap replications; default is {cmd:nboot(1000)}{p_end}
{synopt:{opt nolambda(#)}}number of subsamples for RMS test; default is {cmd:nolambda(5)}{p_end}
{synopt:{opt seed(#)}}random seed for reproducibility{p_end}
{synopt:{opt nodots}}suppress bootstrap progress display{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator: {cmd:ols}, {cmd:robust}, {cmd:hc2}, {cmd:hc3}, {cmd:hac}, {cmd:cluster}, {cmd:cr0}, {cmd:cr1}, or {cmd:hc1_cluster}{p_end}
{synopt:{opt cluster(varname)}}cluster variable; requires {cmd:vce(cluster)}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:equivtest} performs equivalence testing for pre-trends in difference-in-differences
(DiD) estimation. It provides a unified interface to three test types:

{phang2}
{cmd:type(max)} - Maximum absolute placebo coefficient test using the Intersection-Union
approach or bootstrap methods.

{phang2}
{cmd:type(mean)} - Mean placebo coefficient test based on the average of placebo effects.

{phang2}
{cmd:type(rms)} - Root mean squared placebo coefficient test.

{pstd}
The tests provide a framework for equivalence testing
as an alternative to traditional pre-trend tests in DiD designs.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt type(string)} specifies the test type. {cmd:max} performs the maximum absolute
coefficient test, {cmd:mean} performs the mean coefficient test, and {cmd:rms}
performs the root mean squared coefficient test.

{phang}
{opt id(varname)} specifies the panel identifier variable.

{phang}
{opt group(varname)} specifies the treatment group indicator, which must be 0 for
control units and 1 for treated units. Abbreviation {cmd:g()} is accepted.

{phang}
{opt time(varname)} specifies the time period variable. Alias {cmd:period()} is 
also accepted for backward compatibility.

{dlgtab:Model}

{phang}
{opt x(varlist)} specifies control variables to include in the regression.

{phang}
{opt pretreatment(numlist)} specifies which pre-treatment periods to include.
If omitted, all pre-treatment periods are used.

{phang}
{opt baseperiod(#)} specifies the base period for placebo coefficient construction.
If omitted, the last pre-treatment period is used.

{dlgtab:Test}

{phang}
{opt threshold(#)} specifies the equivalence threshold. If omitted, the command
computes the minimum threshold at which equivalence can be established.

{phang}
{opt alpha(#)} specifies the significance level. The default is {cmd:alpha(0.05)}.
For {cmd:type(rms)}, alpha must be one of 0.01, 0.025, 0.05, 0.1, or 0.2.

{phang}
{opt method(string)} specifies the method for {cmd:type(max)}. {cmd:iu} uses the
Intersection-Union approach (default), {cmd:boot} uses parametric bootstrap 
valid under spherical (homoskedastic) errors, and {cmd:wild} uses wild cluster 
bootstrap robust to heteroskedasticity and serial correlation.

{phang}
{opt nboot(#)} specifies the number of bootstrap replications for {cmd:method(boot)}
or {cmd:method(wild)}. The default is {cmd:nboot(1000)}.

{phang}
{opt nolambda(#)} specifies the number of subsamples K for {cmd:type(rms)}.
The default is {cmd:nolambda(5)}.

{phang}
{opt seed(#)} specifies the random seed for reproducibility in bootstrap and
RMS tests.

{phang}
{opt nodots} suppresses the progress display during bootstrap replications.
By default, the {cmd:boot} and {cmd:wild} methods show a dot-based progress 
indicator as bootstrap iterations proceed. This option is useful for batch 
processing or when embedding the command in loops.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies the variance estimator. Not available with 
{cmd:type(rms)} or bootstrap methods. Options are:

{phang2}
{cmd:ols} - Standard OLS variance (homoskedastic). Default.

{phang2}
{cmd:robust} or {cmd:hc1} - HC1 heteroskedasticity-robust (White 1980).

{phang2}
{cmd:hc2} - HC2 leverage-adjusted heteroskedasticity-robust (MacKinnon & White 1985).
Better finite-sample properties than HC1.

{phang2}
{cmd:hc3} - HC3 more conservative leverage adjustment (Davidson & MacKinnon 1993).
Recommended for small samples.

{phang2}
{cmd:hac} - Arellano HAC estimator for panel data.

{phang2}
{cmd:cluster} or {cmd:cr0} - CR0 cluster-robust without small-sample adjustment.
Requires {opt cluster()} option.

{phang2}
{cmd:cr1} - CR1 cluster-robust with G/(G-1) small-sample adjustment (Stata default).
More conservative than CR0. Requires {opt cluster()} option.

{phang2}
{cmd:hc1_cluster} - HC1 cluster-robust with small-sample adjustment.
Requires {opt cluster()} option.

{phang}
{opt cluster(varname)} specifies the cluster variable for cluster-robust
standard errors. Required for {cmd:vce(cluster)}, {cmd:vce(cr0)}, {cmd:vce(cr1)}, and {cmd:vce(hc1_cluster)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. gen treat = (idcode > 2000)}{p_end}

{pstd}Maximum test with IU method{p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) id(idcode) group(treat) time(year)}{p_end}

{pstd}Maximum test with specified threshold{p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) id(idcode) group(treat) time(year) threshold(0.1)}{p_end}

{pstd}Maximum test with bootstrap{p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) id(idcode) group(treat) time(year) method(boot) nboot(500)}{p_end}

{pstd}Mean test{p_end}
{phang2}{cmd:. equivtest ln_wage, type(mean) id(idcode) group(treat) time(year)}{p_end}

{pstd}RMS test{p_end}
{phang2}{cmd:. equivtest ln_wage, type(rms) id(idcode) group(treat) time(year) alpha(0.05)}{p_end}

{pstd}With robust standard errors{p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) id(idcode) group(treat) time(year) vce(robust)}{p_end}

{pstd}Using short option names (g and period are also accepted){p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) id(idcode) g(treat) period(year)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:equivtest} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars (common)}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of individuals{p_end}
{synopt:{cmd:e(no_placebos)}}number of placebo coefficients{p_end}
{synopt:{cmd:e(alpha)}}significance level{p_end}
{synopt:{cmd:e(base_period)}}base period{p_end}
{synopt:{cmd:e(is_balanced)}}1 if balanced panel, 0 otherwise{p_end}
{synopt:{cmd:e(T_pre)}}number of pre-treatment periods (balanced panels only){p_end}
{synopt:{cmd:e(T_min)}}minimum periods per individual (unbalanced panels only){p_end}
{synopt:{cmd:e(T_max)}}maximum periods per individual (unbalanced panels only){p_end}
{synopt:{cmd:e(threshold_specified)}}1 if threshold specified, 0 otherwise{p_end}
{synopt:{cmd:e(threshold)}}equivalence threshold (if specified){p_end}
{synopt:{cmd:e(reject)}}1 if null rejected, 0 otherwise (if threshold specified){p_end}
{synopt:{cmd:e(min_threshold)}}minimum threshold (if threshold not specified){p_end}

{p2col 5 28 32 2: Scalars (type=max)}{p_end}
{synopt:{cmd:e(max_abs_coef)}}maximum absolute placebo coefficient{p_end}
{synopt:{cmd:e(nboot)}}number of bootstrap replications (method=boot or wild){p_end}
{synopt:{cmd:e(wild)}}1 if wild bootstrap, 0 otherwise (method=boot or wild){p_end}
{synopt:{cmd:e(boot_critical)}}bootstrap critical value (method=boot or wild){p_end}

{p2col 5 28 32 2: Scalars (type=mean)}{p_end}
{synopt:{cmd:e(abs_mean_placebo)}}absolute value of mean placebo coefficient{p_end}
{synopt:{cmd:e(var_mean_placebo)}}variance of mean placebo coefficient{p_end}
{synopt:{cmd:e(se_mean_placebo)}}standard error of mean placebo coefficient{p_end}
{synopt:{cmd:e(p_value)}}p-value (if threshold specified){p_end}
{synopt:{cmd:e(mean_critical_value)}}critical value for mean test (if threshold specified){p_end}

{p2col 5 28 32 2: Scalars (type=rms)}{p_end}
{synopt:{cmd:e(rms_placebo_coefs)}}root mean square of placebo coefficients{p_end}
{synopt:{cmd:e(nolambda)}}number of lambda subsamples{p_end}
{synopt:{cmd:e(rms_critical_value)}}RMS critical value (if threshold specified){p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:equivtest}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(type)}}test type{p_end}
{synopt:{cmd:e(method)}}method (for type=max){p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(idvar)}}panel identifier{p_end}
{synopt:{cmd:e(groupvar)}}treatment group variable{p_end}
{synopt:{cmd:e(timevar)}}time variable{p_end}
{synopt:{cmd:e(vce)}}variance estimator{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if cluster-robust VCE){p_end}
{synopt:{cmd:e(placebo_coef_names)}}names of placebo coefficients{p_end}
{synopt:{cmd:e(preperiods)}}specified pre-treatment periods (if pretreatment() used){p_end}
{synopt:{cmd:e(properties)}}{cmd:b}{p_end}

{p2col 5 28 32 2: Matrices}{p_end}
{synopt:{cmd:e(b_placebo)}}placebo coefficient vector (not stored for method=boot or method=wild){p_end}
{synopt:{cmd:e(V_placebo)}}variance-covariance matrix (type=max with method=iu, or type=mean){p_end}
{synopt:{cmd:e(se_placebo)}}standard errors (type=max with method=iu){p_end}
{synopt:{cmd:e(IU_critical_values)}}critical values per coefficient (type=max, method=iu, threshold specified){p_end}
{synopt:{cmd:e(min_equiv_thresholds)}}minimum thresholds per coefficient (type=max, method=iu, threshold not specified){p_end}


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
Cai Xuanyu ({browse "mailto:xuanyuCAI@outlook.com":xuanyuCAI@outlook.com}), 
Xu Wenli ({browse "mailto:wlxu@cityu.edu.mo":wlxu@cityu.edu.mo})


{marker alsosee}{...}
{title:Also see}

{psee}
Online: {helpb maxequivtest}, {helpb meanequivtest}, {helpb rmsequivtest}, 
{helpb equivsim}, {helpb equivtest_plot}
{p_end}
