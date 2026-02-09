{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{viewerjumpto "Syntax" "maxequivtest##syntax"}{...}
{viewerjumpto "Description" "maxequivtest##description"}{...}
{viewerjumpto "Options" "maxequivtest##options"}{...}
{viewerjumpto "Methods" "maxequivtest##methods"}{...}
{viewerjumpto "Examples" "maxequivtest##examples"}{...}
{viewerjumpto "Stored results" "maxequivtest##results"}{...}
{viewerjumpto "References" "maxequivtest##references"}{...}
{viewerjumpto "Also see" "maxequivtest##alsosee"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{cmd:maxequivtest} {hline 2}}Maximum equivalence test for pre-trends in difference-in-differences estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:maxequivtest}
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
{synopt:{opt th:reshold(#)}}equivalence threshold δ{p_end}
{synopt:{opt pre:treatment(numlist)}}pre-treatment periods{p_end}
{synopt:{opt base:period(#)}}base period for placebo construction{p_end}
{synopt:{opt method(string)}}test method: {cmd:iu}, {cmd:boot}, or {cmd:wild}; default is {cmd:iu}{p_end}
{synopt:{opt alpha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt reps(#)}}bootstrap replications; default is {cmd:reps(1000)}{p_end}
{synopt:{opt seed(#)}}random seed for reproducibility{p_end}
{synopt:{opt vce(vcetype)}}variance-covariance estimator (IU method only){p_end}
{synopt:{opt cluster(varname)}}cluster variable; requires {cmd:vce(cluster)}{p_end}
{synopt:{opt nodots}}suppress bootstrap progress display{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:maxequivtest} implements the Maximum Equivalence Test for pre-trends in 
Difference-in-Differences (DiD) estimation.

{pstd}
The test examines whether the maximum absolute value of placebo coefficients 
is within an equivalence bound, providing evidence for the parallel trends 
assumption.

{pstd}
{bf:Hypotheses:}

{pmore}
H0: ||β||∞ ≥ δ  (pre-trends are NOT equivalent to zero)

{pmore}
H1: ||β||∞ < δ  (pre-trends ARE equivalent to zero)

{pstd}
where ||β||∞ = max|β_l| is the maximum absolute value of placebo coefficients.

{pstd}
If {opt threshold()} is specified, the command performs a hypothesis test.
If {opt threshold()} is omitted, the command computes the minimum equivalence 
threshold δ* at which H0 can be rejected.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt id(varname)} specifies the variable that identifies panel individuals.

{phang}
{opt group(varname)} specifies the treatment group indicator variable. This must 
be a binary variable (0/1) where 1 indicates the treatment group. 
Abbreviation {cmd:g()} is accepted.

{phang}
{opt time(varname)} specifies the time period variable. 
Alias {cmd:period()} is also accepted for backward compatibility.

{dlgtab:Optional}

{phang}
{opt x(varlist)} specifies control variables to include in the TWFE regression.

{phang}
{opt threshold(#)} specifies the equivalence threshold δ. Must be strictly positive (> 0).
If omitted, the command computes the minimum equivalence threshold.

{phang}
{opt pretreatment(numlist)} specifies which periods are pre-treatment periods.
If omitted, all periods before the maximum period are used.

{phang}
{opt baseperiod(#)} specifies the base period for placebo construction. 
If omitted, defaults to the last pre-treatment period.

{phang}
{opt method(string)} specifies the test method:

{phang2}
{cmd:iu} - Intersection-Union test (default). Analytical test based on the 
folded normal distribution. Fast and exact.

{phang2}
{cmd:boot} - Parametric bootstrap test. Valid under spherical errors 
(homoskedastic and serially uncorrelated). Uses constrained estimation 
under the null hypothesis.

{phang2}
{cmd:wild} - Wild cluster bootstrap test using Rademacher weights. Robust to 
heteroskedasticity and serial correlation within clusters. Recommended for 
empirical applications.

{phang}
{opt alpha(#)} specifies the significance level. Default is 0.05.

{phang}
{opt reps(#)} specifies the number of bootstrap replications for {cmd:boot} 
and {cmd:wild} methods. Default is 1000. Ignored for {cmd:iu} method.

{phang}
{opt seed(#)} specifies the random seed for reproducibility.

{phang}
{opt vce(vcetype)} specifies the variance-covariance estimator for the 
{cmd:iu} method. Options are:

{phang2}
{cmd:ols} - Standard OLS variance (homoskedastic). Default.

{phang2}
{cmd:robust} or {cmd:hc1} - HC1 heteroskedasticity-robust (White 1980). 
Corresponds to R package {cmd:vcov="HC"}.

{phang2}
{cmd:hc2} - HC2 leverage-adjusted heteroskedasticity-robust (MacKinnon & White 1985).
Better finite-sample properties than HC1. Uses Stata native {cmd:vce(hc2)}.

{phang2}
{cmd:hc3} - HC3 more conservative leverage adjustment (Davidson & MacKinnon 1993).
Recommended for small samples. Uses Stata native {cmd:vce(hc3)}.

{phang2}
{cmd:hac} - Arellano HAC estimator. Corresponds to R package {cmd:vcov="HAC"}.

{phang2}
{cmd:cluster} - CR0 cluster-robust without small-sample adjustment.
Corresponds to R package {cmd:vcov="CL"}. Requires {opt cluster()} option.

{phang2}
{cmd:cr0} - CR0 cluster-robust using Stata native {cmd:vce(cr0)}.
Requires {opt cluster()} option.

{phang2}
{cmd:cr1} - CR1 cluster-robust with G/(G-1) small-sample adjustment (Stata default).
More conservative than CR0. Uses Stata native {cmd:vce(cluster)}. 
Requires {opt cluster()} option.

{phang2}
{cmd:hc1_cluster} - HC1 cluster-robust with small-sample adjustment.
Requires {opt cluster()} option.

{phang}
{opt cluster(varname)} specifies the cluster variable for cluster-robust 
standard errors. Required for {cmd:vce(cluster)}, {cmd:vce(cr0)}, {cmd:vce(cr1)}, and {cmd:vce(hc1_cluster)}.

{phang}
{opt nodots} suppresses the progress display during bootstrap replications.
By default, the {cmd:boot} and {cmd:wild} methods show a dot-based progress 
indicator as bootstrap iterations proceed. This option is useful for batch 
processing or when embedding the command in loops.


{marker methods}{...}
{title:Methods}

{pstd}
{bf:IU (Intersection-Union) Method}

{pstd}
The IU method applies the intersection-union principle (Berger and Hsu, 1996)
using the folded normal distribution to construct critical values.
For each placebo coefficient β̂_l with estimated variance Σ̂_ll/n, the test 
rejects H0 if:

{pmore}
|β̂_l| < Q_{N_F(δ, Σ̂_ll/n)}(α)  for all l = 1,...,T

{pstd}
where Q_{N_F(μ, σ²)}(α) denotes the α-quantile of the folded normal 
distribution with mean μ and variance σ². This test is computationally 
attractive but tends to be conservative, especially when T is large.

{pstd}
{bf:Bootstrap Methods}

{pstd}
The bootstrap methods construct the null distribution of the test statistic 
through resampling. Both methods first compute a constrained estimator β̂̂_c 
that satisfies the null hypothesis ||β||∞ = δ, then generate bootstrap 
samples under this constrained model.

{phang2}
{cmd:boot} - Parametric bootstrap valid under spherical (homoskedastic and 
serially uncorrelated) errors. Generates bootstrap errors from 
N(0, σ̂̂_c), where σ̂̂_c is the constrained variance estimate.

{phang2}
{cmd:wild} - Wild cluster bootstrap using Rademacher weights (±1 with equal 
probability). Robust to heteroskedasticity and serial correlation within 
clusters. Recommended for empirical applications.


{marker examples}{...}
{title:Examples}

{pstd}Setup: Generate simulation data{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(4) beta(0 0 0 0) clear seed(12345)}{p_end}

{pstd}Basic IU test with specified threshold{p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period) threshold(0.5)}{p_end}

{pstd}Compute minimum equivalence threshold{p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period)}{p_end}
{phang2}{cmd:. display "Minimum threshold: " e(min_threshold)}{p_end}

{pstd}Bootstrap test{p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period) threshold(0.5) method(boot) reps(1000) seed(12345)}{p_end}

{pstd}Wild bootstrap test{p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period) threshold(0.5) method(wild) reps(1000) seed(12345)}{p_end}

{pstd}With robust standard errors (IU method){p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period) threshold(0.5) vce(robust)}{p_end}

{pstd}With cluster-robust standard errors{p_end}
{phang2}{cmd:. maxequivtest Y, id(id) group(G) time(period) threshold(0.5) vce(cluster) cluster(id)}{p_end}

{pstd}Using short option names (g and period are also accepted){p_end}
{phang2}{cmd:. maxequivtest Y, id(id) g(G) period(period) threshold(0.5)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:maxequivtest} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}total number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of individuals{p_end}
{synopt:{cmd:e(T)}}number of time periods{p_end}
{synopt:{cmd:e(no_placebos)}}number of placebo coefficients{p_end}
{synopt:{cmd:e(max_abs_coef)}}maximum absolute placebo coefficient{p_end}
{synopt:{cmd:e(alpha)}}significance level{p_end}
{synopt:{cmd:e(base_period)}}base period{p_end}
{synopt:{cmd:e(is_balanced)}}1 if balanced panel, 0 otherwise{p_end}
{synopt:{cmd:e(threshold_specified)}}1 if threshold specified, 0 otherwise{p_end}

{p2col 5 25 29 2: For unbalanced panels (when e(is_balanced)==0):}{p_end}
{synopt:{cmd:e(T_min)}}minimum number of time periods per individual{p_end}
{synopt:{cmd:e(T_max)}}maximum number of time periods per individual{p_end}

{p2col 5 25 29 2: When threshold specified:}{p_end}
{synopt:{cmd:e(threshold)}}equivalence threshold{p_end}
{synopt:{cmd:e(critical_value)}}critical value (IU method){p_end}
{synopt:{cmd:e(reject)}}1 if H0 rejected, 0 otherwise{p_end}

{p2col 5 25 29 2: When threshold not specified:}{p_end}
{synopt:{cmd:e(min_threshold)}}minimum equivalence threshold δ*{p_end}

{p2col 5 25 29 2: For bootstrap methods:}{p_end}
{synopt:{cmd:e(B)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(seed)}}random seed (. if not specified){p_end}
{synopt:{cmd:e(sigma_hathat_c)}}constrained sigma estimate{p_end}

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Matrices}{p_end}
{synopt:{cmd:e(b_placebo)}}placebo coefficient vector (no_placebos × 1){p_end}
{synopt:{cmd:e(V_placebo)}}variance-covariance matrix (IU method only){p_end}
{synopt:{cmd:e(se_placebo)}}standard errors (IU method only){p_end}
{synopt:{cmd:e(crit_values)}}critical values for each placebo (IU method with threshold){p_end}
{synopt:{cmd:e(min_thresholds)}}minimum thresholds for each placebo (IU method without threshold){p_end}

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"maxequivtest"{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(idvar)}}individual variable name{p_end}
{synopt:{cmd:e(groupvar)}}group variable name{p_end}
{synopt:{cmd:e(timevar)}}time variable name{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable name (if specified){p_end}
{synopt:{cmd:e(placebo_names)}}names of placebo variables{p_end}
{synopt:{cmd:e(method)}}test method used{p_end}
{synopt:{cmd:e(vce)}}VCE type{p_end}


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
Online: {helpb equivtest}, {helpb meanequivtest}, {helpb rmsequivtest}, 
{helpb equivsim}, {helpb equivtest_plot}
{p_end}
