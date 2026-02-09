{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "equivsim" "help equivsim"}{...}
{vieweralsosee "equivtest_plot" "help equivtest_plot"}{...}
{viewerjumpto "Syntax" "rmsequivtest##syntax"}{...}
{viewerjumpto "Description" "rmsequivtest##description"}{...}
{viewerjumpto "Options" "rmsequivtest##options"}{...}
{viewerjumpto "Examples" "rmsequivtest##examples"}{...}
{viewerjumpto "Stored results" "rmsequivtest##results"}{...}
{viewerjumpto "References" "rmsequivtest##references"}{...}
{viewerjumpto "Authors" "rmsequivtest##authors"}{...}
{viewerjumpto "Also see" "rmsequivtest##alsosee"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{cmd:rmsequivtest} {hline 2}}RMS equivalence test for pre-trends in difference-in-differences estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:rmsequivtest}
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
{synopt:{opt th:reshold(#)}}equivalence threshold ζ{p_end}
{synopt:{opt pre:treatment(numlist)}}pre-treatment periods{p_end}
{synopt:{opt base:period(#)}}base period for placebo construction{p_end}
{synopt:{opt nol:ambda(#)}}number of subsamples K (default: 5){p_end}
{synopt:{opt alpha(#)}}significance level (default: 0.05){p_end}
{synopt:{opt level(#)}}confidence level for CI (default: 95){p_end}
{synopt:{opt seed(#)}}random seed for reproducibility{p_end}
{synopt:{opt nodis:play}}suppress results display{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rmsequivtest} implements the Root Mean Squared (RMS) Equivalence Test for 
pre-trends in Difference-in-Differences (DiD) estimation.

{pstd}
The test examines whether the root mean square of placebo coefficients is 
within an equivalence bound, using a self-normalized pivot statistic based 
on subsampling.

{pstd}
{bf:Hypotheses:}

{pmore}
H0: β_RMS >= ζ  (pre-trends are NOT equivalent to zero)

{pmore}
H1: β_RMS < ζ   (pre-trends ARE equivalent to zero)

{pstd}
where β_RMS = sqrt((1/T) Σ β_l²) is the root mean square of placebo coefficients.

{pstd}
{bf:Test Statistic:}

{pmore}
The test uses a self-normalized pivot statistic based on subsampling:

{pmore}
M̂_n = (β̂_RMS²(1) - β_RMS²) / V̂_n

{pmore}
where V̂_n is the self-normalized variance computed from K subsamples.

{pstd}
{bf:Rejection Rule:}

{pmore}
Reject H0 if: β̂_RMS² < ζ² + Q_W(α) * V̂_n

{pmore}
where Q_W(α) is the α-quantile of the W distribution.


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
{opt threshold(#)} specifies the equivalence threshold ζ. If omitted, the 
command computes the minimum equivalence threshold that would lead to rejection.

{phang}
{opt pretreatment(numlist)} specifies which periods are pre-treatment periods.
If omitted, all periods are used and the maximum period is treated as the 
base period.

{phang}
{opt baseperiod(#)} specifies the base period for placebo construction. 
If omitted, defaults to the last pre-treatment period.

{phang}
{opt nolambda(#)} specifies the number of subsamples K for computing the 
self-normalized variance. Default is 5. Must be at least 2.

{phang}
{opt alpha(#)} specifies the significance level. Default is 0.05. 
Must be one of: 0.01, 0.025, 0.05, 0.1, or 0.2.

{pstd}
{err:{bf:Important:}} The RMS test only supports the following significance levels:{p_end}

{p 8 8 2}
α ∈ {c -(}0.01, 0.025, 0.05, 0.1, 0.2{c )-}

{pstd}
This restriction exists because W distribution critical values are pre-computed 
lookup tables. Using other alpha values will result in an error.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence 
intervals of the RMS and mean squared placebo coefficients. The default is 
{cmd:level(95)}, meaning 95% confidence intervals. Supported values are 80, 90, 
95, 98, and 99.

{pmore}
The confidence interval is computed using the self-normalized approach described 
in Remark 3(b) of Dette and Schumann (2024):

{pmore}
CI_α(β²_RMS) = [ β̂²_RMS + Q_W(α/2)·V̂_n,  β̂²_RMS + Q_W(1-α/2)·V̂_n ]

{pmore}
where Q_W are quantiles of the limiting W distribution and V̂_n is the 
self-normalized variance estimator. The CI for β_RMS is obtained by taking 
square roots, ensuring the lower bound is non-negative.

{phang}
{opt seed(#)} specifies the random seed for reproducibility. The RMS test 
involves random subsampling, so results may vary between runs without a seed.

{phang}
{opt nodisplay} suppresses the display of results. This option is useful for 
batch processing or when embedding the command in loops where output is not 
needed.


{marker examples}{...}
{title:Examples}

{pstd}Setup: Panel data with treatment and control groups{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. gen treat = (race == 2)}{p_end}

{pstd}Basic RMS equivalence test (compute minimum threshold){p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) seed(12345)}{p_end}

{pstd}Test with specific equivalence threshold{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) threshold(0.5) seed(12345)}{p_end}

{pstd}Test with control variables{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) x(age ttl_exp) seed(12345)}{p_end}

{pstd}Test with different significance level{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) alpha(0.1) seed(12345)}{p_end}

{pstd}Test with more subsamples for better variance estimation{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) nolambda(10) seed(12345)}{p_end}

{pstd}Test with 90% confidence interval{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) level(90) seed(12345)}{p_end}

{pstd}Test with 99% confidence interval{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) level(99) seed(12345)}{p_end}

{pstd}Access confidence interval bounds after estimation{p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) group(treat) time(year) seed(12345)}{p_end}
{phang2}{cmd:. display "RMS 95% CI: [" e(RMS_ci_lower) ", " e(RMS_ci_upper) "]"}{p_end}

{pstd}Using short option names (g and period are also accepted){p_end}
{phang2}{cmd:. rmsequivtest ln_wage, id(idcode) g(treat) period(year) seed(12345)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rmsequivtest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}total number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of individuals{p_end}
{synopt:{cmd:e(T)}}number of time periods{p_end}
{synopt:{cmd:e(no_placebos)}}number of placebo coefficients{p_end}
{synopt:{cmd:e(MS_placebo)}}mean squared placebo coefficient{p_end}
{synopt:{cmd:e(RMS_placebo)}}root mean squared placebo coefficient{p_end}
{synopt:{cmd:e(V_n)}}self-normalized variance{p_end}
{synopt:{cmd:e(Q_W)}}W distribution critical value{p_end}
{synopt:{cmd:e(alpha)}}significance level{p_end}
{synopt:{cmd:e(nolambda)}}number of subsamples{p_end}
{synopt:{cmd:e(base_period)}}base period{p_end}
{synopt:{cmd:e(is_balanced)}}1 if balanced panel, 0 otherwise{p_end}
{synopt:{cmd:e(seed)}}random seed (. if not specified){p_end}
{synopt:{cmd:e(threshold_specified)}}1 if threshold specified, 0 otherwise{p_end}

{p2col 5 20 24 2: Confidence Interval Results:}{p_end}
{synopt:{cmd:e(level)}}confidence level for CI{p_end}
{synopt:{cmd:e(RMS_ci_lower)}}lower bound of RMS confidence interval{p_end}
{synopt:{cmd:e(RMS_ci_upper)}}upper bound of RMS confidence interval{p_end}
{synopt:{cmd:e(MS_ci_lower)}}lower bound of mean squared CI{p_end}
{synopt:{cmd:e(MS_ci_upper)}}upper bound of mean squared CI{p_end}

{p2col 5 20 24 2: For unbalanced panels (when e(is_balanced)==0):}{p_end}
{synopt:{cmd:e(T_min)}}minimum number of time periods per individual{p_end}
{synopt:{cmd:e(T_max)}}maximum number of time periods per individual{p_end}

{p2col 5 20 24 2: When threshold specified:}{p_end}
{synopt:{cmd:e(threshold)}}equivalence threshold{p_end}
{synopt:{cmd:e(MS_critical)}}MS critical value{p_end}
{synopt:{cmd:e(RMS_critical)}}RMS critical value{p_end}
{synopt:{cmd:e(reject)}}1 if H0 rejected, 0 otherwise{p_end}

{p2col 5 20 24 2: When threshold not specified:}{p_end}
{synopt:{cmd:e(min_threshold)}}minimum equivalence threshold{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b_placebo)}}placebo coefficient vector (1 × no_placebos){p_end}
{synopt:{cmd:e(MS_lambda)}}MS values for each subsample (1 × nolambda){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"rmsequivtest"{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(idvar)}}individual variable name{p_end}
{synopt:{cmd:e(groupvar)}}group variable name{p_end}
{synopt:{cmd:e(timevar)}}time variable name{p_end}
{synopt:{cmd:e(placebo_names)}}names of placebo period variables{p_end}
{synopt:{cmd:e(properties)}}"b"{p_end}


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
Online: {helpb equivtest}, {helpb maxequivtest}, {helpb meanequivtest}, 
{helpb equivsim}, {helpb equivtest_plot}
{p_end}
