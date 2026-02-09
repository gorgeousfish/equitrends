{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{vieweralsosee "equivtest_plot" "help equivtest_plot"}{...}
{viewerjumpto "Syntax" "equivsim##syntax"}{...}
{viewerjumpto "Description" "equivsim##description"}{...}
{viewerjumpto "Options" "equivsim##options"}{...}
{viewerjumpto "DGP Specification" "equivsim##dgp"}{...}
{viewerjumpto "Modes" "equivsim##modes"}{...}
{viewerjumpto "Examples" "equivsim##examples"}{...}
{viewerjumpto "Stored results" "equivsim##results"}{...}
{viewerjumpto "References" "equivsim##references"}{...}
{viewerjumpto "Authors" "equivsim##authors"}{...}
{viewerjumpto "Also see" "equivsim##alsosee"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:equivsim} {hline 2}}Generate panel data for equivalence testing simulation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:equivsim}
{cmd:,}
{opt n(#)}
{opt preperiods(#)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt n(#)}}number of individuals{p_end}
{synopt:{opt preperiods(#)}}number of pre-treatment periods (T){p_end}

{syntab:DGP Specification}
{synopt:{opt dgp(string)}}DGP type: {cmd:spherical}, {cmd:ar1}, or {cmd:ar3}; default is {cmd:ar3}{p_end}
{synopt:{opt phi(numlist)}}AR coefficients; default depends on {opt dgp()}{p_end}
{synopt:{opt burnins(#)}}burn-in periods for AR process; default is {cmd:burnins(100)}{p_end}

{syntab:Error Variance}
{synopt:{opt sd(#)}}base standard deviation; default is {cmd:sd(1)}{p_end}
{synopt:{opt het(#)}}heteroskedasticity parameter; default is {cmd:het(1)}{p_end}

{syntab:Treatment Effects}
{synopt:{opt beta(numlist)}}placebo effects for pre-treatment periods{p_end}
{synopt:{opt piatt(#)}}average treatment effect on treated; default is {cmd:piatt(0)}{p_end}

{syntab:Additional Effects}
{synopt:{opt dip(#)}}mean of Ashenfelter's dip shock distribution{p_end}
{synopt:{opt trend(#)}}linear trend coefficient for treated group{p_end}

{syntab:Fixed Effects}
{synopt:{opt eta(numlist)}}individual fixed effects (length n){p_end}
{synopt:{opt lambda(numlist)}}time fixed effects (length T+2){p_end}

{syntab:Mode}
{synopt:{opt rcompat}}R compatibility mode{p_end}

{syntab:Data Management}
{synopt:{opt clear}}clear existing data{p_end}
{synopt:{opt replace}}replace existing data{p_end}
{synopt:{opt seed(#)}}random number seed for reproducibility{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:equivsim} generates balanced panel data for Monte Carlo
simulation studies of equivalence testing methods in difference-in-differences settings.

{pstd}
The generated dataset contains the following variables:

{p2colset 9 20 22 2}{...}
{p2col:{cmd:id}}individual identifier (1 to n){p_end}
{p2col:{cmd:period}}time period (1 to T+2){p_end}
{p2col:{cmd:Y}}outcome variable{p_end}
{p2col:{cmd:G}}treatment group indicator (0=control, 1=treated){p_end}
{p2colreset}{...}

{pstd}
The panel structure consists of T pre-treatment periods, one base period (T+1),
and one post-treatment period (T+2).


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt n(#)} specifies the number of individuals in the panel. Must be a positive integer.

{phang}
{opt preperiods(#)} specifies the number of pre-treatment periods T. Must be a positive integer.
The total number of periods is T+2 (T pre-treatment + 1 base + 1 post-treatment).

{dlgtab:DGP Specification}

{phang}
{opt dgp(string)} specifies the data generating process type. Options are:

{phang2}
{cmd:spherical} - i.i.d. N(0, sigma^2) errors with no serial correlation and 
homoskedastic variance (het is ignored).

{phang2}
{cmd:ar1} - AR(1) process with phi=0.5 and heteroskedasticity.

{phang2}
{cmd:ar3} - AR(3) process with phi=(0.5, 0.3, 0.1) and heteroskedasticity.
This is the default and matches the DGP in Dette and Schumann (2024).

{phang}
{opt phi(numlist)} specifies the AR coefficients for the error process. The length of the
numlist determines the AR order p. Default depends on {opt dgp()}: 0 for spherical,
0.5 for AR(1), and (0.5, 0.3, 0.1) for AR(3). User-specified values override the default.
The coefficients must satisfy the stationarity condition (all eigenvalues of the 
companion matrix have modulus less than 1).

{phang}
{opt burnins(#)} specifies the number of burn-in periods for the AR process to reach
stationarity. Default is {cmd:burnins(100)}.

{dlgtab:Error Variance}

{phang}
{opt sd(#)} specifies the base standard deviation for the error term. Must be positive.
Default is {cmd:sd(1)}. The actual error standard deviation depends on the group:
sigma_i = sd * (1 + het) for treated (G=1) and sigma_i = sd for control (G=0).

{phang}
{opt het(#)} specifies the heteroskedasticity parameter. Must be non-negative.
Default is {cmd:het(1)}. The treated group has standard deviation sd*(1+het), 
while the control group has standard deviation sd. With default values (sd=1, het=1),
the treated group has variance 4 and the control group has variance 1, matching
the paper's specification. For {cmd:dgp(spherical)}, het is ignored (homoskedastic).

{dlgtab:Treatment Effects}

{phang}
{opt beta(numlist)} specifies the placebo effects for pre-treatment periods.

{pstd}
{err:{bf:Important:}} The length of {opt beta()} must equal {opt preperiods()}. 
For example, if {cmd:preperiods(4)} is specified, {cmd:beta()} must contain 
exactly 4 values (e.g., {cmd:beta(0 0 0 0)}).

{pstd}
Default is all zeros (no placebo effects).

{phang}
{opt piatt(#)} specifies the average treatment effect on the treated in the
post-treatment period. Default is {cmd:piatt(0)}.

{dlgtab:Additional Effects}

{phang}
{opt dip(#)} specifies the mean of the Ashenfelter's dip shock distribution.
For each treated individual at the base period (T+1), a shock V_i ~ N(dip, 1)
is added to the error term.

{phang}
{opt trend(#)} specifies the linear trend coefficient. The trend effect equals
psi * t for the treated group, where t is the period number.

{dlgtab:Fixed Effects}

{phang}
{opt eta(numlist)} specifies user-defined individual fixed effects. The length
must equal n. If not specified, eta_i ~ N(0,1) in paper mode or eta_i = 0 in rcompat mode.

{phang}
{opt lambda(numlist)} specifies user-defined time fixed effects. The length must
equal T+2. If not specified, lambda_t ~ N(0,1) in paper mode or lambda_t = 0 in rcompat mode.

{dlgtab:Mode}

{phang}
{opt rcompat} enables R compatibility mode. See {help equivsim##modes:Modes} below.

{dlgtab:Data Management}

{phang}
{opt clear} clears existing data before generating new data.

{phang}
{opt replace} replaces existing data if present.

{phang}
{opt seed(#)} sets the random number seed for reproducibility.


{marker dgp}{...}
{title:DGP Specification}

{pstd}
The outcome variable Y is generated according to:

{p 8 8 2}
Y_it = eta_i + lambda_t + beta_t * G_i * D_t + piatt * G_i * D_{T+2} + dip * G_i * D_{T+1} + trend * G_i * t + u_it

{pstd}
where:

{p2colset 9 25 27 2}{...}
{p2col:{it:eta_i}}individual fixed effect{p_end}
{p2col:{it:lambda_t}}time fixed effect{p_end}
{p2col:{it:beta_t}}placebo effect in period t{p_end}
{p2col:{it:G_i}}treatment group indicator{p_end}
{p2col:{it:D_t}}indicator for pre-treatment period t{p_end}
{p2col:{it:D_{T+1}}}indicator for base period{p_end}
{p2col:{it:D_{T+2}}}indicator for post-treatment period{p_end}
{p2col:{it:piatt}}average treatment effect on treated{p_end}
{p2col:{it:dip}}Ashenfelter's dip mean (V_i ~ N(dip, 1)){p_end}
{p2col:{it:trend}}linear trend coefficient{p_end}
{p2col:{it:u_it}}AR(p) error process{p_end}
{p2colreset}{...}

{pstd}
The error term u_it follows an AR(p) process:

{p 8 8 2}
u_it = phi_1 * u_{i,t-1} + phi_2 * u_{i,t-2} + ... + phi_p * u_{i,t-p} + sigma_i * epsilon_it

{pstd}
where epsilon_it ~ N(0,1) and sigma_i = sd*(1+het) for treated (G=1) and sigma_i = sd for
control (G=0). With default values (sd=1, het=1), the treated group has SD=2 and control 
has SD=1. For {cmd:dgp(spherical)}, errors are homoskedastic (het is forced to 0).


{marker modes}{...}
{title:Modes}

{pstd}
{cmd:equivsim} supports two modes:

{pstd}
{bf:Default Mode (Paper DGP):}

{p2colset 9 30 32 2}{...}
{p2col:DGP Type}AR(3) (can be changed via {opt dgp()}){p_end}
{p2col:AR Process}phi = (0.5, 0.3, 0.1){p_end}
{p2col:Heteroskedasticity}σ_i = sd*(1+het*G_i); with defaults: treated SD=2, control SD=1{p_end}
{p2col:Treatment Assignment}Random: Pr(G=1) = 0.5{p_end}
{p2col:Fixed Effects}η_i ~ N(0,1), λ_t ~ N(0,1){p_end}
{p2colreset}{...}

{pstd}
{bf:R Compatibility Mode ({opt rcompat}):}

{p2colset 9 30 32 2}{...}
{p2col:DGP Type}AR(1) by default (unless {opt dgp()} specified){p_end}
{p2col:AR Process}phi = 0.5{p_end}
{p2col:Heteroskedasticity}Homoskedastic (both groups have SD=sd){p_end}
{p2col:Treatment Assignment}Deterministic: G=1 if mod(id,2)==0{p_end}
{p2col:Fixed Effects}η_i = 0, λ_t = 0{p_end}
{p2colreset}{...}

{pstd}
Use {opt rcompat} mode for cross-validation with the R package EquiTrends.

{pstd}
{bf:DGP Types Summary:}

{p2colset 9 20 22 2}{...}
{p2col:{cmd:spherical}}i.i.d. errors, no serial correlation, homoskedastic{p_end}
{p2col:{cmd:ar1}}AR(1) with phi=0.5, heteroskedastic{p_end}
{p2col:{cmd:ar3}}AR(3) with phi=(0.5, 0.3, 0.1), heteroskedastic (default){p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Examples}

{pstd}Basic usage with default parameters:{p_end}
{phang2}{cmd:. equivsim, n(100) preperiods(4) clear}

{pstd}Generate data with specific placebo effects:{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(4) beta(0.1 0.2 0.3 0.4) clear}

{pstd}Generate data with treatment effect:{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(4) beta(0 0 0 0) piatt(0.5) clear}

{pstd}Generate data with Ashenfelter's dip:{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(4) dip(0.3) clear}

{pstd}Generate data with linear trend:{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(10) trend(0.025) clear}

{pstd}R compatibility mode for cross-validation:{p_end}
{phang2}{cmd:. equivsim, n(500) preperiods(4) rcompat clear seed(12345)}

{pstd}Custom AR process:{p_end}
{phang2}{cmd:. equivsim, n(200) preperiods(20) phi(0.6) clear}

{pstd}Monte Carlo simulation example:{p_end}
{phang2}{cmd:. forvalues i = 1/1000 {c -(}}{p_end}
{phang2}{cmd:.     quietly equivsim, n(500) preperiods(4) clear seed(`i')}{p_end}
{phang2}{cmd:.     // Run estimation and store results}{p_end}
{phang2}{cmd:. {c )-}}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:equivsim} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of individuals{p_end}
{synopt:{cmd:r(preperiods)}}number of pre-treatment periods{p_end}
{synopt:{cmd:r(n_obs)}}total number of observations{p_end}
{synopt:{cmd:r(n_treated)}}number of treated individuals{p_end}
{synopt:{cmd:r(n_control)}}number of control individuals{p_end}
{synopt:{cmd:r(seed)}}random seed used (if specified){p_end}
{p2colreset}{...}


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
{helpb rmsequivtest}, {helpb equivtest_plot}
{p_end}
