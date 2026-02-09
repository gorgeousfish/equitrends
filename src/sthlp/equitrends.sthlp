{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{vieweralsosee "equivsim" "help equivsim"}{...}
{vieweralsosee "equivtest_plot" "help equivtest_plot"}{...}
{vieweralsosee "equitrends_data" "help equitrends_data"}{...}
{viewerjumpto "Description" "equitrends##description"}{...}
{viewerjumpto "Commands" "equitrends##commands"}{...}
{viewerjumpto "References" "equitrends##references"}{...}
{viewerjumpto "Authors" "equitrends##authors"}{...}
{viewerjumpto "Also see" "equitrends##alsosee"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{cmd:equitrends} {hline 2}}Equivalence testing for pre-trends in difference-in-differences estimation{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
The {cmd:equitrends} package implements equivalence tests for assessing the 
plausibility of the parallel trends assumption (PTA) in Difference-in-Differences 
(DiD) estimation. The methods are based on Dette and Schumann (2024).

{pstd}
Traditional pre-tests for parallel trends test the null hypothesis of "no 
differences in trends" between treatment and control groups. However, failure 
to reject this null does not imply the absence of differences. The equivalence 
testing approach reverses the burden of proof by testing the null hypothesis 
that trend differences are {it:non-negligible}, requiring the data to provide 
evidence {it:in favor} of similar trends.

{pstd}
The package provides three types of equivalence tests:

{p 8 12 2}
{bf:Maximum test}: Tests whether the maximum absolute placebo coefficient 
exceeds a threshold delta. This is useful when all pre-treatment periods 
should show similar trends.
{p_end}

{p 8 12 2}
{bf:Mean test}: Tests whether the average of placebo coefficients exceeds 
a threshold tau. This is appropriate when violations of parallel trends 
are expected to be of the same sign (e.g., monotone violations).
{p_end}

{p 8 12 2}
{bf:RMS test}: Tests whether the root mean square of placebo coefficients 
exceeds a threshold zeta. This provides a balance between the maximum and 
mean tests and is robust to cancellation effects.
{p_end}


{marker commands}{...}
{title:Commands}

{pstd}
The package provides the following commands:

{synoptset 18 tabbed}{...}
{synopthdr:Command}
{synoptline}
{syntab:Testing}
{synopt:{helpb equivtest}}Unified interface for all three equivalence tests{p_end}
{synopt:{helpb maxequivtest}}Maximum absolute placebo coefficient test{p_end}
{synopt:{helpb meanequivtest}}Mean placebo coefficient test{p_end}
{synopt:{helpb rmsequivtest}}Root mean square placebo coefficient test{p_end}

{syntab:Utilities}
{synopt:{helpb equivsim}}Monte Carlo simulation for power analysis{p_end}
{synopt:{helpb equivtest_plot}}Visualization of placebo coefficients and thresholds{p_end}
{synopt:{helpb equitrends_data}}Load bundled example datasets{p_end}
{synoptline}


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
{helpb rmsequivtest}, {helpb equivsim}, {helpb equivtest_plot}, {helpb equitrends_data}
{p_end}
