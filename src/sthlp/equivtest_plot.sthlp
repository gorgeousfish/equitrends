{smcl}
{* *! version 0.1.0}{...}
{vieweralsosee "equivtest" "help equivtest"}{...}
{vieweralsosee "maxequivtest" "help maxequivtest"}{...}
{vieweralsosee "meanequivtest" "help meanequivtest"}{...}
{vieweralsosee "rmsequivtest" "help rmsequivtest"}{...}
{vieweralsosee "equivsim" "help equivsim"}{...}
{viewerjumpto "Syntax" "equivtest_plot##syntax"}{...}
{viewerjumpto "Description" "equivtest_plot##description"}{...}
{viewerjumpto "Options" "equivtest_plot##options"}{...}
{viewerjumpto "Examples" "equivtest_plot##examples"}{...}
{viewerjumpto "Remarks" "equivtest_plot##remarks"}{...}
{viewerjumpto "References" "equivtest_plot##references"}{...}
{viewerjumpto "Authors" "equivtest_plot##authors"}{...}
{viewerjumpto "Also see" "equivtest_plot##alsosee"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{cmd:equivtest_plot} {hline 2}}Visualization of equivtest results{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:equivtest_plot}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt thresh:old(#)}}equivalence threshold for horizontal lines{p_end}
{synopt:{opt ci}}display confidence intervals (IU test only){p_end}
{synopt:{opt level(#)}}confidence level; default is {cmd:level(95)}{p_end}

{syntab:Plot elements}
{synopt:{opt con:nect}}connect points with lines{p_end}
{synopt:{opt noline}}suppress zero reference line{p_end}
{synopt:{opt nobase}}suppress base period reference point{p_end}
{synopt:{opt nothreshold}}suppress equivalence threshold lines{p_end}

{syntab:Point markers}
{synopt:{opt msymbol(symbolstyle)}}marker symbol for placebo coefficients{p_end}
{synopt:{opt msize(markersizestyle)}}marker size{p_end}
{synopt:{opt mcolor(colorstyle)}}marker color{p_end}

{syntab:Base period markers}
{synopt:{opt basemsymbol(symbolstyle)}}marker symbol for base period{p_end}
{synopt:{opt basemcolor(colorstyle)}}marker color for base period{p_end}

{syntab:Confidence interval style}
{synopt:{opt cilcolor(colorstyle)}}CI line color{p_end}
{synopt:{opt cilwidth(linewidthstyle)}}CI line width{p_end}

{syntab:Threshold line style}
{synopt:{opt threshlcolor(colorstyle)}}threshold line color{p_end}
{synopt:{opt threshlwidth(linewidthstyle)}}threshold line width{p_end}
{synopt:{opt threshlpattern(linepatternstyle)}}threshold line pattern{p_end}

{syntab:Titles}
{synopt:{opt title(string)}}graph title{p_end}
{synopt:{opt subtitle(string)}}graph subtitle{p_end}
{synopt:{opt xtitle(string)}}x-axis title{p_end}
{synopt:{opt ytitle(string)}}y-axis title{p_end}

{syntab:Axes}
{synopt:{opt xlabel(rule_or_values)}}x-axis labels{p_end}
{synopt:{opt ylabel(rule_or_values)}}y-axis labels{p_end}

{syntab:Other}
{synopt:{opt legend(contents)}}legend options{p_end}
{synopt:{opt note(string)}}graph note{p_end}
{synopt:{opt scheme(schemename)}}graph scheme{p_end}

{syntab:Output}
{synopt:{opt saving(filename)}}save graph to file{p_end}
{synopt:{opt replace}}replace existing file{p_end}
{synopt:{opt name(windowname)}}graph window name{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:equivtest_plot} creates a coefficient plot of pre-trend placebo coefficients
from {cmd:equivtest} results. The plot displays:

{phang2}
{bf:Point estimates} - Placebo coefficients for each pre-treatment period

{phang2}
{bf:Base period reference} - A reference point at (0, 0) for the base period

{phang2}
{bf:Confidence intervals} - Optional confidence intervals (only available for
IU test method)

{phang2}
{bf:Threshold lines} - Horizontal dashed lines at ±threshold

{phang2}
{bf:Zero line} - A horizontal reference line at y=0

{pstd}
The x-axis shows relative time (period - base_period), so the base period
appears at x=0.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt threshold(#)} specifies the equivalence threshold for drawing horizontal
reference lines at ±threshold. If not specified, the command uses the threshold
from {cmd:e(threshold)} if available, or {cmd:e(min_threshold)} as a fallback.
The threshold must be strictly positive (greater than 0).

{phang}
{opt ci} requests confidence intervals to be displayed. This option is only
available when the previous {cmd:equivtest} command used {cmd:type(max)}
{cmd:method(iu)}. For other test types (bootstrap, mean, RMS), this option
is ignored with a warning because period-specific standard errors are not
available.

{phang}
{opt level(#)} specifies the confidence level for confidence intervals.
The default is {cmd:level(95)}. The value must be between 0 and 100 (exclusive).

{dlgtab:Plot elements}

{phang}
{opt connect} connects the point estimates with lines. By default, points are
displayed without connecting lines.

{phang}
{opt noline} suppresses the horizontal reference line at y=0.

{phang}
{opt nobase} suppresses the base period reference point at (0, 0).

{phang}
{opt nothreshold} suppresses the horizontal dashed lines at ±threshold. Use this
option when you do not want to display the equivalence bounds on the plot.

{dlgtab:Point markers}

{phang}
{opt msymbol(symbolstyle)} specifies the marker symbol for placebo coefficients.
The default is {cmd:msymbol(O)} (circle). See {manhelp symbolstyle G-4}.

{phang}
{opt msize(markersizestyle)} specifies the marker size. The default is
{cmd:msize(medium)}. See {manhelp markersizestyle G-4}.

{phang}
{opt mcolor(colorstyle)} specifies the marker color. The default is
{cmd:mcolor(navy)}. See {manhelp colorstyle G-4}.

{dlgtab:Base period markers}

{phang}
{opt basemsymbol(symbolstyle)} specifies the marker symbol for the base period
reference point. The default is {cmd:basemsymbol(S)} (square).

{phang}
{opt basemcolor(colorstyle)} specifies the marker color for the base period.
The default is the same as {opt mcolor()}.

{dlgtab:Confidence interval style}

{phang}
{opt cilcolor(colorstyle)} specifies the color for confidence interval lines.
The default is {cmd:cilcolor(navy)}.

{phang}
{opt cilwidth(linewidthstyle)} specifies the width for confidence interval lines.
The default is {cmd:cilwidth(medium)}.

{dlgtab:Threshold line style}

{phang}
{opt threshlcolor(colorstyle)} specifies the color for threshold lines.
The default is {cmd:threshlcolor(red)}.

{phang}
{opt threshlwidth(linewidthstyle)} specifies the width for threshold lines.
The default is {cmd:threshlwidth(medium)}.

{phang}
{opt threshlpattern(linepatternstyle)} specifies the pattern for threshold lines.
The default is {cmd:threshlpattern(dash)}.

{dlgtab:Titles}

{phang}
{opt title(string)} specifies the graph title. The default is
"Pre-trend Placebo Coefficients".

{phang}
{opt subtitle(string)} specifies the graph subtitle. When not specified and an
equivalence threshold is available, an auto-generated subtitle is displayed that
describes the equivalence bound and its interpretation for the test type. For
example, the maximum test shows "Equivalence bound: ±δ = ±0.15 (all |β_l| < δ)",
while the mean and RMS tests include a "conservative visual" qualifier to
indicate that the ±threshold lines apply to the aggregate statistic, not to
individual coefficients.

{phang}
{opt xtitle(string)} specifies the x-axis title. The default is
"Period relative to treatment".

{phang}
{opt ytitle(string)} specifies the y-axis title. The default is "Coefficient".

{dlgtab:Axes}

{phang}
{opt xlabel(rule_or_values)} specifies x-axis labels. See {manhelp axis_label_options G-3}.

{phang}
{opt ylabel(rule_or_values)} specifies y-axis labels. See {manhelp axis_label_options G-3}.

{dlgtab:Other}

{phang}
{opt legend(contents)} specifies legend options. By default, the legend is
turned off. See {manhelp legend_options G-3}.

{phang}
{opt note(string)} specifies a note to appear at the bottom of the graph.
When not specified, an auto-generated note is displayed that includes:
(1) the test type, null hypothesis, and method;
(2) the test statistic value and rejection decision (when available).
For the mean and RMS tests, the note explicitly states that the equivalence
bound applies to the aggregate statistic (average or root mean square), not
to individual coefficients. Specifying this option overrides the auto-generated
note entirely.

{phang}
{opt scheme(schemename)} specifies the graph scheme. See {manhelp scheme G-4}.

{dlgtab:Output}

{phang}
{opt saving(filename)} saves the graph to the specified file in Stata's
.gph format.

{phang}
{opt replace} allows {opt saving()} to overwrite an existing file.

{phang}
{opt name(windowname)} specifies the name for the graph window.


{marker examples}{...}
{title:Examples}

{pstd}Setup and run equivtest{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. gen treat = (idcode > 2000)}{p_end}
{phang2}{cmd:. equivtest ln_wage, type(max) method(iu) id(idcode) group(treat) time(year) threshold(0.1)}{p_end}

{pstd}Basic plot{p_end}
{phang2}{cmd:. equivtest_plot}{p_end}

{pstd}Plot with confidence intervals{p_end}
{phang2}{cmd:. equivtest_plot, ci}{p_end}

{pstd}Plot with 90% confidence intervals{p_end}
{phang2}{cmd:. equivtest_plot, ci level(90)}{p_end}

{pstd}Plot with connected points{p_end}
{phang2}{cmd:. equivtest_plot, connect}{p_end}

{pstd}Plot with custom threshold{p_end}
{phang2}{cmd:. equivtest_plot, threshold(0.2)}{p_end}

{pstd}Plot with custom styling{p_end}
{phang2}{cmd:. equivtest_plot, ci msymbol(D) mcolor(blue) threshlcolor(green)}{p_end}

{pstd}Publication-quality plot{p_end}
{phang2}{cmd:. equivtest_plot, ci level(95) title("Pre-trend Analysis") ///}{p_end}
{phang2}{cmd:     subtitle("IU Method") xtitle("Relative Time") ///}{p_end}
{phang2}{cmd:     msymbol(O) mcolor(navy) threshlpattern(dash) ///}{p_end}
{phang2}{cmd:     scheme(s2color)}{p_end}

{pstd}Save graph to file{p_end}
{phang2}{cmd:. equivtest_plot, ci saving(pretrend_plot) replace}{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Confidence intervals:} Confidence intervals are only available when the
previous {cmd:equivtest} command used {cmd:type(max)} {cmd:method(iu)}. This is
because only the IU (Intersection-Union) method provides period-specific
standard errors. Bootstrap, mean, and RMS tests do not provide these standard
errors, so confidence intervals cannot be computed.

{pstd}
{bf:Threshold priority:} When determining which threshold to display, the
command uses the following priority:

{phang2}1. User-specified {opt threshold()} option{p_end}
{phang2}2. {cmd:e(threshold)} from the previous {cmd:equivtest} command{p_end}
{phang2}3. {cmd:e(min_threshold)} as a fallback{p_end}

{pstd}
{bf:Relative time:} The x-axis displays relative time, calculated as
(period - base_period). This means the base period always appears at x=0,
and earlier periods have negative values.

{pstd}
{bf:Base period reference:} The base period is displayed as a reference point
at (0, 0) with a different marker symbol (square by default) to distinguish
it from the placebo coefficients.

{pstd}
{bf:Interpreting threshold lines by test type:} The ±threshold horizontal lines
have different interpretations depending on the test type:

{phang2}{bf:Maximum test} ({cmd:type(max)}): The threshold lines represent ±δ.
The visual interpretation is exact: if all plotted coefficient points fall
within the ±δ band, then max|β_l| < δ, which is precisely the alternative
hypothesis. The plot directly corresponds to the test conclusion.{p_end}

{phang2}{bf:Mean test} ({cmd:type(mean)}): The threshold lines represent ±τ.
The test evaluates whether |β̄| < τ, where β̄ is the {it:average} of the
placebo coefficients. The lines provide a {it:conservative} visual check:
if all individual coefficients fall within ±τ, then the average must also
be within ±τ (since |β̄| ≤ max|β_l|). However, the test may reject H₀
even when some individual coefficients exceed ±τ, because positive and
negative coefficients can offset each other in the average.{p_end}

{phang2}{bf:RMS test} ({cmd:type(rms)}): The threshold lines represent ±ζ.
The test evaluates whether β_RMS < ζ, where β_RMS is the root mean square
of the placebo coefficients. As with the mean test, the lines provide a
conservative visual check: all points within ±ζ implies β_RMS < ζ, but
the test may reject H₀ even when some individual coefficients exceed ±ζ.{p_end}

{pstd}
{bf:Auto-generated annotations:} When the {opt subtitle()} and {opt note()}
options are not specified, {cmd:equivtest_plot} automatically generates
informative annotations based on the test type and stored results. The
subtitle describes the equivalence bound and its interpretation, while the
note reports the test statistic value and rejection decision. For the mean
and RMS tests, these annotations explicitly flag that the threshold lines
are a "conservative visual" to prevent misinterpretation. Specifying
{opt subtitle()} or {opt note()} overrides the corresponding auto-generated
text.


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
{helpb rmsequivtest}, {helpb equivsim}
{p_end}
