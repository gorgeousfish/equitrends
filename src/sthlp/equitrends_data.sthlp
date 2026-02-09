{smcl}
{* *! version 2.0.0  08feb2026}{...}
{vieweralsosee "[D] use" "help use"}{...}
{vieweralsosee "equitrends" "help equitrends"}{...}
{viewerjumpto "Syntax" "equitrends_data##syntax"}{...}
{viewerjumpto "Description" "equitrends_data##description"}{...}
{viewerjumpto "Options" "equitrends_data##options"}{...}
{viewerjumpto "Examples" "equitrends_data##examples"}{...}
{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{cmd:equitrends_data} {hline 2}}Load bundled example datasets for equitrends package{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:equitrends_data} [{it:datasetname}] [{cmd:,} {opt clear}]


{marker description}{...}
{title:Description}

{pstd}
{cmd:equitrends_data} loads example datasets bundled with the equitrends package.
The command locates the dataset from the locally installed package files,
so users do not need to know the installation path or change their working
directory.

{pstd}
If {it:datasetname} is not specified, the default dataset {bf:MonthlyPanel} is loaded.

{pstd}
The {cmd:.dta} extension may be included or omitted in {it:datasetname}.


{marker options}{...}
{title:Options}

{phang}
{opt clear} specifies that it is okay to replace the data in memory, even
    though the current data have not been saved to disk.


{title:Available Datasets}

{synoptset 20 tabbed}{...}
{synopt:{opt MonthlyPanel}}Di Tella and Schargrodsky (2004) crime data for empirical examples{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Load the default MonthlyPanel dataset:{p_end}
{phang2}{cmd:. equitrends_data, clear}{p_end}

{pstd}Explicitly specify the dataset name:{p_end}
{phang2}{cmd:. equitrends_data MonthlyPanel, clear}{p_end}


{title:Authors}

{pstd}
Cai Xuanyu ({browse "mailto:xuanyuCAI@outlook.com":xuanyuCAI@outlook.com})
{p_end}

{pstd}
Xu Wenli ({browse "mailto:wlxu@cityu.edu.mo":wlxu@cityu.edu.mo})
{p_end}

{pstd}
{browse "https://github.com/gorgeousfish/equitrends"}
{p_end}


{title:Also see}

{psee}
{space 2}Help:  {manhelp use D}, {help equitrends}
{p_end}
