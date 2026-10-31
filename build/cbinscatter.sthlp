{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "cbinscatter##syntax"}{...}
{viewerjumpto "Description" "cbinscatter##description"}{...}
{viewerjumpto "Options" "cbinscatter##options"}{...}
{viewerjumpto "Examples" "cbinscatter##examples"}{...}
{viewerjumpto "Stored results" "cbinscatter##results"}{...}
{title:Title}

{phang}
{bf:cbinscatter} {hline 2} C-accelerated binned scatter plots


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cbinscatter}
{it:yvar}
{it:xvar}
{ifin}
{weight}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt nq:uantiles(#)}}number of bins (default 20){p_end}
{synopt:{opt c:ontrols(varlist)}}control variables to partial out via OLS{p_end}
{synopt:{opt a:bsorb(varlist)}}fixed effects to absorb via HDFE{p_end}
{synopt:{opt by(varname)}}create separate series by group{p_end}
{synopt:{opt meth:od(method)}}residualization method: {opt classic} (default) or {opt binsreg}{p_end}

{syntab:Line Fitting}
{synopt:{opt line:type(type)}}fit line type: {opt none}, {opt linear} (default), {opt qfit}, {opt cubic}{p_end}

{syntab:Data Options}
{synopt:{opt discrete}}treat x as discrete (one bin per unique value){p_end}

{synopt:{opt save:data(filename)}}save bin data to file{p_end}

{syntab:Graph Options}
{synopt:{opt nog:raph}}suppress graph output{p_end}
{synopt:{opt title(string)}}graph title{p_end}
{synopt:{opt yti:tle(string)}}y-axis title (default: yvar name){p_end}
{synopt:{opt xti:tle(string)}}x-axis title (default: xvar name){p_end}
{synopt:{opt legend(string)}}legend options{p_end}
{synopt:{opt colors(string)}}colors for series{p_end}
{synopt:{opt msymbols(string)}}marker symbols for series{p_end}
{synopt:{opt mlabels(string)}}marker labels for series{p_end}
{synopt:{it:twoway_options}}additional twoway graph options{p_end}

{syntab:Reporting}
{synopt:{opt reportreg}}report underlying regression{p_end}
{synopt:{opt v:erbose}}display progress and elapsed times{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{opt aweight}s, {opt fweight}s, {opt pweight}s, and {opt iweight}s are allowed;
see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cbinscatter} divides x into quantile bins and plots the mean of y against
the mean of x within each bin. It can adjust for controls and fixed effects,
plot separate series by group, and overlay a polynomial fit. Computation runs
in C; Stata draws the graph.

{pstd}
Tied x values stay in the same bin, so bins need not contain equal numbers of
observations and some requested bins may be empty. Weights affect cutpoints,
bin means, and fitted lines.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt nquantiles(#)} specifies the requested number of quantile bins.
The default is 20. Must be between 2 and 1000.

{phang}
{opt controls(varlist)} specifies control variables to partial out from both
x and y before computing bins. The residualization is performed via OLS.

{phang}
{opt absorb(varlist)} specifies categorical variables representing fixed effects
to absorb from both x and y. Can be combined with {opt controls()}. The variables may be
numeric or string; their values only label groups, so negative, fractional and
large codes are allowed.

{phang}
{opt by(varname)} creates separate scatter series for each level of the
specified variable. As in {cmd:binscatter}, the bins are the quantiles of x in
the whole sample, and each group's points are its means within those bins;
with {cmd:method(binsreg)} each group is binned on its own x, as in
{cmd:binsreg}. Each group gets its own optional fit line.

{phang}
{opt method(method)} specifies the residualization method when controls or
fixed effects are specified. This follows the methodology discussion in
Cattaneo, Crump, Farrell, and Feng (2024), "On Binscatter".

{phang2}{opt classic} (the default) residualizes both x and y on the controls/FE,
then creates bins based on the residualized x values. This is the traditional
binscatter approach.{p_end}

{phang2}{opt binsreg} creates bins based on the raw (non-residualized) x values,
then computes the conditional mean of y within each bin after partialling out
controls/FE. Controls and fixed effects enter the regression jointly with bin indicators.{p_end}

{pstd}
Without controls, fixed effects, or {opt by()}, both methods produce the
same bin means.

{pstd}
With {opt method(binsreg)}, singleton fixed-effect groups are retained in the
bins and in the sample used to center adjusted means.


{dlgtab:Line Fitting}

{phang}
{opt linetype(type)} specifies the type of fit line to overlay on the scatter.
Options are:

{phang2}{opt none} - no fit line{p_end}
{phang2}{opt linear} (or {opt lfit} or {opt line}) - linear fit (default){p_end}
{phang2}{opt qfit} or {opt quadratic} - quadratic fit{p_end}
{phang2}{opt cubic} - cubic fit{p_end}
{phang2}{opt connect} - connect bin means with lines (no fit){p_end}

{pstd}
Fit lines use the underlying observations.
With {opt method(classic)} and {opt controls()} or {opt absorb()}, the line is
fit to the residualized x and y after adding their sample means back, as
{cmd:binscatter} does, so line and points share one scale and the coefficients
match {cmd:binscatter}'s {cmd:e(y1_coefs)}. With {opt method(binsreg)}, the
polynomial is fit jointly with the controls and absorbed effects and evaluated
at the controls' means, as {cmd:binsreg}'s {opt polyreg()} does.

{dlgtab:Data Options}

{phang}
{opt discrete} treats x as a discrete variable, creating one bin for each
unique value of x instead of quantile-based bins. With {opt method(binsreg)},
the bin means are adjusted for {opt controls()} and {opt absorb()} as with
quantile bins.

{phang}
{opt savedata(filename)} saves the bin data (bin means, counts, etc.) to
the specified Stata data file, replacing an existing file of that name.
Use this option for aggregate bin data; {opt genxq()} is not supported.

{dlgtab:Graph Options}

{phang}
{opt nograph} suppresses the graph output. Useful when you only need the
computed bin statistics in {cmd:e()}.

{phang}
{opt title(string)}, {opt ytitle(string)}, {opt xtitle(string)} specify
graph and axis titles. Defaults for axis titles are the variable names.

{phang}
{opt colors(string)} specifies colors for the scatter series when using
{opt by()}. Default colors match those used by {cmd:binscatter}.

{phang}
{opt msymbols(string)} specifies marker symbols for the scatter series.
The {it:i}th symbol applies to the {it:i}th series, as in {cmd:binscatter}.
Default is filled circles (O).

{phang}
{opt mlabels(string)} specifies marker labels for the scatter series. The
{it:i}th word labels every bin marker of the {it:i}th series; use quotes for
labels containing spaces. {cmd:binscatter} has no equivalent option.

{phang}
{it:twoway_options} any other options are passed through to the underlying
{cmd:twoway} graph command.

{dlgtab:Reporting}

{phang}
{opt reportreg} displays the fitted polynomial coefficients, R-squared, and
number of observations for each group.

{phang}
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.


{title:Control and weight handling}

{pstd}
Weight expressions are evaluated once. Missing weights are excluded, selected
weights must be positive, and fweights must be integers. Factor-variable controls
and interactions are materialized before adjustment. Redundant control columns
are omitted; a failed adjustment returns an error. Fixed effects are removed by iterative projection.
{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Basic binned scatter plot:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbinscatter price mpg}{p_end}

{pstd}With 10 bins and no fit line:{p_end}
{phang2}{cmd:. cbinscatter price mpg, nquantiles(10) linetype(none)}{p_end}

{pstd}Control for other variables:{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight length)}{p_end}

{pstd}With fixed effects:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. cbinscatter ln_wage tenure, absorb(idcode year)}{p_end}

{pstd}Separate series by group:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. cbinscatter price mpg, by(foreign) nquantiles(10)}{p_end}

{pstd}Quadratic fit with custom titles:{p_end}
{phang2}{cmd:. cbinscatter price mpg, linetype(qfit) title("Price vs Fuel Efficiency") ytitle("Price (USD)") xtitle("Miles per Gallon")}{p_end}

{pstd}With weights:{p_end}
{phang2}{cmd:. cbinscatter price mpg [aw=weight], nquantiles(15)}{p_end}

{pstd}Save bin data without displaying graph:{p_end}
{phang2}{cmd:. cbinscatter price mpg, nograph savedata(mybins)}{p_end}

{pstd}Using the binsreg method (Cattaneo et al.):{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(binsreg)}{p_end}

{pstd}Compare classic vs binsreg methods:{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(classic) savedata(classic_bins)}{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(binsreg) savedata(binsreg_bins)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cbinscatter} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations used{p_end}
{synopt:{cmd:e(N_dropped)}}observations dropped due to missing values{p_end}
{synopt:{cmd:e(nquantiles)}}requested number of quantile bins{p_end}
{synopt:{cmd:e(num_groups)}}number of by-groups{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cbinscatter}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(xvar)}}x variable{p_end}
{synopt:{cmd:e(controls)}}control variables{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(by)}}by variable{p_end}
{synopt:{cmd:e(linetype)}}fit line type{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bindata)}}matrix of bin statistics: by_group, bin_id, x_mean, y_mean, n_obs{p_end}
{synopt:{cmd:e(coefs)}}fit line coefficients, one row per by-group: constant, x, x^2, x^3 (for a polynomial fit){p_end}
{synopt:{cmd:e(fit_stats)}}R-squared and number of observations of each fit{p_end}


{title:References}

{phang}
Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024. On Binscatter.
{it:American Economic Review} 114(5): 1488-1514.
{browse "https://doi.org/10.1257/aer.20221576"}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command is inspired by and provides similar functionality to {cmd:binscatter}
by Michael Stepner. The original {cmd:binscatter} package is available at
{browse "https://github.com/michaelstepner/binscatter"}.


{title:Also see}

{psee}
Online: {browse "https://github.com/michaelstepner/binscatter":binscatter} (if installed),
{help ctools}, {help creghdfe}, {help twoway}
{p_end}
