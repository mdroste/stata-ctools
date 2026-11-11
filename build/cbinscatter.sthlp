{smcl}
{* *! version 1.2.6 20261111}{...}
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
{it:yvar} [{it:yvar} ...]
{it:xvar}
{ifin}
{weight}
[{cmd:,} {it:options}]

{pstd}
As in {cmd:binscatter}, several {it:yvar}s may precede {it:xvar}. They share one
sample (every variable nonmissing), so the bins of {it:xvar} are the same; each
{it:yvar} gets its own binned means and fit line, and all are drawn in one graph:
for each by-group one series per {it:yvar}, the legend naming the {it:yvar}
(and the group). {opt savedata()} then saves every {it:yvar}'s bins with the
variables {cmd:y_index} and {cmd:depvar}, and the results are stored as
{cmd:e(bindata_y}{it:#}{cmd:)} and {cmd:e(coefs_y}{it:#}{cmd:)}, with
{cmd:e(num_depvars)} and {cmd:e(graphcmd)}.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt n:quantiles(#)}}number of bins (default 20){p_end}
{synopt:{opt c:ontrols(varlist)}}control variables to partial out via OLS{p_end}
{synopt:{opt a:bsorb(varlist)}}fixed effects to absorb via HDFE{p_end}
{synopt:{opt by(varname)}}create separate series by group{p_end}
{synopt:{opt meth:od(method)}}residualization method: {opt classic} (default) or {opt binsreg}{p_end}

{syntab:Line Fitting}
{synopt:{opt line:type(type)}}fit line type: {opt none}, {opt linear} (default), {opt qfit}, {opt cubic}{p_end}

{syntab:Data Options}
{synopt:{opt discrete}}treat x as discrete (one bin per unique value){p_end}
{synopt:{opt xq(varname)}}bins given by the values of {it:varname}{p_end}
{synopt:{opt gen:xq(newvar)}}store each observation's bin in {it:newvar}{p_end}
{synopt:{opt med:ians}}plot bin medians instead of means{p_end}
{synopt:{opt noa:ddmean}}do not add the means back to the residualized x and y{p_end}
{synopt:{opt rd(numlist)}}fit lines separately on each side of these x values{p_end}
{synopt:{opt randvar(varname)}}quantile cutpoints from observations with {it:varname} <= {opt randcut()}{p_end}
{synopt:{opt randcut(#)}}share (or {opt randvar()} cutoff) of the random subsample; default 1{p_end}
{synopt:{opt randn(#)}}size of the random subsample{p_end}

{synopt:{opt save:data(filename)}}save bin data to file{p_end}
{synopt:{opt replace}}let {opt savedata()} and {opt savegraph()} overwrite files{p_end}

{syntab:Graph Options}
{synopt:{opt nog:raph}}suppress graph output{p_end}
{synopt:{opt title(string)}}graph title{p_end}
{synopt:{opt yti:tle(string)}}y-axis title (default: yvar name){p_end}
{synopt:{opt xti:tle(string)}}x-axis title (default: xvar name){p_end}
{synopt:{opt legend(string)}}legend options{p_end}
{synopt:{opt colors(string)}}colors for series{p_end}
{synopt:{opt mc:olors(string)}}marker colors for series{p_end}
{synopt:{opt lc:olors(string)}}line colors for series{p_end}
{synopt:{opt savegraph(filename)}}save the graph ({cmd:.gph}) or export it{p_end}
{synopt:{opt m:symbols(string)}}marker symbols for series{p_end}
{synopt:{opt mlabels(string)}}marker labels for series{p_end}
{synopt:{it:twoway_options}}additional twoway graph options{p_end}

{syntab:binscatter legacy options}
{synopt:{opt nbins(#)}, {opt create_xq}, {opt x_q(varname)}, {opt symbols(string)}, {opt unique(string)}, {cmd:linetype(noline)}}accepted as {cmd:binscatter} accepts them; see {help cbinscatter##legacy:Legacy options}{p_end}

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
The default is 20 and the smallest is 2; the bins of all by() groups must fit
in one matrix ({cmd:c(max_matdim)} rows). {cmd:binscatter} stops with r(198)
("32,740 invalid name") when it forms more than 100 quantile bins; {cmd:cbinscatter}
forms them. As with {cmd:binscatter}, when
the estimation sample has no more than {it:#} distinct values of x (with
{opt controls()} or {opt absorb()}, of the residualized x plus its mean as
{cmd:binscatter} stores it, in float precision), x is treated as discrete: each
by-group gets one bin per value it contains. This applies to the default method
only.

{phang}
{opt controls(varlist)} specifies control variables to partial out from both
x and y before computing bins. The residualization is performed via OLS.
Factor variables and time-series operators are allowed.

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
fixed effects are specified. See Cattaneo, Crump, Farrell, and Feng (2024),
"On Binscatter", for a comparison of these methods.

{phang2}{opt classic} (the default) residualizes both x and y on the
controls/FE, then creates bins based on the residualized x values. This is the
method used by {cmd:binscatter}.{p_end}

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
Without controls in the fit, when x has no more distinct values than the
polynomial's order, a lower power is omitted, as {cmd:regress} omits it in
{cmd:binscatter}'s fit: with two values, {cmd:linetype(qfit)} keeps x^2 and
omits x.

{dlgtab:Data Options}

{phang}
{opt discrete} treats x as a discrete variable, creating one bin for each
unique value of x instead of quantile-based bins. With {opt method(binsreg)},
the bin means are adjusted for {opt controls()} and {opt absorb()} as with
quantile bins. As in {cmd:binscatter}, x (the residualized x with
{opt method(classic)} and {opt controls()} or {opt absorb()}, as {cmd:binscatter}
stores it in float precision) may have at most
half as many distinct values as the sample has observations; otherwise the
command exits with r(134). The bins must also fit in a Stata matrix
({cmd:c(max_matdim)} rows, over all by-groups), or the command exits with
r(503).

{phang}
{opt xq(varname)} takes the bins from {it:varname}, which must hold positive
integers (in every observation, as {cmd:binscatter} checks it): one bin per
value in the sample, in ascending order of the values, plotted at the mean of
x (the residualized x with {opt controls()} or {opt absorb()}) and of y within
it. Observations with {it:varname} missing are excluded. {cmd:e(bindata)} numbers
the bins by their values. As in {cmd:binscatter}, {opt xq()} cannot be combined
with {opt discrete}, {opt genxq()} or {opt nquantiles()}, more distinct values
than half the sample exit with r(134), and with {opt controls()} or
{opt absorb()} a warning recalls that the bins should be bins of the
residualized x. Not with {cmd:method(binsreg)}.

{phang}
{opt genxq(newvar)} stores the quantile bin of each observation in the
sample, as {cmd:binscatter} numbers it ({cmd:1} + the number of cutpoints
below the observation's x, so tied cutpoints skip numbers). With {opt by()}
the bins are those of the pooled x. The variable is a {cmd:byte} (an
{cmd:int} above 100 bins), missing outside the sample, labeled
"{it:#} quantiles of {it:xvar}" ({cmd:binscatter} names the temporary
residual variable there under {opt controls()} or {opt absorb()}). When x has
no more distinct values than {opt nquantiles()} it is treated as discrete and,
as in {cmd:binscatter}, a note says {opt genxq()} is ignored and no variable is
created. Not with {opt xq()} or {opt discrete}. With {cmd:method(binsreg)} it
holds the bins of that method (each by-group's own quantiles).

{phang}
{opt medians} plots the median of x and of y within each bin instead of the
means, as {cmd:binscatter} does: {cmd:_pctile}'s 50th percentile (the average
of the two middle values for an even count), weighted with the weights. For
discrete bins x is the bin's value. The fit lines are unchanged. Not with
{cmd:method(binsreg)}.

{phang}
{opt noaddmean} leaves the residualized x and y without their sample means
when {opt controls()} or {opt absorb()} are given, as in {cmd:binscatter}; the
bins, points and fit lines are then on the residuals' scale. Not with
{cmd:method(binsreg)}.

{phang}
{opt rd(numlist)} (ascending) draws a dashed vertical line at each value and
fits the line separately in each interval: x <= the first value, between
consecutive values (above the lower, up to the higher), and above the last,
on the x that is binned (residualized with {opt controls()} or
{opt absorb()}). Each line runs from the first bin (or the break point) to the
break point (or the last bin). The bins themselves are not split, as in
{cmd:binscatter}. An interval without observations has no line (a message says
so); one with a single observation exits with r(2001), as {cmd:regress} does.
{cmd:e(coefs)} and {cmd:e(fit_stats)} then have one row per group and interval
(group by group), and {cmd:e(rdintervals)} holds the intervals' bounds. Not
with {cmd:method(binsreg)}.

{phang}
{opt randvar(varname)}, {opt randcut(#)} and {opt randn(#)} take the quantile
cutpoints from a random subsample, as {cmd:binscatter}'s {cmd:fastxtile} does;
every observation is then binned on them. With {opt randvar()}, the subsample
is the observations with {it:varname} <= {opt randcut()} (default 1). Without
it, {opt randcut()} must be in (0,1), and a uniform draw is made for every
observation as {cmd:binscatter} makes it ({cmd:runiform()} for the data sorted
by x after the sample), so the same {help set seed:seed} gives the same
subsample. {opt randn(#)} sets {opt randcut()} to # divided by the sample size
({it:#} larger than the sample uses the whole sample, with a note);
{opt randvar()} must then lie in [0,1]. More bins than the subsample's
observations plus one exit with r(198). The options do nothing with
{opt discrete} or {opt xq()}, or when x is treated as discrete (no draw is then
made, though, unlike {cmd:binscatter}, invalid values are still refused). With weights and ties in x, which tied observations are drawn may
differ from {cmd:binscatter}. Not with {cmd:method(binsreg)} or
{opt nofastxtile}.

{phang}
{opt savedata(filename)} saves the bin data. Given without an extension (and
with the default method), it writes {cmd:binscatter}'s two files,
{it:filename}{cmd:.csv} with the scatter points ({it:xvar} of each by() group,
then each {it:yvar} for each group, {cmd:binscatter}'s rows and number format)
and {it:filename}{cmd:.do}, which reads the CSV with {cmd:insheet}, labels the
variables and redraws the graph with {cmd:binscatter}'s {cmd:twoway} command,
and also {it:filename}{cmd:.dta} with the bins (by_group, bin_id, x_mean,
y_mean, n_obs; with several y also y_index and depvar). A file name with an
extension (for example {cmd:mybins.dta}) is that Stata data file alone, and so
is any name with {cmd:method(binsreg)}, whose bins are not {cmd:binscatter}'s
({cmd:binsreg}'s own {opt savedata()} layout is not written). An existing file
is replaced only with {opt replace}; otherwise the command exits with r(602)
before computing anything.

{dlgtab:Graph Options}

{phang}
{opt nograph} suppresses the graph. Bin statistics are still stored in
{cmd:e()}.

{phang}
{opt title()}, {opt ytitle()}, {opt xtitle()} specify graph and axis titles.
Defaults for axis titles are the variable names. As in {cmd:binscatter}, they
are passed to {cmd:twoway} as typed, after the defaults, so suboptions such as
{cmd:title("Text", size(small))}, several lines ({cmd:title("One" "Two")}) and
empty titles ({cmd:ytitle("")}) work as in {help title_options}.

{phang}
{opt legend()} specifies legend options, merged after the default legend as in
{cmd:binscatter}: with one series the default is {cmd:legend(off)}, which
{cmd:legend(on)} overrides. With {opt by()}, each key reads
{it:byvar}{cmd:=}{it:value}, the value's label when it has one.

{phang}
{opt colors(string)} specifies colors for the scatter series when using
{opt by()}. Default colors match those used by {cmd:binscatter}: with one
series (and no {cmd:linetype(connect)}) the first color marks the bins and the
second draws the line; otherwise the {it:i}th series takes the {it:i}th color.

{phang}
{opt mcolors(string)} and {opt lcolors(string)} give the marker colors and the
line colors (of the fit lines and of {cmd:connect}) themselves, the {it:i}th
for the {it:i}th series, as in {cmd:binscatter}; either defaults to the rule
above.

{phang}
{opt savegraph(filename)} saves the graph after drawing it, as
{cmd:binscatter} does: with {cmd:graph save} for a {cmd:.gph} file or a name
without an extension, otherwise with {cmd:graph export} (the extension gives
the format). An existing file is replaced only with {opt replace}; otherwise
the command exits with r(602) before computing anything. Not with
{opt nograph}.

{phang}
{opt msymbols(string)} (abbreviation {opt m()}, as in {cmd:binscatter})
specifies marker symbols for the scatter series.
The {it:i}th symbol applies to the {it:i}th series, as in {cmd:binscatter}.
Default is filled circles (O).

{phang}
{opt mlabels(string)} specifies marker labels for the scatter series. The
{it:i}th word labels every bin marker of the {it:i}th series; use quotes for
labels containing spaces. {cmd:binscatter} has no equivalent option.

{phang}
{it:twoway_options} any other options are passed through to the underlying
{cmd:twoway} graph command. With {opt nograph} they are still checked, so an
invalid option is an error, as it is with the graph; {opt name()} and
{opt saving()} then store and write nothing.

{phang}
{opt replace} lets {opt savedata()} and {opt savegraph()} overwrite existing
files, as in {cmd:binscatter}. {opt nofastxtile} is accepted and changes
nothing (quantiles are always exact), except that, as in {cmd:binscatter}, it
cannot be combined with {opt randvar()}, {opt randcut()} or {opt randn()}.
Where {cmd:binscatter} displays an error and stops without an error code (for
example {opt xq()} with {opt discrete}), {cmd:cbinscatter} exits with r(198).

{pstd}
All of these options also work with several {it:yvar}s (the bins, the random
subsample and {opt genxq()} are common to all of them), with {opt by()}, and
with {opt controls()} and {opt absorb()}.

{dlgtab:Reporting}

{phang}
{opt reportreg} displays the regressions behind the fit lines, as
{cmd:binscatter} does, on the same sample as the bins, and stores their
coefficients in {cmd:e(reportreg_b)}. All weight types are allowed; with
{opt iweight}s, absorbed effects enter these regressions as indicators.

{phang}
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.


{title:Control and weight handling}

{pstd}
Weight expressions are evaluated once. Missing weights are excluded, selected
weights must be positive, and fweights must be integers. Factor-variable
controls and interactions are expanded before estimation. Collinear control
variables are omitted. If the adjustment cannot be computed, the command
returns an error. Fixed effects are removed by iterative projection. {p_end}


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

{pstd}Medians, separate fits on each side of 25 mpg, and the bins in a variable:{p_end}
{phang2}{cmd:. cbinscatter price mpg, medians rd(25) genxq(mpg_bin)}{p_end}

{pstd}Save bin data without displaying graph:{p_end}
{phang2}{cmd:. cbinscatter price mpg, nograph savedata(mybins)}{p_end}

{pstd}Using the binsreg method (Cattaneo et al.):{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(binsreg)}{p_end}

{pstd}Compare classic vs binsreg methods:{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(classic) savedata(classic_bins)}{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight) method(binsreg) savedata(binsreg_bins)}{p_end}


{marker legacy}{...}
{dlgtab:Legacy options}

{pstd}
{cmd:binscatter} 7.02's legacy options work as there, each with
{cmd:binscatter}'s note: {opt nbins(#)} is {opt nquantiles(#)}, {opt create_xq}
is {cmd:genxq(q_}{it:xvar}{cmd:)}, {opt x_q(varname)} is {opt xq(varname)},
{opt symbols()} is {opt msymbols()}, {cmd:linetype(noline)} is
{cmd:linetype(none)}, and {opt unique()} is ignored. With both an option and its
legacy form, {cmd:binscatter} shows its message and stops with r(0);
{cmd:cbinscatter} shows the same message and exits with r(198).
({cmd:binscatter}'s legacy {opt method()}, which it ignores, is
{cmd:cbinscatter}'s {opt method()}.)

{marker binsreg}{...}
{dlgtab:Differences from binsreg}

{pstd}
{cmd:method(binsreg)} computes {cmd:binsreg}'s binned means on fixed quantile
bins of the raw x. {cmd:binsreg} was not available to verify against, so its
other features are not provided: the IMSE-optimal number of bins and
{opt binsmethod()}, {opt nbinsrot()}, {opt binspos()}, {opt pselect()} and
{opt sselect()}; {opt dots(p s)} with p > 0, {opt line()}, {opt ci()}, {opt cb()},
{opt deriv()}, {opt at()}, {opt vce()}, {opt level()} and {opt nsims()} (no
confidence intervals, bands or tests); {opt polyreg(p)} above 3 (the fit is at
most cubic, {cmd:linetype(cubic)}); {opt samebinsby}, {opt bycolors()},
{opt bysymbols()} and {opt bylpatterns()} (each by() group is binned on its own,
and the series take {opt colors()}, {opt mcolors()}, {opt lcolors()} and
{opt msymbols()}); {opt masspoints()}, {opt dfcheck()}, {opt noplot} (use
{opt nograph}), {opt plotxrange()} and {opt plotyrange()} (use {cmd:xscale()} or
{cmd:yscale()}), {opt reghdfeopt()} and {cmd:binsreg}'s {opt savedata()} layout.
{cmd:binsreg}'s positional controls ({cmd:binsreg y x w1 w2}) are given in
{opt controls(w1 w2)}: with more than two variables, {cmd:cbinscatter} reads them
as several y variables and one x, as {cmd:binscatter} does. None of
{cmd:binsreg}'s {cmd:e()} results are stored; see below for those of
{cmd:cbinscatter}.

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
{synopt:{cmd:e(discrete)}}1 if the bins are the values of x ({opt discrete}, or fewer values than {opt nquantiles()}){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cbinscatter}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable (variables){p_end}
{synopt:{cmd:e(xvar)}}x variable{p_end}
{synopt:{cmd:e(controls)}}control variables{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(by)}}by variable{p_end}
{synopt:{cmd:e(linetype)}}fit line type{p_end}
{synopt:{cmd:e(rd)}}with {opt rd()}, the break points{p_end}
{synopt:{cmd:e(graphcmd)}}the {cmd:twoway} command that drew the graph{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bindata)}}matrix of bin statistics: by_group, bin_id, x_mean, y_mean, n_obs{p_end}
{synopt:{cmd:e(coefs)}}fit line coefficients, one row per by-group (and {opt rd()} interval): constant, x, x^2, x^3 (for a polynomial fit){p_end}
{synopt:{cmd:e(fit_stats)}}R-squared and number of observations of each fit{p_end}
{synopt:{cmd:e(reportreg_b)}}with {opt reportreg}, the coefficients of each group's fit regression, one row per by-group (and {opt rd()} interval){p_end}
{synopt:{cmd:e(rdintervals)}}with {opt rd()}, each interval's bounds ({cmd:gt}, {cmd:lt_eq}), as in {cmd:binscatter}{p_end}
{synopt:{cmd:e(y}{it:#}{cmd:_coefs)}}with a fit line, {cmd:binscatter}'s coefficients of the {it:#}th y's fits: the x terms from the highest power down and {cmd:_cons}, one row per by() group and {opt rd()} interval, named {cmd:by}{it:#}{cmd::rd}{it:#}{p_end}
{synopt:{cmd:e(byvalues)}}with a numeric {opt by()} variable, its values, one row per group ({cmd:by1}, {cmd:by2}, ...), as in {cmd:binscatter}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the estimation sample, as in {cmd:binscatter}{p_end}


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
