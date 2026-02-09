{smcl}
{* *! version 1.0.1 07Feb2026}{...}
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
{synopt:{opt genxq(varname)}}generate bin assignment variable{p_end}
{synopt:{opt save:data(filename)}}save bin data to file{p_end}

{syntab:Graph Options}
{synopt:{opt nog:raph}}suppress graph output{p_end}
{synopt:{opt title(string)}}graph title{p_end}
{synopt:{opt yt:itle(string)}}y-axis title (default: yvar name){p_end}
{synopt:{opt xt:itle(string)}}x-axis title (default: xvar name){p_end}
{synopt:{opt legend(string)}}legend options{p_end}
{synopt:{opt colors(string)}}colors for series{p_end}
{synopt:{opt msymbols(string)}}marker symbols for series{p_end}
{synopt:{opt mlabels(string)}}marker labels for series{p_end}
{synopt:{it:twoway_options}}additional twoway graph options{p_end}

{syntab:Reporting}
{synopt:{opt reportreg}}report underlying regression{p_end}
{synopt:{opt v:erbose}}display detailed progress information and timing breakdown{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{opt aweight}s, {opt fweight}s, {opt pweight}s, and {opt iweight}s are allowed;
see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:cbinscatter} is a high-performance replacement for {cmd:binscatter} (v7.02) by Michael Stepner that uses
a C plugin for all data computations. It creates binned scatter plots by dividing
the x variable into equal-sized quantile bins and plotting the mean of y within
each bin.

{pstd}
All heavy computation (data loading, residualization, bin computation, and line
fitting) is performed in optimized C code with OpenMP parallelization. This makes
{cmd:cbinscatter} substantially faster than pure Stata implementations, especially
on large datasets.

{pstd}
Key features:

{phang2}{bf:Residualization:} Control variables are partialled out using OLS. Fixed
effects specified with {opt absorb()} are absorbed using the same HDFE algorithm
as {cmd:creghdfe}.{p_end}

{phang2}{bf:Histogram-based binning:} Uses an O(N) histogram algorithm for bin
assignment instead of sorting, providing dramatic speedups on large data.{p_end}

{phang2}{bf:Single-pass line fitting:} Linear and quadratic fits use closed-form
solutions computed in a single pass through the data.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt nquantiles(#)} specifies the number of equal-sized bins to create.
The default is 20. Must be between 2 and 1000.

{phang}
{opt controls(varlist)} specifies control variables to partial out from both
x and y before computing bins. The residualization is performed via OLS.

{phang}
{opt absorb(varlist)} specifies categorical variables representing fixed effects
to absorb from both x and y. Uses the same iterative demeaning algorithm as
{cmd:creghdfe}. Can be combined with {opt controls()}.

{phang}
{opt by(varname)} creates separate scatter series for each level of the
specified variable. Each group gets its own bins and optional fit line.

{phang}
{opt method(method)} specifies the residualization method when controls or
fixed effects are specified. This follows the methodology discussion in
Cattaneo, Crump, Farrell, and Feng (2024), "On Binscatter".

{phang2}{opt classic} (the default) residualizes both x and y on the controls/FE,
then creates bins based on the residualized x values. This is the traditional
binscatter approach.{p_end}

{phang2}{opt binsreg} creates bins based on the raw (non-residualized) x values,
then computes the conditional mean of y within each bin after partialling out
controls/FE. This is the approach recommended by Cattaneo et al. and implemented
in their {cmd:binsreg} package. It avoids potential distortions that can occur
when the conditional expectation function is nonlinear.{p_end}

{pstd}
When no controls or absorb variables are specified, both methods produce
identical results.

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
Fit lines are computed from the underlying microdata, not the bin means.

{dlgtab:Data Options}

{phang}
{opt discrete} treats x as a discrete variable, creating one bin for each
unique value of x instead of quantile-based bins.

{phang}
{opt genxq(varname)} generates a new variable containing the bin assignment
for each observation. (Not yet implemented.)

{phang}
{opt savedata(filename)} saves the bin data (bin means, counts, etc.) to
the specified Stata data file.

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
Default is filled circles (O). (Not yet implemented.)

{phang}
{opt mlabels(string)} specifies marker labels for the scatter series.
(Not yet implemented.)

{phang}
{it:twoway_options} any other options are passed through to the underlying
{cmd:twoway} graph command.

{dlgtab:Reporting}

{phang}
{opt verbose} displays detailed progress information and timing breakdown.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:cbinscatter} uses all available CPU cores.


{marker examples}{...}
{title:Examples}

{pstd}Basic binned scatter plot:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. cbinscatter price mpg}{p_end}

{pstd}With 10 bins and no fit line:{p_end}
{phang2}{cmd:. cbinscatter price mpg, nquantiles(10) linetype(none)}{p_end}

{pstd}Control for other variables:{p_end}
{phang2}{cmd:. cbinscatter price mpg, controls(weight length)}{p_end}

{pstd}With fixed effects:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. cbinscatter ln_wage tenure, absorb(idcode year)}{p_end}

{pstd}Separate series by group:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
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
{synopt:{cmd:e(nquantiles)}}number of bins{p_end}
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
{synopt:{cmd:e(coefs)}}fit line coefficients (if linetype specified){p_end}
{synopt:{cmd:e(fit_stats)}}R-squared and other fit statistics{p_end}


{title:Performance}

{pstd}
{cmd:cbinscatter} achieves significant speedups through:

{phang2}1. {bf:Parallel data loading:} Data is loaded from Stata in parallel with 8-way loop unrolling.{p_end}

{phang2}2. {bf:Histogram-based binning:} Instead of sorting (O(N log N)), uses a histogram approach with O(N) complexity and sequential memory access, providing 40x+ speedups for bin computation.{p_end}

{phang2}3. {bf:Single-pass statistics:} Bin means and fit coefficients are computed in single passes through the data using closed-form solutions.{p_end}

{phang2}4. {bf:Minimal Stata overhead:} Graph generation uses optimized approaches to minimize overhead on large datasets.{p_end}

{pstd}
On a dataset with 25 million observations, {cmd:cbinscatter} completes in under 1 second including graph generation.


{title:References}

{phang}
Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024. On Binscatter.
{it:American Economic Review} 114(5): 1488-1514.
{browse "https://doi.org/10.1257/aer.20221576"}

{pstd}
The {opt method(binsreg)} option implements the approach recommended in this paper,
which avoids potential distortions in the classic binscatter method when the
conditional expectation function is nonlinear.


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
This command is inspired by and provides similar functionality to {cmd:binscatter}
by Michael Stepner. The original {cmd:binscatter} package is available at
{browse "https://github.com/michaelstepner/binscatter"}. We thank Michael Stepner
for his pioneering work on binned scatter plot visualization in Stata.


{title:Also see}

{psee}
Online: {browse "https://github.com/michaelstepner/binscatter":binscatter} (if installed),
{help ctools}, {help creghdfe}, {help twoway}
{p_end}
