{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "cinterpolate##syntax"}{...}
{viewerjumpto "Description" "cinterpolate##description"}{...}
{viewerjumpto "Options" "cinterpolate##options"}{...}
{viewerjumpto "Examples" "cinterpolate##examples"}{...}
{viewerjumpto "Stored results" "cinterpolate##results"}{...}

{title:Title}

{pstd}{bf:cinterpolate} {hline 2} C-accelerated linear and cubic interpolation (ctools)

{marker syntax}{...}

{title:Syntax}

{p 8 12 2}{cmd:cinterpolate} {it:yvar xvar} {ifin}, {opt gen:erate(newvar)}
[{opt by(varlist)} {opt ep:olate} {opt cub:ic} {opt threads(#)} {opt verbose}]

{marker description}{...}

{title:Description}

{pstd}Defaults to the linear interpolation performed by {help ipolate}; {opt cubic}
selects local cubic interpolation as in Nicholas J. Cox's SSC {help cipolate}.
Neither reference command is required. Within each
by-group, nonmissing y values at the same x are averaged. Missing y values between
known points are interpolated. Missing x values and observations outside if/in
produce missing results. The output is double precision and the input observation
order is preserved. Without epolate, results outside the known range are missing.
For cubic interpolation, additional gaps near the ends remain missing unless
two distinct known x points are available on each side.

{marker options}{...}

{title:Options}

{pstd}{opt generate()} names the required new output variable.

{pstd}{opt by()} computes separate interpolations for numeric or string groups.
The {cmd:by:} prefix is also accepted; do not combine it with by(). Extended missing
group values remain distinct. Group string data must fit within 2045 bytes.

{pstd}{opt epolate} extrapolates at both ends using the nearest two distinct
known x points. A group with only one known point cannot extrapolate.

{pstd}{opt cubic} fits a cubic polynomial through the nearest two distinct known
x points on each side of a missing value. If either side has fewer than two
anchors, the result remains missing. This is local Lagrange interpolation, not
a spline, and can overshoot. Duplicate x values are averaged before choosing
distinct anchors; extended missing x and y values are treated as missing.
These rules avoid repeated anchors and extended-missing quirks in legacy SSC
{cmd:cipolate}. The options {opt cubic} and {opt epolate} cannot be combined.

{pstd}{opt threads(#)} sets the maximum number of threads, from 1 to 256.
Omitting it or specifying zero uses the {help ctools##threads:ctools setting}.

{pstd}{opt verbose} reports loading, sorting, interpolation, and storage times.

{marker results}{...}

{title:Stored results}

{pstd}The generated variable is the result; the command does not store results in {cmd:r()} or {cmd:e()}.

{marker examples}{...}

{title:Examples}

{phang2}{cmd:. cinterpolate income year, generate(income_i) by(person)}{p_end}

{phang2}{cmd:. bysort person: cinterpolate income year, generate(income_e) epolate}{p_end}

{phang2}{cmd:. cinterpolate income year, generate(income_c) by(person) cubic}{p_end}

{title:Remarks}

{pstd}The command changes neither the source variables nor observation order.
The output variable must be new. Interpolation uses the selected sample only;
{cmd:if}/{cmd:in} excludes observations from both the known points and the output.


{pstd}For scripts using the former ctools linear command, replace {cmd:cipolate}
with {cmd:cinterpolate}. No ctools alias named {cmd:cipolate} is provided, so the
SSC command can coexist. Remove stale ctools {cmd:cipolate.ado} and
{cmd:cipolate.sthlp} files when updating an existing installation; locate them
with {cmd:which cipolate} and inspect their headers before removal.

{title:Also see}

{pstd}{help ipolate}, {help cipolate}, {help ctools}
