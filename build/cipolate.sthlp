{smcl}
{* *! version 1.2.0 20261029}{...}
{viewerjumpto "Syntax" "cipolate##syntax"}{...}
{viewerjumpto "Description" "cipolate##description"}{...}
{viewerjumpto "Options" "cipolate##options"}{...}
{viewerjumpto "Examples" "cipolate##examples"}{...}
{viewerjumpto "Stored results" "cipolate##results"}{...}

{title:Title}

{pstd}{bf:cipolate} {hline 2} C-accelerated linear interpolation (ctools)

{marker syntax}{...}

{title:Syntax}

{p 8 12 2}{cmd:cipolate} {it:yvar xvar} {ifin}, {opt gen:erate(newvar)}
[{opt by(varlist)} {opt ep:olate} {opt threads(#)} {opt verbose}]

{marker description}{...}

{title:Description}

{pstd}Computes the linear interpolation performed by {help ipolate}. Within each
by-group, nonmissing y values at the same x are averaged. Missing y values between
known points are interpolated. Missing x values and observations outside if/in
produce missing results. The output is double precision and the input observation
order is preserved. Without epolate, results outside the known range are missing.
This is a linear interpolator, not the separately distributed
SSC cubic-interpolation command also named cipolate. Check {cmd:which cipolate}
if both packages are installed.

{marker options}{...}

{title:Options}

{pstd}{opt generate()} names the required new output variable.

{pstd}{opt by()} computes separate interpolations for numeric or string groups.
The {cmd:by:} prefix is also accepted; do not combine it with by(). Extended missing
group values remain distinct. Group string data must fit within 2045 bytes.

{pstd}{opt epolate} extrapolates at both ends using the nearest two distinct
known x points. A group with only one known point cannot extrapolate.

{pstd}{opt threads(#)} sets the maximum number of threads, from 1 to 256.
Omitting it or specifying zero uses the {help ctools##threads:ctools setting}.

{pstd}{opt verbose} reports loading, sorting, interpolation, and storage times.

{marker results}{...}

{title:Stored results}

{pstd}The generated variable is the result; the command does not store results in {cmd:r()} or {cmd:e()}.

{marker examples}{...}

{title:Examples}

{phang2}{cmd:. cipolate income year, generate(income_i) by(person)}{p_end}

{phang2}{cmd:. bysort person: cipolate income year, generate(income_e) epolate}{p_end}

{title:Remarks}

{pstd}The command changes neither the source variables nor observation order.
The output variable must be new. Interpolation uses the selected sample only;
{cmd:if}/{cmd:in} excludes observations from both the known points and the output.


{title:Also see}

{pstd}{help ipolate}, {help ctools}
