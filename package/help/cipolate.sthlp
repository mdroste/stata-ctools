{smcl}
{title:Title}
{pstd}{bf:cipolate} {hline 2} C-accelerated linear interpolation (ctools)

{title:Syntax}
{p 8 12 2}{cmd:cipolate} {it:yvar xvar} {ifin}, {opt gen:erate(newvar)}
[{opt by(varlist)} {opt ep:olate} {opt threads(#)} {opt verbose}]

{title:Description}
{pstd}Computes the linear interpolation performed by {help ipolate}. Within each
by-group, nonmissing y values at the same x are averaged. Missing y values between
known points are interpolated. Missing x values and observations outside if/in
produce missing results. The output is double precision and the input observation
order is preserved. Without epolate, results outside the known range are missing.
This is a linear interpolator, not the separately distributed
SSC cubic-interpolation command also named cipolate. Check {cmd:which cipolate}
if both packages are installed.

{title:Options}
{pstd}{opt generate()} names the required new output variable.
{pstd}{opt by()} computes separate interpolations for numeric or string groups.
The {cmd:by:} prefix is also accepted; do not combine it with by(). Extended missing
group values remain distinct. Group string data must fit within 2045 bytes.
{pstd}{opt epolate} extrapolates at both ends using the nearest two distinct
known x points. A group with only one known point cannot extrapolate.
{pstd}{opt threads(#)} accepts a nonnegative integer and limits C worker threads; zero selects the default.
{opt verbose} reports loading, sorting, interpolation, and storage times.

{title:Stored results}
{pstd}The generated variable is the result; no command-specific stored results
are promised.

{title:Examples}
{phang2}{cmd:. cipolate income year, generate(income_i) by(person)}
{phang2}{cmd:. bysort person: cipolate income year, generate(income_e) epolate}

{title:Implementation}
{pstd}Loads only the selected observations and relevant columns, sorts an index
once, and scans each group. Numeric keys use stable parallel radix passes;
mixed and string keys use parallel merges. Already ordered inputs skip sorting.
Output is staged in a temporary variable until the
calculation succeeds. Performance depends on group sizes, data transfer, and
hardware; no fixed speedup is promised.

{title:Also see}
{pstd}{help ipolate}, {help ctools}
