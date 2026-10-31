{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "cdistinct##syntax"}{...}
{viewerjumpto "Description" "cdistinct##description"}{...}
{viewerjumpto "Options" "cdistinct##options"}{...}
{viewerjumpto "Examples" "cdistinct##examples"}{...}
{viewerjumpto "Stored results" "cdistinct##results"}{...}
{title:Title}

{pstd}{bf:cdistinct} {hline 2} Report distinct values of variables

{marker syntax}{...}
{title:Syntax}

{p 8 12 2}{cmd:cdistinct} [{varlist}] {ifin}
[{cmd:,} {opt miss:ing} {opt j:oint} {opt a:bbrev(#)}
{opt min:imum(#)} {opt max:imum(#)} {opt verbose} {opt threads(#)}]

{pstd}{cmd:by:} is allowed.

{marker description}{...}
{title:Description}

{pstd}{cmd:cdistinct} reports the observation count and distinct-value count for
each variable, following the user-written {cmd:distinct} command. Omitting the
variable list selects all variables. Numeric and string variables are allowed.
Missing values are excluded separately for each variable unless {opt missing}
is specified. The command leaves values and observation order unchanged.

{pstd}The C backend scans ordered keys directly and counts repeated values in
a bounded hash table. It switches to sorting when the table fills. Counts are
exact. Text in {cmd:strL} variables
is limited to 2,045 bytes; longer values produce an error rather than being truncated.

{marker options}{...}
{title:Options}

{pstd}{opt joint} counts combinations of the entire variable list. Without
{opt missing}, a row with any missing key is excluded.

{pstd}{opt missing} includes numeric missing values and empty strings.
Extended numeric missing codes count separately.

{pstd}{opt abbrev(#)} sets the width used for abbreviated variable names in the
display. The default is the longest variable name; the minimum width is five.

{pstd}{opt minimum(#)} and {opt maximum(#)} display only counts within the
specified range. The defaults are zero and no upper bound. Reversed limits
are exchanged. These options do not change the returned results.

{pstd}{opt threads(#)} sets the maximum worker count, from 1 to 256.
Zero or omission uses the {help ctools##threads:ctools setting}.

{pstd}{opt verbose} reports load and count times above the table.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. cdistinct firm year}{p_end}
{phang2}{cmd:. cdistinct firm year, joint missing}{p_end}
{phang2}{cmd:. bysort industry: cdistinct firm}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:r(N)} is the number of eligible observations and {cmd:r(ndistinct)}
is the number of distinct values. Without {opt joint}, both refer to the last
variable examined, even if its row is hidden by a display filter. An empty sample
returns zero.

{title:Remarks}

{pstd}{cmd:in} selects the same original rows for every variable. In {cmd:distinct}
1.2.1, internal sorting can change which rows {cmd:in} selects after the first
variable. Use a sample indicator with that reference command when comparing
several variables under a row restriction.

{title:Also see}
{pstd}{help ctools}, {help distinct}
