{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "crangejoin##syntax"}{...}
{viewerjumpto "Description" "crangejoin##description"}{...}
{viewerjumpto "Options" "crangejoin##options"}{...}
{viewerjumpto "Examples" "crangejoin##examples"}{...}

{title:Title}

{pstd}{bf:crangejoin} {hline 2} C-accelerated inclusive interval joins

{marker syntax}{...}

{title:Syntax}

{p 8 12 2}{cmd:crangejoin} {it:keyvar low high} {cmd:using} {it:filename}
[{cmd:,} {opt by(varlist)} {opt k:eepusing(varlist)} {opt p:refix(string)}
{opt s:uffix(string)} {opt a:ll} {opt threads(#)} {opt verbose}
{opt snap:shot} {opt fast}]

{marker description}{...}

{title:Description}

{pstd}{cmd:crangejoin} joins observations by numeric intervals.
Each master row is paired with every using row whose numeric key is within its
inclusive bounds and whose by() values agree. Unmatched master rows remain with
missing using values. Master order is preserved; matches within a master row are
ordered by using key, with ties preserving their original using order.
The joined dataset replaces the data in memory. There is no if/in syntax;
filter the master dataset beforehand. The command does not store results in {cmd:r()} or {cmd:e()}.

{pstd}A bound can be a numeric master variable (absolute bounds) or a number.
Numeric bounds are offsets from the master key when it exists, and absolute
bounds otherwise. A missing lower/upper bound is unbounded in that direction.
Using rows with missing keys never match. If the master contains keyvar, its
missing keys never match, even with explicit bound variables. Reversed intervals
produce no match. As in rangejoin, no valid master intervals or no nonmissing
using keys returns error 2000.

{marker options}{...}

{title:Options}

{pstd}{opt by()} restricts matches to equal numeric/string groups; missing group
values participate and extended numeric missings remain distinct. keyvar cannot
also be a by() variable.

{pstd}{opt keepusing()} selects using variables. Key and by variables are always
loaded.

{pstd}{opt prefix()} and {opt suffix()} rename using variables sharing names
with master variables. The default is suffix(_U) when neither is specified.
{opt all} renames every using
variable except the by() variables. Name conflicts return an error.

{pstd}{opt threads(#)} sets the maximum number of threads, from 1 to 256.
Omitting it or specifying zero uses the {help ctools##threads:ctools setting}. {opt verbose}
reports loading, index sorting, group/key indexing, matching/counting, and
output writing times.

{pstd}By default, an error or interruption restores the master dataset.
Recovery uses an in-memory copy on Stata 16 or later and {help preserve} on
earlier versions.

{pstd}{opt snapshot} uses {help preserve} for recovery on every Stata version.
Stata may write the copy to disk when it exceeds {cmd:c(max_preservemem)}.

{pstd}{opt fast} skips recovery. An error or interruption can leave a partially
joined dataset. It cannot be combined with {opt snapshot}.


{marker examples}{...}

{title:Examples}

{phang2}{cmd:. crangejoin date -5 5 using events.dta, by(firm)}{p_end}

{phang2}{cmd:. crangejoin price lower upper using houses.dta, keepusing(address)}{p_end}

{phang2}{cmd:. crangejoin date . . using history.dta, by(person) suffix(_history)}{p_end}

{title:Remarks}

{pstd}Every matching pair produces an output observation, so the joined dataset
can be much larger than either input. Inputs and output are limited to
2,147,483,647 observations. Numeric and fixed-width string variables are
supported. Master variables and retained using variables cannot be {cmd:strL}.
{cmd:rangejoin} and {cmd:rangestat} are not required.


{title:Also see}

{pstd}{help cmerge}, {help crangestat}, {help ctools}
