{smcl}
{title:Title}
{pstd}{bf:crangejoin} {hline 2} C-accelerated inclusive interval joins

{title:Syntax}
{p 8 12 2}{cmd:crangejoin} {it:keyvar low high} {cmd:using} {it:filename}
[{cmd:,} {opt by(varlist)} {opt k:eepusing(varlist)} {opt p:refix(string)}
{opt s:uffix(string)} {opt a:ll} {opt threads(#)} {opt verbose}]

{title:Description}
{pstd}Implements the matching interface of Robert Picard's {cmd:rangejoin} 1.1.3.
Each master row is paired with every using row whose numeric key is within its
inclusive bounds and whose by() values agree. Unmatched master rows remain with
missing using values. Master order is preserved; matches within a master row are
ordered by using key, with ties preserving their original using order.
The joined dataset replaces the data in memory. There is no if/in syntax;
filter the master dataset beforehand. No command-specific stored results are
promised.

{pstd}A bound can be a numeric master variable (absolute bounds) or a number.
Numeric bounds are offsets from the master key when it exists, and absolute
bounds otherwise. A missing lower/upper bound is unbounded in that direction.
Using rows with missing keys never match. If the master contains keyvar, its
missing keys never match, even with explicit bound variables. Reversed intervals
produce no match. As in rangejoin, no valid master intervals or no nonmissing
using keys returns error 2000.

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
{pstd}{opt threads(#)} accepts a nonnegative integer and limits C workers; zero uses the default. {opt verbose}
reports loading, index sorting, group/key indexing, matching/counting, and
output writing times.

{title:Examples}
{phang2}{cmd:. crangejoin date -5 5 using events.dta, by(firm)}
{phang2}{cmd:. crangejoin price lower upper using houses.dta, keepusing(address)}
{phang2}{cmd:. crangejoin date . . using history.dta, by(person) suffix(_history)}

{title:Implementation and limits}
{pstd}The C kernel sorts a using index, finds intervals by binary search, checks
output size, and writes matched columns directly. A group index and contiguous
key array avoid repeated multi-column comparisons during interval searches.
Output writes are divided into row tiles so even narrow or highly expanded
joins can use multiple workers. The wrapper preserves variable
metadata and restores the master dataset on failure. Runtime and memory depend
on output size; all combinations may be much larger than either input. The
plugin limits inputs and output to 2^31-1 observations. Numeric and fixed-width
strings are supported. Master variables and selected using variables cannot be
strL; recast them to fixed strings when lossless.
{cmd:rangejoin} and {cmd:rangestat} are not required to run crangejoin.

{title:Also see}
{pstd}{help cmerge}, {help crangestat}, {help ctools}
