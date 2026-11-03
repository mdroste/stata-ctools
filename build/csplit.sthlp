{smcl}
{* *! version 1.2.4 20261103}{...}
{viewerjumpto "Syntax" "csplit##syntax"}{...}
{viewerjumpto "Description" "csplit##description"}{...}
{viewerjumpto "Options" "csplit##options"}{...}
{viewerjumpto "Examples" "csplit##examples"}{...}
{viewerjumpto "Stored results" "csplit##results"}{...}

{title:Title}

{pstd}{bf:csplit} {hline 2} C-accelerated string splitting

{marker syntax}{...}

{title:Syntax}

{p 8 12 2}{cmd:csplit} {it:strvar} {ifin}
[{cmd:,} {opt gen:erate(stub)} {opt parse(parse_strings)} {opt notrim}

{pstd}{opt limit(#)} {opt destring} {opt force} {opt float} {opt ignore(strings)}
{opt percent} {opt threads(#)} {opt verbose}]

{marker description}{...}

{title:Description}

{pstd}{cmd:csplit} splits a string variable into new variables named stub1,
stub2, and so on. It follows {help split}'s rules for matching delimiter
bytes. The default stub is the input variable name. By default, spaces
separate fields, and leading and trailing spaces are removed from the
remaining text after each split. Explicit parse() strings can contain multiple
characters. For delimiters matching at the same position, the last listed
delimiter takes precedence. Adjacent delimiters preserve empty fields.
Trailing separators do not create an extra variable unless another observation
requires it. Unselected observations have empty outputs (missing after numeric
conversion). The source variable and observation order are preserved. An empty
selected sample returns error 2000.

{marker options}{...}

{title:Options}

{pstd}{opt generate(stub)} names the new variables stub1, stub2, and so on.
The default stub is the source variable name. Output names must be new.

{pstd}{opt parse(parse_strings)} specifies one or more literal delimiters, for
example {cmd:parse("|" ";" "::")}. Without {opt parse()}, spaces separate
fields, and leading and trailing spaces are removed from the remaining text
after each split. With explicit delimiters, outer spaces are trimmed initially
but spaces within the remaining text are preserved. This is literal splitting,
without CSV-style quoted-field handling.

{pstd}{opt notrim} retains outer spaces and requires explicit parse() strings.

{pstd}{opt limit(#)} creates at most # fields; the remaining unsplit text is discarded.

{pstd}{opt destring} uses Stata's {help destring} to convert the split fields
to numeric variables. It accepts the options {opt force}, {opt float},
{opt ignore()}, and {opt percent}. Those options require destring. {opt force}
turns unconvertible values into missing; {opt float} requests float storage;
{opt ignore(strings)} removes specified characters; {opt percent} converts
percentage values. Without force, a field containing nonnumeric values can
remain a string, as with Stata's {cmd:destring}.

{pstd}{opt threads(#)} sets the maximum number of threads, from 1 to 256.
Omitting it or specifying zero uses the {help ctools##threads:ctools setting}.

{pstd}{opt verbose} reports the time spent loading strings, splitting them and
determining output lengths, and writing the new variables.

{marker results}{...}

{title:Stored results}

{pstd}For fixed-width input, the macros {cmd:r(varlist)} and {cmd:r(nvars)} contain the generated variable
names and their count. Scalar {cmd:r(k_new)} also contains the count. For {cmd:strL} input, results follow {help split}.

{marker examples}{...}

{title:Examples}

{phang2}{cmd:. csplit address, generate(part) parse(",")}{p_end}

{phang2}{cmd:. csplit codes, parse("::" ":") limit(4)}{p_end}

{phang2}{cmd:. csplit values, parse(",") destring}{p_end}

{title:Remarks}

{pstd}The ctools plugin splits fixed-width string variables. For {cmd:strL}
variables, {cmd:csplit} uses Stata's {help split}. Delimiters match literal
bytes; quoted CSV fields are not interpreted. Use {help cimport} to read a
delimited file with quoted fields.


{title:Also see}

{pstd}{help split}, {help destring}, {help ctools}
