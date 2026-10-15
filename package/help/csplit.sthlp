{smcl}
{title:Title}
{pstd}{bf:csplit} {hline 2} C-accelerated string splitting

{title:Syntax}
{p 8 12 2}{cmd:csplit} {it:strvar} {ifin}
[{cmd:,} {opt gen:erate(stub)} {opt parse(parse_strings)} {opt notrim}
{opt limit(#)} {opt destring} {opt force} {opt float} {opt ignore(strings)}
{opt percent} {opt threads(#)} {opt verbose}]

{title:Description}
{pstd}Splits a string into new variables named stub1, stub2, and so on, using the
byte-string and delimiter semantics of {help split}. The default stub is the input
variable name. Default parsing uses spaces and trims successive remainders.
Explicit parse() strings can contain multiple characters. For delimiters matching
at the same position, the last listed delimiter takes precedence. Adjacent
delimiters preserve empty fields. Trailing
separators do not create an extra variable unless another observation requires it.
Unselected observations have empty outputs (missing after numeric conversion).
The source variable and observation order are preserved. An empty selected
sample returns error 2000.

{title:Options}
{pstd}{opt generate(stub)} names the new variables stub1, stub2, and so on.
The default stub is the source variable name. Output names must be new.
{pstd}{opt parse(parse_strings)} specifies one or more literal delimiters, for
example {cmd:parse("|" ";" "::")}. Without parse(), spaces delimit fields and
successive remainders are trimmed. With explicit delimiters, outer spaces are
trimmed initially but spaces within the remaining text are preserved.
This is literal splitting, without CSV-style quoted-field handling.
{pstd}{opt notrim} retains outer spaces and requires explicit parse() strings.
{opt limit(#)} creates at most # fields; the remaining unsplit text is discarded.
{pstd}{opt destring} applies native {help destring} to the staged fields, with its
force, float, ignore(), and percent options. Those options require destring.
{opt force} turns unconvertible values into missing; {opt float} requests float
storage; {opt ignore(strings)} removes specified characters; {opt percent}
converts percentage values. Without force, a field containing nonnumeric values
can remain a string, as with native destring.
{pstd}{opt threads(#)} accepts a nonnegative integer and limits C workers; zero uses the default.
{opt verbose} reports load, tokenize/size, and write times.

{title:Stored results}
{pstd}The macros {cmd:r(varlist)} and {cmd:r(nvars)} contain the generated variable
names and their count. Scalar {cmd:r(k_new)} also contains the count, as with split.

{title:Examples}
{phang2}{cmd:. csplit address, generate(part) parse(",")}
{phang2}{cmd:. csplit codes, parse("::" ":") limit(4)}
{phang2}{cmd:. csplit values, parse(",") destring}

{title:Implementation and limits}
{pstd}Fixed-width strings up to 2045 bytes use the C tokenizer. The first pass determines the
number and widths of fields. A single delimiter is terminated in place; multiple
delimiters use packed tokens in the loaded buffer. The write pass reuses these
tokens without parsing again. Both passes parallelize
across observations. Naming errors
and conversion failures leave the source data unchanged. String bytes are
preserved; no Unicode normalization is performed. strL input uses native split
to retain long-string behavior, because the plugin cannot write strL output.

{title:Also see}
{pstd}{help split}, {help destring}, {help ctools}
