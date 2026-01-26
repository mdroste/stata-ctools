{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "ctools##syntax"}{...}
{viewerjumpto "Description" "ctools##description"}{...}
{viewerjumpto "Options" "ctools##options"}{...}
{viewerjumpto "Commands" "ctools##commands"}{...}
{viewerjumpto "Examples" "ctools##examples"}{...}
{title:Title}

{phang}
{bf:ctools} {hline 2} High-performance C-based tools for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:ctools}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt help}}display help information{p_end}
{synopt:{opt update}}check for updates (placeholder){p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ctools} is a collection of high-performance C-based tools for Stata.
These tools leverage optimized C code to provide significant speed improvements
over native Stata commands for common data operations.


{marker options}{...}
{title:Options}

{phang}
{opt help} displays this help file.

{phang}
{opt update} checks for available updates. (Not yet implemented.)


{marker commands}{...}
{title:Available Commands}

{phang}
{help csort} {hline 2} High-performance radix sort. Significantly faster than
Stata's native {cmd:sort} command for large datasets.

{phang}
{help cmerge} {hline 2} C-accelerated merge. Drop-in replacement for Stata's
{cmd:merge} command with parallel data loading and optimized sorting.

{phang}
{help cimport} {hline 2} C-accelerated CSV import. High-performance replacement
for {cmd:import delimited} with multi-threaded parallel parsing.

{phang}
{help cexport} {hline 2} C-accelerated CSV export. High-performance replacement
for {cmd:export delimited} with parallel data loading.

{phang}
{help creghdfe} {hline 2} C-accelerated high-dimensional fixed effects regression.
High-performance replacement for {cmd:reghdfe} with optimized fixed effects absorption.


{marker examples}{...}
{title:Examples}

{phang}{cmd:. ctools}{p_end}
{phang}{cmd:. ctools, help}{p_end}
{phang}{cmd:. csort myvar}{p_end}
{phang}{cmd:. csort var1 var2 var3, verbose}{p_end}
{phang}{cmd:. cmerge m:1 id using lookup.dta}{p_end}
{phang}{cmd:. cimport delimited using data.csv, clear}{p_end}
{phang}{cmd:. cexport delimited using output.csv, replace}{p_end}
{phang}{cmd:. creghdfe y x1 x2, absorb(firm year)}{p_end}


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help csort}, {help cmerge}, {help cimport}, {help cexport}, {help creghdfe}
{p_end}
