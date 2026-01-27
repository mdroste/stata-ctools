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
{cmd:ctools} is a suite of high-performance C-accelerated Stata commands.
These commands use optimized C plugins with OpenMP parallelization to provide
significant speed improvements over native Stata commands and popular user-written
packages for common data operations and statistical estimation.

{pstd}
All commands are designed as drop-in replacements with syntax that closely
matches the original commands they replace.


{marker options}{...}
{title:Options}

{phang}
{opt help} displays this help file.

{phang}
{opt update} checks for available updates. (Not yet implemented.)


{marker commands}{...}
{title:Available Commands}

{pstd}
{ul:Data Management}

{p2colset 5 20 22 2}{...}
{p2col:{help csort}}High-performance parallel sorting (replaces {help sort}){p_end}
{p2col:{help cmerge}}C-accelerated merge (replaces {help merge}){p_end}
{p2col:{help cimport}}Multi-threaded CSV import (replaces {help import delimited}){p_end}
{p2col:{help cexport}}Parallel data export to CSV/Excel (replaces {help export delimited}/{help export excel}){p_end}
{p2col:{help cencode}}Parallel string encoding (replaces {help encode}){p_end}
{p2col:{help cdecode}}Parallel numeric decoding (replaces {help decode}){p_end}
{p2col:{help cdestring}}Parallel string-to-numeric conversion (replaces {help destring}){p_end}

{pstd}
{ul:Sampling}

{p2col:{help csample}}Random sampling without replacement (replaces {help sample}){p_end}
{p2col:{help cbsample}}Bootstrap sampling with replacement (replaces {help bsample}){p_end}

{pstd}
{ul:Data Transformation}

{p2col:{help cwinsor}}Parallel winsorization (replaces {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2}){p_end}

{pstd}
{ul:Statistical Estimation}

{p2col:{help creghdfe}}HDFE linear regression (replaces {browse "https://github.com/sergiocorreia/reghdfe":reghdfe}){p_end}
{p2col:{help civreghdfe}}IV/2SLS with HDFE (replaces {browse "https://github.com/sergiocorreia/ivreghdfe":ivreghdfe}){p_end}
{p2col:{help cqreg}}Quantile regression with HDFE (replaces {help qreg}){p_end}

{pstd}
{ul:Visualization}

{p2col:{help cbinscatter}}Binned scatter plots (replaces {browse "https://github.com/michaelstepner/binscatter":binscatter}){p_end}
{p2colreset}{...}


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
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
Several commands in this package are inspired by excellent user-written packages.
We gratefully acknowledge:

{p 8 12 2}- Sergio Correia for {browse "https://github.com/sergiocorreia/reghdfe":reghdfe} and
{browse "https://github.com/sergiocorreia/ivreghdfe":ivreghdfe}{p_end}
{p 8 12 2}- Michael Stepner for {browse "https://github.com/michaelstepner/binscatter":binscatter}{p_end}
{p 8 12 2}- Yujun Lian for {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2}{p_end}


{title:Also see}

{psee}
Online: {help csort}, {help cmerge}, {help cimport}, {help cexport}, {help creghdfe},
{help civreghdfe}, {help cqreg}, {help cbinscatter}, {help cencode}, {help cdecode},
{help cdestring}, {help cwinsor}, {help csample}, {help cbsample}
{p_end}
