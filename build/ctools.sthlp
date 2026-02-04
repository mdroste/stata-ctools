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
[{cmd:,} {opt env:ironment_check} {opt update}]


{marker options}{...}
{title:Options}

{phang}
{opt environment_check} loads the ctools C plugin and reports whether it was
successfully loaded. This is useful for diagnosing installation issues. The
plugin is automatically detected for your platform (macOS ARM/Intel, Windows,
Linux).

{phang}
{opt update} downloads and installs the latest version of ctools from GitHub.
You may need to restart Stata after updating for all changes to take effect.


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


{marker commands}{...}
{title:Available Commands}

{pstd}
{ul:Data Management}

{p2colset 5 20 22 2}{...}
{p2col:{help cimport}}Import text-delimited and Excel data (replaces {help import delimited}){p_end}
{p2col:{help cexport}}Export text-delimited and Excel data (replaces {help export delimited}/{help export excel}){p_end}
{p2col:{help csort}}Sort dataset (replaces {help sort}){p_end}
{p2col:{help cmerge}}Merge (join) datasets (replaces {help merge}){p_end}
{p2col:{help csample}}Resampling without replacement (replaces {help sample}){p_end}
{p2col:{help cbsample}}Resampling with replacement (replaces {help bsample}){p_end}
{p2col:{help cencode}}Recast string as labeled numeric (replaces {help encode}){p_end}
{p2col:{help cdecode}}Recast labeled numeric as string (replaces {help decode}){p_end}
{p2col:{help cdestring}}Recast string as numeric type (replaces {help destring}){p_end}
{p2col:{help cwinsor}}Winsorize variables (replaces {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2}/{help gstats winsor:gstats winsor}){p_end}
{p2col:{help crangestat}}Range statistics of variables (replaces {browse "https://ideas.repec.org/c/boc/bocode/s458118.html":rangestat}){p_end}

{pstd}
{ul:Estimation}

{p2col:{help creghdfe}}OLS with multi-way fixed effects (replaces {browse "https://github.com/sergiocorreia/reghdfe":reghdfe}){p_end}
{p2col:{help civreghdfe}}2SLS/GMM with multi-way fixed effects (replaces {browse "https://github.com/sergiocorreia/ivreghdfe":ivreghdfe}){p_end}
{p2col:{help cqreg}}Quantile regression (replaces {help qreg}){p_end}
{p2col:{help cpsmatch}}Propensity score matching (replaces {browse "https://ideas.repec.org/c/boc/bocode/s457730.html":psmatch2}){p_end}

{pstd}
{ul:Visualization}

{p2col:{help cbinscatter}}Binned scatter plots (replaces {browse "https://github.com/michaelstepner/binscatter":binscatter}){p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Examples}

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
Thanks to:

{p 8 12 2}- Sergio Correia for {browse "https://github.com/sergiocorreia/ftools":ftools},
{browse "https://github.com/sergiocorreia/reghdfe":reghdfe}, and
{browse "https://github.com/sergiocorreia/ivreghdfe":ivreghdfe} (with Lars Vilhuber){p_end}
{p 8 12 2}- Mauricio Caceres Bravo for {browse "https://github.com/mcaceresb/gtools":gtools}{p_end}
{p 8 12 2}- Christopher (Kit) Baum, Mark E Schaffer, and Steven Stillman for
{browse "https://ideas.repec.org/c/boc/bocode/s425401.html":ivreg2}{p_end}
{p 8 12 2}- Robert Picard, Nicholas J. Cox, and Roberto Ferrer for
{browse "https://ideas.repec.org/c/boc/bocode/s458161.html":rangestat}{p_end}
{p 8 12 2}- Michael Stepner for {browse "https://github.com/michaelstepner/binscatter":binscatter}{p_end}
{p 8 12 2}- Yujun Lian for {browse "https://ideas.repec.org/c/boc/bocode/s457765.html":winsor2}{p_end}
{p 8 12 2}- Sascha Witt for the {browse "https://github.com/SaschaWitt/ips4o":IPS4o} sorting algorithm{p_end}


{title:Also see}

{psee}
Online: {help cimport}, {help cexport}, {help csort}, {help cmerge}, {help cencode},
{help cdecode}, {help cdestring}, {help csample}, {help cbsample}, {help cwinsor},
{help crangestat}, {help creghdfe}, {help civreghdfe}, {help cqreg}, {help cpsmatch},
{help cbinscatter}
{p_end}
