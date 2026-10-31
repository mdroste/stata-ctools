{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "ctools##syntax"}{...}
{viewerjumpto "Description" "ctools##description"}{...}
{viewerjumpto "Options" "ctools##options"}{...}
{viewerjumpto "Threads" "ctools##threads"}{...}
{viewerjumpto "Commands" "ctools##commands"}{...}
{viewerjumpto "Examples" "ctools##examples"}{...}
{title:Title}

{phang}
{bf:ctools} {hline 2} C-accelerated tools for Stata


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:ctools}
[{cmd:,} {opt version} {opt env:ironment_check} {opt update}]

{p 8 17 2}
{cmdab:ctools}{cmd:,}
{opt thr:eads(#|default)}
[{opt perm:anently}]


{marker description}{...}
{title:Description}

{pstd}
{cmd:ctools} provides commands for data management, estimation, and binned
scatter plots. Computation runs in C, with parallel processing where available.
The commands and their help files are listed below.

{pstd}
Run {cmd:ctools} without options to list the commands. Use
{cmd:ctools, version} to report the installed ado and plugin versions.
The {opt verbose} option is accepted but does not change the output.


{marker options}{...}
{title:Options}

{phang}
{opt version} displays the ado and plugin versions and the plugin revision.
It stores them in {cmd:r(version)}, {cmd:r(plugin_version)}, and
{cmd:r(plugin_revision)}.

{phang}
{opt environment_check} loads the ctools C plugin and reports whether it was
successfully loaded. This is useful for diagnosing installation issues. The
plugin is automatically detected for your platform (macOS ARM/Intel, Windows,
Linux). It also reports the thread settings described {help ctools##threads:below}.

{phang}
{opt threads(#)} sets the number of threads ctools commands use for the rest of
the Stata session when they are called without their own {opt threads()}
option. {it:#} is an integer from 1 through 256. {opt threads(default)} removes
the session setting. The setting is stored in the global macro
{cmd:CTOOLS_THREADS}.

{phang}
{opt permanently}, specified with {opt threads()}, also saves the setting for
future Stata sessions, or with {opt threads(default)} removes the saved
setting. The setting is saved in a small per-user file outside the ado
directories, so it survives updates and reinstalls of ctools; its location is
shown by {cmd:ctools, environment_check}.

{phang}
{opt update} downloads and installs the latest version of ctools from GitHub.
You may need to restart Stata after updating for all changes to take effect.


{marker threads}{...}
{title:Threads}

{pstd}
By default, ctools uses one thread per CPU it is allowed to run on, but at most
16. The first of the following that is set determines the thread limit:

{p 8 12 2}1. the {opt threads(#)} option of the command;{p_end}
{p 8 12 2}2. the session setting, {cmd:ctools, threads(#)};{p_end}
{p 8 12 2}3. the environment variable {cmd:CTOOLS_NUM_THREADS}, intended for
administrators and batch job scripts;{p_end}
{p 8 12 2}4. the saved setting, {cmd:ctools, threads(#) permanently};{p_end}
{p 8 12 2}5. the environment variable {cmd:OMP_NUM_THREADS} (its first value);{p_end}
{p 8 12 2}6. the number of usable CPUs, at most 16.{p_end}

{pstd}
In every case ctools uses no more threads than there are CPUs it may run on,
and respects {cmd:OMP_THREAD_LIMIT}. On Linux, the usable CPUs follow the
process's CPU affinity (for example {cmd:taskset}, Slurm or PBS allocations)
and container CPU limits. To use more than 16 threads on a large machine, set
a higher value with any of the methods above.

{pstd}
{cmd:ctools, threads()} and {cmd:ctools, environment_check} store
{cmd:r(threads)} (threads used by commands without {opt threads()}),
{cmd:r(threads_default)}, {cmd:r(threads_cap)}, {cmd:r(threads_limit)}
(most threads any setting can request here), {cmd:r(cpus_detected)},
{cmd:r(cpus_available)}, and the locals {cmd:r(threads_source)},
{cmd:r(threads_session)}, {cmd:r(threads_saved)}, and {cmd:r(threads_file)}.

{marker commands}{...}
{title:Available Commands}

{pstd}
{ul:Data Management}

{p2colset 5 20 22 2}{...}
{p2col:{help cimport}}Import text, Excel, dBase, SAS, SPSS, and shapefiles{p_end}
{p2col:{help cexport}}Export text, Excel, dBase, SAS, SPSS, and shapefiles{p_end}
{p2col:{help cdistinct}}Report distinct values{p_end}
{p2col:{help chdfe}}Remove fixed effects from variables{p_end}
{p2col:{help csort}}Sort dataset{p_end}
{p2col:{help cmerge}}Merge (join) datasets{p_end}
{p2col:{help csample}}Resampling without replacement{p_end}
{p2col:{help cbsample}}Resampling with replacement{p_end}
{p2col:{help cencode}}Recast string as labeled numeric{p_end}
{p2col:{help cdecode}}Recast labeled numeric as string{p_end}
{p2col:{help cdestring}}Recast string as numeric type{p_end}
{p2col:{help cwinsor}}Winsorize variables{p_end}
{p2col:{help crangestat}}Range statistics of variables{p_end}

{p2col:{help cinterpolate}}Linear and cubic interpolation{p_end}
{p2col:{help csplit}}Split strings into fields{p_end}
{p2col:{help crangejoin}}Inclusive interval joins{p_end}

{pstd}
{ul:Estimation}

{p2col:{help creghdfe}}OLS with multi-way fixed effects{p_end}
{p2col:{help civreghdfe}}2SLS/GMM with multi-way fixed effects{p_end}
{p2col:{help cppmlhdfe}}PPML with multi-way fixed effects{p_end}
{p2col:{help cqreg}}Quantile regression{p_end}
{p2col:{help cpsmatch}}Propensity score matching{p_end}

{pstd}
{ul:Visualization}

{p2col:{help cbinscatter}}Binned scatter plots{p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Examples}

{phang}{cmd:. csort myvar}{p_end}
{phang}{cmd:. csort var1 var2 var3, verbose}{p_end}
{phang}{cmd:. cmerge m:1 id using lookup.dta}{p_end}
{phang}{cmd:. cimport delimited using data.csv, clear}{p_end}
{phang}{cmd:. cexport delimited using output.csv, replace}{p_end}
{phang}{cmd:. creghdfe y x1 x2, absorb(firm year)}{p_end}

{pstd}Thread settings:{p_end}
{phang}{cmd:. ctools, environment_check}{p_end}
{phang}{cmd:. ctools, threads(8)}{p_end}
{phang}{cmd:. ctools, threads(32) permanently}{p_end}
{phang}{cmd:. ctools, threads(default) permanently}{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
The package draws on the following commands and algorithms:

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
{help cbinscatter}, {help cppmlhdfe}, {help cinterpolate}, {help csplit},
{help crangejoin}
{p_end}
