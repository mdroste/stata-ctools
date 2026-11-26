{smcl}
{* *! version 1.3.1 20261124}{...}
{viewerjumpto "Syntax" "ctools##syntax"}{...}
{viewerjumpto "Description" "ctools##description"}{...}
{viewerjumpto "Options" "ctools##options"}{...}
{viewerjumpto "Threads" "ctools##threads"}{...}
{viewerjumpto "Commands" "ctools##commands"}{...}
{viewerjumpto "Examples" "ctools##examples"}{...}
{viewerjumpto "Acknowledgments" "ctools##acknowledgments"}{...}
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
{opt version} displays the ado and plugin versions, the plugin revision and
the plugin protocol. It stores them in {cmd:r(version)},
{cmd:r(plugin_version)}, {cmd:r(plugin_revision)}, and
{cmd:r(plugin_protocol)}. Every ctools command first checks that the plugin's
protocol matches the one the ado files need; a plugin from another release
stops with an error that asks you to reinstall ctools and restart Stata.

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
Restart Stata after updating, because a loaded plugin is replaced only when
Stata restarts. If GitHub cannot be reached or the installation fails,
{cmd:ctools, update} reports the failure and exits with {cmd:net install}'s
return code.


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
When global {cmd:CTOOLS_THREADS} is invalid, commands fail until it is
cleared, so {cmd:r(threads)} is missing and {cmd:r(threads_session)} is
{cmd:invalid}; {cmd:r(threads_default)} reports the count that applies once
it is cleared.

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


{* ctools-attribution:begin (generated by scripts/sync_attribution.py; edit config/attribution.json)}{...}
{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
ctools reimplements existing Stata commands in C. Each ctools command takes
its syntax, options, output, and stored results from the command it replaces,
and is pinned to one version of that command. The Acknowledgments section of
each command's help file has details and references.

{pstd}
{ul:Community-contributed commands}

{p2colset 5 20 22 2}{...}
{p2col:{help cdistinct##acknowledgments:cdistinct}}{cmd:distinct} 1.4.3
(September 3, 2023) by Gary Longton and Nicholas J. Cox{p_end}
{p2col:{help crangestat##acknowledgments:crangestat}}{cmd:rangestat} 1.1.1
(May 9, 2017) by Robert Picard, Nicholas J. Cox, and Roberto Ferrer{p_end}
{p2col:{help crangejoin##acknowledgments:crangejoin}}{cmd:rangejoin} 1.1.3
(April 13, 2021) by Robert Picard{p_end}
{p2col:{help cinterpolate##acknowledgments:cinterpolate}}Stata's {cmd:ipolate}
({cmd:ipolate.ado} 1.3.5, April 23, 2020); {cmd:cipolate} 1.0.0 (July 4, 2002)
by Nicholas J. Cox{p_end}
{p2col:{help cwinsor##acknowledgments:cwinsor}}{cmd:winsor2} 1.1 (December 16,
2014) by Yujun Lian; {cmd:gstats winsor} from gtools 1.10.1 (December 5, 2022)
by Mauricio Caceres Bravo{p_end}
{p2col:{help cbinscatter##acknowledgments:cbinscatter}}{cmd:binscatter} 7.02
(November 24, 2013) by Michael Stepner; {cmd:binsreg} 2.2 (August 20, 2026) by
Matias D. Cattaneo, Richard K. Crump, Max H. Farrell, and Yingjie Feng{p_end}
{p2col:{help cpsmatch##acknowledgments:cpsmatch}}{cmd:psmatch2} 4.1.0 (June
11, 2026) by Edwin Leuven and Barbara Sianesi{p_end}
{p2col:{help creghdfe##acknowledgments:creghdfe}}{cmd:reghdfe} 6.13.1 (January
10, 2026) by Sergio Correia and Noah Constantine{p_end}
{p2col:{help civreghdfe##acknowledgments:civreghdfe}}{cmd:ivreghdfe} 1.1.4
(November 29, 2025) by Sergio Correia; {cmd:ivreg2} 4.1.12 (August 14, 2024)
by Christopher F Baum, Mark E Schaffer, and Steven Stillman{p_end}
{p2col:{help cppmlhdfe##acknowledgments:cppmlhdfe}}{cmd:ppmlhdfe} 2.3.3
(November 2, 2025) by Sergio Correia, Paulo Guimarães, and Thomas
Zylkin{p_end}
{p2col:{help chdfe##acknowledgments:chdfe}}{cmd:hdfe} 3.2.9 (February 21,
2016; SSC distribution of August 6, 2016) by Sergio Correia{p_end}

{pstd}
{ul:Official Stata commands} (StataNow 19.5, revision 14 Oct 2026)

{p2col:{help csplit##acknowledgments:csplit}}Stata's {cmd:split}
({cmd:split.ado} 2.1.0, November 2, 2021), written by Nicholas J. Cox{p_end}
{p2col:{help cdestring##acknowledgments:cdestring}}Stata's {cmd:destring}
({cmd:destring.ado} 2.4.4, September 15, 2016), originally written by Nicholas
J. Cox{p_end}
{p2col:{help csort##acknowledgments:csort}}Stata's {cmd:sort} (built
in){p_end}
{p2col:{help cmerge##acknowledgments:cmerge}}Stata's {cmd:merge}
({cmd:merge.ado} 3.2.13, September 12, 2023){p_end}
{p2col:{help cimport##acknowledgments:cimport}}Stata's {cmd:import delimited}
({cmd:import_delimited.ado} 1.3.5, April 18, 2022); Stata's {cmd:import excel}
({cmd:import_excel.ado} 1.0.15, September 20, 2019); Stata's other
{cmd:import} commands{p_end}
{p2col:{help cexport##acknowledgments:cexport}}Stata's {cmd:export delimited}
({cmd:export_delimited.ado} 1.0.9, June 13, 2020); Stata's {cmd:export excel}
({cmd:export_excel.ado} 1.1.1, July 30, 2018); Stata's other {cmd:export}
commands{p_end}
{p2col:{help csample##acknowledgments:csample}}Stata's {cmd:sample}
({cmd:sample.ado} 2.4.0, January 20, 2021){p_end}
{p2col:{help cbsample##acknowledgments:cbsample}}Stata's {cmd:bsample}
({cmd:bsample.ado} 2.0.8, January 17, 2025){p_end}
{p2col:{help cencode##acknowledgments:cencode}}Stata's {cmd:encode} (built
in){p_end}
{p2col:{help cdecode##acknowledgments:cdecode}}Stata's {cmd:decode} (built
in){p_end}
{p2col:{help cqreg##acknowledgments:cqreg}}Stata's {cmd:qreg} ({cmd:qreg.ado}
4.4.2, May 8, 2020){p_end}

{pstd}
{ul:Supporting packages}

{p2col:{cmd:ftools}}2.50.0 (January 9, 2026) by Sergio Correia. Required by
{cmd:reghdfe}, {cmd:ivreghdfe}, and {cmd:ppmlhdfe}.{p_end}
{p2col:{cmd:reghdfe}}6.13.1 (January 10, 2026) by Sergio Correia and Noah
Constantine. Fixed-effects solver used by {cmd:ivreghdfe} and
{cmd:ppmlhdfe}.{p_end}
{p2col:{cmd:ranktest}}2.0.04 (September 21, 2020) by Frank Kleibergen, Mark E
Schaffer, and Frank Windmeijer. Rank tests reported by {cmd:ivreg2} and
{cmd:ivreghdfe}.{p_end}
{p2col:{cmd:boottest}}4.5.5 (October 24, 2026) by David Roodman. Reference for
the wild bootstrap p-values in {cmd:creghdfe} and {cmd:civreghdfe}.{p_end}
{p2colreset}{...}

{pstd}
Versions were recorded on November 25, 2026.

{pstd}
Please report bugs and questions about ctools commands at
{browse "https://github.com/mdroste/stata-ctools/issues":github.com/mdroste/stata-ctools/issues}.
The authors of the original commands do not maintain ctools. If you use a
ctools command in published work, please cite the original command.
{browse "https://github.com/mdroste/stata-ctools/blob/main/docs/ATTRIBUTION.md":docs/ATTRIBUTION.md}
has file checksums and all references.

{pstd}
{ul:Other acknowledgments}

{p 8 10 2}- Michael Axtmann, Sascha Witt, Daniel Ferizovic, and Peter Sanders for the
In-place Parallel Super Scalar Samplesort
({browse "https://github.com/SaschaWitt/ips4o":IPS4o}) algorithm used in
{cmd:csort}.{p_end}
{p 8 10 2}- Sergio Correia for {browse "https://github.com/sergiocorreia/ftools":ftools}
and Mauricio Caceres Bravo for
{browse "https://github.com/mcaceresb/gtools":gtools}.{p_end}
{p 8 10 2}- Lars Vilhuber for contributions to {cmd:ivreghdfe}.{p_end}
{p 8 10 2}- The authors of the libraries compiled into the ctools plugins, listed in
{cmd:help ctools_licenses}.{p_end}
{* ctools-attribution:end}{...}


{title:Also see}

{psee}
Online: {help cimport}, {help cexport}, {help csort}, {help cmerge}, {help cencode},
{help cdecode}, {help cdestring}, {help csample}, {help cbsample}, {help cwinsor},
{help crangestat}, {help creghdfe}, {help civreghdfe}, {help cqreg}, {help cpsmatch},
{help cbinscatter}, {help cppmlhdfe}, {help cinterpolate}, {help csplit},
{help crangejoin}, {help cdistinct}, {help chdfe}
{p_end}

{psee}
Licenses: {help ctools_licenses} (ctools and the libraries in its plugins)
{p_end}
