{smcl}
{* *! version 0.9.0 26Jan2026}{...}
{viewerjumpto "Syntax" "csort##syntax"}{...}
{viewerjumpto "Description" "csort##description"}{...}
{viewerjumpto "Options" "csort##options"}{...}
{viewerjumpto "Algorithms" "csort##algorithms"}{...}
{viewerjumpto "Remarks" "csort##remarks"}{...}
{viewerjumpto "Examples" "csort##examples"}{...}
{viewerjumpto "Stored results" "csort##results"}{...}
{title:Title}

{phang}
{bf:csort} {hline 2} C-accelerated parallel sorting for Stata datasets


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:csort}
{varlist}
{ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt alg:orithm(name)}}sorting algorithm (see below){p_end}
{synopt:{opt str:eam(#)}}streaming mode, loading {it:#} variables at a time (1-16){p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt nosort:edby}}do not set Stata's sortedby attribute{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}

{pstd}
where {it:name} is one of: {opt auto} (default), {opt counting}, {opt lsd}, {opt msd},
{opt ips4o}, {opt timsort}, {opt sample}, or {opt merge}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csort} is a high-performance drop-in replacement for Stata's {help sort:sort}
command. It uses a C plugin with parallel algorithms to achieve significant
speed improvements over the native command, especially for large datasets.

{pstd}
{cmd:csort} sorts the observations of the current dataset by the values of the
variables in {varlist}. The sort is {it:stable}, meaning observations with equal
key values maintain their relative order.


{marker options}{...}
{title:Options}

{phang}
{opt algorithm(name)} specifies which sorting algorithm to use. The default is
{opt auto}, which intelligently selects the optimal algorithm based on your data:
string variables use MSD radix sort, while numeric variables use counting sort
(with automatic fallback to LSD radix for large ranges). This provides excellent
performance for common Stata workloads without manual tuning. See
{help csort##algorithms:Algorithms} for details on manual algorithm selection.

{phang}
{opt stream(#)} enables streaming mode, which reduces memory usage for wide datasets.
In streaming mode, only the sort key variables are loaded into C memory. After
sorting, the permutation is applied to non-key variables using sequential Stata
I/O. This is useful when sorting datasets with many columns that would otherwise
exceed available memory.

{pmore}
The argument {it:#} specifies how many variables to process at a time (1-16).
Values above 16 are capped at 16.
Higher values use more memory but may improve performance by allowing better
parallelization. For example, {cmd:stream(4)} processes 4 variables simultaneously,
using approximately 4× the buffer memory of {cmd:stream(1)}.

{phang}
{opt threads(#)} specifies the maximum number of threads to use for parallel
operations. By default, {cmd:csort} uses all available CPU cores as reported by
OpenMP. Use this option to limit parallelism, for example when running multiple
jobs simultaneously or to reduce resource usage.

{phang}
{opt nosortedby} prevents {cmd:csort} from setting Stata's internal {it:sortedby}
attribute after sorting. Normally, {cmd:csort} calls Stata's {cmd:sort} command
at the end to register the sort order. With {opt nosortedby}, this step is
skipped, which can be faster when Stata does not need to recognize the data as
sorted (e.g., for intermediate operations). Note: commands that rely on Stata's
{it:sortedby} attribute (like {cmd:by:}) will not recognize the sort order.

{phang}
{opt verbose} displays a timing breakdown showing time spent in each phase of
the sort: loading data to C, sorting, and storing data back to Stata.


{marker algorithms}{...}
{title:Algorithms}

{pstd}
By default, {cmd:csort} automatically selects the optimal sorting algorithm based
on your data characteristics. This usually provides the best performance without
manual tuning. You can override the default with the {opt algorithm()} option.

{dlgtab:Auto (default)}

{pstd}
{opt algorithm(auto)} intelligently selects the best algorithm:

{p 8 12 2}- {bf:String variables}: Uses MSD radix sort, which handles variable-length
strings efficiently by processing from the most significant byte.{p_end}

{p 8 12 2}- {bf:Numeric variables}: Uses counting sort for O(n+k) performance when
the value range is small (common for dates, years, categories). Automatically
falls back to LSD radix sort for larger ranges.{p_end}

{pstd}
This provides massive speedups for typical Stata workloads. For example, sorting
panel data by date or ID can be 100x faster than comparison-based sorts.

{dlgtab:Manual Algorithm Selection}

{pstd}
For specialized use cases, you can manually select an algorithm:

{p2colset 5 28 30 2}{...}
{p2col :{opt counting}}Integer data with small range (dates, years, categories). O(n+k) performance.
Falls back to LSD radix if range is too large.{p_end}
{p2col :{opt lsd}}Numeric data. Processes bytes from least to most significant.{p_end}
{p2col :{opt msd}}String data. Processes from most significant byte, can skip trailing bytes.{p_end}
{p2col :{opt ips4o}}General-purpose parallel samplesort. Good default for mixed workloads.{p_end}
{p2col :{opt timsort}}Partially sorted data. O(n) for nearly sorted, O(n log n) worst case.{p_end}
{p2col :{opt sample}}Very large datasets with many cores. Near-linear parallel scaling.{p_end}
{p2col :{opt merge}}Guaranteed stable O(n log n). Predictable performance.{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:csort} performs all operations in C using parallelized algorithms. The
implementation uses:

{p 8 12 2}1. Parallel data loading with OpenMP{p_end}
{p 8 12 2}2. Optimized sort algorithm (selected via {opt algorithm()} option){p_end}
{p 8 12 2}3. Parallel data storing back to Stata{p_end}

{pstd}
For large datasets (millions of observations), {cmd:csort} can be 2-10x faster
than the native {cmd:sort} command, depending on data characteristics and the
number of available CPU cores.

{pstd}
{cmd:csort} supports sorting on multiple variables. When sorting on multiple
variables, the first variable in {varlist} is the primary sort key, the second
variable is the secondary key, and so forth.

{pstd}
Both numeric and string variables are supported. Numeric variables are sorted
in ascending numerical order. String variables are sorted in ascending
lexicographic (ASCII) order.

{pstd}
{bf:Limitations:} {cmd:csort} does not support datasets exceeding 2^31 - 1
(2,147,483,647) observations.


{marker examples}{...}
{title:Examples}

{pstd}Sort by a single numeric variable:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csort price}{p_end}

{pstd}Sort by a string variable:{p_end}
{phang2}{cmd:. csort make}{p_end}

{pstd}Sort by multiple variables (primary and secondary keys):{p_end}
{phang2}{cmd:. csort foreign mpg}{p_end}

{pstd}Sort with verbose output showing timing:{p_end}
{phang2}{cmd:. csort price, verbose}{p_end}

{pstd}Force a specific algorithm (normally not needed):{p_end}
{phang2}{cmd:. csort year, algorithm(counting)}{p_end}
{phang2}{cmd:. csort make, algorithm(msd)}{p_end}
{phang2}{cmd:. csort idcode, algorithm(ips4o)}{p_end}

{pstd}Sort with detailed progress information:{p_end}
{phang2}{cmd:. csort price mpg, verbose}{p_end}

{pstd}Use streaming mode for wide datasets with limited memory:{p_end}
{phang2}{cmd:. csort idcode, stream(1)}{p_end}

{pstd}Use streaming mode with 4 variables processed at a time:{p_end}
{phang2}{cmd:. csort idcode, stream(4)}{p_end}

{pstd}Limit parallelism to 4 threads:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. csort price, threads(4)}{p_end}

{pstd}Sort without setting Stata's sortedby attribute (faster for intermediate ops):{p_end}
{phang2}{cmd:. csort price, nosortedby}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:csort} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations sorted{p_end}
{synopt:{cmd:r(time_load)}}time to load data from Stata to C (seconds){p_end}
{synopt:{cmd:r(time_sort)}}time to sort data (seconds){p_end}
{synopt:{cmd:r(time_permute)}}time to apply permutation (seconds){p_end}
{synopt:{cmd:r(time_store)}}time to store data back to Stata (seconds){p_end}
{synopt:{cmd:r(time_stream)}}time to stream non-key variables (streaming mode only){p_end}
{synopt:{cmd:r(time_cleanup)}}time to free C memory (seconds){p_end}
{synopt:{cmd:r(time_total)}}total C plugin time (seconds){p_end}
{synopt:{cmd:r(stream)}}1 if streaming mode was used, 0 otherwise{p_end}
{synopt:{cmd:r(threads_max)}}maximum threads available{p_end}
{synopt:{cmd:r(openmp_enabled)}}1 if OpenMP is enabled, 0 otherwise{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(algorithm)}}sorting algorithm used (auto, counting, lsd, msd, ips4o, timsort, sample, or merge){p_end}
{synopt:{cmd:r(sortvars)}}variables used for sorting{p_end}
{synopt:{cmd:r(cmd)}}{cmd:csort}{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] sort}

{psee}
Online: {help sort}, {help gsort}, {help ctools}, {help cmerge}
{p_end}
