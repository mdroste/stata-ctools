{smcl}
{* *! version 1.1.0 20261016}{...}
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
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt alg:orithm(name)}}sorting algorithm (see below){p_end}
{synopt:{opt str:eam(#)}}streaming mode, loading {it:#} variables at a time (1-16){p_end}
{synopt:{opt nostr:eam}}disable automatic streaming for wide datasets{p_end}
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
{cmd:csort} sorts the observations of the current dataset by the values of the
variables in {varlist}. The sort is {it:stable}, meaning observations with equal
key values maintain their relative order.


{pstd}
Like native {cmd:sort}, {cmd:csort} rejects {cmd:if} and {cmd:in}. Fixed-width
strings are supported; {cmd:strL} variables are rejected.

{marker options}{...}
{title:Options}

{phang}
{opt algorithm(name)} specifies which sorting algorithm to use. The default is
{opt auto}, which selects an algorithm from the key types and ranges:
numeric variables use counting sort (with automatic fallback to LSD radix for
large ranges), string variables of up to 16 bytes, alone or with numeric
variables, use LSD radix sort, and wider string variables use MSD radix sort.
See
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
{opt nostream} disables automatic streaming mode. By default, {cmd:csort}
automatically enables streaming for datasets with many variables to reduce
memory usage. Use {opt nostream} to force all variables to be loaded into C
memory simultaneously.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt nosortedby} prevents {cmd:csort} from setting Stata's internal {it:sortedby}
attribute after sorting. Normally, {cmd:csort} calls Stata's {cmd:sort} command
at the end to register the sort order. With {opt nosortedby}, this step is
skipped, which can be faster when Stata does not need to recognize the data as
sorted (e.g., for intermediate operations). Any previous sort order is cleared,
so commands that rely on Stata's {it:sortedby} attribute (like {cmd:by:}) will
not treat the data as sorted.

{phang}
{opt verbose} displays progress and elapsed times.


{marker algorithms}{...}
{title:Algorithms}

{pstd}
Use {opt algorithm()} to select an algorithm. The default is {opt auto}.

{dlgtab:Auto (default)}

{pstd}
{opt algorithm(auto)} uses the following rules:

{p 8 12 2}- {bf:Numeric variables}: Uses counting sort for O(n+k) performance when
the value range is small (common for dates, years, categories). Automatically
falls back to LSD radix sort for larger ranges.{p_end}

{p 8 12 2}- {bf:String variables of up to 16 bytes} ({cmd:str1}-{cmd:str16}), alone
or with numeric variables: Uses LSD radix sort, which reads each string as one or
two 8-byte words and sorts them together with the numeric keys, comparing later
keys only for observations tied on the earlier ones.{p_end}

{p 8 12 2}- {bf:Wider string variables}: Uses MSD radix sort, which handles long
variable-length strings efficiently by processing from the most significant
byte.{p_end}

{pstd}
The benefit depends on key ranges, dataset shape, and hardware.

{dlgtab:Manual Algorithm Selection}

{pstd}
For specialized use cases, you can manually select an algorithm:

{p2colset 5 28 30 2}{...}
{p2col :{opt counting}}Integer data with small range (dates, years, categories). O(n+k) performance.
Falls back to LSD radix if range is too large.{p_end}
{p2col :{opt lsd}}Numeric data and strings of up to 16 bytes. Processes bytes from least to most significant.{p_end}
{p2col :{opt msd}}Long string data. Processes from most significant byte, can skip trailing bytes.{p_end}
{p2col :{opt ips4o}}General-purpose parallel samplesort.{p_end}
{p2col :{opt timsort}}Comparison sort that uses existing ordered runs; O(n log n) worst case.{p_end}
{p2col :{opt sample}}Parallel sample sort.{p_end}
{p2col :{opt merge}}Guaranteed stable O(n log n). Predictable performance.{p_end}


{marker remarks}{...}
{title:Remarks}


{pstd}
{cmd:csort} supports sorting on multiple variables. When sorting on multiple
variables, the first variable in {varlist} is the primary sort key, the second
variable is the secondary key, and so forth.

{pstd}
Both numeric and string variables are supported. Numeric variables are sorted
in ascending numerical order. String variables are sorted in ascending
byte order. Numeric missing values sort after finite numbers, in the order
{cmd:.}, {cmd:.a}, ..., {cmd:.z}; empty strings sort first.

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
{phang2}{cmd:. csort rep78, algorithm(counting)}{p_end}
{phang2}{cmd:. csort make, algorithm(msd)}{p_end}
{phang2}{cmd:. csort make, algorithm(ips4o)}{p_end}

{pstd}Sort with detailed progress information:{p_end}
{phang2}{cmd:. csort price mpg, verbose}{p_end}

{pstd}Use streaming mode for wide datasets with limited memory:{p_end}
{phang2}{cmd:. csort make, stream(1)}{p_end}

{pstd}Use streaming mode with 4 variables processed at a time:{p_end}
{phang2}{cmd:. csort make, stream(4)}{p_end}

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
{synopt:{cmd:r(algorithm)}}requested algorithm (auto, counting, lsd, msd, ips4o, timsort, sample, or merge){p_end}
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
