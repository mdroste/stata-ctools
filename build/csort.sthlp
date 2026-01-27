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

{pstd}
{cmd:csort} is a high-performance drop-in replacement for Stata's {help sort:sort}
command.


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
{synopt:{opt str:eam}}streaming mode for reduced memory usage{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt nosort:edby}}do not set Stata's sortedby attribute{p_end}
{synopt:{opt v:erbose}}display timing breakdown{p_end}
{synoptline}

{pstd}
where {it:name} is one of: {opt ips4o} (default), {opt lsd}, {opt msd}, {opt timsort},
{opt sample}, {opt counting}, or {opt merge}


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
{opt ips4o} (In-place Parallel Super Scalar Samplesort), a state-of-the-art
parallel sorting algorithm. See {help csort##algorithms:Algorithms} for details
on when to use each algorithm.

{phang}
{opt stream} enables streaming mode, which reduces memory usage for wide datasets.
In streaming mode, only the sort key variables are loaded into C memory. After
sorting, the permutation is applied to non-key variables one at a time using
sequential Stata I/O. This is useful when sorting datasets with many columns
that would otherwise exceed available memory. Performance is comparable to
standard mode for most datasets.

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
{cmd:csort} provides seven sorting algorithms optimized for different data patterns
and parallelization strategies:

{dlgtab:IPS4o (default)}

{pstd}
{opt algorithm(ips4o)} uses IPS4o (In-place Parallel Super Scalar Samplesort),
a state-of-the-art parallel sorting algorithm that combines samplesort with
branchless bucket classification and cache-efficient block-based data movement.

{p 8 12 2}{it:Best for:} General-purpose sorting, large datasets, high core counts{p_end}
{p 8 12 2}{it:Performance:} O(n log n) with excellent cache efficiency and parallel scaling{p_end}
{p 8 12 2}{it:Parallelization:} Fully parallel classification and bucket sorting{p_end}
{p 8 12 2}{it:Key features:} Branchless operations, minimal auxiliary memory, adaptive bucket count{p_end}

{dlgtab:LSD Radix Sort}

{pstd}
{opt algorithm(lsd)} uses Least Significant Digit radix sort. This is a highly
parallelized sorting algorithm that processes data byte-by-byte from the least
significant to most significant byte.

{p 8 12 2}{it:Best for:} Fixed-width numeric keys, uniformly distributed data{p_end}
{p 8 12 2}{it:Performance:} O(n * k) where k is key width; excellent cache utilization{p_end}
{p 8 12 2}{it:Parallelization:} Fully parallel histogram and scatter phases{p_end}

{dlgtab:MSD Radix Sort}

{pstd}
{opt algorithm(msd)} uses Most Significant Digit radix sort. This algorithm
processes data from the most significant byte first, allowing early termination
for strings that differ in their prefixes.

{p 8 12 2}{it:Best for:} Variable-length strings, data with common prefixes, short strings{p_end}
{p 8 12 2}{it:Performance:} O(n * k) but can skip trailing bytes; adaptive to string lengths{p_end}
{p 8 12 2}{it:Parallelization:} Parallel histogram computation and partition scatter{p_end}

{dlgtab:Timsort}

{pstd}
{opt algorithm(timsort)} uses Timsort, an adaptive hybrid algorithm that
combines merge sort and insertion sort. It detects and exploits existing runs
(pre-sorted sequences) in the data.

{p 8 12 2}{it:Best for:} Partially sorted data, data with natural runs, nearly sorted data{p_end}
{p 8 12 2}{it:Performance:} O(n log n) worst case, O(n) for nearly sorted data{p_end}
{p 8 12 2}{it:Parallelization:} Limited; primarily sequential{p_end}

{dlgtab:Sample Sort}

{pstd}
{opt algorithm(sample)} uses parallel sample sort, the gold standard for parallel
sorting. It samples data to select splitters, partitions into buckets, then sorts
each bucket independently in parallel.

{p 8 12 2}{it:Best for:} Very large datasets (millions of observations), high core counts (8+ cores){p_end}
{p 8 12 2}{it:Performance:} O(n log n / p) with p processors; near-linear speedup{p_end}
{p 8 12 2}{it:Parallelization:} Excellent; all phases are parallel{p_end}

{dlgtab:Counting Sort}

{pstd}
{opt algorithm(counting)} uses parallel counting sort, optimal for integer data
with a small range of distinct values. Automatically falls back to LSD radix
if data is not suitable.

{p 8 12 2}{it:Best for:} Integer data with range < 1,000,000 (year, state codes, categorical variables){p_end}
{p 8 12 2}{it:Performance:} O(n + k) where k is range; faster than comparison sorts{p_end}
{p 8 12 2}{it:Parallelization:} Fully parallel histogram and scatter{p_end}
{p 8 12 2}{it:Note:} Not suitable for non-integer or string data{p_end}

{dlgtab:Parallel Merge Sort}

{pstd}
{opt algorithm(merge)} uses parallel merge sort, a stable algorithm with
predictable O(n log n) performance. Divides data into blocks, sorts each
in parallel, then merges.

{p 8 12 2}{it:Best for:} When stable sort is required, predictable performance needed{p_end}
{p 8 12 2}{it:Performance:} O(n log n) guaranteed; stable and predictable{p_end}
{p 8 12 2}{it:Parallelization:} Parallel block sort + parallel merge{p_end}

{dlgtab:Algorithm Selection Guide}

{pstd}
Use the following guidelines to select an algorithm:

{p2colset 5 28 30 2}{...}
{p2col :{opt ips4o}}Default; best for most sorting tasks, excellent parallel scaling{p_end}
{p2col :{opt lsd}}Alternative for numeric data with uniform distribution{p_end}
{p2col :{opt msd}}Better for string variables with common prefixes{p_end}
{p2col :{opt timsort}}When data is already partially sorted{p_end}
{p2col :{opt sample}}Very large datasets with many CPU cores{p_end}
{p2col :{opt counting}}Integer data with small range (year, categories){p_end}
{p2col :{opt merge}}When guaranteed stable O(n log n) is needed{p_end}


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

{pstd}Use MSD radix sort for string sorting:{p_end}
{phang2}{cmd:. csort make, algorithm(msd)}{p_end}

{pstd}Use sample sort for large datasets with many cores:{p_end}
{phang2}{cmd:. csort id, algorithm(sample)}{p_end}

{pstd}Use counting sort for integer data (e.g., year):{p_end}
{phang2}{cmd:. csort year, algorithm(counting)}{p_end}

{pstd}Use parallel merge sort for stable sorting:{p_end}
{phang2}{cmd:. csort id, algorithm(merge)}{p_end}

{pstd}Use IPS4o for large datasets with memory efficiency:{p_end}
{phang2}{cmd:. csort id, algorithm(ips4o)}{p_end}

{pstd}Sort with detailed progress information:{p_end}
{phang2}{cmd:. csort price mpg, verbose}{p_end}

{pstd}Use streaming mode for wide datasets with limited memory:{p_end}
{phang2}{cmd:. csort id, stream}{p_end}

{pstd}Limit parallelism to 4 threads:{p_end}
{phang2}{cmd:. csort price, threads(4)}{p_end}

{pstd}Sort without setting Stata's sortedby attribute (faster for intermediate ops):{p_end}
{phang2}{cmd:. csort id, nosortedby}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
When the {opt verbose} option is specified, {cmd:csort} stores the following
in global macros:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars (globals)}{p_end}
{synopt:{cmd:_csort_time_load}}time to load data from Stata to C (seconds){p_end}
{synopt:{cmd:_csort_time_sort}}time to sort data (seconds){p_end}
{synopt:{cmd:_csort_time_store}}time to store data back to Stata (seconds){p_end}
{synopt:{cmd:_csort_time_stream}}time to stream non-key variables (streaming mode only){p_end}
{synopt:{cmd:_csort_time_total}}total elapsed time (seconds){p_end}


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
