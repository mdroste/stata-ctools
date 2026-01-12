{smcl}
{* *! version 2.0.0}{...}
{viewerjumpto "Syntax" "cmerge##syntax"}{...}
{viewerjumpto "Description" "cmerge##description"}{...}
{viewerjumpto "Options" "cmerge##options"}{...}
{viewerjumpto "Remarks" "cmerge##remarks"}{...}
{viewerjumpto "Examples" "cmerge##examples"}{...}
{viewerjumpto "Stored results" "cmerge##results"}{...}
{title:Title}

{phang}
{bf:cmerge} {hline 2} C-accelerated merge for Stata datasets


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cmerge}
{it:merge_type}
{varlist}
{cmd:using}
{it:filename}
[{cmd:,} {it:options}]

{pstd}
where {it:merge_type} is one of

{p2colset 9 22 24 2}{...}
{p2col :{opt 1:1}}one-to-one merge{p_end}
{p2col :{opt m:1}}many-to-one merge{p_end}
{p2col :{opt 1:m}}one-to-many merge{p_end}
{p2col :{opt m:m}}many-to-many merge (not recommended){p_end}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt keep(results)}}which observations to keep from merged data{p_end}
{synopt:{opt assert(results)}}verify that merge results are as expected{p_end}
{synopt:{opt gen:erate(varname)}}name of variable to mark merge results; default is {cmd:_merge}{p_end}
{synopt:{opt nogen:erate}}do not create {cmd:_merge} variable{p_end}
{synopt:{opt keepus:ing(varlist)}}variables to keep from using data{p_end}
{synopt:{opt sorted}}assert that both datasets are sorted on key variables{p_end}
{synopt:{opt force}}allow string/numeric variable type mismatches{p_end}
{synopt:{opt norep:ort}}do not display result summary{p_end}
{synopt:{opt verbose}}display detailed timing and progress information{p_end}
{synopt:{opt nolabel}}do not copy value labels from using data{p_end}
{synopt:{opt nonotes}}do not copy variable notes from using data{p_end}
{synopt:{opt update}}update missing values of same-named variables with using data{p_end}
{synopt:{opt replace}}replace all values of same-named variables with using data{p_end}
{synoptline}

{pstd}
{it:results} is one or more of: {opt match}, {opt master}, {opt using}, {opt 1}, {opt 2}, {opt 3}


{marker description}{...}
{title:Description}

{pstd}
{cmd:cmerge} is a high-performance drop-in replacement for Stata's {help merge:merge}
command. It uses a C plugin with parallel data loading and optimized sorting
to achieve significant speed improvements over the native command.

{pstd}
{cmd:cmerge} joins the dataset currently in memory (the master data) with
{it:filename} (the using data), matching on the key variables specified in
{it:varlist}. The result replaces the data in memory.


{marker options}{...}
{title:Options}

{phang}
{opt keep(results)} specifies which observations are to be kept from the
merged data. Valid results are {opt match} (or {opt 3}), {opt master} (or {opt 1}),
and {opt using} (or {opt 2}).

{phang}
{opt assert(results)} specifies assertions about the match results that should
be verified. If any assertion fails, an error is raised.

{phang}
{opt generate(varname)} specifies the name of the variable to be created
marking the source of each observation. The default name is {cmd:_merge}.

{phang}
{opt nogenerate} specifies that the merge-result variable should not be created.

{phang}
{opt keepusing(varlist)} specifies which variables from the using dataset
should be kept in the merged result.

{phang}
{opt sorted} specifies that both datasets are already sorted on the key
variables. This skips the internal sorting step and can improve performance
for pre-sorted data.

{phang}
{opt force} allows merging when key variables have different types (string vs.
numeric) in master and using datasets.

{phang}
{opt noreport} suppresses the table showing the merge result summary.

{phang}
{opt verbose} displays detailed timing information and progress messages
during the merge operation.

{phang}
{opt nolabel} prevents value labels from being copied from the using dataset.
By default, value labels attached to variables in the using dataset are copied
to the merged result.

{phang}
{opt nonotes} prevents variable notes (characteristics) from being copied from
the using dataset. By default, variable notes are preserved.

{phang}
{opt update} specifies that for observations that match, missing values of
variables that exist in both datasets should be updated with corresponding
non-missing values from the using dataset. Variables unique to using are
always added.

{phang}
{opt replace} specifies that for observations that match, all values of
variables that exist in both datasets should be replaced with corresponding
values from the using dataset. This is more aggressive than {opt update}
as it replaces non-missing values too.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:cmerge} performs all merge operations entirely in C for maximum speed.
The algorithm uses:

{p 8 12 2}1. Parallel data loading into C memory{p_end}
{p 8 12 2}2. High-performance radix sort on key variables{p_end}
{p 8 12 2}3. Single-pass sorted merge join{p_end}
{p 8 12 2}4. Parallel output construction{p_end}

{pstd}
For large datasets (millions of observations), {cmd:cmerge} can be 3-5x faster
than the native {cmd:merge} command.


{marker examples}{...}
{title:Examples}

{pstd}Simple one-to-one merge:{p_end}
{phang2}{cmd:. cmerge 1:1 id using data2.dta}{p_end}

{pstd}Many-to-one merge keeping only matched observations:{p_end}
{phang2}{cmd:. cmerge m:1 state year using statedata.dta, keep(match)}{p_end}

{pstd}Merge with verbose output and no merge variable:{p_end}
{phang2}{cmd:. cmerge m:1 id using lookup.dta, nogenerate verbose}{p_end}

{pstd}Keep only selected variables from using dataset:{p_end}
{phang2}{cmd:. cmerge 1:1 id using fulldata.dta, keepusing(var1 var2)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cmerge} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations after merge{p_end}
{synopt:{cmd:r(N_1)}}number of observations only in master{p_end}
{synopt:{cmd:r(N_2)}}number of observations only in using{p_end}
{synopt:{cmd:r(N_3)}}number of matched observations{p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(using)}}name of using file{p_end}
{synopt:{cmd:r(keyvars)}}key variable names{p_end}


{title:Author}

{pstd}
ctools package


{title:Also see}

{psee}
{space 2}Help: {help merge}, {help ctools}, {help csort}
{p_end}
