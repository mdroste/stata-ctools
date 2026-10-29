{smcl}
{* *! version 1.1.0 20261016}{...}
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
{p2col :{opt m:m}}match repeated keys sequentially within each key group{p_end}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt keep(results)}}which observations to keep from merged data{p_end}
{synopt:{opt assert(results)}}verify that merge results are as expected{p_end}
{synopt:{opt gen:erate(varname)}}name of variable to mark merge results; default is {cmd:_merge}{p_end}
{synopt:{opt nogen:erate}}do not create {cmd:_merge} variable{p_end}
{synopt:{opt keepus:ing(varlist)}}variables to keep from using data{p_end}
{syntab:Merge behavior}
{synopt:{opt sorted}}assert that both datasets are sorted on key variables{p_end}
{synopt:{opt force}}allow string/numeric mismatches for non-key variables only{p_end}
{synopt:{opt nolabel}}do not copy value labels from using data{p_end}
{synopt:{opt nonotes}}do not copy variable notes from using data{p_end}
{synopt:{opt update}}update missing values of same-named variables with using data{p_end}
{synopt:{opt replace}}replace conflicting nonmissing values; requires {opt update}{p_end}
{synopt:{opt preserve_order(#)}}preserve observation order from master dataset (0=no, 1=yes){p_end}
{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}
{syntab:Reporting}
{synopt:{opt norep:ort}}do not display result summary{p_end}
{synopt:{opt verbose}}display detailed timing and progress information{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{pstd}
{it:results} is one or more of: {opt match}, {opt master}, {opt using}, {opt match_update}, {opt match_conflict}, or codes {opt 1} through {opt 5}


{marker description}{...}
{title:Description}


{pstd}
{cmd:cmerge} joins the dataset currently in memory (the master data) with
{it:filename} (the using data), matching on the key variables specified in
{it:varlist}. The result replaces the data in memory.


{pstd}
For an observation-by-observation merge, use
{cmd:cmerge 1:1 _n using} {it:filename}. For a keyed merge, missing key values
participate in matching. {cmd:m:m} pairs observations sequentially within each
key group and repeats the last observation of the shorter group; it does not
form every pair. Use {help joinby} for all combinations.

{marker options}{...}
{title:Options}

{phang}
{opt keep(results)} specifies which observations are to be kept from the
merged data. As with {help merge}, results are {opt master} (or {opt 1}),
{opt using} (or {opt 2}), {opt match} (or {opt 3}), {opt match_update}
(or {opt 4}), and {opt match_conflict} (or {opt 5}), with the same
abbreviations; the last two arise only with {opt update}.

{phang}
{opt assert(results)} specifies assertions about the match results that should
be verified; results are specified as in {opt keep()}. If any assertion fails,
an error is raised (r(9)) and the merged data remain in memory. The assertion
is evaluated before {opt keep()} restricts the result.

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
for pre-sorted data. As with {help merge}, data that are not sorted give an
error (r(5)).

{phang}
{opt force} allows incompatible string/numeric types in shared non-key variables.
The master storage type is retained, and incompatible using-side values are
missing. Key variables must have matching string/numeric types even with
{opt force}; otherwise the command returns error 106 before modifying data.
{p_end}

{phang}
{opt nolabel} prevents copying value-label definitions from the using data.
Without this option, labels are copied, but an existing master definition
with the same name takes precedence.

{phang}
{opt nonotes} prevents copying notes from the using data.

{phang}
{opt update} fills missing master values of shared non-key variables with
nonmissing using values. Matched observations receive code 4 if any missing
value was filled, or code 5 if any nonmissing values conflict; code 5 takes
precedence. Other matched observations receive code 3.

{phang}
{opt replace}, with {opt update}, also replaces conflicting nonmissing master
values with nonmissing using values. Missing using values do not replace
nonmissing master values.

{phang}
{opt preserve_order(1)} retains master observation order, with using-only
observations appended. The default is {cmd:preserve_order(0)}.

{phang}
{opt noreport} suppresses the table of merge results.

{phang}
{opt verbose} displays progress and elapsed times.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt snapshot} uses {help preserve} to copy the whole master dataset for
recovery after an error or interruption. Stata may write the copy to disk
when it exceeds {cmd:c(max_preservemem)}.

{phang}
{opt fast} skips recovery to reduce memory use. An error or interruption can
leave a partially merged dataset. It cannot be combined with {opt snapshot}.

{marker remarks}{...}
{title:Remarks}

{pstd}
Shared numeric and fixed-width string variables are widened as needed to hold
both inputs; a {cmd:long}/{cmd:float} combination uses {cmd:double}. New variables
retain the using data's formats, labels, notes, and characteristics. If either
dataset contains a retained {cmd:strL} variable, the merge uses Stata's
{cmd:merge}; {opt preserve_order()} does not apply on this path.

{marker failure}{...}
{pstd}
By default, errors and interruptions restore the master data, except that a
failed {opt assert()} leaves the merged data in memory. Recovery uses an
in-memory copy; some merges use {help preserve} on Stata versions before 16.
{opt snapshot} selects {cmd:preserve} for every merge. {opt fast} disables
recovery.


{marker examples}{...}
{title:Examples}

{pstd}Create example datasets for merging:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. keep make price mpg}{p_end}
{phang2}{cmd:. save auto_master, replace}{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. keep make weight length foreign}{p_end}
{phang2}{cmd:. save auto_using, replace}{p_end}

{pstd}Simple one-to-one merge on make:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using}{p_end}

{pstd}Keep only matched observations:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, keep(match)}{p_end}

{pstd}Merge with verbose output and no merge variable:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, nogenerate verbose}{p_end}

{pstd}Keep only selected variables from using dataset:{p_end}
{phang2}{cmd:. use auto_master, clear}{p_end}
{phang2}{cmd:. cmerge 1:1 make using auto_using, keepusing(weight foreign)}{p_end}

{pstd}Many-to-one merge example with panel data:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. keep idcode year ln_wage}{p_end}
{phang2}{cmd:. save nlswork_wages, replace}{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. collapse (first) race, by(idcode)}{p_end}
{phang2}{cmd:. save nlswork_person, replace}{p_end}
{phang2}{cmd:. use nlswork_wages, clear}{p_end}
{phang2}{cmd:. cmerge m:1 idcode using nlswork_person}{p_end}


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
{synopt:{cmd:r(N_4)}}matches with missing values updated (with {opt update}){p_end}
{synopt:{cmd:r(N_5)}}matches with conflicting nonmissing values (with {opt update}){p_end}
{synopt:{cmd:r(time)}}elapsed time in seconds; missing on the native merge path{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(using)}}name of using file{p_end}
{synopt:{cmd:r(keyvars)}}key variable names{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Also see}

{psee}
Manual: {bf:[D] merge}

{psee}
Online: {help merge}, {help joinby}, {help ctools}, {help csort}
{p_end}
